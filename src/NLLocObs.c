/*
 * Copyright (C) 1999-2015 Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser Public License for more details.

 * You should have received a copy of the GNU Lesser Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */


/*   NLLocObs.c

        Observation-reading and arrival-processing for NonLinLoc: checkObs,
        GetObservations, GetNextObs (the multi-format phase reader), and the
        arrival helpers InitializeArrivalFields, isExcluded, lastLegType,
        CalcSimpleElevCorr, EvaluateArrivalAlias, ApplyTimeDelays,
        ApplySurfaceTimeDelay and ExtractFilenameInfo.

        Extracted verbatim from NLLocLib.c during the modernization effort
        (Phase 3: modularization). Pure code movement, no logic change; the
        public prototypes remain in NLLocLib.h, so callers are unaffected.
 */


/* NLLocLib.h is not self-contained; these must be included before it,
   in the same order as NLLocLib.c */
#include "GridLib.h"
#include "ran1/ran1.h"
#include "velmod.h"
#include "GridMemLib.h"
#include "calc_crust_corr.h"
#include "phaseloclist.h"
#include "otime_limit.h"
#include "NLLocLib.h"
#include "json_io.h"

#ifdef CUSTOM_ETH
#include "custom_eth/eth_functions.h"
#endif


/** function to check and initialize an observation */

int checkObs(ArrivalDesc *arrival, int nobs) {

    int istat;

    // 20100506 AJL - added to support preservation of observation index order for calls to NLLoc() function (e.g. from SeisComp3)
    // 20130326 AJL - moved to calling function
    //arrival[nobs].original_obs_index = nobs;  // sets index of observation in order originally read from input

    /* check for aliased arrival label */
    if ((istat = EvaluateArrivalAlias(arrival + nobs)) < 0)
        ;
    /* set some fields */
    InitializeArrivalFields(arrival + nobs);
    /* check some fields */
    if (!isgraph(arrival[nobs].phase[0]))
        snprintf(arrival[nobs].phase, sizeof(arrival[nobs].phase), "%s", ARRIVAL_NULL_STR);
    if (!isgraph(arrival[nobs].comp[0]))
        snprintf(arrival[nobs].comp, sizeof(arrival[nobs].comp), "%s", ARRIVAL_NULL_STR);
    if (!isgraph(arrival[nobs].onset[0]))
        snprintf(arrival[nobs].onset, sizeof(arrival[nobs].onset), "%s", ARRIVAL_NULL_STR);
    if (!isgraph(arrival[nobs].first_mot[0]))
        snprintf(arrival[nobs].first_mot, sizeof(arrival[nobs].first_mot), "%s", ARRIVAL_NULL_STR);
    if (arrival[nobs].coda_dur < VERY_SMALL_DOUBLE)
        arrival[nobs].coda_dur = CODA_DUR_NULL;
    if (arrival[nobs].amplitude < VERY_SMALL_DOUBLE)
        arrival[nobs].amplitude = AMPLITUDE_NULL;
    if (arrival[nobs].period < VERY_SMALL_DOUBLE)
        arrival[nobs].period = PERIOD_NULL;
    if (message_flag >= 3) {
        // display arrival parameters
        snprintf(MsgStr, sizeof(MsgStr),
                "Arrival %d:  %s (%s)  %s %s %s %d  %4.4d %2.2d %2.2d   %2.2d %2.2d %lf  Unc: %s %lf  Amp: %lf  Dur: %lf  Per: %lf",
                nobs,
                arrival[nobs].label,
                arrival[nobs].time_grid_label,
                arrival[nobs].onset,
                arrival[nobs].phase,
                arrival[nobs].first_mot,
                arrival[nobs].quality,
                arrival[nobs].year,
                arrival[nobs].month,
                arrival[nobs].day,
                arrival[nobs].hour,
                arrival[nobs].min,
                arrival[nobs].sec,
                arrival[nobs].error_type,
                arrival[nobs].error,
                arrival[nobs].amplitude,
                arrival[nobs].coda_dur,
                arrival[nobs].period);
        nll_putmsg(3, MsgStr);
    }
    /* remove blanks/whitespace from phase string */
    removeSpace(arrival[nobs].phase);
    /* check for reject phase code */
    if (IsPhaseID(arrival[nobs].phase, "$")) {
        snprintf(MsgStr, sizeof(MsgStr),
                "WARNING: phase code is $, rejecting observation: %s %s", arrival[nobs].label, arrival[nobs].phase);
        nll_putmsg(2, MsgStr);
        return (-1);
    }
    /* check for valid P or S phase code */
    /*INGV		if (!IsPhaseID(arrival[nobs].phase, "P")
            && !IsPhaseID(arrival[nobs].phase, "S")) {
            snprintf(MsgStr, sizeof(MsgStr),
            "WARNING: phase code not in P or S phase id list, rejecting observation: %s %s", arrival[nobs].label, arrival[nobs].phase);
            nll_putmsg(2, MsgStr);
            return(-1);
    }
            INGV*/
    /* check for duplicate arrival */
    if (nll_mode != MODE_DIFFERENTIAL) {
        // AJL 20041201 - algorithm changed, less strict
        //if ( IsDuplicateArrival(arrival, nobs + 1, nobs, NULL) >= 0
        //		|| IsDuplicateArrival(arrival, nobs + 1, nobs, arrival[nobs].phase) >= 0 )
        //		) {
        // AJL 20200131 - iRejectDuplicateArrivals > 1 flags no check, accept all arrivals
        //if (IsDuplicateArrival(arrival, nobs + 1, nobs, !iRejectDuplicateArrivals) >= 0) {
        if (iRejectDuplicateArrivals > -2 && IsDuplicateArrival(arrival, nobs + 1, nobs, !iRejectDuplicateArrivals) >= 0) {

            snprintf(MsgStr, sizeof(MsgStr),
                    "WARNING: duplicate arrival, rejecting observation: %s %s", arrival[nobs].label, arrival[nobs].phase);
            nll_putmsg(2, MsgStr);

            return (-1);
        }
    }

    // all OK
    return (1);

}

/** function read observation file and open station grid files */

int GetObservations(FILE* fp_obs, char* ftype_obs, char fn_ttgrids[MAX_NUM_TIME_GRID_PATHS][FILENAME_MAX],
        ArrivalDesc *arrival, int *pi_end_of_input,
        int* pnignore, int* pnreject, int maxNumArrivals, HypoDesc* phypo,
        int* pMaxArrExceeded, int *pnumSArrivals, int nobs_prev) {

    int nobs, nobs_read, nobs_total;
    int index_obs_original; // 20130326 AJL - added
    int istat, ntry, nLocate, n_compan, n_time_grid;
    char filename[FILENAME_MAX];
    char eval_phase[PHASE_LABEL_LEN];
    int isZeroWeight;
    char arrival_phase[PHASE_LABEL_LEN];
    char tmp_phase[PHASE_LABEL_LEN];
    int read_2d_sheets;
    int i_need_elev_corr, n_phs_try;

    SourceDesc* pstation;

    *pnumSArrivals = 0;

    /* initalize format specific event data */
    HypoInverseArchiveSumHdr[0] = '\0';


    /** read observations to arrival array */

    nobs = nobs_prev;
    index_obs_original = 0; // 20130326 AJL - added
    *pnignore = 0;
    *pnreject = 0;
    nLocate = 0;
    ntry = 0;

    nll_putmsg(4, "Dummy message");

    // 20180907 AJL - added phypo to recover location information (e.g. magntiude) from observation file
    // 20180907 AJL while ((istat = GetNextObs(fp_obs, arrival + nobs, ftype_obs, ntry++ == 0)) != EOF) {
    while ((istat = GetNextObs(phypo, fp_obs, arrival + nobs, ftype_obs, ntry++ == 0)) != EOF) {


        if (istat == OBS_FILE_INTERNAL_ERROR) {
            nll_putmsg(1, "ERROR: internal error reading observation file.");
            break;
        }
        if (istat == OBS_FILE_FORMAT_ERROR) {
            nll_putmsg(1, "ERROR: format error in observation file.");
            break;
        }
        if (istat == OBS_FILE_END_OF_EVENT) {
            if (nobs == 0) {
                ntry = 0; // 20160925 AJL - bug fix, to reset event fields as if new event
                continue;
            } else {
                break;
            }
        }
        if (istat == OBS_FILE_END_OF_INPUT) {
            *pi_end_of_input = 1;
            break;
        }
        if (istat == OBS_FILE_INVALID_PHASE) {
            index_obs_original++; // 20130326 AJL - added
            continue;
        }
        if (istat == OBS_FILE_INVALID_DATE) {
            index_obs_original++; // 20130326 AJL - added
            continue;
        }
        if (istat == OBS_FILE_SKIP_INPUT_LINE)
            continue;
        /* SH if comment line found -> read next obs */
        if (istat == OBS_IS_COMMENT_LINE)
            continue;
        if (nobs == MAX_NUM_ARRIVALS - 1) {
            if (!*pMaxArrExceeded) {
                *pMaxArrExceeded = 1;
                snprintf(MsgStr, sizeof(MsgStr), "WARNING: maximum number of arrivals exceeded, only first %d will be processed.", MAX_NUM_ARRIVALS);
                nll_putmsg(1, MsgStr);
            }
            continue;
        }

        // 20100506 AJL - added to support preservation of observation index order for calls to NLLoc() function (e.g. from SeisComp3)
        arrival[nobs].original_obs_index = index_obs_original; // sets index of observation in order originally read from input    // 20130326 AJL - added
        index_obs_original++; // 20130326 AJL - added
        // check and initialize observation
        if (checkObs(arrival, nobs) > 0)
            nobs++;
        // check and initialize second observation
        // SH if istat = OBS_FILE_TWO_ARRIVALS_READ then S-arrival has been read in -> two observations per line
        if (istat == OBS_FILE_TWO_ARRIVALS_READ) {
            // 20100506 AJL - added to support preservation of observation index order for calls to NLLoc() function (e.g. from SeisComp3)
            arrival[nobs].original_obs_index = index_obs_original; // sets index of observation in order originally read from input    // 20130326 AJL - added
            index_obs_original++; // 20130326 AJL - added
            if (checkObs(arrival, nobs + 1) > 0)
                nobs++;
        }

    }


    // save number of obs read
    nobs_total = nobs;
    nobs_read = nobs - nobs_prev;


    /** check for minimum number of arrivals */

    //printf("*pi_end_of_input %d  nobs_read %d  MinNumArrLoc %d\n", *pi_end_of_input, nobs_read, MinNumArrLoc);

    if (*pi_end_of_input && nobs_total < 1) {
        return (nLocate + *pnignore);
    }



    /** process arrivals, determine if reject or ignore */

    // check for no absolute timing (inst begins with '*')
    CheckAbsoluteTiming(arrival + nobs_prev, nobs_read);

    /* homogenize date/time */
    if ((istat = HomogDateTime(arrival, nobs_total, phypo)) < 0) {
        nll_puterr("ERROR: in arrival date/times.");
        if (istat == OBS_FILE_ARRIVALS_CROSS_YEAR_BOUNDARY)
            nll_puterr("ERROR: arrivals cross year boundary.");
        return (-1);
    }


    /* sort to get arrivals in time order */
    if (nll_mode != MODE_DIFFERENTIAL && (istat = SortArrivalsTime(arrival, nobs_total)) < 0) {
        nll_puterr("ERROR: sorting arrivals by time.");
        return (-1);
    }


    /* check each arrival and initialize */

    Num3DGridReadToMemory = 0;
    for (nobs = nobs_prev; nobs < nobs_total; nobs++) {

        if (message_flag >= 3) {
            snprintf(MsgStr, sizeof(MsgStr), "Checking Arrival %d:  %s (%s)  %s %s %s %d",
                    nobs,
                    arrival[nobs].label,
                    arrival[nobs].time_grid_label,
                    arrival[nobs].onset,
                    arrival[nobs].phase,
                    arrival[nobs].first_mot,
                    arrival[nobs].quality);
            nll_putmsg(3, MsgStr);
        } else if (nobs_read > 1000 && message_flag > 0 && nobs % 100 == 0) {
            fprintf(stdout, "Checking Arrival %d/%d\r", nobs, nobs_read);
            fflush(stdout);
        }

        // make sure station location is null
        arrival[nobs].station.x = arrival[nobs].station.y = arrival[nobs].station.z = -LARGE_DOUBLE;

        // set some flags
        read_2d_sheets = 1;
        i_need_elev_corr = 0;

        snprintf(arrival_phase, sizeof(arrival_phase), "%s", arrival[nobs].phase);

        // set phase groups
        arrival[nobs].isP = IsPhaseID(arrival_phase, "P");
        arrival[nobs].isS = IsPhaseID(arrival_phase, "S");

        // check for zero weight phase
        isZeroWeight = 0;
        if (arrival_phase[0] == '*') {
            snprintf(MsgStr, sizeof(MsgStr),
                    "INFO: arrival is isZeroWeight since first character of phase is *: %s %s",
                    arrival[nobs].label, arrival[nobs].phase);
            nll_putmsg(3, MsgStr);
            isZeroWeight = 1;
            snprintf(tmp_phase, sizeof(tmp_phase), "%s", arrival_phase + 1);
            snprintf(arrival_phase, sizeof(arrival_phase), "%s", tmp_phase);
        } else if (arrival[nobs].apriori_weight < VERY_SMALL_DOUBLE) {
            snprintf(MsgStr, sizeof(MsgStr),
                    "INFO: arrival is isZeroWeight since apriori_weight < VERY_SMALL_DOUBLE: %s %s",
                    arrival[nobs].label, arrival[nobs].phase);
            nll_putmsg(3, MsgStr);
            isZeroWeight = 1;
            /*printf("isZeroWeight Arrival %d:  %s (%s)  %s %s %s %d %lf",
            nobs,
            arrival[nobs].label,
            arrival[nobs].time_grid_label,
            arrival[nobs].onset,
            arrival[nobs].phase,
            arrival[nobs].first_mot,
            arrival[nobs].quality,
            arrival[nobs].apriori_weight);*/
        } else if (arrival[nobs].error >= ARRIVAL_ERROR_NULL_TEST) {
            snprintf(MsgStr, sizeof(MsgStr),
                    "INFO: arrival is isZeroWeight since arrival[nobs].error >= ARRIVAL_ERROR_NULL_TEST: %s %s",
                    arrival[nobs].label, arrival[nobs].phase);
            nll_putmsg(3, MsgStr);
            isZeroWeight = 1;
        }

        snprintf(arrival[nobs].fileroot, sizeof(arrival[nobs].fileroot), "%s", "\0");
        arrival[nobs].n_companion = -1;
        arrival[nobs].n_time_grid = -1;
        arrival[nobs].tfact = 1.0;

        /* if Vp/Vs > 0, check for companion phase, its time grid will be used for times */
        if (VpVsRatio > 0.0) {
            /* try finding previously initialized companion phase */
            if (IsPhaseID(arrival_phase, "S") &&
                    (n_compan =
                    IsSameArrival(arrival, nobs, nobs, "P")) >= 0 &&
                    arrival[n_compan].flag_ignore == 0) {
                arrival[nobs].tfact = VpVsRatio;
                arrival[nobs].gdesc.type = arrival[n_compan].gdesc.type;
                arrival[nobs].station = arrival[n_compan].station;
                arrival[nobs].n_companion = n_compan;
                if (message_flag >= 3) {
                    snprintf(MsgStr, sizeof(MsgStr),
                            "INFO: S phase: %d %s %s using companion phase %d travel time grids.",
                            nobs, arrival[nobs].label, arrival[nobs].phase, arrival[nobs].n_companion);
                    nll_putmsg(3, MsgStr);
                }
                // save filename as grid identifier (needed for GridMemList)
                snprintf(arrival[nobs].gdesc.title, sizeof(arrival[nobs].gdesc.title), "%s", arrival[n_compan].gdesc.title);
            }
        }

        /* if MODE_DIFFERENTIAL, check for same station phase, its time grid will be used for times */
        if (nll_mode == MODE_DIFFERENTIAL && arrival[nobs].n_companion < 0) {
            /* try finding previously initialized companion phase */
            if ((n_compan = IsSameArrival(arrival, nobs, nobs, NULL)) >= 0 &&
                    arrival[n_compan].flag_ignore == 0) {
                arrival[nobs].gdesc.type = arrival[n_compan].gdesc.type;
                arrival[nobs].station = arrival[n_compan].station;
                arrival[nobs].n_companion = n_compan;
                if (message_flag >= 3) {
                    snprintf(MsgStr, sizeof(MsgStr),
                            "INFO: MODE_DIFFERENTIAL: %d %s %s using companion phase %d travel time grids.",
                            nobs, arrival[nobs].label, arrival[nobs].phase, arrival[nobs].n_companion);
                    nll_putmsg(3, MsgStr);
                }
                // save filename as grid identifier (needed for GridMemList)
                snprintf(arrival[nobs].gdesc.title, sizeof(arrival[nobs].gdesc.title), "%s", arrival[n_compan].gdesc.title);
            }
        }

        /* no previously initialized companion phase */
        if (arrival[nobs].n_companion < 0) {

            for (int n_ttgrid = 0; n_ttgrid < NumTimeGridPaths; n_ttgrid++) { // 20251027 add support for alternative travel-time grid path/root

                // DEBUG
                /*snprintf(MsgStr, sizeof(MsgStr), "INFO: n_ttgrid: %d, NumGridBufFilesOpen: %d, NumGridHdrFilesOpen: %d, NumFilesOpen: %d",
                        n_ttgrid, NumGridBufFilesOpen, NumGridHdrFilesOpen, NumFilesOpen);
                nll_putmsg(0, MsgStr);*/

                // try to open time grid file using original phase ID
                snprintf(arrival[nobs].fileroot, sizeof(arrival[nobs].fileroot), "%s.%s.%s", fn_ttgrids[n_ttgrid],
                        arrival_phase, arrival[nobs].time_grid_label);
                snprintf(filename, sizeof(filename), "%s.time", arrival[nobs].fileroot);
                // try opening time grid file for this phase
                istat = OpenGrid3dFile(filename,
                        &(arrival[nobs].fpgrid),
                        &(arrival[nobs].fphdr),
                        &(arrival[nobs].gdesc), "time",
                        &(arrival[nobs].station),
                        arrival[nobs].gdesc.iSwapBytes);

                if (istat < 0) {
                    // try to open time grid file using LOCPHASEID mapped phase ID
                    EvalPhaseID(eval_phase, sizeof(eval_phase), arrival_phase);
                    snprintf(arrival[nobs].fileroot, sizeof(arrival[nobs].fileroot), "%s.%s.%s", fn_ttgrids[n_ttgrid],
                            eval_phase, arrival[nobs].time_grid_label);
                    snprintf(filename, sizeof(filename), "%s.time", arrival[nobs].fileroot);
                    /* try opening time grid file for this phase */
                    istat = OpenGrid3dFile(filename,
                            &(arrival[nobs].fpgrid),
                            &(arrival[nobs].fphdr),
                            &(arrival[nobs].gdesc), "time",
                            &(arrival[nobs].station),
                            arrival[nobs].gdesc.iSwapBytes);
                }

                /* try opening P time grid file for S if no P companion phase */
                if (istat < 0 && VpVsRatio > 0.0 && IsPhaseID(arrival_phase, "S")) {
                    arrival[nobs].tfact = VpVsRatio;
                    snprintf(arrival[nobs].fileroot, sizeof(arrival[nobs].fileroot), "%s.%s.%s", fn_ttgrids[n_ttgrid],
                            "P", arrival[nobs].time_grid_label);
                    snprintf(filename, sizeof(filename), "%s.time", arrival[nobs].fileroot);
                    istat = OpenGrid3dFile(filename,
                            &(arrival[nobs].fpgrid),
                            &(arrival[nobs].fphdr),
                            &(arrival[nobs].gdesc), "time",
                            &(arrival[nobs].station),
                            arrival[nobs].gdesc.iSwapBytes);
                    if (message_flag >= 3) {
                        snprintf(MsgStr, sizeof(MsgStr),
                                "INFO: S phase: using P phase travel time grid file: %s", filename);
                        nll_putmsg(3, MsgStr);
                    }
                }

                // check if station/source in grid hdr file was DEFAULT
                if (istat >= 0 && strcmp(arrival[nobs].station.label, "DEFAULT") == 0) {
                    // get station/source coordinates, etc
                    pstation = FindSource(arrival[nobs].time_grid_label);
                    if (pstation != NULL) {
                        arrival[nobs].station = *pstation;
                        i_need_elev_corr = 1;
                    } else
                        istat = -1;
                }

                // try opening DEFAULT time grid file
                if (istat < 0) {
                    pstation = FindSource(arrival[nobs].time_grid_label);
                    if (pstation != NULL) {
                        // open DEFAULT time grid for this phase
                        n_phs_try = 0;
                        while (istat < 0 && n_phs_try < 2) {
                            n_phs_try++;
                            if (n_phs_try == 1) // try to open time grid file using original phase ID
                                snprintf(eval_phase, sizeof(eval_phase), "%s", arrival_phase);
                            else // try to open time grid file using LOCPHASEID mapped phase ID
                                EvalPhaseID(eval_phase, sizeof(eval_phase), arrival_phase);
                            snprintf(arrival[nobs].fileroot, sizeof(arrival[nobs].fileroot), "%s.%s.%s", fn_ttgrids[n_ttgrid], eval_phase, "DEFAULT");
                            snprintf(filename, sizeof(filename), "%s.time", arrival[nobs].fileroot);

                            //#define LOC2SSST_CLUGE
#ifdef LOC2SSST_CLUGE
                            // 20201022 AJL - Cluge, assume need to byte swap if DEFAULT  // TODO: make this automatic or configurable
                            arrival[nobs].gdesc.iSwapBytes = 1;
                            //
#endif

                            /* check if time grid already read */
                            if ((n_time_grid = FindDuplicateTimeGrid(arrival, nobs, nobs)) >= 0 && arrival[n_time_grid].flag_ignore == 0) {
                                arrival[nobs].gdesc.type = arrival[n_time_grid].gdesc.type;
                                arrival[nobs].sheetdesc = arrival[n_time_grid].sheetdesc;
                                arrival[nobs].station = *pstation;
                                arrival[nobs].n_time_grid = n_time_grid;
                                if (message_flag >= 3) {
                                    snprintf(MsgStr, sizeof(MsgStr),
                                            "INFO: DEFAULT travel time: %d %s %s using previous phase %d travel time grids.",
                                            nobs, arrival[nobs].label, arrival[nobs].phase, arrival[nobs].n_time_grid);
                                    nll_putmsg(3, MsgStr);
                                }
                                // save filename as grid identifier (needed for GridMemList)
                                snprintf(arrival[nobs].gdesc.title, sizeof(arrival[nobs].gdesc.title), "%s", arrival[n_time_grid].gdesc.title);
                                istat = 1;
                            } else {
                                istat = OpenGrid3dFile(filename,
                                        &(arrival[nobs].fpgrid),
                                        &(arrival[nobs].fphdr),
                                        &(arrival[nobs].gdesc), "time",
                                        &(arrival[nobs].station),
                                        arrival[nobs].gdesc.iSwapBytes);
                                if (istat >= 0 && message_flag >= 3) {
                                    snprintf(MsgStr, sizeof(MsgStr),
                                            "INFO: using DEFAULT travel time grid file: %s", filename);
                                    nll_putmsg(3, MsgStr);
                                }
                                arrival[nobs].station = *pstation;
                            }
                        }
                        i_need_elev_corr = 1;
                        //read_2d_sheets = 0; // too slow!
                    }
                }

                // save filename as grid identifier (needed for GridMemList)
                snprintf(arrival[nobs].gdesc.title, sizeof(arrival[nobs].gdesc.title), "%s", filename);

                if (istat < 0) {
                    CloseGrid3dFile(&(Arrival[nobs].gdesc), &(Arrival[nobs].fpgrid), &(arrival[nobs].fphdr));
                    snprintf(arrival[nobs].fileroot, sizeof(arrival[nobs].fileroot), "%s", "\0");
                    if (n_ttgrid < NumTimeGridPaths - 1) { // 20251027 add support for alternative travel-time grid path/root
                        continue;
                    }
                    snprintf(MsgStr, sizeof(MsgStr),
                            "WARNING: cannot open time grid file: %s: rejecting grid for: %s %s",
                            filename, arrival[nobs].label, arrival[nobs].phase);
                    nll_putmsg(2, MsgStr);
                    goto RejectArrival;
                }


                /* check that search grid is inside time grid (3D grids) */

                if (arrival[nobs].gdesc.type == GRID_TIME &&
                        !IsGridInside(LocGrid + 0, &(arrival[nobs].gdesc), 0)) {
                    CloseGrid3dFile(&(Arrival[nobs].gdesc), &(Arrival[nobs].fpgrid), &(arrival[nobs].fphdr));
                    if (n_ttgrid < NumTimeGridPaths - 1) { // 20251027 add support for alternative travel-time grid path/root
                        continue;
                    }
                    snprintf(MsgStr, sizeof(MsgStr),
                            "WARNING: initial location search grid not contained inside arrival time grid, rejecting grid for: %s %s", arrival[nobs].label, arrival[nobs].phase);
                    nll_putmsg(1, MsgStr);
                    goto RejectArrival;
                }

                /* (3D grids) check that
                                         (1) distance from center of search grid to station
                                         is greater than DistStaGridMin (if DistStaGridMin > 0)
                                         is less than DistStaGridMax (if DistStaGridMax > 0) */
                // 20190522 AJL - added to support LOCMETH maxDistStaGrid and minDistStaGrid with 3D grids

                if (arrival[nobs].gdesc.type == GRID_TIME &&
                        (istat = IsDistStaGridOK(LocGrid + 0,
                        &(arrival[nobs].station),
                        DistStaGridMin, DistStaGridMax,
                        LocGrid[0].origx + (LocGrid[0].dx
                        * (double) (LocGrid[0].numx - 1)) / 2.0,
                        LocGrid[0].origy + (LocGrid[0].dy
                        * (double) (LocGrid[0].numy - 1)) / 2.0)
                        ) != 1) {
                    CloseGrid3dFile(&(Arrival[nobs].gdesc), &(Arrival[nobs].fpgrid), &(arrival[nobs].fphdr));
                    if (istat == -2) {
                        if (n_ttgrid < NumTimeGridPaths - 1) { // 20251027 add support for alternative travel-time grid path/root
                            continue;
                        }
                        snprintf(MsgStr, sizeof(MsgStr),
                                "WARNING: distance from grid center to station \n\texceeds maximum station distance, ignoring observation in misfit calculation: %s %s",
                                arrival[nobs].label, arrival[nobs].phase);
                        nll_putmsg(2, MsgStr);
                        arrival[nobs].flag_ignore = 1;
                        goto IgnoreArrival;
                    }
                }


                /* (2D grids) check that
                                         (1) greatest distance from search grid to station
                                         is inside time grid, and
                                         (2) distance from center of search grid to station
                                         is greater than DistStaGridMin (if DistStaGridMin > 0)
                                         is less than DistStaGridMax (if DistStaGridMax > 0) */

                if (arrival[nobs].gdesc.type == GRID_TIME_2D &&
                        (istat = IsGrid2DBigEnough(LocGrid + 0,
                        &(arrival[nobs].gdesc),
                        &(arrival[nobs].station),
                        DistStaGridMin, DistStaGridMax,
                        LocGrid[0].origx + (LocGrid[0].dx
                        * (double) (LocGrid[0].numx - 1)) / 2.0,
                        LocGrid[0].origy + (LocGrid[0].dy
                        * (double) (LocGrid[0].numy - 1)) / 2.0)
                        ) != 1) {
                    CloseGrid3dFile(&(Arrival[nobs].gdesc), &(Arrival[nobs].fpgrid), &(arrival[nobs].fphdr));
                    if (istat == -1) {
                        if (n_ttgrid < NumTimeGridPaths - 1) { // 20251027 add support for alternative travel-time grid path/root
                            continue;
                        }
                        snprintf(MsgStr, sizeof(MsgStr),
                                "WARNING: greatest distance from initial 3D location search grid to station \n\texceeds 2D time grid size, rejecting grid for: %s %s",
                                arrival[nobs].label, arrival[nobs].phase);
                        nll_putmsg(2, MsgStr);
                        goto RejectArrival;
                    } else if (istat == -2) {
                        if (n_ttgrid < NumTimeGridPaths - 1) { // 20251027 add support for alternative travel-time grid path/root
                            continue;
                        }
                        snprintf(MsgStr, sizeof(MsgStr),
                                "WARNING: distance from grid center to station \n\texceeds maximum station distance, ignoring observation in misfit calculation: %s %s",
                                arrival[nobs].label, arrival[nobs].phase);
                        nll_putmsg(2, MsgStr);
                        arrival[nobs].flag_ignore = 1;
                        goto IgnoreArrival;
                    }
                    if (istat == -3) {
                        if (n_ttgrid < NumTimeGridPaths - 1) { // 20251027 add support for alternative travel-time grid path/root
                            continue;
                        }
                        snprintf(MsgStr, sizeof(MsgStr),
                                "WARNING: depth range of initial 3D location search grid exceeds that of 2D time grid size, rejecting grid for: %s %s",
                                arrival[nobs].label, arrival[nobs].phase);
                        nll_putmsg(2, MsgStr);
                        goto RejectArrival;
                    }
                }

                // 20251027 add support for alternative travel-time grid path/root
                // made it this far, so have valid travel-time grid
                break;

            } // for (int n_ttgrid = 0; n_ttgrid < NumTimeGridPaths; n_ttgrid++) {
        } // if (arrival[nobs].n_companion < 0)


        /* check for time delays */
        ApplyTimeDelays(arrival + nobs);
        ;
        // calculate elevation correction if needed
        if (ApplyElevCorrFlag && i_need_elev_corr) {
            arrival[nobs].elev_corr = CalcSimpleElevCorr(arrival, nobs, ElevCorrVelP, ElevCorrVelS);
        }



        /* check if arrival is excluded */

        // zero weight or excluded in control file
        if (isZeroWeight || isExcluded(arrival[nobs].label, arrival[nobs].phase)) {
            snprintf(MsgStr, sizeof(MsgStr),
                    "INFO: arrival is excluded, ignoring observation in misfit calculation: %s %s",
                    arrival[nobs].label, arrival[nobs].phase);
            nll_putmsg(2, MsgStr);
            arrival[nobs].flag_ignore = 1;
            CloseGrid3dFile(&(Arrival[nobs].gdesc), &(Arrival[nobs].fpgrid), &(arrival[nobs].fphdr));
            goto IgnoreArrival;
        }
        // no absolute time and not EDT
        if (!arrival[nobs].abs_time && LocMethod != METH_EDT) {
            snprintf(MsgStr, sizeof(MsgStr),
                    "INFO: arrival does not have absolute timing, ignoring observation in misfit calculation: %s %s %s",
                    arrival[nobs].label, arrival[nobs].phase, arrival[nobs].inst);
            arrival[nobs].flag_ignore = 1;
            CloseGrid3dFile(&(Arrival[nobs].gdesc), &(Arrival[nobs].fpgrid), &(arrival[nobs].fphdr));
            nll_putmsg(2, MsgStr);
            goto IgnoreArrival;
        }
        // EDT_BOX but not BOX error
        if (!strcmp(arrival[nobs].error_type, "BOX") && LocMethod != METH_EDT_BOX) {
            snprintf(MsgStr, sizeof(MsgStr),
                    "INFO: method is EDT_BOX but arrival does not have error type BOX, ignoring observation in misfit calculation: %s %s %s",
                    arrival[nobs].label, arrival[nobs].phase, arrival[nobs].inst);
            arrival[nobs].flag_ignore = 1;
            CloseGrid3dFile(&(Arrival[nobs].gdesc), &(Arrival[nobs].fpgrid), &(arrival[nobs].fphdr));
            nll_putmsg(2, MsgStr);
            goto IgnoreArrival;
        }


        /* check if maximum number of arrivals for location exceeded */

        if (nLocate + nobs_prev >= maxNumArrivals) {
            nll_putmsg(2,
                    "WARNING: maximum number of arrivals for location exceeded, \n\tignoring observation in misfit calculation.");
            arrival[nobs].flag_ignore = 1;
            CloseGrid3dFile(&(Arrival[nobs].gdesc), &(Arrival[nobs].fpgrid), &(arrival[nobs].fphdr));
            goto IgnoreArrival;
        }



        /** arrival to be used for location */

        /* no further processing for arrivals with companion time grids */
        if (arrival[nobs].n_companion >= 0 || arrival[nobs].n_time_grid >= 0)
            goto AcceptArrival;

        /** prepare time grids access in memory or on disk */

        /* construct dual-sheet description
        (2D grids for all search types,
        and 3D grids for grid-search) */

        if (SearchType == SEARCH_GRID
                || (read_2d_sheets && arrival[nobs].gdesc.type == GRID_TIME_2D)) {

            arrival[nobs].sheetdesc = arrival[nobs].gdesc;
            //INGV ??
            //if (arrival[nobs].gdesc.numx > 1)
            arrival[nobs].sheetdesc.numx = 2;
            /* allocate grid */
            arrival[nobs].sheetdesc.buffer =
                    AllocateGrid(&(arrival[nobs].sheetdesc));
            if (arrival[nobs].sheetdesc.buffer == NULL) {
                nll_puterr(
                        "ERROR: allocating memory for arrival sheet buffer.");
                goto RejectArrival;
                //return(EXIT_ERROR_MEMORY);
            }
            /* create array access pointers */
            arrival[nobs].sheetdesc.array =
                    CreateGridArray(&(arrival[nobs].sheetdesc));
            if (arrival[nobs].sheetdesc.array == NULL) {
                nll_puterr(
                        "ERROR: creating array for accessing arrival sheet buffer.");
                goto RejectArrival;
                //return(EXIT_ERROR_MEMORY);
            }
            arrival[nobs].sheetdesc.origx = VERY_LARGE_DOUBLE;

        }


        /* read 3D grid into memory (3D grids for Metropolis or Octtree search) */

        //int XX_last = NumAllocations;
        if ((SearchType == SEARCH_MET || SearchType == SEARCH_OCTTREE)
                && arrival[nobs].gdesc.type == GRID_TIME
                && (MaxNum3DGridMemory < 0 || Num3DGridReadToMemory < MaxNum3DGridMemory)) {

            /* allocate grid */
            arrival[nobs].gdesc.buffer = NLL_AllocateGrid(&(arrival[nobs].gdesc));
            //printf("DEBUG: ALLOCATE: nobs %d Arrival[nobs].n_time_grid %d\n", nobs, Arrival[nobs].n_time_grid);
            if (arrival[nobs].gdesc.buffer == NULL) {
                nll_puterr(
                        "ERROR: allocating memory for arrival time grid buffer; grid will be read from disk.");
            } else {
                /* create array access pointers */
                arrival[nobs].gdesc.array = NLL_CreateGridArray(&(arrival[nobs].gdesc));
                if (arrival[nobs].gdesc.array == NULL) {
                    nll_puterr(
                            "ERROR: creating array for accessing arrival time grid buffer.");
                    goto RejectArrival;
                    //return(EXIT_ERROR_MEMORY);
                }
                /* read time grid */
                if (NLL_ReadGrid3dBuf(&(arrival[nobs].gdesc), arrival[nobs].fpgrid) < 0) {
                    nll_puterr(
                            "ERROR: reading arrival time grid buffer.");
                    goto RejectArrival;
                    //return(EXIT_ERROR_MEMORY);
                }
                CloseGrid3dFile(&(Arrival[nobs].gdesc), &(Arrival[nobs].fpgrid), &(arrival[nobs].fphdr));
                Num3DGridReadToMemory++;
            }
        }
        //printf("XXX: NLLoc try put in memory: NumAllocations %d->%d\n", XX_last, NumAllocations);


        /* read time grid and close file (2D grids)*/

        if (read_2d_sheets && arrival[nobs].gdesc.type == GRID_TIME_2D) {
            istat = ReadArrivalSheets(1, &(arrival[nobs]), 0.0);
            CloseGrid3dFile(&(Arrival[nobs].gdesc), &(Arrival[nobs].fpgrid), &(arrival[nobs].fphdr));
            if (istat < 0) {
                snprintf(MsgStr, sizeof(MsgStr),
                        "ERROR: reading arrival travel time sheets (2D grid), rejecting observation: %s %s",
                        arrival[nobs].label, arrival[nobs].phase);
                nll_puterr(MsgStr);
                goto RejectArrival;
            }
        }

        /* arrival accepted, use for location */
AcceptArrival:
        arrival[nobs].flag_ignore = 0;
        nLocate++;
        if (IsPhaseID(arrival[nobs].phase, "S"))
            (*pnumSArrivals)++;


        /* arrival accepted, ignore for location */
IgnoreArrival:

        *pnignore += arrival[nobs].flag_ignore;

        continue;


        /* arrival rejected */
RejectArrival:

        (*pnreject)++;
        arrival[nobs].flag_ignore = 999;
        if (message_flag >= 3) {
            snprintf(MsgStr, sizeof(MsgStr), "   Rejected Arrival %d:  %s (%s)  %s %s %s %d",
                    nobs,
                    arrival[nobs].label,
                    arrival[nobs].time_grid_label,
                    arrival[nobs].onset,
                    arrival[nobs].phase,
                    arrival[nobs].first_mot,
                    arrival[nobs].quality);
            nll_putmsg(3, MsgStr);
        }

    }


    /* avoid returning 0 if arrivals were read, return 0 indicates end of file */
    if (nLocate + *pnignore == 0 && nobs_read > 0) {
        snprintf(MsgStr, sizeof(MsgStr),
                "WARNING: %d arrivals read, but none accepted for location.", nobs_read);
        nll_putmsg(2, MsgStr);

        return (-1);
    }

    // AJL 20091208 Zero weight phase modification
    // following block is present in NLLoc1.c
    //if (!(*pi_end_of_input) && (nobs_total > 0 && nobs_total < MinNumArrLoc)) {
    /*if (!(*pi_end_of_input) && (nobs_total > 0 && nLocate < MinNumArrLoc)) {
        snprintf(MsgStr, sizeof(MsgStr),
                "WARNING: too few observations to locate (%d available, %d needed), skipping event.", nobs_total, MinNumArrLoc);
        nll_putmsg(1, MsgStr);
        if (message_flag >= 3) {
            snprintf(MsgStr, sizeof(MsgStr),
                    "INFO: %d observations needed (specified in control file entry LOCMETH).",
                    MinNumArrLoc);
            nll_putmsg(3, MsgStr);
        }
        return (-1);
    }*/

    return (nLocate + *pnignore);
}

/** function to initialize arrival fields */

void InitializeArrivalFields(ArrivalDesc * arrival) {

    arrival->abs_time = 1;
    arrival->obs_centered = 0.0;
    arrival->pred_travel_time = 0.0;
    arrival->pred_travel_time_best = -LARGE_DOUBLE;
    arrival->pred_centered = 0.0;
    arrival->cent_resid = 0.0;
    arrival->obs_travel_time = 0.0;
    // 	arrival->delay = 0.0;	// do not initialize here!
    arrival->elev_corr = 0.0;
    arrival->residual = 0.0;
    arrival->dist = 0.0;
    arrival->azim = 0.0;
    // 20110601 AJL - default/unset values set to -1.0 so unset angles can be more easily identified in later processing (e.g. in SeisComp3).
    //arrival->ray_azim = 0.0;
    //arrival->ray_dip = 0.0;
    arrival->ray_azim = -1.0;
    arrival->ray_dip = -1.0;
    arrival->ray_qual = 0;
    arrival->pdf_residual_sum = 0.0;
    arrival->pdf_weight_sum = 0.0;

    arrival->fpgrid = NULL;
    arrival->fphdr = NULL;
    arrival->gdesc.buffer = NULL;
    arrival->gdesc.iSwapBytes = iSwapBytesOnInput;
    arrival->sheetdesc.buffer = NULL;

    arrival->station_weight = 1.0;

    // 20180608 AJL - bug fix
    arrival->tt_error = 0.0;


    //DD
    if (nll_mode != MODE_DIFFERENTIAL) {

        arrival->weight = 0.0;
    }


}

/** function to test if arrival is excluded */

int isExcluded(char *label, char *phase) {

    int slen1 = strlen(label); // 20200727 AJL - added so only prefix (e.g. network) of excluded stations can be provided

    for (int nexclude = 0; nexclude < NumLocExclude; nexclude++) {
        //printf("DEBUG: isExcluded <%s> <%s> ? <%s> <%s>\n", label, phase, LocExclude[nexclude].label, LocExclude[nexclude].phase);
        int slen2 = strlen(LocExclude[nexclude].label);
        int slen = slen1 < slen2 ? slen1 : slen2;
        if (strncmp(label, LocExclude[nexclude].label, slen) == 0
                && (strcmp(phase, LocExclude[nexclude].phase) == 0
                || strcmp("*", LocExclude[nexclude].phase) == 0)) // 20191208 AJL - added phase wild-card matching
            return (1);
        //printf("DEBUG: isExcluded NO MATCH!\n");
    }

    // check if not excluded
    // 20200605 AJL - added
    if (NumLocInclude > 0) { // activate if one or more LOCINCLUDE statements are present
        for (int ninclude = 0; ninclude < NumLocInclude; ninclude++) {
            //printf("DEBUG: isIncluded <%s> <%s> ? <%s> <%s>\n", label, phase, LocInclude[ninclude].label, LocInclude[ninclude].phase);
            int slen2 = strlen(LocInclude[ninclude].label);
            int slen = slen1 < slen2 ? slen1 : slen2;
            if (strncmp(label, LocInclude[ninclude].label, slen) == 0
                    && (strcmp(phase, LocInclude[ninclude].phase) == 0
                    || strcmp("*", LocInclude[ninclude].phase) == 0))

                return (0);
            //printf("DEBUG: isIncluded NO MATCH!\n");
        }
        // not explicitly included, must be excluded
        return (1);
    }

    return (0);

}



/** function to determine type P or S of last leg of phase */

// returns 'P' if P, 'S' if S, ' ' otherwise

char lastLegType(ArrivalDesc * arrival) {
    char *c_last_p, *c_last_P, *c_last_s, *c_last_S;
    int i_last_p, i_last_P, i_last_s, i_last_S;


    // get offsets to last char of types p/P or s/S
    c_last_p = strrchr(arrival->phase, 'p');
    i_last_p = c_last_p == NULL ? -1 : (int) (c_last_p - arrival->phase);
    c_last_P = strrchr(arrival->phase, 'P');
    i_last_P = c_last_P == NULL ? -1 : (int) (c_last_P - arrival->phase);
    c_last_s = strrchr(arrival->phase, 's');
    i_last_s = c_last_s == NULL ? -1 : (int) (c_last_s - arrival->phase);
    c_last_S = strrchr(arrival->phase, 'S');
    i_last_S = c_last_S == NULL ? -1 : (int) (c_last_S - arrival->phase);

    // determine type of last char
    i_last_p = i_last_p > i_last_P ? i_last_p : i_last_P;
    i_last_s = i_last_s > i_last_S ? i_last_s : i_last_S;
    if (i_last_p >= 0 && i_last_p > i_last_s)
        return ('P');
    if (i_last_s >= 0 && i_last_s > i_last_p)

        return ('S');
    return (' ');

}

/** function to calculuate simple elevation correction based on surface velcoity */

double CalcSimpleElevCorr(ArrivalDesc *arrival, int narr, double pvel, double svel) {
    double yval_grid, t_surface, t_elev, elev_corr;
    int n_compan, diagnostic;

    // TODO: procedure could  be simplified by using cell slowness from model files fp_model_grid_P, etc., if available

    // Switch debug messages from this function on/off (1/0).
    diagnostic = message_flag >= 3;
    //diagnostic = 1;

    // check for companion
    if ((n_compan = arrival[narr].n_companion) >= 0) {
        if (diagnostic) {
            snprintf(MsgStr, sizeof(MsgStr), "CalcSimpleElevCorr: n_compan=%d", n_compan);
            nll_putmsg(1, MsgStr);
        }
        if ((elev_corr = arrival[n_compan].elev_corr) < 0.0)
            return (0.0);
        // check for specific P or S vel
    } else if (pvel > 0.0 && lastLegType(arrival + narr) == 'P') {
        elev_corr = -arrival[narr].station.depth / pvel;
    } else if (svel > 0.0 && lastLegType(arrival + narr) == 'S') {
        elev_corr = -arrival[narr].station.depth / svel;
        // else check grid type
    } else {
        if (arrival[narr].gdesc.type == GRID_TIME) {
            if (diagnostic) {
                snprintf(MsgStr, sizeof(MsgStr), "CalcSimpleElevCorr: GRID_TIME");
                nll_putmsg(1, MsgStr);
            }
            /* 3D grid */
            if ((t_surface = (double) ReadAbsInterpGrid3d(arrival[narr].fpgrid,
                    &(arrival[narr].gdesc), 0.0, 0.0, 0.0, 0)) < 0.0)
                return (0.0);
            // read pos along X dir because grid may not extend above surface
            if ((t_elev = (double) ReadAbsInterpGrid3d(arrival[narr].fpgrid,
                    &(arrival[narr].gdesc), fabs(arrival[narr].station.depth), 0.0, 0.0, 0)) < 0.0)
                return (0.0);
        } else {
            if (diagnostic) {
                snprintf(MsgStr, sizeof(MsgStr), "CalcSimpleElevCorr: GRID_TIME_2D");
                nll_putmsg(1, MsgStr);
            }
            /* 2D grid (1D model) */
            if ((t_surface = ReadAbsInterpGrid2d(arrival[narr].fpgrid,
                    &(arrival[narr].gdesc), 0.0, 0.0)) < 0.0)
                return (0.0);
            // read pos along Horiz dir because grid may not extend above surface
            yval_grid = fabs(arrival[narr].station.depth);
            if (GeometryMode == MODE_GLOBAL)
                yval_grid *= KM2DEG;
            if ((t_elev = ReadAbsInterpGrid2d(arrival[narr].fpgrid,
                    &(arrival[narr].gdesc), yval_grid, 0.0)) < 0.0)
                return (0.0);
        }
        if (arrival[narr].station.depth > 0.0) // below surface
            t_elev *= -1.0;
        elev_corr = t_elev - t_surface;
    }

    elev_corr *= arrival[narr].tfact;

    if (diagnostic) {

        snprintf(MsgStr, sizeof(MsgStr), "CalcSimpleElevCorr: lat=%.3f  lon=%.3f  depth=%.3f  elev_corr=%.3f",
                arrival[narr].station.dlat, arrival[narr].station.dlong, arrival[narr].station.depth, elev_corr);
        nll_putmsg(1, MsgStr);
    }

    return (elev_corr);

}

/** function to evaluate arrival label aliases */

int EvaluateArrivalAlias(ArrivalDesc * arrival) {
    int nAlias;
    int checkAgain = 1, icount = 0, aliasApplied = 0;
    char *pchr, tmpLabel[MAXLINE];


    snprintf(tmpLabel, sizeof(tmpLabel), "%s", arrival->label);

    if (message_flag >= 4) {
        snprintf(MsgStr, sizeof(MsgStr), "Checking for station name alias: %s", tmpLabel);
        nll_putmsg(4, MsgStr);
    }


    /* evaluate aliases until no replacement done */

    while (checkAgain && icount < MAX_NUM_LOC_ALIAS_CHECKS) {

        checkAgain = 0;
        icount++;

        for (nAlias = 0; nAlias < NumLocAlias; nAlias++) {

            /* check if alias can be rejected */

            if (strcmp(LocAlias[nAlias].name, tmpLabel) != 0)
                continue;
            if (LocAlias[nAlias].byr > arrival->year)
                continue;
            if (LocAlias[nAlias].byr == arrival->year) {
                if (LocAlias[nAlias].bmo > arrival->month)
                    continue;
                if (LocAlias[nAlias].bmo == arrival->month
                        && LocAlias[nAlias].bday > arrival->day)
                    continue;
            }
            if (LocAlias[nAlias].eyr < arrival->year)
                continue;
            if (LocAlias[nAlias].eyr == arrival->year) {
                if (LocAlias[nAlias].emo < arrival->month)
                    continue;
                if (LocAlias[nAlias].emo == arrival->month
                        && LocAlias[nAlias].eday < arrival->day)
                    continue;
            }

            /* apply alias */

            aliasApplied = 1;

            snprintf(tmpLabel, sizeof(tmpLabel), "%s", LocAlias[nAlias].alias);
            if (message_flag >= 3) {
                snprintf(MsgStr, sizeof(MsgStr), " -> %s", tmpLabel);
                nll_putmsg(4, MsgStr);
            }

            /* if alias label is same as arrival label, end check to avoid infinite recursion */
            if (strcmp(tmpLabel, arrival->label) == 0)
                checkAgain = 0;
            else
                checkAgain = 1;

            break;

        }
    }


    /* check for possible recursion in alias specifications */

    if (icount >= MAX_NUM_LOC_ALIAS_CHECKS) {
        if (message_flag >= 4)
            nll_putmsg(4, "");
        nll_puterr("ERROR: possible infinite recursion in station name alias.");
        return (-1);
    }


    /* update arrival label */

    snprintf(arrival->time_grid_label, sizeof(arrival->time_grid_label), "%s", tmpLabel);
    if ((pchr = strrchr(tmpLabel, '_')) != NULL)
        *pchr = '\0';
    //18JUL2002 AJL (IRSN)
    // to prevent station label change in output file
    //	snprintf(arrival->label, sizeof(arrival->label), "%s", tmpLabel);


    /* return if no alias applied */

    if (!aliasApplied) {
        if (message_flag >= 4)
            nll_putmsg(4, "");
        return (0);
    }

    if (message_flag >= 4)
        nll_putmsg(4, "");

    return (0);
}

/** function to apply time delays */

int ApplyTimeDelays(ArrivalDesc * arrival) {
    int nDelay, i;
    int ifound;
    double tmp_delay;

    double tmp_std_dev;

    //* // 20191210 AJL - Bug fix: re-introduced phase mapping which was commented out
    // added SH 23.10.2007
    char eval_phase[PHASE_LABEL_LEN];
    char arrival_phase[PHASE_LABEL_LEN];

    // SH 23.10.2007 do phase mapping before checking for time delays
    snprintf(arrival_phase, sizeof(arrival_phase), "%s", arrival->phase);
    EvalPhaseID(eval_phase, sizeof(eval_phase), arrival_phase);

    //*/

    if (message_flag >= 4) {
        snprintf(MsgStr, sizeof(MsgStr), "Checking for time delay: %s %s",
                arrival->label, arrival_phase);
        nll_putmsg(4, MsgStr);
    }

    arrival->delay = 0.0;

    /* check time delays for match to label/phase */

    ifound = 0;
    for (nDelay = 0; !ifound && nDelay < NumTimeDelays; nDelay++) {

        //printf("DEBUG: comparing delay %s %s to obs %s %s\n", TimeDelay[nDelay].label, TimeDelay[nDelay].phase, arrival->label, eval_phase);

        // 20200405 AJL - check station, and phase with and without phase mapping
        if (strcmp(TimeDelay[nDelay].label, arrival->label) == 0
                && (strcmp(TimeDelay[nDelay].phase, eval_phase) == 0 || strcmp(TimeDelay[nDelay].phase, arrival->phase) == 0)) {
            //printf("DEBUG: SUCCESS!\n");

            /* apply station delays */

            tmp_delay = TimeDelay[nDelay].delay;
            arrival->delay = 0.0;
            if (fabs(tmp_delay) > VERY_SMALL_DOUBLE) {
                // DELAY_CORR			arrival->sec -= tmp_delay;
                arrival->delay = tmp_delay;
                arrival->obs_time -= (long double) arrival->delay; // DELAY_CORR	- incorporating delay so subtract (Tcorr = Tobs - (O-C))
                if (message_flag >= 4) {
                    snprintf(MsgStr, sizeof(MsgStr),
                            "   delay of %lf sec subtracted from obs time.",
                            tmp_delay);
                    nll_putmsg(4, MsgStr);
                }
                ifound = 1;
                // adjust error based on std-dev of delays
                // AJL 20030826 added to test if pick error should be set based on residual statistics
                if (0) {
                    tmp_std_dev = TimeDelay[nDelay].std_dev;
                    if (tmp_std_dev > arrival->error) {
                        if (message_flag >= 4) {
                            snprintf(MsgStr, sizeof(MsgStr), "   error set from %f to %f sec.", arrival->error, tmp_std_dev);
                            nll_putmsg(4, MsgStr);
                        }
                        arrival->error = tmp_std_dev;
                    } else {
                        if (message_flag >= 4) {
                            snprintf(MsgStr, sizeof(MsgStr), "   error not changed from %f to %f sec.", arrival->error, tmp_std_dev);
                            nll_putmsg(4, MsgStr);
                        }
                    }
                }
            }
            break;
        }
    }
    if (message_flag >= 4)
        nll_putmsg(4, "");


    // if time delay not found, check for delay surface
    if (!ifound && NumTimeDelaySurface) {
        tmp_delay = LARGE_FLOAT;
        for (i = 0; i < NumTimeDelaySurface; i++) {
            if (strcmp(eval_phase, TimeDelaySurfacePhase[i]) == 0) {
                tmp_delay = ApplySurfaceTimeDelay(i, arrival);
                tmp_delay *= TimeDelaySurfaceMultiplier[i];
                break;
            }
        }
        if (i < NumTimeDelaySurface && tmp_delay < LARGE_FLOAT / 2.0) {
            arrival->delay = tmp_delay;
            arrival->obs_time -= (long double) arrival->delay; // DELAY_CORR	- incorporating delay so subtract (Tcorr = Tobs - (O-C))
            printf("%s %s %s, ", arrival->label, eval_phase, TimeDelaySurfacePhase[i]);
            if (message_flag >= 1) {

                snprintf(MsgStr, sizeof(MsgStr),
                        "    %s surface delay of %lf sec at lat %f, long %f subtracted from obs time.",
                        TimeDelaySurfacePhase[i], tmp_delay,
                        arrival->station.dlat, arrival->station.dlong);
                //nll_putmsg(4, MsgStr);
                nll_putmsg(1, MsgStr);
            }
        }
    }


    return (0);
}

/** function to get time delay from time delay surface */

double ApplySurfaceTimeDelay(int nsurface, ArrivalDesc * arrival) {
    double x, y;

    if (arrival->station.is_coord_latlon) {
        x = arrival->station.dlong;
        y = arrival->station.dlat;
    } else {

        return (LARGE_FLOAT);
    }

    return (get_surface_z(nsurface, x, y));
}




/** function to extract arrival information observation file name */

/* returns: 0 if nothing done, 1 if info extracted, -1 if unexpected error */

int ExtractFilenameInfo(char *filename, char *type_obs) {
    int istat;
    char *filepos, *extpos;

    if (strcmp(ftype_obs, "RENASS_DEP") == 0) {
        /* find beginning of filename */
        if ((filepos = strrchr(filename, '/')) == NULL)
            return (-1);
        /* get date/time from filename */
        /* try long format (i.e. NICE199801311202.dep) */
        if ((extpos = strstr(filepos, ".dep")) != NULL &&
                extpos - filepos - 12 >= 0) {
            if ((istat = sscanf(extpos - 12, "%4d%2d%2d%2d%2d",
                    &EventTime.year, &EventTime.month, &EventTime.day,
                    &EventTime.hour, &EventTime.min))
                    != 5)
                return (-1);
        }/* try short format (i.e. g504210802.dep) */
        else if ((extpos = strstr(filepos, ".dep")) != NULL &&
                extpos - filepos - 9 >= 0) {
            if ((istat = sscanf(extpos - 9, "%1d%2d%2d%2d%2d",
                    &EventTime.year, &EventTime.month, &EventTime.day,
                    &EventTime.hour, &EventTime.min))
                    != 5)

                return (-1);
            EventTime.year += 1990;
        }

        return (1);
    }

    return (0);

}

/** function to read arrival from observation file */

int GetNextObs(HypoDesc* phypo, FILE* fp_obs, ArrivalDesc *arrival, char* ftype_obs, int nfirst) {

    int istat = 0, iloop;
    char *cstat;
    char chr, instruction[MAXSTRING];
    char chrtmp[MAXSTRING];
    char eval_phase_tmp[PHASE_LABEL_LEN];

    double psec, ssec;

    double weight, ttime;

    int ioff;

    static char line[MAXLINE_LONG];

    static int date_saved, year_save, month_save, day_save;
    static int check_for_S_arrival;

    static int in_hypocenter_event;

    int ifound;

    int itest, itest2;

    // ETH LOC format
    char eth_line_key[MAXSTRING];
    int eth_use_loc = 0, eth_use_mag = 0;
    double vpvs;

    // NEIC / ISC format
    static char last_label[10];
    char cmonth[4];
    char* pchr;
    static int origin_hour = 0;

    // ISC format
    char isc_time_str[10];

    // DD
    static int hypo_cc_flag;
    static long int dd_event_id_1, dd_event_id_2;
    static double dd_otime_corr;
    double tt_sta1, tt_sta2;

    // HYPOINVERSE_Y2000_ARC
    int idummy;

    // SAFOD
    int yearday;
    double left_uncertainty, right_uncertainty;


    /* if no obs read for this event, set date saved flag to 0 */
    if (nfirst) {
        date_saved = 0;
        check_for_S_arrival = 0;
        in_hypocenter_event = 0;
    }


    /* check for special control instructions */
    iloop = 1;
    while (iloop && !check_for_S_arrival) {
        iloop = 0;
        if ((chr = fgetc(fp_obs)) == '!'/* || chr == '#'*/) {
            ungetc(chr, fp_obs);
            if (fgets(line, MAXLINE_LONG, fp_obs) == NULL)
                return (OBS_FILE_END_OF_INPUT);
            printf("1 %s", line);
            /* read instruction */
            istat = sscanf(line + 1, "%s", instruction);
            if (istat == EOF)
                return (OBS_FILE_END_OF_INPUT);
            if (istat == 1) {
                if (strcmp(instruction, "END_EVENT") == 0) {
                    return (OBS_FILE_END_OF_EVENT);
                } else if (strcmp(instruction, "END_FILE") == 0) {
                    return (OBS_FILE_END_OF_INPUT);
                } else if (strcmp(instruction, "SKIP_NEXT_LINE") == 0) {
                    if (fgets(line, MAXLINE_LONG, fp_obs) == NULL)
                        return (OBS_FILE_END_OF_INPUT);
                    printf("2 %s", line);

                    iloop = 1;
                    continue;
                } else { // skip this line
                    iloop = 1;
                    continue;
                }
            }
            if (istat != 1) {
                nll_puterr2("WARNING: unrecognized control statement", line);
            }
        } else {
            ungetc(chr, fp_obs);
        }
    }


    /* set field defaults */
    snprintf(arrival->label, sizeof(arrival->label), "%s", ARRIVAL_NULL_STR);
    snprintf(arrival->network, sizeof(arrival->network), "%s", ARRIVAL_NULL_STR);
    snprintf(arrival->inst, sizeof(arrival->inst), "%s", ARRIVAL_NULL_STR);
    snprintf(arrival->comp, sizeof(arrival->comp), "%s", ARRIVAL_NULL_STR);
    snprintf(arrival->onset, sizeof(arrival->onset), "%s", ARRIVAL_NULL_STR);
    snprintf(arrival->phase, sizeof(arrival->phase), "%s", ARRIVAL_NULL_STR);
    snprintf(arrival->first_mot, sizeof(arrival->first_mot), "%s", ARRIVAL_NULL_STR);
    arrival->quality = 99;
    arrival->first_mot_quality = 1.0; // 20200829 AJL - is initialized to 1.0 but may be changed (e.g. )
    snprintf(arrival->error_type, sizeof(arrival->error_type), "%s", "GAU");
    arrival->error = ARRIVAL_ERROR_NULL;
    arrival->coda_dur = CODA_DUR_NULL;
    arrival->amplitude = AMPLITUDE_NULL;
    arrival->period = PERIOD_NULL;
    arrival->amp_mag = MAGNITUDE_NULL;
    arrival->dur_mag = MAGNITUDE_NULL;
    arrival->apriori_weight = 1.0;

    arrival->dd_event_id_1 = -1;
    arrival->flag_ignore = 0; // 20150521 AJL

    // 20211211 AJL - Bug fix?
    arrival->sheetdesc.array = NULL;
    arrival->sheetdesc.buffer = NULL;

    /* attempt to read obs based on obs file type */

    if (strncmp(ftype_obs, "NLLOC_OBS", 9) == 0) {

        // *_LOCPHASEID - convert phase name using LOCPHASEID (homogenizes names for LOCDELAY accumulation)

        /* read next line */
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);

        // check for event public_id (assume listed before arrivals, e.g. in NLL Hypocenter-Phase file)
        // 20190823 AJL - added
        if (sscanf(line, "PUBLIC_ID %s", phypo->public_id) == 1) {
            return (OBS_FILE_SKIP_INPUT_LINE);
        }
        // check for event QUALITY (assume listed before arrivals, e.g. in NLL Hypocenter-Phase file)
        // 20191019 AJL - added
        if (sscanf(line,
                "QUALITY %*s %Lf %*s %lf %*s %lf %*s %lf %*s %d %*s %lf %*s %lf %*s %lf %d %*s %lf %d",
                &phypo->probmax, &phypo->misfit,
                &phypo->grid_misfit_max,
                &phypo->rms, &phypo->nreadings, &phypo->gap,
                &phypo->dist,
                &phypo->amp_mag, &phypo->num_amp_mag,
                &phypo->dur_mag, &phypo->num_dur_mag
                ) == 11) {
            return (OBS_FILE_SKIP_INPUT_LINE);
        }
        // check for event FOCALMECH (assume listed before arrivals, e.g. in NLL Hypocenter-Phase file)
        // 20200218 AJL - added
        /* FOCALMECH */
        /* Hyp dlat dlong depth Mech dipDir dipAng rake mf misfit nObs nObs */
        if (sscanf(line,
                "FOCALMECH %*s %lf %lf %lf %*s %lf %lf %lf %*s %lf %*s %d",
                &phypo->focMech.dlat, &phypo->focMech.dlong,
                &phypo->focMech.depth,
                &phypo->focMech.dipDir, &phypo->focMech.dipAng,
                &phypo->focMech.rake,
                &phypo->focMech.misfit, &phypo->focMech.nObs
                ) == 8) {
            return (OBS_FILE_SKIP_INPUT_LINE);
        }

        istat = ReadArrival(line, arrival, IO_ARRIVAL_OBS);
        if (istat < 1) {
            if (in_hypocenter_event) {
                return (OBS_FILE_END_OF_EVENT);
            } else {
                return (OBS_FILE_SKIP_INPUT_LINE);
            }
        }

        in_hypocenter_event = 1;

        /* convert error to quality */
        if ((arrival->quality = Err2Qual(arrival)) < 0)
            arrival->quality = 99;

        // convert phase name using LOCPHASEID if requested (homogenizes names for LOCDELAY accumulation)
        // 20200812 AJL - added
        if (strstr(ftype_obs, "_LOCPHASEID") != NULL) {
            //printf("DEBUG: arrival->phase %s", arrival->phase);
            EvalPhaseID(eval_phase_tmp, sizeof(eval_phase_tmp), arrival->phase);
            snprintf(arrival->phase, sizeof(arrival->phase), "%s", eval_phase_tmp);
            //printf(" -> arrival->phase %s\n", arrival->phase);
        }

        return (istat);

    } else if (strncmp(ftype_obs, "INGV_JSON", 9) == 0) {

        // INGV json format 20211012
        // assumes a single event in each input file

        istat = json_read_next_arrival_INGV(phypo, fp_obs, arrival, nfirst);
        if (istat < 1) {
            return (OBS_FILE_END_OF_INPUT);
        }

        return (istat);

    } else if (strcmp(ftype_obs, "UUSS") == 0) {
        /* SH Seismograph Station University of Utah format (PING)
                                    first line: reftime (yymmddhhmm)
                                    next lines arrival times rel. to reftime incl. duration
         */
        /* interpret line: reftime or phase information or comment line etc */
        chr = fgetc(fp_obs);
        /* fprintf(stderr,"GetNextObs: chr %d\n",chr); */
        if (chr == EOF) {
            return (OBS_FILE_END_OF_INPUT);
        } else if (chr == 99 || chr == 67) {
            /* comment line starts with 'c' */
            /* read until end of line */
            ungetc(chr, fp_obs);
            cstat = fgets(line, MAXLINE_LONG, fp_obs);
            return (OBS_IS_COMMENT_LINE);
        } else if (chr == 35) {
            /* end of event marked by '#' */
            /* read until end of line */
            ungetc(chr, fp_obs);
            cstat = fgets(line, MAXLINE_LONG, fp_obs);
            return (OBS_FILE_END_OF_EVENT);
        } else if (chr == 32) {
            /* line contains reference time starting with a space */
            ungetc(chr, fp_obs);
            cstat = fgets(line, MAXLINE_LONG, fp_obs);
            istat = ReadFortranInt(line, 2, 2, &EventTime.year);
            if (EventTime.year < 20)
                EventTime.year += 100;
            EventTime.year += 1900;
            istat += ReadFortranInt(line, 4, 2, &EventTime.month);
            istat += ReadFortranInt(line, 6, 2, &EventTime.day);
            istat += ReadFortranInt(line, 8, 2, &EventTime.hour);
            istat += ReadFortranInt(line, 10, 2, &EventTime.min);
            arrival->sec = 0.0;
        } else {
            /* line contains phase information */
            ungetc(chr, fp_obs);
        }

        /* read next line */
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);
        /* now read phase info */
        istat = ReadFortranString(line, 1, 4, arrival->label);
        /* check if comment line */
        if ((strncmp(arrival->label, "c", 1) == 0) || (strncmp(arrival->label, "C", 1) == 0)) {
            return (OBS_IS_COMMENT_LINE);
        }
        TrimString(arrival->label);
        snprintf(arrival->phase, sizeof(arrival->phase), "%s", "P");
        if ((istat += ReadFortranString(line, 6, 1, arrival->first_mot)) != 2) {
            /* read next lines until arrival time is found */
            while ((cstat = fgets(line, MAXLINE_LONG, fp_obs)) != NULL) {
                if (strncmp(line, "#", 1) == 0)
                    /* end of event detected */
                    return (OBS_FILE_END_OF_EVENT);
                istat = 0;
                /* fprintf(stderr,"GetNextObs: %s\n",line); */
                /* fprintf(stderr,"GetNextObs: %d %s\n",istat,arrival->first_mot); */
                if ((istat += ReadFortranString(line, 6, 1, arrival->first_mot)) == 1) {
                    istat = ReadFortranString(line, 1, 4, arrival->label);
                    if ((strncmp(arrival->label, "c", 1) == 0) || (strncmp(arrival->label, "C", 1) == 0)) {
                        return (OBS_IS_COMMENT_LINE);
                    }
                    TrimString(arrival->label);
                    break;
                }
            }
        }

        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);

        /* continue with reading of remaining phase information */
        istat += ReadFortranInt(line, 10, 1, &arrival->quality);
        istat += ReadFortranReal(line, 12, 6, &arrival->sec);

        arrival->year = EventTime.year;
        arrival->month = EventTime.month;
        arrival->day = EventTime.day;
        arrival->hour = EventTime.hour;
        arrival->min = EventTime.min;

        /* convert quality to error */
        Qual2Err(arrival);

        /* fprintf(stderr,"GetNextObs: phase %4s %1s %1s %1d %5.2f\n",arrival->label,
        arrival->phase,arrival->first_mot,arrival->quality,arrival->sec); */

        /* check for possible S-arrival */
        istat += ReadFortranReal(line, 21, 6, &ssec);
        if (ssec > 0.0) {
            arrival++;
            istat = ReadFortranString(line, 1, 4, arrival->label);
            TrimString(arrival->label);
            snprintf(arrival->phase, sizeof(arrival->phase), "%s", "S");
            istat += ReadFortranInt(line, 19, 1, &arrival->quality);
            arrival->year = EventTime.year;
            arrival->month = EventTime.month;
            arrival->day = EventTime.day;
            arrival->hour = EventTime.hour;
            arrival->min = EventTime.min;
            arrival->sec = ssec;
            /* convert quality to error */
            Qual2Err(arrival);
            return (OBS_FILE_TWO_ARRIVALS_READ);
        } else
            return (istat);


    } else if (strcmp(ftype_obs, "HYPO71") == 0 ||
            strcmp(ftype_obs, "HYPO71_OV") == 0 ||
            strcmp(ftype_obs, "HYPO71_S_QUAL_PLUS_1") == 0 ||
            strcmp(ftype_obs, "HYPOELLIPSE") == 0) {

        if (check_for_S_arrival) {
            /* check for S phase input in last input line read */

            /* set S read offset to allow correction of incorrect hypo71 format */
            ioff = 0;
            if (strcmp(ftype_obs, "HYPO71_OV") == 0) {
                istat = ReadFortranString(line, 31, 1, chrtmp);
                if (strcmp(chrtmp, " ") != 0)
                    ioff = -1;
            }

            /* check for zero or blank S phase time */
            istat = ReadFortranString(line, 32 + ioff, 5, chrtmp);
            if (istat > 0
                    && strncmp(chrtmp, "     ", 5) != 0
                    && strncmp(chrtmp, "    0", 5) != 0) {
                /* read S phase input in last input line read */
                istat = ReadFortranString(line, 1, 4, arrival->label);
                TrimString(arrival->label);
                istat += ReadFortranInt(line, 10, 2, &arrival->year);
                // temporary fix for HYPO71 Y2K problem
                // 20231017 //if (arrival->year < 20)
                if (arrival->year < 40)
                    arrival->year += 100;
                arrival->year += 1900;
                istat += ReadFortranInt(line, 12, 2, &arrival->month);
                istat += ReadFortranInt(line, 14, 2, &arrival->day);
                istat += ReadFortranInt(line, 16, 2, &arrival->hour);
                istat += ReadFortranInt(line, 18, 2, &arrival->min);
                istat += ReadFortranReal(line, 32 + ioff, 5, &arrival->sec);
                /* 20090126 AJL bug fix - original version fails if integer sec < 1.00
                // check for integer sec format
                                if (arrival->sec > 99.999)
                                arrival->sec /= 100.0;
                 */
                strncpy(chrtmp, line + 31 + ioff, 5);
                chrtmp[5] = '\0';
                //printf("%s", line);
                //printf("chrtmp=|%s| arrival->sec=%f", chrtmp, arrival->sec);
                if (strchr(chrtmp, '.') == NULL)
                    arrival->sec /= 100.0;
                //printf(" -> arrival->sec=%f\n", arrival->sec);
                // END - 20090126 AJL bug fix
                istat += ReadFortranReal(line, 20, 5, &psec);
                // 20140528 AJL - bug fix, check for integer P time
                strncpy(chrtmp, line + 19, 5);
                chrtmp[5] = '\0';
                if (strchr(chrtmp, '.') == NULL)
                    psec /= 100.0;
                // END - 20140528 AJL - bug fix
                /* check for P second >= 60.0 */
                if (psec >= 60.0 && arrival->sec < 60.0)
                    arrival->sec += 60.0;
                arrival->phase[0] = '\0';
                istat += ReadFortranString(line, 38 + ioff, 1, arrival->phase);
                TrimString(arrival->phase);
                istat += ReadFortranInt(line, 40 + ioff, 1, &arrival->quality);
            }
        }

        /* check for S arrival input found */
        if (check_for_S_arrival && istat == 10
                && IsPhaseID(arrival->phase, "S")
                && IsGoodDate(arrival->year,
                arrival->month, arrival->day)) {

            //snprintf(arrival->phase, sizeof(arrival->phase), "%s", "S");

            /* set error fields */
            if (strcmp(ftype_obs, "HYPO71_S_QUAL_PLUS_1") == 0) {
                arrival->quality += 1;
            }
            snprintf(arrival->error_type, sizeof(arrival->error_type), "%s", "GAU");
            if (arrival->quality >= 0 &&
                    arrival->quality < NumQuality2ErrorLevels) {
                arrival->error =
                        Quality2Error[arrival->quality];
            } else {
                arrival->error =
                        Quality2Error[NumQuality2ErrorLevels - 1];
                nll_puterr("WARNING: invalid arrival weight.");
            }

            line[0] = '\0';
            check_for_S_arrival = 0;
            return (istat);
        } else if (check_for_S_arrival) {
            check_for_S_arrival = 0;
            return (OBS_FILE_SKIP_INPUT_LINE);
        }




        /* read next line */
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);


        // 20231021 AJL - Added support for NLL header lines
        // check for event public_id (assume listed before arrivals, e.g. in NLL Hypocenter-Phase file)
        // 20190823 AJL - added
        if (sscanf(line, "PUBLIC_ID %s", phypo->public_id) == 1) {
            return (OBS_FILE_SKIP_INPUT_LINE);
        }
        // check for event QUALITY (assume listed before arrivals, e.g. in NLL Hypocenter-Phase file)
        // 20191019 AJL - added
        if (sscanf(line,
                "QUALITY %*s %Lf %*s %lf %*s %lf %*s %lf %*s %d %*s %lf %*s %lf %*s %lf %d %*s %lf %d",
                &phypo->probmax, &phypo->misfit,
                &phypo->grid_misfit_max,
                &phypo->rms, &phypo->nreadings, &phypo->gap,
                &phypo->dist,
                &phypo->amp_mag, &phypo->num_amp_mag,
                &phypo->dur_mag, &phypo->num_dur_mag
                ) == 11) {
            return (OBS_FILE_SKIP_INPUT_LINE);
        }


        /* read formatted P (or S) arrival input */
        istat = ReadFortranString(line, 1, 4, arrival->label);
        TrimString(arrival->label);
        istat += ReadFortranString(line, 5, 1, arrival->onset);
        TrimString(arrival->onset);
        istat += ReadFortranString(line, 6, 1, arrival->phase);
        TrimString(arrival->phase);
        istat += ReadFortranString(line, 7, 1, arrival->first_mot);
        TrimString(arrival->first_mot);
        istat += ReadFortranInt(line, 8, 1, &arrival->quality);
        istat += ReadFortranInt(line, 10, 2, &arrival->year);
        // temporary fix for HYPO71 Y2K problem
        // 20231017 //if (arrival->year < 20)
        if (arrival->year < 40)
            arrival->year += 100;
        arrival->year += 1900;
        istat += ReadFortranInt(line, 12, 2, &arrival->month);
        istat += ReadFortranInt(line, 14, 2, &arrival->day);
        istat += ReadFortranInt(line, 16, 2, &arrival->hour);
        istat += ReadFortranInt(line, 18, 2, &arrival->min);
        istat += ReadFortranReal(line, 20, 5, &arrival->sec);

        if (istat != 11) {
            line[0] = '\0';
            return (OBS_FILE_END_OF_EVENT);
        }

        /* read optional amplitude/period fields */
        istat += ReadFortranReal(line, 44, 4, &arrival->amplitude);
        //printf("AMPLITUDE: %lf (%s)\n", arrival->amplitude, line);
        istat += ReadFortranReal(line, 48, 3, &arrival->period);
        // 20030826 AJL coda_dur added
        istat += ReadFortranReal(line, 71, 5, &arrival->coda_dur);


        check_for_S_arrival = 1;
        /* check for valid phase code */
        //		if (IsPhaseID(arrival->phase, "P")) {
        //snprintf(arrival->phase, sizeof(arrival->phase), "%s", "P");
        check_for_S_arrival = 1;

        //		} else if (IsPhaseID(arrival->phase, "S")) {
        //snprintf(arrival->phase, sizeof(arrival->phase), "%s", "S");
        //			check_for_S_arrival = 0;
        //		} else
        //			return(OBS_FILE_END_OF_EVENT);

        if (!IsGoodDate(arrival->year, arrival->month, arrival->day))
            return (OBS_FILE_END_OF_EVENT);


        /* 20090126 AJL bug fix - original version fails if integer sec < 1.00
        // check for integer sec format
                          if (arrival->sec > 99.999)
                          arrival->sec /= 100.0;
         */
        strncpy(chrtmp, line + 19, 5);
        chrtmp[5] = '\0';
        //printf("%s", line);
        //printf("chrtmp=|%s| arrival->sec=%f", chrtmp, arrival->sec);
        if (strchr(chrtmp, '.') == NULL)
            arrival->sec /= 100.0;
        //printf(" -> arrival->sec=%f\n", arrival->sec);
        // END - 20090126 AJL bug fix

        /* convert quality to error */
        // check if arrival is S
        if (IsPhaseID(arrival->phase, "S") && strcmp(ftype_obs, "HYPO71_S_QUAL_PLUS_1") == 0) {
            arrival->quality += 1;
        }
        Qual2Err(arrival);

        return (istat);

    } else if (strcmp(ftype_obs, "HYPOINVERSE_Y2000_ARC") == 0 ||
            strcmp(ftype_obs, "HYPOINVERSE_Y2000_ARC_NET") == 0) // concatenate network with station: NN_SSS
    {

        /*
                        Y2000 (station) archive format
                        StartCol. Len. FortranFormat Data (* revised from pre-Y2000 format)
                        1 5 A5 5-letter station site code, left justified. *
                        6 2 A2 2-letter seismic network code. *
                        8 1 1X Blank *
                        9 1 A1 One letter station component code.
                        10 3 A3 3-letter station component code. *
                        13 1 1X Blank *
                        14 2 A2 P remark such as "IP".
                        16 1 A1 P first motion.
                        17 1 I1 Assigned P weight code.
                        18 4 I4 Year. *
                        22 8 4I2 Month, day, hour and minute.
                        30 5 F5.2 Second of P arrival.
                        35 4 F4.2 P travel time residual.
                        39 3 F3.2 P weight actually used.
                        42 5 F5.2 Second of S arrival.
                        47 2 A2 S remark such as "ES".
                        49 1 1X Blank
                        50 1 I1 Assigned S weight code.
                        51 4 F4.2 S travel time residual.
                        55 7 F7.2 Amplitude (Normally peak-to-peak). *
                        62 2 I2 Amp units code. 0=PP mm, 1=0 to peak mm (UCB), 2=digital counts. *
                        64 3 F3.2 S weight actually used.
                        67 4 F4.2 P delay time.
                        71 4 F4.2 S delay time.
                        75 4 F4.1 Epicentral distance (km).
                        79 3 F3.0 Emergence angle at source.
                        82 1 I1 Amplitude magnitude weight code.
                        83 1 I1 Duration magnitude weight code.
                        84 3 F3.2 Period at which the amplitude was measured for this station.
                        87 1 A1 1-letter station remark.
                        88 4 F4.0 Coda duration in seconds.
                        92 3 F3.0 Azimuth to station in degrees E of N.
                        95 3 F3.2 Duration magnitude for this station. *
                        98 3 F3.2 Amplitude magnitude for this station. *
                        101 4 F4.3 Importance of P arrival.
                        105 4 F4.3 Importance of S arrival.
                        109 1 A1 Data source code.
                        110 1 A1 Label code for duration magnitude from FC1 or FC2 command.
                        111 1 A1 Label code for amplitude magnitude from XC1 or XC2 command.
                        112 2 A2 2-letter station location code (component extension).
                        113 is the last filled column.

        Example HYPOINVERSE doc
                :
        19991201 0 8 64735 5775120 3094 355 0 160
        PST NC VVHZ IPD0199912010008 7.56 00 0 4 0 0
        PMM NC VVHZ EPU2199912010008 7.71 700 0 4 19M 0 0
        PVC NC VVHZ IPD1199912010008 7.81 8.98ESD2 00 0 0 6 0
        PPC NC VVHZ IPU0199912010008 8.13 9.19ESU2 1600 0 4 7M 0 0
        PHP NC VVHZ EPU2199912010008 8.39 00 0 4 0 0
        PSM NC VVHZ EPU2199912010008 9.98 00 0 4 0 0
        BMS NC VVHZ P 419991201000899.99 00 0 4 22 0
        21069577

        Example INGV:

        200512132009226042 4387 12E4874   90  0 15137 24 11620767 507 7814 311298    203    0 301 469 12   0 120  0 15Ita WW   16    0  0   0  0    141429 298 120   0   01
        $1
        LNSS IV EBHE P U0200512132009 2792  62181    0   0   0      0 0  0   0   0 235 9200  0   59127315  0 498   0W     0
        $   6 5.20 1.80 6.40 2.65 0.06 NSR0  -47 PHP0  820 151886 19 943 23 497 31 329 39 159 47  91
        LNSS IV EBHE P U2200512132009 2794  64 90    0   0   0      0 0  0   0   0 235 9200  0   72127334  0 124   0W     0
        $   6 5.18 1.80 5.50 2.05 0.04 NSR0  -47 PHP0  954 151220 19 636 23 531 31 289 39 172 47 105
        MNS  IV EBHE P  3200512132009 2982 -75 36    0   0   0      0 0  0   0   0 399 9100  0   49196299  0  40   0W     0
        $   6 4.91 1.80 5.38 2.18 0.07 NSR0  -43 PHP0   38 15 525 19 483 23 258 31 168 39  80 43  66

         */

        if (check_for_S_arrival) {
            /* check for S phase input in last input line read */

            // check for zero or blank S phase time and remark
            // 42 5 F5.2 Second of S arrival.
            // 47 2 A2 S remark such as "ES".
            istat = ReadFortranString(line, 42, 7, chrtmp);
            //printf("chrtmp <%s>\n", chrtmp);
            if (istat > 0
                    && strncmp(chrtmp, "     ", 5) != 0
                    && strncmp(chrtmp, "    0", 5) != 0
                    && strncmp(chrtmp + 5, "  ", 2) != 0) {
                //printf("ACCEPT chrtmp <%s>\n", chrtmp);
                // read S phase input in last input line read
                istat = ReadFortranString(line, 1, 5, arrival->label);
                TrimString(arrival->label);
                istat = ReadFortranString(line, 6, 2, arrival->network); // network code ignored in NLL
                int clen_net = TrimString(arrival->network);
                if (clen_net > 0 && strcmp(ftype_obs, "HYPOINVERSE_Y2000_ARC_NET") == 0) {
                    snprintf(chrtmp, sizeof(chrtmp), "%s", arrival->label);
                    snprintf(arrival->label, sizeof(arrival->label), "%s", arrival->network);
                    strncat(arrival->label, "_", sizeof(arrival->label) - strlen(arrival->label) - 1);
                    strncat(arrival->label, chrtmp, sizeof(arrival->label) - strlen(arrival->label) - 1);
                }
                istat += ReadFortranString(line, 9, 1, arrival->comp);
                TrimString(arrival->comp);
                istat += ReadFortranString(line, 10, 3, arrival->inst);
                TrimString(arrival->inst);
                // try and decode S remark
                arrival->onset[0] = ARRIVAL_NULL_CHR;
                istat += ReadFortranString(line, 47, 2, chrtmp);
                if (chrtmp[0] == 'S' || chrtmp[0] == 's') {
                    snprintf(arrival->phase, sizeof(arrival->phase), "%s", chrtmp);
                    TrimString(arrival->phase);
                } else if (chrtmp[1] == 'S' || chrtmp[1] == 's') {
                    arrival->onset[0] = chrtmp[0];
                    arrival->phase[0] = chrtmp[1];
                } else {
                    arrival->phase[0] = 'S';
                }
                //istat += ReadFortranString(line, 49, 1, arrival->first_mot);
                //TrimString(arrival->first_mot);
                istat += ReadFortranInt(line, 18, 4, &arrival->year);
                istat += ReadFortranInt(line, 22, 2, &arrival->month);
                istat += ReadFortranInt(line, 24, 2, &arrival->day);
                istat += ReadFortranInt(line, 26, 2, &arrival->hour);
                istat += ReadFortranInt(line, 28, 2, &arrival->min);
                istat += ReadFortranReal(line, 42, 5, &arrival->sec);
                /* 20090126 AJL bug fix - original version fails if integer sec < 1.00
                // check for integer sec format
                        if (arrival->sec > 99.999)
                        arrival->sec /= 100.0;
                 */
                strncpy(chrtmp, line + 41, 5);
                chrtmp[5] = '\0';
                //printf("%s", line);
                //printf("chrtmp=|%s| arrival->sec=%f\n", chrtmp, arrival->sec);
                //printf("chrtmp[4]=|%c|\n", chrtmp[4]);
                if (strchr(chrtmp, '.') == NULL) {
                    arrival->sec /= 100.0;
                }
                //printf(" -> arrival->sec=%f\n", arrival->sec);
                // END - 20090126 AJL bug fix
                istat += ReadFortranReal(line, 30, 5, &psec);
                // 2020/08/11 - JMS
                // test for dummy P-phase (psec = 9999), if yes, ignore psec
                ReadFortranInt(line, 30, 5, &itest);
                if (itest == 9999) {
                    psec = 0.0;
                }
                // END - 2020/08/11 - JMS
                /* 20090126 AJL bug fix - original version fails if integer sec < 1.00
                // check for integer sec format
                                if (psec > 99.999)
                                psec /= 100.0;
                 */
                strncpy(chrtmp, line + 29, 5);
                chrtmp[5] = '\0';
                //printf("%s", line);
                //printf("chrtmp=|%s| arrival->sec=%f\n", chrtmp, psec);
                //printf("chrtmp[4]=|%c|\n", chrtmp[4]);
                if (strchr(chrtmp, '.') == NULL) {
                    psec /= 100.0;
                }
                //printf(" -> psec=%f\n", psec);
                // END - 20090126 AJL bug fix
                // check for P second >= 60.0 - not needed (?)
                if (psec >= 60.0 && arrival->sec < 60.0)
                    arrival->sec += 60.0;
                istat += ReadFortranInt(line, 50, 1, &arrival->quality);
                // read optional amplitude/period fields - read with P
                //istat += ReadFortranReal(line, 55, 7, &arrival->amplitude);
                //istat += ReadFortranReal(line, 84, 3, &arrival->period);
                //istat += ReadFortranReal(line, 88, 4, &arrival->coda_dur);
            }

            // check for valid S arrival input found
            if (istat == 12
                    && IsPhaseID(arrival->phase, "S")
                    && IsGoodDate(arrival->year, arrival->month, arrival->day)) {
                // convert quality to error
                Qual2Err(arrival);
                line[0] = '\0';
                check_for_S_arrival = 0;
                return (istat);
            } else { // not found
                check_for_S_arrival = 0;
                return (OBS_FILE_SKIP_INPUT_LINE);
            }

        }


        // read next line
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);

        // check for shawdow card - line starting with $
        if (strncmp(line, "$", 1) == 0)
            return (OBS_FILE_SKIP_INPUT_LINE);

        // check for summary header - line starting with year
        if (ReadFortranInt(line, 1, 4, &idummy) == 1 && idummy > 0) {
            snprintf(HypoInverseArchiveSumHdr, sizeof(HypoInverseArchiveSumHdr), "%s", line);
            //printf("HypoInverseArchiveSumHdr:\n|%s|\n", HypoInverseArchiveSumHdr);
            // 20220129 AJL - added reading of preferred magnitude as NLL amp mag.
            // read magnitude
            // 148 3 F3.2 Preferred magnitude, chosen by the Hypoinverse PRE command.
            // 151 4 F4.1 Total of the preferred mag weights (~ number of readings).
            //  fprintf(fpio, "X%3.0lf%4.0lf", 100.0 * amp_mag, 10.0 * amp_mag_wt);
            sscanf(HypoInverseArchiveSumHdr + 147, "%3lf", &phypo->amp_mag);
            phypo->amp_mag /= 100.0;
            // 20220129 AJL - added reading of Event identification number as NLL public_id.
            // 137 10 I10 Event identification number
            strncpy(chrtmp, HypoInverseArchiveSumHdr + 136, 11);
            chrtmp[10] = '\0';
            sscanf(chrtmp, "%s", phypo->public_id);
            return (OBS_FILE_SKIP_INPUT_LINE);
        }

        // 2020/08/11 - JMS
        // check for dummy P-phase (second 9999 and weight of 4)
        ReadFortranInt(line, 30, 5, &itest);
        ReadFortranInt(line, 17, 1, &itest2);
        // if dummy P-phase, skip reading and try read S-phase
        if (itest == 9999 && itest2 == 4) {
            check_for_S_arrival = 1;
            return (OBS_FILE_SKIP_INPUT_LINE);
        } else {
            // read formatted P arrival input
            istat = ReadFortranString(line, 1, 5, arrival->label);
            TrimString(arrival->label);
            istat = ReadFortranString(line, 6, 2, arrival->network); // network code ignored in NLL
            int clen_net = TrimString(arrival->network);
            if (clen_net > 0 && strcmp(ftype_obs, "HYPOINVERSE_Y2000_ARC_NET") == 0) {
                snprintf(chrtmp, sizeof(chrtmp), "%s", arrival->label);
                snprintf(arrival->label, sizeof(arrival->label), "%s", arrival->network);
                strncat(arrival->label, "_", sizeof(arrival->label) - strlen(arrival->label) - 1);
                strncat(arrival->label, chrtmp, sizeof(arrival->label) - strlen(arrival->label) - 1);
            }
            istat += ReadFortranString(line, 9, 1, arrival->comp);
            TrimString(arrival->comp);
            istat += ReadFortranString(line, 10, 3, arrival->inst);
            TrimString(arrival->inst);
            // try and decode P remark
            arrival->onset[0] = ARRIVAL_NULL_CHR;
            istat += ReadFortranString(line, 14, 2, chrtmp);
            if (chrtmp[0] == 'P' || chrtmp[0] == 'p') {
                snprintf(arrival->phase, sizeof(arrival->phase), "%s", chrtmp);
                TrimString(arrival->phase);
            } else if (chrtmp[1] == 'P' || chrtmp[1] == 'p') {
                arrival->onset[0] = chrtmp[0];
                arrival->phase[0] = chrtmp[1];
            } else {
                arrival->phase[0] = 'P';
            }
            TrimString(arrival->phase);
            istat += ReadFortranString(line, 16, 1, arrival->first_mot);
            TrimString(arrival->first_mot);
            istat += ReadFortranInt(line, 18, 4, &arrival->year);
            istat += ReadFortranInt(line, 22, 2, &arrival->month);
            istat += ReadFortranInt(line, 24, 2, &arrival->day);
            istat += ReadFortranInt(line, 26, 2, &arrival->hour);
            istat += ReadFortranInt(line, 28, 2, &arrival->min);
            istat += ReadFortranReal(line, 30, 5, &arrival->sec);
            /* 20090126 AJL bug fix - original version fails if integer sec < 1.00
            // check for integer sec format
            if (arrival->sec > 99.999)
                    arrival->sec /= 100.0;
             */
            strncpy(chrtmp, line + 29, 5);
            chrtmp[5] = '\0';
            //printf("%s", line);
            //printf("chrtmp=|%s| arrival->sec=%f\n", chrtmp, arrival->sec);
            //printf("chrtmp[4]=|%c|\n", chrtmp[4]);
            if (strchr(chrtmp, '.') == NULL) {
                arrival->sec /= 100.0;
            }
            //printf(" -> arrival->sec=%f\n", arrival->sec);
            // END - 20090126 AJL bug fix
            istat += ReadFortranInt(line, 17, 1, &arrival->quality);

            if (istat != 12) {
                line[0] = '\0';
                return (OBS_FILE_END_OF_EVENT);
            }

            // read optional amplitude/period fields
            istat += ReadFortranReal(line, 55, 7, &arrival->amplitude);
            istat += ReadFortranReal(line, 84, 3, &arrival->period);
            istat += ReadFortranReal(line, 88, 4, &arrival->coda_dur);


            /* check for valid phase code */
            //		if (IsPhaseID(arrival->phase, "P")) {
            //snprintf(arrival->phase, sizeof(arrival->phase), "%s", "P");
            check_for_S_arrival = 1;

            //		} else if (IsPhaseID(arrival->phase, "S")) {
            //snprintf(arrival->phase, sizeof(arrival->phase), "%s", "S");
            //			check_for_S_arrival = 0;
            //		} else
            //			return(OBS_FILE_END_OF_EVENT);

            if (!IsGoodDate(arrival->year, arrival->month, arrival->day))
                return (OBS_FILE_END_OF_EVENT);

            // convert quality to error
            Qual2Err(arrival);

            return (istat);
        }

    } else if (strcmp(ftype_obs, "INGV_BOLL") == 0 ||
            strcmp(ftype_obs, "INGV_BOLL_LOCAL") == 0 ||
            strcmp(ftype_obs, "INGV_ARCH") == 0) {
        /*
                        1 0
                        DOI 1016   EZPG  03243563 EZSG  243882                             23
                        11 0
                        SEI 1016   EZPG  04132382 EZSG  132881
                        VMG 1016   EZPG  04132514 EZSG  133068                             56
                        ZCCA1016   EZPG  04132655                                          47
         */
        if (check_for_S_arrival > 0) {
            /* check for additional phase input in last input line read */

            ioff = 0;
            if (strcmp(ftype_obs, "INGV_ARCH") == 0)
                ioff = 13 * (check_for_S_arrival - 1);
            /* check for zero or blank S phase time */
            istat = ReadFortranString(line, 33 + ioff, 6, chrtmp);
            if (istat > 0
                    && strncmp(chrtmp, "      ", 6) != 0
                    && strncmp(chrtmp, "     0", 6) != 0) {
                /* read S phase input in last input line read */
                istat = ReadFortranString(line, 1, 4, arrival->label);
                TrimString(arrival->label);
                //istat += ReadFortranInt(line, 10, 2, &arrival->year);
                // no year in INGV boll!
                arrival->year = 2003;
                istat += ReadFortranInt(line, 5, 2, &arrival->month);
                istat += ReadFortranInt(line, 7, 2, &arrival->day);
                istat += ReadFortranInt(line, 18, 2, &arrival->hour);
                istat += ReadFortranInt(line, 33 + ioff, 2, &arrival->min);
                istat += ReadFortranReal(line, 35 + ioff, 4, &arrival->sec);
                arrival->sec /= 100.0;
                istat += ReadFortranReal(line, 22, 4, &psec);
                psec /= 100.0;
                // check for P second >= 60.0 (correct for INGV???)
                if (psec >= 60.0 && arrival->sec < 60.0)
                    arrival->sec += 60.0;
                istat += ReadFortranString(line, 27 + ioff, 1, arrival->onset);
                // set weight
                arrival->quality = arrival->onset[0] == 'i' ? 1 : 2;
                istat += ReadFortranString(line, 28 + ioff, 1, arrival->comp);
                if (arrival->comp[0] == ' ')
                    arrival->comp[0] = ARRIVAL_NULL_CHR;
                arrival->phase[0] = '\0';
                istat += ReadFortranString(line, 29 + ioff, 4, arrival->phase);
                TrimString(arrival->phase);
            }
        }

        /* check for S arrival input found */
        if (check_for_S_arrival > 0 && istat == 10
                //				&& IsPhaseID(arrival->phase, "S")
                && IsGoodDate(arrival->year,
                arrival->month, arrival->day)) {

            //snprintf(arrival->phase, sizeof(arrival->phase), "%s", "S");

            /* set error fields */
            snprintf(arrival->error_type, sizeof(arrival->error_type), "%s", "GAU");
            if (arrival->quality >= 0 &&
                    arrival->quality < NumQuality2ErrorLevels) {
                arrival->error =
                        Quality2Error[arrival->quality];
            } else {
                arrival->error =
                        Quality2Error[NumQuality2ErrorLevels - 1];
                nll_puterr("WARNING: invalid arrival weight.");
            }

            // increment additional arrival count
            if (strcmp(ftype_obs, "INGV_ARCH") == 0 && check_for_S_arrival < 3) {
                check_for_S_arrival++;
            } else {
                check_for_S_arrival = 0;
                line[0] = '\0';
            }
            return (istat);
        } else if (check_for_S_arrival > 0) {
            check_for_S_arrival = 0;
            return (OBS_FILE_SKIP_INPUT_LINE);
        }



        /*
                                        1 0
                                        DOI 1016   EZPG  03243563 EZSG  243882                             23
                                        11 0
                                        SEI 1016   EZPG  04132382 EZSG  132881
                                        VMG 1016   EZPG  04132514 EZSG  133068                            156
                                        ZCCA1016   EZPG  04132655                                          47
         */
        // read next line
        cstat = fgets(line, MAXLINE_LONG, fp_obs);

        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);

        /* read formatted P (or S) arrival input */
        istat = ReadFortranString(line, 1, 4, arrival->label);
        if (isspace(arrival->label[0])) {
            // check to skip next event(s) if is teleseism
            if (strcmp(ftype_obs, "INGV_BOLL_LOCAL") == 0) {
                while (1) {
                    TrimString(line);
                    if (strlen(line) > 10) {
                        while (1) { // read to next event
                            cstat = fgets(line, MAXLINE_LONG, fp_obs);
                            if (cstat == NULL)
                                return (OBS_FILE_END_OF_INPUT);
                            if (isspace(line[0]))
                                break;
                        }
                    } else {
                        break;
                    }
                }
            }
            line[0] = '\0';
            return (OBS_FILE_END_OF_EVENT);
        }
        TrimString(arrival->label);
        istat += ReadFortranString(line, 12, 1, arrival->onset);
        TrimString(arrival->onset);
        istat += ReadFortranString(line, 13, 1, arrival->comp);
        if (arrival->comp[0] == ' ')
            arrival->comp[0] = ARRIVAL_NULL_CHR;
        // set weight
        arrival->quality = arrival->onset[0] == 'i' ? 0 : 1;
        istat += ReadFortranString(line, 14, 4, arrival->phase);
        TrimString(arrival->phase);
        istat += ReadFortranString(line, 10, 1, arrival->first_mot);
        TrimString(arrival->first_mot);
        // no year in INGV boll!
        arrival->year = 2003;
        istat += ReadFortranInt(line, 5, 2, &arrival->month);
        istat += ReadFortranInt(line, 7, 2, &arrival->day);
        istat += ReadFortranInt(line, 18, 2, &arrival->hour);
        istat += ReadFortranInt(line, 20, 2, &arrival->min);
        istat += ReadFortranReal(line, 22, 4, &arrival->sec);
        arrival->sec /= 100.0;


        if (istat != 10) {
            line[0] = '\0';
            return (OBS_FILE_END_OF_EVENT);
        }

        /* read optional amplitude/period fields */
        istat += ReadFortranReal(line, 66, 4, &arrival->amplitude);
        //printf("AMPLITUDE: %lf (%s)\n", arrival->amplitude, line);
        //istat += ReadFortranReal(line, 48, 3, &arrival->period);
        // 20030826 AJL coda_dur added
        //istat += ReadFortranReal(line, 71, 5, &arrival->coda_dur);


        check_for_S_arrival = 1;
        /* check for valid phase code */
        //		if (IsPhaseID(arrival->phase, "P")) {
        //snprintf(arrival->phase, sizeof(arrival->phase), "%s", "P");
        //			check_for_S_arrival = 1;

        //		} else if (IsPhaseID(arrival->phase, "S")) {
        //snprintf(arrival->phase, sizeof(arrival->phase), "%s", "S");
        //			check_for_S_arrival = 0;
        //		} else
        //			return(OBS_FILE_END_OF_EVENT);

        if (!IsGoodDate(arrival->year, arrival->month, arrival->day))
            return (OBS_FILE_END_OF_EVENT);


        /* check for integer sec format */
        //if (arrival->sec > 99.999)
        //	arrival->sec /= 100.0;

        /* convert quality to error */
        Qual2Err(arrival);

        return (istat);
    } else if (strcmp(ftype_obs, "NCSN_Y2K_5") == 0) {

        if (check_for_S_arrival) {
            /* check for S phase input in last input line read */

            /* check for zero or blank S phase time */
            istat = ReadFortranString(line, 48, 1, chrtmp);
            if (istat > 0 && (strncmp(chrtmp, "S", 1) == 0 || strncmp(chrtmp, "s", 1) == 0)) {
                /* read S phase input in last input line read */
                istat = ReadFortranString(line, 1, 5, arrival->label);
                TrimString(arrival->label);
                // 20120928 AJL - added
                istat = ReadFortranString(line, 6, 2, arrival->network); // network code ignored in NLL
                TrimString(arrival->network);
                istat += ReadFortranString(line, 9, 1, arrival->comp);
                TrimString(arrival->comp);
                istat += ReadFortranString(line, 10, 3, arrival->inst);
                TrimString(arrival->inst);
                // END - 20120928 AJL - added
                istat += ReadFortranInt(line, 18, 4, &arrival->year);
                istat += ReadFortranInt(line, 22, 2, &arrival->month);
                istat += ReadFortranInt(line, 24, 2, &arrival->day);
                istat += ReadFortranInt(line, 26, 2, &arrival->hour);
                istat += ReadFortranInt(line, 28, 2, &arrival->min);
                istat += ReadFortranReal(line, 42, 5, &arrival->sec);
                /* check for integer sec format */
                //				if (arrival->sec > 99.999)
                arrival->sec /= 100.0;
                istat += ReadFortranReal(line, 30, 5, &psec);
                psec /= 100.0;
                /* check for P second >= 60.0 */
                if (psec >= 60.0 && arrival->sec < 60.0)
                    arrival->sec += 60.0;
                arrival->phase[0] = '\0';
                istat += ReadFortranString(line, 48, 1, arrival->phase);
                TrimString(arrival->phase);
                istat += ReadFortranString(line, 47, 1, arrival->onset);
                TrimString(arrival->onset);
                istat += ReadFortranInt(line, 50, 1, &arrival->quality);
            }
        }

        /* check for S arrival input found */
        if (check_for_S_arrival && istat == 11
                && IsPhaseID(arrival->phase, "S")
                && IsGoodDate(arrival->year, arrival->month, arrival->day)) {

            /* set error fields */
            snprintf(arrival->error_type, sizeof(arrival->error_type), "%s", "GAU");
            if (arrival->quality >= 0 && arrival->quality < NumQuality2ErrorLevels) {
                arrival->error = Quality2Error[arrival->quality];
            } else {
                arrival->error = Quality2Error[NumQuality2ErrorLevels - 1];
                nll_puterr("WARNING: invalid arrival weight.");
            }

            line[0] = '\0';
            check_for_S_arrival = 0;

            // AJL 20070608 - this format may have muliple entries for each phase, need to reject
            // earlier phases with large error so that later entries will be used
            if (arrival->error > 999.0)
                return (OBS_FILE_SKIP_INPUT_LINE);

            return (istat);
        } else if (check_for_S_arrival) {
            check_for_S_arrival = 0;
            return (OBS_FILE_SKIP_INPUT_LINE);
        }


        /* read next line */
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);
        if (strcmp(line, "    ") == 0) {
            line[0] = '\0';
            return (OBS_FILE_END_OF_EVENT);
        }

        /* read formatted P (or S) arrival input */
        istat = ReadFortranString(line, 1, 5, arrival->label);
        TrimString(arrival->label);
        // 20120928 AJL - added
        istat = ReadFortranString(line, 6, 2, arrival->network); // network code ignored in NLL
        //printf("DEBUG: arrival->network %s [%s]", arrival->network, line);
        TrimString(arrival->network);
        // END - 20120928 AJL - added
        istat += ReadFortranString(line, 9, 1, arrival->comp);
        TrimString(arrival->comp);
        istat += ReadFortranString(line, 10, 3, arrival->inst);
        TrimString(arrival->inst);
        istat += ReadFortranString(line, 14, 1, arrival->onset);
        TrimString(arrival->onset);
        if (strpbrk(arrival->onset, "XYZ") != NULL) { // skip RTP (phases with X,Y,Z onset)
            return (OBS_FILE_SKIP_INPUT_LINE);
        }
        istat += ReadFortranString(line, 15, 1, arrival->phase);
        TrimString(arrival->phase);
        istat += ReadFortranString(line, 16, 1, arrival->first_mot);
        TrimString(arrival->first_mot);
        istat += ReadFortranInt(line, 17, 1, &arrival->quality);
        /*
                        if (arrival->quality > 3) {			// skip phases with quality > 4
                        return(OBS_FILE_SKIP_INPUT_LINE);
        }
         */
        istat += ReadFortranInt(line, 18, 4, &arrival->year);
        istat += ReadFortranInt(line, 22, 2, &arrival->month);
        istat += ReadFortranInt(line, 24, 2, &arrival->day);

        if (!IsGoodDate(arrival->year, arrival->month, arrival->day))
            return (OBS_FILE_END_OF_EVENT);

        istat += ReadFortranInt(line, 26, 2, &arrival->hour);
        istat += ReadFortranInt(line, 28, 2, &arrival->min);
        istat += ReadFortranReal(line, 30, 5, &arrival->sec);
        /* check for integer sec format */
        //		if (arrival->sec > 99.999)
        arrival->sec /= 100.0;

        if (istat != 13) {
            line[0] = '\0';
            return (OBS_FILE_END_OF_EVENT);
        }

        /* read optional amplitude/period fields */
        istat += ReadFortranReal(line, 55, 7, &arrival->amplitude);
        arrival->amplitude /= 100.0;
        istat += ReadFortranReal(line, 84, 3, &arrival->period);
        arrival->period /= 100.0;


        /* check for valid phase code */
        if (IsPhaseID(arrival->phase, "P")) {
            //snprintf(arrival->phase, sizeof(arrival->phase), "%s", "P");
            check_for_S_arrival = 1;

        } else if (IsPhaseID(arrival->phase, "S")) {
            //snprintf(arrival->phase, sizeof(arrival->phase), "%s", "S");
            check_for_S_arrival = 0;
        } else
            return (OBS_FILE_SKIP_INPUT_LINE);

        /* convert quality to error */
        Qual2Err(arrival);

        // AJL 20070608 - this format may have muliple entries for each phase, need to reject
        // earlier phases with large error so that later entries will be used
        if (arrival->error > 999.0)
            return (OBS_FILE_SKIP_INPUT_LINE);

        return (istat);

    } else if (strcmp(ftype_obs, "NCSN_Y2K_5_SCSN") == 0) // 20251103 - support 9 char station labels from E Hauksson
    {
        /* read next line */
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);
        if (strcmp(line, "    ") == 0) {
            line[0] = '\0';
            return (OBS_FILE_END_OF_EVENT);
        }

        // check if event line and parse desired fields
        //printf("DEBUG: line %s\n", line);
        istat = ReadFortranString(line, 139, 8, phypo->public_id);
        if (istat > 0) { // success, must be event line (assume all other lines in file are < 139 char
            printf("DEBUG: phypo->public_id %s istat %d\n", phypo->public_id, istat);
            // add additional reads of NLL supported fields here, e.g.
            //istat += ReadFortranReal(line, 5, 4, &phypo->amp_mag); 	// arbitrary columns, not correct
            //printf("DEBUG: phypo->amp_mag %f istat %d\n", phypo->amp_mag, istat);
            return (OBS_FILE_SKIP_INPUT_LINE); // return so next phase file line will be read
        }

        // NCSN_Y2K_5_SCSN has S on separate line, always check
        // check for S phase input

        /* check for zero or blank S phase time */
        istat = ReadFortranString(line, 48, 1, chrtmp);
        //printf("DEBUG: istat %d chrtmp %s\n", istat, chrtmp);
        if (istat > 0 && (strncmp(chrtmp, "S", 1) == 0 || strncmp(chrtmp, "s", 1) == 0)) {
            /* read S phase input in last input line read */
            // 20251103 - support 9 char station labels from E Hauksson
            istat = ReadFortranString(line, 1, 9, arrival->label);
            TrimString(arrival->label);
            istat += 2; // skipping network and comp
            istat += ReadFortranString(line, 10, 3, arrival->inst);
            TrimString(arrival->inst);
            // END - 20120928 AJL - added
            istat += ReadFortranInt(line, 18, 4, &arrival->year);
            istat += ReadFortranInt(line, 22, 2, &arrival->month);
            istat += ReadFortranInt(line, 24, 2, &arrival->day);
            istat += ReadFortranInt(line, 26, 2, &arrival->hour);
            istat += ReadFortranInt(line, 28, 2, &arrival->min);
            istat += ReadFortranReal(line, 42, 5, &arrival->sec);
            /* check for integer sec format */
            //				if (arrival->sec > 99.999)
            arrival->sec /= 100.0;
            istat += ReadFortranReal(line, 30, 5, &psec);
            psec /= 100.0;
            /* check for P second >= 60.0 */
            if (psec >= 60.0 && arrival->sec < 60.0)
                arrival->sec += 60.0;
            arrival->phase[0] = '\0';
            istat += ReadFortranString(line, 48, 1, arrival->phase);
            TrimString(arrival->phase);
            istat += ReadFortranString(line, 47, 1, arrival->onset);
            TrimString(arrival->onset);
            istat += ReadFortranInt(line, 50, 1, &arrival->quality);
            //printf("DEBUG: istat %d  arrival->phase %s\n", istat, arrival->phase);
        }

        /* check for S arrival input found */
        if (istat >= 12
                && IsPhaseID(arrival->phase, "S")
                && IsGoodDate(arrival->year, arrival->month, arrival->day)) {

            /* set error fields */
            snprintf(arrival->error_type, sizeof(arrival->error_type), "%s", "GAU");
            if (arrival->quality >= 0 && arrival->quality < NumQuality2ErrorLevels) {
                arrival->error = Quality2Error[arrival->quality];
            } else {
                arrival->error = Quality2Error[NumQuality2ErrorLevels - 1];
                nll_puterr("WARNING: invalid arrival weight.");
            }
            //printf("DEBUG: arrival->quality %d  arrival->error %f\n", arrival->quality, arrival->error);

            line[0] = '\0';

            // AJL 20070608 - this format may have multiple entries for each phase, need to reject
            // earlier phases with large error so that later entries will be used
            if (arrival->error > 999.0)
                return (OBS_FILE_SKIP_INPUT_LINE);

            return (istat);
        }

        /* read formatted P (or S) arrival input */
        // 20251103 - support 9 char station labels from E Hauksson
        istat = ReadFortranString(line, 1, 9, arrival->label);
        TrimString(arrival->label);
        istat += 2; // skipping network and comp
        istat += ReadFortranString(line, 10, 3, arrival->inst);
        TrimString(arrival->inst);
        istat += ReadFortranString(line, 14, 1, arrival->onset);
        TrimString(arrival->onset);
        if (strpbrk(arrival->onset, "XYZ") != NULL) { // skip RTP (phases with X,Y,Z onset)
            return (OBS_FILE_SKIP_INPUT_LINE);
        }
        istat += ReadFortranString(line, 15, 1, arrival->phase);
        TrimString(arrival->phase);
        istat += ReadFortranString(line, 16, 1, arrival->first_mot);
        TrimString(arrival->first_mot);
        istat += ReadFortranInt(line, 17, 1, &arrival->quality);
        /*
                        if (arrival->quality > 3) {			// skip phases with quality > 4
                        return(OBS_FILE_SKIP_INPUT_LINE);
        }
         */
        istat += ReadFortranInt(line, 18, 4, &arrival->year);
        istat += ReadFortranInt(line, 22, 2, &arrival->month);
        istat += ReadFortranInt(line, 24, 2, &arrival->day);

        if (!IsGoodDate(arrival->year, arrival->month, arrival->day))
            return (OBS_FILE_END_OF_EVENT);

        istat += ReadFortranInt(line, 26, 2, &arrival->hour);
        istat += ReadFortranInt(line, 28, 2, &arrival->min);
        istat += ReadFortranReal(line, 30, 5, &arrival->sec);
        /* check for integer sec format */
        //		if (arrival->sec > 99.999)
        arrival->sec /= 100.0;

        if (istat != 14) {
            line[0] = '\0';
            return (OBS_FILE_END_OF_EVENT);
        }

        /* read optional amplitude/period fields */
        istat += ReadFortranReal(line, 55, 7, &arrival->amplitude);
        arrival->amplitude /= 100.0;
        istat += ReadFortranReal(line, 84, 3, &arrival->period);
        arrival->period /= 100.0;


        /* check for valid phase code */
        if (IsPhaseID(arrival->phase, "P")) {
            //snprintf(arrival->phase, sizeof(arrival->phase), "%s", "P");
            ;

        } else if (IsPhaseID(arrival->phase, "S")) {
            //snprintf(arrival->phase, sizeof(arrival->phase), "%s", "S");
            ;
        } else
            return (OBS_FILE_SKIP_INPUT_LINE);

        /* convert quality to error */
        Qual2Err(arrival);

        // AJL 20070608 - this format may have muliple entries for each phase, need to reject
        // earlier phases with large error so that later entries will be used
        if (arrival->error > 999.0)
            return (OBS_FILE_SKIP_INPUT_LINE);

        return (istat);


    } else if (strcmp(ftype_obs, "PAOLO_OV") == 0 ||
            strcmp(ftype_obs, "ALBERTO_3D") == 0 ||
            strcmp(ftype_obs, "ALBERTO_3D_4") == 0 ||
            strcmp(ftype_obs, "SIMULPS") == 0 ||
            strcmp(ftype_obs, "VELEST") == 0) {
        /*
                        "PAOLO_OV"
                        87 119 23 5 45.00 40 49.17  14 25.80   0.43   0.00
                        BKEP0 45.51 BKWP0 45.61 BKNP0 45.66 BAFP0 45.74 BKES1 45.86 BKWS1 46.04
                        BKNS1 46.06 BAFS1 46.20

                        "ALBERTO_3D"
                        891018  0 4 15.34 37  2.66 121 52.54  16.92   5.80
                        JTGVIPU018.3000JBZVIPU018.7600JECVIPD018.5700JRGVIPU018.4500JHLVIPD018.5500
                        HLTVIPU024.9600BSLVIPU025.5800BPCVEPD025.4900HJSVIPU025.4900

                        "SIMULPS"
                        from Stephan Husen   email:  stephan@tomo.ig.erdw.ethz.ch
                        summary line containing event origin time, location, and magnitude
                        format: a4,a2,1x,a2,i2,1x,f5.2,i3,a1,f5.2,1x,i3,a1,f5.2,f7.2,f7.2)
                        iyrmo, iday, ihr, mino, seco, ltde, cns, eltm, lnde, cew, elnm, dep, mag
                        observed traveltimes to stations, six sets per line
                        each set is in the format: a4,a4,f6.2       sta, rmki, tt        rmki -> (I,E)(P,S,SP)(U,D,-)(0-4)
                        blank terminates event's data
        example:
                        84 3 6  6 7 58.58 46N19.78   7E26.85   3.77   2.00
                        SIE_IP-1  1.45DIX2IP-1  4.94DIX_IP-1  4.95EMS_IP-1  8.14EMV_IP-1  8.54MMK2IP-1  8.58
                        MMK_IP-1  8.67STG_IP-1 16.48SLE_IP-1 27.54WIL_IP-1 28.14SAX_IP-1 30.78

                        84 3 7  4 7 56.92 46N20.89   7E26.59   6.09   1.50
                        ZZB_IP-1  1.03ZZE_IP-1  1.06ZZA_IP-1  1.07ZZC_IP-1  1.11ZZD_IP-1  1.11ZZF_IP-1  1.11

                        "VELEST"
                        phase format of VELEST program; from Stephan Husen email: stephan@tomo.ig.erdw.ethz.ch
                        summary line containing event origin time, location, and magnitude
                        format:a4,a2,1x,a2,i2,1x,f5.2,1x,f7.4,a1,1x,f8.4,a1,f7.2,f7.2
                        iyrmo, iday, ihr, mino, seco, ltde, cns, lnde, cew, dep, mag
                        observed travel times (six sets per line)
                        each set is in format a4,a1,a1,f6.2
                        sta, phaseID, weight, tt

         */
        /* check for end of event (assumes no blanks after last phase) */
        chr = fgetc(fp_obs);
        if (chr == EOF) {
            return (OBS_FILE_END_OF_INPUT);
        } else if (chr == '\n') {
            if ((chr = fgetc(fp_obs)) != EOF && !isdigit(chr)) {
                // another phase follows
                ungetc(chr, fp_obs);
            } else if (chr == '0') {
                // end of event "0" line
                // read to end of line
                while ((chr = fgetc(fp_obs)) != EOF && chr != '\n')
                    ;
                return (OBS_FILE_END_OF_EVENT);
            } else {
                // end of event
                ungetc(chr, fp_obs);
                return (OBS_FILE_END_OF_EVENT);
            }
        } else {
            ungetc(chr, fp_obs);
        }


        /* check for event hypocenter line */
        if ((chr = fgetc(fp_obs)) != EOF && isdigit(chr)) {
            // read hypocenter time
            ungetc(chr, fp_obs);
            /* read next line */
            cstat = fgets(line, MAXLINE_LONG, fp_obs);
            if (cstat == NULL)
                return (OBS_FILE_END_OF_INPUT);
            istat = ReadFortranInt(line, 1, 2, &EventTime.year);
            istat += ReadFortranInt(line, 3, 2, &EventTime.month);
            istat += ReadFortranInt(line, 5, 2, &EventTime.day);
            istat += ReadFortranInt(line, 8, 2, &EventTime.hour);
            istat += ReadFortranInt(line, 10, 2, &EventTime.min);
            if (istat == 5) {
                // temporary fix for HYPO71 Y2K problem
                // 20231017 //if (arrival->year < 20)
                if (arrival->year < 40)
                    EventTime.year += 100;
                EventTime.year += 1900;
                // read to end of line
                //while ((chr = fgetc(fp_obs)) != EOF && chr != '\n')
                //	;
                if (chr == EOF)
                    return (OBS_FILE_END_OF_INPUT);
            }
        } else {
            if (chr == EOF)
                return (OBS_FILE_END_OF_INPUT);
            ungetc(chr, fp_obs);
        }


        if (strcmp(ftype_obs, "PAOLO_OV") == 0) {
            // try to read phase arrival input
            istat = fscanf(fp_obs, " %s %lf", line, &arrival->sec);
            if (istat == EOF)
                return (OBS_FILE_END_OF_INPUT);
            istat += ReadFortranString(line, 1, 3, arrival->label);
            istat += ReadFortranString(line, 4, 1, arrival->phase);
            istat += ReadFortranInt(line, 5, 1, &arrival->quality);
            if (istat != 5)
                return (OBS_FILE_END_OF_EVENT);

        } else if (strcmp(ftype_obs, "SIMULPS") == 0) {
            // try to read phase arrival input
            istat = fscanf(fp_obs, "%14c", line);
            line[14] = '\0';
            if (istat == EOF)
                return (OBS_FILE_END_OF_INPUT);
            // printf("<%s>\n", line);
            istat += ReadFortranString(line, 1, 4, arrival->label);
            istat += ReadFortranString(line, 5, 1, arrival->onset);
            istat += ReadFortranString(line, 6, 1, arrival->phase);
            istat += ReadFortranString(line, 7, 1, arrival->first_mot);
            istat += ReadFortranInt(line, 8, 1, &arrival->quality);
            istat += ReadFortranReal(line, 9, 6, &arrival->sec);

            if (istat != 7)
                return (OBS_FILE_END_OF_EVENT);

        } else if (strcmp(ftype_obs, "VELEST") == 0) {
            /* try to read phase arrival input */
            istat = fscanf(fp_obs, "%12c", line);
            line[12] = '\0';
            if (istat == EOF)
                return (OBS_FILE_END_OF_INPUT);
            /* printf("<%s>\n", line); */
            istat += ReadFortranString(line, 1, 4, arrival->label);
            istat += ReadFortranString(line, 5, 1, arrival->phase);
            istat += ReadFortranInt(line, 6, 1, &arrival->quality);
            istat += ReadFortranReal(line, 7, 6, &arrival->sec);

            if (istat != 5)
                return (OBS_FILE_END_OF_EVENT);

        } else if (strcmp(ftype_obs, "ALBERTO_3D") == 0
                || strcmp(ftype_obs, "ALBERTO_3D_4") == 0) {
            // try to read phase arrival input
            istat = fscanf(fp_obs, "%15c", line);
            line[15] = '\0';
            if (istat == EOF)
                return (OBS_FILE_END_OF_INPUT);
            //printf("<%s>\n", line);
            if (strcmp(ftype_obs, "ALBERTO_3D_4") == 0) {
                istat += ReadFortranString(line, 1, 4, arrival->label);
                istat++;
            } else {
                istat += ReadFortranString(line, 1, 3, arrival->label);
                istat += ReadFortranString(line, 4, 1, arrival->comp);
            }
            istat += ReadFortranString(line, 5, 1, arrival->onset);
            istat += ReadFortranString(line, 6, 1, arrival->phase);
            istat += ReadFortranString(line, 7, 1, arrival->first_mot);
            istat += ReadFortranInt(line, 8, 1, &arrival->quality);
            istat += ReadFortranReal(line, 9, 7, &arrival->sec);

            if (istat != 8)
                return (OBS_FILE_END_OF_EVENT);
        }


        arrival->year = EventTime.year;
        arrival->month = EventTime.month;
        arrival->day = EventTime.day;
        arrival->hour = EventTime.hour;
        arrival->min = EventTime.min;


        /* convert quality to error */
        Qual2Err(arrival);

        return (istat);

    } else if (strcmp(ftype_obs, "PAL") == 0) {
        // Zhou, Y., Yue, H., Fang, L., Zhou, S., Zhao, L. & Ghosh, A., 2021. An Earthquake Detection and Location Architecture for Continuous Seismograms: Phase Picking, Association, Location, and Matched Filter (PALM). Seismological Research Letters. doi:10.1785/0220210111

        /*
2023-02-01T07:02:02.724000Z,38.28,38.68,23,2.37,0.5
KO.SVRC,2023-02-01T07:02:12.920000Z,2023-02-01T07:02:21.690000Z,4.4354171239049296e-07,42808.9
TU.ATAB,2023-02-01T07:02:20.190000Z,2023-02-01T07:02:33.830000Z,2.2707865061916482e-07,36.4
KO.DARE,2023-02-01T07:02:21.750000Z,2023-02-01T07:02:36.290000Z,6.171520362434723e-08,28.4
TU.NARI,2023-02-01T07:02:11.210000Z,2023-02-01T07:02:17.650000Z,1.2742076995164313e-07,16.7
TU.AKCD,2023-02-01T07:02:15.210000Z,2023-02-01T07:02:24.530000Z,7.080782126560284e-07,25167.6
TU.HEKM,2023-02-01T07:02:19.710000Z,2023-02-01T07:02:32.270000Z,6.990056816112568e-07,37.8
2023-02-01T03:38:52.499600Z,38.38,39.36,5,0.96,0.3
TU.MKAM,2023-02-01T03:39:01.900000Z,2023-02-01T03:39:09.390000Z,8.320017592216995e-09,49.8
KO.SVRC,2023-02-01T03:38:53.860000Z,2023-02-01T03:38:54.960000Z,1.1081106036964678e-06,107958.4
         */

        char *cptr;

        //*DEBUG*/printf("TP 00");
        if (check_for_S_arrival) {
            check_for_S_arrival = 0;
            // check for S phase input in last input line read
            // KO.SVRC,2023-02-01T07:02:12.920000Z,2023-02-01T07:02:21.690000Z,4.4354171239049296e-07,42808.9
            // replace ',' with ' '
            while ((cptr = strchr(line, ',')) != NULL) {
                *cptr = ' ';
            }
            istat = sscanf(line, "%s %*d-%*d-%*dT%*d:%*d:%*fZ %d-%d-%dT%d:%d:%lfZ",
                    arrival->label, &arrival->year, &arrival->month, &arrival->day, &arrival->hour, &arrival->min, &arrival->sec
                    );
            //*DEBUG*/printf("\n%d  %s\n", istat, line);
            if (istat != 7) {
                return (OBS_FILE_SKIP_INPUT_LINE);
            }
            // convert periods in NET.STA to _
            while ((cptr = strchr(arrival->label, '.')) != NULL) {
                *cptr = '_';
            }
            strncpy(arrival->phase, "S", 2);
            arrival->quality = 1;
            // convert quality to error
            Qual2Err(arrival);

            return (istat);
        }

        //*DEBUG*/printf(" 01");
        // check for end of event or event hypocenter line (assumes no blanks after last phase)
        chr = fgetc(fp_obs);
        //*DEBUG*/printf("chr %c %d %d\n", chr, isdigit(chr), in_hypocenter_event);
        ungetc(chr, fp_obs);
        if (chr != EOF && isdigit(chr)) {
            if (in_hypocenter_event) {
                // end of event
                in_hypocenter_event = 0;
                return (OBS_FILE_END_OF_EVENT);
            } else {
                // read hypocenter line
                cstat = fgets(line, MAXLINE_LONG, fp_obs);
                if (cstat == NULL)
                    return (OBS_FILE_END_OF_INPUT);

                // 2023-02-01T07:02:02.724000Z,38.28,38.68,23,2.37,0.5
                istat = sscanf(line, "%d-%d-%dT%d:%d:%lfZ,%lf,%lf,%lf,%lf",
                        &EventTime.year, &EventTime.month, &EventTime.day,
                        &EventTime.hour, &EventTime.min, &EventTime.sec,
                        &phypo->dlat, &phypo->dlong, &phypo->depth, &phypo->amp_mag
                        );
                if (istat != 10) {
                    return (OBS_FILE_END_OF_EVENT);
                }
                in_hypocenter_event = 1;
                check_for_S_arrival = 0;
            }
        } else {
            if (chr == EOF)
                return (OBS_FILE_END_OF_INPUT);
        }

        //*DEBUG*/printf(" 03");
        // P phase
        // read next line
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);
        //*DEBUG*/printf(" 04");

        // KO.SVRC,2023-02-01T07:02:12.920000Z,2023-02-01T07:02:21.690000Z,4.4354171239049296e-07,42808.9
        // replace ',' with ' '
        while ((cptr = strchr(line, ',')) != NULL) {
            *cptr = ' ';
        }
        istat = sscanf(line, "%s %d-%d-%dT%d:%d:%lfZ",
                arrival->label, &arrival->year, &arrival->month, &arrival->day, &arrival->hour, &arrival->min, &arrival->sec
                );
        //*DEBUG*/printf(" 05");
        //*DEBUG*/printf("\n%d  %s\n", istat, line);
        if (istat != 7) {
            return (OBS_FILE_END_OF_EVENT);
        }
        //*DEBUG*/printf(" 06");
        // convert periods in NET.STA to _
        while ((cptr = strchr(arrival->label, '.')) != NULL) {
            *cptr = '_';
        }
        strncpy(arrival->phase, "P", 2);
        arrival->quality = 0;
        // convert quality to error
        Qual2Err(arrival);

        check_for_S_arrival = 1;

        //*DEBUG*/printf(" 99\n");
        return (istat);

    } else if (strcmp(ftype_obs, "SED_LOC") == 0 || strcmp(ftype_obs, "SED_LOC_ERR") == 0) {

        /* example:
                        SH 06/11/2005 LOC format has been extended to include uncertainty interval for picking
                        SH 23/11/2005 LOC format has been extended to include observation weight (quality) in last column

                        1 1   200.0   300.0    1.71 1                        INST
                        0.000     0.000  10.00 0000 00 00 00 00  0.00 0    TRIAL
                        nico                                                  AUTOR
                        2002/01/04 00:01                                      Local            20020015
                        FUSIO   P       ID  27.454 1 0.21      186 NWA       0 Pick HHN  29.626  -0.050   0.050 0
                        FUSIO   S       E   29.474 1 0.21      186 NWA       0 Pick HHN  29.626  -0.250   0.250 1
                        VDL     P       ID  33.578 1 0.24       16 NWA       0 Pick HHN  41.014  -0.050   0.050 1
                        VDL     S       Q   39.651 1 0.24       16 NWA       0 Pick HHN  41.014  -0.500   0.500 3
                        LLS     P       IU  34.752 1 0.36       14 NWA       0 Pick HHE  43.193  -0.050   0.050 2
                        LLS     S       I   41.622 1 0.36       14 NWA       0 Pick HHE  43.193  -0.050   0.050 1
                        BNALP   P       Q   35.959 1 0.59        8 NWA       0 Pick HHE  47.417  -0.500   0.500 3
                        BNALP   S       E   43.504 1 0.59        8 NWA       0 Pick HHE  47.417  -0.250   0.250 2
                        HASLI   P       E   36.250 1 0.36       32 NWA       0 Pick HHE  44.465  -0.250   0.250 2
                        HASLI   S       E   44.209 1 0.36       32 NWA       0 Pick HHE  44.465  -0.250   0.250 2
                        MMK     P       E   36.365 1 0.22       32 NWA       0 Pick HHN  46.943  -0.250   0.250 2
                        LKBD    P       E   39.858 1 1.77       18 NWA       0 Pick HHN  53.318  -0.250   0.250 2
                        BERNI   P       E   40.349 1 0.58       10 NWA       0 Pick HHN  56.013  -0.250   0.250 2
                        CHDAW   P       E   40.741 1 0.29       14 NWA       0 Pick HHN  59.776  -0.250   0.250 2
                        FUORN   P       IU  44.069 1 0.25       10 NWA       0 Pick HHE  65.863  -0.050   0.050 0
                        SKIP
                        20020104000124246358N008807E00515Ml1272705135008005022SEDN023156015414         A
                        KP200201040001                                        ARCHIVE
         */

        // read line
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);

        // check for end of event (assumes empty lines between event)
        if (LineIsBlank(line)) {
            // end of event
            return (OBS_FILE_END_OF_EVENT);
        }

        // 20100607 AJL - initialize to garbage value to prevent valgrind uninitialized error
        snprintf(eth_line_key, sizeof(eth_line_key), "%s", "$@GARBAGE");

        // read event hypocenter line
        if (!in_hypocenter_event) {
            // find origin time line
            ifound = 0;
            while ((istat = ReadFortranString(line, 55, 4, eth_line_key)) > 0) {
                /* SH 03/05/2003
                                read control parameters, which are specified in first line of SED_LOC format
                                so far only VpVSRatio is of use
                                the value of VpVsRatio specified in NLLoc control file will be overwritten!
                                SH 11/25/2004
                                VpVsRatio should only be read from LOCfile if NLLoc is used within SNAP */
#ifdef CUSTOM_ETH
                if (strcmp(eth_line_key, "INST") == 0) {
                    istat = ReadFortranReal(line, 23, 6, &vpvs);
                    if (istat != 1) {
                        nll_puterr(
                                "WARNING: could not read VpVsRatio! Use value of control file\n");
                    } else {
                        VpVsRatio = vpvs;
                    }
                }
#endif
                /* SH 07/26/2004
                                                other identiers for SED_LOC are regi, Tele and Unkn  */
                if ((strcmp(eth_line_key, "Loca") == 0) ||
                        (strcmp(eth_line_key, "Regi") == 0) ||
                        (strcmp(eth_line_key, "Tele") == 0) ||
                        (strcmp(eth_line_key, "Unkn") == 0)) {
                    ifound = 1;
                    break;
                }
                // read next line
                cstat = fgets(line, MAXLINE_LONG, fp_obs);
                if (cstat == NULL)
                    return (OBS_FILE_END_OF_INPUT);
            }
            if (!ifound)
                return (OBS_FILE_END_OF_EVENT);

            // read hypocenter time
            //2001/12/04 01:29                                      Local           20010748
            istat = ReadFortranInt(line, 1, 4, &EventTime.year);
            istat += ReadFortranInt(line, 6, 2, &EventTime.month);
            istat += ReadFortranInt(line, 9, 2, &EventTime.day);
            istat += ReadFortranInt(line, 12, 2, &EventTime.hour);
            istat += ReadFortranInt(line, 15, 2, &EventTime.min);
            /* SH 03/05/2004  added event number  */
            istat += ReadFortranInt(line, 71, 8, &EventTime.ev_nr);
            /* SH 23/11/2004  event number is undertermined when reading GSE2 files */
            if (istat == 4 || istat == 6) {
                if (istat == 4)
                    EventTime.ev_nr = 0; /* no event number */
                in_hypocenter_event = 1;
                // read until phase line reached
                while ((istat = ReadFortranString(line, 55, 4, eth_line_key)) > 0
                        && strcmp(line, "    ") != 0
                        /* SH 06/11/2005  change to deal with extended LOC format */
                        && strcmp(eth_line_key, " Pic") != 0) {
                    cstat = fgets(line, MAXLINE_LONG, fp_obs);
                    if (cstat == NULL)
                        return (OBS_FILE_END_OF_INPUT);
                }
            } else {
                return (OBS_FILE_END_OF_EVENT);
            }
        }

        // check if valid phase line: key is blank or strlen < 55
        if ((istat = ReadFortranString(line, 55, 4, eth_line_key)) > 0
                && strcmp(line, "    ") != 0
                /* SH 06/11/2005  change to deal with extended LOC format */
                && strcmp(eth_line_key, " Pic") != 0) {
            return (OBS_FILE_SKIP_INPUT_LINE);
        }

        /* read phase arrival input
                           SH 06/11/2005 LOC format has been extended to include uncertainty interval for picking
                           SH 23/11/2005 LOC format has been extended to include obs weight in the last coluumn
                           FUSIO   P       ID  27.454 1 0.21      186 NWA       0 Pick HHN  29.626  -0.050   0.050 1
                           FUSIO   S       E   29.474 1 0.21      186 NWA       0 Pick HHN  29.626  -0.250   0.250 3
         */
        istat = ReadFortranString(line, 1, 5, arrival->label);
        istat += ReadFortranString(line, 9, 6, arrival->phase);
        istat += ReadFortranString(line, 17, 1, arrival->onset);
        istat += ReadFortranString(line, 18, 1, arrival->first_mot);
        istat += ReadFortranReal(line, 20, 7, &arrival->sec);
        //		istat += ReadFortranReal(line, 21, 6, &arrival->sec);
        istat += ReadFortranInt(line, 28, 1, &eth_use_loc);
        istat += ReadFortranReal(line, 29, 5, &arrival->period);
        istat += ReadFortranReal(line, 34, 9, &arrival->amplitude);
        /* SH 02/25/2004 bug fix
                                           arrival->inst is only of type char[5] and not char[9]!
                                           istat += ReadFortranString(line, 44,9, arrival->inst);  */
        istat += ReadFortranString(line, 44, 4, arrival->inst);
        if (strcmp(arrival->inst, "    ") == 0)
            snprintf(arrival->inst, sizeof(arrival->inst), "%s", ARRIVAL_NULL_STR);
        istat += ReadFortranInt(line, 54, 1, &arrival->clipped);
        /* SH 23/11/2005 read arrival quality from last column (for extended LOC format) */
        if (strcmp(eth_line_key, " Pic") == 0)
            if ((ReadFortranInt(line, 89, 1, &arrival->quality)) > 0) istat++;


        if (istat < 10) {
            return (OBS_FILE_END_OF_EVENT);
        }

        /* remove blanks/whitespace */
        removeSpace(arrival->label);
        removeSpace(arrival->phase);
        removeSpace(arrival->inst);

        /* AJL EvalPhaseID not used here anymore, used later */
        /*		if (EvalPhaseID(arrival->phase) < 0) {
                                           nll_puterr2("WARNING: phase ID not found", arrival->phase);
                                           return(OBS_FILE_INVALID_PHASE);
                }
         */

        arrival->year = EventTime.year;
        arrival->month = EventTime.month;
        arrival->day = EventTime.day;
        arrival->hour = EventTime.hour;
        arrival->min = EventTime.min;

        /* convert onset to error */
        //Uncertainties of phase readings:
        //                              I             E            Q
        //for Key 'Loca':            < +/- 0.05   < +/- 0.2     > +/- 0.2
        //for Key 'Regi' or 'Tele'   < +/- 0.20   < +/- 1.0     > +/- 1.0

        /* SH 23/11/2005 for extended LOC format arrival quality is now in the last column;
                                           arrival onset should no longer be used to estimate arrival error */

        // AJL 20070308 added SED_LOC_ERR
        if (strcmp(ftype_obs, "SED_LOC_ERR") == 0) {
            istat = ReadFortranReal(line, 73, 7, &left_uncertainty);
            istat += ReadFortranReal(line, 81, 7, &right_uncertainty);
            snprintf(arrival->error_type, sizeof(arrival->error_type), "%s", "GAU");
            arrival->error = (right_uncertainty - left_uncertainty) / 2.0;
            ;
        } else if (istat == 11) { /* extended LOC format */
            if (eth_use_loc == 0) {
                arrival->error = ARRIVAL_ERROR_NULL;
            } else {
                Qual2Err(arrival);
            } /* convert quality to uncertainty in s */
        } else { /* old LOC format */
            snprintf(arrival->error_type, sizeof(arrival->error_type), "%s", "GAU");
            if (eth_use_loc == 0) {
                arrival->error = ARRIVAL_ERROR_NULL;
            } else if (strcmp(arrival->onset, "I") == 0) {
                arrival->error = 0.05;
            } else if (strcmp(arrival->onset, "E") == 0) {
                arrival->error = 0.2;
            } else if (strcmp(arrival->onset, "Q") == 0) {
                arrival->error = 0.2;
            } else {
                arrival->error = ARRIVAL_ERROR_NULL;
            }
            arrival->quality = 0;
        }

        return (istat);
    } else if (strcmp(ftype_obs, "SED_LOC_OLD") == 0) {

        /* example:
                        1 1   200.0   300.0    1.71 1                        INST
                        0.000     0.000  10.00 0000 00 00 00 00  0.00 0    TRIAL
                        maraini                                               AUTOR
                        2001/12/04 01:29                                      Local           20010748
                        FUORN   P       E   14.910 1 0.14      448 NWA       0
                        FUORN   S       E   16.140 1 0.14      448 NWA       0
                        OSS2    P       ID  16.428 1 0.09     2554 AHP       0
                        OSS     P       ED  16.524 1 0.11     1672 AHP       0
                        BERNI   P       E   17.930 1 0.18       46 NWA       0
                        BERNI   S       E   21.139 1 0.18       46 NWA       0
                        CHDAW   P       Q   20.016 1 0.56       10 NWA       0
                        VDL2    P       IU  23.704 1 0.20      246 AHP       0
                        VDL     P       E   23.863 1 0.12      188 AHP       0
                        VDL     P       E   23.822 1 0.28        4 NWA       1
                        VDL     S       E   31.177 1 0.28        4 NWA       1
                        SKIP
                        20011204012913146568N010226E00614Ml1239814161001001000SEDG024196006511         B
                        KP200112040128                                        ARCHIVE
         */

        // read line
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);

        // check for end of event (assumes empty lines between event)
        if (LineIsBlank(line)) {
            // end of event
            return (OBS_FILE_END_OF_EVENT);
        }


        // read event hypocenter line
        if (!in_hypocenter_event) {
            // find origin time line
            ifound = 0;
            while ((istat = ReadFortranString(line, 55, 4, eth_line_key)) > 0) {
                /* SH 03/05/2003
                                read control parameters, which are specified in first line of SED_LOC format
                                so far only VpVSRatio is of use
                                the value of VpVsRatio specified in NLLoc control file will be overwritten!  */
                if (strcmp(eth_line_key, "INST") == 0) {
                    istat = ReadFortranReal(line, 23, 6, &vpvs);
                    if (istat != 1) {
                        nll_puterr(
                                "WARNING: could not read VpVsRatio! Use value of control file\n");
                    } else {
                        VpVsRatio = vpvs;
                    }
                }
                if (strcmp(eth_line_key, "Loca") == 0) {
                    ifound = 1;
                    break;
                }
                // read next line
                cstat = fgets(line, MAXLINE_LONG, fp_obs);
                if (cstat == NULL)
                    return (OBS_FILE_END_OF_INPUT);
            }
            if (!ifound)
                return (OBS_FILE_END_OF_EVENT);

            // read hypocenter time
            //2001/12/04 01:29                                      Local           20010748
            istat = ReadFortranInt(line, 1, 4, &EventTime.year);
            istat += ReadFortranInt(line, 6, 2, &EventTime.month);
            istat += ReadFortranInt(line, 9, 2, &EventTime.day);
            istat += ReadFortranInt(line, 12, 2, &EventTime.hour);
            istat += ReadFortranInt(line, 15, 2, &EventTime.min);
            /* SH 03/05/2004  added event number  */
            istat += ReadFortranInt(line, 71, 8, &EventTime.ev_nr);
            if (istat == 6) {
                in_hypocenter_event = 1;
                // read until phase line reached
                while ((istat = ReadFortranString(line, 55, 4, eth_line_key)) > 0
                        && strcmp(line, "    ") != 0) {
                    cstat = fgets(line, MAXLINE_LONG, fp_obs);
                    if (cstat == NULL)
                        return (OBS_FILE_END_OF_INPUT);
                }
            } else {
                return (OBS_FILE_END_OF_EVENT);
            }
        }

        // check if valid phase line: key is blank or strlen < 55
        if ((istat = ReadFortranString(line, 55, 4, eth_line_key)) > 0
                && strcmp(line, "    ") != 0) {
            return (OBS_FILE_SKIP_INPUT_LINE);
        }

        // read phase arrival input
        //CHDAW   P       Q   20.016 1 0.56       10 NWA       0
        //OSS     P       ED  16.524 1 0.11     1672 AHP       0
        //BERNI   P       E   17.930 1 0.18       46 NWA       0

        istat = ReadFortranString(line, 1, 5, arrival->label);
        istat += ReadFortranString(line, 9, 6, arrival->phase);
        istat += ReadFortranString(line, 17, 1, arrival->onset);
        istat += ReadFortranString(line, 18, 1, arrival->first_mot);
        istat += ReadFortranReal(line, 20, 7, &arrival->sec);
        //		istat += ReadFortranReal(line, 21, 6, &arrival->sec);
        istat += ReadFortranInt(line, 28, 1, &eth_use_loc);
        istat += ReadFortranReal(line, 29, 5, &arrival->period);
        istat += ReadFortranReal(line, 34, 9, &arrival->amplitude);
        /* SH 02/25/2004 bug fix
                                           arrival->inst is only of type char[5] and not char[9]!
                                           istat += ReadFortranString(line, 44,9, arrival->inst);  */
        istat += ReadFortranString(line, 44, 4, arrival->inst);
        /* SH 07232004 added */
        //istat += ReadFortranInt(line, 54, 1, &eth_use_mag);
        istat += ReadFortranInt(line, 54, 1, &arrival->clipped);


        if (istat != 10) {
            return (OBS_FILE_END_OF_EVENT);
        }

        /* remove blanks/whitespace */
        removeSpace(arrival->label);
        removeSpace(arrival->phase);
        removeSpace(arrival->inst);

        /* AJL EvalPhaseID not used here anymore, used later */
        /*		if (EvalPhaseID(arrival->phase) < 0) {
                                           nll_puterr2("WARNING: phase ID not found", arrival->phase);
                                           return(OBS_FILE_INVALID_PHASE);
                }
         */

        arrival->year = EventTime.year;
        arrival->month = EventTime.month;
        arrival->day = EventTime.day;
        arrival->hour = EventTime.hour;
        arrival->min = EventTime.min;

        /* convert onset to error */
        //Uncertainties of phase readings:
        //                              I             E            Q
        //for Key 'Loca':            < +/- 0.05   < +/- 0.2     > +/- 0.2
        //for Key 'Regi' or 'Tele'   < +/- 0.20   < +/- 1.0     > +/- 1.0
        snprintf(arrival->error_type, sizeof(arrival->error_type), "%s", "GAU");
        if (eth_use_loc == 0) {
            arrival->error = ARRIVAL_ERROR_NULL;
        } else if (strcmp(arrival->onset, "I") == 0) {
            arrival->error = 0.05;
        } else if (strcmp(arrival->onset, "E") == 0) {
            arrival->error = 0.2;
        } else if (strcmp(arrival->onset, "Q") == 0) {
            arrival->error = 0.2;
        } else {
            arrival->error = ARRIVAL_ERROR_NULL;
        }
        arrival->quality = 0;

        // check if should use amplitude
        if (eth_use_mag != 0)
            arrival->amplitude = AMPLITUDE_NULL;

        return (istat);
    }//  ETH_LOC replaced by SED_LOC
        //else if (strcmp(ftype_obs, "ETH_LOC") == 0)
    else if (0) {

        /* example:
                        1 1   200.0   300.0    1.71 1                        INST
                        0.000     0.000  10.00 0000 00 00 00 00  0.00 0    TRIAL
                        maraini                                               AUTOR
                        2001/12/04 01:29                                      Local           20010748
                        FUORN   P       E   14.910 1 0.14      448 NWA       0
                        FUORN   S       E   16.140 1 0.14      448 NWA       0
                        OSS2    P       ID  16.428 1 0.09     2554 AHP       0
                        OSS     P       ED  16.524 1 0.11     1672 AHP       0
                        BERNI   P       E   17.930 1 0.18       46 NWA       0
                        BERNI   S       E   21.139 1 0.18       46 NWA       0
                        CHDAW   P       Q   20.016 1 0.56       10 NWA       0
                        VDL2    P       IU  23.704 1 0.20      246 AHP       0
                        VDL     P       E   23.863 1 0.12      188 AHP       0
                        VDL     P       E   23.822 1 0.28        4 NWA       1
                        VDL     S       E   31.177 1 0.28        4 NWA       1
                        SKIP
                        20011204012913146568N010226E00614Ml1239814161001001000SEDG024196006511         B
                        KP200112040128                                        ARCHIVE
         */

        // read line
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);

        // check for end of event (assumes empty lines between event)
        if (LineIsBlank(line)) {
            // end of event
            return (OBS_FILE_END_OF_EVENT);
        }


        // read event hypocenter line
        if (!in_hypocenter_event) {
            // find origin time line
            ifound = 0;
            while ((istat = ReadFortranString(line, 55, 4, eth_line_key)) > 0) {
                if (strcmp(eth_line_key, "Loca") == 0
                        || strcmp(eth_line_key, "Regi") == 0
                        || strcmp(eth_line_key, "Tele") == 0) {
                    ifound = 1;
                    break;
                }
                // read next line
                cstat = fgets(line, MAXLINE_LONG, fp_obs);
                if (cstat == NULL)
                    return (OBS_FILE_END_OF_INPUT);
            }
            if (!ifound)
                return (OBS_FILE_END_OF_EVENT);

            // read hypocenter time
            //2001/12/04 01:29                                      Local           20010748
            istat = ReadFortranInt(line, 1, 4, &EventTime.year);
            istat += ReadFortranInt(line, 6, 2, &EventTime.month);
            istat += ReadFortranInt(line, 9, 2, &EventTime.day);
            istat += ReadFortranInt(line, 12, 2, &EventTime.hour);
            istat += ReadFortranInt(line, 15, 2, &EventTime.min);
            //printf("EventTime: %d/%d/%d %d:%d\n",
            //EventTime.year, EventTime.month, EventTime.day,
            //EventTime.hour, EventTime.min);
            if (istat == 5) {
                in_hypocenter_event = 1;
                // read until phase line reached
                while ((istat = ReadFortranString(line, 55, 4, eth_line_key)) > 0
                        && strcmp(line, "    ") != 0) {
                    cstat = fgets(line, MAXLINE_LONG, fp_obs);
                    if (cstat == NULL)
                        return (OBS_FILE_END_OF_INPUT);
                }
            } else {
                return (OBS_FILE_END_OF_EVENT);
            }
        }

        // check if valid phase line: key is blank or strlen < 55
        if ((istat = ReadFortranString(line, 55, 4, eth_line_key)) > 0
                && strcmp(line, "    ") != 0) {
            return (OBS_FILE_SKIP_INPUT_LINE);
        }

        // read phase arrival input
        //CHDAW   P       Q   20.016 1 0.56       10 NWA       0
        //OSS     P       ED  16.524 1 0.11     1672 AHP       0
        //BERNI   P       E   17.930 1 0.18       46 NWA       0

        istat = ReadFortranString(line, 1, 5, arrival->label);
        istat += ReadFortranString(line, 9, 6, arrival->phase);
        istat += ReadFortranString(line, 17, 1, arrival->onset);
        istat += ReadFortranString(line, 18, 1, arrival->first_mot);
        istat += ReadFortranReal(line, 20, 7, &arrival->sec);
        //              istat += ReadFortranReal(line, 21, 6, &arrival->sec);
        istat += ReadFortranInt(line, 28, 1, &eth_use_loc);
        istat += ReadFortranReal(line, 29, 5, &arrival->period);
        istat += ReadFortranReal(line, 34, 9, &arrival->amplitude);
        //20040226 S Husen bug fix  istat += ReadFortranString(line, 44,9, arrival->inst);
        istat += ReadFortranString(line, 44, 4, arrival->inst);
        istat += ReadFortranInt(line, 54, 1, &eth_use_mag);


        if (istat != 10) {
            return (OBS_FILE_END_OF_EVENT);
        }

        /* remove blanks/whitespace */
        removeSpace(arrival->label);
        removeSpace(arrival->phase);
        removeSpace(arrival->inst);

        arrival->year = EventTime.year;
        arrival->month = EventTime.month;
        arrival->day = EventTime.day;
        arrival->hour = EventTime.hour;
        arrival->min = EventTime.min;

        /* convert onset to error */
        //Uncertainties of phase readings:
        //                              I             E            Q
        //for Key 'Loca':            < +/- 0.05   < +/- 0.2     > +/- 0.2
        //for Key 'Regi' or 'Tele'   < +/- 0.20   < +/- 1.0     > +/- 1.0
        snprintf(arrival->error_type, sizeof(arrival->error_type), "%s", "GAU");
        if (eth_use_loc == 0) {
            arrival->error = ARRIVAL_ERROR_NULL;
        } else if (strcmp(arrival->onset, "I") == 0) {
            arrival->error = 0.05;
        } else if (strcmp(arrival->onset, "E") == 0) {
            arrival->error = 0.2;
        } else if (strcmp(arrival->onset, "Q") == 0) {
            arrival->error = 0.2;
        } else {
            arrival->error = ARRIVAL_ERROR_NULL;
        }
        arrival->quality = 0;

        // check if should use amplitude
        if (eth_use_mag != 0)
            arrival->amplitude = AMPLITUDE_NULL;

        return (istat);
    } else if (strcmp(ftype_obs, "NEIC") == 0) {

        /* example:
        20 JUN 2003  (171)

        ot  = 06:19:38.51   +/-   0.17              AMAZONAS, BRAZIL
        lat =      -7.513   +/-    3.5
        lon =     -71.642   +/-    4.6              MAGNITUDE 7.1 (HRV)
        dep =       555.2  (depth phases)

        110 km (70 miles) E of Cruzeiro do Sul, Brazil (pop N/A)
        335 km (205 miles) ENE of Pucallpa, Peru
        410 km (255 miles) SSW of Leticia, Colombia
        2730 km (1700 miles) WNW of BRASILIA, Brazil

        nph =  194 of 401    se = 0.83        FE=113                      A

        error ellipse = ( 65.1,155.1,  0.0;  0.0,  0.0,  0.0;  6.3,  4.0,  0.0)

        mb = 6.4 ( 93)  ML = 0.0 (  0)  mblg = 0.0 (  0)  md = 0.0 (  0)  MS = 0.0 (  0)

        sta  phase     arrival     res   dist azm    amp  per mag     amp  per mag  sta
        SAML iPc     06:21:44.12   2.0    8.5 100                                   SAML
        eS      06:23:25.14  103.X
        LPAZ eP      06:21:51.01  -0.7    9.4 159 L:3.7+0 .93 6.2X g:3.0+0 .93 5.4X LPAZ
        OTAV ePc     06:22:02.21   1.7   10.3 318 g:9.9+0 1.0 6.0X                  OTAV
        eS      06:23:47.88  107.X
         */
        /* read line */
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        /*DEBUG*///printf("%s", line);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);
        /* check for end of event (assumes empty lines between event) */
        if (LineIsBlank(line)) {
            /* end of event */
            return (OBS_FILE_END_OF_EVENT);
        }


        /* read event hypocenter line */
        if (!in_hypocenter_event) {
            //(">>event hypocenter line --- %s\n", line);
            /* read hypocenter time */
            //     20 JUN 2003  (171)
            istat = ReadFortranInt(line, 13, 4, &EventTime.year);
            istat += ReadFortranString(line, 9, 3, cmonth);
            istat += ReadFortranInt(line, 6, 2, &EventTime.day);
            if (istat == 3 && EventTime.year > 0 && EventTime.year < 2199 && EventTime.day > 0 && EventTime.day < 32 && strlen(cmonth) == 3) {
                in_hypocenter_event = 1;
                EventTime.month = Month2Int(cmonth);
                origin_hour = 0;
                // find phs line
                while (strncmp(line, " sta", 4) != 0) {
                    // check for and get origin hour
                    if (strstr(line, "ot  =") != NULL) {
                        ReadFortranInt(line, 12, 2, &origin_hour);
                    }
                    // read next line
                    cstat = fgets(line, MAXLINE_LONG, fp_obs);
                    /*DEBUG*///printf("%s", line);
                    if (cstat == NULL)
                        return (OBS_FILE_END_OF_INPUT);
                }
                /* read next line */
                cstat = fgets(line, MAXLINE_LONG, fp_obs);
                /*DEBUG*///printf("%s", line);
                if (cstat == NULL)
                    return (OBS_FILE_END_OF_INPUT);
            } else {
                return (OBS_FILE_END_OF_EVENT);
            }
        }


        /* read phase arrival input */
        //printf(">>phase arrival line --- %s\n", line);
        /*
         SAML iPc     06:21:44.12   2.0    8.5 100                                   SAML
              eS      06:23:25.14  103.X
         LPAZ eP      06:21:51.01  -0.7    9.4 159 L:3.7+0 .93 6.2X g:3.0+0 .93 5.4X LPAZ
         KDAK ePn     09:34:19.19  -1.8   16.7  57 b:1.2+3 .50 6.3                   KDAK
         SVW2 eP      09:34:20.52   2.8   16.3  44                                   SVW2
         PMR  ePn     09:34:54.42   1.4   19.3  46 b:1.3+3 1.1 6.1                   PMR
         */
        istat = ReadFortranString(line, 2, 4, arrival->label);
        if (strncmp(arrival->label, "    ", 4) == 0)
            snprintf(arrival->label, sizeof(arrival->label), "%s", last_label);
        else
            snprintf(last_label, sizeof(last_label), "%s", arrival->label);
        TrimString(arrival->label);
        TrimString(last_label);
        istat += ReadFortranString(line, 7, 1, arrival->onset);
        istat += ReadFortranString(line, 8, 6, arrival->phase);
        TrimString(arrival->phase);
        // check for valid onset
        pchr = arrival->onset;
        // onset 'x' added by AJL to allow zero weighting of phases
        if (!(*pchr == ' ' || *pchr == 'i' || *pchr == 'e' || *pchr == 'q' || *pchr == 'x')) {
            snprintf(chrtmp, sizeof(chrtmp), "%c%s", *pchr, arrival->phase);
            snprintf(arrival->phase, sizeof(arrival->phase), "%s", chrtmp);
            *pchr = ARRIVAL_NULL_CHR;
        }
        // check for first motion
        pchr = &(arrival->phase[strlen(arrival->phase) - 1]);
        if (*pchr == 'c' || *pchr == 'd') {
            arrival->first_mot[0] = *pchr;
            *pchr = '\0';
        }
        // set weight, check for P - reset weight
        if (arrival->onset[0] == 'i')
            arrival->quality = 0;
        else if (arrival->onset[0] == 'e')
            arrival->quality = 1;
        else if (arrival->onset[0] == 'x')
            arrival->quality = 4;
        else
            arrival->quality = 2;
        if (strncmp(arrival->phase, "P", 1) != 0) {
            arrival->quality += 1;
        }
        istat += ReadFortranInt(line, 15, 2, &arrival->hour);
        // check for day jump
        //printf("   arrival->hour %d  origin_hour %d\n", arrival->hour, origin_hour);
        if (arrival->hour < origin_hour)
            arrival->hour += 24;
        //printf(">>>   arrival->hour %d  origin_hour %d\n", arrival->hour, origin_hour);
        istat += ReadFortranInt(line, 18, 2, &arrival->min);
        istat += ReadFortranReal(line, 21, 5, &arrival->sec);

        if (istat != 6) {
            return (OBS_FILE_END_OF_EVENT);
        }

        arrival->year = EventTime.year;
        arrival->month = EventTime.month;
        arrival->day = EventTime.day;


        /* convert quality to error */
        Qual2Err(arrival);

        return (istat);


    } else if (strcmp(ftype_obs, "NEIC_LATEST_EQS") == 0) {

        // http://earthquake.usgs.gov/earthquakes/map/

        /* example:

        Channel	Distance	Azimuth	Phase	Arrival Time	Status	Residual	Weight
        IM PD31 BHZ --	0.378028	236.911	Pn	2013-09-21T13:16:45.66Z	manual	-0.10	0.0160
        IM PD31 BHN --	0.378028	236.911	Sn	2013-09-21T13:16:55.11Z	manual	0.00	0.0340
        IM PDAR SHZ FB	0.378181	236.836	Pn	2013-09-21T13:16:45.75Z	manual	0.00	0.0160
        IM PDAR SHZ FB	0.378181	236.836	Sn	2013-09-21T13:16:55.35Z	manual	0.20	0.0340
        US BW06 BHZ 00	0.378425	236.86	Pn	2013-09-21T13:16:45.58Z	manual	-0.20	0.0160
        US BW06 BHN 00	0.378425	236.86	Pn	2013-09-21T13:16:45.67Z	manual	0.00	0.0000

         *        */
        /* read line */
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        //printf(">1 %s", line);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);
        /* check for end of event (assumes empty lines between event) */
        if (LineIsBlank(line)) {
            /* end of event */
            return (OBS_FILE_END_OF_EVENT);
        }


        if (0 && !in_hypocenter_event) {
            // find phase line
            while (strncmp(line, "Channel", 7) != 0) {
                // read next line
                cstat = fgets(line, MAXLINE_LONG, fp_obs);
                //printf(">2 %s", line);
                if (cstat == NULL)
                    return (OBS_FILE_END_OF_INPUT);
            }
            in_hypocenter_event = 1;
        }


        /* read phase arrival input */
        /*
        IM PD31 BHZ --	0.378028	236.911	Pn	2013-09-21T13:16:45.66Z	manual	-0.10	0.0160
        IM PD31 BHN --	0.378028	236.911	Sn	2013-09-21T13:16:55.11Z	manual	0.00	0.0340
         */
        istat = sscanf(line, "%s %s %s %*s %*f %*f %s %4d-%2d-%2dT%2d:%2d:%lfZ",
                arrival->network, arrival->label, arrival->inst, arrival->phase, &arrival->year, &arrival->month, &arrival->day, &arrival->hour, &arrival->min, &arrival->sec);
        //printf("%s %s %s - - - %s %4d-%2d-%2dT%2d:%2d:%fZ\n", arrival->network, arrival->label, arrival->inst, arrival->phase, arrival->year, arrival->month, arrival->day, arrival->hour, arrival->min, arrival->sec);

        if (istat != 10) {
            return (OBS_FILE_END_OF_EVENT);
        }

        /* convert quality to error */
        //Qual2Err(arrival);
        snprintf(arrival->error_type, sizeof(arrival->error_type), "%s", "GAU");
        arrival->error = 2.0;
        // error if first arrival P
        if (strstr("P$Pg$Pn$Pb$P0$P1$PKP$PKPdf", arrival->phase) != NULL) {
            arrival->error = 1.0;
        }


        return (istat);

    } else if (strcmp(ftype_obs, "GEM") == 0) {

        /* example:

        132750016  261603062   KSA       P       P  84.13 355.75 1929-07-03 01:05:49 ?
        132750016  261603063   KSA       S       S  84.13 355.75 1929-07-03 01:15:55 ?
        132750017  261603064   LPZ       P         100.36 105.55 1929-07-03 01:45:50 ?

        #     evid      hypid author origin time               lat      lon    depth F  nsta  nass  ndef
    905632  607283695    GEM 1933/06/12 15:23:41.916  61.231 -151.354  15.00 H    21    68    26

        #     rdid       phid   sta  Pphase  Aphase  delta   esaz    arrival time    ch
        132814560  261709084   VIC       P      Pn  20.30 116.11 1933-06-12 15:28:16 ?
        132814560  261709085   VIC      SS      Sn  20.30 116.11 1933-06-12 15:32:19 ?
        132814561  261709086   BZM      SS      SS  28.04 105.18 1933-06-12 15:36:24 ?
        132814562  261709087   UKI       S       S  28.12 128.91 1933-06-12 15:34:30 ?

         *        */
        /* read line */
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        //printf(">1 %s", line);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);
        /* check for end of event (assumes empty lines between event) */
        if (LineIsBlank(line)) {
            /* end of event */
            return (OBS_FILE_END_OF_EVENT);
        }

        if (!in_hypocenter_event) {
            // find phase line
            while (strstr(line, "Pphase") == NULL) {
                // read next line
                cstat = fgets(line, MAXLINE_LONG, fp_obs);
                //printf(">2 %s", line);
                if (cstat == NULL)
                    return (OBS_FILE_END_OF_INPUT);
            }
            in_hypocenter_event = 1;
            cstat = fgets(line, MAXLINE_LONG, fp_obs);
        }


        /* read phase arrival input */
        /*
        132814560  261709084   VIC       P      Pn  20.30 116.11 1933-06-12 15:28:16 ?
         */
        istat = ReadFortranString(line, 25, 9, chrtmp);
        istat = sscanf(chrtmp, "%s", arrival->label);
        istat = ReadFortranString(line, 42, 3, chrtmp);
        istat = sscanf(chrtmp, "%s", arrival->phase);
        istat = ReadFortranString(line, 59, 22, chrtmp);
        istat = sscanf(chrtmp, "%4d-%2d-%2d %2d:%2d:%lf %s",
                &arrival->year, &arrival->month, &arrival->day, &arrival->hour, &arrival->min, &arrival->sec, arrival->comp);
        //printf("%s %s %4.4d-%2.2d-%2.2d %2.2d:%2.2d:%f %s\n", arrival->label, arrival->phase, arrival->year, arrival->month, arrival->day, arrival->hour, arrival->min, arrival->sec, arrival->comp);

        if (istat != 7) {
            return (OBS_FILE_END_OF_EVENT);
        }

        /* convert quality to error */
        arrival->error = Quality2Error[0];

        return (istat);

    } else

        if (strncmp(ftype_obs, "ISC_IMS1.0", 10) == 0) {

        // ISC_IMS1.0 - standard
        // extensions (e.g. ISC_IMS1.0_LOCPHASEID,  ISC_IMS1.0_LOCPHASEID_RESID_WT, etc.)
        // *_LOCPHASEID - convert phase name using LOCPHASEID (homogenizes names for LOCDELAY accumulation)
        // *_RESID_WT - set error to 9999.0 if ISC Tres > 4*nominal_error from LOCQUAL2ERR


        /* example:
                        ...
                        DATA_TYPE BULLETIN IMS1.0:short
                        ISC Comprehensive Bulletin
                        Event  2957613 Sicily

                        Date       Time        Err   RMS Latitude Longitude  Smaj  Smin  Az Depth   Err Ndef Nsta Gap  mdist  Mdist Qual   Author      OrigID
                        2002/04/05 04:52:15.70         0.60  38.7500   14.2800                  14.0               6                          BJI        4180045
                        2002/04/05 04:52:20.23   0.05        37.9440   16.4600  38.9  17.4  -1  10.0   1.1         9 348                      PDG        4681048
                        ...
                        Sta     Dist  EvAz Phase        Time      TRes  Azim AzRes   Slow   SRes Def   SNR       Amp   Per Qual Magnitude    ArrID
                        MSI     0.45 121.8 Pg       04:52:30.8     1.0                           T__                        _e            71903810
                        SCLL    0.55 110.1 Pg       04:52:32.5     1.0                           T__                        _e            71903811
                        MTTG    0.67 131.6 Pg       04:52:34.5     0.6                           T__                        _e            71903812
                        ...
                        SGO     2.12   5.0 PN       04:52:59.0     2.6                           T__                        _e            71903828
                        LVI     2.19 258.8 PN       04:52:58.3     0.8                           T__                        _e            71903826
                        LVI     2.19 258.8                                                       ___           736.0  0.60  __            71903827
                        MRLC    2.33   8.0 PN       04:53:01.3     1.9                           T__                        _e            71903830
                        CSSN    2.41 359.2 PN       04:53:01.0     0.4                           T__                        _e            71903829
                        ...
                        AVF    11.99 317.8 P        04:55:12.9   -1.37                           T__            20.8  1.26  _e            71869701
                        LOR    12.04 320.6 P        04:55:13.3   -1.73                           T__            22.9  1.39  _e            71869702
                        LOR    12.04 320.6 R                                                     ___           126.3 18.75  _e            71869703
                        SSF    12.08 319.1 P        04:55:15.3    -0.3                           T__            12.1  1.04  _e            71869705
                        BGF    12.11 315.9 P        04:55:14.6   -1.37                           T__            34.9  1.38  _e            71869704
                        MOX    12.44 349.8 P        04:55:31.3   10.96                           T__                        _e            75622714
                        MOX    12.44 349.8 LR                                                    ___           700.0 16.00  __        3.8 75622715
                        MOX    12.44 349.8 LR                                                    ___           300.0 16.00  __        3.8 75622717
                        MOX    12.44 349.8 LR                                                    ___           600.0 15.00  __            75622716
                        BRG    12.45 356.7 P        04:55:21.5     1.0                           T__             9.7  1.20  __            75622621
                        BRG    12.45 356.7          04:55:40.5                                   ___                        _i            75622622
                        BRG    12.45 356.7 LR                                                    ___           870.0 13.80  __        3.9 75622623
                        BRG    12.45 356.7 LR                                                    ___            60.0 13.80  __        3.8 75622625
                        BRG    12.45 356.7 LR                                                    ___           690.0 13.80  __            75622624
                        ETSF   12.66 295.5 P        04:55:22.4    -1.0                           T__            31.9  1.41  _e            71869706
                        CLL    12.94 354.2 P        04:55:26     -1.05                           T__                        __            79343809
                        ...


                        STOP

         */

        long int file_pos = 0;

        /* read line */
        do {
            file_pos = ftell(fp_obs);
            cstat = fgets(line, MAXLINE_LONG, fp_obs);
            //printf("IMS-1 %s", line);
            if (cstat == NULL)
                return (OBS_FILE_END_OF_INPUT);
        } while (LineIsBlank(line));
        // check for end of event (assumes Event or STOP at end of event or following DATA_TYPE)
        if (in_hypocenter_event && (strncmp(line, "Event", 5) == 0 || strncmp(line, "EVENT", 5) == 0)) {
            /* end of event */
            fseek(fp_obs, file_pos, SEEK_SET);
            return (OBS_FILE_END_OF_EVENT);
        }
        if (strncmp(line, "STOP", 4) == 0) {
            /* end of event */
            return (OBS_FILE_END_OF_EVENT);
        }
        // 20210917 AJL - added to support DATA_TYPE BULLETIN IMS1.0:short with ISF2.0 extensions Bulletin from IGN
        if (strncmp(line, "DATA_TYPE", 4) == 0) {
            /* end of event */
            return (OBS_FILE_END_OF_EVENT);
        }


        /* not yet in event, find and read event hypocenter line */
        if (!in_hypocenter_event) {
            // assume may already be in "Event" block
            snprintf(phypo->public_id, sizeof(phypo->public_id), "%s", "-1"); // 20191209 AJL - added
            // find date
            while (strncmp(line, "   Date", 7) != 0) {
                // check for Event line to get event id  // 20191209 AJL - added
                if (strncmp(line, "Event", 5) == 0 || strncmp(line, "EVENT", 5) == 0) {
                    sscanf(line, "%*s %s ", phypo->public_id);
                }
                // read next line
                do {
                    file_pos = ftell(fp_obs);
                    cstat = fgets(line, MAXLINE_LONG, fp_obs);
                    //printf("IMS-2 %s", line);
                    if (cstat == NULL)
                        return (OBS_FILE_END_OF_INPUT);
                } while (LineIsBlank(line));
            }
            // read next line
            do {
                file_pos = ftell(fp_obs);
                cstat = fgets(line, MAXLINE_LONG, fp_obs);
                //printf("IMS-3 %s", line);
                if (cstat == NULL)
                    return (OBS_FILE_END_OF_INPUT);
            } while (LineIsBlank(line));
            /* read hypocenter date */
            //2002/04/05 04:52:15.70         0.60  38.7500   14.2800                  14.0               6                          BJI        4180045
            istat = ReadFortranInt(line, 1, 4, &EventTime.year);
            istat += ReadFortranInt(line, 6, 2, &EventTime.month);
            istat += ReadFortranInt(line, 9, 2, &EventTime.day);
            istat += ReadFortranInt(line, 12, 2, &EventTime.hour);
            if (istat == 4) {
                in_hypocenter_event = 1;
                int in_magnitude = 0;
                double amp_mag = -9.9;
                int num_amp_mag = -1;
                // find phases lines
                while (strncmp(line, "Sta     Dist", 12) != 0) {
                    // read next line
                    do {
                        file_pos = ftell(fp_obs);
                        cstat = fgets(line, MAXLINE_LONG, fp_obs);
                        //printf("IMS-4 %s", line);
                        if (cstat == NULL)
                            return (OBS_FILE_END_OF_INPUT);
                    } while (LineIsBlank(line));
                    // check for end of event (assumes Event or STOP at end of event
                    if (strncmp(line, "Event", 5) == 0 || strncmp(line, "EVENT", 5) == 0) {
                        /* end of event */
                        fseek(fp_obs, file_pos, SEEK_SET);
                        return (OBS_FILE_END_OF_EVENT);
                    }
                    if (strncmp(line, "STOP", 4) == 0) {
                        /* end of event */
                        return (OBS_FILE_END_OF_EVENT);
                    }
                    if (in_magnitude) {
                        // read magnitude
                        //Magnitude  Err Nsta Author      OrigID
                        //mb     6.5 0.2  156 ISC        1090877
                        //MS     7.9 0.2   37 ISC        1090877
                        istat = sscanf(line, "%*s %lf %*s %d %*s %*s", &amp_mag, &num_amp_mag);
                        if (istat == 2) {
                            // use largest magnitude read (!!!)
                            if (amp_mag > phypo->amp_mag) {
                                phypo->amp_mag = amp_mag;
                                phypo->num_amp_mag = num_amp_mag;
                            }
                            //printf("DEBUG: ISC_IMS1.0 magnitude: %lf %d\n", phypo->amp_mag, phypo->num_amp_mag);
                        } else {
                            in_magnitude = 0;
                        }
                    } else if (strncmp(line, "Magnitude", 9) == 0) {
                        in_magnitude = 1;
                    }


                }
                /* read next line */
                do {
                    file_pos = ftell(fp_obs);
                    cstat = fgets(line, MAXLINE_LONG, fp_obs);
                    //printf("IMS-5 %s", line);
                    if (cstat == NULL)
                        return (OBS_FILE_END_OF_INPUT);
                } while (LineIsBlank(line));
                // check for end of event (assumes Event or STOP at end of event
                if (strncmp(line, "Event", 5) == 0 || strncmp(line, "EVENT", 5) == 0) {
                    /* end of event */
                    fseek(fp_obs, file_pos, SEEK_SET);
                    return (OBS_FILE_END_OF_EVENT);
                }
                if (strncmp(line, "STOP", 4) == 0) {
                    /* end of event */
                    return (OBS_FILE_END_OF_EVENT);
                }
            } else {
                return (OBS_FILE_END_OF_EVENT); // blank line or ignored line
            }
        }


        /* read phase arrival input */
        /*
                                           LVI     2.19 258.8 PN       04:52:58.3     0.8                           T__                        _e            71903826
                                           LVI     2.19 258.8                                                       ___           736.0  0.60  __            71903827
         */

        // check if is a reading with a seconds time
        //while (ReadFortranString(line, 44, 3, isc_time_str) != 1
        while (ReadFortranString(line, 35, 3, isc_time_str) != 1
                || strncmp(isc_time_str, "   ", 3) == 0) {
            // not a time, check if end
            if (cstat == NULL)
                return (OBS_FILE_END_OF_INPUT);
            // check for end of event (assumes Event or STOP at end of event
            if (strncmp(line, "Event", 5) == 0 || strncmp(line, "EVENT", 5) == 0) {
                /* end of event */
                fseek(fp_obs, file_pos, SEEK_SET);
                return (OBS_FILE_END_OF_EVENT);
            }
            if (strncmp(line, "STOP", 4) == 0) {
                /* end of event */
                return (OBS_FILE_END_OF_EVENT);
            }
            // not a time, not end, read next line
            do {
                file_pos = ftell(fp_obs);
                cstat = fgets(line, MAXLINE_LONG, fp_obs);
                //printf("IMS-6 %s", line);
                if (cstat == NULL)
                    return (OBS_FILE_END_OF_INPUT);
            } while (LineIsBlank(line));
            // check for end of event (assumes Event or STOP at end of event
            if (strncmp(line, "Event", 5) == 0 || strncmp(line, "EVENT", 5) == 0) {
                /* end of event */
                fseek(fp_obs, file_pos, SEEK_SET);
                return (OBS_FILE_END_OF_EVENT);
            }
            if (strncmp(line, "STOP", 4) == 0) {
                /* end of event */
                return (OBS_FILE_END_OF_EVENT);
            }
        }

        // read phase time reading
        istat = ReadFortranString(line, 1, 5, arrival->label);
        if (strncmp(arrival->label, "     ", 5) == 0)
            snprintf(arrival->label, sizeof(arrival->label), "%s", last_label);
        else
            snprintf(last_label, sizeof(last_label), "%s", arrival->label);
        TrimString(arrival->label);
        TrimString(last_label);
        istat += ReadFortranString(line, 20, 8, arrival->phase);
        TrimString(arrival->phase);

        // convert phase name using LOCPHASEID if requested (homogenizes names for LOCDELAY accumulation)
        // 20190822 AJL - added
        if (strstr(ftype_obs, "_LOCPHASEID") != NULL) {
            EvalPhaseID(eval_phase_tmp, sizeof(eval_phase_tmp), arrival->phase);
            snprintf(arrival->phase, sizeof(arrival->phase), "%s", eval_phase_tmp);
        }

        int ioffset = 0;
        if (strstr(ftype_obs, "_GEIN_NOA") != NULL) { // DATA_TYPE BULLETIN IMS1.0:LONG ???
            ioffset = 1;
        }

        istat += ReadFortranInt(line, 29 + ioffset, 2, &arrival->hour);
        // BUG FIX - 20091019 signaled by Edith Korger (Edith.Korger@awi.de)
        // ISC_IMS1.0 uses 00h for phases arriving after 24h of day of event !!!
        if (arrival->hour < EventTime.hour)
            arrival->hour += 24;
        //
        istat += ReadFortranInt(line, 32 + ioffset, 2, &arrival->min);
        istat += ReadFortranReal(line, 35 + ioffset, 6, &arrival->sec);
        //printf("IMS-date %s\n", line);
        //printf("IMS-date %s\n", line + 29 + ioffset);
        //printf("IMS-date %d %d %f\n", arrival->hour, arrival->min, arrival->sec);
        double tres = 0.0;
        istat += ReadFortranReal(line, 40 + ioffset, 7, &tres); // TRes
        istat += ReadFortranString(line, 101, 1, arrival->first_mot);
        istat += ReadFortranString(line, 102, 1, arrival->onset);
        // set weight, check for P - reset weight
        if (arrival->onset[0] == 'i')
            arrival->quality = 0;
        else if (arrival->onset[0] == 'e')
            arrival->quality = 1;
        else if (arrival->onset[0] == 'q')
            arrival->quality = 2;
        else if (arrival->onset[0] == '_')
            arrival->quality = 2;
        else
            arrival->quality = 3;
        // check if "P" arrival
        /* 20160929 AJL - replace with multiplication of error below
        if (if arrival->quality < 3 && strncmp(arrival->phase, "P", 1) != 0) {
            arrival->quality += 1;
        }*/

        // special processing for quality in phase code (e.g. "P_3")
        if (strstr(ftype_obs, "_GEIN_NOA") != NULL) {
            char *csep = 0;
            if ((csep = strstr(arrival->phase, "_")) != NULL) {
                int iqual;
                int istat2 = sscanf(csep + 1, "%d", &iqual);
                if (istat2 > 0) {
                    arrival->quality = iqual;
                    arrival->phase[csep - arrival->phase] = '\0';
                }
            }

        }



        if (istat != 8) {
            printf("DEBUG_IMS OBS_FILE_INVALID_PHASE %d\n", istat);
            return (OBS_FILE_INVALID_PHASE);
        }

        arrival->year = EventTime.year;
        arrival->month = EventTime.month;
        arrival->day = EventTime.day;


        // convert quality to error
        Qual2Err(arrival);
        if (strstr(ftype_obs, "_GEIN_NOA") == NULL) {
            // 20160929 AJL - multiplication of error if not first arrival P
            // 20190815  if (0 && strstr("P$Pg$Pn$Pb$P0$P1$PKP$PKPdf", arrival->phase) == NULL) {
            if (strstr("P$Pg$Pn$Pb$P0$P1$PKP$PKPdf", arrival->phase) == NULL) {
                arrival->error *= 4.0;
            }
        }

        // set error to 9999.0 if ISC Tres > 4*nominal_error from LOCQUAL2ERR
        // 20190913 AJL - added
        if (strstr(ftype_obs, "_RESID_WT") != NULL) {
            if (fabs(tres) > arrival->error * 2.0) {
                arrival->error = 9999.0;
            }
        }


        return (istat);

    } else if (strcmp(ftype_obs, "CSEM_GSE2.0") == 0) {

        /* example:
                        ...
                        BEGIN GSE2.0
                        MSG_TYPE DATA
                        MSG_ID 041203222858 EMSC
                        DATA_TYPE BULLETIN GSE2.0

                        EVENT 20041203222858
                        Date       Time       Latitude Longitude    Depth    Ndef Nsta Gap    Mag1  N    Mag2  N    Mag3  N  Author          ID
                        rms   OT_Error      Smajor Sminor Az        Err   mdist  Mdist     Err        Err        Err     Quality

                        2004/12/03 22:28:58.5     44.2742    7.5419      8.0      66   53 109  ML 3.3 11                        MIX       00000001
                        1.22   +-  0.10       2.6    1.6  123               0.42   7.13   +-0.6                           a i ke

                        ( Data from stations operated by :     LDG,   MDD,   ZUR )

                        NORTHERN ITALY
                        Sta     Dist  EvAz     Phase      Date       Time     TRes  Azim  AzRes  Slow  SRes Def   SNR       Amp   Per   Mag1   Mag2       ID
                        SBF     0.42 190.6 a   Pn      2004/12/03 22:29:09.3   1.3                          T     1.0                                    LDG
                        MBDF    0.71 309.7 a   Pg      2004/12/03 22:29:09.8  -2.4                          T     1.0                                    LDG
                        FRF     0.96 222.5 a   Pg      2004/12/03 22:29:16.4  -0.2                          T     1.0                                    LDG
                        FRF     0.96 222.5 a   Sg      2004/12/03 22:29:27.9  -1.3                          T     1.0     134.3  0.16 ML 4.0             LDG
                        ...
                        BAIF    6.22 339.7 a   Pn      2004/12/03 22:30:30.2   1.9                          T     1.0                                    LDG
                        ERTA    6.27 240.5 m   P       2004/12/03 22:30:31.9  -1.8                          T               0.8  0.16                    MAD
                        LDF     6.83 311.9 a   Pn      2004/12/03 22:30:38.8   2.1                          T     1.0                                    LDG
                        GRR     7.12 308.2 a   Pn      2004/12/03 22:30:41.8   1.1                          T     1.0                                    LDG
                        FLN     7.13 311.9 a   Pn      2004/12/03 22:30:41.8   1.1                          T     1.0                                    LDG

                        STOP

                        ...

         */
        /* read line */
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        printf("IMS %s", line);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);
        /* check for end of event (assumes Event or STOP at end of event) */
        if ((in_hypocenter_event && strncmp(line, "EVENT", 5) == 0)
                || strncmp(line, "STOP", 4) == 0) {
            /* end of event */
            return (OBS_FILE_END_OF_EVENT);
        }


        /* not yet in event, find and read event hypocenter line */
        if (!in_hypocenter_event) {
            // assume may alredy be in "Event" block
            // find phases lines
            // following works for CSEM:
            //while (strncmp(line, "Sta     Dist", 12) != 0) {
            // does not work for INGV which has: "Sta    Dist"
            while (strncmp(line, "Sta    ", 7) != 0) {
                // read next line
                cstat = fgets(line, MAXLINE_LONG, fp_obs);
                //printf("4 %s", line);
                if (cstat == NULL)
                    return (OBS_FILE_END_OF_INPUT);
            }
            // read next line
            cstat = fgets(line, MAXLINE_LONG, fp_obs);
            //printf("5 %s", line);
            if (cstat == NULL)
                return (OBS_FILE_END_OF_INPUT);
            in_hypocenter_event = 1;
        }


        // read phase arrival input



        /* this block needed for ISC_IMSL because of amplitude readings, etc.
        // check if is a reading with a seconds time
        //while (ReadFortranString(line, 44, 3, isc_time_str) != 1
                        while (ReadFortranString(line, 35, 3, isc_time_str) != 1
                        || strncmp(isc_time_str, "   ", 3) == 0) {
        // not a time, check if end
                        if (cstat == NULL || strncmp(line, "Event", 5 || strncmp(line, "EVENT", 5) == 0) == 0
                        || strncmp(line, "STOP", 4) == 0)
                        return(OBS_FILE_END_OF_EVENT);
        // not a time, not end, read next line
                        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        //printf("6 %s", line);
                        if (cstat == NULL)
                        return(OBS_FILE_END_OF_INPUT);
        }
         */

        // read phase time reading
        istat = ReadFortranString(line, 1, 5, arrival->label);
        if (strncmp(arrival->label, "     ", 5) == 0)
            snprintf(arrival->label, sizeof(arrival->label), "%s", last_label);
        else
            snprintf(last_label, sizeof(last_label), "%s", arrival->label);
        TrimString(arrival->label);
        TrimString(last_label);
        istat += ReadFortranString(line, 24, 8, arrival->phase);
        TrimString(arrival->phase);
        istat += ReadFortranInt(line, 32, 4, &arrival->year);
        istat += ReadFortranInt(line, 37, 2, &arrival->month);
        istat += ReadFortranInt(line, 40, 2, &arrival->day);
        istat += ReadFortranInt(line, 43, 2, &arrival->hour);
        istat += ReadFortranInt(line, 46, 2, &arrival->min);
        istat += ReadFortranReal(line, 49, 5, &arrival->sec);
        istat += ReadFortranString(line, 21, 1, arrival->first_mot);
        istat += ReadFortranString(line, 22, 1, arrival->onset);
        // set weight, check for P - reset weight
        if (arrival->onset[0] == 'i')
            arrival->quality = 0;
        else if (arrival->onset[0] == ' ')
            arrival->quality = 0;
        else if (arrival->onset[0] == 'e')
            arrival->quality = 1;
        else if (arrival->onset[0] == 'q')
            arrival->quality = 2;
        else if (arrival->onset[0] == '_')
            arrival->quality = 2;
        else
            arrival->quality = 3;
        // check if "P" arrival
        if (strncmp(arrival->phase, "P", 1) != 0) {
            arrival->quality += 1;
        }

        if (istat != 10) {
            return (OBS_FILE_INVALID_PHASE);
        }

        // convert quality to error
        Qual2Err(arrival);

        return (istat);
    } else if (strcmp(ftype_obs, "CSEM_ALERT") == 0) {

        /* example:

                        ...
                        Sta       Phase     Date     Time       Res     Dist Azm    Net
                        ---------------------------------------------------------------
                        PRNI  m   Pn    2004/02/11 08:15:27.7   0.4     1.41 197    GII
                        PRNI  m   S     2004/02/11 08:15:46.6   1.5     1.41 197    GII
                        KSDI  m   Pn    2004/02/11 08:15:29.9   1.6     1.49   5    GII
                        KSDI  a   P     2004/02/11 08:15:28.8   0.3     1.49   5    ODC
                        MAMC  ad  P     2004/02/11 08:16:03.9   0.4     3.95 331    CYP
                        UZH   m i P     2004/02/11 08:19:33.0   1.1    19.64 333    LVV
         */
        /* read line */
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);
        /* check for end of event (assumes empty lines between event) */
        if (LineIsBlank(line)) {
            /* end of event */
            return (OBS_FILE_END_OF_EVENT);
        }


        if (!in_hypocenter_event) {
            // find phase line
            while (strncmp(line, "Sta       Phase", 15) != 0) {
                // read next line
                cstat = fgets(line, MAXLINE_LONG, fp_obs);
                if (cstat == NULL)
                    return (OBS_FILE_END_OF_INPUT);
            }
            // skip next line */
            cstat = fgets(line, MAXLINE_LONG, fp_obs);
            if (cstat == NULL)
                return (OBS_FILE_END_OF_INPUT);
            in_hypocenter_event = 1;
        }


        /* read phase arrival input */
        /*
                        PRNI  m   Pn    2004/02/11 08:15:27.7   0.4     1.41 197    GII
                        PRNI  m   S     2004/02/11 08:15:46.6   1.5     1.41 197    GII
                        KSDI  m   Pn    2004/02/11 08:15:29.9   1.6     1.49   5    GII
                        KSDI  a   P     2004/02/11 08:15:28.8   0.3     1.49   5    ODC
                        MAMC  ad  P     2004/02/11 08:16:03.9   0.4     3.95 331    CYP
                        UZH   m i P     2004/02/11 08:19:33.0   1.1    19.64 333    LVV
         */
        istat = ReadFortranString(line, 1, 4, arrival->label);
        TrimString(arrival->label);
        istat += ReadFortranString(line, 11, 5, arrival->phase);
        TrimString(arrival->phase);
        istat += ReadFortranString(line, 8, 1, arrival->first_mot);
        if (arrival->first_mot[0] == ' ')
            arrival->first_mot[0] = ARRIVAL_NULL_CHR;
        istat += ReadFortranString(line, 9, 1, arrival->onset);
        if (arrival->onset[0] == ' ')
            arrival->onset[0] = ARRIVAL_NULL_CHR;
        // set weight, check for P - reset weight
        arrival->quality = arrival->onset[0] == 'i' ? 0 : 1;
        if (strncmp(arrival->phase, "P", 1) != 0) {
            arrival->quality += 1;
        }
        istat += ReadFortranInt(line, 17, 4, &arrival->year);
        istat += ReadFortranInt(line, 22, 2, &arrival->month);
        istat += ReadFortranInt(line, 25, 2, &arrival->day);
        istat += ReadFortranInt(line, 28, 2, &arrival->hour);
        istat += ReadFortranInt(line, 31, 2, &arrival->min);
        istat += ReadFortranReal(line, 34, 5, &arrival->sec);

        if (istat != 10) {
            return (OBS_FILE_END_OF_EVENT);
        }


        /* convert quality to error */
        Qual2Err(arrival);

        return (istat);


    } else if (strcmp(ftype_obs, "TEXNET_BULLETIN") == 0 || strcmp(ftype_obs, "TEXNET_BULLETIN__ZERO_WT_X") == 0) {

        /* example:
        Event:
    Public ID              texnet2017qjsf
    Type                   earthquake
    Description
      region name: Western Texas
        Origin:
    Date                   2017-08-21
    Time                   23:59:32.064
    Latitude                31.11088 deg  +/-    0.943 km
    Longitude             -103.21398 deg  +/-    1.163 km
    Depth                     10.065 km   +/-    2.882 km
    Agency                 TXNet
    Mode                   manual
    Status                 final
    Residual RMS               0.306 s
    Azimuthal gap               70.2 deg

        3 Network magnitudes:
    ML        2.04 +/- 0.04   7 preferred
    MLv       2.38 +/- 0.14   9
    M         2.26            9

        20 Phase arrivals:
    sta   net      dist   azi  phase   time             res     wt  sta
    PB08  TX     38.669 129.0  P       23:59:39.396  -0.044 M  1.3  PB08
    PB08  TX     38.669 129.0  S       23:59:46.382   1.408 MX 0.0  PB08
    PB02  TX     44.113 318.6  P       23:59:42.210   1.903 MX 0.0  PB02
    PB02  TX     44.113 318.6  S       23:59:47.113   0.623 M  0.6  PB02
    MNHN  TX     49.214  58.8  P       23:59:41.457   0.318 M  1.1  MNHN
    MNHN  TX     49.214  58.8  S       23:59:47.759  -0.187 M  1.1  MNHN
    PECS  TX     69.203 294.6  P       23:59:44.397  -0.078 M  1.3  PECS
    PECS  TX     69.203 294.6  S       23:59:53.018  -0.767 M  0.4  PECS
    ALPN  TX     90.224 204.9  P       23:59:48.309   0.260 M  1.2  ALPN
    ALPN  TX     90.224 204.9  S       23:59:59.725  -0.315 M  1.0  ALPN
    ODSA  TX    126.621  27.7  P       23:59:54.041   0.088 M  1.2  ODSA
    ODSA  TX    126.621  27.7  S       00:00:09.858  -0.517 M  0.7  ODSA
    CL2B  SC    144.494 332.5  P       23:59:57.516   0.588 M  0.7  CL2B
    CL7   SC    156.814 337.3  P       23:59:58.850  -0.107 M  1.2  CL7
    HTMS  SC    157.709 343.6  P       23:59:59.155   0.013 M  1.2  HTMS
    GDL2  SC    165.061 317.2  P       00:00:00.116  -0.163 M  1.1  GDL2
    SAND  TX    168.584 139.6  P       00:00:01.024   0.439 M  0.9  SAND
    VHRN  TX    171.704 257.9  P       00:00:01.384   0.241 M  1.1  VHRN
    VHRN  TX    171.704 257.9  S       00:00:22.569  -0.389 M  0.8  VHRN
    TX31  IM    200.960 190.8  P       00:00:06.262   1.544 MX 0.0  TX31

        18 Station magnitudes:
    sta   net      dist   azi  type   value   res        amp  per
    PB08  TX     38.669 129.0  ML      2.06  0.02    1.12763
    PB08  TX     38.669 129.0  MLv     2.02 -0.36   0.587801
    MNHN  TX     49.214  58.8  ML      2.09  0.05   0.870497
    MNHN  TX     49.214  58.8  MLv     2.30 -0.08   0.559911
    PECS  TX     69.203 294.6  ML      1.99 -0.05   0.482742
    PECS  TX     69.203 294.6  MLv     2.25 -0.12   0.257526
    ALPN  TX     90.224 204.9  ML      1.94 -0.10   0.304913
    ALPN  TX     90.224 204.9  MLv     2.45  0.07   0.311978
    ODSA  TX    126.621  27.7  ML      2.08  0.04   0.284525
    ODSA  TX    126.621  27.7  MLv     2.15 -0.22   0.102733
    CL2B  SC    144.494 332.5  MLv     2.51  0.14   0.200041
    CL7   SC    156.814 337.3  MLv     2.64  0.26   0.230992
    HTMS  SC    157.709 343.6  MLv     2.66  0.28   0.238286
    GDL2  SC    165.061 317.2  MLv     2.36 -0.02   0.109993
    SAND  TX    168.584 139.6  ML      2.05  0.01   0.196648
    SAND  TX    168.584 139.6  MLv     2.44  0.07   0.130226
    VHRN  TX    171.704 257.9  ML      2.02 -0.02   0.172512
    VHRN  TX    171.704 257.9  MLv     2.29 -0.09   0.0844888


         */

        // read next line
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        //printf("1 %s", line);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);
        if (in_hypocenter_event) {
            if (strncmp(line, "Event:", 6) == 0) {
                // end of event
                return (OBS_FILE_END_OF_EVENT);
            }
        }


        // not yet in event, find and read event hypocenter line
        if (!in_hypocenter_event) {
            // assume may already be in "Event" block
            // find date
            while (strstr(line, "Public ID") == NULL) {
                // read next line
                cstat = fgets(line, MAXLINE_LONG, fp_obs);
                //printf("2 %s", line);
                if (cstat == NULL)
                    return (OBS_FILE_END_OF_INPUT);
                if (strncmp(line, "Event:", 6) == 0) {
                    // end of event
                    return (OBS_FILE_END_OF_EVENT);
                }
            }
            // read hypocenter public id
            //    Date                   2017-08-01
            istat = sscanf(line, " Public ID %s", phypo->public_id);
            // find date
            while (strstr(line, "Date") == NULL) {
                // read next line
                cstat = fgets(line, MAXLINE_LONG, fp_obs);
                //printf("2 %s", line);
                if (cstat == NULL)
                    return (OBS_FILE_END_OF_INPUT);
                if (strncmp(line, "Event:", 6) == 0) {
                    // end of event
                    return (OBS_FILE_END_OF_EVENT);
                }
            }
            // read hypocenter date
            //    Date                   2017-08-01
            istat = sscanf(line, " Date %d-%d-%d", &EventTime.year, &EventTime.month, &EventTime.day);

            if (istat == 3) {
                in_hypocenter_event = 1;

                // find time line
                //     Time                   23:59:32.064
                int ifound = 1;
                while (strstr(line, "Time") == NULL) {
                    // read next line
                    cstat = fgets(line, MAXLINE_LONG, fp_obs);
                    //printf("DEBUG: TEXNET_BULLETIN preferred magnitude: <%s>\n", line);
                    //printf("3 %s", line);
                    if (cstat == NULL)
                        return (OBS_FILE_END_OF_INPUT);
                    if (strncmp(line, "Event:", 6) == 0) {
                        // end of event
                        return (OBS_FILE_END_OF_EVENT);
                    }
                    if (strstr(line, "Phase") != NULL) {
                        ifound = 0;
                        break;
                    }
                }
                if (ifound) {
                    // read time
                    //     Time                   23:59:32.064
                    istat = sscanf(line, " Time %d:%d:%lf", &EventTime.hour, &EventTime.min, &EventTime.sec);
                    //printf("DEBUG: TEXNET_BULLETIN Time: %2.2d:%2.2d:%05.2f\n", EventTime.hour, EventTime.min, EventTime.sec);
                }

                // find preferred magnitude line
                //        ML        0.99 +/- 0.11   4 preferred
                ifound = 1;
                while (strstr(line, "preferred") == NULL) {
                    // read next line
                    cstat = fgets(line, MAXLINE_LONG, fp_obs);
                    //printf("DEBUG: TEXNET_BULLETIN preferred magnitude: <%s>\n", line);
                    //printf("3 %s", line);
                    if (cstat == NULL)
                        return (OBS_FILE_END_OF_INPUT);
                    if (strncmp(line, "Event:", 6) == 0) {
                        // end of event
                        return (OBS_FILE_END_OF_EVENT);
                    }
                    if (strstr(line, "Phase") != NULL) {
                        ifound = 0;
                        break;
                    }
                }
                if (ifound) {
                    // read magnitude
                    //        ML        0.99 +/- 0.11   4 preferred
                    istat = sscanf(line, "%*s %lf %*s %*f %d", &phypo->amp_mag, &phypo->num_amp_mag);
                    //printf("DEBUG: TEXNET_BULLETIN preferred magnitude: %lf %d\n", phypo->amp_mag, phypo->num_amp_mag);

                }

                // find phases lines
                //    sta  net   dist azi  phase   time         res     wt  sta
                while (strncmp(line, "    sta", 7) != 0) {
                    // read next line
                    cstat = fgets(line, MAXLINE_LONG, fp_obs);
                    //printf("3 %s", line);
                    if (cstat == NULL)
                        return (OBS_FILE_END_OF_INPUT);
                    if (strncmp(line, "Event:", 6) == 0) {
                        // end of event
                        return (OBS_FILE_END_OF_EVENT);
                    }
                }
                // read next line/
                cstat = fgets(line, MAXLINE_LONG, fp_obs);
                //printf("4 %s", line);
                if (cstat == NULL)
                    return (OBS_FILE_END_OF_INPUT);
                if (strncmp(line, "Event:", 6) == 0) {
                    // end of event
                    return (OBS_FILE_END_OF_EVENT);
                }
            } else {
                return (OBS_FILE_END_OF_EVENT); // blank line or ignored line
            }
        }


        /* read phase arrival input */
        /*
    sta  net   dist azi  phase   time         res     wt  sta
    PECS  TX    0.2 213  P       15:18:57.657  -0.1 M  1.5  PECS
    PECS  TX    0.2 213  S       15:19:00.866   0.4 M  1.4  PECS
    PB02  TX    0.2 122  S       15:19:01.446  -0.5 M  1.3  PB02
    HTMS  SC    0.9   6  P       15:19:12.251   0.1 M  1.5  HTMS
    ALPN  TX    1.2 175  S       15:19:29.664  -2.9 MX 0.0  ALPN
    ODSA  TX    1.2  60  P       15:19:19.587   3.0 MX 0.1  ODSA
    ODSA  TX    1.2  60  S       15:19:32.332  -0.6 M  1.2  ODSA
    ...
         */

        // read phase time reading
        istat = sscanf(line, "%s %*s %*f %*f %s %d:%d:%lf %*f %s",
                arrival->label, arrival->phase, &arrival->hour, &arrival->min, &arrival->sec, chrtmp);
        //istat = ReadFortranString(line, 5, 5, arrival->label);
        TrimString(arrival->label);
        //istat += ReadFortranString(line, 26, 7, arrival->phase);
        TrimString(arrival->phase);
        //istat += ReadFortranInt(line, 34, 2, &arrival->hour);
        // set weight, check for P - reset weight
        if (strcmp(arrival->phase, "P") == 0)
            arrival->quality = 0;
        else
            arrival->quality = 1;

        //printf("4A %s %s %d:%d:%lf %s\n",
        //        arrival->label, arrival->phase, arrival->hour, arrival->min, arrival->sec, chrtmp);
        if (istat != 6) {
            return (OBS_FILE_INVALID_PHASE);
        }
        //printf("4B %s %s %d:%d:%lf %s\n",
        //        arrival->label, arrival->phase, arrival->hour, arrival->min, arrival->sec, chrtmp);

        arrival->year = EventTime.year;
        arrival->month = EventTime.month;
        arrival->day = EventTime.day;

        // check for arrival time over day boundary relative to event time
        if (arrival->hour == 0 && EventTime.hour == 23) {
            arrival->hour = 24;
            //printf("DEBUG: TEXNET_BULLETIN Time past day bndry -> %2.2d:%2.2d:%05.2f\n", arrival->hour, arrival->min, arrival->sec);
        }


        // convert quality to error
        Qual2Err(arrival);
        // multiplication of error if not first arrival P
        //if (0 && strstr("P$Pg$Pn$Pb$P0$P1$PKP$PKPdf", arrival->phase) == NULL) {
        //    arrival->error *= 2.0;
        //}

        if (strcmp(ftype_obs, "TEXNET_BULLETIN__ZERO_WT_X") == 0) {
            // check if phase is excluded for location
            //istat = ReadFortranString(line, 54, 1, chrtmp);
            if (strlen(chrtmp) > 1 && chrtmp[1] == 'X') {
                arrival->apriori_weight = 0.0;
            }
        }


        return (istat);

    } else if (strcmp(ftype_obs, "INGV_AUTOPICKS") == 0) {

        /* example:
                        ...
                        FUORNHHE 2005/04/30 23:22:15.29
                        FUORNHHZ 2005/04/30 23:23:50.97
                        FUORNHHE 2005/04/30 23:27:59.99
                        FUORNHHE 2005/04/30 23:27:59.99  4715 63.21 238.72
                        TOLF HHE 2005/04/30 23:38:48.56
                        TOLF HHE 2005/04/30 23:38:48.56  1770 0.06 26.75
                        ...
         */

        /* read line */
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);
        /* check for end of event (assumes empty lines between event) */
        if (LineIsBlank(line)) {
            /* end of event */
            return (OBS_FILE_END_OF_EVENT);
        }

        /* read phase arrival input */
        istat = ReadFortranString(line, 1, 5, arrival->label);
        TrimString(arrival->label);
        istat += ReadFortranString(line, 6, 2, arrival->inst);
        istat += ReadFortranString(line, 8, 1, arrival->comp);
        snprintf(arrival->phase, sizeof(arrival->phase), "%s", "P");
        arrival->quality = 0;
        istat += ReadFortranInt(line, 10, 4, &arrival->year);
        istat += ReadFortranInt(line, 15, 2, &arrival->month);
        istat += ReadFortranInt(line, 18, 2, &arrival->day);
        istat += ReadFortranInt(line, 21, 2, &arrival->hour);
        istat += ReadFortranInt(line, 24, 2, &arrival->min);
        istat += ReadFortranReal(line, 27, 5, &arrival->sec);

        if (istat != 9) {
            return (OBS_FILE_END_OF_EVENT);
        }


        /* convert quality to error */
        Qual2Err(arrival);

        return (istat);
    } else if (strcmp(ftype_obs, "NCEDC_UCB") == 0) {

        /* example:
                        ...
                        19990818 010619.0000 +37.8950 -122.7000 010.0000 0.00 000 4.87 008 4.99 040 0.00 000 0.000e+00 000 099 180 009.47 00.0300 00.2000 00.2000 00.4000 00.1700 x F
                        $COM Bolinas, 19 km WSW of San Rafael, CA
                        $COM BRK: Mo=1.23E+23 dyne-cm, Mw=4.7.  Minor damage reported in Bolinas and San Rafael; felt from Santa Rosa in the north to San Jose in the
                        $COM south.
                        $AMP CMSB CL   2 WAS  407.40 001.10 0.00e+00
                        $AMP CMSB CL   3 WAS  596.20 000.77 0.00e+00
                        $PHS CMSB CL   1 i P        d 19990818 010625.9300 0039.535 093.58 000.00 000.00 0.00e+00 Y
                        $PHS CMSB CL   3 e S        x 19990818 010630.8700 0039.535 093.58 000.00 000.00 0.00e+00 N
                        $AMP CRQB CL   2 WAS  972.90 000.76 0.00e+00
                        $AMP CRQB CL   3 WAS  1122.00 000.82 0.00e+00
                        $PHS CRQB CL   1 i P        d 19990818 010627.5600 0045.409 066.72 000.00 000.00 0.00e+00 N
                        $PHS CRQB CL   3 i S        x 19990818 010633.6500 0045.409 066.72 000.00 000.00 0.00e+00 N
         */
        /* read line */
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);
        /* check for end of event (assumes empty lines between event) */
        if (LineIsBlank(line) || strncmp(line, "$END", 4) == 0) {
            /* end of event */
            return (OBS_FILE_END_OF_EVENT);
        }


        // find next phase line
        while (strncmp(line, "$PHS", 4) != 0) {
            // read next line
            cstat = fgets(line, MAXLINE_LONG, fp_obs);
            if (cstat == NULL || strncmp(line, "$END", 4) == 0)
                return (OBS_FILE_END_OF_EVENT);
        }


        /* read phase arrival input */
        /* example:
                        $PHS CRQB CL   1 i P        d 19990818 010627.5600 0045.409 066.72 000.00 000.00 0.00e+00 N
                        $PHS CRQB CL   3 i S        x 19990818 010633.6500 0045.409 066.72 000.00 000.00 0.00e+00 N
                        $PHS MHC  VSP  Z i P        c 19990818 010636.8500 0111.833 123.00 000.00 000.00 0.00e+00 N
                        $PHS MHC  VSP  E i S        x 19990818 010651.9000 0111.833 123.00 000.00 000.00 0.00e+00 N
         */
        /*
                        $PHS ZSP  TS13 Z i P        d 19890212 195041.9000 0032.674 068.58 000.00 000.00 0.00e+00 N
                        $PHS ZSP  TS13 Z i S        x 19890212 195047.4000 0032.674 068.58 000.00 000.00 0.00e+00 N
                        $PHS NFI  x    Z 0 P        c 19890212 195042.2500 0038.265 246.22 000.00 000.00 0.00e+00 Y
                        $PHS JEG  x    Z 0 P        d 19890212 195042.4500 0038.037 161.07 000.00 000.00 0.00e+00 N
                        $PHS PCC  B14  Z i P        c 19890212 195042.8000 0232.545 146.11 000.00 000.00 0.00e+00 Y
                        $PHS PCC  B14  Z i S        x 19890212 195048.7000 0232.545 146.11 000.00 000.00 0.00e+00 N
                        $PHS NPR  x    Z 3 P        c 19890212 195044.1300 0040.492 295.89 000.00 000.00 0.00e+00 N
         */
        istat = ReadFortranString(line, 6, 4, arrival->label);
        TrimString(arrival->label);
        istat += ReadFortranString(line, 11, 4, arrival->inst);
        TrimString(arrival->inst);
        istat += ReadFortranString(line, 16, 1, arrival->comp);
        TrimString(arrival->comp);
        istat += ReadFortranString(line, 18, 1, arrival->onset);
        if (arrival->onset[0] == '?')
            arrival->onset[0] = ARRIVAL_NULL_CHR;
        istat += ReadFortranString(line, 20, 8, arrival->phase);
        TrimString(arrival->phase);
        istat += ReadFortranString(line, 29, 1, arrival->first_mot);
        if (arrival->first_mot[0] == 'x')
            arrival->first_mot[0] = ARRIVAL_NULL_CHR;
        // set weight, check for P - reset weight
        if (arrival->onset[0] == 'x' || arrival->onset[0] == 'X')
            arrival->quality = 1;
        if (arrival->onset[0] == 'i')
            arrival->quality = 0;
        if (arrival->onset[0] == 'e')
            arrival->quality = 1;
        if (arrival->onset[0] == '?')
            arrival->quality = 1;
        if (arrival->onset[0] == '0')
            arrival->quality = 0;
        if (arrival->onset[0] == '1')
            arrival->quality = 1;
        if (arrival->onset[0] == '2')
            arrival->quality = 2;
        if (arrival->onset[0] == '3')
            arrival->quality = 3;
        if (arrival->onset[0] == '4')
            arrival->quality = 4;
        if (strncmp(arrival->phase, "P", 1) != 0) {
            arrival->quality += 1;
        }
        istat += ReadFortranInt(line, 31, 4, &arrival->year);
        istat += ReadFortranInt(line, 35, 2, &arrival->month);
        istat += ReadFortranInt(line, 37, 2, &arrival->day);
        istat += ReadFortranInt(line, 40, 2, &arrival->hour);
        istat += ReadFortranInt(line, 42, 2, &arrival->min);
        istat += ReadFortranReal(line, 44, 7, &arrival->sec);

        if (istat != 12) {
            return (OBS_FILE_END_OF_EVENT);
        }


        /* convert quality to error */
        Qual2Err(arrival);

        return (istat);
    } else if (strcmp(ftype_obs, "HYPOCENTER") == 0) {

        /* example:
                        67  827 1632 12.4 L
                        ZAK0SZ  PG      1632  33.30
                        ZAK0SZ  SG      1632  50.80
                        IRK0SZ  PG      1632  29.50
                        IRK0SZ  SG      1632  43.90
                        KHT0SZ  PG      1632  34.00
         */

        /* read line */
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);

        /* check for end of event (assumes empty lines between event) */
        if (LineIsBlank(line)) {
            /* end of event */
            return (OBS_FILE_END_OF_EVENT);
        }


        /* read event hypocenter line */
        if (!in_hypocenter_event) {
            /* read hypocenter time */
            istat = ReadFortranInt(line, 1, 5, &EventTime.year);
            istat += ReadFortranInt(line, 7, 2, &EventTime.month);
            istat += ReadFortranInt(line, 9, 2, &EventTime.day);
            if (istat == 3) {
                if (EventTime.year < 2000)
                    EventTime.year += 1900;
                /* !! assume any 2 digit year is 1900-1999 */
                in_hypocenter_event = 1;
                /* read next line */
                cstat = fgets(line, MAXLINE_LONG, fp_obs);
                if (cstat == NULL)
                    return (OBS_FILE_END_OF_INPUT);
            } else {
                return (OBS_FILE_END_OF_EVENT);
            }
        }


        /* read phase arrival input */

        istat = ReadFortranString(line, 3, 3, arrival->label);
        istat += ReadFortranInt(line, 6, 1, &arrival->quality); // !!! is this field quality???
        istat += ReadFortranString(line, 7, 2, arrival->comp);
        istat += ReadFortranString(line, 11, 6, arrival->phase);
        istat += ReadFortranInt(line, 19, 2, &arrival->hour);
        istat += ReadFortranInt(line, 21, 2, &arrival->min);
        // format of seconds varies, so cannot use fixed fortran read
        //istat += ReadFortranReal(line, 25, 5, &arrival->sec);
        istat += sscanf(line, "%*s %*s %*s %lf", &arrival->sec);

        if (istat != 7) {
            return (OBS_FILE_END_OF_EVENT);
        }


        arrival->year = EventTime.year;
        arrival->month = EventTime.month;
        arrival->day = EventTime.day;


        /* convert quality to error */
        Qual2Err(arrival);

        return (istat);
    } else if (strcmp(ftype_obs, "HYPODD_PHA") == 0 ||
            strcmp(ftype_obs, "HYPODD_PHA_S_QUAL_PLUS_1") == 0) {

        /* example:
        # 1985  1 24  2 19 58.71  37.8832 -122.2415    9.80 1.40  0.15  0.51  0.02      38542
                        NCCSP       2.850  -1.000   P
                        NCCBW       3.430  -1.000   P
                        NCCMC       2.920   0.200   P
                        NCCAI       3.440  -1.000   P
                        NCCBR       3.940   1.000   P
                        NCCLC       4.610   0.500   P
                        NCJPR       4.470   0.100   P
                        NCNLH       6.000   0.200   P
                        NCCSH       6.130   0.100   P
                        NCNTA       5.630   0.100   P
                        NCJMG       6.240   0.200   P
         */

        if (in_hypocenter_event) {
            /* check for end of event (assumes no blanks after last phase) */
            chr = fgetc(fp_obs);
            if (chr == EOF) {
                return (OBS_FILE_END_OF_INPUT);
            } else if (chr == '#') {
                // end of event
                ungetc(chr, fp_obs);
                return (OBS_FILE_END_OF_EVENT);
            } else {
                ungetc(chr, fp_obs);
            }
        }

        /* read line */
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);


        /* read event hypocenter line */
        if (!in_hypocenter_event) {
            /* read hypocenter time */
            // DATE, TIME, LAT, LON, DEP, MAG, EH, EV, RMS, ID
            // 2025  02  28   23    59    33.212    36.6105    25.6038    0.00     0.768     0.0     0.0    0.0    34423
            // 20250314 AJL - added reading of MAG to hypo
            istat = sscanf(line, "# %d %d %d %d %d %lf %*f %*f %*f %lf %*f %*f %*f %ld",
                    &EventTime.year, &EventTime.month, &EventTime.day,
                    &EventTime.hour, &EventTime.min, &EventTime.sec, &phypo->amp_mag, &EventID);
            if (istat == 8) {
                in_hypocenter_event = 1;
                /* read next line */
                cstat = fgets(line, MAXLINE_LONG, fp_obs);
                if (cstat == NULL)
                    return (OBS_FILE_END_OF_INPUT);
            } else {
                return (OBS_FILE_END_OF_EVENT);
            }
        }


        // read phase arrival input

        istat = sscanf(line, "%s %lf %lf %s", chrtmp, &ttime, &weight, arrival->phase);
        if (istat != 4)
            return (OBS_FILE_END_OF_EVENT);

        // trim label (to agree with Alberto's 3 char convention)
        //		if (strlen(chrtmp) > 4 && strncmp(chrtmp, "NC", 2) == 0)
        //			snprintf(arrival->label, sizeof(arrival->label), "%s", chrtmp + 2);
        //		else
        snprintf(arrival->label, sizeof(arrival->label), "%s", chrtmp);

        arrival->hour = EventTime.hour;
        arrival->min = EventTime.min;
        arrival->sec = EventTime.sec + ttime;

        arrival->year = EventTime.year;
        arrival->month = EventTime.month;
        arrival->day = EventTime.day;

        // 20110620 AJL - preserve event id if available
        arrival->dd_event_id_1 = EventID;
        snprintf(phypo->public_id, sizeof (phypo->public_id), "%ld", EventID); // 20250221 AJL - added

        // convert weight to quality (following Appendix C in hypoDD.pdf)
        weight = fabs(weight);
        if (weight > 0.99)
            arrival->quality = 0;
        else if (weight > 0.49)
            arrival->quality = 1;
        else if (weight > 0.19)
            arrival->quality = 2;
        else if (weight > 0.09)
            arrival->quality = 3;
        else
            arrival->quality = 4;
        // increase S uncertainty
        if (IsPhaseID(arrival->phase, "S") && strcmp(ftype_obs, "HYPODD_PHA_S_QUAL_PLUS_1") == 0) {
            arrival->quality += 1;
        }
        // convert quality to error
        Qual2Err(arrival);

        return (istat);
    }// DD
    else if (strcmp(ftype_obs, "HYPODD_CC") == 0 || strcmp(ftype_obs, "HYPODD_CT") == 0
            || strcmp(ftype_obs, "HYPODD_") == 0) {

        /* HYPODD_CC example:
        #    28136    46442     -0.174000
                        NCCMO    -0.038245201    0.63    P
                        NCCRP    -0.062200069    0.78    P
                        NCJBG    -0.020649910    0.60    P
                        NCJPS    -0.011000156    0.50    P
         */
        /* OTC - Origin time correction relative to the event origin time reported
                        in the catalog data. If not known for all available data,
                        then cross correlation data can not be used in combination
                        with catalog data. Set OTC to -999 if not known for individual
                        observations; Set to 0.0 if cross-correlation and catalog origin
                        times are identical, or if only cross-correlation data is used.
         */
        /* HYPODD_CT example:
        #     38542     38520
                        NCCBW     3.430   3.430 1.0000 P
                        NCCSP     2.850   2.840 1.0000 P
                        NCNLN     8.950   8.950 1.0000 P
                        NCCAI     3.440   3.420 1.0000 P
         */


        /* read line */
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);


        /* read events line */
        if (line[0] == '#') {
            hypo_cc_flag = 1;
            istat = sscanf(line, "# %ld %ld %lf",
                    &dd_event_id_1, &dd_event_id_2, &dd_otime_corr);
            if (dd_event_id_1 == dd_event_id_2)
                printf("ERROR: dd_event_id_1 == &dd_event_id_2  %ld %ld\n", dd_event_id_1, dd_event_id_2);
            if (istat == 2)
                hypo_cc_flag = 0;

            if (istat >= 2) {
                in_hypocenter_event = 1;
                /* read next line */
                cstat = fgets(line, MAXLINE_LONG, fp_obs);
                if (cstat == NULL)
                    return (OBS_FILE_END_OF_INPUT);
            } else {
                printf("ERROR: format error: %s\n", line);
                return (OBS_FILE_FORMAT_ERROR);
            }
        }


        // read phase arrival input

        arrival->xcorr_flag = hypo_cc_flag;

        if (hypo_cc_flag) { // HYPODD_CC

            istat = sscanf(line, "%s %lf %lf %s",
                    arrival->label, &arrival->dd_dtime, &arrival->weight, arrival->phase);
            //printf("%s %lf %lf %s ignore:%d\n", arrival->label, arrival->dd_dtime, arrival->weight, arrival->phase, arrival->flag_ignore);
            if (istat != 4) {
                printf("ERROR: HYPODD_CC format error: %s\n", line);
                return (OBS_FILE_FORMAT_ERROR);
            }

            arrival->dd_event_id_1 = dd_event_id_1;
            arrival->dd_event_id_2 = dd_event_id_2;
            // incorporate OTC into dd_dtime when reading dt file
            arrival->dd_dtime -= dd_otime_corr;
            arrival->error = 0.0;
            snprintf(arrival->error_type, sizeof(arrival->error_type), "%s", "XCC");

        } else { // HYPODD_CT

            istat = sscanf(line, "%s %lf %lf %lf %s",
                    arrival->label, &tt_sta1, &tt_sta2, &arrival->weight, arrival->phase);
            //printf("%s %lf %lf %lf %s ignore:%d\n", arrival->label, tt_sta1, tt_sta2, arrival->weight, arrival->phase, arrival->flag_ignore);
            if (istat != 5) {
                printf("ERROR: HYPODD_CT format error: %s\n", line);
                return (OBS_FILE_FORMAT_ERROR);
            }

            arrival->dd_event_id_1 = dd_event_id_1;
            arrival->dd_event_id_2 = dd_event_id_2;
            arrival->dd_dtime = tt_sta1 - tt_sta2;
            arrival->error = 0.0;
            snprintf(arrival->error_type, sizeof(arrival->error_type), "%s", "CAT");

        }

        // add noise
        //arrival->dd_dtime += 0.02 * get_rand_double(-1.0, 1.0);
        //

        return (istat);
    } else if (strcmp(ftype_obs, "RENASS_WWW") == 0) {

        /* read next line */
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);


        /* read formatted P arrival input */
        istat = sscanf(line, "%s %s %dh %dmn %lfsec",
                arrival->label, arrival->phase, &arrival->hour,
                &arrival->min, &arrival->sec);

        if (istat != 5) {
            line[0] = '\0';
            return (OBS_FILE_SKIP_INPUT_LINE);
        }

        arrival->quality = 0;
        arrival->year = 1900;
        arrival->month = 01;
        arrival->day = 01;


        /* convert quality to error */
        Qual2Err(arrival);

        return (istat);
    } else if (strcmp(ftype_obs, "RENASS_DEP") == 0) {

        /* read next line */
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);


        /* read formatted arrival input */
        istat = ReadFortranString(line, 1, 4, arrival->label);
        TrimString(arrival->label);
        istat += ReadFortranString(line, 20, 4, arrival->phase);
        TrimString(arrival->phase);
        istat += ReadFortranInt(line, 25, 2, &arrival->hour);
        istat += ReadFortranInt(line, 28, 2, &arrival->min);
        istat += ReadFortranReal(line, 30, 5, &arrival->sec);

        if (istat != 5) {
            line[0] = '\0';
            return (OBS_FILE_SKIP_INPUT_LINE);
        }


        arrival->quality = 0;
        arrival->year = EventTime.year;
        arrival->month = EventTime.month;
        arrival->day = EventTime.day;
        if (arrival->hour != EventTime.hour)
            nll_puterr("WARNING: filename and arrival hours do not match.");
        if (arrival->min != EventTime.min)
            nll_puterr("WARNING: filename and arrival minutes do not match.");


        /* convert quality to error */
        Qual2Err(arrival);

        return (istat);
    } else if (strcmp(ftype_obs, "SEISAN") == 0 || strcmp(ftype_obs, "NORDIC") == 0) {

        /* read next line */
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);

        // check for blank line
        if (LineIsBlank(line))
            return (OBS_FILE_END_OF_EVENT);

        // check for line without blank first character
        if (line[0] != ' ')
            return (OBS_FILE_END_OF_EVENT);



        // check for date info
        if (!date_saved) {
            // Check whether it is a hypocenter line
            // 2023 0220 1219 05.0 L  37.359  37.066  5.3  AFD  6 0.1 1.6LAFD                1
            if (line[79] != '1') {
                nll_puterr2("ERROR: bad SEISAN hypocentre line", line);
                return (OBS_FILE_FORMAT_ERROR);
            }
            istat = ReadFortranInt(line, 2, 4, &year_save);
            /* in SEISAN data the year is always fully specified. The correction below is susceptiple to millenium bug */
            /* if (year_save < 500) */
            /* 	year_save += 1900; */
            istat += ReadFortranInt(line, 7, 2, &month_save);
            istat += ReadFortranInt(line, 9, 2, &day_save);
            //printf("SEISAN_DATE Read: istat %d   -  %d %d %d\n", istat, year_save, month_save, day_save);
            if (istat == 3 && IsGoodDate(year_save, month_save, day_save)) {
                date_saved = 1;
                // read magnitude
                istat += ReadFortranReal(line, 57, 3, &phypo->amp_mag);
                phypo->num_amp_mag = 1;
                return (OBS_FILE_SKIP_INPUT_LINE);
            }
        }


        // ignore all lines which are not phase observation lines
        if (line[79] != ' ' && line[79] != '4') {
            return (OBS_FILE_SKIP_INPUT_LINE);
        }

        // ignore amplitude reading lines (AFAD)
        // KAHM HE  IAML    1219 12.00      126.9   0.2                          18.9
        if (strncmp(line + 10, "IAML", 4) == 0) {
            return (OBS_FILE_SKIP_INPUT_LINE);
        }
        // ignore SPEC reading lines with bad format (AFAD)
        // SPEC SABUHH Z T140552 K 0.040 GD625.0 VP 6.20 DE 3.00 Q0200.0 QA 0.60 VS      3
        // SPEC EKARHH Z MO18.47 ST45.73 OM 5.79 f0      R  6.57 AL      WI 104.9 MW  6.2 3
        //printf("len %ld\n", strnlen(line, MAXLINE_LONG));
        if (strnlen(line, MAXLINE_LONG) > 81 && strncmp(line + 1, "SPEC", 4) == 0) {
            return (OBS_FILE_SKIP_INPUT_LINE);
        }

        // read formatted P arrival input
        // KHMR HZ IP       1219 08.60                                           14.2 264
        istat = ReadFortranString(line, 2, 5, arrival->label);
        TrimString(arrival->label);
        istat += ReadFortranString(line, 7, 1, arrival->inst);
        TrimString(arrival->inst);
        if (strlen(arrival->inst) < 1) { // 20231108 AJL - make sure empty fields are set to "?"
            strncpy(arrival->inst, "?", sizeof (arrival->inst));
        }
        istat += ReadFortranString(line, 8, 1, arrival->comp);
        TrimString(arrival->comp);
        if (strlen(arrival->comp) < 1) { // 20231108 AJL - make sure empty fields are set to "?"
            strncpy(arrival->comp, "?", sizeof (arrival->comp));
        }
        istat += ReadFortranString(line, 10, 1, arrival->onset);
        TrimString(arrival->onset);
        if (strlen(arrival->onset) < 1) { // 20231108 AJL - make sure empty fields are set to "?"
            strncpy(arrival->onset, "?", sizeof (arrival->onset));
        }
        istat += ReadFortranString(line, 11, 4, arrival->phase);
        TrimString(arrival->phase);
        if (strlen(arrival->phase) < 1) { // 20231108 AJL - make sure empty fields are set to "?"
            strncpy(arrival->phase, "?", sizeof (arrival->phase));
        }
        istat += ReadFortranInt(line, 15, 1, &arrival->quality);
        istat += ReadFortranString(line, 17, 1, arrival->first_mot);
        TrimString(arrival->first_mot);
        if (strlen(arrival->first_mot) < 1) { // 20231108 AJL - make sure empty fields are set to "?"
            strncpy(arrival->first_mot, "?", sizeof (arrival->first_mot));
        }
        istat += ReadFortranInt(line, 19, 2, &arrival->hour);
        istat += ReadFortranInt(line, 21, 2, &arrival->min);
        istat += ReadFortranReal(line, 23, 6, &arrival->sec);
        istat += ReadFortranReal(line, 30, 4, &arrival->coda_dur);
        istat += ReadFortranReal(line, 34, 7, &arrival->amplitude);
        istat += ReadFortranReal(line, 42, 4, &arrival->period);
        if (date_saved) {
            arrival->year = year_save;
            arrival->month = month_save;
            arrival->day = day_save;
        } else {
            arrival->year = 1900;
            arrival->month = 01;
            arrival->day = 01;
        }

        /*printf(
        "SEISAN Read: istat %d   -  %s %s %s %d %s %d %d %d %d %d %lf %lf\n",
        istat, arrival->label, arrival->onset, arrival->phase, arrival->quality, arrival->first_mot, arrival->year, arrival->month, arrival->day, arrival->hour, arrival->min, arrival->sec, arrival->coda_dur);*/

        if (istat != 13) {
            nll_puterr2("ERROR: bad SEISAN phase line:", line);
            return (OBS_FILE_SKIP_INPUT_LINE);
        }

        /* convert quality to error */
        Qual2Err(arrival);

        return (istat);

    } else if (strcmp(ftype_obs, "SAFOD") == 0) {

        /* example:

                Location Header
                arrival 1
                arrival 2
                ...
                arrival n
                        (blank line signifying end of phases for this event)

        Here's an example:

        1994 259  3  0   0.0000   0.46899    1.72781  -0.7000      994Shot01
        MMCR  1994 259  3  0   0.330 P  -0.020   0.020
        MMEC  1994 259  3  0   0.720 P  -0.724   0.724
        MMEN  1994 259  3  0   0.950 P  -0.020   0.020


        The Header:
        The fields in the location header are: year, julian day, hour, minute, second, local x, local y, local z, event ID.
        The time fields give the current estimate of the origin time, the x,y,z coordinates are in a reference fram local to SAFOD. I will be providing a program to allow you to go between this system and UTM NAD 27 soon.

        The Phases:
        The fields in the phase line are: station code, year, julian day, hour, minute, second, phase designator (P or S), left uncertainty, right uncertainty.

        Often, the left uncertainty field is blank. Any information in colums to the right of the right uncertainty is superfluous (sometimes there will be a residual from a location run, for example).
         */

        // read line
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);
        // check for end of event (assumes empty lines between event)
        if (LineIsBlank(line)) {
            // end of event
            return (OBS_FILE_END_OF_EVENT);
        }

        // read phase arrival input
        istat = ReadFortranString(line, 1, 5, arrival->label);
        istat += ReadFortranInt(line, 7, 4, &arrival->year);
        istat += ReadFortranInt(line, 12, 3, &yearday);
        istat += ReadFortranInt(line, 16, 2, &arrival->hour);
        istat += ReadFortranInt(line, 19, 2, &arrival->min);
        istat += ReadFortranReal(line, 22, 7, &arrival->sec);
        istat += ReadFortranString(line, 30, 1, arrival->phase);
        istat += ReadFortranReal(line, 32, 7, &left_uncertainty);
        istat += ReadFortranReal(line, 40, 7, &right_uncertainty);

        TrimString(arrival->label);
        MonthDay(arrival->year, yearday, &arrival->month, &arrival->day);
        if (strcmp(arrival->phase, "0") == 0)
            snprintf(arrival->phase, sizeof(arrival->phase), "%s", "P");

        if (istat != 9) {
            if (nfirst) // may be header line
                return (OBS_FILE_SKIP_INPUT_LINE);
            nll_puterr2("ERROR: bad SAFOD phase line:", line);
            return (OBS_FILE_END_OF_EVENT);
        }

        // set error fields
        snprintf(arrival->error_type, sizeof(arrival->error_type), "%s", "GAU");
        arrival->error = right_uncertainty;
        /*
        if (strcmp(arrival->phase, "P") == 0)
        arrival->error = 0.1;
        if (strcmp(arrival->phase, "S") == 0)
        arrival->error = 0.2;
         */

        return (istat);

    } else if (strcmp(ftype_obs, "INERIS") == 0) {

        /* example:

        #      051013_07425894-TG13  864088 194914 96    -9999
        TGe3D   P  2005 10 13 07 42  2.718E-01  0.005
        TGe1DB  P  2005 10 13 07 42  2.773E-01  0.005
        TGe1DS  P  2005 10 13 07 42  2.786E-01  0.005
        TLi3D   P  2005 10 13 07 42  1.766E-01  0.005
        TLi1DB  P  2005 10 13 07 42  1.828E-01  0.005

         */

        // read line
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);
        // check for end of event (assumes empty lines between event)
        if (LineIsBlank(line)) {
            // end of event
            return (OBS_FILE_END_OF_EVENT);
        }

        // read phase arrival input
        istat = sscanf(line, "%s %s %d %d %d %d %d %lf %lf",
                arrival->label, arrival->phase, &arrival->year, &arrival->month, &arrival->day,
                &arrival->hour, &arrival->min, &arrival->sec, &arrival->error);

        if (istat != 9) {
            if (nfirst) // may be header line
                return (OBS_FILE_SKIP_INPUT_LINE);
            nll_puterr2("ERROR: bad INERIS phase line:", line);
            return (OBS_FILE_SKIP_INPUT_LINE);
        }

        // set error fields
        snprintf(arrival->error_type, sizeof(arrival->error_type), "%s", "GAU");

        return (istat);

    } else if (strcmp(ftype_obs, "EW_PKFILTER") == 0) {

        /* example:

        Pass: 10  6 73 4771 FG4  IVSHZ U3  20070323133511.42      33      37      27
        Pass: 10  6 73 4772 LAV9 IVSHZ D2  20070323133517.18      50      20      23
        Pass: 10  6 73 4773 LAV9 IVSHZ D3  20070323133533.68      35      43      44
        Pass: 10  6 73 4774 MSAG IVHHZ  2  20070323133538.59     376     465     346
        Pass: 10  6 73 4775 MSC  OVV   D3  20070323133540.00      62     202     211
        Pass: 10  6 73 4776 MSC  OVV    1  20070323133549.24     482    1754     270
        Pass: 10  6 73 4777 NRCA IVSHZ D0  20070323133542.24     783     259     638
        Pass: 10  6 73 4778 BDI  IVHHZ  2  20070323133544.54     575     265      93
        Pass: 10  6 73 4779 ARVD IVEHZ  3  20070323133553.28      94      45     135
        Pass: 10  6 73 4780 AOI  IVBHZ U2  20070323133554.52     187     154     283
        Pass: 10  6 73 4781 SNTG IVSHZ  2  20070323133549.12      40      80      69

         */

        // read line
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);
        // check for end of event (assumes empty lines between event)
        if (LineIsBlank(line)) {
            // end of event
            return (OBS_FILE_END_OF_EVENT);
        }

        // read phase arrival input
        istat = ReadFortranString(line, 21, 5, arrival->label);
        istat += ReadFortranString(line, 26, 2, arrival->network);
        istat += ReadFortranString(line, 28, 3, arrival->comp);
        istat += ReadFortranString(line, 32, 1, arrival->first_mot);
        istat += ReadFortranInt(line, 33, 1, &arrival->quality);
        istat += ReadFortranInt(line, 36, 4, &arrival->year);
        istat += ReadFortranInt(line, 40, 2, &arrival->month);
        istat += ReadFortranInt(line, 42, 2, &arrival->day);
        istat += ReadFortranInt(line, 44, 2, &arrival->hour);
        istat += ReadFortranInt(line, 46, 2, &arrival->min);
        istat += ReadFortranReal(line, 48, 6, &arrival->sec);

        if (istat != 11) {
            nll_puterr2("ERROR: bad EW_PKFILTER phase line:", line);
            return (OBS_FILE_END_OF_EVENT);
        }

        TrimString(arrival->label);
        TrimString(arrival->network);
        TrimString(arrival->comp);

        snprintf(arrival->phase, sizeof(arrival->phase), "%s", "P");

        // convert quality to error
        Qual2Err(arrival);
        // set error fields
        snprintf(arrival->error_type, sizeof(arrival->error_type), "%s", "GAU");

        return (istat);

    } else if (strcmp(ftype_obs, "EW_PTWC_HAWAII") == 0) {

        /* example:

        10 70  9 3888 MLOA PTHNZ  1  20150124165003.00     107     103 68
        10 70  9 3889 SPDD HVEHZ D1  20150124165021.21     141     154 88
        10 70  9 3890 MLOA PTHNZ  2  20150124165039.47      55      71 11
        10 70  9 3891 POLD HVEHZ U1  20150124165049.40     155     363 312
        10 70  9 3892 KLUD HVEHZ U3  20150124165043.35      26      15 27
        10 70  9 3893 HTCD HVEHZ U0  20150124165047.15     754     858 574

         */

        // read line
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);
        // check for end of event (assumes empty lines between event)
        if (LineIsBlank(line)) {
            // end of event
            return (OBS_FILE_END_OF_EVENT);
        }

        // read phase arrival input
        istat = ReadFortranString(line, 16, 5, arrival->label);
        istat += ReadFortranString(line, 21, 2, arrival->network);
        istat += ReadFortranString(line, 23, 3, arrival->comp);
        istat += ReadFortranString(line, 27, 1, arrival->first_mot);
        istat += ReadFortranInt(line, 28, 1, &arrival->quality);
        istat += ReadFortranInt(line, 31, 4, &arrival->year);
        istat += ReadFortranInt(line, 35, 2, &arrival->month);
        istat += ReadFortranInt(line, 37, 2, &arrival->day);
        istat += ReadFortranInt(line, 39, 2, &arrival->hour);
        istat += ReadFortranInt(line, 41, 2, &arrival->min);
        istat += ReadFortranReal(line, 43, 6, &arrival->sec);

        if (istat != 11) {
            nll_puterr2("ERROR: bad EW_PTWC_HAWAII phase line:", line);
            return (OBS_FILE_END_OF_EVENT);
        }

        TrimString(arrival->label);
        TrimString(arrival->network);
        TrimString(arrival->comp);

        snprintf(arrival->phase, sizeof(arrival->phase), "%s", "P");

        // convert quality to error
        Qual2Err(arrival);
        // set error fields
        snprintf(arrival->error_type, sizeof(arrival->error_type), "%s", "GAU");

        return (istat);

    } else if (strcmp(ftype_obs, "GEOFON") == 0) {

        /* example:
          Northern Sumatra, Indonesia  M=6.3  2008/01/22  17:14:57.4  1.10 N  97.49 E   31 km

          Stat  Net   Date       Time          Amp    Per   Res  Dist  Az mb  ML  mB
          GSI   GE  08/01/22  17:15:04.5         0.0  0.0  -0.6   0.2  22 0.0 5.6 0.0
          PPI   IA  08/01/22  17:15:46.5         0.0  0.0  -0.5   3.3 118 0.0 6.1 0.0
          PDSI  IA  08/01/22  17:15:48.8         0.0  0.0  -2.2   3.6 124 0.0 6.2 0.0
          IPM   MY  08/01/22  17:16:10.5         0.0  0.0   1.8   4.9  46 0.0 6.6 0.0
          BSI   IA  08/01/22  17:16:06.1         0.0  0.0  -2.8   4.9 333 0.0 6.4 0.0
          RGRI  IA  08/01/22  17:16:11.7      1867.5  0.8   0.5   5.1 106 6.1 6.5 0.0
          KUM   MY  08/01/22  17:16:14.4      2527.4  1.1   0.9   5.2  37 6.1 6.5 0.0
          NTU   MS  08/01/22  17:16:27.4        34.5  0.9   0.6   6.2  88 4.5 5.6 0.0
          KOM   MY  08/01/22  17:16:30.4       549.7  1.3   0.9   6.4  84 5.5 6.1 6.3
                        ...
         */

        /* read line */
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);
        /* check for end of event (assumes empty lines between event) */
        if (LineIsBlank(line)) {
            /* end of event */
            return (OBS_FILE_END_OF_EVENT);
        }

        /* read phase arrival input */
        istat = ReadFortranString(line, 3, 5, arrival->label);
        TrimString(arrival->label);
        istat += ReadFortranString(line, 9, 2, arrival->network);
        TrimString(arrival->network);
        snprintf(arrival->phase, sizeof(arrival->phase), "%s", "P");
        arrival->quality = 0;
        istat += ReadFortranInt(line, 13, 2, &arrival->year);
        arrival->year += 2000;
        istat += ReadFortranInt(line, 16, 2, &arrival->month);
        istat += ReadFortranInt(line, 19, 2, &arrival->day);
        istat += ReadFortranInt(line, 23, 2, &arrival->hour);
        istat += ReadFortranInt(line, 26, 2, &arrival->min);
        istat += ReadFortranReal(line, 29, 5, &arrival->sec);

        if (istat != 8) {
            return (OBS_FILE_END_OF_EVENT);
        }


        /* convert quality to error */
        Qual2Err(arrival);

        return (istat);

    } else if (strcmp(ftype_obs, "KOERI_EQEVENTS") == 0) {

        // http://www.koeri.boun.edu.tr/sismo/2/latest-earthquakes/automatic-solutions/

        /* example:

        -
        -   Ozet-Sum
        -   Varislar-Arrivals
        -   Kalite-Quality
        -   Istasyonlar-Stations

        Deprem ID-Earthquake ID : koeri2023ctri    

        Stat    Net   Time UTC                  Phase   Res     Distance   Azimuth   Status      Amp (MLv)      Mag (MLv)
        ------- ----- ------------------------- ------- ------- ---------- --------- ----------- -------------- -----------
        SARI    KO    2023-02-09 05:24:33.449   P       -0.2    0.4        2.1       automatic   69.56897307    4.17
        KMRS    KO    2023-02-09 05:24:38.389   P       0.3     0.6        5.0       automatic   50.37548677    4.54
        SARI    KO    2023-02-09 05:24:38.999   S       -0.2    0.4        2.1       automatic   -              -
        DARE    KO    2023-02-09 05:24:41.420   P       -0.5    0.8        0.6       automatic   11.34292752    4.00
        KHMN    KO    2023-02-09 05:24:42.129   P       0.5     0.8        5.2       automatic   42.34370569    4.57
        KOZT    KO    2023-02-09 05:24:43.039   P       -0.4    0.9        3.9       automatic   11.99514082    4.08
        KMRS    KO    2023-02-09 05:24:46.939   S       0.3     0.6        5.0       automatic   -              -
        CEYT    KO    2023-02-09 05:24:50.900   P       -0.1    1.3        4.1       automatic   4.763989123    3.90

         */


        if (!in_hypocenter_event) {
            // find phase line
            while (strncmp(line, "  -----", 7) != 0) {
                // read next line
                cstat = fgets(line, MAXLINE_LONG, fp_obs);
                //printf(">0 %s", line);
                if (cstat == NULL)
                    return (OBS_FILE_END_OF_INPUT);
            }
            in_hypocenter_event = 1;
        }

        /* read line */
        cstat = fgets(line, MAXLINE_LONG, fp_obs);
        //printf(">1 %s", line);
        if (cstat == NULL)
            return (OBS_FILE_END_OF_INPUT);
        /* check for end of event (assumes empty lines between event) */
        if (LineIsBlank(line)) {
            /* end of event */
            return (OBS_FILE_END_OF_EVENT);
        }

        /* read phase arrival input */
        /*
        SARI    KO    2023-02-09 05:24:33.449   P       -0.2    0.4        2.1       automatic   69.56897307    4.17
        KMRS    KO    2023-02-09 05:24:38.389   P       0.3     0.6        5.0       automatic   50.37548677    4.54
         */
        char koeri_sta[16], koeri_net[16];
        double koeri_mag;
        istat = sscanf(line, "%s %s %4d-%2d-%2d %2d:%2d:%lf %s %*f %*f %*f  %*s  %*f  %lf",
                koeri_sta, koeri_net, &arrival->year, &arrival->month, &arrival->day, &arrival->hour, &arrival->min, &arrival->sec, arrival->phase, &koeri_mag);
        //printf("DEBUG: KOERI_EQEVENTS: %s %s %4d-%2d-%2dT%2d:%2d:%f\n", arrival->label, arrival->phase, arrival->year, arrival->month, arrival->day, arrival->hour, arrival->min, arrival->sec);

        if (istat < 9) {
            return (OBS_FILE_END_OF_EVENT);
        }

        // check for valid magnitude
        if (istat > 9) {
            //istat = sscanf(line, "%*s %lf %*s %*f %d", &phypo->amp_mag, &phypo->num_amp_mag);
            if (phypo->num_amp_mag > 0) {
                phypo->amp_mag = (phypo->amp_mag * (double) phypo->num_amp_mag + koeri_mag);
                (phypo->num_amp_mag)++;
                phypo->amp_mag /= (double) phypo->num_amp_mag;
            } else {
                phypo->amp_mag = koeri_mag;
                phypo->num_amp_mag = 1;
            }
        }

        snprintf(arrival->label, sizeof(arrival->label), "%s_%s", koeri_net, koeri_sta);

        // set error */
        snprintf(arrival->error_type, sizeof(arrival->error_type), "%s", "GAU");
        // error if not P
        arrival->error = Quality2Error[1];
        // error if first arrival P
        if (strstr("P$Pg$Pn$Pb$P0$P1$PKP$PKPdf", arrival->phase) != NULL) {
            arrival->error = Quality2Error[0];
        }

        return (istat);


    } else {
        nll_puterr2("ERROR: unrecognized observation file type", ftype_obs);

        return (OBS_FILE_END_OF_INPUT);
    }
}
