/*
 * Copyright (C) 1999-2023 Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
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


/*   loc_combine.c

        Program to merge into a set of primary NLL locations the arrivals form a second set of NLL locations.

 */



/*
        history:

        ver 01    28Apr2023  AJL  Original version



.........1.........2.........3.........4.........5.........6.........7.........8

 */



#include "GridLib.h"
#include "phaseloclist.h"


/* defines */

#define MAX_NUM_INPUT_FILES 50000

// 20221020 AJL - replaced X_MAX_NUM_ARRIVALS with this defined constant so larger numbers of phases can be processed
#define MAX_NUM_ARRIVALS_LOC2SSST 20000

/* globals */

// station list
int NumStationPhases;
SourceDesc StationPhaseList[MAX_NUM_ARRIVALS_LOC2SSST];
// station/phase ssst grids
GridDesc ssst_grid_template;
GridDesc ssst_time_grid_template;
// LocNode array
int NumLocNodes;
LocNode *LocNodeArray[MAX_NUM_INPUT_FILES];
// PhsNode array
int NumPhsNodes;
PhsNode *PhsNodeArray[MAX_NUM_ARRIVALS_LOC2SSST];

// station list input line can be very long
#define MAX_LEN_STATION_LINE 64000

typedef struct {
    double max_distance_horiz; // Maximum epicentral difference between hypocenters to combine (km)
    double max_distance_vert; // Maximum depth difference between hypocenters to combine (km)
    double max_otime_diff; // Maximum otime difference between hypocenters to combine (sec)
    int use_rejected; // Flag to indicate that NLL REJECTED locations should be accepted for processing (default=0)
}
LC_Params;
LC_Params Params;

typedef struct {
    double RMSMax;
    int NRdgsMin;
    double GapMax, PResidualMax, SResidualMax, EllLen3Max;
    double DepthMin, DepthMax; // 20221116 AJL - added to support filtering of events (phases) by hypocenter depth
}
LC_Phstat;
LC_Phstat PhsStat;

char secondary_phase_filter[MAXLINE_LONG];
int have_secondary_phase_filter = 0;

int angle_mode = ANGLE_MODE_NO; /* angle mode - ANGLE_MODE_NO, ANGLE_MODE_YES */

char fn_hypos_primary_in[FILENAME_MAX];
char fn_hypos_secondary_in[FILENAME_MAX];
char fn_hypos_combined_out[FILENAME_MAX];
char fn_time_input[FILENAME_MAX] = "";
double VpVsRatio;
int iSwapBytesOnInput;


/* functions */

int Readloc_combineInput(FILE*);
int Doloc_combine();
int read_input_locations(LocNode **loc_list_head, char* fn_hypos_primary_in);



/*** program to combine readings in two sets of hyp files */

#define PNAME  "loc_combine"
#define NARGS 2

int main(int argc, char** argv) {

    int istat;

    /* set program name */

    // 20210511 AJL - Bug fix: changed to strncpy from strlcpy which is not available in Linux.
    strncpy(prog_name, PNAME, sizeof (prog_name));


    // check command line for correct usage

    if (argc != NARGS) {
        disp_usage(prog_name, "merge into a set of primary NLL locations the arrivals form a second set of NLL locations");
        disp_usage(prog_name, "<control file>");
        exit(EXIT_ERROR_USAGE);
    }

    SetConstants();
    PhaseFormat = FORMAT_PHASE_2;

    // initializations before reading input
    Params.max_distance_horiz = LARGE_FLOAT;
    Params.max_distance_vert = LARGE_FLOAT;
    Params.max_otime_diff = LARGE_FLOAT;
    Params.use_rejected = 0;

    PhsStat.RMSMax = 1.0e6;
    PhsStat.NRdgsMin = 0;
    PhsStat.GapMax = 360.0;
    PhsStat.PResidualMax = 1.0e6;
    PhsStat.SResidualMax = 1.0e6;
    PhsStat.EllLen3Max = 1.0e6;
    PhsStat.DepthMin = -1.0e6;
    PhsStat.DepthMax = 1.0e6;

    strcpy(secondary_phase_filter, "\0");

    // read control file
    strncpy(fn_control, argv[1], sizeof (fn_control) - 1);
    if ((fp_control = fopen(fn_control, "r")) == NULL) {
        nll_puterr("ERROR: opening control file.");
        exit(EXIT_ERROR_FILEIO);
    }

    if ((istat = Readloc_combineInput(fp_control)) < 0) {
        nll_puterr("ERROR: reading loc_combine control file.");
        exit(EXIT_ERROR_FILEIO);
    }

    if ((istat = Doloc_combine()) < 0) {
        nll_puterr("ERROR: doing loc_combine process.");
        exit(-1);
    }



    exit(0);

}

/*** function to read input file name ***/

int get_lc_inpfile(char* line1) {

    sscanf(line1, "%s %s", fn_hypos_primary_in, fn_hypos_secondary_in);

    snprintf(MsgStr, sizeof (MsgStr), "LCLOCFILES:  Input: Primary: %s  Secondary: %s", fn_hypos_primary_in, fn_hypos_secondary_in);
    nll_putmsg(1, MsgStr);

    return (0);
}

/*** function to read output file name ***/

int get_lc_outfile(char* line1) {

    sscanf(line1, "%s", fn_hypos_combined_out);

    snprintf(MsgStr, sizeof (MsgStr), "LCOUT:  Output: %s.*",
            fn_hypos_combined_out);
    nll_putmsg(1, MsgStr);

    return (0);
}

/*** function to read secondary hypos file phase filter ***/

int get_lc_phfilter(char* line1) {

    sscanf(line1, "%s", secondary_phase_filter);

    snprintf(MsgStr, sizeof (MsgStr), "LCPHFILTER: %s", secondary_phase_filter);
    nll_putmsg(1, MsgStr);

    return (0);
}

/*** function to read hypocenter filters ***/

int get_lc_phstat(char* line1) {

    int istat = sscanf(line1, "%lf %d %lf %lf %lf %lf %lf %lf",
            &PhsStat.RMSMax, &PhsStat.NRdgsMin, &PhsStat.GapMax, &PhsStat.PResidualMax, &PhsStat.SResidualMax, &PhsStat.EllLen3Max, &PhsStat.DepthMin, &PhsStat.DepthMax);

    sprintf(MsgStr,
            "LCPHSTAT:  RMSMax: %f  NRdgsMin: %d  GapMax: %.3g  PResidualMax: %.3g SResidualMax: %.3g EllLen3Max %.3g DepthMin %.3g DepthMax %.3g",
            PhsStat.RMSMax, PhsStat.NRdgsMin, PhsStat.GapMax, PhsStat.PResidualMax, PhsStat.SResidualMax, PhsStat.EllLen3Max, PhsStat.DepthMin, PhsStat.DepthMax);
    nll_putmsg(1, MsgStr);

    if (istat < 6)
        return (-1);

    return (0);
}

/*** function to read ls parameters ***/

int get_lc_params(char* line1) {

    int istat = sscanf(line1, "%lf %lf %lf %d",
            &Params.max_distance_horiz, &Params.max_distance_vert, &Params.max_otime_diff, &Params.use_rejected);

    sprintf(MsgStr,
            "LCPARAMS:  max_distance_horiz: %f  max_distance_vert:%f  max_otime_diff:%f  use_rejected:%d",
            Params.max_distance_horiz, Params.max_distance_vert, Params.max_otime_diff, Params.use_rejected);
    nll_putmsg(1, MsgStr);

    if (istat < 1)

        return (-1);

    return (0);
}

/*** function to read input file */

int Readloc_combineInput(FILE * fp_input) {

    int istat, iscan;
    char param[MAXLINE], *pchr;
    char line[MAX_LEN_STATION_LINE], *fgets_return;

    int flag_control = 0, flag_inpfile = 0, flag_outfile = 0, flag_params = 0;
    int flag_include = 1;
    int flag_lcfilter = 0;
    int flag_phstat = 1;

    // read each input line

    while ((fgets_return = fgets(line, MAX_LEN_STATION_LINE, fp_input)) != NULL
            || fp_include != NULL) {


        // check for end of include file

        if (fgets_return == NULL && fp_include != NULL) {
            SwapBackIncludeFP(&fp_input);
            continue;
        }


        istat = -1;

        // read parameter line

        if ((iscan = sscanf(line, "%s", param)) < 0)
            continue;

        // skip comment line or white space

        if (strncmp(param, "#", 1) == 0 || isspace(param[0]))
            istat = 0;


        // read include file params and set input to include file

        if (strcmp(param, "INCLUDE") == 0) {
            if ((istat = GetIncludeFile(strchr(line, ' '),
                    &fp_input)) < 0) {
                nll_puterr("ERROR: processing INCLUDE file.");
                flag_include = 0;
            }
        }


        // read control params

        if (strcmp(param, "CONTROL") == 0) {
            if ((istat = get_control(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading CONTROL parameters.");
            else
                flag_control = 1;
        }


        // read output file name

        if (strcmp(param, "LCLOCFILES") == 0) {
            if ((istat = get_lc_inpfile(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading LCLOCFILES parameters.");
            else
                flag_inpfile = 1;
        }


        // read output file name

        if (strcmp(param, "LCOUT") == 0) {
            if ((istat = get_lc_outfile(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading LCOUT parameters.");
            else
                flag_outfile = 1;
        }


        // read secondary phase filter

        if (strcmp(param, "LCPHFILTER") == 0) {
            if ((istat = get_lc_phfilter(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading LCPHFILTER parameters.");
            else {
                flag_lcfilter = 1;
                have_secondary_phase_filter = 1;
            }
        }


        // read parameters

        if (strcmp(param, "LCPARAMS") == 0) {
            if ((istat = get_lc_params(strchr(line, ' '))) < 0) {
                nll_puterr("ERROR: reading LCPARAMS parameters.");
            } else {
                flag_params = 1;
            }
        }


        // read filters

        if (strcmp(param, "LCPHSTAT") == 0) {
            if ((istat = get_lc_phstat(strchr(line, ' '))) < 0) {
                nll_puterr("ERROR: reading LCPHSTAT parameters.");
                flag_phstat = 0;
            }
        }



        // unrecognized input

        if (istat < 0) {
            if ((pchr = strchr(line, '\n')) != NULL)
                *pchr = '\0';
            snprintf(MsgStr, sizeof (MsgStr), "Skipping input: %s", line);
            nll_putmsg(4, MsgStr);
        }

    }

    // check for missing input

    if (!flag_control)
        nll_puterr("ERROR: no control (CONTROL) params read.");
    if (!flag_inpfile)
        nll_puterr("ERROR: no inputfile (LCLOCFILES) params read.");
    if (!flag_outfile)
        nll_puterr("ERROR: no outputfile (LCOUT) params read.");
    if (!flag_params)
        nll_puterr("ERROR: no parameters (LCPARAMS) read.");
    if (!flag_lcfilter)
        nll_putmsg(2, "INFO: no parameters (LCPHFILTER) read.");

    return (flag_include * flag_control * flag_inpfile * flag_outfile * flag_params * flag_phstat - 1);
}

/** function to calculate distance between 2 XYZ points */

double Dist3D2_loc_combine(double x1, double x2, double y1, double y2, double z1, double z2, int latlon) {

    if (latlon) {
        double hdist = GCDistance(y1, x1, y2, x2);
        double dz = z1 - z2;
        return (hdist * hdist + dz * dz); // approximate for short distances on sphere!
    } else {
        double dx, dy, dz;
        dx = x1 - x2;
        dy = y1 - y2;
        dz = z1 - z2;
        return (dx * dx + dy * dy + dz * dz);
    }
}

/** function to calculate ssst corrections */

int Doloc_combine() {


    // read input locations
    LocNode *loc_list_primary_head = NULL; // root node of primary locations list
    //int n_hypos_primary_read = read_input_locations(&loc_list_primary_head, fn_hypos_primary_in);
    LocNode *loc_list_secondary_head = NULL; // root node of secondary locations list
    //int n_hypos_secondary_read = read_input_locations(&loc_list_secondary_head, fn_hypos_secondary_in);

    // set some variables
    char fn_root_out[FILENAME_MAX], fname[2 * FILENAME_MAX], fn_root_out_last[FILENAME_MAX];
    strcpy(fn_root_out_last, "");
    int n_file_root_count = 1;
    LocNode *loc_list_secondary_curr = loc_list_secondary_head; // current node of secondary locations list

    // loop over primary locations
    int nmatch = 0;
    LocNode* locNode_primary = NULL;
    int loc_id_primary = 0;
    //printf("TP:00\n");
    while ((locNode_primary = getLocationFromLocList(loc_list_primary_head, loc_id_primary)) != NULL) {
        //printf("TP:00\n");
        HypoDesc *phypo_primary = locNode_primary->plocation->phypo;
        // loop over secondary locations
        LocNode* locNode_secondary = NULL;
        int loc_id_secondary = 0;
        int ifound = 0;
        while ((locNode_secondary = getLocationFromLocList(loc_list_secondary_curr, loc_id_secondary)) != NULL) {
            //printf("TP:01\n");
            HypoDesc *phypo_secondary = locNode_secondary->plocation->phypo;
            // check otime match
            double tdiff = getLocTimeValue(phypo_primary) - getLocTimeValue(phypo_secondary);
            if (fabs(tdiff) < Params.max_otime_diff) {
                //printf("TP:02\n");
                // check depth match
                double zdiff = phypo_primary->depth - phypo_secondary->depth;
                if (fabs(zdiff) < Params.max_distance_vert) {
                    //printf("TP:03\n");
                    // check epicenter match
                    double hdiff = GCDistance(
                            phypo_primary->dlat, phypo_primary->dlong,
                            phypo_secondary->dlat, phypo_secondary->dlong);
                    if (fabs(hdiff) < Params.max_distance_horiz) {
                        //printf("TP:04\n");
                        // save locations with merged arrivals
                        printf("Match:\n");
                        printf("   %4.4d %2.2d %2.2d  %2.2d %2.2d %lf  Lat %lf Long %lf Depth %lf  N=%d\n",
                                phypo_primary->year, phypo_primary->month, phypo_primary->day,
                                phypo_primary->hour, phypo_primary->min, (double) phypo_primary->sec,
                                phypo_primary->dlat, phypo_primary->dlong, phypo_primary->depth, loc_id_primary);
                        printf("   %4.4d %2.2d %2.2d  %2.2d %2.2d %lf  Lat %lf Long %lf Depth %lf  N=%d\n",
                                phypo_secondary->year, phypo_secondary->month, phypo_secondary->day,
                                phypo_secondary->hour, phypo_secondary->min, (double) phypo_secondary->sec,
                                phypo_secondary->dlat, phypo_secondary->dlong, phypo_secondary->depth, loc_id_secondary);
                        /*// merge arrivals
                        int NumArrivals = 0;
                        for (int narr = 0; narr < locNode_primary->plocation->narrivals; narr++) {
                            Arrival[NumArrivals] = *(locNode_primary->plocation->parrivals + narr);
                            NumArrivals++;
                        }
                        for (int narr = 0; narr < locNode_secondary->plocation->narrivals; narr++) {
                            Arrival[NumArrivals] = *(locNode_secondary->plocation->parrivals + narr);
                            NumArrivals++;
                        }
                        // write primary location and merged arrivals to hyp file
                        int iSaveDecSec = 0;
                        int iSavePublicID = 0;
                        SetOutName(Arrival + 0, fn_hypos_combined_out, fn_root_out, fn_root_out_last, iSaveDecSec,
                                iSavePublicID, phypo_primary->public_id, &n_file_root_count);
                        sprintf(fname, "%s.grid0.loc.hyp", fn_root_out);
                        int iWriteArrivals = 1;
                        int iWriteEndLoc = 1;
                        int iWriteMinimal = 1;
                        WriteLocation(NULL, phypo_primary, Arrival, NumArrivals, fname,
                                iWriteArrivals, iWriteEndLoc, iWriteMinimal,
                                locNode_primary->plocation->pgrid, -1);*/
                        nmatch++;
                        // break
                        ifound = 1;
                    }

                }
            }
            if (ifound) {
                break;
            }
            loc_id_secondary++;
        }
        // include primary location arrivals to hyp file
        int NumArrivals = 0;
        for (int narr = 0; narr < locNode_primary->plocation->narrivals; narr++) {
            Arrival[NumArrivals] = *(locNode_primary->plocation->parrivals + narr);
            NumArrivals++;
        }
        if (ifound) {
            // match, include secondary location arrivals in hyp file
            for (int narr = 0; narr < locNode_secondary->plocation->narrivals; narr++) {
                // check secondary phase filter
                if (have_secondary_phase_filter) {
                    if (strstr(secondary_phase_filter, (locNode_secondary->plocation->parrivals + narr)->phase) != NULL) {
                        continue;
                    }
                }
                Arrival[NumArrivals] = *(locNode_secondary->plocation->parrivals + narr);
                NumArrivals++;
            }
            // continue search from current secondary hyp (assumes primary and secondary hypos sorted in time!)
            loc_list_secondary_curr = locNode_secondary;
        }
        int iSaveDecSec = 0;
        int iSavePublicID = 0;
        SetOutName(Arrival + 0, fn_hypos_combined_out, fn_root_out, fn_root_out_last, iSaveDecSec,
                iSavePublicID, phypo_primary->public_id, &n_file_root_count);
        sprintf(fname, "%s.grid0.loc.hyp", fn_root_out);
        int iWriteArrivals = 1;
        int iWriteEndLoc = 1;
        int iWriteMinimal = 1;
        WriteLocation(NULL, phypo_primary, Arrival, NumArrivals, fname,
                iWriteArrivals, iWriteEndLoc, iWriteMinimal,
                locNode_primary->plocation->pgrid, -1);

        loc_id_primary++;
    }

    printf("Num Primary: %d, Num Matched: %d\n", loc_id_primary, nmatch);
    printf("Output path: %s\n", fn_hypos_combined_out);


    return (0);

}

/** function to read input NLL Location files to location list */

int read_input_locations(LocNode **loc_list_head, char* fn_hypos_in) {

    int istat;

    double Len3Mean = 0.0, RMSMean = 0.0;
    double NRdgsMean = 0.0, GapMean = 0.0;
    double DepthMean = 0.0;
    int Len3Reject = 0, RMSReject = 0;
    int NRdgsReject = 0, GapReject = 0;
    int DepthMinReject = 0, DepthMaxReject = 0;
    int AbortedReject = 0, RejectedReject = 0;
    int iReject;

    //char *pchr;
    //char sys_string[FILENAME_MAX];
    //char filename[FILENAME_MAX];
    FILE *fp_hypo;

    int nFile, numFiles, nLocRead, nLocAccepted;
    //char fn_hyp_in_list[MAX_NUM_INPUT_FILES][FILENAME_MAX];
    // make fn_hyp_in_list static to avoid stack overflow problem (F Tilmann) //
    static char fn_hyp_in_list[MAX_NUM_INPUT_FILES][FILENAME_MAX];
    static HypoDesc hypo;
    int nHypos;

    // check for wildcards in input file name
    if ((numFiles = ExpandWildCards(fn_hypos_in, fn_hyp_in_list, MAX_NUM_INPUT_FILES)) < 1) {
        nll_puterr("ERROR: no matching .hyp files found.");
        return (-1);
    }
    if (numFiles >= MAX_NUM_INPUT_FILES) {
        snprintf(MsgStr, sizeof (MsgStr), "WARNING: maximum number of event files exceeded, only first %d will be processed.", MAX_NUM_INPUT_FILES);
        nll_puterr(MsgStr);
    }
    nLocRead = 0;
    nLocAccepted = 0;
    for (nFile = 0; nFile < numFiles; nFile++) {

        //fprintf(OUT_LEVEL_1, "Reading location <%s>                       \r", fn_hyp_in_list[nFile]);

        // open hypocenter file
        if ((fp_hypo = fopen(fn_hyp_in_list[nFile], "r")) == NULL) {
            nll_puterr2("ERROR: opening hypocenter file, ignoring event, file",
                    fn_hyp_in_list[nFile]);
            continue;
        }

        // loop over events in file

        while (1) {

            istat = GetHypLoc(fp_hypo, NULL, &hypo, Arrival, &NumArrivals, 1, NULL, 0);
            if (istat == EOF) {
                break;
            }
            if (istat < 0) {
                nll_puterr2("ERROR: reading hypocenter file, ignoring event, file",
                        fn_hyp_in_list[nFile]);
                break;
            }

            if (strcmp(hypo.locStat, "ABORTED") == 0) {
                //nll_puterr("WARNING: location ABORTED, ignoring event");
                AbortedReject++;
                continue;
            } else if (!Params.use_rejected && strcmp(hypo.locStat, "REJECTED") == 0) {
                //nll_puterr("WARNING: location REJECTED, ignoring event");
                RejectedReject++;
                continue;
            }


            Len3Mean += hypo.ellipsoid.len3;
            RMSMean += hypo.rms;
            NRdgsMean += (double) hypo.nreadings;
            GapMean += hypo.gap;
            DepthMean += hypo.z;
            nLocRead++;

            iReject = 0;
            if (hypo.rms > PhsStat.RMSMax) {
                //nll_puterr("WARNING: location RMS is Greater than RMSMax, ignoring event");
                RMSReject++;
                iReject = 1;
            }
            if (hypo.nreadings < PhsStat.NRdgsMin) {
                //nll_puterr("WARNING: location num readings is less than NRdgsMin, ignoring event");
                NRdgsReject++;
                iReject = 1;
            }
            if (hypo.gap > PhsStat.GapMax) {
                //nll_puterr("WARNING: location gap is greater than GapMax, ignoring event");
                GapReject++;
                iReject = 1;
            }
            if (hypo.ellipsoid.len3 > PhsStat.EllLen3Max) {
                //nll_puterr("WARNING: ellipsoid.len3 is greater than EllLen3Max, ignoring event");
                Len3Reject++;
                iReject = 1;
            }
            if (hypo.z < PhsStat.DepthMin) {
                nll_puterr("WARNING: hypo.z is less than DepthMin, ignoring event");
                DepthMinReject++;
                iReject = 1;
            }
            if (hypo.z > PhsStat.DepthMax) {
                nll_puterr("WARNING: hypo.z is greater than DepthMax, ignoring event");
                DepthMaxReject++;
                iReject = 1;
            }
            if (iReject)
                continue;

            // initialize some hypo fields
            latlon2rect(0, hypo.dlat, hypo.dlong, &(hypo.x), &(hypo.y));
            hypo.z = hypo.depth;

            Location* ploc_list_node = newLocation(
                    cloneHypoDesc(&hypo),
                    cloneArrivalDescArray(Arrival, NumArrivals),
                    NumArrivals, NULL,
                    NULL,
                    NULL);
            *loc_list_head = addLocationToLocList(loc_list_head, ploc_list_node, nLocAccepted);

            nLocAccepted++;


        }
        fclose(fp_hypo);

    }
    fprintf(OUT_LEVEL_1, "\n");

    if (nLocRead > 0) {
        // write messages
        fprintf(stdout,
                "%d location files read, %d accepted.\n",
                numFiles, nLocAccepted);
        fprintf(stdout,
                "Loc Aborted: Reject %d\n", AbortedReject);
        fprintf(stdout,
                "Loc Rejected: Reject %d\n", RejectedReject);
        fprintf(stdout,
                "RMS Mean: %lf, Reject %d\n",
                RMSMean / (double) nLocRead, RMSReject);
        fprintf(stdout,
                "NRdgs Mean: %lf, Reject %d\n",
                NRdgsMean / (double) nLocRead, NRdgsReject);
        fprintf(stdout,
                "Gap Mean: %lf, Reject %d\n",
                GapMean / (double) nLocRead, GapReject);
        fprintf(stdout,
                "Len3 Mean: %lf, Reject %d\n",
                Len3Mean / (double) nLocRead, Len3Reject);
        fprintf(stdout,
                "Depth Mean: %lf, DepthMinReject %d\n",
                DepthMean / (double) nLocRead, DepthMinReject);
        fprintf(stdout,
                "Depth Mean: %lf, DepthMaxReject %d\n",
                DepthMean / (double) nLocRead, DepthMaxReject);
    }

    nHypos = nLocAccepted;


    return (nHypos);

}
