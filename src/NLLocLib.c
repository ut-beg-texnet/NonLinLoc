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


/*   NLLoc.c

        Program to do global search earthquake location in 3-D models

 */


/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
Mouans-Sartoux, France
e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


/*
        history:	(see also http://alomax.net/nlloc -> Updates)

        ver 01    26SEP1997  AJL  Original version
        ver 02    08JUN1998  AJL  Metropolis added
        ver  2         2000  AJL  Oct-Tree added
        ver  3      DEC2003  AJL  EDT added
        ver  4    10MAY2004  AJL  Added following changes from S. Husen:
                27MAR2002  *SH   VELEST phase format added (changes in GetNextObservation)
                28AUG2002  *SH   event origin time is now calculated relative to
                                the second (and not minute as before); initial
                                OT seconds are now read in and added to arrival
                                time (changes in GetNextObservation)
                17NOV2002  *SH   UUSS phase format added (changes in GetObservations
                                and GetNextObservation)
                12JUN2003  *SH   Fixed bug with arrival times > 100s in UUSS phase format
                01OCT2003  *SH   added SED format "SED_LOC"
                                code was written by A. Lomax; bug fix by Danijel Schorlemmer
                02MAR2004  *SH   modifications to make NLLoc compatible for routine earthquake
                                location with SNAP (SED):
                                - introduced 2nd argument snap_pid, which is the snap_pid of SNAP; needed
                                to form filename of outputfile hyprint{snap_pid}; usage of NLLoc
                                now is:
                                NLLoc <control file> <snap_pid>
                                - added new subroutine WriteSnapSum: output location results
                                into file hyprint{snap_pid} in format readable by SNAP; format
                                is identical to output format of program grid_search by
                                M. Baer of SED;
                                file hyprint{snap_pid} will be written if control parameter
                                LOCHYPOUT is set to SAVE_SNAP_SUM
                                - added subroutine get_region_names_nr and associated subroutines
                                to convert lat/lon into Swiss coordinates in km and to find
                                region name for local earthquakes in Switzerland
                NOV2004  AJL   Split off NLLocLib and created NLDiffLoc (non-linear double-difference location)
                APR2005  AJL   Added LOCSTAWT, LOCELEVCORR, LOCDELAY_SURFACE
                JAN2006 Frederik Tilmann    Changes to SEISAN reader
                        - 5 character station names now permitted
                        - checks format line code in column 80 to only read phase lines.
                                Previously some 'unlucky' lines
                                of other types got interpreted as phase lines
                        - can now read high precision phase times (accurate to 0.001 s) .
                                Previously these high precision
                                picks drops 10s of seconds, resulting in completely wrong times
                        - now read instrument and component
                APR2006  AJL   Added  LOCMETH-EDT_OT_WT, LOCSEARCH-OCT:useStationsDensity, stopOnMinNodeSize
                ...
                20100506 AJL - added to support preservation of observation index order for calls to NLLoc() function (e.g. from SeisComp3)
                20130627 AJL - add prior pick weighting, change station distribution weighing from sum to product
                201501   AJL - added EW_PTWC_HAWAII obs format
                20150324 AJL - added L1_NORM
                20170811 AJL - added HYPO_TYPE_EXPECTATION: support for expectation hypocenter results output
                20261103 AJL - applied important fixes from ML analysis of code by Mustafa Çomoğlu, see NLL_fixes/NLL_fixes_for_src_20261102




.........1.........2.........3.........4.........5.........6.........7.........8

 */


/* References */
/*
        TV82	Tarantola and Valette,  (1982)
                "Inverse Problems = Quest for Information",
                J Geophys 50, 159-170.
        MEN92	Moser, van Eck and Nolet,  (1992)
                "Hypocenter Determination ... Shortest Path Method",
                JGR 97, B5, 6563-6572.
 */



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


// define globals

char f_outpath[FILENAME_MAX];
GaussLocParams Gauss;
Gauss2LocParams Gauss2;
int iUseGauss2;
ScatterParams Scatter;
int NumEvents;
int NumEventsLocated;
int NumLocationsCompleted;
int NumObsFiles;
int NumArrivalsRead;
int NumArrivalsLocation;
char fn_loc_obs[MAX_NUM_OBS_FILES][FILENAME_MAX];
char ftype_obs[MAXLINE];

// 20251027 add support for alternative travel-time grid path/root
char fn_time_grids[MAX_NUM_TIME_GRID_PATHS][FILENAME_MAX];
int NumTimeGridPaths;

char fn_path_output[FILENAME_MAX];
int iSwapBytesOnInput;
FILE *fp_model_grid_P;
FILE *fp_model_hdr_P;
GridDesc model_grid_P;
FILE *fp_model_grid_S;
FILE *fp_model_hdr_S;
GridDesc model_grid_S;
int SearchType;
SearchPdfGridDesc SearchPrior;
int iUseSearchPrior;
SearchPdfGridDesc SearchPosterior;
int iUseSearchPosterior;
int LocMethod;
int EDT_use_otime_weight;
int EDT_otime_weight_active;
double DistStaGridMin;
double DistStaGridMax;
int MinNumArrLoc;
int MaxNumArrLoc;
int MinNumSArrLoc;
double VpVsRatio;
char LocSignature[MAXLINE_LONG];
GridDesc LocGrid[MAX_NUM_LOCATION_GRIDS];
int NumLocGrids;
int LocGridSave[MAX_NUM_LOCATION_GRIDS]; /* !should be in GridDesc */
//int Num3DGridReadToMemory, MaxNum3DGridMemory;
int iWriteHypHeader[MAX_NUM_LOCATION_GRIDS];
char HypoInverseArchiveSumHdr[MAXLINE_LONG];
int iSaveNLLocEvent, iSaveNLLocSum, iSaveNLLocSumCSV, iSaveNLLocOctree,
iSaveHypo71Event, iSaveHypo71Sum,
iSaveHypoEllEvent, iSaveHypoEllSum,
iSaveHypoInvSum, iSaveHypoInvY2KArc, iSaveAlberto4Sum, iSaveFmamp,
iSaveSnapSum, iCalcSedOrigin, iSaveDecSec, iSavePublicID, iSaveNone;
int iSaveNLLocExpectation;
int iSaveNLLocEvent_JSON;
int iUseArrivalPriorWeights;
int iSetStationDistributionWeights;
double stationDistributionWeightCutoff;
double AveInterStationDistance;
int NumForceOctTreeStaDenWt;
int iRejectDuplicateArrivals;
EventTimeExtract EventTime;
long int EventID;
int NumMagnitudeMethods;
MagDesc Magnitude[MAX_NUM_MAG_METHODS];
CompDesc Component[MAX_NUM_COMP_DESC];
int NumCompDesc;
AliasDesc LocAlias[MAX_NUM_LOC_ALIAS];
int NumLocAlias;
ExcludeDesc LocExclude[MAX_NUM_LOC_EXCLUDE];
int NumLocExclude;
ExcludeDesc LocInclude[MAX_NUM_LOC_INCLUDE];
int NumLocInclude;
TimeDelayDesc TimeDelay[MAX_NUM_STA_DELAYS];
int NumTimeDelays;
char TimeDelaySurfacePhase[MAX_SURFACES][PHASE_LABEL_LEN];
double TimeDelaySurfaceMultiplier[MAX_SURFACES];
int NumTimeDelaySurface;
int ApplyElevCorrFlag;
double ElevCorrVelP;
double ElevCorrVelS;
int ApplyCrustElevCorrFlag;
double MinDistCrustElevCorr;
struct surface *topo_surface;
int topo_surface_index; // topo surface index is velmod.h.MAX_SURFACES-1 so as not to interferce with any TimeDelaySurfaces read in
int NumStationPhases;
SourceDesc StationPhaseList[X_MAX_NUM_ARRIVALS];
int FixOriginTimeFlag;
WalkParams Metrop; /* walk parameters */
int MetNumSamples; /* number of samples to evaluate */
int MetLearn; /* learning length in number of samples for calculation of sample statistics */
int MetEquil; /* number of samples to equil before using */
int MetStartSave; /* number of sample to begin saving */
int MetSkip; /* number of samples to wait between saves */
double MetStepInit; /* initial step size (km) (< 0.0 for auto) */
double MetStepMin; /* minimum step size (km) */
double MetStepMax; /* maximum step size (km) (NLDiffLoc) */
double MetStepFact; /* step size factor */
double MetProbMin; /* minimum likelihood necessary after learn */
double MetVelocity; /* velocity for conversion of distance to time */
double MetInititalTemperature; /* initial temperature */
int MetUse; /* number of samples to use = MetNumSamples - MetEquil */
OcttreeParams octtreeParams; /* Octtree parameters */
Tree3D* octTree; /* Octtree */
ResultTreeNode* resultTreeRoot; /* Octtree likelihood*volume results tree root node */
//ResultTreeNode* resultTreeLikelihoodRoot;	/* Octtree likelihood results tree root node */
int angleMode; /* angle mode - ANGLE_MODE_NO, ANGLE_MODE_YES */
int iAngleQualityMin; /* minimum quality for angles to be used */
OtimeLimit** OtimeLimitList;
int NumOtimeLimit;
int NRdgs_Min;
double RMS_Max, Gap_Max;
double P_ResidualMax;
double S_ResidualMax;
double Ell_Len3_Max;
double Hypo_Depth_Min;
double Hypo_Depth_Max;
double Hypo_Dist_Max;
//char snap_pid[255];

FILE *pSumFileHypNLLoc[MAX_NUM_LOCATION_GRIDS];
FILE *pSumFileHypNLLocCSV[MAX_NUM_LOCATION_GRIDS]; // 20240826 AJL - added
FILE *pSumFileHypo71[MAX_NUM_LOCATION_GRIDS];
FILE *pSumFileHypoEll[MAX_NUM_LOCATION_GRIDS];
FILE *pSumFileHypoInv[MAX_NUM_LOCATION_GRIDS];
FILE *pSumFileHypoInvY2K[MAX_NUM_LOCATION_GRIDS];
FILE *pSumFileAlberto4[MAX_NUM_LOCATION_GRIDS];
FILE *pSumFileFmamp[MAX_NUM_LOCATION_GRIDS];


// AJL - 20080710 (valgrind)
/* locally allocated memory which must be cleaned up */

int clean_memory(int istat);

// EDT_OT_WT_ML allocations
double *ot_ml_arrival = NULL; // array of ot estimate for each arrival
double *ot_ml_arrival_edt_sum = NULL; // array of weight of ot estimate for each arrival
int isize_ot_ml_array = 0;

// ConstWeightMatrix() allocations
MatrixDouble wt_matrix = NULL;
MatrixDouble edt_matrix = NULL;
int last_matrix_alloc_size = -1;

/** function to perform grid search location */

int Locate(int ngrid, char* fn_loc_obs, char* fn_root_out, int numArrivalsReject, int return_locations, int return_oct_tree_grid, int return_scatter_sample, LocNode **ploc_list_head) {

    int istat, n, narr;
    char fnout[4 * MAXLINE];

    FILE *fpio;
    char fname[5 * MAXLINE]; // should be bigger than fnout, to avoid compiler warnings
    float *fdata = NULL;
    float ftemp;
    int iSizeOfFdata;
    double oct_node_value_max, oct_tree_integral = 0.0;
    double oct_tree_prob_integral = 0.0;

    // AJL 20071219
    Location *ploc_list_node;


    Hypocenter.nScatterSaved = -1;



    /* write message */

    nll_putmsg(2, "");
    if (SearchType == SEARCH_GRID)
        snprintf(MsgStr, sizeof(MsgStr), "Searching Grid %d:", ngrid);
    else if (SearchType == SEARCH_MET)
        snprintf(MsgStr, sizeof(MsgStr), "Applying Metropolis within Grid %d:", ngrid);
    else if (SearchType == SEARCH_OCTTREE)
        snprintf(MsgStr, sizeof(MsgStr), "Applying Octtree search within Grid %d:", ngrid);
    nll_putmsg(2, MsgStr);
    if (message_flag >= 3)
        display_grid_param(LocGrid + ngrid);



    /* set output name */
    snprintf(fnout, sizeof(fnout), "%s.grid%d", fn_root_out, ngrid);
    snprintf(Hypocenter.fileroot, sizeof(Hypocenter.fileroot), "%s", fnout);

    /* initialize hypocenter fields */
    snprintf(Hypocenter.locStat, sizeof(Hypocenter.locStat), "LOCATED");
    snprintf(Hypocenter.locStatComm, sizeof(Hypocenter.locStatComm), "Location completed.");
    Hypocenter.x = Hypocenter.y = Hypocenter.z = 0.0;
    Hypocenter.ix = Hypocenter.iy = Hypocenter.iz = -1;
    // 20110620 AJL - preserve event id if available
    if (NumArrivalsLocation > 0 && Arrival[0].dd_event_id_1 >= 0)
        Hypocenter.event_id = Arrival[0].dd_event_id_1;
    else
        Hypocenter.event_id = -1;
    // 20170811 AJL - support for expectation hypocenter results output
    if (iSaveNLLocExpectation) {
        snprintf(Hypocenter.type, sizeof(Hypocenter.type), "%s", HYPO_TYPE_EXPECTATION);
    } else {
        snprintf(Hypocenter.type, sizeof(Hypocenter.type), "%s", HYPO_TYPE_MAXIMUM_LIKELIHOOD);
    }


    /* search type dependent initializations */

    if (SearchType == SEARCH_GRID) {

        /* check that current grid is contained within first grid */

        if (!IsGridInside(LocGrid + ngrid, LocGrid, 0)) {
            nll_puterr(
                    "WARNING: this grid not entirely contained inside 0th grid, ending search for this event.");
            return (clean_memory(GRID_NOT_INSIDE));
        }

        /* initialize 3D location grid */

        /* allocate location grid */
        LocGrid[ngrid].buffer = AllocateGrid(LocGrid + ngrid);
        if (LocGrid[ngrid].buffer == NULL) {
            nll_puterr(
                    "ERROR: allocating memory for 3D location grid buffer.");
            return (clean_memory(EXIT_ERROR_MEMORY));
        }
        /* create array access pointers */
        LocGrid[ngrid].array = CreateGridArray(LocGrid + ngrid);
        if (LocGrid[ngrid].array == NULL) {
            nll_puterr("ERROR: creating array for accessing 3D location grid buffer.");
            return (clean_memory(EXIT_ERROR_MEMORY));
        }
        LocGrid[ngrid].sum = 0.0;


        /* reset y-z dual-sheet grids (3D time grids) */

        for (narr = 0; narr < NumArrivalsLocation; narr++) {
            if (Arrival[narr].sheetdesc.type == GRID_TIME)
                Arrival[narr].sheetdesc.origx =
                    VERY_LARGE_DOUBLE;
        }



    } else if (SearchType == SEARCH_MET) {

        /* test change 17JAN2000 AJL */
        /*		InitializeMetropolisWalk(LocGrid + ngrid ,
                                Arrival, NumArrivalsLocation, &Metrop,
                                MetNumSamples, MetStepInit);
         */

        InitializeMetropolisWalk(LocGrid + ngrid,
                Arrival, NumArrivalsLocation, &Metrop,
                MetLearn + MetEquil, MetStepInit);

        /* allocate scatter array for saved samples */
        iSizeOfFdata = (1 + MetUse / MetSkip) * 4 * sizeof (float);
        if ((fdata = (float *) malloc(iSizeOfFdata)) == NULL) {
            nll_puterr("ERROR: creating array for scatter samples.");
            return (clean_memory(EXIT_ERROR_LOCATE));
        }
        //NumAllocations++;

    } else if (SearchType == SEARCH_OCTTREE) {

        if (0 && LocMethod == METH_OT_STACK && octtreeParams.use_stations_density) {
            snprintf(MsgStr, sizeof(MsgStr), "WARNING: LOCSEARCH use_stations_density disabled with LOCMETHOD OT_STACK.");
            nll_putmsg(1, MsgStr);
            octtreeParams.use_stations_density = 0;
        }

        // station density weighting
        if (octtreeParams.use_stations_density) {
            AveInterStationDistance = calcAveInterStationDistance(StationPhaseList, NumStationPhases);
            snprintf(MsgStr, sizeof(MsgStr), "Station Density Weight:  Ave Station Distance: %lf", AveInterStationDistance);
            nll_putmsg(1, MsgStr);
            if (AveInterStationDistance < SMALL_DOUBLE) { // should not get here
                nll_puterr("ERROR: cannot apply OctTree Station Density Weight: Ave Station Distance is zero!");
            }
            NumForceOctTreeStaDenWt = 0;
        }

        // initialize memory/arrays for regular, initial oct-tree search grid
        // this is an x, y, z array of oct-tree root nodes,
        // a true oct-tree is created at each of these roots
        octTree = InitializeOcttree(LocGrid + ngrid, &octtreeParams);
        //NumAllocations++;

        // allocate scatter array for saved samples
        iSizeOfFdata = octtreeParams.num_scatter * 4 * sizeof (float);
        iSizeOfFdata = (12 * iSizeOfFdata) / 10; // sample may be slightly larger than requested
        if ((fdata = (float *) malloc(iSizeOfFdata)) == NULL) {
            nll_puterr("ERROR: creating array for scatter samples.");
            return (clean_memory(EXIT_ERROR_LOCATE));
        }
        //NumAllocations++;

    }


    /* since sorted, reset companion indices */
    if (VpVsRatio > 0.0) {
        //for (narr = 0; narr < NumArrivalsLocation; narr++) {
        for (narr = 0; narr < NumArrivals; narr++) {
            if (Arrival[narr].n_companion < 0)
                continue;
            int n_companion_save = Arrival[narr].n_companion;
            // 20160805 AJL - bug fix  narr -> NumArrivals
            //             if (IsPhaseID(Arrival[narr].phase, "S") &&
            //        (Arrival[narr].n_companion = IsSameArrival(Arrival, narr, narr, "P")) < 0) {
            if (IsPhaseID(Arrival[narr].phase, "S") &&
                    (Arrival[narr].n_companion = IsSameArrival(Arrival, NumArrivals, narr, "P")) < 0) {
                //
                snprintf(MsgStr, sizeof(MsgStr), "ERROR: cannot find companion arrival: %s %s n_companion %d->%d", Arrival[narr].label, Arrival[narr].phase, n_companion_save, Arrival[narr].n_companion);
                nll_puterr(MsgStr);
                // DEBUG
                if (1) {
                    snprintf(MsgStr, sizeof(MsgStr), "Target:   narr %d %s label %s  time_grid_label %s", narr, Arrival[narr].phase, Arrival[narr].label, Arrival[narr].time_grid_label);
                    nll_puterr(MsgStr);
                    for (n = 0; n < NumArrivals; n++) {
                        snprintf(MsgStr, sizeof(MsgStr), "      narr %d %s label %s  time_grid_label %s", n, Arrival[n].phase, Arrival[n].label, Arrival[n].time_grid_label);
                        nll_puterr(MsgStr);
                    }
                }
                return (clean_memory(EXIT_ERROR_LOCATE));
            }
        }
    }


    /* do search */

    if (SearchType == SEARCH_GRID) {

        /* grid-search location (fill location grid) */
        if ((istat =
                LocGridSearch(ngrid, NumArrivals, NumArrivalsLocation,
                Arrival, LocGrid + ngrid, &Gauss, &Hypocenter)) < 0) {
            nll_puterr("ERROR: in grid search location.");
            return (clean_memory(EXIT_ERROR_LOCATE));
        }

    } else if (SearchType == SEARCH_MET) {

        /* Metropolis location (random walk) */
        if ((Hypocenter.nScatterSaved =
                LocMetropolis(ngrid, NumArrivals, NumArrivalsLocation,
                Arrival, LocGrid + ngrid,
                &Gauss, &Hypocenter, &Metrop, fdata)) < 0) {
            nll_puterr("ERROR: in Metropolis location.");
            return (clean_memory(EXIT_ERROR_LOCATE));
        }

    } else if (SearchType == SEARCH_OCTTREE) {

        /* do Octree location (importance sampling) */
        if ((Hypocenter.nScatterSaved =
                LocOctree(ngrid, NumArrivals, NumArrivalsLocation,
                Arrival, LocGrid + ngrid,
                &Gauss, &Hypocenter, &octtreeParams,
                octTree, fdata, &oct_node_value_max, &oct_tree_integral)) < 0) {
            nll_puterr("ERROR: in Octree location.");
            return (clean_memory(EXIT_ERROR_LOCATE));
        }

    }

    /* 20170911 moved below

        // clean up dates, calculate rms
        StdDateTime(Arrival, NumArrivals, &Hypocenter);

        // determine azimuth gaps
        double gap_secondary;
        Hypocenter.gap = CalcAzimuthGap(Arrival, NumArrivalsLocation, &gap_secondary);
        Hypocenter.gap_secondary = gap_secondary;

        // re-sort arrivals by distance
        if ((istat = SortArrivalsDist(Arrival, NumArrivals)) < 0) {
            nll_puterr("ERROR: sorting arrivals by distance.");
            return (clean_memory(EXIT_ERROR_LOCATE));
        }
     */

    /* search type dependent processing */

    if (SearchType == SEARCH_GRID && LocGridSave[ngrid]) {

        /* calculate confidence intervals and save to disk */

        if (LocGrid[ngrid].type == GRID_PROB_DENSITY) {
            if ((istat = CalcConfidenceIntrvl(LocGrid + ngrid, &Hypocenter, fnout)) < 0) {
                nll_puterr("ERROR: calculating confidence intervals.");
                return (clean_memory(EXIT_ERROR_LOCATE));
            }

            /* generate probabilistic scatter of events */
            char fnscatout[5 * MAXLINE]; // 20240129 AJL
            snprintf(fnscatout, sizeof(fnscatout), "%s.loc", fnout);
            if ((istat = GenEventScatterGrid(LocGrid + ngrid, &Hypocenter, &Scatter, fnscatout)) < 0) {
                nll_puterr("ERROR: calculating event scatter.");
            }

            /* calculate "traditional" statistics */
            Hypocenter.expect = CalcExpectation(LocGrid + ngrid, NULL);
            istat = rect2latlon(0, Hypocenter.expect.x, Hypocenter.expect.y, &(Hypocenter.expect_dlat), &(Hypocenter.expect_dlong));
            Hypocenter.cov = CalcCovariance(LocGrid + ngrid, &Hypocenter.expect, NULL);
            Hypocenter.ellipsoid = CalcErrorEllipsoid(&Hypocenter.cov, DELTA_CHI_SQR_68_3);
            Hypocenter.ellipse = CalcHorizontalErrorEllipse(&Hypocenter.cov, DELTA_CHI_SQR_68_2);

        } else {
            Hypocenter.probmax = -1.0;
        }

    } else if ((SearchType == SEARCH_MET || SearchType == SEARCH_OCTTREE) && LocGridSave[ngrid]) {

        if (octtreeParams.use_stations_density) {
            snprintf(MsgStr, sizeof(MsgStr), "Station Density Weight:  Number Force Divide: %d  max_num_nodes: %d", NumForceOctTreeStaDenWt, octtreeParams.max_num_nodes);
            nll_putmsg(1, MsgStr);
            if (NumForceOctTreeStaDenWt >= octtreeParams.max_num_nodes) {
                nll_puterr("ERROR: Number Force Divide > max_num_nodes !  Must reduce LOCSEARCH use_stations_density level.");
            } else if (NumForceOctTreeStaDenWt > (9 * octtreeParams.max_num_nodes) / 10) {
                nll_puterr("WARNING: Number Force Divide > 90% max_num_nodes !  Should reduce LOCSEARCH use_stations_density level.");
            }
        }


        if (SearchType == SEARCH_OCTTREE) {

            /*
            // determine integral of all oct-tree leaf node pdf values
            oct_tree_integral = integrateResultTree(resultTreeRoot, 0.0, oct_node_value_max);
            snprintf(MsgStr, sizeof(MsgStr), "Octree oct_node_value_max= %le oct_tree_integral= %le", oct_node_value_max, oct_tree_integral);
            nll_putmsg(1, MsgStr);*/

            // generate scatter sample
            if (Hypocenter.nScatterSaved == 0) // not saved during search
                Hypocenter.nScatterSaved = GenEventScatterOcttree(&octtreeParams, oct_node_value_max, fdata, oct_tree_integral, &Hypocenter);

        }

        /* write scatter file */
        if (iSaveNLLocEvent) {
            snprintf(fname, sizeof(fname), "%s.loc.scat", fnout);
            if ((fpio = fopen(fname, "w")) != NULL) {
                /* write scatter file header information */
                fseek(fpio, 0, SEEK_SET);
                fwrite(&(Hypocenter.nScatterSaved), sizeof (int), 1, fpio);
                ftemp = (float) Hypocenter.probmax;
                fwrite(&ftemp, sizeof (float), 1, fpio);
                /* skip header record */
                fseek(fpio, 4 * sizeof (float), SEEK_SET);
                /* write scatter samples */
                fwrite(fdata, 4 * sizeof (float), Hypocenter.nScatterSaved, fpio);
                fclose(fpio);
            } else {
                nll_puterr("ERROR: opening scatter output file.");
                return (clean_memory(EXIT_ERROR_IO));
            }
        }

        if (SearchType == SEARCH_OCTTREE && (iSaveNLLocOctree || return_oct_tree_grid)) {

            if (LocGrid[ngrid].type == GRID_PROB_DENSITY) {
                // convert oct tree values to likelihood
                oct_tree_prob_integral = convertOcttreeValuesToProbabilityDensity(
                        resultTreeRoot, VALUE_IS_LOG_PROB_DENSITY_IN_NODE, 0.0, oct_node_value_max);
                octTree->data_code = GRID_LIKELIHOOD;
                octTree->integral = oct_tree_prob_integral;
                //printf("DEBUG: oct_node_value_max %f  oct_tree_prob_integral %f\n", oct_node_value_max, oct_tree_prob_integral);
                // norm     // 20190626 AJL - added
                octTree->integral = normalizeProbabilityDensityOcttree(resultTreeRoot, 0.0, octTree->integral);
                //printf("DEBUG: octTree->integral %f\n", octTree->integral);
                // create new result tree sorted by node values only, without multiplication by volume
                //				resultTreeLikelihoodRoot = NULL;
                //				resultTreeLikelihoodRoot = createResultTree(resultTreeRoot, resultTreeLikelihoodRoot);
                snprintf(MsgStr, sizeof(MsgStr), "Oct tree structure converted to probability.");
                nll_putmsg(1, MsgStr);
                // convert oct tree values to confidence
                //convertOcttreeValuesToConfidence(resultTreeRoot, 0.0);
            }

            if (iSaveNLLocOctree) {
                // write oct tree structure to file
                snprintf(fname, sizeof(fname), "%s.loc.octree", fnout);
                if ((fpio = fopen(fname, "w")) != NULL) {
                    istat = writeTree3D(fpio, octTree);
                    //printf("DEBUG: write output oct tree file: %s\n", fname);
                    fclose(fpio);
                    snprintf(MsgStr, sizeof(MsgStr), "Oct tree structure written to file : %d nodes", istat);
                    nll_putmsg(1, MsgStr);
                } else {
                    nll_puterr("ERROR: opening oct tree structure output file.");
                    return (clean_memory(EXIT_ERROR_IO));
                }
            }
        }

        /* calculate "traditional" statistics */
        Hypocenter.expect = CalcExpectationSamples(fdata, Hypocenter.nScatterSaved);
        istat = rect2latlon(0, Hypocenter.expect.x, Hypocenter.expect.y, &(Hypocenter.expect_dlat), &(Hypocenter.expect_dlong));
        Hypocenter.cov = CalcCovarianceSamples(fdata, Hypocenter.nScatterSaved, &Hypocenter.expect);
        if (Hypocenter.nScatterSaved) {
            Hypocenter.ellipsoid = CalcErrorEllipsoid(&Hypocenter.cov, DELTA_CHI_SQR_68_3);
            Hypocenter.ellipse = CalcHorizontalErrorEllipse(&Hypocenter.cov, DELTA_CHI_SQR_68_2);
            snprintf(MsgStr, sizeof(MsgStr), "ellipsoid_volume = %le", (4.0 / 3.0) * cPI * 8.0 * Hypocenter.ellipsoid.len1 * Hypocenter.ellipsoid.len2 * Hypocenter.ellipsoid.len3);
            nll_putmsg(2, MsgStr);
        }
    }


    /* search type independent processing */

    // re-calculate solution and arrival statistics for expectation hypocenter in case that expectation results are to be saved
    // 20170811 AJL - added to allow saving of expectation hypocenter results instead of maximum likelihood
    if (iSaveNLLocExpectation) {

        // set hypocenter x,y,z to expectation
        Hypocenter.max_like.x = Hypocenter.x;
        Hypocenter.max_like.y = Hypocenter.y;
        Hypocenter.max_like.z = Hypocenter.z;
        istat = rect2latlon(0, Hypocenter.max_like.x, Hypocenter.max_like.y, &(Hypocenter.max_like_dlat), &(Hypocenter.max_like_dlong));
        int hypo_hour, hypo_min;
        hypotime2hrminsec(Hypocenter.time, &hypo_hour, &hypo_min, &(Hypocenter.max_like_sec));
        Hypocenter.x = Hypocenter.expect.x;
        Hypocenter.y = Hypocenter.expect.y;
        Hypocenter.z = Hypocenter.expect.z;

        double cell_diagonal_time_var_best = 0.0; // TODO: add to Grid Search ?
        double cell_diagonal_best = 0.0; // TODO: add to Grid Search ?
        double cell_volume_best = 0.0; // TODO: add to Grid Search ?
        double misfit_max = Hypocenter.grid_misfit_max; // maximum likelihood misfit_max set in previous call to SaveBestLocation
        SaveBestLocation(NULL, NumArrivals, NumArrivalsLocation, Arrival, LocGrid + ngrid,
                &Gauss, &Hypocenter, misfit_max, LocGrid[ngrid].type, 1, cell_diagonal_time_var_best, cell_diagonal_best, cell_volume_best);

    }

    // clean up dates, calculate rms
    StdDateTime(Arrival, NumArrivals, &Hypocenter);

    // determine azimuth gaps
    double gap_secondary;
    Hypocenter.gap = CalcAzimuthGap(Arrival, NumArrivalsLocation, &gap_secondary);
    Hypocenter.gap_secondary = gap_secondary;

    // re-sort arrivals by distance
    if ((istat = SortArrivalsDist(Arrival, NumArrivals)) < 0) {
        nll_puterr("ERROR: sorting arrivals by distance.");
        return (clean_memory(EXIT_ERROR_LOCATE));
    }


    // QML fields added for compatibility with QuakeML OriginQuality attributes (AJL 201005)
    int usedPhaseCount; // QML - Number of defining phases, i. e., phase observations that were actually used for computing
    // the origin. Note that there may be more than one defining phase per station.
    int associatedPhaseCount; // QML - Number of associated phases, regardless of their use for origin computation.
    int associatedStationCount; // QML - Number of stations at which the event was observed.
    int usedStationCount; // QML - Number of stations from which data was used for origin computation.
    int depthPhaseCount; // QML - Number of depth phases (typically pP, sometimes sP) used in depth computation.
    usedPhaseCount = CalcArrivalCounts(Arrival, NumArrivals, NumArrivalsRead, &associatedPhaseCount, &associatedStationCount, &usedStationCount, &depthPhaseCount);
    if (usedPhaseCount != Hypocenter.nreadings) {
        snprintf(MsgStr, sizeof(MsgStr), "ERROR: usedPhaseCount %d != Hypocenter.nreadings %d: this should not happen!\n", usedPhaseCount, Hypocenter.nreadings);
        nll_puterr(MsgStr);
    }
    //printf("DEBUG: usedPhaseCount %d  Hypocenter.nreadings %d\n", usedPhaseCount, Hypocenter.nreadings);
    //printf("DEBUG: associatedPhaseCount %d  associatedStationCount %d  usedStationCount %d  depthPhaseCount %d\n",
    //        associatedPhaseCount, associatedStationCount, usedStationCount, depthPhaseCount);
    Hypocenter.associatedPhaseCount = associatedPhaseCount;
    Hypocenter.associatedStationCount = associatedStationCount;
    Hypocenter.usedStationCount = usedStationCount;
    Hypocenter.depthPhaseCount = depthPhaseCount;

    // save distances
    // QML fields added for compatibility with QuakeML OriginQuality attributes (AJL 201005)
    double minimumDistance; // QML - Epicentral distance of station closest to the epicenter. Unit: km
    double maximumDistance; // QML - Epicentral distance of station farthest from the epicenter. Unit: km
    double medianDistance; // QML - Median epicentral distance of used stations. Unit: km
    minimumDistance = CalcArrivalDistances(Arrival, NumArrivalsLocation, &maximumDistance, &medianDistance, usedStationCount);
    //printf("DEBUG: minimumDistance %.1f  maximumDistance %.1f  medianDistance %.1f\n",
    //        minimumDistance, maximumDistance, medianDistance);
    Hypocenter.dist = minimumDistance;
    Hypocenter.minimumDistance = minimumDistance;
    Hypocenter.maximumDistance = maximumDistance;
    Hypocenter.medianDistance = medianDistance;

    // mist QML fields
    snprintf(Hypocenter.groundTruthLevel, sizeof(Hypocenter.groundTruthLevel), "%s", "-");

    // SED-ETH fields added for compatibility with legacy SED location quality indicators (AJL 201006)
    // algorithm from SH 29JUL2004
    if (iCalcSedOrigin) {
        /* determine quality factor
               A: RMS < 0.5 s; diff < 0.5 km; errh < 2.0 km & errz < 2.0 km
               B: RMS < 0.5 s; diff < 0.5 km; errh >= 2.0 km & errz >= 2.0 km
               C: RMS < 0.5 s; diff >= 0.5 km;
               C: RMS >= 0.5 s
         */
        // difference between maximum likelihood and expectation hypocenter locations
        double diff = sqrt(
                (Hypocenter.expect.x - Hypocenter.x) * (Hypocenter.expect.x - Hypocenter.x) +
                (Hypocenter.expect.y - Hypocenter.y) * (Hypocenter.expect.y - Hypocenter.y) +
                (Hypocenter.expect.z - Hypocenter.z) * (Hypocenter.expect.z - Hypocenter.z)
                );
        //printf("\nDEBUG: WriteSnapSum: diff maximum likelihood and expectation hypo: %f\n", diff);
        // std err x, y, z
        double errx = sqrt(Hypocenter.cov.xx);
        double erry = sqrt(Hypocenter.cov.yy);
        double errz = sqrt(Hypocenter.cov.zz);
        if (GeometryMode == MODE_GLOBAL) { //  GLOBAL - convert err to deg
            errx *= KM2DEG;
            erry *= KM2DEG;
            errz *= KM2DEG;
        }
        char qual = '-';
        if (Hypocenter.rms >= 0.5) {
            qual = 'D';
        } else if (diff > 0.5) {
            qual = 'C';
        } else if ((errx > 2.0 || erry > 2.0) && errz > 2.0) {
            qual = 'B';
        } else {
            qual = 'A';
        }
        Hypocenter.diffMaxLikeExpect = diff;
        Hypocenter.qualitySED = qual; // flags SED fields available - will be written to NLL Hypocenter .
    } else {
        Hypocenter.qualitySED = '\0'; // flags SED fields not available.
    }




    /* search type dependent results saving */

    if (SearchType == SEARCH_GRID) {

        /* save location grid to disk */

        if (!iSaveNone && LocGridSave[ngrid])
            if ((istat = WriteGrid3dBuf(LocGrid + ngrid, NULL,
                    fnout, "loc")) < 0) {
                nll_puterr("ERROR: writing location grid to disk.");
                return (clean_memory(EXIT_ERROR_IO));
            }
    } else if (SearchType == SEARCH_MET || SearchType == SEARCH_OCTTREE) {

        /* save location grid header to disk */

        if (!iSaveNone && LocGridSave[ngrid])
            if ((istat = WriteGrid3dHdr(LocGrid + ngrid, NULL, fnout, "loc")) < 0) {
                nll_puterr("ERROR: writing grid header to disk.");
                return (clean_memory(EXIT_ERROR_IO));
            }
    }



    /* display and save minimum misfit location to file */

    if (LocGridSave[ngrid]) {
        /* calculate magnitudes */
        // 20180907 AJL - following 4 lines moved to NLLoc() since may be modified when reading observations
        /*Hypocenter.amp_mag = MAGNITUDE_NULL;
        Hypocenter.num_amp_mag = 0;
        Hypocenter.dur_mag = MAGNITUDE_NULL;
        Hypocenter.num_dur_mag = 0;*/
        // 20121015 AJL - bug fix
        //for (n = 0; n < MAX_NUM_MAG_METHODS; n++)
        for (n = 0; n < NumMagnitudeMethods; n++)
            CalculateMagnitude(&Hypocenter, Arrival, NumArrivals,
                Component, NumCompDesc, Magnitude + n);
        /* calculate estimated VpVs ratio */
        CalculateVpVsEstimate(&Hypocenter, Arrival, NumArrivals);
        /* save location */
        if ((istat = SaveLocation(&Hypocenter, ngrid, fn_loc_obs, fnout, numArrivalsReject, "grid", 1, &Gauss)) < 0) {
            nll_puterr("ERROR: saving location.");
            return (clean_memory(istat));
        }
        /* add location to loclist */
        if (return_locations) {
            ploc_list_node = newLocation(
                    cloneHypoDesc(&Hypocenter),
                    cloneArrivalDescArray(Arrival, NumArrivals),
                    NumArrivals, cloneGridDesc(&(LocGrid[ngrid])),
                    return_oct_tree_grid ? octTree : NULL,
                    return_scatter_sample ? fdata : NULL);
            *ploc_list_head = addLocationToLocList(ploc_list_head, ploc_list_node, NumEventsLocated);
        }
        /* update station statistics table */
        if (
                ((LocGrid[ngrid].numz == 2 && strncmp(Hypocenter.locStat, "ABORTED", 7) != 0) // 20200812 AJL - Added so that station statistics will be accumulated when depth fixed (i.e. numz = 1)
                || strncmp(Hypocenter.locStat, "LOCATED", 7) == 0)
                && Hypocenter.rms <= RMS_Max
                && Hypocenter.nreadings >= NRdgs_Min
                && Hypocenter.gap <= Gap_Max
                && Hypocenter.ellipsoid.len3 <= Ell_Len3_Max
                && Hypocenter.z >= Hypo_Depth_Min
                && Hypocenter.z <= Hypo_Depth_Max) {
            UpdateStaStat(ngrid, Arrival, NumArrivals, P_ResidualMax, S_ResidualMax, Hypo_Dist_Max, 1.0);
            //printf("INSTALLED in Stat Table: ");
        } else {
            //printf("NOT INSTALLED in Stat Table: ");
        }
        //printf("Hypo: %s %f %d %d %f %f\n", Hypocenter.locStat, Hypocenter.rms, Hypocenter.nreadings, Hypocenter.gap, Hypocenter.ellipsoid.len3, Hypocenter.z);
    }



    /* search type dependent cleanup */

    if (SearchType == SEARCH_GRID) {

        /* free grid memory */

        DestroyGridArray(LocGrid + ngrid);
        FreeGrid(LocGrid + ngrid);

        /* intialize next grid origin location */

        if (ngrid < NumLocGrids - 1) {
            if (LocGrid[ngrid + 1].autox)
                LocGrid[ngrid + 1].origx = Hypocenter.x
                    - 0.5 * (double) (LocGrid[ngrid + 1].numx - 1)
                * LocGrid[ngrid + 1].dx;
            if (LocGrid[ngrid + 1].autoy)
                LocGrid[ngrid + 1].origy = Hypocenter.y
                    - 0.5 * (double) (LocGrid[ngrid + 1].numy - 1)
                * LocGrid[ngrid + 1].dy;
            if (LocGrid[ngrid + 1].autoz)
                LocGrid[ngrid + 1].origz = Hypocenter.z
                    - 0.5 * (double) (LocGrid[ngrid + 1].numz - 1)
                * LocGrid[ngrid + 1].dz;

            /* try to make sure new grid is inside initial grid */
            if (!IsGridInside(LocGrid + ngrid + 1, LocGrid, 1))
                nll_puterr(
                    "WARNING: cannot get next grid entirely contained inside 0th grid.");
        }

    } else if (SearchType == SEARCH_MET) {

        /* free saved samples memory */
        if (!return_scatter_sample) {
            free(fdata);
            fdata = NULL;
            //NumAllocations--;
        }

    } else if (SearchType == SEARCH_OCTTREE) {

        // free results tree - IMPORTANT!
        freeResultTree(resultTreeRoot);

        /* free oct-tree memory */
        if (!return_oct_tree_grid) {
            freeTree3D(octTree, 1);
            //NumAllocations--;
        }

        /* free saved samples memory */
        if (!return_scatter_sample) {
            free(fdata);
            fdata = NULL;
            //NumAllocations--;
        }

    }

    // GetNLLoc_PdfGrid
    if (SearchPrior.coherence != NULL) {
        free(SearchPrior.coherence);
        SearchPrior.coherence = NULL;
    }
    if (SearchPrior.weight != NULL) {
        free(SearchPrior.weight);
    }
    if (SearchPosterior.coherence != NULL) {
        free(SearchPosterior.coherence);
        SearchPosterior.coherence = NULL;
    }
    if (SearchPosterior.weight != NULL) {
        free(SearchPosterior.weight);
        SearchPosterior.weight = NULL;
    }
    if (SearchPosterior.nfirst_motion_arrivals != NULL) {
        free(SearchPosterior.nfirst_motion_arrivals);
        SearchPosterior.nfirst_motion_arrivals = NULL;
    }
    if (SearchPosterior.first_motion_arrivals != NULL) {
        for (int i = 0; i < SearchPosterior.nGrids; i++) {
            free(SearchPosterior.first_motion_arrivals[i]);
        }
        free(SearchPosterior.first_motion_arrivals);
        SearchPosterior.first_motion_arrivals = NULL;
    }

    if (SearchPosterior.tree3D != NULL) {
        for (int i = 0; i < SearchPosterior.nGrids; i++) {
            freeTree3D(SearchPosterior.tree3D[i], 1);
        }
        free(SearchPosterior.tree3D);
        SearchPosterior.tree3D = NULL;
    }



    clean_memory(0);


    /* re-sort to get location arrivals in time order */

    if ((istat =
            SortArrivalsIgnore(Arrival, NumArrivals)) < 0) {
        nll_puterr(
                "ERROR: sorting arrivals by ignore flag.");
        return (clean_memory(EXIT_ERROR_LOCATE));
    }
    if ((istat = SortArrivalsTime(Arrival, NumArrivalsLocation)) < 0) {
        nll_puterr("ERROR: sorting arrivals by time.");
        return (EXIT_ERROR_LOCATE);
    }



    /* search type dependent return */

    if (SearchType == SEARCH_GRID) {

        return (0);

    } else if (SearchType == SEARCH_MET || SearchType == SEARCH_OCTTREE) {

        if (Hypocenter.nScatterSaved == 0)
            return (1);

        return (0);
    }


    return (0);


}

/** function to do misc memory cleanup */

int clean_memory(int istat) {

    // AJL - 20080710 (valgrind)
    // free EDT_OT_WT memory
    if (ot_ml_arrival != NULL)
        free(ot_ml_arrival);
    ot_ml_arrival = NULL;
    if (ot_ml_arrival_edt_sum != NULL)
        free(ot_ml_arrival_edt_sum);
    ot_ml_arrival_edt_sum = NULL;
    isize_ot_ml_array = 0;

    return (istat);

}

/** function to initialize Metropolis walk */

void InitializeMetropolisWalk(GridDesc* ptgrid, ArrivalDesc* parrivals, int
        numArrLoc, WalkParams* pMetrop, int numSamples, double initStep) {
    int narr;
    double xlen, ylen, zlen, dminlen;
    double xmin, xmax, ymin, ymax;
    SourceDesc* pstation;


    /* set walk limits equal to grid limits */
    xmin = ptgrid->origx;
    xmax = xmin + (double) (ptgrid->numx - 1) * ptgrid->dx;
    ymin = ptgrid->origy;
    ymax = ymin + (double) (ptgrid->numy - 1) * ptgrid->dy;

    /* find station with earliest arrival and non-zero weight */
    narr = 0;
    while (narr < numArrLoc && parrivals[narr].weight < 0.001)
        narr++;

    /* initialize walk location */
    if (narr < numArrLoc)
        pstation = &(parrivals[narr].station);
    if (narr < numArrLoc &&
            pstation->x >= xmin && pstation->x <= xmax
            && pstation->y >= ymin && pstation->y <= ymax) {
        /* start walk at location of station with earliest arrival */
        pMetrop->x = pstation->x;
        pMetrop->y = pstation->y;
    } else {
        /* start walk at grid center */
        pMetrop->x = ptgrid->origx
                + (double) (ptgrid->numx - 1) * ptgrid->dx / 2.0;
        pMetrop->y = ptgrid->origy
                + (double) (ptgrid->numy - 1) * ptgrid->dy / 2.0;
    }
    /* start walk at grid center depth */
    pMetrop->z = ptgrid->origz
            + (double) (ptgrid->numz - 1) * ptgrid->dz / 2.0;

    /* calculate initial step size */
    if (initStep < 0.0) {
        xlen = (double) ptgrid->numx * ptgrid->dx / 2.0;
        ylen = (double) ptgrid->numy * ptgrid->dy / 2.0;
        zlen = (double) ptgrid->numz * ptgrid->dz / 2.0;
        dminlen = xlen < ylen ? xlen : ylen;
        dminlen = dminlen < zlen ? dminlen : zlen;
        /* step is size that tiles plane parallel to max len sides */
        pMetrop->dx = sqrt((xlen * ylen * zlen / dminlen) / (double) numSamples);
        /* step is size that tiles search volume */
        /*pMetrop->dx = pow(xlen * ylen * zlen / (double) numSamples, 1.0/3.0);*/
    } else {
        pMetrop->dx = initStep;
    }

    if (message_flag >= 4) {
        snprintf(MsgStr, sizeof(MsgStr),
                "INFO: Metropolis initial step size: %lf", pMetrop->dx);
        nll_putmsg(4, MsgStr);
    }

    /* set likelihood */
    pMetrop->likelihood = -1.0;

}

static int save_location_count = 0;

/** function to display and save minimum misfit location to file */

int SaveLocation(HypoDesc* hypo, int ngrid, char* fnobs, char *fnout, int numArrivalsReject,
        char* loctypename, int isave_phases, GaussLocParams * gauss_par) {
    int istat;
    char *pchr;
    //char sys_command[2 * FILENAME_MAX];
    char sourcefname[FILENAME_MAX], targetfname[2 * FILENAME_MAX];
    char fname[3 * FILENAME_MAX], frootname[2 * FILENAME_MAX];
    FILE *fp_tmp;

    /* set signature string */
    snprintf(hypo->signature, sizeof(hypo->signature), "%s   obs:%s   %s:v%s(%s)   run:%s",
            LocSignature, fnobs, prog_name, PVER, PDATE, CurrTimeStr());
    while ((pchr = strchr(hypo->signature, '\n')))
        *pchr = ' ';

    /* display hypocenter to std out */
    if (message_flag >= 3)
        WriteLocation(stdout, hypo, Arrival,
            NumArrivals + numArrivalsReject, fnout,
            isave_phases, 1, 0, LocGrid + ngrid, 0);

    /*  save requested hypocenter/phase formats */

#ifdef CUSTOM_ETH
    /* SH 02/26/2004
            added new routine WriteSnapSum to write hypocenter summary to file (SNAP format) */
    // must call WriteSnapSum before other output because it sets ETH magnitudes
    if (iSaveSnapSum) {
        WriteSnapSum(NULL, hypo, Arrival, NumArrivals);
    }
#endif

    if (iSaveNLLocEvent) {
        /* write NLLoc hypocenter to event file */
        snprintf(frootname, sizeof(frootname), "%s.loc", fnout);
        snprintf(fname, sizeof(fname), "%s.hyp", frootname);
        if ((istat = WriteLocation(NULL, hypo, Arrival,
                NumArrivals + numArrivalsReject, fname, isave_phases, 1, 0,
                LocGrid + ngrid, 0)) < 0) {
            nll_puterr("ERROR: writing location to event file.");
            return (EXIT_ERROR_IO);
        }
        /* copy event file to last.hyp */
        /* Modified Jan Wiszniowski 2022-02-02
        snprintf(sys_command, sizeof(sys_command), "cp %s %slast.hyp", fname, f_outpath);
        system(sys_command);
         */
        snprintf(targetfname, sizeof(targetfname), "%slast.hyp", f_outpath);
        copy_file(fname, targetfname);
        /**/
        snprintf(fname, sizeof(fname), "%s.hdr", frootname);
        /* Modified Jan Wiszniowski 2022-02-02
        snprintf(sys_command, sizeof(sys_command), "cp %s %slast.hdr", fname, f_outpath);
        system(sys_command);
         */
        snprintf(targetfname, sizeof(targetfname), "%slast.hdr", f_outpath);
        copy_file(fname, targetfname);
        /**/ snprintf(fname, sizeof(fname), "%s.scat", frootname);
        if ((fp_tmp = fopen(fname, "rb")) != NULL) {
            fclose(fp_tmp);
            /* Modified Jan Wiszniowski 2022-02-02
            snprintf(sys_command, sizeof(sys_command), "cp %s %slast.scat", fname, f_outpath);
            system(sys_command);
             */
            snprintf(targetfname, sizeof(targetfname), "%slast.scat", f_outpath);
            copy_file(fname, targetfname);
            /**/
        }
    }


    // 20220131 AJL - added to support JSON output of location results
    if (iSaveNLLocEvent_JSON) {

#ifdef _GNU_SOURCE

        printf("DEBUG: iSaveNLLocEvent_JSON\n");
        // write NLLoc hypocenter to memory stream, convert stream to JSON

        // GNU C library extensions to support memory streams (function open_memstream).
        char *bp_memory_stream = NULL;
        size_t memory_stream_size;
        FILE *fp_memory_stream = NULL;
        fp_memory_stream = open_memstream(&bp_memory_stream, &memory_stream_size);
        if (fp_memory_stream == NULL) {
            nll_puterr("ERROR: Cannot write NLLoc hypocenter-phase file to JSON: GNU C library extensions needed to support memory streams (function open_memstream).");
        } else {
            // write location to memory stream
            if ((istat = WriteLocation(fp_memory_stream, hypo, Arrival,
                    NumArrivals + numArrivalsReject, NULL, isave_phases, 1, 0,
                    LocGrid + ngrid, 0)) < 0) {
                nll_puterr("ERROR: writing location to memory stream, Cannot write NLLoc hypocenter-phase file to JSON.");
            } else {
                printf("DEBUG: iSaveNLLocEvent_JSON: wrote location to memory stream\n");
                // convert to JSON
                snprintf(frootname, sizeof(frootname), "%s.loc", fnout);
                snprintf(fname, sizeof(fname), "%s.hyp.json", frootname);
                FILE *fp_json_out = NULL;
                if ((fp_json_out = fopen(fname, "w")) == NULL) {
                    nll_puterr2("ERROR: opening hypocenter JSON output file", fname);
                } else {
                    fflush(fp_memory_stream);
                    rewind(fp_memory_stream);
                    printf("DEBUG: iSaveNLLocEvent_JSON: json_write_NLL_location\n");
                    json_write_NLL_location(bp_memory_stream, memory_stream_size, fp_json_out);
                    fclose(fp_json_out);
                }

                // send message

            }

            // cleanup
            fclose(fp_memory_stream);

        }

#else
        nll_puterr("ERROR: Cannot write NLLoc hypocenter-phase file to JSON: GNU C library extensions needed to support memory streams (function open_memstream(); see compiler define _GNU_SOURCE).");
#endif

    }

    if (iSaveNLLocSum) {
        /* write NLLoc hypocenter to summary file */
        if ((istat = WriteLocation(pSumFileHypNLLoc[ngrid],
                hypo,
                Arrival, NumArrivals, fnout, 0, 1, 0,
                LocGrid + ngrid, 0)) < 0) {
            nll_puterr("ERROR: writing location to summary file.");
            return (EXIT_ERROR_IO);
        }
        fflush(pSumFileHypNLLoc[ngrid]);
        /* copy event grid header to .sum header */
        /* Modified Jan Wiszniowski 2022-02-02
        snprintf(sys_command, sizeof(sys_command),
                "cp %s.loc.hdr %s.sum.%s%d.loc.hdr",
                fnout, fn_path_output, loctypename, ngrid);
        system(sys_command);
         */
        snprintf(sourcefname, sizeof(sourcefname), "%s.loc.hdr", fnout);
        snprintf(targetfname, sizeof(targetfname), "%s.sum.%s%d.loc.hdr", fn_path_output, loctypename, ngrid);
        copy_file(sourcefname, targetfname);
        /**/
    }

    if (iSaveNLLocSumCSV) {
        /* write NLLoc hypocenter to CSV summary file */
        if ((istat = WriteLocationCSV(pSumFileHypNLLocCSV[ngrid], hypo, fnout)) < 0) {
            nll_puterr("ERROR: writing location CSV to summary file.");
            return (EXIT_ERROR_IO);
        }
        fflush(pSumFileHypNLLocCSV[ngrid]);
    }

    if (iSaveHypo71Event) {
        /* write HYPO71 hypocenter to event file */
        WriteHypo71(NULL, hypo, Arrival, NumArrivals, fnout, 1, 1);
    }
    if (iSaveHypo71Sum) {
        /* write HYPO71 hypocenter to summary file */
        WriteHypo71(pSumFileHypo71[ngrid], hypo,
                Arrival, NumArrivals, fnout, iWriteHypHeader[ngrid], 0);
    }

    if (iSaveHypoEllEvent) {
        /* write pseudo-HypoEllipse hypo to event file */
        WriteHypoEll(NULL, hypo, Arrival, NumArrivals, fnout, 1, 1);
    }
    if (iSaveHypoEllSum) {
        /* write pseudo-HypoEllipse hypo to summary file */
        WriteHypoEll(pSumFileHypoEll[ngrid], hypo,
                Arrival, NumArrivals, fnout, iWriteHypHeader[ngrid], 0);
    }

    if (iSaveHypoInvSum) {
        /* write HypoInverseArchive hypocenter to summary file */
        WriteHypoInverseArchive(pSumFileHypoInv[ngrid], hypo, Arrival, NumArrivals,
                fnout, 0, 1, gauss_par->arrivalWeightMax);
        /* also write to last.hypo_inv */
        snprintf(fname, sizeof(fname), "%slast.hypo_inv", f_outpath);
        if ((fp_tmp = fopen(fname, "w")) != NULL) {
            WriteHypoInverseArchive(fp_tmp, hypo, Arrival, NumArrivals,
                    fnout, 0, 1, gauss_par->arrivalWeightMax);
            fclose(fp_tmp);
        }
    }

    if (iSaveHypoInvY2KArc) {
        /* write HypoInverseArchive hypocenter to summary file */
        WriteHypoInverseArchive(pSumFileHypoInvY2K[ngrid], hypo, Arrival, NumArrivals,
                fnout, 1, 1, gauss_par->arrivalWeightMax);
        /* also write to last.arc */
        snprintf(fname, sizeof(fname), "%slast.arc", f_outpath);
        if ((fp_tmp = fopen(fname, "w")) != NULL) {
            WriteHypoInverseArchive(fp_tmp, hypo, Arrival, NumArrivals,
                    fnout, 1, 1, gauss_par->arrivalWeightMax);
            fclose(fp_tmp);
        }
    }

    if (iSaveAlberto4Sum) {
        /* write Alberto 4 SIMULPS format */
        WriteHypoAlberto4(pSumFileAlberto4[ngrid], hypo, Arrival, NumArrivals, fnout);
    }

    if (iSaveFmamp) {
        // check for special processing of arrivals for fmamp output
        // 20200829 AJL - added to support posterior location with multiple events and observed fm polarities for the same station+phase
        if (iUseSearchPosterior && SearchPosterior.first_motion_arrivals != NULL && SearchPosterior.nfirst_motion_arrivals != NULL) {
            SearchPdfGridDesc *searchPdfGrid = &SearchPosterior;
            // write fmamp format with combined event arrivals
            WriteHypoFmampSearchPosterior(searchPdfGrid, pSumFileFmamp[ngrid], hypo, fnout, save_location_count < 1);
        } else {
            // write fmamp format with event arrivals
            WriteHypoFmamp(pSumFileFmamp[ngrid], hypo, Arrival, NumArrivals, fnout, save_location_count < 1);
        }
    }

    iWriteHypHeader[ngrid] = 0;

    save_location_count++;

    return (0);

}

/** function to combine arrivals and accumulate weighted first motion for multiple SearchPosterior events
 *
 * 20200829 AJL - added to support posterior location with multiple events and observed fm polarities for the same station+phase
 */

int WriteHypoFmampSearchPosterior(SearchPdfGridDesc *searchPdfGrid, FILE *fpio, HypoDesc* phypo, char* filename, int write_header) {

    ArrivalDesc *pfmarrivals;
    if ((pfmarrivals = (ArrivalDesc *) calloc(MAX_NUM_ARRIVALS, sizeof (ArrivalDesc))) == NULL) {
        nll_puterr("ERROR: allocating memory for temporary first_motion_arrivals for writing fmamp.");
        return (-1);
    }
    int nfmarrivals = 0;


    // allocate weight_sum
    double *weight_sum;
    if ((weight_sum = (double *) malloc(MAX_NUM_ARRIVALS * sizeof (double))) == NULL) {
        nll_puterr("ERROR: allocating memory for weight_sum for writing fmamp.");
        return (-1);
    }
    // allocate nweight
    double *fm_weight_sum;
    if ((fm_weight_sum = (double *) malloc(MAX_NUM_ARRIVALS * sizeof (double))) == NULL) {
        nll_puterr("ERROR: allocating memory for nweight for writing fmamp.");
        return (-1);
    }

    // write fmamp format with combined event arrivals
    ArrivalDesc *parr, *pfmarr;
    for (int nFile = 0; nFile < searchPdfGrid->nGrids; nFile++) { // loop over events
        for (int narr = 0; narr < searchPdfGrid->nfirst_motion_arrivals[nFile]; narr++) { // loop over event arrivals
            parr = searchPdfGrid->first_motion_arrivals[nFile] + narr;
            int ifound = -1;
            // check if station+phase is already in fm arrivals array
            for (int nfmarr = 0; nfmarr < nfmarrivals; nfmarr++) {
                pfmarr = pfmarrivals + nfmarr;
                if (!strcmp(parr->label, pfmarr->label) && !strcmp(parr->phase, pfmarr->phase)) { // compare on station+phase
                    ifound = nfmarr;
                    break;
                }
            }
            if (ifound < 0) { // station+phase not yet in array, initialize
                pfmarrivals[nfmarrivals] = *parr;
                weight_sum[nfmarrivals] = 0.0;
                fm_weight_sum[nfmarrivals] = 0.0;
                ifound = nfmarrivals;
                nfmarrivals++;
            }
            // get first motion
            int ifm = 0;
            if (strstr("CcUu+", parr->first_mot)) { // follows fmamp conventions in fmamp/read_input.c
                ifm = 1;
            } else if (strstr("DdRr-", parr->first_mot)) { // follows fmamp conventions in fmamp/read_input.c
                ifm = -1;
            } else {

                continue; // no first motion, skip arrival
            }
            // update stats for this station+phase
            // weight is polarity * grid weight
            weight_sum[ifound] += searchPdfGrid->weight[nFile];
            fm_weight_sum[ifound] += (double) ifm * searchPdfGrid->weight[nFile];
            //printf("DEBUG: ifm %d  searchPdfGrid->weight[nFile] %f\n", ifm, searchPdfGrid->weight[nFile]);
        }
    }

    // set final arrival first motion, fm weight and take-off angles for each station+phase
    char fileroot[4 * MAXLINE];
    for (int nfmarr = 0; nfmarr < nfmarrivals; nfmarr++) {
        double first_motion = 0.0;
        if (weight_sum[nfmarr] > FLT_MIN) {
            first_motion = fm_weight_sum[nfmarr] / weight_sum[nfmarr];
        }
        pfmarr = pfmarrivals + nfmarr;
        if (first_motion >= 0.0) {
            snprintf(pfmarr->first_mot, sizeof(pfmarr->first_mot), "%s", "+");
        } else {
            snprintf(pfmarr->first_mot, sizeof(pfmarr->first_mot), "%s", "-");
        }
        pfmarr->first_mot_quality = fabs(first_motion);
        //printf("DEBUG: first_motion %f  first_mot_quality %f\n", first_motion, pfmarr->first_mot_quality);
        // set take-off angles at search posterior hypocenter
        // try to open time grid file using original phase ID
        EvaluateArrivalAlias(pfmarr);
        snprintf(fileroot, sizeof(fileroot), "%s.%s.%s.angle", fn_time_grids[0], pfmarr->phase, pfmarr->time_grid_label);
        int iavailable;
        // need to get grid type from file on disk  TODO: this could be integrated into ReadTakeOffAnglesFile function
        FILE *fp_grid, *fp_hdr;
        GridDesc gdesc;
        if (OpenGrid3dFile(fileroot, &fp_grid, &fp_hdr, &gdesc, "angle", NULL, iSwapBytesOnInput) < 0) {
            if (message_flag >= 3) {
                snprintf(MsgStr, sizeof(MsgStr), "WARNING: cannot open angle grid file, ignoring angles: %s", fileroot);
                nll_putmsg(3, MsgStr);
            }
            //angles = SetTakeOffAngles(0.0, 0.0, 0);
            //GetTakeOffAngles(&angles, &(pfmarr->ray_azim), &(pfmarr->ray_dip), &(pfmarr->ray_qual));
            iavailable = -1;
        } else {
            //printf("DEBUG: pfmarr->gdesc.type %d  GRID_ANGLE %d\n", gdesc.type, GRID_ANGLE);
            int gdesc_type = gdesc.type;
            CloseGrid3dFile(&gdesc, &fp_grid, &fp_hdr);
            if (gdesc_type == GRID_ANGLE) {
                // 3D grid
                iavailable = ReadTakeOffAnglesFile(fileroot,
                        phypo->x, phypo->y, phypo->z,
                        &(pfmarr->ray_azim),
                        &(pfmarr->ray_dip),
                        &(pfmarr->ray_qual), -1.0, iSwapBytesOnInput);
            } else {
                // 2D grid (1D model)
                iavailable = ReadTakeOffAnglesFile(fileroot,
                        0.0,
                        GeometryMode == MODE_GLOBAL ? pfmarr->dist * KM2DEG : pfmarr->dist,
                        phypo->z,
                        &(pfmarr->ray_azim),
                        &(pfmarr->ray_dip),
                        &(pfmarr->ray_qual), pfmarr->azim, iSwapBytesOnInput);
            }
        }
        if (iavailable < 0) { // angles grids not available
            pfmarr->first_mot_quality = 0.0;
        }
        // check some things
        if (pfmarr->ray_azim < 0.0 || pfmarr->ray_azim > 360.0 || pfmarr->ray_dip < 0.0 || pfmarr->ray_dip > 180.0) {
            pfmarr->first_mot_quality = 0.0;
        }
        //printf("DEBUG: ReadTakeOffAnglesFile: fn_loc_grids %s\n", fn_loc_grids);
        //printf("DEBUG: ReadTakeOffAnglesFile: pfmarr->time_grid_label %s\n", pfmarr->time_grid_label);
        //printf("DEBUG: ReadTakeOffAnglesFile: fileroot %s\n", fileroot);
        //printf("DEBUG: ReadTakeOffAnglesFile: iavailable %d  iAngleQualityMin %d\n", iavailable, iAngleQualityMin);
        //printf("DEBUG: ReadTakeOffAnglesFile: x y z %f %f %f  az %f  raz %f  rdip %f rq  %d\n", phypo->x, phypo->y, phypo->z, pfmarr->azim, pfmarr->ray_azim, pfmarr->ray_dip, pfmarr->ray_qual);

    }

    WriteHypoFmamp(fpio, phypo, pfmarrivals, nfmarrivals, filename, write_header);

    free(pfmarrivals);
    free(weight_sum);
    free(fm_weight_sum);

    return (nfmarrivals);

}


/** function to check if a date is reasonable */

int IsGoodDate(int iyear, int imonth, int iday) {
    if (iyear >= SMALLEST_EVENT_YEAR && iyear <= LARGEST_EVENT_YEAR
            && imonth > 0 && imonth < 13
            && iday > 0 && iday < 32)

        return (1);

    return (0);
}

/** function to homogenize date / time of arrivals */

int HomogDateTime(ArrivalDesc *arrival, int num_arrivals, HypoDesc * phypo) {
    int narr;
    int dofymin = 10000, yearmin = 10000;
    int test_month, test_day;

    for (narr = 0; narr < num_arrivals; narr++) {
        if (arrival[narr].year < yearmin)
            yearmin = arrival[narr].year;
        // AJL 20060615 - now allow crossing of year boundary (requires that first reading is earlier year, etc...)
        if (arrival[narr].year != yearmin) {
            // ok if Dec 31 -> Jan 01
            if ((arrival[narr].year == yearmin + 1)
                    && (arrival[narr].month == 1) && (arrival[narr].day == 1)) {
                arrival[narr].year = yearmin;
                arrival[narr].month = 12;
                arrival[narr].day = 31;
                arrival[narr].hour += 24;
            } else {
                return (OBS_FILE_ARRIVALS_CROSS_YEAR_BOUNDARY);
            }
        }
        arrival[narr].day_of_year =
                DayOfYear(arrival[narr].year, arrival[narr].month, arrival[narr].day);
        if (arrival[narr].day_of_year < dofymin)
            dofymin = arrival[narr].day_of_year;
    }

    for (narr = 0; narr < num_arrivals; narr++) {
        if (arrival[narr].day_of_year > dofymin) {
            arrival[narr].day_of_year--;
            arrival[narr].day--;
            arrival[narr].hour += 24;
        }
    }

    for (narr = 0; narr < num_arrivals; narr++)
        arrival[narr].obs_time = (long double) arrival[narr].sec
            //			- (long double) arrival[narr].delay	// DELAY_CORR	- incorporating delay so subtract (Tcorr = Tobs - (O-C))
            + 60.0L * ((long double) arrival[narr].min
            + 60.0L * (long double) arrival[narr].hour);

    if (!FixOriginTimeFlag) {
        /* initialize hypocenter year/month/day if origin time not fixed */
        phypo->year = yearmin;
        MonthDay(yearmin, dofymin, &(phypo->month), &(phypo->day));
    } else {
        /* homogenize hypocenter otime if origin time fixed */
        MonthDay(yearmin, dofymin, &test_month, &test_day);
        if (phypo->year != yearmin || test_month != phypo->month
                || test_day != phypo->day) {
            nll_puterr(
                    "ERROR: earliest arrivals year/month/day does not match fixed origin time year/month/day, ignoring observation set.");

            return (OBS_FILE_ARRIVALS_CROSS_YEAR_BOUNDARY);
        }
        phypo->time = (long double) phypo->sec
                + 60.0L * ((long double) phypo->min
                + 60.0L * (long double) phypo->hour);
        phypo->min = 0;
        phypo->hour = 0;
    }

    return (0);


}

/** function to check for arrivals with no absolute timing */

int CheckAbsoluteTiming(ArrivalDesc *arrival, int num_arrivals) {
    int narr;
    int nNoAbs = 0;

    for (narr = 0; narr < num_arrivals; narr++) {
        if (arrival[narr].inst[0] == '*') {
            arrival[narr].abs_time = 0;
            nNoAbs++;
        } else {

            arrival[narr].abs_time = 1;
        }

    }

    return (nNoAbs);


}

int hypotime2hrminsec(long double phypo_time, int *phypo_hour, int *phypo_min, double *phypo_sec) {

    long double hyp_time_tmp = phypo_time;
    *phypo_hour = (int) (hyp_time_tmp / 3600.0L);
    hyp_time_tmp -= (long double) *phypo_hour * 3600.0L;
    *phypo_min = (int) (hyp_time_tmp / 60.0L);
    hyp_time_tmp -= (long double) *phypo_min * 60.0L;
    *phypo_sec = (double) hyp_time_tmp;

    // TODO: check for case of origin time in previous day and correct date and time (currently get negative otime min, sec)  20150716 AJL

    return (0);

}

/** function to regenerate hypo date and time fields, robust to negative sec in hypo time */

int reset_hypodatetime(long double phypo_time, HypoDesc *phypo) {

    long double hyp_time_tmp = phypo_time;
    phypo->hour = (int) (hyp_time_tmp / 3600.0L);
    hyp_time_tmp -= (long double) phypo->hour * 3600.0L;
    phypo->min = (int) (hyp_time_tmp / 60.0L);
    hyp_time_tmp -= (long double) phypo->min * 60.0L;
    phypo->sec = (double) hyp_time_tmp;

    if (phypo->sec < 0.0) {

        printf("DEBUG: reset_hypodatetime: phypo time in: %4.4d-%2.2d-%2.2dT%2.2d:%2.2d:%lf\n", phypo->year, phypo->month, phypo->day, phypo->hour, phypo->min, phypo->sec);

        double hypo_sec = -phypo->sec;
        int int_sec_offset = (int) hypo_sec + 1;
        double dec_sec_corr = (double) int_sec_offset - hypo_sec;

        struct tm hypo_time;
        hypo_time.tm_year = phypo->year - 1900;
        hypo_time.tm_mon = phypo->month - 1;
        hypo_time.tm_mday = phypo->day;
        hypo_time.tm_hour = phypo->hour;
        hypo_time.tm_min = phypo->min;
        hypo_time.tm_sec = 0;

        // 20260619 AJL - Bug fix: do time conversion in UTC, not local timezone
        //20260619 time_t time_seconds = mktime(&hypo_time);
        time_t time_seconds = timegm(&hypo_time);
        // 20260619 AJL - Bug fix: do not add 1 hour
        //20260619time_seconds += 3600 - int_sec_offset;
        time_seconds -= int_sec_offset;

        struct tm *new_hypo_time = gmtime(&time_seconds);

        phypo->year = new_hypo_time->tm_year + 1900;
        phypo->month = new_hypo_time->tm_mon + 1;
        phypo->day = new_hypo_time->tm_mday;
        phypo->hour = new_hypo_time->tm_hour;
        phypo->min = new_hypo_time->tm_min;
        phypo->sec = (double) new_hypo_time->tm_sec + dec_sec_corr;

        printf("DEBUG: reset_hypodatetime: phypo time out: %4.4d-%2.2d-%2.2dT%2.2d:%2.2d:%lf\n", phypo->year, phypo->month, phypo->day, phypo->hour, phypo->min, phypo->sec);

    }

    return (0);

}

/** function to standardize date / time of arrivals, calculate rms */

int StdDateTime(ArrivalDesc *arrival, int num_arrivals, HypoDesc * phypo) {
    int narr;
    double rms_resid = 0.0, weight_sum = 0.0;
    long double sec_tmp;


    for (narr = 0; narr < num_arrivals; narr++) {
        // calc obs travel time and residual for arrivals with finite travel time and abs timing
        // AJL 20050926 bug fix - do not update rms with arrivals with no travel time!
        // if (arrival[narr].abs_time) {
        //printf("DEBUG: abs_time: %d, tt_Pred: %9.4le\n", arrival[narr].abs_time, arrival[narr].pred_travel_time);
        if (arrival[narr].abs_time && arrival[narr].pred_travel_time > 0.0) {
            arrival[narr].obs_travel_time =
                    arrival[narr].obs_time - phypo->time;
            arrival[narr].residual = arrival[narr].obs_travel_time -
                    arrival[narr].pred_travel_time;
            //printf("DEBUG: residual: %9.4le, tt_obs: %9.4le, tt_Pred: %9.4le, ", arrival[narr].residual, arrival[narr].obs_travel_time, arrival[narr].pred_travel_time);
            rms_resid += arrival[narr].weight *
                    arrival[narr].residual * arrival[narr].residual;
            weight_sum += arrival[narr].weight;
        } else {
            arrival[narr].obs_travel_time = 0.0;
            arrival[narr].residual = 0.0;
        }
        /* convert time to year/month/day/hour/min */
        // DELAY_CORR		sec_tmp = arrival[narr].obs_time;
        // removing delay so add (Tobs = Tcorr + (O-C))
        sec_tmp = arrival[narr].obs_time + arrival[narr].delay;
        arrival[narr].hour = (int) (sec_tmp / 3600.0L);
        sec_tmp -= (long double) arrival[narr].hour * 3600.0L;
        arrival[narr].min = (int) (sec_tmp / 60.0L);
        sec_tmp -= (long double) arrival[narr].min * 60.0L;
        arrival[narr].sec = (double) sec_tmp;
        MonthDay(arrival[narr].year, arrival[narr].day_of_year,
                &(arrival[narr].month), &(arrival[narr].day));
    }

    // set rms if not set earlier in SaveBestLocation
    if (phypo->rms < 0.0) {
        phypo->rms = 999.99;
        if (weight_sum > 0.0)
            phypo->rms = sqrt(rms_resid / weight_sum); // 20150324 AJL - TODO: should be mean of residuals for METH_L1_NORM ???
    }

    hypotime2hrminsec(phypo->time, &(phypo->hour), &(phypo->min), &(phypo->sec));
    //printf("DEBUG: StdDateTime: phypo->time %4.4d-%2.2d-%2.2dT%2.2d:%2.2d:%lf\n", phypo->year, phypo->month, phypo->day, phypo->hour, phypo->min, phypo->sec);
    // 20250207 AJL - Bug fix: if negative time correct data/time
    if (phypo->sec < 0.0) {

        reset_hypodatetime(phypo->time, phypo);
    }

    /*hyp_time_tmp = phypo->time;
                                                            phypo->hour = (int) (hyp_time_tmp / 3600.0L);
                                                            hyp_time_tmp -= (long double) phypo->hour * 3600.0L;
                                                            phypo->min = (int) (hyp_time_tmp / 60.0L);
                                                            hyp_time_tmp -= (long double) phypo->min * 60.0L;
                                                            phypo->sec = (double) hyp_time_tmp;*/

    return (0);

}

/** function to check for duplicate label and phase (and time) in arrival */

int IsDuplicateArrival(ArrivalDesc *arrival, int num_arrivals, int ntest, int rejectOnlyForExactTimeMatch) {
    int narr;

    for (narr = 0; narr < num_arrivals; narr++) {
        if (narr != ntest
                && !strcmp(arrival[narr].time_grid_label, arrival[ntest].time_grid_label)
                && !strcmp(arrival[narr].phase, arrival[ntest].phase)) {
            if (rejectOnlyForExactTimeMatch) {
                if (fabs(arrival[narr].sec - arrival[ntest].sec) <=
                        ((arrival[narr].error + arrival[ntest].error) / 2.0)) {
                    if (arrival[narr].min == arrival[ntest].min &&
                            arrival[narr].hour == arrival[ntest].hour &&
                            arrival[narr].day == arrival[ntest].day &&
                            arrival[narr].month == arrival[ntest].month &&
                            arrival[narr].year == arrival[ntest].year
                            )
                        return (narr);
                }
            } else {

                return (narr);
            }
        }
    }

    return (-1);

}

/** function to check for duplicate label and phase in arrival */

int IsSameArrival(ArrivalDesc *arrival, int num_arrivals, int ntest, char *phase_test) {
    int narr;

    if (phase_test == NULL) {
        for (narr = 0; narr < num_arrivals; narr++) {
            if (narr != ntest
                    && ((IsPhaseID(arrival[narr].phase, "P") &&
                    IsPhaseID(arrival[ntest].phase, "P"))
                    || (IsPhaseID(arrival[narr].phase, "S") &&
                    IsPhaseID(arrival[ntest].phase, "S")))
                    && !strcmp(arrival[narr].time_grid_label, arrival[ntest].time_grid_label))
                return (narr);
        }
    } else {
        for (narr = 0; narr < num_arrivals; narr++) {
            if (narr != ntest
                    && !strcmp(arrival[narr].time_grid_label, arrival[ntest].time_grid_label)
                    && IsPhaseID(arrival[narr].phase, phase_test))

                return (narr);
        }
    }

    return (-1);

}

/** function to check for duplicate label and phase in arrival */

int FindDuplicateTimeGrid(ArrivalDesc *arrival, int num_arrivals, int ntest) {
    int narr;

    for (narr = 0; narr < num_arrivals; narr++) {
        if (narr != ntest
                && !strcmp(arrival[narr].fileroot, arrival[ntest].fileroot)
                && arrival[narr].flag_ignore == 0
                )

            return (narr);
    }

    return (-1);

}


/** function to read y-z travel time sheet from disk for each arrival */

int ReadArrivalSheets(int num_arrivals, ArrivalDesc *arrival, double xsheet) {

    int istat, narr, ixsheet;
    void **array_tmp;
    double sheet_origx, sheet_dx;


    /* loop over arrivals */

    for (narr = 0; narr < num_arrivals; narr++) {

        /* skip sheet read if arrival has companion */
        if (arrival[narr].n_companion >= 0)
            continue;

        /* skip sheet read or set xsheet to zero for 2D grid */
        if (arrival[narr].gdesc.type == GRID_TIME_2D) {
            if (arrival[narr].sheetdesc.origx < LARGE_DOUBLE)
                continue;
            xsheet = 0.0;
        }

        sheet_origx = arrival[narr].sheetdesc.origx;
        sheet_dx = arrival[narr].sheetdesc.dx;


        /* check which sheets are required from disc */

        /* both required sheets already read */
        if (sheet_origx <= xsheet && xsheet < sheet_origx + sheet_dx)
            continue;

        /* find x index in disk grid of lower plane of dual sheet */
        if (arrival[narr].gdesc.numx > 1)
            ixsheet = (int) ((xsheet - arrival[narr].gdesc.origx)
                / arrival[narr].gdesc.dx);
        else
            ixsheet = 0;
        if (ixsheet < 0 || ixsheet > arrival[narr].gdesc.numx - 1) {
            nll_puterr("WARNING: invalid ixsheet value:");
            snprintf(MsgStr, sizeof(MsgStr), "  Arr: %d  ixsheet: %d", narr, ixsheet);
            nll_puterr(MsgStr);
        }

        /* one required sheet already read */
        if (sheet_origx + sheet_dx <= xsheet &&
                xsheet < sheet_origx + 2.0 * sheet_dx) {
            /* exchange sheet pointers */
            array_tmp = arrival[narr].sheetdesc.array[0];
            arrival[narr].sheetdesc.array[0] =
                    arrival[narr].sheetdesc.array[1];
            arrival[narr].sheetdesc.array[1] = array_tmp;

            /* read next sheet if xsheet not exactly on last sheet */
            /*			if (fabs(xsheet - (sheet_origx + sheet_dx)) */
            /*					> VERY_SMALL_DOUBLE) { */
            /* read new sheet */
            if ((istat =
                    ReadGrid3dBufSheet(
                    arrival[narr].sheetdesc.array[1][0],
                    &(arrival[narr].gdesc),
                    arrival[narr].fpgrid, ixsheet + 1)) < 0)
                nll_puterr(
                    "ERROR: reading new arrival travel time sheet.");
            /*			} */

            /* set dual-sheet origin */
            arrival[narr].sheetdesc.origx += sheet_dx;
        }/* no required sheets already read */
        else {

            /* read lower sheet */
            if ((istat =
                    ReadGrid3dBufSheet(
                    arrival[narr].sheetdesc.array[0][0],
                    &(arrival[narr].gdesc),
                    arrival[narr].fpgrid, ixsheet)) < 0)
                nll_puterr(
                    "ERROR: reading lower arrival travel time sheet.");

            /* read upper sheet if not at last sheet */
            if (ixsheet + 1 < arrival[narr].gdesc.numx) {

                if ((istat =
                        ReadGrid3dBufSheet(
                        arrival[narr].sheetdesc.array[1][0],
                        &(arrival[narr].gdesc),
                        arrival[narr].fpgrid, ixsheet + 1)) < 0)
                    nll_puterr(
                        "ERROR: reading upper arrival travel time sheet.");
            }

            /* set dual-sheet origin */
            arrival[narr].sheetdesc.origx =
                    (double) ixsheet * sheet_dx
                    + arrival[narr].gdesc.origx;
        }

        /*Narr %d O %lf %lf %lf  N %d %d %d  dx %lf %lf %lf\n", narr, arrival[narr].sheetdesc.origx, arrival[narr].sheetdesc.origy, arrival[narr].sheetdesc.origz, arrival[narr].sheetdesc.numx, arrival[narr].sheetdesc.numy, arrival[narr].sheetdesc.numz, arrival[narr].sheetdesc.dx, arrival[narr].sheetdesc.dy, arrival[narr].sheetdesc.dz);*/
    }

    return (0);

}

/** function to construct weight matrix (inverse of covariance matrix) */

int ConstWeightMatrix(int num_arrivals, ArrivalDesc *arrival, GaussLocParams * gauss_par) {

    //printf("DEBUG: ConstWeightMatrix: num_arrivals %d\n", num_arrivals);

    int istat, nrow, ncol;
    double sigmaT2, corr_len2;
    int corr_len_nonzero = 1;
    double dx, dy, dz, dist2;
    double weight_sum;
    double sta_wt, prior_wt;
    SourceDesc *sta1, *sta2;
    double arrivalWeightMax = -1.0;

    double sigmaT, corr_len, dist; // 20150324 AJL - METH_L1_NORM

    // free old matrices
    if (last_matrix_alloc_size > 0) {
        free_matrix_double(edt_matrix, last_matrix_alloc_size, last_matrix_alloc_size);
        free_matrix_double(wt_matrix, last_matrix_alloc_size, last_matrix_alloc_size);
    }
    last_matrix_alloc_size = num_arrivals;
    // allocate square matrices
    edt_matrix = matrix_double(num_arrivals, num_arrivals);
    wt_matrix = matrix_double(num_arrivals, num_arrivals);


    /* set constants */

    sigmaT2 = gauss_par->SigmaT * gauss_par->SigmaT;
    corr_len2 = gauss_par->CorrLen * gauss_par->CorrLen;
    // 20150324 AJL - METH_L1_NORM
    sigmaT = gauss_par->SigmaT;
    corr_len = gauss_par->CorrLen;

    // AJL 20041201 - corr_len_nonzero flag added, before corr_len2 was set to 1.0 (was bug?)
    if (corr_len2 < VERY_SMALL_DOUBLE || gauss_par->CorrLen < 0.0) {
        corr_len_nonzero = 0;
        snprintf(MsgStr, sizeof(MsgStr), "LOCGAU param CorrLen is zero, will not be used: %lf", gauss_par->CorrLen);
        nll_putmsg(2, MsgStr);
    } else {
        corr_len_nonzero = 1;
        snprintf(MsgStr, sizeof(MsgStr), "LOCGAU param CorrLen is non-zero, will be used: %lf", gauss_par->CorrLen);
        nll_putmsg(2, MsgStr);
    }


    /* load covariances */

    for (nrow = 0; nrow < num_arrivals; nrow++) {
        sta1 = &(arrival[nrow].station);
        arrival[nrow].tt_error = gauss_par->SigmaT;
        for (ncol = 0; ncol <= nrow; ncol++) {
            sta2 = &(arrival[ncol].station);

            /* travel time error (TV82, eq. 10-14; MEN92, eq. 22) */
            if (strcmp(arrival[nrow].phase, arrival[ncol].phase) == 0) {
                // same phase types, include spatial correlation
                dx = sta1->x - sta2->x;
                dy = sta1->y - sta2->y;
                dz = sta1->z - sta2->z;
                dist2 = dx * dx + dy * dy + dz * dz;
                if (GeometryMode == MODE_GLOBAL) {
                    dist2 *= DEG2KM * DEG2KM;
                }
                // 20150324 AJL - METH_L1_NORM
                dist = sqrt(dist2);
                // EDT
                if (ncol == nrow) { // diagonal of EDT gets gaussian model time error
                    edt_matrix[nrow][ncol] = sigmaT2;
                } else { // off-diagonal of EDT gets gaussian model weight
                    if (corr_len_nonzero)
                        edt_matrix[nrow][ncol] = edt_matrix[ncol][nrow] = exp(-0.5 * dist2 / corr_len2);
                    else
                        edt_matrix[nrow][ncol] = edt_matrix[ncol][nrow] = 0.0;
                }
                // LS/L2 gaussian model time error
                // AJL 20050914 - bug?  added ncol == nrow case so error is non-zero
                if (ncol == nrow) { // diagonal gets gaussian model time error
                    wt_matrix[nrow][ncol] =
                            (LocMethod == METH_L1_NORM) ?
                            sigmaT // L1  METH_L1_NORM
                            : sigmaT2; // L2  METH_GAU_ANALYTIC
                } else { // off-diagonal
                    if (corr_len_nonzero) {
                        wt_matrix[nrow][ncol] = wt_matrix[ncol][nrow] =
                                (LocMethod == METH_L1_NORM) ?
                                sigmaT * exp(-1.0 * dist / corr_len) // L1  METH_L1_NORM
                                : sigmaT2 * exp(-0.5 * dist2 / corr_len2); // L2  METH_GAU_ANALYTIC
                    } else {
                        wt_matrix[nrow][ncol] = wt_matrix[ncol][nrow] = 0.0;
                    }
                }
            } else {
                // different phase types, assumed no spatial correlation
                edt_matrix[nrow][ncol] = edt_matrix[ncol][nrow] = 0.0;
                wt_matrix[nrow][ncol] = wt_matrix[ncol][nrow] = 0.0;
            }

            /* obs time error */
            if (ncol == nrow) {
                edt_matrix[nrow][ncol] += arrival[nrow].error * arrival[nrow].error;
                wt_matrix[nrow][ncol] +=
                        (LocMethod == METH_L1_NORM) ?
                        arrival[nrow].error // L1  METH_L1_NORM
                        : arrival[nrow].error * arrival[nrow].error; // L2  METH_GAU_ANALYTIC
            }

        }
    }

    if (message_flag >= 5)
        display_matrix_double("Covariance", wt_matrix, num_arrivals, num_arrivals);


    /* invert covariance matrix to obtain weight matrix */

    //if ((istat = nll_dgaussj(wt_matrix, num_arrivals, null_mtrx, 0)) < 0) {
    if ((istat = matrix_double_inverse(wt_matrix, num_arrivals, num_arrivals)) < 0) {
        nll_puterr("ERROR: inverting covariance matrix.");
        return (-1);
    }

    if (message_flag >= 5)
        display_matrix_double("Weight", wt_matrix, num_arrivals, num_arrivals);


    // station distance weighting
    if (iSetStationDistributionWeights) {
        for (nrow = 0; nrow < num_arrivals; nrow++) {
            //printf("station weight: %s %s %s weight: %lf\n", arrival[nrow].label, arrival[nrow].inst, arrival[nrow].comp, arrival[nrow].station_weight);
            for (ncol = 0; ncol <= nrow; ncol++) {
                // 20130627 AJL - change weighing from sum to product
                //sta_wt = (arrival[nrow].station_weight + arrival[ncol].station_weight) / 2.0;
                sta_wt = sqrt(arrival[nrow].station_weight * arrival[ncol].station_weight);
                wt_matrix[nrow][ncol] *= sta_wt;
                if (ncol != nrow) // 20130627 AJL - bug fix
                    wt_matrix[ncol][nrow] *= sta_wt;
            }
        }
    }

    // prior arrival weighting
    // 20130627 AJL - add prior weighting
    if (iUseArrivalPriorWeights) {
        for (nrow = 0; nrow < num_arrivals; nrow++) {
            //printf("station weight: %s %s %s weight: %lf\n", arrival[nrow].label, arrival[nrow].inst, arrival[nrow].comp, arrival[nrow].station_weight);
            for (ncol = 0; ncol <= nrow; ncol++) {
                if (iUseArrivalPriorWeights && arrival[nrow].apriori_weight >= -VERY_SMALL_DOUBLE && arrival[ncol].apriori_weight >= -VERY_SMALL_DOUBLE) {
                    prior_wt = sqrt(arrival[nrow].apriori_weight * arrival[ncol].apriori_weight);
                    wt_matrix[nrow][ncol] *= prior_wt;
                    if (ncol != nrow)
                        wt_matrix[ncol][nrow] *= prior_wt;
                }
            }
        }
    }

    /* get row weights & sum of weights */

    weight_sum = 0.0;
    for (nrow = 0; nrow < num_arrivals; nrow++) {
        arrival[nrow].weight = 0.0;
        for (ncol = 0; ncol < num_arrivals; ncol++) {
            arrival[nrow].weight += wt_matrix[nrow][ncol];
            weight_sum += wt_matrix[nrow][ncol];
            //printf("row %d col %d: wt_tx(r,c) %f  arr(row)_wt %f   wt_sum %lf\n", nrow, ncol, wt_matrix[nrow][ncol], arrival[nrow].weight, weight_sum);
        }
    }
    for (nrow = 0; nrow < num_arrivals; nrow++) {
        arrival[nrow].weight = (double) num_arrivals * arrival[nrow].weight / weight_sum;
        //printf("observation weight: %s %s %s weight: %lf\n", arrival[nrow].label, arrival[nrow].inst, arrival[nrow].comp, arrival[nrow].weight);
        if (arrival[nrow].weight < 0.0) {
            snprintf(MsgStr, sizeof(MsgStr),
                    "ERROR: negative observation weight: %s %s %s weight: %lf",
                    arrival[nrow].label, arrival[nrow].inst,
                    arrival[nrow].comp, arrival[nrow].weight);
            nll_puterr(MsgStr);
            nll_puterr("   Gaussian model error (see LOCGAU) may be too large relative to obs uncertainty (see LOCQUAL2ERR, or NLL-Phase format ErrMag).");
        }
        if (arrival[nrow].weight > arrivalWeightMax)
            arrivalWeightMax = arrival[nrow].weight;
    }
    if (message_flag >= 4) {
        snprintf(MsgStr, sizeof(MsgStr), "Weight Matrix sum: %lf", weight_sum);
        nll_putmsg(4, MsgStr);
    }


    // set global variables
    gauss_par->EDTMtrx = edt_matrix;
    gauss_par->WtMtrx = wt_matrix;
    gauss_par->WtMtrxSum = weight_sum;
    gauss_par->arrivalWeightMax = arrivalWeightMax;

    return (0);

}

/** function to do weight matrix memory cleanup */

int CleanWeightMatrix() {

    // AJL - 20080710 (valgrind)
    // free EDT_OT_WT memory
    if (edt_matrix != NULL)
        free_matrix_double(edt_matrix, last_matrix_alloc_size, last_matrix_alloc_size);
    edt_matrix = NULL;
    if (wt_matrix != NULL)
        free_matrix_double(wt_matrix, last_matrix_alloc_size, last_matrix_alloc_size);
    wt_matrix = NULL;
    last_matrix_alloc_size = -1;

    return (0);

}



/** function to calculate weighted mean of observed arrival times */

/*		(TV82, eq. A-38) */

void CalcCenteredTimesObs(int num_arrivals, ArrivalDesc *arrival,
        GaussLocParams* gauss_par, HypoDesc * phypo) {

    int nrow, ncol, narr;
    long double sum, weighted_mean;
    MatrixDouble wtmtx;
    double *wtmtxrow, wt_sum;


    if (!FixOriginTimeFlag) {

        /* calculate weighted mean of observed times */

        wtmtx = gauss_par->WtMtrx;
        sum = 0.0L;
        wt_sum = 0.0;

        for (nrow = 0; nrow < num_arrivals; nrow++) {
            if (!arrival[nrow].abs_time)
                continue; // ignore obs without absolute timing
            wtmtxrow = wtmtx[nrow];
            for (ncol = 0; ncol < num_arrivals; ncol++) {
                if (!arrival[ncol].abs_time)
                    continue; // ignore obs without absolute timing
                sum += (long double) *(wtmtxrow + ncol) * arrival[ncol].obs_time;
                wt_sum += (double) *(wtmtxrow + ncol);
            }
        }
        if (wt_sum > 0.0)
            weighted_mean = sum / (long double) wt_sum;
        else
            weighted_mean = (long double) arrival[0].obs_time;

    } else {

        /* use fixed origin time as reference */

        weighted_mean = phypo->time;
    }


    /* set centered observed times */

    if (message_flag >= 3) {
        nll_putmsg(3, "");
        nll_putmsg(3, "Delayed, Sorted, Centered Observations:");
    }
    for (narr = 0; narr < num_arrivals; narr++) {
        arrival[narr].obs_centered =
                (double) (arrival[narr].obs_time - weighted_mean);
        if (message_flag >= 3) {

            snprintf(MsgStr, sizeof(MsgStr),
                    "  %3d  %-12s %-6s %2.2d:%2.2d:%7.4lf - %7.4lfs -> %8.4lf (%10.4lf)",
                    narr, arrival[narr].label, arrival[narr].phase,
                    arrival[narr].hour, arrival[narr].min,
                    arrival[narr].sec, arrival[narr].delay, arrival[narr].obs_centered,
                    ((double) arrival[narr].obs_time));
            nll_putmsg(3, MsgStr);
        }
    }

    gauss_par->meanObs = weighted_mean;

}


/** function to calculate weighted mean of predicted travel times */

/*		(TV82, eq. A-38) */

void CalcCenteredTimesPred(int num_arrivals, ArrivalDesc *arrival, GaussLocParams * gauss_par) {

    int nrow, ncol, narr;
    double sum, weighted_mean, pred_time_row;
    MatrixDouble wtmtx;
    double *wtmtxrow, wt_sum;


    if (!FixOriginTimeFlag) {

        wtmtx = gauss_par->WtMtrx;
        sum = 0.0;
        wt_sum = 0.0;

        for (nrow = 0; nrow < num_arrivals; nrow++) {
            // AJL 20041115 bug fix!
            if (arrival[nrow].pred_travel_time <= 0.0)
                continue; // ignore obs without predicted times
            // END
            if (!arrival[nrow].abs_time)
                continue; // ignore obs without absolute timing
            wtmtxrow = wtmtx[nrow];
            pred_time_row = arrival[nrow].pred_travel_time;
            for (ncol = 0; ncol < num_arrivals; ncol++) {
                // AJL 20041115 bug fix!
                if (arrival[ncol].pred_travel_time <= 0.0)
                    continue; // ignore obs without predicted times
                // END
                if (!arrival[ncol].abs_time)
                    continue; // ignore obs without absolute timing
                sum += (double) *(wtmtxrow + ncol) * pred_time_row;
                wt_sum += (double) *(wtmtxrow + ncol);
            }
        }

        if (wt_sum > 0.0)
            weighted_mean = sum / wt_sum;
        else
            weighted_mean = (long double) arrival[0].pred_travel_time;

    } else {

        // for fixed origin time use travel time directly
        weighted_mean = 0.0;

    }



    /* set centered predicted times */

    for (narr = 0; narr < num_arrivals; narr++) {
        // AJL 20041115 bug fix!
        if (arrival[narr].pred_travel_time <= 0.0)

            continue; // ignore obs without predicted times
        // END
        arrival[narr].pred_centered = arrival[narr].pred_travel_time - weighted_mean;
    }


    gauss_par->meanPred = (double) weighted_mean;

}

//static double maxvalue = -1.0;






/** function to open summary output files */

int OpenSummaryFiles(char *path_output, char* loctypename) {

    int ngrid;
    char fname[FILENAME_MAX];



    for (ngrid = 0; ngrid < NumLocGrids; ngrid++) {

        if (!LocGridSave[ngrid])
            continue;

        /* Grid Hyp format */

        pSumFileHypNLLoc[ngrid] = NULL;
        snprintf(fname, sizeof(fname), "%s.sum.%s%d.loc.hyp", path_output, loctypename, ngrid);
        if ((pSumFileHypNLLoc[ngrid] = fopen(fname, "w")) == NULL) {
            nll_puterr2("ERROR: opening summary output file", fname);
            return (-1);
        } else {
            NumFilesOpen++;
        }


        /* CSV Hyp format */

        pSumFileHypNLLocCSV[ngrid] = NULL;
        snprintf(fname, sizeof(fname), "%s.sum.%s%d.loc.csv", path_output, loctypename, ngrid);
        if ((pSumFileHypNLLocCSV[ngrid] = fopen(fname, "w")) == NULL) {
            nll_puterr2("ERROR: opening CSV summary output file", fname);
            return (-1);
        } else {
            NumFilesOpen++;
        }
        WriteLocationCSVheader(pSumFileHypNLLocCSV[ngrid]);

        iWriteHypHeader[ngrid] = 1;

        /* Hypo71 format */
        pSumFileHypo71[ngrid] = NULL;
        if (iSaveHypo71Sum) {
            snprintf(fname, sizeof(fname), "%s.sum.%s%d.loc.hypo_71", path_output, loctypename, ngrid);
            if ((pSumFileHypo71[ngrid] = fopen(fname, "w"))
                    == NULL) {
                nll_puterr2(
                        "ERROR: opening HYPO71 summary output file",
                        fname);
                return (-1);
            } else {
                NumFilesOpen++;
            }
            fprintf(pSumFileHypo71[ngrid], "%s\n",
                    Hypocenter.comment);
        }


        /* HypoEllipse format */
        pSumFileHypoEll[ngrid] = NULL;
        if (iSaveHypoEllSum) {
            snprintf(fname, sizeof(fname), "%s.sum.%s%d.loc.hypo_ell", path_output, loctypename, ngrid);
            if ((pSumFileHypoEll[ngrid] = fopen(fname, "w"))
                    == NULL) {
                nll_puterr2(
                        "ERROR: opening HypoEllipse summary output file",
                        fname);
                return (-1);
            } else {
                NumFilesOpen++;
            }
            fprintf(pSumFileHypoEll[ngrid], "%s\n",
                    Hypocenter.comment);
        }


        /* HypoInverse Archive format */
        pSumFileHypoInv[ngrid] = NULL;
        if (iSaveHypoInvSum) {
            snprintf(fname, sizeof(fname), "%s.sum.%s%d.loc.hypo_inv", path_output, loctypename, ngrid);
            if ((pSumFileHypoInv[ngrid] = fopen(fname, "w"))
                    == NULL) {
                nll_puterr2(
                        "ERROR: opening HypoInverse Archive summary output file",
                        fname);
                return (-1);
            } else {
                NumFilesOpen++;
            }
        }

        /* HypoInverse Archive Y2000 format */
        pSumFileHypoInvY2K[ngrid] = NULL;
        if (iSaveHypoInvY2KArc) {
            snprintf(fname, sizeof(fname), "%s.sum.%s%d.loc.arc", path_output, loctypename, ngrid);
            if ((pSumFileHypoInvY2K[ngrid] = fopen(fname, "w"))
                    == NULL) {
                nll_puterr2(
                        "ERROR: opening HypoInverse Archive Y2000 summary output file",
                        fname);
                return (-1);
            } else {
                NumFilesOpen++;
            }
        }

        /* Alberto 3D 4 chr sta SIMULPS format */
        pSumFileAlberto4[ngrid] = NULL;
        if (iSaveAlberto4Sum) {
            snprintf(fname, sizeof(fname), "%s.sum.%s%d.loc.sim", path_output, loctypename, ngrid);
            if ((pSumFileAlberto4[ngrid] = fopen(fname, "w"))
                    == NULL) {
                nll_puterr2(
                        "ERROR: opening Alberto 3D, 4 chr sta, SIMULPS output file",
                        fname);
                return (-1);
            } else {

                NumFilesOpen++;
            }
        }

        // fmamp hypocenter-phase format  // 20160920 AJL - added
        pSumFileFmamp[ngrid] = NULL;
        if (iSaveFmamp) {
            snprintf(fname, sizeof(fname), "%s.sum.%s%d.loc.fmamp", path_output, loctypename, ngrid);
            if ((pSumFileFmamp[ngrid] = fopen(fname, "w"))
                    == NULL) {
                nll_puterr2(
                        "ERROR: opening Fmamp output file",
                        fname);
                return (-1);
            } else {

                NumFilesOpen++;
            }
        }

    }


    return (0);

}

/** function to close summary output files */

int CloseSummaryFiles() {
    int ngrid;


    for (ngrid = 0; ngrid < NumLocGrids; ngrid++) {

        if (!LocGridSave[ngrid])
            continue;

        /* Grid Hyp format */

        if (pSumFileHypNLLoc[ngrid] != NULL) {
            fclose(pSumFileHypNLLoc[ngrid]);
            pSumFileHypNLLoc[ngrid] = NULL;
            NumFilesOpen--;
        }

        /* CSV Hyp format */

        if (pSumFileHypNLLocCSV[ngrid] != NULL) {
            fclose(pSumFileHypNLLocCSV[ngrid]);
            pSumFileHypNLLocCSV[ngrid] = NULL;
            NumFilesOpen--;
        }

        /* Hypo71 format */

        if (pSumFileHypo71[ngrid] != NULL) {
            fclose(pSumFileHypo71[ngrid]);
            NumFilesOpen--;
        }

        /* HypoEll format */

        if (pSumFileHypoEll[ngrid] != NULL) {
            fclose(pSumFileHypoEll[ngrid]);
            NumFilesOpen--;
        }

        /* HypoInv format */

        if (pSumFileHypoInv[ngrid] != NULL) {
            fclose(pSumFileHypoInv[ngrid]);
            NumFilesOpen--;
        }

        /* HypoInv Y2000 format */

        if (pSumFileHypoInvY2K[ngrid] != NULL) {
            fclose(pSumFileHypoInvY2K[ngrid]);
            NumFilesOpen--;
        }

        /* Alberto 3D 4 chr sta SIMULPS format */
        if (pSumFileAlberto4[ngrid] != NULL) {
            fclose(pSumFileAlberto4[ngrid]);
            NumFilesOpen--;
        }

        // fmamp hypocenter-phase format */
        if (pSumFileFmamp[ngrid] != NULL) {

            fclose(pSumFileFmamp[ngrid]);
            NumFilesOpen--;
        }


    }

    return (0);

}









/** function to set station distribution weights */

int setStationDistributionWeights(SourceDesc *stations, int numStations, ArrivalDesc *arrival, int nArrivals) {

    int i, nsta, m, nError = 0;
    double x, y, dist, station_weight, weight;
    double dist_ave, cutoff2;
    double station_weight_sum;
    ArrivalDesc *arr;


    // get cutoff distance
    if (stationDistributionWeightCutoff > 0.0) {
        cutoff2 = stationDistributionWeightCutoff * stationDistributionWeightCutoff;
    } else {
        // use average distance
        dist_ave = calcAveInterStationDistance(stations, numStations);
        if (message_flag >= 2) {
            snprintf(MsgStr, sizeof(MsgStr), "Station Dist Weight:  Ave Station Distance: %lf", dist_ave);
            nll_putmsg(2, MsgStr);
        }
        if (dist_ave <= 0.0)
            return (-1);
        cutoff2 = dist_ave * dist_ave;
    }



    // calculate station weights

    station_weight_sum = 0.0;
    nsta = 0;
    for (i = 0; i < nArrivals; i++) {
        arr = arrival + i;
        station_weight = 0.0;
        x = arr->station.x;
        y = arr->station.y;
        if (x == 0.0 && y == 0.0) // station location not known
            continue;
        for (m = 0; m < numStations; m++) {
            if ((stations + m)->ignored)
                continue;
            // 20240207 AJL - Bug fix: skip stations without coordinates
            if ((stations + m)->x < -1.0e6)
                continue;
            dist = GetEpiDist((stations + m), x, y);
            weight = exp(-(dist * dist) / cutoff2); // 0.0 for farthest -> 1.0 for same position
            //printf("other dist/weigh: %s d=%lf w=%lf (%lf,%lf,%lf)\n",
            //        (stations + m)->label, dist, weight, (stations + m)->x, (stations + m)->y, (stations + m)->z);
            station_weight += weight; // count of number close stations
        }
        nsta++;
        station_weight = 1.0 / station_weight; // weight inversely prop to number of close stations
        arr->station_weight = station_weight;
        //printf("Station Weight: %s %lf (%lf,%lf,%lf) cutoff2 %lf\n", arr->label, arr->station_weight, arr->station.x, arr->station.y, arr->station.z, sqrt(cutoff2));
        station_weight_sum += station_weight;
    }
    if (nsta > 0) {
        // normalize
        station_weight_sum /= (double) nsta;
        for (i = 0; i < nArrivals; i++) {
            arr = arrival + i;
            arr->station_weight /= station_weight_sum;
            if (message_flag >= 2) {

                snprintf(MsgStr, sizeof(MsgStr), "Station Dist Weight: %s %lf (%lf,%lf,%lf)",
                        arr->label, arr->station_weight, arr->station.x, arr->station.y, arr->station.z);
                nll_putmsg(2, MsgStr);
            }
        }
    }


    return (nError);


}

/** function to get travel times for all observed arrivals */

int getTravelTimes(ArrivalDesc *arrival, int num_arr_loc, double xval, double yval, double zval) {

    int nReject;
    int narr, n_compan;
    FILE* fp_grid;
    double yval_grid = 0.0;
    GridDesc* ptgrid;

    // 20101005 AJL - added calculation of mean slowness
    double slowness_P = -1.0;
    double slowness_S = -1.0;
    if (LocMethod == METH_OT_STACK) {
        if (fp_model_grid_P != NULL) {
            if (model_grid_P.numx > 2) {
                // 3D grid
                slowness_P = (double) ReadAbsInterpGrid3d(fp_model_grid_P, &model_grid_P, xval, yval, zval, 0);
            } else {
                // 2D grid (1D model)
                yval_grid = model_grid_P.dy; // aribitrary, small y grid value
                slowness_P = ReadAbsInterpGrid2d(fp_model_grid_P, &model_grid_P, yval_grid, zval);
            }
            if (GeometryMode != MODE_GLOBAL)
                slowness_P /= model_grid_P.dy; // value in model file is slowness * ds
        }
        if (fp_model_grid_S != NULL) {
            if (model_grid_S.numx > 2) {
                // 3D grid
                slowness_S = (double) ReadAbsInterpGrid3d(fp_model_grid_S, &model_grid_S, xval, yval, zval, 0);
            } else {
                // 2D grid (1D model)
                yval_grid = model_grid_S.dy; // aribitrary, small y grid value
                slowness_S = ReadAbsInterpGrid2d(fp_model_grid_S, &model_grid_S, yval_grid, zval);
            }
            if (GeometryMode != MODE_GLOBAL)
                slowness_S /= model_grid_S.dy; // value in model file is slowness * ds
        }
        if (slowness_P <= SMALL_FLOAT)
            slowness_P = -1.0;
        if (slowness_S < 0.0 && VpVsRatio > 0.0)
            slowness_S = slowness_P * VpVsRatio;
        if (slowness_S <= SMALL_FLOAT)
            slowness_S = -1.0;
        /*static int icount = 0;
        if (icount++ % 100 == 0) {
            //printf("depth=%f  slowness_P=%f  slowness_S=%f\n", zval, slowness_P, slowness_S);
            printf("depth=%f  1/slowness_P=%f  1/slowness_S=%f\n", zval, 1.0 / slowness_P, 1.0 / slowness_S);
            //printf("model_grid_P.numx=%d  model_grid_P.dy=%f  yval_grid=%f\n", model_grid_P.numx, model_grid_P.dy, yval_grid);
        }*/
    }

    /* loop over observed arrivals */

    nReject = 0;
    for (narr = 0; narr < num_arr_loc; narr++) {
        /* check for companion */
        if ((n_compan = arrival[narr].n_companion) >= 0) {
            if ((arrival[narr].pred_travel_time = arrival[n_compan].pred_travel_time) < 0.0)
                nReject++;
            arrival[narr].pred_travel_time *= arrival[narr].tfact;
            /* else check grid type */
        } else {
            if (arrival[narr].gdesc.type == GRID_TIME) {
                /* 3D grid */
                if (arrival[narr].gdesc.buffer == NULL) {
                    /* read time grid from disk */
                    fp_grid = arrival[narr].fpgrid;
                } else {
                    /* read time grid from memory buffer */
                    fp_grid = NULL;
                }
                if ((arrival[narr].pred_travel_time = (double) ReadAbsInterpGrid3d(fp_grid, &(arrival[narr].gdesc),
                        xval, yval, zval, 0)) < 0.0)
                    nReject++;
            } else {
                /* 2D grid (1D model) */
                yval_grid = GetEpiDist(&(arrival[narr].station), xval, yval);
                if (GeometryMode == MODE_GLOBAL)
                    yval_grid *= KM2DEG;
                if (arrival[narr].sheetdesc.buffer == NULL) {
                    /* read time grid from disk */
                    fp_grid = arrival[narr].fpgrid;
                    ptgrid = &(arrival[narr].gdesc);
                } else {
                    /* read time grid from memory buffer */
                    fp_grid = NULL;
                    ptgrid = &(arrival[narr].sheetdesc);
                }
                if ((arrival[narr].pred_travel_time = ReadAbsInterpGrid2d(fp_grid, ptgrid, yval_grid, zval)) < 0.0)
                    nReject++;
                //printf("DEBUG: getTT:  xval %lf yval %lf yval_grid %lf zval %lf t %lf \n", xval, yval, yval_grid, zval, arrival[narr].pred_travel_time);
                //display_grid_param(&(arrival[narr].sheetdesc));
            }
            arrival[narr].pred_travel_time *= arrival[narr].tfact;
            // apply crustal correction
            if (ApplyCrustElevCorrFlag && GeometryMode == MODE_GLOBAL
                    && arrival[narr].pred_travel_time > 0.0) {
                //printf("arrival[narr].pred_travel_time before: %f", arrival[narr].pred_travel_time);
                if (yval_grid > MinDistCrustElevCorr)
                    arrival[narr].pred_travel_time += applyCrustElevCorrection(arrival + narr, xval, yval, zval);
                //printf(" -> after: %f\n", arrival[narr].pred_travel_time);
            } else if (ApplyElevCorrFlag) {

                if (arrival[narr].pred_travel_time > 0.0) // ignore arrivals with no pred tt
                    arrival[narr].pred_travel_time += arrival[narr].elev_corr;
            }
        }

        // set slowness
        if (arrival[narr].isS)
            arrival[narr].slowness = slowness_S;

        else
            arrival[narr].slowness = slowness_P;

    }


    return (nReject);

}



/** function to apply crustal correction and elevation correction */

// assumes vertical ray (dtdd = 0.0) !!!

double applyCrustElevCorrection(ArrivalDesc* parrival, double xval, double yval, double zval) {

    double correction = 0.0;
    double dtdd = 0.0;
    char cphase;

    if (IsPhaseID(parrival->phase, "P"))
        cphase = 'P';
    else if (IsPhaseID(parrival->phase, "S"))
        cphase = 'S';
    else
        return (0.0);


    // source
    correction = calc_crust_corr(cphase, yval, xval, zval, VERY_LARGE_DOUBLE, dtdd);
    // receiver
    correction +=
            calc_crust_corr(cphase, parrival->station.dlat,
            parrival->station.dlong, 0.0, -1000.0 * parrival->station.depth, dtdd);

    return (correction);

}

/** function to check if xyz location is above topography */

int isAboveTopo(double xval, double yval, double zval) {

    double xlon, ylat, elev;
    double topo_elev;

    // check topo if available
    if (topo_surface_index >= 0) {
        if (model_surface[topo_surface_index].is_latlon) {
            // convert xyz to lat/long/elev
            rect2latlon(0, xval, yval, &ylat, &xlon);
        } else {
            xlon = xval;
            ylat = yval;
        }
        if (map_itype[0] != MAP_TRANS_NONE) // 20110105 AJL
            elev = -zval * 1000.0; // elevation in meters
        else
            elev = -zval; // 20110105 AJL
        // get elevation of topo at location
        topo_elev = get_surface_z(topo_surface_index, xlon, ylat);
        int iabove = elev > topo_elev;

        /*
            if (iabove) {
                printf("xyz %f %f %f  lon/lat/elev %f %f %f  topo_elev %f  above %d\n", xval, yval, zval, ylat, xlon, elev, topo_elev, iabove);
            }//*/

        return (iabove);
    }

    return (0);

}





