/*
 * Copyright (C) 1999-2020 Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
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


/*   Loc2ssst.c

        Program to generate SSST station correction grids from a set of NLL locations.

 */



/*
        history:

        ver 01    21Jan2020  AJL  Original version



.........1.........2.........3.........4.........5.........6.........7.........8

 */



#include "GridLib.h"
#include "phaseloclist.h"


/* defines */

//#define MAX_NUM_INPUT_FILES 50000
/* MAX_NUM_INPUT_FILES values of > 215 give a segemtation violation when jumping into the subroutine
    (gcc (GCC) 3.2 20020903 (Red Hat Linux 8.0 3.2-7) Frederik Tilmann */
// reduce array sizes to see if this is cause of segmentation error on some systems.
//#define MAX_NUM_INPUT_FILES 50000
/* MAX_NUM_INPUT_FILES values of > 12254 give a segmentation violation when jumping into the subroutine Sum_Locations
      (gcc (GCC) 3.2 20020903 (Red Hat Linux 8.0 3.2-7)
also: (gcc (GCC) 4.4.4 20100630 (Red Hat 4.4.4-10)   Kernel: 2.6.32.26-175.fc12.x86_64
 // reduce array sizes to see if this is cause of segmentation error on some systems.
 Frederik Tilmann */
//#define MAX_NUM_INPUT_FILES 32000
// 20170210 AJL - reset to 50k
#define MAX_NUM_INPUT_FILES 50000


/* globals */

// station list
int NumStationPhases;
SourceDesc StationPhaseList[X_MAX_NUM_ARRIVALS];
// station/phase ssst grids
GridDesc ssst_grid_template;
GridDesc ssst_time_grid_template;
// LocNode array
int NumLocNodes;
LocNode *LocNodeArray[MAX_NUM_INPUT_FILES];
// PhsNode array
int NumPhsNodes;
PhsNode *PhsNodeArray[X_MAX_NUM_ARRIVALS];

typedef struct {
    double char_dist; // Characteristic event-station distance in km for weighting contribution of an event to SSST correction for a station
    double weight_floor; // (0.0-1.0); small value added to events-station weights so ssst values at large event-station distance remain non-zero (station static).
    char stations[MAXLINE]; // List of stations to process, separated by non-character, no spaces
    int use_rejected; // Flag to indicate that NLL REJECTED locations should be accepted for SSST processing (default=0)
}
LS_Params;
LS_Params Params;

typedef struct {
    double RMSMax;
    int NRdgsMin;
    double GapMax, PResidualMax, SResidualMax, EllLen3Max;
}
LS_Phstat;
LS_Phstat PhsStat;

int angle_mode = ANGLE_MODE_NO; /* angle mode - ANGLE_MODE_NO, ANGLE_MODE_YES */

char fn_hypos_in[FILENAME_MAX];
char fn_ls_output[FILENAME_MAX];
char fn_time_input[FILENAME_MAX] = "";
int ihave_time_input_grids = 0;
double VpVsRatio;
int iSwapBytesOnInput;


/* functions */

int GenSSST(int argc, char** argv);
int ReadLoc2ssstInput(FILE*);
int DoLoc2ssst();
double get_ssst_value(double xval, double yval, double zval,
        PhsNode **phs_node_array, int num_phs_nodes, LocNode **loc_node_array, int num_loc_nodes, LS_Params *pparams);
int open_traveltime_grid(ArrivalDesc* parr, char *fn_time_grid_input, char *stacode, char *phasecode, double vp_vs_ratio, double *ptfact);
int add_ssst_to_traveltime_grid(char *phasecode, char *stacode, GridDesc *pssst_grid, GridDesc *ptraveltime_grid, GridDesc *pssst_time_grid, SourceDesc* psrce, double tfact);
int GenAngleGrid(GridDesc* ptgrid, SourceDesc* psource, char *filename, GridDesc* pagrid, int angle_mode);


/*** program to sum event scatter files */

#define PNAME  "Loc2ssst"
#define NARGS 2

int main(int argc, char** argv) {

    int istat;

    /* set program name */

    // 20210511 AJL - Bug fix: changed to strncpy from strlcpy which is not available in Linux.
    strncpy(prog_name, PNAME, sizeof (prog_name));


    // check command line for correct usage

    if (argc != NARGS) {
        disp_usage(prog_name, "<control file>");
        exit(EXIT_ERROR_USAGE);
    }

    SetConstants();

    // initializations before reading input
    Params.char_dist = LARGE_FLOAT;
    Params.weight_floor = 0.0;
    strcpy(Params.stations, "");
    Params.use_rejected = 0;

    PhsStat.RMSMax = 1.0e6;
    PhsStat.NRdgsMin = 0;
    PhsStat.GapMax = 360.0;
    PhsStat.PResidualMax = 1.0e6;
    PhsStat.SResidualMax = 1.0e6;
    PhsStat.EllLen3Max = 1.0e6;

    // read control file
    // 20210511 AJL - Bug fix: changed to strncpy from strlcpy which is not available in Linux.
    strncpy(fn_control, argv[1], sizeof (fn_control));
    if ((fp_control = fopen(fn_control, "r")) == NULL) {
        nll_puterr("ERROR: opening control file.");
        exit(EXIT_ERROR_FILEIO);
    }

    if ((istat = ReadLoc2ssstInput(fp_control)) < 0) {
        nll_puterr("ERROR: reading Loc2ssst control file.");
        exit(EXIT_ERROR_FILEIO);
    }

    if ((istat = DoLoc2ssst()) < 0) {
        nll_puterr("ERROR: doing Loc2ssst process.");
        exit(-1);
    }



    exit(0);

}

/*** function to read input file name ***/

int get_ls_inpfile(char* line1) {

    sscanf(line1, "%s", fn_hypos_in);

    snprintf(MsgStr, sizeof (MsgStr), "LSLOCFILES:  Input: %s", fn_hypos_in);
    nll_putmsg(1, MsgStr);

    return (0);
}

/*** function to read output file name ***/

int get_ls_outfile(char* line1) {

    sscanf(line1, "%s", fn_ls_output);

    snprintf(MsgStr, sizeof (MsgStr), "LSOUT:  Output: %s.*",
            fn_ls_output);
    nll_putmsg(1, MsgStr);

    return (0);
}

/*** function to read hypocenter filters ***/

int get_ls_phstat(char* line1) {

    int istat = sscanf(line1, "%lf %d %lf %lf %lf %lf",
            &PhsStat.RMSMax, &PhsStat.NRdgsMin, &PhsStat.GapMax, &PhsStat.PResidualMax, &PhsStat.SResidualMax, &PhsStat.EllLen3Max);

    sprintf(MsgStr,
            "LSPHSTAT:  RMSMax: %f  NRdgsMin: %d  GapMax: %.3g  PResidualMax: %.3g SResidualMax: %.3g EllLen3Max %.3g",
            PhsStat.RMSMax, PhsStat.NRdgsMin, PhsStat.GapMax, PhsStat.PResidualMax, PhsStat.SResidualMax, PhsStat.EllLen3Max);
    nll_putmsg(1, MsgStr);

    if (istat < 6)

        return (-1);

    return (0);
}

/*** function to read ls parameters ***/

int get_ls_params(char* line1) {

    int istat = sscanf(line1, "%lf %lf %d",
            &Params.char_dist, &Params.weight_floor, &Params.use_rejected);

    if (Params.weight_floor < 0.0) {
        Params.weight_floor = 0.0;
        snprintf(MsgStr, sizeof (MsgStr), "WARNING: WeightFloor < 0.0, reset to 0.0.");
        nll_puterr(MsgStr);
    }

    sprintf(MsgStr,
            "LSPARAMS:  CharDist: %f  WeightFloor:%f  UseRejected:%d",
            Params.char_dist, Params.weight_floor, Params.use_rejected);
    nll_putmsg(1, MsgStr);

    if (istat < 1)

        return (-1);

    return (0);
}

/*** function to read ls stations ***/

int get_ls_stations(char* line1) {

    int istat = sscanf(line1, "%s", Params.stations);

    sprintf(MsgStr, "LSSTATIONS:  %s", Params.stations);
    nll_putmsg(1, MsgStr);

    if (istat < 1)

        return (-1);

    return (0);
}

/*** function to read grid mode params ***/

int get_ls_mode(char* line1) {

    char str_angle_mode[MAXLINE];


    sscanf(line1, "%s", str_angle_mode);

    sprintf(MsgStr, "LSMODE:  %s", str_angle_mode);
    nll_putmsg(1, MsgStr);

    if (strcmp(str_angle_mode, "ANGLES_YES") == 0)
        angle_mode = ANGLE_MODE_YES;
    else if (strcmp(str_angle_mode, "ANGLES_NO") == 0)
        angle_mode = ANGLE_MODE_NO;
    else if (strcmp(str_angle_mode, "ANGLES_INCLINATION") == 0)
        angle_mode = ANGLE_MODE_INCLINATION;
    else {
        angle_mode = ANGLE_MODE_UNDEF;
        nll_puterr("ERROR: unrecognized angle mode");
        return (-1);
    }

    return (0);

}

/** function to read output file name
 *
 * NOTE: Do not change format, see NLlocLib.c->GetNLLoc_Method()
 *
 */


int GetNLLoc_Files(char* line1) {

    int istat;

    //NLLocLib.c  istat = sscanf(line1, "%s %s %s %s %d", fnobs, ftype_obs, fn_loc_grids, fn_path_output, &iSwapBytesOnInput);
    istat = sscanf(line1, "%*s %*s %s %*s %d", fn_time_input, &iSwapBytesOnInput);
    if (istat < 2)
        iSwapBytesOnInput = 0;

    sprintf(MsgStr, "LOCFILES:  InputTimeGrids: %s.* iSwapBytesOnInput: %d", fn_time_input, iSwapBytesOnInput);
    nll_putmsg(1, MsgStr);

    return (0);
}

/** function to read method
 *
 * NOTE: Do not change format, see NLlocLib.c->GetNLLoc_Method()
 *
 */

int GetNLLoc_Method(char* line1) {

    int istat;

    //     istat = sscanf(line1, "%s %lf %d %d %d %lf %d %lf %d", loc_method,
    //        &DistStaGridMax, &MinNumArrLoc, &MaxNumArrLoc, &MinNumSArrLoc,
    //        &VpVsRatio, &MaxNum3DGridMemory, &DistStaGridMin, &iRejectDuplicateArrivals);
    istat = sscanf(line1, "%*s %*lf %*d %*d %*d %lf", &VpVsRatio);

    sprintf(MsgStr, "LOCMETHOD:  VpVsRatio: %lf", VpVsRatio);
    nll_putmsg(1, MsgStr);

    // 20200203 AJL - not sure if this is correct, may be OK?  TODO:
    /*if (VpVsRatio > 0.0 && GeometryMode == MODE_GLOBAL) {
        nll_puterr("ERROR: cannot use VpVsRatio>0 with TRANSFORM GLOBAL.");
        return (-1);
    }*/

    return (0);
}

/*** function to read input file */

int ReadLoc2ssstInput(FILE * fp_input) {

    int istat, iscan;
    char param[MAXLINE], *pchr;
    char line[2 * MAXLINE], *fgets_return;

    int flag_control = 0, flag_inpfile = 0, flag_outfile = 0, flag_grid = 0,
            flag_trans = 0, flag_params = 0;
    int flag_include = 1;
    int flag_method = 1;
    int flag_phase_id = 1;
    int flag_phstat = 1;
    int flag_mode = 1;
    int flag_stations = 1;

    int flag_source = 0;
    int flag_nlloc_outfile = 0;
    int flag_out_grid = 0;

    // read each input line

    while ((fgets_return = fgets(line, 2 * MAXLINE, fp_input)) != NULL
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


        // read transform params

        if (strcmp(param, "TRANS") == 0) {
            if ((istat = get_transform(0, strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading TRANS parameters.");
            else
                flag_trans = 1;
        }


        /* read file names */

        if (strcmp(param, "LOCFILES") == 0) {
            if ((istat = GetNLLoc_Files(strchr(line, ' '))) < 0) {
                nll_puterr("ERROR: reading NLLoc output file name.");
            } else {
                flag_nlloc_outfile = 1;
            }
        }


        /* read method */

        if (strcmp(param, "LOCMETH") == 0) {
            if ((istat = GetNLLoc_Method(strchr(line, ' '))) < 0) {
                nll_puterr("ERROR: reading NLLoc method.");
                flag_method = 0;
            }
        }


        /* read phase identifier values */

        if (strcmp(param, "LOCPHASEID") == 0) {
            if ((istat = GetPhaseID(strchr(line, ' '))) < 0) {
                nll_puterr("ERROR: reading phase identifier values.");
                flag_phase_id = 0;
            }
        }


        /* read source params */

        if (strcmp(param, "LOCSRCE") == 0 || strcmp(param, "GTSRCE") == 0) {
            if ((istat = GetNextSource(strchr(line, ' '))) < 0) {
                nll_puterr("ERROR: reading source params:");
                nll_puterr(line);
            } else
                flag_source = 1;
        }



        // read input file name

        if (strcmp(param, "LSLOCFILES") == 0) {
            if ((istat = get_ls_inpfile(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading LSLOCFILES parameters.");
            else
                flag_inpfile = 1;
        }


        // read output file name

        if (strcmp(param, "LSOUT") == 0) {
            if ((istat = get_ls_outfile(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading LSOUT parameters.");
            else
                flag_outfile = 1;
        }


        // read grid params

        if (strcmp(param, "LSGRID") == 0) {
            if ((istat = get_grid(strchr(line, ' '))) < 0) {
                nll_puterr("ERROR: reading LSGRID parameters.");
            } else {
                flag_grid = 1;
                ssst_grid_template = grid_in;
            }
        }


        // read grid params

        if (strcmp(param, "LSOUTGRID") == 0) {
            if ((istat = get_grid(strchr(line, ' '))) < 0) {
                nll_puterr("ERROR: reading LSOUTGRID parameters.");
            } else {
                flag_out_grid = 1;
                ssst_time_grid_template = grid_in;
            }
        }


        // read parameters

        if (strcmp(param, "LSPARAMS") == 0) {
            if ((istat = get_ls_params(strchr(line, ' '))) < 0) {
                nll_puterr("ERROR: reading LSPARAMS parameters.");
            } else {
                flag_params = 1;
            }
        }


        // read filters

        if (strcmp(param, "LSPHSTAT") == 0) {
            if ((istat = get_ls_phstat(strchr(line, ' '))) < 0) {
                nll_puterr("ERROR: reading LSPHSTAT parameters.");
                flag_phstat = 0;
            }
        }


        /* read stations */

        if (strcmp(param, "LSSTATIONS") == 0) {
            if ((istat = get_ls_stations(strchr(line, ' '))) < 0) {
                nll_puterr("ERROR: reading LSSTATIONS parameters.");
                flag_stations = 0;
            }
        }

        /* read grid mode names */

        if (strcmp(param, "LSMODE") == 0) {
            if ((istat = get_ls_mode(strchr(line, ' '))) < 0) {
                nll_puterr("ERROR: reading LSMODE mode.");
                flag_mode = 0;
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
        nll_puterr("ERROR: no inputfile (LSINP) params read.");
    if (!flag_outfile)
        nll_puterr("ERROR: no outputfile (LSOUT) params read.");
    if (!flag_grid)
        nll_puterr("ERROR: no grid (LSGRID) params read.");
    if (!flag_params)
        nll_puterr("ERROR: no parameters (LSPARAMS) read.");

    if (!flag_trans) {
        snprintf(MsgStr, sizeof (MsgStr), "INFO: no transformation (TRANS) params read.");
        nll_putmsg(1, MsgStr);
        Hypocenter.comment[0] = '\0';
    }

    ihave_time_input_grids = flag_out_grid * flag_nlloc_outfile;

    //printf("DEBUG: %d %d %d %d %d %d %d %d %d %d\n", flag_include, flag_control, flag_method, flag_phase_id, flag_inpfile, flag_outfile, flag_grid, flag_params, flag_phstat, flag_stations);
    return (flag_include * flag_control * flag_method * flag_phase_id * flag_inpfile * flag_outfile * flag_grid * flag_params * flag_phstat * flag_mode * flag_stations - 1);
}

/** function to calculate distance between 2 XYZ points */

double Dist3D2_Loc2ssst(double x1, double x2, double y1, double y2, double z1, double z2, int latlon) {

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

int DoLoc2ssst() {

    int istat;

    double Len3Mean = 0.0, RMSMean = 0.0;
    double NRdgsMean = 0.0, GapMean = 0.0;
    int Len3Reject = 0, RMSReject = 0;
    int NRdgsReject = 0, GapReject = 0;
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

    LocNode *loc_list_head = NULL; // root node of location list
    PhsNode *phs_list_head = NULL; // root node of station/phase arrivals list

    char arrival_phase[PHASE_LABEL_LEN];

    //SourceDesc* Srce = NULL;


    // read input NLL Location files

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
                continue;
            } else if (!Params.use_rejected && strcmp(hypo.locStat, "REJECTED") == 0) {
                //nll_puterr("WARNING: location REJECTED, ignoring event");
                continue;
            }


            Len3Mean += hypo.ellipsoid.len3;
            RMSMean += hypo.rms;
            NRdgsMean += (double) hypo.nreadings;
            GapMean += hypo.gap;
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
            if (iReject)
                continue;

            // initialize some hypo fields
            latlon2rect(0, hypo.dlat, hypo.dlong, &(hypo.x), &(hypo.y));
            hypo.z = hypo.depth;

            // do phase mapping before saving location arrivals and before accumulating time corrections
            for (int i = 0; i < NumArrivals; i++) {
                // DEBUG
                /*if (strcmp(Arrival[i].label, "PB07") == 0 && strcmp(Arrival[i].phase, "P") == 0) {
                if (strcmp(Arrival[i].label, "PB07") == 0) {
                    printf("DEBUG: PB07 %s found in arrival list res=%f hyp=%s\n", Arrival[i].phase, Arrival[i].residual, fn_hyp_in_list[nFile]);
                }
                 */
                strcpy(arrival_phase, Arrival[i].phase);
                EvalPhaseID(Arrival[i].phase, arrival_phase);
                strcpy(Arrival[i].time_grid_label, Arrival[i].label);
                // update swap bytes
                Arrival[i].gdesc.iSwapBytes = iSwapBytesOnInput;
            }
            int iuse_phaseid_in_label = 1;
            int i_check_station_has_XYZ_coords = 1; // 20210903 AJL - Bug fix: do not add stations for arrivals without coordinates
            NumStationPhases = addToStationList(StationPhaseList, NumStationPhases, Arrival, NumArrivals, iuse_phaseid_in_label, i_check_station_has_XYZ_coords);

            Location* ploc_list_node = newLocation(
                    cloneHypoDesc(&hypo),
                    cloneArrivalDescArray(Arrival, NumArrivals),
                    NumArrivals, NULL,
                    NULL,
                    NULL);
            loc_list_head = addLocationToLocList(&loc_list_head, ploc_list_node, nLocAccepted);

            nLocAccepted++;


        }
        fclose(fp_hypo);

    }
    fprintf(OUT_LEVEL_1, "\n");

    // write message
    fprintf(stdout,
            "%d location files read, %d accepted.\n",
            numFiles, nLocAccepted);
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

    nHypos = nLocAccepted;



    // loop over station/phase, accumulate corrections in output grid
    fprintf(stdout, "\nAccumulating SSST corrections in output grid for %d station/phases:\n", NumStationPhases);

    // allocate grids
    GridDesc ssst_grid = ssst_grid_template;
    GridDesc *pssst_grid = &ssst_grid;
    DuplicateGrid(pssst_grid, &ssst_grid_template, ssst_grid_template.chr_type);
    GridDesc ssst_time_grid = ssst_time_grid_template;
    GridDesc *pssst_time_grid = &ssst_time_grid;
    DuplicateGrid(pssst_time_grid, &ssst_time_grid_template, ssst_time_grid_template.chr_type);

    // convert source location coordinates
    istat = ConvertSourceLoc(0, StationPhaseList, NumStationPhases, 1, 1);

    SourceDesc *station_phase;
    char stacode[SOURCE_LABEL_LEN];
    char phasecode[SOURCE_LABEL_LEN];
    // 20210622 AJL - Bug fix: replaced possibly non-existent arrival pointer (PhsNodeArray[0]->parrival) with temp arrival pointer.
    ArrivalDesc arrival_tmp;
    ArrivalDesc* parr_tmp = &arrival_tmp;
    for (int n = 0; n < NumStationPhases; n++) {

        // initialize ssst grid for this station/phase
        station_phase = StationPhaseList + n;
        strcpy(stacode, station_phase->label);
        *strchr(stacode, '#') = '\0';
        strcpy(phasecode, strchr(station_phase->label, '#') + 1);
        // rename to station code only
        strcpy(station_phase->label, stacode);
        if (0) {
            // DEBUG
            if (strcmp(stacode, "PRSN") != 0 || strcmp(phasecode, "P") != 0) {
                continue;
            }
            // END - DEBUG
        }
        // check if station requested
        //printf("DEBUG: <%s> %s\n", Params.stations, stacode);
        if (strlen(Params.stations) > 0 && strstr(Params.stations, stacode) == NULL) {
            continue;
        }
        fprintf(stdout, "----------\n");
        fprintf(stdout, "%d/%d %s %s  x %lf y %lf z %lf  lat %lf lon %lf  xyz %d  latlon %d\n",
                n, NumStationPhases, stacode, phasecode, station_phase->x, station_phase->y, station_phase->z,
                station_phase->dlat, station_phase->dlong, station_phase->is_coord_xyz, station_phase->is_coord_latlon);
        if (station_phase->x < -LARGE_DOUBLE / 2.0 || station_phase->y < -LARGE_DOUBLE / 2.0 || station_phase->z < -LARGE_DOUBLE / 2.0) {
            snprintf(MsgStr, sizeof (MsgStr), "WARNING: No station coordinates available: %s %s", stacode, phasecode);
            nll_puterr(MsgStr);
            continue;
        }

        // loop over locations, look for arrival for this stacode and phasecode
        double residualmax = LARGE_FLOAT;
        if (IsPhaseID(phasecode, "P")) {
            residualmax = PhsStat.PResidualMax;
        } else if (IsPhaseID(phasecode, "S")) {
            residualmax = PhsStat.SResidualMax;
        }
        LocNode* locNode = NULL;
        int loc_id = 0;
        int phs_id = 0;
        while ((locNode = getLocationFromLocList(loc_list_head, loc_id)) != NULL) {
            ArrivalDesc* parrivals = locNode->plocation->parrivals;
            int narrivals = locNode->plocation->narrivals;
            // loop over arrivals, look for arrival for this stacode and phasecode
            ArrivalDesc* parr;
            for (int narr = 0; narr < narrivals; narr++) {
                parr = parrivals + narr;
                if (fabs(parr->residual) > residualmax || strcmp(parr->label, stacode) != 0 || strcmp(parr->phase, phasecode) != 0) {
                    //nll_puterr("WARNING: residual is greater than residualmax, ignoring arrival");
                    continue;
                }
                int addDuplicates = 1;
                PhsNode *phsNode = addArrivalToPhaseList(&phs_list_head, parr, phs_id, addDuplicates);
                addRemoveLocationInAssocLocationsList(phsNode, loc_id, 1);
                //printf("DEBUG: addRemoveLocationInAssocLocationsList phs_id %d loc_id%d\n", phs_id, phsNode->passoc_locations[0]);
                phs_id++;
            }
            loc_id++;
        }
        // put locations in array for efficiency
        NumLocNodes = 0;
        while ((locNode = getLocationFromLocList(loc_list_head, NumLocNodes)) != NULL) {
            if (NumLocNodes >= MAX_NUM_INPUT_FILES) {
                snprintf(MsgStr, sizeof (MsgStr), "ERROR: size LocNodeArray exceeded, only first phases %d will be processed.", MAX_NUM_INPUT_FILES);
                nll_puterr(MsgStr);
                break;
            }
            LocNodeArray[NumLocNodes++] = locNode;
        }
        // put phases in array for efficiency
        PhsNode *phsNode;
        NumPhsNodes = 0;
        while ((phsNode = getPhsNodeFromPhaseList(phs_list_head, NumPhsNodes)) != NULL) {
            if (NumPhsNodes >= X_MAX_NUM_ARRIVALS) {
                snprintf(MsgStr, sizeof (MsgStr), "ERROR: size PhsNodeArray exceeded, only first phases %d will be processed.", X_MAX_NUM_ARRIVALS);
                nll_puterr(MsgStr);
                break;
            }
            PhsNodeArray[NumPhsNodes++] = phsNode;
        }
        // clean up
        // 20201104 AJL - Buf fix: moved this free to end of block, PhsNodeArray points to allocations made in addArrivalToPhaseList()
        //freePhaseList(phs_list_head, 0);
        //phs_list_head = NULL;

        //display_grid_param(pssst_grid);

        // get ssst corrections for this station/phase
        int ix, iy, iz;
        GRID_FLOAT_TYPE xval, yval, zval, ssstval, ssstval_sum = 0.0;
        xval = pssst_grid->origx;
        for (ix = 0; ix < pssst_grid->numx; ix++) {
            fprintf(stdout, "\r%d%%  ", (int) (0.5 + 100.0 * ((float) ix / (pssst_grid->numx - 1))));
            fflush(stdout);
            yval = pssst_grid->origy;
            for (iy = 0; iy < pssst_grid->numy; iy++) {
                //fprintf(stdout, "\rx%d%% y%d%%  ", (int) (0.5 + 100.0 * ((float) ix / pssst_grid->numx)), (int) (0.5 + 100.0 * ((float) iy / pssst_grid->numy)));
                //fflush(stdout);
                zval = pssst_grid->origz;
                for (iz = 0; iz < pssst_grid->numz; iz++) {
                    ssstval = (GRID_FLOAT_TYPE) get_ssst_value(xval, yval, zval, PhsNodeArray, NumPhsNodes, LocNodeArray, NumLocNodes, &Params);
                    if (ssstval > -LARGE_FLOAT) {
                        ((GRID_FLOAT_TYPE ***) pssst_grid->array)[ix][iy][iz] = ssstval;
                    }
                    ssstval_sum += ssstval;
                    zval += pssst_grid->dz;
                }
                yval += pssst_grid->dy;
            }
            xval += pssst_grid->dx;
        }
        fprintf(stdout, "ssstval: mean: %f last: %f\n", ssstval_sum / (double) (pssst_grid->numx * pssst_grid->numy * pssst_grid->numz), ssstval);

        // write ssst correction grid to disk
        char filename[MAXLINE_LONG];
        sprintf(filename, "%s.%s.%s", fn_ls_output, phasecode, stacode);
        sprintf(MsgStr, "Finished calculation, SSST grid output files: %s.*", filename);
        nll_putmsg(1, MsgStr);
        istat = WriteGrid3dBuf(pssst_grid, station_phase, filename, "ssst");
        if (istat < 0) {
            nll_puterr2("ERROR: writing SSST grid to disk", filename);
            return (-1);
        }
        // save station coordinates to file
        char fn_stations[FILENAME_MAX];
        FILE* fp_stations;
        sprintf(fn_stations, "%s.ssst.stations", filename);
        if ((fp_stations = fopen(fn_stations, "w")) != NULL) {
            fprintf(fp_stations, "%s_%s %f %f %f", stacode, phasecode, station_phase->x, station_phase->y, station_phase->z);
            fclose(fp_stations);
        } else {
            if (message_flag >= 1)
                nll_putmsg2(1, "INFO: cannot open station list file", fn_stations);
            // will use stations found in hypocenter-phase files
        }

        // update input time grids to create ssst corrected time grids for location
        if (ihave_time_input_grids) {
            nll_putmsg2(1, "INFO: Opening existing time grid", fn_time_input);
            double tfact = 1.0;
            ArrivalDesc* parr;
            // DEBUG!!
            //int message_flag_save = message_flag;
            //message_flag = 99;
            // END DEBUG
            istat = open_traveltime_grid(parr_tmp, fn_time_input, stacode, phasecode, VpVsRatio, &tfact);
            // DEBUG!!
            //message_flag = message_flag_save;
            // END DEBUG
            if (istat < 0) {
                nll_puterr2("ERROR: Opening existing time grid", parr_tmp->fileroot);
                return (-1);
            }
            // reset ssst grid type to TIME
            add_ssst_to_traveltime_grid(phasecode, stacode, pssst_grid, &(parr_tmp->gdesc), pssst_time_grid, station_phase, tfact);
            //sprintf(filename, "%s_ssst_corr.%s.%s", fn_ls_output, phasecode, stacode); // 20201010 AJL - bug fix: adding _ssst_corr to file root makes file management difficult
            sprintf(filename, "%s.%s.%s", fn_ls_output, phasecode, stacode);
            // remove existing files since may be links  // 20201010 AJL - added
            char fname_remove[MAXLINE_LONG];
            sprintf(fname_remove, "%s.%s.hdr", filename, "time");
            remove(fname_remove);
            sprintf(fname_remove, "%s.%s.buf", filename, "time");
            remove(fname_remove);
            istat = WriteGrid3dBuf(pssst_time_grid, station_phase, filename, "time");
            if (istat < 0) {
                nll_puterr2("ERROR: writing SSST corrected time grid to disk", filename);
                return (-1);
            }
            nll_putmsg2(1, "INFO: SSST corrected time grid written to disk", filename);

            // angles
            if (angle_mode == ANGLE_MODE_YES || angle_mode == ANGLE_MODE_INCLINATION) {
                GridDesc angle_grid;
                if (angle_mode == ANGLE_MODE_YES) {
                    DuplicateGrid(&angle_grid, pssst_time_grid, "ANGLE");
                    if ((istat = GenAngleGrid(pssst_time_grid, station_phase, filename,
                            &angle_grid, angle_mode)) < 0)
                        nll_puterr("ERROR: calculating take-off angles.");
                } else if (angle_mode == ANGLE_MODE_INCLINATION) {
                    GridDesc angle_grid;
                    DuplicateGrid(&angle_grid, pssst_time_grid, "INCLINATION");
                    if ((istat = GenAngleGrid(pssst_time_grid, station_phase, filename,
                            &angle_grid, angle_mode)) < 0)
                        nll_puterr("ERROR: calculating inclination angles.");
                }
                DestroyGridArray(&angle_grid);
                FreeGrid(&angle_grid);
            }

            // free grid
            //printf("DEBUG: DestroyGridArray %s %s\n", phasecode, stacode);
            DestroyGridArray(&(parr_tmp->gdesc));
            //printf("DEBUG: FreeGrid %s %s\n", phasecode, stacode);
            FreeGrid(&(parr_tmp->gdesc));
        }

        // make Grid2GMT cpt file
        double cpt_start_stop = PhsStat.PResidualMax / 2.0;
        double cpt_increment;
        if (cpt_start_stop <= 1.0 + FLT_MIN) {
            cpt_start_stop = 1.0;
            cpt_increment = 0.2;
        } else if (cpt_start_stop <= 2.0) {
            cpt_start_stop = 2.0;
            cpt_increment = 0.5;
        } else if (cpt_start_stop <= 5.0) {
            cpt_start_stop = 5.0;
            cpt_increment = 1.0;
        } else if (cpt_start_stop <= 10.0) {
            cpt_start_stop = 10.0;
            cpt_increment = 2.0;
        } else {

            cpt_start_stop = 10.0;
            cpt_increment = 2.0;
        }
        char system_str[MAXLINE_LONG];
        sprintf(system_str, "makecpt -Z -Cred2green -T-%f/%f/%f > Grid2GMT_SSST.cpt", cpt_start_stop, cpt_start_stop, cpt_increment);
        system(system_str);
        sprintf(MsgStr, "Grid2GMT cpt file written to: %s", "Grid2GMT_SSST.cpt");
        nll_putmsg(1, MsgStr);

        // clean up
        // 20201104 AJL - Buf fix: moved this free from above
        freePhaseList(phs_list_head, 0);
        phs_list_head = NULL;

    }

    return (0);

}

/** function to calculate ssst correction value at an xyz point for a specified station and phase */

double get_ssst_value(double xval, double yval, double zval,
        PhsNode **phs_node_array, int num_phs_nodes, LocNode **loc_node_array, int num_loc_nodes, LS_Params *pparams) {

    double char_dist2 = pparams->char_dist * pparams->char_dist;
    double weight_floor = pparams->weight_floor;

    double ssst_corr_sum = 0.0;
    double ssst_wt_sum = 0.0;
    int nssst_corr = 0;

    // loop over phase node arrivals, look for arrival for this stacode and phasecode

    double xloc, yloc, zloc;
    int latlon = GeometryMode == MODE_GLOBAL;
    LocNode *locNode = NULL;
    PhsNode *phsNode = NULL;
    for (int nphs = 0; nphs < num_phs_nodes; nphs++) {

        phsNode = *(phs_node_array + nphs);

        //printf("\nDEBUG: phsNode %d ->passoc_locations[0] %d\n", phs_id, phsNode->passoc_locations[0]);
        locNode = *(loc_node_array + phsNode->passoc_locations[0]);
        if (latlon) {
            xloc = locNode->plocation->phypo->dlong;
            yloc = locNode->plocation->phypo->dlat;
            zloc = locNode->plocation->phypo->z;
        } else {
            xloc = locNode->plocation->phypo->x;
            yloc = locNode->plocation->phypo->y;
            zloc = locNode->plocation->phypo->z;
        }
        // get distance from event to center of node
        double event_node_dist2 = Dist3D2_Loc2ssst(xloc, xval, yloc, yval, zloc, zval, latlon);
        //event_node_dist = Dist2D(xloc, xval, yloc, yval); // TEST ONLY!
        // get weight
        double weight = exp(-(event_node_dist2 / char_dist2));
        weight += weight_floor;
        //if (weight < weight_floor) {
        //    weight = weight_floor;
        //}
        ArrivalDesc* parr = phsNode->parrival;
        ssst_corr_sum += (parr->residual + parr->delay) * weight;
        ssst_wt_sum += weight;
        nssst_corr++;
        //if (fabs(xval - 0.0) < 0.01 && fabs(yval - 0.0) < 0.01 && fabs(zval - 5.0) < 0.01) {
        //if (fabs(xval - -103.5) < 0.05 && fabs(yval - 31.4) < 0.01 && fabs(zval - 5.0) < 0.05) {
        //printf("\nDEBUG: phs_id %d xval %f yval %f zval %f xloc %f yloc %f zloc %f parr->residual %f parr->delay %f dist %f weight %f ssst_corr_sum %f \n", nphs, xval, yval, zval, xloc, yloc, zloc, parr->residual, parr->delay, sqrt(event_node_dist2), weight, ssst_corr_sum);
        //}

    }

    double ssst_corr = 0.0;
    if (ssst_wt_sum > FLT_MIN) {

        ssst_corr = (ssst_corr_sum / ssst_wt_sum);
    }

    //printf("\nDEBUG: ============= ssst_corr_sum %f ssst_wt_sum %f ssst_corr %f \n", ssst_corr_sum, ssst_wt_sum, ssst_corr);
    return (ssst_corr);

}

/** function to find and open existing travel-time grids */

int open_traveltime_grid(ArrivalDesc* parr, char *fn_time_grid_input, char *stacode, char *phasecode, double vp_vs_ratio, double *ptfact) {

    int istat;
    int DEBUG = 0;

    char filename[MAXLINE_LONG];
    //sprintf(filename, "%s.%s.%s", fn_time_grid_input, stacode, phasecode);

    //char arrival_phase[PHASE_LABEL_LEN];
    //strcpy(arrival_phase, parr->phase);

    *ptfact = 1.0;

    // try to open time grid file using original phase ID
    //sprintf(parr->fileroot, "%s.%s.%s", fn_time_grid_input, phasecode, parr->time_grid_label);
    sprintf(parr->fileroot, "%s.%s.%s", fn_time_grid_input, phasecode, stacode);
    sprintf(filename, "%s.time", parr->fileroot);
    if (DEBUG) {
        printf("DEBUG: parr->fileroot <%s>\n", parr->fileroot);
        printf("DEBUG: filename <%s>\n", filename);
        nll_puterr2("INFO: Opening existing time grid: ", filename);
    }
    // try opening time grid file for this phase
    istat = OpenGrid3dFile(filename,
            &(parr->fpgrid),
            &(parr->fphdr),
            &(parr->gdesc), "time",
            &(parr->station),
            parr->gdesc.iSwapBytes);
    if (DEBUG && istat < 0) {
        nll_puterr2("WARNING: Cannot open existing time grid: ", filename);
    }

    /* already mapped earlier
    if (istat < 0) {
        // try to open time grid file using LOCPHASEID mapped phase ID
        EvalPhaseID(eval_phase, arrival_phase);
        sprintf(parr->fileroot, "%s.%s.%s", fn_grids, eval_phase, parr->time_grid_label);
        sprintf(filename, "%s.time", parr->fileroot);
        // try opening time grid file for this phase
        istat = OpenGrid3dFile(filename,
                &(parr->fpgrid),
                &(parr->fphdr),
                &(parr->gdesc), "time",
                &(parr->station),
                parr->gdesc.iSwapBytes);
    }*/

    /* try opening P time grid file for S if no P companion phase */
    if (istat < 0 && vp_vs_ratio > 0.0 && IsPhaseID(phasecode, "S")) {
        *ptfact = vp_vs_ratio;
        sprintf(parr->fileroot, "%s.%s.%s", fn_time_grid_input, "P", stacode);
        sprintf(filename, "%s.time", parr->fileroot);
        if (DEBUG && istat < 0) {
            nll_puterr2("INFO: Opening existing time grid: ", filename);
        }
        istat = OpenGrid3dFile(filename,
                &(parr->fpgrid),
                &(parr->fphdr),
                &(parr->gdesc), "time",
                &(parr->station),
                parr->gdesc.iSwapBytes);
        if (DEBUG && istat < 0) {
            nll_puterr2("WARNING: Cannot open existing time grid: ", filename);
        }
        if (message_flag >= 3) {
            sprintf(MsgStr,
                    "INFO: S phase: using P phase travel time grid file: %s", filename);
            nll_putmsg(3, MsgStr);
        }
    }


    // try opening DEFAULT time grid file
    int i_need_elev_corr = 0; // TODO: implement elevation corr
    if (istat < 0) {
        SourceDesc* pstation = FindSource(stacode);
        if (pstation != NULL) {
            // open DEFAULT time grid for this phase
            int iSwapBytes = parr->gdesc.iSwapBytes;
            // 20201022 AJL - Cluge, assume need to byte swap if DEFAULT  // TODO: make this automatic or configurable
            iSwapBytes = 1;
            //
            sprintf(filename, "%s.%s.%s.time", fn_time_grid_input, phasecode, "DEFAULT");
            if (DEBUG && istat < 0) {
                nll_puterr2("INFO: Opening existing time grid: ", filename);
            }
            istat = OpenGrid3dFile(filename,
                    &(parr->fpgrid),
                    &(parr->fphdr),
                    &(parr->gdesc), "time",
                    &(parr->station),
                    iSwapBytes);
            if (DEBUG && istat < 0) {
                nll_puterr2("WARNING: Cannot open existing time grid: ", filename);
            }
            if (istat >= 0 && message_flag >= 3) {
                sprintf(MsgStr,
                        "INFO: using DEFAULT travel time grid file: %s", filename);
                nll_putmsg(3, MsgStr);
            }
            parr->station = *pstation;
            i_need_elev_corr = 1;
        }
    }

    if (istat < 0) {
        return (istat);
    }

    // allocate grids
    GridDesc *ptime_grid = &(parr->gdesc);
    ptime_grid->buffer = AllocateGrid(ptime_grid);
    if (ptime_grid->buffer == NULL) {
        nll_puterr(
                "ERROR: allocating memory for search PDF grid buffer.");
        return (EXIT_ERROR_MEMORY);
    }
    // create grid array access pointers
    ptime_grid->array = CreateGridArray(ptime_grid);
    if (ptime_grid->array == NULL) {
        nll_puterr(
                "ERROR: creating array for accessing search PDF grid buffer.");
        return (EXIT_ERROR_MEMORY);
    }
    // read grid
    if ((istat =
            ReadGrid3dBuf(ptime_grid, parr->fpgrid)) < 0) {
        nll_puterr("ERROR: reading search PDF grid from disk.");
        return (EXIT_ERROR_FILEIO);
    }
    CloseGrid3dFile(ptime_grid, &parr->fpgrid, &parr->fphdr);

    return (istat);
}

/** add original traveltimes to ssst corrections for a station/phase ssst grid
 */
int add_ssst_to_traveltime_grid(char *phasecode, char *stacode, GridDesc *pssst_grid, GridDesc *ptraveltime_grid, GridDesc *pssst_time_grid, SourceDesc* psrce, double tfact) {


    int is3D = 0;
    if (ptraveltime_grid->type == GRID_TIME) {
        is3D = 1;
    }

    int debug_count = 0;
    int ix, iy, iz;
    double xval, yval, zval, ttval, arrival_dist, ssst_corr;
    xval = pssst_time_grid->origx;
    for (ix = 0; ix < pssst_time_grid->numx; ix++) {
        yval = pssst_time_grid->origy;
        for (iy = 0; iy < pssst_time_grid->numy; iy++) {
            //fprintf(stdout, "\rx%d/%d y%d/%d  ", ix, pssst_time_grid->numx, iy, pssst_time_grid->numy);
            //fflush(stdout);
            zval = pssst_time_grid->origz;
            for (iz = 0; iz < pssst_time_grid->numz; iz++) {
                if (is3D) {
                    ttval = ReadAbsInterpGrid3d(NULL, ptraveltime_grid, xval, yval, zval, 1);
                } else {
                    arrival_dist = GetEpiDist(psrce, xval, yval);
                    ttval = ReadAbsInterpGrid2d(NULL, ptraveltime_grid,
                            GeometryMode == MODE_GLOBAL ? arrival_dist * KM2DEG : arrival_dist, zval);
                }
                if (ttval < -LARGE_FLOAT) {
                    if (debug_count++ < 10) {
                        printf("DEBUG ERROR: ttval < -LARGE_FLOAT: %s %s  is3D %d, xval %f, yval %f, zval %f\n",
                                phasecode, stacode, is3D, xval, yval, zval);
                    }
                    ttval = 0.0;
                    ssst_corr = 0.0;
                } else {
                    ssst_corr = ReadAbsInterpGrid3d(NULL, pssst_grid, xval, yval, zval, 1);
                    if (ssst_corr < -LARGE_FLOAT) {
                        if (debug_count++ < 10) {
                            printf("DEBUG ERROR: ssst_corr < -LARGE_FLOAT: %s %s  is3D %d, xval %f, yval %f, zval %f\n",
                                    phasecode, stacode, is3D, xval, yval, zval);
                        }
                        ssst_corr = 0.0;
                    }
                }
                if (fabs(ttval * tfact + ssst_corr) < SMALL_FLOAT) {
                    if (debug_count++ < 10) {
                        printf("DEBUG ERROR: fabs(ttval * tfact + ssst_corr) < SMALL_FLOAT: %s %s  is3D %d, xval %f, yval %f, zval %f\n",
                                phasecode, stacode, is3D, xval, yval, zval);
                    }
                }
                ((GRID_FLOAT_TYPE ***) pssst_time_grid->array)[ix][iy][iz] = (GRID_FLOAT_TYPE) (ttval * tfact + ssst_corr);
                zval += pssst_time_grid->dz;
            }
            yval += pssst_time_grid->dy;
        }
        xval += pssst_time_grid->dx;
    }
    if (debug_count > 0) {
                      printf("DEBUG ERROR: %d total errors\n", debug_count);
    }

    return (0);

}

/*** function to generate take-off angle grid */

int GenAngleGrid(GridDesc* ptgrid, SourceDesc* psource, char *filename, GridDesc* pagrid, int angle_mode) {

    int istat;

    double xsource, ysource, zsource;

    int grid_mode = GRID_TIME; // must be GRID_TIME for SSST time grid

    /* check grid mode, make appropriate adjustments */

    xsource = psource->x;
    ysource = psource->y;
    zsource = psource->z;


    /* generate angle grid */

    /*if (angle_calc_meth == ANGLE_METHOD_GRADIENT)*/
    if (1) {
        /* check things */
        if (ptgrid->type != GRID_TIME && ptgrid->type != GRID_TIME_2D) {
            nll_puterr(
                    "ERROR: Gradient take-off angle algorithm requires TIME grid.");
            return (-1);
        }
        if (ptgrid->dx != ptgrid->dy || ptgrid->dx != ptgrid->dz) {
            nll_puterr(
                    "ERROR: Gradient take-off angle algorithm requires cubic grid, i.e. dx=dy=dz.");
            return (-1);
        }

        /* run gradient take-off angle algorithm */
        if ((istat = CalcAnglesGradient(ptgrid, pagrid, angle_mode, grid_mode)) < 0)
            return (-1);

    }


    /* save angle grid to disk */

    sprintf(MsgStr,
            "Finished calculation, take-off angles grid output files: %s.*",
            filename);
    nll_putmsg(1, MsgStr);
    /* need only ix=0 sheet for 2D grids */
    if (angle_mode == ANGLE_MODE_YES)
        istat = WriteGrid3dBuf(pagrid, psource, filename, "angle");
    else if (angle_mode == ANGLE_MODE_INCLINATION)
        istat = WriteGrid3dBuf(pagrid, psource, filename, "inclination");
    if (istat < 0) {
        nll_puterr("ERROR: writing take-off angles grid to disk.");
        return (-1);
    }


    return (0);

}
