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


/*   NLLocInput.c

        NLLoc control-file (input) parsing for NonLinLoc: is_nll_control_json,
        ReadNLLoc_Input and the GetNLLoc_* / Get* statement parsers.

        Extracted verbatim from NLLocLib.c during the modernization effort
        (Phase 3: modularization). Pure code movement, no logic change; all
        function prototypes remain in NLLocLib.h, so callers are unaffected.
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


/** function to check if nll input file is nll-control JSON format */

// AJL 20211014 - added

int is_nll_control_json(FILE* fp_input) {

    char line_buf[MAXLINE_LONG];

    // read each input line and check for nll-control JSON name tag

    while (fgets(line_buf, MAXLINE_LONG, fp_input) != NULL) { // read next line

        // skip comment line
        if (strncmp(line_buf, "#", 1) == 0) {
            continue;
        }

        // check for JSON name tag
        if (strstr(line_buf, "nll-control") != NULL) {
            rewind(fp_input);
            return (1);
        }
    }

    rewind(fp_input);

    return (0);

}

/** function to read input control file */

// AJL 20071217 - added char** passing of parameters

int ReadNLLoc_Input(FILE* fp_input, char** param_line_array, int n_param_lines) {
    int istat, iscan;
    char param[MAXLINE] = "\0", *pchr;
    char line_buf[MAXLINE_LONG];
    char *line;
    char *fgets_return;

    int flag_control = 0, flag_outfile = 0, flag_grid = 0, flag_search = 0, flag_prior = 0, flag_posterior = 0,
            flag_method = 0, flag_comment = 0, flag_signature = 0,
            flag_hyptype = 0, flag_gauss = 0, flag_gauss2 = 0, flag_trans = 0, flag_comp = 0,
            flag_phstat = 0, flag_phase_id = 0, flag_sta_wt = 0, flag_qual2err = 0,
            flag_mag = 0, flag_alias = 0, flag_exclude = 0, flag_include = 0, flag_time_delay = 0,
            flag_topo_surface = 0, flag_time_delay_surface = 0, flag_elev_corr = 0,
            flag_otime = 0, flag_angles = 0, flag_source = 0;
    int flag_include_file = 1;

    int ok_search_pdf = 1;


    char **param_lines;


    // param lines not needed by NLDiffLoc
    if (nll_mode == MODE_DIFFERENTIAL) {
        flag_search = 1;
    }


    // flag some memory that may be allocated in reading input control file
    SearchPrior.coherence = NULL;
    SearchPrior.weight = NULL;
    SearchPosterior.coherence = NULL;
    SearchPosterior.weight = NULL;
    SearchPosterior.first_motion_arrivals = NULL;
    SearchPosterior.nfirst_motion_arrivals = NULL;


    /* read each input line */

    line = line_buf; // for reading from file case
    param_lines = param_line_array;
    fgets_return = NULL;

    while (
            // reading from file (control file or included file
            (fp_input != NULL && (fgets_return = fgets(line, MAXLINE_LONG, fp_input)) != NULL)
            ||
            fp_include != NULL
            ||
            // reading from string
            (param_line_array != NULL && ((param_lines - param_line_array) < n_param_lines)
            && (fgets_return = line = *(param_lines++)) != NULL)
            ) {

        /* check for end of include file */

        if (fgets_return == NULL && fp_include != NULL) {
            SwapBackIncludeFP(&fp_input);
            continue;
        }

        //printf(line);

        istat = -1;

        /*read parameter line */

        if ((iscan = sscanf(line, "%s", param)) < 0)
            continue;

        /* skip comment line or white space */

        if (strncmp(param, "#", 1) == 0 || isspace(param[0]))
            istat = 0;



        /* read include file params and set input to include file */

        /*
           if (strcmp(param, "INCLUDE") == 0)
            if ((istat = GetIncludeFile(strchr(line, ' '), &fp_input)) < 0) {
                nll_puterr("ERROR: processing include file.");
                flag_include_file = 0;
            }
         */


        /*
         * 20100913 Jan Becker - bug fix
         a bugfix for NLL that fixes a segfault when the include statement
within the control buffer is used along with the library call. The (const)
input parameter array was used to load the file contents what was wrong. The
patch below fixes this problem.
         */

        if (strcmp(param, "INCLUDE") == 0) {
            if ((istat = GetIncludeFile(strchr(line, ' '), &fp_input)) < 0) {
                nll_puterr("ERROR: processing include file.");
                flag_include_file = 0;
            } else
                line = line_buf;
        }




        /* read control params */

        if (strcmp(param, "CONTROL") == 0) {
            if ((istat = get_control(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading control params.");
            else
                flag_control = 1;
        }

        /* read transform params */

        if (strcmp(param, "TRANS") == 0) {
            if ((istat = get_transform(0, strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading transformation parameters.");
            else
                flag_trans = 1;
        }


        /* read grid params */

        if (strcmp(param, "LOCGRID") == 0) {
            if ((istat = GetNLLoc_Grid(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading grid parameters.");
            else
                flag_grid = 1;
        }


        /* read file names */

        if (strcmp(param, "LOCFILES") == 0) {
            if ((istat = GetNLLoc_Files(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading NLLoc output file name.");
            else
                flag_outfile = 1;
        }


        /* read output file types names */

        if (strcmp(param, "LOCHYPOUT") == 0) {
            if ((istat = GetNLLoc_HypOutTypes(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading NLLoc hyp output file types.");
            else
                flag_hyptype = 1;
        }


        /* read search type */

        if (strcmp(param, "LOCSEARCH") == 0) {
            if ((istat = GetNLLoc_SearchType(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading NLLoc search type.");
            else
                flag_search = 1;
        }


        /* read search prior */
        // 20190510 AJL - added

        if (strcmp(param, "LOCPRIOR") == 0) {
            if ((istat = GetNLLoc_PdfGrid(strchr(line, ' '), PDF_GRID_PRIOR)) < 0)
                nll_puterr("ERROR: reading NLLoc search prior PDF.");
            else
                flag_prior = 1;
        }


        /* read search posterior */
        // 201900610 AJL - added

        if (strcmp(param, "LOCPOSTERIOR") == 0) {
            if ((istat = GetNLLoc_PdfGrid(strchr(line, ' '), PDF_GRID_POSTERIOR)) < 0)
                if (istat == -9) {
                    ok_search_pdf = 0;
                } else {
                    nll_puterr("ERROR: reading NLLoc search posterior PDF.");
                } else
                flag_posterior = 1;
        }


        /* read method */

        if (strcmp(param, "LOCMETH") == 0) {
            if ((istat = GetNLLoc_Method(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading NLLoc method.");
            else
                flag_method = 1;
        }


        /* read fixed origin time parameters */

        if (strcmp(param, "LOCFIXOTIME") == 0) {
            if ((istat = GetNLLoc_FixOriginTime(
                    strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading NLLoc fixed origin time params.");
            else
                flag_otime = 1;
        }


        /* read phase identifier values */

        if (strcmp(param, "LOCPHASEID") == 0) {
            if ((istat = GetPhaseID(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading phase identifier values.");
            else
                flag_phase_id = 1;
        }


        /* read station distance weighting values */

        if (strcmp(param, "LOCSTAWT") == 0) {
            if ((istat = GetStaWeight(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading station distance weighting values.");
            else
                flag_sta_wt = 1;
        }


        /* read quality2error values */

        if (strcmp(param, "LOCQUAL2ERR") == 0) {
            if ((istat = GetQuality2Err(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading quality2error values.");
            else
                flag_qual2err = 1;
        }


        /* read comment */

        if (strcmp(param, "LOCCOM") == 0) {
            snprintf(Hypocenter.comment, sizeof(Hypocenter.comment), "%s", strchr(line, ' ') + 1);
            //*(strchr(Hypocenter.comment, '\n')) = '\0';
            TrimString(Hypocenter.comment);
            snprintf(MsgStr, sizeof(MsgStr), "LOCCOMMENT:  %s\n", Hypocenter.comment);
            nll_putmsg(3, MsgStr);
            flag_comment = 1;
        }


        /* read signature */

        if (strcmp(param, "LOCSIG") == 0) {
            snprintf(LocSignature, sizeof(LocSignature), "%s", strchr(line, ' ') + 1);
            //*(strchr(LocSignature, '\n')) = '\0';
            TrimString(LocSignature);
            snprintf(MsgStr, sizeof(MsgStr), "LOCSIGNATURE:  %s\n",
                    LocSignature);
            nll_putmsg(3, MsgStr);
            flag_signature = 1;
        }


        /* read gauss params */

        if (strcmp(param, "LOCGAU2") == 0) {
            if ((istat = GetNLLoc_Gaussian2(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading Gaussian2 parameters.");
            else
                flag_gauss2 = 1;
        }


        /* read gauss params */

        if (strcmp(param, "LOCGAU") == 0) {
            if ((istat = GetNLLoc_Gaussian(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading Gaussian parameters.");
            else
                flag_gauss = 1;
        }



        /* read phase statistics params */

        if (strcmp(param, "LOCPHSTAT") == 0) {
            if ((istat = GetNLLoc_PhaseStats(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading Phase Statistics parameters.");
            else
                flag_phstat = 1;
        }


        /* read take-off angles params */

        if (strcmp(param, "LOCANGLES") == 0) {
            if ((istat = GetNLLoc_Angles(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading Take-off Angles parameters.");
            else
                flag_angles = 1;
        }


        /* read magnitude calculation params */

        if (strcmp(param, "LOCMAG") == 0) {
            if ((istat = GetNLLoc_Magnitude(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading Magnitude Calculation parameters.");
            else
                flag_mag = 1;
        }


        /* read component params */

        if (strcmp(param, "LOCCMP") == 0) {
            if ((istat = GetCompDesc(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading Component description parameters.");
            else
                flag_comp = 1;
        }


        /* read alias params */

        if (strcmp(param, "LOCALIAS") == 0) {
            if ((istat = GetLocAlias(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading Alias parameters.");
            else
                flag_alias = 1;
        }


        /* read exclude params */

        if (strcmp(param, "LOCEXCLUDE") == 0) {
            if ((istat = GetLocExclude(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading Exclude parameters.");
            else
                flag_exclude = 1;
        }


        /* read exclude params */

        if (strcmp(param, "LOCINCLUDE") == 0) {
            if ((istat = GetLocInclude(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading Include parameters.");
            else
                flag_include = 1;
        }


        /* read station time delay params */

        if (strcmp(param, "LOCDELAY") == 0) {
            if ((istat = GetTimeDelays(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading Time Delay parameters.");
            else
                flag_time_delay = 1;
        }


        /* read surface time delay params */

        if (strcmp(param, "LOCTOPO_SURFACE") == 0) {
            if ((istat = GetTopoSurface(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading Topo Surface parameters.");
            else
                flag_topo_surface = 1;
        }


        /* read surface time delay params */

        if (strcmp(param, "LOCDELAY_SURFACE") == 0) {
            if ((istat = GetTimeDelaySurface(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading Time Delay Surface parameters.");
            else
                flag_time_delay_surface = 1;
        }


        /* read station elevation correction params */

        if (strcmp(param, "LOCELEVCORR") == 0) {
            if ((istat = GetElevCorr(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading Elevation Correction parameters.");
            else
                flag_elev_corr = 1;
        }


        /* read source params */

        if (strcmp(param, "LOCSRCE") == 0 || strcmp(param, "GTSRCE") == 0) {
            if ((istat = GetNextSource(strchr(line, ' '))) < 0) {
                nll_puterr("ERROR: reading source params:");
                nll_puterr(line);
            } else
                flag_source = 1;
        }



        /* unrecognized input */

        if (istat < 0) {
            if ((pchr = strchr(line, '\n')) != NULL)
                *pchr = '\0';
            snprintf(MsgStr, sizeof(MsgStr), "Skipping input: %s", line);
            nll_putmsg(5, MsgStr);
        }

    }


    // test for invalid configuration
    int test_search_pdf = 1;
    if (iUseSearchPrior || iUseSearchPosterior) {
        for (int ngrid = 0; ngrid < NumLocGrids; ngrid++) {
            if (LocGrid[ngrid].type != GRID_PROB_DENSITY) {
                nll_puterr("ERROR: cannot define NLLoc search PDF (LOCPRIOR or LOCPOSTERIOR) if any location grid (LOCGRID) is not type GRID_PROB_DENSITY.");
                test_search_pdf = 0;
            }

        }
    }

    /* check for missing required input */

    if (!flag_control)
        nll_puterr("ERROR: reading control (CONTROL) params.");
    if (!flag_outfile)
        nll_puterr("ERROR: reading i/o file (LOCFILES) params.");
    if (!flag_trans)
        nll_puterr("ERROR: reading transformation (TRANS) params.");
    if (!flag_grid)
        nll_puterr("ERROR: reading grid (LOCGRID) params.");
    if (!flag_search)
        nll_puterr("ERROR: reading search type (LOCSEARCH) params.");
    if (!flag_method)
        nll_puterr("ERROR: reading method (LOCMETH) params.");
    if (!flag_gauss)
        nll_puterr("ERROR: reading Gaussian (LOCGAU) params.");
    if (!flag_qual2err)
        nll_puterr("ERROR: reading Quality2Error (LOCQUAL2ERR) params.");


    /* check for missing optional input */

    if (!flag_gauss2) {
        snprintf(MsgStr, sizeof(MsgStr), "INFO: no Gaussian2 (LOCGAU2) params read.");
        nll_putmsg(2, MsgStr);
    }
    if (!flag_comment) {
        snprintf(MsgStr, sizeof(MsgStr), "INFO: no comment (LOCCOM) params read.");
        nll_putmsg(2, MsgStr);
        Hypocenter.comment[0] = '\0';
    }
    if (!flag_signature) {
        snprintf(MsgStr, sizeof(MsgStr), "INFO: no signature (LOCSIG) params read.");
        nll_putmsg(2, MsgStr);
        LocSignature[0] = '\0';
    }
    if (!flag_hyptype) {
        snprintf(MsgStr, sizeof(MsgStr),
                "INFO: no hypocenter output file type (LOCHYPOUT) params read.");
        nll_putmsg(2, MsgStr);
        snprintf(MsgStr, sizeof(MsgStr),
                "INFO: DEFAULT: \"LOCHYPOUT SAVE_NLLOC_ALL SAVE_HYPOINVERSE_Y2000_ARC SAVE_FMAMP\"");
        nll_putmsg(2, MsgStr);
        iSaveNLLocEvent = iSaveNLLocSum = 1;
        iSaveNLLocSumCSV = 1; // 20240826 AJL - added
        iSaveNLLocExpectation = 0;
        iSaveHypo71Sum = iSaveHypoEllSum = 0;
        iSaveHypo71Event = iSaveHypoEllEvent = 0;
        iSaveHypoInvSum = 0;
        iSaveHypoInvY2KArc = 1;
        iSaveAlberto4Sum = 0;
        iSaveFmamp = 1;
        iSaveSnapSum = 0;
        iCalcSedOrigin = 0;
        iSaveDecSec = 0;
        iSavePublicID = 0; // 20211208 AJL - added
        iSaveNone = 0;
    }
    if (!flag_phase_id) {
        snprintf(MsgStr, sizeof(MsgStr), "INFO: no phase identifier (LOCPHASEID) values read.");
        nll_putmsg(2, MsgStr);
    }
    if (!flag_sta_wt) {
        snprintf(MsgStr, sizeof(MsgStr), "INFO: no station distance weighting (LOCSTAWT) values read.");
        nll_putmsg(2, MsgStr);
    }
    if (!flag_prior) {
        snprintf(MsgStr, sizeof(MsgStr), "INFO: no Search Prior (LOCPRIOR) values read.");
        nll_putmsg(2, MsgStr);
    }
    if (!flag_posterior) {
        snprintf(MsgStr, sizeof(MsgStr), "INFO: no Search Prior (LOCPOSTERIOR) values read.");
        nll_putmsg(2, MsgStr);
    }
    if (!flag_mag) {
        snprintf(MsgStr, sizeof(MsgStr), "INFO: no Magnitude Calculation (LOCMAG) params read.");
        nll_putmsg(2, MsgStr);
    }
    if (!flag_phstat) {
        snprintf(MsgStr, sizeof(MsgStr),
                "INFO: no PhaseStatistics (LOCPHSTAT) params read.");
        nll_putmsg(2, MsgStr);
        RMS_Max = VERY_LARGE_DOUBLE;
        NRdgs_Min = -1;
        Gap_Max = VERY_LARGE_DOUBLE;
        P_ResidualMax = VERY_LARGE_DOUBLE;
        S_ResidualMax = VERY_LARGE_DOUBLE;
        Ell_Len3_Max = VERY_LARGE_DOUBLE;
        Hypo_Depth_Min = -VERY_LARGE_DOUBLE;
        Hypo_Depth_Max = VERY_LARGE_DOUBLE;
        Hypo_Dist_Max = VERY_LARGE_DOUBLE; // 20190812 AJL - bug fix
    }
    if (!flag_comp) {
        snprintf(MsgStr, sizeof(MsgStr), "INFO: no Component Descirption (LOCCMP) params read.");
        nll_putmsg(2, MsgStr);
    }
    if (!flag_alias) {
        snprintf(MsgStr, sizeof(MsgStr), "INFO: no Alias (LOCALIAS) params read.");
        nll_putmsg(2, MsgStr);
    }
    if (!flag_exclude) {
        snprintf(MsgStr, sizeof(MsgStr), "INFO: no Exclude (LOCEXCLUDE) params read.");
        nll_putmsg(2, MsgStr);
    }
    if (!flag_include) {
        snprintf(MsgStr, sizeof(MsgStr), "INFO: no Include (LOCINCLUDE) params read.");
        nll_putmsg(2, MsgStr);
    }
    if (!flag_time_delay) {
        snprintf(MsgStr, sizeof(MsgStr), "INFO: no Time Delay (LOCDELAY) params read.");
        nll_putmsg(2, MsgStr);
    }
    if (!flag_topo_surface) {
        snprintf(MsgStr, sizeof(MsgStr), "INFO: no Topo Surface (LOCTOPO_SURFACE) params read.");
        nll_putmsg(2, MsgStr);
    }
    if (!flag_time_delay_surface) {
        snprintf(MsgStr, sizeof(MsgStr), "INFO: no Time Delay Surface (LOCDELAY_SURFACE) params read.");
        nll_putmsg(2, MsgStr);
    }
    if (!flag_elev_corr) {
        snprintf(MsgStr, sizeof(MsgStr), "INFO: no Elevation Correction (LOCELEVCORR) params read.");
        nll_putmsg(2, MsgStr);
    }
    if (!flag_otime) {
        snprintf(MsgStr, sizeof(MsgStr), "INFO: no Fixed Origin Time (LOCFIXOTIME) params read.");
        nll_putmsg(2, MsgStr);
    }
    if (!flag_angles) {
        snprintf(MsgStr, sizeof(MsgStr), "INFO: no Take-off Angles (LOCANGLES) params read,");
        nll_putmsg(2, MsgStr);
        snprintf(MsgStr, sizeof(MsgStr), "      default is angleMode=ANGLES_NO, qualtiyMin=5 .");
        nll_putmsg(2, MsgStr);
        angleMode = ANGLE_MODE_NO;
        iAngleQualityMin = 5;
    }
    if (!flag_source) {

        snprintf(MsgStr, sizeof(MsgStr), "INFO: no Station (LOCSRCE or GTSRCE) params read.");
        nll_putmsg(2, MsgStr);
    }



    return (test_search_pdf * ok_search_pdf *
            flag_include_file * flag_control * flag_outfile * flag_grid *
            flag_search *
            flag_method * flag_gauss * flag_qual2err * flag_trans - 1);
}

/** function to read output file name
 *
 * NOTE: if the format of this control statement is changed, also update in Loc2ssst.c->GetNLLoc_Files()
 *
 */

int GetNLLoc_Files(char* line1) {

    int istat, nObsFile;
    char fnobs[FILENAME_MAX];
    char fn_time_grids_list[MAX_NUM_TIME_GRID_PATHS * (FILENAME_MAX + 2)]; // 20251027 add support for alternative travel-time grid path/root

    istat = sscanf(line1, "%s %s %s %s %d", fnobs, ftype_obs, fn_time_grids_list,
            fn_path_output, &iSwapBytesOnInput);
    if (istat < 5)
        iSwapBytesOnInput = 0;

    //printf("TEST!!! --> line1: %s\n", line1);
    //printf("TEST!!! --> fn_path_output: %s\n", fn_path_output);

    /* check for wildcards in observation file name */
    NumObsFiles = ExpandWildCards(fnobs, fn_loc_obs, MAX_NUM_OBS_FILES);

    if (message_flag >= 3) {
        snprintf(MsgStr, sizeof(MsgStr),
                "LOCFILES:  ObsType: %s  InGrids: %s.*  OutPut: %s.* iSwapBytesOnInput: %d",
                ftype_obs, fn_time_grids_list, fn_path_output, iSwapBytesOnInput);
        nll_putmsg(3, MsgStr);
        for (nObsFile = 0; nObsFile < NumObsFiles; nObsFile++) {
            snprintf(MsgStr, sizeof (MsgStr), "   Obs File: %3d  %s", nObsFile, fn_loc_obs[nObsFile]);
            nll_putmsg(3, MsgStr);
        }
    }

    if (NumObsFiles == MAX_NUM_OBS_FILES)
        nll_putmsg(1, "LOCFILES: WARNING: maximum number of files/events reached");

    // parse time grid path/roots
    NumTimeGridPaths = 0;
    char* token = strtok(fn_time_grids_list, ",");
    while (token != NULL) {
        strncpy(fn_time_grids[NumTimeGridPaths++], token, FILENAME_MAX);
        token = strtok(NULL, ",");
    }

    return (0);
}

/** function to read search type ***/

int GetNLLoc_SearchType(char* line1) {
    int istat, ierr;

    char search_type[MAXLINE];


    istat = sscanf(line1, "%s", search_type);

    if (istat != 1)
        return (-1);

    if (strcmp(search_type, "GRID") == 0) {

        SearchType = SEARCH_GRID;
        istat = sscanf(line1, "%s %d", search_type, &(Scatter.npts));
        if (istat != 2)
            return (-1);

        snprintf(MsgStr, sizeof(MsgStr), "LOCSEARCH:  Type: %s NumScatter %d",
                search_type, Scatter.npts);
        nll_putmsg(3, MsgStr);

    } else if (strcmp(search_type, "MET") == 0) {

        SearchType = SEARCH_MET;
        istat = sscanf(line1, "%s %d %d %d %d %d %lf %lf %lf %lf",
                search_type, &MetNumSamples, &MetLearn, &MetEquil,
                &MetStartSave, &MetSkip,
                &MetStepInit, &MetStepMin, &MetStepFact, &MetProbMin);
        ierr = 0;

        snprintf(MsgStr, sizeof(MsgStr),
                "LOCSEARCH:  Type: %s  numSamples %d  numLearn %d  numEquilibrate %d  startSave %d  numSkip %d  stepInit %lf  stepMin %lf  stepFact %lf  probMin %lf",
                search_type, MetNumSamples, MetLearn, MetEquil,
                MetStartSave, MetSkip,
                MetStepInit, MetStepMin, MetStepFact, MetProbMin);
        nll_putmsg(3, MsgStr);

        if (checkRangeInt("LOCSEARCH", "numSamples", MetNumSamples, 1, 0, 0, 0) != 0)
            ierr = -1;
        if (checkRangeInt("LOCSEARCH", "numLearn", MetLearn, 1, 0, 0, 0) != 0)
            ierr = -1;
        if (checkRangeInt("LOCSEARCH", "numEquilibrate", MetEquil, 1, 0, 0, 0) != 0)
            ierr = -1;
        if (checkRangeInt("LOCSEARCH", "startSave", MetStartSave, 1, 0, 0, 0) != 0)
            ierr = -1;
        if (checkRangeInt("LOCSEARCH", "numSkip", MetSkip, 1, 1, 0, 0) != 0)
            ierr = -1;
        if (checkRangeDouble("LOCSEARCH", "stepMin", MetStepMin, 1, 0.0, 0, 0.0) != 0)
            ierr = -1;
        if (ierr < 0)
            return (-1);
        if (istat != 10)
            return (-1);

        //?? AJL 17JAN2000 MetUse = MetNumSamples - MetEquil;
        MetUse = MetNumSamples - MetStartSave;

        /* check for "normal" StartSave value*/
        if (MetStartSave < MetLearn + MetEquil) {
            snprintf(MsgStr, sizeof(MsgStr),
                    "LOCSEARCH:  WARNING: Metropolis StartSave < NumLearn + NumEquilibrate.");
            nll_putmsg(1, MsgStr);
        }

    } else if (strcmp(search_type, "OCT") == 0) {

        SearchType = SEARCH_OCTTREE;
        istat = sscanf(line1, "%s %d %d %d %lf %d %d %d %d %lf",
                search_type, &octtreeParams.init_num_cells_x,
                &octtreeParams.init_num_cells_y, &octtreeParams.init_num_cells_z,
                &octtreeParams.min_node_size, &octtreeParams.max_num_nodes,
                &octtreeParams.num_scatter, &octtreeParams.use_stations_density,
                &octtreeParams.stop_on_min_node_size,
                &octtreeParams.mean_cell_velocity);

        if (istat < 8)
            octtreeParams.use_stations_density = 0;
        if (octtreeParams.use_stations_density < 0)
            octtreeParams.use_stations_density = 0;

        if (istat < 9)
            octtreeParams.stop_on_min_node_size = 1;

        if (istat < 10)
            octtreeParams.mean_cell_velocity = -1.0;


        snprintf(MsgStr, sizeof(MsgStr),
                "LOCSEARCH:  Type: %s  init_num_cells_x %d  init_num_cells_y %d  init_num_cells_z %d  min_node_size %f  max_num_nodes %d  num_scatter %d  use_stations_density %d  stop_on_min_node_size %d  octtreeParams.mean_cell_velocity %f",
                search_type, octtreeParams.init_num_cells_x, octtreeParams.init_num_cells_y,
                octtreeParams.init_num_cells_z,
                octtreeParams.min_node_size, octtreeParams.max_num_nodes,
                octtreeParams.num_scatter, octtreeParams.use_stations_density,
                octtreeParams.stop_on_min_node_size, octtreeParams.mean_cell_velocity);
        nll_putmsg(3, MsgStr);


        // check for valid input values
        ierr = 0;
        if (checkRangeInt("LOCSEARCH", "init_num_cells_x",
                octtreeParams.init_num_cells_x, 1, 0, 0, 0) != 0)
            ierr = -1;
        if (checkRangeInt("LOCSEARCH", "init_num_cells_y",
                octtreeParams.init_num_cells_y, 1, 0, 0, 0) != 0)
            ierr = -1;
        if (checkRangeInt("LOCSEARCH", "init_num_cells_z",
                octtreeParams.init_num_cells_z, 1, 0, 0, 0) != 0)
            ierr = -1;
        if (checkRangeDouble("LOCSEARCH", "min_node_size",
                octtreeParams.min_node_size, 1, 0.0, 0, 0.0) != 0)
            ierr = -1;
        if (checkRangeInt("LOCSEARCH", "max_num_nodes",
                octtreeParams.max_num_nodes, 1, 0, 0, 0) != 0)
            ierr = -1;
        if (checkRangeInt("LOCSEARCH", "num_scatter",
                octtreeParams.num_scatter, 1, 0, 0, 0) != 0)
            ierr = -1;

        // check for valid OctTree values
        int init_n_cells = octtreeParams.init_num_cells_x * octtreeParams.init_num_cells_y * octtreeParams.init_num_cells_z;
        if (init_n_cells >= octtreeParams.max_num_nodes) {
            snprintf(MsgStr, sizeof(MsgStr), "ERROR: LOCSEARCH OCT: OctTree init_num_cells (%d) >= max_num_nodes (%d): no oct-tree subdivision can be performed.",
                    init_n_cells, octtreeParams.max_num_nodes);
            nll_putmsg(1, MsgStr);
            ierr = -1;
        } else if (octtreeParams.max_num_nodes - init_n_cells < 10000) {
            snprintf(MsgStr, sizeof(MsgStr), "WARNING: LOCSEARCH OCT: OctTree max_num_nodes - init_num_cells (%d) < 10000: very few oct-tree subdivisions can be performed.",
                    octtreeParams.max_num_nodes - init_n_cells);
            nll_putmsg(1, MsgStr);
        }

        if (ierr < 0)
            return (-1);
        if (istat < 7)

            return (-1);

    }


    return (0);
}

/** function to read search prior parameters
 *  20190510 AJL - added
 **/

int GetNLLoc_PdfGrid(char* line1, int prior_type) {

    int istat, ierr;

    SearchPdfGridDesc *searchPdfGrid = NULL;
    if (prior_type == PDF_GRID_PRIOR) {
        searchPdfGrid = &SearchPrior;
    } else if (prior_type == PDF_GRID_POSTERIOR) {
        searchPdfGrid = &SearchPosterior;
    }

    char grid_type[MAXLINE];
    static char file_line[MAXLINE_LONG];

    istat = sscanf(line1, "%s", grid_type);

    if (istat != 1)
        return (-1);

    if (strcmp(grid_type, "OCT_TREE") == 0) {

        searchPdfGrid->gridType = PDF_GRID_OCT_TREE;
        searchPdfGrid->max_total_other_weight = -1.0; // default
        searchPdfGrid->max_count_other = -1; // default   // 20211031 AJL - added
        searchPdfGrid->max_se3 = -1.0; // default
        searchPdfGrid->max_mag_diff = -1.0; // default
        searchPdfGrid->min_mag = -999; // default
        istat = sscanf(line1, "%*s %s %lf %lf %lf %lf %d %lf %lf",
                searchPdfGrid->grid_file_path, &(searchPdfGrid->default_value),
                &(searchPdfGrid->coherence_min), &(searchPdfGrid->max_total_other_weight), &(searchPdfGrid->max_se3), &(searchPdfGrid->max_count_other),
                &(searchPdfGrid->max_mag_diff), &(searchPdfGrid->min_mag));
        snprintf(MsgStr, sizeof(MsgStr), "LOCPRIOR/LOCPOSTERIOR:  Type: %s  GridFile: %s  DefaultValue: %e  CoherenceMin: %f  MaxOtherWeight: %f  MaxSE3: %f  MaxNother: %d  MaxMagDiff: %f  MinMag: %f",
                grid_type, searchPdfGrid->grid_file_path, searchPdfGrid->default_value,
                searchPdfGrid->coherence_min, searchPdfGrid->max_total_other_weight, searchPdfGrid->max_se3, searchPdfGrid->max_count_other,
                searchPdfGrid->max_mag_diff, searchPdfGrid->min_mag);
        nll_putmsg(3, MsgStr);
        snprintf(MsgStr, sizeof(MsgStr), "Dummy Message"); // 20210527 AJL - Bug Fix: TODO: somewhere this message is printed!
        ierr = 0;
        if (checkRangeDouble("LOCPRIOR/LOCPOSTERIOR", "DefaultValue", searchPdfGrid->default_value, 1, 0.0, 0, 0.0) != 0)
            ierr = -1;
        if (ierr < 0 || istat < 3)
            return (-1);

        // read oct-tree grids
        // check for wildcards in observation file name
        static char fn_pdf_grid[MAX_NUM_PDF_GRID_FILES][FILENAME_MAX];
        static double coherence[MAX_NUM_PDF_GRID_FILES];
        int numPdfGridFiles = ExpandWildCards(searchPdfGrid->grid_file_path, fn_pdf_grid, MAX_NUM_PDF_GRID_FILES);
        if (numPdfGridFiles >= MAX_NUM_PDF_GRID_FILES) {
            snprintf(MsgStr, sizeof(MsgStr), "WARNING: maximum number of pdf grid files files exceeded, only first %d will be processed.", MAX_NUM_PDF_GRID_FILES);
            nll_puterr(MsgStr);
        }
        // if single file, check if is *.stream_coherences file and read file names and coherence
        int found_valid_stream_coherences = 0;
        if (numPdfGridFiles == 1) {
            FILE *fp_coherence_test;
            if ((fp_coherence_test = fopen(fn_pdf_grid[0], "r")) != NULL) {
                // test if stream_coherences file
                if (fscanf(fp_coherence_test, "%s", file_line) > 0 && strcmp(file_line, "STREAM_COHERENCES") == 0) {
                    // second line is self location
                    double file_coherence = 1.0;
                    // self line may or may not have coherence
                    if (fscanf(fp_coherence_test, "%lf %s", &file_coherence, file_line) == 2) {
                        found_valid_stream_coherences = 1;
                        // with coherence value
                    } else if (fscanf(fp_coherence_test, "%s", file_line) > 0) {
                        // no coherence value
                        found_valid_stream_coherences = 1;
                    } else {
                        found_valid_stream_coherences = 0;
                    }
                    // check min_mag
                    static HypoDesc hypo_self;
                    if (ReadHypoDesc(file_line, &hypo_self) < -1) {
                        nll_puterr2("ERROR: opening or reading self event hypo file", file_line);
                    }
                    //printf("DEBUG: searchPdfGrid->min_mag %f hypo_self.amp_mag %f\n", searchPdfGrid->min_mag, hypo_self.amp_mag);
                    if (searchPdfGrid->min_mag > -99 && fabs(hypo_self.amp_mag - MAGNITUDE_NULL) > 0.01) {
                        if (hypo_self.amp_mag < searchPdfGrid->min_mag) {
                            snprintf(MsgStr, sizeof(MsgStr), "INFO: GetNLLoc_PdfGrid: amp_mag: %f < searchPdfGrid->min_mag %f  %s : EVENT IGNORED",
                                    hypo_self.amp_mag, searchPdfGrid->min_mag, file_line);
                            nll_putmsg(1, MsgStr);
                            return (-9);
                        }
                    }
                    if (found_valid_stream_coherences) {
                        numPdfGridFiles = 0;
                        // include self if posterior
                        if (prior_type == PDF_GRID_POSTERIOR && file_coherence >= searchPdfGrid->coherence_min) {
                            snprintf(fn_pdf_grid[numPdfGridFiles], sizeof(fn_pdf_grid[numPdfGridFiles]), "%s", file_line);
                            strncat(fn_pdf_grid[numPdfGridFiles], ".octree", sizeof(fn_pdf_grid[numPdfGridFiles]) - strlen(fn_pdf_grid[numPdfGridFiles]) - 1);
                            coherence[numPdfGridFiles] = file_coherence;
                            numPdfGridFiles++;
                        }
                        // next lines are coherence and oct-tree file root for each child
                        while (fscanf(fp_coherence_test, "%lf %s", &(coherence[numPdfGridFiles]), file_line) > 1) {
                            // check hypo filters
                            static HypoDesc hypo_other;
                            if (searchPdfGrid->max_se3 > 0.0 || searchPdfGrid->max_mag_diff > 0.0) {
                                if (ReadHypoDesc(file_line, &hypo_other) < -1) {
                                    nll_puterr2("ERROR: opening or reading other event hypo file", file_line);
                                }
                                // check se3
                                if (searchPdfGrid->max_se3 > 0.0) {
                                    if (hypo_other.ellipsoid.len3 > searchPdfGrid->max_se3) {
                                        snprintf(MsgStr, sizeof(MsgStr), "INFO: GetNLLoc_PdfGrid: se3: %f > searchPdfGrid->max_se3 %f  %s : IGNORED OTHER",
                                                hypo_other.ellipsoid.len3, searchPdfGrid->max_se3, file_line);
                                        nll_putmsg(3, MsgStr);
                                        continue;
                                    }
                                }
                                // check max_mag_diff
                                //printf("DEBUG: searchPdfGrid->max_mag_diff %f hypo_self.amp_mag %f\n", searchPdfGrid->max_mag_diff, hypo_self.amp_mag);
                                if (searchPdfGrid->max_mag_diff > 0.0
                                        && fabs(hypo_other.amp_mag - MAGNITUDE_NULL) > 0.01 && fabs(hypo_self.amp_mag - MAGNITUDE_NULL) > 0.01) {
                                    if (fabs(hypo_other.amp_mag - hypo_self.amp_mag) > searchPdfGrid->max_mag_diff) {
                                        snprintf(MsgStr, sizeof(MsgStr), "INFO: GetNLLoc_PdfGrid: amp_mag diff: %f-%f > searchPdfGrid->max_mag_diff %f  %s : IGNORED OTHER",
                                                hypo_self.amp_mag, hypo_other.amp_mag, searchPdfGrid->min_mag, file_line);
                                        nll_putmsg(3, MsgStr);
                                        continue;
                                    }
                                }
                            }
                            // check coherence min
                            if (coherence[numPdfGridFiles] >= searchPdfGrid->coherence_min) {
                                snprintf(fn_pdf_grid[numPdfGridFiles], sizeof(fn_pdf_grid[numPdfGridFiles]), "%s", file_line);
                                strncat(fn_pdf_grid[numPdfGridFiles], ".octree", sizeof(fn_pdf_grid[numPdfGridFiles]) - strlen(fn_pdf_grid[numPdfGridFiles]) - 1);
                                numPdfGridFiles++;
                                if (numPdfGridFiles >= MAX_NUM_PDF_GRID_FILES) {
                                    snprintf(MsgStr, sizeof(MsgStr),
                                            "WARNING: maximum number of coherence pdf grid files files reached, only first %d will be processed.",
                                            MAX_NUM_PDF_GRID_FILES);
                                    nll_putmsg(1, MsgStr);
                                    break;
                                }
                                if (searchPdfGrid->max_count_other >= 0
                                        && numPdfGridFiles > searchPdfGrid->max_count_other) {
                                    snprintf(MsgStr, sizeof(MsgStr),
                                            "WARNING: maximum number of other coherence pdf grid files files reached, only first %d will be processed.",
                                            searchPdfGrid->max_count_other);
                                    nll_putmsg(1, MsgStr);
                                    break;
                                }
                            }
                        }
                    }
                }
            } else {
                // file does not exist, do normal location
                if (message_flag >= 1)
                    fprintf(stdout, "INFO: Ignoring LOCPRIOR/LOCPOSTERIOR: File does not exist: %s\n", fn_pdf_grid[0]);
                return (0);
            }
            fclose(fp_coherence_test);
        }
        if (!found_valid_stream_coherences) {
            for (int nFile = 0; nFile < numPdfGridFiles; nFile++) {
                coherence[nFile] = 1.0;
            }
        }
        // allocate oct-tree grids
        if ((searchPdfGrid->tree3D = (Tree3D **) malloc(numPdfGridFiles * sizeof (Tree3D*))) == NULL) {
            nll_puterr("ERROR: allocating memory for search PDF oct-tree grid.");
            return (-1);
        }
        // allocate coherence
        if ((searchPdfGrid->coherence = (double *) malloc(numPdfGridFiles * sizeof (double))) == NULL) {
            nll_puterr("ERROR: allocating memory for search PDF coherence.");
            return (-1);
        }
        // allocate weight
        if ((searchPdfGrid->weight = (double *) malloc(numPdfGridFiles * sizeof (double))) == NULL) {
            nll_puterr("ERROR: allocating memory for search PDF weight.");
            return (-1);
        }
        // allocate first_motion_arrivals
        if ((searchPdfGrid->first_motion_arrivals = (ArrivalDesc **) malloc(numPdfGridFiles * sizeof (ArrivalDesc *))) == NULL) {
            nll_puterr("ERROR: allocating memory for search PDF first_motion_arrivals array.");
            return (-1);
        }
        // allocate nfirst_motion_arrivals
        if ((searchPdfGrid->nfirst_motion_arrivals = (int *) malloc(numPdfGridFiles * sizeof (int))) == NULL) {
            nll_puterr("ERROR: allocating memory for search PDF nfirst_motion_arrivals.");
            return (-1);
        }
        // read grid files
        FILE *fp_oct_in;
        searchPdfGrid->nGrids = 0;
        double tot_other_wt = 0.0; // 20200727 AJL - added
        // temporary arrival array
        ArrivalDesc* arrival_tmp;
        if ((arrival_tmp = (ArrivalDesc *) calloc(MAX_NUM_ARRIVALS, sizeof (ArrivalDesc))) == NULL) {
            nll_puterr("ERROR: allocating memory for search PDF temporary first_motion_arrivals.");
            return (-1);
        }
        for (int nFile = 0; nFile < numPdfGridFiles; nFile++) {
            // open input grid file
            if ((fp_oct_in = fopen(fn_pdf_grid[nFile], "r")) == NULL) {
                nll_puterr2("ERROR: opening input oct tree file", fn_pdf_grid[nFile]);
                return (-1);
            }
            searchPdfGrid->tree3D[nFile] = readTree3D(fp_oct_in);
            searchPdfGrid->coherence[nFile] = coherence[nFile];
            // weight is zero at coherence_min and 1.0 at coherence=1.0
            //searchPdfGrid->weight[nFile] = (searchPdfGrid->coherence[nFile] - searchPdfGrid->coherence_min) / (1.0 - searchPdfGrid->coherence_min);
            // TEST 20210126 AJL - 0 -> 1 cosine taper weighting
            //if (1) {
            // weight is zero at coherence_min and 1.0 at coherence=0.9
            double wt_tmp = (searchPdfGrid->coherence[nFile] - searchPdfGrid->coherence_min) / (0.9 - searchPdfGrid->coherence_min);
            if (wt_tmp >= 1.0) {
                wt_tmp = 1.0;
            } else if (wt_tmp <= 0.0) {
                wt_tmp = 0.0;
            } else {
                //printf("DEBUG: read input oct tree file: coherence: %f  weight %f  %s", coherence[nFile], searchPdfGrid->weight[nFile], fn_pdf_grid[nFile]);
                // use cos instead of sin   wt_tmp = cPI * (wt_tmp - 0.5); // -PI/2 -> PI/2
                wt_tmp = cPI * (1.0 - wt_tmp); // PI -> 0
                //printf(" -> %f", wt_tmp);
                // use cos instead of sin   wt_tmp = 0.5 * (sin(wt_tmp) + 1.0); // 0 -> 1 sine
                wt_tmp = 0.5 * cos(wt_tmp) + 0.5; // 0 -> 1 cos
                //printf(" -> %f", wt_tmp);
                //printf("\n");
            }
            searchPdfGrid->weight[nFile] = wt_tmp;
            //}
            // END TEST
            // TEST 20210110 AJL - 0 -> 1 sine weighting
            if (0) {
                //printf("DEBUG: read input oct tree file: coherence: %f  weight %f  %s", coherence[nFile], searchPdfGrid->weight[nFile], fn_pdf_grid[nFile]);
                double wt_tmp = cPI * (searchPdfGrid->weight[nFile] - 0.5); // -PI/2 -> PI/2
                //printf(" -> %f", wt_tmp);
                wt_tmp = 0.5 * (sin(wt_tmp) + 1.0); // 0 -> 1 sine
                //printf(" -> %f", wt_tmp);
                // sqrt(sin)
                if (1) {
                    // =IF(D2>=0.5,0.5+0.5*POWER(2*(D2-0.5),1/2),0.5-0.5*POWER(-2*(D2-0.5),1/2))
                    if (wt_tmp >= 0.5) {
                        wt_tmp = 0.5 + 0.5 * sqrt(2.0 * (wt_tmp - 0.5));
                    } else {
                        wt_tmp = 0.5 - 0.5 * sqrt(-2.0 * (wt_tmp - 0.5));
                    }
                    //printf(" -> %f", wt_tmp);
                }
                //printf("\n");
                searchPdfGrid->weight[nFile] = wt_tmp;
            }
            // END TEST
            // Limit total weight of other events  // 20200727 AJL - added
            if (nFile > 0) {
                tot_other_wt += searchPdfGrid->weight[nFile];
            }
            fclose(fp_oct_in);
            searchPdfGrid->nGrids++;
            // read arrivals with first motions
            char fn_hypo_root[FILENAME_MAX];
            snprintf(fn_hypo_root, sizeof(fn_hypo_root), "%s", fn_pdf_grid[nFile]);
            *(strrchr(fn_hypo_root, '.')) = '\0'; // get hypo root name
            FILE *fpio_tmp = NULL;
            ReadFirstMotionArrivals(&fpio_tmp, fn_hypo_root, arrival_tmp, &(searchPdfGrid->nfirst_motion_arrivals[nFile]));
            // allocate and load arrivals to first_motion_arrivals
            if ((searchPdfGrid->first_motion_arrivals[nFile]
                    = (ArrivalDesc *) malloc(searchPdfGrid->nfirst_motion_arrivals[nFile] * sizeof (ArrivalDesc))) == NULL) {
                nll_puterr("ERROR: allocating memory for search PDF first_motion_arrivals.");
                return (-1);
            }
            for (int narr = 0; narr < searchPdfGrid->nfirst_motion_arrivals[nFile]; narr++) {
                searchPdfGrid->first_motion_arrivals[nFile][narr] = arrival_tmp[narr];
            }
        }
        // Limit total weight of other events  // 20200727 AJL - added
        if (searchPdfGrid->max_total_other_weight > 0.0 && tot_other_wt > searchPdfGrid->max_total_other_weight) {
            for (int nFile = 1; nFile < numPdfGridFiles; nFile++) {
                //searchPdfGrid->weight[nFile] *= searchPdfGrid->max_total_other_weight / tot_other_wt;
                // 20231114 AJL - Bug fix: //searchPdfGrid->weight[nFile] /= tot_other_wt;
                searchPdfGrid->weight[nFile] *= searchPdfGrid->max_total_other_weight / tot_other_wt; // 20231114 AJL - Bug fix
            }
        }
        free(arrival_tmp);

    } else if (strcmp(grid_type, "GRID") == 0) {

        searchPdfGrid->gridType = PDF_GRID_GRID;
        int iswap_bytes;
        istat = sscanf(line1, "%*s %s %d %lf",
                searchPdfGrid->grid_file_path, &iswap_bytes, &(searchPdfGrid->default_value));
        snprintf(MsgStr, sizeof(MsgStr), "LOCPRIOR/LOCPOSTERIOR:  Type: %s  GridFile: %s  SwapBytes: %d  DefaultValue: %e",
                grid_type, searchPdfGrid->grid_file_path, iswap_bytes, searchPdfGrid->default_value);
        nll_putmsg(3, MsgStr);
        ierr = 0;
        if (checkRangeDouble("LOCPRIOR/LOCPOSTERIOR", "DefaultValue", searchPdfGrid->default_value, 1, 0.0, 0, 0.0) != 0)
            ierr = -1;
        if (ierr < 0 || istat < 3)
            return (-1);

        // read and initialize grid
        // open grid file and read header
        FILE * fp_prior_grid, *fp_prior_hdr;
        if ((istat = OpenGrid3dFile(searchPdfGrid->grid_file_path, &fp_prior_grid, &fp_prior_hdr,
                &searchPdfGrid->grid, " ", NULL, iswap_bytes)) < 0) {
            CloseGrid3dFile(&searchPdfGrid->grid, &fp_prior_grid, &fp_prior_hdr);
            nll_puterr2("ERROR: cannot open PDF grid", searchPdfGrid->grid_file_path);
            return (EXIT_ERROR_FILEIO);
        }
        if (message_flag >= 3)
            display_grid_param(&searchPdfGrid->grid);
        // allocate grids
        searchPdfGrid->grid.buffer = AllocateGrid(&searchPdfGrid->grid);
        if (searchPdfGrid->grid.buffer == NULL) {
            nll_puterr(
                    "ERROR: allocating memory for search PDF grid buffer.");
            return (EXIT_ERROR_MEMORY);
        }
        // create grid array access pointers
        searchPdfGrid->grid.array = CreateGridArray(&searchPdfGrid->grid);
        if (searchPdfGrid->grid.array == NULL) {
            nll_puterr(
                    "ERROR: creating array for accessing search PDF grid buffer.");
            return (EXIT_ERROR_MEMORY);
        }
        // read grid
        if ((istat =
                ReadGrid3dBuf(&searchPdfGrid->grid, fp_prior_grid)) < 0) {
            nll_puterr("ERROR: reading search PDF grid from disk.");
            return (EXIT_ERROR_FILEIO);
        }
        CloseGrid3dFile(&searchPdfGrid->grid, &fp_prior_grid, &fp_prior_hdr);

    } else {
        searchPdfGrid->gridType = PDF_GRID_UNDEF;
        nll_puterr2("ERROR: unrecognized search PDF grid type:", grid_type);
        return (-1);
    }

    if (prior_type == PDF_GRID_PRIOR) {
        iUseSearchPrior = 1;
    } else if (prior_type == PDF_GRID_POSTERIOR) {

        iUseSearchPosterior = 1;
        /* 20220107 AJL - Revert Bug fix: faster to put grids in memory, maybe too much disk I/O otherwise
        // 20211026 AJL - Bug fix: do not put 3D grids in memory: LOCMETH maximum_number_3D_grids, not needed and can use much memory
        if (MaxNum3DGridMemory != 0) {
            MaxNum3DGridMemory = 0;
            snprintf(MsgStr, sizeof(MsgStr), "INFO: LOCPOSTERIOR is active: LOCMETH maximum_number_3D_grids reset to 0");
            nll_putmsg(1, MsgStr);
        }
         */
    }

    return (0);
}

/** function to read requested hypocenter output file types ***/

int GetNLLoc_HypOutTypes(char* line1) {
    int istat;

    char *pchr, hyp_type[MAXLINE];


    snprintf(MsgStr, sizeof(MsgStr), "LOCHYPOUT:  ");

    pchr = line1;
    do {

        /* check for blank line */
        while (*pchr == ' ')
            pchr++;
        if (isspace(*pchr))
            break;

        if ((istat = sscanf(pchr, "%s", hyp_type)) != 1)
            return (-1);

        if (strcmp(hyp_type, "SAVE_NLLOC_ALL") == 0) {
            iSaveNLLocEvent = iSaveNLLocSum = 1;
            iSaveNLLocSumCSV = 1; // 20240826 AJL - added
        } else if (strcmp(hyp_type, "SAVE_NLLOC_SUM") == 0) {
            iSaveNLLocSum = 1;
            iSaveNLLocSumCSV = 1; // 20240826 AJL - added
        } else if (strcmp(hyp_type, "SAVE_NLLOC_EXPECTATION") == 0) // 20170811 AJL - added
            iSaveNLLocExpectation = 1;
        else if (strcmp(hyp_type, "SAVE_NLLOC_OCTREE") == 0)
            iSaveNLLocOctree = 1;
        else if (strcmp(hyp_type, "SAVE_NLLOC_JSON") == 0) // 20220131 AJL - added to support JSON output of location results
            iSaveNLLocEvent_JSON = 1;
        else if (strcmp(hyp_type, "SAVE_HYPO71_ALL") == 0)
            iSaveHypo71Event = iSaveHypo71Sum = 1;
        else if (strcmp(hyp_type, "SAVE_HYPO71_SUM") == 0)
            iSaveHypo71Sum = 1;
        else if (strcmp(hyp_type, "SAVE_HYPOELL_ALL") == 0)
            iSaveHypoEllEvent = iSaveHypoEllSum = 1;
        else if (strcmp(hyp_type, "SAVE_HYPOELL_SUM") == 0)
            iSaveHypoEllSum = 1;
        else if (strcmp(hyp_type, "SAVE_HYPOINV_SUM") == 0)
            iSaveHypoInvSum = 1;
        else if (strcmp(hyp_type, "SAVE_HYPOINVERSE_Y2000_ARC") == 0)
            iSaveHypoInvY2KArc = 1;
        else if (strcmp(hyp_type, "SAVE_ALBERTO_3D_4") == 0)
            iSaveAlberto4Sum = 1;
        else if (strcmp(hyp_type, "SAVE_FMAMP") == 0)
            iSaveFmamp = 1;
        else if (strcmp(hyp_type, "SAVE_SNAP_SUM") == 0)
            /* SH 02/26/2004 added SNAP summary file */
            iSaveSnapSum = 1;
            /* filename format, int sec or 5.2 decimal seconds */
            // 20100617 AJL -  added: calculate and report SED origin, e.g. SED location quality indicators
        else if (strcmp(hyp_type, "CALC_SED_ORIGIN") == 0)
            iCalcSedOrigin = 1;
            /* filename format, int sec or 5.2 decimal seconds */
        else if (strcmp(hyp_type, "FILENAME_DEC_SEC") == 0)
            iSaveDecSec = 1;
        else if (strcmp(hyp_type, "FILENAME_PUBLIC_ID") == 0) // 20211208 AJL - added
            iSavePublicID = 1;
            /* new values NLL PHASE_2 format*/
            /* 20060629 AJL - Added */
        else if (strcmp(hyp_type, "NLL_FORMAT_VER_2") == 0)
            PhaseFormat = FORMAT_PHASE_2;
        else if (strcmp(hyp_type, "NONE") == 0) {
            iSaveNone = 1;
            iSaveNLLocEvent = iSaveNLLocSum = iSaveHypo71Sum = iSaveHypoEllSum =
                    iSaveHypo71Event = iSaveHypoEllEvent = iSaveHypoInvSum = iSaveHypoInvY2KArc =
                    iSaveAlberto4Sum = iSaveFmamp = iSaveSnapSum = iCalcSedOrigin = iSaveDecSec = iSavePublicID = 0;
            iSaveNLLocExpectation = 0; // 20170811 AJL - added
            iSaveNLLocSumCSV = 0; // 20240826 AJL - added
        } else
            return (-1);

        strncat(MsgStr, hyp_type, sizeof(MsgStr) - strlen(MsgStr) - 1);
        strncat(MsgStr, " ", sizeof(MsgStr) - strlen(MsgStr) - 1);

    } while ((pchr = strchr(pchr + 1, ' ')) != NULL);

    nll_putmsg(3, MsgStr);

    return (0);
}

/** function to read method
 *
 * NOTE: if the format of this control statement is changed, also update in Loc2ssst.c->GetNLLoc_Method()
 *
 */

int GetNLLoc_Method(char* line1) {
    int istat, ierr;

    char loc_method[MAXLINE];


    istat = sscanf(line1, "%s %lf %d %d %d %lf %d %lf %d", loc_method,
            &DistStaGridMax, &MinNumArrLoc, &MaxNumArrLoc, &MinNumSArrLoc,
            &VpVsRatio, &MaxNum3DGridMemory, &DistStaGridMin, &iRejectDuplicateArrivals);
    if (istat < 8)
        DistStaGridMin = -1.0;
    if (istat < 9)
        iRejectDuplicateArrivals = 1;

    snprintf(MsgStr, sizeof(MsgStr),
            "LOCMETH:  method: %s  minDistStaGrid: %lf  maxDistStaGrid: %lf  minNumberPhases: %d  maxNumberPhases: %d  minNumberSphases: %d  VpVsRatio: %lf  max3DGridMemory: %d  DistStaGridMin: %f  iRejectDuplicateArrivals: %d",
            loc_method, DistStaGridMin, DistStaGridMax, MinNumArrLoc, MaxNumArrLoc,
            MinNumSArrLoc, VpVsRatio, MaxNum3DGridMemory, DistStaGridMin, iRejectDuplicateArrivals);
    nll_putmsg(3, MsgStr);

    /* 20220107 AJL - Revert Bug fix: faster to put grids in memory, maybe too much disk I/O otherwise
    // 20211026 AJL - Bug fix: do not put 3D grids in memory: LOCMETH maximum_number_3D_grids, not needed and can use much memory
    if (iUseSearchPosterior == 1) {
        // 20211026 AJL - Bug fix: do not put 3D grids in memory LOCMETH maximum_number_3D_grids, not needed and can use much memory
        if (MaxNum3DGridMemory != 0) {
            MaxNum3DGridMemory = 0;
            snprintf(MsgStr, sizeof(MsgStr), "INFO: LOCPOSTERIOR is active: LOCMETH maximum_number_3D_grids reset to 0");
            nll_putmsg(1, MsgStr);
        }
    }*/

    // 20170922 AJL - bug fix, since MaxNum3DGridMemory used in GridMemLib.c with values assumed >=0, just set very large value here if <0
    if (MaxNum3DGridMemory < 0) {
        MaxNum3DGridMemory = INT_MAX;
    }

    ierr = 0;

    if (ierr < 0 || istat < 7)
        return (-1);

    EDT_use_otime_weight = 0;

    if (strcmp(loc_method, "GAU_ANALYTIC") == 0) {
        LocMethod = METH_GAU_ANALYTIC;
    } else if (strcmp(loc_method, "GAU_TEST") == 0) {
        LocMethod = METH_GAU_TEST;
    } else if (strcmp(loc_method, "OT_STACK") == 0) {
        LocMethod = METH_OT_STACK;
    } else if (strcmp(loc_method, "ML_OT") == 0) {
        LocMethod = METH_ML_OT;
        EDT_use_otime_weight = 2;
    } else if (strcmp(loc_method, "EDT") == 0 || strcmp(loc_method, "EDT_TEST") == 0) {
        LocMethod = METH_EDT;
    } else if (strcmp(loc_method, "EDT_OT_WT") == 0) {
        LocMethod = METH_EDT;
        EDT_use_otime_weight = 1;
    } else if (strcmp(loc_method, "EDT_OT_WT_ML") == 0) {
        LocMethod = METH_EDT;
        EDT_use_otime_weight = 2;
    } else if (strcmp(loc_method, "EDT_BOX") == 0) {
        LocMethod = METH_EDT_BOX;
    } else if (strcmp(loc_method, "L1_NORM") == 0) { // 20140515 AJL - added for NLDiffLoc    // 20150324 AJL - added for NLLoc
        LocMethod = METH_L1_NORM;
    } else {
        LocMethod = METH_UNDEF;
        nll_puterr2("ERROR: unrecognized location method:", loc_method);
        return (EXIT_ERROR_LOCATE);
    }

    if (MaxNumArrLoc < 1)
        MaxNumArrLoc = MAX_NUM_ARRIVALS;

    // 20200203 AJL - not sure if this is correct, may be OK?  TODO:

    /*if (VpVsRatio > 0.0 && GeometryMode == MODE_GLOBAL) {
                                            nll_puterr("ERROR: cannot use VpVsRatio>0 with TRANSFORM GLOBAL.");

                                            return (EXIT_ERROR_LOCATE);
                                        }*/

    return (0);
}

/** function to read fixed origin time parameters ***/

int GetNLLoc_FixOriginTime(char* line1) {
    int istat;


    istat = sscanf(line1, "%d %d %d %d %d %lf",
            &Hypocenter.year, &Hypocenter.month, &Hypocenter.day,
            &Hypocenter.hour, &Hypocenter.min, &Hypocenter.sec);

    snprintf(MsgStr, sizeof(MsgStr),
            "LOCFIXOTIME:  %4.4d%2.2d%2.2d %2.2d%2.2d %5.2lf",
            Hypocenter.year, Hypocenter.month, Hypocenter.day,
            Hypocenter.hour, Hypocenter.min, Hypocenter.sec);
    nll_putmsg(3, MsgStr);

    if (istat != 6)
        return (-1);

    FixOriginTimeFlag = 1;

    return (0);
}

/** function to read grid params */

int GetNLLoc_Grid(char* input_line) {
    int istat;
    char str_save[20];

    istat = sscanf(input_line, "%d %d %d %lf %lf %lf %lf %lf %lf %s %s",
            &(grid_in.numx), &(grid_in.numy), &(grid_in.numz),
            &(grid_in.origx), &(grid_in.origy), &(grid_in.origz),
            &(grid_in.dx), &(grid_in.dy), &(grid_in.dz), grid_in.chr_type,
            str_save);

    convert_grid_type(&grid_in, 1);
    if (message_flag >= 2)
        display_grid_param(&grid_in);
    snprintf(MsgStr, sizeof(MsgStr), "LOCGRID: Save: %s", str_save);
    nll_putmsg(3, MsgStr);

    if (istat != 11)
        return (-1);

    if (NumLocGrids < MAX_NUM_LOCATION_GRIDS) {
        LocGrid[NumLocGrids] = grid_in;
        LocGrid[NumLocGrids].autox = 0;
        LocGrid[NumLocGrids].autoy = 0;
        LocGrid[NumLocGrids].autoz = 0;
        if (LocGrid[NumLocGrids].origx < -LARGE_DOUBLE)
            LocGrid[NumLocGrids].autox = 1;
        if (LocGrid[NumLocGrids].origy < -LARGE_DOUBLE)
            LocGrid[NumLocGrids].autoy = 1;
        if (LocGrid[NumLocGrids].origz < -LARGE_DOUBLE)
            LocGrid[NumLocGrids].autoz = 1;
        if (strcmp(str_save, "SAVE") == 0)
            LocGridSave[NumLocGrids] = 1;
        else
            LocGridSave[NumLocGrids] = 0;
        NumLocGrids++;
    } else
        nll_puterr("WARNING: maximum number of location grids exceeded.");

    return (0);
}

/** function to read station distance weighting params ***/

int GetStaWeight(char* line1) {
    int istat, ierr;


    istat = sscanf(line1, "%d %lf", &iSetStationDistributionWeights, &stationDistributionWeightCutoff);

    snprintf(MsgStr, sizeof(MsgStr), "LOCSTAWT:  flag: %d  CutoffDist: %f",
            iSetStationDistributionWeights, stationDistributionWeightCutoff);
    nll_putmsg(3, MsgStr);

    ierr = 0;
    //if (checkRangeDouble("LOCSTAWT", "Station distribution weight cutoff distance",
    //		stationDistributionWeightCutoff, 1, 0.0, 0, 0.0) != 0)
    //	ierr = -1;

    if (ierr < 0 || istat != 2)

        return (-1);

    return (0);

}

/** function to read gaussian params ***/

int GetNLLoc_Gaussian2(char* line1) {
    int istat, ierr;


    istat = sscanf(line1, "%lf %lf %lf", &(Gauss2.SigmaTfraction), &(Gauss2.SigmaTmin), &(Gauss2.SigmaTmax));

    snprintf(MsgStr, sizeof(MsgStr), "LOCGAUSS2:  SigmaTfraction: %lf  SigmaTmin: %lf  SigmaTmax: %lf",
            Gauss2.SigmaTfraction, Gauss2.SigmaTmin, Gauss2.SigmaTmax);
    //nll_putmsg(1, MsgStr);
    nll_putmsg(3, MsgStr);

    ierr = 0;
    if (checkRangeDouble("LOCGAU2", "SigmaTfraction",
            Gauss2.SigmaTfraction, 1, 0.0, 1, 1.0) != 0)
        ierr = -1;
    if (checkRangeDouble("LOCGAU2", "SigmaTmin",
            Gauss2.SigmaTmin, 1, 0.0, 0, 0.0) != 0)
        ierr = -1;
    if (checkRangeDouble("LOCGAU2", "SigmaTmax",
            Gauss2.SigmaTmax, 1, 0.0, 0, 0.0) != 0)
        ierr = -1;

    if (ierr < 0 || istat != 3)
        return (-1);

    iUseGauss2 = 1;

    return (0);

}

/** function to read gaussian params ***/

int GetNLLoc_Gaussian(char* line1) {
    int istat, ierr;


    istat = sscanf(line1, "%lf %lf", &(Gauss.SigmaT), &(Gauss.CorrLen));

    snprintf(MsgStr, sizeof(MsgStr), "LOCGAUSS:  SigmaT: %lf  CorrLen: %lf",
            Gauss.SigmaT, Gauss.CorrLen);
    nll_putmsg(3, MsgStr);

    ierr = 0;
    if (checkRangeDouble("LOCGAU", "SigmaT",
            Gauss.SigmaT, 1, 0.0, 0, 0.0) != 0)
        ierr = -1;
    if (checkRangeDouble("LOCGAU", "CorrLen",
            Gauss.CorrLen, 1, 0.0, 0, 0.0) != 0)
        ierr = -1;

    if (ierr < 0 || istat != 2)

        return (-1);

    return (0);
}

/** function to read magnitude calculation type ***/

int GetNLLoc_Magnitude(char* line1) {
    int istat, ierr;

    char mag_type[MAXLINE];

    if (NumMagnitudeMethods >= MAX_NUM_MAG_METHODS) {
        nll_puterr2("ERROR: maximum number of LOCMAG statements read: ignoring: ", line1);
        return (-1);
    }

    istat = sscanf(line1, "%s", mag_type);

    if (istat != 1)
        return (-1);

    if (strcmp(mag_type, "ML_HB") == 0) {

        // default values
        Magnitude[NumMagnitudeMethods].hb_Ro = 100.0;
        Magnitude[NumMagnitudeMethods].hb_Mo = 3.0;

        Magnitude[NumMagnitudeMethods].type = MAG_ML_HB;
        istat = sscanf(line1, "%s %lf %lf %lf %lf %lf",
                mag_type, &(Magnitude[NumMagnitudeMethods].amp_fact_ml_hb),
                &(Magnitude[NumMagnitudeMethods].hb_n), &(Magnitude[NumMagnitudeMethods].hb_K),
                &(Magnitude[NumMagnitudeMethods].hb_Ro), &(Magnitude[NumMagnitudeMethods].hb_Mo));
        snprintf(MsgStr, sizeof(MsgStr), "LOCMAGNITUDE:  Type: %s  f %e  n %f  K %f  Ro %f  Mo %f",
                mag_type, Magnitude[NumMagnitudeMethods].amp_fact_ml_hb, Magnitude[NumMagnitudeMethods].hb_n,
                Magnitude[NumMagnitudeMethods].hb_K,
                Magnitude[NumMagnitudeMethods].hb_Ro, Magnitude[NumMagnitudeMethods].hb_Mo);
        nll_putmsg(3, MsgStr);

        ierr = 0;
        if (checkRangeDouble("LOCMAG", "f", Magnitude[NumMagnitudeMethods].amp_fact_ml_hb, 1, 0.0, 0, 0.0) != 0)
            ierr = -1;

        if (ierr < 0 || istat < 4)
            return (-1);

    } else if (strcmp(mag_type, "MD_FMAG") == 0) {

        Magnitude[NumMagnitudeMethods].type = MAG_MD_FMAG;
        istat = sscanf(line1, "%s %lf %lf %lf %lf %lf",
                mag_type, &(Magnitude[NumMagnitudeMethods].fmag_c1), &(Magnitude[NumMagnitudeMethods].fmag_c2),
                &(Magnitude[NumMagnitudeMethods].fmag_c3), &(Magnitude[NumMagnitudeMethods].fmag_c4), &(Magnitude[NumMagnitudeMethods].fmag_c5));
        snprintf(MsgStr, sizeof(MsgStr), "LOCMAGNITUDE:  Type: %s  C1 %lf  C2 %lf  C3 %lf  C4 %lf  C5 %lf",
                mag_type, Magnitude[NumMagnitudeMethods].fmag_c1, Magnitude[NumMagnitudeMethods].fmag_c2, Magnitude[NumMagnitudeMethods].fmag_c3,
                Magnitude[NumMagnitudeMethods].fmag_c4, Magnitude[NumMagnitudeMethods].fmag_c5);
        nll_putmsg(3, MsgStr);

        if (istat != 6)
            return (-1);

    } else {
        Magnitude[NumMagnitudeMethods].type = MAG_UNDEF;
        nll_puterr2("ERROR: unrecognized magnitude calculation type:", mag_type);
    }

    NumMagnitudeMethods++;

    return (0);
}

/** function to read phase statistics params ***/

int GetNLLoc_PhaseStats(char* line1) {
    int istat;


    istat = sscanf(line1, "%lf %d %lf %lf %lf %lf %lf %lf %lf",
            &RMS_Max, &NRdgs_Min, &Gap_Max, &P_ResidualMax, &S_ResidualMax, &Ell_Len3_Max, &Hypo_Depth_Min, &Hypo_Depth_Max, &Hypo_Dist_Max);

    if (istat < 6)
        Ell_Len3_Max = VERY_LARGE_DOUBLE;
    if (istat < 7)
        Hypo_Depth_Min = -VERY_LARGE_DOUBLE;
    if (istat < 8)
        Hypo_Depth_Max = VERY_LARGE_DOUBLE;
    if (istat < 9)
        Hypo_Dist_Max = VERY_LARGE_DOUBLE;

    snprintf(MsgStr, sizeof(MsgStr),
            "LOCPHSTAT:  RMS_Max: %f  NRdgs_Min: %d  Gap_Max: %.3g  P_ResidualMax: %.3g S_ResidualMax: %.3g Ell_Len3_Max %.3g Hypo_Depth_min %.3g Hypo_Depth_max %.3g Hypo_Dist_Max %.3g",
            RMS_Max, NRdgs_Min, Gap_Max,
            P_ResidualMax, S_ResidualMax, Ell_Len3_Max, Hypo_Depth_Min, Hypo_Depth_Max, Hypo_Dist_Max);
    nll_putmsg(3, MsgStr);

    if (istat < 5)

        return (-1);

    return (0);
}

/** function to read angles mode params ***/

int GetNLLoc_Angles(char* line1) {
    char strAngleMode[MAXLINE];


    sscanf(line1, "%s %d", strAngleMode, &iAngleQualityMin);

    snprintf(MsgStr, sizeof(MsgStr), "LOCANGLES:  %s  %d", strAngleMode, iAngleQualityMin);
    nll_putmsg(4, MsgStr);

    if (strcmp(strAngleMode, "ANGLES_YES") == 0)
        angleMode = ANGLE_MODE_YES;
    else if (strcmp(strAngleMode, "ANGLES_NO") == 0)
        angleMode = ANGLE_MODE_NO;
    else {
        angleMode = ANGLE_MODE_UNDEF;
        nll_puterr("ERROR: unrecognized angle mode");

        return (-1);
    }

    return (0);

}

/** function to read component description ***/

int GetCompDesc(char* line1) {
    int istat, ierr;

    if (NumCompDesc >= MAX_NUM_COMP_DESC) {
        snprintf(MsgStr, sizeof(MsgStr), "%s", line1);
        nll_putmsg(1, MsgStr);
        snprintf(MsgStr, sizeof(MsgStr),
                "WARNING: maximum number of component descriptions reached, ignoring description.");
        nll_putmsg(1, MsgStr);
        return (-1);
    }

    Component[NumCompDesc].sta_corr_md_fmag = 1.0; // fmag default

    istat = sscanf(line1, "%s %s %s %lf %lf %lf",
            Component[NumCompDesc].label, Component[NumCompDesc].inst,
            Component[NumCompDesc].comp, &(Component[NumCompDesc].amp_fact_ml_hb),
            &(Component[NumCompDesc].sta_corr_ml_hb),
            &(Component[NumCompDesc].sta_corr_md_fmag));

    snprintf(MsgStr, sizeof(MsgStr),
            "LOCCMP:  Label: %s  Inst: %s  Comp: %s  Afact: %lf  StaCorr_ML_HB: %lf  StaCorr_MD_FMAG: %lf",
            Component[NumCompDesc].label, Component[NumCompDesc].inst,
            Component[NumCompDesc].comp, Component[NumCompDesc].amp_fact_ml_hb,
            Component[NumCompDesc].sta_corr_ml_hb,
            Component[NumCompDesc].sta_corr_md_fmag);
    nll_putmsg(3, MsgStr);

    ierr = 0;
    if (checkRangeDouble("LOCCMP", "amp_fact_ml_hb",
            Component[NumCompDesc].amp_fact_ml_hb, 1, 0.0, 0, 0.0) != 0)
        ierr = -1;

    if (ierr < 0 || istat < 5)
        return (-1);

    NumCompDesc++;

    return (0);
}

/** function to read arrival label alias ***/

int GetLocAlias(char* line1) {

    if (NumLocAlias >= MAX_NUM_LOC_ALIAS) {
        snprintf(MsgStr, sizeof(MsgStr), "%s", line1);
        nll_putmsg(1, MsgStr);
        snprintf(MsgStr, sizeof(MsgStr),
                "WARNING: maximum number of aliases reached, ignoring alias.");
        nll_putmsg(1, MsgStr);
        return (-1);
    }

    sscanf(line1, "%s %s  %d %d %d  %d %d %d",
            LocAlias[NumLocAlias].name, LocAlias[NumLocAlias].alias,
            &(LocAlias[NumLocAlias].byr), &(LocAlias[NumLocAlias].bmo),
            &(LocAlias[NumLocAlias].bday),
            &(LocAlias[NumLocAlias].eyr), &(LocAlias[NumLocAlias].emo),
            &(LocAlias[NumLocAlias].eday));

    snprintf(MsgStr, sizeof(MsgStr),
            "LOCALIAS:  Name: %s  Alias: %s  Valid: %4.4d %2.2d %2.2d -> %4.4d %2.2d %2.2d",
            LocAlias[NumLocAlias].name, LocAlias[NumLocAlias].alias,
            LocAlias[NumLocAlias].byr, LocAlias[NumLocAlias].bmo,
            LocAlias[NumLocAlias].bday,
            LocAlias[NumLocAlias].eyr, LocAlias[NumLocAlias].emo,
            LocAlias[NumLocAlias].eday);
    nll_putmsg(3, MsgStr);

    NumLocAlias++;

    return (0);
}

/** function to read exclude arrival label and phase ***/

int GetLocExclude(char* line1) {

    if (NumLocExclude >= MAX_NUM_LOC_EXCLUDE) {
        snprintf(MsgStr, sizeof(MsgStr), "%s", line1);
        nll_putmsg(1, MsgStr);
        snprintf(MsgStr, sizeof(MsgStr),
                "WARNING: maximum number of LOCEXCLUDE phases reached, ignoring exclude.");
        nll_putmsg(1, MsgStr);
        return (-1);
    }

    sscanf(line1, "%s %s",
            LocExclude[NumLocExclude].label, LocExclude[NumLocExclude].phase);

    if (message_flag >= 3) {
        snprintf(MsgStr, sizeof(MsgStr), "LOCEXCLUDE:  Name: %s  Phase: %s",
                LocExclude[NumLocExclude].label, LocExclude[NumLocExclude].phase);
        nll_putmsg(3, MsgStr);
    }

    NumLocExclude++;

    return (0);
}

/** function to read exclude arrival label and phase ***/

int GetLocInclude(char* line1) {

    if (NumLocInclude >= MAX_NUM_LOC_INCLUDE) {
        snprintf(MsgStr, sizeof(MsgStr), "%s", line1);
        nll_putmsg(1, MsgStr);
        snprintf(MsgStr, sizeof(MsgStr),
                "WARNING: maximum number of LOCINCLUDE phases reached, ignoring include.");
        nll_putmsg(1, MsgStr);
        return (-1);
    }

    sscanf(line1, "%s %s",
            LocInclude[NumLocInclude].label, LocInclude[NumLocInclude].phase);

    if (message_flag >= 3) {
        snprintf(MsgStr, sizeof(MsgStr), "LOCINCLUDE:  Name: %s  Phase: %s",
                LocInclude[NumLocInclude].label, LocInclude[NumLocInclude].phase);
        nll_putmsg(3, MsgStr);
    }

    NumLocInclude++;

    return (0);
}

/** function to read station phase time delays ***/

int GetTimeDelays(char* line1) {

    if (NumTimeDelays >= MAX_NUM_STA_DELAYS) {
        snprintf(MsgStr, sizeof(MsgStr), "%s", line1);
        nll_putmsg(3, MsgStr);
        snprintf(MsgStr, sizeof(MsgStr),
                "WARNING: maximum number of station delays reached, ignoring alias.");
        nll_putmsg(2, MsgStr);
        return (-1);
    }

    sscanf(line1, "%s %s %d %lf %lf",
            TimeDelay[NumTimeDelays].label, TimeDelay[NumTimeDelays].phase,
            &(TimeDelay[NumTimeDelays].n_residuals),
            &(TimeDelay[NumTimeDelays].delay),
            &(TimeDelay[NumTimeDelays].std_dev));

    if (message_flag >= 3) {
        snprintf(MsgStr, sizeof(MsgStr),
                "LOCDELAY:  Label: %s  Phase: %s  NumResiduals: %d  TimeDelay: %lf  StdDev: %lf",
                TimeDelay[NumTimeDelays].label, TimeDelay[NumTimeDelays].phase,
                TimeDelay[NumTimeDelays].n_residuals,
                TimeDelay[NumTimeDelays].delay,
                TimeDelay[NumTimeDelays].std_dev);
        nll_putmsg(3, MsgStr);
    }

    NumTimeDelays++;

    return (0);
}

/** function to read topo surface (GMT GRD file with x=long, y=lat, z=delay in sec ***/

int GetTopoSurface(char* line1) {

    int idump_decimation = 0;
    char dump_file[FILENAME_MAX];


    // initialize topo surface fields
    topo_surface = model_surface + (MAX_SURFACES - 1);
    topo_surface_index = MAX_SURFACES - 1;

    sscanf(line1, "%s %d", topo_surface->grd_file, &idump_decimation);

    snprintf(MsgStr, sizeof(MsgStr), "LOCTOPO_SURFACE:  GMT GRD File: %s  Dump to file decimation: %d", topo_surface->grd_file, idump_decimation);
    nll_putmsg(3, MsgStr);
    //nll_putmsg(0, MsgStr);

    if (read_grd(topo_surface, message_flag >= 2) < 0) {
        nll_puterr2("ERROR: reading Topo Surface GMT GRD File: ",
                topo_surface->grd_file);
        return (-1);
    }

    // print grid limits info (for seismicitydefaults)
    double lat_ul, lon_ul, lat_ur, lon_ur, lat_lr, lon_lr, lat_ll, lon_ll;
    if (!topo_surface->is_latlon) {
        rect2latlon(0, topo_surface->hdr->x_min, topo_surface->hdr->y_max, &lat_ul, &lon_ul);
        rect2latlon(0, topo_surface->hdr->x_max, topo_surface->hdr->y_max, &lat_ur, &lon_ur);
        rect2latlon(0, topo_surface->hdr->x_max, topo_surface->hdr->y_min, &lat_lr, &lon_lr);
        rect2latlon(0, topo_surface->hdr->x_min, topo_surface->hdr->y_min, &lat_ll, &lon_ll);
        snprintf(MsgStr, sizeof(MsgStr), "LOCTOPO_SURFACE:  FileURL; lat, long upper left; lat, long upper right; lat, long lower right; lat, long lower left;");
        nll_putmsg(1, MsgStr);
        snprintf(MsgStr, sizeof(MsgStr), "LOCTOPO_SURFACE:  %s; %f,%f; %f,%f; %f,%f; %f,%f;",
                topo_surface->grd_file, lat_ul, lon_ul, lat_ur, lon_ur, lat_lr, lon_lr, lat_ll, lon_ll);
        nll_putmsg(1, MsgStr);
    }


    if (idump_decimation) {

        snprintf(dump_file, sizeof(dump_file), "%s", topo_surface->grd_file);
        strncat(dump_file, ".bin", sizeof(dump_file) - strlen(dump_file) - 1);
        dump_grd(topo_surface_index, idump_decimation, 1.0, 1.0, -0.001, dump_file);
        snprintf(MsgStr, sizeof(MsgStr), "LOCTOPO_SURFACE:  Grid dumped to: %s", dump_file);
        nll_putmsg(1, MsgStr);
    }

    return (0);
}

/** function to read time delay surface (GMT GRD file with x=long, y=lat, z=delay in sec ***/

int GetTimeDelaySurface(char* line1) {

    sscanf(line1, "%s %lf %s",
            TimeDelaySurfacePhase[NumTimeDelaySurface],
            &TimeDelaySurfaceMultiplier[NumTimeDelaySurface],
            model_surface[NumTimeDelaySurface].grd_file);

    if (message_flag >= 1) {
        snprintf(MsgStr, sizeof(MsgStr), "LOCDELAY_SURFACE:  Phase: %s  Mult: %f  GMT GRD File: %s",
                TimeDelaySurfacePhase[NumTimeDelaySurface],
                TimeDelaySurfaceMultiplier[NumTimeDelaySurface],
                model_surface[NumTimeDelaySurface].grd_file);
        //nll_putmsg(3, MsgStr);
        nll_putmsg(1, MsgStr);
    }

    if (read_grd(&model_surface[NumTimeDelaySurface], message_flag > 2) < 0) {
        nll_puterr2("ERROR: reading Surface Delay GMT GRD File: ",
                model_surface[NumTimeDelaySurface].grd_file);
        return (-1);
    }

    NumTimeDelaySurface++;

    return (0);
}

/** function to read elevation correction params ***/

int GetElevCorr(char* line1) {

    int istat;

    istat = sscanf(line1, "%d %lf %lf",
            &ApplyElevCorrFlag, &ElevCorrVelP, &ElevCorrVelS);

    snprintf(MsgStr, sizeof(MsgStr), "LOCELEVCORR:  Flag: %d  VelP: %lf  VelS: %lf",
            ApplyElevCorrFlag, ElevCorrVelP, ElevCorrVelS);
    //nll_putmsg(3, MsgStr);
    nll_putmsg(1, MsgStr);

    if (istat != 3)

        return (-1);

    return (0);
}
