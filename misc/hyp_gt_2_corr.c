/*
 * Copyright (C) 1999-2018 Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
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


/*   hyp_gt_2_corr.c

        Program to combine hypocenter and probabilistic ground truth information into station corrections

 */



/*
        history:

        ver 01    04Oct2018  AJL  Original version


 */



#include "../src/GridLib.h"
#include "../src/ran1/ran1.h"
#include "../src/velmod.h"
#include "../src/GridMemLib.h"
#include "../src/calc_crust_corr.h"
#include "../src/phaseloclist.h"
#include "../src/otime_limit.h"
#include "../src/NLLocLib.h"


/* defines */


/* globals */


/* functions */

int HypGt2Corr(int argc, char** argv);



/*** program to sum event scatter files */

#define PNAME  "hyp_gt_2_corr"

int main(int argc, char** argv) {

    int istat, narg;


    // set program name
    strcpy(prog_name, PNAME);


    // check command line for correct usage
    fprintf(stdout, "%s Arguments: ", prog_name);
    for (narg = 0; narg < argc; narg++)
        fprintf(stdout, "<%s> ", argv[narg]);
    fprintf(stdout, "\n");

    if (argc < 2) {
        nll_puterr("ERROR wrong number of command line arguments.");
        disp_usage(PNAME,
                "<output_file_root> <hyp_gt_file> [<locphstat_statement> [latCent,longCent,radius [fixDepth]]]"
                );
        exit(-1);
    }

    //SetConstants();
    message_flag = 1;
    DispProgInfo();
    message_flag = 0;

    if ((istat = HypGt2Corr(argc, argv)) < 0) {
        nll_puterr("ERROR doing HypGt2Corr process.");
        exit(-1);
    }



    exit(0);

}

/** function to calculate distance between 2 XYZ points */

double Dist2D_LocSum(double x1, double x2, double y1, double y2, int latlon) {

    if (latlon) {
        return (GCDistance(y1, x1, y2, x2));
        //return(EllipsoidDistance(yval, xval, psta->y, psta->x));
    } else {
        double dx, dy;
        dx = x1 - x2;
        dy = y1 - y2;
        return (sqrt(dx * dx + dy * dy));
    }
}

int HypGt2Corr(int argc, char** argv) {

    int istat;
    char *cstat, *cstat2;
    char line[MAXLINE_LONG];
    char line2[MAXLINE_LONG];

    char locphstat_statement[MAXLINE_LONG];

    char fn_control[FILENAME_MAX];
    FILE *fp_control;
    char fn_control_temp[FILENAME_MAX];
    FILE *fp_control_temp;
    char fn_sta_corr[FILENAME_MAX];
    FILE *fp_sta_corr;
    char fn_root_out[FILENAME_MAX];
    char fn_hyp_gt[FILENAME_MAX];
    FILE *fp_hyp_gt;

    char sys_command[2 * FILENAME_MAX];

    char event_id[FILENAME_MAX];
    char fn_hyp[FILENAME_MAX];
    FILE *fp_hypo;
    char fn_ttime_grid_root[FILENAME_MAX];

    GridDesc locgrid;
    HypoDesc Hypo;

    LocNode *loc_list_head = NULL; // root node of location list

    int index;
    char uwi_api[32];
    double mean_lat, mean_lon, mean_depth, prob;


    int NRdgs_Min_local;
    double RMS_Max_local, Gap_Max_local;
    double P_ResidualMax_local;
    double S_ResidualMax_local;
    double Ell_Len3_Max_local;
    double Hypo_Depth_Min_local;
    double Hypo_Depth_Max_local;
    double Hypo_Dist_Max_local;

    // get command line parameters
    sprintf(fn_root_out, "%s", argv[1]);
    strcpy(fn_hyp_gt, argv[2]);
    int ihave_locphstat_statement = 0;
    if (argc > 3) {
        strcpy(locphstat_statement, argv[3]);
        ihave_locphstat_statement = 1;
        istat = sscanf(locphstat_statement, "LOCPHSTAT %lf %d %lf %lf %lf %lf %lf %lf %lf",
                &RMS_Max_local, &NRdgs_Min_local, &Gap_Max_local, &P_ResidualMax_local, &S_ResidualMax_local,
                &Ell_Len3_Max_local, &Hypo_Depth_Min_local, &Hypo_Depth_Max_local, &Hypo_Dist_Max_local);
        if (istat < 9) {
            printf("ERROR: reading LOCPHSTAT parameters: %s\n", locphstat_statement);
        }
    }

    double latRad = 0.0, longRad = 0.0, radius = -1.0;
    if (argc > 4) {
        istat = sscanf(argv[4], "%lf,%lf,%lf", &latRad, &longRad, &radius);
        if (istat < 3) {
            printf("ERROR: reading latCent,longCent,radius parameters: %s\n", argv[4]);
        }
        printf("Geog Radius Cut: Lat: %f, Long: %f, Rad: %f\n", latRad, longRad, radius);
    }

    double fixDepth = FLT_MAX;
    if (argc > 5) {
        istat = sscanf(argv[5], "%lf", &fixDepth);
        if (istat < 1) {
            printf("ERROR: reading fixDepth parameter: %s\n", argv[4]);
        }
        printf("Fix Depth: fixDepth: %f\n", fixDepth);
    }

    SetConstants();
    int nHypRead = 0;
    int nHypAccepted = 0;
    int nHypProcessed = 0;
    int nHypGtProcessed = 0;
    int nWriteStaStatTable = 0;

    int nSta_corr_grid = MAX_NUM_LOCATION_GRIDS - 1; // use last grid so that will not be freed in FreeStaStatTable in NLLoc()

    // open hypocenter - ground truth file
    //sprintf(strstr(fn_hyp_in_list[nFile], test_str), "\0");
    if ((fp_hyp_gt = fopen(fn_hyp_gt, "r")) == NULL) {
        printf("ERROR: opening hypocenter - ground truth file: %s\n", fn_hyp_gt);
        return (-1);
    }

    // loop over events in file

    // read first line
    cstat = fgets(line, MAXLINE_LONG, fp_hyp_gt);
    //printf("0 %s", line);
    while (1) {
        // TEST while (nHypProcessed < 10) {

        if (cstat == NULL)
            break;

        while (strncmp(line, "Event", 5) != 0) {
            // read next line
            cstat = fgets(line, MAXLINE_LONG, fp_hyp_gt);
            //printf("1 %s", line);
            if (cstat == NULL)
                break;
        }
        if (cstat == NULL)
            break;

        // parse event line
        // Event: 2017cnti /temp/TexasNetwork/WestTexas/westTX_3Dmodel_global/loc_20180917/loc/westTX_3Dmodel_global.20170205.233532.grid0.loc.hyp
        istat = sscanf(line, "Event: %s %s", event_id, fn_hyp);
        //printf("fn_hyp %s\n", fn_hyp);
        // read hypocenter
        if ((fp_hypo = fopen(fn_hyp, "r")) == NULL) {
            printf("ERROR: opening hypocenter file, ignoring event, file: %s\n", fn_hyp);
            strcpy(line, "XXX"); // force read to next event
            continue;
        }
        istat = GetHypLoc(fp_hypo, NULL, &Hypo, Arrival, &NumArrivals, 0, &locgrid, 0);
        fclose(fp_hypo);
        if (istat == EOF) {
            strcpy(line, "XXX"); // force read to next event
            continue;
        }
        if (istat < 0) {
            printf("ERROR: reading hypocenter file, ignoring event, file: %s\n", fn_hyp);
            strcpy(line, "XXX"); // force read to next event
            continue;
        }
        nHypRead++;
        if (strncmp(Hypo.locStat, "ABORTED", 7) == 0) {
            printf("WARNING: location ABORTED, ignoring event: %s\n", fn_hyp);
            strcpy(line, "XXX"); // force read to next event
            continue;
        } else if (strncmp(Hypo.locStat, "REJECTED", 8) == 0) {
            printf("WARNING: location REJECTED, ignoring event: %s\n", fn_hyp);
            strcpy(line, "XXX"); // force read to next event
            continue;
        }

        HypoDesc* phypo_test = &Hypo;
        int ireject = 0;

        double distance;
        if (radius > 0.0 && (distance = Dist2D_LocSum(phypo_test->dlong, longRad, phypo_test->dlat, latRad, 1)) > radius) {
            printf(">>> Reject stats: phypo distance %f > radius %f\n", distance, radius);
            ireject = 1;
        }
        if (phypo_test->rms > RMS_Max_local) {
            printf(">>> Reject stats: phypo_test->rms %f > RMS_Max_local %f\n", phypo_test->rms, RMS_Max_local);
            ireject = 1;
        }
        if (phypo_test->nreadings < NRdgs_Min_local) {
            printf(">>> Reject stats: hypo->nreadings %d < NRdgs_Min_local %d\n", phypo_test->nreadings, NRdgs_Min_local);
            ireject = 1;
        }
        if (phypo_test->gap > Gap_Max_local) {
            printf(">>> Reject stats: phypo_test->gap %f > Gap_Max_local %f\n", phypo_test->gap, Gap_Max_local);
            ireject = 1;
        }
        if (phypo_test->ellipsoid.len3 > Ell_Len3_Max_local) {
            printf(">>> Reject stats: phypo_test->ellipsoid.len3 %f > Ell_Len3_Max_local %f\n", phypo_test->ellipsoid.len3, Ell_Len3_Max_local);
            ireject = 1;
        }
        if (phypo_test->z < Hypo_Depth_Min_local) {
            printf(">>> Reject stats: phypo_test->z %f < Hypo_Depth_Min_local %f\n", phypo_test->z, Hypo_Depth_Min_local);
            ireject = 1;
        }
        if (phypo_test->z > Hypo_Depth_Max_local) {
            printf(">>> Reject stats: phypo_test->z %f > Hypo_Depth_Max_local %f\n", phypo_test->z, Hypo_Depth_Max_local);
            ireject = 1;
        }
        if (ireject) {
            printf("WARNING: location REJECTED STATS, ignoring event: %s\n", fn_hyp);
            strcpy(line, "XXX"); // force read to next event
            continue;
        }

        nHypAccepted++;

        // construct name and read last.in control file, and get travel-time root
        // construct name
        strcpy(fn_control, fn_hyp);
        char *chr = strrchr(fn_control, '/');
        if (chr != NULL) {
            *(chr + 1) = '\0';
            strcat(fn_control, "last.in");
        } else {
            strcpy(fn_control, "last.in");
        }
        //printf("fn_control %s\n", fn_control);
        // open and read last.in control file
        strcpy(fn_ttime_grid_root, "");
        if ((fp_control = fopen(fn_control, "r")) == NULL) {
            printf("ERROR: opening control file %s, ignoring event, file: %s\n", fn_control, fn_hyp);
            strcpy(line, "XXX"); // force read to next event
            continue;
        }
        int iSwapBytesOnInput = 0;
        int flag_locfiles = 0;
        int flag_trans = 0;
        int flag_locphstat = ihave_locphstat_statement;
        cstat2 = fgets(line2, MAXLINE_LONG, fp_control);
        //printf("0 %s", line2);
        while (!(flag_locfiles * flag_trans * flag_locphstat)) {
            if (strncmp(line2, "LOCFILES", 8) == 0) {
                // read LOCFILES
                istat = sscanf(line2, "LOCFILES %*s %*s %s %*s %d", fn_ttime_grid_root, &iSwapBytesOnInput);
                if (istat < 1) {
                    printf("ERROR: reading LOCFILES parameters: %s\n", line2);
                } else {
                    flag_locfiles = 1;
                    //printf("fn_ttime_grid_root %s\n", fn_ttime_grid_root);
                }
                if (istat < 2)
                    iSwapBytesOnInput = 0;
            } else if (strncmp(line2, "TRANS", 5) == 0) {
                if ((istat = get_transform(0, strchr(line2, ' '))) < 0) {
                    printf("ERROR: reading TRANS parameters: %s\n", line2);
                } else {
                    flag_trans = 1;
                    //printf("%s", line2);
                }
            } else if (!ihave_locphstat_statement && strncmp(line2, "LOCPHSTAT", 9) == 0) {
                istat = sscanf(line2, "LOCPHSTAT %lf %d %lf %lf %lf %lf %lf %lf %lf",
                        &RMS_Max_local, &NRdgs_Min_local, &Gap_Max_local, &P_ResidualMax_local, &S_ResidualMax_local, &Ell_Len3_Max_local, &Hypo_Depth_Min_local, &Hypo_Depth_Max_local, &Hypo_Dist_Max_local);
                if (istat < 9) {
                    printf("ERROR: reading LOCPHSTAT parameters: %s\n", line2);
                } else {
                    flag_locphstat = 1;
                    //printf("%s", line2);
                }
            }
            // read next line
            cstat2 = fgets(line2, MAXLINE_LONG, fp_control);
            //printf("1 %s", line2);
            if (cstat2 == NULL)
                break;
        }
        fclose(fp_control);
        if (!flag_locfiles) {
            printf("ERROR: finding LOCFILES line in %s, ignoring event: %s\n", fn_control, fn_hyp);
            strcpy(line, "XXX"); // force read to next event
            continue;
        }
        if (!flag_trans) {
            printf("ERROR: finding TRANS in %s, ignoring event: %s\n", fn_control, fn_hyp);
            strcpy(line, "XXX"); // force read to next event
            continue;
        }
        if (!flag_locphstat) {
            printf("ERROR: finding LOCPHSTAT in %s, ignoring event: %s\n", fn_control, fn_hyp);
            strcpy(line, "XXX"); // force read to next event
            continue;
        }


        int first = 1;
        while (1) {

            // read next line
            cstat = fgets(line, MAXLINE_LONG, fp_hyp_gt);
            //printf("2 %s", line);
            if (cstat == NULL)
                break;
            if (strncmp(line, "Event", 5) == 0) {
                // no more gt for this event
                break;
            }

            if (first) {
                printf("fn_hyp %s\n", fn_hyp);
                nHypProcessed++;
                first = 0;
            }
            // parse hypo gt line
            istat = sscanf(line, "%d %s %lf %lf %lf %lf", &index, uwi_api, &mean_lat, &mean_lon, &mean_depth, &prob);
            nHypGtProcessed++;

            // copy control file to temp control file
            strcpy(fn_control_temp, "nll_temp.in");
            sprintf(sys_command, "cp %s %s", fn_control, fn_control_temp);
            system(sys_command);

            // comment out existing lines in temp control file
            sprintf(sys_command, "sed -i '' 's/LOCFILES/#LOCFILES/g' %s", fn_control_temp);
            system(sys_command);
            sprintf(sys_command, "sed -i '' 's/LOCGRID/#LOCGRID/g' %s", fn_control_temp);
            system(sys_command);
            sprintf(sys_command, "sed -i '' 's/LOCSEARCH/#LOCSEARCH/g' %s", fn_control_temp);
            system(sys_command);
            // remove any previous station corrections applied
            sprintf(sys_command, "sed -i '' 's/INCLUDE.*stat_totcorr/#INCLUDE*stat_totcorr/g' %s", fn_control_temp);
            system(sys_command);

            // open temp control file for append
            if ((fp_control_temp = fopen(fn_control_temp, "a")) == NULL) {
                printf("ERROR: opening temp control file %s, ignoring event, file: %s\n", fn_control_temp, fn_hyp);
                strcpy(line, "XXX"); // force read to next event
                continue;
            }

            // append modified control file lines
            //
            fprintf(fp_control_temp, "LOCFILES %s %s %s %s %d\n", fn_hyp, "NLLOC_OBS", fn_ttime_grid_root, fn_root_out, iSwapBytesOnInput);
            //
            double mean_x = mean_lon;
            double mean_y = mean_lat;
            double ds = 1.0 / 100.0; // 10 meters
            double ds_xy = ds;
            if (GeometryMode == MODE_GLOBAL) {
                ds_xy *= KM2DEG;
            } else {
                // int latlon2rect(int n_proj, double dlat, double dlong, double* pxrect, double* pyrect);
                //printf("latlon2rect: mean_x %f, mean_y %f ->", mean_x, mean_y);
                latlon2rect(0, mean_lat, mean_lon, &mean_x, &mean_y);
                //printf("latlon2rect: mean_x %f, mean_y %f\n", mean_x, mean_y);
            }
            double xmin = mean_x - ds_xy / 2.0;
            double ymin = mean_y - ds_xy / 2.0;
            //
            // Edit here for: /Users/anthony/work/TexasNetwork/WestTexas/20180903-FRAC/python_TexNetFrac RUN -> HORIZ_ONLY
            double zmin = mean_depth - ds / 2.0; // depth is frac well depth
            if (fixDepth != FLT_MAX) {
                zmin = fixDepth - ds / 2.0; // depth is fixed
            }
            // TEST  ONLY! double zmin = phypo_test->z - ds / 2.0;     // depth is hypo depth TEST ONLY!
            //
            //LOCGRID 401 361 53  -108.0 28.0249 -4.0  0.021774249999999995 0.017986388888888892 2.0 PROB_DENSITY SAVE
            fprintf(fp_control_temp, "LOCGRID  2 2 2  %f %f %f  %f %f %f PROB_DENSITY SAVE\n", xmin, ymin, zmin, ds_xy, ds_xy, ds);
            //
            //LOCSEARCH OCT  20 18 2 0.01 25000 5000 0 1
            fprintf(fp_control_temp, "LOCSEARCH OCT  1 1 1 0.01 110 10 0 1\n");

            fclose(fp_control_temp);

            // run NLLoc
            int nll_stat;
            char pid_main[255]; // string process id
            int return_locations = 1;
            int return_oct_tree_grid = 0;
            int return_scatter_sample = 0;
            nll_stat = NLLoc(pid_main, fn_control_temp, NULL, -1, NULL, -1,
                    return_locations, return_oct_tree_grid, return_scatter_sample, &loc_list_head);

            printf("Finished NLLoc: fn_root_out %s\n", fn_root_out);

            // loop over returned location results
            int id = 0;
            LocNode* locNode = NULL;
            while ((locNode = getLocationFromLocList(loc_list_head, id)) != NULL) {

                //HypoDesc* phypo = locNode->plocation->phypo;
                int narrivals = locNode->plocation->narrivals;
                ArrivalDesc* parrivals = locNode->plocation->parrivals;
                //ArrivalDesc* parr;

                printf("Processing location %d, %d arrivals\n", id, narrivals);

                //sprintf(frootname, "out/%3.3d", id_filename);
                //sprintf(fname, "%s.loc.hyp", frootname);

                // NOTE: the angles in locNode->plocation->ellipsoid, ellipse and
                //    phase/arrival angles in locNode->plocation->parrivals are in NLL internal coordinates
                //    and must be transformed to geographic directions if there is a non-zero rotAngle in TRANS.
                //    For example:
                //          ellipsoid_az1 = rect2latlonAngle(0, phypo->ellipsoid.az1);
                //          ellipsoid_az2 = rect2latlonAngle(0, phypo->ellipsoid.az2);
                //    These two transformations ARE NOT applied in WriteLocation() used below!
                //          sta_azim = rect2latlonAngle(0, parr->azim);
                //          ray_azim = rect2latlonAngle(0, parr->ray_azim);
                //    These two transformations ARE applied in WriteLocation() used below.

                // write NLLoc Hypocenter-Phase file to disk
                //if ((istat = WriteLocation(NULL, locNode->plocation->phypo, locNode->plocation->parrivals,
                //        locNode->plocation->narrivals, fname, 1, 1, 0, locNode->plocation->pgrid, 0)) < 0) {
                //    nll_puterr2("ERROR: writing location to event file: %s", fname);
                //}

                // loop over arrivals
                /*for (int narr = 0; narr < narrivals; narr++) {
                    parr = parrivals + narr;
                    printf("   >> %s %s: hyp_tt %f, sta_res %f, sta_delay %f\n", parr->phase, parr->label, parr->pred_travel_time,
                            parr->residual, parr->delay);
                }*/
                // update station statistics table
                double weight = prob;
                HypoDesc* phypo_test = &Hypo;
                // following no longer needed, tests done above
                /*if (phypo_test->rms > RMS_Max_local)
                    printf(">>> Reject stats: phypo_test->rms %f > RMS_Max_local %f\n", phypo_test->rms, RMS_Max_local);
                if (phypo_test->nreadings < NRdgs_Min_local)
                    printf(">>> Reject stats: hypo->nreadings %d < NRdgs_Min_local %d\n", phypo_test->nreadings, NRdgs_Min_local);
                if (phypo_test->gap > Gap_Max_local)
                    printf(">>> Reject stats: phypo_test->gap %f > Gap_Max_local %f\n", phypo_test->gap, Gap_Max_local);
                if (phypo_test->ellipsoid.len3 > Ell_Len3_Max_local)
                    printf(">>> Reject stats: phypo_test->ellipsoid.len3 %f > Ell_Len3_Max_local %f\n", phypo_test->ellipsoid.len3, Ell_Len3_Max_local);
                if (phypo_test->z < Hypo_Depth_Min_local)
                    printf(">>> Reject stats: phypo_test->z %f < Hypo_Depth_Min_local %f\n", phypo_test->z, Hypo_Depth_Min_local);
                if (phypo_test->z > Hypo_Depth_Max_local)
                    printf(">>> Reject stats: phypo_test->z %f > Hypo_Depth_Max_local %f\n", phypo_test->z, Hypo_Depth_Max_local);
                 */
                int iAcceptRadius = 1;
                if (radius > 0.0 && Dist2D_LocSum(Hypo.dlong, longRad, Hypo.dlat, latRad, 1) > radius) {
                    iAcceptRadius = 0;
                }
                if (
                        iAcceptRadius
                        // following no longer needed, tests done above
                        //strncmp(phypo_test->locStat, "LOCATED", 7) == 0 &&  // highly constrained LOCGRID always rejected!
                        && phypo_test->rms <= RMS_Max_local
                        && phypo_test->nreadings >= NRdgs_Min_local
                        && phypo_test->gap <= Gap_Max_local
                        && phypo_test->ellipsoid.len3 <= Ell_Len3_Max_local
                        && phypo_test->z >= Hypo_Depth_Min_local
                        && phypo_test->z <= Hypo_Depth_Max_local) {
                    UpdateStaStat(nSta_corr_grid, parrivals, narrivals, P_ResidualMax_local, S_ResidualMax_local, Hypo_Dist_Max_local, weight);
                    nWriteStaStatTable++;
                }

                // clean up
                freeLocList(loc_list_head, 1);
                loc_list_head = NULL;

                id++;

                printf("\n");

            }

        }

    }
    fclose(fp_hyp_gt);

    // write delays only
    sprintf(fn_sta_corr, "%s.sum.grid0.loc.hyp_gt_2_corr.stat_totcorr", fn_root_out);
    if ((fp_sta_corr = fopen(fn_sta_corr, "w")) == NULL) {
        printf("ERROR: opening total phase corrections output file: %s\n", fn_sta_corr);
    } else {
        WriteStaStatTable(nSta_corr_grid, fp_sta_corr,
                RMS_Max_local, NRdgs_Min_local, Gap_Max_local,
                P_ResidualMax_local, S_ResidualMax_local, Ell_Len3_Max_local,
                Hypo_Depth_Min_local, Hypo_Depth_Max_local, Hypo_Dist_Max_local, WRITE_RES_DELAYS);
    }
    fclose(fp_sta_corr);


    // write messages
    printf("%d location files read, %d accepted, %d processed, %d hyp-gt processed, %d WriteStaStatTable.\n",
            nHypRead, nHypAccepted, nHypProcessed, nHypGtProcessed, nWriteStaStatTable);
    printf("Total phase corrections output file: %s\n", fn_sta_corr);
    //fprintf(stdout,
    //        "Station corrections written to <%s>\n", num_points_written, nHypRead, fn_hyp_scat_out);


    return (0);

}

