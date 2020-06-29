/*
 * Copyright (C) 2004 Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */


/*   Loc2ddct.c

        Program to sum location files and convert to hypoDD catalog *.ct format.

 */



/*
        history:

        ver 01    05Aug2004  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

 */

#define PNAME  "Loc2ddct"

#include "GridLib.h"


// defines

#define MAX_NUM_INPUT_FILES 4096
//#define SMALL_FILENAME_MAX 256

#define VERBOSE 0


// globals

//char fn_hyp_in_list[MAX_NUM_INPUT_FILES][FILENAME_MAX];
char fn_hyp_in_list[MAX_NUM_INPUT_FILES][FILENAME_MAX];


// functions

int ConvertLocToCT(int argc, char *argv[]);

/*** program to sum event scatter files */

int main(int argc, char *argv[]) {

    int istat, narg;

    // set program name

    strcpy(prog_name, PNAME);


    // check command line for correct usage

    fprintf(stdout, "%s Arguments: ", prog_name);
    for (narg = 0; narg < argc; narg++)
        fprintf(stdout, "<%s> ", argv[narg]);
    fprintf(stdout, "\n");

    if (argc < 5) {
        nll_puterr("ERROR wrong number of command line arguments.");
        disp_usage(PNAME, "<add_file_list> <output_file_root> max_event_dist arrival_weight_min [max_station_dist]");
        printf("Example: Loc2ddct \"/temp/nlloc_tmp/loc/layer.*.*.grid0.loc\" Loc2ddct_layer 10 0.1\n");
        exit(-1);
    }

    SetConstants();

    // re-allocate arrivals arrays
    // allocate arrivals array
    MAX_NUM_STATIONS = X_MAX_NUM_STATIONS_DIFF;
    MAX_NUM_ARRIVALS = MAX_NUM_ARRIVALS_STA * MAX_NUM_STATIONS;
    if ((Arrival = (ArrivalDesc *) malloc(MAX_NUM_ARRIVALS * sizeof (ArrivalDesc))) == NULL) {
        nll_puterr("ERROR: re-allocating Arrival array.");
        return (EXIT_ERROR_MEMORY);
    }

    message_flag = 1;
    DispProgInfo();
    message_flag = 0;

    if ((istat = ConvertLocToCT(argc, argv)) < 0) {
        nll_puterr("ERROR converting to hypoDD catalog-differences format.");
        exit(-1);
    }



    exit(0);

}


#define  MAX_NUM_DIFF_HYPOCENTERS 100000
HypoDesc diffHypos[MAX_NUM_DIFF_HYPOCENTERS];

int ConvertLocToCT(int argc, char *argv[]) {

    int istat, i;

    int nArrivals_read;

    double event_dist_max, epi_dist;

    char fn_hyp_in[FILENAME_MAX], fn_root_out[FILENAME_MAX];
    char fn_hyp_out[FILENAME_MAX], fn_diff_out[FILENAME_MAX], fn_xyz_out[FILENAME_MAX];
    FILE *fp_hyp_in, *fp_hyp_out, *fp_diff_out, *fp_xyz_out;

    int nHypo, numFiles, nLocWritten, nLocAccepted;
    char test_str[10];

    GridDesc locgrid;
    HypoDesc *phyp1, *phyp2;

    double weight, weight_min;

    int nwritten, narr1, narr2;
    ArrivalDesc *parr1, *parr2;
    long dd_event_id_1, dd_event_id_2;

    strcpy(test_str, ".hyp");

    // get command line parameters
    strcpy(fn_hyp_in, argv[1]);
    strcpy(fn_root_out, argv[2]);
    sscanf(argv[3], "%lf", &event_dist_max);
    sscanf(argv[4], "%lf", &weight_min);
    double station_dist_max = -1.0;
    if (argc > 5) {
        sscanf(argv[5], "%lf", &station_dist_max);
    }


    // open output files
    sprintf(fn_hyp_out, "%s.hyp", fn_root_out);
    if ((fp_hyp_out = fopen(fn_hyp_out, "w")) == NULL) {
        nll_puterr("ERROR: opening hyp output file.");
        return (-1);
    }
    sprintf(fn_diff_out, "%s.ct", fn_root_out);
    if ((fp_diff_out = fopen(fn_diff_out, "w")) == NULL) {
        nll_puterr("ERROR: opening ct output file.");
        return (-1);
    }
    sprintf(fn_xyz_out, "%s.xyz", fn_root_out);
    if ((fp_xyz_out = fopen(fn_xyz_out, "w")) == NULL) {
        nll_puterr("ERROR: opening xyz output file.");
        return (-1);
    }

    // sum requested loc files into output grid

    // check for wildcards in input file name
    strcat(fn_hyp_in, test_str);
    if ((numFiles = ExpandWildCards(fn_hyp_in,
            fn_hyp_in_list, MAX_NUM_INPUT_FILES)) < 1) {
        nll_puterr("ERROR: no matching .hyp files found.");
        return (-1);
    }
    if (numFiles >= MAX_NUM_INPUT_FILES) {
        sprintf(MsgStr, "WARNING: maximum number of event files exceeded, only first %d will be processed.", MAX_NUM_INPUT_FILES);
        nll_puterr(MsgStr);
    }


    nLocWritten = 0;
    nLocAccepted = 0;
    for (nHypo = 0; nHypo < numFiles; nHypo++) {

        // open hypocenter file
        // 20160621 AJL - bug (compiler warning) fix  sprintf(strstr(fn_hyp_in_list[nHypo], test_str), "\0");
        sprintf(strstr(fn_hyp_in_list[nHypo], test_str), "%c", '\0');
        sprintf(fn_hyp_in, "%s.hyp", fn_hyp_in_list[nHypo]);
        if (VERBOSE > 0) fprintf(stdout, "Opening: %s\n", fn_hyp_in);
        if ((fp_hyp_in = fopen(fn_hyp_in, "r")) == NULL) {
            nll_puterr2("ERROR: opening hypocenter file, ignoring event, file",
                    fn_hyp_in_list[nHypo]);
            continue;
        }

        while (1) {

            istat = GetHypLoc(fp_hyp_in, fn_hyp_in_list[nHypo], diffHypos + nLocWritten, Arrival + NumArrivals,
                    &nArrivals_read, 1, &locgrid, 0);
            if (istat == EOF) {
                fclose(fp_hyp_in);
                break;
            }
            if (VERBOSE > 1) fprintf(stdout, "Reading: location %d\n", nLocWritten);
            // 20180612 AJL - added
            if (strstr(MapProjStr[0], "GLOBAL") != NULL) {
                GeometryMode = MODE_GLOBAL;
            }

            if (strcmp(diffHypos[nLocWritten].locStat, "ABORTED") == 0) {
                //nll_puterr("WARNING: location ABORTED, ignoring event");
                if (VERBOSE > 0) fprintf(stdout, "Skip: diffHypos[nLocWritten].locStat, \"ABORTED\" %d\n", nLocWritten);
                continue;
            } else if (strcmp(diffHypos[nLocWritten].locStat, "REJECTED") == 0) {
                //nll_puterr("WARNING: location REJECTED, ignoring event");
                if (VERBOSE > 0) fprintf(stdout, "Skip: diffHypos[nLocWritten].locStat, \"REJECTED\" %d\n", nLocWritten);
                continue;
            }
            nLocAccepted++;


            diffHypos[nLocWritten].event_id = nLocWritten;
            weight = 0.0;
            for (i = NumArrivals; i < NumArrivals + nArrivals_read; i++) {
                Arrival[i].dd_event_id_1 = nLocWritten;
                weight = Arrival[i].weight > weight ? Arrival[i].weight : weight;
            }
            // normalize weight
            for (i = NumArrivals; i < NumArrivals + nArrivals_read; i++)
                Arrival[i].weight /= weight;

            WriteLocation(fp_hyp_out, diffHypos + nLocWritten,
                    Arrival, nArrivals_read, fn_hyp_out, 0, 0, 0, &locgrid, 0);
            nLocWritten++;

            // write end line and blank line
            fprintf(fp_hyp_out, "END_NLLOC\n\n");

            NumArrivals += nArrivals_read;

            if (nLocWritten >= MAX_NUM_DIFF_HYPOCENTERS) {
                sprintf(MsgStr, "WARNING: maximum number of input hypocenters read, only first %d will be processed.", MAX_NUM_DIFF_HYPOCENTERS);
                nll_puterr(MsgStr);
                break;
            }

        }

    }


    fclose(fp_hyp_out);

    // write message
    fprintf(stdout,
            "%d location files read, %d events accepted.\n", numFiles, nLocAccepted);
    fprintf(stdout, "%d locations written to ascii sumfile <%s>\n", nLocWritten, fn_hyp_out);

    // difference and write arrivals

    /* example:
    #     38542     38520
    NCCBW     3.430   3.430 1.0000 P
    NCCSP     2.850   2.840 1.0000 P
    NCNLN     8.950   8.950 1.0000 P
    NCCAI     3.440   3.420 1.0000 P
     */

    dd_event_id_1 = dd_event_id_2 = -1;
    nwritten = 0;

    for (narr1 = 0; narr1 < NumArrivals; narr1++) {

        fprintf(stdout, "Arr: %d/%d        \r", narr1, NumArrivals);
        if (VERBOSE > 1) fprintf(stdout, "\n");

        parr1 = Arrival + narr1;

        if (parr1->weight < SMALL_DOUBLE) {
            if (VERBOSE > 1) fprintf(stdout, "Skip: parr1->weight < SMALL_DOUBLE %d %f\n", narr1, parr1->weight);
            continue; // not used for location
        }

        for (narr2 = narr1 + 1; narr2 < NumArrivals; narr2++) {

            parr2 = Arrival + narr2;

            if (parr2->weight < SMALL_DOUBLE) {
                if (VERBOSE > 1) fprintf(stdout, "Skip: narr2->weight < SMALL_DOUBLE %d %f\n", narr2, parr2->weight);
                continue; // not used for location
            }

            if (parr1->dd_event_id_1 == parr2->dd_event_id_1) {
                if (VERBOSE > 1) fprintf(stdout, "Skip: parr1->dd_event_id_1 == parr2->dd_event_id_1 %d %d\n", narr1, narr2);
                continue; // same event
            }

            if (strcmp(parr1->label, parr2->label)) {
                if (VERBOSE > 1) fprintf(stdout, "Skip: strcmp(parr1->label, parr2->label) %d %d\n", narr1, narr2);
                continue; // not same station
            }

            if (strcmp(parr1->phase, parr2->phase)) {
                if (VERBOSE > 1) fprintf(stdout, "Skip: strcmp(parr1->phase, parr2->phase) %d %d\n", narr1, narr2);
                continue; // not same phase
            }

            phyp1 = &(diffHypos[parr1->dd_event_id_1]);
            phyp2 = &(diffHypos[parr2->dd_event_id_1]);

            // check distance between events
            if (event_dist_max > 0.0) {
                // 20180612 AJL  epi_dist = Dist3D(phyp1->x, phyp2->x, phyp1->y, phyp2->y, 0.0, 0.0); // epi
                epi_dist = Dist2D(phyp1->x, phyp2->x, phyp1->y, phyp2->y); // epi
                //phyp1->z, phyp2->z);  // hypo
                if (epi_dist > event_dist_max) {
                    if (VERBOSE > 1) fprintf(stdout, "Skip: epi_dist > event_dist_max %d %d\n", narr1, narr2);
                    continue;
                }
            }

            // check distance between station and events
            if (station_dist_max > 0.0) {
                if (parr1->dist * DEG2KM > station_dist_max) {
                    if (VERBOSE > 1) fprintf(stdout, "Skip: parr1->dist > station_dist_max %d\n", narr1);
                    continue;
                }
                if (parr2->dist * DEG2KM > station_dist_max) {
                    if (VERBOSE > 1) fprintf(stdout, "Skip: parr1->dist > station_dist_max %d\n", narr2);
                    continue;
                }
            }

            // check weight
            //printf("%d/%d  parr1->weight %f  parr2->weight %f -> ", narr1, narr2, parr1->weight, parr2->weight);
            weight = parr1->weight < parr2->weight ? parr1->weight : parr2->weight;
            if (weight < weight_min) {
                //printf("REJECT\n");
                if (VERBOSE > 1) fprintf(stdout, "Skip: weight < weight_min %d %f %d %f\n", narr1, parr1->weight, narr2, parr2->weight);
                continue;
            }
            //printf("Accept\n");
            if (VERBOSE > 1) fprintf(stdout, "Accept: %d w %f  %d w %f\n", narr1, parr1->weight, narr2, parr2->weight);

            // write event line if needed
            if (!(parr1->dd_event_id_1 == dd_event_id_1
                    && parr2->dd_event_id_1 == dd_event_id_2)) {
                dd_event_id_1 = parr1->dd_event_id_1;
                dd_event_id_2 = parr2->dd_event_id_1;
                fprintf(fp_diff_out, "# %8ld %8ld\n", dd_event_id_1, dd_event_id_2);
            }

            // write phase line
            fprintf(fp_diff_out, "%-7s %7.3lf %7.3lf %6.4lf %4s\n",
                    parr1->label,
                    parr1->pred_travel_time + parr1->residual,
                    parr2->pred_travel_time + parr2->residual,
                    weight, parr1->phase);

            fprintf(fp_xyz_out, "> GMT_LATLONDEPTH\n");
            fprintf(fp_xyz_out, "%lf %lf %lf\n", phyp1->dlat, phyp1->dlong, phyp1->depth);
            fprintf(fp_xyz_out, "%lf %lf %lf\n", phyp2->dlat, phyp2->dlong, phyp2->depth);

            nwritten++;

        }

    }
    fprintf(stdout, "\n");


    fclose(fp_xyz_out);
    fclose(fp_diff_out);

    fprintf(stdout, "%d arrivals differenced, %d dd written to ct file <%s>\n",
            NumArrivals, nwritten, fn_diff_out);

    return (0);

}



