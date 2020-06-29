/*
 * Copyright (C) 1999 Anthony Lomax <lomax@faille.unice.fr>
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


/*   LocSum.c

        Program to sum (concatinate) location and binary event scatter files

 */

/*------------------------------------------------------------/ */
/* Anthony Lomax           | email: lomax@faille.unice.fr     / */
/* UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax / */
/* 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        / */
/* 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        / */
/*------------------------------------------------------------/ */


/*
        history:

        ver 01    04Nov1997  AJL  Original version

                   Unregistered bug fixes Jan 2006 Frederik Tilmann
                        MAX_NUM_INPUT_FILES reduced to 32000 (larger values result in segmentation violation on jumping into subroutine Sum_Locations)
                   20060213 AJL MAX_NUM_INPUT_FILES reduced to 5000
                   20110727 AJL Incorporated changes to MAX_NUM_INPUT_FILES etc. proposed by Frederik Tilmann :
                       Unregistered bug fixes Frederik Tilmann
                       2011-07-27
                        // make fn_hyp_in_list static to avoid stack overflow problem (F Tilmann) - previously //
                        // Add fclose(fp_hypo) to for loop over all events as otherwise run out of file pointers after about
                            1000 events

.........1.........2.........3.........4.........5.........6.........7.........8

 */



#include "GridLib.h"


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
#define MAX_NUM_INPUT_FILES 32000


/* globals */



/* functions */

int SumLocations(int argc, char** argv);



/*** program to sum event scatter files */

#define PNAME  "LocSum"

int main(int argc, char** argv) {

    int istat, narg;


    /* set program name */

    strcpy(prog_name, PNAME);


    /* check command line for correct usage */

    fprintf(stdout, "%s Arguments: ", prog_name);
    for (narg = 0; narg < argc; narg++)
        fprintf(stdout, "<%s> ", argv[narg]);
    fprintf(stdout, "\n");

    if (argc < 5) {
        nll_puterr("ERROR wrong number of command line arguments.");
        disp_usage(PNAME,
                "<size_gridfile> <scat_decim> <output_file_root> <add_file_list> [Len3Max [ProbMin [RMSMax [NRdgsMin [GapMax [latMin,latMax,longMin,longMax]]]]]]");
        exit(-1);
    }

    SetConstants();
    message_flag = 1;
    DispProgInfo();
    message_flag = 0;

    if ((istat = SumLocations(argc, argv)) < 0) {
        nll_puterr("ERROR doing ADD process.");
        exit(-1);
    }



    exit(0);

}

int SumLocations(int argc, char** argv) {

    int istat;
    int num_decim;
    int npt, npoints;
    long num_points_tot = 0L, num_points_read = 0L, num_points_written = 0L;

    char *pchr;
    char sys_string[FILENAME_MAX];
    char filename[FILENAME_MAX];
    char fn_grid_size[FILENAME_MAX];
    char fn_scatter[FILENAME_MAX];
    char fn_hyp_scat_out[FILENAME_MAX];
    char fn_root_out[FILENAME_MAX], fn_hypos_in[FILENAME_MAX],
            fn_scat_out[FILENAME_MAX];
    FILE *fp_hypo, *fp_dummy, *fp_hyp_scat_out, *fp_scat_out,
            *fp_scat_in, *fp_grid, *fp_hdr;
    float fdata[4], probmax = -VERY_LARGE_FLOAT;

    int nFile, numFiles, nLocWritten, nLocAccepted;
    //char fn_hyp_in_list[MAX_NUM_INPUT_FILES][FILENAME_MAX];
    // make fn_hyp_in_list static to avoid stack overflow problem (F Tilmann) //
    static char fn_hyp_in_list[MAX_NUM_INPUT_FILES][FILENAME_MAX];
    char test_str[10];

    double xmin, xmax, ymin, ymax, zmin, zmax;

    double Len3Max, ProbMin, RMSMax;
    int NRdgsMin, GapMax;
    double Len3Mean = 0.0, ProbMean = 0.0, RMSMean = 0.0;
    double NRdgsMean = 0.0, GapMean = 0.0;
    int Len3Reject = 0, ProbReject = 0, RMSReject = 0;
    int NRdgsReject = 0, GapReject = 0, CutReject = 0;
    int iReject;

    double latMin, latMax, longMin, longMax;
    int icut = 0;

    GridDesc Grid, locgrid;
    SourceDesc* Srce = NULL;
    HypoDesc Hypo;


    strcpy(test_str, ".hyp");

    /* get command line parameters */
    strcpy(fn_grid_size, argv[1]);
    sscanf(argv[2], "%d", &num_decim);
    if (num_decim < 1)
        fprintf(stdout,
            "  Scatter Decimation < 1, no scatter points will be written to ouptut.\n");
    else
        fprintf(stdout, "  Decimation: %d\n", num_decim);

    sprintf(fn_root_out, "%s", argv[3]);
    strcpy(fn_hypos_in, argv[4]);

    Len3Max = 1.0e6;
    if (argc > 5) {
        sscanf(argv[5], "%lf", &Len3Max);
    }
    fprintf(stdout, "  Ellipsoid Len3 Maximum: %lf\n", Len3Max);
    ProbMin = -1.0;
    if (argc > 6) {
        sscanf(argv[6], "%lf", &ProbMin);
    }
    fprintf(stdout, "  Probablity Minimum: %lf\n", ProbMin);
    RMSMax = 1.0e6;
    if (argc > 7) {
        sscanf(argv[7], "%lf", &RMSMax);
    }
    fprintf(stdout, "  RMS Maximum: %lf\n", RMSMax);
    NRdgsMin = 0;
    if (argc > 8) {
        sscanf(argv[8], "%d", &NRdgsMin);
    }
    fprintf(stdout, "  Num Readings Minimum: %d\n", NRdgsMin);
    GapMax = 360;
    if (argc > 9) {
        sscanf(argv[9], "%d", &GapMax);
    }
    fprintf(stdout, "  Gap Maximum: %d\n", GapMax);
    if (argc > 10) {
        sscanf(argv[10], "%lf,%lf,%lf,%lf", &latMin, &latMax, &longMin, &longMax);
        icut = 1;
        fprintf(stdout, "  Geog Cut Limits: Lat: %f -> %f, Long: %f -> %f\n", latMin, latMax, longMin, longMax);
    } else {
        icut = 0;
        fprintf(stdout, "  No Geog Cut\n");
    }


    /* duplicate size grid files to make dummy output grid files */

    sprintf(sys_string, "cp %s.hdr %s.hdr", fn_grid_size, fn_root_out);
    if ((istat = system(sys_string)) != 0) {
        sprintf(MsgStr, "system return value = %d", istat);
        nll_puterr2("ERROR: copying header file", MsgStr);
        return (-1);
    }
    sprintf(filename, "%s.buf", fn_root_out);
    fp_dummy = fopen(filename, "w");
    fclose(fp_dummy);

    /* read grid header file and get grid bounds */

    Grid.iSwapBytes = 0;
    if ((istat = OpenGrid3dFile(fn_root_out, &fp_grid, &fp_hdr, &Grid, "scat", Srce, Grid.iSwapBytes)) < 0) {
        sprintf(MsgStr, "%s.hdr", fn_root_out);
        nll_puterr2("ERROR: open grid header file", MsgStr);
        return (-1);
    }
    xmin = Grid.origx;
    xmax = xmin + (double) (Grid.numx - 1) * Grid.dx;
    ymin = Grid.origy;
    ymax = ymin + (double) (Grid.numy - 1) * Grid.dy;
    zmin = Grid.origz;
    zmax = zmin + (double) (Grid.numz - 1) * Grid.dz;



    /* open ascii hypocenter/scatter file */

    sprintf(fn_hyp_scat_out, "%s.hyp", fn_root_out);
    if ((fp_hyp_scat_out = fopen(fn_hyp_scat_out, "w")) == NULL) {
        nll_puterr("ERROR: opening scatter ascii output file.");
        return (-1);
    }

    /* open scatter file */

    sprintf(fn_scat_out, "%s.scat", fn_root_out);
    if ((fp_scat_out = fopen(fn_scat_out, "w")) == NULL) {
        nll_puterr("ERROR: opening scatter output file.");
        return (-1);
    }
    /* skip header record */
    fseek(fp_scat_out, 4 * sizeof (float), SEEK_SET);


    /* sum requested grid files into output grid */

    /* check for wildcards in input file name */
    strcat(fn_hypos_in, test_str);
    if ((numFiles = ExpandWildCards(fn_hypos_in, fn_hyp_in_list, MAX_NUM_INPUT_FILES)) < 1) {
        nll_puterr("ERROR: no matching .hyp files found.");
        return (-1);
    }
    if (numFiles >= MAX_NUM_INPUT_FILES) {
        sprintf(MsgStr, "WARNING: maximum number of event files exceeded, only first %d will be processed.", MAX_NUM_INPUT_FILES);
        nll_puterr(MsgStr);
    }


    nLocWritten = 0;
    nLocAccepted = 0;
    for (nFile = 0; nFile < numFiles; nFile++) {

        fprintf(OUT_LEVEL_1, "Adding location <%s>.\n", fn_hyp_in_list[nFile]);

        /* open hypocenter file */
        //sprintf(strstr(fn_hyp_in_list[nFile], test_str), "\0");
        if ((fp_hypo = fopen(fn_hyp_in_list[nFile], "r")) == NULL) {
            nll_puterr2("ERROR: opening hypocenter file, ignoring event, file",
                    fn_hyp_in_list[nFile]);
            continue;
        }

        // loop over events in file

        while (1) {

            istat = GetHypLoc(fp_hypo, NULL, &Hypo, Arrival, &NumArrivals, 1, &locgrid, 0);
            if (istat == EOF) {
                break;
            }
            if (istat < 0) {
                nll_puterr2("ERROR: reading hypocenter file, ignoring event, file",
                        fn_hyp_in_list[nFile]);
                break;
            }

            if (strcmp(Hypo.locStat, "ABORTED") == 0) {
                //nll_puterr("WARNING: location ABORTED, ignoring event");
                continue;
            } else if (strcmp(Hypo.locStat, "REJECTED") == 0) {
                //nll_puterr("WARNING: location REJECTED, ignoring event");
                continue;
            }
            nLocAccepted++;
            Len3Mean += Hypo.ellipsoid.len3;
            ProbMean += Hypo.probmax;
            RMSMean += Hypo.rms;
            NRdgsMean += (double) Hypo.nreadings;
            GapMean += Hypo.gap;
            iReject = 0;
            if (Hypo.ellipsoid.len3 > Len3Max) {
                //nll_puterr("WARNING: location ellipsoid Len3 is greater than Len3Max, ignoring event");
                Len3Reject++;
                iReject = 1;
            }
            if (Hypo.probmax < ProbMin) {
                //nll_puterr("WARNING: location Prob max is less than ProbMin, ignoring event");
                ProbReject++;
                iReject = 1;
            }
            if (Hypo.rms > RMSMax) {
                //nll_puterr("WARNING: location RMS is Greater than RMSMax, ignoring event");
                RMSReject++;
                iReject = 1;
            }
            if (Hypo.nreadings < NRdgsMin) {
                //nll_puterr("WARNING: location num readings is less than NRdgsMin, ignoring event");
                NRdgsReject++;
                iReject = 1;
            }
            if (Hypo.gap > GapMax) {
                //nll_puterr("WARNING: location gap is greater than GapMax, ignoring event");
                GapReject++;
                iReject = 1;
            }
            if (icut && (Hypo.dlat < latMin || Hypo.dlat > latMax || Hypo.dlong < longMin || Hypo.dlong > longMax)) {
                CutReject++;
                iReject = 1;
            }

            if (iReject)
                continue;


            PhaseFormat = FORMAT_PHASE_2; // 20110105 AJL - to allow long station names
            WriteLocation(fp_hyp_scat_out, &Hypo, Arrival, NumArrivals, fn_hyp_scat_out, 1, 0, 0, &locgrid, 0);

            nLocWritten++;

            if (num_decim > 0) {

                /* open scatter file */
                strcpy(fn_scatter, fn_hyp_in_list[nFile]);
                pchr = strstr(fn_scatter, test_str);
                if (pchr != NULL)
                    *pchr = '\0';
                strcat(fn_scatter, ".scat\0");
                if ((fp_scat_in = fopen(fn_scatter, "r")) == NULL) {
                    nll_puterr2("ERROR: opening scatter file", fn_scatter);
                    fprintf(fp_hyp_scat_out, "SCATTER Nsamples %d\n", 0);
                    fprintf(fp_hyp_scat_out, "END_SCATTER\n\n");
                    continue;
                }


                /* read header record */
                fseek(fp_scat_in, 0, SEEK_SET);
                fread(&npoints, sizeof (int), 1, fp_scat_in);

                fprintf(fp_hyp_scat_out, "SCATTER Nsamples %d\n", npoints / num_decim);

                /* skip header record */
                fseek(fp_scat_in, 4 * sizeof (float), SEEK_SET);

                /* copy date records */
                /*fprintf(stdout, "  Summing %d samples...\n", npoints);*/
                for (npt = 0; npt < npoints; npt += num_decim) {

                    if ((istat = fread(fdata, sizeof (float), 4, fp_scat_in)) != 4) {
                        sprintf(MsgStr, "ERROR: freed = %d != 4!!! (%d/%d)",
                                istat, npt, npoints);
                        nll_puterr(MsgStr);
                    }
                    num_points_read++;
                    fseek(fp_scat_in, (num_decim - 1) * 4 * sizeof (float),
                            SEEK_CUR);

                    /* clip - check that sample is within grid */
                    /*if (fdata[0] < xmin || fdata[0] > xmax
                            || fdata[1] < ymin || fdata[1] > ymax
                            || fdata[2] < zmin || fdata[2] > zmax)
                            continue;*/

                    fwrite(fdata, sizeof (float), 4, fp_scat_out);
                    num_points_written++;

                    if (fdata[3] > probmax)
                        probmax = fdata[3];

                    fprintf(fp_hyp_scat_out, "%9.4lf %9.4lf %9.4lf %9.2le\n",
                            fdata[0], fdata[1], fdata[2], fdata[3]);

                }

                num_points_tot += npoints;

                fprintf(fp_hyp_scat_out, "END_SCATTER\n");

                fclose(fp_scat_in);
            }


            /* write end line and blank line */
            fprintf(fp_hyp_scat_out, "END_NLLOC\n\n");



        }
        fclose(fp_hypo);

    }


    /* write header informaion */
    fseek(fp_scat_out, 0, SEEK_SET);
    fwrite(&num_points_written, sizeof (int), 1, fp_scat_out);
    fdata[0] = (float) probmax;
    fwrite(fdata, sizeof (float), 1, fp_scat_out);

    fclose(fp_scat_out);
    fclose(fp_hyp_scat_out);

    /* write message */
    fprintf(stdout,
            "%d location files read, %d accepted.\n",
            numFiles, nLocAccepted);
    fprintf(stdout,
            "Len3 Mean: %lf, Reject %d\n",
            Len3Mean / (double) nLocAccepted, Len3Reject);
    fprintf(stdout,
            "Prob Mean: %lf, Reject %d\n",
            ProbMean / (double) nLocAccepted, ProbReject);
    fprintf(stdout,
            "RMS Mean: %lf, Reject %d\n",
            RMSMean / (double) nLocAccepted, RMSReject);
    fprintf(stdout,
            "NRdgs Mean: %lf, Reject %d\n",
            NRdgsMean / (double) nLocAccepted, NRdgsReject);
    fprintf(stdout,
            "Gap Mean: %lf, Reject %d\n",
            GapMean / (double) nLocAccepted, GapReject);
    fprintf(stdout,
            "Cut Reject %d\n", CutReject);
    fprintf(stdout,
            "%ld samples total in input scatter files.\n", num_points_tot);
    fprintf(stdout,
            "%ld decimated samples read.\n", num_points_read);
    fprintf(stdout, "%ld samples written to binary sumfile <%s>\n",
            num_points_written, fn_scat_out);
    fprintf(stdout,
            "%ld samples and %d locations written to ascii sumfile <%s>\n",
            num_points_written, nLocWritten, fn_hyp_scat_out);


    return (0);

}



/*------------------------------------------------------------/ */
/* Anthony Lomax           | email: lomax@faille.unice.fr     / */
/* UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax / */
/* 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        / */
/* 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        / */
/*------------------------------------------------------------/ */

