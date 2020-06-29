/*
 * Copyright (C) 1999-2010 Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
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



/*   scat2latlon.c

        Program to convert binary event scatter files to lat/lon ascii files

 */

/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
161 Allee du Micocoulier, 06370 Mouans-Sartoux, France
tel: +33(0)493752502  e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


/*
        history:

        ver 01    08Jan2013    AJL    Original version

 */



#include "GridLib.h"


/* defines */

#define MAX_NUM_INPUT_FILES 32000


/* globals */



/* functions */

int ConvertScat(int argc, char** argv);



/*** program to sum event scatter files */

#define PNAME  "scat2latlon"

int main(int argc, char** argv) {

    int istat, narg;


    /* set program name */

    strcpy(prog_name, PNAME);


    /* check command line for correct usage */

    fprintf(stdout, "%s Arguments: ", prog_name);
    for (narg = 0; narg < argc; narg++)
        fprintf(stdout, "<%s> ", argv[narg]);
    fprintf(stdout, "\n");

    if (argc < 4) {
        nll_puterr("ERROR wrong number of command line arguments.");
        disp_usage(PNAME,
                "<decim_factor> <output_dir> <hyp_file_list>"
                "\n   decim_factor   decimation factor for skipping scatter samples (1+)"
                "\n   output_dir     existing directory for writing output ascii scatter files"
                "\n   hyp_file_list  path and filename to NLL *.hyp files to process, can be wildcarded"
                );
        exit(-1);
    }

    SetConstants();
    message_flag = 1;
    DispProgInfo();
    message_flag = 0;

    if ((istat = ConvertScat(argc, argv)) < 0) {
        nll_puterr("ERROR doing ADD process.");
        exit(-1);
    }



    exit(0);

}

int ConvertScat(int argc, char** argv) {

    int istat;
    int num_decim;
    int n, npt, npoints;
    long num_points_tot = 0L, num_points_read = 0L, num_points_written = 0L;

    char *pchr;
    char filename[FILENAME_MAX];
    char fn_scatter[FILENAME_MAX];
    char fn_root_out[FILENAME_MAX], fn_hypos_in[FILENAME_MAX], fn_scat_out[FILENAME_MAX];
    FILE *fp_hypo, *fp_scat_out, *fp_scat_in;
    float fdata[4], probmax = -VERY_LARGE_FLOAT;
    double dlat, dlon;

    int nFile, numFiles, nLocWritten, nLocAccepted;
    // make fn_hyp_in_list static to avoid stack overflow problem (F Tilmann) //
    static char fn_hyp_in_list[MAX_NUM_INPUT_FILES][FILENAME_MAX];
    char test_str[10];

    GridDesc locgrid;
    HypoDesc Hypo;
    char map_trans_str[2 * MAXLINE];
    int n_proj = 0;


    strcpy(test_str, ".hyp");

    // get command line parameters
    sscanf(argv[1], "%d", &num_decim);
    if (num_decim < 1)
        fprintf(stdout,
            "  Decimation < 1, no scatter points will be written to ouptut.\n");
    else
        fprintf(stdout, "  Decimation: %d\n", num_decim);

    sprintf(fn_root_out, "%s", argv[2]);

    // check if input file list given on command line
    if (argc > 4) { // assume args are list of hyp file to process
        numFiles = 0;
        for (n = 3; n < argc && numFiles < MAX_NUM_INPUT_FILES; n++) {
            strcpy(fn_hyp_in_list[numFiles], argv[n]);
            numFiles++;
        }
    } else { // assume argv[3] is wildcarded input file root
        strcpy(fn_hypos_in, argv[3]);
        // make sure input file name ends in test_str
        strcat(fn_hypos_in, test_str);
        // check for wildcards in input file name
        if ((numFiles = ExpandWildCards(fn_hypos_in, fn_hyp_in_list, MAX_NUM_INPUT_FILES)) < 1) {
            nll_puterr("ERROR: no matching .hyp files found or other error.");
            return (-1);
        }
    }

    if (numFiles >= MAX_NUM_INPUT_FILES) {
        sprintf(MsgStr, "WARNING: maximum number of event files reached, only first %d will be processed.", MAX_NUM_INPUT_FILES);
        nll_puterr(MsgStr);
    }


    nLocWritten = 0;
    nLocAccepted = 0;
    for (nFile = 0; nFile < numFiles; nFile++) {

        fprintf(OUT_LEVEL_1, "Converting scatter file for location <%s>.\n", fn_hyp_in_list[nFile]);

        // open hypocenter file
        //sprintf(strstr(fn_hyp_in_list[nFile], test_str), "\0");
        if ((fp_hypo = fopen(fn_hyp_in_list[nFile], "r")) == NULL) {
            nll_puterr2("ERROR: opening hypocenter file, ignoring event, file",
                    fn_hyp_in_list[nFile]);
            continue;
        }

        // open scatter file
        char *cptr = strrchr(fn_hyp_in_list[nFile], '/');
        if (cptr != NULL)
            strcpy(filename, cptr + 1);
        else
            strcpy(filename, fn_hyp_in_list[nFile]);
        sprintf(fn_scat_out, "%s/%s.scat.xyz", fn_root_out, filename);
        if ((fp_scat_out = fopen(fn_scat_out, "w")) == NULL) {
            nll_puterr2("ERROR: opening scatter output file:", fn_scat_out);
            return (-1);
        }
        fprintf(fp_scat_out, "> GMT_LATLONDEPTH\n");
        fprintf(fp_scat_out, "> lat(deg) lon(deg) depth(km) likelihood\n");

        num_points_tot = 0L;
        num_points_read = 0L;
        num_points_written = 0L;

        // loop over events in file

        while (1) {

            istat = GetHypLoc(fp_hypo, NULL, &Hypo, Arrival, &NumArrivals, 1, &locgrid, n_proj);
            if (istat == EOF) {
                break;
            }
            if (istat < 0) {
                nll_puterr2("ERROR: reading hypocenter file, ignoring event, file",
                        fn_hyp_in_list[nFile]);
                break;
            }

            if (strcmp(Hypo.locStat, "ABORTED") == 0) {
                nll_puterr("WARNING: location ABORTED, ignoring event");
                continue;
            } else if (strcmp(Hypo.locStat, "REJECTED") == 0) {
                nll_puterr("WARNING: location REJECTED, ignoring event");
                continue;
            }
            nLocAccepted++;

            if (num_decim > 0) {

                // set map transform
                printf("MapProjStr[n_proj] <%s>\n", MapProjStr[n_proj]);
                projection_str2transform_str(map_trans_str, MapProjStr[n_proj]);
                printf("map_trans_str <%s>\n", map_trans_str);
                if ((istat = get_transform(n_proj, map_trans_str)) < 0) {
                    nll_puterr2("ERROR: reading transformation parameters:", map_trans_str);
                    continue;
                }

                // open scatter file
                strcpy(fn_scatter, fn_hyp_in_list[nFile]);
                pchr = strstr(fn_scatter, test_str);
                if (pchr != NULL)
                    *pchr = '\0';
                strcat(fn_scatter, ".scat\0");
                if ((fp_scat_in = fopen(fn_scatter, "r")) == NULL) {
                    nll_puterr2("ERROR: opening scatter file", fn_scatter);
                    continue;
                }

                fprintf(fp_scat_out, "> %s\n", Hypo.fileroot);


                // read header record
                fseek(fp_scat_in, 0, SEEK_SET);
                fread(&npoints, sizeof (int), 1, fp_scat_in);

                // skip header record
                fseek(fp_scat_in, 4 * sizeof (float), SEEK_SET);

                // copy date records
                /*fprintf(stdout, "  Summing %d samples...\n", npoints);*/
                for (npt = 0; npt < npoints; npt += num_decim) {

                    if ((istat = fread(fdata, sizeof (float), 4, fp_scat_in)) != 4) {
                        sprintf(MsgStr, "ERROR: freed = %d != 4!!! (%d/%d)", istat, npt, npoints);
                        nll_puterr(MsgStr);
                    }
                    num_points_read++;
                    fseek(fp_scat_in, (num_decim - 1) * 4 * sizeof (float),
                            SEEK_CUR);

                    if (fdata[3] > probmax)
                        probmax = fdata[3];

                    istat = rect2latlon(n_proj, fdata[0], fdata[1], &dlat, &dlon);
                    //fprintf(fp_scat_out, "%9.4lf %9.4lf %9.4lf  0 %9.2le\n", dlat, dlon, fdata[2], fdata[3]);
                    fprintf(fp_scat_out, "%f %f %f  0.001  %e\n", dlat, dlon, fdata[2], fdata[3]);

                    num_points_written++;
                }

                num_points_tot += npoints;

                fclose(fp_scat_in);
            }


        }

        fclose(fp_scat_out);
        fclose(fp_hypo);

        fprintf(stdout,
                "   %ld samples total in input scatter files.\n", num_points_tot);
        fprintf(stdout,
                "   %ld decimated samples read.\n", num_points_read);
        fprintf(stdout, "   %ld samples written to binary sumfile <%s>\n",
                num_points_written, fn_scat_out);

    }



    /* write message */
    fprintf(stdout,
            "%d location files read, %d accepted.\n",
            numFiles, nLocAccepted);


    return (0);

}


