/*
 * Copyright (C) 2016- Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
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



/*   MergeFmamp.c

        Program to merge fmamp mechanism results into  NonLinLoc hyp files

 */



/*
        history:

        ver 01    20160920  AJL  Original version adapted from MergeFmamp.c


.........1.........2.........3.........4.........5.........6.........7.........8

 */



#include "../src/include/GridLib.h"


/* defines */

#define MAX_NUM_INPUT_FILES 4096


/* globals */

int DEBUG = 0;

/* functions */

int MergeFmamp(int argc, char** argv);
int compareTimes(int iy1, int iy2, int im1, int im2, int id1, int id2, int ih1, int ih2,
        int imn1, int imn2, double s1, double s2, double tolerance);


/*** program to merge fmamp results into NLL hyp file */

#define PNAME  "MergeFmamp"

/** function to read fpfit summary record to HypoDesc structure */

int ReadFmampSum(FILE *fp_in, HypoDesc *phypo, char *phypo_fmamp_quality) {

    static char line[MAXLINE_LONG];

    // currently unused fields
    long unique_id;
    double errh;
    double errz;
    double dist_max;
    double ampAttenPower;
    char mag_type[64];
    char variant_name[64];
    int nreadings;
    double stat_dist_ratio;
    double sum_misfit_weight;
    double mean_dist_P;
    double mean_dist_T;

    // NLL FocalMech uses FPFIT conventions for fault planes:
    //	double dipDir; // Dip direction (downdip azimuth in degrees,clockwise from north)
    //	double dipAng; // Dip angle in degrees down from horizontal
    //	double rake; // Rake in degrees: 0=left lateral, 90=reverse, +-180=right lateral, -90=normal
    // FMAMP uses Aki&Richards conventions for fault planes:
    //  double strike;  // Strike direction (along strike azimuth in degrees,clockwise from north) (0->360deg))
    //  double dip;     // Dip angle in degrees down from horizontal (0->90deg)
    //  double rake;    // Rake in degrees (-180->180deg) : 0=left lateral, 90=reverse, +-180=right lateral, -90=normal
    //
    // so must convert fmamp strike to NLL FocalMech dipDir:
    // dipDir = strike +90
    double strike;


    // read next line
    char *cstat = fgets(line, MAXLINE_LONG, fp_in);
    if (cstat == NULL)
        return (EOF);

    int istat = sscanf(line, "%ld %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %d %lf %lf %lf %lf %lf %lf %s %s %lf %lf %lf %d %d %lf %lf %lf %lf %lf %c ",
            // event information
            &unique_id,
            &phypo->year, &phypo->month, &phypo->day, &phypo->hour, &phypo->min, &phypo->sec, &phypo->rms,
            &phypo->dlat, &phypo->dlong, &errh, &phypo->depth, &errz,
            &phypo->nreadings, &phypo->dist, &dist_max, &phypo->gap, &phypo->gap_secondary,
            &ampAttenPower, &phypo->amp_mag, mag_type,
            // mechanism information
            variant_name,
            &strike, &phypo->focMech.dipAng, &phypo->focMech.rake,
            &nreadings, &phypo->focMech.nObs,
            &stat_dist_ratio,
            &phypo->focMech.misfit, &sum_misfit_weight,
            &mean_dist_P, &mean_dist_T,
            phypo_fmamp_quality
            );
    phypo->focMech.dipDir = strike + 90;
    if (phypo->focMech.dipDir >= 360.0) {
        phypo->focMech.dipDir -= 360.0;
    }

    if (istat != 33) {
        return (-2);
    }
    return (istat);



}

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
                "<NLLhyp_files [in quotes if wildcards used]> <fmamp_file> <out_hyp_directory> time_tolerance [MechQualityAccept [NObsUsedMin [GapMax [MechMisfitMax]]]]");
        exit(-1);
    }

    SetConstants();
    message_flag = 1;
    DispProgInfo();
    message_flag = 0;

    if ((istat = MergeFmamp(argc, argv)) < 0) {
        nll_puterr("ERROR doing fmamp merge process.");
        exit(-1);
    }



    exit(0);

}

int MergeFmamp(int argc, char** argv) {

    int istat1, istat2;
    int matched;
    int numNLLFiles;
    int nHypMerged, nFmampRead, nNLLRead, nNLLFilesRead;

    char fn_hyp_in_list[MAX_NUM_INPUT_FILES][FILENAME_MAX];

    char fn_nllhyp_in[FILENAME_MAX];
    char fn_fmamp_in[FILENAME_MAX];
    char fdir_hyp_out[FILENAME_MAX], fn_out[FILENAME_MAX], fn_out_fopen[FILENAME_MAX], filename[FILENAME_MAX];
    FILE *fp_fmamp_in;


    double time_tolerance;
    double idiff;


    GridDesc locgrid;
    HypoDesc HypoFmamp, HypoNLL;

    char *pchr;
    char test_str[] = ".hyp";

    /* get command line parameters */
    strcpy(fn_nllhyp_in, argv[1]);
    strcpy(fn_fmamp_in, argv[2]);
    strcpy(fdir_hyp_out, argv[3]);
    sscanf(argv[4], "%lf", &time_tolerance);

    fprintf(stdout, "DEBUG: argc %d\n", argc);
    char MechQualityAccept[64] = "";
    if (argc > 5) {
        sscanf(argv[5], "%s", MechQualityAccept);
    }
    fprintf(stdout, "  MechQualityAccept: %s\n", MechQualityAccept);
    int NObsUsedMin = 0;
    if (argc > 6) {
        sscanf(argv[6], "%d", &NObsUsedMin);
    }
    fprintf(stdout, "  NObsUsedMin mechanism: %d\n", NObsUsedMin);
    int GapMax = 999;
    if (argc > 7) {
        sscanf(argv[7], "%d", &GapMax);
    }
    fprintf(stdout, "  Gap Maximum hypo: %d\n", GapMax);
    double MechMisfitMax = 1.0e6;
    if (argc > 8) {
        sscanf(argv[8], "%lf", &MechMisfitMax);
    }
    fprintf(stdout, "  MechMisfitMax: %lf\n", MechMisfitMax);


    /* open fmamp summary input files */

    if ((fp_fmamp_in = fopen(fn_fmamp_in, "r")) == NULL) {
        nll_puterr2("ERROR: opening fmamp summary file: ", fn_fmamp_in);
        return (-1);
    }


    /* check for wildcards in input file name */
    if ((pchr = strstr(fn_nllhyp_in, test_str)) == NULL) {
        strcat(fn_nllhyp_in, test_str);
    }
    if ((numNLLFiles = ExpandWildCards(fn_nllhyp_in, fn_hyp_in_list, MAX_NUM_INPUT_FILES)) < 1) {
        nll_puterr("ERROR: no matching .hyp files found.");
        return (-1);
    }
    if (numNLLFiles >= MAX_NUM_INPUT_FILES) {
        sprintf(MsgStr, "WARNING: maximum number of event files exceeded, only first %d will be processed.", MAX_NUM_INPUT_FILES);
        nll_puterr(MsgStr);
    }


    /* merge cooresponding hypocenters */

    char hypo_fmamp_quality;

    nHypMerged = 0;
    nFmampRead = 0;
    nNLLRead = 0;
    nNLLFilesRead = 0;
    FILE *fpHypo = NULL;
    FILE *fp_out = NULL;

    // read fmamp header
    if ((istat2 = ReadFmampSum(fp_fmamp_in, &HypoFmamp, &hypo_fmamp_quality)) == EOF) {
        nll_puterr("ERROR reading fmamp summary file header.");
        return (-1);
    }

    while (1) {

        // read next fmamp hypocenter
        if ((istat2 = ReadFmampSum(fp_fmamp_in, &HypoFmamp, &hypo_fmamp_quality)) == EOF) {
            nll_puterr2("INFO: EOF reading fmamp summary file", fn_fmamp_in);
            break;
        } else if (istat2 < 0) {
            nll_puterr2("ERROR: reading fmamp summary file", fn_fmamp_in);
            continue;
        }
        nFmampRead++;

        while (istat2 >= 0 && istat2 != EOF && nNLLFilesRead < numNLLFiles) {
            //printf("DEBUG: istat2 %d >= 0 && istat2 != EOF && nNLLFilesRead %d < numNLLFiles %d, nNLLRead %d\n", istat2, nNLLFilesRead, numNLLFiles, nNLLRead);

            // read next NLL hypo
            while (1) {
                if (fpHypo == NULL) {
                    if ((fpHypo = fopen(fn_hyp_in_list[nNLLFilesRead], "r")) == NULL) {
                        nll_puterr2("ERROR: opening hypocenter file", fn_hyp_in_list[nNLLFilesRead]);
                        return (-1);
                    }
                    if ((pchr = strstr(fn_hyp_in_list[nNLLFilesRead], test_str)) != NULL)
                        sprintf(pchr, "");
                    //
                    strcpy(fn_out, fdir_hyp_out);
                    strcat(fn_out, "/");
                    if (strrchr(fn_hyp_in_list[nNLLFilesRead], '/') != NULL)
                        strcpy(filename, strrchr(fn_hyp_in_list[nNLLFilesRead], '/'));
                    else
                        strcpy(filename, fn_hyp_in_list[nNLLFilesRead]);
                    strcat(fn_out, filename);
                    strcpy(fn_out_fopen, fn_out);
                    strcat(fn_out_fopen, ".hyp");
                    if ((fp_out = fopen(fn_out_fopen, "w")) == NULL) {
                        nll_puterr2("ERROR: opening output file", fn_out);
                        return (-1);
                    }
                }
                istat1 = GetHypLoc(fpHypo, fn_hyp_in_list[nNLLFilesRead], &HypoNLL, Arrival,
                        &NumArrivals, 1, &locgrid, 0);
                //printf("DEBUG: fpHypo %ld, istat1 %d == EOF &&  && nNLLFilesRead %d < numNLLFiles - 1 %d, nNLLRead %d\n", fpHypo, istat2, nNLLFilesRead, numNLLFiles - 1, nNLLRead);
                if (istat1 == EOF && nNLLFilesRead < numNLLFiles - 1) {
                    // try next file
                    fclose(fpHypo);
                    fpHypo = NULL;
                    fclose(fp_out);
                    fp_out = NULL;
                    nNLLFilesRead++;
                    continue;
                }
                break;
            }
            nNLLRead++;


            if (DEBUG) {
                printf("\nFmamp:\n");
                WriteLocation(stdout, &HypoFmamp, Arrival, 0, fn_nllhyp_in, 0, 0, 1, &locgrid, 0);
                printf("NLL (%s):\n", fn_hyp_in_list[nNLLFilesRead - 1]);
                WriteLocation(stdout, &HypoNLL, Arrival, 0, fn_fmamp_in, 0, 0, 1, &locgrid, 0);
            }

            if (istat1 >= 0 && istat1 != EOF) {

                matched = 0;
                while (1) {

                    idiff = compareTimes(HypoFmamp.year, HypoNLL.year, HypoFmamp.month, HypoNLL.month,
                            HypoFmamp.day, HypoNLL.day, HypoFmamp.hour, HypoNLL.hour,
                            HypoFmamp.min, HypoNLL.min, HypoFmamp.sec, HypoNLL.sec,
                            time_tolerance);

                    if (idiff == 0) { // match
                        matched = 1;
                        break;
                    } else if (idiff < 0) {
                        // fmamp hyp earlier, read next fmamp hyp
                        if ((istat2 = ReadFmampSum(fp_fmamp_in, &HypoFmamp, &hypo_fmamp_quality)) == EOF)
                            break;
                        else if (istat2 < 0) {
                            nll_puterr2("ERROR: reading fmamp summary file", fn_fmamp_in);
                            continue;
                        }
                        nFmampRead++;
                    } else {
                        // fmamp hyp later, read next NLL hyp
                        break;
                    }
                }
                // 20200818 AJL - write all hypos to output, event if no fmamp match
                //if (!matched)
                //continue;

                // check optional limits

                if (strlen(MechQualityAccept) > 0 && strchr(MechQualityAccept, hypo_fmamp_quality) == NULL) {
                    if (DEBUG) {
                        sprintf(MsgStr, "WARNING: location quality %c is not contained in MechQualityAccept %s, ignoring event", hypo_fmamp_quality, MechQualityAccept);
                        nll_puterr(MsgStr);
                    }
                    continue;
                } else if (HypoFmamp.focMech.nObs < NObsUsedMin) {
                    if (DEBUG) nll_puterr(
                            "WARNING: mechanims num readings is less than NObsUsedMin, ignoring event");
                    continue;
                } else if (HypoFmamp.gap > GapMax) {
                    if (DEBUG) nll_puterr(
                            "WARNING: location gap is greater than GapMax, ignoring event");
                    continue;
                } else if (HypoFmamp.focMech.misfit > MechMisfitMax) {
                    if (DEBUG) nll_puterr(
                            "WARNING: solution misfit is greater than MFMax, ignoring event");
                    continue;
                }

                // merge fmamp info into NLL hyp
                // 20200818 AJL - write all hypos to output, event if no fmamp match
                if (matched) {
                    HypoNLL.focMech = HypoFmamp.focMech;
                }

                // output merged hyp
                WriteGrid3dHdr(&locgrid, NULL, fn_out, NULL);
                WriteLocation(fp_out, &HypoNLL, Arrival, NumArrivals, fn_out, 1, 1, 0, &locgrid, 0);

                nHypMerged++;
                //printf("MERGED!!!!!!!!!!! %d\n", nHypMerged);

            }

            if (matched) {
                break; // next fmamp hypo
            }

        }

    }
    fclose(fpHypo);
    fclose(fp_out);

    printf("%d fmamp and %d NLL hypocenters read, %d merged, output in: %s\n:", nFmampRead, nNLLRead, nHypMerged, fdir_hyp_out);

    return (0);

}

int compareTimes(int iy1, int iy2, int im1, int im2, int id1, int id2, int ih1, int ih2,
        int imn1, int imn2, double s1, double s2, double tolerance) {
    double time1, time2;

    if (iy1 != iy2)
        return (iy1 < iy2 ? -1 : 1);
    if (im1 != im2)
        return (im1 < im2 ? -1 : 1);
    if (id1 != id2)
        return (id1 < id2 ? -1 : 1);

    time1 = 60.0 * (60.0 * (double) ih1 + (double) imn1) + s1;
    time2 = 60.0 * (60.0 * (double) ih2 + (double) imn2) + s2;
    if (fabs(time1 - time2) > tolerance)
        return (time1 < time2 ? -1 : 1);

    return (0);

}


