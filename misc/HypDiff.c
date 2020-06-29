/*
 * Copyright (C) 1999-2003 Anthony Lomax <anthony@alomax.net>
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


/*   HypDiff.c

        Program to find the difference in GEOGRAPHIC hypocenters in 2 NonLinLoc hyp files

 */



/*
        history:

        ver 01    23Apr2001  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

 */



#include "../src/GridLib.h"


/* defines */

#define MAX_NUM_INPUT_FILES 4096


/* globals */



/* functions */

int DiffHypocenters(int argc, char** argv);
int compareTimes(int iy1, int iy2, int im1, int im2, int id1, int id2, int ih1, int ih2,
        int imn1, int imn2, double s1, double s2, double tolerance, double *ptdiff);


#define PNAME  "HypDiff"

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
                "hypfile1 hypfile2 time_tolerance nobs_tolerance");
        exit(-1);
    }

    SetConstants();
    message_flag = 1;
    DispProgInfo();
    message_flag = 0;

    if ((istat = DiffHypocenters(argc, argv)) < 0) {
        nll_puterr("ERROR doing difference process.");
        exit(-1);
    }



    exit(0);

}

int DiffHypocenters(int argc, char** argv) {

    int istat1, istat2;
    char fn_hyp1[FILENAME_MAX];
    char fn_hyp2[FILENAME_MAX];
    char fn_out[FILENAME_MAX];
    char fn_out1[FILENAME_MAX];
    char fn_out2[FILENAME_MAX];
    FILE *fp_hyp1, *fp_hyp2;
    FILE *fp_hyp1_out, *fp_hyp2_out, *fp_gmt_out;
    FILE *fp_gmt_dist_out, *fp_gmt_az_out, *fp_gmt_depth_out, *fp_gmt_horiz_out, *fp_gmt_all_out;

    char *pchr, *pstmp1, *pstmp2;
    double time_tolerance;
    int nobs_tolerance;
    int nLocMatched, nHypo1 = 0, nHypo2 = 0;
    int idiff, ok;

    HypoDesc Hypo1, Hypo2;

    GridDesc locgrid;

    double angle, x1, x2, y1, y2, dx, dy, dz, azim, tdiff;


    /* get command line parameters */
    strcpy(fn_hyp1, argv[1]);
    strcpy(fn_hyp2, argv[2]);
    sscanf(argv[3], "%lf", &time_tolerance);
    sscanf(argv[4], "%d", &nobs_tolerance);

    pstmp1 = strrchr(fn_hyp1, '/');
    if (pstmp1 == NULL)
        pstmp1 = fn_hyp1;
    else
        pstmp1++;
    pstmp2 = strrchr(fn_hyp2, '/');
    if (pstmp2 == NULL)
        pstmp2 = fn_hyp2;
    else
        pstmp2++;

    /* open hyp input files */

    if ((fp_hyp1 = fopen(fn_hyp1, "r")) == NULL) {
        nll_puterr2("ERROR: opening hyp file: ", fn_hyp1);
        return (-1);
    }
    if ((fp_hyp2 = fopen(fn_hyp2, "r")) == NULL) {
        nll_puterr2("ERROR: opening hyp file: ", fn_hyp2);
        return (-1);
    }
    /* open hyp output files */

    sprintf(fn_out1, "diff_%s", pstmp1);
    if ((fp_hyp1_out = fopen(fn_out1, "w")) == NULL) {
        nll_puterr2("ERROR: opening hyp output file: ", fn_out1);
        return (-1);
    }
    sprintf(fn_out2, "diff_%s", pstmp2);
    if ((fp_hyp2_out = fopen(fn_out2, "w")) == NULL) {
        nll_puterr2("ERROR: opening hyp output file: ", fn_out2);
        return (-1);
    }

    /* open gmt output files */

    sprintf(fn_out, "diff_%s_%s.xyz", pstmp1, pstmp2);
    if ((fp_gmt_out = fopen(fn_out, "w")) == NULL) {
        nll_puterr2("ERROR: opening gmt output file:", fn_out);
        return (-1);
    }
    sprintf(fn_out, "diff_%s_%s.dist", pstmp1, pstmp2);
    if ((fp_gmt_dist_out = fopen(fn_out, "w")) == NULL) {
        nll_puterr2("ERROR: opening gmt dist output file:", fn_out);
        return (-1);
    }
    sprintf(fn_out, "diff_%s_%s.az", pstmp1, pstmp2);
    if ((fp_gmt_az_out = fopen(fn_out, "w")) == NULL) {
        nll_puterr2("ERROR: opening gmt az output file:", fn_out);
        return (-1);
    }
    sprintf(fn_out, "diff_%s_%s.depth", pstmp1, pstmp2);
    if ((fp_gmt_depth_out = fopen(fn_out, "w")) == NULL) {
        nll_puterr2("ERROR: opening gmt depth output file:", fn_out);
        return (-1);
    }
    sprintf(fn_out, "diff_%s_%s.horiz", pstmp1, pstmp2);
    if ((fp_gmt_horiz_out = fopen(fn_out, "w")) == NULL) {
        nll_puterr2("ERROR: opening gmt horiz output file:", fn_out);
        return (-1);
    }

    sprintf(fn_out, "diff_%s_%s.all", pstmp1, pstmp2);
    if ((fp_gmt_all_out = fopen(fn_out, "w")) == NULL) {
        nll_puterr2("ERROR: opening gmt horiz output file:", fn_out);
        return (-1);
    }
    fprintf(fp_gmt_all_out,
            "year1 month1 day1  hour1 min1 sec1  lat1 long1 depth1 nphs1   "
            "year2 month2 day2  hour2 min2 sec2  lat2 long2 depth2 nphs2   "
            "dt ds dh dz daz"
            "\n");


    /* difference corresponding hypocenters */

    nLocMatched = 0;
    if ((istat2 = GetHypLoc(fp_hyp2, fn_hyp2, &Hypo2, NULL, &NumArrivals, 0, &locgrid, 0))
            >= 0 && istat2 != EOF)
        nHypo2++;

    while (istat2 >= 0 && istat2 != EOF &&
            (istat1 = GetHypLoc(fp_hyp1, fn_hyp1, &Hypo1, NULL, &NumArrivals, 0, &locgrid, 0))
            >= 0 && istat1 != EOF) {

        nHypo1++;

        ok = 0;
        while (1) {
            idiff = compareTimes(Hypo1.year, Hypo2.year, Hypo1.month, Hypo2.month,
                    Hypo1.day, Hypo2.day, Hypo1.hour, Hypo2.hour,
                    Hypo1.min, Hypo2.min, Hypo1.sec, Hypo2.sec,
                    time_tolerance, &tdiff);

            //WriteLocation(stdout, &Hypo1, Arrival, 0, fn_hyp1, 0, 0, 1, &locgrid, 0);
            //WriteLocation(stdout, &Hypo2, Arrival, 0, fn_hyp2, 0, 0, 1, &locgrid, 0);

            if (idiff == 0) { // match
                ok = 1;
                break;
            } else if (idiff > 0) { // hypo2 earlier, read next hypo2
                if ((istat2 = GetHypLoc(fp_hyp2, fn_hyp2, &Hypo2, NULL,
                        &NumArrivals, 0, &locgrid, 0))
                        >= 0 && istat2 != EOF) {
                    nHypo2++;
                    continue;
                } else // error or EOF hypo2
                    break;
            } else { // hypo2 later, read next hypo1
                break;
            }
        }
        if (!ok)
            continue;

        if (fabs((double) (Hypo1.nreadings - Hypo2.nreadings)) > (double) (nobs_tolerance) + 0.1) {
            //fprintf(stdout, "nobs diff:\n");
            //WriteLocation(stdout, &Hypo1, Arrival, 0, fn_hyp1, 0, 0, 1, &locgrid, 0);
            //WriteLocation(stdout, &Hypo2, Arrival, 0, fn_hyp2, 0, 0, 1, &locgrid, 0);
            continue;
        }
        // output difference

        WriteLocation(fp_hyp1_out, &Hypo1, Arrival, 0, fn_hyp1, 0, 0, 0, &locgrid, 0);
        fprintf(fp_hyp1_out, "END_NLLOC\n\n");
        WriteLocation(fp_hyp2_out, &Hypo2, Arrival, 0, fn_hyp2, 0, 0, 0, &locgrid, 0);
        fprintf(fp_hyp2_out, "END_NLLOC\n\n");

        strcpy(fn_out, fn_out1);
        pchr = strstr(fn_out, ".hyp");
        *pchr = '\0';
        WriteGrid3dHdr(&locgrid, NULL, fn_out, NULL);

        strcpy(fn_out, fn_out2);
        pchr = strstr(fn_out, ".hyp");
        *pchr = '\0';
        WriteGrid3dHdr(&locgrid, NULL, fn_out, NULL);

        fprintf(fp_gmt_out, "> GMT_LATLONDEPTH\n");
        fprintf(fp_gmt_out, "%lf %lf %lf\n", Hypo1.dlat, Hypo1.dlong, Hypo1.depth);
        fprintf(fp_gmt_out, "%lf %lf %lf\n", Hypo2.dlat, Hypo2.dlong, Hypo2.depth);

        // set simple transform (assume no rotation!)
        map_itype[0] = MAP_TRANS_SIMPLE;
        strcpy(map_trans_type[0], "SIMPLE");
        map_orig_lat[0] = Hypo1.dlat;
        map_orig_long[0] = Hypo2.dlong;
        map_rot[0] = 0.0;
        angle = -cRPD * map_rot[0];
        map_cosang[0] = cos(angle);
        map_sinang[0] = sin(angle);

        latlon2rect(0, Hypo1.dlat, Hypo1.dlong, &x1, &y1);
        latlon2rect(0, Hypo2.dlat, Hypo2.dlong, &x2, &y2);
        dx = x1 - x2;
        dy = y1 - y2;
        dz = Hypo1.depth - Hypo2.depth;
        fprintf(fp_gmt_dist_out, "%lf\n", sqrt(dx * dx + dy * dy + dz * dz));
        fprintf(fp_gmt_depth_out, "%lf\n", dz);
        fprintf(fp_gmt_horiz_out, "%lf\n", sqrt(dx * dx + dy * dy));
        azim = atan2(dx, dy) / cRPD;
        if (azim < 0.0)
            azim += 360.0;
        fprintf(fp_gmt_az_out, "%lf\n", azim);

        // all parameters and measures
        fprintf(fp_gmt_all_out, "%4.4d %2.2d %2.2d  %2.2d %2.2d %lf  %lf %lf %lf %d    ",
                Hypo1.year, Hypo2.month, Hypo1.day,
                Hypo1.hour, Hypo1.min, (double) Hypo1.sec,
                Hypo1.dlat, Hypo1.dlong, Hypo1.depth, Hypo1.nreadings
                );
        fprintf(fp_gmt_all_out, "%4.4d %2.2d %2.2d  %2.2d %2.2d %lf  %lf %lf %lf %d    ",
                Hypo2.year, Hypo2.month, Hypo2.day,
                Hypo2.hour, Hypo2.min, (double) Hypo2.sec,
                Hypo2.dlat, Hypo2.dlong, Hypo2.depth, Hypo2.nreadings
                );
        fprintf(fp_gmt_all_out, "%lf %lf %lf %lf %lf ",
                tdiff,
                sqrt(dx * dx + dy * dy + dz * dz),
                sqrt(dx * dx + dy * dy),
                dz,
                azim
                );
        fprintf(fp_gmt_all_out, "\n");


        nLocMatched++;

        if ((istat2 = GetHypLoc(fp_hyp2, fn_hyp2, &Hypo2, NULL, &NumArrivals, 0, &locgrid, 0))
                >= 0 && istat2 != EOF)
            nHypo2++;

    }

    printf("\n%d/%d read, %d matched.\n:", nHypo1, nHypo2, nLocMatched);

    return (0);

}

int compareTimes(int iy1, int iy2, int im1, int im2, int id1, int id2, int ih1, int ih2,
        int imn1, int imn2, double s1, double s2, double tolerance, double *ptdiff) {

    double time1, time2;

    *ptdiff = 0.0;

    if (iy1 != iy2)
        return (iy1 < iy2 ? -1 : 1);
    if (im1 != im2)
        return (im1 < im2 ? -1 : 1);
    if (id1 != id2)
        return (id1 < id2 ? -1 : 1);

    time1 = 60.0 * (60.0 * (double) ih1 + (double) imn1) + s1;
    time2 = 60.0 * (60.0 * (double) ih2 + (double) imn2) + s2;
    *ptdiff = time1 - time2;
    if (fabs(*ptdiff) > tolerance)
        return (time1 < time2 ? -1 : 1);

    return (0);

}


