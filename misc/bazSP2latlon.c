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


/*   bazSP2latlon.c

        Program to read a 2D grid at specified distance and depth.


 */


/*
        history:	(see also http://alomax.net/nlloc -> Updates)

        ver 01    28FEB2018  AJL  Original version

        see GridLib.c


.........1.........2.........3.........4.........5.........6.........7.........8

 */





#define PNAME  "bazSP2latlon"

#include "GridLib.h"




/** program to demonstrate getting travel times from NLL travel time grids through a function call
 *
 *  See Makefile->bazSP2latlon for compile/link requirements
 *
 *  This demonstration program reads all input from a disk file.
 *
 *  The required travel-time grids must be present as *.buf and *.hyp disk files.
 *
 */


#define NARGS_MIN 8
#define ARG_DESC "<input P time grid 2D root> <input S time grid 2D root> <back azimuth (deg)> <S-P time (s)> <depth( km)> <lat0> <lon0>"

#define TOLERANCE 0.01  // sec

int main(int argc, char *argv[]) {

    int istat;

    // set program name
    strcpy(prog_name, PNAME);

    // check command line for correct usage
    if (argc < NARGS_MIN) {
        disp_usage(prog_name, ARG_DESC);
        printf("Output:\n");
        printf("baz sp_time lat lon depth\n");
        printf("Example:\n");
        printf("bazSP2latlon time/testModel-r4-s3-t3LV_ALL.P.TX31.time time/testModel-r4-s3-t3LV_ALL.S.TX31.time 0 20 5 29.334200 -103.667800\n");
        printf("0.000000 20.000000 30.719041 -103.667800 5.000000\n");

        return (EXIT_ERROR_USAGE);
    }

    // grid file name
    char fn_P_grid_root[MAXLINE];
    strcpy(fn_P_grid_root, argv[1]);
    FILE *fp_grid_P, *fp_hdr_P;
    GridDesc grid_P;
    if ((istat = OpenGrid3dFile(fn_P_grid_root, &fp_grid_P, &fp_hdr_P, &grid_P, "", NULL, 0)) < 0) {
        nll_puterr2("ERROR opening grid file:", fn_P_grid_root);
        return (-1);
    }

    char fn_S_grid_root[MAXLINE];
    strcpy(fn_S_grid_root, argv[2]);
    FILE *fp_grid_S, *fp_hdr_S;
    GridDesc grid_S;
    if ((istat = OpenGrid3dFile(fn_S_grid_root, &fp_grid_S, &fp_hdr_S, &grid_S, "", NULL, 0)) < 0) {
        nll_puterr2("ERROR opening grid file:", fn_S_grid_root);
        return (-1);
    }

    // distance, depth
    double baz, sp_time, depth, lat0, lon0;
    sscanf(argv[3], "%lf", &baz);
    sscanf(argv[4], "%lf", &sp_time);
    sscanf(argv[5], "%lf", &depth);
    sscanf(argv[6], "%lf", &lat0);
    sscanf(argv[7], "%lf", &lon0);

    //printf("DEBUG: fn_P_grid_root=%s fn_S_grid_root=%s baz=%f s-p=%f depth=%f lat0=%f lon0=%f\n", fn_P_grid_root, fn_S_grid_root, baz, sp_time, depth, lat0, lon0);

    // convert baz to forward azimuth
    double azimuth = baz + 180.0;
    if (azimuth >= 360.0) {
        azimuth -= 360.0;
    }
    SetConstants();

    double ymin = grid_P.origy; // should be 0.0 for grid 2D
    double ymax = ymin + grid_P.dy * (double) (grid_P.numy - 1);
    double yloc = (ymin + ymax) / 2.0;
    double yloc_max = ymax;
    double yloc_min = ymin;
    int niter = 0;
    while (niter++ < 1000) {
        double ptime = ReadAbsInterpGrid2d(fp_grid_P, &grid_P, yloc, depth);
        double stime = ReadAbsInterpGrid2d(fp_grid_S, &grid_S, yloc, depth);
        double sp_time_calc = stime - ptime;
        double diff = sp_time - sp_time_calc;

        //printf("n %d  ymin %f  ymax %f  y %f  calc %f  diff %f\n", niter, yloc_min, yloc_max, yloc, sp_time_calc, diff);

        if (fabs(diff) < TOLERANCE) {
            break;
        }
        // binary search
        if (diff > 0.0) {
            yloc_min = yloc;
            yloc = (yloc + yloc_max) / 2.0;
        } else {
            yloc_max = yloc;
            yloc = (yloc_min + yloc) / 2.0;
        }
    }

    double lat, lon;
    PointAtGCDistanceAzimuth(lat0, lon0, yloc * KM2DEG, azimuth, &lat, &lon);

    printf("%f %f %f %f %f\n", baz, sp_time, lat, lon, depth);

    return (0);

}




