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


/*   oct2grid.c

        Program to convert OctTree file to Grid file



 */

/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
161 Allee du Micocoulier, 06370 Mouans-Sartoux, France
tel: +33(0)493752502  e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


/*
        history:

        ver 01    25Nov2005  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

 */



#include "GridLib.h"


// defines


// globals


// functions

int apply_oct2grid(int, char **);
int CalcConfidenceIntrvl(GridDesc* ptgrid, char* filename);



/*** Program to convert an octtree grid to a regular grid
 *
 * output grid is optionally normalized by its integral
 *
 */

#define PNAME  "oct2grid"

int main(int argc, char *argv[]) {

    int narg;


    // set program name

    strcpy(prog_name, PNAME);


    // check command line for correct usage

    fprintf(stdout, "\n%s Arguments: ", prog_name);
    for (narg = 0; narg < argc; narg++)
        fprintf(stdout, "<%s> ", argv[narg]);
    fprintf(stdout, "\n");

    if (argc < 5) {
        disp_usage(PNAME, "<input octree> <output grid> <integral_norm_flag> <dx> [<dy> <dz>]");
        exit(-1);
    }

    message_flag = 99;

    apply_oct2grid(argc, argv);

    exit(0);

}


#define N_STEPS_CONF 11

int apply_oct2grid(int argc, char *argv[]) {
    int istat;

    char fn_oct_in[FILENAME_MAX];
    char fn_grid_out[FILENAME_MAX];
    FILE *fp_oct_in;

    GridDesc grid_out;
    Tree3D* ptree;

    double dx, dy, dz;
    //char grid_type[MAXLINE];

    //char fn_conf_out[FILENAME_MAX];
    //FILE *fpio;
    //int ilevel;
    //double conf_level, conf_incr;


    // open input grid file
    strcpy(fn_oct_in, argv[1]);
    if ((fp_oct_in = fopen(fn_oct_in, "r")) == NULL) {
        nll_puterr("ERROR opening input oct tree file.");
        return (-1);
    }

    // normalize flag
    int integral_norm_flag = 0;
    sscanf(argv[3], "%d", &integral_norm_flag);

    // create output grid
    sscanf(argv[4], "%lf", &dx);
    if (argc > 5) {
        sscanf(argv[5], "%lf", &dy);
        sscanf(argv[6], "%lf", &dz);
    } else {
        dz = dy = dx;
    }
    ptree = readTree3D(fp_oct_in);
    fclose(fp_oct_in);
    //strcpy(grid_type, "LIKELIHOOD");
    //strcpy(grid_type, "PROB_DENSITY");
    //strcpy(grid_type, "MISFIT");

    // check for global projection
    if (ptree->isSpherical) {
        strcpy(MapProjStr[0], "TRANSFORM  GLOBAL");
    }


    printf("Convert oct-tree to grid...\n");
    ConvertOctTree2Grid(ptree, dx, dy, dz, NULL, &grid_out);

    grid_out.type = GRID_PROB_DENSITY; // TODO: assume input oct-tree grid is relative prob density
    convert_grid_type(&grid_out, 0);

    if (integral_norm_flag) {
        int flag_normalize = 1;
        double integral = IntegrateGrid(&grid_out, flag_normalize);
        printf("Normalized to integral=%f\n", integral);
    }

    // output file name
    strcpy(fn_grid_out, argv[2]);
    // save grid to disk
    printf("Write grid to file %s.octree ...\n", fn_grid_out);
    if ((istat = WriteGrid3dBuf(&grid_out, NULL, fn_grid_out, "octree")) < 0) {
        nll_puterr("ERROR: writing oct tree grid to disk.\n");
        return (-1);
    }

    // calculate and write confidence intervals
    CalcConfidenceIntrvl(&grid_out, fn_grid_out);


    return (0);

}



/** function to calculate confidence intervals and save to file */
/*		(MEN92, eq. 25ff) */

#define N_STEPS_SRCH 101
#define N_STEPS_CONF 11

int CalcConfidenceIntrvl(GridDesc* ptgrid, char* filename) {
    FILE *fpio;
    char fname[FILENAME_MAX];

    int ix, iy, iz;
    int iconf, isrch;
    double srch_level, srch_incr, conf_level, conf_incr, prob_den;
    double srch_sum[N_STEPS_SRCH];
    double contour[N_STEPS_CONF];


    /* write message */
    printf("Calculating confidence intervals over grid...\n");

    // find maximum value in grid
    printf("   find maximum value in grid...");
    double value_max = -1.0;
    for (ix = 0; ix < ptgrid->numx; ix++) {
        for (iy = 0; iy < ptgrid->numy; iy++) {
            for (iz = 0; iz < ptgrid->numz; iz++) {
                prob_den = ((GRID_FLOAT_TYPE ***) ptgrid->array)[ix][iy][iz];
                if (prob_den >= value_max)
                    value_max = prob_den;
            }
        }
    }
    printf(" = %f\n", value_max);

    for (isrch = 0; isrch < N_STEPS_SRCH; isrch++)
        srch_sum[isrch] = 0.0;


    /* accumulate approx integral of probability density in search bins */

    printf("   accumulate approx integral of pdf in bins...\n");
    srch_incr = value_max / (N_STEPS_SRCH - 1);
    for (ix = 0; ix < ptgrid->numx; ix++) {
        for (iy = 0; iy < ptgrid->numy; iy++) {
            for (iz = 0; iz < ptgrid->numz; iz++) {
                prob_den = ((GRID_FLOAT_TYPE ***) ptgrid->array)[ix][iy][iz];
                srch_level = 0.0;
                for (isrch = 0; isrch < N_STEPS_SRCH; isrch++) {
                    if (prob_den >= srch_level)
                        srch_sum[isrch] += prob_den;
                    srch_level += srch_incr;
                }

            }
        }
    }
    ptgrid->sum = 1.0;

    /* normalize by 100% confidence level sum */

    for (isrch = 1; isrch < N_STEPS_SRCH; isrch++)
        srch_sum[isrch] /= srch_sum[0];
    srch_sum[0] = 1.0;


    /* open confidence interval file */

    sprintf(fname, "%s.octree.conf", filename);
    if ((fpio = fopen(fname, "w")) == NULL) {
        nll_puterr("ERROR: opening confidence interval output file.");
        return (-1);
    }

    /* find confidence levels and write to file */

    printf("   write confidence levels to file %s ...\n", fname);
    conf_incr = 1.0 / (N_STEPS_CONF - 1);
    conf_level = 1.0;
    iconf = N_STEPS_CONF - 1;
    for (isrch = 0; isrch < N_STEPS_SRCH; isrch++) {
        if (srch_sum[isrch] <= conf_level) {
            contour[iconf] = (double) isrch * srch_incr;
            fprintf(fpio, "%lf C %.2lf\n",
                    contour[iconf], conf_level);
            if (--iconf < 0)
                break;
            conf_level -= conf_incr;
        }
    }


    fclose(fpio);

    return (0);

}
