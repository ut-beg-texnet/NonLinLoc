/*
 * Copyright (C) 2019 Anthony Lomax <anthony@alomax.net>
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


/*   grid2scat.c

        Program to convert a grid file to binary scatter file

 */

/*
        history:

        ver 01    09May2019  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

 */



#include "GridLib.h"


// defines


// globals


// functions

int GenScatter(int, char **);



/*** Program to  convert a grid file to a binary scatter file */

#define PNAME  "grid2scat"

int main(int argc, char *argv[]) {

    int narg;


    // set program name

    strcpy(prog_name, PNAME);


    // check command line for correct usage

    fprintf(stdout, "%s Arguments: ", prog_name);
    for (narg = 0; narg < argc; narg++)
        fprintf(stdout, "<%s> ", argv[narg]);
    fprintf(stdout, "\n");

    if (argc < 3) {
        disp_usage(PNAME, "<input grid> <num_scatter_points>");
        exit(-1);
    }

    GenScatter(argc, argv);

    exit(0);

}

int GenScatter(int argc, char *argv[]) {

    int istat;

    char fn_grid_in[FILENAME_MAX];
    FILE *fp_grid_in;
    FILE *fp_grid_in_hdr;

    GridDesc grid_in;
    SourceDesc srce;
    HypoDesc hypo;
    ScatterParams scatter;


    // set NLL constants
    SetConstants();
    // initialize random number generator
    SRAND_FUNC(RandomNumSeed);


    // open input grid file

    strcpy(fn_grid_in, argv[1]);
    if ((istat = OpenGrid3dFile(fn_grid_in, &fp_grid_in, &fp_grid_in_hdr,
            &grid_in, "", &srce, 0)) < 0) {
        nll_puterr("ERROR opening input grid file.");
        return (-1);
    }
    // allocate grid
    grid_in.buffer = AllocateGrid(&grid_in);
    if (grid_in.buffer == NULL) {
        nll_puterr(
                "ERROR: allocating memory for input grid buffer.\n");
        return (-1);
    }
    // create grid array access pointers
    grid_in.array = CreateGridArray(&grid_in);
    if (grid_in.array == NULL) {
        nll_puterr(
                "ERROR: creating array for accessing input grid buffer.\n");
        return (-1);
    }

    if ((istat = ReadGrid3dBuf(&grid_in, fp_grid_in)) < 0) {
        nll_puterr("ERROR: reading PDF grid from disk.");
        exit(EXIT_ERROR_IO);
    }


    // set number of scatter points
    sscanf(argv[2], "%d", &scatter.npts);


    // dump grid header and info
    printf("Grid Header:\n");
    printf("%d %d %d  %lf %lf %lf  %lf %lf %lf %s\n",
            grid_in.numx, grid_in.numy, grid_in.numz,
            grid_in.origx, grid_in.origy, grid_in.origz,
            grid_in.dx, grid_in.dy, grid_in.dz, grid_in.chr_type);

    if (grid_in.type == GRID_TIME || grid_in.type == GRID_TIME_2D
            || grid_in.type == GRID_ANGLE || grid_in.type == GRID_ANGLE_2D)
        printf("%s %lf %lf %lf\n", srce.label, srce.x, srce.y, srce.z);

    printf("Map projection: %s", grid_in.mapProjStr);

    char fnout[FILENAME_MAX];
    sprintf(fnout, "%s.loc", fn_grid_in);
    hypo.probmax = 1.0; //  TODO: not used?
    message_flag = 999;
    printf("Generating event scatter: num_points=%d\n", scatter.npts);
    int flag_normalize = 1;
    IntegrateGrid(&grid_in, flag_normalize);
    if ((istat = GenEventScatterGrid(&grid_in, &hypo, &scatter, fnout)) < 0) {
        nll_puterr("ERROR: generating event scatter.");
    }

    fclose(fp_grid_in);
    fclose(fp_grid_in_hdr);

    return (0);

}



