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


/*   read_grid2d.c

        Program to read a 2D grid at specified distance and depth.


 */


/*
        history:	(see also http://alomax.net/nlloc -> Updates)

        ver 01    28FEB2018  AJL  Original version

        see GridLib.c


.........1.........2.........3.........4.........5.........6.........7.........8

 */





#define PNAME  "read_grid2d"

#include "GridLib.h"




/** program to demonstrate getting travel times from NLL travel time grids through a function call
 *
 *  See Makefile->read_grid2d for compile/link requirements
 *
 *  This demonstration program reads all input from a disk file.
 *
 *  The required travel-time grids must be present as *.buf and *.hyp disk files.
 *
 */


#define NARGS_MIN 5
#define ARG_DESC "<input grid root> <type time|angle> <distance(km)> <depth<km>"

int main(int argc, char *argv[]) {

    // set program name
    strcpy(prog_name, PNAME);

    // check command line for correct usage
    if (argc < NARGS_MIN) {
        disp_usage(prog_name, ARG_DESC);
        return (EXIT_ERROR_USAGE);
    }

    // grid file name
    char fn_grid_root[MAXLINE];
    strcpy(fn_grid_root, argv[1]);
    char file_type[MAXLINE];
    strcpy(file_type, argv[2]);

    // distance, depth
    double xloc[1], yloc[1], zloc[1];
    xloc[0] = 0.0;
    sscanf(argv[3], "%lf", yloc);
    sscanf(argv[4], "%lf", zloc);
    //printf("DEBUG: fn_grid_root=%s x= %f dist=%f depth=%f\n", fn_grid_root, xloc[0], yloc[0], zloc[0]);

    SetConstants();


    /** ===========================================================================
     *  call function that gets values from grid (see GridLib.c ReadGridFile())
     */
    GRID_FLOAT_TYPE grid_value[1];
    int iSwapBytes = 0;
    ReadGridFile(grid_value, fn_grid_root, file_type, xloc, yloc, zloc, 1, iSwapBytes);

    // write grid value to stdout
    //printf("dist=%f depth=%f value=%f\n", yloc[0], zloc[0], grid_value[0]);
    printf("%f\n", grid_value[0]);


return (0);

}




