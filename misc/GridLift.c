/*
 * Copyright (C) 2002 Anthony Lomax <anthony@alomax.net>
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


/*   GridLift.c

	Program to lift cells in a 3D grid file
	
	Replace value in cell X with value in cell below cell X if original cell X value is > cutoff.

*/

//------------------------------------------------------------/
// Anthony Lomax           | email: anthony@alomax.net        /
// UMR Geosciences Azur    | web: www.alomax.net              /
// 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        /
// 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        /
//------------------------------------------------------------/


/*
	history:

	ver 01    18Mar2002  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

*/



#include "../src/GridLib.h"


// defines


// globals


// functions

int DoLiftProcess(int , char ** , double );



/*** Program to process (add, ) 3D grid files */

#define PNAME  "GridLift"


main(int argc, char *argv[])
{

	int istat, narg, nchar;
	double cutoff;


	// set program name

	strcpy(prog_name, PNAME);


	// check command line for correct usage

	fprintf(stdout, "\n%s Arguments: ", prog_name);
	for (narg = 0; narg < argc; narg++)
		fprintf(stdout, "<%s> ", argv[narg]);
	fprintf(stdout, "\n");

	if (argc > 1) {
		sscanf(argv[1], "%lf", &cutoff);
	} else {
		sprintf(MsgStr,
			"ERROR bad cutoff value - <%s>.", argv[1]);
		nll_puterr(MsgStr);
		disp_usage(PNAME, "<cutoff> <input grid> <output grid>");
		exit(-1);
	}

	DoLiftProcess(argc, argv, cutoff);

	exit(0);

}



int DoLiftProcess(int argc, char *argv[], double cutoff)
{
	int istat, narg;

	char fn_grid_in[FILENAME_MAX];
	char fn_grid_out[FILENAME_MAX];
	FILE *fp_grid_in;
	FILE *fp_grid_in_hdr;

	GridDesc grid_in, grid_out;
	
	int ix, iy, iz;
	float val, val_below;
	int nlift = 0;

	

	// open input grid file

	strcpy(fn_grid_in, argv[2]);
	if ((istat = OpenGrid3dFile(fn_grid_in, &fp_grid_in, &fp_grid_in_hdr,
			&grid_in, "", NULL, 0)) < 0)
	{
		nll_puterr("ERROR opening input grid file.");
		return(-1);
	}

	// cread and initialize output grid

	// output file name
	strcpy(fn_grid_out, argv[3]);

	// create output grid description
	grid_out = grid_in;
grid_in.iSwapBytes = 1;
grid_out.iSwapBytes = 0;

	// allocate grid
	grid_out.buffer = AllocateGrid(&grid_out);
	if (grid_out.buffer == NULL) {
		nll_puterr("ERROR: allocating memory for output grid buffer.\n");
		return(-1);
	}

	// create grid array access pointers
	grid_out.array = CreateGridArray(&grid_out);
	if (grid_out.array == NULL) {
		nll_puterr("ERROR: creating array for accessing output grid buffer.\n");
		return(-1);
	}

	// lift input grid file into output grid

	nlift = 0;
	for (ix = 0; ix < grid_in.numx; ix++) {
		printf("ix = %d/%d\r", ix, grid_in.numx);
		for (iy = 0; iy < grid_in.numy; iy++) {
			for (iz = 0; iz < grid_in.numz; iz++) {
				val = ReadGrid3dValue(fp_grid_in, ix, iy, iz, &grid_in);
				if (val >= cutoff && iz < grid_in.numz - 1
					&& (val_below =
						ReadGrid3dValue(fp_grid_in, ix, iy, iz + 1, &grid_in))
						< cutoff) {
					grid_out.array[ix][iy][iz] = val_below;
					nlift++;
				} else {
					grid_out.array[ix][iy][iz] = val;
				}
			}
		}
	}
	printf("\n");
	printf("%d cells lifted, Nx*Ny=%d, cutoff=%lf\n", nlift, grid_in.numx *grid_in.numy, cutoff);


	// save sum grid to disk

	if ((istat = WriteGrid3dBuf(&grid_out, NULL, fn_grid_out, "lift")) < 0) {
		nll_puterr("ERROR: writing lifted grid to disk.\n");
		return(-1);
	}


	close(fp_grid_in);
	close(fp_grid_in_hdr);

	return(0);

}



//------------------------------------------------------------/
// Anthony Lomax           | email: anthony@alomax.net        /
// UMR Geosciences Azur    | web: www.alomax.net              /
// 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        /
// 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        /
//------------------------------------------------------------/

