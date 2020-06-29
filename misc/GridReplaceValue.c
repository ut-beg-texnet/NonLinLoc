/*
 * Copyright (C) 2007 Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
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


/*   GridReplaceValue.c

	Program to replace a value in a 3D grid file

	Replace value in cell X with new value if original X value is within a tolerance of a specified old value.

*/



/*
	history:

	ver 01    15May2007  AJL  Original version from GridReplaceValue.c


.........1.........2.........3.........4.........5.........6.........7.........8

*/



#include "../src/GridLib.h"


// defines


// globals


// functions

int DoLiftProcess(int , char ** , double );



/*** Program to process (add, ) 3D grid files */

#define PNAME  "GridReplaceValue"


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

	if (argc < 6){
		disp_usage(PNAME, "<old value> <new value> <tolerance> <input grid> <output grid>");
		fprintf(stderr, "   Replaces value in cell X with new value if original X value is within a tolerance of a specified old value.\n");
		exit(-1);
	}

	DoGridReplaceValueProcess(argc, argv);

	exit(0);

}



int DoGridReplaceValueProcess(int argc, char *argv[])
{
	int istat, narg;

	char fn_grid_in[FILENAME_MAX];
	char fn_grid_out[FILENAME_MAX];
	FILE *fp_grid_in;
	FILE *fp_grid_in_hdr;

	float old_value, new_value, tolerance;
	GridDesc grid_in, grid_out;

	int ix, iy, iz;
	float val, old_value_minus_tolerance, old_value_plus_tolerance;
	int nreplaced = 0;


	// get numeric parameters
	sscanf(argv[1], "%f", &old_value);
	sscanf(argv[2], "%f", &new_value);
	sscanf(argv[3], "%f", &tolerance);
	old_value_minus_tolerance = old_value - tolerance;
	old_value_plus_tolerance = old_value + tolerance;

	// open input grid file

	strcpy(fn_grid_in, argv[4]);
	if ((istat = OpenGrid3dFile(fn_grid_in, &fp_grid_in, &fp_grid_in_hdr,
			&grid_in, "", NULL, 0)) < 0)
	{
		nll_puterr("ERROR opening input grid file.");
		return(-1);
	}

	// cread and initialize output grid

	// output file name
	strcpy(fn_grid_out, argv[5]);

	// create output grid description
	grid_out = grid_in;
grid_in.iSwapBytes = 0;
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

	// copy input grid values into output grid, replacing where necessary

	nreplaced = 0;
	for (ix = 0; ix < grid_in.numx; ix++) {
		printf("ix = %d/%d\r", ix, grid_in.numx);
		for (iy = 0; iy < grid_in.numy; iy++) {
			for (iz = 0; iz < grid_in.numz; iz++) {
				val = ReadGrid3dValue(fp_grid_in, ix, iy, iz, &grid_in);
				if (val >= old_value_minus_tolerance && val <= old_value_plus_tolerance) {
					grid_out.array[ix][iy][iz] = new_value;
					nreplaced++;
				} else {
					grid_out.array[ix][iy][iz] = val;
				}
			}
		}
	}
	printf("\n");
	printf("%d/%d cell values replaced\n", nreplaced, grid_in.numx * grid_in.numy * grid_in.numz);


	// save sum grid to disk

	if ((istat = WriteGrid3dBuf(&grid_out, NULL, fn_grid_out, "replace")) < 0) {
		nll_puterr("ERROR: writing output grid to disk.\n");
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

