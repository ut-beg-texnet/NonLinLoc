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


/*   GridXY2XYZ.c

	Program to ocnvert an GRID3D xyz grid with 1 z layer to a GRID3D xyz grid
	iwth 2 z layers valif for time3D
	(copy bottom z sheet of a grid to a new z sheet)
	

*/

//------------------------------------------------------------/
// Anthony Lomax           | email: anthony@alomax.net        /
// UMR Geosciences Azur    | web: www.alomax.net              /
// 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        /
// 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        /
//------------------------------------------------------------/


/*
	history:

	ver 01    25Aug2004  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

*/



#include "../src/GridLib.h"


// defines


// globals


// functions

int AddZ(int , char ** );



/*** Program to process (add, ) 3D grid files */

#define PNAME  "GridXY2XYZ"


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


	if (argc < 5) {
		disp_usage(PNAME, "<input grid> <output grid> <minval> <maxval>");
		exit(-1);
	}

	AddZ(argc, argv);

	exit(0);

}



int AddZ(int argc, char *argv[])
{
	int istat, narg;

	char fn_grid_in[FILENAME_MAX];
	char fn_grid_out[FILENAME_MAX];
	FILE *fp_grid_in;
	FILE *fp_grid_in_hdr;

	GridDesc grid_in, grid_out;
	
	int ix, iy, iz;
	float val, val_min, val_max;

	int iSwapBytes = 0;

	

	// open input grid file

	strcpy(fn_grid_in, argv[1]);
	iSwapBytes = 0;
	if ((istat = OpenGrid3dFile(fn_grid_in, &fp_grid_in, &fp_grid_in_hdr,
			&grid_in, "", NULL, iSwapBytes)) < 0)
	{
		nll_puterr("ERROR opening input grid file.");
		return(-1);
	}

	// cread and initialize output grid

	// output file name
	strcpy(fn_grid_out, argv[2]);

	// min. max grid values
	sscanf(argv[3], "%f", &val_min);
	// min. max grid values
	sscanf(argv[4], "%f", &val_max);

	// create output grid description
	grid_out = grid_in;
	grid_out.numz++;

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

	for (ix = 0; ix < grid_in.numx; ix++) {
		printf("ix = %d/%d\r", ix, grid_in.numx);
		for (iy = 0; iy < grid_in.numy; iy++) {
			for (iz = 0; iz < grid_in.numz; iz++) {
				val = ReadGrid3dValue(fp_grid_in, ix, iy, iz, &grid_in);
				if (val < val_min)
					val = val_min;
				if (val > val_max)
					val = val_max;
				grid_out.array[ix][iy][iz] = val;
				if (iz == grid_in.numz - 1)
					grid_out.array[ix][iy][iz + 1] = val;
			}
		}
	}
	printf("\n");


	// save sum grid to disk

	if ((istat = WriteGrid3dBuf(&grid_out, NULL, fn_grid_out, "zadd")) < 0) {
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

