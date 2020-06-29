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


/*   GridRotate.c

	Program to rotate a 3D grid file
	1) 90CW =  90 deg clockwise (map Y to X, X to -Y

*/

//------------------------------------------------------------/
// Anthony Lomax           | email: lomax@faille.unice.fr     /
// UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax /
// 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        /
// 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        /
//------------------------------------------------------------/


/*
	history:

	ver 01    14Mar2002  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

*/



#include "../src/GridLib.h"


// defines

#define MAX_NUM_INPUT_FILES 500


// globals


// functions

int DoRotateProcess(int , char **);



/*** Program to process (add, ) 3D grid files */

#define PNAME  "GridRotate"


main(int argc, char *argv[])
{

	int istat, narg, nchar;
	char process[MAXLINE];


	// set program name

	strcpy(prog_name, PNAME);


	// check command line for correct usage

	fprintf(stdout, "\n%s Arguments: ", prog_name);
	for (narg = 0; narg < argc; narg++)
		fprintf(stdout, "<%s> ", argv[narg]);
	fprintf(stdout, "\n");

	if (argc > 1)
		sscanf(argv[1], "%s", process);
	else
		strcpy(process, "");
	nchar = 0;
	while ( nchar < MAXLINE && (process[nchar] = toupper(process[nchar])) )
		nchar++;

	if (strcmp(process, "90CW") == 0) {
		if (argc < 4) {
			nll_puterr("ERROR wrong number of command line arguments.");
			disp_usage(PNAME, "90CW <input_gridfile> <output_gridfile>");
			exit(-1);
		}
		if ((istat = Rotate90CW(argc, argv)) < 0) {
			nll_puterr("ERROR doing Rotate 90CW process.");
			exit(-1);
		}
	} else {
		sprintf(MsgStr,
			"ERROR unrcognized process - <%s>.", process);
		nll_puterr(MsgStr);
		disp_usage(PNAME, "90CW/... <arguments>");
		exit(-1);
	}


	exit(0);

}



int Rotate90CW(int argc, char *argv[])
{
	int istat, narg;

	char fn_grid_in[FILENAME_MAX];
	char fn_grid_out[FILENAME_MAX];
	FILE *fp_grid_in;
	FILE *fp_grid_in_hdr;

	GridDesc grid_in, grid_out;
	
	int ix, iy, iz;
	float val;

	

	// open input grid file

	strcpy(fn_grid_in, argv[2]);
	if ((istat = OpenGrid3dFile(fn_grid_in, &fp_grid_in, &fp_grid_in_hdr,
			&grid_in, "", NULL)) < 0)
	{
		nll_puterr("ERROR opening input grid file.");
		return(-1);
	}

	// cread and initialize output grid

	// output file name
	strcpy(fn_grid_out, argv[3]);

	// create output grid description
	grid_out = grid_in;
	grid_out.numx = grid_in.numy;
	grid_out.numy = grid_in.numx;

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

	// rotate input grid file into output grid

	for (ix = 0; ix < grid_in.numx; ix++) {
		printf("ix = %d/%d\r", ix, grid_in.numx);
		for (iy = 0; iy < grid_in.numy; iy++) {
			for (iz = 0; iz < grid_in.numz; iz++) {
				val = ReadGrid3dValue(fp_grid_in, ix, iy, iz, &grid_in);
				grid_out.array[iy][grid_out.numy - ix - 1][iz] = val;
			}
		}
	}
	printf("\n");


	// save sum grid to disk

	if ((istat = WriteGrid3dBuf(&grid_out, NULL, fn_grid_out, "rot90CW")) < 0) {
		nll_puterr("ERROR: writing rotated grid to disk.\n");
		return(-1);
	}


	close(fp_grid_in);
	close(fp_grid_in_hdr);

	return(0);

}



//------------------------------------------------------------/
// Anthony Lomax           | email: lomax@faille.unice.fr     /
// UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax /
// 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        /
// 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        /
//------------------------------------------------------------/

