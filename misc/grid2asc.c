/*
 * Copyright (C) 2006 Anthony Lomax <anthony@alomax.net>
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


/*   grid2asc.c

	Program to dump contents of grid file to ascii

*/

/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
161 Allee du Micocoulier, 06370 Mouans-Sartoux, France
tel: +33(0)493752502  e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/

/*
	history:

	ver 01    18Mar2002  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

*/



#include "../src/GridLib.h"


// defines


// globals


// functions

int DumpToASCII(int , char ** );



/*** Program to  convert grid files to ascii */

#define PNAME  "grid2asc"


int main(int argc, char *argv[])
{

	int istat, narg, nchar;


	// set program name

	strcpy(prog_name, PNAME);


	// check command line for correct usage

	fprintf(stdout, "%s Arguments: ", prog_name);
	for (narg = 0; narg < argc; narg++)
		fprintf(stdout, "<%s> ", argv[narg]);
	fprintf(stdout, "\n");

	if (argc < 2) {
		disp_usage(PNAME, "<input grid>");
		exit(-1);
	}

	DumpToASCII(argc, argv);

	exit(0);

}



int DumpToASCII(int argc, char *argv[])
{
	int istat, narg;

	char fn_grid_in[FILENAME_MAX];
	FILE *fp_grid_in;
	FILE *fp_grid_in_hdr;

	GridDesc grid_in;
	SourceDesc srce;

	int ix, iy, iz;
	float val;



	// open input grid file

	strcpy(fn_grid_in, argv[1]);
	if ((istat = OpenGrid3dFile(fn_grid_in, &fp_grid_in, &fp_grid_in_hdr,
	     &grid_in, "", &srce, 0)) < 0)
	{
		nll_puterr("ERROR opening input grid file.");
		return(-1);
	}


	// dump grid header and info

	printf("Dump of NLL grid file: %s\n", fn_grid_in);
	printf("Order: (((vector in Z) for each row in Y) for each column in X)\n");
	printf("i.e., v(x,y,z):  v(0,0,0),v(0,0,1),...,v(numx-1,numy-1,numz-2),v(numx-1,numy-1,numz-1)\n");
	printf("Each vector in Z is terminated by end-of-line <cr><lf>\n");
	printf("Grid Header:\n");
	printf("%d %d %d  %lf %lf %lf  %lf %lf %lf %s\n",
		grid_in.numx, grid_in.numy, grid_in.numz,
		grid_in.origx, grid_in.origy, grid_in.origz,
		grid_in.dx, grid_in.dy, grid_in.dz, grid_in.chr_type);

	if (grid_in.type == GRID_TIME || grid_in.type == GRID_TIME_2D
		   || grid_in.type == GRID_ANGLE || grid_in.type == GRID_ANGLE_2D)
		printf("%s %lf %lf %lf\n",
			srce.label, srce.x, srce.y, srce.z);

	printf("%s\n", MapProjStr[0]);


	// dump input grid file to ASCII

	for (ix = 0; ix < grid_in.numx; ix++) {
		for (iy = 0; iy < grid_in.numy; iy++) {
			for (iz = 0; iz < grid_in.numz; iz++) {
				val = ReadGrid3dValue(fp_grid_in, ix, iy, iz, &grid_in);
				printf("%f ", val);
			}
			printf("\n");
		}
	}
	printf("\n");


	fclose(fp_grid_in);
	fclose(fp_grid_in_hdr);

	return(0);

}



