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

int apply_oct2grid(int , char ** );



/*** Program to process (add, ) 3D grid files */

#define PNAME  "oct2grid"


int main(int argc, char *argv[])
{

	int narg;


	// set program name

	strcpy(prog_name, PNAME);


	// check command line for correct usage

	fprintf(stdout, "\n%s Arguments: ", prog_name);
	for (narg = 0; narg < argc; narg++)
		fprintf(stdout, "<%s> ", argv[narg]);
	fprintf(stdout, "\n");

	if (argc < 4) {
		disp_usage(PNAME, "<input octtree> <output grid> <dx> [<dy> <dz>]");
		exit(-1);
	}

	apply_oct2grid(argc, argv);

	exit(0);

}


#define N_STEPS_CONF	11

int apply_oct2grid(int argc, char *argv[])
{
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
	if ((fp_oct_in = fopen(fn_oct_in, "r")) == NULL)
	{
		nll_puterr("ERROR opening input oct tree file.");
		return(-1);
	}

	// cread output grid
	sscanf(argv[3], "%lf", &dx);
	if (argc > 4) {
		sscanf(argv[4], "%lf", &dy);
		sscanf(argv[5], "%lf", &dz);
	} else {
		dz = dy = dx;
	}
	ptree = readTree3D(fp_oct_in);
	fclose(fp_oct_in);
	//strcpy(grid_type, "LIKELIHOOD");
	//strcpy(grid_type, "PROB_DENSITY");
	//strcpy(grid_type, "MISFIT");
	ConvertOctTree2Grid(ptree, dx, dy, dz, NULL, &grid_out);

	// output file name
	strcpy(fn_grid_out, argv[2]);
	// save grid to disk
	if ((istat = WriteGrid3dBuf(&grid_out, NULL, fn_grid_out, "octree")) < 0) {
		nll_puterr("ERROR: writing oct tree grid to disk.\n");
		return(-1);
	}


	// open confidence interval file
	/*
	sprintf(fn_conf_out, "%s.octree.conf", fn_grid_out);
	if ((fpio = fopen(fn_conf_out, "w")) == NULL) {
		nll_puterr("ERROR: opening confidence interval output file.");
		return(-1);
	}
	// write confidence levels to file
	conf_incr = 1.0 / (N_STEPS_CONF - 1);
	conf_level = 1.0;
	for (ilevel = 0; ilevel < N_STEPS_CONF; ilevel++) {
		fprintf(fpio, "%lf C %.2lf\n", conf_level, conf_level);
		conf_level -= conf_incr;
	}
	fclose(fpio);
	*/


	return(0);

}



