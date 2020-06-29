/* 
 * Copyright (C) 1999 Anthony Lomax <lomax@faille.unice.fr>
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


/*   Grid2Stat.c

	Program to calculate "traditional" measures (mean, covariance, ...) 
	from 3D grid PROB_DENSITY files 

	output to stdout

*/

//------------------------------------------------------------/
// Anthony Lomax           | email: lomax@faille.unice.fr     /
// UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax /
// 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        /
// 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        /
//------------------------------------------------------------/


/*
	history:

	ver 01    08JAN1998  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

*/



#include "../src/GridLib.h"
#include "GridGraphLib.h"


// defines



// globals

char fninput[FILENAME_MAX];
#define MAX_NUM_INPUT_FILES 500
char fnroot_input_list[MAX_NUM_INPUT_FILES][FILENAME_MAX];
char fnroot_input[FILENAME_MAX];


// functions

int GenStats(char* );
void GenErrorBars(FILE* , Vect3D* , Mtrx3D* );
void GenCross(FILE* , Vect3D* );


/** program to generate "traditional" statistical measures from PDF grid files 
	and write output to stdout */

#define PNAME  "Grid2Stat"


main(int argc, char *argv[])
{

	int istat, narg;
	char test_str[] = ".buf";
	char *chrpos;
	int NumFiles, nFile;


	// set program name

	strcpy(prog_name, PNAME);


	// check command line for correct usage

	fprintf(stdout, "\n%s Arguments: ", prog_name);
	for (narg = 0; narg < argc; narg++)
		fprintf(stdout, "<%s> ", argv[narg]);
	fprintf(stdout, "\n");

	if (argc != 2)
	{
		nll_puterr("ERROR wrong number of command line arguments.");
		disp_usage(PNAME, "<grdfile>");
		exit(-1);
	}


	// read command line arguments

	strcpy(fninput, argv[1]);



	// set constants

	SetConstants();
	message_flag = 0;


	// check for wildcards in input file name
	strcat(fninput, test_str);
	NumFiles = ExpandWildCards(fninput, 
			fnroot_input_list, MAX_NUM_INPUT_FILES);


	for (nFile = 0; nFile < NumFiles; nFile++) {
		chrpos = strstr(fnroot_input_list[nFile], test_str);
		*chrpos = '\0';
		strcpy(fnroot_input, fnroot_input_list[nFile]);
		if ((chrpos = strrchr(fnroot_input, '/')) == NULL)
			chrpos = fnroot_input;
		else
			chrpos++;
		fprintf(stdout, "GRID_FILE: %s.*\n", fnroot_input);
		if ((istat = GenStats(fnroot_input)) < 0) {
			exit(EXIT_ERROR_MISC);
		}
	}


	exit(0);

}



/** function to generate "traditional" statistics from grid file */

int GenStats(char* fnroot)
{


	int istat;

	FILE *fp_grid, *fp_hdr;
	GridDesc pdf_grid;

	Vect3D expectation;
	Mtrx3D covariance;

	int npts_ellipse;
	Vect2D* ellipse_array;



	// open grid file and header file and read header file

	if ((istat = OpenGrid3dFile(fnroot, &fp_grid, &fp_hdr,
			&pdf_grid, "", NULL)) < 0)
	{
		nll_puterr("ERROR opening grid file.");
		return(-1);
	}
	display_grid_param(&pdf_grid);


	// allocate grid buffer

	pdf_grid.buffer = AllocateGrid(&pdf_grid);
	if (pdf_grid.buffer == NULL) {
		nll_puterr(
"ERROR: allocating memory for 3D PDF grid buffer.\n");
		exit(EXIT_ERROR_MEMORY);
	}

	// create grid array access pointers

	pdf_grid.array = CreateGridArray(&pdf_grid);
	if (pdf_grid.array == NULL) {
		nll_puterr(
"ERROR: creating array for accessing 3D PDF grid buffer.\n");
		exit(EXIT_ERROR_MEMORY);
	}


	// read PDF grid

	if (istat = 
		ReadGrid3dBuf(&pdf_grid, fp_grid)
				< 0) {
		nll_puterr("ERROR: reading PDF grid from disk.\n");
		exit(EXIT_ERROR_IO);
	}


	// calculate expectation

	expectation = CalcExpectation(&pdf_grid);
	fprintf(stdout, "EXPECTATION: { x: %lf  y:%lf  z:%lf }\n",
		expectation.x, expectation.y, expectation.z);

	// calculate covariance matrix

	covariance = CalcCovariance(&pdf_grid, &expectation);
	fprintf(stdout, "COVARIANCE: {\n");
	fprintf(stdout, "   xx: %lf  xy:%lf  xz:%lf\n",
		covariance.xx, covariance.xy, covariance.xz);
	fprintf(stdout, "   yx: %lf  yy:%lf  yz:%lf\n",
		covariance.yx, covariance.yy, covariance.yz);
	fprintf(stdout, "   zx: %lf  zy:%lf  zz:%lf\n",
		covariance.zx, covariance.zy, covariance.zz);
	fprintf(stdout, "}\n");


	// generate error plots

	//GenErrorBars(stdout, &expectation, &covariance);
	GenCross(stdout, &expectation);


	// ellipses
	npts_ellipse = 100;
	ellipse_array = 
		(Vect2D *) malloc((size_t) npts_ellipse * sizeof(Vect2D));
	ellipse_array = Cov2Ellipse(covariance.xx, covariance.zz, 
			covariance.xz, expectation.x, expectation.z, 
		ellipse_array, npts_ellipse);
	Vect2DArray2GMT(stdout, ellipse_array, npts_ellipse, "XZ");
	ellipse_array = Cov2Ellipse(covariance.xx, covariance.yy, 
			covariance.xy, expectation.x, expectation.y, 
		ellipse_array, npts_ellipse);
	Vect2DArray2GMT(stdout, ellipse_array, npts_ellipse, "XY");
	ellipse_array = Cov2Ellipse(covariance.zz, covariance.yy, 
			covariance.zy, expectation.z, expectation.y, 
		ellipse_array, npts_ellipse);
	Vect2DArray2GMT(stdout, ellipse_array, npts_ellipse, "ZY");




	fprintf(stdout, "\n");


	// clean up

	FreeGrid(&pdf_grid);
	DestroyGridArray(&pdf_grid);
	CloseGrid3dFile(&fp_grid, &fp_hdr);
	free(ellipse_array);

	return(0);

}





/** function to generate error graphics calls */

void GenErrorBars(FILE* fp_io, Vect3D* pexpect, Mtrx3D* pcov)
{
	double barlen = 0.75;


	fprintf(fp_io, 
		"#XZ\npsxy $JVAL $RVAL -W1/0/0/0 -M -K -O << END >> %s.ps\n",
		"${POSTSCRIPT_NAME}");
	Err2GMT(fp_io, pexpect->x, pexpect->z, 
		2.0 * sqrt(pcov->xx), barlen, 2.0 * sqrt(pcov->zz), barlen);
	fprintf(fp_io, "END\n\n");

	fprintf(fp_io, 
		"#XY\npsxy $JVAL $RVAL -W1/0/0/0 -M -K -O << END >> %s.ps\n",
		"${POSTSCRIPT_NAME}");
	Err2GMT(fp_io, pexpect->x, pexpect->y, 
		2.0 * sqrt(pcov->xx), barlen, 2.0 * sqrt(pcov->yy), barlen);
	fprintf(fp_io, "END\n\n");

	fprintf(fp_io, 
		"#ZY\npsxy $JVAL $RVAL -W1/0/0/0 -M -K -O << END >> %s.ps\n",
		"${POSTSCRIPT_NAME}");
	Err2GMT(fp_io, pexpect->z, pexpect->y, 
		2.0 * sqrt(pcov->zz), barlen, 2.0 * sqrt(pcov->yy), barlen);
	fprintf(fp_io, "END\n\n");


}




/** function to generate error graphics calls */

void GenCross(FILE* fp_io, Vect3D* pexpect)
{
	double barlen = 0.75;


	fprintf(fp_io, 
		"#XZ\npsxy $JVAL $RVAL -W1/0/0/0 -M -K -O << END >> %s.ps\n",
		"${POSTSCRIPT_NAME}");
	Err2GMT(fp_io, pexpect->x, pexpect->z, 
		barlen, 0.0, barlen, 0.0);
	fprintf(fp_io, "END\n\n");

	fprintf(fp_io, 
		"#XY\npsxy $JVAL $RVAL -W1/0/0/0 -M -K -O << END >> %s.ps\n",
		"${POSTSCRIPT_NAME}");
	Err2GMT(fp_io, pexpect->x, pexpect->y, 
		barlen, 0.0, barlen, 0.0);
	fprintf(fp_io, "END\n\n");

	fprintf(fp_io, 
		"#ZY\npsxy $JVAL $RVAL -W1/0/0/0 -M -K -O << END >> %s.ps\n",
		"${POSTSCRIPT_NAME}");
	Err2GMT(fp_io, pexpect->z, pexpect->y, 
		barlen, 0.0, barlen, 0.0);
	fprintf(fp_io, "END\n\n");


}




//------------------------------------------------------------/
// Anthony Lomax           | email: lomax@faille.unice.fr     /
// UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax /
// 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        /
// 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        /
//------------------------------------------------------------/

