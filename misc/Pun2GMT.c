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


/*   Pun2GMT.c

	Program to generate GMT command and psxy data from hypo71 .pun file

	output in GMT readable format

*/

//------------------------------------------------------------/
// Anthony Lomax           | email: lomax@faille.unice.fr     /
// UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax /
// 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        /
// 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        /
//------------------------------------------------------------/


/*
	history:

	ver 01    10Nov1997  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

*/



#define PNAME  "Pun2GMT"

#include "../src/GridLib.h"
#include "GridGraphLib.h"

#define PLOT_WIDTH 20.0
#define PLOT_HEIGHT 28.8
#define TITLE_FONT_SIZE 18
#define TITLE_FONT 4
#define HYPO_FONT_SIZE 14
#define HYPO_FONT 4
#define STA_FONT_SIZE 10
#define STA_FONT 4
#define ANNOTATION_FONT_SIZE 7
#define ANNOTATION_FONT 4


// globals


// functions

int ReadPun2GMT_Input(FILE* );
int ConvertPun2GMT(char* , char* , char* , char* );


/*** program to generate  GTM command files for plotting grid data */


main(int argc, char *argv[])
{

	int istat, narg;
	char fngrid_control[MAXLINE], fn_input[MAXLINE], fn_outroot[MAXLINE];
	FILE *fp_grid_control, *fp_outroot, *fp_gmt;
	char filename[MAXLINE], fn_xy_out[MAXLINE];


	// set program name

	strcpy(prog_name, PNAME);


	// check command line for correct usage

	fprintf(stdout, "\n%s Arguments: ", prog_name);
	for (narg = 0; narg < argc; narg++)
		fprintf(stdout, "<%s> ", argv[narg]);
	fprintf(stdout, "\n");

	if (argc != 4)
	{
		nll_puterr("ERROR wrong number of command line arguments.");
		disp_usage(PNAME, 
		    "<controlfile> <punfile> <outroot>");
		exit(-1);
	}


	// read command line arguments

	strcpy(fngrid_control, argv[1]);
	strcpy(fn_input, argv[2]);
	strcpy(fn_outroot, argv[3]);



	// set constants

	SetConstants();


	// read control file

	if ((fp_grid_control = fopen(fngrid_control, "r")) == NULL) {
		nll_puterr("ERROR opening control file.");
		exit(EXIT_ERROR_FILEIO);
	}


	if ((istat = ReadPun2GMT_Input(fp_grid_control)) < 0) {
		nll_puterr("ERROR reading control file.");
		exit(EXIT_ERROR_FILEIO);
	}


	// open GMT command file

	sprintf(filename, "%s.gmt", fn_outroot);
	if ((fp_gmt = fopen(filename, "w")) == NULL) {
		nll_puterr("ERROR opening output GMT script file.");
		exit(EXIT_ERROR_FILEIO);
	}


	// convert .pun to xy

	ConvertPun2GMT(fn_input, "XZ", fn_outroot, fn_xy_out);
	fprintf(fp_gmt, 
		"psxy %s $JVAL $RVAL -W1/0/0/0 -M -K -O >> %s.ps\n\n",
		fn_xy_out, "${POSTSCRIPT_NAME}");
	ConvertPun2GMT(fn_input, "XY", fn_outroot, fn_xy_out);
	fprintf(fp_gmt, 
		"psxy %s $JVAL $RVAL -W1/0/0/0 -M -K -O >> %s.ps\n\n",
		fn_xy_out, "${POSTSCRIPT_NAME}");
	ConvertPun2GMT(fn_input, "ZY", fn_outroot, fn_xy_out);
	fprintf(fp_gmt, 
		"psxy %s $JVAL $RVAL -W1/0/0/0 -M -K -O >> %s.ps\n\n",
		fn_xy_out, "${POSTSCRIPT_NAME}");


	exit(0);

}



/** function to convert hypo71 .pun file to gmt psxyz file and command */

int ConvertPun2GMT(char* fn_input, char* orientation, char* fn_outroot, 
	char* fn_xy_out)
{
	int istat;
	FILE *fp_input, *fp_xy_out, *fp_lonlat_out;
	char  line_in[2*MAXLINE], fn_lonlat_out[MAXLINE];
	int or_xy, or_xz, or_yz, or_zy;
	double vlat, vlong, dlat, dlong, vdepth;
	int ilat, ilong;
	double xgrid, ygrid;
	double errx, erry, errz;
	double barlen = 0.75;


	if ((fp_input = fopen(fn_input, "r")) == NULL) {
		fprintf(stderr, "ERROR: Cannot open .pun file <%s>.\n",
				fn_input);
		return(-1);
    	}

	sprintf(fn_xy_out, "%s.%s", fn_outroot, orientation);
	if ((fp_xy_out = fopen(fn_xy_out, "w")) == NULL) {
		fprintf(stderr, 
			"ERROR: Cannot open xy output file <%s>.\n",
				fn_xy_out);
		return(-1);
    	}
	sprintf(fn_lonlat_out, "%s.lonlat.%s", fn_outroot, orientation);
	if ((fp_lonlat_out = fopen(fn_lonlat_out, "w")) == NULL) {
		fprintf(stderr, 
			"ERROR: Cannot open lonlat output file <%s>.\n",
				fn_lonlat_out);
		return(-1);
    	}

	or_xy = or_xz = or_yz = or_zy = 0;
	if (strcmp(orientation, "XY") == 0)
		or_xy = 1;
	else if (strcmp(orientation, "XZ") == 0)
		or_xz = 1;
	else if (strcmp(orientation, "YZ") == 0)
		or_yz = 1;
	else if (strcmp(orientation, "ZY") == 0)
		or_zy = 1;


	// skip header line

	if (fgets(line_in, 2*MAXLINE, fp_input) == NULL)
		nll_puterr("ERROR: reading header line in .pun file");


	// convert each hypocenter

	fprintf(fp_xy_out, ">\n");
	fprintf(fp_lonlat_out, ">\n");

	while(fgets(line_in, 2*MAXLINE, fp_input) != NULL) {
		dlat = dlong = vdepth = errx = erry = errz = 2.0e2;
		ilat = ilong = 0;
		if ((istat = sscanf(line_in, 
"%*d %*d %*f %d-%lf %d-%lf %lf %*d %*d %*f %*f %*f %*f %*s %lf %lf %lf", 
				&ilat, &dlat, &ilong, &dlong, &vdepth, 
				&errx, &erry, &errz)) != 8) {
			nll_puterr("WARNING: Error reading line in .pun file:");
			nll_puterr(line_in);
		}
		vlat = (double) ilat + dlat / 60.0;
		vlong = ilong + dlong / 60.0;
		latlon2rect(vlat, vlong, &xgrid, &ygrid);
		if (or_xy) {
			Err2GMT(fp_xy_out, xgrid, ygrid, 
				errx, barlen, erry, barlen);
			Err2GMT(fp_lonlat_out, vlong, vlat, 
				errx / (c111 * cos(rpd * vlat)), 
				barlen / (c111 * cos(rpd * vlat)),
				erry / c111, barlen / c111);
		} else if (or_xz) {
			Err2GMT(fp_xy_out, xgrid, vdepth, 
				errx, barlen, errz, barlen);
			Err2GMT(fp_lonlat_out, vlong, vdepth, 
				errx / (c111 * cos(rpd * vlat)), 
				barlen / (c111 * cos(rpd * vlat)),
				errz, barlen);
		} else if (or_yz) {
			Err2GMT(fp_xy_out, ygrid, vdepth, 
				erry, barlen, errz, barlen);
			Err2GMT(fp_lonlat_out, vlat, vdepth, 
				erry / c111, barlen / c111, errz, barlen);
		} else if (or_zy) {
			Err2GMT(fp_xy_out, vdepth, ygrid, 
				errz, barlen, erry, barlen);
			Err2GMT(fp_lonlat_out, vdepth, vlat, 
				errz, barlen, erry / c111, barlen / c111);
		}
	}

	fclose(fp_input);
	fclose(fp_xy_out);
	fclose(fp_lonlat_out);

	return(0);

}






/*** function to read input file */

int ReadPun2GMT_Input(FILE* fp_input)
{
	int istat, iscan;
	char param[MAXLINE];
	char line[MAXLINE];

	int flag_control = 0, flag_trans = 0;



	// read each input line

	while (fgets(line, MAXLINE, fp_input) != NULL) { 

		istat = -1;

		//read parameter line

		if ((iscan = sscanf(line, "%s", param)) < 0 )
			continue;

		// skip comment line or white space

		if (strncmp(param, "#", 1) == 0 || isspace(param[0]))
			istat = 0;


		// read control params

		if (strncmp(param, "CONTROL", 6) == 0)
			if ((istat = get_control(strchr(line, ' '))) < 0) 
				nll_puterr("Error reading control params.");
			else
				flag_control = 1;


		//read transform params

		if (strncmp(param, "TRANS", 5) == 0)
    			if ((istat = get_transform(strchr(line, ' '))) < 0)
			    nll_puterr("ERROR reading transformation parameters.");
			else
				flag_trans = 1;


		// unrecognized input

		if (istat < 0 && message_flag > 1) {
			fprintf(stdout, "Skipping input: %s", line);
		}

	}


	// check for missing input

	if (!flag_control) 
		nll_puterr("ERROR no control (CONTROL) params read.");
	if (!flag_trans) 
		nll_puterr("ERROR no transformation (TRANS) params read.");

	
	return (flag_control + flag_trans - 1);
}





//------------------------------------------------------------/
// Anthony Lomax           | email: lomax@faille.unice.fr     /
// UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax /
// 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        /
// 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        /
//------------------------------------------------------------/

