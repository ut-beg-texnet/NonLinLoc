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


/*   Resid2GMT.c

	Program to generate GMT psxy data from NLLoc .stat file

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

	ver 01    18Dec1998  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

*/



#define PNAME  "Resid2GMT"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


#define MAXLINE 100


// structures

typedef struct
{
	double x, y;
	char name[10];
} StaLoc;


// globals

// stations with vertical comp only commented out for S red plot

StaLoc staList[] = {
	{0.910400, -0.310300, "BKE"},
	{-2.564900, 0.634200, "OVO"},
	{-1.626800, -2.754700, "CSM"},
	{0.419700, 3.652700, "SMC"},
//	{-0.617500, -4.384300, "CPV"},
	{0.658100, -0.291700, "NBK"},
	{0.504000, -1.051000, "OBK"},
//	{-4.371500, 2.356400, "SSB"},
//	{9.864500, 11.060100, "NL9"},
//	{4.811400, -7.884300, "PS9"},
//	{-6.953300, -1.865800, "HR9"},
	{-0.589300, -1.606600, "BKS"},
//	{3.882100, -1.532500, "TRZ"},
	{0.924700, -2.736200, "FTC"},
	{-3.000200, -1.551000, "TDG"},
	{0.111500, 0.949000, "BKN"},
	{-1.191800, -0.365800, "SGV"},
	{3.221800, 1.856400, "OTV"},
	{11.868900, 10.649000, "VIS"},
	{-0.507900, 0.026800, "CRT"},
	{-6.613000, -0.695400, "POR"},
//	{3.501600, 7.757600, "SCI"},
	{-1.388100, -0.995400, "BAF"},
	{0.0, 0.0, "$$$"}
};


// functions

int Convert2MapGMT(char* fn_nlloc_stat, char* fn_out);
StaLoc *findStaLoc(char *staName, StaLoc stalist[]);


// Program to generate GMT psxy data from NLLoc .stat file


main(int argc, char *argv[])
{

	int istat;
	int c;
	char fn_nlloc_stat[MAXLINE];
	char fn_out[MAXLINE];


	/* read each set of residuals */

	while (1)
	{
		/* get stat file name */
		if ((istat = fscanf(stdin, "%s", fn_nlloc_stat)) != 1 || istat == EOF)
			break;
		/* get output file name */
		if ((istat = fscanf(stdin, "%s", fn_out)) != 1 || istat == EOF)
			break;

		printf("Processing file: %s\n", fn_nlloc_stat);

		if (Convert2MapGMT(fn_nlloc_stat, fn_out) < 0)
			break;
	}

	exit(0);

}



/** function to convert residuals list to gmt psxy (residual at sta x,y) and pstext file  */


#define RESID_MIN 0.02

int Convert2MapGMT(char* fn_nlloc_stat, char* fn_out)
{
	int istat, nRdg, nresid;
	int c = 0;
	int ifound_readings;
	char label[MAXLINE];
	char fname[MAXLINE];
	FILE *fp_in, *fp_xy_out_P, *fp_xy_out_S, *fp_out, *fp_text_out;

	char staName[20];
	char phase[20];
	double resid, resid_plot;

	StaLoc *staloc;


	if ((fp_in = fopen(fn_nlloc_stat, "r")) == NULL) {
		fprintf(stderr,
			"ERROR: Cannot open input file <%s>.\n",
				fn_nlloc_stat);
		return(-1);
    	}
	sprintf(fname, "%s.P.resid.map.xy", fn_out);
	if ((fp_xy_out_P = fopen(fname, "w")) == NULL) {
		fprintf(stderr,
			"ERROR: Cannot open xy output file <%s>.\n",
				fname);
		return(-1);
    	}
	sprintf(fname, "%s.S.resid.map.xy", fn_out);
	if ((fp_xy_out_S = fopen(fname, "w")) == NULL) {
		fprintf(stderr,
			"ERROR: Cannot open xy output file <%s>.\n",
				fname);
		return(-1);
    	}
	sprintf(fname, "%s.resid.map.text", fn_out);
	if ((fp_text_out = fopen(fname, "w")) == NULL) {
		fprintf(stderr,
			"ERROR: Cannot open text output file <%s>.\n",
				fname);
		return(-1);
    	}


	// find beginning of Total Phase Corrections list
	do {
		istat = fscanf(fp_in, "%s", label);
	} while (istat != EOF && strcmp(label, "Total") != 0);


	/* read residuals */

	nRdg = 0;
	ifound_readings = 0;
	do {

		istat = fscanf(fp_in, "%s", label);
		if (istat == EOF)
			break;
		if (strcmp(label, "LOCDELAY") != 0) {
			if (ifound_readings)	// do not pass to next set of readings
				break;
			while((c = fgetc(fp_in)) != '\n' && c !=EOF)
				;
		} else {
			ifound_readings = 1;
			istat = fscanf(fp_in, "%s %s %d %lf\n", staName, phase, &nresid, &resid);
			if (istat == EOF)
				break;
			printf("%s %s %d %lf\n", staName, phase, nresid, resid);

			if ((staloc = findStaLoc(staName, staList)) != NULL) {
				if (strcmp(phase, "P") == 0)
					fp_out = fp_xy_out_P;
				else
					fp_out = fp_xy_out_S;
				if (resid > 0.0) {
					resid_plot = resid < RESID_MIN ? RESID_MIN : resid;
					fprintf(fp_out,
"psxy $JVAL $RVAL -Sc -W8/0/0/255 -K -O  << END >> ${POSTSCRIPT_NAME}\n%lf %lf %lf\nEND\n",
						staloc->x, staloc->y, fabs(resid_plot * 5.0));
				} else {
					resid_plot = resid > -RESID_MIN ? -RESID_MIN : resid;
					fprintf(fp_out,
"psxy $JVAL $RVAL -Sd -W8/255/0/0 -K -O  << END >> ${POSTSCRIPT_NAME}\n%lf %lf %lf\nEND\n",
						staloc->x, staloc->y, fabs(resid_plot * 5.0));
				}

				fprintf(fp_text_out, "%lf %lf 10 0 4 6 %s\n",
					staloc->x, staloc->y, staName);
				nRdg++;
			}

		}

	} while (c != EOF);


	fclose(fp_xy_out_P);
	fclose(fp_xy_out_S);
	fclose(fp_text_out);

	return(0);

}


/** function to find a station name in a StaList */

StaLoc *findStaLoc(char *staName, StaLoc stalist[])
{

	int n;

	n = 0;
	while (strcmp(stalist[n].name, "$$$") != 0) {
		if (strcmp(stalist[n].name, staName) == 0)
			return(&(stalist[n]));
		n++;
	}

	return(NULL);

}




//------------------------------------------------------------/
// Anthony Lomax           | email: lomax@faille.unice.fr     /
// UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax /
// 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        /
// 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        /
//------------------------------------------------------------/

