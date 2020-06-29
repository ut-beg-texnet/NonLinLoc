/*
 * Copyright (C) 2004 Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
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



/*   Corr2MapGMT.c

	Program to generate GMT psxy data from NLLoc .stat_totcorr file

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



#define PNAME  "ConvertCorr2MapGMT2GMT"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "../src/GridLib.h"




// globals



// functions

int ConvertCorr2MapGMT(char* fn_nlloc_stat, char* fn_stations, char* fn_out, double res_max);
SourceDesc* findStaLoc(char *staName);


// Program to generate GMT psxy data from NLLoc .stat file


main(int argc, char *argv[])
{

	int istat, narg;
	int c;
	char fn_nlloc_stat[MAXLINE_LONG];
	char fn_stations[MAXLINE_LONG];
	char fn_out[MAXLINE_LONG];
	double res_max;


	fprintf(stdout, "%s Arguments: ", prog_name);
	for (narg = 0; narg < argc; narg++)
		fprintf(stdout, "<%s> ", argv[narg]);
	fprintf(stdout, "\n");

	if (argc < 5) {
		nll_puterr("ERROR wrong number of command line arguments.");
		disp_usage(PNAME, "<stat_corr_file> <station_list> <output_file_root> <res_max>");
		exit(-1);
	}

	/* read each set of residuals */

	/* get stat file name */
	strcpy(fn_nlloc_stat, argv[1]);
	/* get station file name */
	strcpy(fn_stations, argv[2]);
	/* get output file name */
	strcpy(fn_out, argv[3]);
	/* get output file name */
	sscanf(argv[4], "%lf", &res_max);

	printf("Processing file: %s\n", fn_nlloc_stat);

	ConvertCorr2MapGMT(fn_nlloc_stat, fn_stations, fn_out, res_max);

	exit(0);

}



/** function to convert residuals list to gmt psxy (residual at sta x,y) and pstext file  */


#define RESID_MIN 0.02

int ConvertCorr2MapGMT(char* fn_nlloc_stat, char* fn_stations, char* fn_out, double res_max)
{
	int istat, nRdg, nresid;
	int c = 0;
	int ifound_readings;
	char line[MAXLINE_LONG];
	char label[MAXLINE_LONG];
	char fname[MAXLINE_LONG];
	FILE *fp_in, *fp_sta, *fp_xy_out_Ppos, *fp_xy_out_Spos, *fp_xy_out_Pneg, *fp_xy_out_Sneg, *fp_out, *fp_text_out;
	FILE *fp_xy_out_P, *fp_xy_out_S;
	FILE *fp_xy_out_Pg;

	char staName[20];
	char phase[20];
	double resid, resid_plot;
	double stddev;      // 20190822 AJL - added

	SourceDesc* staloc;

	// read stations
	printf("Reading station file <%s>.\n", fn_stations);
	if ((fp_sta = fopen(fn_stations, "r")) == NULL) {
		fprintf(stderr,
			"ERROR: Cannot open station file <%s>.\n", fn_stations);
		return(-1);
    	}
	while (fgets(line, MAXLINE_LONG, fp_sta) != NULL)
		if ((istat = GetNextSource(strchr(line, ' '))) < 0) {
			nll_puterr("ERROR: reading source params:");
			nll_puterr(line);
		}
	fclose(fp_sta);


	// open files
	printf("Opening sta_corr file <%s>.\n", fn_nlloc_stat);
	if ((fp_in = fopen(fn_nlloc_stat, "r")) == NULL) {
		fprintf(stderr,
			"ERROR: Cannot open input file <%s>.\n", fn_nlloc_stat);
		return(-1);
    	}

	// open output files
	printf("Opening output files <%s>.\n", fn_out);

	// all res (for further GMT processing)
	sprintf(fname, "%s.P.resid.map.xyc", fn_out);
	if ((fp_xy_out_P = fopen(fname, "w")) == NULL) {
		fprintf(stderr,
			"ERROR: Cannot open xy output file <%s>.\n", fname);
		return(-1);
    	}
	sprintf(fname, "%s.S.resid.map.xyc", fn_out);
	if ((fp_xy_out_S = fopen(fname, "w")) == NULL) {
		fprintf(stderr,
			"ERROR: Cannot open xy output file <%s>.\n", fname);
		return(-1);
    	}
	sprintf(fname, "%s.Pg.resid.map.xyc", fn_out);
	if ((fp_xy_out_Pg = fopen(fname, "w")) == NULL) {
		fprintf(stderr,
			"ERROR: Cannot open xy output file <%s>.\n", fname);
		return(-1);
    	}

	// pos or neg P and S (for SeismicityViewer)
	sprintf(fname, "%s.Ppos.resid.map.xyz", fn_out);
	if ((fp_xy_out_Ppos = fopen(fname, "w")) == NULL) {
		fprintf(stderr,
			"ERROR: Cannot open xy output file <%s>.\n", fname);
		return(-1);
    	}
	fprintf(fp_xy_out_Ppos, "> GMT_LATLONDEPTH\n");
	sprintf(fname, "%s.Spos.resid.map.xyz", fn_out);
	if ((fp_xy_out_Spos = fopen(fname, "w")) == NULL) {
		fprintf(stderr,
			"ERROR: Cannot open xy output file <%s>.\n", fname);
		return(-1);
    	}
	fprintf(fp_xy_out_Spos, "> GMT_LATLONDEPTH\n");
	sprintf(fname, "%s.Pneg.resid.map.xyz", fn_out);
	if ((fp_xy_out_Pneg = fopen(fname, "w")) == NULL) {
		fprintf(stderr,
			"ERROR: Cannot open xy output file <%s>.\n", fname);
		return(-1);
    	}
	fprintf(fp_xy_out_Pneg, "> GMT_LATLONDEPTH\n");
	sprintf(fname, "%s.Sneg.resid.map.xyz", fn_out);
	if ((fp_xy_out_Sneg = fopen(fname, "w")) == NULL) {
		fprintf(stderr,
			"ERROR: Cannot open xy output file <%s>.\n", fname);
		return(-1);
    	}
	fprintf(fp_xy_out_Sneg, "> GMT_LATLONDEPTH\n");
	// text (for SeismicityViewer)
	sprintf(fname, "%s.resid.map.text", fn_out);
	if ((fp_text_out = fopen(fname, "w")) == NULL) {
		fprintf(stderr,
			"ERROR: Cannot open text output file <%s>.\n", fname);
		return(-1);
    	}


	// find beginning of Total Phase Corrections list
	//do {
	//	istat = fscanf(fp_in, "%s", label);
	//} while (istat != EOF && strcmp(label, "Total") != 0);


	/* read residuals */

	printf("Reading and converting corrections.\n");
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
			istat = fscanf(fp_in, "%s %s %d %lf %lf\n", staName, phase, &nresid, &resid, &stddev);
			if (istat == EOF)
				break;
			printf("%s %s %d %lf %lf\n", staName, phase, nresid, resid, stddev);

			if (fabs(resid) > res_max)
				continue;

			if ((staloc = findStaLoc(staName)) != NULL) {
				if (resid > 0.0) {
					if (strcmp(phase, "P") == 0)
						fp_out = fp_xy_out_Ppos;
					else
						fp_out = fp_xy_out_Spos;
					resid_plot = resid < RESID_MIN ? RESID_MIN : resid;
					fprintf(fp_out, "%lf %lf %lf %lf\n",
						staloc->dlat, staloc->dlong, staloc->depth, fabs(resid_plot * 20.0));
				} else {
					if (strcmp(phase, "P") == 0)
						fp_out = fp_xy_out_Pneg;
					else
						fp_out = fp_xy_out_Sneg;
					resid_plot = resid > -RESID_MIN ? -RESID_MIN : resid;
					fprintf(fp_out, "%lf %lf %lf %lf\n",
						staloc->dlat, staloc->dlong, staloc->depth, fabs(resid_plot * 20.0));
				}

				fprintf(fp_text_out, "%lf %lf 10 0 4 6 %s\n",
					staloc->dlat, staloc->dlong, staName);

				if (strcmp(phase, "P") == 0)
					fprintf(fp_xy_out_P, "%lf %lf %lf %d\n",
						staloc->dlat, staloc->dlong, resid, nresid);
				if (strcmp(phase, "S") == 0)
					fprintf(fp_xy_out_S, "%lf %lf %lf %d\n",
						staloc->dlat, staloc->dlong, resid, nresid);
				// additional phases
				if (strcmp(phase, "Pg") == 0)
					fprintf(fp_xy_out_Pg, "%lf %lf %lf %d\n",
						staloc->dlat, staloc->dlong, resid, nresid);

				nRdg++;
			}

		}

	} while (c != EOF);


	fclose(fp_xy_out_Ppos);
	fclose(fp_xy_out_Spos);
	fclose(fp_xy_out_Pneg);
	fclose(fp_xy_out_Sneg);
	fclose(fp_text_out);
	fclose(fp_xy_out_P);
	fclose(fp_xy_out_S);
	fclose(fp_xy_out_Pg);

	return(0);

}


/** function to find a station name in a StaList */

SourceDesc* findStaLoc(char *staName)
{

	SourceDesc* pstation;

	pstation = FindSource(staName);

	return(pstation);

}




//------------------------------------------------------------/
// Anthony Lomax           | email: lomax@faille.unice.fr     /
// UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax /
// 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        /
// 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        /
//------------------------------------------------------------/

