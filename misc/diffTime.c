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


/*   diffTime.c

	Program to calculate average "static" by 
		differencing travel times in 2 observation files for a single
		station and outputing average differences

*/

/*------------------------------------------------------------/ */
/* Anthony Lomax           | email: lomax@faille.unice.fr     / */
/* UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax / */
/* 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 26 25        / */
/* 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 26 10        / */
/*------------------------------------------------------------/ */


/*
	history:

	ver 01    26JAN1999  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

*/



#include "../src/GridLib.h"


/* defines */



/* globals */



/* functions */

int DiffTime(int , char **);



/*** program to difference 2 observatin files */

#define PNAME  "diffTime"


main(int argc, char *argv[])
{

	int istat, narg;


	/* set program name */

	strcpy(prog_name, PNAME);


	/* check command line for correct usage */

	fprintf(stderr, "\n%s Arguments: ", prog_name);
	for (narg = 0; narg < argc; narg++)
		fprintf(stderr, "<%s> ", argv[narg]);
	fprintf(stderr, "\n");

	if (argc != 6) {
		nll_puterr("ERROR wrong number of command line arguments.");
		disp_usage(PNAME, 
"<obsFile_1> <obsFile_2> staName phaseID staticShift\n (difference: obsFile_1 - obsFile_2 + staticShift)");
		exit(-1);
	}

	if ((istat = DiffTime(argc, argv)) < 0) {
		nll_puterr("ERROR differencing observation files.");
		exit(-1);
	}



	exit(0);

}



int DiffTime(int argc, char *argv[])
{
	int istat, narg;
	int nObsRead;
	char *cstat1, *cstat2;
	char staName[ARRIVAL_LABEL_LEN];
	char phaseID[PHASE_LABEL_LEN];
	char fn_obs1_in[FILENAME_MAX],  fn_obs2_in[FILENAME_MAX];
	FILE *fp_obs1_in, *fp_obs2_in;
	char line1[MAXLINE_LONG], line2[MAXLINE_LONG];
	double staticShift, time1, time2, difference;

	ArrivalDesc arr1, arr2;



	/* get command line parameters */
	strcpy(fn_obs1_in, argv[1]);
	strcpy(fn_obs2_in, argv[2]);

	strcpy(staName, argv[3]);
	strcpy(phaseID, argv[4]);

	sscanf(argv[5], "%lf", &staticShift);


	/* open observation files */

	if ((fp_obs1_in = fopen(fn_obs1_in, "r")) == NULL) {
		nll_puterr("ERROR: opening observation file 1.");
		return(-1);
	}

	if ((fp_obs2_in = fopen(fn_obs2_in, "r")) == NULL) {
		nll_puterr("ERROR: opening observation file 2.");
		return(-1);
	}


	nObsRead = 0;
	difference = 0.0;
	while (1) {

		/* read next observation */

		/* file 1 read next line */
		cstat1 = fgets(line1, MAXLINE_LONG, fp_obs1_in);
		/* file 2 read next line */
		cstat2 = fgets(line2, MAXLINE_LONG, fp_obs2_in);

		if (cstat1 == NULL)
			break;
	    	if ((istat = ReadArrival(line1, &arr1, IO_ARRIVAL_OBS)) == EOF)
			continue;
		else if (istat < 0) {
			continue;
		}
		if (cstat2 == NULL)
			break;
	    	if ((istat = ReadArrival(line2, &arr2, IO_ARRIVAL_OBS)) == EOF)
			continue;
		else if (istat < 0) {
			continue;
		}

		/* check for requested station */

		if (strcmp(arr1.label, staName) != 0 
				|| strcmp(arr2.label, staName) != 0)
			continue;
		if (strcmp(arr2.phase, phaseID) != 0 
				|| strcmp(arr2.phase, phaseID) != 0)
			continue;

		nObsRead++;
//fprintf(stderr, "1: %s", line1);
//fprintf(stderr, "2: %s", line2);
//fprintf(stderr, "%s %s ", staName, phaseID);

		/* update difference */

		time1 = arr1.sec 
			+ 60.0 * ((double) arr1.min + (60.0 * (double) arr1.hour));
		time2 = arr2.sec 
			+ 60.0 * ((double) arr2.min + (60.0 * (double) arr2.hour));
		difference += time1 - time2;
//fprintf(stderr, "  t1-t2: %lf ", time1 - time2);
		difference += staticShift;
//fprintf(stderr, "  t1-t2+static: %lf\n", time1 - time2 + staticShift);

	}


	fclose(fp_obs1_in);
	fclose(fp_obs2_in);

	difference /= (double) nObsRead;

	/* write summary */
	fprintf(stdout, "diffTime %s %s %d %lf\n", 
		staName, phaseID, nObsRead, difference);

	return(0);



}


/*------------------------------------------------------------/ */
/* Anthony Lomax           | email: lomax@faille.unice.fr     / */
/* UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax / */
/* 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        / */
/* 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        / */
/*------------------------------------------------------------/ */

