/*
* Copyright (C) 1999-2005 Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
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


/*   hyp2bull.c

	Program to convert NLLoc .hyp to a simple bulletin format

*/

/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
161 Allee du Micocoulier, 06370 Mouans-Sartoux, France
tel: +33(0)493752502  e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


/*
	history:

	ver 01    29Jul2008  AJL  Original version


*/



#include "../src/GridLib.h"


/* defines */

#define MAX_NUM_INPUT_FILES 5000


/* globals */



/* functions */

int hyp2bull(int , char **);



/*** program to sum event scatter files */

#define PNAME  "hyp2bull"


int main(int argc, char *argv[])
{

	int istat, narg;


	/* set program name */

	strcpy(prog_name, PNAME);


	/* check command line for correct usage */

	fprintf(stdout, "\n%s Arguments: ", prog_name);
	for (narg = 0; narg < argc; narg++)
		fprintf(stdout, "<%s> ", argv[narg]);
	fprintf(stdout, "\n");

	if (argc < 3) {
		nll_puterr("ERROR wrong number of command line arguments.");
		disp_usage(PNAME,
			   "<hyp_file_list> <output_file> ");
		exit(-1);
	}

	if ((istat = hyp2bull(argc, argv)) < 0) {
		nll_puterr("ERROR converting hyp file.");
		exit(-1);
	}



	exit(0);

}



int hyp2bull(int argc, char *argv[])
{
	int istat;
	int nLocWritten, nLocRead;
	int numFiles, nFile;
	char fn_out[FILENAME_MAX];
	char fn_hyp_in[FILENAME_MAX];
	char label_last[PHASE_LABEL_LEN];
	FILE *fp_out, *fp_hypo;

	GridDesc locgrid;
	HypoDesc Hypo;
	ArrivalDesc *parr;
	int narr;	

	char fn_hyp_in_list[MAX_NUM_INPUT_FILES][FILENAME_MAX];

	
	SetConstants();

	/* get command line parameters */
	strcpy(fn_hyp_in, argv[1]);
	sprintf(fn_out, "%s", argv[2]);



	/* check for wildcards in input file name */
	strcat(fn_hyp_in, ".hyp");
	if ((numFiles = ExpandWildCards(fn_hyp_in, fn_hyp_in_list, MAX_NUM_INPUT_FILES)) < 1) {
		nll_puterr("ERROR: no matching .hyp files found.");
		return(-1);
	}
	if (numFiles >= MAX_NUM_INPUT_FILES) {
		sprintf(MsgStr, "WARNING: maximum number of event files exceeded, only first %d will be processed.", MAX_NUM_INPUT_FILES);
		nll_puterr(MsgStr);
	}

	/* open ascii hypocenter output file */
	if ((fp_out = fopen(fn_out, "w")) == NULL) {
		nll_puterr("ERROR: opening bulletin output file.");
		return(-1);
	}



	nLocWritten = 0;
	nLocRead = 0;
	for (nFile = 0; nFile < numFiles; nFile++) {

		fprintf(OUT_LEVEL_1, "Adding location <%s>.\n", fn_hyp_in_list[nFile]);

		/* open hypocenter file */
		//sprintf(strstr(fn_hyp_in_list[nFile], test_str), "\0");
		if ((fp_hypo = fopen(fn_hyp_in_list[nFile], "r")) == NULL)
		{
			nll_puterr2("ERROR: opening hypocenter file, ignoring event, file",
				fn_hyp_in_list[nFile]);
			continue;
		}
		nLocRead++;
		// loop over events in file

		while (1) {

			istat = GetHypLoc(fp_hypo, NULL, &Hypo, Arrival, &NumArrivals, 1, &locgrid, 0);
			if (istat == EOF) {
				break;
			}
			if (istat < 0) {
				nll_puterr2("ERROR: reading hypocenter file, ignoring event, file",
					fn_hyp_in_list[nFile]);
				break;
			}

			if (strcmp(Hypo.locStat, "ABORTED") == 0) {
			//nll_puterr("WARNING: location ABORTED, ignoring event");
				continue;
			} else if (strcmp(Hypo.locStat, "REJECTED") == 0) {
			//nll_puterr("WARNING: location REJECTED, ignoring event");
				continue;
			}

			nLocWritten++;

			// write hypocenter
			fprintf(fp_out, 
				"evenement      date        heure           lat.         long.       prof.     rms\n");
//    3       11/02/2006    20:26:30.16     43.62 N       5.60 E      7 km     0.10 sec

			fprintf(fp_out, "  %3d", nLocWritten);
			fprintf(fp_out, "       %2.2d/%2.2d/%4.4d", Hypo.day, Hypo.month, Hypo.year);
			fprintf(fp_out, "    %2.2d:%2.2d:%5.2f", Hypo.hour, Hypo.min, Hypo.sec);
			fprintf(fp_out, "    %6.2f N", Hypo.dlat);
			fprintf(fp_out, "    %7.2f E", Hypo.dlong);
			fprintf(fp_out, "    %3.0f km", Hypo.depth);
			fprintf(fp_out, "    %6.3f sec", Hypo.rms);
			
			fprintf(fp_out, "\n");
			
			for (narr = 0; narr < NumArrivals; narr++) {
				Arrival[narr].obs_time = (long double) Arrival[narr].sec
						+ 60.0L * ((long double) Arrival[narr].min
						+ 60.0L * (long double) Arrival[narr].hour);
			}
			if (SortArrivalsTime(Arrival, NumArrivals) < 0)
				nll_puterr("ERROR: sorting arrivals by time.");
			if (SortArrivalsDist(Arrival, NumArrivals) < 0)
				nll_puterr("ERROR: sorting arrivals by distance.");

			strcpy(label_last, "XXX");
			for (narr = 0; narr < NumArrivals; narr++) {
				parr = Arrival + narr;
				if (strcmp(parr->phase, "S") == 0) {	// S
					if (strcmp(parr->label, label_last) == 0) {	// same statinon, S after P
						fprintf(fp_out, "      S:");
						strcpy(label_last, "$$$");
					} else {	// lone S
						fprintf(fp_out, "\n");
						fprintf(fp_out, "    %s", parr->label);
						fprintf(fp_out, "                                            S:");
						strcpy(label_last, "$$$");
					}
					fprintf(fp_out, "  %2.2d:%2.2d:%5.2f", parr->hour, parr->min, parr->sec);
				} else {	// P
					fprintf(fp_out, "\n");
					fprintf(fp_out, "    %s", parr->label);
					fprintf(fp_out, "        %s%s:", parr->phase, strcmp(parr->first_mot, "?") == 0 ? " " : parr->first_mot);
					strcpy(label_last, parr->label);
					fprintf(fp_out, "  %2.2d:%2.2d:%6.3f", parr->hour, parr->min, parr->sec);
				}
				fprintf(fp_out, " +/-%5.2f sec", parr->error);
			}

			if (strcmp(label_last, "$$$") == 0)	// S after P
				fprintf(fp_out, "\n");
			fprintf(fp_out, "\n");
			fprintf(fp_out, "\n");

		}

	}


	fclose(fp_out);

	/* write message */
	fprintf(stdout,
		"%d locations read, %d written to ascii bulletin file <%s>\n",
		nLocRead, nLocWritten, fn_out);


	return(0);

}



