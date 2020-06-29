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


/*   mag_func_test.c

	Program to demonstrate calculating magnitudes through a function call.


*/

/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
161 Allee du Micocoulier, 06370 Mouans-Sartoux, France
tel: +33(0)493752502  e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


/*
	history:	(see also http://alomax.net/nlloc -> Updates)

	ver 01    05AUG2008  AJL  Original version

	see GridLib.c


.........1.........2.........3.........4.........5.........6.........7.........8

*/





#define PNAME  "mag_func_test"

#include "GridLib.h"
#include "ran1/ran1.h"
#include "velmod.h"
#include "GridMemLib.h"
#include "calc_crust_corr.h"
#include "phaseloclist.h"
#include "otime_limit.h"
#include "NLLocLib.h"




/** program to demonstrate getting magnitudes through a function call
 *
 *  See Makefile->mag_func_test for compile/link requirements
 *
 *  This demonstration program reads all input from a disk file.
 *
 */


#define NARGS_MIN 2
#define ARG_DESC "<input file>"

int main(int argc, char *argv[])
{

	char event_id[64];
	double event_lat, event_long, event_depth;
	char comp_id[64];
	double comp_lat, comp_long, comp_depth, comp_amp, comp_dur, comp_corr, comp_amp_fact;


	// set program name
	strcpy(prog_name, PNAME);

	// check command line for correct usage
	if (argc < NARGS_MIN) {
		disp_usage(prog_name, ARG_DESC);
		return(EXIT_ERROR_USAGE);
	}

	// set file names
	// control file name
	char fn_input[MAXLINE];
	strcpy(fn_input, argv[1]);
	char in_line[4*MAXLINE];


	SetConstants();

	// read input
	FILE* fp_input;
	if ((fp_input = fopen(fn_input, "r")) == NULL) {
		nll_puterr("FATAL ERROR: opening input file.");
		return(EXIT_ERROR_FILEIO);
	}

	// get magnitude parameters, sets Magnitude[0]
	// LOCMAG  ML_HB  1.82e7  1.678  -0.00514  20  1.8
	if ((fgets(in_line, 4*MAXLINE, fp_input) == NULL) || (GetNLLoc_Magnitude(strchr(in_line, ' ')) < 0)) {
		nll_puterr("ERROR: reading magnitude parameters.");
		return(EXIT_ERROR_FILEIO);
	}

	// loop over events
	while (1) {

		if ((fgets(in_line, 4*MAXLINE, fp_input) == NULL) || (sscanf(in_line, "%s %lf %lf %lf", event_id, &event_lat, &event_long, &event_depth) < 4))
			break;

		// read compenents
		int nvalues = 0;
		double mag_mean = 0.0;
		while ((fgets(in_line, 4*MAXLINE, fp_input) != NULL) && (sscanf(in_line, "%s %lf %lf %lf %lf %lf %lf %lf", comp_id, &comp_lat, &comp_long, &comp_depth, &comp_amp, &comp_dur, &comp_corr, &comp_amp_fact) == 8)) {

			//printf("%s %f %f %f %f %f %f %f", comp_id, comp_lat, comp_long, comp_depth, comp_amp, comp_dur, comp_corr, comp_amp_fact);

			if (strncmp(comp_id, "XXX", 3) == 0)
				break;

			// get epicentral distance
			double dist = GCDistance(event_lat, event_long, comp_lat, comp_long);

			/** ===========================================================================
			 *  call functions that get magnitude (see NLLocLib.c)
			*/
			double mag = 0.0;
			if (Magnitude[0].type == MAG_ML_HB) {
				mag = Calc_ML_HuttonBoore(comp_amp * comp_amp_fact * Magnitude[0].amp_fact_ml_hb, dist, event_depth, comp_corr, Magnitude[0].hb_n, Magnitude[0].hb_K, Magnitude[0].hb_Ro, Magnitude[0].hb_Mo);

			}
			else if (Magnitude[0].type == MAG_MD_FMAG) {
				mag = Calc_MD_FMAG(comp_dur, dist, event_depth, comp_corr, Magnitude[0].fmag_c1, Magnitude[0].fmag_c2, Magnitude[0].fmag_c3, Magnitude[0].fmag_c4, Magnitude[0].fmag_c5);
			}


			// write magnitude to stdout
			printf("%s %s  mag=%.3f\n", event_id, comp_id, mag);

			mag_mean += mag;
			nvalues++;
		}

		// write mean magnitude to stdout
		printf("%s  n_mag=%d  mean_mag=%.3f\n\n", event_id, nvalues, mag_mean / (double) nvalues);

	}


	// clean up
	fclose(fp_input);


	return(0);

}




