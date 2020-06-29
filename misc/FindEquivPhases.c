/*
 * Copyright (C) 2004 Anthony Lomax <www.alomax.net, anthony@alomax.net>
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


/*   FindEquivPhases.c

	Program to find and write a (sub-)set of phases from a hypocenter that most closely match
	the phases of another hypocenter.  Designed to find sub-set of phases in a modern event that most
	closely match the phases in an early-instrumental event.

*/


/*
	history:

	ver 01    09Nov2004  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

*/



#include "../src/GridLib.h"


/* defines */


/* globals */



/* functions */

int doFindEquivPhases(int argc, char** argv);
int IsEquivPhase(char *target_phase, char *source_phase);
char lastLegType(char *phase);



/*** program to sum event scatter files */

#define PNAME  "FindEquivPhases"


int main(int argc, char** argv)
{

	int narg;
	int nPhases;


	/* set program name */

	strcpy(prog_name, PNAME);


	/* check command line for correct usage */

	fprintf(stdout, "%s Arguments: ", prog_name);
	for (narg = 0; narg < argc; narg++)
		fprintf(stdout, "<%s> ", argv[narg]);
	fprintf(stdout, "\n");

	if (argc < 5) {
		nll_puterr("ERROR wrong number of command line arguments.");
		disp_usage(PNAME,
"<nll_hyp_target> <nll_hyp_phs_source> <nll_hyp_out> dist_max_fract");
		exit(0);
	}

	SetConstants();
	message_flag = 1;
	DispProgInfo();
	message_flag = 0;

	if ((nPhases = doFindEquivPhases(argc, argv)) < 0) {
		nll_puterr("ERROR doing Phase Equivalence process.");
		exit(0);
	}



	exit(nPhases);

}



int doFindEquivPhases(int argc, char** argv)
{


	int istat, narg;
	int narr0, narr1;;
	int nPhaseMatched, nAttempt;

	char fn_hyp_target[FILENAME_MAX];
	char fn_hyp_phs_source[FILENAME_MAX];
	char fn_hyp_out[FILENAME_MAX];

	FILE *fp_hyp_target, *fp_hyp_phs_source, *fp_hyp_out;

	HypoDesc Hypo_target, Hypo_phs_source, Hypo_out;
	GridDesc Grid_target, Grid_phs_source, Grid_out;

	int NumArrivals_target;
	ArrivalDesc Arrival_target[X_MAX_NUM_ARRIVALS];
	int NumArrivals_phs_source;
	ArrivalDesc Arrival_phs_source[X_MAX_NUM_ARRIVALS];
	int NumArrivals_out;
	ArrivalDesc Arrival_out[X_MAX_NUM_ARRIVALS];

	int isSelected[X_MAX_NUM_ARRIVALS];


	double dist_max_fract, xval, yval, dist, dist_min, dist_epi;
	int index_min;
	char *target_phase, *source_phase;


	// get command line parameters

	narg = 1;

	strcpy(fn_hyp_target, argv[narg]);
	narg++;
	strcpy(fn_hyp_phs_source, argv[narg]);
	narg++;
	strcpy(fn_hyp_out, argv[narg]);
	narg++;
	dist_max_fract = VERY_LARGE_DOUBLE;
	sscanf(argv[narg], "%lf", &dist_max_fract);
	narg++;


	// associate each NLL hypo file

	nAttempt = 0;
	nPhaseMatched = 0;


	// open "target" NLL hypocenter file
	istat = GetHypLoc(NULL, fn_hyp_target, &Hypo_target, Arrival_target, &NumArrivals_target, 1, &Grid_target, 0);
	if (istat < 0) {
		nll_puterr2("ERROR: opening target NLL hypocenter-phase file.", fn_hyp_target);
		return(-1);
	}
	fprintf(stdout, "NLL target hypocenter-phase file: %s\n", fn_hyp_target);
	// open "phs_source" NLL hypocenter file
	istat = GetHypLoc(NULL, fn_hyp_phs_source, &Hypo_phs_source, Arrival_phs_source, &NumArrivals_phs_source, 1, &Grid_phs_source, 0);
	if (istat < 0) {
		nll_puterr2("ERROR: opening phs_source NLL hypocenter-phase file.", fn_hyp_target);
		return(-1);
	}
	fprintf(stdout, "NLL phs_source hypocenter-phase file: %s\n", fn_hyp_target);
	// open copy of "target" NLL hypocenter file for output
	istat = GetHypLoc(NULL, fn_hyp_target, &Hypo_out, Arrival_out, &NumArrivals_out, 1, &Grid_out, 0);

	// set projection
printf("MapProjStr[0]=%s\n", MapProjStr[0]);
	if (strstr(MapProjStr[0], "GLOBAL") != NULL)
		GeometryMode = MODE_GLOBAL;
printf("GeometryMode=%d  MODE_GLOBAL=%d\n", GeometryMode, MODE_GLOBAL);

	// initialize selected array
	for (narr1 = 0; narr1 < NumArrivals_phs_source; narr1++)
		isSelected[narr1] = 0;


	// loop over target arrivals
	for (narr0 = 0; narr0 < NumArrivals_target; narr0++) {

		if (Arrival_target[narr0].weight < SMALL_DOUBLE)
			continue;
		nAttempt++;
		target_phase = Arrival_target[narr0].phase;
		xval = Arrival_target[narr0].station.x;
		yval = Arrival_target[narr0].station.y;
		dist_min = VERY_LARGE_DOUBLE;
		index_min = -1;

		// loop over phase source arrivals, find closest
		for (narr1 = 0; narr1 < NumArrivals_phs_source; narr1++) {
//printf("%d/%d %s/%s", narr0, narr1, Arrival_phs_source[narr0].phase, Arrival_phs_source[narr1].phase);
			if (isSelected[narr1])
				continue;
			// check for match of phase
			source_phase = Arrival_phs_source[narr1].phase;
			if (IsEquivPhase(target_phase, source_phase)) {
				// compare distances
				dist = GetEpiDist(&(Arrival_phs_source[narr1].station), xval, yval);
				if (dist < dist_min) {
					dist_min = dist;
					index_min = narr1;
				}
//printf(" d=%lf dmin=%lf imin=%d\n", dist, dist_min, index_min);
			}
		}
		dist_epi = (Arrival_phs_source[index_min].dist + Arrival_target[narr0].dist) / 2.0;
		if (GeometryMode == MODE_GLOBAL)
			dist_epi /= KM2DEG;
//printf("dist_epi %f (%f) source %f  target %f\n", dist_epi, dist_max_fract * dist_epi,
//Arrival_phs_source[index_min].dist, Arrival_target[narr0].dist);

		if (index_min >= 0 && dist_min <= dist_max_fract * dist_epi) {
//printf("MATCHED: dist_min %f <= dist_max_fract * dist_epi %f\n", dist_min, dist_max_fract * dist_epi);
			Arrival_out[nPhaseMatched] = Arrival_phs_source[index_min];
			isSelected[index_min] = 1;
			nPhaseMatched++;
		}

	}


	// write output location
	if (nPhaseMatched) {
		// open output NLL hypocenter file
		if ((fp_hyp_out = fopen(fn_hyp_out, "w")) == NULL) {
			nll_puterr2("ERROR: opening output NLL phase file", fn_hyp_out);
			return(-1);
		}
		WritePhases(fp_hyp_out, &Hypo_out,
			Arrival_out, nPhaseMatched, NULL, 1, 0, 1, &Grid_out, 0, IO_ARRIVAL_ALL);
		fclose(fp_hyp_out);
	}


	// write message
	fprintf(stdout,
		"\n%d target phases, %d attempt, %d source phases, %d matched, %d written to %s.\n",
			NumArrivals_target, nAttempt, NumArrivals_phs_source,
			nPhaseMatched, nPhaseMatched, fn_hyp_out);

	return(nPhaseMatched);

}


int IsEquivPhase(char *target_phase, char *source_phase) {

//printf("source %s  target %s  match %d\n", target_phase, source_phase,
//lastLegType(target_phase) == lastLegType(source_phase));

	if (lastLegType(target_phase) == lastLegType(source_phase))
		return(1);

	return(0);


}


/*** function to determine type P or S of last leg of phase */

// returns 'P' if P, 'S' if S, ' ' otherwise

char lastLegType(char *phase)
{
	char *c_last_p, *c_last_P, *c_last_s, *c_last_S;
	int i_last_p, i_last_P, i_last_s, i_last_S;

	// get offsets to last char of types p/P or s/S
	c_last_p = strrchr(phase, 'p');
	i_last_p = c_last_p == NULL ? -1 : (int) (c_last_p - phase);
	c_last_P = strrchr(phase, 'P');
	i_last_P = c_last_P == NULL ? -1 : (int) (c_last_P - phase);
	c_last_s = strrchr(phase, 's');
	i_last_s = c_last_s == NULL ? -1 : (int) (c_last_s - phase);
	c_last_S = strrchr(phase, 'S');
	i_last_S = c_last_S == NULL ? -1 : (int) (c_last_S - phase);

	// determine type of last char
	i_last_p = i_last_p > i_last_P ? i_last_p : i_last_P;
	i_last_s = i_last_s > i_last_S ? i_last_s : i_last_S;
	if (i_last_p >= 0 && i_last_p > i_last_s)
		return('P');
	if (i_last_s >= 0 && i_last_s > i_last_p)
		return('S');
	return(' ');

}


