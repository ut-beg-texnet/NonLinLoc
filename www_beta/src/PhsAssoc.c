/*
 * Copyright (C) 2003 Anthony Lomax <www.alomax.net, anthony@alomax.net>
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


/*   PhsAssoc.c

	Program to associate phase picks with phase codes using the NLLoc maximum likelihood hypocenter
	and travel time table for each possible phase code.

*/


/*
	history:

	ver 01    12Dec2003  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

*/



#include "GridLib.h"


/* defines */

#define MAX_NUM_INPUT_FILES 4096

#define VERY_LARGE_RESIDUAL 999999.9
#define TIME_TOLERANCE 0.001


/* globals */



/* functions */

int AssociatePhases(int argc, char** argv);
int CalcObsTravelTimes(ArrivalDesc *arrival, int num_arrivals, HypoDesc* phypo);


/*** program to sum event scatter files */

#define PNAME  "PhsAssoc"


int main(int argc, char** argv)
{

	int narg;
	int nPhaseIDChanged;


	/* set program name */

	strcpy(prog_name, PNAME);


	/* check command line for correct usage */

	fprintf(stdout, "%s Arguments: ", prog_name);
	for (narg = 0; narg < argc; narg++)
		fprintf(stdout, "<%s> ", argv[narg]);
	fprintf(stdout, "\n");

	if (argc < 7) {
		nll_puterr("ERROR wrong number of command line arguments.");
		disp_usage(PNAME,
"max_residual i_swap_bytes i_look_for_sta_grids i_write_all <time_grid_root> <nll_hyp_files...>");
		exit(0);
	}

	SetConstants();
	message_flag = 1;
	DispProgInfo();
	message_flag = 0;

	if ((nPhaseIDChanged = AssociatePhases(argc, argv)) < 0) {
		nll_puterr("ERROR doing Phase Association process.");
		exit(0);
	}



	exit(nPhaseIDChanged);

}



int AssociatePhases(int argc, char** argv)
{

	int narg;

	int iTryStationGrids, iWriteAll;
	double residual_max;

	int istat;
	int nArr, nTimeGrid;
	int nPhaseIDChanged, nPhaseIDChanged_total, nZeroWeight;
	int nEventsRead, nEventsChanged, nEventsOutside;
	int numFiles1, numFiles2;

	char fn_hyp[FILENAME_MAX];
	char fn_loc_grids[FILENAME_MAX], fn_loc_grids1[FILENAME_MAX], fn_loc_grids2[FILENAME_MAX];
	char fn_hyp_out[FILENAME_MAX];

	char fn_loc_grids_list[MAX_NUM_INPUT_FILES][FILENAME_MAX];

	char *pchr0, *pchr1;
	char tmp_str[FILENAME_MAX];
	char best_phase_id[PHASE_LABEL_LEN];
	char orig_phase_id[PHASE_LABEL_LEN];

	FILE *fp_hyp_out;
	FILE *fp_grid, *fp_hdr;

	double obs_travel_time, pred_travel_time, best_pred_travel_time = -1.0;
	double residual, orig_residual, best_residual;
	double yval_grid;

	GridDesc Grid, locgrid;
	SourceDesc Srce;
	HypoDesc Hypo;


	// get command line parameters

	narg = 1;

	residual_max = VERY_LARGE_RESIDUAL;
	sscanf(argv[narg], "%lf", &residual_max);
	narg++;

	Grid.iSwapBytes = 0;
	sscanf(argv[narg], "%d", &(Grid.iSwapBytes));
	narg++;

	iTryStationGrids = 1;
	sscanf(argv[narg], "%d", &iTryStationGrids);
	narg++;

	iWriteAll = 1;
	sscanf(argv[narg], "%d", &iWriteAll);
	narg++;

	strcpy(fn_loc_grids, argv[narg]);
	narg++;


	// associate each NLL hypo file

	nPhaseIDChanged_total = 0;
	nEventsRead = 0;
	nEventsChanged = 0;
	nEventsOutside = 0;

	for (; narg < argc; narg++) {

		strcpy(fn_hyp, argv[narg]);

		// open input NLL hypocenter file
		pchr0 = strstr(fn_hyp, ".hyp");
		if (pchr0 != NULL)
			*pchr0 = '\0';
		istat = GetHypLoc(NULL, fn_hyp, &Hypo, Arrival, &NumArrivals, 1, &locgrid, 0);

		if (istat < 0) {
			nll_puterr2("ERROR: opening NLL hypocenter-phase file.", fn_hyp);
			continue;
		}
		fprintf(stdout,
			"\nNLL hypocenter-phase file: %s\n", fn_hyp);

		nEventsRead++;

		if (!IsPointInsideGrid(&locgrid, Hypo.x, Hypo.y, Hypo.z)) {
			fprintf(stdout,
				"Outside grid, ignoring event: %s\n", fn_hyp);
			nEventsOutside++;
			continue;
		}

		if (isOnGridBoundary(Hypo.x, Hypo.y, Hypo.z, &locgrid, locgrid.dx, locgrid.dz, 0)) {
			fprintf(stdout,
				"WARNING: Event on grid boundary: %s\n", fn_hyp);
		}

		if(CalcObsTravelTimes(Arrival, NumArrivals, &Hypo) < 0)
			continue;

		// check transform type
		if (strstr(MapProjStr[0], "GLOBAL") != NULL)
			GeometryMode = MODE_GLOBAL;


		nPhaseIDChanged = 0;
		nZeroWeight = 0;
		for (nArr = 0; nArr < NumArrivals; nArr++) {

			//WriteArrival(stdout,  Arrival + nArr, IO_ARRIVAL_ALL);

			// check for fixed phase type (inst begins with '+')
			// do not change arrival
			if (Arrival[nArr].inst[0] == '+')
				continue;

			// check for no absolute timing (inst begins with '*')
			// do not change arrival
			if (Arrival[nArr].inst[0] == '*')
				continue;

			obs_travel_time = Arrival[nArr].obs_travel_time;
			orig_residual = Arrival[nArr].residual;
			best_residual = VERY_LARGE_RESIDUAL;
			strcpy(orig_phase_id, Arrival[nArr].phase);
			strcpy(best_phase_id, Arrival[nArr].phase);
			// remove zero weight phase tag
			if (best_phase_id[0] == '*') {
				strcpy(best_phase_id, best_phase_id + 1);
				best_residual = VERY_LARGE_RESIDUAL;
			}

			// get list of time grid files for this station
			numFiles1 = 0;
			if (iTryStationGrids) {
				sprintf(fn_loc_grids1, "%s.*.%s.time.buf", fn_loc_grids, Arrival[nArr].label);
				if ((numFiles1 = ExpandWildCards(fn_loc_grids1,
						fn_loc_grids_list, MAX_NUM_INPUT_FILES)) < 0) {
					nll_puterr("ERROR: getting list of station specifc time grid files.");
						numFiles1 = 0;
				}
			}
			if (numFiles1 == 0) {
				// if no station grid, get list of default time grid files
				sprintf(fn_loc_grids2, "%s.*.DEFAULT.time.buf", fn_loc_grids);
				if ((numFiles2 = ExpandWildCards(fn_loc_grids2,
						fn_loc_grids_list + numFiles1, MAX_NUM_INPUT_FILES - numFiles1)) < 0) {
					nll_puterr("ERROR: getting list of DEFAULT time grid files.");
						numFiles2 = 0;
				}
			} else {
				numFiles2 = 0;
			}
			if (numFiles1 + numFiles2 >= MAX_NUM_INPUT_FILES) {
				sprintf(MsgStr, "WARNING: maximum number of travel-time files exceeded, only first %d will be processed.", MAX_NUM_INPUT_FILES);
				nll_puterr(MsgStr);
			}

			for (nTimeGrid = 0; nTimeGrid < numFiles1 + numFiles2; nTimeGrid++) {
				// open next time grid file
				pchr0 = strstr(fn_loc_grids_list[nTimeGrid], ".buf");
				if (pchr0 != NULL)
					*pchr0 = '\0';
				if ((istat = OpenGrid3dFile(fn_loc_grids_list[nTimeGrid], &fp_grid, &fp_hdr,
							&Grid, "time", &Srce, Grid.iSwapBytes)) < 0) {
						sprintf(MsgStr, "%s.time.*", fn_loc_grids_list[nTimeGrid]);
					nll_puterr2("ERROR: opening grid file", MsgStr);
					continue;
				}
				if (Arrival[nArr].gdesc.type == GRID_TIME) {
					pred_travel_time =
						(double) ReadAbsInterpGrid3d(
							fp_grid, &Grid, Hypo.x, Hypo.y, Hypo.z);
				} else {
					yval_grid = GetEpiDist(&(Arrival[nArr].station), Hypo.x, Hypo.y);
					if (GeometryMode == MODE_GLOBAL)
						yval_grid *= KM2DEG;
					pred_travel_time =
						ReadAbsInterpGrid2d(fp_grid, &Grid, yval_grid, Hypo.z);
				}
				CloseGrid3dFile(&fp_grid, &fp_hdr);
				if (pred_travel_time < 0.0)
					continue;
				// check if pred time is closest to obs time
				// iasp91.ScP.DEFAULT.time
				residual = obs_travel_time - pred_travel_time;
				if (fabs(residual) < fabs(best_residual) - TIME_TOLERANCE) {
					best_residual = residual;
					best_pred_travel_time = pred_travel_time;
					pchr0 = fn_loc_grids_list[nTimeGrid];
					pchr1 = strrchr(pchr0, '.');
					strncpy(tmp_str, pchr0, pchr1 - pchr0);
					pchr1 = strrchr(tmp_str, '.');
					*pchr1 = '\0';
					pchr1 = strrchr(tmp_str, '.');
					strcpy(best_phase_id, pchr1 + 1);
				}
/*if (strcmp(Arrival[nArr].label, "Berkel") == 0) printf("%d %s (%fdeg) %s %s (%fs) -> %s  (%fs)  current: %s  (%fs)\n",
	nTimeGrid, Arrival[nArr].label, yval_grid, Arrival[nArr].comp, orig_phase_id,
	orig_residual, best_phase_id, best_residual, fn_loc_grids_list[nTimeGrid], residual);
*/
			}
/*					fprintf(stdout, "!!! %s %s (%fs) -> %s  (%fs)\n",
						Arrival[nArr].label, orig_phase_id,
						orig_residual, best_phase_id, best_residual);
*/
   // check if best phase is different from original phase
			if (strcmp(Arrival[nArr].phase, best_phase_id) != 0) {
				// update arrival
				strcpy(Arrival[nArr].phase, best_phase_id);
				Arrival[nArr].pred_travel_time = best_pred_travel_time;
				Arrival[nArr].residual = best_residual;
			}
			// check residual
			if (fabs(Arrival[nArr].residual) > residual_max) {
				// update arrival
				sprintf(tmp_str, "*%s", Arrival[nArr].phase);
				strcpy(Arrival[nArr].phase, tmp_str);
			}
			// check if changed
			if (strcmp(orig_phase_id, Arrival[nArr].phase) != 0) {
				// phase ID changed
				nPhaseIDChanged++;
				fprintf(stdout, "%s %s %s (%fs) -> %s  (%fs)",
					Arrival[nArr].label, Arrival[nArr].comp, orig_phase_id,
					orig_residual, best_phase_id, best_residual);
				if (Arrival[nArr].phase[0] == '*') {
					// -> zero weight
					nZeroWeight++;
					fprintf(stdout, " -> %s (zero weight)", Arrival[nArr].phase);
				}
				fprintf(stdout, "\n");

			}


		}


		if (nPhaseIDChanged > 0)
			nEventsChanged++;

		// write output location
		if (iWriteAll || nPhaseIDChanged > 0) {
			// open output NLL hypocenter file
			strcpy(fn_hyp_out, fn_hyp);
			pchr0 = strstr(fn_hyp_out, ".hyp");
			if (pchr0 != NULL)
				*pchr0 = '\0';
			strcat(fn_hyp_out, ".phs_assoc");
			if ((fp_hyp_out = fopen(fn_hyp_out, "w")) == NULL) {
				nll_puterr2("ERROR: opening output NLL associated phase file", fn_hyp_out);
				continue;
			}
			WritePhases(fp_hyp_out, &Hypo,
				Arrival, NumArrivals, NULL, 1, 0, 1, &locgrid, 0, IO_ARRIVAL_OBS);
			fclose(fp_hyp_out);
		}

		// write message
		fprintf(stdout,
			"%d phases read, %d IDs changed, %d given 0 weight.\n",
				NumArrivals, nPhaseIDChanged, nZeroWeight);

		nPhaseIDChanged_total += nPhaseIDChanged;

	}

	// write message
	fprintf(stdout,
		"\n%d events read, %d outside grid, %d changed, %d phases changed.\n",
			nEventsRead, nEventsOutside, nEventsChanged, nPhaseIDChanged_total);

	return(nPhaseIDChanged_total);

}



/*** function to complete and verify obs_travel_time of arrivals */

int CalcObsTravelTimes(ArrivalDesc *arrival, int num_arrivals, HypoDesc* phypo)
{
	int narr;
	int dofymin = 10000, yearmin = 10000;

	double obs_time, obs_travel_time, residual;



	// calc and check each arrival

	// get day of year and min doy, year
	for (narr = 0; narr < num_arrivals; narr++) {
		arrival[narr].day_of_year =
			DayOfYear(arrival[narr].year, arrival[narr].month, arrival[narr].day);
		if (arrival[narr].day_of_year < dofymin)
			dofymin = arrival[narr].day_of_year;
		if (arrival[narr].year < yearmin)
			yearmin = arrival[narr].year;
		if (arrival[narr].year != yearmin) {
printf("ERROR: arrivals cross year boundary, ignoring observation set.");
			return(-1);
		}
	}

	// calc hypocenter otime in sec from beginning of day of year;
	if (phypo->year != yearmin || DayOfYear(phypo->year, phypo->month, phypo->day) != dofymin) {
		nll_puterr(
"ERROR: earliest arrivals year/month/day does not match fixed origin time year/month/day, ignoring observation set.");
		return(OBS_FILE_ARRIVALS_CROSS_YEAR_BOUNDARY);
	}
	phypo->time = (long double) phypo->sec
		+ 60.0L * ((long double) phypo->min
		+ 60.0L * (long double) phypo->hour);


	// process arrivals

	for (narr = 0; narr < num_arrivals; narr++) {

		// correct to min day of year
		if (arrival[narr].day_of_year > dofymin) {
			arrival[narr].day_of_year--;
			arrival[narr].day--;
			arrival[narr].hour += 24;
		}

		// calc  observed time in sec from beginning of day of year;
		obs_time = arrival[narr].sec
		// DELAY_CORR
		//	- (long double) arrival[narr].delay
			+ 60.0 * ((double) arrival[narr].min
				+ 60.0 * (double) arrival[narr].hour);
/*printf("%s %s obs_time (%f) (%d %d %f)\n",
arrival[narr].label, arrival[narr].phase, obs_time, arrival[narr].hour,
arrival[narr].min, (float) arrival[narr].sec);*/
			arrival[narr].obs_time = obs_time;

		// calc  observed travel time obs_time - hypo.time from beginning  of day of year;
		obs_travel_time = arrival[narr].obs_time - (double) phypo->time;
		arrival[narr].obs_travel_time =
			arrival[narr].pred_travel_time + arrival[narr].residual + arrival[narr].delay;
		if (fabs(arrival[narr].obs_travel_time) > TIME_TOLERANCE &&
			fabs(arrival[narr].obs_travel_time - obs_travel_time) > TIME_TOLERANCE) {
printf("WARNING: %s %s arrival.obs_travel_time (%f) != obs_travel_time (%f)\n",
arrival[narr].label, arrival[narr].phase, arrival[narr].obs_travel_time, obs_travel_time);
			arrival[narr].obs_travel_time = obs_travel_time;
		}

		// calc residual;
		if (arrival[narr].pred_travel_time < TIME_TOLERANCE) {
			arrival[narr].residual = VERY_LARGE_RESIDUAL;
		} else {
			residual = arrival[narr].obs_travel_time
				- arrival[narr].pred_travel_time - arrival[narr].delay;
			if (fabs(arrival[narr].residual) > TIME_TOLERANCE &&
				fabs(arrival[narr].residual - residual) > TIME_TOLERANCE) {
printf("WARNING: %s %s arrival.residual (%f) != residual (%f)\n",
arrival[narr].label, arrival[narr].phase, arrival[narr].residual, residual);
				arrival[narr].residual = residual;
			}
		}

	}


	return(0);


}





