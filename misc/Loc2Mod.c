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


/*   Loc2Mod.c

	Program to find optimal velcoity model based on location residual statistics 

*/

//------------------------------------------------------------/
// Anthony Lomax           | email: lomax@faille.unice.fr     /
// UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax /
// 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        /
// 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        /
//------------------------------------------------------------/


/*
	history:

	ver 01    05NOV1999  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

*/



#include "../src/GridLib.h"
#include "ran1.h"

// defines



// globals

char fninput[FILENAME_MAX];
#define MAX_NUM_INPUT_FILES 500
char fnroot_input_list[MAX_NUM_INPUT_FILES][FILENAME_MAX];
char fnroot_input[FILENAME_MAX];


// functions



/** 	Program to find optimal velcoity model based on location residual statistics 
 */

#define PNAME  "Loc2Mod"


main(int argc, char *argv[])
{

	int istat, narg;


	// set program name

	strcpy(prog_name, PNAME);


	// check command line for correct usage

	fprintf(stdout, "\n%s Arguments: ", prog_name);
	for (narg = 0; narg < argc; narg++)
		fprintf(stdout, "<%s> ", argv[narg]);
	fprintf(stdout, "\n");

	if (argc != 1)
	{
		nll_puterr("ERROR wrong number of command line arguments.");
		disp_usage(PNAME, "< >");
		exit(-1);
	}


	// read command line arguments

	//strcpy(fninput, argv[1]);


	// find optimal velocity model

	doSearchGradModel();

	exit(0);

}


#define RANDOM_SEED 4321

#define BIN_PATH "/u/thermie/export/home/lomax/bin/"

#define FILENAME_CONTROL "loc2mod.in"
#define FILENAME_VEL_INCLUDE "loc2model_vel_include.in"
#define FILENAME_LOC_INCLUDE "loc2model_loc_include.in"
#define FILENAME_LOC_STAT "/teraves/lomax/work/Loc2Mod/loc/loc2mod.sum.grid0.loc.stat"
#define FILENAME_RESULTS "/teraves/lomax/work/Loc2Mod/loc2model.out"

#define VEL0_MIN 2.0
#define VEL0_MAX 4.0

#define VEL_GRAD_MIN 0.0
#define VEL_GRAD_MAX 1.0

#define TOP_OF_LAYER -2.0
#define VEL_MAX_DEPTH 2.0
#define VEL_MAX_AT_DEPTH 6.0

#define VP_VS_MIN 1.6
#define VP_VS_MAX 2.0001
#define VP_VS_STEP 0.05


/** function to find optimal gradient model */

int doSearchGradModel(char* fnroot)
{

	int istat;
	char command_string[MAXLINE];
	char file_string[MAXLINE];
	FILE *fp_stat, *fp_include, *fp_out;

	double vel0, vel_grad, vp_vs;

	char *cstat;
	char line[MAXLINE_LONG];
	int i_found_delays;
	char keyword[20], label[12], phase[12];
	int num_residuals;
	double residual, residual_std, residual_min, residual_max;

	double ave_res_sum, std_dev_sum;
	int ave_res_weight_sum, std_dev_weight_sum;
	double ave_res_sum_min = VERY_LARGE_DOUBLE, std_dev_sum_min = VERY_LARGE_DOUBLE;

	int i_minimum;
	int ntest = 0;


	/* initialize random number generator */
	SRAND_FUNC(RANDOM_SEED);


	remove(FILENAME_RESULTS);


	while (1) {

		// Monte Carlo search over models

		vel0 = get_rand_double(VEL0_MIN, VEL0_MAX);
		vel_grad = get_rand_double(VEL_GRAD_MIN, VEL_GRAD_MAX);

		// check that velocity at depth is reasonable
		if (vel0 + (VEL_MAX_DEPTH - TOP_OF_LAYER) * vel_grad > VEL_MAX_AT_DEPTH) {
printf("REJECT: vel0 %lf  vel_grad %lf  vel_at_depth %lf\n", vel0, vel_grad, vel0 + (VEL_MAX_DEPTH - TOP_OF_LAYER) * vel_grad);
			continue;
		}

		fp_include = fopen(FILENAME_VEL_INCLUDE, "w");
		// gradient layer
		sprintf(file_string, "LAYER  %lf %lf %lf  9.99 0.00  9.9 0.0\n", 
			TOP_OF_LAYER, vel0, vel_grad);
		fprintf(fp_include, "%s", file_string);
		printf(file_string);
		// constant velocity basement
		sprintf(file_string, "LAYER  %lf %lf %lf  9.99 0.00  9.9 0.0\n", 
			VEL_MAX_DEPTH, VEL_MAX_AT_DEPTH, 0.0);
		fprintf(fp_include, "%s", file_string);
		printf(file_string);
		fclose(fp_include);

		// run Vel2Grid
		sprintf(command_string, "%sVel2Grid %s\n", BIN_PATH, FILENAME_CONTROL);
		//printf(command_string);
		system(command_string);
		// run Grid2Time
		sprintf(command_string, "%sGrid2Time %s\n", BIN_PATH, FILENAME_CONTROL);
		//printf(command_string);
		system(command_string);


		// grid search over Vp/Vs

		for (vp_vs = VP_VS_MIN; vp_vs <= VP_VS_MAX; vp_vs += VP_VS_STEP) {

			fp_include = fopen(FILENAME_LOC_INCLUDE, "w");
/* method (GAU_ANALYTIC), maximum_dist_sta_to_grid, minimum_number_phases, maximum_number_phases, minimum_number_S_phases, Vp/Vs, maximum_number_3D_grids */
			sprintf(file_string, "LOCMETH GAU_ANALYTIC 9999.0 11 -1 -1 %lf -1\n", vp_vs);
			fprintf(fp_include, "%s", file_string);
			//printf(file_string);
			fclose(fp_include);

			// run NLLoc
			sprintf(command_string, "%sNLLoc %s\n", BIN_PATH, FILENAME_CONTROL);
			//printf(command_string);
			system(command_string);

			ntest++;

			// read stat file
			fp_stat = fopen(FILENAME_LOC_STAT, "r");

			i_found_delays = 0;
			std_dev_weight_sum = 0;
			while((cstat = fgets(line, MAXLINE_LONG, fp_stat)) != NULL) {

				istat = sscanf(line, "%s %s %s %d %lf %lf %lf %lf", 
						keyword,
						label, phase, &num_residuals, &residual, &residual_std, 
						&residual_min, &residual_max);

				if (istat != 8 || strcmp(keyword, "LOCDELAY") != 0) {

					// have read all delays
					if (i_found_delays)
						break;

					// have not started reading delays
					continue;
				}


				// have not started reading delays
				if (!i_found_delays) {

					i_found_delays = 1;

					ave_res_sum = (double) num_residuals * fabs(residual);
					ave_res_weight_sum = num_residuals;
					std_dev_sum = (double) (num_residuals - 1) * residual_std;
					std_dev_weight_sum = (num_residuals - 1);

				// are reading delays
				} else  {

					ave_res_sum += (double) num_residuals * fabs(residual);
					ave_res_weight_sum += num_residuals;
					std_dev_sum += (double) (num_residuals - 1) * residual_std;
					std_dev_weight_sum += (num_residuals - 1);

				}

			}

			fclose(fp_stat);

			if (!i_found_delays || std_dev_weight_sum < 1)
				continue;

			ave_res_sum /= (double) ave_res_weight_sum;
			std_dev_sum /= (double) std_dev_weight_sum;

			// check for minimum residuals
			i_minimum = 0;
			if (ave_res_sum <= ave_res_sum_min) {
				i_minimum = 1;
				ave_res_sum_min = ave_res_sum;
			}
			if (std_dev_sum <= std_dev_sum_min) {
				i_minimum = 1;
				std_dev_sum_min = std_dev_sum;
			}
			if (i_minimum) {
				fp_out = fopen(FILENAME_RESULTS, "a+");
				fseek(fp_out, 0L, SEEK_END);
				fprintf(fp_out, 
"ntest %d  ave_res %lf (%d)  ave_std_dev %lf (%d)  vel0 %lf  vel_grad %lf  vp_vs  %lf\n",
					ntest, ave_res_sum, ave_res_weight_sum, std_dev_sum, std_dev_weight_sum, 
					vel0, vel_grad, vp_vs);
				fclose(fp_out);

			}

		}

	}

	return(0);

}






//------------------------------------------------------------/
// Anthony Lomax           | email: lomax@faille.unice.fr     /
// UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax /
// 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        /
// 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        /
//------------------------------------------------------------/

