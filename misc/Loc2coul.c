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


/*   Loc2coul.c

	Program to sum location files into a Coulomb grid.

*/



/*
	history:

	ver 01    25Oct2004  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

*/



#include "../src/GridLib.h"
#include "ran1.h"
#include "../okada/okada.h"


/* defines */

#define MAX_NUM_INPUT_FILES 4096


/* globals */

double mech_dipDir, mech_dipAng, mech_rake;


/* functions */

int ConvertToCoulomb(int argc, char** argv);
double calc_rad(double ray_vert_ang, double ray_horiz_ang, char orig_wave);
double faultShearStress(const double *slip, const double *fault_normal, double *stress);
void PointStressResponse(double c, double dip,
		      const double *slip, const double *obspoint,
		      double *stress);
double pointFaultShearStressStress(double c, double dip,
		      const double *slip, const double *obspoint, double *stress);
int ReadCoul_Input(FILE* fp_input);
int saveSamples(int n_samples, GridDesc* ptgrid, char* fileroot, double cutoff, char* code);



/*** program to sum event scatter files */

#define PNAME  "Loc2coul"


int main(int argc, char** argv)
{

	int istat, narg;


	/* set program name */

	strcpy(prog_name, PNAME);


	/* check command line for correct usage */

	fprintf(stdout, "%s Arguments: ", prog_name);
	for (narg = 0; narg < argc; narg++)
		fprintf(stdout, "<%s> ", argv[narg]);
	fprintf(stdout, "\n");

	if (argc < 9) {
		nll_puterr("ERROR wrong number of command line arguments.");
		disp_usage(PNAME, "<control_file> <add_file_list> <output_file_root> <output_grid_hdr> nobs_min misfit_max n_smooth(+y,-x) n_samples");
		exit(-1);
	}


	SetConstants();
	message_flag = 1;
	DispProgInfo();
	message_flag = 3;
	RandomNumSeed = 1234;

	/* initialize random number generator */
	SRAND_FUNC(RandomNumSeed);
	if (message_flag >= 4)
		test_rand_int();

	if ((istat = ConvertToCoulomb(argc, argv)) < 0) {
		nll_puterr("ERROR converting to Coulomb grid.");
		exit(-1);
	}



	exit(0);

}




int ConvertToCoulomb(int argc, char** argv)
{

	int istat, i;
	int nArrivals;

	int n_obs_min;
	double misfit_max;
	int n_smooth, n_samples;

	char fn_hyp_in[FILENAME_MAX], fn_root_out[FILENAME_MAX];
	char grid_hdr_filename[FILENAME_MAX], fngrid_control[FILENAME_MAX];
	FILE *fp_hyp_in, *fp_grid_control;

	int nHypo, numFiles, nLocWritten, nLocAccepted;
	char fn_hyp_in_list[MAX_NUM_INPUT_FILES][FILENAME_MAX_SMALL];
	char test_str[10];

	GridDesc loc_grid, coul_grid, wt_grid, smoothed_grid;
	HypoDesc Hypo;

	ArrivalDesc Arrival[1];

	int ix, iy, iz;
	double dx, dy, dz;
	double xval, yval, xloc, yloc, zdepth;
	double mech_strikeDir;
	double azimuth, ray_dip, dist_xyz, dist_horiz;
	double coulomb, weight;
	double value_max_abs;

	int ny, nx;
	double smooth_sum;

	// okada
	double slip[3], obspoint[3], stress[9];
	double okada_az;




	strcpy(test_str, ".hyp");

	/* get command line parameters */
	strcpy(fn_hyp_in, argv[2]);
	strcpy(fn_root_out, argv[3]);
	strcpy(grid_hdr_filename, argv[4]);
	sscanf(argv[5], "%d", &n_obs_min);
	sscanf(argv[6], "%lf", &misfit_max);
	sscanf(argv[7], "%d", &n_smooth);
	sscanf(argv[8], "%d", &n_samples);

	/* read control file */

	strcpy(fngrid_control, argv[1]);
	if ((fp_grid_control = fopen(fngrid_control, "r")) == NULL) {
		nll_puterr("ERROR opening control file.");
		exit(EXIT_ERROR_FILEIO);
	}
	if ((istat = ReadCoul_Input(fp_grid_control)) < 0) {
		exit(EXIT_ERROR_FILEIO);
	}
	fclose(fp_grid_control);


	/* initialize 3D grid */
	ReadGrid3dHdr(&coul_grid, NULL, grid_hdr_filename, "coul");
	/* allocate grid */
	coul_grid.buffer = AllocateGrid(&coul_grid);
	if (coul_grid.buffer == NULL) {
		nll_puterr("ERROR: allocating memory for 3D Coulomb grid buffer.");
		exit(EXIT_ERROR_MEMORY);
	}
	/* create array access pointers */
	coul_grid.array = CreateGridArray(&coul_grid);
	if (coul_grid.array == NULL) {
		nll_puterr("ERROR: creating array for accessing 3D Coulomb grid buffer.");
		exit(EXIT_ERROR_MEMORY);
	}
	// create weights grid
	DuplicateGrid(&wt_grid , &coul_grid , "GRID_UNDEF");

	// initialize output grid
	for (ix = 0; ix <  coul_grid.numx; ix++)
		for (iy = 0; iy <  coul_grid.numy; iy++)
			for (iz = 0; iz <  coul_grid.numz; iz++) {
				coul_grid.array[ix][iy][iz] = 0.0f;
				wt_grid.array[ix][iy][iz] = 0.0f;
			}


	/* sum requested loc files into output grid */

	/* check for wildcards in input file name */
	strcat(fn_hyp_in, test_str);
	if ((numFiles = ExpandWildCards(fn_hyp_in,
			fn_hyp_in_list, MAX_NUM_INPUT_FILES)) < 1) {
		nll_puterr("ERROR: no matching .hyp files found.");
		return(-1);
	}
	if (numFiles >= MAX_NUM_INPUT_FILES) {
		sprintf(MsgStr, "WARNING: maximum number of event files exceeded, only first %d will be processed.", MAX_NUM_INPUT_FILES);
		nll_puterr(MsgStr);
	}


	nLocWritten = 0;
	nLocAccepted = 0;
	for (nHypo = 0; nHypo < numFiles; nHypo++) {

		/* open hypocenter file */
		sprintf(strstr(fn_hyp_in_list[nHypo], test_str), "\0");
		sprintf(fn_hyp_in, "%s.hyp", fn_hyp_in_list[nHypo]);
		if ((fp_hyp_in = fopen(fn_hyp_in, "r")) == NULL) {
			nll_puterr2("ERROR: opening hypocenter file, ignoring event, file",
				fn_hyp_in_list[nHypo]);
			continue;
		}

		while (1) {

			istat = GetHypLoc(fp_hyp_in, fn_hyp_in_list[nHypo], &Hypo, Arrival,
					&nArrivals, 0, &loc_grid, 0);
			if (istat == EOF) {
				fclose(fp_hyp_in);
				break;
			}
			if (strcmp(Hypo.locStat, "ABORTED") == 0) {
				//nll_puterr("WARNING: location ABORTED, ignoring event");
				continue;
			} else if (strcmp(Hypo.locStat, "REJECTED") == 0) {
				//nll_puterr("WARNING: location REJECTED, ignoring event");
				continue;
			}
			nLocAccepted++;


			// do coulomb calculation
			if (Hypo.focMech.nObs >= n_obs_min && Hypo.focMech.misfit <= misfit_max) {
fprintf(stdout,
"%4.4d%2.2d%2.2d %2.2d%2.2d%2.2d  Mech: dipDir %f  dipAng %f  rake %f  nObs %d  misfit %f\n",
Hypo.year, Hypo.month, Hypo.day, Hypo.hour, Hypo.min, (int) Hypo.sec,
Hypo.focMech.dipDir, Hypo.focMech.dipAng, Hypo.focMech.rake, Hypo.focMech.nObs, Hypo.focMech.misfit);

// dipDir;	Dip direction (downdip azimuth in degrees,clockwise from north)
// dipAng;	Dip angle in degrees down from horizontal
// rake;	Rake in degrees: 0=left lateral, 0 to 180=thrust, +-180=right lateral, -180 to 0=normal */

				mech_dipDir = Hypo.focMech.dipDir;
				// adjust azimuth to A+R convention
				mech_strikeDir = Hypo.focMech.dipDir - 90.0;
				// adjust azimuth for projection azimuth
				mech_dipDir = latlon2rectAngle(0, mech_dipDir);
				mech_strikeDir = latlon2rectAngle(0, mech_strikeDir);

				mech_dipDir *= cRPD;
				mech_strikeDir *= cRPD;    // radians CW from Y axis
				mech_dipAng = Hypo.focMech.dipAng * cRPD;
				mech_rake = Hypo.focMech.rake * cRPD;
//fprintf(stdout, "mech_strikeDir %f  mech_dipAng %f  mech_rake %f\n", mech_strikeDir / cRPD, mech_dipAng / cRPD, mech_rake / cRPD);

				// okada
				slip[0] = cos(mech_rake);	// strike-slip
				slip[1] = sin(mech_rake);	// dip-slip

				slip[2] = 0.0;	// tensile

				// ajdust for right/left lateral
				//if (slip[0] < 0) {
				//	slip[0] *= -1.0;
				//	slip[1] *= -1.0;
				//}


				// loop over cells

				xval = coul_grid.origx + coul_grid.dx / 2.0;
				for (ix = 0; ix <  coul_grid.numx; ix++) {

					yval = coul_grid.origy + coul_grid.dy / 2.0;
					for (iy = 0; iy <  coul_grid.numy; iy++) {

						//if (ModelCoordsMode == COORDS_LATLON) {
						//	rect2latlon(0, xval, yval, &yloc, &xloc);
						//} else {
							xloc = xval;
							yloc = yval;
						//}

						dx = xloc - Hypo.x;
						dy = yloc - Hypo.y;
						azimuth = atan2(dx, dy);  // radians CW from Y axis
						okada_az = -(azimuth - mech_strikeDir);
//fprintf(stdout,"azimuth %f  mech_strikeDir %f  okada_az %f\n", azimuth / cRPD, mech_strikeDir / cRPD, okada_az / cRPD);
						dist_horiz = sqrt(dx * dx + dy * dy);
						// observation point with x axis along
						obspoint[0] = dist_horiz * cos(okada_az);
						obspoint[1] = dist_horiz * sin(okada_az);

//						if (dist_horiz > 0.0) {

							zdepth = coul_grid.origz + coul_grid.dz / 2.0;
							for (iz = 0; iz <  coul_grid.numz; iz++) {

								dz = zdepth - Hypo.z;

								dist_xyz = sqrt(dx * dx + dy * dy + dz * dz);

								if (dist_xyz > 2.0 * coul_grid.dz) {

									obspoint[2] = -zdepth;

									// update Coulomb value
									if (1) {	// point shear dislocation
										coulomb = pointFaultShearStressStress(
										-Hypo.z, mech_dipAng,
										slip, obspoint, stress);
										//weight = exp(-dist_xyz / 5.0);
										weight = 1.0 / pow(dist_xyz, 2.5);
										//weight = 1.0;
										coulomb *= weight;
										//if (coulomb > 0.0) {
											//if (dist_xyz > 0.0)
											coul_grid.array[ix][iy][iz] += coulomb;
											//wt_grid.array[ix][iy][iz] += weight;
										//}
									} else {	// SH radiation
										ray_dip = acos(dz / dist_xyz);
										coulomb = calc_rad(ray_dip, azimuth, 'H');
										//coulomb = fabs(coulomb) - 0.31831;  // 1/PI
										coulomb = fabs(coulomb);
										coul_grid.array[ix][iy][iz] += coulomb;
									}

								}

								zdepth += coul_grid.dz;
							}
//						}
						yval += coul_grid.dy;
					}
					xval += coul_grid.dx;
				}


			nLocWritten++;

			}

		}

	}

	// normaloze output grid
	value_max_abs = 0.0;
	for (ix = 0; ix < coul_grid.numx; ix++)
		for (iy = 0; iy < coul_grid.numy; iy++)
			for (iz = 0; iz < coul_grid.numz; iz++) {
				//if (wt_grid.array[ix][iy][iz] > 1.0)
				//	coul_grid.array[ix][iy][iz] /= wt_grid.array[ix][iy][iz];
				if (coul_grid.array[ix][iy][iz] > value_max_abs)
					value_max_abs = coul_grid.array[ix][iy][iz];
			}
	for (ix = 0; ix < coul_grid.numx; ix++)
		for (iy = 0; iy < coul_grid.numy; iy++)
			for (iz = 0; iz < coul_grid.numz; iz++)
				coul_grid.array[ix][iy][iz] /= value_max_abs;


	// save samples if requested (before smoothing)
	if (n_samples > 0) {
		saveSamples(n_samples, &coul_grid, fn_root_out, 0.5, "50");
		saveSamples(n_samples, &coul_grid, fn_root_out, 0.7, "70");
		saveSamples(n_samples, &coul_grid, fn_root_out, 0.8, "80");
	}

	// smooth if requested
	if (fabs(n_smooth) > 0) {

		ReadGrid3dHdr(&smoothed_grid, NULL, grid_hdr_filename, "coul");
		/* allocate grid */
		smoothed_grid.buffer = AllocateGrid(&smoothed_grid);
		if (smoothed_grid.buffer == NULL) {
			nll_puterr("ERROR: allocating memory for 3D Coulomb smooth grid buffer.");
			exit(EXIT_ERROR_MEMORY);
		}
		/* create array access pointers */
		smoothed_grid.array = CreateGridArray(&smoothed_grid);
		if (smoothed_grid.array == NULL) {
			nll_puterr("ERROR: creating array for accessing 3D Coulomb smooth grid buffer.");
			exit(EXIT_ERROR_MEMORY);
		}

		if (n_smooth > 0) {

			value_max_abs = 0.0;
			for (ix = 0; ix < smoothed_grid.numx; ix++)
				for (iy = 0; iy < smoothed_grid.numy; iy++)
					for (iz = 0; iz < smoothed_grid.numz; iz++) {
						smooth_sum = 0.0;
						for (ny = -n_smooth; ny <= n_smooth; ny++) {
							if ((iy + ny) >= 0 && (iy + ny) < coul_grid.numy)
								smooth_sum += coul_grid.array[ix][iy + ny][iz]
									* (1.0 -
									((double) fabs(ny)) / (double) (n_smooth + 1));
						}
						smoothed_grid.array[ix][iy][iz] = smooth_sum;
						// normalization constant
						if (smoothed_grid.array[ix][iy][iz] > value_max_abs)
							value_max_abs = smoothed_grid.array[ix][iy][iz];
					}

		} else {

			n_smooth = -n_smooth;
			value_max_abs = 0.0;
			for (ix = 0; ix < smoothed_grid.numx; ix++)
				for (iy = 0; iy < smoothed_grid.numy; iy++)
					for (iz = 0; iz < smoothed_grid.numz; iz++) {
						smooth_sum = 0.0;
						for (nx = -n_smooth; nx <= n_smooth; nx++) {
							if ((ix + nx) >= 0 && (ix + nx) < coul_grid.numx)
								smooth_sum += coul_grid.array[ix + nx][iy][iz]
									* (1.0 -
									((double) fabs(nx)) / (double) (n_smooth + 1));
						}
						smoothed_grid.array[ix][iy][iz] = smooth_sum;
						// normalization constant
						if (smoothed_grid.array[ix][iy][iz] > value_max_abs)
							value_max_abs = smoothed_grid.array[ix][iy][iz];
					}

		}

		coul_grid = smoothed_grid;

	}

	// normaloze output grid
	for (ix = 0; ix < coul_grid.numx; ix++)
		for (iy = 0; iy < coul_grid.numy; iy++)
			for (iz = 0; iz < coul_grid.numz; iz++)
				coul_grid.array[ix][iy][iz] /= value_max_abs;


	// save grid to disk
	if (istat =
		WriteGrid3dBuf(&coul_grid, NULL, fn_root_out, "coul") < 0) {
		nll_puterr("ERROR: writing slowness grid to disk.");
		exit(EXIT_ERROR_IO);
	}

	// write message
	fprintf(stdout,
		"%d location files read, %d events accepted.\n", numFiles, nLocAccepted);
	fprintf(stdout, "%d Coulomb patters written to grid <%s.*>\n", nLocWritten, fn_root_out);




	return(0);

}



/*** function to calculate maxium shear stress on fault plane */

double pointFaultShearStressStress(double c, double dip,
		      const double *slip, const double *obspoint,
		      double *stress) {

	//double fault_normal[3] = {0.0, 0.0, 1.0};
	//double slip0[3] = {1.0, 0.0, 0.0};

	PointStressResponse(c, dip, slip, obspoint, stress);
	//return(stress[8]);	// normal only
	//return(stress[2]);	// shear only
	return(stress[2] - 0.5 * stress[8]);	// Coulomb

	//return(faultShearStress(slip0, fault_normal, stress));

}



/*** function to calculate maxium shear stress on fault plane */

// from http://richter.colorado.edu/~sethmc/thesis/node8.html

// for calc fault plane is horizonatal (normal in z-dir)
// slip is in x-dir

// returns stress tensor in the fault plane coordinate frame

#define ALPHA 5.0
#define BETA 3.0

void PointStressResponse(double c, double dip,
		      const double *slip, const double *obspoint,
		      double *stress)
{

/*
     s = 0 1 2
         3 4 5
         6 7 8 ) */

	double x0, y0, y1, z0, r0;
	double x, y, z, r;
	double x2, y2, z2, r2, c1, a2, inv_r5;
	double cosdip, sindip, temp;
	double cosslip, sinslip;

	x0 = obspoint[0];
	y0 = obspoint[1];
	z0 = obspoint[2] - c;

	// unit sphere
	//r0 = sqrt(x0 * x0 + y0 * y0 + z0 * z0);
	//temp = r0;
	//x0 /= temp;
	//y0 /= temp;
	//z0 /= temp;

	// rotate to horizontal fault plane
	cosdip = cos(dip);
	sindip = sin(dip);
	y1 = y0 * cosdip + z0 * sindip;
	z = -y0 * sindip + z0 * cosdip;

	// rotate to slip direction
	temp = sqrt(slip[0] * slip[0] + slip[1] * slip[1]);
	cosslip = slip[0] / temp;
	sinslip = slip[1] / temp;
	y = y1 * cosslip - x0 * sinslip;
	x = y1 * sinslip + x0 * cosslip;


	x2 = x * x;
	y2 = y * y;
	z2 = z * z;
	r2 = x2 + y2 + z2;
	r = sqrt(r2);
	//inv_r5 = 1.0 / pow(r, 5);
	//inv_r5 = 1.0 / pow(r, 3);
	//inv_r5 = 1.0 / r0;
	//inv_r5 = 1.0 / sqrt(r0);
	inv_r5 = 1.0;

	c1 = 1.0 / (BETA * BETA) - 1.0 / (ALPHA * ALPHA);
	a2 = ALPHA * ALPHA;

	//stress[0] = inv_r5 * 2.0 * x * z * ( c1 * (2.0 - 5.0 * x2 / r2) - 2.0 / a2 );
	//stress[4] = inv_r5 * 2.0 * x * z * ( c1 * (1.0 - 5.0 * y2 / r2) - 1.0 / a2 );
	stress[8] = inv_r5 * 2.0 * x * z * ( c1 * (2.0 - 5.0 * z2 / r2) - 2.0 / a2 );
	//stress[1] = stress[3] = inv_r5 * y * z * ( c1 * (1.0 - 10.0 * x2 / r2) - 1.0 / a2 );
	//stress[5] = stress[7] = inv_r5 * x * y * ( c1 * (1.0 - 10.0 * z2 / r2) - 1.0 / a2 );
	stress[2] = stress[6] = inv_r5 * ( c1 * (x2 + z2 - 10.0 * x2 * z2 / r2) + (3 * y2 - r2) / (3 * a2) );

}



/*** function to read input file */

int ReadCoul_Input(FILE* fp_input)
{
	int istat, iscan;
	char param[MAXLINE];
	char line[MAXLINE];

	int flag_control = 0, flag_trans = 0, flag_grid = 0;


	/* read each input line */

	while (fgets(line, MAXLINE, fp_input) != NULL) {

		istat = -1;

		/*read parameter line */

		if ((iscan = sscanf(line, "%s", param)) < 0 )
			continue;

		/* skip comment line or white space */

		if (strncmp(param, "#", 1) == 0 || isspace(param[0]))
			istat = 0;


		/* read control params */

		if (strncmp(param, "CONTROL", 6) == 0)
			if ((istat = get_control(strchr(line, ' '))) < 0)
				nll_puterr("Error reading control params.");
			else
				flag_control = 1;


		/*read transform params */

		// input location trans
		if (strncmp(param, "TRANS", 5) == 0)
			if ((istat = get_transform(0, strchr(line, ' '))) < 0)
			    nll_puterr("ERROR reading transformation parameters.");
			else {
				flag_trans = 1;
			}



		/* unrecognized input */

		if (istat < 0 && message_flag > 3) {
			fprintf(stdout, "Skipping input: %s", line);
		}

	}


	/* check for missing input */

	if (!flag_control)
		nll_puterr("ERROR no control (CONTROL) params read.");
	if (!flag_trans)
		nll_puterr("ERROR no transformation (TRANS) params read.");


	return (flag_control + flag_trans - 1);
}


/*** function to generate and save samples */

int saveSamples(int n_samples, GridDesc* ptgrid, char* fileroot, double cutoff, char* code) {

	int istat;
	int ix, iy, iz;
	float ftemp;
	double sum;
	char filename[FILENAME_MAX];
	GridDesc samples_grid;


	// create new grid
	DuplicateGrid(&samples_grid , ptgrid , "COULOMB");

	// convert output grid to prob density - sum=1
	sum = 0.0;
	for (ix = 0; ix < ptgrid->numx; ix++)
		for (iy = 0; iy < ptgrid->numy; iy++)
			for (iz = 0; iz < ptgrid->numz; iz++) {
				ftemp = ptgrid->array[ix][iy][iz];
				ftemp = ftemp > cutoff ? ftemp - cutoff : 0.0;
				ftemp = pow(ftemp, 2);
				samples_grid.array[ix][iy][iz] = ftemp;
				sum += samples_grid.array[ix][iy][iz];
			}
	sum *= (samples_grid.dx * samples_grid.dy * samples_grid.dz);
	for (ix = 0; ix < samples_grid.numx; ix++)
		for (iy = 0; iy < samples_grid.numy; iy++)
			for (iz = 0; iz < samples_grid.numz; iz++)
				samples_grid.array[ix][iy][iz] /= sum;

	strcpy(filename, fileroot);
	strcat(filename, ".");
	strcat(filename, code);
	strcat(filename, ".coul");

	istat = GenEventScatterGrid1(&samples_grid, 1.0, n_samples, filename);

	return(istat);

}



/*** function to generate sample (scatter) of location PDF */

int GenEventScatterGrid1(GridDesc* ptgrid, double probmax, int n_samples,
	char* filename)
{
	FILE *fpio;
	char fname[FILENAME_MAX];

	int ix, iy, iz;
	double origx, origy, origz;
	double dx, dy, dz, dvol;
	double xval, yval, zval;
	float fdata[4];
	int tot_npoints = 0;
	double xnpt, xnpoints;
	double prob_den, prob;



	/* return if no scatter samples requested */
	if (n_samples < 1)
		return(0);

	/* write message */
	putmsg(3, "");
	putmsg(3, "Generating event scatter file...");

	/* open scatter file */

	sprintf(fname, "%s.scat", filename);
//printf("%s\n", fname);
	if ((fpio = fopen(fname, "w")) == NULL) {
		nll_puterr("ERROR: opening scatter output file.");
		return(-1);
	} else {
		NumFilesOpen++;
	}
	/* skip header record (used later to store number of samples taken) */
	fseek(fpio, 4 * sizeof(float), SEEK_SET);


	/* generate N=Scatter->npts events with prob P=prob_den/probmax at  */
	/*	uniformly-randomly chosen grid locations */

	origx = ptgrid->origx;
	origy = ptgrid->origy;
	origz = ptgrid->origz;
	dx = ptgrid->dx;
	dy = ptgrid->dy;
	dz = ptgrid->dz;
	dvol = dx * dy * dz;
//printf("%lf %lf %lf %lf %lf %lf\n", origx, origy, origz, dx, dy, dz);


	for (ix = 0; ix <  ptgrid->numx; ix++) {
	    for (iy = 0; iy <  ptgrid->numy; iy++) {
		for (iz = 0; iz <  ptgrid->numz; iz++) {

			prob_den =  ptgrid->array[ix][iy][iz];

			xnpoints = (double) n_samples * dvol * prob_den;
//printf("%d %lf %lf\n", n_samples, dvol, prob_den);

			xval = origx + (double) ix * dx + dx / 2.0;
			yval = origy + (double) iy * dy + dy / 2.0;
			zval = origz + (double) iz * dz + dz / 2.0;
//printf("%lf %lf %lf %lf\n", xval, yval, zval, xnpoints);
//printf("rand xval = %lf\n", get_rand_double(-dx / 2.0, dx / 2.0));
//printf("rand yval = %lf\n", get_rand_double(-dy / 2.0, dy / 2.0));
//printf("rand zval = %lf\n", get_rand_double(-dz / 2.0, dz / 2.0));

			while (xnpoints > 0.0) {

			    if (xnpoints > 1.0 ||
					xnpoints - (double) ((int) xnpoints)
						> get_rand_double(0.0, 1.0)) {
				fdata[0] = (float) (xval + get_rand_double(-dx / 2.0, dx / 2.0));
				fdata[1] = (float) (yval + get_rand_double(-dy / 2.0, dy / 2.0));
				fdata[2] = (float) (zval + get_rand_double(-dz / 2.0, dz / 2.0));
				fdata[3] = prob_den;
//printf("%f\n", (float) (xval + get_rand_double(-dx / 2.0, dx / 2.0)));
//printf("%f %f %f %f\n", fdata[0], fdata[1], fdata[2], fdata[3]);
				fwrite(fdata, sizeof(float), 4, fpio);

				tot_npoints++;
			    }

			    xnpoints -= 1.0;

			}

		}
	    }
	}


	/* write header informaion */
	fseek(fpio, 0, SEEK_SET);
	fwrite(&tot_npoints, sizeof(int), 1, fpio);
	fdata[0] = (float) probmax;
	fwrite(fdata, sizeof(float), 1, fpio);

	fclose(fpio);
	NumFilesOpen--;

	/* write message */
	sprintf(MsgStr, "  %d points generated.", tot_npoints);
	putmsg(3, MsgStr);
	sprintf(MsgStr, "  (%d points requested, dvol= %lf, probmax=%lf)",
			n_samples, dvol, probmax);
	putmsg(3, MsgStr);

	return(0);

}



/*** function to calculate radiation amplitude */

/* A & R figs 4.20 & 5.5 */

double calc_rad(double ray_vert_ang, double ray_horiz_ang, char orig_wave)
{
	double radamp;

		/* cal radiation pattern (from Aki & Richards eqs. 4.84 - 4.86) */

	if (orig_wave == 'P') {
		radamp = cos(mech_rake) * sin(mech_dipAng) * pow(sin(ray_vert_ang),2.0)
				* sin(2.0 * (ray_horiz_ang - mech_dipDir))
			- cos(mech_rake) * cos(mech_dipAng) * sin(2.0 * ray_vert_ang)
				* cos(ray_horiz_ang - mech_dipDir)
			+ sin(mech_rake) * sin(2.0 * mech_dipAng)
				* (pow(cos(ray_vert_ang),2.0) - pow(sin(ray_vert_ang)
				* sin(ray_horiz_ang - mech_dipDir), 2.0))
			+ sin(mech_rake) * cos(2.0 * mech_dipAng) * sin(2.0 * ray_vert_ang)
				* sin(ray_horiz_ang - mech_dipDir)
		;
	} else if (orig_wave == 'V') {
		radamp = sin(mech_rake) * cos(2.0 * mech_dipAng) * cos(2.0 * ray_vert_ang)
				* sin(ray_horiz_ang - mech_dipDir)
			- cos(mech_rake) * cos(mech_dipAng) * cos(2.0 * ray_vert_ang)
				* cos(ray_horiz_ang - mech_dipDir)
			+ 0.5 * cos(mech_rake) * sin(mech_dipAng)
				* sin(2.0 * ray_vert_ang) * sin(2.0 * (ray_horiz_ang - mech_dipDir))
			- 0.5 *  sin(mech_rake) * sin(2.0 * mech_dipAng) * sin(2.0 * ray_vert_ang)
				* (1.0 + pow(sin(ray_horiz_ang - mech_dipDir), 2.0))
		;
		radamp *= -1.0;
			/* mult by -1 since A & R change dir of + in figs 4.20 & 5.5 */
	} else if (orig_wave == 'H') {
		radamp = cos(mech_rake) * cos(mech_dipAng) * cos(ray_vert_ang)
				* sin(ray_horiz_ang - mech_dipDir)
			+ cos(mech_rake) * sin(mech_dipAng) * sin(ray_vert_ang)
				* cos(2.0 * (ray_horiz_ang - mech_dipDir))
			+ sin(mech_rake) * cos(2.0 * mech_dipAng)
				* cos(ray_vert_ang) * cos(ray_horiz_ang - mech_dipDir)
			- 0.5 * sin(mech_rake) * sin(2.0 * mech_dipAng) * sin(ray_vert_ang)
				* sin(2.0 * (ray_horiz_ang - mech_dipDir))
		;
	} else {
		sprintf(MsgStr,
			"WARNING: Mechanism: unrecognized original wave type: %c", orig_wave);
		putmsg(2, MsgStr);
	}
	return(radamp);
}








/*** function to calculate maxium shear stress on fault plane */

double faultShearStress(const double *slip, const double *fault_normal, double *stress)
{
/*
     s = 0 1 2
         3 4 5
         6 7 8 ) */

	double traction[3];
	double mag_traction_2;
	double mag_normal_traction_2;
	double mag_shear_traction_2;
	double sign;

	// mulitply stress tensor by fault_normal to get traction across fault_normal
	traction[0] = stress[0] * fault_normal[0] + stress[1] * fault_normal[1] + stress[2] * fault_normal[2];
	traction[1] = stress[3] * fault_normal[0] + stress[4] * fault_normal[1] + stress[5] * fault_normal[2];
	traction[2] = stress[6] * fault_normal[0] + stress[7] * fault_normal[1] + stress[8] * fault_normal[2];

	// mag traction
	mag_traction_2 = traction[0] * traction[0]
		+ traction[1] * traction[1]
		+ traction[2] * traction[2];

	//  mag normal traction = traction DOT normal
	mag_normal_traction_2 = traction[0] * fault_normal[0]
		+ traction[1] * fault_normal[1]
		+ traction[2] * fault_normal[2];
	mag_normal_traction_2 *= mag_normal_traction_2;

	// mag shear traction
	mag_shear_traction_2 = mag_traction_2 - mag_normal_traction_2;

	// sign = traction DOT slip
	sign = traction[0] * slip[0]
		+ traction[1] * slip[1]
		+ traction[2] * slip[2];
	sign /= fabs(sign);

	return(sign * sqrt(mag_shear_traction_2));

}


