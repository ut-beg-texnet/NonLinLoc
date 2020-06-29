/*
 * Copyright (C) 2002 Anthony Lomax <anthony@alomax.net>
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


/*   calcpd.c

	Program to calculate partial derivatives at the source for a list of stations

*/

/*------------------------------------------------------------*/
/* Anthony Lomax           | email: anthony@alomax.net        */
/* Scientific Software     | web: www.alomax.net              */
/* Mouans-Sartoux, FRANCE  |                                  */
/*------------------------------------------------------------*/


/*
	history:

	ver 01    15FEB2002  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

*/



#define PNAME  "calcpd"

#include "../src/GridLib.h"


/* defines */

/* globals  */
char fn_time_grids[MAXLINE], fn_model_grid[MAXLINE];
int NumStations, NumSources;
double VpVsRatio;



/* function declarations */

int ReadCalcpd_Input(FILE* );
int GetCalcpd_Files(char* );
double ReadGridValue(FILE* , GridDesc* , HypoDesc* , SourceDesc* );
Vect3D CalcPartialDerivs(char *filename, GridDesc* ptgrid,
			HypoDesc*  pevent, SourceDesc* psta,
			double vel_source, ArrivalDesc* parr);
int WritePhaseArrival(FILE* fp_out, ArrivalDesc* parr,
		double travel_time, double residual, double distance,
		Vect3D deriv, SourceDesc* psta);
int WriteHypocenter(FILE* fp_out, HypoDesc* phypo, double vel_p, double vel_s);



/*** Program to calculate partial derivatives at the source for a list of stations */

#define NARGS 2

main(int argc, char *argv[])
{

	int istat;

	FILE *fp_in;
	FILE *fp_out;

	/* set program name */
	strcpy(prog_name, PNAME);

	/* check command line for correct usage */

	if (argc != NARGS) {
		disp_usage(prog_name, "<control file>");
		exit(EXIT_ERROR_USAGE);
	}



	/* set constants */

	SetConstants();
	NumStations = 0;
	NumSources = 0;
	VpVsRatio = -1.0;



	/* read control file */

	strcpy(fn_control, argv[1]);
	if ((fp_control = fopen(fn_control, "r")) == NULL) {
		nll_puterr("ERROR: opening control file.");
		exit(EXIT_ERROR_FILEIO);
	}

	if ((istat = ReadCalcpd_Input(fp_control)) < 0) {
		exit(EXIT_ERROR_FILEIO);
	}
	fclose(fp_control);
	
	// write format to stderr
	fprintf(stderr, "Format: Hypocenter line:\n");
	fprintf(stderr, "year month day  hour min sec  lat long depth  mag  Vp Vs\n");
	fprintf(stderr, "Format: Phase lines:\n");
	fprintf(stderr, "station onset phase first_mot quality sec  pred_travel_time residual distance  ");
	fprintf(stderr, "dT/dx dT/dy dT/dz  ray_azim ray_dip ray_qual  sta.x sta.y sta.z\n");
		


	/* set input and output files */

	fp_in = stdin;
	fp_out = stdout;
	
	calcpd_alberto(fp_in, fp_out);
	
}


/*** function to calculate partial derivtives */

int calcpd_alberto(FILE* fp_in, FILE* fp_out)
{

	int istat;
	char line_in[MAXLINE_LONG];
	char *cstat;
	
	int lat_deg;
	char lat_str[2];
	char lat_chr;
	double lat_min;
	int lon_deg;
	char lon_str[2];
	char lon_chr;
	double lon_min;
	char sign;
	
	char chr;
	char phase_code[12];
	
	char filename[MAXLINE], filename_model[MAXLINE], filename_angle[MAXLINE];
	FILE *fp_grid, *fp_hdr;
	GridDesc grid;
	
	HypoDesc hypo;
	SourceDesc station;
	ArrivalDesc arrival;
	
	double residual, travel_time, vel_source;
	double ray_azim, ray_dip;
	int ray_qual;
	Vect3D deriv;
	double value, ratio;


	/* generate arrival times for each source */

	while (1) {

		// read hypocenter line
		cstat = fgets(line_in, MAXLINE_LONG, fp_in);
//printf(">%s", line_in);
		if (cstat == NULL)
			return(0);
		istat = ReadFortranInt(line_in, 1, 2, &hypo.year);
		hypo.year += 1900;
		istat += ReadFortranInt(line_in, 3, 2, &hypo.month);
		istat += ReadFortranInt(line_in, 5, 2, &hypo.day);
		istat += ReadFortranInt(line_in, 8, 2, &hypo.hour);
		istat += ReadFortranInt(line_in, 10, 2, &hypo.min);
		istat += ReadFortranReal(line_in, 13, 5, &hypo.sec);
		istat += ReadFortranInt(line_in, 19, 2, &lat_deg);
		istat += ReadFortranString(line_in, 21, 1, lat_str);
		istat += ReadFortranReal(line_in, 22, 5, &lat_min);
		istat += ReadFortranInt(line_in, 28, 3, &lon_deg);
		istat += ReadFortranString(line_in, 31, 1, lon_str);
		istat += ReadFortranReal(line_in, 32, 5, &lon_min);
		istat += ReadFortranReal(line_in, 38, 6, &hypo.depth);
		istat += ReadFortranReal(line_in, 47, 4, &hypo.amp_mag);
		
		if (istat != 14) {
			nll_puterr2("ERROR: reading hypocenter line:", line_in);
			return(-1);
		}
		
		// decode lat/long
		lat_chr = lat_str[0];
		if (lat_chr == ' ')
			lat_chr = 'N';
		lon_chr = lon_str[0];
		if (lon_chr == ' ')
			lon_chr = 'W';
		sign = toupper(lat_chr) == 'N' ? 1.0 : -1.0;
		hypo.dlat = sign * ((double) lat_deg + lat_min / 60.0);
		sign = toupper(lon_chr) == 'E' ? 1.0 : -1.0;
		hypo.dlong = sign * ((double) lon_deg + lon_min / 60.0);
		// convert to x/y
		istat = latlon2rect(0, hypo.dlat, hypo.dlong, &(hypo.x), &(hypo.y));
		hypo.z = hypo.depth;

		sprintf(filename_model, "%s.%s.mod", fn_model_grid, "P");
		if ((istat = OpenGrid3dFile(filename_model, &fp_grid,
				&fp_hdr, &grid, "model", &station)) < 0) {
			CloseGrid3dFile(&fp_grid, &fp_hdr);
			nll_puterr2("ERROR: opening model grid files: ", filename_model);
			continue;
		}
		
		/* get velocity at source */
		value = ReadGridValue(fp_grid, &grid, &hypo, NULL);
		if (grid.type == GRID_VELOCITY)
			vel_source = value;
		else if (grid.type == GRID_VELOCITY_METERS)
			vel_source = value / 1000.0;
		else if (grid.type == GRID_SLOWNESS)
			vel_source = 1.0 / (value);
		else if (grid.type == GRID_SLOW_LEN)
			vel_source = 1.0 / (value / grid.dx);
		else if (grid.type == GRID_VEL2)
			vel_source = sqrt(value);
		else if (grid.type == GRID_SLOW2)
			vel_source = sqrt(1.0 / value);
		else if (grid.type == GRID_SLOW2_METERS)
			vel_source = sqrt(1.0 / value) / 1000.0;
		CloseGrid3dFile(&fp_grid, &fp_hdr);

		if ((istat = WriteHypocenter(fp_out, &hypo, vel_source, vel_source / VpVsRatio)) < 0)
			nll_puterr("ERROR: writing hypocenter.");
		
		
		// read phases
		
		while (1) {
			
			// check for end of event (assumes no blanks after last phase)
			chr = fgetc(fp_in);
			if (chr == EOF) {
				return(0);
			} else if (chr == '\n') {
				if ((chr = fgetc(fp_in)) != EOF && !isdigit(chr)) {
					// another phase follows
					ungetc(chr, fp_in);
				} else if (chr == '0') {
					// end of event "0" line
					// read to end of line
					while ((chr = fgetc(fp_in)) != EOF && chr != '\n')
						;
					break;
				} else {
					// end of event
					ungetc(chr, fp_in);
					break;
				}
			} else {
				ungetc(chr, fp_in);
			}


			// read phase arrival input
			istat = fscanf(fp_in, "%15c", line_in);
			line_in[15] = '\0';
//printf("<%s>\n", line_in);
			line_in[15]= '\0';
			if (istat == EOF)
				return(0);
			istat = ReadFortranString(line_in, 1, 4, arrival.label);
			istat += ReadFortranString(line_in, 5, 1, arrival.onset);
			if (arrival.onset[0] == ' ')
				arrival.onset[0] = '?';
			istat += ReadFortranString(line_in, 6, 1, arrival.phase);
			if (arrival.phase[0] == ' ')
				arrival.phase[0] = '?';
			istat += ReadFortranString(line_in, 7, 1, arrival.first_mot);
			if (arrival.first_mot[0] == ' ')
				arrival.first_mot[0] = '?';
			istat += ReadFortranInt(line_in, 8, 1, &arrival.quality);
			istat += ReadFortranReal(line_in, 9, 7, &arrival.sec);

			if (istat != 6) {
				nll_puterr2("ERROR: reading phase:", line_in);
				return(-1);
			}

		
			// check for S and valid VpVs
			if (strcmp(arrival.phase, "S") == 0 && VpVsRatio > 0.0) {
				strcpy(phase_code, "P");
				ratio = VpVsRatio;
			} else {
				strcpy(phase_code, arrival.phase);
				ratio = 1.0;
			}
			
			// get travel time
			sprintf(filename, "%s.%s.%s.time", fn_time_grids,
				phase_code, arrival.label);
			if ((istat = OpenGrid3dFile(filename, &fp_grid,
					&fp_hdr, &grid, "time", &station)) < 0) {
				CloseGrid3dFile(&fp_grid, &fp_hdr);
				nll_puterr2("ERROR: opening time grid files: ", filename);
				continue;
			}
			travel_time = ReadGridValue(fp_grid, &grid, &hypo, &station);
			CloseGrid3dFile(&fp_grid, &fp_hdr);
			sprintf(MsgStr,
"Calculating travel time for phase: %s  %s  X %.2lf  Y %.2lf  Z %.2lf",
					arrival.label, arrival.phase,
					station.x, station.y, station.z);
			//nll_puterr(MsgStr);
			if (travel_time < 0.0) {
				sprintf(MsgStr, "ERROR: calculating travel time: t = %lf", travel_time);
				nll_puterr(MsgStr);
				continue;
			}

			// calculate residual
			residual = arrival.sec - hypo.sec - travel_time;
		
			// calculate angles and partial derivatives
			sprintf(filename_angle, "%s.%s.%s.angle", fn_time_grids,
				phase_code, arrival.label);
			deriv = CalcPartialDerivs(filename_angle, &grid, &hypo, &station,
				vel_source, &arrival);
				
			if ((istat = WritePhaseArrival(fp_out, &arrival,
					travel_time, residual,
					GetEpiDist(&station, hypo.x, hypo.y), deriv, &station)) < 0)
				nll_puterr("ERROR: writing phase arrival.");

		}
		
		// next hypocenter
		fprintf(fp_out, "\n");
	}


	return(0);

}




/*** function to read a value from  a time grid file */

double ReadGridValue(FILE* fpgrid, GridDesc* ptgrid,
			HypoDesc*  pevent, SourceDesc* psta)
{

	int istat;
	double value, yval_grid;

	/* get travel time */

	if (ptgrid->type == GRID_TIME) {
		/* 3D grid */
		value = ReadAbsInterpGrid3d(fpgrid, ptgrid,
			pevent->x, pevent->y, pevent->z);
	} else {
		/* 2D grid (1D model) */
		if (psta == NULL)
			yval_grid = 0.0;
		else
			yval_grid = GetEpiDist(psta, pevent->x, pevent->y);
		value =
			ReadAbsInterpGrid2d(fpgrid, ptgrid, yval_grid, pevent->z);
	}


	return(value);

}


/*** function to calc first motion from angle grid file */

Vect3D CalcPartialDerivs(char *filename, GridDesc* ptgrid,
			HypoDesc*  pevent, SourceDesc* psta,
			double vel_source, ArrivalDesc* parr)
{

	int istat;
	double yval_grid, azim, radamp;
	int ipolarity = 0;
	double sind;
	
	Vect3D deriv;
	
	deriv.x = deriv.y = deriv.z = -LARGE_DOUBLE;

	/* get take-off angles */

	if (ptgrid->type == GRID_TIME) {
		/* 3D grid */
		if (ReadTakeOffAnglesFile(filename,
				pevent->x, pevent->y, pevent->z,
				&(parr->ray_azim), &(parr->ray_dip), &(parr->ray_qual), -1.0) < 0)
			return(deriv);
	} else {
		/* 2D grid (1D model) */
		yval_grid = GetEpiDist(psta, pevent->x, pevent->y);
		azim = GetEpiAzim(psta, pevent->x, pevent->y);
		if (ReadTakeOffAnglesFile(filename,
				0.0, yval_grid, pevent->z,
				&(parr->ray_azim), &(parr->ray_dip), &(parr->ray_qual), azim) < 0)
			return(deriv);
	}
	
	/* calc partial derivatives */
//vel_source = 1.0;
	deriv.z = -cos(rpd * parr->ray_dip) / vel_source;
	sind = sin(rpd * parr->ray_dip);
	deriv.x = -sin(rpd * parr->ray_azim) * sind / vel_source;
	deriv.y = -cos(rpd * parr->ray_azim) * sind / vel_source;

	return(deriv);

}



/*** function to write arrival to file */

int WritePhaseArrival(FILE* fp_out, ArrivalDesc* parr,
		double travel_time, double residual, double distance, Vect3D deriv, SourceDesc* psta)
{


	/* write arrival to output file */

	fprintf(fp_out, " %s %s %s %s %d %lf   %lf %lf %lf   %lf %lf %lf   %lf %lf %d    %lf %lf %lf\n",
		parr->label, parr->onset, parr->phase, parr->first_mot, parr->quality, parr->sec,
		travel_time, residual, distance, deriv.x, deriv.y, deriv.z,
		parr->ray_azim, parr->ray_dip, parr->ray_qual, psta->x, psta->y, psta->z
	);


	return(0);

}



/*** function to write hypocenter to file */

int WriteHypocenter(FILE* fp_out, HypoDesc* phypo, double vel_p, double vel_s)
{


	/* write hypocenter to output file */

	fprintf(fp_out, "%d %d %d  %d %d %lf  %lf %lf %lf  %lf  %lf %lf\n",
		phypo->year, phypo->month, phypo->day, phypo->hour, phypo->min, phypo->sec,
		phypo->dlat, phypo->dlong, phypo->depth, phypo->amp_mag,
		vel_p, vel_s
	);


	return(0);

}





/*** function to read input file */

int ReadCalcpd_Input(FILE* fp_input)
{
	int istat, iscan;
	char param[MAXLINE];
	char line_in[4*MAXLINE], *fgets_return;

	int flag_control = 0, flag_files = 0,
		flag_trans = 0, flag_vp_vs = 0;

	int flag_include = 1;


	/* read each input line */

	while ((fgets_return = fgets(line_in, 4*MAXLINE, fp_input)) != NULL
			|| fp_include != NULL) {


		/* check for end of include file */

		if (fgets_return == NULL && fp_include != NULL) {
			SwapBackIncludeFP(&fp_input);
			continue;
		}


		istat = -1;

		/*read parameter line */

		if ((iscan = sscanf(line_in, "%s", param)) < 0 )
			continue;

		/* skip comment line or white space */

		if (strncmp(param, "#", 1) == 0 || isspace(param[0]))
			istat = 0;


		/* read include file params and set input to include file */

		if (strcmp(param, "INCLUDE") == 0)
			if ((istat = GetIncludeFile(strchr(line_in, ' '),
							&fp_input)) < 0) {
				nll_puterr("ERROR: processing include file.");
				flag_include = 0;
			}


		/* read control params */

		if (strcmp(param, "CONTROL") == 0)
			if ((istat = get_control_proxy(strchr(line_in, ' '))) < 0)
				nll_puterr("ERROR: reading control params.");
			else
				flag_control = 1;


		/*read transform params */

		if (strcmp(param, "TRANS") == 0)
    			if ((istat = get_transform(0, strchr(line_in, ' '))) < 0)
			    nll_puterr("ERROR: reading transformation parameters.");
			else
				flag_trans = 1;


		/* read file names */

		if (strcmp(param, "CPDFILES") == 0)
			if ((istat = GetCalcpd_Files(strchr(line_in, ' '))) < 0)
			  nll_puterr("ERROR: reading output file name.");
			else
				flag_files = 1;


		/* read VpVs params */

		if (strcmp(param, "CPDVPVS") == 0)
			if ((istat = get_vp_vs(strchr(line_in, ' '))) < 0)
				nll_puterr("ERROR: reading calcpd Vp/Vs params.");
			else
				flag_vp_vs = 1;


		/* unrecognized input */

		if (istat < 0) {
			sprintf(MsgStr, "Skipping input: %s", line_in);
			putmsg(4, MsgStr);
		}

	}


	/* check for missing input */

	if (!flag_control)
		nll_puterr("ERROR: no control (CONTROL) params read.");
	if (!flag_trans)
		nll_puterr("INFO: no transformation (TRANS) params read.");
	if (!flag_files)
		nll_puterr("ERROR: no i/o file (CPDFILES) params read.");
	
	if (!flag_vp_vs)
		nll_puterr("INFO: no Vp/Vs (CPDVPVS) params read.");


	return (flag_include * flag_control * flag_files - 1);
}



/*** function to read file names ***/

int GetCalcpd_Files(char* line1)
{

	sscanf(line1, "%s %s", fn_model_grid, fn_time_grids);

	fprintf(stderr, "CPDFILES:  Model grids: %s.*  Time grid: %s\n", fn_model_grid, fn_time_grids);

	return(0);
}



/*** function to read Vp / Vs ratio ***/

int get_vp_vs(char* line1)
{
	int istat;

	istat = sscanf(line1, "%lf", &VpVsRatio);

	fprintf(stderr,"CPDVPVS: VpVsRatio=%lf\n", VpVsRatio);

	if (istat != 1)
		return(-1);

	return(0);

}




/*** function to read control params ***/

int get_control_proxy(char* line1)
{
	int istat;

	istat = sscanf(line1, "%d", &message_flag, &RandomNumSeed);

	if (istat == 1)
		RandomNumSeed = 837465;

	/* display program information */
	//DispProgInfo();

	sprintf(MsgStr, "CONTROL:  MessageFlag: %d  RandomNumSeed: %d",
		message_flag, RandomNumSeed);
	putmsg(3, MsgStr);

	if (checkRangeInt("CONTROL", "MessageFlag", message_flag, 1, 0, 0, 0) != 0)
		return(-1);

	if (istat != 1 && istat != 2)
		return(-1);

	return(0);

}



/*------------------------------------------------------------/ */
/* Anthony Lomax           | email: lomax@faille.unice.fr     / */
/* UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax / */
/* 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        / */
/* 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        / */
/*------------------------------------------------------------/ */

