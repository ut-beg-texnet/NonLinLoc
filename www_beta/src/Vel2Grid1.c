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


/*   Vel2Grid.c

	Program to generate 3-D vel/slowness model grid from velocity model description

*/

/*------------------------------------------------------------/ */
/* Anthony Lomax           | email: lomax@faille.unice.fr     / */
/* UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax / */
/* 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        / */
/* 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        / */
/*------------------------------------------------------------/ */


/*
	history:

	ver 01    22SEP1997  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

*/



#define PNAME  "Vel2Grid"

#include "GridLib.h"
#include "velmod.h"



/* globals  */

char fn_vg_output[MAXLINE];

/* wave type (P, S, ...) for vel grids */
#define MAX_NUM_WAVE_TYPES 10
char WaveType[MAX_NUM_WAVE_TYPES][12];
int NumWaveTypes;


/* function declarations */

int ReadVel2GridInput(FILE* );
int VelModToGrid3d(GridDesc* , char * );
int get_vg_outfile(char* );
int get_vg_type(char* );



/*** program to generate  3-D vel/slowness grid */


#define NARGS 2

int main(int argc, char *argv[])
{

	int istat;
	int nWaveType;
	char fileRoot[MAXLINE];

	GridDesc mod_grid;	/* model grid */



	/* set program name */
	strcpy(prog_name, PNAME);

	/* check command line for correct usage */

	if (argc != NARGS) {
		disp_usage(prog_name, "<control file>");
		exit(EXIT_ERROR_USAGE);
	}



	/* set constants */

	prog_mode_3d = 1;
	prog_mode_Mod2D3D = 0;
	NumWaveTypes = 0;
	SetConstants();



	/* read control file */

	strcpy(fn_control, argv[1]);
	if ((fp_control = fopen(fn_control, "r")) == NULL) {
		nll_puterr("ERROR: opening control file.");
		exit(EXIT_ERROR_FILEIO);
	}

	if ((istat = ReadVel2GridInput(fp_control)) < 0) {
		exit(EXIT_ERROR_FILEIO);
	}
	mod_grid = grid_in;

	/* determine model coordinates mode - rect or latlon */
	SetModelCoordsMode(num_surfaces);


	/* initialize 3D grid */

	/* allocate model grid */
	mod_grid.buffer = AllocateGrid(&mod_grid);
	if (mod_grid.buffer == NULL) {
		nll_puterr(
"ERROR: allocating memory for 3D slowness grid buffer.");
		exit(EXIT_ERROR_MEMORY);
	}
	/* create array access pointers */
	mod_grid.array = CreateGridArray(&mod_grid);
	if (mod_grid.array == NULL) {
		nll_puterr(
"ERROR: creating array for accessing 3D vel/slowness grid buffer.");
		exit(EXIT_ERROR_MEMORY);
	}



	/* create model for each wave type */

	for (nWaveType = 0; nWaveType < NumWaveTypes; nWaveType++) {

		sprintf(fileRoot, "%s.mod", WaveType[nWaveType]);
		sprintf(MsgStr, "Creating model grid files: %s.%s.*",
			fn_vg_output, fileRoot);
		nll_putmsg(1, MsgStr);

		/* load vel model to grid */

		if ((istat = VelModToGrid3d(&mod_grid, WaveType[nWaveType])) < 0) {
			nll_puterr("ERROR: loading velocity model to grid.");
			exit(EXIT_ERROR_MODEL);
		}

		/* save grid to disk */

		if ((istat = WriteGrid3dBuf(&mod_grid, NULL, fn_vg_output, fileRoot)) < 0) {
			nll_puterr("ERROR: writing slowness grid to disk.");
			exit(EXIT_ERROR_IO);
		}

	}


	exit(EXIT_NORMAL);

}




/*** function to read input file */

int ReadVel2GridInput(FILE* fp_input)
{
	int istat, iscan;
	char param[MAXLINE], *pchr;
	char line[2*MAXLINE], *fgets_return;

	int flag_control = 0, flag_outfile = 0, flag_grid = 0, flag_type = 0,
		flag_trans = 0;
	int flag_include = 1;



	/* read each input line */

	while ((fgets_return = fgets(line, 2*MAXLINE, fp_input)) != NULL
			|| fp_include != NULL) {


		/* check for end of include file */

		if (fgets_return == NULL && fp_include != NULL) {
			SwapBackIncludeFP(&fp_input);
			continue;
		}


		istat = -1;

		/*read parmeter line */

		if ((iscan = sscanf(line, "%s", param)) < 0 )
			continue;

		/* skip comment line or white space */

		if (strncmp(param, "#", 1) == 0 || isspace(param[0]))
			istat = 0;



		/* read include file params and set input to include file */

		if (strcmp(param, "INCLUDE") == 0)
			if ((istat = GetIncludeFile(strchr(line, ' '), 
							&fp_input)) < 0) {
				nll_puterr("ERROR: processing include file.");
				flag_include = 0;
			}


		/* read control params */

		if (strcmp(param, "CONTROL") == 0) {
			if ((istat = get_control(strchr(line, ' '))) < 0) 
				nll_puterr("ERROR: reading control params.");
			else
				flag_control = 1;
		}


		/*read transform params */

		if (strcmp(param, "TRANS") == 0) {
    			if ((istat = get_transform(0, strchr(line, ' '))) < 0)
			    nll_puterr("ERROR: reading transformation parameters.");
			else
				flag_trans = 1;
		}


		/* read output file name (OUTfile) */

		if (strcmp(param, "VGOUT") == 0) {
			if ((istat = get_vg_outfile(strchr(line, ' '))) < 0)
				nll_puterr(
				"ERROR: reading Vel2Grid output file name.");
			else
				flag_outfile = 1;
		}

		/* read grid params */

		if (strcmp(param, "VGGRID") == 0) {
    			if ((istat = get_grid(strchr(line, ' '))) < 0)
				nll_puterr("ERROR: reading grid parameters.");
			else
				flag_grid = 1;
		}



		/* read grid type (VGTYPE) */

		if (strcmp(param, "VGTYPE") == 0) {
			if ((istat = get_vg_type(strchr(line, ' '))) < 0)
				nll_puterr("ERROR: reading Vel2Grid grid type.");
			else
				flag_type = 1;
		}

		/* check for velocity model input */

		istat = read_vel_mod_input(fp_input, param, 
			line, istat, message_flag >= 2);
			//strchr(line, ' '), istat, message_flag >= 2);



		/* unrecognized input */

		if (istat < 0) {
			if ((pchr = strchr(line, '\n')) != NULL)
				*pchr = '\0';
			sprintf(MsgStr, "Skipping input: %s", line);
			nll_putmsg(4, MsgStr);
		}

	}
	


	/* check for missing input */

	if (!flag_control) 
		nll_puterr("ERROR: no control (CONTROL) params read.");
	if (!flag_outfile) 
		nll_puterr("ERROR: no outputfile (VGOUT) params read.");
	if (!flag_type) 
		nll_puterr("ERROR: no type (VGTYPE) params read.");
	if (!flag_grid) 
		nll_puterr("ERROR: no grid (VGGRID) params read.");

	if (!flag_trans) {
		sprintf(MsgStr, "INFO: no transformation (TRANS) params read.");
		nll_putmsg(1, MsgStr);
		Hypocenter.comment[0] = '\0';
	}
	
	return (flag_include * flag_control * flag_outfile * flag_grid * 
			flag_type - 1);
}



/*** function to read output file name ***/

int get_vg_outfile(char* line1)
{

	sscanf(line1, "%s", fn_vg_output);

	sprintf(MsgStr, "Vel2Grid files:  Output: %s.*",
		 fn_vg_output);
	nll_putmsg(3, MsgStr);

	return(0);
}





/*** function to read output file name ***/

int get_vg_type(char* line1)		     
{

	if (NumWaveTypes >= MAX_NUM_WAVE_TYPES) {
		nll_puterr("WARNING: maximum number of wave types reached, ignoring wave type.");
		return(-1);
	}


	sscanf(line1, " %s", WaveType[NumWaveTypes]);

	sprintf(MsgStr, "Vel2Grid wave type:  %s", WaveType[NumWaveTypes]);
	nll_putmsg(3, MsgStr);

	NumWaveTypes++;


	return(0);
}





/*** function to load 3D velocity model to model grid ***/
/*  Notes:  

	(1)	Implicit staggered grids used -
	model space and travel time grid is numx X numy X numz while 
	slowness grid is numx-1 X numy-1 X numz-1.
	Slowness values are for mid-points of travel time grid cells,
	i.e. slowness grid is shifted (+dx/2,+dy/2,+dz/2) in space relative to 
	travel time grid.

	(2)	Podvin Lecomte FFD uses cubic cells, i.e. dx=dy=dz.

*/

#define DUMP_LAT_LON_TO_FILE 0

int VelModToGrid3d(GridDesc* grid, char *waveType)
{

	int ix, iy, iz;
	int imodel;
	char cWaveType;
	double xval, yval, xloc, yloc, zdepth;
	double vel, den, vel1;
	
	double dlat, dlon;
	FILE *fpfile;


	/* check wavetype */

	if (strcmp(waveType, "P") == 0)
		cWaveType = 'P';
	else if (strcmp(waveType, "S") == 0)
		cWaveType = 'S';
	else {
		nll_puterr2( "ERROR: unrecognized wave type", waveType);
		return(-1);
	}
	
	
	if (DUMP_LAT_LON_TO_FILE) {
		fpfile = fopen("lon_lat.txt", "w");
	}


	/* generate grid values */
	/* Note:  staggered grid assumed, thus vel lookup is shifted +dx/2, etc. */

	xval = grid->origx + grid->dx / 2.0;
	for (ix = 0; ix <  grid->numx; ix++) {
		
		yval = grid->origy + grid->dy / 2.0;
		for (iy = 0; iy <  grid->numy; iy++) {
			
			if (ModelCoordsMode == COORDS_LATLON) {
				rect2latlon(0, xval, yval, &yloc, &xloc);
/*printf("rect2latlon(0, xval %lf yval %lf yloc %lf xloc %lf\n", xval, yval, yloc, xloc);*/
			} else {
				xloc = xval;
				yloc = yval;
			}
			
			zdepth = grid->origz + grid->dz / 2.0;
			for (iz = 0; iz <  grid->numz; iz++) {
				
				if (DUMP_LAT_LON_TO_FILE) {
					rect2latlon(0, xval, yval, &dlat, &dlon);
					if (zdepth >= grid->dz / 2.0)
						fprintf(fpfile, "%f %f %f\n", dlon, dlat, -zdepth * 1000.0);
					else
						fprintf(fpfile, "%f %f %f\n", dlon, dlat, -1000.0 * grid->dz / 2.0);
				}


				/* check for non-lat/lon and non-layer 
							vel mod element */
				vel = get_vel(xval, yval, zdepth, 
					cWaveType, &den, 0, &imodel);
				if (imodel >= LAYEROFFSET || imodel < 0) {
					/* check for surface */
					vel1 = get_surface_vel(
						xloc, yloc, zdepth,
					    cWaveType, model_surface,
					    num_surfaces, &den, 0);
					vel = vel1 > 0.0 ? vel1 : vel;
/*if (ix == 5 && iy == 5)
printf("xloc %lf yloc %lf zdepth %lf cWaveType %c imodel %d vel %lg\n", xloc, yloc, zdepth, cWaveType, imodel, vel);*/
				}

				if (vel < 0.0) {
					nll_puterr("ERROR: cannot get velocity.");
					return(-1);
				}

				switch (grid->type) {

				case GRID_VELOCITY:
				    ((GRID_FLOAT_TYPE***) grid->array)[ix][iy][iz] = vel;
				    break;

				case GRID_VELOCITY_METERS:
				    ((GRID_FLOAT_TYPE***) grid->array)[ix][iy][iz] = 1000.0 * vel;
				    break;

				case GRID_SLOW_LEN:
				    ((GRID_FLOAT_TYPE***) grid->array)[ix][iy][iz] = grid->dx / vel;
				    break;

				case GRID_SLOW2_METERS:
				    ((GRID_FLOAT_TYPE***) grid->array)[ix][iy][iz] =
					(1.0e-3 / vel) * (1.0e-3 / vel);
				    break;

				default:
				    nll_puterr("ERROR: unrecognized grid type.");
					return(-1);

				}


				zdepth += grid->dz;
			}
			yval += grid->dy;
		}
		xval += grid->dx; 
	}

	if (DUMP_LAT_LON_TO_FILE) {
		fclose(fpfile);
	}


	return (0);

}




/*------------------------------------------------------------/ */
/* Anthony Lomax           | email: lomax@faille.unice.fr     / */
/* UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax / */
/* 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        / */
/* 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        / */
/*------------------------------------------------------------/ */

