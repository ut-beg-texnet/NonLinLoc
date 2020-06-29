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


/*   Sec2Grid.c

	Program to generate 3-D vel/slowness model grid from 2D velocity sections

*/

/*------------------------------------------------------------/ */
/* Anthony Lomax           | email: lomax@faille.unice.fr     / */
/* UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax / */
/* 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        / */
/* 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        / */
/*------------------------------------------------------------/ */


/*
	history:

	ver 01    18AUG1999  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

*/



#define PNAME  "Sec2Grid"

#include "../src/GridLib.h"
#include "ran1.h"
#include "velmod.h"

// Linux fix
#define xmalloc malloc


// do_spline allocations

#define DO_SPLINE_LEN 100
#define DO_SPLINE_OUT_LEN 1000
#define YDIMENSION 1
double *do_spline_t;
double **do_spline_y;
double **do_spline_z;
double *do_spline_x_out;
double **do_spline_yy_out;

#define bool int

inline void do_spline (int used, int len, double **t, int ydimension, double **y, double **z, 
	   double tension, bool periodic, bool spec_boundary_condition,
	   double k, int precision, double first_t, double last_t, 
	   double spacing_t, int no_of_intervals, bool spec_first_t, 
	   bool spec_last_t, bool spec_spacing_t, 
	   bool spec_no_of_intervals, bool suppress_abscissa,
	   double *do_spline_x, double **do_spline_yy
	);


/* structures  */

enum { OBLIQUE, X_CONST_NORTH, X_CONST_SOUTH, Y_CONST_EAST, Y_CONST_WEST };

typedef struct {

	struct surface grid2d;
	double x0, y0, z0;		/* coordinates of origin corner in 3D coord system */
	double azimuth;			/* azimuth of section CW from North */
	double scale2km;		/* scale to venvert grid length units to km */
	double x_orig_sect;		/* x coord on 2D sect to place at origin in 3D sys */

	double extrap_min, extrap_max;	/* distance to extrapolate ends of grid */

	double cosAz, sinAz;		/* trig values */
	double xmin, ymin, zmin;	/* section limits in 3D coord system */
	double xmax, ymax, zmax;	/* section limits in 3D coord system */
	double m_const, b_const;	/* const for eq of line through section in 3D coord system
						y = mx + b */
	int line_type; 		/* type of line OBLIQUE, ... */

} Sect2D;

//XXXX


typedef struct
{
 	double value; 		/* data value */
	double azimuth;		/* azimuth in 3D coord system */
}
SectValue;



/* velocity mapping parameters */

typedef struct
{
 	double vOld; 		/* original velocity */
	double vNew;		/* new velocity */
}
VelocityMapping;



/* isoline cutoff parameters */

typedef struct {

	double vCut;			/* cutoff velocity (km/sec) */
	double vReplaceRefLevel;	/* replacement velocity reference level*/
	double vReplace;		/* replacement velocity  (km/sec) at ref level */
	double vReplaceGrad;		/* replacement velocity gradient  (km/sec) */
	double depthForceReplace; 	/* depth below which to force vel replacement */
	int cutAtFirstMax;		/* start cutoff if max vel found */
	double velMinCutAtMax;		/* min vel to apply cutoff if max vel found */
	double velMax;			/* maximum velocity cutoff for replace */

} IsoCut;




/* globals  */

char fn_sg_output[MAXLINE];

double splineTension;	// spline tension parameter (see GNU spline doc)
double defaultValue;	// default value when no interpo value available
int nDefaultValues;
double valueCutoffMin;	// min cutoff value 
int nCutoffMinValues;
double valueCutoffMax;	// max cutoff value
int nCutoffMaxValues;
double interpDistMax;	// max distance along circumference from section to interpolate
int nInterpDistMax;
double radiusMin;	// min radius to apply spline interpolation (otherwise use mean)
int nRadiusMin;
int iExtrapolateEnd;	// extrapolate grid edge values if off grid

int splineError_NotEnoughNodes;
int splineError_OpeningFile;
int splineError_ReadingFile;

int spline_argc;
#define NUM_ARGS 20
char *spline_argv[NUM_ARGS];

SectValue *secValue;

/* wave type (P, S, ...) for vel grids */
#define MAX_NUM_WAVE_TYPES 10
char WaveType[MAX_NUM_WAVE_TYPES][12];
int NumWaveTypes;

double xCent, yCent;	/* center of polar coord system */


#define MAX_SECTIONS 100
Sect2D sections[MAX_SECTIONS];
Sect2D sections_depth[MAX_SECTIONS];
GridDesc  sections_depth_grid[MAX_SECTIONS];
int numSections;	/* number of section 2D grids */

#define MAX_VEL_MAPPINGS 100
VelocityMapping velMapping[MAX_VEL_MAPPINGS];
int numVelMapping;	/* number of velocity mappings */

#define MAX_ISOLINE_CUT 1
IsoCut isolineCut[MAX_ISOLINE_CUT];
int numIsolineCut;	/* number of isoline cutoffs */


/* function declarations */

void spline_main(int argc, char *argv[]);

int ReadSec2GridInput(FILE* );
int get_section2D(Sect2D *psur, int nsurface, char *input_line, int imessage);
int get_sg_spline(char* line1);	     
int SetLineParam(Sect2D* psection);
int SetLimits(Sect2D* psection);
inline double section2value3D(double xval, double yval, double zdepth, 
	Sect2D* psection, int numSections, double xCent, double yCent, 
	double radMin, double interpDMax, int iExtrapolateEnd);
int SecToGrid3dPolar(GridDesc* grid, Sect2D* psection, int numSections, 
	double xCent, double yCent, double radMin, double interpDMax, int iExtrapolateEnd,
	int no_negative, double defaultVal, int *pnDefaultVals,
	double valCutoffMin, int *pnCutoffMinVals,
	double valCutoffMax, int *pnCutoffMaxVals, GridDesc* depth_cutoff_grid, IsoCut *pisocut
	);
inline int getSectValueAtIntersection(SectValue secVal[2], double xCent, double yCent, 
	double radius, Sect2D *psection, double zdepth, int iExtrapolateEnd);
inline double getGrid2DValue(Sect2D *psection, double xval, double yval, double zdepth,
			int iExtrapolateEnd);
inline int findIntersectionLineAndCircle(double m, double b, 
	double xCent, double yCent, double radius, double xIntrsct[2], double yIntrsct[2]);
inline int calculateRoots(double a, double b, double c, double roots[]);

inline int SortSectValueAzimuth(SectValue* secValue, int num_val);
inline int CmpSectValueAzimuth(const SectValue *keyval, const SectValue *datum );

int applyIsolineCutoffUpwards(IsoCut *pisocut, GridDesc* grid, GridDesc* depth_grid);
int surfaceToGrid(struct surface *psurf, GridDesc* pgrid_desc, char *type_str);
inline double getValueBelowSection(double zdepth);

GridDesc* grid3D_to_gridXY(GridDesc* pgrid, char *type);
int Sect2DToSect1D(Sect2D *psect2d, Sect2D *psect1d);

int get_velocity_mapping(VelocityMapping *velMap, int *pnVelMap, char *line1);
int applyVelocityMapping(VelocityMapping *velMapping, int numVelMapping, GridDesc* grid);
double mapVelocity(double vel, VelocityMapping *velMapping, int numVelMapping);

int mapIsolineCutoff(IsoCut *pisocut, GridDesc* grid, char *fileroot);


/*** program to generate  3-D vel/slowness grid */


#define NARGS 2

int main(int argc, char *argv[])
{

	int istat, n, i;
	int nWaveType, nsec;
	int nisocut;
	char fileRoot[MAXLINE], secFileRoot[MAXLINE], gmtFileRoot[MAXLINE];

	GridDesc mod_grid;	/* model 3D grid */
	GridDesc tmp_2D_grid;	/* temp 2D grid */
	GridDesc* pdepth_grid2D;	/* depth 2D grid */

	int nx, nz, idummy;
	long ioffset;




	/* set program name */
	strcpy(prog_name, PNAME);

	/* check command line for correct usage */

	if (argc != NARGS) {
		disp_usage(prog_name, "<control file>");
		exit(EXIT_ERROR_USAGE);
	}



	do_spline_t = (double *)xmalloc (sizeof(double) * DO_SPLINE_LEN);
	do_spline_y = (double **)xmalloc (sizeof(double *) * YDIMENSION);
	do_spline_z = (double **)xmalloc (sizeof(double *) * YDIMENSION);
	do_spline_yy_out = (double **)xmalloc (sizeof(double *) * YDIMENSION);
	do_spline_x_out = (double *)xmalloc (sizeof(double) * DO_SPLINE_OUT_LEN);
	for (i = 0; i < YDIMENSION; i++) {
		do_spline_y[i] = (double *)xmalloc (sizeof(double) * DO_SPLINE_LEN);
		do_spline_z[i] = (double *)xmalloc (sizeof(double) * DO_SPLINE_LEN);
		do_spline_yy_out[i] = (double *)xmalloc (sizeof(double) * DO_SPLINE_OUT_LEN);
	}


	/* set constants */

	prog_mode_3d = 1;
	NumWaveTypes = 0;
	numIsolineCut = 0;
	numSections = 0;
	splineTension = 0.0;
	SetConstants();

	/* read control file */

	strcpy(fn_control, argv[1]);
	if ((fp_control = fopen(fn_control, "r")) == NULL) {
		nll_puterr("ERROR: opening control file.");
		exit(EXIT_ERROR_FILEIO);
	}

	if ((istat = ReadSec2GridInput(fp_control)) < 0) {
		exit(EXIT_ERROR_FILEIO);
	}
	mod_grid = grid_in;



	/* set spline function arguments */
	for (n = 0; n < NUM_ARGS; n++)
		spline_argv[n] = (char *) malloc((size_t) 100);
	spline_argc = 0;
	strcpy(spline_argv[spline_argc++], "spline_main");
	strcpy(spline_argv[spline_argc++], "-I");
	strcpy(spline_argv[spline_argc++], "d");
	strcpy(spline_argv[spline_argc++], "-O");
	strcpy(spline_argv[spline_argc++], "d");
	strcpy(spline_argv[spline_argc++], "-p");
	strcpy(spline_argv[spline_argc++], "-n");
	strcpy(spline_argv[spline_argc++], "360");
	strcpy(spline_argv[spline_argc++], "-T");
	sprintf(spline_argv[spline_argc++], "%lf", splineTension);
	strcpy(spline_argv[spline_argc++], "-o");
	strcpy(spline_argv[spline_argc++], "spline.out");
	strcpy(spline_argv[spline_argc++], "spline.in");



	/* initialize random number generator */
	SRAND_FUNC(RandomNumSeed);

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

	/* allocate secValue */
	secValue = (SectValue *) malloc((size_t) (2 * numSections) * sizeof(SectValue));

	/* calculate line segment parameterization of section in 3D coords */
	for (nsec = 0; nsec < numSections; nsec++) {
		SetLimits(sections + nsec);
		SetLineParam(sections + nsec);
	}


	/* apply isoline cutoffs to sections */
	for (nsec = 0; nsec < numSections; nsec++) {
		surfaceToGrid(&((sections + nsec)->grid2d), &tmp_2D_grid, "VELOCITY");
		// save file before applying cutoff
		sprintf(secFileRoot, "%s.raw", (sections + nsec)->grid2d.grd_file);
		WriteGrid3dBuf(&tmp_2D_grid, NULL, secFileRoot, "vel");
		// create 1D grid for storing cutoff depths
		Sect2DToSect1D(sections + nsec, sections_depth + nsec);
		surfaceToGrid(&((sections_depth + nsec)->grid2d), sections_depth_grid + nsec, "DEPTH");
		// apply cutoff
		applyIsolineCutoffUpwards(isolineCut, &tmp_2D_grid, sections_depth_grid + nsec);
		/* load 2D section grid back to section */
		for (nz = 0; nz < tmp_2D_grid.numz; nz++) {
			ioffset = nz * tmp_2D_grid.numx;
			for (nx = 0; nx < tmp_2D_grid.numx; nx++) {
				*((sections + nsec)->grid2d.zdata + ioffset++) 
					= tmp_2D_grid.array[nx][0][nz];
			}
		}
		/* save file after applying cutoff */
		sprintf(secFileRoot, "%s.cut", (sections + nsec)->grid2d.grd_file);
		WriteGrid3dBuf(&tmp_2D_grid, NULL, secFileRoot, "vel");
		FreeGrid(&tmp_2D_grid);
		DestroyGridArray(&tmp_2D_grid);

		/* load 1D depth grid back to section */
		for (nx = 0; nx < (sections_depth_grid + nsec)->numx; nx++) {
			*((sections_depth + nsec)->grid2d.zdata + nx) 
				= (sections_depth_grid + nsec)->array[nx][0][0];
		}
	}

	/* interpolate section depths to 2D depth grid */
	sprintf(fileRoot, "depth0");
	sprintf(MsgStr, "Creating 2D depth grid file: %s.%s.*", fn_sg_output, fileRoot);
	putmsg(1, MsgStr);
	pdepth_grid2D = grid3D_to_gridXY(&mod_grid, "VELOCITY");
	if ((istat = 
		SecToGrid3dPolar(pdepth_grid2D, sections_depth, numSections,
			xCent, yCent, radiusMin, interpDistMax, iExtrapolateEnd,
			0, 0.0, &idummy, -LARGE_DOUBLE, &idummy, 
			LARGE_DOUBLE, &idummy, NULL, NULL)) < 0) {
		nll_puterr("ERROR: interpolating section depths to 2D grid.");
		exit(EXIT_ERROR_MODEL);
	}
	/* save depth grid to disk */
	if ((istat = 
		WriteGrid3dBuf(pdepth_grid2D, NULL, fn_sg_output, fileRoot)) < 0) {
		nll_puterr("ERROR: writing section depth grid to disk.");
		exit(EXIT_ERROR_IO);
	}


	/* create model for each wave type */

	for (nWaveType = 0; nWaveType < NumWaveTypes; nWaveType++) {

		nDefaultValues = 0;
		nCutoffMinValues = 0;
		nCutoffMaxValues = 0;
		nInterpDistMax = 0;
		nRadiusMin = 0;
		splineError_NotEnoughNodes = 0;
		splineError_OpeningFile = 0;
		splineError_ReadingFile = 0;

		sprintf(fileRoot, "%s.mod", WaveType[nWaveType]);
		sprintf(MsgStr, "Creating model grid files: %s.%s.*",
			fn_sg_output, fileRoot);
		putmsg(1, MsgStr);



		/* load vel model to grid */

		if ((istat = 
			SecToGrid3dPolar(&mod_grid, sections, numSections,
				xCent, yCent, radiusMin, interpDistMax, iExtrapolateEnd,
				1, defaultValue, &nDefaultValues,
				valueCutoffMin, &nCutoffMinValues,
				valueCutoffMax, &nCutoffMaxValues, pdepth_grid2D, isolineCut)) < 0) {
			nll_puterr("ERROR: loading velocity model to grid.");
			exit(EXIT_ERROR_MODEL);
		}


		/* generate depth for isoline cutoffs 3D grid */
		sprintf(gmtFileRoot, "%s.%s", fn_sg_output, WaveType[nWaveType]);
		mapIsolineCutoff(isolineCut, &mod_grid, gmtFileRoot);


		/* apply velocity mapping to 3D grid */
		if (numVelMapping > 2) {
			applyVelocityMapping(velMapping, numVelMapping, &mod_grid);
		}


		/* save grid to disk */

		if ((istat = 
			WriteGrid3dBuf(&mod_grid, NULL, fn_sg_output, fileRoot))
				< 0) {
			nll_puterr("ERROR: writing slowness grid to disk.");
			exit(EXIT_ERROR_IO);
		}

		sprintf(MsgStr, 
			"Number of Values set to: Default (interp error) = %d, MinCutoffValue = %d, MaxCutoffValue = %d", 
			nDefaultValues, nCutoffMinValues, nCutoffMaxValues);
		putmsg(1, MsgStr);

		sprintf(MsgStr, 
			"Number of Values for: MaxInterpDist = %d, UsedMean = %d", 
			nInterpDistMax, nRadiusMin);
		putmsg(1, MsgStr);

		sprintf(MsgStr, 
			"Spline errors: NotEnoughNodes = %d, OpeningFile = %d, ReadingFile = %d", 
				splineError_NotEnoughNodes, splineError_OpeningFile, 
				splineError_ReadingFile);
		putmsg(1, MsgStr);

	}


	FreeGrid(pdepth_grid2D);
	DestroyGridArray(pdepth_grid2D);


	exit(EXIT_NORMAL);

}




/*** function to convert 2D surface to grid */

int surfaceToGrid(struct surface *pgrid2d, GridDesc* pgrid_desc, char *type_str)
{
	int nx, nz;
	long ioffset;

	/* setup grid description */

	pgrid_desc->numx = pgrid2d->hdr->nx;
	pgrid_desc->numy = 1;
	pgrid_desc->numz = pgrid2d->hdr->ny;
 	pgrid_desc->origx = pgrid2d->hdr->x_min;
	pgrid_desc->origy = 0.0;
	pgrid_desc->origz = pgrid2d->hdr->y_min;
 	pgrid_desc->autox = pgrid_desc->autoy = pgrid_desc->autoz = 0; 
	pgrid_desc->dx = pgrid2d->hdr->x_inc;
	pgrid_desc->dy = pgrid2d->hdr->x_inc;
	pgrid_desc->dz = pgrid2d->hdr->y_inc;
	strcpy(pgrid_desc->chr_type, type_str);
	pgrid_desc->type = convert_grid_type(pgrid_desc);
	strcpy(pgrid_desc->title, pgrid2d->hdr->title);
	pgrid_desc->sum = 0.0;


	/* allocate model grid */
	pgrid_desc->buffer = AllocateGrid(pgrid_desc);
	if (pgrid_desc->buffer == NULL) {
		nll_puterr("ERROR: allocating memory for 2D temp grid buffer.");
		exit(EXIT_ERROR_MEMORY);
	}
	/* create array access pointers */
	pgrid_desc->array = CreateGridArray(pgrid_desc);
	if (pgrid_desc->array == NULL) {
	nll_puterr("ERROR: creating array for accessing 2D temp grid buffer.");
		exit(EXIT_ERROR_MEMORY);
	}

	/* load section data to grid */
	for (nz = 0; nz < pgrid_desc->numz; nz++) {
		ioffset = nz * pgrid_desc->numx;
		for (nx = 0; nx < pgrid_desc->numx; nx++) {
			pgrid_desc->array[nx][0][nz] = *(pgrid2d->zdata + ioffset++);
		}
	}
}


/*** function to create 1D grid corresponding to 2D surface */

int Sect2DToSect1D(Sect2D *psect2d, Sect2D *psect1d)
{
	int nx;


	/* copy section fields */
	*psect1d = *psect2d;
	psect1d->grid2d.hdr = malloc(sizeof(struct grd_hdr));
	*(psect1d->grid2d.hdr) = *(psect2d->grid2d.hdr);
//XXXX

	/* set surface fields */
	psect1d->grid2d.hdr->ny = 1;
	psect1d->grid2d.hdr->y_min = psect1d->grid2d.hdr->y_max = 0.0;
	psect1d->grid2d.hdr->y_inc = LARGE_DOUBLE;
	psect1d->grid2d.hdr->z_scale_factor = 1.0;
	psect1d->grid2d.hdr->z_add_offset = 0.0;

	psect1d->grid2d.zdata = malloc(sizeof(double) * psect1d->grid2d.hdr->nx);
	for (nx = 0; nx < psect1d->grid2d.hdr->nx; nx++)
		*(psect1d->grid2d.zdata + nx) = 0.0;

	return(0);

}


/*** function to read input file */

int ReadSec2GridInput(FILE* fp_input)
{
	int istat, iscan;
	char param[MAXLINE], *pchr;
	char line[2*MAXLINE], *fgets_return;

	int flag_control = 0, flag_outfile = 0, flag_grid = 0, flag_type = 0,
		flag_trans = 0, flag_spline = 0, flag_sections = 0, flag_isocut = 0,
		flag_velmap = 0;
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

		if (strcmp(param, "CONTROL") == 0)
			if ((istat = get_control(strchr(line, ' '))) < 0) 
				nll_puterr("ERROR: reading control params.");
			else
				flag_control = 1;


		/*read transform params */

		if (strcmp(param, "TRANS") == 0)
    			if ((istat = get_transform(strchr(line, ' '))) < 0)
			    nll_puterr("ERROR: reading transformation parameters.");
			else
				flag_trans = 1;


		/* read output file name (OUTfile) */

		if (strcmp(param, "SGOUT") == 0)
			if ((istat = get_sg_outfile(strchr(line, ' '))) < 0) 
				nll_puterr(
				"ERROR: reading Sec2Grid output file name.");
			else
				flag_outfile = 1;

		/* read grid params */

		if (strcmp(param, "SGGRID") == 0)
    			if ((istat = get_grid(strchr(line, ' '))) < 0)
				nll_puterr("ERROR: reading grid parameters.");
			else
				flag_grid = 1;



		/* read grid type */

		if (strcmp(param, "SGTYPE") == 0)
			if ((istat = get_sg_type(strchr(line, ' '))) < 0) 
				nll_puterr("ERROR: reading Sec2Grid grid type.");
			else
				flag_type = 1;


		/* read spline params */

		if (strcmp(param, "SGSPLINE") == 0)
			if ((istat = get_sg_spline(strchr(line, ' '))) < 0) 
				nll_puterr("ERROR: reading Sec2Grid spline params.");
			else
				flag_spline = 1;


		/* read spline params */

		if (strcmp(param, "SGDISTAZ") == 0)
			if ((istat = get_sg_distaz(strchr(line, ' '))) < 0) 
				nll_puterr("ERROR: reading Sec2Grid distance/azimuth params.");


		/* read 2D section desc */

		if (strcmp(param, "SG2DSECT") == 0) {	
    			if ((istat = get_section2D(sections, numSections, 
					strchr(line, ' '), message_flag > 1)) < 0)
				fprintf(stderr, "ERROR: reading 2D section.\n");
			else {
				numSections++;
				flag_sections = 1;
			}
		}


		/* read velocity mapping params */

		if (strcmp(param, "SGVELMAP") == 0) {	
    			if ((istat = get_velocity_mapping(velMapping, &numVelMapping, 
					strchr(line, ' '))) < 0)
				fprintf(stderr, "ERROR: reading velocity mapping params.\n");
			else {
				flag_velmap = 1;
			}
		}


		/* read isoline cutoff params */

		if (strcmp(param, "SGISOCUT") == 0) {	
    			if ((istat = get_isoline_cut(isolineCut, numIsolineCut, 
					strchr(line, ' '), message_flag > 1)) < 0)
				fprintf(stderr, "ERROR: reading isoline cutoff params.\n");
			else {
				numIsolineCut++;
				flag_isocut = 1;
			}
		}


		/* check for velocity model input */

/*!!!! REMOVE	istat = read_vel_mod_input(fp_input, param, 
			line, istat, message_flag >= 2);
			//strchr(line, ' '), istat, message_flag >= 2);
*/


		/* unrecognized input */

		if (istat < 0) {
			if ((pchr = strchr(line, '\n')) != NULL)
				*pchr = '\0';
			sprintf(MsgStr, "Skipping input: %s", line);
			putmsg(4, MsgStr);
		}

	}
	


	/* check for missing input */

	if (!flag_control) 
		nll_puterr("ERROR: no control (CONTROL) params read.");
	if (!flag_outfile) 
		nll_puterr("ERROR: no outputfile (SGOUT) params read.");
	if (!flag_type) 
		nll_puterr("ERROR: no type (SGTYPE) params read.");
	if (!flag_spline) 
		nll_puterr("ERROR: no spline (SGSPLINE) params read.");
	if (!flag_grid) 
		nll_puterr("ERROR: no grid (SGGRID) params read.");
	if (!flag_sections) 
		nll_puterr("ERROR: no 2D section (SG2DSECT) params read.");

	if (!flag_velmap) {
		sprintf(MsgStr, "INFO: no velocity mapping (SGVELMAP) params read.");
		putmsg(1, MsgStr);
	}

	if (!flag_isocut) {
		sprintf(MsgStr, "INFO: no isoline cutoff (SGISOCUT) params read.");
		putmsg(1, MsgStr);
	}

	if (!flag_trans) {
		sprintf(MsgStr, "INFO: no transformation (TRANS) params read.");
		putmsg(1, MsgStr);
		Hypocenter.comment[0] = '\0';
	}
	
	return (flag_include * flag_control * flag_outfile * flag_grid * 
			flag_type * flag_spline * flag_sections - 1);
}



/*** function to read output file name ***/

int get_sg_outfile(char* line1)		     
{

	int istat;

	if ((istat = sscanf(line1, "%s", fn_sg_output)) != 1)
		return(-1);


	sprintf(MsgStr, "Sec2Grid files:  Output: %s.*",
		 fn_sg_output);
	putmsg(1, MsgStr);

	return(0);
}





/*** function to read wave type ***/

int get_sg_spline(char* line1)		     
{
	int istat;
	char coord_type[MAXLINE];
	double dlat, dlong;


	if ((istat = sscanf(line1, "%s %lf %lf %lf %lf %lf %lf %lf %lf %d", coord_type,
			&xCent, &yCent, &splineTension, 
			&valueCutoffMin, &valueCutoffMax, &defaultValue,
			&interpDistMax, &radiusMin, &iExtrapolateEnd)) != 10)
		return(-1);


	if (strncmp(coord_type, "XYZ", 3) == 0){
		;
	}
	else if (strcmp(coord_type, "LATLON") == 0) {
		dlat = xCent;
		dlong = yCent;
		istat = latlon2rect(dlat, dlong, &xCent, &yCent);
	}


	sprintf(MsgStr, "Sec2Grid spline:  center: (%lf,%lf)  tension: %lf  value limits: %lf -> %lf  value default: %lf  interpDistMax: %lf  radiusMin: %lf  extrapolateEnds: %d", 
		xCent, yCent, splineTension, valueCutoffMin, valueCutoffMax, defaultValue,
		interpDistMax, radiusMin, iExtrapolateEnd);
	putmsg(1, MsgStr);

	return(0);
}



/*** function to read and calculate distance/azimuths between 2 points ***/

int get_sg_distaz(char* line1)		     
{
	int istat;
	char name1[MAXLINE], coord_type1[MAXLINE];
	char name2[MAXLINE], coord_type2[MAXLINE];
	double x1, y1, x2, y2;
	double dlat, dlong;
	double dist, az;


	if ((istat = sscanf(line1, "%s %s %lf %lf %s %s %lf %lf", 
			name1, coord_type1, &x1, &y1,
			name2, coord_type2, &x2, &y2
			)) != 8)
		return(-1);


	if (strncmp(coord_type1, "XYZ", 3) == 0){
		;
	}
	else if (strcmp(coord_type1, "LATLON") == 0) {
		dlat = x1;
		dlong = y1;
		istat = latlon2rect(dlat, dlong, &x1, &y1);
	}
	if (strncmp(coord_type2, "XYZ", 3) == 0){
		;
	}
	else if (strcmp(coord_type2, "LATLON") == 0) {
		dlat = x2;
		dlong = y2;
		istat = latlon2rect(dlat, dlong, &x2, &y2);
	}

	dist = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
	az = atan2(x2 - x1, y2 - y1) / rpd;
	if (az < 0.0)
		az += 360.0;


	sprintf(MsgStr, 
		"Distance/azimuth calc:  %s to %s  distance: %lf km  azimuth: %lf deg", 
		name1, name2, dist, az);
	putmsg(1, MsgStr);
	sprintf(MsgStr, 
		"   %s is at (%lf,%lf),  %s is at (%lf,%lf)", 
		name1, x1, y1, name2, x2, y2);
	putmsg(1, MsgStr);

	return(0);
}



/*** function to read limits ***/

int get_sg_type(char* line1)		     
{

	int istat;

	if (NumWaveTypes >= MAX_NUM_WAVE_TYPES) {
		nll_puterr("WARNING: maximum number of wave types reached, ignoring wave type.");
		return(-1);
	}


	if ((istat = sscanf(line1, " %s", WaveType[NumWaveTypes])) != 1)
		return(-1);


	sprintf(MsgStr, "Sec2Grid wave type:  %s", WaveType[NumWaveTypes]);
	putmsg(1, MsgStr);

	NumWaveTypes++;


	return(0);
}





/*** function to read isoline cutoff params */

/*	comments

*/

int get_isoline_cut(IsoCut *pcut, int nisocut, char *input_line, int imessage)
{

	int istat;
	IsoCut *ps;

	if (nisocut >= MAX_ISOLINE_CUT) {
		fprintf(stderr, "ERROR: Maximum number of isoline cutoffs read: %d\n", nisocut);
		return(-1);
	}

	ps = pcut + nisocut;

	if ((istat = sscanf(input_line, "%lf %lf %lf %lf %lf %d %lf %lf",
			&ps->vCut, 
			&ps->vReplaceRefLevel, &ps->vReplace, &ps->vReplaceGrad, 
			&ps->depthForceReplace, &ps->cutAtFirstMax, &ps->velMinCutAtMax, &ps->velMax))
				!= 8)
		return(-1);

	sprintf(MsgStr, "Sec2Grid isoline cutoff:  vCut: %lf  vReplaceRefLevel: %lf  vReplace: %lf  vReplaceGrad: %lf  depthForceReplace: %lf  cutAtFirstMax: %d  velMinCutAtMax: %lf  velMax: %lf", 
		ps->vCut, ps->vReplaceRefLevel, ps->vReplace, 
		ps->vReplaceGrad, ps->depthForceReplace, ps->cutAtFirstMax, ps->velMinCutAtMax, ps->velMax);
	putmsg(1, MsgStr);

	return (1);
}





/*** function to read model surface from input file */

/*	surface is defined by a GMT grd file
	surface array is in order of reading of surfaces
	velocity at point x,y,z is determined by first surface found, 
		working backwards from end of surface array, that is above point
*/

int get_section2D(Sect2D *psur, int nsurface, char *input_line, int imessage)
{

	int istat;
	char coord_type[MAXLINE];
	char str_ref_type[MAXLINE];
	Sect2D *ps;
	double dlat, dlong;


	if (nsurface >= MAX_SECTIONS) {
		fprintf(stderr, "ERROR: Maximum number of sections read: %d\n", nsurface);
		return(-1);
	}

	ps = psur + nsurface;

	if ((istat = sscanf(input_line, "%s %s %lf %lf %lf %lf %lf %lf %lf %lf",
			ps->grid2d.grd_file, coord_type,
			&ps->x0, &ps->y0, &ps->z0,
			&ps->azimuth, &ps->scale2km, &ps->x_orig_sect,
			&ps->extrap_min, &ps->extrap_max))
				!= 10)
		return(-1);

	if (strncmp(coord_type, "XYZ", 3) == 0){
		;
	}
	else if (strcmp(coord_type, "LATLON") == 0) {
		dlat = ps->x0;
		dlong = ps->y0;
		istat = latlon2rect(dlat, dlong, &(ps->x0), &(ps->y0));
	}

	if (read_grd(&(ps->grid2d), imessage) < 0)
		return(-1);

	if (convGridTokm(&(ps->grid2d), imessage) < 0)
		return(-1);

	return (1);
}



int get_velocity_mapping(VelocityMapping *velMap, int *pnVelMap, char *line1)
{

	int istat, nvelmap;
	double vel1, vel2;
	char frmt[MAXLINE] = "%lf %lf";
	char frmttmp[MAXLINE];

	*pnVelMap = 0;

	while (*pnVelMap < MAX_VEL_MAPPINGS && (istat = sscanf(line1, frmt, &vel1, &vel2)) == 2) {
		velMapping[*pnVelMap].vOld =  vel1;
		velMapping[*pnVelMap].vNew =  vel2;
		(*pnVelMap)++;
		sprintf(frmttmp, "%%*lf %%*lf %s", frmt);
		strcpy(frmt, frmttmp);
	}


	sprintf(MsgStr, "Sec2Grid SGVELMAP: 946");
	putmsg(2, MsgStr);
	for (nvelmap = 0; nvelmap < *pnVelMap; nvelmap++) {
		sprintf(MsgStr, " %lf->%lf ", 
			velMapping[nvelmap].vOld, velMapping[nvelmap].vNew);
		putmsg(1, MsgStr);
	}

	return(0);
}



/*** function to calculate line parameterization of section in 3D coords */

int SetLineParam(Sect2D* psection)
{

	// section parallel to Y axis
	if (fabs(psection->azimuth) < SMALL_DOUBLE 
				|| fabs(psection->azimuth - 180.0) < SMALL_DOUBLE) {
		/* line parallel to Y axis */
		psection->m_const = VERY_LARGE_DOUBLE;
		psection->b_const = VERY_LARGE_DOUBLE;
		if (fabs(psection->azimuth) < SMALL_DOUBLE)
			psection->line_type = X_CONST_NORTH;
		else
			psection->line_type = X_CONST_SOUTH;
		sprintf(MsgStr, "WARNING: Section is parallel to Y axis, will be ingnored: %s\n", 
			psection->grid2d.grd_file);
		nll_puterr(MsgStr);
		return(psection->line_type);
	}

	// section parallel to X axis
	if (fabs(psection->azimuth - 90.0) < SMALL_DOUBLE 
				|| fabs(psection->azimuth - 270.0) < SMALL_DOUBLE) {
		/* line parallel to X axis */
		psection->m_const = 0.0;
		psection->b_const = psection->y0;
		if (fabs(psection->azimuth) - 90.0 < SMALL_DOUBLE)
			psection->line_type = Y_CONST_EAST;
		else
			psection->line_type = Y_CONST_WEST;
//printf(" Y_CONST line: az %lf, y = %f x + %lf\n", psection->azimuth, psection->m_const, psection->b_const);
		return(psection->line_type);
	}

	// oblique section

	// get slope
	psection->m_const = cos(psection->azimuth * rpd) 
			/ sin(psection->azimuth * rpd);  // azimuth deg CW from N
	psection->b_const = psection->y0 - psection->m_const * psection->x0;
//printf("OBLIQUE line: az %lf, y = %f x + %lf\n", psection->azimuth, psection->m_const, psection->b_const);

	psection->line_type = OBLIQUE;
	return(psection->line_type);

}



/*** function to calculate limits section in 3D coords and set trig values */

int SetLimits(Sect2D* psection)
{

	char fn_mapline[MAXLINE];
	FILE *fp_mapline;

	double xmax, ymax, zmax;
	double xmin, ymin, zmin;

	int istep;
	double xloc, yloc;
	double xpos, ypos, xstep, ystep;


	// set trig constants
	psection->cosAz = cos(psection->azimuth * rpd);
	psection->sinAz = sin(psection->azimuth * rpd);

	// set limits    (azimuth deg CW from N; 2D->3D x->x, y->z)
	xmin = psection->x0 - psection->sinAz * psection->x_orig_sect;
	xmax = psection->x0 + psection->sinAz * 
		((double) (psection->grid2d.hdr->nx - 1) 
			* psection->grid2d.hdr->x_inc - psection->x_orig_sect);
	ymin = psection->y0 - psection->cosAz * psection->x_orig_sect;
	ymax = psection->y0 + psection->cosAz * 
		((double) (psection->grid2d.hdr->nx - 1) 
			* psection->grid2d.hdr->x_inc - psection->x_orig_sect);
	zmin = psection->z0;
	zmax = psection->z0 + 
		((double) (psection->grid2d.hdr->ny - 1) 
			* psection->grid2d.hdr->y_inc);

	/* write limits to plot file */



	sprintf(fn_mapline, "%s.limits.xy", psection->grid2d.grd_file);
	if ((fp_mapline = fopen(fn_mapline, "w")) == NULL) {
		nll_puterr("ERROR: opening mapline file.");
	} else {
		fprintf(fp_mapline, "> GMT_LONLAT\n");
		xpos = xmin;
		ypos = ymin;
		xstep = (xmax - xmin) / (double) 99;
		ystep = (ymax - ymin) / (double) 99;
		for (istep = 0; istep < 100; istep++) {
			xpos += xstep;
			ypos += ystep;
			rect2latlon(xpos, ypos, &yloc, &xloc);
			fprintf(fp_mapline, "%lf %lf\n", xloc, yloc);
		}
		fprintf(fp_mapline, ">\n");
		fclose(fp_mapline);
	}


	// sort
	psection->xmin = xmin < xmax ? xmin : xmax;
	psection->xmax = xmin < xmax ? xmax : xmin;
	psection->ymin = ymin < ymax ? ymin : ymax;
	psection->ymax = ymin < ymax ? ymax : ymin;
	psection->zmin = zmin < zmax ? zmin : zmax;
	psection->zmax = zmin < zmax ? zmax : zmin;

//printf("LIMITS line: cosAz %lf  sinAz %f  x(%lf,%lf)  y(%lf,%lf)\n", psection->cosAz, psection->sinAz, psection->xmin, psection->xmax, psection->ymin, psection->ymax, psection->zmin, psection->zmax);

	return(1);

}



/*** function to calculate real roots of a 2nd deg polynomial */

inline int calculateRoots(double a, double b, double c, double root[])
{

	double disciminant, sqrt_disc;

	root[0] = root[1] = -VERY_LARGE_DOUBLE;

	// set trig constants
	disciminant = (b * b) - 4 * a * c;
//printf("calculateRoots: disciminant = %lf \n", disciminant);

	// no real roots
	if (disciminant < 0.0) {
		return(0);

	// two real roots
	} else if (disciminant > 0.0) {
		sqrt_disc = sqrt(disciminant);
		root[0] = (-b + sqrt_disc) / (2 * a);
		root[1] = (-b - sqrt_disc) / (2 * a);
		return(2);

	// one real roots
	} else {
		root[0] = root[1] = -b / (2 * a);
		return(1);
	}


}



/*** function to apply isoline cutoff to model grid ***/
/*  Notes:  

	(1)	xxxx

*/

#define THRESHOLD 0
#define PEAK 1

int applyIsolineCutoffUpwards(IsoCut *pisocut, GridDesc* grid, GridDesc* depth_grid)
{

	int n, ix, iy, iz;
	int n_cut, iz_cut, iz_cut_test[10], type[10];
	double zdepth, zdepth_cut, zdepth_cut_test[10];
	double zdepth_cut_last, depth_diff, depth_diff_closest;
	double vel;
	double vel_last;
	int isIncreasing, iBeginCheck;



	// set zdepth_cut_last so that shallowest cut will be chosen first
	zdepth_cut_last = grid->origz;

	/* find isoline level in each grid colums */
	for (ix = 0; ix <  grid->numx; ix++) {

		// set zdepth_cut_last so that shallowest cut will be chosen first
		if (grid->numy > 1)
			zdepth_cut_last = grid->origz;

		for (iy = 0; iy <  grid->numy; iy++) {

			// find isoline level
			zdepth = grid->origz + grid->dz / 2.0 
					+ (double) (grid->numz - 1) * grid->dz;
			// save default max depth to depth grid
			depth_grid->array[ix][iy][0] = (float) zdepth;
			vel_last = -VERY_LARGE_DOUBLE;
			isIncreasing = 0;
			iBeginCheck = 0;
			n_cut = 0;
			for (iz = grid->numz - 1; iz >=0; iz--) {

				// get velocity at grid node

				switch (grid->type) {
					case GRID_VELOCITY:
					    vel = grid->array[ix][iy][iz];
					    break;
					case GRID_VELOCITY_METERS:
					    vel = grid->array[ix][iy][iz] / 1000.0;
					    break;
					case GRID_SLOW_LEN:
					    vel = grid->dx / grid->array[ix][iy][iz];
					    break;
					case GRID_SLOW2_METERS:
					    vel = 1.0e-3 / sqrt(grid->array[ix][iy][iz]);
					    break;
					default:
					     nll_puterr2("ERROR: unrecognized grid type", 
							grid->chr_type);
						return(-1);
				}

				// check for cutoff

				if (vel <= pisocut->vCut && vel_last > pisocut->vCut) {
					type[n_cut] = THRESHOLD;
					iz_cut_test[n_cut] = iz;
					zdepth_cut_test[n_cut++] = zdepth;
				}
				else if (pisocut->cutAtFirstMax && isIncreasing
							&& vel >= pisocut->velMinCutAtMax
							&& vel < vel_last) {
					type[n_cut] = PEAK;
					iz_cut_test[n_cut] = iz;
					zdepth_cut_test[n_cut++] = zdepth;
				}

				isIncreasing = iBeginCheck && vel >= vel_last;
				// avoid starting cutoff while below model
				if (fabs(vel - defaultValue) > 1.0e-6) {
					iBeginCheck = 1;
				}

				vel_last = vel;
				zdepth -= grid->dz;

			}


			// choose cut depth closest to last cut
			depth_diff_closest = VERY_LARGE_DOUBLE;
			for (n = n_cut - 1; n >= 0; n--) {
				if (type[n] == THRESHOLD) {
					iz_cut = iz_cut_test[n];
					zdepth_cut = zdepth_cut_test[n];
					depth_diff_closest = depth_diff;
					break;
				}
				depth_diff = fabs(zdepth_cut_test[n] - zdepth_cut_last);
				if (depth_diff < depth_diff_closest) {
					iz_cut = iz_cut_test[n];
					zdepth_cut = zdepth_cut_test[n];
					depth_diff_closest = depth_diff;
				}
			}

			if (n_cut == 0 || zdepth_cut > pisocut->depthForceReplace) {
				iz_cut = (pisocut->depthForceReplace 
					- grid->origz - grid->dz / 2.0) / grid->dz;
				zdepth_cut = pisocut->depthForceReplace;
			}

			zdepth_cut_last = zdepth_cut;

			if (iz_cut >= 0) {  // found - replace velocities

				zdepth = grid->origz + grid->dz / 2.0 
					+ (double) iz_cut * grid->dz;
				// save depth to depth grid
				depth_grid->array[ix][iy][0] = (float) zdepth;

				for (iz = iz_cut; iz <  grid->numz; iz++) {

					vel = pisocut->vReplace 
						+ (zdepth - pisocut->vReplaceRefLevel) 
							* pisocut->vReplaceGrad;
					vel = vel > pisocut->velMax ?  pisocut->velMax : vel;
					switch (grid->type) {
						case GRID_VELOCITY:
						    grid->array[ix][iy][iz] = vel;
						    break;
						case GRID_VELOCITY_METERS:
						    grid->array[ix][iy][iz] = 1000.0 * vel;
						    break;
						case GRID_SLOW_LEN:
						    grid->array[ix][iy][iz] = grid->dx / vel;
						    break;
						case GRID_SLOW2_METERS:
						    grid->array[ix][iy][iz] = 
							(1.0e-3 / vel) * (1.0e-3 / vel);
						    break;
						default:
						     nll_puterr2("ERROR: unrecognized grid type", 
								grid->chr_type);
							return(-1);
					}

					zdepth += grid->dz;
				}

			}
		}
	}

	return (0);

}



/*** function to create and allocate XY 2D grid corresponding to an XYZ 3D grid */

GridDesc* grid3D_to_gridXY(GridDesc* pgrid, char *type)
{

	GridDesc *pdepth_grid;

	pdepth_grid = (GridDesc *) malloc(sizeof(GridDesc));


	/* copy-initialize grid */
	*pdepth_grid = *pgrid;
	pdepth_grid->numz = 1;
	pdepth_grid->origz = 0.0;
	strcpy(pdepth_grid->chr_type, type);
	convert_grid_type(pdepth_grid);


	/* allocate grid */
	pdepth_grid->buffer = AllocateGrid(pdepth_grid);
	if (pdepth_grid->buffer == NULL) {
		nll_puterr(
"ERROR: allocating memory for isoline depth grid buffer.");
		return(NULL);
	}
	/* create array access pointers */
	pdepth_grid->array = CreateGridArray(pdepth_grid);
	if (pdepth_grid->array == NULL) {
		nll_puterr(
"ERROR: creating array for accessing isoline depth grid buffer.");
		return(NULL);
	}


	return(pdepth_grid);

}



/*** function to generate map of isoline cutoff surface ***/

int mapIsolineCutoff(IsoCut *pisocut, GridDesc* grid, char *fileroot)
{

	int ix, iy, iz;
	int iSearching;
	double zdepth;
	double vel, velrep;
	double vel_last;

	GridDesc* pdepth_grid;


	/* initialize 2D grid */
	pdepth_grid = grid3D_to_gridXY(grid, "DEPTH");


	/* find isoline replacement level in each grid colums */

	for (ix = 0; ix <  grid->numx; ix++) {
		for (iy = 0; iy <  grid->numy; iy++) {
			zdepth = grid->origz + grid->dz / 2.0;
			for (iz = 0; iz <  grid->numz; iz++) {

				/* get velocity at grid node */

				switch (grid->type) {
					case GRID_VELOCITY:
					    vel = grid->array[ix][iy][iz];
					    break;
					case GRID_VELOCITY_METERS:
					    vel = grid->array[ix][iy][iz] / 1000.0;
					    break;
					case GRID_SLOW_LEN:
					    vel = grid->dx / grid->array[ix][iy][iz];
					    break;
					case GRID_SLOW2_METERS:
					    vel = 1.0e-3 / sqrt(grid->array[ix][iy][iz]);
					    break;
					default:
					     nll_puterr2("ERROR: unrecognized grid type", 
							grid->chr_type);
						return(-1);
				}

	
				// check if exceeds replacement velocity

				velrep = pisocut->vReplace 
						+ (zdepth - pisocut->vReplaceRefLevel) 
							* pisocut->vReplaceGrad;
				velrep = velrep > pisocut->velMax ?  pisocut->velMax : velrep;

				if (vel >= velrep) {
					pdepth_grid->array[ix][iy][0] = (float) zdepth;
					break;
				}

				zdepth += grid->dz;
			}

			if (iz == grid->numz) {
				/* not found */
				pdepth_grid->array[ix][iy][0] = (float) zdepth;
			}


		}
	}



	if (WriteGrid3dBuf(pdepth_grid, NULL, fileroot, "depth") < 0) {
		nll_puterr("ERROR: writing isoline depth grid to disk.");
		return(-1);
	}


	FreeGrid(pdepth_grid);
	DestroyGridArray(pdepth_grid);


fprintf(stderr, "\nMapIsolineCutoff - Done looping through 3D grid!\n");
	return (0);

}




/*** function to apply velocity mapping to model grid ***/

int applyVelocityMapping(VelocityMapping *velMapping, int numVelMapping, GridDesc* grid)
{

	int ix, iy, iz;
	double zdepth;
	double vel;
	double vel_last;


	/* loop through grid */
	for (ix = 0; ix <  grid->numx; ix++) {
		for (iy = 0; iy <  grid->numy; iy++) {
			for (iz = 0; iz <  grid->numz; iz++) {

					// get velocity at grid node

					switch (grid->type) {
						case GRID_VELOCITY:
						    vel = grid->array[ix][iy][iz];
						    break;
						case GRID_VELOCITY_METERS:
						    vel = grid->array[ix][iy][iz] / 1000.0;
						    break;
						case GRID_SLOW_LEN:
						    vel = grid->dx / grid->array[ix][iy][iz];
						    break;
						case GRID_SLOW2_METERS:
						    vel = 1.0e-3 / sqrt(grid->array[ix][iy][iz]);
						    break;
						default:
						     nll_puterr2("ERROR: unrecognized grid type", 
								grid->chr_type);
							return(-1);
					}

					// replace velocity

					vel = mapVelocity(vel, velMapping, numVelMapping);

					switch (grid->type) {
						case GRID_VELOCITY:
						    grid->array[ix][iy][iz] = vel;
						    break;
						case GRID_VELOCITY_METERS:
						    grid->array[ix][iy][iz] = 1000.0 * vel;
						    break;
						case GRID_SLOW_LEN:
						    grid->array[ix][iy][iz] = grid->dx / vel;
						    break;
						case GRID_SLOW2_METERS:
						    grid->array[ix][iy][iz] = 
							(1.0e-3 / vel) * (1.0e-3 / vel);
						    break;
						default:
						     nll_puterr2("ERROR: unrecognized grid type", 
								grid->chr_type);
							return(-1);
					}
	

			}
		}
	}


fprintf(stderr, "\nVelocity mapping - Done looping through 3D grid!\n");
	return (0);

}



double mapVelocity(double vel, VelocityMapping *velMapping, int numVelMapping)
{
	int n1, n2;
	double newVel;


	// must be at least 2 mappings
	if (numVelMapping < 2)
		return(vel);

	// vel must be > than first mapping
	if (vel < velMapping[0].vOld)
		return(vel);

	// find relevant mapping interval
	n2 = 1;
	while (n2 < numVelMapping && velMapping[n2].vOld <= vel)
		n2++;
	if (n2 == numVelMapping)
		return(vel);
	n1 = n2 - 1;

	// calculate mapped velocity
	if (velMapping[n2].vOld - velMapping[n1].vOld > SMALL_FLOAT)
		newVel = velMapping[n1].vNew + (vel - velMapping[n1].vOld) 
			* (velMapping[n2].vNew - velMapping[n1].vNew)
			/ (velMapping[n2].vOld - velMapping[n1].vOld);
	else
		newVel = velMapping[n1].vNew; 

	return(newVel);


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

int SecToGrid3dPolar(GridDesc* grid, Sect2D* psection, int numSections, 
	double xCent, double yCent, double radMin, double interpDMax, int iExtrapolateEnd,
	int no_negative, double defaultVal, int *pnDefaultVals,
	double valCutoffMin, int *pnCutoffMinVals,
	double valCutoffMax, int *pnCutoffMaxVals, 
	GridDesc* depth_cutoff_grid, IsoCut *pisocut)
{

	int ix, iy, iz;
	int imodel;
	char cWaveType;
	double xval, yval, xloc, yloc, zdepth;
	double vel, den, vel1;



	/* generate grid values */
	/* Note:  staggered grid assumed, thus vel lookup is shifted +dx/2, etc. */
	xval = grid->origx + grid->dx / 2.0;
	for (ix = 0; ix <  grid->numx; ix++) {
fprintf(stderr, "-> ix = %d, iy = 0, iz = 0        \r", ix);
		yval = grid->origy + grid->dy / 2.0;
		for (iy = 0; iy <  grid->numy; iy++) {

if (grid->numx == 1)
fprintf(stderr, "-> ix = %d, iy = %d, iz = 0        \r", ix, iy);

			if (ModelCoordsMode == COORDS_LATLON) {
				rect2latlon(xval, yval, &yloc, &xloc);
/*printf("rect2latlon(xval %lf yval %lf yloc %lf xloc %lf\n", xval, yval, yloc, xloc);*/
			} else {
				xloc = xval;
				yloc = yval;
			}
			zdepth = grid->origz + grid->dz / 2.0;
			for (iz = 0; iz <  grid->numz; iz++) {
//printf("xval %lf yval %lf zdepth %lf psection %ld numSections %d xCent  %lf yCent %lf\n", xval, yval, zdepth, psection, numSections, xCent, yCent);

				/* check for depth below cutoff */
				if (depth_cutoff_grid != NULL 
						&& zdepth >= depth_cutoff_grid->array[ix][iy][0]) { 
					vel= pisocut->vReplace 
						+ (zdepth - pisocut->vReplaceRefLevel) 
							* pisocut->vReplaceGrad;
					vel = vel > pisocut->velMax ?  pisocut->velMax : vel;
				} else {
					vel = section2value3D(xval, yval, zdepth, 
						psection, numSections, xCent, yCent, 
						radMin, interpDMax, iExtrapolateEnd);
				}
				if (no_negative && vel < 0.0) {
					vel = defaultVal;
					(*pnDefaultVals)++;
				} else if (vel < valCutoffMin) {
					vel = valCutoffMin;
					(*pnCutoffMinVals)++;
				} else if (vel > valCutoffMax) {
					vel = valCutoffMax;
					(*pnCutoffMaxVals)++;
				}

//printf("xval %lf yval %lf zdepth %lf vel %lf\n", xval, yval, zdepth, vel);

				switch (grid->type) {
					case GRID_VELOCITY:
					    grid->array[ix][iy][iz] = vel;
					    break;
					case GRID_VELOCITY_METERS:
					    grid->array[ix][iy][iz] = 1000.0 * vel;
					    break;
					case GRID_SLOW_LEN:
					    grid->array[ix][iy][iz] = grid->dx / vel;
					    break;
					case GRID_SLOW2_METERS:
					    grid->array[ix][iy][iz] = 
						(1.0e-3 / vel) * (1.0e-3 / vel);
					    break;
					default:
					     nll_puterr2("ERROR: unrecognized grid type", 
							grid->chr_type);
						return(-1);
				}


				zdepth += grid->dz;
			}
			yval += grid->dy;
		}
		xval += grid->dx; 
	}


fprintf(stderr, "\nDone looping through 3D grid!\n");
	return (0);

}



/*** function to find velocity at 3G Grid location from 2D sections
	by cubic spline interpolation through sections around a circumference 
	centered at xCent, yCent ***/

inline double section2value3D(double xval, double yval, double zdepth, 
	Sect2D* psection, int numSections, double xCent, double yCent, 
	double radMin, double interpDMax, int iExtrapolateEnd)
{

	int nsec, nval, nvalue, n;
	int numSplinePoints, nIntersections;
	static char command_string[MAXLINE];
	FILE *fp_in, *fp_out;
	double radius, azimuth, rotate_tmp, vmean;
	double angleMaxDeg, aztest, aztmp;
	double az0, az1, val0, val1, diff, value3D;
	size_t iread;
	double dbl_max = DBL_MAX;

	double weightSum, weight, angDist;

	SectValue sectval[2];

	// do_spline params
	int do_spline_dim;


	/* get section values at intersection of circumference */

	radius = sqrt((xval - xCent) * (xval - xCent) 
			+ (yval - yCent) * (yval - yCent));
	nvalue = 0;
	for (nsec = 0; nsec < numSections; nsec++) {
		nIntersections = 
			getSectValueAtIntersection(sectval, xCent, yCent, radius, 
			psection + nsec, zdepth, iExtrapolateEnd);
		// load intersection values
		for (n = 0; n < nIntersections; n++) {
			secValue[nvalue++] = sectval[n];
		}
	}


	// calc azimuth of 3D grid location
	azimuth = atan2((xval - xCent), (yval - yCent)) / rpd;
	if (azimuth < 0.0)
		azimuth += 360.0;

	// check that closest section value is at most interpDMax from 3D location
	angleMaxDeg = (interpDMax / radius) / rpd;
	if (angleMaxDeg < 180.0) {
		for (nval = 0; nval < nvalue; nval++) {
			aztmp = azimuth;
			aztest = aztmp - secValue[nval].azimuth;
			while (aztest > 180.0) {
				aztmp -= 360.0;
				aztest = aztmp - secValue[nval].azimuth;
			}
			while (aztest < -180.0) {
				aztmp += 360.0;
				aztest = aztmp - secValue[nval].azimuth;
			}
			if (fabs(aztest) <= angleMaxDeg)  // ok, section is close enough
				break;
		}
		if (nval == nvalue) {
			nInterpDistMax++;
			return(-1.0);
		}
	}


	// check that we have enough values to make spline
	if (nvalue < 2) {
		splineError_NotEnoughNodes++;
	}

	// if not enough nodes, return average of values
	if (nvalue == 1) {
		vmean = 0.0;
		for (nval = 0; nval < nvalue; nval++)
			vmean += secValue[nval].value;
		nRadiusMin++;
		return(vmean / (double) nvalue);
	} else if (nvalue == 0)
		return(-1.0);


	/* sort by azimuth */
	SortSectValueAzimuth(secValue, nvalue);


	/* interpolate */

	if (splineTension < 99.0) {

		/* run spline to interpolate value at 3D point */

		/* rotate by minimum azimuth to ensure 0->360 deg spline*/
		do_spline_dim = YDIMENSION - 1;
		rotate_tmp = secValue[0].azimuth;
		/* load values */
		for (nval = 0; nval < nvalue; nval++) {
			diff = secValue[nval].azimuth - rotate_tmp;
			do_spline_t[nval] = diff;
			do_spline_y[do_spline_dim][nval] = secValue[nval].value;
		}
		/* close circuit at 360 deg */
		diff = 360.0;
		do_spline_t[nvalue] = diff;
		do_spline_y[do_spline_dim][nvalue] = secValue[0].value;

		/* call spline interpolation */
		numSplinePoints = 360;

		// function call of do_spline function only 
		 do_spline (nvalue, DO_SPLINE_LEN, 
		       &do_spline_t, 1, do_spline_y, do_spline_z, splineTension, 1,
		       0, 1.0, 6,
		       0.0, 0.0, 0.0, numSplinePoints,
		       0, 0, 0, 
		       1, 0, do_spline_x_out, do_spline_yy_out);

		/* find spline values that bracket azimuth */
		azimuth -= rotate_tmp;
		if (azimuth < 0.0)
			azimuth += 360.0;

	     	az0 = do_spline_x_out[0];
	     	val0 = do_spline_yy_out[do_spline_dim][0];
		nval = 1;
		while (nval < numSplinePoints + 1) {
 	    		az1 = do_spline_x_out[nval];
 	    		val1 = do_spline_yy_out[do_spline_dim][nval];
			if (az1 >= azimuth)
				break;
			az0 = az1;
			val0 = val1;
			nval++;
		}

		/* linear interpolation to get value3D */
		value3D = val0 + ((azimuth - az0) / (az1 - az0)) * (val1 - val0);

	} else if (splineTension < 9999.0) {

		/* do linear interpolation to get value at 3D point */
		for (nval = 0; nval < nvalue; nval++) {
			if (secValue[nval].azimuth >= azimuth) {
				val1 = secValue[nval].value;
				az1 = secValue[nval].azimuth;
				break;
			}
		}
		if (nval == 0) {
			/* azimuth is less than first section azimuth
				or greater than last section azimuth */
			val0 = secValue[nvalue - 1].value;
			az0 = secValue[nvalue - 1].azimuth - 360.0;
			val1 = secValue[0].value;
			az1 = secValue[0].azimuth;
		} else if (nval == nvalue) {
			/* azimuth is greater than last section azimuth */
			val0 = secValue[nvalue - 1].value;
			az0 = secValue[nvalue - 1].azimuth;
			val1 = secValue[0].value;
			az1 = secValue[0].azimuth + 360.0;
		} else {
			/* azimuth is between first section azimuth
				and last section azimuth */
			val0 = secValue[nval - 1].value;
			az0 = secValue[nval - 1].azimuth;
		}
	weight = fabs((azimuth - az0) / (az1 - az0));
	//weight *= weight;
	//weight = sin(weight * pi / 2.0);
	value3D = (1.0 - weight) * val0 + weight * val1;
		//value3D = val0 + ((azimuth - az0) / (az1 - az0)) * (val1 - val0);
//printf("az0 %lf azimuth %lf az1 %lf  val0 %lf value3D %lf val1 %lf weight %lf\n", az0, azimuth, az1,  val0, value3D, val1, weight);
	} else  {

		/* do weighted cosine interpolation to get value at 3D point */
		value3D = 0.0;
		weightSum = 0.0;
		for (nval = 0; nval < nvalue; nval++) {
			angDist = fabs(secValue[nval].azimuth - azimuth);
			if (angDist > 180.0)
				angDist = 360.0 - angDist;
			weight = (90.0 - angDist) / 90.0;
			if (weight > 0.0) {
				weight = sin(weight * weight * pi / 2.0);
				value3D += weight * secValue[nval].value;
				weightSum += weight;
			}
//if (nvalue < 5) printf("nval %d nvalue %d az %lf azval %lf angDist %lf weight %lf weightSum %lf val %lf value3D  %lf\n", nval, nvalue, azimuth, secValue[nval].azimuth, angDist, weight, weightSum, secValue[nval].value, value3D);
		}
		value3D /= weightSum;
	}
//if (nvalue < 5) printf("value3D %lf\n", value3D);


	/* if close to center, return weighted average of mean value 
		 and interpolated value */
	if (radius < radMin) {
		vmean = 0.0;
		for (nval = 0; nval < nvalue; nval++)
			vmean += secValue[nval].value;
		nRadiusMin++;
		return((vmean / (double) nvalue) * (radMin - radius) / radMin 
			+ value3D * radius / radMin);
	}

	return(value3D);
	
}


/*** function to find value at intersection(s) of circle with section ***/

inline int getSectValueAtIntersection(SectValue secVal[2], double xCent, double yCent, 
	double radius, Sect2D *psection, double zdepth, int iExtrapolateEnd)
{

	int nIntersections, nValues;
	double xIntrsct[2], yIntrsct[2];
	double value;

	if (psection->line_type == X_CONST_NORTH || psection->line_type == X_CONST_SOUTH) {
// !!!! ADD CODE HERE TO DEAL WITH CASE OF EXACTLY NORTH-SOUTH SECTION ???
		return(0);
	} else {
		nIntersections = findIntersectionLineAndCircle(
			psection->m_const, psection->b_const, 
			xCent, yCent, radius, xIntrsct, yIntrsct);
	}

	// no intersection
	if (nIntersections == 0) {
//printf("getSectValueAtIntersection: NO INTERSECTIONS\n");
		return(0);

	// one or two intersections
	} else {
		nValues = 0;
		secVal[nValues].azimuth = atan2(xIntrsct[0] - xCent, yIntrsct[0] - yCent) / rpd;
		if (secVal[nValues].azimuth < 0.0)
			secVal[nValues].azimuth += 360.0;
		if ((value = getGrid2DValue(psection, xIntrsct[0], yIntrsct[0], 
					zdepth, iExtrapolateEnd)) 
				> -LARGE_DOUBLE)
			secVal[nValues++].value = value;
//printf("getSectValueAtIntersection: intsct 1  xIntrsct[0] %lf yIntrsct[0] %lf -> az %lf value %lf\n", xIntrsct[0], yIntrsct[0], secVal[nValues - 1].azimuth, value);
		if (nIntersections == 2) {
			secVal[nValues].azimuth = 
				atan2(xIntrsct[1] - xCent, yIntrsct[1] - yCent) / rpd;
			if (secVal[nValues].azimuth < 0.0)
				secVal[nValues].azimuth += 360.0;
			if ((value = getGrid2DValue(
				psection, xIntrsct[1], yIntrsct[1], zdepth, iExtrapolateEnd))
					> -LARGE_DOUBLE)
				secVal[nValues++].value = value;
//printf("getSectValueAtIntersection: intsct 1  xIntrsct[1] %lf yIntrsct[1] %lf -> az %lf value %lf\n", xIntrsct[1], yIntrsct[1], secVal[nValues - 1].azimuth, value);
		}
		return(nValues);

	}

}

/*** function to find value in Sec2D at  ***/

inline double getGrid2DValue(Sect2D *psection, double xval, double yval, double zdepth,
			int iExtrapolateEnd)
{

	int ix, iy, iz;
	long ioffset;
	double dist_along_sect, xdiff, ydiff, aztest, sign;
	double value;

	// check that point is in rect containing section
	if (!iExtrapolateEnd && (xval < psection->xmin || xval > psection->xmax ||
			yval < psection->ymin || yval > psection->ymax))
		return(-VERY_LARGE_DOUBLE);

	// find 2D section node indexes of 3D point
	xdiff =  xval - psection->x0;
	ydiff =  yval - psection->y0;
	// get direction along section from sect origin
	aztest = atan2(xdiff, ydiff) / rpd;
	if (aztest < 0.0)
		aztest += 360.0;
	sign = fabs(aztest - psection->azimuth) < 90.0 ? 1.0 : -1.0;
	// calc distance along section
	dist_along_sect = sign * sqrt(xdiff * xdiff + ydiff * ydiff) + psection->x_orig_sect;
	ix = (int) (dist_along_sect / psection->grid2d.hdr->x_inc);
	iy = (int) ((zdepth - psection->grid2d.hdr->y_min) / psection->grid2d.hdr->y_inc);
//fprintf(stderr, "getGrid2DValue: iy %d, zdepth %lf, y_min %lf, y_inc %lf\n",  iy, zdepth, psection->grid2d.hdr->y_min, psection->grid2d.hdr->y_inc);

//fprintf(stderr, "getGrid2DValue: xval %lf, yval %lf, ix %d, iy %d\n",  xval, yval, ix, iy);
	/* check that node is in section */
	if (ix < 0 || ix >= psection->grid2d.hdr->nx) {
		if (!iExtrapolateEnd)
			return(-VERY_LARGE_DOUBLE);
		if (ix < 0 && -dist_along_sect < psection->extrap_min)
			ix = 0;
		else if (ix >= psection->grid2d.hdr->nx  
				&& dist_along_sect - (double) (psection->grid2d.hdr->nx - 1) 
				* psection->grid2d.hdr->x_inc < psection->extrap_max)
			ix = psection->grid2d.hdr->nx - 1;
		else
			return(-VERY_LARGE_DOUBLE);
	}
	if (iy >= psection->grid2d.hdr->ny) {
		return(getValueBelowSection(zdepth));
	}
	if (iy < 0) {
		return(-VERY_LARGE_DOUBLE);
	}

	/* get value at node */
	ioffset = iy * psection->grid2d.hdr->nx + ix;
//fprintf(stderr, "getGrid2DValue: ix %d  iy %d ioffset %ld\n", ix, iy, ioffset);
	value = *(psection->grid2d.zdata + ioffset);

//fprintf(stderr, "getGrid2DValue: ioffset %ld, value %lf\n", ioffset, value);
	return(value);

}



/*** function to find intersection points of a line with a circle ***/

/* line: y = mx +b
	circle: (x - xCent)**2 + (y - yCent)**2 = radius**2

	assumed that line is not vertical (x const, m undef)
*/

inline int findIntersectionLineAndCircle(double m, double b, 
		double xCent, double yCent, double radius, double xIntrsct[2], double yIntrsct[2])
{

	int nroots;
int n;
	double aCoeff, bCoeff, cCoeff, root[2], B;


//printf("findIntersect: m %lf b %lf xCent %lf yCent %lf radius %lf\n", m, b, xCent, yCent, radius);
	B = b - yCent;

	aCoeff = m * m + 1.0;
	bCoeff = 2.0 * (B * m - xCent);
	cCoeff = B * B - radius * radius + xCent * xCent;
	
	// calculate x val of intersection(s)
	nroots = calculateRoots(aCoeff, bCoeff, cCoeff, root);
//for (n = 0; n < nroots; n++)
//printf("findIntersect: root %d = %lf\n", n, root[n]);

	// no real roots
	if (nroots == 0) {
		return(0);

	// two real roots
	} else if (nroots == 2) {
		xIntrsct[0] = root[0];
		yIntrsct[0] = m * root[0] + b;
		xIntrsct[1] = root[1];
		yIntrsct[1] = m * root[1] + b;
		return(2);

	// one real roots
	} else {
		xIntrsct[0] = root[0];
		yIntrsct[0] = m * root[0] + b;
		return(1);
	}


}



/*** function to sort SectValue by azimuth field */

inline int SortSectValueAzimuth(SectValue* secValue, int num_val)
{

	qsort((void *) secValue, (size_t) num_val, sizeof(SectValue), 
		(int (*)(const void *, const void * )) CmpSectValueAzimuth);


	return(0);

}

/*** function to compare SectValue by azimuth field */

inline int CmpSectValueAzimuth(const SectValue *keyval, const SectValue *datum )
{
	if (keyval->azimuth < datum->azimuth)
		return(-1);
	if (keyval->azimuth > datum->azimuth)
		return(1);

	return(0);
}




/*** function to get velocity below section if isoline cut defined */

inline double getValueBelowSection(double zdepth)
{
	double vel;

	if (numIsolineCut < 1)
		return(-VERY_LARGE_DOUBLE);

	vel = isolineCut->vReplace 
		+ (zdepth - isolineCut->vReplaceRefLevel) * isolineCut->vReplaceGrad;
	vel = vel > isolineCut->velMax ?  isolineCut->velMax : vel;

	return(vel);
}


/*------------------------------------------------------------/ */
/* Anthony Lomax           | email: lomax@faille.unice.fr     / */
/* UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax / */
/* 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        / */
/* 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        / */
/*------------------------------------------------------------/ */

