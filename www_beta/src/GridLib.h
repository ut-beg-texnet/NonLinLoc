/*
 * Copyright (C) 1999-2015 Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
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


/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
161 Allee du Micocoulier, 06370 Mouans-Sartoux, France
tel: +33(0)493752502  e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


/*  GridLib.h

        include file for grid library functions

 */




#define PACKAGE  "NonLinLoc"
#define PVER  "6.04.00"
#define PDATE "26Mar2015"
/*#define PCOPYRIGHT "\nCopyright (C) 1999-2010 Anthony Lomax\n"*/
#define PCOPYRIGHT "\0"


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <time.h>
#include <dirent.h>
#include <fnmatch.h>
/* SH 07222004 added */
#include <ctype.h>
/*
#ifdef __APPLE__
#include <sys/types.h>
#endif
 */

#ifdef EXTERN_MODE
#define	EXTERN_TXT extern
#else
#define EXTERN_TXT
#endif

#include "geometry/geometry.h"
#include "alomax_matrix/alomax_matrix.h"
#include "alomax_matrix/alomax_matrix_svd.h"
#include "matrix_statistics/matrix_statistics.h"
#include "util.h"
#include "geo.h"
#include "ran1/ran1.h"
#include "map_project.h"
#include "octtree/octtree.h"

/*
// following inline probably does nothing - functions body must be defined in header file
#ifndef INLINE
#ifndef __GNUC__
#define inline
#else
#define __inline__
#endif
#endif
*/

// the following defines needed for old cc version
#ifndef SEEK_SET
#define SEEK_SET 0
#endif
#ifndef SEEK_CUR
#define SEEK_CUR 1
#endif


// 20100614 AJL -
// the following sets globally the float size (float or double) for NLL grids - modify with care!
// normal NonLinLoc usage currently require float type
// take-off angle and byte-swapping support requires float type
#ifdef GRID_FLOAT_TYPE_DOUBLE
#define GRID_FLOAT_TYPE double
#else
#define GRID_FLOAT_TYPE float
#endif


/* misc defines */

#define DOUBLE double

#define OUT_LEVEL_0 stderr
#define OUT_LEVEL_1 stdout
#define OUT_LEVEL_2 fp_null

#ifndef MAXLINE
#define MAXLINE 101
#endif
#define MAXLINE_LONG 4*MAXLINE
#define MAXSTRING 21

//#ifndef FILENAME_MAX
#undef FILENAME_MAX 	// FILENAME_MAX may be too small: was 14 on ETH HP!
#define FILENAME_MAX 1024
//#endif



#ifndef SMALL_DOUBLE
#define SMALL_DOUBLE 1.0e-20
#endif
#ifndef LARGE_DOUBLE
#define LARGE_DOUBLE 1.0e20
#endif
#ifndef VERY_SMALL_DOUBLE
#define VERY_SMALL_DOUBLE 1.0e-30
#endif
#ifndef VERY_LARGE_DOUBLE
#define VERY_LARGE_DOUBLE 1.0e30
#endif
#ifndef SMALL_FLOAT
#define SMALL_FLOAT 1.0e-20F
#endif
#ifndef LARGE_FLOAT
#define LARGE_FLOAT 1.0e20F
#endif
#ifndef VERY_LARGE_FLOAT
#define VERY_LARGE_FLOAT 1.0e30F
#endif
#ifndef VERY_SMALL_FLOAT
#define VERY_SMALL_FLOAT 1.0e-30F
#endif



/* NLL */

// DD
/* program mode */
#define MODE_ABSOLUTE  			0
#define MODE_DIFFERENTIAL  		1
EXTERN_TXT int nll_mode;


#define EXIT_NORMAL 0
#define EXIT_ERROR_MEMORY -1
#define EXIT_ERROR_USAGE -2
#define EXIT_ERROR_FILEIO -3
#define EXIT_ERROR_MODEL -4
#define EXIT_ERROR_TTIME -5
#define EXIT_ERROR_LOCATE -6
#define EXIT_ERROR_DATES -7
#define EXIT_ERROR_IO -10
#define EXIT_ERROR_MISC -100


#define GRID_UNDEF		0		/* undefined */

#define GRID_VELOCITY		1		/* velocity (km/sec) */
#define GRID_VELOCITY_METERS	11		/* velocity (m/sec) */
#define GRID_SLOWNESS		2		/* slowness (sec/km) */
#define GRID_VEL2		3		/* velocity**2 (km/sec)**2 */
#define GRID_SLOW2		4		/* slowness**2 (sec/km)**2 */
#define GRID_SLOW2_METERS	44		/* slowness**2 (sec/m)**2 */
#define GRID_SLOW_LEN		5		/* slowness*length (sec) */

#define GRID_TIME		1000		/* time (sec) 3D grid */
#define GRID_TIME_2D		1001		/* time (sec)
							2D grid / 1D model */
#define GRID_PROB_DENSITY	2001		/* probability density (may be confidence levels) */
#define GRID_MISFIT		2002		/* misfit (sec) */
#define GRID_LIKELIHOOD		2003		/* relative likelihood normalised 0-1 (may be same as probability density) */

#define GRID_ANGLE		3000		/* take-off angles 3D grid */
#define GRID_ANGLE_2D		3001		/* take-off angles 2D grid / 1D model */
#define GRID_INCLINATION	3100		/* take-off angles 3D grid */
#define GRID_INCLINATION_2D	3101		/* take-off angles 2D grid / 1D model */

#define GRID_DEPTH		4000		/* depth (km) 3D grid */
#define GRID_LENGTH		4001		/* distance or length (km) 3D grid */

#define GRID_COULOMB		5000		/* Coulomb */

#define ANGLE_MODE_NO	0
#define ANGLE_MODE_YES	1
#define ANGLE_MODE_INCLINATION	2
#define ANGLE_MODE_UNDEF	-1

/* error codes (-55000) */
#define OBS_FILE_SKIP_INPUT_LINE 		-55011
#define OBS_FILE_ARRIVALS_CROSS_YEAR_BOUNDARY 	-55022
#define OBS_FILE_INVALID_PHASE			-55033
#define OBS_FILE_INVALID_DATE			-55044
#define OBS_FILE_END_OF_EVENT			-55055
#define OBS_FILE_END_OF_INPUT			-55066
#define OBS_FILE_TWO_ARRIVALS_READ			-55077
#define OBS_FILE_FORMAT_ERROR		-55099

#define GRID_NOT_INSIDE				-55111

/* SH found comment line in UUSS format; starts with 'c' */
#define OBS_IS_COMMENT_LINE                     -55222

#define SMALLEST_EVENT_YEAR	1800
#define LARGEST_EVENT_YEAR	2100


/*------------------------------------------------------------/ */
/* macros */
/*------------------------------------------------------------/ */


/* */
/*------------------------------------------------------------/ */





/*------------------------------------------------------------/ */
/* structures */
/*------------------------------------------------------------/ */

/* grid  description */

typedef struct {
    void *buffer; /* buffer (contiguous floats) */
    void ***array; /* array access to buffer */
    int numx, numy, numz; /* grid size */
    double origx, origy, origz; /* orig (km) */
    int autox, autoy, autoz; /* automatic position of orig (1=yes) */
    double dx, dy, dz; /* len side (km) */
    int type; /* grid type (VELOCITY, SLOWNESS, SLOW2, etc) */
    char chr_type[MAXLINE]; /* char desc of grid type */
    char title[FILENAME_MAX]; /* description */
    double sum; /* sum of grid values */
    int iSwapBytes; /* flag to specify if hi/lo bytes should be swapped when reading grid from disk files */
    char float_type[MAXLINE]; /* 20100616 AJL - char desc of float type: FLOAT or DOUBLE or empty (backwards compatibility for type FLOAT */
}
GridDesc;


 // 20101209 AJL#define X_MAX_NUM_STATIONS	500
#define X_MAX_NUM_STATIONS	5000
#define X_MAX_NUM_STATIONS_DIFF	500000
#define MAX_NUM_ARRIVALS_STA	2
#define X_MAX_NUM_ARRIVALS 	MAX_NUM_ARRIVALS_STA*X_MAX_NUM_STATIONS
#define ARRIVAL_LABEL_LEN	64  // 20100409 AJL changed from 7 to 64, also changed code in GridLib.c->FindSource()

/* phase (sythetic arrival observation) */

typedef struct {
    char label[ARRIVAL_LABEL_LEN]; /* char label */
    double time; /* arrival time */
    double error; /* error in arrival time to calc */
    char error_type[MAXLINE]; /* error type */
    double error_report; /* error in arrival time  to report*/
    char error_report_type[MAXLINE]; /* error type */
}
PhaseDesc;


/* source or station */

#define MAX_NUM_SOURCES 	20000
#define SOURCE_LABEL_LEN	ARRIVAL_LABEL_LEN

typedef struct {
    int is_coord_xyz; /* xyz coord flag */
    double x, y, z; /* xyz loc (km) */
    int is_coord_latlon; /* lat/long/depth coord flag */
    double dlat, dlong, depth; /* loc (lat/long (deg), depth (km)) */
    double otime; /* origin time */
    char label[SOURCE_LABEL_LEN]; /* char label */
    int ignored; /* 1 = ignored - not used for location */
    double station_weight; /* station specific weight */
}
SourceDesc;

/* station */

typedef struct {
    char label[ARRIVAL_LABEL_LEN]; /* char label */
    double x, y, z; /* loc (km) */
    PhaseDesc phs[MAX_NUM_ARRIVALS_STA];
    double prob_active; /* probability station is active */
}
StationDesc;


/* arrival observation */

#define ARRIVAL_NULL_STR "?"
#define ARRIVAL_NULL_CHR '?'
#define FORMAT_PHASE_1 0
#define FORMAT_PHASE_2 1
#define IO_ARRIVAL_OBS 0
#define IO_ARRIVAL_ALL 1
// AJL 20091208 Zero weight phase modification
#define ARRIVAL_ERROR_NULL 9999.0
#define ARRIVAL_ERROR_NULL_TEST 9998.0
#define CODA_DUR_NULL -1.0
#define AMPLITUDE_NULL -1.0
#define PERIOD_NULL -1.0
// 20060623 AJL MAG_NULL changed to MAGNITUDE_NULL to resolve potential conflict with Earthworm defines
#define MAGNITUDE_NULL -9.9
#define INST_LABEL_LEN 32
#define COMP_LABEL_LEN 32
//INGV #define PHASE_LABEL_LEN	7
#define PHASE_LABEL_LEN	32

typedef struct {
    // 20100506 AJL - added to support preservation of observation index order for calls to NLLoc() function (e.g. from SeisComp3)
    int original_obs_index; // index of observation in order originally read from input

    char label[ARRIVAL_LABEL_LEN]; /* char label (i.e. station or site code) */
    char network[ARRIVAL_LABEL_LEN]; /* char network name */
    char time_grid_label[ARRIVAL_LABEL_LEN]; /* char label for time grid */
    char inst[INST_LABEL_LEN]; /* instrument code */
    char comp[COMP_LABEL_LEN]; /* component (ie Z N 128) */

    /* observed values */

    char phase[PHASE_LABEL_LEN]; /* char phase id */
    char onset[2]; /* char onset (ie E I) */
    char first_mot[2]; /* char first motion id */
    int quality; /* pick quality (ie weight 0 1 2 3 4) */
    int year, month, day; /* observed arrival date */
    int hour, min; /* observed arrival hour/min */
    double sec; /* observed arrival seconds */

    double error; /* error in arrival time */
    char error_type[MAXLINE]; /* error type */

    double coda_dur; /* coda duration */
    double amplitude; /* amplitude */
    double period; /* period */
    /* SH 07212004 added */
    int clipped; /* 0 = amplitude is ok
                               1 = don't use amplitude for magnitude calculation
                              -1 = amplitude is clipped */

    /* new values NLL PHASE_2 format*/
    /* 20060629 AJL - Added */
    double apriori_weight; /* a priori weight of datum */


    /* calculated values */

    double tt_error; /* travel time error (=LOCGAU if no LOCGAU2, otherwise calculated
						from travel time using LOCGAU2 params) */
    double delay; /* time delay (is subtracted from arival
						seconds when phase read */
    double elev_corr; /* elevation correction (is added to arival
						seconds when phase read */
    int day_of_year; /* day of year (of earliest arrival) */
    long double obs_time; /* corrected observed time; secs from beginning of day of year */

    /* flags */

    int flag_ignore; /* ignore arrival for location */
    int abs_time; /* absolute timing flag: 1 = has absolute time, 0 =  no */

    /* calculated values */

    double obs_centered; /* centered observed time */

    double pred_travel_time; /* predicted travel time */
    double pred_centered; /* centered predicted travel time */
    double pred_travel_time_best; /* predicted travel time from best solution */

    double cent_resid; /* residual (centered obs - centered calc) */

    double obs_travel_time; /* obs travel time (obs_time - orig time) */
    double residual; /* residual (obs - calc) */
    double weight; /* effective weight of datum */

    double dist; /* distance from hypocenter (km) */
    double azim; /* azimuth from hypocenter (deg CW from N) */

    double ray_azim; /* ray take-off azimuth (deg CW from N) */
    double ray_dip; /* ray take-off dip (0 (down) - 180 (up) deg) */
    int ray_qual; /* ray take-off quality (0 to 10) */

    double amp_mag; /* amplitude magnitude */
    double dur_mag; /* duration magnitude */

    /* statistical residuals */
    double pdf_residual_sum; /* residual (PDF weighted residual) */
    double pdf_weight_sum; /* cumulative PDF weight */

    /* travel time grid file */
    int n_companion; /* companion phase for time grids
						= nphase or -1 for none */
    int n_time_grid; /* existing phase for time grids (DEFAULT 2D grids)
						= nphase or -1 for none */
    double tfact; /* factor to multiply by time grid values
						= 1.0 or 1/VpVs */
    char fileroot[FILENAME_MAX]; /* root name for grid file */
    FILE* fpgrid; /* 3D travel time grid file */
    FILE* fphdr; /* 3D grid file header */
    GridDesc gdesc; /* description for grid in disk file */
    GridDesc sheetdesc; /* description for dual-sheet in memory */

    SourceDesc station; /* station description */

    double station_weight; /* station specific weight */

    // 20101005 AJL - slowness value, for OT_STACK, etc.
    double slowness; // wave slowness at current hypocenter
    int isP;    // is a P phase
    int isS;    // is an S phase

    // DD
    // see hypoDD doc (Open File Report 01-113)
    int xcorr_flag; /* 1 = xcorr data, 0 = catalog data */
    long int dd_event_id_1, dd_event_id_2; /* ID of event 1, event 2 */
    int dd_event_index_1, dd_event_index_2; /* HypoDesc index of event 1, event 2 */
    double dd_dtime; /* Differential time (s) between event 1 and event 2 at station. DT = T1-T2. */

}
ArrivalDesc;


/* focal mechanism (taken from FPFIT file fpfit.doc) */

typedef struct {
    double dlat, dlong, depth; // lat, lon, depth of hypocenter

    double dipDir; // Dip direction (downdip azimuth in
    //	degrees,clockwise from north)
    double dipAng; // Dip angle in degrees down from horizontal
    double rake; // Rake in degrees: 0=left lateral, 90=reverse,
    //	+-180=right lateral, -90=normal
    double misfit; // Solution misfit value. 0=perfect fit,
    //	1=perfect misfit (never exceeds 0.5 in reality).
    int nObs; // Number of first motion observations used in solution

    /* FPFIT parameters unused by NonLinLoc --------------------------- */
    double misfit90; // Solution misfit value + 90% confidence estimate.
    //	Useful for testing whether restricted solutions
    //	are within uncertainty region of free solution
    double staDist; // Station distribution ratio (0-1).  Lower numbers indicate
    //	data lie near nodal planes.
    double ratioMH; // (# of machine picks)/(# of hand picks)
    double conf90strike; // maximum half-width of 90% confidence range of strike
    double conf90dip; // maximum half-width of 90% confidence range of dip
    double conf90rake; // maximum half-width of 90% confidence range of rake
    char convFlag[2]; // Convergence flag: C=no convergence; otherwise blank
    char multSolFlag[2]; // Multiple solution flag: *=multiple; otherwise blank

} FocalMech;

/* hypocenter
 *
 * QML indicates fields added for compatibility with QuakeML OriginQuality attributes (AJL 201005)

 */

typedef struct {
    double x, y, z; /* loc (km) */
    int ix, iy, iz; /* loc (grid) */
    double dlat, dlong, depth; /* loc (geographic) */
    int year, month, day; /* date */
    int hour, min; /* hour/min */
    double sec; /* origin time (s) */
    long double time; /* secs from beginning of day */

    int nreadings; /* number arrivals used in solution */
    double gap; /* largest azimuth separation between stations as seen from epicenter (deg) */
    double gap_secondary; // QML - largest azimuth separation between stations filled by a single station as seen from epicenter (deg)
    double dist; /* epi distance to closest sta (km or deg (MODE_GLOBAL)) */
    double rms; /* rms error of arrival times */

    // QML fields added for compatibility with QuakeML OriginQuality attributes (AJL 201005)
    int associatedPhaseCount; // QML - Number of associated phases, regardless of their use for origin computation.
    // [->nreadings] int usedPhaseCount; // QML - Number of defining phases, i. e., phase observations that were actually used for computing
    // the origin. Note that there may be more than one defining phase per station.
    int associatedStationCount; // QML - Number of stations at which the event was observed.
    int usedStationCount; // QML - Number of stations from which data was used for origin computation.
    int depthPhaseCount; // QML - Number of depth phases (typically pP, sometimes sP) used in depth computation.
    // [->rms] double standardError; // QML - RMS of the travel time residuals of the arrivals used for the origin computation. Unit: s
    // [->gap] double azimuthalGap; // QML - Largest azimuthal gap in station distribution as seen from epicenter. Unit: deg
    // [->gap_secondary] double secondaryAzimuthalGap; // QML - Secondary azimuthal gap in station distribution, i. e., the largest azimuthal gap a station closes. Unit: deg
    char groundTruthLevel[8]; // QML - String describing ground-truth level, e. g. GT0, GT5, etc.
    double minimumDistance; // QML - Epicentral distance of station closest to the epicenter. Unit: km
    double maximumDistance; // QML - Epicentral distance of station farthest from the epicenter. Unit: km
    double medianDistance; // QML - Median epicentral distance of used stations. Unit: km

    // SED-ETH fields added for compatibility with legacy SED location quality indicators (AJL 201006)
    // SED_Origin errx 33.0775  erry 7.5298  errz 27.5298  diffMaxLikeExpect 4.9584 quality B
    double diffMaxLikeExpect; // Euclidean distance between XYZ maximum likelihood and expectation hypocenters
    char qualitySED; // A-D SED quality indicator

    double amp_mag; /* amplitude magnitude */
    int num_amp_mag; /* number of readings */
    double dur_mag; /* duration magnitude */
    int num_dur_mag; /* number of readings */
    /* SH 07232004 added */
    double mag_err; /* error in magnitude */

    long double probmax; /* probability density at hypo (maximum in grid) */ // 20130314 C Satriano, AJL - changed to long double

    double misfit; /* misfitat hypo (minimum in grid) */
    double grid_misfit_max; /* max misfit in loc grid */

    Vect3D expect; /* "traditional" expectation */
    Mtrx3D cov; /* "traditional" covariance matrix - units always in km, including for MODE_GLOBAL, see CalcCovarianceSamplesGlobal() */
    double expect_dlat, expect_dlong; /* lat / long trans of expectation */

    Ellipse2D ellipse; /* horizontal (X/Y) error ellipse description */
    Ellipsoid3D ellipsoid; /* error ellipsoid description */
    FocalMech focMech; /* focal mechanism description */

    double VpVs; /* Vp / Vs ratio */
    int nVpVs; /* num P-S pairs used to calc VpVs */
    double tsp_min_max_diff; /* max P-S - min P-S (sec) */

    char label[MAXLINE]; /* char label */
    char fileroot[2 * MAXLINE]; /* fileroot for original location */
    char comment[2 * MAXLINE]; /* char comment */
    char signature[6 * MAXLINE]; /* char signature/program/date/etc */
    char searchInfo[2 * MAXLINE]; /* char search type dependent info */

    char locStat[32]; /* char location status LOCATED, ABORTED, IGNORED, REJECTED */
    char locStatComm[2 * MAXLINE]; /* char location status comment */

    // DD
    // see hypoDD doc (Open File Report 01-113)
    long int event_id;
    // DD Metropolis
    int nSamples;
    int ipos;
    int nScatterSaved;
    int numClipped;
    int flag_ignore; /* ignore hypocenter for location */
    double dotime; /* sec correction to otime */
    //
}
HypoDesc;


/* take-off angles */

#define ANGLE_QUALITY_CUTOFF 5		/* value to determine "bad" angles */
#define ANGLES_OFFSET 16

typedef union {
    unsigned short ival[2]; /* unsigned short values:
		(ival[0] bits 0-3 = quality: 0 to 10;
			bits 4-15 = dip: 0 (down) to 1800 (up) in tenths deg)
		(ival[1] = azimuth: 0 to 3600 in tenths deg) */

    float fval; /* float value (dummy) */
}
TakeOffAngles;



/* */
/*------------------------------------------------------------/ */



/*------------------------------------------------------------/ */
/* globals  */
/*------------------------------------------------------------/ */

EXTERN_TXT char fn_control[MAXLINE]; /* control file name */
EXTERN_TXT FILE *fp_control; /* control file pointer */
EXTERN_TXT char fn_output[MAXLINE]; /* output file name */

/* miscellaneous */
EXTERN_TXT int RandomNumSeed;
EXTERN_TXT int NumFilesOpen;
EXTERN_TXT int NumGridBufFilesOpen, NumGridHdrFilesOpen;
EXTERN_TXT int NumAllocations;

/* algorithm constants */
EXTERN_TXT int prog_mode_3d;
EXTERN_TXT int prog_mode_2dto3d;

// mode
#define MODE_RECT 		0	// rectangular cartesian x(km),y(km),z:depth(km)
#define MODE_GLOBAL 		1	// spherical x:longitdue(deg),y:latittude(deg),z:depth(km)
EXTERN_TXT int GeometryMode;

/* 3D grid description */
EXTERN_TXT int grid_type; /* grid type (VELOCITY, SLOWNESS, SLOW2, etc) */
EXTERN_TXT GridDesc grid_in;

/* source */
EXTERN_TXT int NumSources;
EXTERN_TXT SourceDesc Source[MAX_NUM_SOURCES];

/* stations */
//EXTERN_TXT int NumStations;
EXTERN_TXT StationDesc Station[MAX_NUM_SOURCES];

/* arrivals */
EXTERN_TXT int PhaseFormat;
EXTERN_TXT int MAX_NUM_STATIONS;
EXTERN_TXT int MAX_NUM_ARRIVALS;
EXTERN_TXT int NumArrivals;
EXTERN_TXT ArrivalDesc* Arrival;

/* hypocenter */
EXTERN_TXT HypoDesc Hypocenter;

/* geographic transformations (lat/long <=> x/y) */
#define NUM_PROJ_MAX 		3
#define MAP_TRANS_UNDEF 	-1
#define MAP_TRANS_NONE	 	0
#define MAP_TRANS_GLOBAL 	1
#define MAP_TRANS_SIMPLE 	2
#define MAP_TRANS_LAMBERT 	3
#define MAP_TRANS_SDC 		4
EXTERN_TXT char map_trans_type[NUM_PROJ_MAX][MAXLINE]; /* name of projection */
EXTERN_TXT int map_itype[NUM_PROJ_MAX]; /* int id of projection */
EXTERN_TXT char MapProjStr[NUM_PROJ_MAX][2 * MAXLINE]; /* string description of proj params */
EXTERN_TXT char map_ref_ellipsoid[NUM_PROJ_MAX][MAXLINE]; /* name of reference ellipsoid */
/* general map parameters */
EXTERN_TXT double map_orig_lat[NUM_PROJ_MAX], map_orig_long[NUM_PROJ_MAX], map_rot[NUM_PROJ_MAX];
EXTERN_TXT double map_cosang[NUM_PROJ_MAX], map_sinang[NUM_PROJ_MAX]; /* rotation */
/* LAMBERT projection parameters */
EXTERN_TXT double map_lambert_1st_std_paral[NUM_PROJ_MAX], map_lambert_2nd_std_paral[NUM_PROJ_MAX];
/* SDC Short Distance Coversion projection parameters */
EXTERN_TXT double map_sdc_xltkm[NUM_PROJ_MAX], map_sdc_xlnkm[NUM_PROJ_MAX];
#define MAP_TRANS_SDC_DRLT 0.99330647

/* constants */
EXTERN_TXT double cPI;
EXTERN_TXT double cRPD;
EXTERN_TXT double c111;

/* include file */
EXTERN_TXT char fn_include[FILENAME_MAX];
EXTERN_TXT FILE* fp_include;
EXTERN_TXT FILE* fp_input_save;

/* take-off angle */
EXTERN_TXT TakeOffAngles AnglesNULL;

/* quality to error mapping (hypo71, etc) */
#define MAX_NUM_QUALITY_LEVELS 50
EXTERN_TXT double Quality2Error[MAX_NUM_QUALITY_LEVELS];
EXTERN_TXT int NumQuality2ErrorLevels;

/* model coordinates */
#define COORDS_RECT	0
#define COORDS_LATLON	1
int ModelCoordsMode;

/* */
/*------------------------------------------------------------/ */




/*------------------------------------------------------------/ */
/* function declarations */
/*------------------------------------------------------------/ */

/*void disp_usage(char* , char*);
void puterr(char *pm);
void puterr2(char *, char *);
void putmsg(int , char *);*/
void SetConstants(void);

int get_control(char*);
int get_grid(char*);
int convert_grid_type(GridDesc*, int);
int display_grid_param(GridDesc*);
int get_mcsyn(char*);
int get_path_method(char*);
int GetNextSource(char*);
int GetSource(char*, SourceDesc*, int);
SourceDesc* FindSource(char* label);
char* projection_str2transform_str(char* trans_str, char* proj_str);
int get_transform(int, char*);

void get_velfile(char*);
int convertCoordsRect(int, int, double, double, double *, double *);
int latlon2rect(int, double, double, double*, double*);
int rect2latlon(int, double, double, double*, double*);
double rect2latlonAngle(int, double);
double latlon2rectAngle(int, double);
double getGMTJVAL(int, char*, double, double, double, double, double, double);
int ConvertSourceLoc(int, SourceDesc *, int, int, int);
int ConvertASourceLocation(int n_proj, SourceDesc *srce_in, int toXY, int toLatLon);


/* grid functions */
void* AllocateGrid(GridDesc*);
void FreeGrid(GridDesc*);
int InitializeGrid(GridDesc*, GRID_FLOAT_TYPE);
void*** CreateGridArray(GridDesc*);
void DestroyGridArray(GridDesc*);
void DuplicateGrid(GridDesc*, GridDesc*, char *);
int CheckGridArray(GridDesc*, double, double, double, double);
int SumGrids(GridDesc* pgrid_sum, GridDesc* pgrid_new, FILE* fp_grid_new, double factor);
int WriteGrid3dBuf(GridDesc*, SourceDesc*, char*, char*);
int WriteGrid3dHdr(GridDesc*, SourceDesc*, char*, char*);
int ReadGrid3dBuf(GridDesc*, FILE*);
int ReadGrid3dHdr(GridDesc*, SourceDesc*, char*, char*);
int ReadGrid3dHdr_grid_description(FILE *fpio, GridDesc* pgrid);
int ReadGrid3dBufSheet(GRID_FLOAT_TYPE *, GridDesc*, FILE*, int);
GRID_FLOAT_TYPE ReadAbsGrid3dValue(FILE*, GridDesc*, double, double,
        double, int);
int SwapBytes(float *buffer, long int bufsize);
int OpenGrid3dFile(char *, FILE **, FILE **, GridDesc*,
        char*, SourceDesc*, int);
void CloseGrid3dFile(FILE **, FILE **);
GRID_FLOAT_TYPE* ReadGridFile(GRID_FLOAT_TYPE* values, char *fname, char* file_type, double* xloc, double* yloc, double* zloc, int nvalues, int iSwapBytes);
GRID_FLOAT_TYPE ReadGrid3dValue(FILE *, int, int, int, GridDesc*);
DOUBLE InterpCubeLagrange(DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE,
        DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE);
float InterpCubeAngles(DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE,
        DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE);
GRID_FLOAT_TYPE ReadAbsInterpGrid3d(FILE *, GridDesc*, double, double,
        double);
DOUBLE InterpSquareLagrange(DOUBLE, DOUBLE,
        DOUBLE, DOUBLE, DOUBLE, DOUBLE);
DOUBLE ReadAbsInterpGrid2d(FILE *, GridDesc*,
        double, double);
int isOnGridBoundary(double, double, double, GridDesc*, double, double, int);
int IsPointInsideGrid(GridDesc*, double, double, double);
int IsGridInside(GridDesc*, GridDesc*, int);
int IsGrid2DBigEnough(GridDesc*, GridDesc*, SourceDesc*,
        double, double, double, double);

/* statistics functions */
double GaussDev();
void TestGaussDev();
int GenTraditionStats(GridDesc*, Vect3D*, Mtrx3D*, FILE*);
Vect3D CalcExpectation(GridDesc*, FILE*);
Mtrx3D CalcCovariance(GridDesc*, Vect3D*, FILE*);
Mtrx3D CalcCovarianceSamples(float*, int, Vect3D*);

/* hypocenter functions */
int WriteLocation(FILE*, HypoDesc*, ArrivalDesc*, int, char*,
        int, int, int, GridDesc*, int);
int WritePhases(FILE *fpio, HypoDesc* phypo, ArrivalDesc* parrivals,
        int narrivals, char* filename,
        int iWriteArrivals, int iWriteEndLoc, int iWriteMinimal,
        GridDesc* pgrid, int n_proj, int io_arrival_mode);
int GetHypLoc(FILE*, const char*, HypoDesc*, ArrivalDesc*, int*, int, GridDesc*, int);
int ReadArrival(char*, ArrivalDesc*, int);
int WriteArrival(FILE*, ArrivalDesc*, int);
int WriteArrivalHypo(FILE*, ArrivalDesc*, int);
int ReadHypStatistics(FILE **, char*, Vect3D*, Vect3D*,
        Mtrx3D*, Ellipsoid3D*, ArrivalDesc*, int*);
int ReadFocalMech(FILE **pfpio, char* fnroot_in, FocalMech* pfocalMech,
        ArrivalDesc* parrivals, int *pnarrivals);
void Qual2Err(ArrivalDesc *);
int Err2Qual(ArrivalDesc *);
int GetQuality2Err(char*);

/* strucutre utility functions */
HypoDesc* cloneHypoDesc(HypoDesc *phypo_orig);
ArrivalDesc* cloneArrivalDescArray(ArrivalDesc* parrivals_orig, int narrivals);
GridDesc* cloneGridDesc(GridDesc* pgrid_orig);

/* source/station functions */
double GetEpiDist(SourceDesc*, double, double);
double GetEpiAzim(SourceDesc*, double, double);
double GetEpiDistSta(StationDesc*, double, double);
double GetEpiAzimSta(StationDesc*, double, double);
double Dist3D(double, double, double, double, double, double);
double calcAveInterStationDistance(SourceDesc *stations, int numStations);
int stationLocationIsKnown(double x, double y);

/* date functions */
int DayOfYear(int, int, int);
void MonthDay(int, int, int*, int*);
int Month2Int(char* cmonth);

/* time functions */
char* CurrTimeStr(void);

/* file list functions */
int ExpandWildCards(char*, char[][FILENAME_MAX], int);
int fnmatch_wrapper(const struct dirent* entry);
EXTERN_TXT char ExpandWildCards_pattern[FILENAME_MAX];

/* string / char functions */
int TrimString(char*);

/* arrival sorting functions */
int SortArrivalsTime(ArrivalDesc*, int);
int CmpArrivalsTime(const ArrivalDesc *, const ArrivalDesc *);
int SortArrivalsIgnore(ArrivalDesc*, int);
int CmpArrivalsIgnore(const ArrivalDesc *, const ArrivalDesc *);
int SortArrivalsDist(ArrivalDesc*, int);
int CmpArrivalsDist(const ArrivalDesc *, const ArrivalDesc *);

/* generic sorting functions */
int SortDoubles(double *, int);
int CmpDoubles(const double *, const double *);

/* include file functions */
int GetIncludeFile(char*, FILE **);
void SwapBackIncludeFP(FILE **);

/* source take-off angle functions */
TakeOffAngles SetTakeOffAngles(double, double, int);
void SetAnglesFloat(TakeOffAngles*, float);
int GetTakeOffAngles(TakeOffAngles *, double *, double *, int *);
int ReadTakeOffAnglesFile(char *, double, double, double,
        double *, double *, int *, double, int);
int CalcAnglesGradient(GridDesc* ptgrid, GridDesc* pagrid, int angle_mode, int grid_mode);
TakeOffAngles GetGradientAngles(double vcent, double xlow, double xhigh,
        double ylow, double yhigh, double zlow, double zhigh,
        double dx, double dy, double dz, int iflag2D,
        double *pazim, double *pdip, int *piqual);
int CalcAnglesQuality(double grad_low, double grad_high);

int SetModelCoordsMode(int);

/* FORTRAN like read statements */
int LineIsBlank(char *);
int ReadFortranString(char*, int, int, char*);
int ReadFortranInt(char*, int, int, int*);
int ReadFortranReal(char*, int, int, double*);


/* OctTree support */
int ConvertOctTree2Grid(Tree3D* tree, double dx, double dy, double dz, char *grid_type, GridDesc *pgrid_out);

/* misc functions */
int ReadFpfitSum(FILE *fp_in, HypoDesc *phypo);

// DD
int ReadHypoDDInitHypo(FILE *fp_in, HypoDesc *phypo, int n_proj);
int WriteHypoDDInitHypo(FILE *fp_out, HypoDesc *phypo);
int WriteHypoDDXCorrDiff(FILE *fp_out, int num_arrivals, ArrivalDesc *arrival, HypoDesc *phypo);
int WriteDiffArrival(FILE* fpio, HypoDesc* hypos, ArrivalDesc* parr, int iWriteType);


/* */
/*------------------------------------------------------------/ */


