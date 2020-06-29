/*
 * Copyright (C) 1999-2012 Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
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


/*   fmm2grid.c
 * Program to convert arrtimes.dat ASCII grid to Grid3D format grid
 * arrtimes.dat ASCII grid are output by the ANU-FMM multi-stage 3D
 * fast marching code  (http://rses.anu.edu.au/seismology/soft/fmmcode/)
 * when the save_timefields_mode is set to True
 */


/*
        history:

        ver 01    01May2012  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

 */



#define PNAME  "fmm2grid"

#include "GridLib.h"
#include "velmod.h"


/* globals  */

char fn_fg_input[MAXLINE];
char fn_fg_output[MAXLINE];


/* wave type (P, S, PmP...) for each time grids */
#define MAX_NUM_PATH_PHASE_CODES 999
char PathPhaseCodes[MAX_NUM_PATH_PHASE_CODES][64];
int NumPathPhaseCodes;


/* function declarations */

int ReadFmm2GridInput(FILE*);
int ReadFmmTimesHeader(FILE* fp_fmmtimes, GridDesc* fmm_times);
int ReadFmmTimes(FILE* fp_fmmtimes, GridDesc* fmm_times);
int FmmTimesToGrid3d(GridDesc*, GridDesc*);
int AllocateFmmTimes(GridDesc*);
void FreeFmmTimes(GridDesc*);
int FmmTimesToFMM(GridDesc* fmm_times, GridDesc* grid, FILE* fp_fmmVgridFile, FILE* fp_fmmPropGridFile);

/*** program to generate  3-D vel/slowness grid */


#define NARGS 2
#define LENDELTA 6

#define LENVEL 5

int main(int argc, char *argv[]) {

    int istat;
    int nPathGrid;
    char fileRoot[MAXLINE];

    GridDesc time_grid; /* model grid */
    GridDesc fmm_times; /* Velocity model */



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
    NumPathPhaseCodes = 0;
    SetConstants();


    /* read control file */

    strcpy(fn_control, argv[1]);
    if ((fp_control = fopen(fn_control, "r")) == NULL) {
        nll_puterr("ERROR: opening control file.");
        exit(EXIT_ERROR_FILEIO);
    }

    if ((istat = ReadFmm2GridInput(fp_control)) < 0) {
        exit(EXIT_ERROR_FILEIO);
    }
    time_grid = grid_in;

    /* convert source location coordinates  */
    istat = ConvertSourceLoc(0, Source, NumSources, 1, 1);

    /* determine model coordinates mode - rect or latlon */
    SetModelCoordsMode(num_surfaces);


    /* read FMM Times input file */

    FILE* fp_fmmtimes = NULL;
    if ((fp_fmmtimes = fopen(fn_fg_input, "r")) == NULL) {
        nll_puterr("ERROR: opening FMM Times input file.");
        return (-1);
    }

    if ((istat = ReadFmmTimesHeader(fp_fmmtimes, &fmm_times)) < 0) {
        nll_puterr("ERROR: Reading FMM Times input file header lines.");
        exit(EXIT_ERROR_FILEIO);
    }

    /* initialize 3D grids */

    /* allocate FMM time grid */
    fmm_times.buffer = AllocateGrid(&fmm_times);
    if (fmm_times.buffer == NULL) {
        nll_puterr("ERROR: allocating memory for FMM time grid buffer.");
        exit(EXIT_ERROR_MEMORY);
    }
    /* create array access pointers */
    fmm_times.array = CreateGridArray(&fmm_times);
    if (fmm_times.array == NULL) {
        nll_puterr("ERROR: creating array for accessing FMM time grid buffer.");
        exit(EXIT_ERROR_MEMORY);
    }
    fmm_times.type = GRID_TIME;

    /* allocate NLL time grid */
    time_grid.buffer = AllocateGrid(&time_grid);
    if (time_grid.buffer == NULL) {
        nll_puterr("ERROR: allocating memory for NLL time grid buffer.");
        exit(EXIT_ERROR_MEMORY);
    }
    /* create array access pointers */
    time_grid.array = CreateGridArray(&time_grid);
    if (time_grid.array == NULL) {
        nll_puterr("ERROR: creating array for accessing NLL time grid buffer.");
        exit(EXIT_ERROR_MEMORY);
    }


    /* create time grid for each set of FMM times */

    int n_values_expected = fmm_times.numx * fmm_times.numy * fmm_times.numz;
    for (nPathGrid = 0; nPathGrid < NumPathPhaseCodes; nPathGrid++) {

        sprintf(fileRoot, "%s.%s.%s", fn_fg_output, PathPhaseCodes[nPathGrid], (Source + 0)->label);
        sprintf(MsgStr, "Creating time grid files: %s.time.*", fileRoot);
        nll_putmsg(1, MsgStr);

        int nread = ReadFmmTimes(fp_fmmtimes, &fmm_times);
        if (nread != n_values_expected) {
            sprintf(MsgStr, "ERROR: number values read for time grid %d != number expected (nx*ny*nz %d)", nread, n_values_expected);
            nll_puterr(MsgStr);
            exit(EXIT_ERROR_FILEIO);
        }

        // load vel model to grid

        if ((istat = FmmTimesToGrid3d(&fmm_times, &time_grid)) < 0) {
            nll_puterr("ERROR: loading velocity model to grid.");
            exit(EXIT_ERROR_MODEL);
        }

        // save grid to disk

        if ((istat = WriteGrid3dBuf(&time_grid, (Source + 0), fileRoot, "time")) < 0) {
            nll_puterr("ERROR: writing time grid to disk.");
            exit(EXIT_ERROR_IO);
        }

    }

    exit(EXIT_NORMAL);

}

/*** function to read FMM Times input file */

int ReadFmmTimesHeader(FILE* fp_fmmtimes, GridDesc* fmm_times) {

    char line[2 * MAXLINE];

    /*
     arrtimes.dat
    ------------
    The output file arrtimes.dat is generated if you request that some grids
    of arrival times be saved by setting save_timefields_mode in the file
    mode_set.in to true (T). The specific paths for which arrival time grids
    are saved are specified in the input file gridsave.in.

    The format of arrtimes.dat is as follows:

    First line : nr of nodes in radial, latitude and longitude direction

    Second line : grid spacing in r, lat and long (km and degrees respectively)

    Third line : origin of the grid in radius(km) lat and long (degrees)

    Fourth line : the number of sets of arrival time values grids that follow

    For each set of arrival times a first line is written that specifies the
    source  and the path from this source from which these arrival times on
    the grid are calculated. After this the arrival times on the grid nodes
    follow, one per line, with the r-index varying fastest, then the
    latitude-index, and the longitude-index varying slowest.  Missing values
    are possible in general, and are indicated by a value of -1.0.
     */


    // First line : nr of nodes in radial, latitude and longitude direction
    if ((fgets(line, 2 * MAXLINE, fp_fmmtimes)) &&
            (sscanf(line, "%d %d %d", &fmm_times->numz, &fmm_times->numy, &fmm_times->numx) == 3)) {
        fprintf(stderr, "ReadFmmTimesHeader: no. grid nodes in rad/lat/lon -> %d %d %d \n", fmm_times->numz, fmm_times->numy, fmm_times->numx);
    } else {
        nll_puterr2("ERROR: reading first line of FMM Times input file: %s", line);
        return (-1);
    }

    // Second line : grid spacing in r, lat and long (km and degrees respectively)
    if ((fgets(line, 2 * MAXLINE, fp_fmmtimes)) &&
            (sscanf(line, "%lf %lf %lf", &fmm_times->dz, &fmm_times->dy, &fmm_times->dx) == 3)) {
        fprintf(stderr, "ReadFmmTimesHeader: grid spacing in rad/lat/lon -> %f %f %f \n", fmm_times->dz, fmm_times->dy, fmm_times->dx);
    } else {
        nll_puterr2("ERROR: reading second line of FMM Times input file: %s", line);
        return (-1);
    }

    // Third line : origin of the grid in radius(km) lat and long (degrees)
    if ((fgets(line, 2 * MAXLINE, fp_fmmtimes)) &&
            (sscanf(line, "%lf %lf %lf", &fmm_times->origz, &fmm_times->origy, &fmm_times->origx) == 3)) {
        // convert z orig from FMM (min radius) to NLL (top of grid) convention
        fmm_times->origz = AVG_ERAD - fmm_times->origz - fmm_times->dz * (double) (fmm_times->numz - 1);
        fprintf(stderr, "ReadFmmTimesHeader: origin in depth/lat/lon -> %f %f %f \n", fmm_times->origz, fmm_times->origy, fmm_times->origx);
    } else {
        nll_puterr2("ERROR: reading third line of FMM Times input file: %s", line);
        return (-1);
    }

    // Fourth line : the number of sets of arrival time values grids that follow
    int n_grids = -1;
    if ((fgets(line, 2 * MAXLINE, fp_fmmtimes)) &&
            (sscanf(line, "%d", &n_grids) == 1)) {
        fprintf(stderr, "ReadFmmTimesHeader: number of sets of arrival time values grids -> %d\n", n_grids);
        if (n_grids != NumPathPhaseCodes) {
            sprintf(MsgStr, "ERROR: number of sets of arrival time values grids is not equal to number of path phase codes (FGPATHPHASECODE): %d", NumPathPhaseCodes);
            nll_puterr(MsgStr);
            return (-1);
        }
    } else {
        nll_puterr2("ERROR: reading fourth line of FMM Times input file: %s", line);
        return (-1);
    }

    return (0);

}

int ReadFmmTimes(FILE* fp_fmmtimes, GridDesc* fmm_times) {

    /*
     For each set of arrival times a first line is written that specifies the
    source  and the path from this source from which these arrival times on
    the grid are calculated. After this the arrival times on the grid nodes
    follow, one per line, with the r-index varying fastest, then the
    latitude-index, and the longitude-index varying slowest.  Missing values
    are possible in general, and are indicated by a value of -1.0.
     */

    int i1, i2, i3;
    char line[2 * MAXLINE];
    fgets(line, 2 * MAXLINE, fp_fmmtimes); // source line (ignored)
    if (sscanf(line, "%d %d %d", &i1, &i2, &i3) != 3) {
        sprintf(MsgStr, "ERROR: reading arrival times first line (source/path): %s\n", line);
        nll_puterr(MsgStr);
        return (0);
    }

    double value;
    int nread = 0;
    int i, j, k;
    for (i = 0; i < fmm_times->numx; i++) {
        for (j = 0; j < fmm_times->numy; j++) {
            for (k = fmm_times->numz - 1; k >= 0; k--) { // FMM starts from bottom of grid
                if (fscanf(fp_fmmtimes, "%lf ", &value) != 1) {
                    sprintf(MsgStr, "ERROR: reading time value in FMM Times input file: nread %d: x %d  y %d  z %d ", nread, i, j, k);
                    nll_puterr(MsgStr);
                    return (nread);
                }
                ((GRID_FLOAT_TYPE ***) fmm_times->array)[i][j][k] = value;
                nread++;
            }
        }
    }

    return (nread);
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

int FmmTimesToGrid3d(GridDesc* fmm_times, GridDesc* grid) {

    int ix, iy, iz;
    double xval, yval, xloc, yloc, zdepth;
    double time;

    int n_time_error = 0;


    /* generate grid values */

    xloc = yloc = 0.0;
    xval = grid->origx;
    for (ix = 0; ix < grid->numx; ix++) {
        yval = grid->origy;
        for (iy = 0; iy < grid->numy; iy++) {
            rect2latlon(0, xval, yval, &yloc, &xloc);
            zdepth = grid->origz;
            for (iz = 0; iz < grid->numz; iz++) {
                time = (GRID_FLOAT_TYPE) ReadAbsInterpGrid3d(NULL, fmm_times, xloc, yloc, zdepth, 0);
                // test nearest neighbour read to avoid no time data near interface (e.g. no times at base of crust)
                if (1 && time < 0.0)    // try nearest neighbour
                    time = (GRID_FLOAT_TYPE) ReadAbsGrid3dValue(NULL, fmm_times, xloc, yloc, zdepth, 0);
                //
                if (time < 0.0) {
                    n_time_error++;
                    time = -1.0;
                    //return (-1);
                }

                switch (grid->type) {

                    case GRID_TIME:
                        ((GRID_FLOAT_TYPE ***) grid->array)[ix][iy][iz] = time;
                        break;

                    default:
                        nll_puterr("ERROR: unrecognized grid type.");
                        return (-1);

                }
                zdepth += grid->dz;
            }
            yval += grid->dy;
        }
        xval += grid->dx;
    }

    if (n_time_error > 0) {
        sprintf(MsgStr, "WARNING: could not get FMM time for %d grid points.", n_time_error);
        nll_puterr(MsgStr);
    }

    return (0);

}

/*** function to read output file name ***/

int get_fg_outfile(char* line1) {

    sscanf(line1, "%s", fn_fg_output);

    sprintf(MsgStr, "fmm2grid files:  Output: %s.*",
            fn_fg_output);
    nll_putmsg(3, MsgStr);

    return (0);
}

/*** function to read wave type ***/

int get_fg_type(char* line1) {

    if (NumPathPhaseCodes >= MAX_NUM_PATH_PHASE_CODES) {
        nll_puterr("WARNING: maximum number of path phase codes reached, ignoring path phase code.");
        return (-1);
    }


    sscanf(line1, " %s", PathPhaseCodes[NumPathPhaseCodes]);

    sprintf(MsgStr, "fmm2grid path phase code:  %s", PathPhaseCodes[NumPathPhaseCodes]);
    nll_putmsg(3, MsgStr);

    NumPathPhaseCodes++;


    return (0);
}

/*** function to read input file name ***/

int get_fg_inpfile(char* line1) {

    char type[MAXLINE];

    sscanf(line1, "%s %s", fn_fg_input, type);

    sprintf(MsgStr, "fmm2grid files:  Input: %s   Type: %s", fn_fg_input, type);
    nll_putmsg(3, MsgStr);

    return (0);
}

/*** function to read input file */

int ReadFmm2GridInput(FILE* fp_input) {
    int istat, iscan;
    char param[MAXLINE], *pchr;
    char line[2 * MAXLINE], *fgets_return;

    int flag_control = 0, flag_inpfile = 0, flag_outfile = 0, flag_grid = 0,
            flag_source = 0, flag_path_phase_code = 0,
            flag_trans = 0;
    int flag_include = 1;

    /* read each input line */


    while ((fgets_return = fgets(line, 2 * MAXLINE, fp_input)) != NULL
            || fp_include != NULL) {


        /* check for end of include file */

        if (fgets_return == NULL && fp_include != NULL) {
            SwapBackIncludeFP(&fp_input);
            continue;
        }


        istat = -1;

        /*read parmeter line */

        if ((iscan = sscanf(line, "%s", param)) < 0)
            continue;

        /* skip comment line or white space */

        if (strncmp(param, "#", 1) == 0 || isspace(param[0]))
            istat = 0;



        /* read include file params and set input to include file */

        if (strcmp(param, "INCLUDE") == 0) {
            if ((istat = GetIncludeFile(strchr(line, ' '),
                    &fp_input)) < 0) {
                nll_puterr("ERROR: processing include file.");
                flag_include = 0;
            }
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


        /* read Input file name (INPfile) */

        if (strcmp(param, "FGINP") == 0) {
            if ((istat = get_fg_inpfile(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading 3D Velocity input parameters.");
            else
                flag_inpfile = 1;
        }


        /* read output file name  */

        if (strcmp(param, "FGOUT") == 0) {
            if ((istat = get_fg_outfile(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading fmm2grid output file name.");
            else
                flag_outfile = 1;
        }

        /* read source params */

        if (strcmp(param, "FGSRCE") == 0) {
            if ((istat = GetNextSource(strchr(line, ' '))) < 0) {
                nll_puterr("ERROR: reading source params:");
                nll_puterr(line);
            } else
                flag_source = 1;
        }


        /* read grid params */

        if (strcmp(param, "FGGRID") == 0) {
            if ((istat = get_grid(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading grid parameters.");
            else
                flag_grid = 1;
        }



        /* read path phase code */

        if (strcmp(param, "FGPATHPHASECODE") == 0) {
            if ((istat = get_fg_type(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading fmm2grid path phase code.");
            else
                flag_path_phase_code = 1;
        }



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
    if (!flag_inpfile)
        nll_puterr("ERROR: no inputfile (FGINP) params read.");
    if (!flag_outfile)
        nll_puterr("ERROR: no outputfile (FGOUT) params read.");
    if (!flag_source)
        nll_puterr("ERROR: no source (FGSRCE) params read.");
    if (!flag_path_phase_code)
        nll_puterr("ERROR: no type (FGPATHPHASECODE) params read.");
    if (!flag_grid)
        nll_puterr("ERROR: no grid (FGGRID) params read.");

    if (!flag_trans) {
        sprintf(MsgStr, "INFO: no transformation (TRANS) params read.");
        nll_putmsg(1, MsgStr);
        Hypocenter.comment[0] = '\0';
    }

    return (flag_include * flag_control * flag_inpfile * flag_outfile * flag_grid * flag_path_phase_code - 1);
}

