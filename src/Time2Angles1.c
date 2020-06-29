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




/*   Time2Angles.c

        Program to calculate take-off angle grids from travel-time grids

 */

/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
161 Allee du Micocoulier, 06370 Mouans-Sartoux, France
tel: +33(0)493752502  e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


/*
        history:

        ver 01    24Jul2012  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

 */



/*

# =============================================================================
# =============================================================================
# Time2Angles control file statements
# =============================================================================
#
#

# input, output filename root
# (TAFILES <input time file root> <output angles file root> wave_type (P, S))
#
TAFILES  /temp/nlloc_tmp/eth_topo/time_3d/ch  /temp/nlloc_tmp/eth_topo/angles_3d/ch P

# angles grid mode
# (TAMODE angle_mode)
#    (char[])   angle_mode (ANGLES_YES, ANGLES_NO, ANGLES_INCLINATION)
#
TAMODE ANGLES_YES

# source description (multiple sources can be specified)
# (TASRCE (see #GTSRCE)
#
# no entries - use GTSRCE sources

#
#
# =============================================================================
# END of Time2Angles control file statements
# =============================================================================
# =============================================================================

 */



#include "GridLib.h"



/*------------------------------------------------------------/ */
/* declarations for grid modes */
/*------------------------------------------------------------/ */

int angle_mode; /* angle mode - ANGLE_MODE_NO, ANGLE_MODE_YES */



/* globals  */


char fn_at_input[MAXLINE], fn_at_output[MAXLINE];
int iSwapBytesOnInput;


/* function declarations */

int ReadTime2AnglesInput(FILE*);
int get_ta_files(char*);
int get_grid_mode(char*);
void InitAnglesGrid(GridDesc* pangle_grid, GridDesc* ptime_grid);
int GenAngleGrid(GridDesc*, SourceDesc*, GridDesc*, int);



/*** program to generate  3-D travel time grid */

#define PNAME  "Time2Angles"

#define NARGS 2

int main(int argc, char *argv[]) {

    int istat;
    int nsrce;

    int ix, iy, iz, iymax, izmax, iystep, izstep;

    char fn_time[MAXLINE];
    FILE *fp_time_grid, *fp_time_hdr;

    GridDesc time_grid, angle_grid;


    /* set program name */
    strcpy(prog_name, PNAME);

    /* check command line for correct usage */

    if (argc != NARGS) {
        disp_usage(prog_name, "<control file>");
        exit(EXIT_ERROR_USAGE);
    }



    /* set constants */

    SetConstants();
    prog_mode_3d = 1;
    NumSources = 0;



    /* read control file */

    strcpy(fn_control, argv[1]);
    if ((fp_control = fopen(fn_control, "r")) == NULL) {
        nll_puterr("ERROR: opening control file.");
        exit(EXIT_ERROR_FILEIO);
    }

    if ((istat = ReadTime2AnglesInput(fp_control)) < 0) {
        exit(EXIT_ERROR_FILEIO);
    }

    /* convert source location coordinates  */

    istat = ConvertSourceLoc(0, Source, NumSources, 1, 1);


    /* generate travel time and take-off angle grids for each source */

    for (nsrce = 0; nsrce < NumSources; nsrce++) {
        sprintf(MsgStr,
                "\nCalculating angles for source: %s  X %.4lf  Y %.4lf  Z %.4lf (lat/lon/depth  %f  %f  %f) ...",
                (Source + nsrce)->label, (Source + nsrce)->x,
                (Source + nsrce)->y, (Source + nsrce)->z,
                (Source + nsrce)->dlat, (Source + nsrce)->dlong, (Source + nsrce)->depth
                );
        nll_putmsg(1, MsgStr);

        /* open time file and read header */

        sprintf(fn_time, "%s.%s.time", fn_at_input, (Source + nsrce)->label);
        if ((istat = OpenGrid3dFile(fn_time, &fp_time_grid, &fp_time_hdr,
                &time_grid, " ", NULL, iSwapBytesOnInput)) < 0) {
            CloseGrid3dFile(&time_grid, &fp_time_grid, &fp_time_hdr);
            nll_puterr2("ERROR: cannot open time grid", fn_time);
            exit(EXIT_ERROR_FILEIO);
        }

        /* initialize 3D grids */

        InitAnglesGrid(&angle_grid, &time_grid);

        if (angle_mode == ANGLE_MODE_YES) {
            if (time_grid.type == GRID_TIME_2D)
                DuplicateGrid(&angle_grid, &time_grid, "ANGLE2D");
            else
                DuplicateGrid(&angle_grid, &time_grid, "ANGLE");
        } else if (angle_mode == ANGLE_MODE_INCLINATION) {
            if (time_grid.type == GRID_TIME_2D)
                DuplicateGrid(&angle_grid, &time_grid, "INCLINATION2D");
            else
                DuplicateGrid(&angle_grid, &time_grid, "INCLINATION");
        }

        /* allocate time grids */
        time_grid.buffer = AllocateGrid(&time_grid);
        if (time_grid.buffer == NULL) {
            nll_puterr(
                    "ERROR: allocating memory for 3D slowness grid buffer.");
            exit(EXIT_ERROR_MEMORY);
        }

        /* create grid array access pointers */
        time_grid.array = CreateGridArray(&time_grid);
        if (time_grid.array == NULL) {
            nll_puterr(
                    "ERROR: creating array for accessing 3D slowness grid buffer.");
            exit(EXIT_ERROR_MEMORY);
        }



        /* read time grid */

        if ((istat =
                ReadGrid3dBuf(&time_grid, fp_time_grid)) < 0) {
            nll_puterr("ERROR: reading time grid from disk.");
            exit(EXIT_ERROR_IO);
        }


        /* check time grid */
        if (1) {
            ix = time_grid.numx / 2;
            iymax = time_grid.numy - 1;
            iystep = iymax / 10;
            if (iystep < 1) iystep = 1;
            izmax = time_grid.numz - 1;
            izstep = izmax / 10;
            if (izstep < 1) izstep = 1;
            sprintf(MsgStr, "Sample of Time Grid: X=%d  Y=0,%d,%d  Z=0,%d,%d",
                    ix, iymax, iystep, izmax, izstep);
            nll_putmsg(2, MsgStr);
            if (message_flag >= 2) {
                for (iz = 0; iz < izmax; iz += izstep) {
                    for (iy = 0; iy < iymax; iy += iystep)
                        fprintf(stdout, "%.2e ", ((GRID_FLOAT_TYPE ***) time_grid.array)[0][iy][iz]);
                }
                fprintf(stdout, "\n");
            }
            if ((istat = CheckGridArray(&time_grid,
                    VERY_LARGE_FLOAT, VERY_LARGE_FLOAT,
                    -VERY_LARGE_FLOAT, -VERY_LARGE_FLOAT)) < 0) {
                nll_puterr("WARNING: time grid contains invalid or missing values.");
                //exit(EXIT_ERROR_MODEL);
            }
        }



        if (angle_mode == ANGLE_MODE_YES) {
            if ((istat = GenAngleGrid(&time_grid, Source + nsrce,
                    &angle_grid, angle_mode)) < 0)
                nll_puterr("ERROR: calculating take-off angles.");
        } else if (angle_mode == ANGLE_MODE_INCLINATION) {
            if ((istat = GenAngleGrid(&time_grid, Source + nsrce,
                    &angle_grid, angle_mode)) < 0)
                nll_puterr("ERROR: calculating inclination angles.");
        }
    }




    exit(EXIT_NORMAL);

}

/*** function to initialize travel time grid description */

void InitAnglesGrid(GridDesc* pangle_grid, GridDesc* ptime_grid) {

    char chr_type[MAXLINE];

    /* set grid type */
    if (ptime_grid->type == GRID_TIME_2D)
        strcpy(chr_type, "ANGLE2D");
    else
        strcpy(chr_type, "ANGLE");

    /* duplicate grid and allocate memory */
    DuplicateGrid(pangle_grid, ptime_grid, chr_type);


}

/*** function to read input file */

int ReadTime2AnglesInput(FILE* fp_input) {
    int istat, iscan;
    char param[MAXLINE], *pchr;
    char line[4 * MAXLINE], *fgets_return;

    int flag_control = 0, flag_grid_mode = 0, flag_outfile = 0, flag_source = 0, flag_trans = 0;
    int flag_include = 1;


    /* read each input line */

    while ((fgets_return = fgets(line, 4 * MAXLINE, fp_input)) != NULL
            || fp_include != NULL) {


        /* check for end of include file */

        if (fgets_return == NULL && fp_include != NULL) {
            SwapBackIncludeFP(&fp_input);
            continue;
        }


        istat = -1;

        /*read parameter line */

        if ((iscan = sscanf(line, "%s", param)) < 0)
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


        /* read grid mode names */

        if (strcmp(param, "TAMODE") == 0) {
            if ((istat = get_grid_mode(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading Time2Angles grid mode.");
            else
                flag_grid_mode = 1;
        }


        /* read file names */

        if (strcmp(param, "TAFILES") == 0) {
            if ((istat = get_ta_files(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading Time2Angles file names.");
            else
                flag_outfile = 1;
        }


        /* read source params */

        if (strcmp(param, "TASRCE") == 0 || strcmp(param, "GTSRCE") == 0) {
            if ((istat = GetNextSource(strchr(line, ' '))) < 0) {
                nll_puterr("ERROR: reading source params:");
                nll_puterr(line);
            } else
                flag_source = 1;
        }

        /*read transform params */

        if (strcmp(param, "TRANS") == 0) {
            if ((istat = get_transform(0, strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading transformation parameters.");
            else
                flag_trans = 1;
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
    if (!flag_grid_mode)
        nll_puterr("ERROR: no grid mode (TAMODE) params read.");
    if (!flag_outfile)
        nll_puterr("ERROR: no file (TAFILES) params read.");
    if (!flag_source)
        nll_puterr("ERROR: no source (TASRCE or GTSRCE) params read.");
    if (!flag_trans)
        nll_puterr("ERROR: no transformation (TRANS) params read.");


    return (flag_include * flag_control * flag_grid_mode * flag_outfile * flag_source * flag_trans - 1);
}

/*** function to read output file name ***/

int get_ta_files(char* line1) {
    int istat;
    char waveType[12];

    istat = sscanf(line1, "%s %s %s %d", fn_at_input, fn_at_output, waveType, &iSwapBytesOnInput);
    if (istat < 4)
        iSwapBytesOnInput = 0;

    strcat(strcat(fn_at_input, "."), waveType);
    strcat(strcat(fn_at_output, "."), waveType);

    sprintf(MsgStr,
            "Time2Angles TAFILES:  Input: %s.*  Output: %s.*  wavetype: %s.*  iSwapBytesOnInput: %d",
            fn_at_input, fn_at_output, waveType, iSwapBytesOnInput);
    nll_putmsg(3, MsgStr);

    return (0);
}

/*** function to read grid mode params ***/

int get_grid_mode(char* line1) {

    char str_angle_mode[MAXLINE];

    sscanf(line1, "%s", str_angle_mode);

    sprintf(MsgStr, "Time2Angles TAMODE:  %s", str_angle_mode);
    nll_putmsg(3, MsgStr);

    if (strcmp(str_angle_mode, "ANGLES_YES") == 0)
        angle_mode = ANGLE_MODE_YES;
    else if (strcmp(str_angle_mode, "ANGLES_NO") == 0)
        angle_mode = ANGLE_MODE_NO;
    else if (strcmp(str_angle_mode, "ANGLES_INCLINATION") == 0)
        angle_mode = ANGLE_MODE_INCLINATION;
    else {
        angle_mode = ANGLE_MODE_UNDEF;
        nll_puterr("ERROR: unrecognized angle mode");
        return (-1);
    }

    return (0);

}

/*** function to generate take-off angle grid */

int GenAngleGrid(GridDesc* ptgrid, SourceDesc* psource, GridDesc* pagrid, int angle_mode) {

    int istat, itemp = 0;
    char filename[MAXLINE];

    double xsource, ysource, zsource;



    /* check grid mode, make appropriate adjustments */

    if (ptgrid->type == GRID_TIME_2D) {
        /* set horiz source location to grid origin */
        xsource = pagrid->origx;
        ysource = pagrid->origy;
        zsource = psource->z;
    } else {
        xsource = psource->x;
        ysource = psource->y;
        zsource = psource->z;
    }



    /* generate angle grid */

    /*if (angle_calc_meth == ANGLE_METHOD_GRADIENT)*/
    if (1) {
        // DEBUG
        //display_grid_param(ptgrid);

        /* check things */
        if (ptgrid->type != GRID_TIME && ptgrid->type != GRID_TIME_2D) {
            nll_puterr(
                    "ERROR: Gradient take-off angle algorithm requires TIME grid.");
            return (-1);
        }
        if (ptgrid->dx != ptgrid->dy || ptgrid->dx != ptgrid->dz) {
            nll_puterr(
                    "ERROR: Gradient take-off angle algorithm requires cubic grid, i.e. dx=dy=dz.");
            return (-1);
        }

        /* run gradient take-off angle algorithm */
        if ((istat = CalcAnglesGradient(ptgrid, pagrid, angle_mode, ptgrid->type)) < 0)
            return (-1);

    }




    /* save angle grid to disk */

    sprintf(filename, "%s.%s", fn_at_output, psource->label);
    sprintf(MsgStr,
            "Finished calculation, take-off angles grid output files: %s.*",
            filename);
    nll_putmsg(1, MsgStr);
    /* need only ix=0 sheet for 2D grids */
    if (ptgrid->type == GRID_TIME_2D) {
        itemp = pagrid->numx;
        pagrid->numx = 1;
    }
    if (angle_mode == ANGLE_MODE_YES)
        istat = WriteGrid3dBuf(pagrid, psource, filename, "angle");
    else if (angle_mode == ANGLE_MODE_INCLINATION)
        istat = WriteGrid3dBuf(pagrid, psource, filename, "inclination");
    if (ptgrid->type == GRID_TIME_2D)
        pagrid->numx = itemp;
    if (istat < 0) {
        nll_puterr("ERROR: writing take-off angles grid to disk.");
        return (-1);
    }


    return (0);

}


