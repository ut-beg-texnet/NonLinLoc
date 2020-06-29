/*
 * Copyright (C) 1999-2008 Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
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




/*   Grid2Time.c

        Program to calculate travel times for 3-D grids

 */

/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
161 Allee du Micocoulier, 06370 Mouans-Sartoux, France
tel: +33(0)493752502  e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


/*
        history:

        ver 01    22SEP1997  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

 */



/*#define EXTERN_MODE 1 */

#include "GridLib.h"

#define DEBUG_GRID2TIME 0



/*------------------------------------------------------------/ */
/* declarations for grid modes */
/*------------------------------------------------------------/ */

int grid_mode; /* grid mode - GRID_TIME, GRID_TIME_2D */

int angle_mode; /* angle mode - ANGLE_MODE_NO, ANGLE_MODE_YES */



/*------------------------------------------------------------/ */
/* declarations for various travel time calculation methods */
/*------------------------------------------------------------/ */

int tt_calc_meth; /*	travel time calc method */


/*------------------------------------------------------------/ */
/* METHOD_PODLECFD = Podvin & Lecomte Finite Diff  (TIme_3d_fs.c) */
/*			Podvin & Lecomte, Geophys.J.Intnl. 105, 271-284, 1991. */
#define METHOD_PODLECFD 1
int time_3d(GRID_FLOAT_TYPE *HS, GRID_FLOAT_TYPE *T, int NX, int NY, int NZ, GRID_FLOAT_TYPE XS, GRID_FLOAT_TYPE YS, GRID_FLOAT_TYPE ZS, GRID_FLOAT_TYPE HS_EPS_INIT, int MSG);
/*
int time_3d(HS,T,NX,NY,NZ,XS,YS,ZS,HS_EPS_INIT,MSG)
GRID_FLOAT_TYPE *HS,*T,HS_EPS_INIT,XS,YS,ZS;
int   NX,NY,NZ,MSG;
 */

/* Podvin & Lecomte Finite Diff parameters */

double plfd_hs_eps_init;
int plfd_message;

/*------------------------------------------------------------/ */



/*------------------------------------------------------------/ */
/* METHOD_WAVEFRONT_RAY = wavefront ray tracing (green3d)  */
/*		Programed by p.S.Lucio and G.C.Lambare  */
/*		GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP) */
/*		Lambare, Lucio and Hanyga, 1996, GJI, 125, 584-598. */
#define METHOD_WAVEFRONT_RAY 2
int get_gt_wavefront(char* line1);

/* wavefront ray tracing (green3d) parameters */

int wvfrnt_nir; /* maximum number of stored arrival */
/* (by default it always keep the first arrival  */
/* traveltime see in ajuste.f for any modification of  */
/* this default) */

int wvfrnt_npr; /* number of parameters for which to compute maps */
/* it should be less than 27 and greater or equal to 3 */
/* see in ajuste.f for any modification */

float wvfrnt_targ_orient[7][3];
/* comments from modelrese.f:  (Note: FORTRAN/C indexing) */
/* origine point of target x,y,z */
/* 	targ_orient(1,1),targ_orient(2,1),targ_orient(3,1) */
/* 	targ_orient[0][0],targ_orient[0][1],targ_orient[0][2] */
/* vector for first increment in target */
/* 	targ_orient(1,2),targ_orient(2,2),targ_orient(3,2) */
/* 	targ_orient[1][0],targ_orient[1][1],targ_orient[1][2] */
/* vector for second increment in target */
/* 	targ_orient(1,3),targ_orient(2,3),targ_orient(3,3) */
/* 	targ_orient[2][0],targ_orient[2][1],targ_orient[2][2] */
/* vector for third increment in target */
/* 	targ_orient(1,4),targ_orient(2,4),targ_orient(3,4) */
/* 	targ_orient[3][0],targ_orient[3][1],targ_orient[3][2] */

int *wvfrnt_imap;
/* imap(nzr,nxr,nyr) */
/*	total number of arrivals by points */

float wvfrnt_fi1min, wvfrnt_fi2min, wvfrnt_fi1max, wvfrnt_fi2max;
/* initial angular aperture (in radian) */
/* 	see (Lucio et al, 1996) */

float wvfrnt_dxmin2, wvfrnt_dpmin2;
/* dxmax**2,dpmax**2  see (Lucio et al, 1996) */
/* 	it gives some idea of the precision of the ray */
/* 	field sampling (in m**2 and (s/m)**2) */
/* 	example (5**2 , (5E-E06)**2) */
/*	increase this number in order to go faster */

float wvfrnt_dtemps;
/* travel time step (in s)  */

/*------------------------------------------------------------/ */


/*------------------------------------------------------------/ */
/* END declarations for various travel time calculation methods */
/*------------------------------------------------------------/ */



/* globals  */


char fn_gt_input[MAXLINE], fn_gt_output[MAXLINE];
int iSwapBytesOnInput;


/* function declarations */

int ReadGrid2TimeInput(FILE*);
int get_gt_files(char*);
int get_grid_mode(char*);
int get_gt_plfd(char*);
int GenTimeGrid(GridDesc*, SourceDesc*, GridDesc*, char*);
int GenAngleGrid(GridDesc*, SourceDesc*, GridDesc*, int);
void InitTimeGrid(GridDesc*, GridDesc*);
int RunGreen3d(GridDesc*, SourceDesc*, GridDesc*, char*);



/*** program to generate  3-D travel time grid */

#define PNAME  "Grid2Time"

#define NARGS 2

int main(int argc, char *argv[]) {

    int istat;
    int nsrce;

    int ix, iy, iz, iymax, izmax, iystep, izstep;

    char fn_model[MAXLINE];
    FILE *fp_model_grid, *fp_model_hdr;

    GridDesc mod_grid, time_grid, angle_grid;


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

    if ((istat = ReadGrid2TimeInput(fp_control)) < 0) {
        exit(EXIT_ERROR_FILEIO);
    }



    /* convert source location coordinates  */

    istat = ConvertSourceLoc(0, Source, NumSources, 1, 1);


    /* open model file and read header */

    sprintf(fn_model, "%s.mod", fn_gt_input);
    if ((istat = OpenGrid3dFile(fn_model, &fp_model_grid, &fp_model_hdr,
            &mod_grid, " ", NULL, iSwapBytesOnInput)) < 0) {
        CloseGrid3dFile(&fp_model_grid, &fp_model_hdr);
        nll_puterr2("ERROR: cannot open model grid", fn_model);
        exit(EXIT_ERROR_FILEIO);
    }

    /* check grid x size */

    if (grid_mode == GRID_TIME_2D && mod_grid.numx != 2) {
        nll_puterr(
                "ERROR: grid xNum must be 2 for gridMode GRID2D");
        exit(EXIT_ERROR_TTIME);
    } else if (grid_mode == GRID_TIME && mod_grid.numx <= 2) {
        sprintf(MsgStr,
                "WARNING: grid xNum = %d is very small for gridMode GRID3D",
                mod_grid.numx);
        nll_putmsg(1, MsgStr);
    }


    /* initialize 3D grids */

    InitTimeGrid(&time_grid, &mod_grid);

    if (angle_mode == ANGLE_MODE_YES) {
        if (grid_mode == GRID_TIME_2D)
            DuplicateGrid(&angle_grid, &time_grid, "ANGLE2D");
        else
            DuplicateGrid(&angle_grid, &time_grid, "ANGLE");
    } else if (angle_mode == ANGLE_MODE_INCLINATION) {
        if (grid_mode == GRID_TIME_2D)
            DuplicateGrid(&angle_grid, &time_grid, "INCLINATION2D");
        else
            DuplicateGrid(&angle_grid, &time_grid, "INCLINATION");
    }

    /* allocate model grids */
    mod_grid.buffer = AllocateGrid(&mod_grid);
    if (mod_grid.buffer == NULL) {
        nll_puterr(
                "ERROR: allocating memory for 3D slowness grid buffer.");
        exit(EXIT_ERROR_MEMORY);
    }

    /* create grid array access pointers */
    mod_grid.array = CreateGridArray(&mod_grid);
    if (mod_grid.array == NULL) {
        nll_puterr(
                "ERROR: creating array for accessing 3D slowness grid buffer.");
        exit(EXIT_ERROR_MEMORY);
    }



    /* read vel/slowness model grid */

    if ((istat =
            ReadGrid3dBuf(&mod_grid, fp_model_grid)) < 0) {
        nll_puterr("ERROR: reading vel/slowness model grid from disk.");
        exit(EXIT_ERROR_IO);
    }


    /* check model grid */
    if (1) {
        ix = mod_grid.numx / 2;
        iymax = mod_grid.numy - 1;
        iystep = iymax / 10;
        if (iystep < 1) iystep = 1;
        izmax = mod_grid.numz - 1;
        izstep = izmax / 10;
        if (izstep < 1) izstep = 1;
        sprintf(MsgStr, "Sample of Model Grid: X=%d  Y=0,%d,%d  Z=0,%d,%d",
                ix, iymax, iystep, izmax, izstep);
        nll_putmsg(2, MsgStr);
        if (message_flag >= 2) {
            for (iz = 0; iz < izmax; iz += izstep) {
                for (iy = 0; iy < iymax; iy += iystep)
                    fprintf(stdout, "%.2e ", ((GRID_FLOAT_TYPE ***) mod_grid.array)[0][iy][iz]);
            }
            fprintf(stdout, "\n");
        }
        if ((istat = CheckGridArray(&mod_grid,
                VERY_LARGE_FLOAT, VERY_LARGE_FLOAT,
                -VERY_LARGE_FLOAT, -VERY_LARGE_FLOAT)) < 0) {
            nll_puterr("ERROR: invalid vel/slowness model grid.");
            exit(EXIT_ERROR_MODEL);
        }
    }

    // 20101005 AJL - added to allow retrieval of model velocity at a point (for NLL OT_STACK etc.)
    // save copy of model file in time output directory
    if ((istat = WriteGrid3dBuf(&mod_grid, NULL, fn_gt_output, "mod")) < 0) {
        nll_puterr2("ERROR: writing slowness grid to disk", fn_gt_output);
        exit(EXIT_ERROR_IO);
    }


    /* generate travel time and take-off angle grids for each source */

    for (nsrce = 0; nsrce < NumSources; nsrce++) {
        sprintf(MsgStr,
                "\nCalculating travel times for source: %s  X %.4lf  Y %.4lf  Z %.4lf (lat/lon/depth  %f  %f  %f) ...",
                (Source + nsrce)->label, (Source + nsrce)->x,
                (Source + nsrce)->y, (Source + nsrce)->z,
                (Source + nsrce)->dlat, (Source + nsrce)->dlong, (Source + nsrce)->depth
                );
        nll_putmsg(1, MsgStr);
        if ((istat = GenTimeGrid(&mod_grid, Source + nsrce, &time_grid, fn_model)) < 0)
            nll_puterr("ERROR: calculating travel times.");
        else if (angle_mode == ANGLE_MODE_YES) {
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

void InitTimeGrid(GridDesc* ptime_grid, GridDesc* pmod_grid) {
    char chr_type[MAXLINE];


    /* set grid type */
    if (grid_mode == GRID_TIME_2D)
        strcpy(chr_type, "TIME2D");
    else
        strcpy(chr_type, "TIME");



    /* duplicate grid and allocate memory */
    if (tt_calc_meth == METHOD_PODLECFD) {
        DuplicateGrid(ptime_grid, pmod_grid, chr_type);
    } else if (tt_calc_meth == METHOD_WAVEFRONT_RAY) {
        DuplicateGrid(ptime_grid, &grid_in, chr_type);
    }



    /* allocate additional grids */
    if (tt_calc_meth == METHOD_WAVEFRONT_RAY) {
        if ((wvfrnt_imap = (int *) malloc((size_t)
                (ptime_grid->numx * ptime_grid->numy * ptime_grid->numz
                * sizeof (int)))) == NULL) {
            nll_puterr(
                    "ERROR: allocating memory for 3D time grid buffer.");
            exit(EXIT_ERROR_MEMORY);
        }

    }

}

/*** function to generate travel time grid */

int GenTimeGrid(GridDesc* pmgrid, SourceDesc* psource, GridDesc* ptt_grid, char* fn_model) {

    int istat, itemp = 0;
    char filename[MAXLINE];
    double xsource, ysource, zsource;
    double vel_source;


    // Podvin & Lecomte Finite Diff parameters - use GRID_FLOAT_TYPE
    GRID_FLOAT_TYPE xsource_igrid, ysource_igrid, zsource_igrid;


    /* check grid mode, make appropriate adjustments */

    if (grid_mode == GRID_TIME_2D) {
        /* set horiz source location to grid origin */
        xsource = ptt_grid->origx;
        ysource = ptt_grid->origy;
        zsource = psource->z;
        vel_source = ReadAbsInterpGrid2d(NULL, pmgrid, ysource, zsource);
    } else {
        xsource = psource->x;
        ysource = psource->y;
        zsource = psource->z;
        /*vel_source = ReadAbsInterpGrid3d(NULL, pmgrid,
                xsource, ysource, zsource);*/
        vel_source = ReadAbsGrid3dValue(NULL, pmgrid, xsource, ysource, zsource, 1);
    }


    /* calculate source grid location */

    xsource_igrid = (xsource - ptt_grid->origx) / ptt_grid->dx;
    ysource_igrid = (ysource - ptt_grid->origy) / ptt_grid->dy;
    zsource_igrid = (zsource - ptt_grid->origz) / ptt_grid->dz;


    /* check that source in inside model grid */

    if (!IsPointInsideGrid(pmgrid, xsource, ysource, zsource)) {
        nll_puterr(
                "ERROR: Source point is not inside model grid.");
        sprintf(MsgStr,
                "Source:  GridLoc: ix=%lf iy=%lf iz=%lf",
                xsource_igrid, ysource_igrid, zsource_igrid);
        nll_putmsg(0, MsgStr);
        return (-1);
    }


    /* display velocity at source */

    if (pmgrid->type == GRID_VELOCITY)
        ;
    else if (pmgrid->type == GRID_VELOCITY_METERS)
        vel_source = vel_source / 1000.0;
    else if (pmgrid->type == GRID_SLOWNESS)
        vel_source = 1.0 / (vel_source);
    else if (pmgrid->type == GRID_SLOW_LEN)
        vel_source = 1.0 / (vel_source / pmgrid->dx);
    else if (pmgrid->type == GRID_VEL2)
        vel_source = sqrt(vel_source);
    else if (pmgrid->type == GRID_SLOW2)
        vel_source = sqrt(1.0 / vel_source);
    else if (pmgrid->type == GRID_SLOW2_METERS)
        vel_source = sqrt(1.0 / vel_source) / 1000.0;
    sprintf(MsgStr,
            "Source:  Velocity: %lf km/sec  GridLoc: ix=%lf iy=%lf iz=%lf",
            vel_source, xsource_igrid, ysource_igrid, zsource_igrid);
    nll_putmsg(1, MsgStr);



    /* generate travel time grid */

    if (tt_calc_meth == METHOD_PODLECFD) {
        /* check things */
        if (pmgrid->type != GRID_SLOW_LEN) {
            nll_puterr(
                    "ERROR: Podvin-Lecomte algorithm requires SLOW_LEN grid.");
            return (-1);
        }
        if (pmgrid->dx != pmgrid->dy || pmgrid->dx != pmgrid->dz) {
            nll_puterr(
                    "ERROR: Podvin-Lecomte algorithm requires cubic grid, i.e. dx=dy=dz.");
            return (-1);
        }

        if (DEBUG_GRID2TIME) {
            fprintf(stdout, "ptt_grid->numx %d  ptt_grid->numy %d  ptt_grid->numz %d  xsource_igrid %.3f  ysource_igrid %.3f  zsource_igrid %.3f  plfd_hs_eps_init %.2e\n",
                    ptt_grid->numx, ptt_grid->numy, ptt_grid->numz, xsource_igrid, ysource_igrid, zsource_igrid, plfd_hs_eps_init);
        }

        /* run Podvin-Lecomte algorithm */
        istat = time_3d(pmgrid->buffer, ptt_grid->buffer,
                ptt_grid->numx, ptt_grid->numy, ptt_grid->numz,
                // AJL 2010 - following seems to avoid all zero travel time grids in some cases ! (GTSRCE  TARGET_0   XYZ  0.0 0.0 2.200000 0.0 ??  Mac OS X ??)
                //xsource_igrid + ptt_grid->dx / 2.0, ysource_igrid, zsource_igrid,
                (GRID_FLOAT_TYPE) xsource_igrid, (GRID_FLOAT_TYPE) ysource_igrid, (GRID_FLOAT_TYPE) zsource_igrid,
                (GRID_FLOAT_TYPE) plfd_hs_eps_init, plfd_message);
        if (DEBUG_GRID2TIME) {
            fprintf(stdout, "time_3d returned: %d\n", istat);
        }
        if (istat)
            return (-1);


    } else if (tt_calc_meth == METHOD_WAVEFRONT_RAY) {
        /* check things */
        if (pmgrid->type != GRID_VELOCITY_METERS) {
            nll_puterr(
                    "ERROR: Wavefront ray tracing (green3d) algorithm requires VELOCITY_METERS (meters/sec) grid.");
            return (-1);
        }
        if (grid_mode != GRID_TIME) {
            nll_puterr(
                    "ERROR: Wavefront ray tracing (green3d) algorithm requires grid mode GRID3D.");
            return (-1);
        }

        /* run wavefront algorithm */
        if ((istat = RunGreen3d(pmgrid, psource, ptt_grid, fn_model))
                < 0)
            return (-1);

    } else {
        nll_puterr("ERROR: unrecognized travel time calculation method");
        return (-1);
    }


    /* check time grid */
    if (DEBUG_GRID2TIME) {
        int ix, iy, iz, iymax, izmax, iystep, izstep;
        ix = ptt_grid->numx / 2;
        iymax = ptt_grid->numy - 1;
        iystep = iymax / 10;
        if (iystep < 1) iystep = 1;
        izmax = ptt_grid->numz - 1;
        izstep = izmax / 10;
        if (izstep < 1) izstep = 1;
        sprintf(MsgStr, "Sample of Time Grid: X=%d  Y=0,%d,%d  Z=0,%d,%d",
                ix, iymax, iystep, izmax, izstep);
        nll_putmsg(2, MsgStr);
        if (message_flag >= 2) {
            for (iz = 0; iz < izmax; iz += izstep) {
                for (iy = 0; iy < iymax; iy += iystep)
                    fprintf(stdout, "%.2e ", ((GRID_FLOAT_TYPE ***) ptt_grid->array)[0][iy][iz]);
            }
            fprintf(stdout, "\n");
        }
        if ((istat = CheckGridArray(ptt_grid,
                VERY_LARGE_FLOAT, VERY_LARGE_FLOAT,
                -VERY_LARGE_FLOAT, -VERY_LARGE_FLOAT)) < 0) {
            nll_puterr("ERROR: invalid vel/slowness model grid.");
            exit(EXIT_ERROR_MODEL);
        }
    }

    // output ascii table of grid
#define WRITE_ASCII_GRID 0
#define USE_COORDS_LATLON 1
    if (WRITE_ASCII_GRID) {
        sprintf(filename, "%s.%s.time.csv", fn_gt_output, psource->label);
        FILE *fp_ascii_grid = NULL;
        if ((fp_ascii_grid = fopen(filename, "w")) == NULL) {
            nll_puterr2("ERROR: opening ascii grid output file: %s", filename);
        } else {
            int ix, iy, iz, ixmax, iymax, izmax, ixstep, iystep, izstep;
            ix = iy = iz = 0;
            ixmax = ptt_grid->numx - 1;
            ixstep = 1;
            iymax = ptt_grid->numy - 1;
            iystep = 1;
            izmax = ptt_grid->numz - 1;
            izstep = 1;
            double xloc, yloc;
            double xval = ptt_grid->origx;
            for (ix = 0; ix < ixmax; ix += ixstep) {
                double yval = ptt_grid->origy;
                for (iy = 0; iy < iymax; iy += iystep) {
                    if (USE_COORDS_LATLON) {
                        rect2latlon(0, xval, yval, &yloc, &xloc);
                    } else {
                        xloc = xval;
                        yloc = yval;
                    }
                    double zloc = ptt_grid->origz;
                    for (iz = 0; iz < izmax; iz += izstep) {
                        fprintf(fp_ascii_grid, "%g %g %g %g\n", xloc, yloc, zloc, ((GRID_FLOAT_TYPE ***) ptt_grid->array)[ix][iy][iz]);
                        zloc += ptt_grid->dz;
                    }
                    yval += ptt_grid->dy;
                }
                xval += ptt_grid->dx;
            }
            fclose(fp_ascii_grid);
            sprintf(MsgStr, "Ascii grid output written to file: %s", filename);
            nll_putmsg(1, MsgStr);
        }
    }


    /* save time grid to disk */

    sprintf(filename, "%s.%s", fn_gt_output, psource->label);
    sprintf(MsgStr,
            "Finished calculation, time grid output files: %s.*", filename);
    nll_putmsg(1, MsgStr);
    /* need only ix=0 sheet for 2D grids */
    if (grid_mode == GRID_TIME_2D) {
        itemp = ptt_grid->numx;
        ptt_grid->numx = 1;
    }
    istat = WriteGrid3dBuf(ptt_grid, psource, filename, "time");
    if (grid_mode == GRID_TIME_2D)
        ptt_grid->numx = itemp;
    if (istat < 0) {
        nll_puterr("ERROR: writing slowness grid to disk.");
        return (-1);
    }


    return (0);

}

/* function to run wavefront-ray algorithm green3d using a system call */

int RunGreen3d(GridDesc* pmgrid, SourceDesc* psource, GridDesc* ptt_grid,
        char* fn_model) {

    int istat;

    FILE *fp_green_in, *fp_green_io;
    char fn_green_in[] = "green3d.input";
    char fn_numarrival[] = "green3d.numarrival";
    char fn_temps[] = "green3d.temps";
    char fn_ampl[] = "green3d.ampl";

    char system_str[MAXLINE];

    double km2m = 1000.0;
    float dxm, dym, dzm;
    float xsrc[3];


    /* set parameters */

    dxm = km2m * pmgrid->dx;
    dym = km2m * pmgrid->dy;
    dzm = km2m * pmgrid->dz;

    wvfrnt_targ_orient[0][0] = km2m * ptt_grid->origy;
    wvfrnt_targ_orient[0][1] = km2m * ptt_grid->origx;
    wvfrnt_targ_orient[0][2] = km2m * ptt_grid->origz;
    wvfrnt_targ_orient[1][0] = km2m * ptt_grid->dy;
    wvfrnt_targ_orient[1][1] = 0.0;
    wvfrnt_targ_orient[1][2] = 0.0;
    wvfrnt_targ_orient[2][0] = 0.0;
    wvfrnt_targ_orient[2][1] = km2m * ptt_grid->dx;
    wvfrnt_targ_orient[2][2] = 0.0;
    wvfrnt_targ_orient[3][0] = 0.0;
    wvfrnt_targ_orient[3][1] = 0.0;
    wvfrnt_targ_orient[3][2] = km2m * ptt_grid->dz;

    xsrc[0] = km2m * psource->y;
    xsrc[1] = km2m * psource->x;
    xsrc[2] = km2m * psource->z;



    /* write green3d input file */
    /* Note: green3a map storage is FORTRAN (z,x,y) = C (y,x,z) */
    /*   but we use C (x,y,z), thus we must exchange some */
    /*   x and y arguments to green3a. */

    if ((fp_green_in = fopen(fn_green_in, "w")) == NULL) {
        nll_puterr("ERROR: opening green3d input file.");
        return (-1);
    }

    fprintf(fp_green_in,
            "%d %d %d dimension of the velocity model nx,ny,nz\n",
            pmgrid->numy, pmgrid->numx, pmgrid->numz);
    fprintf(fp_green_in,
            "%d %d %d %d dimension of the reservoir nxr,nyr,nzr nber of arrivals nir\n",
            ptt_grid->numy, ptt_grid->numx, ptt_grid->numz, wvfrnt_nir);
    fprintf(fp_green_in, "%d number_parameter_maps\n", wvfrnt_npr);
    fprintf(fp_green_in, "%f %f xmin, dx for the velocity field\n",
            km2m * pmgrid->origy, dym);
    fprintf(fp_green_in, "%f %f ymin, dy for the velocity field\n",
            km2m * pmgrid->origx, dxm);
    fprintf(fp_green_in, "%f %f zmin, dz for the velocity field\n",
            km2m * pmgrid->origz, dzm);
    fprintf(fp_green_in, "%s.buf\n", fn_model);
    fprintf(fp_green_in, "%f %f %f xmin,ymin,zmin (target)\n",
            wvfrnt_targ_orient[0][0], wvfrnt_targ_orient[0][1],
            wvfrnt_targ_orient[0][2]);
    fprintf(fp_green_in, "%f %f %f dx1,dy1,dz1    (target)\n",
            wvfrnt_targ_orient[1][0], wvfrnt_targ_orient[1][1],
            wvfrnt_targ_orient[1][2]);
    fprintf(fp_green_in, "%f %f %f dx2,dy2,dz2    (target)\n",
            wvfrnt_targ_orient[2][0], wvfrnt_targ_orient[2][1],
            wvfrnt_targ_orient[2][2]);
    fprintf(fp_green_in, "%f %f %f dx3,dy3,dz3    (target)\n",
            wvfrnt_targ_orient[3][0], wvfrnt_targ_orient[3][1],
            wvfrnt_targ_orient[3][2]);
    fprintf(fp_green_in, "%f %f %f xs, ys, zs position of the source\n",
            xsrc[0], xsrc[1], xsrc[2]);
    fprintf(fp_green_in,
            "%f %f angle / vertical direction x fi1min,fi1max\n",
            wvfrnt_fi1min, wvfrnt_fi2min);
    fprintf(fp_green_in,
            "%f %f angle / vertical direction x fi1min,fi1max\n",
            wvfrnt_fi1max, wvfrnt_fi2max);
    fprintf(fp_green_in, "%f precision in x dxmin\n", wvfrnt_dxmin2);
    fprintf(fp_green_in, "%f precision in p dpmin\n", wvfrnt_dpmin2);
    fprintf(fp_green_in, "%f step of the Runge-Kutta scheme\n",
            wvfrnt_dtemps);
    fprintf(fp_green_in, "%s\n", fn_numarrival);
    fprintf(fp_green_in, "%s\n", fn_temps);
    fprintf(fp_green_in, "%s\n", fn_ampl);
    fprintf(fp_green_in, "%s\n", "green3d.px");
    fprintf(fp_green_in, "%s\n", "green3d.py");
    fprintf(fp_green_in, "%s\n", "green3d.pz");
    fprintf(fp_green_in, "%s\n", "green3d.phi1");
    fprintf(fp_green_in, "%s\n", "green3d.phi2");
    fprintf(fp_green_in, "%s\n", "green3d.dx-dfi1");
    fprintf(fp_green_in, "%s\n", "green3d.dy-dfi1");
    fprintf(fp_green_in, "%s\n", "green3d.dz-dfi1");
    fprintf(fp_green_in, "%s\n", "green3d.dpx-dfi1");
    fprintf(fp_green_in, "%s\n", "green3d.dpy-dfi1");
    fprintf(fp_green_in, "%s\n", "green3d.dpz-dfi1");
    fprintf(fp_green_in, "%s\n", "green3d.dpx-dfi2");
    fprintf(fp_green_in, "%s\n", "green3d.dpy-dfi2");
    fprintf(fp_green_in, "%s\n", "green3d.dpz-dfi2");
    fprintf(fp_green_in, "%s\n", "green3d.dpx-dxs");
    fprintf(fp_green_in, "%s\n", "green3d.dpy-dxs");
    fprintf(fp_green_in, "%s\n", "green3d.dpz-dxs");
    fprintf(fp_green_in, "%s\n", "green3d.dpx-dys");
    fprintf(fp_green_in, "%s\n", "green3d.dpy-dys");
    fprintf(fp_green_in, "%s\n", "green3d.dpz-dys");

    fclose(fp_green_in);


    /* run green3d */

    sprintf(system_str, "green3 < %s", fn_green_in);
    sprintf(MsgStr, "Calling green3d [%s] ...", system_str);
    nll_putmsg(2, MsgStr);
    istat = system(system_str);
    sprintf(MsgStr, "Return from green3d - return val is %d", istat);
    nll_putmsg(2, MsgStr);


    /* read wavefront/ray travel time grid and write new grid */
    if ((fp_green_io = fopen(fn_temps, "r")) == NULL) {
        nll_puterr("ERROR: opening wavefront/ray travel time grid file.");
        return (-1);
    }
    if ((istat =
            ReadGrid3dBuf(ptt_grid, fp_green_io)) < 0) {
        nll_puterr(
                "ERROR: reading wavefront/ray travel time grid from disk.");
        return (-1);
    }
    fclose(fp_green_io);
    if ((istat = CheckGridArray(ptt_grid, 1.0e8, -1.0,
            -VERY_LARGE_FLOAT, -VERY_LARGE_FLOAT)) < 0) {
        nll_putmsg(1,
                "WARNING: invalid or incomplete wavefront/ray grid.");
    }


    return (0);


}

/*** function to read input file */

int ReadGrid2TimeInput(FILE* fp_input) {
    int istat, iscan;
    char param[MAXLINE], *pchr;
    char line[4 * MAXLINE], *fgets_return;

    int flag_control = 0, flag_outfile = 0, flag_source = 0,
            flag_grid = 0, flag_plfd = 0, flag_wavefront = 0,
            flag_trans = 0, flag_grid_mode = 0;
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

        if (strcmp(param, "GTMODE") == 0) {
            if ((istat = get_grid_mode(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading Grid2Time grid mode.");
            else
                flag_grid_mode = 1;
        }

        /* read file names */

        if (strcmp(param, "GTFILES") == 0) {
            if ((istat = get_gt_files(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading Grid2Time file names.");
            else
                flag_outfile = 1;
        }

        /* read source params */

        if (strcmp(param, "GTSRCE") == 0) {
            if ((istat = GetNextSource(strchr(line, ' '))) < 0) {
                nll_puterr("ERROR: reading source params:");
                nll_puterr(line);
            } else
                flag_source = 1;
        }


        if (strcmp(param, "GTGRID") == 0) {
            if ((istat = get_grid(strchr(line, ' '))) < 0)
                fprintf(stderr,
                    "ERROR: reading grid parameters.");
            else
                flag_grid = 1;
        }


        /* read PodLec FD params */

        if (strcmp(param, "GT_PLFD") == 0) {
            if ((istat = get_gt_plfd(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading Podvin-Lecomte params.");
            else
                flag_plfd = 1;
        }


        /* read Wavefront params */

        if (strcmp(param, "GT_WAVEFRONT_RAY") == 0) {
            if ((istat = get_gt_wavefront(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading wavefront-ray params.");
            else
                flag_wavefront = 1;
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
    if (!flag_outfile)
        nll_puterr("ERROR: no file (GTFILES) params read.");
    if (!flag_grid_mode)
        nll_puterr("ERROR: no grid mode (GTMODE) params read.");
    if (!flag_source)
        nll_puterr("ERROR: no source (GTSRCE) params read.");
    if (!flag_plfd && !flag_wavefront)
        nll_puterr(
            "ERROR: no Travel Time method (GT_PLFD, GT_WAVEFRONT_RAY) params read.");
    if (flag_plfd + flag_wavefront > 1)
        nll_puterr(
            "ERROR: too many Travel Time methods (GT_PLFD, GT_WAVEFRONT_RAY) read.");
    if (flag_wavefront && !flag_grid)
        nll_puterr("ERROR: no grid (GTGRID) params read.");
    if (flag_plfd && flag_grid)
        nll_putmsg(2,
            "WARNING: grid (GTGRID) params ignored with Podvin-Lecompte FD method.  Podvin-Lecompte FD method reproduces model grid dimensions");
    if (!flag_trans)
        nll_puterr("ERROR: no transformation (TRANS) params read.");


    return (flag_include * flag_control * flag_outfile * flag_grid_mode
            * flag_source * flag_trans
            * (flag_plfd || (flag_wavefront && flag_grid))
            - 1);
}

/*** function to read output file name ***/

int get_gt_files(char* line1) {
    int istat;
    char waveType[12];

    istat = sscanf(line1, "%s %s %s %d", fn_gt_input, fn_gt_output,
            waveType, &iSwapBytesOnInput);
    if (istat < 4)
        iSwapBytesOnInput = 0;

    strcat(strcat(fn_gt_input, "."), waveType);
    strcat(strcat(fn_gt_output, "."), waveType);

    sprintf(MsgStr,
            "Grid2Time GTFILES:  Input: %s.*  Output: %s.*  wavetype: %s.*  iSwapBytesOnInput: %d",
            fn_gt_input, fn_gt_output, waveType, iSwapBytesOnInput);
    nll_putmsg(3, MsgStr);

    return (0);
}

/*** function to read grid mode params ***/

int get_grid_mode(char* line1) {
    char str_grid_mode[MAXLINE];
    char str_angle_mode[MAXLINE];


    sscanf(line1, "%s %s", str_grid_mode, str_angle_mode);

    sprintf(MsgStr, "Grid2Time GTMODE:  %s  %s",
            str_grid_mode, str_angle_mode);
    nll_putmsg(3, MsgStr);

    if (strcmp(str_grid_mode, "GRID3D") == 0)
        grid_mode = GRID_TIME;
    else if (strcmp(str_grid_mode, "GRID2D") == 0)
        grid_mode = GRID_TIME_2D;
    else {
        grid_mode = GRID_UNDEF;
        nll_puterr("ERROR: unrecognized grid mode");
        return (-1);
    }

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

/*** function to read Podvin-Lecompte FD params ***/

int get_gt_plfd(char* line1) {
    int istat, ierr;

    istat = sscanf(line1, "%lf %d", &plfd_hs_eps_init, &plfd_message);

    sprintf(MsgStr, "Grid2Time GT_PLFD: hs_eps_init %f  message_flag %d",
            plfd_hs_eps_init, plfd_message);
    nll_putmsg(3, MsgStr);

    ierr = 0;
    if (checkRangeInt("GT_PLFD", "message_flag", plfd_message, 1, 0, 1, 2) != 0)
        ierr = -1;
    if (checkRangeInt("GT_PLFD", "hs_eps_init", plfd_hs_eps_init, 1, 0.0, 0, 0.0) != 0)
        ierr = -1;

    tt_calc_meth = METHOD_PODLECFD;
    sprintf(MsgStr,
            "  (Method is Podvin-Lecompte Finite-Differences)");
    nll_putmsg(3, MsgStr);

    if (ierr < 0 || istat != 2)
        return (-1);

    return (0);

}

/*** function to read Wavefront params ***/

int get_gt_wavefront(char* line1) {
    int istat;


    istat = sscanf(line1, "%d %d  %f %f %f %f  %f %f  %f",
            &wvfrnt_nir, &wvfrnt_npr,
            &wvfrnt_fi1min, &wvfrnt_fi2min, &wvfrnt_fi1max, &wvfrnt_fi2max,
            &wvfrnt_dxmin2, &wvfrnt_dpmin2, &wvfrnt_dtemps
            );

    if (wvfrnt_nir < 1) {
        wvfrnt_nir = 1;
        nll_putmsg(1,
                "WARNING: max_number_stored_arrivals invalid, reset to 1");
    }
    if (wvfrnt_npr < 3) {
        wvfrnt_npr = 3;
        nll_putmsg(1,
                "WARNING: number_parameter_maps invalid, reset to 3");
    }

    sprintf(MsgStr,
            "Grid2Time GT_WAVEFRONT_RAY: Max_number_stored_arrivals %d  Number_parameter_maps %d",
            wvfrnt_nir, wvfrnt_npr);
    nll_putmsg(3, MsgStr);

    sprintf(MsgStr,
            "  Initial_angular_aperture:  fi1min %f  fi2min %f  fi1max %f fi2max %lf",
            wvfrnt_fi1min, wvfrnt_fi2min, wvfrnt_fi1max, wvfrnt_fi2max);
    nll_putmsg(3, MsgStr);
    sprintf(MsgStr,
            "  Precision of ray field sampling:  dxmax**2 %f  dpmax**2 %f  Travel_time_step %f",
            wvfrnt_dxmin2, wvfrnt_dpmin2, wvfrnt_dtemps);
    nll_putmsg(3, MsgStr);


    tt_calc_meth = METHOD_WAVEFRONT_RAY;
    sprintf(MsgStr,
            "  (Method is wavefront ray tracing - green3d)");
    nll_putmsg(3, MsgStr);

    if (istat != 9)
        return (-1);

    return (0);

}

/*** function to generate take-off angle grid */

int GenAngleGrid(GridDesc* ptgrid, SourceDesc* psource, GridDesc* pagrid, int angle_mode) {

    int istat, itemp = 0;
    char filename[MAXLINE];

    double xsource, ysource, zsource;



    /* check grid mode, make appropriate adjustments */

    if (grid_mode == GRID_TIME_2D) {
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
        if ((istat = CalcAnglesGradient(ptgrid, pagrid, angle_mode, grid_mode)) < 0)
            return (-1);

    }




    /* save angle grid to disk */

    sprintf(filename, "%s.%s", fn_gt_output, psource->label);
    sprintf(MsgStr,
            "Finished calculation, take-off angles grid output files: %s.*",
            filename);
    nll_putmsg(1, MsgStr);
    /* need only ix=0 sheet for 2D grids */
    if (grid_mode == GRID_TIME_2D) {
        itemp = pagrid->numx;
        pagrid->numx = 1;
    }
    if (angle_mode == ANGLE_MODE_YES)
        istat = WriteGrid3dBuf(pagrid, psource, filename, "angle");
    else if (angle_mode == ANGLE_MODE_INCLINATION)
        istat = WriteGrid3dBuf(pagrid, psource, filename, "inclination");
    if (grid_mode == GRID_TIME_2D)
        pagrid->numx = itemp;
    if (istat < 0) {
        nll_puterr("ERROR: writing take-off angles grid to disk.");
        return (-1);
    }


    return (0);

}

/*------------------------------------------------------------/ */
/* Anthony Lomax           | email: lomax@faille.unice.fr     / */
/* UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax / */
/* 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        / */
/* 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        / */
/*------------------------------------------------------------/ */

