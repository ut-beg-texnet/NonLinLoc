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


/*   interface2fmm.c
 * Program to convert pseudo-SimulPS format interface depth data files to
 * the ANU-FMM format interfaces.in file.
 * See the ANU-FMM multi-stage 3D fast marching code
 * (http://rses.anu.edu.au/seismology/soft/fmmcode/)
 */


/*
        history:

        ver 01    16May2012  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

 */




#define PNAME  "interface2fmm.c"

#include "GridLib.h"
#include "velmod.h"

typedef struct {
    float *buffer; /* buffer (contiguous floats) */
    float ***array; /* array access to buffer */
    int numx; /* Number of nodes in x coord */
    int numy; /* Number of nodes in y coord */
    int numz; /* Number of nodes in z coord */
    double origx; /* X orig (km) */
    double origy; /* Y orig (km) */
    double origz; /* V orig (km) */
    float dx; /* distance of nodes in x coord (km) */
    float dy; /* distance of nodes in y coord (km) */
    float dz; /* distance of nodes in z coord (km) */
    float *deltax; /* array with x coords of nodes (km) */
    float *deltay; /* array with y coords of nodes (km) */
    float *deltaz; /* array with z coords of nodes (km) */
    int type; /* type of velocity model [TYPESIMUL/TYPEFDTOMO] */
} VelModel;

#define TYPESIMUL  0


/* globals  */

char fn_vg_input[MAXLINE];
char fn_if_output[MAXLINE];


/* wave type (P, S, ...) for vel grids */
#define MAX_NUM_WAVE_TYPES 10
char WaveType[MAX_NUM_WAVE_TYPES][12];


/* function declarations */

int ReadInterface2FmmInput(FILE*, VelModel*);
int ReadVelModel(VelModel*, GridDesc*);
int AllocateVelModel(VelModel*);
void FreeVelModel(VelModel*);
int VelModToInterface(VelModel* vel_model, GridDesc* grid, FILE* fp_fmmInterfacesFile);
int ModelToGrid3d(VelModel*, GridDesc*);

/*** program to generate  3-D vel/slowness grid */


#define NARGS 2
#define LENDELTA 6

int main(int argc, char *argv[]) {

    int istat;
    char fn_fmmInterfacesFile[MAXLINE];
    FILE* fp_fmmInterfacesFile = NULL;

    GridDesc mod_grid; /* model grid */
    VelModel vel_model; /* Velocity model */

    char fileRoot[MAXLINE];


    vel_model.buffer = NULL;
    vel_model.array = NULL;
    vel_model.deltax = vel_model.deltay = vel_model.deltaz = NULL;


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
    SetConstants();


    /* read control file */

    strcpy(fn_control, argv[1]);
    if ((fp_control = fopen(fn_control, "r")) == NULL) {
        nll_puterr("ERROR: opening control file.");
        exit(EXIT_ERROR_FILEIO);
    }

    if ((istat = ReadInterface2FmmInput(fp_control, &vel_model)) < 0) {
        exit(EXIT_ERROR_FILEIO);
    }
    mod_grid = grid_in;


    /* determine model coordinates mode - rect or latlon */
    SetModelCoordsMode(num_surfaces);


    /* read Interface Model input file */

    if ((istat = ReadVelModel(&vel_model, &mod_grid)) < 0) {
        if (istat == -1) nll_puterr("ERROR: reading first line of Velocity Model input file.");
        if (istat == -2) nll_puterr("ERROR: allocating memory for Velocity Model input file.");
        if (istat == -3) nll_puterr("ERROR: reading pos. of x-nodes of Velocity Model input file.");
        if (istat == -4) nll_puterr("ERROR: reading pos. of y-nodes of Velocity Model input file.");
        if (istat == -5) nll_puterr("ERROR: reading pos. of z-nodes of Velocity Model input file.");
        if (istat == -6) nll_puterr("ERROR: reading velocities of Velocity Model input file.");
        /*	nll_puterr("ERROR: reading Velocity Model input file."); */
        exit(EXIT_ERROR_FILEIO);
    }

    /* initialize 3D grid */

    /* allocate model grid */
    mod_grid.buffer = AllocateGrid(&mod_grid);
    if (mod_grid.buffer == NULL) {
        nll_puterr("ERROR: allocating memory for 3D slowness grid buffer.");
        exit(EXIT_ERROR_MEMORY);
    }
    /* create array access pointers */
    mod_grid.array = CreateGridArray(&mod_grid);
    if (mod_grid.array == NULL) {
        nll_puterr("ERROR: creating array for accessing 3D vel/slowness grid buffer.");
        exit(EXIT_ERROR_MEMORY);
    }

    sprintf(fileRoot, "interface_depth.mod");
    sprintf(MsgStr, "Creating interface depth grid files: %s.%s.*",
            fn_if_output, fileRoot);
    nll_putmsg(1, MsgStr);

    /* load vel model to grid */

    if ((istat = ModelToGrid3d(&vel_model, &mod_grid)) < 0) {
        nll_puterr("ERROR: loading velocity model to grid.");
        exit(EXIT_ERROR_MODEL);
    }

    /* save grid to disk */

    if ((istat = WriteGrid3dBuf(&mod_grid, NULL, fn_if_output, fileRoot)) < 0) {
        nll_puterr("ERROR: writing slowness grid to disk.");
        exit(EXIT_ERROR_IO);
    }


    // open ANU-FMM interfaces.in in current working directory
    sprintf(fn_fmmInterfacesFile, "%s_interfaces.in", fn_if_output);
    if ((fp_fmmInterfacesFile = fopen(fn_fmmInterfacesFile, "w")) == NULL) {
        nll_puterr("ERROR: opening ANU-FMM interfaces.in file.");
        return (-1);
    }

    /* write ANU-FMM interface file interfaces.in */

    if ((istat = VelModToInterface(&vel_model, &mod_grid, fp_fmmInterfacesFile)) < 0) {
        nll_puterr("ERROR: loading velocity model to grid.");
        exit(EXIT_ERROR_MODEL);
    }


    fclose(fp_fmmInterfacesFile);

    FreeVelModel(&vel_model);
    exit(EXIT_NORMAL);

}


#define VEL_FORMAT "%f"

/*** function to read Velocity Model input file */

int ReadVelModel(VelModel* vel_model, GridDesc* mod_grid) {

    int istat;
    int i, j, k;
    double dunit;
    float unit, tempfloat;
    char line[2 * MAXLINE];
    FILE* inpfile;


    if ((inpfile = fopen(fn_vg_input, "r")) == NULL) {
        nll_puterr("ERROR: opening Velocity Model input file.");
        return (-1);
    }

    if (fgets(line, 2 * MAXLINE, inpfile)) {

        // "%4f%3d%3d%3d%3d"
        istat = ReadFortranReal(line, 1, 4, &dunit);
        unit = (float) dunit;
        istat += ReadFortranInt(line, 5, 3, &vel_model->numx);
        istat += ReadFortranInt(line, 8, 3, &vel_model->numy);
        vel_model->numz = 1;
        fprintf(stderr, "\nReadVelModel: no. grid nodes in x,y,z -> %d %d %d \n", vel_model->numx, vel_model->numy, vel_model->numz);
        if (istat < 3) {
            fclose(inpfile);
            return (-1);
        }
        if (AllocateVelModel(vel_model) < 0) {
            fclose(inpfile);
            return (-2);
        }
        //if ((fgets(line, 2 * MAXLINE, inpfile)) && (strlen(line) >= (LENDELTA * vel_model->numx))) {
        fprintf(stderr, "ReadVelModel: grid nodes in x-dir\n\t");
        int nread = 0;
        for (i = 0; i < vel_model->numx; i++) {
            fscanf(inpfile, "%f", &tempfloat);
            /*  vel_model->deltax[i]=vel_model->origx+unit*tempfloat; */
            // AJL vel_model->deltax[i] = -unit*tempfloat;
            vel_model->deltax[i] = unit*tempfloat;
            fprintf(stderr, " %5.1f ", vel_model->deltax[i]);
            nread++;
        }
        fprintf(stderr, "\n");
        //} else {
        if (nread < vel_model->numx) {
            fclose(inpfile);
            return (-3);
        }

        //if ((fgets(line, 2 * MAXLINE, inpfile)) && (strlen(line) >= (LENDELTA * vel_model->numy))) {
        fprintf(stderr, "ReadVelModel: grid nodes in y-dir\n\t");
        nread = 0;
        for (i = 0; i < vel_model->numy; i++) {
            fscanf(inpfile, "%f", &tempfloat);
            /* vel_model->deltay[i]=vel_model->origy+unit*tempfloat; */
            fprintf(stderr, " %5.1f ", tempfloat);
            vel_model->deltay[i] = unit*tempfloat;
            /* fprintf(stderr,"%f\n",vel_model->deltay[i]);  */
            nread++;
        }
        fprintf(stderr, "\n");
        //} else {
        if (nread < vel_model->numy) {
            fclose(inpfile);
            return (-4);
        }

        fscanf(inpfile, "%*f %*f %*f");

        for (k = 0; k < vel_model->numz; k++) {
            /*  	for(j=vel_model->numy;j>0;j--) */
            for (j = 0; j < vel_model->numy; j++) {
                //if ((fgets(line, 2 * MAXLINE, inpfile)) && (strlen(line) >= (LENVEL * vel_model->numx))) {
                nread = 0;
                for (i = 0; i < vel_model->numx; i++) {
                    fscanf(inpfile, VEL_FORMAT, &vel_model->array[i][j][k]);
                    nread++;
                }
                if (nread < vel_model->numx) {
                    fclose(inpfile);
                    return (-6);
                }
            }
        }
    } else {
        fclose(inpfile);
        return (-1);
    }

    fclose(inpfile);
    return (0);
}

/*** function to read output file name ***/

int get_if_outfile(char* line1) {

    sscanf(line1, "%s", fn_if_output);

    sprintf(MsgStr, "interface2fmm files:  Output: %s.*", fn_if_output);
    nll_putmsg(3, MsgStr);

    return (0);
}

/*** returns the velocity in a definite point of the space ***/

float get2Dvalue(VelModel* vel_model, float posx, float posy) {
    float vel;
    int i, j;
    float t, u;

    //printf("get2Dvalue(vel_model, xval %f, yval %f\n", posx, posy);

    // AJL for (i = 0; (i < vel_model->numx) && (vel_model->deltax[i] > posx); i++);
    for (i = 0; (i < vel_model->numx) && (vel_model->deltax[i] < posx); i++);
    for (j = 0; (j < vel_model->numy) && (vel_model->deltay[j] < posy); j++);

    if (i == 0) {
        t = 0;
        i++;
    } else if (i == vel_model->numx) {
        t = 1;
        i--;
    } else
        t = (posx - vel_model->deltax[i - 1]) / (vel_model->deltax[i] - vel_model->deltax[i - 1]);

    if (j == 0) {
        u = 0;
        j++;
    } else if (j == vel_model->numy) {
        u = 1;
        j--;
    } else
        u = (posy - vel_model->deltay[j - 1]) / (vel_model->deltay[j] - vel_model->deltay[j - 1]);

    //printf("get2Dvalue(vel_model, i %d, j %d, t %f, u %f\n", i, j, t, u);

    vel = vel_model->array[i - 1][j - 1][0]*(1 - t)*(1 - u) +
            vel_model->array[i][j - 1][0]*(t)*(1 - u) +
            vel_model->array[i][j][0]*(t)*(u) +
            vel_model->array[i - 1][j][0]*(1 - t)*(u);

    return (vel);
}

/*** function to read input file name ***/

int get_if_inpfile(VelModel* vel_model, char* line1) {
    char type[MAXLINE];

    sscanf(line1, "%s %s", fn_vg_input, type);

    if (sscanf(line1, "%*s %*s %lf %lf %lf", &vel_model->origx, &vel_model->origy, &vel_model->origz) != 3) {
        nll_puterr("ERROR: reading SIMUL2K 3D Velocity input parameters");
        return (-1);
    }
    vel_model->type = TYPESIMUL;

    sprintf(MsgStr, "interface2fmm files:  Input: %s   Type: %s", fn_vg_input, type);
    nll_putmsg(3, MsgStr);

    return (0);
}

/*** function to read input file */

int ReadInterface2FmmInput(FILE* fp_input, VelModel * vel_model) {
    int istat, iscan;
    char param[MAXLINE], *pchr;
    char line[2 * MAXLINE], *fgets_return;

    int flag_control = 0, flag_inpfile = 0, flag_outfile = 0, flag_grid = 0,
            flag_trans = 0;
    int flag_include = 1;

    /*	vel_model->origx = vel_model->origy = vel_model->origz = 0; */

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

        if (strcmp(param, "IFINP") == 0) {
            if ((istat = get_if_inpfile(vel_model, strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading interface2fmm input parameters.");
            else
                flag_inpfile = 1;
        }


        /* read output file name (OUTfile) */

        if (strcmp(param, "IFOUT") == 0) {
            if ((istat = get_if_outfile(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading interface2fmm output file name.");
            else
                flag_outfile = 1;
        }


        /* read grid params */

        if (strcmp(param, "IFGRID") == 0) {
            if ((istat = get_grid(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading grid parameters.");
            else
                flag_grid = 1;
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
        nll_puterr("ERROR: no inputfile (IFINP) params read.");
    if (!flag_outfile)
        nll_puterr("ERROR: no outputfile (IFOUT) params read.");
    if (!flag_grid)
        nll_puterr("ERROR: no grid (IFGRID) params read.");

    if (!flag_trans) {
        sprintf(MsgStr, "INFO: no transformation (TRANS) params read.");
        nll_putmsg(1, MsgStr);
        Hypocenter.comment[0] = '\0';
    }

    return (flag_include * flag_control * flag_inpfile * flag_outfile * flag_grid - 1);
}

/*** function to load model to model grid ***/

int ModelToGrid3d(VelModel* vel_model, GridDesc* grid) {

    int ix, iy;
    //int imodel;
    double xval, yval, xloc, yloc;
    double value;


    /* generate grid values */
    /* Note:  staggered grid assumed, thus value lookup is shifted +dx/2, etc. */

    xval = grid->origx;
    for (ix = 0; ix < grid->numx; ix++) {
        yval = grid->origy;
        for (iy = 0; iy < grid->numy; iy++) {
            if (ModelCoordsMode == COORDS_LATLON) {
                rect2latlon(0, xval, yval, &yloc, &xloc);
            } else {
                xloc = xval;
                yloc = yval;
            }
            value = get2Dvalue(vel_model, xval, yval);
            //printf("DEBUG: value %f = get2Dvalue(vel_model, xval %f, yval %f)\n", value, xval, yval);

            if (value < 0.0) {
                //nll_puterr("ERROR: cannot get grid value.");
                //return (-1);
            }

            switch (grid->type) {

                case GRID_DEPTH:
                    ((GRID_FLOAT_TYPE ***) grid->array)[ix][iy][0] = value;
                    break;

                default:
                    nll_puterr("ERROR: unrecognized grid type.");
                    return (-1);

            }

            yval += grid->dy;
        }
        xval += grid->dx;
    }


    return (0);

}

int AllocateVelModel(VelModel * vel_model) {
    int ix, iy, numyz;
    float ***garray;


    /* Allocates memory for 3D Velocity Model node distance arrays */

    vel_model->deltax = (float *) malloc((size_t) (vel_model->numx * sizeof (float)));
    vel_model->deltay = (float *) malloc((size_t) (vel_model->numy * sizeof (float)));
    vel_model->deltaz = (float *) malloc((size_t) (vel_model->numz * sizeof (float)));

    if ((vel_model->deltax == NULL) || (vel_model->deltay == NULL) || (vel_model->deltaz == NULL)) {
        nll_puterr("ERROR: allocating memory for 3D Velocity Model node distance arrays.");
        return (-1);
    }
    NumAllocations += 3;



    /* allocate Velocity model buffer */

    vel_model->buffer = (float *) malloc((size_t)
            (vel_model->numx * vel_model->numy * vel_model->numz * sizeof (float)));
    if (vel_model->buffer == NULL) {
        nll_puterr("ERROR: allocating memory for 3D Velocity Model buffer.");
        return (-1);
    }
    NumAllocations++;



    /* creates array for accessing 3D velocity model */

    garray = (float ***) malloc((size_t) vel_model->numx * sizeof (float **));
    if (garray == NULL) {
        nll_puterr("ERROR: creating array for accessing 3D Velocity Model buffer.");
        return (-1);
    }
    NumAllocations++;
    numyz = vel_model->numy * vel_model->numz;
    for (ix = 0; ix < vel_model->numx; ix++) {
        garray[ix] = (float **) malloc((size_t) vel_model->numy * sizeof (float *));
        if (garray[ix] == NULL) {
            nll_puterr("ERROR: creating array for accessing 3D Velocity Model buffer.");
            return (-1);
        }
        NumAllocations++;
        for (iy = 0; iy < vel_model->numy; iy++)
            garray[ix][iy] = vel_model->buffer + ix * numyz + iy * vel_model->numz;
    }
    vel_model->array = garray;
    return (0);
}

/*** function to free the allocations in the velocity model structure ***/
void FreeVelModel(VelModel * vel_model) {
    int ix;

    /*** frees the array for accessing 3D velocity model ***/
    if (vel_model->array) {
        for (ix = 0; ix < vel_model->numx; ix++) {
            free(vel_model->array[ix]);
            NumAllocations--;
        }

        free(vel_model->array);
        NumAllocations--;

        vel_model->array = NULL;
    }


    /*** frees the buffer for 3D velocity model ***/
    if (vel_model->buffer) {
        free(vel_model->buffer);
        NumAllocations--;
        vel_model->buffer = NULL;
    }


    /* frees the arrays for 3D Velocity Model node distance */
    if (vel_model->deltax)
        free(vel_model->deltax);
    if (vel_model->deltay)
        free(vel_model->deltay);
    if (vel_model->deltaz)
        free(vel_model->deltaz);

}

/*** function to write ANU-FMM interface file interfaces.in */

int VelModToInterface(VelModel* vel_model, GridDesc* grid, FILE* fp_fmmInterfacesFile) {


    int pad_for_propagation_grid = 2; // propagation grid is extended at each boudary:
    int pad_for_cublic_spline = 2 + pad_for_propagation_grid; // velocity grid is extended 2 cell beyond pad_for_propagation_grid at each boudary:
    // cubic spline interpolation requires that all the nodes of the FM propagation grid must be
    // INSIDE the SECOND layer of velocity grid nodes as counted from the grid boundaries

    // estimate lat and lon step from lat/lon of upper right and lower left corners of grid
    double xoff = 0.0;
    double yoff = 0.0;
    double lat_ll, lon_ll;
    rect2latlon(0, grid->origx + grid->dx * xoff, grid->origy + grid->dy * yoff, &lat_ll, &lon_ll);
    xoff = (double) (grid->numx - 1);
    yoff = (double) (grid->numy - 1);
    double lat_ur, lon_ur;
    rect2latlon(0, grid->origx + grid->dx * xoff, grid->origy + grid->dy * yoff, &lat_ur, &lon_ur);
    double lat_step = (lat_ur - lat_ll) / (double) (grid->numy - 2);
    double lat_ul, lon_ul;
    rect2latlon(0, grid->origx, grid->origy + grid->dy * yoff, &lat_ul, &lon_ul);
    double lat_lr, lon_lr;
    rect2latlon(0, grid->origx + grid->dx * xoff, grid->origy, &lat_lr, &lon_lr);

    // get lat/lon of grid center
    xoff = (double) (grid->numx / 2);
    yoff = (double) (grid->numy / 2);
    double lat_cent, lon_cent;
    rect2latlon(0, grid->origx + grid->dx * xoff, grid->origy + grid->dy * yoff, &lat_cent, &lon_cent);

    // to guarantee that FMM grid conaints velocity model grid,
    // estimate long step from uppoer (lower) boundary for Northern (Southern) hemisphere
    double lon_step = 0.0;
    if (lat_cent > 0.0) { // Northern hemisphere
        lon_step = (lon_ur - lon_ul) / (double) (grid->numx - 2);
    } else {
        lon_step = (lon_lr - lon_ll) / (double) (grid->numx - 2);
    }

    // get origin for FMM
    double lat_orig = lat_cent - yoff * lat_step; // orig FMM is calculated from grid center, not NLL grid origin.
    double lat_orig_vel = lat_orig - lat_step * pad_for_cublic_spline;
    double lon_orig = lon_cent - xoff * lon_step; // orig FMM is calculated from grid center, not NLL grid origin.
    double lon_orig_vel = lon_orig - lon_step * pad_for_cublic_spline;
    // next    : origin of velocity grid in radius (km) and in lat and long (in radians) for grid 1
    double z_orig = grid->origz + grid->dz * (double) (grid->numz - 1);
    double z_orig_vel = z_orig + grid->dz * pad_for_cublic_spline;

    // write interfaces.in header
    int n_interfaces_to_write = 3;
    fprintf(fp_fmmInterfacesFile, "%d    # number of interfaces\n", n_interfaces_to_write);
    fprintf(fp_fmmInterfacesFile, "%d %d      # number of interface grid nodes in lat and in long\n",
            grid->numy + 2 * pad_for_cublic_spline, grid->numx + 2 * pad_for_cublic_spline);
    fprintf(fp_fmmInterfacesFile, "%lf %lf      # grid node spacing in lat and in long (in radians)  = %fdeg %fdeg\n",
            lat_step * cRPD, lon_step * cRPD, lat_step, lon_step);
    fprintf(fp_fmmInterfacesFile, "%lf %lf      # origin of interface grid in lat and long (in radians)  = %fdeg %fdeg\n",
            lat_orig_vel * cRPD, lon_orig_vel *cRPD, lat_orig_vel, lon_orig_vel);


    // next    : values of the velocity at each velocity grid node for grid 1, in km/sec,
    //           one per line with longitude index varying the fastest, latitude
    //           index the second fastest and radius index the slowest


    // generate grid values
    double xval, yval;
    int ix, iy, iz;
    double interface_depth;
    double radius_top = AVG_ERAD - (grid->origz - grid->dz * pad_for_propagation_grid);
    int n_interfaces_written = 0;
    double zdepth = z_orig_vel;
    int iz_start = -1;
    int iz_end = 1;
    for (iz = iz_start; iz <= iz_end; iz++) {
        double ylat = lat_orig_vel;
        for (iy = -pad_for_cublic_spline; iy < grid->numy + pad_for_cublic_spline; iy++) {
            double xlon = lon_orig_vel;
            for (ix = -pad_for_cublic_spline; ix < grid->numx + pad_for_cublic_spline; ix++) {
                // vgrids.in
                if (ModelCoordsMode == COORDS_LATLON) {
                    xval = xlon;
                    yval = ylat;
                } else {
                    latlon2rect(0, ylat, xlon, &xval, &yval);
                }
                interface_depth = get2Dvalue(vel_model, xval, yval);
                if (interface_depth < 0.0) {
                    interface_depth = radius_top;
                }
                // interface.in
                if (iz == -1)
                    fprintf(fp_fmmInterfacesFile, "%lf\n", radius_top); // top of prop grid
                else if (iz == 0)
                    fprintf(fp_fmmInterfacesFile, "%lf\n", AVG_ERAD - interface_depth); // interface
                else if (iz == 1)
                    fprintf(fp_fmmInterfacesFile, "%lf\n", AVG_ERAD - (grid->origz + grid->dz * (grid->numz - 1 + pad_for_propagation_grid))); // bottom of prop grid

                xlon += lon_step;
            }
            ylat += lat_step;
        }
        n_interfaces_written++;
        zdepth -= grid->dz;
    }


    return (0);

}

