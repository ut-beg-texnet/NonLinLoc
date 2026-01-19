/*
 * Copyright (C) 2025 Anthony Lomax <anthony@alomax.net>
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


/*   GridCut.c

        Program to extract a sub-region from a 3D grid and save to disk

 */


/*
        history:

        ver 01    20250903  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

 */



#include "GridLib.h"


// defines


// globals


// functions

int DoGridCutProcess(int argc, char *argv[]);



/*** Program to process 3D grid files */

#define PNAME  "GridCut"

int main(int argc, char *argv[]) {

    int narg;


    // set program name

    strcpy(prog_name, PNAME);


    // check command line for correct usage

    fprintf(stdout, "\n%s Arguments: ", prog_name);
    for (narg = 0; narg < argc; narg++)
        fprintf(stdout, "<%s> ", argv[narg]);
    fprintf(stdout, "\n");

    disp_usage(PNAME,
            "<input grid path/name> <output grid root> \"<output grid description>\" [\"source description\"]\n"
            "   output grid description - quote delimited string, with:\n"
            "       num_grid_x/y/z orig_grid_x/y/z d_grid_x/y/z\n"
            "       e.g. \"129 105 87  -32.0 -28.0 -3.0  0.5 0.5 0.5\"\n"
            "   source description - optional quote delimited NLL source string.\n"
            "       will supersede orig_grid_x/y in grid description, shifting grid center to source lat/lon.\n"
            "       e.g. \"GTSRCE TX_MB21 LATLON 32.343507 -101.727966 0 0.7880\"\n"
            );

    if (argc < 4) {
        nll_puterr("ERROR: wrong number of command line arguments.");
        exit(-1);
    }

    // set constants
    SetConstants();
    prog_mode_3d = 1;
    NumSources = 0;

    DoGridCutProcess(argc, argv);

    exit(0);

}

int DoGridCutProcess(int argc, char *argv[]) {

    int istat;

    // input file name
    char fn_grid_in[FILENAME_MAX];
    strcpy(fn_grid_in, argv[1]);

    // output file root
    char fn_grid_out[FILENAME_MAX];
    strcpy(fn_grid_out, argv[2]);

    // output file grid description
    char grid_out_desc[MAXLINE];
    strcpy(grid_out_desc, argv[3]);

    int flag_source = 0;
    if (argc > 4) {
        char grid_out_source[MAXLINE];
        strcpy(grid_out_source, argv[4]);
        if ((istat = GetNextSource(strchr(grid_out_source, ' '))) < 0) {
            nll_puterr("ERROR: reading source description:");
            exit(-1);
        } else
            flag_source = 1;
    }


    FILE *fp_grid_in;
    FILE *fp_grid_in_hdr;
    SourceDesc sourceDesc;
    char file_type[FILENAME_MAX];
    // clean up grid filenames
    if (strstr(fn_grid_in, ".buf") != NULL) {
        *strrchr(fn_grid_in, '.') = '\0'; // remove extension from input filename
    } else if (strstr(fn_grid_in, ".hdr") != NULL) {
        *strrchr(fn_grid_in, '.') = '\0'; // remove extension from input filename
    }
    if (strstr(fn_grid_out, ".buf") != NULL) {
        *strrchr(fn_grid_out, '.') = '\0'; // remove extension from output filename
    } else if (strstr(fn_grid_out, ".hdr") != NULL) {
        *strrchr(fn_grid_out, '.') = '\0'; // remove extension from output filename
    }

    // set grid file type (e.g. time, ...)
    if (strrchr(fn_grid_in, '.') != NULL) {
        strcpy(file_type, strrchr(fn_grid_in, '.') + 1);
        *strrchr(fn_grid_out, '.') = '\0'; // remove type from output filename
    } else {
        strcpy(file_type, "");
    }

    // if source specified, append source name
    /*if (flag_source) {
        strcat(fn_grid_out, "_");
        strcat(fn_grid_out, Source[0].label);
    }*/

    printf("Processing grid: %s -> %s\n", fn_grid_in, fn_grid_out);

    // open input grid file
    GridDesc grid_input;
    if ((istat = OpenGrid3dFile(fn_grid_in, &fp_grid_in, &fp_grid_in_hdr,
            &grid_input, file_type, &sourceDesc, 0)) < 0) {
        nll_puterr("ERROR opening input grid file.");
        return (-1);
    }
    // set map projection if available, so will be written to output grid header
    int have_projection = 0;
    if (strlen(grid_input.mapProjStr) > 0) {
        strcpy(MapProjStr[0], grid_input.mapProjStr);
        have_projection = 1;
    }
    // allocate input grid
    grid_input.buffer = AllocateGrid(&grid_input);
    if (grid_input.buffer == NULL) {
        nll_puterr("ERROR: allocating memory for input grid buffer.\n");
        return (-1);
    }
    // create grid array access pointers
    grid_input.array = CreateGridArray(&grid_input);
    if (grid_input.array == NULL) {
        nll_puterr("ERROR: creating array for accessing input grid buffer.\n");
        return (-1);
    }


    // read grid
    printf("Reading grid: %s\n", fn_grid_in);
    if ((istat = ReadGrid3dBuf(&grid_input, fp_grid_in)) < 0) {
        nll_puterr("ERROR: reading input grid grid from disk.");
        return (EXIT_ERROR_FILEIO);
    }
    CloseGrid3dFile(&grid_input, &fp_grid_in, &fp_grid_in_hdr);
    printf("Input grid:\n");
    display_grid_param(&grid_input);

    // create and initialize output grid

    // complete output grid description
    printf("DEBUG: grid_out_desc |%s|\n", grid_out_desc);
    strcat(grid_out_desc, "  ");
    strcat(grid_out_desc, grid_input.chr_type);
    printf("DEBUG: grid_out_desc |%s|\n", grid_out_desc);

    // parse grid description
    get_grid(grid_out_desc);
    GridDesc grid_out = grid_in; // grid_in initialized in GridLib.c
    printf("Output grid:\n");
    display_grid_param(&grid_out);

    // if source specified, shift grid
    if (flag_source) {
        double xrect, yrect;
        if (Source[0].is_coord_xyz) {
            xrect = Source[0].x;
            yrect = Source[0].y;
        } else {
            if (!have_projection) {
                nll_puterr("ERROR: Cannot process lat/lon source, no coordinate transform available.");
                return (-1);
            }
            // set transform
            char trans_str[2 * MAXLINE];
            projection_str2transform_str(trans_str, MapProjStr[0]);
            get_transform(0, trans_str);
            // transform lat/lon to x/y
            latlon2rect(0, Source[0].dlat, Source[0].dlong, &xrect, &yrect);
            printf("DEBUG: lat/lon -> x/y:  %f %f -> %f %f\n", Source[0].dlat, Source[0].dlong, xrect, yrect);
        }
        grid_out.origx = xrect - ((double) grid_out.numx * grid_out.dx / 2.0);
        grid_out.origy = yrect - ((double) grid_out.numy * grid_out.dy / 2.0);
        printf("DEBUG: origx/origy:  %f %f\n", grid_out.origx, grid_out.origy);
        printf("Shifted output grid center to lat/lon: %f %f\n", Source[0].dlat, Source[0].dlong);
        printf("Output grid:\n");
        display_grid_param(&grid_out);
    }

    // check for overlap between in and out grids
    int icount_out = 0;
    if (grid_out.origx > grid_input.origx + (double) grid_input.numx * grid_input.dx) {
        nll_puterr("ERROR: output grid xmin > input grid xmax.");
        icount_out++;
    }
    if (grid_out.origy > grid_input.origy + (double) grid_input.numy * grid_input.dy) {
        nll_puterr("ERROR: output grid ymin > input grid ymax.");
        icount_out++;
    }
    if (grid_out.origz > grid_input.origz + (double) grid_input.numz * grid_input.dz) {
        nll_puterr("ERROR: output grid zmin > input grid zmax.");
        icount_out++;
    }
    if (grid_out.origx + (double) grid_out.numx * grid_out.dx < grid_input.origx) {
        nll_puterr("ERROR: output grid xmax < input grid xmin.");
        icount_out++;
    }
    if (grid_out.origy + (double) grid_out.numy * grid_out.dy < grid_input.origy) {
        nll_puterr("ERROR: output grid ymax < input grid ymin.");
        icount_out++;
    }
    if (grid_out.origz + (double) grid_out.numz * grid_out.dz < grid_input.origz) {
        nll_puterr("ERROR: output grid zmax < input grid zmin.");
        icount_out++;
    }
    if (icount_out > 0) {
        nll_puterr("ERROR: output grid does not overlap input grid.\n");
        return (-1);
    }

    // check for overlap between in and out grids
    int icount_overlap = 0;
    if (grid_out.origx < grid_input.origx) {
        nll_puterr("WARNING: output grid xmin < input grid xmin, clipping output grid to input grid");
        grid_out.numx -= (int) ((grid_input.origx - grid_out.origx) / grid_out.dx);
        grid_out.origx = grid_input.origx;
        icount_overlap++;
    }
    if (grid_out.origy < grid_input.origy) {
        nll_puterr("WARNING: output grid ymin < input grid ymin, clipping output grid to input grid.");
        grid_out.numy -= (int) ((grid_input.origy - grid_out.origy) / grid_out.dy);
        grid_out.origy = grid_input.origy;
        icount_overlap++;
    }
    if (grid_out.origz < grid_input.origz) {
        nll_puterr("WARNING: output grid zmin < input grid zmin, clipping output grid to input grid.");
        grid_out.numz -= (int) ((grid_input.origz - grid_out.origz) / grid_out.dz);
        grid_out.origz = grid_input.origz;
        icount_overlap++;
    }
    double d_origx = grid_out.origx + (double) grid_out.numx * grid_out.dx - (grid_input.origx + (double) grid_input.numx * grid_input.dx);
    if (d_origx > 0.0) {
        nll_puterr("WARNING: output grid xmax > input grid xmax, clipping output grid to input grid.");
        grid_out.numx -= (int) (d_origx / grid_out.dx);
        icount_overlap++;
    }
    double d_origy = grid_out.origy + (double) grid_out.numy * grid_out.dy - (grid_input.origy + (double) grid_input.numy * grid_input.dy);
    if (d_origy > 0.0) {
        nll_puterr("WARNING: output grid ymax > input grid ymax, clipping output grid to input grid.");
        grid_out.numy -= (int) (d_origy / grid_out.dy);
        icount_overlap++;
    }
    double d_origz = grid_out.origz + (double) grid_out.numz * grid_out.dz - (grid_input.origz + (double) grid_input.numz * grid_input.dz);
    if (d_origz > 0.0) {
        nll_puterr("WARNING: output grid zmax > input grid zmax, clipping output grid to input grid.");
        grid_out.numz -= (int) (d_origz / grid_out.dz);
        icount_overlap++;
    }
    if (icount_overlap > 0) {
        printf("INFO: output grid clipped to input grid.\n");
        printf("Output grid:\n");
        display_grid_param(&grid_out);
    }

    // allocate output grid
    grid_out.buffer = AllocateGrid(&grid_out);
    if (grid_out.buffer == NULL) {
        nll_puterr("ERROR: allocating memory for output grid buffer.\n");
        return (-1);
    }
    // create grid array access pointers
    grid_out.array = CreateGridArray(&grid_out);
    if (grid_out.array == NULL) {
        nll_puterr("ERROR: creating array for accessing output grid buffer.\n");
        return (-1);
    }

    // map input grid file into output cut grid
    float val;
    int ix, iy, iz;
    double x_coord = grid_out.origx;
    printf("Processing:\n");
    for (ix = 0; ix < grid_out.numx; ix++) {
        printf("  x=%d/%d\r", ix, grid_out.numx);
        double y_coord = grid_out.origy;
        for (iy = 0; iy < grid_out.numy; iy++) {
            double z_coord = grid_out.origz;
            for (iz = 0; iz < grid_out.numz; iz++) {
                val = ReadAbsInterpGrid3d(fp_grid_in, &grid_input, x_coord, y_coord, z_coord, 0);
                ((GRID_FLOAT_TYPE***) grid_out.array)[ix][iy][iz] = val;
                z_coord += grid_out.dz;
            }
            y_coord += grid_out.dy;
        }
        x_coord += grid_out.dx;
    }
    printf("\n");


    // save processed grid to disk
    printf("Writing output cut grid to disk: %s.*\n", fn_grid_out);
    if ((istat = WriteGrid3dBuf(&grid_out, &sourceDesc, fn_grid_out, file_type)) < 0) {
        nll_puterr("ERROR: writing output cut grid to disk.\n");
        return (-1);
    }

    return (0);

}



