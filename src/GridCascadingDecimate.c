/*
 * Copyright (C) 2016 Anthony Lomax <anthony@alomax.net>
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


/*   GridCascadingDecimate.c

        Program to convert a regular 3D grid file to a cascading decimated file

 */


/*
        history:

        ver 01    20161019  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

 */



#include "GridLib.h"


// defines


// globals


// functions

int DoCascadingDecimateProcess(int argc, char *argv[], double depths[], int ndepths);



/*** Program to process (add, ) 3D grid files */

#define PNAME  "GridCascadingDecimate"

int main(int argc, char *argv[]) {

    int narg;
    char doubling_depths[1024];
    double depths[MAX_NUM_Z_MERGE_DEPTHS];


    // set program name

    strcpy(prog_name, PNAME);


    // check command line for correct usage

    fprintf(stdout, "\n%s Arguments: ", prog_name);
    for (narg = 0; narg < argc; narg++)
        fprintf(stdout, "<%s> ", argv[narg]);
    fprintf(stdout, "\n");

    int ndepths = 0;
    if (argc > 3) {
        sscanf(argv[1], "%s", doubling_depths);
        char *str_pos = strtok(doubling_depths, ",");
        while (str_pos != NULL) {
            depths[ndepths] = atof(str_pos);
            printf("INFO: depths added: %s %f\n", str_pos, depths[ndepths]);
            ndepths++;
            str_pos = strtok(NULL, ",");
        }
    } else {
        disp_usage(PNAME,
                "<doubling_depths> <input grid(s)> <output grid path> [<flag_use_mean_value_in_casc_cell>]\n"
                "   doubling_depths - comma separated, increasing (approx) depths to double cell size"
                );
        exit(-1);
    }

    DoCascadingDecimateProcess(argc, argv, depths, ndepths);

    exit(0);

}

int DoCascadingDecimateProcess(int argc, char *argv[], double depths[], int ndepths) {
    int istat;


#define MAX_NUM_INPUT_FILES 4096
    char fn_grid_in_list[MAX_NUM_INPUT_FILES][FILENAME_MAX];
    char fn_grid_in_base[FILENAME_MAX];
    strcpy(fn_grid_in_base, argv[2]);

    // check for wildcards in input file name
    int numFiles;
    ;
    if ((numFiles = ExpandWildCards(fn_grid_in_base, fn_grid_in_list, MAX_NUM_INPUT_FILES)) < 1) {
        nll_puterr2("ERROR: no matching grid files found: ", fn_grid_in_base);
        return (-1);
    }
    if (numFiles >= MAX_NUM_INPUT_FILES) {
        sprintf(MsgStr, "WARNING: maximum number of grid files exceeded, only first %d will be processed.", MAX_NUM_INPUT_FILES);
        nll_puterr(MsgStr);
    }


    // output file root
    char grid_out_path[FILENAME_MAX];
    strcpy(grid_out_path, argv[3]);

    // flag_use_mean_value_in_casc_cell
    int flag_use_mean_value_in_casc_cell = 0;
    if (argc > 4 && atoi(argv[4]) == 1) {
        flag_use_mean_value_in_casc_cell = 1;
    }

    char fn_grid_in[FILENAME_MAX];
    char fn_grid_out[FILENAME_MAX];
    FILE *fp_grid_in;
    FILE *fp_grid_in_hdr;
    GridDesc grid_in, grid_out;
    SourceDesc sourceDesc;
    char file_type[FILENAME_MAX];


    for (int nFile = 0; nFile < numFiles; nFile++) {

        // input file name
        strcpy(fn_grid_in, fn_grid_in_list[nFile]);
        if (strstr(fn_grid_in, ".buf") != NULL) {
            *strrchr(fn_grid_in, '.') = '\0'; // remove extension from output filename
        } else {
            continue; // not grid buffer file
        }
        //if (strrchr(fn_grid_in, '.') != NULL) {
        //    *strrchr(fn_grid_in, '.') = '\0'; // remove type from output filename
        //}

        // output file name
        strcpy(fn_grid_out, grid_out_path);
        // find input file name (strip path)
        char *cpos;
        if ((cpos = strrchr(fn_grid_in, '/')) != NULL && strlen(fn_grid_in) > cpos - fn_grid_in + 1) {
            ;
        } else {
            cpos = 0;
        }
        strcat(fn_grid_out, "/casc.");
        strcat(fn_grid_out, cpos + 1);

        // set grid file type (e.g. time, ...)
        if (strrchr(fn_grid_in, '.') != NULL) {
            strcpy(file_type, strrchr(fn_grid_in, '.') + 1);
            *strrchr(fn_grid_out, '.') = '\0'; // remove type from output filename
        } else {
            strcpy(file_type, "");
        }

        sprintf(MsgStr, "Processing grid: %s -> %s", fn_grid_in, fn_grid_out);
        nll_putmsg(0, MsgStr);

        // open input grid file
        if ((istat = OpenGrid3dFile(fn_grid_in, &fp_grid_in, &fp_grid_in_hdr,
                &grid_in, file_type, &sourceDesc, 0)) < 0) {
            nll_puterr("ERROR opening input grid file.");
            return (-1);
        }
        // set map projection if available, so will be written to output grid header
        if (strlen(grid_in.mapProjStr) > 0) {
            strcpy(MapProjStr[0], grid_in.mapProjStr);
        }




        // find file root phase.type
        char fileRoot[MAXLINE];
        if (strchr(fn_grid_in, '.') != NULL) {
            strcpy(fileRoot, strchr(fn_grid_in, '.'));
        }

        // create and initialize output grid

        // create output grid description
        grid_out = grid_in;
        grid_in.iSwapBytes = 0;
        grid_out.iSwapBytes = 0;
        setCascadingGrid(&grid_out);
        for (int n = 0; n < ndepths; n++) {
            grid_out.gridDesc_Cascading.z_merge_depths[n] = depths[n];
        }
        grid_out.gridDesc_Cascading.num_z_merge_depths = ndepths;

        // allocate grid
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

        // determine if should use mean slownesses or simple mean
        int iuse_inverse_mean = grid_out.type == GRID_VELOCITY || grid_out.type == GRID_VELOCITY_METERS;
        // create array for storing count of slownesses used for mean
        int nx = 0, ny = 0, nz = 0;
        if (flag_use_mean_value_in_casc_cell) {
            nx = grid_out.numx;
            ny = grid_out.numy;
            nz = grid_out.gridDesc_Cascading.zindex[grid_out.numz - 1];
            printf("INFO: flag_use_mean_value_in_casc_cell: %d\n", flag_use_mean_value_in_casc_cell);
        } else {
            printf("INFO: NO flag_use_mean_value_in_casc_cell: %d\n", flag_use_mean_value_in_casc_cell);

        }
        int (*imean_count)[ny][nz] = malloc(sizeof (int[nx][ny][nz]));
        int ix, iy, iz;
        for (ix = 0; ix < nx; ix++) {
            for (iy = 0; iy < ny; iy++) {
                for (iz = 0; iz < nz; iz++) {
                    imean_count[ix][iy][iz] = 0;
                }
            }
        }

        // map input grid file into output cascaded grid

        int narray = 0;
        int nz_casc;
        int nx_casc;
        int ny_casc;
        int nz_casc_last = -1;
        float val;
        double val_mean;
        GRID_FLOAT_TYPE *pval_casc;
        int icount;
        for (iz = 0; iz < grid_in.numz; iz++) {
            // skip to regular grid z index that corresponds to next cascading grid node
            nz_casc = grid_out.gridDesc_Cascading.zindex[iz];
            if (!flag_use_mean_value_in_casc_cell && nz_casc == nz_casc_last) {
                continue; // in same cascading z cell
            }
            nz_casc_last = nz_casc;
            int xyz_scale = grid_out.gridDesc_Cascading.xyz_scale[iz];
            int nx_casc_last = -1;
            for (ix = 0; ix < grid_in.numx; ix++) {
                nx_casc = ix / xyz_scale;
                // add fractional final casc grid node if not exact alignment with final reg grid node
                if (ix == grid_in.numx - 1 && ix % xyz_scale != 0) {
                    nx_casc++;
                }
                if (!flag_use_mean_value_in_casc_cell && nx_casc == nx_casc_last) {
                    continue; // in same cascading x cell
                }
                nx_casc_last = nx_casc;
                //printf("ix = %d/%d\r", ix, grid_in.numx,);
                int ny_casc_last = -1;
                for (iy = 0; iy < grid_in.numy; iy++) {
                    ny_casc = iy / xyz_scale;
                    // add fractional final casc grid node if not exact alignment with final reg grid node
                    if (iy == grid_in.numy - 1 && iy % xyz_scale != 0) {
                        ny_casc++;
                    }
                    if (!flag_use_mean_value_in_casc_cell && ny_casc == ny_casc_last) {
                        continue; // in same cascading y cell
                    }
                    ny_casc_last = ny_casc;
                    val = ReadGrid3dValue(fp_grid_in, ix, iy, iz, &grid_in, 0);
                    //if ((ix == 0 || ix == grid_in.numx - 1) && (iy == 0 || iy == grid_in.numy - 1)) {
                    //    printf("ixyz %d %d %d, nx_casc %d, ny_casc %d, nz_casc %d, val %f, ptr %ld, size %d\n", ix, iy, iz, nx_casc, ny_casc, nz_casc, val,
                    //            (((GRID_FLOAT_TYPE ***) grid_out.array)[nx_casc][ny_casc] + nz_casc - ((GRID_FLOAT_TYPE ***) grid_out.array)[0][0]),
                    //            sizeof (GRID_FLOAT_TYPE));
                    //}
                    pval_casc = ((GRID_FLOAT_TYPE ***) grid_out.array)[nx_casc][ny_casc] + nz_casc;
                    // update mean slowness in this casc grid cell
                    if (flag_use_mean_value_in_casc_cell) {
                        if ((icount = imean_count[nx_casc][ny_casc][nz_casc]) > 0) {
                            if (iuse_inverse_mean) {
                                val_mean = 1.0 / *pval_casc;
                                val = 1.0 / val;
                            } else {
                                val_mean = *pval_casc;
                            }
                            val_mean = ((val_mean * (double) icount) + val) / (double) (icount + 1);
                            if (iuse_inverse_mean) {
                                val = 1.0 / val_mean;
                            } else {
                                val = val_mean;
                            }
                        }
                        imean_count[nx_casc][ny_casc][nz_casc]++;
                    }
                    // set new value for cascading grid cell
                    *pval_casc = val;
                    narray++;
                }
            }
        }
        //printf("\n");


        // save sum grid to disk

        if ((istat = WriteGrid3dBuf(&grid_out, &sourceDesc, fn_grid_out, file_type)) < 0) {
            nll_puterr("ERROR: writing output cascading grid to disk.\n");
            return (-1);
        }


        fclose(fp_grid_in);
        fclose(fp_grid_in_hdr);

    }

    return (0);

}



//------------------------------------------------------------/
// Anthony Lomax           | email: anthony@alomax.net        /
// UMR Geosciences Azur    | web: www.alomax.net              /
// 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        /
// 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        /
//------------------------------------------------------------/

