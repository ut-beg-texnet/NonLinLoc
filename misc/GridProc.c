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


/*   GridProc.c

        Program to process (add, ) 3D grid files

 */



/*
        history:

        ver 01    04Nov1997  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

 */



#include "../src/GridLib.h"


// defines

#define MAX_NUM_INPUT_FILES 500



// globals


// functions

int DoProcess(int argc, char *argv[], int process_type);



/*** Program to process (add, ) 3D grid files */

#define PNAME  "GridProc"

#define PROC_ADD 0
#define PROC_SUB 1
#define PROC_MUL_CONST 2

int main(int argc, char *argv[]) {

    int istat, narg, nchar;
    char process[MAXLINE];


    // set program name

    strcpy(prog_name, PNAME);

    char usage[] =
            "ADD/SUB/MUL_CONST,... <arguments>"
            "\n\tADD <size_gridfile> <output_gridfile> <add_gridfile_list> [<add_gridfile> ...]"
            "\n\tSUB <size_gridfile> <output_gridfile> <left_gridfile> <right_gridfile>"
            "\n\tMUL_CONST <size_gridfile> <output_gridfile> <in_gridfile> <constant>";


    // check command line for correct usage

    fprintf(stdout, "%s Arguments: ", prog_name);
    for (narg = 0; narg < argc; narg++)
        fprintf(stdout, "<%s> ", argv[narg]);
    fprintf(stdout, "\n");

    if (argc > 1)
        sscanf(argv[1], "%s", process);
    else
        strcpy(process, "");
    nchar = 0;
    while (nchar < MAXLINE && (process[nchar] = toupper(process[nchar])))
        nchar++;

    if (strcmp(process, "ADD") == 0) {
        if (argc < 4) {
            nll_puterr("ERROR wrong number of command line arguments.");
            disp_usage(PNAME, usage);
            exit(-1);
        }
        if ((istat = DoProcess(argc, argv, PROC_ADD)) < 0) {
            nll_puterr("ERROR doing ADD process.");
            disp_usage(PNAME, usage);
            exit(-1);
        }
    } else if (strcmp(process, "SUB") == 0) {
        if (argc < 5) {
            nll_puterr("ERROR wrong number of command line arguments.");
            disp_usage(PNAME, usage);
            exit(-1);
        }
        if ((istat = DoProcess(argc, argv, PROC_SUB)) < 0) {
            nll_puterr("ERROR doing SUB process.");
            disp_usage(PNAME, usage);
            exit(-1);
        }
    } else if (strcmp(process, "MUL_CONST") == 0) {
        if (argc < 5) {
            nll_puterr("ERROR wrong number of command line arguments.");
            disp_usage(PNAME, usage);
            exit(-1);
        }
        if ((istat = DoProcess(argc, argv, PROC_MUL_CONST)) < 0) {
            nll_puterr("ERROR doing MUL_CONST process.");
            disp_usage(PNAME, usage);
            exit(-1);
        }
    } else {
        sprintf(MsgStr,
                "ERROR unrecognized process - <%s>.", process);
        nll_puterr(MsgStr);
        disp_usage(PNAME, usage);
        exit(-1);
    }



    exit(0);

}

int DoProcess(int argc, char *argv[], int process_type) {
    int istat;

    char fn_grid_size[FILENAME_MAX];
    char fn_grid_out[FILENAME_MAX], fn_grid_in[FILENAME_MAX];
    FILE *fp_grid_size, *fp_grid_in;
    FILE *fp_grid_size_hdr, *fp_grid_in_hdr;

    GridDesc grid_size, grid_out, grid_in;

    int nFile, numFiles;
    char fn_grid_in_list[MAX_NUM_INPUT_FILES][FILENAME_MAX];
    char *chrpos, test_str[] = ".buf";

    SourceDesc srce;

    // open size grid file

    strcpy(fn_grid_size, argv[2]);
    if ((istat = OpenGrid3dFile(fn_grid_size, &fp_grid_size,
            &fp_grid_size_hdr,
            &grid_size, "", &srce, 0)) < 0) {
        nll_puterr("ERROR opening size grid file.");
        return (-1);
    }



    // create and initialize output grid

    // output file name
    strcpy(fn_grid_out, argv[3]);

    // output has same dimensions as size grid
    grid_out = grid_size;

    // allocate grid
    grid_out.buffer = AllocateGrid(&grid_out);
    if (grid_out.buffer == NULL) {
        nll_puterr(
                "ERROR: allocating memory for output grid buffer.\n");
        return (-1);
    }

    // create grid array access pointers
    grid_out.array = CreateGridArray(&grid_out);
    if (grid_out.array == NULL) {
        nll_puterr(
                "ERROR: creating array for accessing output grid buffer.\n");
        return (-1);
    }

    // initialize output grid to all zeros

    InitializeGrid(&grid_out, 0.0F);



    // sum requested grid files into output grid

    // input file name
    strcpy(fn_grid_in, argv[4]);

    // check for wildcards in input file name
    strcat(fn_grid_in, test_str);
    numFiles = ExpandWildCards(fn_grid_in, fn_grid_in_list, MAX_NUM_INPUT_FILES);

    // get additional input file (e.g. for ADD multiple arguments or SUB difference)
    // 20190509 AJL  if (argc > 5) {
    int narg = 5;
    if (process_type != PROC_MUL_CONST && numFiles == 1 && argc > 5) {
        while (numFiles < MAX_NUM_INPUT_FILES && narg < argc) {
            strcpy(fn_grid_in_list[numFiles], argv[narg]);
            narg++;
            numFiles++;
        }
    }

    if (process_type == PROC_SUB && numFiles > 2) {
        sprintf(MsgStr, "ERROR: cannot difference more than 2 grids.");
        nll_puterr(MsgStr);
    }

    // get factor (e.g. for SUB difference)
    // 20190509 AJL  if (argc > 5) {
    double mul_const;
    if (process_type == PROC_MUL_CONST && numFiles == 1 && argc > 5) {
        if (sscanf(argv[5], "%lf", &mul_const) < 1) {
            nll_puterr("ERROR opening size grid file.");
            return (-1);
        }
    }


    if (process_type == PROC_ADD || process_type == PROC_SUB) {
        double factor = 1.0;
        for (nFile = 0; nFile < numFiles; nFile++) {
            chrpos = strstr(fn_grid_in_list[nFile], test_str);
            if (chrpos != NULL)
                *chrpos = '\0';
            if (factor > 0.0) {
                fprintf(OUT_LEVEL_1, "Adding grid file <%s>\n", fn_grid_in_list[nFile]);
            } else {
                fprintf(OUT_LEVEL_1, "Subtracting grid file <%s>\n", fn_grid_in_list[nFile]);
            }
            if ((istat = OpenGrid3dFile(fn_grid_in_list[nFile],
                    &fp_grid_in, &fp_grid_in_hdr,
                    &grid_in, "", &srce, 0)) < 0) {
                sprintf(MsgStr, "ERROR: opening grid file <%s>", fn_grid_in_list[nFile]);
                nll_puterr(MsgStr);
            } else {
                if (testIdentical(&grid_out, &grid_in)) {
                    // allocate grid
                    grid_in.buffer = AllocateGrid(&grid_in);
                    if (grid_in.buffer == NULL) {
                        nll_puterr(
                                "ERROR: allocating memory for input grid buffer.\n");
                        return (-1);
                    }
                    // create grid array access pointers
                    grid_in.array = CreateGridArray(&grid_in);
                    if (grid_in.array == NULL) {
                        nll_puterr(
                                "ERROR: creating array for accessing input grid buffer.\n");
                        return (-1);
                    }
                    // read grid
                    if ((istat =
                            ReadGrid3dBuf(&grid_in, fp_grid_in)) < 0) {
                        nll_puterr("ERROR: reading input grid from disk.");
                        return (-1);
                    }
                }
                if ((istat = SumGrids(&grid_out, &grid_in, fp_grid_in, factor)) < 0) {
                    sprintf(MsgStr, "ERROR: processing grid file <%s>.", fn_grid_in_list[nFile]);
                    nll_puterr(MsgStr);
                }
                if (nFile == 0) {
                    grid_out.type = grid_in.type;
                    strcpy(grid_out.chr_type, grid_in.chr_type);
                }
                // now we have added left grid, can subtract right grid
                if (process_type == PROC_SUB) {
                    factor = -1.0;
                }
            }
            fclose(fp_grid_in);
            DestroyGridArray(&grid_in);
            FreeGrid(&grid_in);

        }
    } else if (process_type == PROC_MUL_CONST) {
        nFile = 0;
        chrpos = strstr(fn_grid_in_list[nFile], test_str);
        if (chrpos != NULL)
            *chrpos = '\0';
        fprintf(OUT_LEVEL_1, "Multiplying grid file <%s> by %f\n", fn_grid_in_list[nFile], mul_const);
        if ((istat = OpenGrid3dFile(fn_grid_in_list[nFile],
                &fp_grid_in, &fp_grid_in_hdr,
                &grid_in, "", &srce, 0)) < 0) {
            sprintf(MsgStr, "ERROR: opening grid file <%s>", fn_grid_in_list[nFile]);
            nll_puterr(MsgStr);
        } else {
            if (testIdentical(&grid_out, &grid_in)) {
                // allocate grid
                grid_in.buffer = AllocateGrid(&grid_in);
                if (grid_in.buffer == NULL) {
                    nll_puterr(
                            "ERROR: allocating memory for input grid buffer.\n");
                    return (-1);
                }
                // create grid array access pointers
                grid_in.array = CreateGridArray(&grid_in);
                if (grid_in.array == NULL) {
                    nll_puterr(
                            "ERROR: creating array for accessing input grid buffer.\n");
                    return (-1);
                }
                // read grid
                if ((istat =
                        ReadGrid3dBuf(&grid_in, fp_grid_in)) < 0) {
                    nll_puterr("ERROR: reading input grid from disk.");
                    return (-1);
                }
            }
            if ((istat = MulConstGrid(&grid_out, &grid_in, fp_grid_in, mul_const)) < 0) {
                sprintf(MsgStr, "ERROR: processing grid file <%s>.", fn_grid_in_list[nFile]);
                nll_puterr(MsgStr);
            }
            if (nFile == 0) {
                grid_out.type = grid_in.type;
                strcpy(grid_out.chr_type, grid_in.chr_type);
            }
        }
        fclose(fp_grid_in);
        DestroyGridArray(&grid_in);
        FreeGrid(&grid_in);
    }


    // save sum grid to disk

    char filename_type[MAXLINE] = "ERROR";
    if (process_type == PROC_ADD)
        strcpy(filename_type, "sum");
    else if (process_type == PROC_SUB)
        strcpy(filename_type, "diff");
    else if (process_type == PROC_MUL_CONST)
        strcpy(filename_type, "mul_const");
    if ((istat = WriteGrid3dBuf(&grid_out, &srce, fn_grid_out, filename_type)) < 0) {
        sprintf(MsgStr, "ERROR: writing sum grid to disk: <%s>.", fn_grid_out);
        nll_puterr(MsgStr);
        return (-1);
    }

    fprintf(OUT_LEVEL_1, "Number grids processed: %d\n", numFiles);
    fprintf(OUT_LEVEL_1, "Output grid: <%s>\n", fn_grid_out);

    fclose(fp_grid_size);
    fclose(fp_grid_size_hdr);
    fclose(fp_grid_in_hdr);

    return (0);

}

