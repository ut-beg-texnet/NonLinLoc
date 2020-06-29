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




/*   Grid2GMT.c

        Program to extract cross sections from 3D grid files

        output in GMT readable format

 */

/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
161 Allee du Micocoulier, 06370 Mouans-Sartoux, France
tel: +33(0)493752502  e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


/*
        history:

        ver 01    16SEP1997  AJL  Original version
        ver 02      OCT2000  AJL  Added functionality:
        ver 03      APR2001  AJL  Added functionality:


.........1.........2.........3.........4.........5.........6.........7.........8

 */



#include "GridLib.h"
#include "GridGraphLib.h"

#define GMT_VER_3_3_4	1



#define NUM_COLORS 8
#define MISFIT_CONTOUR_INTERVAL_XY 1.0
#define PLOT_WIDTH 20.0
// 8.5x11 inches
//#define PLOT_HEIGHT 27.0
// A4
#define PLOT_HEIGHT 28.7
#define TITLE_FONT_SIZE 17
#define TITLE_FONT 4
#define HYPO_FONT_SIZE 13
#define HYPO_FONT 4
#define STA_FONT_SIZE 8
#define STA_FONT 4
#define ANNOTATION_FONT_SIZE 7
#define ANNOTATION_FONT 4

#define MAX_GRID_VALUE 			1.0e8
#define SUBSTITUTE_MAX_GRID_VALUE 	-1.0

#define DEFAULT_LENGTH_UNITS "km"
#define DEFAULT_LABEL_CONTOURS " -A- "


// !!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!
// custom user options
// set defaults
#define PLOT_LEGEND_LINE 1





/* globals */

// command line argument flags
char lengthUnits[64] = DEFAULT_LENGTH_UNITS;
char title[MAXLINE] = "\0";
int reverse_xy = 0;
int plot_lat_lon = 0;
char label_contours[64] = DEFAULT_LABEL_CONTOURS;

char fninput[FILENAME_MAX], fnoutput[FILENAME_MAX];
#define MAX_NUM_INPUT_FILES 500
char fnroot_input_list[MAX_NUM_INPUT_FILES][FILENAME_MAX];
char fnroot_input[FILENAME_MAX];
char fngrid_control[FILENAME_MAX], fn_root_output[FILENAME_MAX], fn_ps_output[FILENAME_MAX];
char fn_gmt[FILENAME_MAX], fn_gmtgrd[FILENAME_MAX], fn_cont[FILENAME_MAX];
FILE *fp_grid, *fp_grid_control, *fp_hdr, *fp_gmt;
GridDesc grid0;
int iFirstPlot;
float grid_value_max, grid_value_min;
Vect3D Expectation, MaxLike;
Ellipsoid3D Ellipsoid;
/* stations */
int NumStations;
SourceDesc Stations[X_MAX_NUM_ARRIVALS];
// GMT JVAL, RVAL for lat/long region
char gmt_JVAL_latlong_string[MAXLINE];
char gmt_RVAL_latlong_string[MAXLINE];
// flag to indicate if map grid read from control file.  map grid overrides hypo-file grids
int mapGridRead;

int proj_index_input = -1;
int proj_index_output = -1;
int proj_index_dummy = 1;


/* functions */

int ReadGrid2GMT_Input(FILE*);
double GetContourInterval(double, double, int, int*);
int MakeConfCPT(char*, char*);
int MakeMisfitCPT(char*, double, double);
int GenGMTCommands(char, char, char[][20], int, char*, int, int, int, int, int, GridDesc *);
int MakeTopoCPT(char* fileout);
int GenGridViewGMT(GridDesc*, char, char, char[][20], int,
        int, int, int, int, int,
        double, char *, char*, char*, double *, double *, double *,
        double, double, double, double, int);
int MapFiles2GMT(double, double, double, double, FILE*, char*, int, int);
int MapLines2GMT(int, double, double, double, double, FILE*, char*, int);
int grd2GMT(int, double, double, double, double, FILE*, char*, int);
int GenMapFileName(char *, char *, char *, char *, char *, char *);
int Scat2GMT(char*, char*, int, char*);
double CalcAngleValue(double, int, int);
int PlotTraditionStats(char, char, double, Vect3D*, Vect3D*, Ellipsoid3D *,
        FILE*, char*, char*);
int PlotFocalMechanism(char view_type, double scale, double magnitude, FocalMech *pfocMech,
        FILE* fp_io, char *plt_code, char *linecolor, char *fillcolor,
        double horiz_min, double horiz_max, double vert_min, double vert_max, double plot_scale);
int addStations(SourceDesc *stations, int numStations, ArrivalDesc *arrival, int nArrivals);
int ConvertResids2MapGMT(char* fn_nlloc_stat, char* phaseID, FILE* fp_out,
        SourceDesc* staList, int nstations, double, int);
void genResidualGMT(FILE* fp_out, char* xtra_args, double resid,
        double x, double y, double scale);
SourceDesc *findStaLoc(char *staName, SourceDesc* stalist, int nstations);
int parseCommandArgument(char* arg, char* pcommandchr, char arg_elements[][20]);
void usage();
int parameter_proc(int argcount, char **argvec);
int ReadStationList(FILE* fpio, SourceDesc *stations, int numStations);



/*** program to generate  GMT command files for plotting grid data */

#define PNAME  "Grid2GMT"

int main(int argc, char *argv[]) {

    int istat, narg, n;
    char cplotmode, clatlongmode, cdatatype;
    char arg_elements[10][20];
    int num_arg_elements;
    char args_str[MAXLINE] = "\0";
    int ix1, iy1, ix2, iy2, izlevel;
    char test_str[] = ".hdr";
    char *chrpos;
    int NumFiles, nFile;
    double vlat1, vlong1, vlat2, vlong2;
    double xrect, yrect;

    GridDesc *pgrid_plot;


    /* set program name */

    strcpy(prog_name, PNAME);


    /* check command line for correct usage */

    if (message_flag > 0) {
        fprintf(stdout, "\n%s Arguments: ", prog_name);
        for (narg = 0; narg < argc; narg++)
            fprintf(stdout, "<%s> ", argv[narg]);
        fprintf(stdout, "\n");
    }

    clatlongmode = '\0';
    if (argc > 4)
        sscanf(argv[4], "%c%c", &cplotmode, &clatlongmode);
    else
        cplotmode = '\0';
    cplotmode = toupper(cplotmode);
    if (argc < 6 || cplotmode == '\0' || (cplotmode == 'V' && argc < 10)
            || (cplotmode == 'H' && argc < 7)
            || (cplotmode == 'L' && argc < 6)) {
        nll_puterr("ERROR wrong number of command line arguments.");
        usage();
        exit(-1);
    }


    /* read command line arguments */

    narg = 0;
    strcpy(fngrid_control, argv[++narg]);
    strcpy(fninput, argv[++narg]);
    strcpy(fnoutput, argv[++narg]);
    ++narg;
    // parse plottype
    num_arg_elements = parseCommandArgument(argv[++narg], &cdatatype, arg_elements);

    if (cplotmode == 'V') {
        if (clatlongmode == 'L') { // lat/long specification of ends
            istat = 0;
            istat += sscanf(argv[++narg], "%lf", &vlat1);
            istat += sscanf(argv[++narg], "%lf", &vlong1);
            istat += sscanf(argv[++narg], "%lf", &vlat2);
            istat += sscanf(argv[++narg], "%lf", &vlong2);
        } else { // grid index specification of ends
            istat = 0;
            istat += sscanf(argv[++narg], "%d", &ix1);
            istat += sscanf(argv[++narg], "%d", &iy1);
            istat += sscanf(argv[++narg], "%d", &ix2);
            istat += sscanf(argv[++narg], "%d", &iy2);
        }
        sprintf(args_str, "%s_%s_->_%s_%s", argv[6], argv[7], argv[8], argv[9]);
    } else if (cplotmode == 'H') {
        sscanf(argv[++narg], "%d", &izlevel);
        sprintf(args_str, "%s", argv[6]);
    } else if (cplotmode == 'L') {
        if (argc == 9) {
            ix1 = iy1 = izlevel = 0;
            sscanf(argv[++narg], "%d", &ix1);
            sscanf(argv[++narg], "%d", &iy1);
            sscanf(argv[++narg], "%d", &izlevel);
            sprintf(args_str, "x%s_y%s_z%s", argv[6], argv[7], argv[8]);
        } else if (argc == 7) {
            sprintf(args_str, "z%s", argv[6]);
        }
    }
    // check for length units
    if (argc - narg > 0)
        parameter_proc(argc - narg, argv + narg);




    /* set constants */

    SetConstants();
    SetGraphicConstants();
    NumStations = 0;

    /* read control file */

    strcpy(fngrid_control, argv[1]);
    if ((fp_grid_control = fopen(fngrid_control, "r")) == NULL) {
        nll_puterr("ERROR opening control file.");
        exit(EXIT_ERROR_FILEIO);
    }

    proj_index_input = -1;
    proj_index_output = -1;
    if ((istat = ReadGrid2GMT_Input(fp_grid_control)) < 0) {
        exit(EXIT_ERROR_FILEIO);
    }
    fclose(fp_grid_control);
    if (proj_index_input < 0)
        proj_index_input = proj_index_output;
    else if (proj_index_output < 0)
        proj_index_output = proj_index_input;


    /* check for wildcards in input file name */
    strcat(fninput, test_str);
    NumFiles = ExpandWildCards(fninput, fnroot_input_list, MAX_NUM_INPUT_FILES);


    for (nFile = 0; nFile < NumFiles; nFile++) {
        chrpos = strstr(fnroot_input_list[nFile], test_str);
        *chrpos = '\0';
        strcpy(fnroot_input, fnroot_input_list[nFile]);
        if ((chrpos = strrchr(fnroot_input, '/')) == NULL)
            chrpos = fnroot_input;
        else
            chrpos++;
        sprintf(fn_root_output, "%s%s.%c%c", fnoutput, chrpos,
                cplotmode, cdatatype);
        for (n = 0; n < num_arg_elements; n++) {
            strcat(fn_root_output, "_");
            strcat(fn_root_output, arg_elements[n]);
        }

        if (message_flag > 0) {
            fprintf(stdout, "\nProcessing files %s.*\n",
                    fnroot_input);
            fprintf(stdout, "Output files %s.*\n", fn_root_output);
        }

        /* open grid file and header file and read header file */
        grid0.iSwapBytes = 0;
        // swap bytes if byte order problem
        //grid0.iSwapBytes = 1;
        if ((istat = OpenGrid3dFile(fnroot_input, &fp_grid, &fp_hdr, &grid0, "", NULL,
                grid0.iSwapBytes)) < 0) {
            nll_puterr2("ERROR opening grid file: ", fnroot_input);
            exit(EXIT_ERROR_FILEIO);
        }
        /* convert lat/long specification of ends to grid index */
        if (clatlongmode == 'L') {
            latlon2rect(proj_index_output, vlat1, vlong1, &xrect, &yrect);
            ix1 = (int) (xrect - grid0.origx) / grid0.dx;
            iy1 = (int) (yrect - grid0.origy) / grid0.dy;
            latlon2rect(proj_index_output, vlat2, vlong2, &xrect, &yrect);
            ix2 = (int) (xrect - grid0.origx) / grid0.dx;
            iy2 = (int) (yrect - grid0.origy) / grid0.dy;
            if (ix1 < 0 || ix1 >= grid0.numx || iy1 < 0 || iy1 >= grid0.numy ||
                    ix2 < 0 || ix2 >= grid0.numx || iy2 < 0 || iy2 >= grid0.numy) {
                nll_puterr("ERROR lat long points not inside grid.");
                exit(EXIT_ERROR_MISC);
            }
        }

        // set grid for plotting
        if (mapGridRead) {
            pgrid_plot = &grid_in;
        } else {
            pgrid_plot = &grid0;
        }

        if ((istat = GenGMTCommands(cplotmode, cdatatype,
                arg_elements, num_arg_elements, args_str,
                ix1, iy1, ix2, iy2, izlevel, pgrid_plot)) < 0) {
            exit(EXIT_ERROR_MISC);
        }

        CloseGrid3dFile(&fp_grid, &fp_hdr);
    }


    exit(0);

}

/** function to parse a GMT style program argument */

int parseCommandArgument(char* arg, char* pcommandchr, char arg_elements[][20]) {

    int narg, n;
    char* chr;


    printf("args: <%s>\n", arg);

    // commmand name is first character of program argument
    *pcommandchr = toupper(arg[0]);

    // copy commmand arguments part of program argument to first argument element
    strcpy(arg_elements[0], arg + 1);

    if (strlen(arg_elements[0]) == 0)
        return (0);

    // get remaining argument elements
    narg = 1;
    while ((chr = strchr(arg_elements[narg - 1], '/')) != NULL) {
        *chr = '\0';
        strcpy(arg_elements[narg++], chr + 1);
    }
    for (n = 0; n < narg; n++)
        printf("narg %d arg: <%s>\n", n, arg_elements[n]);


    return (narg);

}

/* prgram usage messages
 */
void usage() {
            disp_usage(PNAME,
                "<controlfile> <gridroot> <outroot> V [G][S][Ennndx][M] ix1 iy1 ix2 iy2 [length_units]");
        disp_usage(PNAME,
                "<controlfile> <gridroot> <outroot> VL [G][S][Ennndx][M] lat1 long1 lat2 long2 [length_units]");
        disp_usage(PNAME,
                "<controlfile> <gridroot> <outroot> H [G][S][Ennndx][M][Rphases/scale] iz [length_units]");
        disp_usage(PNAME,
                "<controlfile> <gridroot> <outroot> L [G][S][Ennndx][M][Rphases/scale] [ix iy iz] [length_units]");
}

/** *************************************************************************
 * parameter_proc:
 *
 * Process the command line parameter flags.
 *
 * Returns 0 on success, and -1 on failure
 ************************************************************************* **/

int parameter_proc(int argcount, char **argvec) {

    int optind;
    int error = 0;

    if (argcount <= 1)
        error++;

    // Process all command line arguments
    for (optind = 1; optind < argcount; optind++) {
        if (strcmp(argvec[optind], "-h") == 0) {
            usage();
            exit(0);
        } else if (strcmp(argvec[optind], "-title") == 0) {
            strcpy(title, argvec[++optind]);
        } else if (strcmp(argvec[optind], "-length-units") == 0) {
            strcpy(lengthUnits, argvec[++optind]);
        } else if (strcmp(argvec[optind], "-reverse-xy") == 0) {
            reverse_xy = 1;
        } else if (strcmp(argvec[optind], "-plot-lat-lon") == 0) {
            plot_lat_lon = 1;
        } else if (strcmp(argvec[optind], "-label_contours") == 0) {
            strcpy(label_contours, " -A ");
        } else if (strncmp(argvec[optind], "-", 1) == 0) {
            fprintf(stderr, "Unknown option: %s\n", argvec[optind]);
            exit(1);
        } else {
            fprintf(stderr, "Unknown option: %s\n", argvec[optind]);
            exit(1);
        }
    }


    return 0;
} // End of parameter_proc()

/** function to generate GMT commands for grid plotting */

int GenGMTCommands(char cplotmode, char cdatatype,
        char arg_elements[][20], int num_arg_elements, char *args_str,
        int ix1, int iy1, int ix2, int iy2, int izlevel, GridDesc *pgrid0) {

    int istat, ihypo;
    char shift_str[MAXLINE], title_str[MAXLINE] = "\0";
    char signature_str[2 * MAXLINE];
    char cpt_str[MAXLINE], scale_label_str[MAXLINE];
    double scale, horiz_width, xlen, ylen, xshift, yshift;
    double xshift_cum = 0.0, yshift_cum = 0.0;
    double scaleshift = 0.0, sxshift = 0.0, syshift = 0.0, plot_width;

    HypoDesc Hypo, hypoDummy;
    int iMultEvent;
    char fn_hypo[FILENAME_MAX];
    FILE* fp_hypo;

    char hypotext[MAXLINE];
    char sys_string[MAXLINE];

    double res_scale;
    int res_min_num_readings;

    double xres, yres, xres_shift, res_legend_mag;
    char res_legend_mag_string[20];


    /* open gmt files */

    sprintf(fn_gmt, "%s.gmt", fn_root_output);
    if ((fp_gmt = fopen(fn_gmt, "w")) == NULL) {
        nll_puterr("ERROR opening gmt output file.");
        return (-1);
    }


    /* write gmt script file */

    fprintf(fp_gmt, "#!/bin/csh\n#\n#\n\n");

    // gmtdefaults
    fprintf(fp_gmt, "gmtset  PAGE_ORIENTATION portrait  X_ORIGIN 0.5  Y_ORIGIN 0.5 \n\n");
    fprintf(fp_gmt, "gmtset  ANNOT_FONT_SIZE_PRIMARY 8  ANNOT_FONT_SIZE_SECONDARY 6  HEADER_FONT_SIZE 12 LABEL_FONT_SIZE 10\n\n");
    fprintf(fp_gmt, "gmtset  LABEL_OFFSET 0.1c  ANNOT_OFFSET_PRIMARY 0.1c ANNOT_OFFSET_SECONDARY 0.1c\n\n");

    fprintf(fp_gmt, "set POSTSCRIPT_NAME = %s\n\n", fn_root_output);
    fprintf(fp_gmt, "unalias rm\n\n");
    fprintf(fp_gmt, "rm -f $POSTSCRIPT_NAME.ps\n\n");
    sprintf(fn_ps_output, "$POSTSCRIPT_NAME");


    fprintf(fp_gmt,
            "# ========================================================\n");
    fprintf(fp_gmt,
            "# Optional X/Y or LAT/LONG plotting\n");
    fprintf(fp_gmt,
            "# NOTE: LAT/LONG plotting works only for unrotated, Horizontal section (H) plots\n");
    fprintf(fp_gmt,
            "#     for X/Y plot, uncomment the folowing line:\n");
    if (!plot_lat_lon)
        fprintf(fp_gmt, "set PLOT_LAT_LONG = 0\n");
    else
        fprintf(fp_gmt, "#set PLOT_LAT_LONG = 0\n");
    //if (doLatLong) {
    fprintf(fp_gmt,
            "#     for LAT/LONG plot, uncomment the folowing line:\n");
    if (plot_lat_lon)
        fprintf(fp_gmt, "set PLOT_LAT_LONG = 1\n");
    else
        fprintf(fp_gmt, "#set PLOT_LAT_LONG = 1\n");
    //}
    fprintf(fp_gmt,
            "# ========================================================\n");
    fprintf(fp_gmt, "\n");


    /* begin plot */

    fprintf(fp_gmt,
            "psbasemap -R0.5/1.0/0.5/1.0 -Bf10 -JX%lf/%lf -K > %s.ps\n\n",
            PLOT_WIDTH, PLOT_HEIGHT, fn_ps_output);


    /* plot view(s) */

    grid_value_max = -VERY_LARGE_FLOAT, grid_value_min = VERY_LARGE_FLOAT;
    iFirstPlot = 1;

    Hypo.ix = ix1;
    Hypo.iy = iy1;
    Hypo.iz = izlevel;
    Hypo.x = Hypo.y = Hypo.z = 0.0;
    Hypo.grid_misfit_max = 1.0;
    strcpy(Hypo.signature, "");

    if (cdatatype == 'G') {
        if (strlen(title) < 1)
            sprintf(title_str, "%s__(%s)", fn_root_output, args_str);
        else
            sprintf(title_str, "%s", title);
    }

    //} else {

    /* open hypocenter file */

    sprintf(fn_hypo, "%s.hyp", fnroot_input);
    if ((fp_hypo = fopen(fn_hypo, "r")) == NULL) {
        if (message_flag >= 1)
            nll_putmsg2(1, "INFO: cannot open hypocenter file", fn_hypo);
    }

    strcpy(Hypo.comment, "---");


    if (fp_hypo != NULL && ((ihypo = GetHypLoc(fp_hypo, NULL, &Hypo, Arrival,
            &NumArrivals, 1, NULL, proj_index_dummy)) == 0 || ihypo == EOF)) {

        if (message_flag >= 2)
            WriteLocation(stdout, &Hypo, Arrival, NumArrivals,
                NULL, 1, 1, 0, NULL, proj_index_output);

        fprintf(fp_gmt,
                "pstext -R0.5/1.0/0.5/1.0 -Bf10 -JX%lf/%lf -K -O << END >> %s.ps\n%lf %lf %d %d %d %d %s\nEND\n\n",
                PLOT_WIDTH, PLOT_HEIGHT, fn_ps_output,
                0.75, 0.98, TITLE_FONT_SIZE,
                0, TITLE_FONT, 2, Hypo.comment);

        /* check for multiple events */
        ihypo = GetHypLoc(fp_hypo, NULL, &hypoDummy, Arrival,
                &NumArrivals, 1, NULL, proj_index_dummy);
        // determine if multi event file
        iMultEvent = ihypo == EOF ? 0 : 1;
        if (!iMultEvent) {
            sprintf(hypotext,
                    "%4.4d %2.2d %2.2d  %2.2d:%2.2d:%7.4lf  Lat: %.3lf  Long: %.3lf  Z: %.2lfkm  RMS: %.2lfs  Mag: %.1lf",
                    Hypo.year, Hypo.month, Hypo.day, Hypo.hour, Hypo.min,
                    (double) Hypo.sec, Hypo.dlat, Hypo.dlong, Hypo.depth, Hypo.rms, Hypo.amp_mag);
            MaxLike.x = Hypo.x;
            MaxLike.y = Hypo.y;
            MaxLike.z = Hypo.z;
            Expectation = Hypo.expect;
            Ellipsoid = Hypo.ellipsoid;
        } else {
            sprintf(hypotext, ".");
            Hypo.x = Hypo.y = Hypo.z = -LARGE_DOUBLE;
            MaxLike.x = Hypo.x;
            MaxLike.y = Hypo.y;
            MaxLike.z = Hypo.z;
            Hypo.ix = Hypo.iy = Hypo.iz = -1;
            Expectation.x = Expectation.y = Expectation.z
                    = -LARGE_DOUBLE;
            Hypo.grid_misfit_max = 999.0;
        }

        fprintf(fp_gmt,
                "pstext -R0.5/1.0/0.5/1.0 -Bf10 -JX%lf/%lf -K -O << END >> %s.ps\n%lf %lf %d %d %d %d %s\nEND\n\n",
                PLOT_WIDTH, PLOT_HEIGHT, fn_ps_output,
                0.75, 0.97, HYPO_FONT_SIZE,
                0, HYPO_FONT, 2, hypotext);

        fclose(fp_hypo);
    }        // no hypocenter file
    else {
        if (strlen(title) < 1)
            sprintf(title_str, "%s__(%s)", fn_root_output, args_str);
        else
            sprintf(title_str, "%s", title);
    }


    // get residual scale
    if (cdatatype == 'R') {
        istat = 0;
        if (num_arg_elements <= 1)
            res_scale = 0.5;
        else
            istat = sscanf(arg_elements[1], "%lf", &res_scale);
        printf("arg_elements[1] <%s>  RES SCALE: %lf\n", arg_elements[1], res_scale);
        if (istat < 0)
            nll_puterr("ERROR: Reading 2nd argument of 'R' datatype.");
        istat = 0;
        if (num_arg_elements <= 2)
            res_min_num_readings = 1;
        else
            istat = sscanf(arg_elements[2], "%d", &res_min_num_readings);
        printf("arg_elements[2] <%s>  RES MIN NUM RADINGS: %d\n", arg_elements[2], res_min_num_readings);
        if (istat < 0)
            nll_puterr("ERROR: Reading 3rd argument of 'R' datatype.");

    }

    if (PLOT_LEGEND_LINE) {

        if (cdatatype == 'S') {
            fprintf(fp_gmt,
                    "pstext -R0.5/1.0/0.5/1.0 -Bf10 -JX%lf/%lf -K -O << END >> %s.ps\n%lf %lf %d %d %d %d %s\nEND\n\n",
                    PLOT_WIDTH, PLOT_HEIGHT, fn_ps_output,
                    0.75, 0.96, HYPO_FONT_SIZE,
                    0, TITLE_FONT, 2, "PDF scatter sample");
        } else if (cdatatype == 'E') {
            fprintf(fp_gmt,
                    "pstext -R0.5/1.0/0.5/1.0 -Bf10 -JX%lf/%lf -K -O << END >> %s.ps\n%lf %lf %d %d %d %d %s\nEND\n\n",
                    PLOT_WIDTH, PLOT_HEIGHT, fn_ps_output,
                    0.75, 0.96, HYPO_FONT_SIZE,
                    0, TITLE_FONT, 2,
                    "Error Ellipsoid, Expectation (dot) and Maximum Likelihood (star)");
        } else if (cdatatype == 'M') {
            fprintf(fp_gmt,
                    "pstext -R0.5/1.0/0.5/1.0 -Bf10 -JX%lf/%lf -K -O << END >> %s.ps\n%lf %lf %d %d %d %d %s\nEND\n\n",
                    PLOT_WIDTH, PLOT_HEIGHT, fn_ps_output,
                    0.75, 0.96, HYPO_FONT_SIZE,
                    0, TITLE_FONT, 2, "Double-couple focal mechanisms");
        } else if (cdatatype == 'R') {
            fprintf(fp_gmt,
                    "pstext -R0.5/1.0/0.5/1.0 -Bf10 -JX%lf/%lf -K -O << END >> %s.ps\n%lf %lf %d %d %d %d %s %s\nEND\n\n",
                    PLOT_WIDTH, PLOT_HEIGHT, fn_ps_output,
                    0.75, 0.96, HYPO_FONT_SIZE,
                    0, TITLE_FONT, 2, arg_elements[0], "Residuals");
        }

    } // PLOT_LEGEND_LINE


    /* draw views */

    if (cplotmode == 'L') {

        /* calculate plot scale */

        // calculate scale with room for horiz view and section
        horiz_width = (double) (pgrid0->numx - 1) *
                pgrid0->dx + (double) (pgrid0->numz - 1) * pgrid0->dz;
        scale = (25.0 / 32.0) * PLOT_WIDTH / horiz_width;
        // calculate horiz width of horiz view
        plot_width = scale * (double) (pgrid0->numx - 1) * pgrid0->dx;


        xshift = 1.75;
        yshift = 3.75;
        sprintf(shift_str, "-X%lf -Y%lf", xshift, yshift);
        //if (cdatatype != 'R')
        GenGridViewGMT(pgrid0, 'L', cdatatype, arg_elements, num_arg_elements,
                0, Hypo.iy, pgrid0->numx,
                Hypo.iy, 0,
                plot_width, "WEnS", "", shift_str,
                &scale, &xlen, &ylen, Hypo.z, Hypo.x,
                Hypo.grid_misfit_max, res_scale, res_min_num_readings);
        xshift_cum += xshift;
        yshift_cum += yshift;

        xshift = 0.0;
        yshift = ylen + PLOT_WIDTH / 28.0;
        sprintf(shift_str, "-X%lf -Y%lf", xshift, yshift);
        GenGridViewGMT(pgrid0, 'H', cdatatype, arg_elements, num_arg_elements,
                0, 0, 0, 0, Hypo.iz,
                plot_width, "WeNs", title_str, shift_str,
                &scale, &xlen, &ylen, Hypo.y, Hypo.x,
                Hypo.grid_misfit_max, res_scale, res_min_num_readings);
        xshift_cum += xshift;
        yshift_cum += yshift;

        xshift = xlen + PLOT_WIDTH / 28.0;
        sprintf(shift_str, "-X%lf -Y%lf", xshift, 0.0);
        //if (cdatatype != 'R')
        GenGridViewGMT(pgrid0, 'L', cdatatype, arg_elements, num_arg_elements,
                Hypo.ix, 0, Hypo.ix,
                pgrid0->numy, 0,
                plot_width, "wENS", "", shift_str,
                &scale, &xlen, &ylen, Hypo.y, Hypo.z,
                Hypo.grid_misfit_max, res_scale, res_min_num_readings);
        xshift_cum += xshift;

        scaleshift = -1.75;
        sxshift = -xshift;
        syshift = -yshift - 2.5;



    } else if (cplotmode == 'H' || cplotmode == 'V') {

        scale = -1.0;
        xshift = 1.75;
        yshift = 3.75;
        sprintf(shift_str, "-X%lf -Y%lf", xshift, yshift);
        //sprintf(title_str, "%s__(%s)", fn_root_output, args_str);
        GenGridViewGMT(pgrid0, cplotmode, cdatatype, arg_elements, num_arg_elements,
                ix1, iy1, ix2, iy2,
                izlevel,
                (26.0 / 32.0) * PLOT_WIDTH, "WESN", title_str, shift_str,
                &scale, &xlen, &ylen, -1.0, -1.0,
                Hypo.grid_misfit_max, res_scale, res_min_num_readings);
        xshift_cum += xshift;
        yshift_cum += yshift;

        scaleshift = 0.0;
        sxshift = -xshift;
        syshift = -yshift + 1.5;
    }


    /* plot scales */

    if (cdatatype != 'S' && cdatatype != 'E' && cdatatype != 'M' && cdatatype != 'R') {

        if (pgrid0->type == GRID_PROB_DENSITY) {
            strcpy(cpt_str, "conflev.");
            strcpy(scale_label_str, "Confidence Level");
        } else if (pgrid0->type == GRID_MISFIT) {
            strcpy(cpt_str, "misfit.");
            strcpy(scale_label_str, "RMS Misfit (sec)");
        } else if (pgrid0->type == GRID_LIKELIHOOD) { // relative likelihood normalised 0-1
            strcpy(cpt_str, "");
            strcpy(scale_label_str, "Relative Likelihood");
        } else if (pgrid0->type == GRID_VELOCITY) {
            strcpy(cpt_str, "");
            strcpy(scale_label_str, "Velocity (km/sec)");
        } else if (pgrid0->type == GRID_VELOCITY_METERS) {
            strcpy(cpt_str, "");
            strcpy(scale_label_str, "Velocity (m/sec)");
        } else if (pgrid0->type == GRID_SLOW2_METERS) {
            strcpy(cpt_str, "");
            strcpy(scale_label_str, "Slowness**2 (sec/m)**2");
        } else if (pgrid0->type == GRID_SLOW_LEN) {
            strcpy(cpt_str, "");
            strcpy(scale_label_str, "Cell Transit Time (sec)");
        } else if (pgrid0->type == GRID_TIME ||
                pgrid0->type == GRID_TIME_2D) {
            strcpy(cpt_str, "");
            strcpy(scale_label_str, "Time (sec)");
        } else if (pgrid0->type == GRID_DEPTH) {
            strcpy(cpt_str, "");
            sprintf(scale_label_str, "Depth (%s)", lengthUnits);
        } else if (pgrid0->type == GRID_LENGTH) {
            strcpy(cpt_str, "");
            sprintf(scale_label_str, "Length (%s)", lengthUnits);
        } else if (pgrid0->type == GRID_ANGLE ||
                pgrid0->type == GRID_ANGLE_2D) {
            strcpy(cpt_str, "");
            if (cplotmode == 'V')
                strcpy(scale_label_str,
                    "Dip Angle(deg) (>180 = unstable)");
            if (cplotmode == 'H')
                strcpy(scale_label_str,
                    "Strike Angle(deg) (<0 = unstable)");
        } else if (pgrid0->type == GRID_INCLINATION ||
                pgrid0->type == GRID_INCLINATION_2D) {
            strcpy(cpt_str, "");
                strcpy(scale_label_str,
                    "Dip Angle(deg) (>180 = unstable)");
        } else {
            strcpy(cpt_str, "");
            strcpy(scale_label_str, "UNKNOWN SCALE");
        }

        fprintf(fp_gmt,
                "psscale -C%s.%scpt -D%lf/%lf/14.0/0.5h -B:\"%s\": -X%lf -Y%lf $SCALE_FLAG -K -O >> %s.ps\n\n",
                fn_root_output, cpt_str,
                (1.0 / 2.0) * PLOT_WIDTH + scaleshift, 0.75,
                scale_label_str, sxshift, syshift, fn_ps_output);
        xshift_cum += sxshift;
        yshift_cum += syshift;
    }

    if (cdatatype == 'R') {

        // plot legend

        res_legend_mag = 0.4 / res_scale;
        //yres = vert_min - 35.0 * plot_scale;
        yres = grid0.origy + (-yshift_cum + 1.25) / scale;
        xres_shift = 2.0 * res_legend_mag * res_scale / scale;
        // positive residual
        xres = grid0.origx + (-xshift_cum + PLOT_WIDTH / 2.0) / scale;
        genResidualGMT(fp_gmt, "-N", res_legend_mag, xres + xres_shift, yres, res_scale);
        sprintf(res_legend_mag_string, "+%.2f sec", res_legend_mag);
        fprintf(fp_gmt,
                "pstext $JVAL $RVAL $BVAL -N -K -O << END >> %s.ps\n%lf %lf %d %d %d %d %s Residuals      %s\nEND\n\n",
                fn_ps_output, xres - 0.25, yres,
                HYPO_FONT_SIZE, 0, STA_FONT, 7, arg_elements[0], res_legend_mag_string);
        // negative residual
        xres += (PLOT_WIDTH / 5.0) / scale;
        genResidualGMT(fp_gmt, "-N", -res_legend_mag, xres + xres_shift, yres, res_scale);
        sprintf(res_legend_mag_string, "%.2f sec", -res_legend_mag);
        fprintf(fp_gmt,
                "pstext $JVAL $RVAL $BVAL -N -K -O << END >> %s.ps\n%lf %lf %d %d %d %d %s\nEND\n\n",
                fn_ps_output, xres - 0.25, yres,
                HYPO_FONT_SIZE, 0, STA_FONT, 7, res_legend_mag_string);

    }



    /* end plot */

    sprintf(signature_str, "%s   %s:v%s %s",
            Hypo.signature, PNAME, PVER, CurrTimeStr());

    fprintf(fp_gmt,
            "pstext -R0.5/1.0/0.5/1.0 -Bf10 -JX%lf/%lf -X%lf -Y%lf -O << END >> %s.ps\n%lf %lf %d %d %d %d %s\nEND\n\n",
            PLOT_WIDTH, PLOT_HEIGHT, -xshift_cum, -yshift_cum,
            fn_ps_output,
            0.501, 0.501, ANNOTATION_FONT_SIZE, 0, ANNOTATION_FONT,
            1, signature_str);

    char *env_name;
    env_name = getenv("PS_VIEWER");
    if (env_name != NULL) {
        fprintf(fp_gmt, "echo \"Finished!  View output:\"\n");
        fprintf(fp_gmt, "echo \"   ${PS_VIEWER} %s.ps\"\n", fn_ps_output);
    } else {
        fprintf(fp_gmt, "echo \"Finished!  Output is:\"\n");
        fprintf(fp_gmt, "echo \"   %s.ps\"\n", fn_ps_output);
    }



    /* clean up */

    fclose(fp_gmt);

    /* run gmt script */

    fprintf(stdout, "\n\nRunning GMT script %s ...\n", fn_gmt);
    sprintf(sys_string, "chmod a+x %s", fn_gmt);
    system(sys_string);
    system(fn_gmt);


    return (0);

}

/*** function to generate GMT commands for a grid view */

int GenGridViewGMT(GridDesc* pgrid, char cviewmode, char cdatatype,
        char arg_elements[][20], int num_arg_elements,
        int ix1, int iy1, int ix2, int iy2, int izlevel,
        double plot_width, char* chr_bounds, char* chr_title,
        char* str_shift, double *pscale, double *pxlen, double *pylen,
        double horiz_line, double vert_line, double misfit_max, double res_scale, int res_min_num_readings) {

    int istat;
    int nx, ny, nz;
    double horiz_min = 0.0, horiz_max = 0.0, horiz_dgrid = 0.0;
    double vert_min = 0.0, vert_max = 0.0, vert_dgrid = 0.0, yscalefact = 0.0, zscalefact = 0.0;
    char horiz_label[MAXLINE], vert_label[MAXLINE];
    char horiz_label_deg[MAXLINE], vert_label_deg[MAXLINE];
    char fnscat[FILENAME_MAX];
    char file_id;
    int nstep;
    double htick_int, vtick_int;
    double vxmin, vxmax, vymin, vymax, vdummy, vdgridx = 0.0, vdgridy = 0.0, vlmean;
    float value;

    /*	union
            {
                    long ival;
                    float fval;
            }
            byteval;
     */


    double plot_scale, gmt_scale;
    char gmt_JVAL[FILENAME_MAX];

    double cos_theta, sin_theta;
    int hypotenouse, nhypot;

    int doLatLong = 0;

    int iAzAngle = 0, iDipAngle = 0;

    int nsta;
    SourceDesc* psta;

    FILE *fp_gmtgrd = NULL;

    static double contour_int;

    int nevents, nplotted;
    FILE *fpio_tmp;
    Vect3D max_like, expect;
    Mtrx3D cov;
    Ellipsoid3D ellipsoid;
    FocalMech focalMech;

    double magnitude;

    char fn_stations[FILENAME_MAX];
    FILE* fp_stations;

    int nresiduals;
    char fn_nlloc_stat[FILENAME_MAX];

    double stax, stay;

    int plotExpScatter = 0;

    /* set plot expectation, ellipsoid and scatter */
    // AJL 20060829
    //if (cdatatype == 'M' || cdatatype == 'R')
    // AJL 20070103
    //if (cdatatype == 'R')
    //	plotExpScatter = 1;

    /* convert S plot type to E to enable E arguments with S */
    // AJL 20060829
    /*if (cdatatype == 'S') {
            cdatatype = 'E';
            plotExpScatter = 1;
    }*/


    /* set lat long alternate */
    if (cviewmode == 'H')
        doLatLong = 1;

    /* set view type identifiers */

    if (cviewmode == 'H')
        file_id = 'H';
    else if (ix1 == ix2)
        file_id = 'X';
    else if (iy1 == iy2)
        file_id = 'Y';
    else
        file_id = 'V';

    /* open gmt grid file */

    if (cdatatype == 'G') {
        sprintf(fn_gmtgrd, "%s.%c.grd", fn_root_output, file_id);
        if ((fp_gmtgrd = fopen(fn_gmtgrd, "w")) == NULL) {
            nll_puterr("ERROR opening gmt grid output file.");
            return (-1);
        }
    }


    if (pgrid->type == GRID_ANGLE || pgrid->type == GRID_ANGLE_2D) {
        if (cviewmode == 'H')
            iAzAngle = 1;
        if (cviewmode == 'V')
            iDipAngle = 1;
    }


    /* write requested grid values */

    // horizontal view
    if (cviewmode == 'H') {
        if (message_flag > 0)
            fprintf(stdout, "Generating horizontal view (izlevel=%d)...\n", izlevel);
        if (cdatatype == 'G') {
            for (ny = pgrid->numy - 1; ny >= 0; ny--) {
                for (nx = 0; nx < pgrid->numx; nx++) {
                    value = ReadGrid3dValue(fp_grid, nx, ny, izlevel, pgrid);
                    /*
                    if (nx % 20 == 0 && ny % 20 == 0) {
                    byteval.fval = value;
                    fprintf(stderr, "%d %d %d value %lf (%lf,%d)\n", nx, ny, izlevel, value, byteval.fval, byteval.ival);
                    }
                    if (value < 0.0)
                    value = 0.0;
                     */
                    if (iAzAngle || iDipAngle)
                        value = CalcAngleValue(value, iAzAngle, iDipAngle);
                    if (value > MAX_GRID_VALUE)
                        value = SUBSTITUTE_MAX_GRID_VALUE;
                    fwrite(&value, sizeof (float), 1, fp_gmtgrd);
                    if (value > grid_value_max)
                        grid_value_max = value;
                    if (value < grid_value_min)
                        grid_value_min = value;

                }
            }
            fclose(fp_gmtgrd);
            if (message_flag > 0)
                fprintf(stdout, "%d rows and %d columns written\n",
                    pgrid->numy, pgrid->numx);
        } else if (cdatatype == 'S' || plotExpScatter) {
            istat = Scat2GMT(fnroot_input, "XY", 1, fnscat);
            istat = Scat2GMT(fnroot_input, "XY", 0, fnscat);
        }
        horiz_min = pgrid->origx;
        horiz_max = horiz_min + (double) (pgrid->numx - 1) *
                pgrid->dx;
        horiz_dgrid = pgrid->dx;
        vert_min = pgrid->origy;
        vert_max = vert_min + (double) (pgrid->numy - 1) * pgrid->dy;
        vert_dgrid = pgrid->dy;
        sprintf(horiz_label, "%s(%s)", reverse_xy ? "Y": "X", lengthUnits);
        sprintf(vert_label, "%s(%s)", reverse_xy ? "X": "Y", lengthUnits);
        strcpy(horiz_label_deg, "Long(deg)");
        strcpy(vert_label_deg, "Lat(deg)");
        yscalefact = 1.0;
        zscalefact = 1.0;
        rect2latlon(proj_index_output, horiz_min, vert_min, &vymin, &vxmin);
        rect2latlon(proj_index_output, horiz_max, vert_max, &vymax, &vxmax);
        vdgridx = (vxmax - vxmin) / (double) (pgrid->numx - 1);
        vdgridy = (vymax - vymin) / (double) (pgrid->numy - 1);
    }// location view NS section
    else if (cviewmode == 'L' && ix1 == ix2) {
        if (message_flag > 0)
            fprintf(stdout, "Generating Y section view...\n");
        if (cdatatype == 'G') {
            for (ny = iy2 - 1; ny >= iy1; ny--) {
                for (nz = 0; nz < pgrid->numz; nz++) {
                    value = ReadGrid3dValue(fp_grid, ix1, ny, nz, pgrid);
                    if (iAzAngle || iDipAngle)
                        value = CalcAngleValue(value, iAzAngle, iDipAngle);
                    if (value > MAX_GRID_VALUE)
                        value = SUBSTITUTE_MAX_GRID_VALUE;
                    fwrite(&value, sizeof (float), 1, fp_gmtgrd);
                    if (value > grid_value_max)
                        grid_value_max = value;
                    if (value < grid_value_min)
                        grid_value_min = value;
                }
            }
            fclose(fp_gmtgrd);
            if (message_flag > 0)
                fprintf(stdout, "%d rows and %d columns written\n",
                    pgrid->numz, iy2 - iy1);
        } else if (cdatatype == 'S' || plotExpScatter) {
            istat = Scat2GMT(fnroot_input, "ZY", 1, fnscat);
            istat = Scat2GMT(fnroot_input, "ZY", 0, fnscat);
        }
        vert_min = pgrid->origy + (double) iy1 * pgrid->dy;
        vert_max = vert_min + (double) (iy2 - iy1 - 1) *
                pgrid->dy;
        vert_dgrid = pgrid->dy;
        horiz_min = pgrid->origz;
        horiz_max = horiz_min + (double) (pgrid->numz - 1) * pgrid->dz;
        horiz_dgrid = pgrid->dz;
        sprintf(horiz_label, "Z(%s)", lengthUnits);
        sprintf(vert_label, "%s(%s)", reverse_xy ? "X": "Y", lengthUnits);
        sprintf(horiz_label_deg, "Z(%s)", lengthUnits);
        strcpy(vert_label_deg, "Lat(deg)");
        yscalefact = 1.0;
        zscalefact = 1.0;
        rect2latlon(proj_index_output, pgrid->origx + (double) ix1 * pgrid->dx,
                (vert_min + vert_max) / 2.0, &vdummy, &vymin);
        if (message_flag > 0)
            fprintf(stdout, "LongitudeXX Section: %lf\n", vymin);
        rect2latlon(proj_index_output, horiz_min, vert_min, &vymin, &vdummy);
        rect2latlon(proj_index_output, horiz_max, vert_max, &vymax, &vdummy);
        vxmin = horiz_min;
        vxmax = horiz_max;
        vdgridx = horiz_dgrid;
        vdgridy = (vymax - vymin) / (double) (pgrid->numy - 1);
    }// NS view
    else if (ix1 == ix2) {
        if (message_flag > 0)
            fprintf(stdout, "Generating Y section view...\n");
        if (cdatatype == 'G') {
            for (nz = pgrid->numz - 1; nz >= 0; nz--)
                /*		for (nz = 0; nz < pgrid->numz; nz++)*/ {
                for (ny = iy1; ny < iy2; ny++) {
                    value = ReadGrid3dValue(fp_grid, ix1, ny, nz, pgrid);
                    if (iAzAngle || iDipAngle)
                        value = CalcAngleValue(value, iAzAngle, iDipAngle);
                    if (value > MAX_GRID_VALUE)
                        value = SUBSTITUTE_MAX_GRID_VALUE;
                    fwrite(&value, sizeof (float), 1, fp_gmtgrd);
                    if (value > grid_value_max)
                        grid_value_max = value;
                    if (value < grid_value_min)
                        grid_value_min = value;
                }
            }
            fclose(fp_gmtgrd);
            if (message_flag > 0)
                fprintf(stdout, "%d rows and %d columns written\n",
                    pgrid->numz, iy2 - iy1);
        } else if (cdatatype == 'S' || plotExpScatter) {
            istat = Scat2GMT(fnroot_input, "YZ", 1, fnscat);
            istat = Scat2GMT(fnroot_input, "YZ", 0, fnscat);
        }
        /*horiz_min = pgrid->origy;*/
        horiz_min = pgrid->origy + (double) iy1 * pgrid->dy;
        /*horiz_max = horiz_min + (double) (pgrid->numy - 1) *
                pgrid->dy;*/
        horiz_max = horiz_min + (double) (iy2 - iy1 - 1) *
                pgrid->dy;
        horiz_dgrid = pgrid->dy;
        vert_min = pgrid->origz;
        vert_max = vert_min + (double) (pgrid->numz - 1) * pgrid->dz;
        vert_dgrid = pgrid->dz;
        sprintf(horiz_label, "%s(%s)", reverse_xy ? "X": "Y", lengthUnits);
        sprintf(vert_label, "Z(%s)", lengthUnits);
        strcpy(horiz_label_deg, "Lat(deg)");
        sprintf(vert_label_deg, "Z(%s)", lengthUnits);
        yscalefact = -1.0;
        zscalefact = 1.0;
        rect2latlon(proj_index_output, pgrid->origx + (double) ix1 * pgrid->dx,
                (horiz_min + horiz_max) / 2.0, &vdummy, &vxmin);
        if (message_flag > 0)
            fprintf(stdout, "LongitudeXX Section: %lf\n", vxmin);
        rect2latlon(proj_index_output, horiz_min, vert_min, &vxmin, &vdummy);
        rect2latlon(proj_index_output, horiz_max, vert_max, &vxmax, &vdummy);
        vymin = vert_min;
        vymax = vert_max;
        vdgridy = horiz_dgrid;
        vdgridx = (vxmax - vxmin) / (double) (pgrid->numx - 1);
    }// EW view
    else if (iy1 == iy2) {
        if (message_flag > 0)
            fprintf(stdout, "Generating X section view...\n");
        if (cdatatype == 'G') {
            for (nz = pgrid->numz - 1; nz >= 0; nz--) {
                for (nx = ix1; nx < ix2; nx++) {
                    value = ReadGrid3dValue(fp_grid, nx, iy1, nz, pgrid);
                    if (iAzAngle || iDipAngle)
                        value = CalcAngleValue(value, iAzAngle, iDipAngle);
                    if (value > MAX_GRID_VALUE)
                        value = SUBSTITUTE_MAX_GRID_VALUE;
                    fwrite(&value, sizeof (float), 1, fp_gmtgrd);
                    if (value > grid_value_max)
                        grid_value_max = value;
                    if (value < grid_value_min)
                        grid_value_min = value;
                }
            }
            fclose(fp_gmtgrd);
            if (message_flag > 0)
                fprintf(stdout, "%d rows and %d columns written\n",
                    pgrid->numz, ix2 - ix1);
        } else if (cdatatype == 'S' || plotExpScatter) {
            istat = Scat2GMT(fnroot_input, "XZ", 1, fnscat);
            istat = Scat2GMT(fnroot_input, "XZ", 0, fnscat);
        }
        horiz_min = pgrid->origx + (double) ix1 * pgrid->dx;
        horiz_max = horiz_min + (double) (ix2 - ix1 - 1) *
                pgrid->dx;
        horiz_dgrid = pgrid->dx;
        vert_min = pgrid->origz;
        vert_max = vert_min + (double) (pgrid->numz - 1) * pgrid->dz;
        vert_dgrid = pgrid->dz;
        sprintf(horiz_label, "%s(%s)", reverse_xy ? "Y": "X", lengthUnits);
        sprintf(vert_label, "Z(%s)", lengthUnits);
        strcpy(horiz_label_deg, "Long(deg)");
        sprintf(vert_label_deg, "Z(%s)", lengthUnits);
        yscalefact = -1.0;
        zscalefact = -1.0;
        rect2latlon(proj_index_output, (horiz_min + horiz_max) / 2.0,
                pgrid->origy + (double) iy1 * pgrid->dy,
                &vlmean, &vdummy);
        if (message_flag > 0)
            fprintf(stdout, "Latitude Section: %lf\n", vlmean);
        rect2latlon(proj_index_output, horiz_min, vert_min, &vdummy, &vxmin);
        rect2latlon(proj_index_output, horiz_max, vert_max, &vdummy, &vxmax);
        vymin = vert_min;
        vymax = vert_max;
        vdgridx = (vxmax - vxmin) / (double) (pgrid->numx - 1);
        vdgridy = vert_dgrid;
    }// oblique vertical view
    else if (cviewmode == 'V') {
        if (message_flag > 0)
            fprintf(stdout, "Generating oblique vertical section view...\n");

        // assume dx = dy
        if (pgrid->dx != pgrid->dy) {
            nll_puterr("ERROR cannot plot oblique vertical section when grid dx != dy.");
            exit(-1);
        }

        hypotenouse = 1 + (int) (sqrt((double) (ix2 - ix1) * (double) (ix2 - ix1)
                + (double) (iy2 - iy1) * (double) (iy2 - iy1)));
        cos_theta = (double) (ix2 - ix1) / (double) (hypotenouse);
        sin_theta = (double) (iy2 - iy1) / (double) (hypotenouse);

        if (cdatatype == 'G') {
            for (nz = pgrid->numz - 1; nz >= 0; nz--) {
                for (nhypot = 0; nhypot < hypotenouse; nhypot++) {
                    nx = ix1 + (int) ((double) nhypot * cos_theta);
                    ny = iy1 + (int) ((double) nhypot * sin_theta);

                    if (nx >= pgrid->numx) nx = pgrid->numx - 1;
                    if (ny >= pgrid->numy) ny = pgrid->numy - 1;
                    if (nx < 0) nx = 0;
                    if (ny < 0) ny = 0;

                    value = ReadGrid3dValue(fp_grid, nx, ny, nz, pgrid);
                    if (iAzAngle || iDipAngle)
                        value = CalcAngleValue(value, iAzAngle, iDipAngle);
                    if (value > MAX_GRID_VALUE)
                        value = SUBSTITUTE_MAX_GRID_VALUE;
                    fwrite(&value, sizeof (float), 1, fp_gmtgrd);
                    if (value > grid_value_max)
                        grid_value_max = value;
                    if (value < grid_value_min)
                        grid_value_min = value;
                }
            }
            fclose(fp_gmtgrd);
            if (message_flag > 0)
                fprintf(stdout, "%d rows and %d columns written\n",
                    pgrid->numz, ix2 - ix1);
        } else if (cdatatype == 'S') {
            nll_puterr("ERROR cannot have datatype = S with oblique Vertical section.");
            exit(-1);
        } else if (cdatatype == 'E') {
            nll_puterr("ERROR cannot have datatype = E with oblique Vertical section.");
            exit(-1);
        } else if (cdatatype == 'M') {
            nll_puterr("ERROR cannot have datatype = M with oblique Vertical section.");
            exit(-1);
        } else if (cdatatype == 'R') {
            nll_puterr("ERROR cannot have datatype = R with oblique Vertical section.");
            exit(-1);
        }
        horiz_dgrid = pgrid->dx; // assume dx = dy
        horiz_min = 0.0;
        horiz_max = horiz_min + (double) (hypotenouse - 1) * horiz_dgrid;
        vert_min = pgrid->origz;
        vert_max = vert_min + (double) (pgrid->numz - 1) * pgrid->dz;
        vert_dgrid = pgrid->dz;
        sprintf(horiz_label, "Dist(%s)", lengthUnits);
        sprintf(vert_label, "Z(%s)", lengthUnits);
        strcpy(horiz_label_deg, "---");
        sprintf(vert_label_deg, "Z(%s)", lengthUnits);
        yscalefact = -1.0;
        zscalefact = -1.0;
        /*		rect2latlon(proj_index_output, (horiz_min + horiz_max) / 2.0,
                                pgrid->origy + (double) iy1 * pgrid->dy,
                                &vlmean, &vdummy);
                        if (message_flag > 0)
                                fprintf(stdout, "Latitude Section: %lf\n", vlmean);
                        rect2latlon(proj_index_output, horiz_min, vert_min, &vdummy, &vxmin);
                        rect2latlon(proj_index_output, horiz_max, vert_max, &vdummy, &vxmax);
         */
        vymin = vert_min;
        vymax = vert_max;
        // 20140102 AJL - bug fix
        //vdgridx = (vxmax - vxmin) / (double) (pgrid->numx - 1);
        vdgridx = (horiz_max - horiz_min) / (double) (pgrid->numx - 1);
        vdgridy = vert_dgrid;
    }



    if (message_flag > 0)
        fprintf(stdout, "Grid value min = %e,  value max = %e\n",
            grid_value_min, grid_value_max);


    /* try to open and read station list file */

    sprintf(fn_stations, "%s.stations", fnroot_input);
    if ((fp_stations = fopen(fn_stations, "r")) != NULL) {
        NumStations = ReadStationList(fp_stations, Stations, 0);
        fclose(fp_stations);
    } else {
        if (message_flag >= 1)
            nll_putmsg2(1, "INFO: cannot open station list file", fn_stations);
        // will use stations found in hypocenter-phase files
    }



    /* write gmt script file */


    /* GMT RVAL, BVAL, JVAL */

    /* RVAL rectangular x/y */
    fprintf(fp_gmt, "# Rect x/y in km\n");
    fprintf(fp_gmt, "set RVAL = \'-R%lf/%lf/%lf/%lf/-9999/9999\'\n",
            horiz_min, horiz_max, vert_min, vert_max);
    if (message_flag > 0)
        fprintf(stdout, "RVAL: RECT/X: %lf/%lf RECTY:%lf/%lf\n",
            horiz_min, horiz_max, vert_min, vert_max);
    /* BVAL rectangular x/y */
    printf("GetContourInterval htick_int R: ");
    htick_int = GetContourInterval(horiz_min, horiz_max, 3, &nstep);
    printf("GetContourInterval vtick_int R: ");
    vtick_int = GetContourInterval(vert_min, vert_max, 3, &nstep);
    fprintf(fp_gmt, "# Rect x/y in km\n");
    fprintf(fp_gmt, "set BVAL = \'-B%lf:%s:/%lf:%s::.%s:%s\'\n",
            htick_int, horiz_label, vtick_int, vert_label,
            chr_title, chr_bounds);
    /* JVAL rectangular x/y */
    // set scale based on plot width
    if (*pscale < 0.0) {
        *pscale = plot_width / (horiz_max - horiz_min);
        // check that plot is not too high
        if ((vert_max - vert_min) * *pscale > 0.7 * PLOT_HEIGHT)
            *pscale = (0.7 * PLOT_HEIGHT) / (vert_max - vert_min);
    }
    plot_scale = *pscale;

    *pxlen = (horiz_max - horiz_min) * plot_scale;
    *pylen = (vert_max - vert_min) * plot_scale;

    fprintf(fp_gmt, "# Rect x/y in km\n");
    fprintf(fp_gmt, "set JVAL = \'-Jx%lf/%lf -Jz%lf\'\n", plot_scale,
            plot_scale * yscalefact, plot_scale * zscalefact);
    fprintf(fp_gmt, "\n");

    if (doLatLong)
        fprintf(fp_gmt, "if (! $PLOT_LAT_LONG) then\n");
    fprintf(fp_gmt,
            "psbasemap $JVAL $RVAL $BVAL %s -K -O >> %s.ps\n", str_shift, fn_ps_output);
    if (doLatLong)
        fprintf(fp_gmt, "endif\n");


    if (doLatLong) {
        fprintf(fp_gmt, "if ($PLOT_LAT_LONG) then\n");
        /* RVAL geographic version */
        fprintf(fp_gmt, "# Latitude/Longitude in degrees\n");
        sprintf(gmt_RVAL_latlong_string, "\'-R%lf/%lf/%lf/%lf\'",
                vxmin, vxmax, vymin, vymax);
        fprintf(fp_gmt, "set RVAL = %s\n", gmt_RVAL_latlong_string);
        if (message_flag > 0)
            fprintf(stdout, "   => LONG/X: %lf/%lf  LAT/Y:%lf/%lf\n",
                vxmin, vxmax, vymin, vymax);
        /* BVAL geographic version */
        printf("GetContourInterval htick_int G: ");
        htick_int = GetContourInterval(vxmin, vxmax, 3, &nstep);
        printf("GetContourInterval vtick_int G: ");
        vtick_int = GetContourInterval(vymin, vymax, 3, &nstep);
        fprintf(fp_gmt, "# Latitude/Longitude in degrees\n");
        fprintf(fp_gmt, "set BVAL = \'-B%lf:%s:/%lf:%s::.%s:%s\'\n",
                htick_int, horiz_label_deg, vtick_int, vert_label_deg,
                chr_title, chr_bounds);
        /* JVAL geographic version */
        gmt_scale = getGMTJVAL(proj_index_output, gmt_JVAL_latlong_string, *pxlen, vxmax, vxmin, *pylen, vymax, vymin);
        fprintf(fp_gmt, "# Latitude/Longitude in degrees\n");
        sprintf(gmt_JVAL, "set JVAL = \'%s -Jz%lf\'",
                gmt_JVAL_latlong_string, gmt_scale * zscalefact);
        fprintf(fp_gmt, "%s\n", gmt_JVAL);

        fprintf(fp_gmt,
                "psbasemap ${JVAL} ${RVAL} ${BVAL} %s -K -O >> %s.ps\n", str_shift, fn_ps_output);

        fprintf(fp_gmt, "endif\n");
    }

    fprintf(fp_gmt, "\n");




    /* plot GMT_GRID geographic features */

    if (cviewmode == 'H') {
        MapFiles2GMT(horiz_min, vert_min, horiz_max, vert_max,
                fp_gmt, fn_ps_output, doLatLong, 1);
    }




    //printf("cdatatype = <%c>\n", cdatatype);

    if (cdatatype == 'G') {

        fprintf(fp_gmt, "# Rect x/y in km\n");
        if (doLatLong)
            fprintf(fp_gmt, "if (! $PLOT_LAT_LONG) then\n");
        if (GMT_VER_3_3_4) {
            fprintf(fp_gmt,
                    "xyz2grd %s -G%sgmt -I%lf/%lf $RVAL -Dkm/km/=/0.0/0.0/%s/remark -V -Zf\n",
                    fn_gmtgrd, fn_gmtgrd, horiz_dgrid, vert_dgrid, fn_root_output);
        } else {
            fprintf(fp_gmt,
                    "xyz2grd %s -G%sgmt -I%lf/%lf $RVAL -Dkm/km/=/0.0/0.0/%s/remark -V -Z -b\n",
                    fn_gmtgrd, fn_gmtgrd, horiz_dgrid, vert_dgrid, fn_root_output);
        }
        if (doLatLong)
            fprintf(fp_gmt, "endif\n\n");

        /* commented out line giving geographic version of xyz2grd */
        if (doLatLong) {
            fprintf(fp_gmt, "# Latitude/Longitude in degrees\n");
            fprintf(fp_gmt, "if ($PLOT_LAT_LONG) then\n");
            if (GMT_VER_3_3_4) {
                fprintf(fp_gmt,
                        "xyz2grd %s -G%sgmt -I%lf/%lf ${RVAL} -Ddeg/deg/=/0.0/0.0/%s/remark -V -Zf\n",
                        fn_gmtgrd, fn_gmtgrd, vdgridx, vdgridy, fn_root_output);
            } else {
                fprintf(fp_gmt,
                        "xyz2grd %s -G%sgmt -I%lf/%lf ${RVAL} -Ddeg/deg/=/0.0/0.0/%s/remark -V -Z -b\n",
                        fn_gmtgrd, fn_gmtgrd, vdgridx, vdgridy, fn_root_output);
            }
            fprintf(fp_gmt, "endif\n\n");
        }

        fprintf(fp_gmt, "set SCALE_FLAG = \n\n");

        if (pgrid->type == GRID_PROB_DENSITY) {

            sprintf(fn_cont, "%s.conf", fnroot_input);
            MakeConfCPT(fn_cont, fn_root_output);
            fprintf(fp_gmt,
                    "grdimage -S-n %sgmt -C%s.conf.cpt $JVAL $RVAL $BVAL -K -O >> %s.ps\n\n",
                    fn_gmtgrd, fn_root_output, fn_ps_output);
            fprintf(fp_gmt,
                    "grdcontour %sgmt %s -C%s.conf.cpt $JVAL $RVAL $BVAL -K -O >> %s.ps\n\n",
                    fn_gmtgrd, label_contours, fn_root_output, fn_ps_output);

        } else if (pgrid->type == GRID_MISFIT) {

            if (iFirstPlot) {
                if (misfit_max > 4.999 * MISFIT_CONTOUR_INTERVAL_XY &&
                        misfit_max < 50.0 * MISFIT_CONTOUR_INTERVAL_XY)
                    contour_int = MISFIT_CONTOUR_INTERVAL_XY;
                else {
                    printf("GetContourInterval contour_int GRID_MISFIT: ");
                    contour_int = GetContourInterval(0.0, (double) grid_value_max, 11, &nstep);
                }
                MakeMisfitCPT(fn_root_output, misfit_max, contour_int);
                iFirstPlot = 0;
            }
            fprintf(fp_gmt,
                    "grdimage -S-n %sgmt -C%s.misfit.cpt $JVAL $RVAL $BVAL -K -O >> %s.ps\n\n",
                    fn_gmtgrd, fn_root_output, fn_ps_output);

            fprintf(fp_gmt,
                    "grdcontour %sgmt %s -C%s.misfit.cpt $JVAL $RVAL $BVAL -K -O >> %s.ps\n\n",
                    fn_gmtgrd, label_contours, fn_root_output, fn_ps_output);

        } else {

            if (iFirstPlot) {
                printf("GetContourInterval contour_int iFirstPlot: ");
                contour_int = GetContourInterval((double) grid_value_min,
                        (double) grid_value_max, NUM_COLORS, &nstep);
                fprintf(fp_gmt,
                        "if (-e Grid2GMT.cpt) then\n");
                fprintf(fp_gmt,
                        "   unalias cp\n");
                fprintf(fp_gmt,
                        "   cp  Grid2GMT.cpt %s.cpt\n", fn_root_output);
                fprintf(fp_gmt,
                        "   set SCALE_FLAG = \n");
                fprintf(fp_gmt,
                        "else\n");
                if (GMT_VER_3_3_4) {
                    if (pgrid->type == GRID_LIKELIHOOD) { // likelihood color table
                        fprintf(fp_gmt,
                                "   makecpt -Z -Chot -I -T0/1/0.1 > %s.cpt\n",
                                fn_root_output);
                    } else { // rainbow
                        double contour_int_cpt = contour_int;
                        double value_min = contour_int * floor(grid_value_min / contour_int);
                        double value_max = contour_int * (1.0 + ceil(grid_value_max / contour_int));
                        if (value_max - value_min < value_min / 1000.0) {
                            value_min -= value_min / 100.0;
                            value_max += value_max / 100.0;
                            contour_int_cpt = (value_max - value_min) / 3.0;
                        }
                        char cpt_colortable[MAXLINE];
                        strcpy(cpt_colortable, "rainbow");
                        /*if (value_min < 0.0 && value_max > 0.0) {
                            // value range straddles zero, set min/max for color tablel equal.
                            strcpy(cpt_colortable, "seis");
                            if (value_max < contour_int)
                                value_max = contour_int;
                            else if (value_min > -contour_int)
                                value_min = -contour_int;
                        } else {
                            // value range positive
                            strcpy(cpt_colortable, "rainbow");
                        }*/
                        char cpt_command[10 * MAXLINE];
                        sprintf(cpt_command, "makecpt -Z -C%s -T%g/%g/%g > %s.cpt",
                                cpt_colortable, value_min, value_max, contour_int_cpt, fn_root_output);
                        fprintf(fp_gmt, "   %s\n", cpt_command);
                        //if (message_flag > 0)
                        nll_putmsg2(1, "INFO:", cpt_command);
                    }
                } else {
                    fprintf(fp_gmt,
                            "   makecpt -C%.1le -S%dc -m%f > %s.cpt\n",
                            contour_int, nstep + 2, contour_int * (double) ((int)
                            ((grid_value_max + grid_value_min)
                            / (2.0 * contour_int))), fn_root_output);
                }
                fprintf(fp_gmt,
                        "   set SCALE_FLAG = -L\n");
                fprintf(fp_gmt,
                        "endif\n\n");
                iFirstPlot = 0;
            }
            fprintf(fp_gmt, "grdimage -S-n %sgmt -C%s.cpt $JVAL $RVAL $BVAL -K -O >> %s.ps\n\n",
                    fn_gmtgrd, fn_root_output, fn_ps_output);

            fprintf(fp_gmt,
                    "grdcontour %sgmt %s -C%s.cpt $JVAL $RVAL $BVAL -K -O >> %s.ps\n\n",
                    fn_gmtgrd, label_contours, fn_root_output, fn_ps_output);

        }

        NumStations += addStations(Stations, NumStations, Arrival, NumArrivals);

    }

    if (cdatatype == 'S' || plotExpScatter) {

        /* draw scatter locations */

        fprintf(fp_gmt, "psxy %s $JVAL $RVAL -W1/255/0/0 -Sp -m -K -O >> %s.ps\n",
                fnscat, fn_ps_output);

        NumStations +=
                addStations(Stations, NumStations, Arrival, NumArrivals);

    }

    if (cdatatype == 'E' || cdatatype == 'S' || plotExpScatter) {

        /* draw "traditional" expectation and error ellipses */

        nevents = nplotted = 0;
        fpio_tmp = NULL;
        while (ReadHypStatistics(&fpio_tmp, fnroot_input,
                &max_like, &expect, &cov, &ellipsoid,
                Arrival, &NumArrivals) != EOF) {
            nevents++;
            NumStations +=
                    addStations(Stations, NumStations, Arrival, NumArrivals);

            if ((istat = PlotTraditionStats(cdatatype, file_id,
                    PLOT_WIDTH / 25.0, &max_like, &expect,
                    &ellipsoid, fp_gmt, arg_elements[0], "0/0/255")) < 0)
                break;
            nplotted += istat;
        }
        if (cdatatype == 'E' && num_arg_elements < 1 && message_flag > 0)
            nll_puterr(
                "ERROR: 'E' datatype string does not have all required arguments.");
        if (message_flag > 0)
            fprintf(stdout,
                "PlotTraditionStats: %d events read, %d plotted\n", nevents, nplotted);


    }


    /* draw lines through maximum likelihood location */

    if (!(cdatatype == 'S' || cdatatype == 'E' || plotExpScatter)
            && (pgrid->type == GRID_PROB_DENSITY
            || pgrid->type == GRID_MISFIT)) {
        fprintf(fp_gmt,
                "psxy $JVAL $RVAL -W1/0/0/0/dotted -m -K -O << END >> %s.ps\n",
                fn_ps_output);
        fprintf(fp_gmt,
                ">\n%lf %lf\n%lf %lf\n>\n%lf %lf\n%lf %lf\nEND\n\n",
                horiz_min, horiz_line, horiz_max, horiz_line,
                vert_line, vert_min, vert_line, vert_max);
    }


    /* draw "traditional" expectation and error ellipses */

    if (!(cdatatype == 'S' || cdatatype == 'E' || plotExpScatter)
            && pgrid->type == GRID_PROB_DENSITY) {
        PlotTraditionStats(cdatatype, file_id,
                PLOT_WIDTH / 25.0, &MaxLike, &Expectation,
                &Ellipsoid, fp_gmt, "111", "0/0/0");

    }


    if (cdatatype == 'M') {

        /* draw focal mechanisms */

        nevents = nplotted = 0;
        fpio_tmp = NULL;
        while (ReadFocalMech(&fpio_tmp, fnroot_input, &focalMech, Arrival, &NumArrivals) != EOF) {

            nevents++;

            // check if null focal mechanism
            if (focalMech.nObs < 1)
                continue;
            NumStations += addStations(Stations, NumStations, Arrival, NumArrivals);
            magnitude = 2.0;
            if ((istat = PlotFocalMechanism(file_id,
                    PLOT_WIDTH / 10.0, magnitude, &focalMech,
                    fp_gmt, arg_elements[0], "-W1/255/0/0 ", "-G255/0/0 ",
                    horiz_min, horiz_max, vert_min, vert_max, plot_scale)) < 0)
                break;
            nplotted += istat;
        }
        if (message_flag > 0)
            fprintf(stdout,
                "PlotFocalMechanisms: %d events read, %d plotted\n", nevents, nplotted);


    }

    if (cdatatype == 'R') {

        /* draw residuals */

        if (cviewmode == 'H') {
            sprintf(fn_nlloc_stat, "%s.stat", fnroot_input);
            nresiduals = ConvertResids2MapGMT(
                    fn_nlloc_stat, arg_elements[0], fp_gmt, Stations, NumStations, res_scale, res_min_num_readings);
            if (num_arg_elements < 1 && message_flag > 0) {
                nll_puterr(
                        "ERROR: 'R' datatype string does not have all required arguments.");
                exit(-1);
            }
            if (message_flag > 0)
                fprintf(stdout,
                    "PlotResiduals: %d plotted\n", nresiduals);


        }

    }


    /* plot stations */

    if (/*cdatatype != 'E' &&*/ cviewmode == 'H') {

        for (nsta = 0; nsta < NumStations; nsta++) {
            psta = Stations + nsta;
            convertCoordsRect(proj_index_input, proj_index_output,
                    psta->x, psta->y, &stax, &stay);

            fprintf(fp_gmt, "# Station\npstext $JVAL $RVAL -S4,0 -G255 -K -O << END >> %s.ps\n%lf %lf %d %d %d %d %s\nEND\n\n", fn_ps_output,
                    stax, stay, STA_FONT_SIZE, 0, STA_FONT, 6, psta->label);
        }
    }

    /* plot line geographic features */

    if (cviewmode == 'H') {
        MapFiles2GMT(horiz_min, vert_min, horiz_max, vert_max,
                fp_gmt, fn_ps_output, doLatLong, 0);
    }



    /* run auxilliary gmt script */

    fprintf(fp_gmt, "if (-e Grid2GMT.%c.gmt) then\n", file_id);
    fprintf(fp_gmt, "   echo  'Running auxilliary GMT script: Grid2GMT.%c.gmt'\n",
            file_id);
    fprintf(fp_gmt, "   source  Grid2GMT.%c.gmt\n", file_id);
    fprintf(fp_gmt, "endif\n\n");


    /* redraw axes */

    fprintf(fp_gmt,
            "psbasemap $JVAL $RVAL $BVAL -K -O >> %s.ps\n", fn_ps_output);
    fprintf(fp_gmt, "\n");
    fprintf(fp_gmt, "\n");


    return (0);

}

/*** function to add station arrival to station list */

int addStations(SourceDesc *stations, int numStations, ArrivalDesc *arrival, int nArrivals) {

    int i, n, nAdded = 0;
    ;


    for (i = 0; i < nArrivals; i++) {

        if (numStations >= MAX_NUM_ARRIVALS)
            return (0);

        for (n = 0; n < numStations; n++) {
            if (strcmp((stations + n)->label, (arrival + i)->label) == 0)
                break;
        }

        if (n == numStations) {
            *(stations + n) = (arrival + i)->station;
            strcpy((stations + n)->label, (arrival + i)->label);
            nAdded++;
            numStations++;
        }

    }

    return (nAdded);


}

/*** function to read station list from file */

int ReadStationList(FILE* fpio, SourceDesc *stations, int numStations) {

    int n;

    n = numStations;
    while (fscanf(fpio, "%s %lf %lf %lf",
            (stations + n)->label, &((stations + n)->x),
            &((stations + n)->y), &((stations + n)->z)) != EOF) {
        n++;
    }

    return (n);
}

/*** function to decode approprate angle value */

double CalcAngleValue(double avalue, int iAzAngle, int iDipAngle) {

    int iqual;
    double azim, dip;
    double value = 0.0;
    TakeOffAngles angles;
    SetAnglesFloat(&angles, avalue);
    iqual = GetTakeOffAngles(&angles, &azim, &dip, &iqual);
    if (iDipAngle) {
        if (iqual < 5)
            value = 181.0;
        else
            value = dip;
    } else if (iAzAngle) {
        if (iqual < 5)
            value = -1.0;
        else
            value = azim;
    }

    return (value);
}

/*** function to select a tick or contour interval */

double GetContourInterval(double value_min, double value_max, int nstep_min, int* pnstep) {
    double contour_int = 1.0e10;


    if (fabs(value_min - value_max) < SMALL_DOUBLE) {
        contour_int = 1.0;
        *pnstep = 3;
    } else {

        while (contour_int > (value_max - value_min) / (double) nstep_min
                && contour_int > SMALL_DOUBLE) {
            contour_int /= 2.0;
            if (contour_int < (value_max - value_min) / (float) nstep_min)
                break;
            contour_int /= 2.5;
            if (contour_int < (value_max - value_min) / (float) nstep_min)
                break;
            contour_int /= 2.0;
        }

        *pnstep = 1 + (int) ((value_max - value_min) / contour_int);
    }

    if (message_flag > 0)
        fprintf(stdout,
            "GetContourInterval: vmin %.2le vmax %.2le  contour_int %.2le  nstep %d \n",
            value_min, value_max, contour_int, *pnstep);

    return (contour_int);

}

/*** function to make a GMT cpt file from confidence contour file */

int MakeConfCPT(char* fn_in, char* fileout) {
    int ndx, nconf, n;
    char fn_out[FILENAME_MAX];
    double contour, conf_level, conf_level_last, cont_low, cont_high;
    FILE *fp_in, *fp_out, *fp_outlevel;

    char conf_lev_str[20][MAXLINE];

    int rgb[11][3] = {
        {255, 0, 0},
        {255, 127, 0},
        {255, 191, 0},
        {255, 223, 0},
        {255, 255, 0},
        {191, 255, 0},
        {0, 255, 0},
        {0, 255, 127},
        {127, 255, 191},
        {254, 254, 254},
        {254, 254, 254}
    };


    if ((fp_in = fopen(fn_in, "r")) == NULL) {
        nll_puterr2("ERROR opening confidence level file", fn_in);
        return (-1);
    }
    sprintf(fn_out, "%s.conf.cpt", fileout);
    if ((fp_out = fopen(fn_out, "w")) == NULL) {
        nll_puterr2("ERROR opening confidence level cpt file", fn_out);
        return (-1);
    }
    sprintf(fn_out, "%s.conflev.cpt", fileout);
    if ((fp_outlevel = fopen(fn_out, "w")) == NULL) {
        nll_puterr2("ERROR opening confidence level cpt file", fn_out);
        return (-1);
    }


    cont_low = -LARGE_DOUBLE;
    conf_level_last = 1.0;
    nconf = 0;
    while (fscanf(fp_in, "%lf C %lf", &contour, &conf_level) != EOF) {
        cont_high = contour;
        ndx = (int) (0.1 + 10.0 * conf_level);
        cont_high = contour;
        fprintf(fp_out, "%le %d %d %d %le %d %d %d\n",
                cont_low, rgb[ndx][0], rgb[ndx][1], rgb[ndx][2],
                cont_high, rgb[ndx][0], rgb[ndx][1], rgb[ndx][2]);
        sprintf(conf_lev_str[nconf], "%le %d %d %d %le %d %d %d\n",
                conf_level, rgb[ndx][0], rgb[ndx][1], rgb[ndx][2],
                conf_level_last, rgb[ndx][0], rgb[ndx][1], rgb[ndx][2]);
        cont_low = cont_high;
        conf_level_last = conf_level;
        nconf++;
    }
    cont_high = LARGE_DOUBLE;
    fprintf(fp_out, "%le %d %d %d %le %d %d %d\n",
            cont_low, rgb[0][0], rgb[0][1], rgb[0][2],
            cont_high, rgb[0][0], rgb[0][1], rgb[0][2]);
    if (conf_level_last > SMALL_DOUBLE)
        sprintf(conf_lev_str[nconf], "%le %d %d %d %le %d %d %d\n",
            0.0, rgb[0][0], rgb[0][1], rgb[0][2],
            conf_level_last, rgb[0][0], rgb[0][1], rgb[0][2]);
    else
        nconf--;

    for (n = nconf; n > 0; n--)
        fprintf(fp_outlevel, "%s", conf_lev_str[n]);


    fclose(fp_in);
    fclose(fp_out);
    fclose(fp_outlevel);

    return (0);


}

/*** function to make a GMT cpt file from confidence contour file */

int MakeMisfitCPT(char* fileout, double value_max, double contour_interval) {
    int ndx;
    char fn_out[FILENAME_MAX];
    double dcont, cont_low, cont_high;
    FILE *fp_out;

    int rgb[11][3] = {
        {255, 0, 0},
        {255, 127, 0},
        {255, 191, 0},
        {255, 223, 0},
        {255, 255, 0},
        {191, 255, 0},
        {0, 255, 0},
        {0, 255, 127},
        {127, 255, 191},
        {254, 254, 254},
        {254, 254, 254}
    };


    sprintf(fn_out, "%s.misfit.cpt", fileout);
    if ((fp_out = fopen(fn_out, "w")) == NULL) {
        nll_puterr("ERROR opening misfit cpt file.");
        return (-1);
    }

    dcont = contour_interval;
    cont_low = 0.0;
    for (ndx = 0; ndx < 10; ndx++) {
        //if (ndx < 9)
        cont_high = cont_low + dcont;
        //else
        //	cont_high = (int) (value_max + contour_interval);
        fprintf(fp_out, "%le %d %d %d %le %d %d %d\n",
                cont_low, rgb[ndx][0], rgb[ndx][1], rgb[ndx][2],
                cont_high, rgb[ndx][0], rgb[ndx][1], rgb[ndx][2]);
        cont_low = cont_high;
    }

    fclose(fp_out);

    return (0);

}

int MakeTopoCPT(char* fileout) {
    int ndx, nstep;
    FILE *fp_out;
    int irgb_min, irgb_max;
    int irgb_low, irgb_high;
    double drgb;
    double cont_low, cont_high, dcont;

    if ((fp_out = fopen(fileout, "w")) == NULL) {
        nll_puterr("ERROR opening topo cpt file.");
        return (-1);
    }

    nstep = 10;
    cont_low = 0.0;
    cont_high = 1.0;
    dcont = (cont_high - cont_low) / (double) (nstep);
    irgb_min = 96;
    irgb_max = 255;
    drgb = (double) (irgb_max - irgb_min) / (double) (nstep);
    for (ndx = 0; ndx < nstep; ndx++) {
        cont_high = cont_low + dcont;
        irgb_high = irgb_min + (int) ((double) ndx + 1) * drgb;
        if (irgb_high > irgb_max)
            irgb_high = irgb_max;
        irgb_low = irgb_min + (int) ((double) ndx) * drgb;
        if (irgb_low > irgb_max)
            irgb_low = irgb_max;
        fprintf(fp_out, "%le %d %d %d %le %d %d %d\n",
                cont_low, 15 * irgb_low / 16, irgb_low, irgb_low,
                cont_high, 15 * irgb_high / 16, irgb_high, irgb_high);
        cont_low = cont_high;
    }

    fclose(fp_out);

    return (0);

}

/*** function to read input file */

int ReadGrid2GMT_Input(FILE* fp_input) {
    int istat, iscan;
    char param[MAXLINE];
    char line[MAXLINE];

    int flag_control = 0, flag_trans = 0, flag_grid = 0;


    mapGridRead = 0;

    /* read each input line */

    while (fgets(line, MAXLINE, fp_input) != NULL) {

        istat = -1;

        /*read parameter line */

        if ((iscan = sscanf(line, "%s", param)) < 0)
            continue;

        /* skip comment line or white space */

        if (strncmp(param, "#", 1) == 0 || isspace(param[0]))
            istat = 0;


        /* read control params */

        if (strncmp(param, "CONTROL", 6) == 0) {
            if ((istat = get_control(strchr(line, ' '))) < 0)
                nll_puterr("Error reading control params.");
            else
                flag_control = 1;
        }


        /* check for graphics input */

        istat = ReadGraphicsInput(fp_input, param, line, istat);


        /*read transform params */

        // output map trans
        if (strncmp(param, "MAPTRANS", 8) == 0) {
            if ((istat = get_transform(1, strchr(line, ' '))) < 0)
                nll_puterr("ERROR reading map transformation parameters.");
            else {
                flag_trans = 1;
                proj_index_output = 1;
            }
        }
        // input location trans
        if (strncmp(param, "TRANS", 5) == 0) {
            if ((istat = get_transform(0, strchr(line, ' '))) < 0)
                nll_puterr("ERROR reading transformation parameters.");
            else {
                flag_trans = 1;
                proj_index_input = 0;
            }
        }


        /* read grid params */

        if (strcmp(param, "MAPGRID") == 0) {
            if ((istat = get_grid(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading grid parameters.");
            else {
                flag_grid = 1;
                mapGridRead = 1;
            }
        }



        /* unrecognized input */

        if (istat < 0 && message_flag > 3) {
            fprintf(stdout, "Skipping input: %s", line);
        }

    }


    /* check for missing input */

    if (!flag_control)
        nll_puterr("ERROR no control (CONTROL) params read.");
    if (!flag_trans)
        nll_puterr("ERROR no transformation (TRANS) params read.");


    return (flag_control + flag_trans - 1);
}

/*** function to draw map file data */

int MapFiles2GMT(double xmin0, double ymin0, double xmax0, double ymax0,
        FILE* fp_gmt, char* fn_ps_output, int doLatLong, int plotGMT_GRD) {
    int nmapfile;

    for (nmapfile = 0; nmapfile < gr_num_map_files; nmapfile++) {
        if (strcmp(mapfile[nmapfile].format, "GMT_GRD") == 0) {
            /* draw grd file */
            if (plotGMT_GRD)
                grd2GMT(nmapfile, xmin0, ymin0, xmax0, ymax0, fp_gmt, fn_ps_output, doLatLong);
        } else {
            /* draw map lines */
            if (!plotGMT_GRD)
                MapLines2GMT(nmapfile, xmin0, ymin0, xmax0, ymax0, fp_gmt, fn_ps_output, doLatLong);
        }
    }


    return (0);

}

/*** function to generate grd command */

int grd2GMT(int nmapfile, double xmin0, double ymin0, double xmax0, double ymax0,
        FILE* fp_gmt, char* fn_ps_output, int doLatLong) {

    /*
     grdimage -S-n grdfile -Ccptfile -Jparameters [  -Btickinfo ]  [
         -Edpi ] [ -G[f|b]rgb ] [ -Iintensfile] [ -K ] [ -M ] [ -O ]
         [ -P ] [ -Rwest/east/south/north[r] ] [ -Ssearch_radius ]  [
         -T[s]  ] [ -U[/dx/dy/][label] ] [ -V ] [ -Xx-shift ] [ -Yy-
         shift ] [ -ccopies ]cpt

    makecpt -Cgebco > my_gebco.cpt
     */

    FILE *fp_tmp;
    char fname_cpt[FILENAME_MAX], fname_temp_cpt[FILENAME_MAX];
    char fname_int[FILENAME_MAX], int_string[FILENAME_MAX + 10];


    // check for cpt file with same name as grid file
    sprintf(fname_cpt, "%s.cpt", mapfile[nmapfile].name);
    if ((fp_tmp = fopen(fname_cpt, "r")) == NULL) {
        // make default cpt file
        sprintf(fname_temp_cpt, "%sgrd.temp.cpt", fnoutput);
        MakeTopoCPT(fname_temp_cpt);
        sprintf(fname_cpt, "%sgrd.cpt", fnoutput);
        fprintf(fp_gmt, "grd2cpt %s -C%s -Z > %s\n", mapfile[nmapfile].name, fname_temp_cpt, fname_cpt);
    }
    if (fp_tmp != NULL)
        fclose(fp_tmp);


    // check for intensity file with same name as grid file
    sprintf(fname_int, "%s.int", mapfile[nmapfile].name);
    if ((fp_tmp = fopen(fname_int, "r")) != NULL) {
        sprintf(int_string, " -I%s ", fname_int);
    } else {
        sprintf(int_string, " ");
    }
    if (fp_tmp != NULL)
        fclose(fp_tmp);


    /* write gmt command */

    //	if (doLatLong)
    //		fprintf(fp_gmt, "if (! $PLOT_LAT_LONG) then\n");
    fprintf(fp_gmt, "grdimage -S-n %s %s -C%s %s %s %s -K -O >> %s.ps\n",
            mapfile[nmapfile].name, int_string, fname_cpt, int_string, gmt_JVAL_latlong_string,
            gmt_RVAL_latlong_string, fn_ps_output);
    //	if (doLatLong)
    //		fprintf(fp_gmt, "endif\n");
    /*
            if (doLatLong) {
                    fprintf(fp_gmt, "if ($PLOT_LAT_LONG) then\n");
                    fprintf(fp_gmt,
    "psxy %s ${JVAL} ${RVAL} -W2/%d/%d/%d -m -K -O >> %s.ps\n",
                            fn_gmtlatlon, ired, igreen, iblue,
                            //texture,
                            fn_ps_output);
                    fprintf(fp_gmt, "endif\n");
            }
     */

    fprintf(fp_gmt, "\n");

    return (0);

}

/*** function to draw map lines */

int MapLines2GMT(int nmapfile, double xmin0, double ymin0, double xmax0, double ymax0,
        FILE* fp_gmt, char* fn_ps_output, int doLatLong) {
    int inside;
    int ired, igreen, iblue;
    char *lstat;
    char line[MAXLINE], texture[MAXLINE],
            fn_gmtxy[FILENAME_MAX], fn_gmtxz[FILENAME_MAX],
            fn_gmtzy[FILENAME_MAX], fn_gmtlatlon[FILENAME_MAX];
    double depth, maplong, maplat, xtemp, ytemp;
    double xmin, xmax, ymin, ymax;
    FILE *fp_map;
    FILE *fp_gmtxy, *fp_gmtxz, *fp_gmtzy, *fp_gmtlatlon;


    /* enlarge bounds sliightly */

    xmin = xmin0 - (xmax0 - xmin0) / 10.0;
    xmax = xmax0 + (xmax0 - xmin0) / 10.0;
    ymin = ymin0 - (ymax0 - ymin0) / 10.0;
    ymax = ymax0 + (ymax0 - ymin0) / 10.0;



    /* draw map lines */

    /*set_line_color(mapfile[nmapfile].rgb, ); */
    /*set_line_style(mapfile[nmapfile].line_style); */

    // 20120504 AJL /* open output files (if exists, do not re-generate) */
    /* open output files */

    GenMapFileName(fn_gmtxy, fn_gmtxz, fn_gmtzy, fn_gmtlatlon, fnoutput,
            mapfile[nmapfile].name);

    // 20120504 AJL if ((fp_gmtxy = fopen(fn_gmtxy, "r")) == NULL || (fp_gmtlatlon = fopen(fn_gmtlatlon, "r")) == NULL) {

        if ((fp_gmtxy = fopen(fn_gmtxy, "w")) == NULL)
            nll_puterr2("ERROR: cannot open map xy output file", fn_gmtxy);
        if ((fp_gmtlatlon = fopen(fn_gmtlatlon, "w")) == NULL)
            nll_puterr2("ERROR: cannot open map lat/lon output file", fn_gmtxy);
        if ((fp_gmtxz = fopen(fn_gmtxz, "w")) == NULL)
            nll_puterr2("ERROR: cannot open map xy output file", fn_gmtxy);
        if ((fp_gmtzy = fopen(fn_gmtzy, "w")) == NULL)
            nll_puterr2("ERROR: cannot open map xy output file", fn_gmtxy);


        /* open input file */

        if ((fp_map = fopen(mapfile[nmapfile].name, "r")) == NULL) {
            nll_puterr2("ERROR: cannot open map file",
                    mapfile[nmapfile].name);
            return (-1);
        }


        if (message_flag > 0)
            fprintf(stdout, "Generating map lines for file %s\n",
                mapfile[nmapfile].name);

        /* draw each segment */

        fprintf(fp_gmtxy, "> XY_LONLAT\n");
        fprintf(fp_gmtxz, "> XY_XZ\n");
        fprintf(fp_gmtzy, "> XY_ZY\n");
        fprintf(fp_gmtlatlon, "> GMT_LONLAT\n");
        depth = 0.0;
        do {

            inside = 0;
            for (;;) {

                if ((lstat = fgets(line, MAXLINE, fp_map))
                        == NULL)
                    break;

                if (strcmp(mapfile[nmapfile].format,
                        "XY_LONLAT") == 0) {
                    sscanf(line, "%lf %lf",
                            &maplong, &maplat);
                    if (maplong < -900.0) {
                        lstat = fgets(line,
                                MAXLINE, fp_map);
                        break;
                    }
                } else if (strncmp(mapfile[nmapfile].format,
                        "GMT_LATLONELEV_M", 17) == 0) {
                    if (sscanf(line, "%lf %lf %lf",
                            &maplat, &maplong, &depth) != 3)
                        break;
                    depth *= -0.001; // convert to detph in km
                } else if (strncmp(mapfile[nmapfile].format,
                        "GMT_LONLATELEV_M", 17) == 0) {
                    if (sscanf(line, "%lf %lf %lf",
                            &maplong, &maplat, &depth) != 3)
                        break;
                    depth *= -0.001; // convert to detph in km
                } else if (strncmp(mapfile[nmapfile].format,
                        "GMT_LATLON", 10) == 0) {
                    if (sscanf(line, "%lf %lf",
                            &maplat, &maplong) != 2)
                        break;
                } else if (strncmp(mapfile[nmapfile].format,
                        "GMT_LONLAT", 10) == 0) {
                    if (sscanf(line, "%lf %lf",
                            &maplong, &maplat) != 2)
                        break;
                } else {
                    nll_puterr2(
                            "ERROR: unrecognized map line type",
                            mapfile[nmapfile].format);
                    lstat = NULL;
                    break;
                }

                latlon2rect(proj_index_output, maplat, maplong, &xtemp, &ytemp);

                if (xtemp < xmin || xtemp > xmax ||
                        ytemp < ymin || ytemp > ymax) {
                    if (inside) {
                        fprintf(fp_gmtxy, ">\n");
                        fprintf(fp_gmtxz, ">\n");
                        fprintf(fp_gmtzy, ">\n");
                        fprintf(fp_gmtlatlon, ">\n");
                    }
                    inside = 0;
                } else {
                    fprintf(fp_gmtxy, "%lf %lf\n",
                            xtemp, ytemp);
                    fprintf(fp_gmtxz, "%lf %lf\n",
                            xtemp, depth);
                    fprintf(fp_gmtzy, "%lf %lf\n",
                            depth, ytemp);
                    fprintf(fp_gmtlatlon, "%lf %lf\n",
                            maplong, maplat);
                    inside = 1;
                }
            }
            if (inside) {
                fprintf(fp_gmtxy, ">\n");
                fprintf(fp_gmtxz, ">\n");
                fprintf(fp_gmtzy, ">\n");
                fprintf(fp_gmtlatlon, ">\n");
            }

        } while (lstat != NULL);


        fclose(fp_map);
        fclose(fp_gmtxy);
        fclose(fp_gmtxz);
        fclose(fp_gmtzy);
        fclose(fp_gmtlatlon);

    // 20120504 AJL }

    // ensure that files are closed
    //	if (fp_gmtxy != NULL)
    //		fclose(fp_gmtxy);
    //	if (fp_gmtlatlon != NULL)
    //	fclose(fp_gmtlatlon);


    /* write gmt command */

    ired = (int) (0.5 + 255.0 * mapfile[nmapfile].rgb.r);
    igreen = (int) (0.5 + 255.0 * mapfile[nmapfile].rgb.g);
    iblue = (int) (0.5 + 255.0 * mapfile[nmapfile].rgb.b);
    if (strcmp(mapfile[nmapfile].line_style, "SOLID") == 0)
        strcpy(texture, "solid");
    else
        strcpy(texture, "solid");

    if (doLatLong)
        fprintf(fp_gmt, "if (! $PLOT_LAT_LONG) then\n");
    fprintf(fp_gmt,
            "psxy %s $JVAL $RVAL -W2/%d/%d/%d -m -K -O >> %s.ps\n",
            fn_gmtxy, ired, igreen, iblue,
            /*texture,*/ fn_ps_output);
    if (doLatLong)
        fprintf(fp_gmt, "endif\n");

    if (doLatLong)
        fprintf(fp_gmt, "if ($PLOT_LAT_LONG) then\n");
    fprintf(fp_gmt,
            "psxy %s $JVAL $RVAL -W2/%d/%d/%d -m -K -O >> %s.ps\n",
            fn_gmtlatlon, ired, igreen, iblue,
            /*texture,*/ fn_ps_output);
    if (doLatLong)
        fprintf(fp_gmt, "endif\n");


    return (0);

}

/*** function to generate a coded name for xy an lat/lon mapfile */

int GenMapFileName(char *fname_xy, char *fname_xz, char *fname_zy,
        char *fname_latlon, char *fnoutput, char *mapfile_name) {
    char *pstring;


    if ((pstring = strrchr(mapfile_name, '/')) != NULL)
        pstring++;
    else
        pstring = mapfile_name;
    sprintf(fname_xy, "%smap.%s.xy", fnoutput, pstring);
    sprintf(fname_xz, "%smap.%s.xz", fnoutput, pstring);
    sprintf(fname_zy, "%smap.%s.zy", fnoutput, pstring);
    sprintf(fname_latlon, "%smap.%s.latlon", fnoutput, pstring);

    return (0);
}

/*** function to convert 3-D event scatter file to 2-D GMT psxy file */

int Scat2GMT(char* fnroot_in, char* orientation, int ilonglat, char* fnscat_out) {

    char fnscat_in[FILENAME_MAX];
    FILE *fp_scat_in, *fp_scat_out;
    int npt, tot_npoints;
    float fdata[4];
    double fdata0, fdata1;
    int or_xy, or_xz, or_yz, or_zy;

    double hypox, hypoy;


    sprintf(fnscat_in, "%s.scat", fnroot_in);
    if ((fp_scat_in = fopen(fnscat_in, "r")) == NULL) {
        fprintf(stderr, "ERROR: Cannot open scatter file `%s'.\n",
                fnscat_in);
        return (-1);
    }

    if (ilonglat)
        sprintf(fnscat_out, "%s.scat.lonlat.%s",
            fnroot_in, orientation);
    else
        sprintf(fnscat_out, "%s.scat.%s", fnroot_in, orientation);
    if ((fp_scat_out = fopen(fnscat_out, "w")) == NULL) {
        fprintf(stderr,
                "ERROR: Cannot open scatter output file `%s'.\n",
                fnscat_out);
        return (-1);
    }

    or_xy = or_xz = or_yz = or_zy = 0;
    if (strcmp(orientation, "XY") == 0)
        or_xy = 1;
    else if (strcmp(orientation, "XZ") == 0)
        or_xz = 1;
    else if (strcmp(orientation, "YZ") == 0)
        or_yz = 1;
    else if (strcmp(orientation, "ZY") == 0)
        or_zy = 1;

    /* read header informaion */
    fseek(fp_scat_in, 0, SEEK_SET);
    fread(&tot_npoints, sizeof (int), 1, fp_scat_in);

    /* skip header record */
    fseek(fp_scat_in, 4 * sizeof (float), SEEK_SET);
    for (npt = 0; npt < tot_npoints; npt++) {
        fread(fdata, sizeof (float), 4, fp_scat_in);
        convertCoordsRect(proj_index_input, proj_index_output,
                fdata[0], fdata[1], &hypox, &hypoy);
        fdata[0] = hypox;
        fdata[1] = hypoy;
        if (ilonglat) {
            rect2latlon(proj_index_output, fdata[0], fdata[1], &fdata1, &fdata0);
        } else {
            fdata0 = (double) fdata[0];
            fdata1 = (double) fdata[1];
        }
        if (or_xy)
            fprintf(fp_scat_out, "%f %f\n", fdata0, fdata1);
        else if (or_xz)
            fprintf(fp_scat_out, "%f %f\n", fdata0, fdata[2]);
        else if (or_yz)
            fprintf(fp_scat_out, "%f %f\n", fdata1, fdata[2]);
        else if (or_zy)
            fprintf(fp_scat_out, "%f %f\n", fdata[2], fdata1);
    }

    fclose(fp_scat_in);
    fclose(fp_scat_out);

    return (0);

}



/** function to generate GMT plots of "traditional" statistics */

/* returns 1 if plotted, 0 if not, -1 if error */

int PlotTraditionStats(char cdatatype, char view_type, double barlen,
        Vect3D* pmax_like, Vect3D* pexpect, Ellipsoid3D *pellipsoid,
        FILE* fp_io, char *plt_code_in, char *GMTcolor) {

    double dist_cutoff, dist;

    int npts_ellipse;
    Vect3D axis1, axis2, axis3;
    Vect3D *ellArray12, *ellArray13, *ellArray23;
    Vect2D* ellipse_array;

    int i;
    double max_like_x, max_like_y, max_like_z;
    double expect_x, expect_y, expect_z;
    double xtmp, ytmp;

    char plt_code[100];

    strcpy(plt_code, plt_code_in);


    /* check plot code string */

    if ((plt_code[0] != '0' && plt_code[0] != '1')
            || (plt_code[1] != '0' && plt_code[1] != '1')
            || (plt_code[2] != '0' && plt_code[2] != '1')) {
        if (cdatatype == 'E')
            nll_puterr2("WARNING invalid entry in ellipse/statistic plot code", plt_code);
        plt_code[0] = '1';
        plt_code[1] = plt_code[2] = '0';
    }


    /* check for |MaxLike - Expect| cutoff */

    if (sscanf(plt_code, "%*c%*c%*c%lf", &dist_cutoff) == 1) {
        dist_cutoff /= 1000.0;
        dist = Dist3D(pmax_like->x, pexpect->x,
                pmax_like->y, pexpect->y,
                pmax_like->z, pexpect->z);
        if (dist > dist_cutoff)
            return (0);
    }



    // convert rect hypo coords to map rect coords
    convertCoordsRect(proj_index_input, proj_index_output,
            pmax_like->x, pmax_like->y, &max_like_x, &max_like_y);
    max_like_z = pmax_like->z;
    convertCoordsRect(proj_index_input, proj_index_output,
            pexpect->x, pexpect->y, &expect_x, &expect_y);
    expect_z = pexpect->z;


    if (plt_code[0] == '1') {

        /* plot maxlimum likelihood star */

        fprintf(fp_io,
                "# Maximum Likelihood\npsxy $JVAL $RVAL -W1/%s -Sa%lf -G%s -K -O << END >> %s.ps\n", GMTcolor, 0.4 * barlen / 2.54, GMTcolor, "${POSTSCRIPT_NAME}");
        // 20110112 AJL  "# Maximum Likelihood\npsxy $JVAL $RVAL -W1/%s -Sa%lf -K -O << END >> %s.ps\n", GMTcolor, 0.25 * barlen / 2.54, "${POSTSCRIPT_NAME}");
        if (view_type == 'Y')
            fprintf(fp_io, "%lf %lf\n", max_like_x, max_like_z);
        else if (view_type == 'H')
            fprintf(fp_io, "%lf %lf\n", max_like_x, max_like_y);
        else if (view_type == 'X')
            fprintf(fp_io, "%lf %lf\n", max_like_z, max_like_y);
        fprintf(fp_io, "END\n\n");
    }



    if (plt_code[1] == '1') {

        /* plot expectation circle */

        fprintf(fp_io,
                "# Expectation\npsxy $JVAL $RVAL -W1/%s -G0 -Sc%lf -K -O << END >> %s.ps\n",
                GMTcolor, 0.4 * barlen / 2.54, "${POSTSCRIPT_NAME}");
        if (view_type == 'Y')
            fprintf(fp_io, "%lf %lf\n", expect_x, expect_z);
        else if (view_type == 'H')
            fprintf(fp_io, "%lf %lf\n", expect_x, expect_y);
        else if (view_type == 'X')
            fprintf(fp_io, "%lf %lf\n", expect_z, expect_y);
        fprintf(fp_io, "END\n\n");
    }



    if (plt_code[2] == '1') {

        /* plot projections of error ellipsoid */

        npts_ellipse = 50;
        npts_ellipse = 20;
        ellipse_array =
                (Vect2D *) malloc((size_t) npts_ellipse
                * sizeof (Vect2D));

        fprintf(fp_io,
                "# Error Ellipsoid\npsxy $JVAL $RVAL -W1/%s -m -K -O << END >> %s.ps\n>\n",
                GMTcolor, "${POSTSCRIPT_NAME}");

        /* convert ellipsoid to 3 3D error ellipses */
        ellipsiod2Axes(pellipsoid, &axis1, &axis2, &axis3);
        ellArray12 = toEllipsoid3D(axis1, axis2, *pexpect, npts_ellipse);
        for (i = 0; i < npts_ellipse; i++) {
            convertCoordsRect(proj_index_input, proj_index_output,
                    ellArray12[i].x, ellArray12[i].y, &xtmp, &ytmp);
            ellArray12[i].x = xtmp;
            ellArray12[i].y = ytmp;
        }

        ellArray13 = toEllipsoid3D(axis1, axis3, *pexpect, npts_ellipse);
        for (i = 0; i < npts_ellipse; i++) {
            convertCoordsRect(proj_index_input, proj_index_output,
                    ellArray13[i].x, ellArray13[i].y, &xtmp, &ytmp);
            ellArray13[i].x = xtmp;
            ellArray13[i].y = ytmp;
        }
        ellArray23 = toEllipsoid3D(axis2, axis3, *pexpect, npts_ellipse);
        for (i = 0; i < npts_ellipse; i++) {
            convertCoordsRect(proj_index_input, proj_index_output,
                    ellArray23[i].x, ellArray23[i].y, &xtmp, &ytmp);
            ellArray23[i].x = xtmp;
            ellArray23[i].y = ytmp;
        }

        if (view_type == 'Y') {
            ellipse_array = Vect3D2To2D(
                    ellArray12, ellipse_array, npts_ellipse, 13);
            Vect2DArray2GMT(fp_io, ellipse_array, npts_ellipse);
            ellipse_array = Vect3D2To2D(
                    ellArray13, ellipse_array, npts_ellipse, 13);
            Vect2DArray2GMT(fp_io, ellipse_array, npts_ellipse);
            ellipse_array = Vect3D2To2D(
                    ellArray23, ellipse_array, npts_ellipse, 13);
            Vect2DArray2GMT(fp_io, ellipse_array, npts_ellipse);
        } else if (view_type == 'H') {
            ellipse_array = Vect3D2To2D(
                    ellArray12, ellipse_array, npts_ellipse, 12);
            Vect2DArray2GMT(fp_io, ellipse_array, npts_ellipse);
            ellipse_array = Vect3D2To2D(
                    ellArray13, ellipse_array, npts_ellipse, 12);
            Vect2DArray2GMT(fp_io, ellipse_array, npts_ellipse);
            ellipse_array = Vect3D2To2D(
                    ellArray23, ellipse_array, npts_ellipse, 12);
            Vect2DArray2GMT(fp_io, ellipse_array, npts_ellipse);
        } else if (view_type == 'X') {
            ellipse_array = Vect3D2To2D(
                    ellArray12, ellipse_array, npts_ellipse, 32);
            Vect2DArray2GMT(fp_io, ellipse_array, npts_ellipse);
            ellipse_array = Vect3D2To2D(
                    ellArray13, ellipse_array, npts_ellipse, 32);
            Vect2DArray2GMT(fp_io, ellipse_array, npts_ellipse);
            ellipse_array = Vect3D2To2D(
                    ellArray23, ellipse_array, npts_ellipse, 32);
            Vect2DArray2GMT(fp_io, ellipse_array, npts_ellipse);
        }
        fprintf(fp_io, "END\n\n");

        /* clean up */
        free(ellipse_array);
        free(ellArray12);
        free(ellArray13);
        free(ellArray23);

    }



    fprintf(fp_io, "\n");




    return (1);

}




/** function to generate GMT plots of focal mechanisms */

/* returns 1 if plotted, 0 if not, -1 if error */

int PlotFocalMechanism(char view_type, double scale, double magnitude, FocalMech *pfocMech,
        FILE* fp_io, char *plt_code, char *linecolor, char *fillcolor,
        double horiz_min, double horiz_max, double vert_min, double vert_max, double plot_scale) {
    double hypox, hypoy;
    char mechAval_string[MAXLINE];
    char R_string[MAXLINE];
    char J_string[MAXLINE];


    /* check plot code string */

    /*	if (plt_code[0] != '0' && plt_code[0] != '1'
                            || plt_code[1] != '0' && plt_code[1] != '1'
                            || plt_code[2] != '0' && plt_code[2] != '1') {
                    nll_puterr2("ERROR invalid entry in ellipse/statistic plot code",
                            plt_code);
                    return(-1);
            }
     */


    //	if (plt_code[0] == '1') {

    /* plot filled mechanism */

    latlon2rect(proj_index_output, pfocMech->dlat, pfocMech->dlong,
            &hypox, &hypoy);
    if (view_type == 'Y') {
        sprintf(mechAval_string, "-Ac%lf/0.0/%lf/0.0/90.0/1.0e6/%lf/%lf",
                horiz_min, horiz_max, vert_max, vert_min);
        sprintf(R_string, "-R%lf/%lf/%lf/%lf",
                0.0, horiz_max - horiz_min, -vert_max, -vert_min);
        sprintf(J_string, "-Jx%lf/%lf", plot_scale, plot_scale);
        fprintf(fp_io,
                "# Focal Mechanism\npscoupe %s %s %s %s %s -Sa%lf -N -K -O << END >> %s.ps\n",
                R_string, J_string, mechAval_string, linecolor, fillcolor,
                scale, "${POSTSCRIPT_NAME}");
        fprintf(fp_io, "%lf %lf  %lf  %lf %lf %lf  %lf 0.0 0.0\n",
                hypox, hypoy, -pfocMech->depth,
                pfocMech->dipDir - 90.0, pfocMech->dipAng, pfocMech->rake,
                magnitude
                );
        fprintf(fp_io, "END\n\n");
    } else if (view_type == 'H') {
        fprintf(fp_io,
                "# Focal Mechanism\npsmeca $JVAL $RVAL %s %s -Sa%lf -K -O << END >> %s.ps\n",
                linecolor, fillcolor, scale, "${POSTSCRIPT_NAME}");
        fprintf(fp_io, "%lf %lf  %lf  %lf %lf %lf  %lf 0.0 0.0\n",
                hypox, hypoy, 0.0,
                pfocMech->dipDir - 90.0, pfocMech->dipAng, pfocMech->rake,
                magnitude
                );
        fprintf(fp_io, "END\n\n");
    } else if (view_type == 'X') {
    }

    //	}


    printf("FOCALMECH Hyp %lf %lf %lf  Mech %lf %lf %lf mf %lf nObs %d\n",
            pfocMech->dlat, pfocMech->dlong, pfocMech->depth,
            pfocMech->dipDir, pfocMech->dipAng, pfocMech->rake,
            pfocMech->misfit, pfocMech->nObs
            );


    fprintf(fp_io, "\n");




    return (1);

}

/** function to convert residuals list to gmt psxy (residual at sta x,y) and pstext file  */


int ConvertResids2MapGMT(char* fn_nlloc_stat, char* phaseID, FILE* fp_out,
        SourceDesc* staList, int nstations, double scale, int min_num_readings) {
    int istat, nRdg, nresid;
    int c = 0;
    int ifound_readings;
    char label[MAXLINE];
    FILE *fp_in;

    char staName[20];
    char phase[20];
    double resid;

    double xsta, ysta;

    SourceDesc *staloc;


    if ((fp_in = fopen(fn_nlloc_stat, "r")) == NULL) {
        nll_puterr2("Cannot open phase statistics file", fn_nlloc_stat);
        return (-1);
    }


    // find beginning of Total Phase Corrections list
    do {
        istat = fscanf(fp_in, "%s", label);
    } while (istat != EOF && strstr(label, "Total") == NULL);


    /* read residuals */

    nRdg = 0;
    ifound_readings = 0;
    do {

        istat = fscanf(fp_in, "%s", label);
        if (istat == EOF)
            break;
        if (strcmp(label, "LOCDELAY") != 0) {
            if (ifound_readings) // do not pass to next set of readings
                break;
            while ((c = fgetc(fp_in)) != '\n' && c != EOF)
                ;
        } else {
            ifound_readings = 1;
            istat = fscanf(fp_in, "%s %s %d %lf", staName, phase, &nresid, &resid);
            if (istat == EOF)
                break;
            //printf("%s %s %d %lf\n", staName, phase, nresid, resid);

            if ((staloc = findStaLoc(staName, staList, nstations)) != NULL) {
                //printf("resid: %s %s ==? %s\n", staName, phase, phaseID);
                if (nresid >= min_num_readings && strstr(phaseID, phase) != NULL) {
                    convertCoordsRect(proj_index_input, proj_index_output,
                            staloc->x, staloc->y, &xsta, &ysta);
                    genResidualGMT(fp_out, "", resid, xsta, ysta, scale);
                    nRdg++;
                }
                //fprintf(fp_text_out, "%lf %lf 10 0 4 6 %s\n", staloc->x, staloc->y, staName);
            }
            while ((c = fgetc(fp_in)) != '\n' && c != EOF)
                ;
        }

    } while (c != EOF);


    fclose(fp_in);

    return (nRdg);

}


/** function to generate GMT code for a positive station residual */

#define RESID_MIN 0.01	// minimum size of residual symbol in inches or cm

void genResidualGMT(FILE* fp_out, char* xtra_args, double resid, double x, double y, double scale) {

    double resid_scaled, resid_plot;

    resid_scaled = resid * scale;

    if (resid_scaled >= 0.0) {
        resid_plot = resid_scaled < RESID_MIN ? RESID_MIN : resid_scaled;
        fprintf(fp_out,
                "psxy $JVAL $RVAL -Sc -W8/0/0/255 %s -K -O  << END >> ${POSTSCRIPT_NAME}.ps\n%lf %lf %lf\nEND\n",
                xtra_args, x, y, fabs(resid_plot));
    } else {
        resid_plot = resid_scaled > -RESID_MIN ? -RESID_MIN : resid_scaled;
        fprintf(fp_out,
                "psxy $JVAL $RVAL -St -W8/255/0/0 %s -K -O  << END >> ${POSTSCRIPT_NAME}.ps\n%lf %lf %lf\nEND\n",
                xtra_args, x, y, fabs(resid_plot));
    }
}

/** function to find a station name in a StaList */

SourceDesc *findStaLoc(char *staName, SourceDesc* stalist, int nstations) {

    int n;

    n = 0;
    while (n < nstations) {
        if (strcmp(stalist[n].label, staName) == 0)
            return (&(stalist[n]));
        n++;
    }

    return (NULL);

}









/*------------------------------------------------------------/ */
/* Anthony Lomax           | email: lomax@faille.unice.fr     / */
/* UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax / */
/* 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        / */
/* 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        / */
/*------------------------------------------------------------/ */

