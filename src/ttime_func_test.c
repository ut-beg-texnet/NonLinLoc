/*
 * Copyright (C) 1999-2010 Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
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


/*   ttime_func_test.c

        Program to demonstrate getting travel times through a function call.


 */

/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
161 Allee du Micocoulier, 06370 Mouans-Sartoux, France
tel: +33(0)493752502  e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


/*
        history:	(see also http://alomax.net/nlloc -> Updates)

        ver 01    24JUL2008  AJL  Original version

        see GridLib.c


.........1.........2.........3.........4.........5.........6.........7.........8

 */





#define PNAME  "ttime_func_test"

#include "GridLib.h"




/** program to demonstrate getting travel times from NLL travel time grids through a function call
 *
 *  See Makefile->ttime_func_test for compile/link requirements
 *
 *  This demonstration program reads all input from a disk file.
 *
 *  The required travel-time grids must be present as *.buf and *.hdr disk files.
 *
 */


#define NARGS_MIN 2
#define ARG_DESC "<input file>"

int main(int argc, char *argv[]) {

    int i;
    SourceDesc sourceDesc[1000];

    /** ===========================================================================
     *  variables and arrays needed by GridLib.c ReadGridFile function
     */
    float travel_times[1000];
    char filename_root[MAXLINE];
    double xloc[1000], yloc[1000], zloc[1000];
    int nvalues;
    int iSwapBytes = 0;
    SourceDesc* psrceIn = NULL; // station source to use instead of srce in grid file (use only for station=DEFAULT 2D grids, use NULL otherwise)


    // set program name
    strcpy(prog_name, PNAME);

    // check command line for correct usage
    if (argc < NARGS_MIN) {
        disp_usage(prog_name, ARG_DESC);
        return (EXIT_ERROR_USAGE);
    }

    // set file names
    // control file name
    char fn_input[MAXLINE];
    strcpy(fn_input, argv[1]);
    char in_line[4 * MAXLINE];


    SetConstants();

    // read input
    FILE* fp_input;
    if ((fp_input = fopen(fn_input, "r")) == NULL) {
        nll_puterr("FATAL ERROR: opening input file.");
        return (EXIT_ERROR_FILEIO);
    }

    // get geographic transform
    if ((fgets(in_line, 4 * MAXLINE, fp_input) == NULL) || (get_transform(0, strchr(in_line, ' '))) < 0) {
        nll_puterr("ERROR: reading transformation parameters.");
        return (EXIT_ERROR_FILEIO);
    }
    //printf("in_line 1: <%s>\n", in_line);

    // loop over station/phase
    char station[ARRIVAL_LABEL_LEN];
    char phase[PHASE_LABEL_LEN];
    double lat, lon, elev; // lat, lon, elev (km) for station, use and required only if station == "DEFAULT" 2D grids.  // 20200724 AJL - added
    SourceDesc srceIn;
    while (1) {

        // travel-time grid files path/root, station, phase
        char fn_grid_root[MAXLINE];

        if (fgets(in_line, 4 * MAXLINE, fp_input) == NULL) {
            // end of file
            break;
        }
        //printf("in_line 2: <%s>\n", in_line);
        int istat;
        if ((istat = sscanf(in_line, "%s %s %s  %lf %lf %lf", fn_grid_root, station, phase, &lat, &lon, &elev)) < 3)
            break;
        sprintf(filename_root, "%s.%s.%s.time", fn_grid_root, phase, station);
        //printf("filename_root: <%s>\n", filename_root);

        // station == "DEFAULT" 2D grids  // 20200724 AJL - added
        if (strcmp(station, "DEFAULT") == 0 && istat >= 6) {
            //    int is_coord_xyz; /* xyz coord flag */
            //    double x, y, z; /* xyz loc (km) */
            //    int is_coord_latlon; /* lat/long/depth coord flag */
            //    double dlat, dlong, depth; /* loc (lat/long (deg), depth (km)) */
            //    double otime; /* origin time */
            //    char label[SOURCE_LABEL_LEN]; /* char label */
            //    int ignored; /* 1 = ignored - not used for location */
            //    double station_weight; /* station specific weight */
            strcpy(srceIn.label, station);
            srceIn.is_coord_xyz = 0;
            srceIn.is_coord_latlon = 1;
            srceIn.dlat = lat;
            srceIn.dlong = lon;
            srceIn.depth = -elev;
            psrceIn = &srceIn;
            ConvertASourceLocation(0, psrceIn, 1, 0);
            //printf("srceIn: lat=%f lon=%f depth=%f\n", srceIn.dlat, srceIn.dlong, srceIn.z);
            // assume Global, TauP type grid byte order
            iSwapBytes = 1;
        }
        // read source locations
        nvalues = 0;
        while (fgets(in_line, 4 * MAXLINE, fp_input) != NULL) {
            //printf("in_line 3: <%s>\n", in_line);
            if (strncmp(in_line, "X", 1) == 0)
                break;
            if (GetSource(in_line, sourceDesc + nvalues, nvalues) < 0)
                break;
            nvalues++;
        }
        // convert source location coordinates
        if (ConvertSourceLoc(0, sourceDesc, nvalues, 1, 1) < 0)
            nll_puterr("ERROR: converting source locations to x/y and lat/lon.");
        for (i = 0; i < nvalues; i++) {
            xloc[i] = sourceDesc[i].x;
            yloc[i] = sourceDesc[i].y;
            zloc[i] = sourceDesc[i].z;
        }

        /** ===========================================================================
         *  call function that gets values from grid (see GridLib.c ReadGridFile())
         */
        ReadGridFile(travel_times, filename_root, "time", xloc, yloc, zloc, nvalues, iSwapBytes, psrceIn);

        // write travel times to stdout
        for (i = 0; i < nvalues; i++)
            printf("%s %s  <-  x=%f y=%f z=%f (lat=%f lon=%f depth=%f)  ttime=%f\n",
                station, phase, xloc[i], yloc[i], zloc[i], sourceDesc[i].dlat, sourceDesc[i].dlong, sourceDesc[i].depth, travel_times[i]);



    }


    // clean up
    fclose(fp_input);


    return (0);

}




