/*
 * Copyright (C) 2017 Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
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


/*   lonlat_2_xy.c

   Program to convert longitude/latitude values to x/y in given transformation


 */

/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


/*
        history:	(see also http://alomax.net/nlloc -> Updates)

        ver 01    30JUN2022  TM   Original version

        see GridLib.c


.........1.........2.........3.........4.........5.........6.........7.........8

 */





#define PNAME  "lonlat_2_xy"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../src/include/GridLib.h"


/** program to test geo.h function calls
 *
 */


#define NARGS_MIN 0
#define ARG_DESC ""

int main(int argc, char *argv[]) {

    // check command line for correct usage
    if (argc < NARGS_MIN) {
        fprintf(stderr, "Usage: %s %s\n", PNAME, ARG_DESC);
        return (-1);
    }

    SetConstants();

// emulate https://epsg.io/25832
// TRANS_MERC GRS-80 0.0 9.0 0.0 1 500000 0.9996
//#define EARTH_ELLIPSOID "GRS-80"
//#define MAP_ORIG_LAT 0.0
//#define MAP_ORIG_LON 9.0
//#define ROTATION_ANGLE 0.0
//#define USE_FALSE_EASTING 1
//#define FALSE_EASTING 500000
//#define SCALE_FACTOR 0.9996

    char trans[999];

    printf("Enter transformation string (the part after 'TRANS '): ");
    scanf("%[^\n]", &trans);
    //sprintf(trans[n++], "TRANS_MERC %s %f  %f  %f %d %d %f ",
    //        EARTH_ELLIPSOID, MAP_ORIG_LAT, MAP_ORIG_LON, ROTATION_ANGLE, USE_FALSE_EASTING, FALSE_EASTING, SCALE_FACTOR);

    get_transform(0, trans);

    float TEST_ORIG_LAT;
    float TEST_ORIG_LON;

    printf("Enter longitude and latitude, separated by space: ");
    scanf("%f %f", &TEST_ORIG_LON, &TEST_ORIG_LAT);
    printf("Longitude: %f\n", TEST_ORIG_LON);
    printf("Latitude: %f\n", TEST_ORIG_LAT);

    if (strcmp(map_trans_type[0], "GLOBAL") == 0) {
        GeometryMode = MODE_GLOBAL;
    } else {
        GeometryMode = MODE_RECT;
    }
    double xrect, yrect;
    latlon2rect(0, TEST_ORIG_LAT, TEST_ORIG_LON, &xrect, &yrect);
    printf("X (in m): %.2f\n", xrect * 1000);
    printf("Y (in m): %.2f\n", yrect * 1000);
    double dlat, dlon;
    rect2latlon(0, xrect, yrect, &dlat, &dlon);
    printf("delta lon (lon/lat -> X/Y -> lon/lat): %g\n", TEST_ORIG_LON - dlon);
    printf("delta lat (lon/lat -> X/Y -> lon/lat): %g\n", TEST_ORIG_LAT - dlat);
    printf("\n");

    return (0);

}

