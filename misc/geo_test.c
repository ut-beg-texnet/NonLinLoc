/*
 * Copyright (C) 2013 Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
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


/*   mag_func_test.c

        Program to demonstrate calculating magnitudes through a function call.


 */

/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
161 Allee du Micocoulier, 06370 Mouans-Sartoux, France
tel: +33(0)493752502  e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


/*
        history:	(see also http://alomax.net/nlloc -> Updates)

        ver 01    27MAR1023  AJL  Original version

        see GridLib.c


.........1.........2.........3.........4.........5.........6.........7.........8

 */





#define PNAME  "geo_test"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "geo.h"


/** program to test geo.h function calls
 *
 */


#define NARGS_MIN 5
#define ARG_DESC "lat1 lon1 lat2 lon2"

int main(int argc, char *argv[]) {


    // check command line for correct usage
    if (argc < NARGS_MIN) {
        fprintf(stderr, "Usage: %s %s\n", PNAME, ARG_DESC);
        return (-1);
    }

    double from_latitude = atof(argv[1]);
    double from_longitude = atof(argv[2]);
    double to_latitude = atof(argv[3]);
    double to_longitude = atof(argv[4]);

    printf("Lat: %f  Lon: %f -> ", from_latitude, from_longitude);
    printf("Lat: %f  Lon: %f\n", to_latitude, to_longitude);
    printf("\tdist: %f\n", GCDistance(from_latitude, from_longitude, to_latitude, to_longitude));
    printf("\tazimuth: %f\n", GCAzimuth(from_latitude, from_longitude, to_latitude, to_longitude));

    return (0);

}




