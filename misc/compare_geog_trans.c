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


/*   compare_geog_trans.c

        Program to compare distances and azimuths from a point with different geographic transformations.


 */

/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


/*
        history:	(see also http://alomax.net/nlloc -> Updates)

        ver 01    22NOV2017  AJL  Original version

        see GridLib.c


.........1.........2.........3.........4.........5.........6.........7.........8

 */





#define PNAME  "compare_geog_trans"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "gridlib.h"


/** program to test geo.h function calls
 *
 */


#define NARGS_MIN 0
#define ARG_DESC ""

double getAzim(double x_from, double y_from, double x_to, double y_to);
double getDist(double x_from, double y_from, double x_to, double y_to);

int main(int argc, char *argv[]) {

    // check command line for correct usage
    if (argc < NARGS_MIN) {
        fprintf(stderr, "Usage: %s %s\n", PNAME, ARG_DESC);
        return (-1);
    }

    SetConstants();

#define EARTH_ELLIPOID "Sphere"
    //#define EARTH_ELLIPOID "WGS-72"
    //#define EARTH_ELLIPOID "WGS-84"

#define MAP_ORIG_LAT 31.31
#define MAP_ORIG_LON -103.56

//#define DO_MAP_ORIG
#ifdef DO_MAP_ORIG
#define TEST_ORIG_LAT MAP_ORIG_LAT
#define TEST_ORIG_LON MAP_ORIG_LON
#else
// 10kmNE
//#define TEST_ORIG_LAT 31.41
//#define TEST_ORIG_LON -103.46
// 100kmNE
#define TEST_ORIG_LAT 32.31
#define TEST_ORIG_LON -102.56
#endif


#define DISTANCE 100.0      // fixed distance (km) for range of azimuth
#define AZIMUTH 45.0      // fixed azimuth (deg) for range of distance

    int num_projections = 5;
    char trans[num_projections][999];

    int n = 0;
    sprintf(trans[n++], "GLOBAL ");
    sprintf(trans[n++], "SIMPLE  %f  %f  0 ", MAP_ORIG_LAT, MAP_ORIG_LON);
    sprintf(trans[n++], "LAMBERT %s %f  %f  30.31 32.31 0 ", EARTH_ELLIPOID, MAP_ORIG_LAT, MAP_ORIG_LON);
    sprintf(trans[n++], "TRANS_MERC %s %f  %f  0 ", EARTH_ELLIPOID, MAP_ORIG_LAT, MAP_ORIG_LON);
    sprintf(trans[n++], "AZIMUTHAL_EQUIDIST %s %f  %f  0 ", EARTH_ELLIPOID, MAP_ORIG_LAT, MAP_ORIG_LON);

    for (int n_proj = 0; n_proj < num_projections; n_proj++) {
        get_transform(n_proj, trans[n_proj]);
    }

    SourceDesc from_point;
    from_point.is_coord_latlon = 1;
    from_point.is_coord_xyz = 0;
    from_point.dlat = TEST_ORIG_LAT;
    from_point.dlong = TEST_ORIG_LON;


    printf("lat0\t lon0\t lat1\t lon1\t dist\t az\t lat2\t lon2");
    for (int n_proj = 0; n_proj < num_projections; n_proj++) {
        printf("\t %s_dist\t %s_az", map_trans_type[n_proj], map_trans_type[n_proj]);
    }
    printf("\n");

//#define DO_AZIMUTH
#ifdef DO_AZIMUTH
    double distance = DISTANCE;
    double distance_step = 0.0;
    double azimuth = 0.0;
    double azimuth_step = 2.0;
#else
    double distance = 2.0;
    double distance_step = 2.0;
    double azimuth = AZIMUTH;
    double azimuth_step = 0.0;
#endif
    // test
    //azimuth = 0.0;
    //azimuth_step = 45.0;

    while ((azimuth_step < FLT_MIN && distance < 200.0) || (distance_step < FLT_MIN && azimuth < 360.0)) {

        double to_latitude;
        double to_longitude;

        PointAtGCDistanceAzimuth(TEST_ORIG_LAT, TEST_ORIG_LON, distance * KM2DEG, azimuth, &to_latitude, &to_longitude);

        printf("%.3f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f",
                TEST_ORIG_LAT, TEST_ORIG_LON,
                from_point.dlat, from_point.dlong,
                distance, azimuth,
                to_latitude, to_longitude
                );

        for (int n_proj = 0; n_proj < num_projections; n_proj++) {
            if (strcmp(map_trans_type[n_proj], "GLOBAL") == 0) {
                GeometryMode = MODE_GLOBAL;
            } else {
                GeometryMode = MODE_RECT;
            }
            ConvertASourceLocation(n_proj, &from_point, 1, 0);
            double xrect, yrect;
            latlon2rect(n_proj, to_latitude, to_longitude, &xrect, &yrect);
            double dist = getDist(from_point.x, from_point.y, xrect, yrect);
            double az = getAzim(from_point.x, from_point.y, xrect, yrect);
            if (az >= 359.999) {
                az = 0.0;
            }
            double az_error = az - azimuth;
            if (az_error > 180.0) {
                az_error -= 360.0;
            } else if (az_error < -180.0) {
                az_error += 360.0;
            }
            printf("\t %.3f\t %.3f", dist - distance, az_error);
            /*printf(" [%.3f, %.3f]", from_point.dlat, from_point.dlong);
            printf(" [%.3f, %.3f]", from_point.x, from_point.y);
            printf(" [%.3f, %.3f]", xrect, yrect);
            printf(" [%d]", GeometryMode);*/
        }
        printf("\n");

        azimuth += azimuth_step;
        distance += distance_step;

    }

    return (0);

}

/** function to calculate epicentral (horizontal) azimuth */

double getAzim(double x_from, double y_from, double x_to, double y_to) {
    double xtmp, ytmp, azim;

    if (GeometryMode == MODE_GLOBAL) {
        return (GCAzimuth(y_from, x_from, y_to, x_to));
    } else {
        xtmp = x_to - x_from;
        ytmp = y_to - y_from;
        azim = atan2(xtmp, ytmp) / cRPD;
        if (azim < 0.0) {
            azim += 360.0;
        } else if (azim >= 360.0) {
            azim -= 360.0;
        }

        return (azim);
    }
}

/** function to calculate epicentral (horizontal) distance */

double getDist(double x_from, double y_from, double x_to, double y_to) {
    double xtmp, ytmp;

    if (GeometryMode == MODE_GLOBAL) {
        return (GCDistance(y_from, x_from, y_to, x_to));
        //return(EllipsoidDistance(y_to, x_to, psrce->y_from, psrce->x_from));
    } else {
        xtmp = x_to - x_from;
        ytmp = y_to - y_from;

        return (sqrt(xtmp * xtmp + ytmp * ytmp));
    }
}




