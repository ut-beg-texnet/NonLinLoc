//	geo.cpp
//
// Adapted from:
//    Geographic Distance and Azimuth Calculations
//    by Andy McGovern
//    http://www.codeguru.com/algorithms/GeoCalc.html

//#include "stdafx.h"
#include "geo.h"

#include <math.h>
#include <float.h>

double GCDistance(double lat1, double lon1, double lat2, double lon2) {
    double d;

    if (lat1 == lat2 && lon1 == lon2)
        return (0.0);

    lat1 *= DE2RA;
    lon1 *= DE2RA;
    lat2 *= DE2RA;
    lon2 *= DE2RA;

    d = sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2);
    return (AVG_ERAD * acos(d));
}

double GCAzimuth(double lat1, double lon1, double lat2, double lon2) {

    double lonA = lon1 * DE2RA;
    double latA = lat1 * DE2RA;
    double lonB = lon2 * DE2RA;
    double latB = lat2 * DE2RA;

    // distance
    double dist =
            acos(
            sin(latA) * sin(latB)
            + cos(latA) * cos(latB) * cos((lonB - lonA)));

    // 20141106 AJL - added following check, to prevent div by 0 of sin(dist))
    if (dist < FLT_MIN) {
        return (0.0);
    }

    // azimuth
    double cosAzimuth =
            (cos(latA) * sin(latB)
            - sin(latA) * cos(latB)
            * cos((lonB - lonA)))
            / sin(dist);
    double sinAzimuth =
            cos(latB) * sin((lonB - lonA)) / sin(dist);
    double az = atan2(sinAzimuth, cosAzimuth) / DE2RA;

    if (isnan(az) && fabs(lon2 - lon1) < 0.000001) {
        if (lat1 > lat2) {
            az = 180.0;
        } else {
            az = 0.0;
        }
    }

    if (az < 0.0) {
        az += 360.0;
    }

    return (az);
}

double GCAzimuth_ERROR(double lat1, double lon1, double lat2, double lon2) {
    double result = 0.0;

    double c, A;

    long ilat1 = (long) (0.50 + lat1 * 360000.0);
    long ilat2 = (long) (0.50 + lat2 * 360000.0);
    long ilon1 = (long) (0.50 + lon1 * 360000.0);
    long ilon2 = (long) (0.50 + lon2 * 360000.0);

    lat1 *= DE2RA;
    lon1 *= DE2RA;
    lat2 *= DE2RA;
    lon2 *= DE2RA;

    if ((ilat1 == ilat2) && (ilon1 == ilon2)) {
        return result;
    } else if (ilon1 == ilon2) {
        if (ilat1 > ilat2)
            result = 180.0;
    } else {
        c = acos(sin(lat2) * sin(lat1) + cos(lat2) * cos(lat1) * cos((lon2 - lon1)));
        A = asin(cos(lat2) * sin((lon2 - lon1)) / sin(c));
        result = (A * RA2DE);

        if ((ilat2 > ilat1) && (ilon2 > ilon1)) {
        } else if ((ilat2 < ilat1) && (ilon2 < ilon1)) {
            result = 180.0 - result;
        } else if ((ilat2 < ilat1) && (ilon2 > ilon1)) {
            result = 180.0 - result;
        } else if ((ilat2 > ilat1) && (ilon2 < ilon1)) {
            result += 360.0;
        }
    }

    return result;
}

double ApproxDistance(double lat1, double lon1, double lat2, double lon2) {

    double F, G, L, sing, cosl, cosf, sinl, sinf, cosg, S, C, W, R, H1, H2, D;

    lat1 = DE2RA * lat1;
    lon1 = -DE2RA * lon1;
    lat2 = DE2RA * lat2;
    lon2 = -DE2RA * lon2;

    F = (lat1 + lat2) / 2.0;
    G = (lat1 - lat2) / 2.0;
    L = (lon1 - lon2) / 2.0;

    sing = sin(G);
    cosl = cos(L);
    cosf = cos(F);
    sinl = sin(L);
    sinf = sin(F);
    cosg = cos(G);

    S = sing * sing * cosl * cosl + cosf * cosf * sinl*sinl;
    C = cosg * cosg * cosl * cosl + sinf * sinf * sinl*sinl;
    W = atan2(sqrt(S), sqrt(C));
    R = sqrt((S * C)) / W;
    H1 = (3 * R - 1.0) / (2.0 * C);
    H2 = (3 * R + 1.0) / (2.0 * S);
    D = 2 * W * ERAD;
    return (D * (1 + FLATTENING * H1 * sinf * sinf * cosg * cosg -
            FLATTENING * H2 * cosf * cosf * sing * sing));
}

double EllipsoidDistance(double lat1, double lon1, double lat2, double lon2) {
    double distance = 0.0;
    double faz, baz;
    double r = 1.0 - FLATTENING;
    double tu1, tu2, cu1, su1, cu2, x, sx, cx, sy, cy, y, sa, c2a, cz, e, c, d;
    double cosy1, cosy2;
    distance = 0.0;

    if ((lon1 == lon2) && (lat1 == lat2)) return distance;
    lon1 *= DE2RA;
    lon2 *= DE2RA;
    lat1 *= DE2RA;
    lat2 *= DE2RA;

    cosy1 = cos(lat1);
    cosy2 = cos(lat2);

    if (cosy1 == 0.0) cosy1 = 0.0000000001;
    if (cosy2 == 0.0) cosy2 = 0.0000000001;

    tu1 = r * sin(lat1) / cosy1;
    tu2 = r * sin(lat2) / cosy2;
    cu1 = 1.0 / sqrt(tu1 * tu1 + 1.0);
    su1 = cu1 * tu1;
    cu2 = 1.0 / sqrt(tu2 * tu2 + 1.0);
    x = lon2 - lon1;

    distance = cu1 * cu2;
    baz = distance * tu2;
    faz = baz * tu1;

    do {
        sx = sin(x);
        cx = cos(x);
        tu1 = cu2 * sx;
        tu2 = baz - su1 * cu2 * cx;
        sy = sqrt(tu1 * tu1 + tu2 * tu2);
        cy = distance * cx + faz;
        y = atan2(sy, cy);
        sa = distance * sx / sy;
        c2a = -sa * sa + 1.0;
        cz = faz + faz;
        if (c2a > 0.0) cz = -cz / c2a + cy;
        e = cz * cz * 2. - 1.0;
        c = ((-3.0 * c2a + 4.0) * FLATTENING + 4.0) * c2a * FLATTENING / 16.0;
        d = x;
        x = ((e * cy * c + cz) * sy * c + y) * sa;
        x = (1.0 - c) * x * FLATTENING + lon2 - lon1;
    } while (fabs(d - x) > EPS);

    x = sqrt((1.0 / r / r - 1.0) * c2a + 1.0) + 1.0;
    x = (x - 2.0) / x;
    c = 1.0 - x;
    c = (x * x / 4.0 + 1.0) / c;
    d = (0.375 * x * x - 1.0) * x;
    x = e * cy;
    distance = 1.0 - e - e;
    distance = ((((sy * sy * 4.0 - 3.0) *
            distance * cz * d / 6.0 - x) * d / 4.0 + cz) * sy * d + y) * c * ERAD * r;

    return distance;
}


/** calculate end point (latitude/longitude) given a starting point, distance (deg), and azimuth (deg).
 *
 * Adapted from:
 *    http://openmap.bbn.com/doc/api/com/bbn/openmap/proj/GreatCircle.html
 *
 * Calculate point at azimuth and distance from another point.
 * <p>
 * Returns a LatLonPoint at arc distance `c' in direction `Az'
 * from start point.
 * <p>
 *
 * @param phi1 latitude in radians of start point
 * @param lambda0 longitude in radians of start point
 * @param c arc radius in radians (0 &lt; c &lt;= PI)
 * @param Az azimuth (direction) east of north (-PI &lt;= Az &lt;
 *        PI)
 * @return LatLonPoint
 *
 */
void PointAtGCDistanceAzimuth(double lat1, double lon1, double dist, double az, double* lat2, double* lon2) {

    double phi1 = DE2RA * lat1;
    double lambda0 = DE2RA * lon1;
    double c = DE2RA * dist;
    double Az = DE2RA * az;

    double cosphi1 = cos(phi1);
    double sinphi1 = sin(phi1);
    double cosAz = cos(Az);
    double sinAz = sin(Az);
    double sinc = sin(c);
    double cosc = cos(c);

    *lat2 = RA2DE * asin(sinphi1 * cosc + cosphi1 * sinc * cosAz);
    *lon2 = RA2DE * (atan2(sinc * sinAz, cosphi1 * cosc - sinphi1 * sinc * cosAz) + lambda0);

}
