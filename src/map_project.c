
/*--------------------------------------------------------------------
    WARNING - this file is not gmt_map.c, it contains a subset of gmt_map.c

    subset of functions selected by A. Lomax, June 1998

    changes or additions indicated by                                 */

/*AJL ... AJL*/

/*AJL*/ /*END AJL*/

/*--------------------------------------------------------------------*/





/*--------------------------------------------------------------------
 *    The GMT-system:	@(#)gmt_map.c	2.56  09 Aug 1995
 *
 *    Copyright (c) 1991-1995 by P. Wessel and W. H. F. Smith
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*
 *
 *			G M T _ M A P . C
 *
 *- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * gmt_map.c contains code related to coordinate transformation
 *
 * 	Map Transformation Setup Routines
 *	These routines initializes the selected map transformation
 *	The names and main function are listed below
 *	NB! Note that the transformation function does not check that they are
 *	passed valid lon,lat numbers. I.e asking for log10 scaling using values
 *	<= 0 results in problems.
 *
 * Map_projections include functions that will set up the transformation
 * between xy and latlon for several map projections.
 *
 * A few of the core coordinate transformation functions are based on similar
 * FORTRAN routines written by Pat Manley, Doug Shearer, and Bill Haxby, and
 * have been rewritten in C and subsequently streamlined.  The rest is based
 * on P. Snyder, "Map Projections - a working manual", USGS Prof paper 1395.
 *
 * Transformations supported (both forward and inverse):
 *	Linear x/y[/z] scaling
 *	Polar (theta, r) scaling
 *	Mercator
 *	Stereographic
 *	Albers Equal-Area Conic
 *	Lambert Conformal Conic
 *	Cassini Cylindrical
 *	Oblique Mercator
 *	TM Transverse Mercator
 *	UTM Universal Transverse Mercator
 *	Lambert Azimuthal Equal-Area
 *	Mollweide Equal-Area
 *	Hammer-Aitoff Equal-Area
 *	Sinusoidal Equal-Area
 *	Winkel Tripel
 *	Orthographic
 *	Azimuthal Equidistant
 *	Robinson
 *	Eckert IV
 *	Cylindrical Equal-area (e.g., Peters, Gall, Behrmann)
 *	Cylindrical Equidistant (Plate Carree)
 *
 * The ellipsoid used is selectable by editing the .gmtdefaults in your
 * home directory.  If no such file, create one by running gmtdefaults.
 *
 * Usage: Initialize system by calling map_setup (separate module), and
 * then just use geo_to_xy() and xy_to_geo() functions.
 *
 * Author:	Paul Wessel
 * Date:	27-JUL-1992
 * Version:	v2.1
 *
 *
 * Functions include:
 *
...

 */


/*AJL*/

/*#include "gmt.h"*/

#include <string.h>
#include <math.h>

#include "map_project.h"

#define PI_2 (2.0*M_PI)
#define D2R (M_PI/180.0)
#define R2D (180.0/M_PI)

#define SMALL 1.0e-10

/*#define d_log(x) ((x) <= 0.0 ? gmt_NaN : log (x))*/
#define d_log(x) ((x) <= 0.0 ? -1.0e10 : log (x))
#define d_sqrt(x) ((x) < 0.0 ? 0.0 : sqrt (x))

/*double gmt_NaN;*/

typedef int BOOLEAN; /* BOOLEAN used for logical variables */

struct ELLIPSOID {
    char name[20];
    int date;
    double eq_radius;
    double pol_radius;
    double flattening;
};

/* Information about a particular ellipsoid */
/* Table taken from Snyder "Map projection - a working manual",
                p 12 Table 1 */

#define N_ELLIPSOIDS 15

struct ELLIPSOID ellipse[N_ELLIPSOIDS] = {
    // name, date, eq_radius, pol_radius, flattening
    { "WGS-84", 1984, 6378137.0, 6356752.1, 1.0 / 298.254},
    { "GRS-80", 1980, 6378137.0, 6356752.3, 1.0 / 298.257},
    { "WGS-72", 1972, 6378135.0, 6356750.5, 1.0 / 298.26},
    { "Australian", 1965, 6378160.0, 6356774.7, 1.0 / 298.25},
    { "Krasovsky", 1940, 6378245.0, 6356863.0, 1.0 / 298.3},
    { "International", 1924, 6378388.0, 6356911.9, 1.0 / 297.0},
    { "Hayford-1909", 1909, 6378388.0, 6356911.9, 1.0 / 297.0},
    { "Clarke-1880", 1880, 6378249.1, 6356514.9, 1.0 / 293.46},
    { "Clarke-1866", 1866, 6378206.4, 6356583.8, 1.0 / 294.98},
    { "Airy", 1830, 6377563.4, 6356256.9, 1.0 / 299.32},
    { "Bessel", 1841, 6377397.2, 6356079.0, 1.0 / 299.15},
    { "Hayford-1830", 1830, 6377276.3, 6356075.4, 1.0 / 300.80},
    { "Sphere", 1980, 6371008.7714, 6371008.7714, 0.0},
        /* https://gisgeography.com/geodetic-datums-nad27-nad83-wgs84/
        NAD27 Datum vs NAD83 Datum
        The NAD27 datum was based on the Clarke Ellipsoid of 1866:
        Semi-major axis: 6,378,206.4 m
        Semi-minor axis: 6,356,583.8 m
        Inverse flattening: 294.98
        The NAD83 datum was based on the Geodetic Reference System (GRS80) Ellipsoid:
        Semi-major axis: 6,378,137.0 m
        Semi-minor axis: 6,356,752.3 m
        Inverse flattening: 298.26
         */
    { "NAD-27", 1927, 6378206.4, 6356583.8, 1.0 / 294.98},
    { "NAD-83", 1983, 6378137.4, 6356752.3, 1.0 / 294.26}

};


// number of projections supported
#define NUM_PROJ_MAX 10

double EQ_RAD[NUM_PROJ_MAX];
double ECC[NUM_PROJ_MAX], ECC2[NUM_PROJ_MAX], ECC4[NUM_PROJ_MAX], ECC6[NUM_PROJ_MAX];
//double M_PR_DEG;

/* fields from struct MAP_PROJECTIONS taken from gmt_project.h,
        and converted to globals.
        WARNING - many fields removed! */

BOOLEAN NorthPole[NUM_PROJ_MAX]; /* TRUE if projection is on northern
					  hermisphere, FALSE on southern */
double CentralMeridian[NUM_PROJ_MAX]; /* Central meridian for projection */
double Pole[NUM_PROJ_MAX]; /* +90 pr -90, depending on which pole */


/* Lambert conformal conic parameters.
                (See Snyder for details on all parameters) */

double LambertConfConic_N[NUM_PROJ_MAX];
double LambertConfConic_F[NUM_PROJ_MAX];
double LambertConfConic_rho0[NUM_PROJ_MAX];




/* set constants that were set in GMT function map_setup */

/* use values from gmt_defaults.h: */

int map_setup_proxy(int n_proj, char* ellipsoid_name) {

    int num_ellipsoid;
    double f;

    /*mknan (gmt_NaN);*/

    /* determine ellipsoid */
    for (num_ellipsoid = 0;
            num_ellipsoid < N_ELLIPSOIDS
            && strcmp(ellipsoid_name, ellipse[num_ellipsoid].name);
            num_ellipsoid++);
    if (num_ellipsoid == N_ELLIPSOIDS)
        return (-1);

    EQ_RAD[n_proj] = ellipse[num_ellipsoid].eq_radius;
    f = ellipse[num_ellipsoid].flattening;
    ECC2[n_proj] = 2 * f - f * f;
    ECC4[n_proj] = ECC2[n_proj] * ECC2[n_proj];
    ECC6[n_proj] = ECC2[n_proj] * ECC4[n_proj];
    ECC[n_proj] = d_sqrt(ECC2[n_proj]);

    return (0);

}

/*END AJL*/




/*
 *	TRANSFORMATION ROUTINES FOR THE LAMBERT CONFORMAL CONIC PROJECTION
 */

/*** function to set up a Lambert Conformal Conic projection */

int vlamb(n_proj, rlong0, rlat0, pha, phb)
int n_proj;
double rlong0, rlat0, pha, phb;
{

    double t_pha, m_pha, t_phb, m_phb, t_rlat0;

    NorthPole[n_proj] = (rlat0 > 0.0);
    // AJL 20090812 BUG FIX
    Pole[n_proj] = (NorthPole[n_proj]) ? 90.0 : -90.0;
    //Pole[n_proj] = (NorthPole) ? 90.0 : -90.0;
    pha *= D2R;
    phb *= D2R;

    t_pha = tan(45.0 * D2R - 0.5 * pha) / pow((1.0 - ECC[n_proj] *
            sin(pha)) / (1.0 + ECC[n_proj] * sin(pha)), 0.5 * ECC[n_proj]);
    m_pha = cos(pha) / d_sqrt(1.0 - ECC2[n_proj]
            * pow(sin(pha), 2.0));
    t_phb = tan(45.0 * D2R - 0.5 * phb) / pow((1.0 - ECC[n_proj] *
            sin(phb)) / (1.0 + ECC[n_proj] * sin(phb)), 0.5 * ECC[n_proj]);
    m_phb = cos(phb) / d_sqrt(1.0 - ECC2[n_proj]
            * pow(sin(phb), 2.0));
    t_rlat0 = tan(45.0 * D2R - 0.5 * rlat0 * D2R) /
            pow((1.0 - ECC[n_proj] * sin(rlat0 * D2R)) /
            (1.0 + ECC[n_proj] * sin(rlat0 * D2R)), 0.5 * ECC[n_proj]);

    if (pha != phb)
        LambertConfConic_N[n_proj] = (d_log(m_pha) - d_log(m_phb))
        / (d_log(t_pha) - d_log(t_phb));
    else
        LambertConfConic_N[n_proj] = sin(pha);

    LambertConfConic_F[n_proj] = m_pha / (LambertConfConic_N[n_proj] *
            pow(t_pha, LambertConfConic_N[n_proj]));
    CentralMeridian[n_proj] = rlong0;
    LambertConfConic_rho0[n_proj] = EQ_RAD[n_proj] * LambertConfConic_F[n_proj] *
            pow(t_rlat0, LambertConfConic_N[n_proj]);

    return (0);
}

/*** function to do x,y to lat,long Lambert Conformal Conic projection */

int lamb(n_proj, lon, lat, x, y)
int n_proj;
double lon, lat, *x, *y;
{
    double rho, theta, hold1, hold2, hold3;

    while ((lon - CentralMeridian[n_proj]) < -180.0) lon += 360.0;
    while ((lon - CentralMeridian[n_proj]) > 180.0) lon -= 360.0;
    lat *= D2R;

    hold2 = pow(((1.0 - ECC[n_proj] * sin(lat))
            / (1.0 + ECC[n_proj] * sin(lat))), 0.5 * ECC[n_proj]);
    hold3 = tan(45.0 * D2R - 0.5 * lat);
    if (fabs(hold3) < SMALL) hold3 = 0.0;
    hold1 = (hold3 == 0.0) ? 0.0 : pow(hold3 / hold2, LambertConfConic_N[n_proj]);
    rho = EQ_RAD[n_proj] * LambertConfConic_F[n_proj] * hold1;
    theta = LambertConfConic_N[n_proj] * (lon - CentralMeridian[n_proj]) * D2R;

    *x = rho * sin(theta);
    *y = LambertConfConic_rho0[n_proj] - rho * cos(theta);

    return (0);
}

/*** function to do lat,long to x,y inverse
                        Lambert Conformal Conic projection */

int ilamb(n_proj, lon, lat, x, y)
int n_proj;
double *lon, *lat, x, y;
{
    int i;
    double theta, temp, rho, t, tphi, phi, delta;

    theta = atan(x / (LambertConfConic_rho0[n_proj] - y));
    *lon = (theta / LambertConfConic_N[n_proj]) * R2D + CentralMeridian[n_proj];

    temp = x * x + (LambertConfConic_rho0[n_proj] - y)
            * (LambertConfConic_rho0[n_proj] - y);
    rho = copysign(d_sqrt(temp), LambertConfConic_N[n_proj]);
    t = pow((rho / (EQ_RAD[n_proj] * LambertConfConic_F[n_proj])), 1. / LambertConfConic_N[n_proj]);
    tphi = 90.0 * D2R - 2.0 * atan(t);
    delta = 1.0;
    for (i = 0; i < 100 && delta > 1.0e-8; i++) {
        temp = (1.0 - ECC[n_proj] * sin(tphi)) / (1.0 + ECC[n_proj] * sin(tphi));
        phi = 90.0 * D2R - 2.0 * atan(t * pow(temp, 0.5 * ECC[n_proj]));
        delta = fabs(fabs(tphi) - fabs(phi));
        tphi = phi;
    }
    *lat = phi * R2D;

    return (0);
}

/* Transverse Mercator Projection (TM) */

struct TRANS_MERCATOR {
    BOOLEAN north_pole; /* TRUE if projection is on northern hemisphere, FALSE on southern */

    int use_false_easting;  // flag to apply false easting (500km added to X when converting geog->UTM, and v.v.)
    double central_meridian; // Central meridian for projection
    double y_central_parralel; // y offset of central parallel
    double t_e2;
    double t_c1, t_c2, t_c3, t_c4;
    double t_ic1, t_ic2, t_ic3, t_ic4;

};
struct TRANS_MERCATOR TransverseMercator[NUM_PROJ_MAX];

double map_scale_factor = 1.0;

/*
 *	TRANSFORMATION ROUTINES FOR THE Transverse Mercator Projection (TM)
 */

/*
int map_init_tm(int n_proj) {
    BOOLEAN search;
    double xmin, xmax, ymin, ymax;

    vtm(TransverseMercator[n_proj].pars[0]);
    if (TransverseMercator[n_proj].units_pr_degree) TransverseMercator[n_proj].pars[1] /= M_PR_DEG;
    TransverseMercator[n_proj].x_scale = TransverseMercator[n_proj].y_scale = TransverseMercator[n_proj].pars[1];
    forward = (PFI) tm;
    inverse = (PFI) itm;

    if (TransverseMercator[n_proj].region) {
        xy_search(&xmin, &xmax, &ymin, &ymax);
        outside = (PFI) wesn_outside;
        crossing = (PFI) wesn_crossing;
        overlap = (PFI) wesn_overlap;
        map_clip = (PFI) wesn_clip;
        left_edge = (PFD) left_rect;
        right_edge = (PFD) right_rect;
        search = FALSE;
    } else { // Find min values
        (*forward) (TransverseMercator[n_proj].w, TransverseMercator[n_proj].s, &xmin, &ymin);
        (*forward) (TransverseMercator[n_proj].e, TransverseMercator[n_proj].n, &xmax, &ymax);
        outside = (PFI) rect_outside;
        crossing = (PFI) rect_crossing;
        overlap = (PFI) rect_overlap;
        map_clip = (PFI) rect_clip;
        left_edge = (PFD) left_rect;
        right_edge = (PFD) right_rect;
        frame_info.check_side = TRUE;
        search = TRUE;
    }

    frame_info.horizontal = TRUE;
    map_setinfo(xmin, xmax, ymin, ymax, TransverseMercator[n_proj].pars[1]);

    gmtdefs.basemap_type = 1;

    return (search);
}*/

void vtm(int n_proj, double lon0, double lat0, int use_false_easting) {
    /* Set up an TM projection */
    double e1;

    e1 = (1.0 - d_sqrt(1.0 - ECC2[n_proj])) / (1.0 + d_sqrt(1.0 - ECC2[n_proj]));
    TransverseMercator[n_proj].t_e2 = ECC2[n_proj] / (1.0 - ECC2[n_proj]);
    TransverseMercator[n_proj].t_c1 = (1.0 - 0.25 * ECC2[n_proj] - 3.0 * ECC4[n_proj] / 64.0 - 5.0 * ECC6[n_proj] / 256.0);
    TransverseMercator[n_proj].t_c2 = (3.0 * ECC2[n_proj] / 8.0 + 3.0 * ECC4[n_proj] / 32.0 + 45.0 * ECC6[n_proj] / 1024.0);
    TransverseMercator[n_proj].t_c3 = (15.0 * ECC4[n_proj] / 256.0 + 45.0 * ECC6[n_proj] / 1024.0);
    TransverseMercator[n_proj].t_c4 = (35.0 * ECC6[n_proj] / 3072.0);
    TransverseMercator[n_proj].t_ic1 = (1.5 * e1 - 27.0 * pow(e1, 3.0) / 32.0);
    TransverseMercator[n_proj].t_ic2 = (21.0 * e1 * e1 / 16.0 - 55.0 * pow(e1, 4.0) / 32.0);
    TransverseMercator[n_proj].t_ic3 = (151.0 * pow(e1, 3.0) / 96.0);
    TransverseMercator[n_proj].t_ic4 = (1097.0 * pow(e1, 4.0) / 512.0);
    TransverseMercator[n_proj].central_meridian = lon0;
    TransverseMercator[n_proj].use_false_easting = use_false_easting;

    // get y offset of central parallel (use lat0 = 0.0 for standard TM w/o offset)
    double lon, lat, x, y;
    lon = lon0;
    lat = lat0;
    tm(n_proj, lon, lat, &x, &y);
    TransverseMercator[n_proj].y_central_parralel = y;
}

void tm(n_proj, lon, lat, x, y)
int n_proj;
double lon, lat, *x, *y;
{
    /* Convert lon/lat to TM x/y */
    double N, T, T2, C, A, M, dlon, tan_lat, cos_lat, A2, A3, A5;

    dlon = lon - TransverseMercator[n_proj].central_meridian;
    if (fabs(dlon) > 360.0) dlon += copysign(360.0, -dlon);
    if (fabs(dlon) > 180.0) dlon = copysign(360.0 - fabs(dlon), -dlon);
    lat *= D2R;
    M = EQ_RAD[n_proj] * (TransverseMercator[n_proj].t_c1 * lat - TransverseMercator[n_proj].t_c2 * sin(2.0 * lat)
            + TransverseMercator[n_proj].t_c3 * sin(4.0 * lat) - TransverseMercator[n_proj].t_c4 * sin(6.0 * lat));
    if (fabs(lat) == M_PI_2) {
        *x = 0.0;
        *y = map_scale_factor * M;
    } else {
        N = EQ_RAD[n_proj] / d_sqrt(1.0 - ECC2[n_proj] * pow(sin(lat), 2.0));
        tan_lat = tan(lat);
        cos_lat = cos(lat);
        T = tan_lat * tan_lat;
        T2 = T * T;
        C = TransverseMercator[n_proj].t_e2 * cos_lat * cos_lat;
        A = dlon * D2R * cos_lat;
        A2 = A * A;
        A3 = A2 * A;
        A5 = A3 * A2;
        *x = map_scale_factor * N * (A + (1.0 - T + C) * (A3 * 0.16666666666666666667)
                + (5.0 - 18.0 * T + T2 + 72.0 * C - 58.0 * TransverseMercator[n_proj].t_e2) * (A5 * 0.00833333333333333333));
        A3 *= A;
        A5 *= A;
        *y = map_scale_factor * (M + N * tan(lat) * (0.5 * A2 + (5.0 - T + 9.0 * C + 4.0 * C * C) * (A3 * 0.04166666666666666667)
                + (61.0 - 58.0 * T + T2 + 600.0 * C - 330.0 * TransverseMercator[n_proj].t_e2) * (A5 * 0.00138888888888888889)));
    }

    // correct for TM x offset of central parallel
    *y -= TransverseMercator[n_proj].y_central_parralel;

    // correct for false easting
    if (TransverseMercator[n_proj].use_false_easting) {
        *x += 500000.0;
    }
}

void itm(n_proj, lon, lat, x, y)
int n_proj;
double *lon, *lat, x, y;
{

    // correct for TM x offset of central parallel
    y += TransverseMercator[n_proj].y_central_parralel;

    // correct for false easting
    if (TransverseMercator[n_proj].use_false_easting) {
        x -= 500000.0;
    }

    /* Convert TM x/y to lon/lat */
    double M, mu, phi1, C1, C12, T1, T12, tmp, tmp2, N1, R1, D, D2, D3, D5, cos_phi1, tan_phi1;

    M = y / map_scale_factor;
    mu = M / (EQ_RAD[n_proj] * TransverseMercator[n_proj].t_c1);
    phi1 = mu + TransverseMercator[n_proj].t_ic1 * sin(2.0 * mu) + TransverseMercator[n_proj].t_ic2 * sin(4.0 * mu)
            + TransverseMercator[n_proj].t_ic3 * sin(6.0 * mu) + TransverseMercator[n_proj].t_ic4 * sin(8.0 * mu);
    cos_phi1 = cos(phi1);
    tan_phi1 = tan(phi1);
    C1 = TransverseMercator[n_proj].t_e2 * cos_phi1 * cos_phi1;
    C12 = C1 * C1;
    T1 = tan_phi1 * tan_phi1;
    T12 = T1 * T1;
    tmp = 1.0 - ECC2[n_proj] * (1.0 - cos_phi1 * cos_phi1);
    tmp2 = d_sqrt(tmp);
    N1 = EQ_RAD[n_proj] / tmp2;
    R1 = EQ_RAD[n_proj] * (1.0 - ECC2[n_proj]) / (tmp * tmp2);
    D = x / (N1 * map_scale_factor);
    D2 = D * D;
    D3 = D2 * D;
    D5 = D3 * D2;

    *lon = TransverseMercator[n_proj].central_meridian + R2D * (D - (1.0 + 2.0 * T1 + C1) * (D3 * 0.16666666666666666667)
            + (5.0 - 2.0 * C1 + 28.0 * T1 - 3.0 * C12 + 8.0 * TransverseMercator[n_proj].t_e2 + 24.0 * T12)
            * (D5 * 0.00833333333333333333)) / cos(phi1);
    D3 *= D;
    D5 *= D;
    *lat = phi1 - (N1 * tan(phi1) / R1) * (0.5 * D2 -
            (5.0 + 3.0 * T1 + 10.0 * C1 - 4.0 * C12 - 9.0 * TransverseMercator[n_proj].t_e2) * (D3 * 0.04166666666666666667)
            + (61.0 + 90.0 * T1 + 298 * C1 + 45.0 * T12 - 252.0 * TransverseMercator[n_proj].t_e2 - 3.0 * C12) * (D5 * 0.00138888888888888889));
    (*lat) *= R2D;
}

/*
 *	TRANSFORMATION ROUTINES FOR THE Universal Transverse Mercator Projection (UTM)
 */

/*
int map_init_utm(int n_proj) {
    BOOLEAN search;
    double xmin, xmax, ymin, ymax, lon0;

    lon0 = 180.0 + 6.0 * TransverseMercator[n_proj].pars[0] - 3.0;
    if (lon0 >= 360.0) lon0 -= 360.0;
    vtm(lon0); // Central meridian for this zone
    if (TransverseMercator[n_proj].units_pr_degree) TransverseMercator[n_proj].pars[1] /= M_PR_DEG;
    TransverseMercator[n_proj].x_scale = TransverseMercator[n_proj].y_scale = TransverseMercator[n_proj].pars[1];
    forward = (PFI) utm;
    inverse = (PFI) iutm;

    if (fabs(TransverseMercator[n_proj].w - TransverseMercator[n_proj].e) > 360.0) { // -R in UTM meters
        iutm(&TransverseMercator[n_proj].w, &TransverseMercator[n_proj].s, TransverseMercator[n_proj].w, TransverseMercator[n_proj].s);
        iutm(&TransverseMercator[n_proj].e, &TransverseMercator[n_proj].n, TransverseMercator[n_proj].e, TransverseMercator[n_proj].n);
        TransverseMercator[n_proj].region = FALSE;
    }
    if (TransverseMercator[n_proj].region) {
        xy_search(&xmin, &xmax, &ymin, &ymax);
        outside = (PFI) wesn_outside;
        crossing = (PFI) wesn_crossing;
        overlap = (PFI) wesn_overlap;
        map_clip = (PFI) wesn_clip;
        left_edge = (PFD) left_rect;
        right_edge = (PFD) right_rect;
        search = FALSE;
    } else {
        utm(TransverseMercator[n_proj].w, TransverseMercator[n_proj].s, &xmin, &ymin);
        utm(TransverseMercator[n_proj].e, TransverseMercator[n_proj].n, &xmax, &ymax);
        outside = (PFI) rect_outside;
        crossing = (PFI) rect_crossing;
        overlap = (PFI) rect_overlap;
        map_clip = (PFI) rect_clip;
        left_edge = (PFD) left_rect;
        right_edge = (PFD) right_rect;
        frame_info.check_side = TRUE;
        search = TRUE;
    }

    frame_info.horizontal = TRUE;
    map_setxy(xmin, xmax, ymin, ymax);

    gmtdefs.basemap_type = 1;

    return (search);
}*/

void utm(n_proj, lon, lat, x, y)
int n_proj;
double lon, lat, *x, *y;
{
    /* Convert lon/lat to UTM x/y */

    if (lon < 0.0) lon += 360.0;
    tm(n_proj, lon, lat, x, y);
    (*x) += 500000.0;
    if (!TransverseMercator[n_proj].north_pole)
        (*y) += 10000000.0; /* For S hemisphere, add 10^6 m */
}

void iutm(n_proj, lon, lat, x, y)
int n_proj;
double *lon, *lat, x, y;
{
    /* Convert UTM x/y to lon/lat */

    x -= 500000.0;
    if (!TransverseMercator[n_proj].north_pole)
        y -= 10000000.0;
    itm(n_proj, lon, lat, x, y);
}

// utm init function
// created by ALomax to enable setting of north_pole flag

void vutm(int n_proj, double lon0, int lat0, int use_false_easting) {
    vtm(n_proj, lon0, lat0, use_false_easting);
    TransverseMercator[n_proj].north_pole = lat0 >= 0.0 ? 1 : 0;
}

/* Azimuthal Equidistant Projection (AE) */

struct AZIMUTHAL_EQUIDIST {
    BOOLEAN north_pole; /* TRUE if projection is on northern hemisphere, FALSE on southern */

    double central_meridian; // Central meridian (longitude) for projection
    double pole; // Central latitude for projection
    double sinp;
    double cosp;

};
struct AZIMUTHAL_EQUIDIST AzimuthalEquidistant[NUM_PROJ_MAX];


/*
 *	TRANSFORMATION ROUTINES FOR THE AZIMUTHAL EQUIDISTANT PROJECTION
 */

/*
int map_init_azeqdist () {
        BOOLEAN search;
        double xmin, xmax, ymin, ymax, dummy, radius;

        gmt_set_spherical ();	// PW: Force spherical for now

        if (project_info.units_pr_degree) {
                vazeqdist (0.0, 90.0);
                azeqdist (0.0, project_info.pars[3], &dummy, &radius);
                if (radius == 0.0) radius = M_PI * EQ_RAD;
                project_info.x_scale = project_info.y_scale = fabs (project_info.pars[2] / radius);
        }
        else
                project_info.x_scale = project_info.y_scale = project_info.pars[2];

        vazeqdist (project_info.pars[0], project_info.pars[1]);
        forward = (PFI)azeqdist;		inverse = (PFI)iazeqdist;

        if (fabs (project_info.pars[1]) == 90.0) {
                project_info.polar = TRUE;
                project_info.north_pole	= (project_info.pars[1] == 90.0);
        }

        if (!project_info.region) {	// Rectangular box given
                (*forward) (project_info.w, project_info.s, &xmin, &ymin);
                (*forward) (project_info.e, project_info.n, &xmax, &ymax);

                outside = (PFI) rect_outside;
                crossing = (PFI) rect_crossing;
                overlap = (PFI) rect_overlap;
                map_clip = (PFI) rect_clip;
                left_edge = (PFD) left_rect;
                right_edge = (PFD) right_rect;
                frame_info.check_side = !gmtdefs.oblique_anotation;
                frame_info.horizontal = (fabs (project_info.pars[1]) < 60.0 && fabs (project_info.n - project_info.s) < 30.0);
                search = TRUE;
        }
        else {
                if (project_info.polar && (project_info.n - project_info.s) < 180.0) {	// Polar aspect
                        if (!project_info.north_pole && project_info.s == -90.0) project_info.edge[0] = FALSE;
                        if (project_info.north_pole && project_info.n == 90.0) project_info.edge[2] = FALSE;
                        if ((fabs (project_info.w - project_info.e) == 360.0 || project_info.w == project_info.e)) project_info.edge[1] = project_info.edge[3] = FALSE;
                        outside = (PFI) polar_outside;
                        crossing = (PFI) wesn_crossing;
                        overlap = (PFI) wesn_overlap;
                        map_clip = (PFI) wesn_clip;
                        frame_info.horizontal = TRUE;
                        gmtdefs.n_lat_nodes = 2;
                        xy_search (&xmin, &xmax, &ymin, &ymax);
                }
                else {	// Global view only, force wesn = 0/360/-90/90
                        frame_info.anot_int[0] = frame_info.anot_int[1] = 0.0;		// No annotations for global mode
                        frame_info.frame_int[0] = frame_info.frame_int[1] = 0.0;	// No tickmarks for global mode
                        project_info.w = 0.0;
                        project_info.e = 360.0;
                        project_info.s = -90.0;
                        project_info.n = 90.0;
                        xmin = ymin = -M_PI * EQ_RAD;
                        xmax = ymax = -xmin;
                        outside = (PFI) eqdist_outside;
                        crossing = (PFI) eqdist_crossing;
                        overlap = (PFI) radial_overlap;
                        map_clip = (PFI) radial_clip;
                        gmtdefs.basemap_type = 1;
                }
                search = FALSE;
                left_edge = (PFD) left_circle;
                right_edge = (PFD) right_circle;
        }

        map_setinfo (xmin, xmax, ymin, ymax, project_info.pars[2]);
        project_info.r = 0.5 * project_info.xmax;

        return (search);
}*/

void vazeqdist(int n_proj, double lon0, double lat0) {
    /* Set up azimuthal equidistant projection */

    AzimuthalEquidistant[n_proj].central_meridian = lon0;
    AzimuthalEquidistant[n_proj].pole = lat0;
    AzimuthalEquidistant[n_proj].sinp = sin(lat0 * D2R);
    AzimuthalEquidistant[n_proj].cosp = cos(lat0 * D2R);
}

void azeqdist(int n_proj, double lon, double lat, double *x, double *y) {
    /* Convert lon/lat to azimuthal equidistant x/y */
    double k, dlon, cc, c, clat, clon, slat;

    while ((lon - AzimuthalEquidistant[n_proj].central_meridian) < -180.0) lon += 360.0;
    while ((lon - AzimuthalEquidistant[n_proj].central_meridian) > 180.0) lon -= 360.0;
    dlon = (lon - AzimuthalEquidistant[n_proj].central_meridian) * D2R;
    lat *= D2R;
    slat = sin(lat);
    clat = cos(lat);
    clon = cos(dlon);

    cc = AzimuthalEquidistant[n_proj].sinp * slat + AzimuthalEquidistant[n_proj].cosp * clat * clon;
    if (fabs(cc) >= 1.0)
        *x = *y = 0.0;
    else {
        c = acos(cc);
        k = EQ_RAD[n_proj] * c / sin(c);
        *x = k * clat * sin(dlon);
        *y = k * (AzimuthalEquidistant[n_proj].cosp * slat - AzimuthalEquidistant[n_proj].sinp * clat * clon);
    }
}

void iazeqdist(int n_proj, double *lon, double *lat, double x, double y) {
    /* Convert azimuthal equidistant x/yto lon/lat */
    double rho, c, sin_c, cos_c;

    rho = hypot(x, y);

    if (rho == 0.0) {
        *lat = AzimuthalEquidistant[n_proj].pole;
        *lon = AzimuthalEquidistant[n_proj].central_meridian;
    } else {
        c = rho / EQ_RAD[n_proj];
        sin_c = sin(c);
        cos_c = cos(c);
        *lat = asin(cos_c * AzimuthalEquidistant[n_proj].sinp + (y * sin_c * AzimuthalEquidistant[n_proj].cosp / rho)) * R2D;
        if (AzimuthalEquidistant[n_proj].pole == 90.0)
            *lon = AzimuthalEquidistant[n_proj].central_meridian + R2D * atan2(x, -y);
        else if (AzimuthalEquidistant[n_proj].pole == -90.0)
            *lon = AzimuthalEquidistant[n_proj].central_meridian + R2D * atan2(x, y);
        else
            *lon = AzimuthalEquidistant[n_proj].central_meridian +
                R2D * atan2(x * sin_c, (rho * AzimuthalEquidistant[n_proj].cosp * cos_c - y * AzimuthalEquidistant[n_proj].sinp * sin_c));
        if ((*lon) <= -180) (*lon) += 360.0;
    }
}
