
/*--------------------------------------------------------------------
    WARNING - this file is not gmt_map.c, it contains a subset of gmt_map.c

    subset of functions selected by A. Lomax, June 1998

    changes or additions indicated by                                 */

/*AJL ... AJL*/

/*AJL*/    /*END AJL*/

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

#define PI_2 (2.0 * M_PI)
#define D2R (M_PI / 180.0)
#define R2D (180.0/M_PI)

#define SMALL 1.0e-10

/*#define d_log(x) ((x) <= 0.0 ? gmt_NaN : log (x))*/
#define d_log(x) ((x) <= 0.0 ? -1.0e10 : log (x))
#define d_sqrt(x) ((x) < 0.0 ? 0.0 : sqrt (x))

/*double gmt_NaN;*/

typedef int BOOLEAN;              /* BOOLEAN used for logical variables */


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

#define N_ELLIPSOIDS 13

struct ELLIPSOID ellipse[N_ELLIPSOIDS] = {
		{ "WGS-84", 1984, 6378137.0, 6356752.1, 1.0/298.254 },
		{ "GRS-80", 1980, 6378137.0, 6356752.3, 1.0/298.257 },
		{ "WGS-72", 1972, 6378135.0, 6356750.5, 1.0/298.26 },
		{ "Australian", 1965, 6378160.0, 6356774.7, 1.0/298.25 },
		{ "Krasovsky", 1940, 6378245.0, 6356863.0, 1.0/298.3 },
		{ "International", 1924, 6378388.0, 6356911.9, 1.0/297.0 },
		{ "Hayford-1909", 1909, 6378388.0, 6356911.9, 1.0/297.0 },
		{ "Clarke-1880", 1880, 6378249.1, 6356514.9, 1.0/293.46 },
		{ "Clarke-1866", 1866, 6378206.4, 6356583.8, 1.0/294.98 },
		{ "Airy", 1830, 6377563.4, 6356256.9, 1.0/299.32 },
		{ "Bessel", 1841, 6377397.2, 6356079.0, 1.0/299.15 },
		{ "Hayford-1830", 1830, 6377276.3, 6356075.4, 1.0/300.80 },
		{ "Sphere", 1980, 6371008.7714, 6371008.7714, 0.0 }
};


// number of projections supported
#define NUM_PROJ_MAX 2

double EQ_RAD[NUM_PROJ_MAX];
double ECC[NUM_PROJ_MAX], ECC2[NUM_PROJ_MAX], ECC4[NUM_PROJ_MAX], ECC6[NUM_PROJ_MAX];
//double M_PR_DEG;

/* fields from struct MAP_PROJECTIONS taken from gmt_project.h,
	and converted to globals.
	WARNING - many fields removed! */

BOOLEAN NorthPole[NUM_PROJ_MAX];		/* TRUE if projection is on northern
					  hermisphere, FALSE on southern */
double CentralMeridian[NUM_PROJ_MAX];		/* Central meridian for projection */
double Pole[NUM_PROJ_MAX];			/* +90 pr -90, depending on which pole */


/* Lambert conformal conic parameters.
		(See Snyder for details on all parameters) */

double LambertConfConic_N[NUM_PROJ_MAX];
double LambertConfConic_F[NUM_PROJ_MAX];
double LambertConfConic_rho0[NUM_PROJ_MAX];




/* set constants that were set in GMT function map_setup */

/* use values from gmt_defaults.h: */

int map_setup_proxy (int n_proj, char* ellipsoid_name) {

	int num_ellipsoid;
	double f;

	/*mknan (gmt_NaN);*/

	/* determine ellipsoid */
	for (num_ellipsoid = 0;
		num_ellipsoid < N_ELLIPSOIDS
		&& strcmp (ellipsoid_name, ellipse[num_ellipsoid].name);
		num_ellipsoid++);
	if (num_ellipsoid == N_ELLIPSOIDS)
		return(-1);

	EQ_RAD[n_proj] = ellipse[num_ellipsoid].eq_radius;
	f = ellipse[num_ellipsoid].flattening;
	ECC2[n_proj] = 2 * f - f * f;
	ECC4[n_proj] = ECC2[n_proj] * ECC2[n_proj];
	ECC6[n_proj] = ECC2[n_proj] * ECC4[n_proj];
	ECC[n_proj] = d_sqrt (ECC2[n_proj]);

	return(0);

}

/*END AJL*/




/*
 *	TRANSFORMATION ROUTINES FOR THE LAMBERT CONFORMAL CONIC PROJECTION
 */

/*** function to set up a Lambert Conformal Conic projection */

int vlamb (n_proj, rlong0, rlat0, pha, phb)
int n_proj;
double rlong0, rlat0, pha, phb; {

	double t_pha, m_pha, t_phb, m_phb, t_rlat0;

	NorthPole[n_proj] = (rlat0 > 0.0);
        // AJL 20090812 BUG FIX
	Pole[n_proj] = (NorthPole[n_proj]) ? 90.0 : -90.0;
	//Pole[n_proj] = (NorthPole) ? 90.0 : -90.0;
	pha *= D2R;
	phb *= D2R;

	t_pha = tan (45.0 * D2R - 0.5 * pha) / pow ((1.0 - ECC[n_proj] *
		sin (pha)) / (1.0 + ECC[n_proj] * sin (pha)), 0.5 * ECC[n_proj]);
	m_pha = cos (pha) / d_sqrt (1.0 - ECC2[n_proj]
		* pow (sin (pha), 2.0));
	t_phb = tan (45.0 * D2R - 0.5 * phb) / pow ((1.0 - ECC[n_proj] *
		sin (phb)) / (1.0 + ECC[n_proj] * sin (phb)), 0.5 * ECC[n_proj]);
	m_phb = cos (phb) / d_sqrt (1.0 - ECC2[n_proj]
		* pow (sin (phb), 2.0));
	t_rlat0 = tan (45.0 * D2R - 0.5 * rlat0 * D2R) /
		pow ((1.0 - ECC[n_proj] * sin (rlat0 * D2R)) /
		(1.0 + ECC[n_proj] * sin (rlat0 * D2R)), 0.5 * ECC[n_proj]);

	if(pha != phb)
		LambertConfConic_N[n_proj] = (d_log(m_pha) - d_log(m_phb))
					/(d_log(t_pha) - d_log(t_phb));
	else
		LambertConfConic_N[n_proj] = sin(pha);

	LambertConfConic_F[n_proj] = m_pha/(LambertConfConic_N[n_proj] *
		pow(t_pha,LambertConfConic_N[n_proj]));
	CentralMeridian[n_proj] = rlong0;
	LambertConfConic_rho0[n_proj] = EQ_RAD[n_proj] * LambertConfConic_F[n_proj] *
		pow(t_rlat0,LambertConfConic_N[n_proj]);

	return(0);
}


/*** function to do x,y to lat,long Lambert Conformal Conic projection */

int lamb (n_proj, lon, lat, x, y)
int n_proj;
double lon, lat, *x, *y; {
	double rho,theta,hold1,hold2,hold3;

	while ((lon - CentralMeridian[n_proj]) < -180.0) lon += 360.0;
	while ((lon - CentralMeridian[n_proj]) > 180.0) lon -= 360.0;
	lat *= D2R;

	hold2 = pow (((1.0 - ECC[n_proj] * sin (lat))
		/ (1.0 + ECC[n_proj] * sin (lat))), 0.5 * ECC[n_proj]);
	hold3 = tan (45.0 * D2R - 0.5 * lat);
	if (fabs (hold3) < SMALL) hold3 = 0.0;
	hold1 = (hold3 == 0.0) ? 0.0 : pow (hold3 / hold2, LambertConfConic_N[n_proj]);
	rho = EQ_RAD[n_proj] * LambertConfConic_F[n_proj] * hold1;
	theta = LambertConfConic_N[n_proj] * (lon - CentralMeridian[n_proj]) * D2R;

	*x = rho * sin (theta);
	*y = LambertConfConic_rho0[n_proj] - rho * cos (theta);

	return(0);
}


/*** function to do lat,long to x,y inverse
			Lambert Conformal Conic projection */

int ilamb (n_proj, lon, lat, x, y)
int n_proj;
double *lon, *lat, x, y; {
	int i;
	double theta, temp, rho, t, tphi, phi, delta;

	theta = atan (x / (LambertConfConic_rho0[n_proj] - y));
	*lon = (theta /LambertConfConic_N[n_proj]) * R2D + CentralMeridian[n_proj];

	temp = x * x + (LambertConfConic_rho0[n_proj] - y)
		* (LambertConfConic_rho0[n_proj] - y);
	rho = copysign (d_sqrt (temp), LambertConfConic_N[n_proj]);
	t = pow ((rho / (EQ_RAD[n_proj] * LambertConfConic_F[n_proj])), 1./LambertConfConic_N[n_proj]);
	tphi = 90.0 * D2R - 2.0 * atan (t);
	delta = 1.0;
	for (i = 0; i < 100 && delta > 1.0e-8; i++) {
		temp = (1.0 - ECC[n_proj] * sin (tphi)) / (1.0 + ECC[n_proj] * sin (tphi));
		phi = 90.0 * D2R - 2.0 * atan (t * pow (temp, 0.5 * ECC[n_proj]));
		delta = fabs (fabs (tphi) - fabs (phi));
		tphi = phi;
	}
	*lat = phi * R2D;

	return(0);
}

