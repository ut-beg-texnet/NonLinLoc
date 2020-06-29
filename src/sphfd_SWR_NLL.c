/*
        NLLoc/src/sphfd_SWR_NLL.c

        Program to calculate FD times in a spherical coordinate system
        and output time grid in NonLinLoc 3D grid file format
        (http://alomax.free.fr/nlloc/soft7.00/formats.html#_grid_)

        Author:  Anthony Lomax, ALomax Scientific
        Version:  4.1  2019.04.05

        Heritage:  This code is an extension of Steve Roecker's sphfd.  The significant
        contributions of Roecker, Hole and Vidale to this code are gratefully acknowledged.
        See below for more details.

 */
/*
        sphfd.c

        Program to calculate FD times in a spherical coordinate system

        Author:  Steve Roecker, RPI
        Version:  4.0  July 2008

        Heritage:  This code is a rewrite of John Hole's cartesian FD code
        "punch.c" and several of the procedures and variable names are retained.
        punch.c was, in turn, derived from John Vidale's original
        FD code.  While this code is a significant rewrite to use spherical
        coordinates, the soundness of the basic framework of those previous
        codes made the writing of this code much easier.  The significant
        contributions of Hole and Vidale to this code are gratefully acknowledged.

        Compilation:  All functions and subroutines are included in this file, and
        a simple command like the following will suffice to compile:

        cc sphfd.c -o sphfd -lm

        Running the Code:  As with the Hole and Vidale codes, a parameter file
        is required.  Some parameters are required, others are optional.  The  file
        is read in using a routine written by Robert Clayton.  To specify a parameter
        file, enter an executable line like the following:

        sphfd par=my.parfile

        The meaning of the variables that the user can specify in the parameter file are
        discussed in the comments below.  Here is an example:

        fxs=74.49440
        fys=42.63330
        fzs=300.0
        nx=101
        ny=101
        nz=65
        x0=70.0
        y0=48.0
        z0=-35.00
        h=10.00
        timefile=AAK.ptimes
        velfile=VP.mod
        reverse=50
        NCUBE=2

        Note that some of these parameters MUST be specified, while others are optional.  In
        the example above, nx, ny, nz, x0, y0, z0, and h specify the grid and are always required.
        fxs, fys, fzs is the location of the startpoint and must be specified as well.   timefile
        and velfile are the principal output (travel time table) and input (wavespeeds) files and
        must be specified.  reverse and NCUBE are optional.

        Directional Conventions:

        This program works in a longitude, latitude, depth (radius) reference frame. Input
        coordinates for sources (e.g. in the parameter file) are geographic latitude (North
        positive), longitude (East positive), and depth (down positive).  Within the program,
        these are stored as:

                 x = longitude; positive east
                 y =  geocentric colatitude: positive south
                 z =  radius, positive away from center (i.e., up)

        Normally, the origin of the grid is the top, NW corner of the model, with elements increasing
        first in longitude (west to east), then in latitude (north to south) and finally in depth
        (top to bottom).

        The size of the cell is specified by

                h  =  the change in radius (nominally km, but need only be consistent with wavespeed units).
                dq = the change in geocentric latitude (degrees)
                df = the change in longitude (degrees)

        Note the actual dimensions of the cell will be (h x rdq x rsin(f)dq) where r is the
        local radius and f is the local latitude.  If dq and df are not specified, the program
        will set them so that the spaital dimensions are approsimately equal at the surface of the
        Earth (r = 6375 km) at the model mid-latitude.

        In this version we upgraded the geocentric conversion routine to take account of changes in elevation.
        Normally this is a very minor correction but we make it anyway so as not to worry about approximations.
        Note that we presume elevations to be specified with respect to the Reference Ellipsoid, not the usual
        Geoidal elevation (w.r.t. MSL).  To convert from one to the the other, you will need to know the local
        Geoidal elevation w.r.t the Reference Ellipsoid and adjust accordingly.

        Note that in the paramter file, fys and y0 are both Geographic/Geodetic latitudes, and we convert to
        gencentric interally.  fzs and z0 are Reference Ellipsoid elevations in km, but note that because these
        are treated as depths, the sign is positive down.  So, for example, a station at 200 meters elevation
        will have fzs = -0.2000.

        On output, the travel time table will have in its header the original x0, y0, z0 as the parameter file,
        BUT the origin (fxs, fys, fzs) in the geocentric frame used internally.  Hence, fxs will be longitude
        in radians, fys will be geocentric colatitude in radians, and fzs will be the distance from the center
        of the earth in km.

        BYTE SWAPPING

        This version allows for byte swapping so that travel time tables and
        wavespeed files can be read in from either little-endian or big-endian
        machines.  The program automatically detects if the platform is little- or
        big-endian.

        By default, the program expects that the input binary files will be in the native format:
        e.g., big-endian on platforms like on pre-Intel Macs or SUNs, little-endian on PCs or Linus-boxes.
        The default action is thus not to swap the binary files on input and ouput.

        The user can override this default using the optional parameter "swab".

        If swab = 0 (default), the program will not swap.  It will read and write these files in their
        native format without regard to machine type.   An example of where this would be
        used is if an Intel machine (i.e., little-endian) is being used and the user wishes
        all the binary files to be little-endian.

        If swab = 1, the program will swap the input files only.  The output file will be in the native order.

        If swab = 2, the program will swap the output files only.  The input file will be in the native order.

        If swab = 3, the program will swap both input and output files.

        The following table summarizes the use of "swab":

        Platform	Input	Output	Action		 Swab
        Big		Big	Big	No Swaps	  0
        Big		Big	Little	Swap Out Only	  2
        Big		Little	Big	Swap In  Only	  1
        Big		Little	Little	Swap In and Out	  3
        Little		Big	Big	Swap In and Out   3
        Little		Big	Little	Swap In Only 	  1
        Little		Little	Big	Swap Out Only 	  2
        Little		Little	Little	No Swaps	  0

 */

#include    <stdio.h>
// following needed for Linux?
#include    <string.h>
#include    <math.h>
#include    <fcntl.h>
/* file header structure */
// 20190405 AJL  #include  "vhead.h"
#include  "sphfd_SWR_NLL.h" // 20190405 AJL

#define PI  3.141592654
#define HPI 1.570796327
#define SQR2 1.414213562
#define SQR3 1.732050808
#define SQR6 2.449489743
#define rearth 6371.0
#define degrad PI/180.0
#define tc(x,y,z)   time0[nxy*(z) + nx*(y) + (x)]
#define sc(x,y,z)   slow0[nxy*(z) + nx*(y) + (x)]
// olddefine rc(z)   (rearth - (z0 + h*(z)))
// z0r is the distance from the center of the earth to the grid origin.   We take
// this to be the outer radius of the grid.
#define rc(z)   (z0r - h*(z))
#define qc(y)   (y0 + dq*(y))
#define fc(x)   (x0 + df*(x))
#define SQR(x) ((x) * (x))
#define DIST(x,y,z,x1,y1,z1) sqrt(SQR(x-(x1))+SQR(y-(y1)) + SQR(z-(z1)))

#define SWAP_2(x) ( (((x) & 0xff) << 8) | ((unsigned short)(x) >> 8) )
#define SWAP_4(x) ( ((x) << 24) | \
         (((x) << 8) & 0x00ff0000) | \
         (((x) >> 8) & 0x0000ff00) | \
         ((x) >> 24) )
#define FIX_SHORT(x) (*(unsigned short *)&(x) = SWAP_2(*(unsigned short *)&(x)))
#define FIX_INT(x)   (*(unsigned int *)&(x)   = SWAP_4(*(unsigned int *)&(x)))
#define FIX_FLOAT(x) FIX_INT(x)

struct sorted {
    float time;
    int i1, i2;
};

/* FUNCTION DECLARATIONS	*/
int
compar();
double
fdsph3d(), fdsphne(), fdsph2d(), fdsphnf(); /*STENCILS */
double
glat(), glath(), rcent, z0r;
void gp_close_dump(); // 20190405 ALomax

int endian();
int litend;

/* -------------------------------------------------------------------------- */

double glat(hlat)

double hlat;
/*
        Convert geographic latitude to geocentric latitude
        hlat (input) = geographic latitude in radians (north positive)
        glat (output)= geocentric latitude in radians (north positive)
 */
{
    double halfpi = 1.570796327;
    double x;

    if ((halfpi - fabs(hlat)) >= 0.05) {
        x = atan(0.993277 * sin(hlat) / cos(hlat));
    } else {
        /* Special formula near pole */
        if (hlat > 0)
            x = hlat / 0.993277 - 0.010632;
        else
            x = hlat / 0.993277 + 0.010632;
    }

    return (x);
}

/* -------------------------------------------------------------------------- */

double glath(hlat, h, r)

double hlat, h, *r;
/*
        Convert geographic latitude to geocentric latitude
       This is a version of glat that takes into account
       changes in elevation.

       This version uses the WGS84 ellipsoid

       Input
               hlat    Geodetic/Geogrphic latitude in radians (north positive)
               h       Ellipsoidal Height in kilometers.  Note this is the usual GPS height
                       rather than the often used Geoidal elevation (MSL).

       Output
               x       Geocentric latitude in radians (north positive) (on return)
               r       Distance from center of the Earth in kilometers

       All operations are double precision

       Author:  S. Roecker, RPI.   July, 2008

 */
{

    double x;

    /*	ap is semi major axis, bp is semiminor axis, f is inverse flattening, esq is the square of
            the ellipticity, which we compute from fi*(2.d0-fi) where fi is 1/f.
     */
    double ap = 6378137.0;
    double f = 298.257223563;
    //        double esq    =  6.69437978616733379e-03;
    double esq = 0.00669437978616733379;
    //        double degrad = 1.74532930056254081e-02;
    //        double degrad = 0.0174532930056254081;

    double anu;
    double sinxl;
    double xp, zp;

    // convert to meters
    h = h * 1000.0;

    sinxl = sin(hlat);

    //       anu is the ellipsoidal radius of curvature at the current geographic latitude
    anu = ap / sqrt(1.0 - esq * sinxl * sinxl);

    xp = (anu + h) * cos(hlat);
    zp = ((1.0 - esq) * anu + h) * sinxl;

    x = atan2(zp, xp);
    *r = sqrt(xp * xp + zp * zp) / 1000.0;

    return (x);
}

compar(a, b)
struct sorted *a, *b;
{
    if (a->time > b->time) return (1);
    if (b->time > a->time) return (-1);
    else return (0);
}

/* 3D TRANSMISSION STENCIL */
double fdsph3d(t, s, r, q, f, g, n, m, h, dq, df)

double t[], s[], r[], q[], f[];
double h, dq, df;
int g[], n[], m[];

/*	Given times, slowness, and related info for points 0-6 on
        a spherical element, solve for time at point 7.
        t   times at points; 0-6 known, t[7] unknown
        s   slownesses at points 0-7
        r   radii
        q   latitudes
        f   longitudes
      g,n,m		signs
        h   radial step
      dq  latitude step
      df  longitude step
 */
{
    double sqrt();

    double x, slo;
    double hsq, dqsq, dfsq, dfsinq7, df2sinq7, dfsinqi, df2sinqi;
    double a, b, c, c2;
    double rad, x1, x2;

    int i, j;
    int notx1, notx2;

    slo = 0.;
    for (i = 0; i < 8; i++) {
        slo += s[i];
    }
    slo = .125 * slo;

    hsq = h*h;
    dqsq = dq*dq;
    dfsq = df*df;
    dfsinq7 = df * sin(q[7]);
    df2sinq7 = df*dfsinq7;

    a = 1.0 / hsq + (1.0 / dqsq + 1. / (dfsinq7 * dfsinq7)) / (r[7] * r[7]);
    b = 0;
    for (i = 0; i < 7; i++) {
        b += t[i]*(g[7] * g[i] / hsq + (n[7] * n[i] / dqsq + m[7] * m[i] / (df2sinq7 * sin(q[i]))) / (r[7] * r[i]));
    }
    b = 2.0 * b;
    c = 0.;
    for (i = 0; i < 7; i++) {
        dfsinqi = df * sin(q[i]);
        c += t[i] * t[i]*(1. / hsq + (1. / dqsq + 1. / (dfsinqi * dfsinqi)) / (r[i] * r[i]));
    }
    for (i = 0; i < 6; i++) {
        c2 = 0.;
        df2sinqi = dfsq * sin(q[i]);
        for (j = i + 1; j < 7; j++) {
            c2 += t[j]*(g[i] * g[j] / hsq + (n[i] * n[j] / dqsq + m[i] * m[j] / (df2sinqi * sin(q[j]))) / (r[i] * r[j]));
        }
        c += 2.0 * c2 * t[i];
    }
    c = c - 16.0 * slo*slo;

    rad = b * b - 4.0 * a*c;
    notx1 = 0;
    notx2 = 0;
    if (rad >= 0.) {
        x1 = (-b + sqrt(rad)) / (2.0 * a);
        if ((x1 < t[0]) || (x1 < t[1]) || (x1 < t[2]) || (x1 < t[3]) || (x1 < t[4]) || (x1 < t[5]) || (x1 < t[6]))
            notx1 = 1;
        x2 = (-b - sqrt(rad)) / (2.0 * a);
        if ((x2 < t[0]) || (x2 < t[1]) || (x2 < t[2]) || (x2 < t[3]) || (x2 < t[4]) || (x2 < t[5]) || (x2 < t[6]))
            notx2 = 1;

        if (notx1) {
            if (notx2) {
                x = 1.e11; /* ACAUSAL; ABORT */
            } else {
                x = x2;
            }
        } else {
            if (notx2) {
                x = x1;
            } else {
                x = x2;
                if (x1 < x2) x = x1;
            }
        }
    } else {
        x = 1.e11; /* SQRT IMAGINARY; ABORT */
    }

    return (x);
}

/* 3D STENCIL FOR NEW EDGE */
double fdsphne(t, s, d01, d12, d14, d25, d45)

double t[], s[];
double d01, d12, d14, d25, d45;

/*  Given times at points 0-4, find time at point 5 using the
    edge stencil */
{
    double x, slo;
    double d1, d2, dt1, dt2, dt3, d1sq;
    double a, b, c;
    double x1, x2, rad;

    double sqrt();

    int notx1, notx2;

    slo = .25 * (s[1] + s[2] + s[4] + s[5]);

    dt1 = (t[4] - t[1]) / d14 - t[2] / d25;
    dt2 = (t[2] - t[1]) / d12 - t[4] / d45;
    dt3 = (t[0] - t[3]) / d01;

    a = 1. / (d25 * d25) + 1. / (d45 * d45);
    b = 2 * (dt1 / d25 + dt2 / d45);
    c = dt1 * dt1 + dt2 * dt2 + dt3 * dt3 - 4 * slo*slo;

    rad = b * b - 4.0 * a*c;
    notx1 = 0;
    notx2 = 0;
    if (rad >= 0.) {
        x1 = (-b + sqrt(rad)) / (2 * a);
        if ((x1 < t[0]) || (x1 < t[1]) || (x1 < t[2]) || (x1 < t[3]) || (x1 < t[4])) notx1 = 1;
        x2 = (-b - sqrt(rad)) / (2 * a);
        if ((x2 < t[0]) || (x2 < t[1]) || (x2 < t[2]) || (x2 < t[3]) || (x2 < t[4])) notx2 = 1;
        if (notx1) {
            if (notx2) {
                x = 1.e11; /* ACAUSAL; ABORT */
            } else {
                x = x2;
            }
        } else {
            if (notx2) {
                x = x1;
            } else {
                x = x1;
                if (x2 < x1) x = x2;
            }
        }
    } else {
        x = 1.e11; /* SQRT IMAGINARY; ABORT */
    }
    return (x);
}

/* 2D TRANSMISSION STENCIL (FOR HEAD WAVES ON FACES OF GRID CELLS) */
double fdsph2d(t, s, d01, d02, d13, d23)
double t[], s[];
double d01, d02, d13, d23;
/* s[0] at newpoint; s[1],t[1] & s[2],t[2] adjacent; s[3],t[3] diagonal
     d01	distance between points 0 and 1
     d02	distance between points 0 and 2
     d13	distance between points 1 and 3
 */
{
    double x, slo;
    double d01sq, d02sq;
    double t1mt3, t2mt3;
    double t1t3t2, t2t3t1;
    double a, b, c;
    double x1, x2, rad;

    int notx1, notx2;

    double sqrt();

    slo = .25 * (s[0] + s[1] + s[2] + s[3]);

    d01sq = d01*d01;
    d02sq = d02*d02;
    t1mt3 = t[1] - t[3];
    t2mt3 = t[2] - t[3];
    t1t3t2 = t1mt3 / d13 - t[2] / d02;
    t2t3t1 = t2mt3 / d23 - t[1] / d01;

    a = (1.0 / d01sq + 1. / d02sq);
    b = 2.0 * (t1t3t2 / d02 + t2t3t1 / d01);
    c = t1t3t2 * t1t3t2 + t2t3t1 * t2t3t1 - 4 * slo*slo;

    rad = b * b - 4.0 * a*c;

    notx1 = 0;
    notx2 = 0;
    if (rad >= 0.) {
        x1 = (-b + sqrt(rad)) / (2 * a);
        if ((x1 < t[1]) || (x1 < t[2]) || (x1 < t[3]) || (x1 < t[4]) || (x1 < t[0])) notx1 = 1;
        x2 = (-b - sqrt(rad)) / (2 * a);
        if ((x2 < t[1]) || (x2 < t[2]) || (x2 < t[3]) || (x2 < t[4]) || (x2 < t[0])) notx2 = 1;
        if (notx1) {
            if (notx2) {
                x = 1.e11; /* ACAUSAL; ABORT */
            } else {
                x = x2;
            }
        } else {
            if (notx2) {
                x = x1;
            } else {
                x = x1;
                if (x2 < x1) x = x2;
            }
        }
    } else {
        x = 1.e11; /* SQRT IMAGINARY; ABORT */
    }

    return (x);
}

/* 3D STENCIL FOR NEW FACE */
double fdsphnf(t, s, d02, d12, d25)
double t[], s[];
double d02, d12, d25;
/* s5 at newpoint (t5 or x); s2,t2 adjacent on old face
     t0,t4 beside t2 on old face and opposite each other;
     t1,t3 beside t1 on old face and opposite each other:

                     4
                     |
                1---2/5---3
                     |
                     0

        Side	dr	dq	df
        Top	5-2	0-2	2-1
        North	2-0	5-2	1-2
        West	2-1	0-2	5-2

 */
{
    double x, slo;
    double t4t0, t1t3;
    double sqrt();

    slo = .5 * (s[5] + s[2]);

    t4t0 = (t[4] - t[0]) / d02;
    t1t3 = (t[1] - t[3]) / d12;
    x = slo * slo - 0.25 * (t4t0 * t4t0 + t1t3 * t1t3);

    if (x >= 0.) {
        x = t[2] + d25 * sqrt(x);
        if ((x < t[0]) || (x < t[1]) || (x < t[3]) || (x < t[4])) /* ACAUSAL; ABORT */
            x = 1.e11;
    } else x = 1.e11; /* SQRT IMAGINARY; ABORT */

    return (x);
}


/* ====================================================================
   ==================================================================== */

/* copyright (c) Robert W. Clayton
 *		 Seismological Laboratory
 *		 Caltech
 *		 Pasadena, CA 91125
 *
 * Getpar routines:
 *
 * Externally visable routines:
 *
 *		setpar(argc,argv)
 *		getpar(name,type,valptr)
 *		mstpar(name,type,valptr)
 *		endpar()
 *
 * To get C-version:
 *		cc -c getpar.c
 *
 * To get F77-version:
 *		cp getpar.c fgetpar.c
 *		cc -c -DFORTRAN fgetpar.c
 *		rm fgetpar.c
 *
 * To get the environment processing stuff add the flag
 *-DENVIRONMENT to each of the cc's above.
 */
#include <stdio.h>

#define MAXLINE  1024 /* max length of line in par file */
#define MAXNAME  64 /* max length of name */
#define MAXVALUE 1024 /* max length of value */
#define MAXFILENAME 64 /* max length of par file name */
#define MAXVECTOR 10 /* max # of elements for unspecified vectors */
#define GETPAR_ERROR 100 /* exit status for getpar error */
#define GETPAR_STOP 101 /* exit status for STOP or mstpar */
#define MAXPARLEVEL 4 /* max recurrsion level for par files */

#ifdef FORTRAN
#define GETPAR getpar_
#define MSTPAR mstpar_
#define ENDPAR endpar_
#else
#define GETPAR getpar
#define MSTPAR mstpar
#define ENDPAR endpar
#endif

#define INIT  1 /* bits for FLAGS (ext_par.argflags) */
#define STOP  2
#define LIST  4
#define END_PAR  8
#define VERBOSE 16

#define LISTINC  32 /* increment size for arglist */
#define BUFINC  1024 /* increment size for argbuf */

struct arglist /* structure of list set up by setpar */ {
    char *argname;
    char *argval;
    int hash;
};

struct ext_par /* global variables for getpar */ {
    char *progname;
    int argflags;
    struct arglist *arglist;
    struct arglist *arghead;
    char *argbuf;
    int nlist;
    int nbuf;
    int listmax;
    int bufmax;
    FILE *listout;
} ext_par;

/* abbreviations: */
#define AL   struct arglist
#define PROGNAME ext_par.progname
#define FLAGS  ext_par.argflags
#define ARGLIST  ext_par.arglist
#define ARGHEAD  ext_par.arghead
#define ARGBUF  ext_par.argbuf
#define NLIST  ext_par.nlist
#define NBUF  ext_par.nbuf
#define LISTMAX  ext_par.listmax
#define BUFMAX  ext_par.bufmax
#define LISTFILE ext_par.listout

#ifdef FORTRAN
setpar_()
#else

setpar(ac, av) /* set up arglist & process INPUT command */
int ac;
char **av;
#endif
{
    register char *pl, *pn, *pv;
    char t, name[MAXNAME], value[MAXVALUE];
    FILE *file, *gp_create_dump();
    int i, addflags, nevlist, testav, testae;
    struct arglist *alptr;
#ifdef FORTRAN
    int ac;
    char **av;
    extern int xargc;
    extern char **xargv;
    ac = xargc;
    av = xargv;
#endif

    PROGNAME = *av;
    FLAGS = INIT;
    LISTFILE = stderr;

    ARGLIST = NULL;
    ARGBUF = NULL;
    NLIST = NBUF = LISTMAX = BUFMAX = 0;
#ifdef ENVIRONMENT
    gp_do_environment(ac, av);
#endif
    nevlist = NLIST;

    while (--ac > 0) {
        av++;
        pl = *av;
        while (*pl == ' ' || *pl == '\t') pl++;
        /* get name */
        pn = name;
        while (*pl != '=' && *pl != '\0') *pn++ = *pl++;
        *pn++ = '\0';
        /* get value */
        if (*pl == '=') pl++;
        pv = value;
        if (*pl == '"' || *pl == '\'') {
            t = *pl++;
            while (*pl != '\0') {
                if (*pl == t) {
                    if (pl[-1] != '\\') break;
                    pv[-1] = t;
                    pl++;
                } else *pv++ = *pl++;
            }
        } else while (*pl) *pv++ = *pl++;
        *pv = '\0';
        if (name[0] == '-') gp_add_entry("SWITCH", &name[1]);
        else gp_add_entry(name, value);
        if (strcmp("par", name) == 0) /* par file */
            gp_do_par_file(value, 1);
    }

    /* do not internally call getpar before this point because
       ARGHEAD is not set. The search will have no stopping point */
    ARGHEAD = ARGLIST;
#ifdef ENVIRONMENT
    *value = '\0';
    if (GETPAR("NOENV", "s", value)) ARGHEAD = ARGLIST + nevlist;
#endif
    addflags = 0;
    *value = '\0';
    if (GETPAR("STOP", "s", value)) addflags |= STOP;
    *value = '\0';
    if (GETPAR("VERBOSE", "s", value)) addflags |= VERBOSE;
    *value = '\0';
    if (GETPAR("LIST", "s", value)) {
        addflags |= LIST;
        LISTFILE = gp_create_dump(value, "list");
    }
    *value = '\0';
    if (GETPAR("INPUT", "s", value)) {
        file = gp_create_dump(value, "list input");
        fprintf(file, "%s: getpar input listing\n", PROGNAME);
        for (i = 0, alptr = ARGLIST; i < NLIST; i++, alptr++) {
            fprintf(file, "%3d: %16s = %s\n",
                    i, alptr->argname, alptr->argval);
        }
        gp_close_dump(file);
    }
    FLAGS |= addflags;
}

gp_add_entry(name, value) /* add an entry to arglist, expanding memory */
register char *name, *value; /* if necessary */
{
    struct arglist *alptr;
    int len;
    register char *ptr;

    /* check arglist memory */
    if (NLIST >= LISTMAX) {
        LISTMAX += LISTINC;
        if (ARGLIST == NULL)
            ARGLIST = (AL *) malloc(LISTMAX * sizeof (AL));
        else ARGLIST = (AL *) realloc(ARGLIST, LISTMAX * sizeof (AL));
    }
    /* check argbuf memory */
    len = strlen(name) + strlen(value) + 2; /* +2 for terminating nulls */
    if (NBUF + len >= BUFMAX) {
        BUFMAX += BUFINC;
        if (ARGBUF == NULL)
            ARGBUF = (char *) malloc(BUFMAX);
        else ARGBUF = (char *) realloc(ARGBUF, BUFMAX);
    }
    if (ARGBUF == NULL || ARGLIST == NULL)
        gp_getpar_err("setpar", "cannot allocate memory");

    /* add name */
    alptr = ARGLIST + NLIST;
    alptr->hash = gp_compute_hash(name);
    ptr = alptr->argname = ARGBUF + NBUF;
    do *ptr++ = *name; while (*name++);

    /* add value */
    NBUF += len;
    alptr->argval = ptr;
    do *ptr++ = *value; while (*value++);
    NLIST++;
}
#ifdef ENVIRONMENT

gp_do_environment(ac, av)
int ac;
char **av;
{
    char **ae;
    register char *pl, *pn, *pv;
    char name[MAXNAME], value[MAXVALUE], t;

    /* The environ pointer ae, is assumed to have a specific relation
       to the arg pointer av. This may not be portable. */
    ae = av + (ac + 1);
    if (ae == NULL) return;

    while (*ae != NULL) {
        pl = *ae++;
        while (*pl == ' ' || *pl == '\t') pl++;
        /* get name */
        pn = name;
        while (*pl != '=' && *pl != '\0') *pn++ = *pl++;
        *pn = '\0';
        if (strcmp("NOENV", pn) == 0) return;

        /* get value */
        if (*pl == '=') pl++;
        pv = value;
        if (*pl == '"' || *pl == '\'') {
            t = *pl++;
            while (*pl != '\0') {
                if (*pl == t) {
                    if (pl[-1] != '\\') break;
                    pv[-1] = t;
                    pl++;
                } else *pv++ = *pl++;
            }
        } else while (*pl) *pv++ = *pl++;
        *pv = '\0';
        gp_add_entry(name, value);
    }
}
#endif

ENDPAR() /* free arglist & argbuf memory, & process STOP command */ {
    if (ARGLIST != NULL) free(ARGLIST);
    if (ARGBUF != NULL) free(ARGBUF);
    ARGBUF = NULL;
    ARGLIST = NULL;
    if (FLAGS & STOP) {
        fprintf(stderr, "%s[endpar]: stop due to STOP in input\n",
                PROGNAME);
        exit(GETPAR_STOP);
    }
    FLAGS = END_PAR; /* this stops further getpar calls */
}

#ifdef FORTRAN
mstpar_(name, type, val, dum1, dum2)
int dum1, dum2; /* dum1 & dum2 are extra args that fortran puts in */
#else

mstpar(name, type, val)
#endif
char *name, *type;
int *val;
{
    int cnt;
    char *typemess;

    if ((cnt = GETPAR(name, type, val)) > 0) return (cnt);

    /* The following line corrects a common input error */
    if (type[1] == 'v') {
        type[1] = type[0];
        type[0] = 'v';
    }

    switch (*type) {
        case 'd': typemess = "an int";
            break;
        case 'f': typemess = "a float";
            break;
        case 'F': typemess = "a double";
            break;
        case 's': typemess = "a string";
            break;
        case 'b': typemess = "a boolean";
            break;
        case 'v': switch (type[1]) {
                case 'd': typemess = "an integer vector";
                    break;
                case 'f': typemess = "a float vector";
                    break;
                case 'F': typemess = "a double vector";
                    break;
                default: typemess = "unknow vectorn (error)";
                    break;
            }
            break;
        default: typemess = "unknown (error)";
            break;
    }
    gp_getpar_err("mstpar", "must specify value for '%s', expecting %s",
            name, typemess);
}

#ifdef FORTRAN
getpar_(name, type, val, dum1, dum2)
int dum1, dum2; /* dum1 & dum2 are extra args that fortran puts in */
#else

getpar(name, type, val)
#endif
char *name, *type;
int *val;
{
    register char *sptr;
    register struct arglist *alptr;
    register int i;
    double atof(), *dbl;
    float *flt;
    int h, hno, hyes, found;
    char line[MAXLINE], *str, *noname;

    if (FLAGS & END_PAR)
        gp_getpar_err("getpar", "called after endpar");
    if ((FLAGS & INIT) == 0)
        gp_getpar_err("getpar", "not initialized with setpar");
    if (FLAGS & VERBOSE)
        fprintf(stderr, "getpar: looking for %s\n", name);

    /* The following line corrects a common input error */
    if (type[1] == 'v') {
        type[1] = type[0];
        type[0] = 'v';
    }


    if (*type == 'b') goto boolean;

    h = gp_compute_hash(name);
    found = 0;
    /* search list backwards, stopping at first find */
    for (alptr = ARGLIST + (NLIST - 1); alptr >= ARGHEAD; alptr--) {
        if (alptr->hash != h) continue;
        if (strcmp(alptr->argname, name) != 0) continue;
        str = alptr->argval;
        switch (*type) {
            case 'd':
                *val = atoi(str);
                found = 1;
                break;
            case 'f':
                flt = (float *) val;
                *flt = atof(str);
                found = 1;
                break;
            case 'F':
                dbl = (double *) val;
                *dbl = atof(str);
                found = 1;
                break;
            case 's':
                sptr = (char *) val;
                while (*str) *sptr++ = *str++;
                *sptr = '\0';
                found = 1;
                break;
            case 'v':
                found = gp_getvector(str, type, val);
                break;
            default:
                gp_getpar_err("getpar",
                        "unknown conversion type %s", type);
                break;
        }
        break;
    }
    goto list;
boolean:
    noname = line;
    sprintf(noname, "no%s", name);
    hno = gp_compute_hash(noname);
    hyes = gp_compute_hash(name);
    found = 0;
    /* search list backwards, stopping at first find */
    for (alptr = ARGLIST + (NLIST - 1); alptr >= ARGHEAD; alptr--) {
        if (alptr->hash != hno && alptr->hash != hyes) continue;
        if (strcmp(alptr->argname, name) == 0) {
            if (alptr->argval[0] == '\0') *val = 1;
            else *val = atol(alptr->argval);
            found++;
            break;
        }
        if (strcmp(alptr->argname, noname) == 0) {
            *val = 0;
            found++;
            break;
        }
    }
list:
    if (FLAGS & LIST) {
        switch (*type) {
            case 'd': sprintf(line, "(int) = %d", *val);
                break;
            case 'f': flt = (float *) val;
                sprintf(line, "(flt) = %14.6e", *flt);
                break;
            case 'F': dbl = (double *) val;
                sprintf(line, "(dbl) = %14.6e", *dbl);
                break;
            case 's': sprintf(line, "(str) = %s", val);
                break;
            case 'b': sprintf(line, "(boo) = %d", *val);
                break;
            case 'v': switch (type[1]) {
                        /* should list these out */
                    case 'd': sprintf(line, "(int vec)");
                        break;
                    case 'f': sprintf(line, "(flt vec)");
                        break;
                    case 'F': sprintf(line, "(dbl vec)");
                        break;
                    default: sprintf(line, " vec type error");
                        break;
                }
                break;
            default: sprintf(line, " type error");
                break;
        }
        fprintf(LISTFILE, "%16s (%s) %s \n", name,
                (found ? "set" : "def"), line);
    }
    return (found);
}

FILE *gp_create_dump(fname, filetype)
char *fname;
char *filetype;
{
    FILE *temp;

    if (*fname == '\0') return (stderr);
    if (strcmp(fname, "stderr") == 0) return (stderr);
    if (strcmp(fname, "stdout") == 0) return (stdout);
    if ((temp = fopen(fname, "w")) != NULL) return (temp);
    fprintf(stderr, "%s[setpar]: cannot create %s file %s\n",
            PROGNAME, filetype, fname);
    return (stderr);
}

//gp_close_dump(file)		// 20190405 ALomax

void gp_close_dump(file) // 20190405 ALomax
FILE *file;
{
    if (file == stderr || file == stdout) return;
    fclose(file);
}

gp_compute_hash(s)
register char *s;
{
    register int h;
    h = s[0];
    if (s[1]) h |= (s[1]) << 8;
    else return (h);
    if (s[2]) h |= (s[2]) << 16;
    else return (h);
    if (s[3]) h |= (s[3]) << 24;
    return (h);
}

gp_do_par_file(fname, level)
char *fname;
int level;
{
    register char *pl, *pn, *pv;
    char t1, t2, line[MAXLINE], name[MAXNAME], value[MAXVALUE];
    FILE *file, *fopen();

    if (level > MAXPARLEVEL)
        gp_getpar_err("setpar", "%d (too many) recursive par file", level);

    if ((file = fopen(fname, "r")) == NULL)
        gp_getpar_err("setpar", "cannot open par file %s", fname);

    while (fgets(line, MAXLINE, file) != NULL) {
        pl = line;
        /* loop over entries on each line */
loop:
        while (*pl == ' ' || *pl == '\t') pl++;
        if (*pl == '\0' || *pl == '\n') continue;
        if (*pl == '#') continue; /* comments on rest of line */

        /* get name */
        pn = name;
        while (*pl != '=' && *pl != '\0' && *pl != ' '
                && *pl != '\t') *pn++ = *pl++;
        *pn = '\0';
        if (*pl == '=') pl++;

        /* get value */
        *value = '\0';
        pv = value;
        if (*pl == '"' || *pl == '\'') {
            t1 = t2 = *pl++;
        } else {
            t1 = ' ';
            t2 = '\t';
        }
        while (*pl != t1 && *pl != t2 &&
                *pl != '\0' && *pl != '\n') *pv++ = *pl++;
        *pv = '\0';
        if (*pl == '"' || *pl == '\'') pl++;
        gp_add_entry(name, value);
        if (strcmp("par", name) == 0)
            gp_do_par_file(value, level + 1);
        goto loop;
    }
    fclose(file);
}

gp_getpar_err(subname, mess, a1, a2, a3, a4)
char *subname, *mess;
int a1, a2, a3, a4;
{
    fprintf(stderr, "\n***** ERROR in %s[%s] *****\n\t",
            (PROGNAME == NULL ? "(unknown)" : PROGNAME), subname);
    fprintf(stderr, mess, a1, a2, a3, a4);
    fprintf(stderr, "\n");
    exit(GETPAR_ERROR);
}

gp_getvector(list, type, val)
char *list, *type;
int *val;
{
    register char *p;
    register int index, cnt;
    char *valptr;
    int limit;
    int ival, *iptr;
    float fval, *fptr;
    double dval, *dptr, atof();

    limit = MAXVECTOR;
    if (type[2] == '(' || type[2] == '[') limit = atol(&type[3]);
    if (limit <= 0)
        gp_getpar_err("getpar", "bad limit=%d specified", limit);
    index = 0;
    p = list;
    while (*p != '\0' && index < limit) {
        cnt = 1;
backup: /* return to here if we find a repetition factor */
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '\0') return (index);
        valptr = p;
        while (*p != ',' && *p != '*' && *p != 'x' && *p != 'X' &&
                *p != '\0') p++;
        if (*p == '*' || *p == 'x' || *p == 'X') {
            cnt = atol(valptr);
            if (cnt <= 0)
                gp_getpar_err("getpar",
                    "bad repetition factor=%d specified",
                    cnt);
            if (index + cnt > limit) cnt = limit - index;
            p++;
            goto backup;
        }
        switch (type[1]) {
            case 'd':
                iptr = (int *) val;
                ival = atol(valptr);
                while (cnt--) iptr[index++] = ival;
                break;
            case 'f':
                fptr = (float *) val;
                fval = atof(valptr);
                while (cnt--) fptr[index++] = fval;
                break;
            case 'F':
                dptr = (double *) val;
                dval = atof(valptr);
                while (cnt--) dptr[index++] = dval;
                break;
            default:
                gp_getpar_err("getpar",
                        "bad vector type=%c specified", type[1]);
                break;
        }
        if (*p != '\0') p++;
    }
    return (index);
}

/* Routine to determine if the machine is big or little endian */
int endian() {
    int *T;

    T = (int *) "\01\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0";
    if (*T == 1) {
        printf("This machine is little-endian.\n");
        return (1);
    } else {
        printf("This machine is big-endian.\n");
        return (0);
    }

}

/*
   Swap bytes for a double precision number
 */
int FixDouble(double *n) {
    unsigned char *cptr, tmp;

    cptr = (unsigned char *) n;
    tmp = cptr[0];
    cptr[0] = cptr[7];
    cptr[7] = tmp;
    tmp = cptr[1];
    cptr[1] = cptr[6];
    cptr[6] = tmp;
    tmp = cptr[2];
    cptr[2] = cptr[5];
    cptr[5] = tmp;
    tmp = cptr[3];
    cptr[3] = cptr[4];
    cptr[4] = tmp;

    return (1);
}

main(ac, av)
int ac;
char **av;
{
    /* NOTE THAT SEVERAL VARIABLES MUST BE SPECIFIED IN par=xxx FILE,
       WHILE OTHERS ARE OPTIONAL:  IF A mstpar STATEMENT READS THE
       VARIABLE BELOW, THEN THE VARIABLE IS REQUIRED;  IF A getpar
       STATEMENT IS USED BELOW, THEN THE VARIABLE IS OPTIONAL */
    int
    nx, /* x-dimension of mesh (LEFT-TO-RIGHT = West to East) */
            ny, /* y-dimension of mesh (FRONT-TO-BACK = North to South) */
            nz, /* z-dimension of mesh  (TOP-TO-BOTTOM) */
            iplus = 1, /* rate of expansion of "cell" in the */
            iminus = 1, /*    plus/minus x/y/z direction */
            jplus = 1,
            jminus = 1,
            kplus = 1,
            kminus = 1,
            igrow, /* counter for "cell" growth */
            srctype = 1, /* if 1, source is a point;
	                              2, source is on the walls of the data volume;
				      3, source on wall, time field known; */
            floatsrc = 1, /* if 0, source must be on a grid node
				      1, source can lie between grid nodes */
            srcwall, /* if 1, source on x=0 wall,  if 2, on x=nx-1 wall
	                           if 3, source on y=0 wall,  if 4, on y=ny-1 wall
	                           if 5, source on z=0 wall,  if 6, on z=nz-1 wall */
            xs, /* shot x position (in grid points) */
            ys, /* shot y position */
            zs, /* shot depth */
            xx, yy, zz, /* Used to loop around xs, ys, zs coordinates	*/
            X1, X2, index, ii, i, j, k, radius,
            tfint, vfint, wfint, ofint, bfint,
            nxy, nyz, nxz, nxyz, nwall, nbox,
            nxp, nyp, nzp,
            /* counters for the position of the sides of current cell */
            x1, x2, y1, y2, z1, z2,
            /* flags set to 1 until a side has been reached */
            dx1 = 1, dx2 = 1, dy1 = 1, dy2 = 1, dz1 = 1, dz2 = 1, rad0 = 1,
            maxrad, /* maximum radius to compute */
            xxx, yyy, zzz, /* used in linear velocity gradient cell source */
            NCUBE = 2, /* size of source cell     1 for 3x3x3, 2 for 5x5x5...	*/
            parse = 1, /* parse the output by this amount */
            fill = 0, /* =0 to fill in starting box with constant wavespeed times  */
            /* =1 to fill in starting box with precomputed times from a file */
            invh = 0,
            invq = 0,
            invf = 0, /* If invh, invq or invf are =1, then the values of h, dq, and df are inverted */
            invx0 = 0,
            invy0 = 0,
            invz0 = 0, /* If invx0, invy0 or invz0 are =1, then the values of h, dq, and df are inverted */
            swab = 0, /* If swab = 0, the program will read and write these files in their
                             native format without regard to machine type. */
            savsrc = 0, /* When we cascade, at some stages we may reference an integer source
                               but we want to preserve the true source location for future refernce.
                               Set savsrc to 1 to preserve the true source in the header, and to  0
                               to redefine it as a grid point */

            reverse = 1, /* will automatically do up to this number of
			           reverse propagation steps to fix waves that travel
			           back into expanding cell */
            headpref = 6, /* if headpref starts > 0, will determine
			           model wall closest to source and will prefer to start
			           reverse calculations on opposite wall */
            head, headw[7]; /* counters for detecting head waves on sides of current cell */

    int
    g[8], /* Stores radial signs for stencils */
            n[8], /* Stores latitude signs for 3D stencils */
            m[8]; /* Stores longitude signs for 3D stencils */

    double
    r[8], /* Stores radial values for 3D stencils */
            q[8], /* Stores latitude values for 3D stencils */
            f[8], /* Stores longitude values for 3D stencils */
            t[8], /* Stores times for stencils */
            s[8]; /* Stores wavespeeds for stencils */

    double d01, d02, d12, d13, d23, d14, d25, d45; /* Distances between stencil points */

    double rz1, rz2, rX2;
    double rzp1, rzm1;
    double rXp1, rXm1;
    double qy1, qy1p1;
    double qy2, qy2m1;
    double qX1, qX1p1, qX1m1;
    double qX2, qX2p1, qX2m1;
    double fx1, fx1p1, fx1m1;
    double fx2, fx2m1;
    double sinqX2, sinqX1;
    double sinqy2, sinqy1;
    double x1p, y1p, z1p;
    double x2p, y2p, z2p;
    double fy1, fz1;

    /* test variables */
    double delt, ddelt, timec, dtc, rtol = 0.1;
    double T0, T1, T2, T3, T4, T5, T6, T7;
    double dtdr, dtdq, dtdf, snew;

    double
    h, /* radial mesh interval (radial length units; e.g, km) */
            dq = 0, /* latitudinal mesh interval (degrees) */
            df = 0, /* longitudinal mesh interval (degrees) */
            y0 = 90, /* latitude  of model origin (degrees) */
            x0 = 0, /* longitude of model origin (degrees) */
            z0 = 0, /* depth of model origin (same units as h) */
            x0_save,
            y0_save,
            ym, /* latitude of model midpoint */
            fxs, /* latitude of shot position (degrees) */
            fys, /* longitude of shot position (degrees) */
            fzs, /* depth of shot position (same units as h) */
            fxss, fyss, fzss,
            /* used to detect head waves:  if headwave operator decreases
               the previously-computed traveltime by at least
               headtest*<~time_across_cell> then the headwave counter is
               triggered */
            fhead, headtest = 1.e-3;

    double try, guess;
    double
    /* used in linear velocity gradient cell source */
    rx, ry, rz, dvx, dvy, dvz, dv, v0,
            rzc, rxyc, rxy1, rho, theta1, theta2;

    float
    *slow0, *time0, *wall, *box0,
            s000, slo,
            maxoff = -1.; /* maximum offset (real units) to compute */

    char
    velfile[80], /* file though which velocity structure is input */
            oldtfile[80], /* file through which old travel times are input */
            timefile[1024], /* file in which travel times appear at the end */
            wallfile[80], /* file containing input wall values of traveltimes */
            boxfile[80]; /* file containing precomputed traveltimes to fill in the start box */

    /* ARRAY TO ORDER SIDE FOR SOLUTION IN THE RIGHT ORDER */
    struct sorted *sort;

    /* These variables are for the header */
    struct vhead headin, headout, headbox;

    fprintf(stderr, "Starting sphfd: by S. Roecker 2003, RPI\n");
    /* INITIALIZE PARAMETERS AND ARRAYS */
    setpar(ac, av);
    mstpar("nx", "d", &nx);
    mstpar("ny", "d", &ny);
    mstpar("nz", "d", &nz);
    mstpar("h", "F", &h);
    getpar("dq", "F", &dq);
    getpar("df", "F", &df);
    getpar("iminus", "d", &iminus);
    getpar("iplus", "d", &iplus);
    getpar("jminus", "d", &jminus);
    getpar("jplus", "d", &jplus);
    getpar("kminus", "d", &kminus);
    getpar("kplus", "d", &kplus);
    getpar("floatsrc", "d", &floatsrc);
    getpar("srctype", "d", &srctype);
    getpar("NCUBE", "d", &NCUBE);
    getpar("reverse", "d", &reverse);
    getpar("headpref", "d", &headpref);
    getpar("maxoff", "f", &maxoff);
    getpar("x0", "F", &x0);
    getpar("y0", "F", &y0);
    getpar("z0", "F", &z0);
    getpar("parse", "d", &parse);
    getpar("fill", "d", &fill);
    getpar("invh", "d", &invh);
    getpar("invf", "d", &invf);
    getpar("invq", "d", &invq);
    getpar("invx0", "d", &invx0);
    getpar("invy0", "d", &invy0);
    getpar("invz0", "d", &invz0);
    getpar("swab", "d", &swab);

    if (invh == 1) h = 1. / h;
    if (invq == 1 && dq != 0) dq = 1. / dq;
    if (invf == 1 && df != 0) df = 1. / df;
    if (invx0 == 1 && x0 != 0) x0 = 1. / x0;
    if (invy0 == 1 && y0 != 0) y0 = 1. / y0;
    if (invz0 == 1 && z0 != 0) z0 = 1. / z0;

    /*
       x0 and y0 are geographic latitude and longitude in degrees.   x0 will be converted to radians, y0 to
       geocentric colatitude.
     */
    x0_save = x0;
    y0_save = y0;

    /* if dq and df have not been specified, then make them so that the
       interval at the surface is equal to h */
    y0 *= degrad;
    /* Convert geographic latitude to geocentric colatitude */
    //	y0 = HPI - glat(y0);
    y0 = HPI - glath(y0, -z0, &z0r);
    fprintf(stderr, " y0, z0, z0r = %g %g %g \n", (HPI - y0) / degrad, z0, z0r);
    /* ignore conversion to compare with Hiajiang tests
            y0 = HPI - y0;
     */
    x0 *= degrad;
    dq *= degrad;
    df *= degrad;
    if (dq == 0) dq = h / rearth;
    /*
            ym = y0 + 0.5*(y0 + dq*ny);
            if (df == 0)  df = fabs(h/(rearth*sin(ym)));
     */
    if (df == 0) df = fabs(h / (rearth * sin(y0)));

    fprintf(stderr, " h, dq, df = %g %g %g \n", h, dq, df);

    /* NB: In this version, only srctype=1 has been tested */
    if (srctype == 1) {
        if (floatsrc == 0) {
            if (savsrc == 1) {
                mstpar("fxs", "f", &fxs);
                mstpar("fys", "f", &fys);
                mstpar("fzs", "f", &fzs);
                fxss = fxs;
                fyss = fys;
                fzss = fzs;
                fprintf(stderr, " Saving Source fxss, fyss, fzss = %g %g %g \n", fxss, fyss, fzss);
            }
            mstpar("xs", "d", &xs);
            mstpar("ys", "d", &ys);
            mstpar("zs", "d", &zs);
            fxs = (double) xs;
            fys = (double) ys;
            fzs = (double) zs;
            fprintf(stderr, " Integer Source xs, ys, zs = %g %g %g \n", fc(xs), qc(ys), rc(zs));
            if (savsrc == 0) {
                fxss = fc(xs);
                fyss = qc(ys);
                fzss = rc(zs);
            }
        } else {
            mstpar("fxs", "F", &fxs);
            mstpar("fys", "F", &fys);
            mstpar("fzs", "F", &fzs);
            fxs *= degrad;
            fys *= degrad;
            /*  Use this to ignore the geocentic converstion to test with Haijiang
                            fys = HPI - fys;
             */
            //		fys = HPI - glat(fys);
            fys = HPI - glath(fys, -fzs, &rcent);
            fxss = fxs;
            fyss = fys;
            //		fzss = rearth - fzs;
            fzss = rcent;
            fxs = (fxs - x0) / df;
            fys = (fys - y0) / dq;
            // NB: fzs is the z distance from the grid origin to the the start point.  In the revised version, rcent
            // is the distance from the earth's center to the start point, and z0r is the distance to the outermost grid,
            // so fzs gives us the distance from the outermost grid to the start point; divided by h is the number of
            // nodes (fractional).
            //		fzs = (fzs-z0)/h;
            fzs = (z0r - rcent) / h;
            xs = (int) (fxs + 0.5);
            ys = (int) (fys + 0.5);
            zs = (int) (fzs + 0.5);
            fprintf(stderr, " Float Source xs, ys, zs = %g %g %g \n", fc(xs), qc(ys), rc(zs));
        }
    } else if (srctype == 2) {
        mstpar("srcwall", "d", &srcwall);
        mstpar("wallfile", "s", wallfile);
    } else if (srctype == 3) {
        mstpar("srcwall", "d", &srcwall);
        mstpar("oldtfile", "s", oldtfile);
    } else {
        fprintf(stderr, "ERROR: incorrect value of srctype\n");
        exit(-1);
    }

    if (fill == 1) {
        mstpar("boxfile", "s", boxfile);
    }

    mstpar("timefile", "s", timefile);
    mstpar("velfile", "s", velfile);
    // 20190405 AJL  endpar();

    if (xs < 0 || ys < 0 || zs < 0 || xs > nx - 1 || ys > ny - 1 || zs > nz - 1) {
        fprintf(stderr, "Error: Source does not appear to be in the model;\n");
        fprintf(stderr, "Please check the values in the parameter file.\n");
        exit(-1);
    }

    if (xs < 2 || ys < 2 || zs < 2 || xs > nx - 3 || ys > ny - 3 || zs > nz - 3) {
        fprintf(stderr, "Source near an edge, beware of traveltime errors\n");
        fprintf(stderr, "for raypaths that travel parallel to edge \n");
        fprintf(stderr, "while wavefronts are strongly curved\n");
    }

    /* SET MAXIMUM RADIUS TO COMPUTE */
    if (maxoff > 0.) {
        maxrad = maxoff / h + 1;
        fprintf(stderr, "WARNING: Computing only to max radius = %d\n", maxrad);
    } else maxrad = 99999999;

    /* Check for endiness */
    litend = endian();

    nxy = nx * ny;
    nyz = ny * nz;
    nxz = nx * nz;
    nxyz = nx * ny * nz;

    /* FORM AND FILL TT AND SLOWNESS ARRAYS */
    /* 20190405 AJL
    if ((tfint = open(timefile, O_CREAT | O_WRONLY | O_TRUNC, 0664)) <= 1) {
        fprintf(stderr, "cannot open %s\n", timefile);
        exit(-1);
    }*/
    if ((vfint = open(velfile, O_RDONLY, 0664)) <= 1) {
        fprintf(stderr, "cannot open %s\n", velfile);
        exit(-1);
    }
    if (fill == 1) {
        if ((bfint = open(boxfile, O_RDONLY, 0664)) <= 1) {
            fprintf(stderr, "cannot open %s\n", boxfile);
            exit(-1);
        }
    }
    if (srctype == 2) {
        if ((wfint = open(wallfile, O_RDONLY, 0664)) <= 1) {
            fprintf(stderr, "cannot open %s\n", wallfile);
            exit(-1);
        }
    }
    if (srctype == 3) {
        if ((ofint = open(oldtfile, O_RDONLY, 0664)) <= 1) {
            fprintf(stderr, "cannot open %s\n", oldtfile);
            exit(-1);
        }
    }

    /* ALLOCATE MAIN AND ALTERNATE GRID FOR SLOWNESSES AND TIMES */
    slow0 = (float *) malloc(4 * nxyz);
    time0 = (float *) malloc(4 * nxyz);

    /* MAKE ARRAY SORT LARGE ENOUGH FOR ANY SIDE */
    if (nx <= ny && nx <= nz) {
        sort = (struct sorted *) malloc(sizeof (struct sorted)*ny * nz);
        nwall = nyz;
    } else if (ny <= nx && ny <= nz) {
        sort = (struct sorted *) malloc(sizeof (struct sorted)*nx * nz);
        nwall = nxz;
    } else {
        sort = (struct sorted *) malloc(sizeof (struct sorted)*nx * ny);
        nwall = nxy;
    }
    wall = (float *) malloc(4 * nwall);
    if (slow0 == NULL || time0 == NULL || sort == NULL || wall == NULL) {
        fprintf(stderr, "cannot allocate memory\n");
        exit(-1);
    }
    /* READ IN VELOCITY FILE */
    read(vfint, &headin, 232);
    if (!strncmp(headin.header, "HEAD", 4)) {
        fprintf(stdout, "File Header Identified\n");
        fprintf(stderr, "Original Header  = %s \n", &headin.header[0]);
        fprintf(stderr, " nx, ny, nz = %d %d %d \n", headin.nx, headin.ny, headin.nz);
        fprintf(stderr, " x0, y0, z0 = %g %g %g \n", headin.x0, headin.y0, headin.z0);
        fprintf(stderr, " dx, dy, dz = %g %g %g \n", headin.dx, headin.dy, headin.dz);
    } else {
        fprintf(stdout, "****WARNING: Cannot identify File Header****\n");
        headin.x0 = x0_save;
        headin.y0 = y0_save;
        headin.z0 = z0;
        headin.dx = df;
        headin.dy = dq;
        headin.dz = h;
        headin.nx = nx;
        headin.ny = ny;
        headin.nz = nz;
        headin.az = 0.;
        headin.clat = y0_save;
        headin.clon = x0_save;
        headin.cz = z0;
        for (i = 0; i < 144; i++) strcpy(&headin.header[i], " ");
        /*	  for (i=0; i<124; i++) headin.header[i] = '\0'; */
        strcpy(&headin.header[0], "H");
        strcpy(&headin.header[1], "E");
        strcpy(&headin.header[2], "A");
        strcpy(&headin.header[3], "D");
        strcpy(&headin.header[4], "F");
        strcpy(&headin.header[5], "I");
        strcpy(&headin.header[6], "N");
        strcpy(&headin.header[7], "E");
        strcpy(&headin.header[8], "S");
        strcpy(&headin.header[9], "P");
        strcpy(&headin.header[10], "H");
        strcpy(&headin.header[11], "R");
        strcpy(&headin.header[12], "V");
        strcpy(&headin.header[13], "P");
        strcpy(&headin.header[14], "M");
        strcpy(&headin.header[15], "D");
        strcpy(&headin.header[16], "N");
        strcpy(&headin.header[17], "O");
        strcpy(&headin.header[18], "F");
        strcpy(&headin.header[19], "L");
        strcpy(&headin.header[20], "s");
        strcpy(&headin.header[21], "p");
        strcpy(&headin.header[22], "h");
        strcpy(&headin.header[23], "f");
        strcpy(&headin.header[24], "d");
        close(vfint);
        if ((vfint = open(velfile, O_RDONLY, 0664)) <= 1) {
            fprintf(stderr, "cannot open %s\n", velfile);
            exit(-1);
        }
    }

    read(vfint, slow0, nxyz * 4);

    /* swap bytes on input if necessary */
    /*        if ((litend && swab==1) || swab==2) { */
    if (swab == 1 || swab == 3) {
        fprintf(stdout, "Swapping bytes on header and wavespeeds ..\n");
        FIX_FLOAT(headin.az);
        FIX_INT(headin.nx);
        FIX_INT(headin.ny);
        FIX_INT(headin.nz);
        FixDouble(&headin.clat);
        FixDouble(&headin.clon);
        FixDouble(&headin.cz);
        FixDouble(&headin.x0);
        FixDouble(&headin.y0);
        FixDouble(&headin.z0);
        FixDouble(&headin.dx);
        FixDouble(&headin.dy);
        FixDouble(&headin.dz);
        for (i = 0; i < nxyz; i++) FIX_FLOAT(slow0[i]);
        fprintf(stdout, ".. Done\n");
    }

    /* CONVERT TO SLOWNESS */
    for (i = 0; i < nxyz; i++) slow0[i] = 1 / slow0[i];

    /* SET TIMES TO DUMMY VALUE */
    for (i = 0; i < nxyz; i++) time0[i] = 1.0e10;

    /* Read in precomputed box times */
    if (fill == 1) {
        x1 = xs - NCUBE;
        if (x1 < 0) x1 = 0;
        x2 = xs + NCUBE;
        if (x2 > nx - 1) x2 = nx - 1;
        y1 = ys - NCUBE;
        if (y1 < 0) y1 = 0;
        y2 = ys + NCUBE;
        if (y2 > ny - 1) y2 = ny - 1;
        z1 = zs - NCUBE;
        if (z1 < 0) z1 = 0;
        z2 = zs + NCUBE;
        if (z2 > nz - 1) z2 = nz - 1;
        nbox = (x2 - x1 + 1)*(y2 - y1 + 1)*(z2 - z1 + 1);
        box0 = (float *) malloc(4 * nbox);
        read(bfint, &headbox, 232);
        read(bfint, box0, nbox * 4);
        fprintf(stderr, " Number of points read from boxfile = %d\n", nbox);
    }

    if (srctype == 1) { /* POINT SOURCE */

        /* FILL IN CELL AROUND SOURCE POINT */
        x1p = fzss * cos(fyss);
        y1p = fzss * sin(fyss) * cos(fxss);
        z1p = fzss * sin(fyss) * sin(fxss);
        s000 = sc(xs, ys, zs);
        ii = 0;
        for (zz = zs - NCUBE; zz <= zs + NCUBE; zz++) {
            if (zz < 0 || zz >= nz) continue;
            for (yy = ys - NCUBE; yy <= ys + NCUBE; yy++) {
                if (yy < 0 || yy >= ny) continue;
                for (xx = xs - NCUBE; xx <= xs + NCUBE; xx++) {
                    if (xx < 0 || xx >= nx) continue;
                    /* Constant Wavespeed */
                    if (fill == 0) {
                        fx1 = fc(xx);
                        qy1 = qc(yy);
                        rz1 = rc(zz);
                        x2p = rz1 * cos(qy1);
                        y2p = rz1 * sin(qy1) * cos(fx1);
                        z2p = rz1 * sin(qy1) * sin(fx1);
                        tc(xx, yy, zz) = s000 * DIST(x1p, y1p, z1p, x2p, y2p, z2p);
                        /*
                                fprintf(stderr, " xx, yy, zz, = %d %d %d \n",  xx, yy, zz );
                                fprintf(stderr, " fxss, fyss, fzss = %g %g %g \n", fxss, fyss, fzss);
                                fprintf(stderr, " fx1, qy1, rz1 = %g %g %g \n", fx1, qy1, rz1);
                                fprintf(stderr, " x1p, y1p, z1p = %g %g %g \n", x1p, y1p, z1p);
                                fprintf(stderr, " x2p, y2p, z2p = %g %g %g \n", x2p, y2p, z2p);
                                fprintf(stderr, " DIST = %g \n", DIST(x1p, y1p, z1p, x2p, y2p, z2p));
                                fprintf(stderr, " Vel = %g \n", 1./s000);
                                fprintf(stderr, " Time = %g \n", tc(xx, yy, zz));
                         */
                    } else {
                        /* Assign precomputed times */
                        tc(xx, yy, zz) = box0[ii];
                        ii++;
                        /* Test lines */
                        /*
                                                                        fx1 = fc(xx);
                                                                        qy1 = qc(yy);
                                                                        rz1 = rc(zz);
                                                                        x2p = rz1*cos(qy1);
                                                                        y2p = rz1*sin(qy1)*cos(fx1);
                                                                        z2p = rz1*sin(qy1)*sin(fx1);
                                                                        timec  = s000 * DIST(x1p, y1p, z1p, x2p, y2p, z2p);
                                fprintf(stderr, " i, j, k, x, y, z, timec, box, dt = %d %d %d %g %g %g %g %g %g\n",
                                        xx, yy, zz, fx1, qy1, rz1, timec, tc(xx, yy, zz), (timec-tc(xx, yy, zz))/timec);
                                exit(-1);
                         */
                        /* End test lines */

                    }
                }
            }
        }
        if (fill == 1) fprintf(stderr, " Number of points assigned from boxfile = %d\n", ii);

        /* SETS LOCATION OF THE SIDES OF THE CELL */
        radius = NCUBE;
        if (xs > NCUBE) x1 = xs - (NCUBE + 1);
        else {
            x1 = -1;
            dx1 = 0;
        }
        if (xs < nx - (NCUBE + 1)) x2 = xs + (NCUBE + 1);
        else {
            x2 = nx;
            dx2 = 0;
        }
        if (ys > NCUBE) y1 = ys - (NCUBE + 1);
        else {
            y1 = -1;
            dy1 = 0;
        }
        if (ys < ny - (NCUBE + 1)) y2 = ys + (NCUBE + 1);
        else {
            y2 = ny;
            dy2 = 0;
        }
        if (zs > NCUBE) z1 = zs - (NCUBE + 1);
        else {
            z1 = -1;
            dz1 = 0;
        }
        if (zs < nz - (NCUBE + 1)) z2 = zs + (NCUBE + 1);
        else {
            z2 = nz;
            dz2 = 0;
        }
    } else if (srctype == 2) { /*  EXTERNAL SOURCE */

        /* FILL IN WALLS' TIMES FROM EXTERNAL DATAFILE */
        read(wfint, wall, 4 * nwall); /* READ X=0 WALL */
        if (wall[0]>-1.e-20) {
            ii = 0;
            for (k = 0; k < nz; k++) {
                for (j = 0; j < ny; j++) {
                    tc(0, j, k) = wall[ii];
                    ii++;
                }
            }
        }
        read(wfint, wall, 4 * nwall); /* READ X=NX-1 WALL */
        if (wall[0]>-1.e-20) {
            ii = 0;
            for (k = 0; k < nz; k++) {
                for (j = 0; j < ny; j++) {
                    tc(nx - 1, j, k) = wall[ii];
                    ii++;
                }
            }
        }
        read(wfint, wall, 4 * nwall); /* READ Y=0 WALL */
        if (wall[0]>-1.e-20) {
            ii = 0;
            for (k = 0; k < nz; k++) {
                for (i = 0; i < nx; i++) {
                    tc(i, 0, k) = wall[ii];
                    ii++;
                }
            }
        }
        read(wfint, wall, 4 * nwall); /* READ Y=NY-1 WALL */
        if (wall[0]>-1.e-20) {
            ii = 0;
            for (k = 0; k < nz; k++) {
                for (i = 0; i < nx; i++) {
                    tc(i, ny - 1, k) = wall[ii];
                    ii++;
                }
            }
        }
        read(wfint, wall, 4 * nwall); /* READ Z=0 WALL */
        if (wall[0]>-1.e-20) {
            ii = 0;
            for (j = 0; j < ny; j++) {
                for (i = 0; i < nx; i++) {
                    tc(i, j, 0) = wall[ii];
                    ii++;
                }
            }
        }
        read(wfint, wall, 4 * nwall); /* READ Z=NZ-1 WALL */
        if (wall[0]>-1.e-20) {
            ii = 0;
            for (j = 0; j < ny; j++) {
                for (i = 0; i < nx; i++) {
                    tc(i, j, nz - 1) = wall[ii];
                    ii++;
                }
            }
        }

        /* SET LOCATIONS OF SIDES OF THE CELL SO THAT CELL IS A FACE  */
        radius = 1;
        if (srcwall == 1) x2 = 1;
        else {
            x2 = nx;
            dx2 = 0;
        }
        if (srcwall == 2) x1 = nx - 2;
        else {
            x1 = -1;
            dx1 = 0;
        }
        if (srcwall == 3) y2 = 1;
        else {
            y2 = ny;
            dy2 = 0;
        }
        if (srcwall == 4) y1 = ny - 2;
        else {
            y1 = -1;
            dy1 = 0;
        }
        if (srcwall == 5) z2 = 1;
        else {
            z2 = nz;
            dz2 = 0;
        }
        if (srcwall == 6) z1 = nz - 2;
        else {
            z1 = -1;
            dz1 = 0;
        }
    } else if (srctype == 3) { /*  REDO OLD TIMES */
        /* READ IN OLD TIME FILE */
        if (srctype == 3) read(ofint, time0, nxyz * 4);
        /* SET LOCATIONS OF SIDES OF THE CELL SO THAT CELL IS A FACE */
        radius = 1;
        if (srcwall == 1) x2 = 1;
        else {
            x2 = nx;
            dx2 = 0;
        }
        if (srcwall == 2) x1 = nx - 2;
        else {
            x1 = -1;
            dx1 = 0;
        }
        if (srcwall == 3) y2 = 1;
        else {
            y2 = ny;
            dy2 = 0;
        }
        if (srcwall == 4) y1 = ny - 2;
        else {
            y1 = -1;
            dy1 = 0;
        }
        if (srcwall == 5) z2 = 1;
        else {
            z2 = nz;
            dz2 = 0;
        }
        if (srcwall == 6) z1 = nz - 2;
        else {
            z1 = -1;
            dz1 = 0;
        }
    } else {
        fprintf(stderr, "incorrect value of srctype = %d\n", srctype);
        exit(-1);
    }

    if (headpref > 0) { /* PREFERRED REVERSE DIRECTION */
        head = nx * ny*nz;
        if (nx > 5 && x2 <= head) {
            headpref = 2;
            head = x2;
        }
        if (nx > 5 && (nx - 1 - x1) <= head) {
            headpref = 1;
            head = nx - 1 - x1;
        }
        if (ny > 5 && y2 <= head) {
            headpref = 4;
            head = y2;
        }
        if (ny > 5 && (ny - 1 - y1) <= head) {
            headpref = 3;
            head = ny - 1 - y1;
        }
        if (nz > 5 && z2 <= head) {
            headpref = 6;
            head = z2;
        }
        if (nz > 5 && (nz - 1 - z1) <= head) {
            headpref = 5;
            head = nz - 1 - z1;
        }
    }

    /* BIGGER LOOP - ALLOWS AUTOMATIC REVERSE PROPAGATION IF
            HEAD WAVES ARE ENCOUNTERED ON FACES OF EXPANDING CELL,
            ALLOWING WAVES TO TRAVEL BACK INTO THE CELL */
    while (reverse > -1) {

        headw[1] = 0;
        headw[2] = 0;
        headw[3] = 0;
        headw[4] = 0;
        headw[5] = 0;
        headw[6] = 0;

        /* BIG LOOP */
        while (rad0 && (dx1 || dx2 || dy1 || dy2 || dz1 || dz2)) {

            /* CALCULATE ON PRIMARY (time0) GRID */
            /* TOP SIDE */
            for (igrow = 1; igrow <= kminus; igrow++) {
                if (dz1) {
                    ii = 0;
                    for (j = y1 + 1; j <= y2 - 1; j++) {
                        for (i = x1 + 1; i <= x2 - 1; i++) {
                            sort[ii].time = tc(i, j, z1 + 1);
                            sort[ii].i1 = i;
                            sort[ii].i2 = j;
                            ii++;
                        }
                    }
                    qsort((char *) sort, ii, sizeof (struct sorted), compar);
                    for (i = 0; i < ii; i++) {
                        X1 = sort[i].i1;
                        X2 = sort[i].i2;
                        index = z1 * nxy + X2 * nx + X1;
                        fhead = 0.;
                        guess = time0[index];

                        rz1 = rc(z1);
                        rzp1 = rc(z1 + 1);
                        qX2 = qc(X2);
                        qX2p1 = qc(X2 + 1);
                        qX2m1 = qc(X2 - 1);
                        fx1 = fc(X1);
                        fx1p1 = fc(X1 + 1);
                        fx1m1 = fc(X1 - 1);

                        sinqX2 = sin(qX2);

                        /*** 3D Trials ***/
                        /* NW corner (point 6) */
                        if (X2 < ny - 1 && X1 < nx - 1 && time0[index + 1] < 1.e9 && time0[index + nx + 1] < 1.e9
                                && time0[index + nx] < 1.e9) {
                            /* Point 0 */
                            t[0] = tc(X1 + 1, X2 + 1, z1 + 1);
                            s[0] = sc(X1 + 1, X2 + 1, z1 + 1);
                            r[0] = rzp1;
                            q[0] = qX2p1;
                            f[0] = fx1p1;
                            g[0] = -1;
                            n[0] = 1;
                            m[0] = 1;
                            /* Point 1 */
                            t[1] = tc(X1, X2 + 1, z1 + 1);
                            s[1] = sc(X1, X2 + 1, z1 + 1);
                            r[1] = rzp1;
                            q[1] = qX2p1;
                            f[1] = fx1;
                            g[1] = -1;
                            n[1] = 1;
                            m[1] = -1;
                            /* Point 2 */
                            t[2] = tc(X1, X2, z1 + 1);
                            s[2] = sc(X1, X2, z1 + 1);
                            r[2] = rzp1;
                            q[2] = qX2;
                            f[2] = fx1;
                            g[2] = -1;
                            n[2] = -1;
                            m[2] = -1;
                            /* Point 3 */
                            t[3] = tc(X1 + 1, X2, z1 + 1);
                            s[3] = sc(X1 + 1, X2, z1 + 1);
                            r[3] = rzp1;
                            q[3] = qX2;
                            f[3] = fx1p1;
                            g[3] = -1;
                            n[3] = -1;
                            m[3] = 1;
                            /* Point 4 */
                            t[4] = tc(X1 + 1, X2 + 1, z1);
                            s[4] = sc(X1 + 1, X2 + 1, z1);
                            r[4] = rz1;
                            q[4] = qX2p1;
                            f[4] = fx1p1;
                            g[4] = 1;
                            n[4] = 1;
                            m[4] = 1;
                            /* Point 5 */
                            t[5] = tc(X1, X2 + 1, z1);
                            s[5] = sc(X1, X2 + 1, z1);
                            r[5] = rz1;
                            q[5] = qX2p1;
                            f[5] = fx1;
                            g[5] = 1;
                            n[5] = 1;
                            m[5] = -1;
                            /* Point 7 */
                            t[6] = tc(X1 + 1, X2, z1);
                            s[6] = sc(X1 + 1, X2, z1);
                            r[6] = rz1;
                            q[6] = qX2;
                            f[6] = fx1p1;
                            g[6] = 1;
                            n[6] = -1;
                            m[6] = 1;
                            /* Point 6 */
                            s[7] = sc(X1, X2, z1);
                            r[7] = rz1;
                            q[7] = qX2;
                            f[7] = fx1;
                            g[7] = 1;
                            n[7] = -1;
                            m[7] = -1;
                            try = fdsph3d(t, s, r, q, f, g, n, m, h, dq, df);
                            if (try < guess) guess = try;
                        }

                        /* NE corner (point 7) */
                        if (X2 < ny - 1 && X1 > 0 && time0[index - 1] < 1.e9 && time0[index + nx - 1] < 1.e9
                                && time0[index + nx] < 1.e9) {
                            /* Point 0 */
                            t[0] = tc(X1, X2 + 1, z1 + 1);
                            s[0] = sc(X1, X2 + 1, z1 + 1);
                            r[0] = rzp1;
                            q[0] = qX2p1;
                            f[0] = fx1;
                            g[0] = -1;
                            n[0] = 1;
                            m[0] = 1;
                            /* Point 1 */
                            t[1] = tc(X1 - 1, X2 + 1, z1 + 1);
                            s[1] = sc(X1 - 1, X2 + 1, z1 + 1);
                            r[1] = rzp1;
                            q[1] = qX2p1;
                            f[1] = fx1m1;
                            g[1] = -1;
                            n[1] = 1;
                            m[1] = -1;
                            /* Point 2 */
                            t[2] = tc(X1 - 1, X2, z1 + 1);
                            s[2] = sc(X1 - 1, X2, z1 + 1);
                            r[2] = rzp1;
                            q[2] = qX2;
                            f[2] = fx1m1;
                            g[2] = -1;
                            n[2] = -1;
                            m[2] = -1;
                            /* Point 3 */
                            t[3] = tc(X1, X2, z1 + 1);
                            s[3] = sc(X1, X2, z1 + 1);
                            r[3] = rzp1;
                            q[3] = qX2;
                            f[3] = fx1;
                            g[3] = -1;
                            n[3] = -1;
                            m[3] = 1;
                            /* Point 4 */
                            t[4] = tc(X1, X2 + 1, z1);
                            s[4] = sc(X1, X2 + 1, z1);
                            r[4] = rz1;
                            q[4] = qX2p1;
                            f[4] = fx1;
                            g[4] = 1;
                            n[4] = 1;
                            m[4] = 1;
                            /* Point 5 */
                            t[5] = tc(X1 - 1, X2 + 1, z1);
                            s[5] = sc(X1 - 1, X2 + 1, z1);
                            r[5] = rz1;
                            q[5] = qX2p1;
                            f[5] = fx1m1;
                            g[5] = 1;
                            n[5] = 1;
                            m[5] = -1;
                            /* Point 6 */
                            t[6] = tc(X1 - 1, X2, z1);
                            s[6] = sc(X1 - 1, X2, z1);
                            r[6] = rz1;
                            q[6] = qX2;
                            f[6] = fx1m1;
                            g[6] = 1;
                            n[6] = -1;
                            m[6] = -1;
                            /* Point 7 */
                            s[7] = sc(X1, X2, z1);
                            r[7] = rz1;
                            q[7] = qX2;
                            f[7] = fx1;
                            g[7] = 1;
                            n[7] = -1;
                            m[7] = 1;

                            try = fdsph3d(t, s, r, q, f, g, n, m, h, dq, df);
                            if (try < guess) guess = try;
                        }

                        /* SW corner (point 5) */
                        if (X2 > 0 && X1 < nx - 1 && time0[index + 1] < 1.e9 && time0[index - nx + 1] < 1.e9
                                && time0[index - nx] < 1.e9) {
                            /* Point 0 */
                            t[0] = tc(X1 + 1, X2, z1 + 1);
                            s[0] = sc(X1 + 1, X2, z1 + 1);
                            r[0] = rzp1;
                            q[0] = qX2;
                            f[0] = fx1p1;
                            g[0] = -1;
                            n[0] = 1;
                            m[0] = 1;
                            /* Point 1 */
                            t[1] = tc(X1, X2, z1 + 1);
                            s[1] = sc(X1, X2, z1 + 1);
                            r[1] = rzp1;
                            q[1] = qX2;
                            f[1] = fx1;
                            g[1] = -1;
                            n[1] = 1;
                            m[1] = -1;
                            /* Point 2 */
                            t[2] = tc(X1, X2 - 1, z1 + 1);
                            s[2] = sc(X1, X2 - 1, z1 + 1);
                            r[2] = rzp1;
                            q[2] = qX2m1;
                            f[2] = fx1;
                            g[2] = -1;
                            n[2] = -1;
                            m[2] = -1;
                            /* Point 3 */
                            t[3] = tc(X1 + 1, X2 - 1, z1 + 1);
                            s[3] = sc(X1 + 1, X2 - 1, z1 + 1);
                            r[3] = rzp1;
                            q[3] = qX2m1;
                            f[3] = fx1p1;
                            g[3] = -1;
                            n[3] = -1;
                            m[3] = 1;
                            /* Point 4 */
                            t[4] = tc(X1 + 1, X2, z1);
                            s[4] = sc(X1 + 1, X2, z1);
                            r[4] = rz1;
                            q[4] = qX2;
                            f[4] = fx1p1;
                            g[4] = 1;
                            n[4] = 1;
                            m[4] = 1;
                            /* Point 6 */
                            t[5] = tc(X1, X2 - 1, z1);
                            s[5] = sc(X1, X2 - 1, z1);
                            r[5] = rz1;
                            q[5] = qX2m1;
                            f[5] = fx1;
                            g[5] = 1;
                            n[5] = -1;
                            m[5] = -1;
                            /* Point 7 */
                            t[6] = tc(X1 + 1, X2 - 1, z1);
                            s[6] = sc(X1 + 1, X2 - 1, z1);
                            r[6] = rz1;
                            q[6] = qX2m1;
                            f[6] = fx1p1;
                            g[6] = 1;
                            n[6] = -1;
                            m[6] = 1;
                            /* Point 5 */
                            s[7] = sc(X1, X2, z1);
                            r[7] = rz1;
                            q[7] = qX2;
                            f[7] = fx1;
                            g[7] = 1;
                            n[7] = 1;
                            m[7] = -1;

                            try = fdsph3d(t, s, r, q, f, g, n, m, h, dq, df);
                            if (try < guess) guess = try;
                        }

                        /* SE corner (point 4) */
                        if (X2 > 0 && X1 > 0 && time0[index - 1] < 1.e9 && time0[index - nx - 1] < 1.e9
                                && time0[index - nx] < 1.e9) {
                            /* Point 0 */
                            t[0] = tc(X1, X2, z1 + 1);
                            s[0] = sc(X1, X2, z1 + 1);
                            r[0] = rzp1;
                            q[0] = qX2;
                            f[0] = fx1;
                            g[0] = -1;
                            n[0] = 1;
                            m[0] = 1;
                            /* Point 1 */
                            t[1] = tc(X1 - 1, X2, z1 + 1);
                            s[1] = sc(X1 - 1, X2, z1 + 1);
                            r[1] = rzp1;
                            q[1] = qX2;
                            f[1] = fx1m1;
                            g[1] = -1;
                            n[1] = 1;
                            m[1] = -1;
                            /* Point 2 */
                            t[2] = tc(X1 - 1, X2 - 1, z1 + 1);
                            s[2] = sc(X1 - 1, X2 - 1, z1 + 1);
                            r[2] = rzp1;
                            q[2] = qX2m1;
                            f[2] = fx1m1;
                            g[2] = -1;
                            n[2] = -1;
                            m[2] = -1;
                            /* Point 3 */
                            t[3] = tc(X1, X2 - 1, z1 + 1);
                            s[3] = sc(X1, X2 - 1, z1 + 1);
                            r[3] = rzp1;
                            q[3] = qX2m1;
                            f[3] = fx1;
                            g[3] = -1;
                            n[3] = -1;
                            m[3] = 1;
                            /* Point 5 */
                            t[4] = tc(X1 - 1, X2, z1);
                            s[4] = sc(X1 - 1, X2, z1);
                            r[4] = rz1;
                            q[4] = qX2;
                            f[4] = fx1m1;
                            g[4] = 1;
                            n[4] = 1;
                            m[4] = -1;
                            /* Point 6 */
                            t[5] = tc(X1 - 1, X2 - 1, z1);
                            s[5] = sc(X1 - 1, X2 - 1, z1);
                            r[5] = rz1;
                            q[5] = qX2m1;
                            f[5] = fx1m1;
                            g[5] = 1;
                            n[5] = -1;
                            m[5] = -1;
                            /* Point 7 */
                            t[6] = tc(X1, X2 - 1, z1);
                            s[6] = sc(X1, X2 - 1, z1);
                            r[6] = rz1;
                            q[6] = qX2m1;
                            f[6] = fx1;
                            g[6] = 1;
                            n[6] = -1;
                            m[6] = 1;
                            /* Point 4 */
                            s[7] = sc(X1, X2, z1);
                            r[7] = rz1;
                            q[7] = qX2;
                            f[7] = fx1;
                            g[7] = 1;
                            n[7] = 1;
                            m[7] = 1;
                            try = fdsph3d(t, s, r, q, f, g, n, m, h, dq, df);
                            if (try < guess) guess = try;
                        }

                        /*** New Edge Trials ***/
                        if (guess > 1.0e9) {
                            /* Eastern edge */
                            if (X1 > 0 && X2 > y1 + 1 && X2 < y2 - 1 && time0[index - 1] < 1.e9) {
                                /* Point 0 */
                                t[0] = tc(X1 - 1, X2 + 1, z1 + 1);
                                s[0] = sc(X1 - 1, X2 + 1, z1 + 1);
                                /* Point 1 */
                                t[1] = tc(X1 - 1, X2, z1 + 1);
                                s[1] = sc(X1 - 1, X2, z1 + 1);
                                /* Point 2 */
                                t[2] = tc(X1, X2, z1 + 1);
                                s[2] = sc(X1, X2, z1 + 1);
                                /* Point 3 */
                                t[3] = tc(X1 - 1, X2 - 1, z1 + 1);
                                s[3] = sc(X1 - 1, X2 - 1, z1 + 1);
                                /* Point 4 */
                                t[4] = tc(X1 - 1, X2, z1);
                                s[4] = sc(X1 - 1, X2, z1);
                                /* Point 5 */
                                s[5] = sc(X1, X2, z1);

                                d01 = rzp1*dq;
                                d12 = rzp1 * sinqX2*df;
                                d14 = h;
                                d25 = h;
                                d45 = rz1 * sinqX2*df;

                                try = fdsphne(t, s, d01, d12, d14, d25, d45);
                                if (try < guess) guess = try;
                            }
                            /* Western Edge */
                            if (X1 < nx - 1 && X2 > y1 + 1 && X2 < y2 - 1 && time0[index + 1] < 1.e9) {
                                /* Point 0 */
                                t[0] = tc(X1 + 1, X2 - 1, z1 + 1);
                                s[0] = sc(X1 + 1, X2 - 1, z1 + 1);
                                /* Point 1 */
                                t[1] = tc(X1 + 1, X2, z1 + 1);
                                s[1] = sc(X1 + 1, X2, z1 + 1);
                                /* Point 2 */
                                t[2] = tc(X1, X2, z1 + 1);
                                s[2] = sc(X1, X2, z1 + 1);
                                /* Point 3 */
                                t[3] = tc(X1 + 1, X2 + 1, z1 + 1);
                                s[3] = sc(X1 + 1, X2 + 1, z1 + 1);
                                /* Point 4 */
                                t[4] = tc(X1 + 1, X2, z1);
                                s[4] = sc(X1 + 1, X2, z1);
                                /* Point 5 */
                                s[5] = sc(X1, X2, z1);

                                d01 = rzp1*dq;
                                d12 = rzp1 * sinqX2*df;
                                d14 = h;
                                d25 = h;
                                d45 = rz1 * sinqX2*df;

                                try = fdsphne(t, s, d01, d12, d14, d25, d45);
                                if (try < guess) guess = try;
                            }
                            /* Northern Edge */
                            if (X2 < ny - 1 && X1 > x1 + 1 && X1 < x2 - 1 && time0[index + nx] < 1.e9) {
                                /* Point 0 */
                                t[0] = tc(X1 + 1, X2 + 1, z1 + 1);
                                s[0] = sc(X1 + 1, X2 + 1, z1 + 1);
                                /* Point 1 */
                                t[1] = tc(X1, X2 + 1, z1 + 1);
                                s[1] = sc(X1, X2 + 1, z1 + 1);
                                /* Point 2 */
                                t[2] = tc(X1, X2, z1 + 1);
                                s[2] = sc(X1, X2, z1 + 1);
                                /* Point 3 */
                                t[3] = tc(X1 - 1, X2 + 1, z1 + 1);
                                s[3] = sc(X1 - 1, X2 + 1, z1 + 1);
                                /* Point 4 */
                                t[4] = tc(X1, X2 + 1, z1);
                                s[4] = sc(X1, X2 + 1, z1);
                                /* Point 5 */
                                s[5] = sc(X1, X2, z1);

                                d01 = rzp1 * sin(qX2p1) * df;
                                d12 = rzp1*dq;
                                d14 = h;
                                d25 = h;
                                d45 = rz1*dq;

                                try = fdsphne(t, s, d01, d12, d14, d25, d45);
                                if (try < guess) guess = try;
                            }
                            /* Southern Edge */
                            if (X2 > 0 && X1 > x1 + 1 && X1 < x2 - 1 && time0[index - nx]) {
                                /* Point 0 */
                                t[0] = tc(X1 - 1, X2 - 1, z1 + 1);
                                s[0] = sc(X1 - 1, X2 - 1, z1 + 1);
                                /* Point 1 */
                                t[1] = tc(X1, X2 - 1, z1 + 1);
                                s[1] = sc(X1, X2 - 1, z1 + 1);
                                /* Point 2 */
                                t[2] = tc(X1, X2, z1 + 1);
                                s[2] = sc(X1, X2, z1 + 1);
                                /* Point 3 */
                                t[3] = tc(X1 + 1, X2 - 1, z1 + 1);
                                s[3] = sc(X1 + 1, X2 - 1, z1 + 1);
                                /* Point 4 */
                                t[4] = tc(X1, X2 - 1, z1);
                                s[4] = sc(X1, X2 - 1, z1);
                                /* Point 5 */
                                s[5] = sc(X1, X2, z1);

                                d01 = rzp1 * sin(qX2m1) * df;
                                d12 = rzp1*dq;
                                d14 = h;
                                d25 = h;
                                d45 = rz1*dq;

                                try = fdsphne(t, s, d01, d12, d14, d25, d45);
                                if (try < guess) guess = try;
                            }
                        }

                        /*** Two-D Trials ***/
                        /* Vertical-East */
                        if (X1 < nx - 1 && time0[index + 1] < 1.e9) {
                            s[0] = sc(X1, X2, z1);
                            /* Point 1 */
                            t[1] = tc(X1 + 1, X2, z1);
                            s[1] = sc(X1 + 1, X2, z1);
                            /* Point 2 */
                            t[2] = tc(X1, X2, z1 + 1);
                            s[2] = sc(X1, X2, z1 + 1);
                            /* Point 3 */
                            t[3] = tc(X1 + 1, X2, z1 + 1);
                            s[3] = sc(X1 + 1, X2, z1 + 1);

                            d01 = rz1 * sinqX2*df;
                            d02 = h;
                            d13 = h;
                            d23 = rzp1 * sinqX2*df;
                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) guess = try;
                        }
                        /* Vertical-North */
                        if (X2 > 0 && time0[index - nx] < 1.e9) {
                            s[0] = sc(X1, X2, z1);
                            /* Point 1 */
                            t[1] = tc(X1, X2 - 1, z1);
                            s[1] = sc(X1, X2 - 1, z1);
                            /* Point 2 */
                            t[2] = tc(X1, X2, z1 + 1);
                            s[2] = sc(X1, X2, z1 + 1);
                            /* Point 3 */
                            t[3] = tc(X1, X2 - 1, z1 + 1);
                            s[3] = sc(X1, X2 - 1, z1 + 1);

                            d01 = rz1*dq;
                            d02 = h;
                            d13 = h;
                            d23 = rzp1*dq;

                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) guess = try;
                        }
                        /* Vertical-West */
                        if (X1 > 0 && time0[index - 1] < 1.e9) {
                            s[0] = sc(X1, X2, z1);
                            /* Point 1 */
                            t[1] = tc(X1 - 1, X2, z1);
                            s[1] = sc(X1 - 1, X2, z1);
                            /* Point 2 */
                            t[2] = tc(X1, X2, z1 + 1);
                            s[2] = sc(X1, X2, z1 + 1);
                            /* Point 3 */
                            t[3] = tc(X1 - 1, X2, z1 + 1);
                            s[3] = sc(X1 - 1, X2, z1 + 1);

                            d01 = rz1 * sinqX2*df;
                            d02 = h;
                            d13 = h;
                            d23 = rzp1 * sinqX2*df;

                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) guess = try;

                        }

                        /* Vertical-South */
                        if (X2 < ny - 1 && time0[index + nx] < 1.e9) {
                            s[0] = sc(X1, X2, z1);
                            /* Point 1 */
                            t[1] = tc(X1, X2 + 1, z1);
                            s[1] = sc(X1, X2 + 1, z1);
                            /* Point 2 */
                            t[2] = tc(X1, X2, z1 + 1);
                            s[2] = sc(X1, X2, z1 + 1);
                            /* Point 3 */
                            t[3] = tc(X1, X2 + 1, z1 + 1);
                            s[3] = sc(X1, X2 + 1, z1 + 1);

                            d01 = rz1*dq;
                            d02 = h;
                            d13 = h;
                            d23 = rzp1*dq;

                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) guess = try;
                        }

                        /* Horizontal - SE */
                        if (X2 < ny - 1 && X1 < nx - 1 && time0[index + 1] < 1.e9 && time0[index + nx + 1] < 1.e9
                                && time0[index + nx] < 1.e9) {
                            s[0] = sc(X1, X2, z1);
                            /* Point 1 */
                            t[1] = tc(X1 + 1, X2, z1);
                            s[1] = sc(X1 + 1, X2, z1);
                            /* Point 2 */
                            t[2] = tc(X1, X2 + 1, z1);
                            s[2] = sc(X1, X2 + 1, z1);
                            /* Point 3 */
                            t[3] = tc(X1 + 1, X2 + 1, z1);
                            s[3] = sc(X1 + 1, X2 + 1, z1);

                            d01 = rz1 * sinqX2*df;
                            d02 = rz1*dq;
                            d13 = rz1*dq;
                            d23 = rz1 * sin(qX2p1) * df;

                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) {
                                fhead = (guess - try) / (d01 * slow0[index]);
                                guess = try;
                            }

                        }

                        /* Horizontal - NE */
                        if (X2 > 0 && X1 < nx - 1 && time0[index + 1] < 1.e9 && time0[index - nx + 1] < 1.e9
                                && time0[index - nx] < 1.e9) {
                            s[0] = sc(X1, X2, z1);
                            /* Point 1 */
                            t[1] = tc(X1 + 1, X2, z1);
                            s[1] = sc(X1 + 1, X2, z1);
                            /* Point 2 */
                            t[2] = tc(X1, X2 - 1, z1);
                            s[2] = sc(X1, X2 - 1, z1);
                            /* Point 3 */
                            t[3] = tc(X1 + 1, X2 - 1, z1);
                            s[3] = sc(X1 + 1, X2 - 1, z1);

                            d01 = rz1 * sinqX2*df;
                            d02 = rz1*dq;
                            d13 = rz1*dq;
                            d23 = rz1 * sin(qX2m1) * df;

                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) {
                                fhead = (guess - try) / (d01 * slow0[index]);
                                guess = try;
                            }
                        }

                        /* Horizontal - SW */
                        if (X2 < ny - 1 && X1 > 0 && time0[index - 1] < 1.e9 && time0[index + nx - 1] < 1.e9
                                && time0[index + nx] < 1.e9) {
                            s[0] = sc(X1, X2, z1);
                            /* Point 1 */
                            t[1] = tc(X1 - 1, X2, z1);
                            s[1] = sc(X1 - 1, X2, z1);
                            /* Point 2 */
                            t[2] = tc(X1, X2 + 1, z1);
                            s[2] = sc(X1, X2 + 1, z1);
                            /* Point 3 */
                            t[3] = tc(X1 - 1, X2 + 1, z1);
                            s[3] = sc(X1 - 1, X2 + 1, z1);

                            d01 = rz1 * sinqX2*df;
                            d02 = rz1*dq;
                            d13 = rz1*dq;
                            d23 = rz1 * sin(qX2p1) * df;

                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) {
                                fhead = (guess - try) / (d01 * slow0[index]);
                                guess = try;
                            }
                        }

                        /* Horizontal - NW */
                        if (X2 > 0 && X1 > 0 && time0[index - 1] < 1.e9 && time0[index - nx - 1] < 1.e9
                                && time0[index - nx] < 1.e9) {
                            s[0] = sc(X1, X2, z1);
                            /* Point 1 */
                            t[1] = tc(X1 - 1, X2, z1);
                            s[1] = sc(X1 - 1, X2, z1);
                            /* Point 2 */
                            t[2] = tc(X1, X2 - 1, z1);
                            s[2] = sc(X1, X2 - 1, z1);
                            /* Point 3 */
                            t[3] = tc(X1 - 1, X2 - 1, z1);
                            s[3] = sc(X1 - 1, X2 - 1, z1);

                            d01 = rz1 * sinqX2*df;
                            d02 = rz1*dq;
                            d13 = rz1*dq;
                            d23 = rz1 * sin(qX2m1) * df;

                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) {
                                fhead = (guess - try) / (d01 * slow0[index]);
                                guess = try;
                            }
                        }

                        /*** New Face Trials ***/
                        if (guess > 1.0e9) {
                            if (X1 > x1 + 1 && X1 < x2 - 1 && X2 > y1 + 1 && X2 < y2 - 1) {
                                /* Point 0 */
                                t[0] = tc(X1, X2 + 1, z1 + 1);
                                s[0] = sc(X1, X2 + 1, z1 + 1);
                                /* Point 1 */
                                t[1] = tc(X1 - 1, X2, z1 + 1);
                                s[1] = sc(X1 - 1, X2, z1 + 1);
                                /* Point 2 */
                                t[2] = tc(X1, X2, z1 + 1);
                                s[2] = sc(X1, X2, z1 + 1);
                                /* Point 3 */
                                t[3] = tc(X1 + 1, X2, z1 + 1);
                                s[3] = sc(X1 + 1, X2, z1 + 1);
                                /* Point 4 */
                                t[4] = tc(X1, X2 - 1, z1 + 1);
                                s[4] = sc(X1, X2 - 1, z1 + 1);
                                /* Point 5 */
                                s[5] = sc(X1, X2, z1);
                                d02 = rzp1*dq;
                                d12 = rzp1 * sinqX2*df;
                                d25 = h;
                                try = fdsphnf(t, s, d02, d12, d25);
                                if (try < guess) guess = try;
                            }
                        }

                        /*** Edge Trials ***/
                        /* Change in Depth */
                        try = tc(X1, X2, z1 + 1) + .5 * (sc(X1, X2, z1) + sc(X1, X2, z1 + 1)) * h;
                        if (try < guess) guess = try;
                        /* Change in East Longitude */
                        if (X1 < nx - 1 && time0[index + 1] < 1.e9) {
                            try = tc(X1 + 1, X2, z1) + .5 * (sc(X1, X2, z1) + sc(X1 + 1, X2, z1)) * rz1 * df*sinqX2;
                            if (try < guess) {
                                fhead = (guess - try) / (rz1 * df * sinqX2 * slow0[index]);
                                guess = try;
                            }
                        }
                        /* Change in West Longitude */
                        if (X1 > 0 && time0[index - 1] < 1.e9) {
                            try = tc(X1 - 1, X2, z1) + .5 * (sc(X1, X2, z1) + sc(X1 - 1, X2, z1)) * rz1 * df*sinqX2;
                            if (try < guess) {
                                fhead = (guess - try) / (rz1 * df * sinqX2 * slow0[index]);
                                guess = try;
                            }
                        }
                        /* Change in South Latitude */
                        if (X2 < ny - 1 & time0[index + nx] < 1.e9) {
                            try = tc(X1, X2 + 1, z1) + .5 * (sc(X1, X2, z1) + sc(X1, X2 + 1, z1)) * rz1*dq;
                            if (try < guess) {
                                fhead = (guess - try) / (rz1 * dq * slow0[index]);
                                guess = try;
                            }
                        }
                        /* Change in North Latitude */
                        if (X2 > 0 && time0[index - nx] < 1.e9) {
                            try = tc(X1, X2 - 1, z1) + .5 * (sc(X1, X2, z1) + sc(X1, X2 - 1, z1)) * rz1*dq;
                            if (try < guess) {
                                fhead = (guess - try) / (rz1 * dq * slow0[index]);
                                guess = try;
                            }
                        }
                        if (guess < time0[index]) {
                            time0[index] = guess;
                            if (fhead > headtest) headw[5]++;
                        }
                    }
                    if (z1 == 0) dz1 = 0;
                    z1--;
                }
                /* End z1 if condition */
            }
            /* End growth loop */

            /* BOTTOM SIDE */
            for (igrow = 1; igrow <= kplus; igrow++) {
                if (dz2) {
                    ii = 0;
                    for (j = y1 + 1; j <= y2 - 1; j++) {
                        for (i = x1 + 1; i <= x2 - 1; i++) {
                            sort[ii].time = tc(i, j, z2 - 1);
                            sort[ii].i1 = i;
                            sort[ii].i2 = j;
                            ii++;
                        }
                    }
                    qsort((char *) sort, ii, sizeof (struct sorted), compar);
                    for (i = 0; i < ii; i++) {
                        X1 = sort[i].i1;
                        X2 = sort[i].i2;
                        index = z2 * nxy + X2 * nx + X1;
                        fhead = 0.;
                        guess = time0[index];

                        rz2 = rc(z2);
                        rzm1 = rc(z2 - 1);
                        qX2 = qc(X2);
                        qX2p1 = qc(X2 + 1);
                        qX2m1 = qc(X2 - 1);
                        fx1 = fc(X1);
                        fx1p1 = fc(X1 + 1);
                        fx1m1 = fc(X1 - 1);

                        sinqX2 = sin(qX2);

                        /*** 3D Trials ***/
                        /* NW corner (point 2) */
                        if (X2 < ny - 1 && X1 < nx - 1 && time0[index + 1] < 1.e9 && time0[index + nx + 1] < 1.e9
                                && time0[index + nx] < 1.e9) {
                            /* Point 4 */
                            t[0] = tc(X1 + 1, X2 + 1, z2 - 1);
                            s[0] = sc(X1 + 1, X2 + 1, z2 - 1);
                            r[0] = rzm1;
                            q[0] = qX2p1;
                            f[0] = fx1p1;
                            g[0] = 1;
                            n[0] = 1;
                            m[0] = 1;
                            /* Point 5 */
                            t[1] = tc(X1, X2 + 1, z2 - 1);
                            s[1] = sc(X1, X2 + 1, z2 - 1);
                            r[1] = rzm1;
                            q[1] = qX2p1;
                            f[1] = fx1;
                            g[1] = 1;
                            n[1] = 1;
                            m[1] = -1;
                            /* Point 6 */
                            t[2] = tc(X1, X2, z2 - 1);
                            s[2] = sc(X1, X2, z2 - 1);
                            r[2] = rzm1;
                            q[2] = qX2;
                            f[2] = fx1;
                            g[2] = 1;
                            n[2] = -1;
                            m[2] = -1;
                            /* Point 7 */
                            t[3] = tc(X1 + 1, X2, z2 - 1);
                            s[3] = sc(X1 + 1, X2, z2 - 1);
                            r[3] = rzm1;
                            q[3] = qX2;
                            f[3] = fx1p1;
                            g[3] = 1;
                            n[3] = -1;
                            m[3] = 1;
                            /* Point 0 */
                            t[4] = tc(X1 + 1, X2 + 1, z2);
                            s[4] = sc(X1 + 1, X2 + 1, z2);
                            r[4] = rz2;
                            q[4] = qX2p1;
                            f[4] = fx1p1;
                            g[4] = -1;
                            n[4] = 1;
                            m[4] = 1;
                            /* Point 1 */
                            t[5] = tc(X1, X2 + 1, z2);
                            s[5] = sc(X1, X2 + 1, z2);
                            r[5] = rz2;
                            q[5] = qX2p1;
                            f[5] = fx1;
                            g[5] = -1;
                            n[5] = 1;
                            m[5] = -1;
                            /* Point 3 */
                            t[6] = tc(X1 + 1, X2, z2);
                            s[6] = sc(X1 + 1, X2, z2);
                            r[6] = rz2;
                            q[6] = qX2;
                            f[6] = fx1p1;
                            g[6] = -1;
                            n[6] = -1;
                            m[6] = 1;
                            /* Point 2 */
                            s[7] = sc(X1, X2, z2);
                            r[7] = rz2;
                            q[7] = qX2;
                            f[7] = fx1;
                            g[7] = -1;
                            n[7] = -1;
                            m[7] = -1;

                            try = fdsph3d(t, s, r, q, f, g, n, m, h, dq, df);
                            if (try < guess) guess = try;

                        }

                        /* NE corner (point 3) */
                        if (X2 < ny - 1 && X1 > 0 && time0[index - 1] < 1.e9 && time0[index + nx - 1] < 1.e9
                                && time0[index + nx] < 1.e9) {
                            /* Point 4 */
                            t[0] = tc(X1, X2 + 1, z2 - 1);
                            s[0] = sc(X1, X2 + 1, z2 - 1);
                            r[0] = rzm1;
                            q[0] = qX2p1;
                            f[0] = fx1;
                            g[0] = 1;
                            n[0] = 1;
                            m[0] = 1;
                            /* Point 5 */
                            t[1] = tc(X1 - 1, X2 + 1, z2 - 1);
                            s[1] = sc(X1 - 1, X2 + 1, z2 - 1);
                            r[1] = rzm1;
                            q[1] = qX2p1;
                            f[1] = fx1m1;
                            g[1] = 1;
                            n[1] = 1;
                            m[1] = -1;
                            /* Point 6 */
                            t[2] = tc(X1 - 1, X2, z2 - 1);
                            s[2] = sc(X1 - 1, X2, z2 - 1);
                            r[2] = rzm1;
                            q[2] = qX2;
                            f[2] = fx1m1;
                            g[2] = 1;
                            n[2] = -1;
                            m[2] = -1;
                            /* Point 7 */
                            t[3] = tc(X1, X2, z2 - 1);
                            s[3] = sc(X1, X2, z2 - 1);
                            r[3] = rzm1;
                            q[3] = qX2;
                            f[3] = fx1;
                            g[3] = 1;
                            n[3] = -1;
                            m[3] = 1;
                            /* Point 0 */
                            t[4] = tc(X1, X2 + 1, z2);
                            s[4] = sc(X1, X2 + 1, z2);
                            r[4] = rz2;
                            q[4] = qX2p1;
                            f[4] = fx1;
                            g[4] = -1;
                            n[4] = 1;
                            m[4] = 1;
                            /* Point 1 */
                            t[5] = tc(X1 - 1, X2 + 1, z2);
                            s[5] = sc(X1 - 1, X2 + 1, z2);
                            r[5] = rz2;
                            q[5] = qX2p1;
                            f[5] = fx1m1;
                            g[5] = -1;
                            n[5] = 1;
                            m[5] = -1;
                            /* Point 2 */
                            t[6] = tc(X1 - 1, X2, z2);
                            s[6] = sc(X1 - 1, X2, z2);
                            r[6] = rz2;
                            q[6] = qX2;
                            f[6] = fx1m1;
                            g[6] = -1;
                            n[6] = -1;
                            m[6] = -1;
                            /* Point 3 */
                            s[7] = sc(X1, X2, z2);
                            r[7] = rz2;
                            q[7] = qX2;
                            f[7] = fx1;
                            g[7] = -1;
                            n[7] = -1;
                            m[7] = 1;

                            try = fdsph3d(t, s, r, q, f, g, n, m, h, dq, df);
                            if (try < guess) guess = try;

                        }

                        /* SW corner (point 1) */
                        if (X2 > 0 && X1 < nx - 1 && time0[index + 1] < 1.e9 && time0[index - nx + 1] < 1.e9
                                && time0[index - nx] < 1.e9) {
                            /* Point 4 */
                            t[0] = tc(X1 + 1, X2, z2 - 1);
                            s[0] = sc(X1 + 1, X2, z2 - 1);
                            r[0] = rzm1;
                            q[0] = qX2;
                            f[0] = fx1p1;
                            g[0] = 1;
                            n[0] = 1;
                            m[0] = 1;
                            /* Point 5 */
                            t[1] = tc(X1, X2, z2 - 1);
                            s[1] = sc(X1, X2, z2 - 1);
                            r[1] = rzm1;
                            q[1] = qX2;
                            f[1] = fx1;
                            g[1] = 1;
                            n[1] = 1;
                            m[1] = -1;
                            /* Point 6 */
                            t[2] = tc(X1, X2 - 1, z2 - 1);
                            s[2] = sc(X1, X2 - 1, z2 - 1);
                            r[2] = rzm1;
                            q[2] = qX2m1;
                            f[2] = fx1;
                            g[2] = 1;
                            n[2] = -1;
                            m[2] = -1;
                            /* Point 7 */
                            t[3] = tc(X1 + 1, X2 - 1, z2 - 1);
                            s[3] = sc(X1 + 1, X2 - 1, z2 - 1);
                            r[3] = rzm1;
                            q[3] = qX2m1;
                            f[3] = fx1p1;
                            g[3] = 1;
                            n[3] = -1;
                            m[3] = 1;
                            /* Point 0 */
                            t[4] = tc(X1 + 1, X2, z2);
                            s[4] = sc(X1 + 1, X2, z2);
                            r[4] = rz2;
                            q[4] = qX2;
                            f[4] = fx1p1;
                            g[4] = -1;
                            n[4] = 1;
                            m[4] = 1;
                            /* Point 2 */
                            t[5] = tc(X1, X2 - 1, z2);
                            s[5] = sc(X1, X2 - 1, z2);
                            r[5] = rz2;
                            q[5] = qX2m1;
                            f[5] = fx1;
                            g[5] = -1;
                            n[5] = -1;
                            m[5] = -1;
                            /* Point 3 */
                            t[6] = tc(X1 + 1, X2 - 1, z2);
                            s[6] = sc(X1 + 1, X2 - 1, z2);
                            r[6] = rz2;
                            q[6] = qX2m1;
                            f[6] = fx1p1;
                            g[6] = -1;
                            n[6] = -1;
                            m[6] = 1;
                            /* Point 1 */
                            s[7] = sc(X1, X2, z2);
                            r[7] = rz2;
                            q[7] = qX2;
                            f[7] = fx1;
                            g[7] = -1;
                            n[7] = 1;
                            m[7] = -1;

                            try = fdsph3d(t, s, r, q, f, g, n, m, h, dq, df);
                            if (try < guess) guess = try;

                        }

                        /* SE corner (point 0) */
                        if (X2 > 0 && X1 > 0 && time0[index - 1] < 1.e9 && time0[index - nx - 1] < 1.e9
                                && time0[index - nx] < 1.e9) {
                            /* Point 4 */
                            t[0] = tc(X1, X2, z2 - 1);
                            s[0] = sc(X1, X2, z2 - 1);
                            r[0] = rzm1;
                            q[0] = qX2;
                            f[0] = fx1;
                            g[0] = 1;
                            n[0] = 1;
                            m[0] = 1;
                            /* Point 5 */
                            t[1] = tc(X1 - 1, X2, z2 - 1);
                            s[1] = sc(X1 - 1, X2, z2 - 1);
                            r[1] = rzm1;
                            q[1] = qX2;
                            f[1] = fx1m1;
                            g[1] = 1;
                            n[1] = 1;
                            m[1] = -1;
                            /* Point 6 */
                            t[2] = tc(X1 - 1, X2 - 1, z2 - 1);
                            s[2] = sc(X1 - 1, X2 - 1, z2 - 1);
                            r[2] = rzm1;
                            q[2] = qX2m1;
                            f[2] = fx1m1;
                            g[2] = 1;
                            n[2] = -1;
                            m[2] = -1;
                            /* Point 7 */
                            t[3] = tc(X1, X2 - 1, z2 - 1);
                            s[3] = sc(X1, X2 - 1, z2 - 1);
                            r[3] = rzm1;
                            q[3] = qX2m1;
                            f[3] = fx1;
                            g[3] = 1;
                            n[3] = -1;
                            m[3] = 1;
                            /* Point 1 */
                            t[4] = tc(X1 - 1, X2, z2);
                            s[4] = sc(X1 - 1, X2, z2);
                            r[4] = rz2;
                            q[4] = qX2;
                            f[4] = fx1m1;
                            g[4] = -1;
                            n[4] = 1;
                            m[4] = -1;
                            /* Point 2 */
                            t[5] = tc(X1 - 1, X2 - 1, z2);
                            s[5] = sc(X1 - 1, X2 - 1, z2);
                            r[5] = rz2;
                            q[5] = qX2m1;
                            f[5] = fx1m1;
                            g[5] = -1;
                            n[5] = -1;
                            m[5] = -1;
                            /* Point 3 */
                            t[6] = tc(X1, X2 - 1, z2);
                            s[6] = sc(X1, X2 - 1, z2);
                            r[6] = rz2;
                            q[6] = qX2m1;
                            f[6] = fx1;
                            g[6] = -1;
                            n[6] = -1;
                            m[6] = 1;
                            /* Point 0 */
                            s[7] = sc(X1, X2, z2);
                            r[7] = rz2;
                            q[7] = qX2;
                            f[7] = fx1;
                            g[7] = -1;
                            n[7] = 1;
                            m[7] = 1;
                            try = fdsph3d(t, s, r, q, f, g, n, m, h, dq, df);
                            if (try < guess) guess = try;

                        }

                        /*** New Edge Trials ***/
                        if (guess > 1.0e9) {
                            /* Eastern edge */
                            if (X1 > 0 && X2 > y1 + 1 && X2 < y2 - 1 && time0[index - 1] < 1.e9) {
                                /* Point 0 */
                                t[0] = tc(X1 - 1, X2 - 1, z2 - 1);
                                s[0] = sc(X1 - 1, X2 - 1, z2 - 1);
                                /* Point 1 */
                                t[1] = tc(X1 - 1, X2, z2 - 1);
                                s[1] = sc(X1 - 1, X2, z2 - 1);
                                /* Point 2 */
                                t[2] = tc(X1, X2, z2 - 1);
                                s[2] = sc(X1, X2, z2 - 1);
                                /* Point 3 */
                                t[3] = tc(X1 - 1, X2 + 1, z2 - 1);
                                s[3] = sc(X1 - 1, X2 + 1, z2 - 1);
                                /* Point 4 */
                                t[4] = tc(X1 - 1, X2, z2);
                                s[4] = sc(X1 - 1, X2, z2);
                                /* Point 5 */
                                s[5] = sc(X1, X2, z2);

                                d01 = rzm1*dq;
                                d12 = rzm1 * sinqX2*df;
                                d14 = h;
                                d25 = h;
                                d45 = rz2 * sinqX2*df;

                                try = fdsphne(t, s, d01, d12, d14, d25, d45);
                                if (try < guess) guess = try;

                            }
                            /* Western Edge */
                            if (X1 < nx - 1 && X2 > y1 + 1 && X2 < y2 - 1 && time0[index + 1] < 1.e9) {
                                /* Point 0 */
                                t[0] = tc(X1 + 1, X2 + 1, z2 - 1);
                                s[0] = sc(X1 + 1, X2 + 1, z2 - 1);
                                /* Point 1 */
                                t[1] = tc(X1 + 1, X2, z2 - 1);
                                s[1] = sc(X1 + 1, X2, z2 - 1);
                                /* Point 2 */
                                t[2] = tc(X1, X2, z2 - 1);
                                s[2] = sc(X1, X2, z2 - 1);
                                /* Point 3 */
                                t[3] = tc(X1 + 1, X2 - 1, z2 - 1);
                                s[3] = sc(X1 + 1, X2 - 1, z2 - 1);
                                /* Point 4 */
                                t[4] = tc(X1 + 1, X2, z2);
                                s[4] = sc(X1 + 1, X2, z2);
                                /* Point 5 */
                                s[5] = sc(X1, X2, z2);

                                d01 = rzm1*dq;
                                d12 = rzm1 * sinqX2*df;
                                d14 = h;
                                d25 = h;
                                d45 = rz2 * sinqX2*df;

                                try = fdsphne(t, s, d01, d12, d14, d25, d45);
                                if (try < guess) guess = try;

                            }
                            /* Northern Edge */
                            if (X2 < ny - 1 && X1 > x1 + 1 && X1 < x2 - 1 && time0[index + nx] < 1.e9) {
                                /* Point 0 */
                                t[0] = tc(X1 - 1, X2 + 1, z2 - 1);
                                s[0] = sc(X1 - 1, X2 + 1, z2 - 1);
                                /* Point 1 */
                                t[1] = tc(X1, X2 + 1, z2 - 1);
                                s[1] = sc(X1, X2 + 1, z2 - 1);
                                /* Point 2 */
                                t[2] = tc(X1, X2, z2 - 1);
                                s[2] = sc(X1, X2, z2 - 1);
                                /* Point 3 */
                                t[3] = tc(X1 + 1, X2 + 1, z2 - 1);
                                s[3] = sc(X1 + 1, X2 + 1, z2 - 1);
                                /* Point 4 */
                                t[4] = tc(X1, X2 + 1, z2);
                                s[4] = sc(X1, X2 + 1, z2);
                                /* Point 5 */
                                s[5] = sc(X1, X2, z2);

                                d01 = rzm1 * sin(qX2p1) * df;
                                d12 = rzm1*dq;
                                d14 = h;
                                d25 = h;
                                d45 = rz2*dq;

                                try = fdsphne(t, s, d01, d12, d14, d25, d45);
                                if (try < guess) guess = try;

                            }
                            /* Southern Edge */
                            if (X2 > 0 && X1 > x1 + 1 && X1 < x2 - 1 && time0[index - nx] < 1.e9) {
                                /* Point 0 */
                                t[0] = tc(X1 + 1, X2 - 1, z2 - 1);
                                s[0] = sc(X1 + 1, X2 - 1, z2 - 1);
                                /* Point 1 */
                                t[1] = tc(X1, X2 - 1, z2 - 1);
                                s[1] = sc(X1, X2 - 1, z2 - 1);
                                /* Point 2 */
                                t[2] = tc(X1, X2, z2 - 1);
                                s[2] = sc(X1, X2, z2 - 1);
                                /* Point 3 */
                                t[3] = tc(X1 - 1, X2 - 1, z2 - 1);
                                s[3] = sc(X1 - 1, X2 - 1, z2 - 1);
                                /* Point 4 */
                                t[4] = tc(X1, X2 - 1, z2);
                                s[4] = sc(X1, X2 - 1, z2);
                                /* Point 5 */
                                s[5] = sc(X1, X2, z2);

                                d01 = rzm1 * sin(qX2m1) * df;
                                d12 = rzm1*dq;
                                d14 = h;
                                d25 = h;
                                d45 = rz2*dq;

                                try = fdsphne(t, s, d01, d12, d14, d25, d45);
                                if (try < guess) guess = try;

                            }
                        }

                        /*** Two-D Trials ***/
                        /* Vertical-East */
                        if (X1 < nx - 1 && time0[index + 1] < 1.e9) {
                            s[0] = sc(X1, X2, z2);
                            /* Point 1 */
                            t[1] = tc(X1 + 1, X2, z2);
                            s[1] = sc(X1 + 1, X2, z2);
                            /* Point 2 */
                            t[2] = tc(X1, X2, z2 - 1);
                            s[2] = sc(X1, X2, z2 - 1);
                            /* Point 3 */
                            t[3] = tc(X1 + 1, X2, z2 - 1);
                            s[3] = sc(X1 + 1, X2, z2 - 1);

                            d01 = rz2 * sinqX2*df;
                            d02 = h;
                            d13 = h;
                            d23 = rzm1 * sinqX2*df;
                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) guess = try;
                        }
                        /* Vertical-North */
                        if (X2 > 0 && time0[index - nx] < 1.e9) {
                            s[0] = sc(X1, X2, z2);
                            /* Point 1 */
                            t[1] = tc(X1, X2 - 1, z2);
                            s[1] = sc(X1, X2 - 1, z2);
                            /* Point 2 */
                            t[2] = tc(X1, X2, z2 - 1);
                            s[2] = sc(X1, X2, z2 - 1);
                            /* Point 3 */
                            t[3] = tc(X1, X2 - 1, z2 - 1);
                            s[3] = sc(X1, X2 - 1, z2 - 1);

                            d01 = rz2*dq;
                            d02 = h;
                            d13 = h;
                            d23 = rzm1*dq;

                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) guess = try;
                        }
                        /* Vertical-West */
                        if (X1 > 0 && time0[index - 1] < 1.e9) {
                            s[0] = sc(X1, X2, z2);
                            /* Point 1 */
                            t[1] = tc(X1 - 1, X2, z2);
                            s[1] = sc(X1 - 1, X2, z2);
                            /* Point 2 */
                            t[2] = tc(X1, X2, z2 - 1);
                            s[2] = sc(X1, X2, z2 - 1);
                            /* Point 3 */
                            t[3] = tc(X1 - 1, X2, z2 - 1);
                            s[3] = sc(X1 - 1, X2, z2 - 1);

                            d01 = rz2 * sinqX2*df;
                            d02 = h;
                            d13 = h;
                            d23 = rzm1 * sinqX2*df;

                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) guess = try;
                        }
                        /* Vertical-South */
                        if (X2 < ny - 1 && time0[index + nx] < 1.e9) {
                            s[0] = sc(X1, X2, z2);
                            /* Point 1 */
                            t[1] = tc(X1, X2 + 1, z2);
                            s[1] = sc(X1, X2 + 1, z2);
                            /* Point 2 */
                            t[2] = tc(X1, X2, z2 - 1);
                            s[2] = sc(X1, X2, z2 - 1);
                            /* Point 3 */
                            t[3] = tc(X1, X2 + 1, z2 - 1);
                            s[3] = sc(X1, X2 + 1, z2 - 1);

                            d01 = rz2*dq;
                            d02 = h;
                            d13 = h;
                            d23 = rzm1*dq;

                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) guess = try;
                        }
                        /* Horizontal - SE */
                        if (X2 < ny - 1 && X1 < nx - 1 && time0[index + 1] < 1.e9 && time0[index + nx + 1] < 1.e9
                                && time0[index + nx] < 1.e9) {
                            s[0] = sc(X1, X2, z2);
                            /* Point 1 */
                            t[1] = tc(X1 + 1, X2, z2);
                            s[1] = sc(X1 + 1, X2, z2);
                            /* Point 2 */
                            t[2] = tc(X1, X2 + 1, z2);
                            s[2] = sc(X1, X2 + 1, z2);
                            /* Point 3 */
                            t[3] = tc(X1 + 1, X2 + 1, z2);
                            s[3] = sc(X1 + 1, X2 + 1, z2);
                            d01 = rz2 * sinqX2*df;
                            d02 = rz2*dq;
                            d13 = rz2*dq;
                            d23 = rz2 * sin(qX2p1) * df;
                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) {
                                fhead = (guess - try) / (d01 * slow0[index]);
                                guess = try;
                            }
                        }
                        /* Horizontal - NE */
                        if (X2 > 0 && X1 < nx - 1 && time0[index + 1] < 1.e9 && time0[index - nx + 1] < 1.e9
                                && time0[index - nx] < 1.e9) {
                            s[0] = sc(X1, X2, z2);
                            /* Point 1 */
                            t[1] = tc(X1 + 1, X2, z2);
                            s[1] = sc(X1 + 1, X2, z2);
                            /* Point 2 */
                            t[2] = tc(X1, X2 - 1, z2);
                            s[2] = sc(X1, X2 - 1, z2);
                            /* Point 3 */
                            t[3] = tc(X1 + 1, X2 - 1, z2);
                            s[3] = sc(X1 + 1, X2 - 1, z2);
                            d01 = rz2 * sinqX2*df;
                            d02 = rz2*dq;
                            d13 = rz2*dq;
                            d23 = rz2 * sin(qX2m1) * df;
                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) {
                                fhead = (guess - try) / (d01 * slow0[index]);
                                guess = try;
                            }
                        }
                        /* Horizontal - SW */
                        if (X2 < ny - 1 && X1 > 0 && time0[index - 1] < 1.e9 && time0[index + nx - 1] < 1.e9
                                && time0[index + nx] < 1.e9) {
                            s[0] = sc(X1, X2, z2);
                            /* Point 1 */
                            t[1] = tc(X1 - 1, X2, z2);
                            s[1] = sc(X1 - 1, X2, z2);
                            /* Point 2 */
                            t[2] = tc(X1, X2 + 1, z2);
                            s[2] = sc(X1, X2 + 1, z2);
                            /* Point 3 */
                            t[3] = tc(X1 - 1, X2 + 1, z2);
                            s[3] = sc(X1 - 1, X2 + 1, z2);
                            d01 = rz2 * sinqX2*df;
                            d02 = rz2*dq;
                            d13 = rz2*dq;
                            d23 = rz2 * sin(qX2p1) * df;
                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) {
                                fhead = (guess - try) / (d01 * slow0[index]);
                                guess = try;
                            }
                        }
                        /* Horizontal - NW */
                        if (X2 > 0 && X1 > 0 && time0[index - 1] < 1.e9 && time0[index - nx - 1] < 1.e9
                                && time0[index - nx] < 1.e9) {
                            s[0] = sc(X1, X2, z2);
                            /* Point 1 */
                            t[1] = tc(X1 - 1, X2, z2);
                            s[1] = sc(X1 - 1, X2, z2);
                            /* Point 2 */
                            t[2] = tc(X1, X2 - 1, z2);
                            s[2] = sc(X1, X2 - 1, z2);
                            /* Point 3 */
                            t[3] = tc(X1 - 1, X2 - 1, z2);
                            s[3] = sc(X1 - 1, X2 - 1, z2);
                            d01 = rz2 * sinqX2*df;
                            d02 = rz2*dq;
                            d13 = rz2*dq;
                            d23 = rz2 * sin(qX2m1) * df;
                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) {
                                fhead = (guess - try) / (d01 * slow0[index]);
                                guess = try;
                            }
                        }


                        /*** New Face Trials ***/
                        if (guess > 1.0e9) {
                            if (X1 > x1 + 1 && X1 < x2 - 1 && X2 > y1 + 1 && X2 < y2 - 1) {
                                /* Point 0 */
                                t[0] = tc(X1, X2 + 1, z2 - 1);
                                s[0] = sc(X1, X2 + 1, z2 - 1);
                                /* Point 1 */
                                t[1] = tc(X1 - 1, X2, z2 - 1);
                                s[1] = sc(X1 - 1, X2, z2 - 1);
                                /* Point 2 */
                                t[2] = tc(X1, X2, z2 - 1);
                                s[2] = sc(X1, X2, z2 - 1);
                                /* Point 3 */
                                t[3] = tc(X1 + 1, X2, z2 - 1);
                                s[3] = sc(X1 + 1, X2, z2 - 1);
                                /* Point 4 */
                                t[4] = tc(X1, X2 - 1, z2 - 1);
                                s[4] = sc(X1, X2 - 1, z2 - 1);
                                /* Point 5 */
                                s[5] = sc(X1, X2, z2);
                                d02 = rzm1*dq;
                                d12 = rzm1 * sinqX2*df;
                                d25 = h;
                                try = fdsphnf(t, s, d02, d12, d25);
                                if (try < guess) guess = try;
                            }
                        }

                        /*** Edge Trials ***/
                        /* Change in Depth */
                        try = tc(X1, X2, z2 - 1) + .5 * (sc(X1, X2, z2) + sc(X1, X2, z2 - 1)) * h;
                        if (try < guess) guess = try;
                        /* Change in East Longitude */
                        if (X1 < nx - 1 && time0[index + 1] < 1.e9) {
                            try = tc(X1 + 1, X2, z2) + .5 * (sc(X1, X2, z2) + sc(X1 + 1, X2, z2)) * rz2 * df*sinqX2;
                            if (try < guess) {
                                fhead = (guess - try) / (rz2 * df * sinqX2 * slow0[index]);
                                guess = try;
                            }
                        }
                        /* Change in West Longitude */
                        if (X1 > 0 && time0[index - 1] < 1.e9) {
                            try = tc(X1 - 1, X2, z2) + .5 * (sc(X1, X2, z2) + sc(X1 - 1, X2, z2)) * rz2 * df*sinqX2;
                            if (try < guess) {
                                fhead = (guess - try) / (rz2 * df * sinqX2 * slow0[index]);
                                guess = try;
                            }
                        }
                        /* Change in South Latitude */
                        if (X2 < ny - 1 && time0[index + nx] < 1.e9) {
                            try = tc(X1, X2 + 1, z2) + .5 * (sc(X1, X2, z2) + sc(X1, X2 + 1, z2)) * rz2*dq;
                            if (try < guess) {
                                fhead = (guess - try) / (rz2 * dq * slow0[index]);
                                guess = try;
                            }
                        }
                        /* Change in North Latitude */
                        if (X2 > 0 && time0[index - nx] < 1.e9) {
                            try = tc(X1, X2 - 1, z2) + .5 * (sc(X1, X2, z2) + sc(X1, X2 - 1, z2)) * rz2*dq;
                            if (try < guess) {
                                fhead = (guess - try) / (rz2 * dq * slow0[index]);
                                guess = try;
                            }
                        }
                        if (guess < time0[index]) {
                            time0[index] = guess;
                            if (fhead > headtest) headw[6]++;
                        }
                    }
                    if (z2 == nz - 1) dz2 = 0;
                    z2++;
                }
                /* End z2 if condition */
            }
            /* End growth loop */


            /* EAST (RIGHT) SIDE */
            for (igrow = 1; igrow <= iplus; igrow++) {
                if (dx2) {
                    ii = 0;
                    for (k = z1 + 1; k <= z2 - 1; k++) {
                        for (j = y1 + 1; j <= y2 - 1; j++) {
                            sort[ii].time = tc(x2 - 1, j, k);
                            sort[ii].i1 = j;
                            sort[ii].i2 = k;
                            ii++;
                        }
                    }
                    qsort((char *) sort, ii, sizeof (struct sorted), compar);

                    for (i = 0; i < ii; i++) {
                        X1 = sort[i].i1;
                        X2 = sort[i].i2;
                        index = X2 * nxy + X1 * nx + x2;
                        fhead = 0.;
                        guess = time0[index];

                        rX2 = rc(X2);
                        rXp1 = rc(X2 + 1);
                        rXm1 = rc(X2 - 1);
                        qX1 = qc(X1);
                        qX1p1 = qc(X1 + 1);
                        qX1m1 = qc(X1 - 1);
                        fx2 = fc(x2);
                        fx2m1 = fc(x2 - 1);

                        sinqX1 = sin(qX1);

                        /*** 3D Trials ***/
                        /* Top North corner (point 7) */
                        /*                       if(time0[index+nx] < 1.e9 && time0[index+nxy+nx] < 1.e9
                                                   && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<ny-1 ) { */
                        if (X2 < nz - 1 && X1 < ny - 1) {
                            if (time0[index + nx] < 1.e9 && time0[index + nxy + nx] < 1.e9
                                    && time0[index + nxy]) {
                                /* Point 1 */
                                t[0] = tc(x2 - 1, X1 + 1, X2 + 1);
                                s[0] = sc(x2 - 1, X1 + 1, X2 + 1);
                                r[0] = rXp1;
                                q[0] = qX1p1;
                                f[0] = fx2m1;
                                g[0] = -1;
                                n[0] = 1;
                                m[0] = -1;
                                /* Point 0 */
                                t[1] = tc(x2, X1 + 1, X2 + 1);
                                s[1] = sc(x2, X1 + 1, X2 + 1);
                                r[1] = rXp1;
                                q[1] = qX1p1;
                                f[1] = fx2;
                                g[1] = -1;
                                n[1] = 1;
                                m[1] = 1;
                                /* Point 3 */
                                t[2] = tc(x2, X1, X2 + 1);
                                s[2] = sc(x2, X1, X2 + 1);
                                r[2] = rXp1;
                                q[2] = qX1;
                                f[2] = fx2;
                                g[2] = -1;
                                n[2] = -1;
                                m[2] = 1;
                                /* Point 2 */
                                t[3] = tc(x2 - 1, X1, X2 + 1);
                                s[3] = sc(x2 - 1, X1, X2 + 1);
                                r[3] = rXp1;
                                q[3] = qX1;
                                f[3] = fx2m1;
                                g[3] = -1;
                                n[3] = -1;
                                m[3] = -1;
                                /* Point 5 */
                                t[4] = tc(x2 - 1, X1 + 1, X2);
                                s[4] = sc(x2 - 1, X1 + 1, X2);
                                r[4] = rX2;
                                q[4] = qX1p1;
                                f[4] = fx2m1;
                                g[4] = 1;
                                n[4] = 1;
                                m[4] = -1;
                                /* Point 4 */
                                t[5] = tc(x2, X1 + 1, X2);
                                s[5] = sc(x2, X1 + 1, X2);
                                r[5] = rX2;
                                q[5] = qX1p1;
                                f[5] = fx2;
                                g[5] = 1;
                                n[5] = 1;
                                m[5] = 1;
                                /* Point 6 */
                                t[6] = tc(x2 - 1, X1, X2);
                                s[6] = sc(x2 - 1, X1, X2);
                                r[6] = rX2;
                                q[6] = qX1;
                                f[6] = fx2m1;
                                g[6] = 1;
                                n[6] = -1;
                                m[6] = -1;
                                /* Point 7 */
                                s[7] = sc(x2, X1, X2);
                                r[7] = rX2;
                                q[7] = qX1;
                                f[7] = fx2;
                                g[7] = 1;
                                n[7] = -1;
                                m[7] = 1;

                                try = fdsph3d(t, s, r, q, f, g, n, m, h, dq, df);
                                if (try < guess) guess = try;
                            }

                        }

                        /* Top South corner (point 4) */
                        /*if(time0[index-nx] < 1.e9 && time0[index+nxy-nx] < 1.e9
                            && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) { */
                        if (X2 < nz - 1 && X1 > 0) {
                            if (time0[index - nx] < 1.e9 && time0[index + nxy - nx] < 1.e9
                                    && time0[index + nxy] < 1.e9) {
                                /* Point 1 */
                                t[0] = tc(x2 - 1, X1, X2 + 1);
                                s[0] = sc(x2 - 1, X1, X2 + 1);
                                r[0] = rXp1;
                                q[0] = qX1;
                                f[0] = fx2m1;
                                g[0] = -1;
                                n[0] = 1;
                                m[0] = -1;
                                /* Point 0 */
                                t[1] = tc(x2, X1, X2 + 1);
                                s[1] = sc(x2, X1, X2 + 1);
                                r[1] = rXp1;
                                q[1] = qX1;
                                f[1] = fx2;
                                g[1] = -1;
                                n[1] = 1;
                                m[1] = 1;
                                /* Point 3 */
                                t[2] = tc(x2, X1 - 1, X2 + 1);
                                s[2] = sc(x2, X1 - 1, X2 + 1);
                                r[2] = rXp1;
                                q[2] = qX1m1;
                                f[2] = fx2;
                                g[2] = -1;
                                n[2] = -1;
                                m[2] = 1;
                                /* Point 2 */
                                t[3] = tc(x2 - 1, X1 - 1, X2 + 1);
                                s[3] = sc(x2 - 1, X1 - 1, X2 + 1);
                                r[3] = rXp1;
                                q[3] = qX1m1;
                                f[3] = fx2m1;
                                g[3] = -1;
                                n[3] = -1;
                                m[3] = -1;
                                /* Point 5 */
                                t[4] = tc(x2 - 1, X1, X2);
                                s[4] = sc(x2 - 1, X1, X2);
                                r[4] = rX2;
                                q[4] = qX1;
                                f[4] = fx2m1;
                                g[4] = 1;
                                n[4] = 1;
                                m[4] = -1;
                                /* Point 7 */
                                t[5] = tc(x2, X1 - 1, X2);
                                s[5] = sc(x2, X1 - 1, X2);
                                r[5] = rX2;
                                q[5] = qX1m1;
                                f[5] = fx2;
                                g[5] = 1;
                                n[5] = -1;
                                m[5] = 1;
                                /* Point 6 */
                                t[6] = tc(x2 - 1, X1 - 1, X2);
                                s[6] = sc(x2 - 1, X1 - 1, X2);
                                r[6] = rX2;
                                q[6] = qX1m1;
                                f[6] = fx2m1;
                                g[6] = 1;
                                n[6] = -1;
                                m[6] = -1;
                                /* Point 4 */
                                s[7] = sc(x2, X1, X2);
                                r[7] = rX2;
                                q[7] = qX1;
                                f[7] = fx2;
                                g[7] = 1;
                                n[7] = 1;
                                m[7] = 1;

                                try = fdsph3d(t, s, r, q, f, g, n, m, h, dq, df);
                                if (try < guess) guess = try;
                            }
                        }

                        /*  Bottom North corner (point 3) */
                        if (X2 > 0 && X1 < ny - 1 && time0[index + nx] < 1.e9 && time0[index - nxy + nx] < 1.e9
                                && time0[index - nxy] < 1.e9) {
                            /* Point 1 */
                            t[0] = tc(x2 - 1, X1 + 1, X2);
                            s[0] = sc(x2 - 1, X1 + 1, X2);
                            r[0] = rX2;
                            q[0] = qX1p1;
                            f[0] = fx2m1;
                            g[0] = -1;
                            n[0] = 1;
                            m[0] = -1;
                            /* Point 0 */
                            t[1] = tc(x2, X1 + 1, X2);
                            s[1] = sc(x2, X1 + 1, X2);
                            r[1] = rX2;
                            q[1] = qX1p1;
                            f[1] = fx2;
                            g[1] = -1;
                            n[1] = 1;
                            m[1] = 1;
                            /* Point 2 */
                            t[2] = tc(x2 - 1, X1, X2);
                            s[2] = sc(x2 - 1, X1, X2);
                            r[2] = rX2;
                            q[2] = qX1;
                            f[2] = fx2m1;
                            g[2] = -1;
                            n[2] = -1;
                            m[2] = -1;
                            /* Point 5 */
                            t[3] = tc(x2 - 1, X1 + 1, X2 - 1);
                            s[3] = sc(x2 - 1, X1 + 1, X2 - 1);
                            r[3] = rXm1;
                            q[3] = qX1p1;
                            f[3] = fx2m1;
                            g[3] = 1;
                            n[3] = 1;
                            m[3] = -1;
                            /* Point 4 */
                            t[4] = tc(x2, X1 + 1, X2 - 1);
                            s[4] = sc(x2, X1 + 1, X2 - 1);
                            r[4] = rXm1;
                            q[4] = qX1p1;
                            f[4] = fx2;
                            g[4] = 1;
                            n[4] = 1;
                            m[4] = 1;
                            /* Point 7 */
                            t[5] = tc(x2, X1, X2 - 1);
                            s[5] = sc(x2, X1, X2 - 1);
                            r[5] = rXm1;
                            q[5] = qX1;
                            f[5] = fx2;
                            g[5] = 1;
                            n[5] = -1;
                            m[5] = 1;
                            /* Point 6 */
                            t[6] = tc(x2 - 1, X1, X2 - 1);
                            s[6] = sc(x2 - 1, X1, X2 - 1);
                            r[6] = rXm1;
                            q[6] = qX1;
                            f[6] = fx2m1;
                            g[6] = 1;
                            n[6] = -1;
                            m[6] = -1;
                            /* Point 3 */
                            s[7] = sc(x2, X1, X2);
                            r[7] = rX2;
                            q[7] = qX1;
                            f[7] = fx2;
                            g[7] = -1;
                            n[7] = -1;
                            m[7] = 1;

                            try = fdsph3d(t, s, r, q, f, g, n, m, h, dq, df);
                            if (try < guess) guess = try;

                        }

                        /* Bottom South corner (point 0) */
                        if (X2 > 0 && X1 > 0 && time0[index - nx] < 1.e9 && time0[index - nxy - nx] < 1.e9
                                && time0[index - nxy] < 1.e9) {
                            /* Point 1 */
                            t[0] = tc(x2 - 1, X1, X2);
                            s[0] = sc(x2 - 1, X1, X2);
                            r[0] = rX2;
                            q[0] = qX1;
                            f[0] = fx2m1;
                            g[0] = -1;
                            n[0] = 1;
                            m[0] = -1;
                            /* Point 3 */
                            t[1] = tc(x2, X1 - 1, X2);
                            s[1] = sc(x2, X1 - 1, X2);
                            r[1] = rX2;
                            q[1] = qX1m1;
                            f[1] = fx2;
                            g[1] = -1;
                            n[1] = -1;
                            m[1] = 1;
                            /* Point 2 */
                            t[2] = tc(x2 - 1, X1 - 1, X2);
                            s[2] = sc(x2 - 1, X1 - 1, X2);
                            r[2] = rX2;
                            q[2] = qX1m1;
                            f[2] = fx2m1;
                            g[2] = -1;
                            n[2] = -1;
                            m[2] = -1;
                            /* Point 5 */
                            t[3] = tc(x2 - 1, X1, X2 - 1);
                            s[3] = sc(x2 - 1, X1, X2 - 1);
                            r[3] = rXm1;
                            q[3] = qX1;
                            f[3] = fx2m1;
                            g[3] = 1;
                            n[3] = 1;
                            m[3] = -1;
                            /* Point 4 */
                            t[4] = tc(x2, X1, X2 - 1);
                            s[4] = sc(x2, X1, X2 - 1);
                            r[4] = rXm1;
                            q[4] = qX1;
                            f[4] = fx2;
                            g[4] = 1;
                            n[4] = 1;
                            m[4] = 1;
                            /* Point 7 */
                            t[5] = tc(x2, X1 - 1, X2 - 1);
                            s[5] = sc(x2, X1 - 1, X2 - 1);
                            r[5] = rXm1;
                            q[5] = qX1m1;
                            f[5] = fx2;
                            g[5] = 1;
                            n[5] = -1;
                            m[5] = 1;
                            /* Point 6 */
                            t[6] = tc(x2 - 1, X1 - 1, X2 - 1);
                            s[6] = sc(x2 - 1, X1 - 1, X2 - 1);
                            r[6] = rXm1;
                            q[6] = qX1m1;
                            f[6] = fx2m1;
                            g[6] = 1;
                            n[6] = -1;
                            m[6] = -1;
                            /* Point 0 */
                            s[7] = sc(x2, X1, X2);
                            r[7] = rX2;
                            q[7] = qX1;
                            f[7] = fx2;
                            g[7] = -1;
                            n[7] = 1;
                            m[7] = 1;
                            try = fdsph3d(t, s, r, q, f, g, n, m, h, dq, df);
                            if (try < guess) guess = try;

                        }

                        /*** New Edge Trials ***/
                        if (guess > 1.0e9) {
                            /* Top edge */
                            /*if(time0[index+nxy] < 1.e9 && X2<nz-1 && X1>y1+1 && X1<y2-1 )  {*/
                            if (X2 < nz - 1 && X1 > y1 + 1 && X1 < y2 - 1) {
                                if (time0[index + nxy] < 1.e9) {
                                    /* Point 0 */
                                    t[0] = tc(x2 - 1, X1 - 1, X2 + 1);
                                    s[0] = sc(x2 - 1, X1 - 1, X2 + 1);
                                    /* Point 1 */
                                    t[1] = tc(x2 - 1, X1, X2 + 1);
                                    s[1] = sc(x2 - 1, X1, X2 + 1);
                                    /* Point 2 */
                                    t[2] = tc(x2 - 1, X1, X2);
                                    s[2] = sc(x2 - 1, X1, X2);
                                    /* Point 3 */
                                    t[3] = tc(x2 - 1, X1 + 1, X2 + 1);
                                    s[3] = sc(x2 - 1, X1 + 1, X2 + 1);
                                    /* Point 4 */
                                    t[4] = tc(x2, X1, X2 + 1);
                                    s[4] = sc(x2, X1, X2 + 1);
                                    /* Point 5 */
                                    s[5] = sc(x2, X1, X2);

                                    d01 = rXp1*dq;
                                    d12 = h;
                                    d14 = rXp1 * sin(qX1) * df;
                                    d25 = rX2 * sin(qX1) * df;
                                    d45 = h;

                                    try = fdsphne(t, s, d01, d12, d14, d25, d45);
                                    if (try < guess) guess = try;
                                }
                            }
                            /* Bottom Edge */
                            if (X2 > 0 && X1 > y1 + 1 && X1 < y2 - 1 && time0[index - nxy] < 1.e9) {
                                /* Point 0 */
                                t[0] = tc(x2 - 1, X1 + 1, X2 - 1);
                                s[0] = sc(x2 - 1, X1 + 1, X2 - 1);
                                /* Point 1 */
                                t[1] = tc(x2 - 1, X1, X2 - 1);
                                s[1] = sc(x2 - 1, X1, X2 - 1);
                                /* Point 2 */
                                t[2] = tc(x2 - 1, X1, X2);
                                s[2] = sc(x2 - 1, X1, X2);
                                /* Point 3 */
                                t[3] = tc(x2 - 1, X1 - 1, X2 - 1);
                                s[3] = sc(x2 - 1, X1 - 1, X2 - 1);
                                /* Point 4 */
                                t[4] = tc(x2, X1, X2 - 1);
                                s[4] = sc(x2, X1, X2 - 1);
                                /* Point 5 */
                                s[5] = sc(x2, X1, X2);

                                d01 = rXm1*dq;
                                d12 = h;
                                d14 = rXm1 * sin(qX1) * df;
                                d25 = rX2 * sin(qX1) * df;
                                d45 = h;

                                try = fdsphne(t, s, d01, d12, d14, d25, d45);
                                if (try < guess) guess = try;

                            }
                            /* Northern Edge */
                            if (X1 < ny - 1 && X2 > z1 + 1 && X2 < z2 - 1 && time0[index + nx] < 1.e9) {
                                /* Point 0 */
                                t[0] = tc(x2 - 1, X1 + 1, X2 + 1);
                                s[0] = sc(x2 - 1, X1 + 1, X2 + 1);
                                /* Point 1 */
                                t[1] = tc(x2 - 1, X1 + 1, X2);
                                s[1] = sc(x2 - 1, X1 + 1, X2);
                                /* Point 2 */
                                t[2] = tc(x2 - 1, X1, X2);
                                s[2] = sc(x2 - 1, X1, X2);
                                /* Point 3 */
                                t[3] = tc(x2 - 1, X1 + 1, X2 - 1);
                                s[3] = sc(x2 - 1, X1 + 1, X2 - 1);
                                /* Point 4 */
                                t[4] = tc(x2, X1 + 1, X2);
                                s[4] = sc(x2, X1 + 1, X2);
                                /* Point 5 */
                                s[5] = sc(x2, X1, X2);

                                d01 = h;
                                d12 = rX2*dq;
                                d14 = rX2 * sin(qX1p1) * df;
                                d25 = rX2 * sin(qX1) * df;
                                d45 = rX2*dq;

                                try = fdsphne(t, s, d01, d12, d14, d25, d45);
                                if (try < guess) guess = try;
                            }
                            /* Southern Edge */
                            if (X1 > 0 && X2 > z1 + 1 && X2 < z2 - 1 && time0[index - nx] < 1.e9) {
                                /* Point 0 */
                                t[0] = tc(x2 - 1, X1 - 1, X2 - 1);
                                s[0] = sc(x2 - 1, X1 - 1, X2 - 1);
                                /* Point 1 */
                                t[1] = tc(x2 - 1, X1 - 1, X2);
                                s[1] = sc(x2 - 1, X1 - 1, X2);
                                /* Point 2 */
                                t[2] = tc(x2 - 1, X1, X2);
                                s[2] = sc(x2 - 1, X1, X2);
                                /* Point 3 */
                                t[3] = tc(x2 - 1, X1 - 1, X2 + 1);
                                s[3] = sc(x2 - 1, X1 - 1, X2 + 1);
                                /* Point 4 */
                                t[4] = tc(x2, X1 - 1, X2);
                                s[4] = sc(x2, X1 - 1, X2);
                                /* Point 5 */
                                s[5] = sc(x2, X1, X2);

                                d01 = h;
                                d12 = rX2*dq;
                                d14 = rX2 * sin(qX1m1) * df;
                                d25 = rX2 * sin(qX1) * df;
                                d45 = rX2*dq;

                                try = fdsphne(t, s, d01, d12, d14, d25, d45);
                                if (try < guess) guess = try;

                            }
                        }

                        /*** Two-D Trials ***/
                        /* East-North */
                        if (X1 > 0 && time0[index - nx] < 1.e9) {
                            s[0] = sc(x2, X1, X2);
                            /* Point 1 */
                            t[1] = tc(x2, X1 - 1, X2);
                            s[1] = sc(x2, X1 - 1, X2);
                            /* Point 2 */
                            t[2] = tc(x2 - 1, X1, X2);
                            s[2] = sc(x2 - 1, X1, X2);
                            /* Point 3 */
                            t[3] = tc(x2 - 1, X1 - 1, X2);
                            s[3] = sc(x2 - 1, X1 - 1, X2);

                            d01 = rX2*dq;
                            d02 = rX2 * sinqX1*df;
                            d13 = rX2 * sin(qX1m1) * df;
                            d23 = rX2*dq;

                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) guess = try;
                        }
                        /* East-South */
                        if (X1 < ny - 1 && time0[index + nx] < 1.e9) {
                            s[0] = sc(x2, X1, X2);
                            /* Point 1 */
                            t[1] = tc(x2, X1 + 1, X2);
                            s[1] = sc(x2, X1 + 1, X2);
                            /* Point 2 */
                            t[2] = tc(x2 - 1, X1, X2);
                            s[2] = sc(x2 - 1, X1, X2);
                            /* Point 3 */
                            t[3] = tc(x2 - 1, X1 + 1, X2);
                            s[3] = sc(x2 - 1, X1 + 1, X2);

                            d01 = rX2*dq;
                            d02 = rX2 * sinqX1*df;
                            d13 = rX2 * sin(qX1p1) * df;
                            d23 = rX2*dq;

                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) guess = try;
                        }
                        /* East-Top */
                        if (X2 > 0 && time0[index - nxy] < 1.e9) {
                            s[0] = sc(x2, X1, X2);
                            /* Point 1 */
                            t[1] = tc(x2, X1, X2 - 1);
                            s[1] = sc(x2, X1, X2 - 1);
                            /* Point 2 */
                            t[2] = tc(x2 - 1, X1, X2);
                            s[2] = sc(x2 - 1, X1, X2);
                            /* Point 3 */
                            t[3] = tc(x2 - 1, X1, X2 - 1);
                            s[3] = sc(x2 - 1, X1, X2 - 1);

                            d01 = h;
                            d02 = rX2 * sinqX1*df;
                            d13 = rXm1 * sinqX1*df;
                            d23 = h;

                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) guess = try;
                        }
                        /* East-Bottom */
                        /*if(time0[index+nxy] < 1.e9 && X2<nz-1 )  {*/
                        if (X2 < nz - 1) {
                            if (time0[index + nxy] < 1.e9) {
                                s[0] = sc(x2, X1, X2);
                                /* Point 1 */
                                t[1] = tc(x2, X1, X2 + 1);
                                s[1] = sc(x2, X1, X2 + 1);
                                /* Point 2 */
                                t[2] = tc(x2 - 1, X1, X2);
                                s[2] = sc(x2 - 1, X1, X2);
                                /* Point 3 */
                                t[3] = tc(x2 - 1, X1, X2 + 1);
                                s[3] = sc(x2 - 1, X1, X2 + 1);

                                d01 = h;
                                d02 = rX2 * sinqX1*df;
                                d13 = rXp1 * sinqX1*df;
                                d23 = h;

                                try = fdsph2d(t, s, d01, d02, d13, d23);
                                if (try < guess) guess = try;
                            }
                        }
                        /* Parallel - North Bottom */
                        /*if(time0[index-nx] < 1.e9 && time0[index+nxy-nx] < 1.e9
                           && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) { */
                        if (X2 < nz - 1 && X1 > 0) {
                            if (time0[index - nx] < 1.e9 && time0[index + nxy - nx] < 1.e9
                                    && time0[index + nxy] < 1.e9) {
                                s[0] = sc(x2, X1, X2);
                                /* Point 1 */
                                t[1] = tc(x2, X1, X2 + 1);
                                s[1] = sc(x2, X1, X2 + 1);
                                /* Point 2 */
                                t[2] = tc(x2, X1 - 1, X2);
                                s[2] = sc(x2, X1 - 1, X2);
                                /* Point 3 */
                                t[3] = tc(x2, X1 - 1, X2 + 1);
                                s[3] = sc(x2, X1 - 1, X2 + 1);
                                d01 = h;
                                d02 = rX2*dq;
                                d13 = rXp1*dq;
                                d23 = h;
                                try = fdsph2d(t, s, d01, d02, d13, d23);
                                if (try < guess) {
                                    fhead = (guess - try) / (d01 * slow0[index]);
                                    guess = try;
                                }
                            }
                        }
                        /* Parallel - North Top */
                        if (X2 > 0 && X1 > 0 && time0[index - nx] < 1.e9 && time0[index - nxy - nx] < 1.e9
                                && time0[index - nxy] < 1.e9) {
                            s[0] = sc(x2, X1, X2);
                            /* Point 1 */
                            t[1] = tc(x2, X1, X2 - 1);
                            s[1] = sc(x2, X1, X2 - 1);
                            /* Point 2 */
                            t[2] = tc(x2, X1 - 1, X2);
                            s[2] = sc(x2, X1 - 1, X2);
                            /* Point 3 */
                            t[3] = tc(x2, X1 - 1, X2 - 1);
                            s[3] = sc(x2, X1 - 1, X2 - 1);
                            d01 = h;
                            d02 = rX2*dq;
                            d13 = rXm1*dq;
                            d23 = h;
                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) {
                                fhead = (guess - try) / (d01 * slow0[index]);
                                guess = try;
                            }
                        }
                        /* Parallel - South Bottom */
                        /*if(time0[index+nx] < 1.e9 && time0[index+nxy+nx] < 1.e9
                            && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<ny-1 ) { */
                        if (X2 < nz - 1 && X1 < ny - 1) {
                            if (time0[index + nx] < 1.e9 && time0[index + nxy + nx] < 1.e9
                                    && time0[index + nxy] < 1.e9) {
                                s[0] = sc(x2, X1, X2);
                                /* Point 1 */
                                t[1] = tc(x2, X1, X2 + 1);
                                s[1] = sc(x2, X1, X2 + 1);
                                /* Point 2 */
                                t[2] = tc(x2, X1 + 1, X2);
                                s[2] = sc(x2, X1 + 1, X2);
                                /* Point 3 */
                                t[3] = tc(x2, X1 + 1, X2 + 1);
                                s[3] = sc(x2, X1 + 1, X2 + 1);
                                d01 = h;
                                d02 = rX2*dq;
                                d13 = rXp1*dq;
                                d23 = h;
                                try = fdsph2d(t, s, d01, d02, d13, d23);
                                if (try < guess) {
                                    fhead = (guess - try) / (d01 * slow0[index]);
                                    guess = try;
                                }
                            }
                        }
                        /* Parallel - South Top */
                        if (X2 > 0 && X1 < ny - 1 && time0[index + nx] < 1.e9 && time0[index - nxy + nx] < 1.e9
                                && time0[index - nxy] < 1.e9) {
                            s[0] = sc(x2, X1, X2);
                            /* Point 1 */
                            t[1] = tc(x2, X1, X2 - 1);
                            s[1] = sc(x2, X1, X2 - 1);
                            /* Point 2 */
                            t[2] = tc(x2, X1 + 1, X2);
                            s[2] = sc(x2, X1 + 1, X2);
                            /* Point 3 */
                            t[3] = tc(x2, X1 + 1, X2 - 1);
                            s[3] = sc(x2, X1 + 1, X2 - 1);
                            d01 = h;
                            d02 = rX2*dq;
                            d13 = rXm1*dq;
                            d23 = h;
                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) {
                                fhead = (guess - try) / (d01 * slow0[index]);
                                guess = try;
                            }
                        }

                        /*** New Face Trials ***/
                        if (guess > 1.0e9) {
                            if (X1 > y1 + 1 && X1 < y2 - 1 && X2 > z1 + 1 && X2 < z2 - 1) {
                                /* Point 0 */
                                t[0] = tc(x2 - 1, X1 + 1, X2);
                                s[0] = sc(x2 - 1, X1 + 1, X2);
                                /* Point 1 */
                                t[1] = tc(x2 - 1, X1, X2 - 1);
                                s[1] = sc(x2 - 1, X1, X2 - 1);
                                /* Point 2 */
                                t[2] = tc(x2 - 1, X1, X2);
                                s[2] = sc(x2 - 1, X1, X2);
                                /* Point 3 */
                                t[3] = tc(x2 - 1, X1, X2 + 1);
                                s[3] = sc(x2 - 1, X1, X2 + 1);
                                /* Point 4 */
                                t[4] = tc(x2 - 1, X1 - 1, X2);
                                s[4] = sc(x2 - 1, X1 - 1, X2);
                                /* Point 5 */
                                s[5] = sc(x2, X1, X2);
                                d02 = rX2*dq;
                                d12 = h;
                                d25 = rX2 * sinqX1*df;
                                try = fdsphnf(t, s, d02, d12, d25);
                                if (try < guess) guess = try;

                            }
                        }

                        /*** Edge Trials ***/
                        /* Change in Longitude */
                        try = tc(x2 - 1, X1, X2) + .5 * (sc(x2, X1, X2) + sc(x2 - 1, X1, X2)) * rX2 * sinqX1*df;
                        if (try < guess) guess = try;
                        /* Change in South Latitude */
                        if (X1 < ny - 1 && time0[index + nx] < 1.e9) {
                            try = tc(x2, X1 + 1, X2) + .5 * (sc(x2, X1, X2) + sc(x2, X1 + 1, X2)) * rX2*dq;
                            if (try < guess) {
                                fhead = (guess - try) / (rX2 * dq * slow0[index]);
                                guess = try;
                            }
                        }
                        /* Change in North Latitude */
                        if (X1 > 0 && time0[index - nx] < 1.e9) {
                            try = tc(x2, X1 - 1, X2) + .5 * (sc(x2, X1, X2) + sc(x2, X1 - 1, X2)) * rX2*dq;
                            if (try < guess) {
                                fhead = (guess - try) / (rX2 * dq * slow0[index]);
                                guess = try;
                            }
                        }
                        /* Change in depth (increasing) */
                        /*if ( time0[index+nxy]<1.e9 && X2<nz-1 )  */
                        if (X2 < nz - 1) {
                            if (time0[index + nxy] < 1.e9) {
                                try = tc(x2, X1, X2 + 1) + .5 * (sc(x2, X1, X2) + sc(x2, X1, X2 + 1)) * h;
                                if (try < guess) {
                                    fhead = (guess - try) / (h * slow0[index]);
                                    guess = try;
                                }
                            }
                        }
                        /* Change in depth (decreasing) */
                        /*if ( time0[index-nxy]<1.e9 && X2>0 )  { */
                        if (X2 > 0) {
                            if (time0[index - nxy] < 1.e9) {
                                try = tc(x2, X1, X2 - 1) + .5 * (sc(x2, X1, X2) + sc(x2, X1, X2 - 1)) * h;
                                if (try < guess) {
                                    fhead = (guess - try) / (h * slow0[index]);
                                    guess = try;
                                }
                            }
                        }
                        if (guess < time0[index]) {
                            time0[index] = guess;
                            if (fhead > headtest) headw[2]++;
                        }
                    }
                    if (x2 == nx - 1) dx2 = 0;
                    x2++;
                }
                /* End x2 if condition */
            }
            /* End growth loop */


            /* WEST (LEFT) SIDE */
            for (igrow = 1; igrow <= iminus; igrow++) {
                if (dx1) {
                    ii = 0;
                    for (k = z1 + 1; k <= z2 - 1; k++) {
                        for (j = y1 + 1; j <= y2 - 1; j++) {
                            sort[ii].time = tc(x1 + 1, j, k);
                            sort[ii].i1 = j;
                            sort[ii].i2 = k;
                            ii++;
                        }
                    }
                    qsort((char *) sort, ii, sizeof (struct sorted), compar);
                    for (i = 0; i < ii; i++) {
                        X1 = sort[i].i1;
                        X2 = sort[i].i2;
                        index = X2 * nxy + X1 * nx + x1;
                        fhead = 0.;
                        guess = time0[index];

                        rX2 = rc(X2);
                        rXp1 = rc(X2 + 1);
                        rXm1 = rc(X2 - 1);
                        qX1 = qc(X1);
                        qX1p1 = qc(X1 + 1);
                        qX1m1 = qc(X1 - 1);
                        fx1 = fc(x1);
                        fx1p1 = fc(x1 + 1);
                        fx1m1 = fc(x1 - 1);

                        sinqX1 = sin(qX1);

                        /*** 3D Trials ***/
                        /* Top North corner (point 6) */
                        /*if(time0[index+nx] < 1.e9 && time0[index+nxy+nx] < 1.e9
                           && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<ny-1 ) { */
                        if (X2 < nz - 1 && X1 < ny - 1) {
                            if (time0[index + nx] < 1.e9 && time0[index + nxy + nx] < 1.e9
                                    && time0[index + nxy] < 1.e9) {
                                /* Point 0 */
                                t[0] = tc(x1 + 1, X1 + 1, X2 + 1);
                                s[0] = sc(x1 + 1, X1 + 1, X2 + 1);
                                r[0] = rXp1;
                                q[0] = qX1p1;
                                f[0] = fx1p1;
                                g[0] = -1;
                                n[0] = 1;
                                m[0] = 1;
                                /* Point 1 */
                                t[1] = tc(x1, X1 + 1, X2 + 1);
                                s[1] = sc(x1, X1 + 1, X2 + 1);
                                r[1] = rXp1;
                                q[1] = qX1p1;
                                f[1] = fx1;
                                g[1] = -1;
                                n[1] = 1;
                                m[1] = -1;
                                /* Point 2 */
                                t[2] = tc(x1, X1, X2 + 1);
                                s[2] = sc(x1, X1, X2 + 1);
                                r[2] = rXp1;
                                q[2] = qX1;
                                f[2] = fx1;
                                g[2] = -1;
                                n[2] = -1;
                                m[2] = -1;
                                /* Point 3 */
                                t[3] = tc(x1 + 1, X1, X2 + 1);
                                s[3] = sc(x1 + 1, X1, X2 + 1);
                                r[3] = rXp1;
                                q[3] = qX1;
                                f[3] = fx1p1;
                                g[3] = -1;
                                n[3] = -1;
                                m[3] = 1;
                                /* Point 4 */
                                t[4] = tc(x1 + 1, X1 + 1, X2);
                                s[4] = sc(x1 + 1, X1 + 1, X2);
                                r[4] = rX2;
                                q[4] = qX1p1;
                                f[4] = fx1p1;
                                g[4] = 1;
                                n[4] = 1;
                                m[4] = 1;
                                /* Point 5 */
                                t[5] = tc(x1, X1 + 1, X2);
                                s[5] = sc(x1, X1 + 1, X2);
                                r[5] = rX2;
                                q[5] = qX1p1;
                                f[5] = fx1;
                                g[5] = 1;
                                n[5] = 1;
                                m[5] = -1;
                                /* Point 7 */
                                t[6] = tc(x1 + 1, X1, X2);
                                s[6] = sc(x1 + 1, X1, X2);
                                r[6] = rX2;
                                q[6] = qX1;
                                f[6] = fx1p1;
                                g[6] = 1;
                                n[6] = -1;
                                m[6] = 1;
                                /* Point 6 */
                                s[7] = sc(x1, X1, X2);
                                r[7] = rX2;
                                q[7] = qX1;
                                f[7] = fx1;
                                g[7] = 1;
                                n[7] = -1;
                                m[7] = -1;

                                try = fdsph3d(t, s, r, q, f, g, n, m, h, dq, df);
                                if (try < guess) guess = try;
                            }

                        }

                        /* Top South corner (point 5) */
                        /*if(time0[index-nx] < 1.e9 && time0[index+nxy-nx] < 1.e9
                           && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) { */
                        if (X2 < nz - 1 && X1 > 0) {
                            if (time0[index - nx] < 1.e9 && time0[index + nxy - nx] < 1.e9
                                    && time0[index + nxy] < 1.e9) {
                                /* Point 0 */
                                t[0] = tc(x1 + 1, X1, X2 + 1);
                                s[0] = sc(x1 + 1, X1, X2 + 1);
                                r[0] = rXp1;
                                q[0] = qX1;
                                f[0] = fx1p1;
                                g[0] = -1;
                                n[0] = 1;
                                m[0] = 1;
                                /* Point 1 */
                                t[1] = tc(x1, X1, X2 + 1);
                                s[1] = sc(x1, X1, X2 + 1);
                                r[1] = rXp1;
                                q[1] = qX1;
                                f[1] = fx1;
                                g[1] = -1;
                                n[1] = 1;
                                m[1] = -1;
                                /* Point 2 */
                                t[2] = tc(x1, X1 - 1, X2 + 1);
                                s[2] = sc(x1, X1 - 1, X2 + 1);
                                r[2] = rXp1;
                                q[2] = qX1m1;
                                f[2] = fx1;
                                g[2] = -1;
                                n[2] = -1;
                                m[2] = -1;
                                /* Point 3 */
                                t[3] = tc(x1 + 1, X1 - 1, X2 + 1);
                                s[3] = sc(x1 + 1, X1 - 1, X2 + 1);
                                r[3] = rXp1;
                                q[3] = qX1m1;
                                f[3] = fx1p1;
                                g[3] = -1;
                                n[3] = -1;
                                m[3] = 1;
                                /* Point 4 */
                                t[4] = tc(x1 + 1, X1, X2);
                                s[4] = sc(x1 + 1, X1, X2);
                                r[4] = rX2;
                                q[4] = qX1;
                                f[4] = fx1p1;
                                g[4] = 1;
                                n[4] = 1;
                                m[4] = 1;
                                /* Point 6 */
                                t[5] = tc(x1, X1 - 1, X2);
                                s[5] = sc(x1, X1 - 1, X2);
                                r[5] = rX2;
                                q[5] = qX1m1;
                                f[5] = fx1;
                                g[5] = 1;
                                n[5] = -1;
                                m[5] = -1;
                                /* Point 7 */
                                t[6] = tc(x1 + 1, X1 - 1, X2);
                                s[6] = sc(x1 + 1, X1 - 1, X2);
                                r[6] = rX2;
                                q[6] = qX1m1;
                                f[6] = fx1p1;
                                g[6] = 1;
                                n[6] = -1;
                                m[6] = 1;
                                /* Point 5 */
                                s[7] = sc(x1, X1, X2);
                                r[7] = rX2;
                                q[7] = qX1;
                                f[7] = fx1;
                                g[7] = 1;
                                n[7] = 1;
                                m[7] = -1;

                                try = fdsph3d(t, s, r, q, f, g, n, m, h, dq, df);
                                if (try < guess) guess = try;
                            }
                        }

                        /*  Bottom North corner (point 2) */
                        if (X2 > 0 && X1 < ny - 1 && time0[index + nx] < 1.e9 && time0[index - nxy + nx] < 1.e9
                                && time0[index - nxy] < 1.e9) {
                            /* Point 0 */
                            t[0] = tc(x1 + 1, X1 + 1, X2);
                            s[0] = sc(x1 + 1, X1 + 1, X2);
                            r[0] = rX2;
                            q[0] = qX1p1;
                            f[0] = fx1p1;
                            g[0] = -1;
                            n[0] = 1;
                            m[0] = 1;
                            /* Point 1 */
                            t[1] = tc(x1, X1 + 1, X2);
                            s[1] = sc(x1, X1 + 1, X2);
                            r[1] = rX2;
                            q[1] = qX1p1;
                            f[1] = fx1;
                            g[1] = -1;
                            n[1] = 1;
                            m[1] = -1;
                            /* Point 3 */
                            t[2] = tc(x1 + 1, X1, X2);
                            s[2] = sc(x1 + 1, X1, X2);
                            r[2] = rX2;
                            q[2] = qX1;
                            f[2] = fx1p1;
                            g[2] = -1;
                            n[2] = -1;
                            m[2] = 1;
                            /* Point 4 */
                            t[3] = tc(x1 + 1, X1 + 1, X2 - 1);
                            s[3] = sc(x1 + 1, X1 + 1, X2 - 1);
                            r[3] = rXm1;
                            q[3] = qX1p1;
                            f[3] = fx1p1;
                            g[3] = 1;
                            n[3] = 1;
                            m[3] = 1;
                            /* Point 5 */
                            t[4] = tc(x1, X1 + 1, X2 - 1);
                            s[4] = sc(x1, X1 + 1, X2 - 1);
                            r[4] = rXm1;
                            q[4] = qX1p1;
                            f[4] = fx1;
                            g[4] = 1;
                            n[4] = 1;
                            m[4] = -1;
                            /* Point 6 */
                            t[5] = tc(x1, X1, X2 - 1);
                            s[5] = sc(x1, X1, X2 - 1);
                            r[5] = rXm1;
                            q[5] = qX1;
                            f[5] = fx1;
                            g[5] = 1;
                            n[5] = -1;
                            m[5] = -1;
                            /* Point 7 */
                            t[6] = tc(x1 + 1, X1, X2 - 1);
                            s[6] = sc(x1 + 1, X1, X2 - 1);
                            r[6] = rXm1;
                            q[6] = qX1;
                            f[6] = fx1p1;
                            g[6] = 1;
                            n[6] = -1;
                            m[6] = 1;
                            /* Point 2 */
                            s[7] = sc(x1, X1, X2);
                            r[7] = rX2;
                            q[7] = qX1;
                            f[7] = fx1;
                            g[7] = -1;
                            n[7] = -1;
                            m[7] = -1;

                            try = fdsph3d(t, s, r, q, f, g, n, m, h, dq, df);
                            if (try < guess) guess = try;
                        }

                        /* Bottom South corner (point 1) */
                        if (X2 > 0 && X1 > 0 && time0[index - nx] < 1.e9 && time0[index - nxy - nx] < 1.e9
                                && time0[index - nxy] < 1.e9) {
                            /* Point 0 */
                            t[0] = tc(x1 + 1, X1, X2);
                            s[0] = sc(x1 + 1, X1, X2);
                            r[0] = rX2;
                            q[0] = qX1;
                            f[0] = fx1p1;
                            g[0] = -1;
                            n[0] = 1;
                            m[0] = 1;
                            /* Point 2 */
                            t[1] = tc(x1, X1 - 1, X2);
                            s[1] = sc(x1, X1 - 1, X2);
                            r[1] = rX2;
                            q[1] = qX1m1;
                            f[1] = fx1;
                            g[1] = -1;
                            n[1] = -1;
                            m[1] = -1;
                            /* Point 3 */
                            t[2] = tc(x1 + 1, X1 - 1, X2);
                            s[2] = sc(x1 + 1, X1 - 1, X2);
                            r[2] = rX2;
                            q[2] = qX1m1;
                            f[2] = fx1p1;
                            g[2] = -1;
                            n[2] = -1;
                            m[2] = 1;
                            /* Point 4 */
                            t[3] = tc(x1 + 1, X1, X2 - 1);
                            s[3] = sc(x1 + 1, X1, X2 - 1);
                            r[3] = rXm1;
                            q[3] = qX1;
                            f[3] = fx1p1;
                            g[3] = 1;
                            n[3] = 1;
                            m[3] = 1;
                            /* Point 5 */
                            t[4] = tc(x1, X1, X2 - 1);
                            s[4] = sc(x1, X1, X2 - 1);
                            r[4] = rXm1;
                            q[4] = qX1;
                            f[4] = fx1;
                            g[4] = 1;
                            n[4] = 1;
                            m[4] = -1;
                            /* Point 6 */
                            t[5] = tc(x1, X1 - 1, X2 - 1);
                            s[5] = sc(x1, X1 - 1, X2 - 1);
                            r[5] = rXm1;
                            q[5] = qX1m1;
                            f[5] = fx1;
                            g[5] = 1;
                            n[5] = -1;
                            m[5] = -1;
                            /* Point 7 */
                            t[6] = tc(x1 + 1, X1 - 1, X2 - 1);
                            s[6] = sc(x1 + 1, X1 - 1, X2 - 1);
                            r[6] = rXm1;
                            q[6] = qX1m1;
                            f[6] = fx1p1;
                            g[6] = 1;
                            n[6] = -1;
                            m[6] = 1;
                            /* Point 1 */
                            s[7] = sc(x1, X1, X2);
                            r[7] = rX2;
                            q[7] = qX1;
                            f[7] = fx1;
                            g[7] = -1;
                            n[7] = 1;
                            m[7] = -1;
                            try = fdsph3d(t, s, r, q, f, g, n, m, h, dq, df);
                            if (try < guess) guess = try;
                        }

                        /*** New Edge Trials ***/
                        if (guess > 1.0e9) {
                            /* Top edge */
                            /*if(time0[index+nxy] < 1.e9 && X2<nz-1 && X1>y1+1 && X1<y2-1 )  {*/
                            if (X2 < nz - 1 && X1 > y1 + 1 && X1 < y2 - 1) {
                                if (time0[index + nxy] < 1.e9) {
                                    /* Point 0 */
                                    t[0] = tc(x1 + 1, X1 + 1, X2 + 1);
                                    s[0] = sc(x1 + 1, X1 + 1, X2 + 1);
                                    /* Point 1 */
                                    t[1] = tc(x1 + 1, X1, X2 + 1);
                                    s[1] = sc(x1 + 1, X1, X2 + 1);
                                    /* Point 2 */
                                    t[2] = tc(x1 + 1, X1, X2);
                                    s[2] = sc(x1 + 1, X1, X2);
                                    /* Point 3 */
                                    t[3] = tc(x1 + 1, X1 - 1, X2 + 1);
                                    s[3] = sc(x1 + 1, X1 - 1, X2 + 1);
                                    /* Point 4 */
                                    t[4] = tc(x1, X1, X2 + 1);
                                    s[4] = sc(x1, X1, X2 + 1);
                                    /* Point 5 */
                                    s[5] = sc(x1, X1, X2);

                                    d01 = rXp1*dq;
                                    d12 = h;
                                    d14 = rXp1 * sin(qX1) * df;
                                    d25 = rX2 * sin(qX1) * df;
                                    d45 = h;

                                    try = fdsphne(t, s, d01, d12, d14, d25, d45);
                                    if (try < guess) guess = try;
                                }
                            }
                            /* Bottom Edge */
                            if (X2 > 0 && X1 > y1 + 1 && X1 < y2 - 1 && time0[index - nxy] < 1.e9) {
                                /* Point 0 */
                                t[0] = tc(x1 + 1, X1 - 1, X2 - 1);
                                s[0] = sc(x1 + 1, X1 - 1, X2 - 1);
                                /* Point 1 */
                                t[1] = tc(x1 + 1, X1, X2 - 1);
                                s[1] = sc(x1 + 1, X1, X2 - 1);
                                /* Point 2 */
                                t[2] = tc(x1 + 1, X1, X2);
                                s[2] = sc(x1 + 1, X1, X2);
                                /* Point 3 */
                                t[3] = tc(x1 + 1, X1 + 1, X2 - 1);
                                s[3] = sc(x1 + 1, X1 + 1, X2 - 1);
                                /* Point 4 */
                                t[4] = tc(x1, X1, X2 - 1);
                                s[4] = sc(x1, X1, X2 - 1);
                                /* Point 5 */
                                s[5] = sc(x1, X1, X2);

                                d01 = rXm1*dq;
                                d12 = h;
                                d14 = rXm1 * sin(qX1) * df;
                                d25 = rX2 * sin(qX1) * df;
                                d45 = h;

                                try = fdsphne(t, s, d01, d12, d14, d25, d45);
                                if (try < guess) guess = try;
                            }
                            /* Northern Edge */
                            if (X1 < ny - 1 && X2 > z1 + 1 && X2 < z2 - 1 && time0[index + nx] < 1.e9) {
                                /* Point 0 */
                                t[0] = tc(x1 + 1, X1 + 1, X2 - 1);
                                s[0] = sc(x1 + 1, X1 + 1, X2 - 1);
                                /* Point 1 */
                                t[1] = tc(x1 + 1, X1 + 1, X2);
                                s[1] = sc(x1 + 1, X1 + 1, X2);
                                /* Point 2 */
                                t[2] = tc(x1 + 1, X1, X2);
                                s[2] = sc(x1 + 1, X1, X2);
                                /* Point 3 */
                                t[3] = tc(x1 + 1, X1 + 1, X2 + 1);
                                s[3] = sc(x1 + 1, X1 + 1, X2 + 1);
                                /* Point 4 */
                                t[4] = tc(x1, X1 + 1, X2);
                                s[4] = sc(x1, X1 + 1, X2);
                                /* Point 5 */
                                s[5] = sc(x1, X1, X2);

                                d01 = h;
                                d12 = rX2*dq;
                                d14 = rX2 * sin(qX1p1) * df;
                                d25 = rX2 * sin(qX1) * df;
                                d45 = rX2*dq;

                                try = fdsphne(t, s, d01, d12, d14, d25, d45);
                                if (try < guess) guess = try;
                            }
                            /* Southern Edge */
                            if (X1 > 0 && X2 > z1 + 1 && X2 < z2 - 1 && time0[index - nx] < 1.e9) {
                                /* Point 0 */
                                t[0] = tc(x1 + 1, X1 - 1, X2 + 1);
                                s[0] = sc(x1 + 1, X1 - 1, X2 + 1);
                                /* Point 1 */
                                t[1] = tc(x1 + 1, X1 - 1, X2);
                                s[1] = sc(x1 + 1, X1 - 1, X2);
                                /* Point 2 */
                                t[2] = tc(x1 + 1, X1, X2);
                                s[2] = sc(x1 + 1, X1, X2);
                                /* Point 3 */
                                t[3] = tc(x1 + 1, X1 - 1, X2 - 1);
                                s[3] = sc(x1 + 1, X1 - 1, X2 - 1);
                                /* Point 4 */
                                t[4] = tc(x1, X1 - 1, X2);
                                s[4] = sc(x1, X1 - 1, X2);
                                /* Point 5 */
                                s[5] = sc(x1, X1, X2);

                                d01 = h;
                                d12 = rX2*dq;
                                d14 = rX2 * sin(qX1m1) * df;
                                d25 = rX2 * sin(qX1) * df;
                                d45 = rX2*dq;

                                try = fdsphne(t, s, d01, d12, d14, d25, d45);
                                if (try < guess) guess = try;
                            }
                        }


                        /*** Two-D Trials ***/
                        /* West-North */
                        if (X1 > 0 && time0[index - nx] < 1.e9) {
                            s[0] = sc(x1, X1, X2);
                            /* Point 1 */
                            t[1] = tc(x1, X1 - 1, X2);
                            s[1] = sc(x1, X1 - 1, X2);
                            /* Point 2 */
                            t[2] = tc(x1 + 1, X1, X2);
                            s[2] = sc(x1 + 1, X1, X2);
                            /* Point 3 */
                            t[3] = tc(x1 + 1, X1 - 1, X2);
                            s[3] = sc(x1 + 1, X1 - 1, X2);

                            d01 = rX2*dq;
                            d02 = rX2 * sinqX1*df;
                            d13 = rX2 * sin(qX1m1) * df;
                            d23 = rX2*dq;

                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) guess = try;
                        }
                        /* West-South */
                        if (X1 < ny - 1 && time0[index + nx] < 1.e9) {
                            s[0] = sc(x1, X1, X2);
                            /* Point 1 */
                            t[1] = tc(x1, X1 + 1, X2);
                            s[1] = sc(x1, X1 + 1, X2);
                            /* Point 2 */
                            t[2] = tc(x1 + 1, X1, X2);
                            s[2] = sc(x1 + 1, X1, X2);
                            /* Point 3 */
                            t[3] = tc(x1 + 1, X1 + 1, X2);
                            s[3] = sc(x1 + 1, X1 + 1, X2);

                            d01 = rX2*dq;
                            d02 = rX2 * sinqX1*df;
                            d13 = rX2 * sin(qX1p1) * df;
                            d23 = rX2*dq;

                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) guess = try;
                        }
                        /* West-Top */
                        if (X2 > 0 && time0[index - nxy] < 1.e9) {
                            s[0] = sc(x1, X1, X2);
                            /* Point 1 */
                            t[1] = tc(x1, X1, X2 - 1);
                            s[1] = sc(x1, X1, X2 - 1);
                            /* Point 2 */
                            t[2] = tc(x1 + 1, X1, X2);
                            s[2] = sc(x1 + 1, X1, X2);
                            /* Point 3 */
                            t[3] = tc(x1 + 1, X1, X2 - 1);
                            s[3] = sc(x1 + 1, X1, X2 - 1);

                            d01 = h;
                            d02 = rX2 * sinqX1*df;
                            d13 = rXm1 * sinqX1*df;
                            d23 = h;

                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) guess = try;
                        }
                        /* West-Bottom */
                        /*if(time0[index+nxy] < 1.e9 && X2<nz-1 )  { */
                        if (X2 < nz - 1) {
                            if (time0[index + nxy] < 1.e9) {
                                s[0] = sc(x1, X1, X2);
                                /* Point 1 */
                                t[1] = tc(x1, X1, X2 + 1);
                                s[1] = sc(x1, X1, X2 + 1);
                                /* Point 2 */
                                t[2] = tc(x1 + 1, X1, X2);
                                s[2] = sc(x1 + 1, X1, X2);
                                /* Point 3 */
                                t[3] = tc(x1 + 1, X1, X2 + 1);
                                s[3] = sc(x1 + 1, X1, X2 + 1);

                                d01 = h;
                                d02 = rX2 * sinqX1*df;
                                d13 = rXp1 * sinqX1*df;
                                d23 = h;

                                try = fdsph2d(t, s, d01, d02, d13, d23);
                                if (try < guess) guess = try;
                            }
                        }
                        /* Parallel - North Bottom */
                        /*if(time0[index-nx] < 1.e9 && time0[index+nxy-nx] < 1.e9
                           && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) { */
                        if (X2 < nz - 1 && X1 > 0) {
                            if (time0[index - nx] < 1.e9 && time0[index + nxy - nx] < 1.e9
                                    && time0[index + nxy] < 1.e9) {
                                s[0] = sc(x1, X1, X2);
                                /* Point 1 */
                                t[1] = tc(x1, X1, X2 + 1);
                                s[1] = sc(x1, X1, X2 + 1);
                                /* Point 2 */
                                t[2] = tc(x1, X1 - 1, X2);
                                s[2] = sc(x1, X1 - 1, X2);
                                /* Point 3 */
                                t[3] = tc(x1, X1 - 1, X2 + 1);
                                s[3] = sc(x1, X1 - 1, X2 + 1);
                                d01 = h;
                                d02 = rX2*dq;
                                d13 = rXp1*dq;
                                d23 = h;
                                try = fdsph2d(t, s, d01, d02, d13, d23);
                                if (try < guess) {
                                    fhead = (guess - try) / (d01 * slow0[index]);
                                    guess = try;
                                }
                            }
                        }
                        /* Parallel - North Top */
                        if (X2 > 0 && X1 > 0 && time0[index - nx] < 1.e9 && time0[index - nxy - nx] < 1.e9
                                && time0[index - nxy] < 1.e9) {
                            s[0] = sc(x1, X1, X2);
                            /* Point 1 */
                            t[1] = tc(x1, X1, X2 - 1);
                            s[1] = sc(x1, X1, X2 - 1);
                            /* Point 2 */
                            t[2] = tc(x1, X1 - 1, X2);
                            s[2] = sc(x1, X1 - 1, X2);
                            /* Point 3 */
                            t[3] = tc(x1, X1 - 1, X2 - 1);
                            s[3] = sc(x1, X1 - 1, X2 - 1);
                            d01 = h;
                            d02 = rX2*dq;
                            d13 = rXm1*dq;
                            d23 = h;
                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) {
                                fhead = (guess - try) / (d01 * slow0[index]);
                                guess = try;
                            }
                        }
                        /* Parallel - South Bottom */
                        /*if(time0[index+nx] < 1.e9 && time0[index+nxy+nx] < 1.e9
                            && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<ny-1 ) { */
                        if (X2 < nz - 1 && X1 < ny - 1) {
                            if (time0[index + nx] < 1.e9 && time0[index + nxy + nx] < 1.e9
                                    && time0[index + nxy] < 1.e9) {
                                s[0] = sc(x1, X1, X2);
                                /* Point 1 */
                                t[1] = tc(x1, X1, X2 + 1);
                                s[1] = sc(x1, X1, X2 + 1);
                                /* Point 2 */
                                t[2] = tc(x1, X1 + 1, X2);
                                s[2] = sc(x1, X1 + 1, X2);
                                /* Point 3 */
                                t[3] = tc(x1, X1 + 1, X2 + 1);
                                s[3] = sc(x1, X1 + 1, X2 + 1);
                                d01 = h;
                                d02 = rX2*dq;
                                d13 = rXp1*dq;
                                d23 = h;
                                try = fdsph2d(t, s, d01, d02, d13, d23);
                                if (try < guess) {
                                    fhead = (guess - try) / (d01 * slow0[index]);
                                    guess = try;
                                }
                            }
                        }
                        /* Parallel - South Top */
                        if (X2 > 0 && X1 < ny - 1 && time0[index + nx] < 1.e9 && time0[index - nxy + nx] < 1.e9
                                && time0[index - nxy] < 1.e9) {
                            s[0] = sc(x1, X1, X2);
                            /* Point 1 */
                            t[1] = tc(x1, X1, X2 - 1);
                            s[1] = sc(x1, X1, X2 - 1);
                            /* Point 2 */
                            t[2] = tc(x1, X1 + 1, X2);
                            s[2] = sc(x1, X1 + 1, X2);
                            /* Point 3 */
                            t[3] = tc(x1, X1 + 1, X2 - 1);
                            s[3] = sc(x1, X1 + 1, X2 - 1);
                            d01 = h;
                            d02 = rX2*dq;
                            d13 = rXm1*dq;
                            d23 = h;
                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) {
                                fhead = (guess - try) / (d01 * slow0[index]);
                                guess = try;
                            }
                        }
                        /*** New Face Trials ***/
                        if (guess > 1.0e9) {
                            if (X1 > y1 + 1 && X1 < y2 - 1 && X2 > z1 + 1 && X2 < z2 - 1) {
                                /* Point 0 */
                                t[0] = tc(x1 + 1, X1 + 1, X2);
                                s[0] = sc(x1 + 1, X1 + 1, X2);
                                /* Point 1 */
                                t[1] = tc(x1 + 1, X1, X2 + 1);
                                s[1] = sc(x1 + 1, X1, X2 + 1);
                                /* Point 2 */
                                t[2] = tc(x1 + 1, X1, X2);
                                s[2] = sc(x1 + 1, X1, X2);
                                /* Point 3 */
                                t[3] = tc(x1 + 1, X1, X2 - 1);
                                s[3] = sc(x1 + 1, X1, X2 - 1);
                                /* Point 4 */
                                t[4] = tc(x1 + 1, X1 - 1, X2);
                                s[4] = sc(x1 + 1, X1 - 1, X2);
                                /* Point 5 */
                                s[5] = sc(x1, X1, X2);
                                d02 = rX2*dq;
                                d12 = h;
                                d25 = rX2 * sinqX1*df;
                                try = fdsphnf(t, s, d02, d12, d25);
                                if (try < guess) guess = try;
                            }
                        }

                        /*** Edge Trials ***/
                        /* Change in Longitude */
                        try = tc(x1 + 1, X1, X2) + .5 * (sc(x1, X1, X2) + sc(x1 + 1, X1, X2)) * rX2 * sinqX1*df;
                        if (try < guess) guess = try;
                        /* Change in South Latitude */
                        if (X1 < ny - 1 && time0[index + nx] < 1.e9) {
                            try = tc(x1, X1 + 1, X2) + .5 * (sc(x1, X1, X2) + sc(x1, X1 + 1, X2)) * rX2*dq;
                            if (try < guess) {
                                fhead = (guess - try) / (rX2 * dq * slow0[index]);
                                guess = try;
                            }
                        }
                        /* Change in North Latitude */
                        if (X1 > 0 && time0[index - nx] < 1.e9) {
                            try = tc(x1, X1 - 1, X2) + .5 * (sc(x1, X1, X2) + sc(x1, X1 - 1, X2)) * rX2*dq;
                            if (try < guess) {
                                fhead = (guess - try) / (rX2 * dq * slow0[index]);
                                guess = try;
                            }
                        }
                        /* Change in depth (increasing) */
                        /*if ( time0[index+nxy]<1.e9 && X2<nz-1 )  { */
                        if (X2 < nz - 1) {
                            if (time0[index + nxy] < 1.e9) {
                                try = tc(x1, X1, X2 + 1) + .5 * (sc(x1, X1, X2) + sc(x1, X1, X2 + 1)) * h;
                                if (try < guess) {
                                    fhead = (guess - try) / (h * slow0[index]);
                                    guess = try;
                                }
                            }
                        }
                        /* Change in depth (decreasing) */
                        if (X2 > 0 && time0[index - nxy] < 1.e9) {
                            try = tc(x1, X1, X2 - 1) + .5 * (sc(x1, X1, X2) + sc(x1, X1, X2 - 1)) * h;
                            if (try < guess) {
                                fhead = (guess - try) / (h * slow0[index]);
                                guess = try;
                            }
                        }
                        if (guess < time0[index]) {
                            time0[index] = guess;
                            if (fhead > headtest) headw[1]++;
                        }
                    }
                    if (x1 == 0) dx1 = 0;
                    x1--;
                }
                /* End x1 if condition */
            }
            /* End growth loop */

            /* NORTH (FRONT) SIDE */
            for (igrow = 1; igrow <= jminus; igrow++) {
                if (dy1) {
                    ii = 0;
                    for (k = z1 + 1; k <= z2 - 1; k++) {
                        for (i = x1 + 1; i <= x2 - 1; i++) {
                            sort[ii].time = tc(i, y1 + 1, k);
                            sort[ii].i1 = i;
                            sort[ii].i2 = k;
                            ii++;
                        }
                    }
                    qsort((char *) sort, ii, sizeof (struct sorted), compar);
                    for (i = 0; i < ii; i++) {
                        X1 = sort[i].i1;
                        X2 = sort[i].i2;
                        index = X2 * nxy + y1 * nx + X1;
                        fhead = 0.;
                        guess = time0[index];

                        rX2 = rc(X2);
                        rXp1 = rc(X2 + 1);
                        rXm1 = rc(X2 - 1);
                        qy1 = qc(y1);
                        qy1p1 = qc(y1 + 1);
                        fx1 = fc(X1);
                        fx1p1 = fc(X1 + 1);
                        fx1m1 = fc(X1 - 1);

                        sinqy1 = sin(qy1);

                        /*** 3D Trials ***/
                        /* Top West corner (point 6) */
                        /*if(time0[index+1] < 1.e9 && time0[index+nxy+1] < 1.e9
                           && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<nx-1 ) { */
                        if (X2 < nz - 1 && X1 < nx - 1) {
                            if (time0[index + 1] < 1.e9 && time0[index + nxy + 1] < 1.e9
                                    && time0[index + nxy] < 1.e9) {
                                /* Point 0 */
                                t[0] = tc(X1 + 1, y1 + 1, X2 + 1);
                                s[0] = sc(X1 + 1, y1 + 1, X2 + 1);
                                r[0] = rXp1;
                                q[0] = qy1p1;
                                f[0] = fx1p1;
                                g[0] = -1;
                                n[0] = 1;
                                m[0] = 1;
                                /* Point 1 */
                                t[1] = tc(X1, y1 + 1, X2 + 1);
                                s[1] = sc(X1, y1 + 1, X2 + 1);
                                r[1] = rXp1;
                                q[1] = qy1p1;
                                f[1] = fx1;
                                g[1] = -1;
                                n[1] = 1;
                                m[1] = -1;
                                /* Point 2 */
                                t[2] = tc(X1, y1, X2 + 1);
                                s[2] = sc(X1, y1, X2 + 1);
                                r[2] = rXp1;
                                q[2] = qy1;
                                f[2] = fx1;
                                g[2] = -1;
                                n[2] = -1;
                                m[2] = -1;
                                /* Point 3 */
                                t[3] = tc(X1 + 1, y1, X2 + 1);
                                s[3] = sc(X1 + 1, y1, X2 + 1);
                                r[3] = rXp1;
                                q[3] = qy1;
                                f[3] = fx1p1;
                                g[3] = -1;
                                n[3] = -1;
                                m[3] = 1;
                                /* Point 4 */
                                t[4] = tc(X1 + 1, y1 + 1, X2);
                                s[4] = sc(X1 + 1, y1 + 1, X2);
                                r[4] = rX2;
                                q[4] = qy1p1;
                                f[4] = fx1p1;
                                g[4] = 1;
                                n[4] = 1;
                                m[4] = 1;
                                /* Point 5 */
                                t[5] = tc(X1, y1 + 1, X2);
                                s[5] = sc(X1, y1 + 1, X2);
                                r[5] = rX2;
                                q[5] = qy1p1;
                                f[5] = fx1;
                                g[5] = 1;
                                n[5] = 1;
                                m[5] = -1;
                                /* Point 7 */
                                t[6] = tc(X1 + 1, y1, X2);
                                s[6] = sc(X1 + 1, y1, X2);
                                r[6] = rX2;
                                q[6] = qy1;
                                f[6] = fx1p1;
                                g[6] = 1;
                                n[6] = -1;
                                m[6] = 1;
                                /* Point 6 */
                                s[7] = sc(X1, y1, X2);
                                r[7] = rX2;
                                q[7] = qy1;
                                f[7] = fx1;
                                g[7] = 1;
                                n[7] = -1;
                                m[7] = -1;

                                try = fdsph3d(t, s, r, q, f, g, n, m, h, dq, df);
                                if (try < guess) guess = try;
                            }
                        }

                        /* Top East corner (point 7) */
                        /*if(time0[index-1] < 1.e9 && time0[index+nxy-1] < 1.e9
                           && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) { */
                        if (X2 < nz - 1 && X1 > 0) {
                            if (time0[index - 1] < 1.e9 && time0[index + nxy - 1] < 1.e9
                                    && time0[index + nxy] < 1.e9) {
                                /* Point 0 */
                                t[0] = tc(X1, y1 + 1, X2 + 1);
                                s[0] = sc(X1, y1 + 1, X2 + 1);
                                r[0] = rXp1;
                                q[0] = qy1p1;
                                f[0] = fx1;
                                g[0] = -1;
                                n[0] = 1;
                                m[0] = 1;
                                /* Point 1 */
                                t[1] = tc(X1 - 1, y1 + 1, X2 + 1);
                                s[1] = sc(X1 - 1, y1 + 1, X2 + 1);
                                r[1] = rXp1;
                                q[1] = qy1p1;
                                f[1] = fx1m1;
                                g[1] = -1;
                                n[1] = 1;
                                m[1] = -1;
                                /* Point 2 */
                                t[2] = tc(X1 - 1, y1, X2 + 1);
                                s[2] = sc(X1 - 1, y1, X2 + 1);
                                r[2] = rXp1;
                                q[2] = qy1;
                                f[2] = fx1m1;
                                g[2] = -1;
                                n[2] = -1;
                                m[2] = -1;
                                /* Point 3 */
                                t[3] = tc(X1, y1, X2 + 1);
                                s[3] = sc(X1, y1, X2 + 1);
                                r[3] = rXp1;
                                q[3] = qy1;
                                f[3] = fx1;
                                g[3] = -1;
                                n[3] = -1;
                                m[3] = 1;
                                /* Point 4 */
                                t[4] = tc(X1, y1 + 1, X2);
                                s[4] = sc(X1, y1 + 1, X2);
                                r[4] = rX2;
                                q[4] = qy1p1;
                                f[4] = fx1;
                                g[4] = 1;
                                n[4] = 1;
                                m[4] = 1;
                                /* Point 5 */
                                t[5] = tc(X1 - 1, y1 + 1, X2);
                                s[5] = sc(X1 - 1, y1 + 1, X2);
                                r[5] = rX2;
                                q[5] = qy1p1;
                                f[5] = fx1m1;
                                g[5] = 1;
                                n[5] = 1;
                                m[5] = -1;
                                /* Point 6 */
                                t[6] = tc(X1 - 1, y1, X2);
                                s[6] = sc(X1 - 1, y1, X2);
                                r[6] = rX2;
                                q[6] = qy1;
                                f[6] = fx1m1;
                                g[6] = 1;
                                n[6] = -1;
                                m[6] = -1;
                                /* Point 7 */
                                s[7] = sc(X1, y1, X2);
                                r[7] = rX2;
                                q[7] = qy1;
                                f[7] = fx1;
                                g[7] = 1;
                                n[7] = -1;
                                m[7] = 1;

                                try = fdsph3d(t, s, r, q, f, g, n, m, h, dq, df);
                                if (try < guess) guess = try;
                            }
                        }

                        /*  Bottom West corner (point 2) */
                        if (X2 > 0 && X1 < nx - 1 && time0[index + 1] < 1.e9 && time0[index - nxy + 1] < 1.e9
                                && time0[index - nxy] < 1.e9) {
                            /* Point 0 */
                            t[0] = tc(X1 + 1, y1 + 1, X2);
                            s[0] = sc(X1 + 1, y1 + 1, X2);
                            r[0] = rX2;
                            q[0] = qy1p1;
                            f[0] = fx1p1;
                            g[0] = -1;
                            n[0] = 1;
                            m[0] = 1;
                            /* Point 1 */
                            t[1] = tc(X1, y1 + 1, X2);
                            s[1] = sc(X1, y1 + 1, X2);
                            r[1] = rX2;
                            q[1] = qy1p1;
                            f[1] = fx1;
                            g[1] = -1;
                            n[1] = 1;
                            m[1] = -1;
                            /* Point 3 */
                            t[2] = tc(X1 + 1, y1, X2);
                            s[2] = sc(X1 + 1, y1, X2);
                            r[2] = rX2;
                            q[2] = qy1;
                            f[2] = fx1p1;
                            g[2] = -1;
                            n[2] = -1;
                            m[2] = 1;
                            /* Point 4 */
                            t[3] = tc(X1 + 1, y1 + 1, X2 - 1);
                            s[3] = sc(X1 + 1, y1 + 1, X2 - 1);
                            r[3] = rXm1;
                            q[3] = qy1p1;
                            f[3] = fx1p1;
                            g[3] = 1;
                            n[3] = 1;
                            m[3] = 1;
                            /* Point 5 */
                            t[4] = tc(X1, y1 + 1, X2 - 1);
                            s[4] = sc(X1, y1 + 1, X2 - 1);
                            r[4] = rXm1;
                            q[4] = qy1p1;
                            f[4] = fx1;
                            g[4] = 1;
                            n[4] = 1;
                            m[4] = -1;
                            /* Point 6 */
                            t[5] = tc(X1, y1, X2 - 1);
                            s[5] = sc(X1, y1, X2 - 1);
                            r[5] = rXm1;
                            q[5] = qy1;
                            f[5] = fx1;
                            g[5] = 1;
                            n[5] = -1;
                            m[5] = -1;
                            /* Point 7 */
                            t[6] = tc(X1 + 1, y1, X2 - 1);
                            s[6] = sc(X1 + 1, y1, X2 - 1);
                            r[6] = rXm1;
                            q[6] = qy1;
                            f[6] = fx1p1;
                            g[6] = 1;
                            n[6] = -1;
                            m[6] = 1;
                            /* Point 2 */
                            s[7] = sc(X1, y1, X2);
                            r[7] = rX2;
                            q[7] = qy1;
                            f[7] = fx1;
                            g[7] = -1;
                            n[7] = -1;
                            m[7] = -1;

                            try = fdsph3d(t, s, r, q, f, g, n, m, h, dq, df);
                            if (try < guess) guess = try;
                        }

                        /* Bottom East corner (point 3) */
                        if (X2 > 0 && X1 > 0 && time0[index - 1] < 1.e9 && time0[index - nxy - 1] < 1.e9
                                && time0[index - nxy] < 1.e9) {
                            /* Point 0 */
                            t[0] = tc(X1, y1 + 1, X2);
                            s[0] = sc(X1, y1 + 1, X2);
                            r[0] = rX2;
                            q[0] = qy1p1;
                            f[0] = fx1;
                            g[0] = -1;
                            n[0] = 1;
                            m[0] = 1;
                            /* Point 1 */
                            t[1] = tc(X1 - 1, y1 + 1, X2);
                            s[1] = sc(X1 - 1, y1 + 1, X2);
                            r[1] = rX2;
                            q[1] = qy1p1;
                            f[1] = fx1m1;
                            g[1] = -1;
                            n[1] = 1;
                            m[1] = -1;
                            /* Point 2 */
                            t[2] = tc(X1 - 1, y1, X2);
                            s[2] = sc(X1 - 1, y1, X2);
                            r[2] = rX2;
                            q[2] = qy1;
                            f[2] = fx1m1;
                            g[2] = -1;
                            n[2] = -1;
                            m[2] = -1;
                            /* Point 4 */
                            t[3] = tc(X1, y1 + 1, X2 - 1);
                            s[3] = sc(X1, y1 + 1, X2 - 1);
                            r[3] = rXm1;
                            q[3] = qy1p1;
                            f[3] = fx1;
                            g[3] = 1;
                            n[3] = 1;
                            m[3] = 1;
                            /* Point 5 */
                            t[4] = tc(X1 - 1, y1 + 1, X2 - 1);
                            s[4] = sc(X1 - 1, y1 + 1, X2 - 1);
                            r[4] = rXm1;
                            q[4] = qy1p1;
                            f[4] = fx1m1;
                            g[4] = 1;
                            n[4] = 1;
                            m[4] = -1;
                            /* Point 6 */
                            t[5] = tc(X1 - 1, y1, X2 - 1);
                            s[5] = sc(X1 - 1, y1, X2 - 1);
                            r[5] = rXm1;
                            q[5] = qy1;
                            f[5] = fx1m1;
                            g[5] = 1;
                            n[5] = -1;
                            m[5] = -1;
                            /* Point 7 */
                            t[6] = tc(X1, y1, X2 - 1);
                            s[6] = sc(X1, y1, X2 - 1);
                            r[6] = rXm1;
                            q[6] = qy1;
                            f[6] = fx1;
                            g[6] = 1;
                            n[6] = -1;
                            m[6] = 1;
                            /* Point 3 */
                            s[7] = sc(X1, y1, X2);
                            r[7] = rX2;
                            q[7] = qy1;
                            f[7] = fx1;
                            g[7] = -1;
                            n[7] = -1;
                            m[7] = 1;
                            try = fdsph3d(t, s, r, q, f, g, n, m, h, dq, df);
                            if (try < guess) guess = try;
                        }

                        /*** New Edge Trials ***/
                        if (guess > 1.0e9) {
                            /* Top edge */
                            /*if(time0[index+nxy] < 1.e9 && X2<nz-1 && X1>x1+1 && X1<x2-1 )  { */
                            if (X2 < nz - 1 && X1 > x1 + 1 && X1 < x2 - 1) {
                                if (time0[index + nxy] < 1.e9) {
                                    /* Point 0 */
                                    t[0] = tc(X1 - 1, y1 + 1, X2 + 1);
                                    s[0] = sc(X1 - 1, y1 + 1, X2 + 1);
                                    /* Point 1 */
                                    t[1] = tc(X1, y1 + 1, X2 + 1);
                                    s[1] = sc(X1, y1 + 1, X2 + 1);
                                    /* Point 2 */
                                    t[2] = tc(X1, y1 + 1, X2);
                                    s[2] = sc(X1, y1 + 1, X2);
                                    /* Point 3 */
                                    t[3] = tc(X1 + 1, y1 + 1, X2 + 1);
                                    s[3] = sc(X1 + 1, y1 + 1, X2 + 1);
                                    /* Point 4 */
                                    t[4] = tc(X1, y1, X2 + 1);
                                    s[4] = sc(X1, y1, X2 + 1);
                                    /* Point 5 */
                                    s[5] = sc(X1, y1, X2);

                                    d01 = rXp1 * sin(qy1p1) * df;
                                    d12 = h;
                                    d14 = rXp1*dq;
                                    d25 = rX2*dq;
                                    d45 = h;

                                    try = fdsphne(t, s, d01, d12, d14, d25, d45);
                                    if (try < guess) guess = try;
                                }
                            }
                            /* Bottom Edge */
                            if (X2 > 0 && X1 > x1 + 1 && X1 < x2 - 1 && time0[index - nxy] < 1.e9) {
                                /* Point 0 */
                                t[0] = tc(X1 + 1, y1 + 1, X2 - 1);
                                s[0] = sc(X1 + 1, y1 + 1, X2 - 1);
                                /* Point 1 */
                                t[1] = tc(X1, y1 + 1, X2 - 1);
                                s[1] = sc(X1, y1 + 1, X2 - 1);
                                /* Point 2 */
                                t[2] = tc(X1, y1 + 1, X2);
                                s[2] = sc(X1, y1 + 1, X2);
                                /* Point 3 */
                                t[3] = tc(X1 - 1, y1 + 1, X2 - 1);
                                s[3] = sc(X1 - 1, y1 + 1, X2 - 1);
                                /* Point 4 */
                                t[4] = tc(X1, y1, X2 - 1);
                                s[4] = sc(X1, y1, X2 - 1);
                                /* Point 5 */
                                s[5] = sc(X1, y1, X2);

                                d01 = rXm1 * sin(qy1p1) * df;
                                d12 = h;
                                d14 = rXm1*dq;
                                d25 = rX2*dq;
                                d45 = h;

                                try = fdsphne(t, s, d01, d12, d14, d25, d45);
                                if (try < guess) guess = try;
                            }
                            /* Eastern Edge */
                            if (X1 > 0 && X2 > z1 + 1 && X2 < z2 - 1 & time0[index - 1] < 1.e9) {
                                /* Point 0 */
                                t[0] = tc(X1 - 1, y1 + 1, X2 - 1);
                                s[0] = sc(X1 - 1, y1 + 1, X2 - 1);
                                /* Point 1 */
                                t[1] = tc(X1 - 1, y1 + 1, X2);
                                s[1] = sc(X1 - 1, y1 + 1, X2);
                                /* Point 2 */
                                t[2] = tc(X1, y1 + 1, X2);
                                s[2] = sc(X1, y1 + 1, X2);
                                /* Point 3 */
                                t[3] = tc(X1 - 1, y1 + 1, X2 + 1);
                                s[3] = sc(X1 - 1, y1 + 1, X2 + 1);
                                /* Point 4 */
                                t[4] = tc(X1 - 1, y1, X2);
                                s[4] = sc(X1 - 1, y1, X2);
                                /* Point 5 */
                                s[5] = sc(X1, y1, X2);

                                d01 = h;
                                d12 = rX2 * sin(qy1p1) * df;
                                d14 = rX2*dq;
                                d25 = rX2*dq;
                                d45 = rX2 * sin(qy1) * df;

                                try = fdsphne(t, s, d01, d12, d14, d25, d45);
                                if (try < guess) guess = try;
                            }
                            /* Western Edge */
                            if (X1 < nx - 1 && X2 > z1 + 1 && X2 < z2 - 1 & time0[index + 1] < 1.e9) {
                                /* Point 0 */
                                t[0] = tc(X1 + 1, y1 + 1, X2 + 1);
                                s[0] = sc(X1 + 1, y1 + 1, X2 + 1);
                                /* Point 1 */
                                t[1] = tc(X1 + 1, y1 + 1, X2);
                                s[1] = sc(X1 + 1, y1 + 1, X2);
                                /* Point 2 */
                                t[2] = tc(X1, y1 + 1, X2);
                                s[2] = sc(X1, y1 + 1, X2);
                                /* Point 3 */
                                t[3] = tc(X1 + 1, y1 + 1, X2 - 1);
                                s[3] = sc(X1 + 1, y1 + 1, X2 - 1);
                                /* Point 4 */
                                t[4] = tc(X1 + 1, y1, X2);
                                s[4] = sc(X1 + 1, y1, X2);
                                /* Point 5 */
                                s[5] = sc(X1, y1, X2);

                                d01 = h;
                                d12 = rX2 * sin(qy1p1) * df;
                                d14 = rX2*dq;
                                d25 = rX2*dq;
                                d45 = rX2 * sin(qy1) * df;

                                try = fdsphne(t, s, d01, d12, d14, d25, d45);
                                if (try < guess) guess = try;
                            }
                        }

                        /*** Two-D Trials ***/
                        /* North-East */
                        if (X1 < nx - 1 && time0[index + 1] < 1.e9) {
                            s[0] = sc(X1, y1, X2);
                            /* Point 1 */
                            t[1] = tc(X1 + 1, y1, X2);
                            s[1] = sc(X1 + 1, y1, X2);
                            /* Point 2 */
                            t[2] = tc(X1, y1 + 1, X2);
                            s[2] = sc(X1, y1 + 1, X2);
                            /* Point 3 */
                            t[3] = tc(X1 + 1, y1 + 1, X2);
                            s[3] = sc(X1 + 1, y1 + 1, X2);

                            d01 = rX2 * sinqy1*df;
                            d02 = rX2*dq;
                            d13 = rX2*dq;
                            d23 = rX2 * sin(qy1p1) * df;

                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) guess = try;
                        }
                        /* North-West */
                        if (X1 > 0 && time0[index - 1] < 1.e9) {
                            s[0] = sc(X1, y1, X2);
                            /* Point 1 */
                            t[1] = tc(X1 - 1, y1, X2);
                            s[1] = sc(X1 - 1, y1, X2);
                            /* Point 2 */
                            t[2] = tc(X1, y1 + 1, X2);
                            s[2] = sc(X1, y1 + 1, X2);
                            /* Point 3 */
                            t[3] = tc(X1 - 1, y1 + 1, X2);
                            s[3] = sc(X1 - 1, y1 + 1, X2);

                            d01 = rX2 * sinqy1*df;
                            d02 = rX2*dq;
                            d13 = rX2*dq;
                            d23 = rX2 * sin(qy1p1) * df;

                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) guess = try;
                        }
                        /* North-Top */
                        if (X2 > 0 && time0[index - nxy] < 1.e9) {
                            s[0] = sc(X1, y1, X2);
                            /* Point 1 */
                            t[1] = tc(X1, y1, X2 - 1);
                            s[1] = sc(X1, y1, X2 - 1);
                            /* Point 2 */
                            t[2] = tc(X1, y1 + 1, X2);
                            s[2] = sc(X1, y1 + 1, X2);
                            /* Point 3 */
                            t[3] = tc(X1, y1 + 1, X2 - 1);
                            s[3] = sc(X1, y1 + 1, X2 - 1);

                            d01 = h;
                            d02 = rX2*dq;
                            d13 = rXm1*dq;
                            d23 = h;

                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) guess = try;
                        }
                        /* North-Bottom */
                        /*if(time0[index+nxy] < 1.e9 && X2<nz-1 )  { */
                        if (X2 < nz - 1) {
                            if (time0[index + nxy] < 1.e9) {
                                s[0] = sc(X1, y1, X2);
                                /* Point 1 */
                                t[1] = tc(X1, y1, X2 + 1);
                                s[1] = sc(X1, y1, X2 + 1);
                                /* Point 2 */
                                t[2] = tc(X1, y1 + 1, X2);
                                s[2] = sc(X1, y1 + 1, X2);
                                /* Point 3 */
                                t[3] = tc(X1, y1 + 1, X2 + 1);
                                s[3] = sc(X1, y1 + 1, X2 + 1);

                                d01 = h;
                                d02 = rX2*dq;
                                d13 = rXp1*dq;
                                d23 = h;

                                try = fdsph2d(t, s, d01, d02, d13, d23);
                                if (try < guess) guess = try;
                            }
                        }
                        /* Parallel - East Bottom */
                        /*if(time0[index+1] < 1.e9 && time0[index+nxy+1] < 1.e9
                           && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<nx-1 ) { */
                        if (X2 < nz - 1 && X1 < nx - 1) {
                            if (time0[index + 1] < 1.e9 && time0[index + nxy + 1] < 1.e9
                                    && time0[index + nxy] < 1.e9) {
                                s[0] = sc(X1, y1, X2);
                                /* Point 1 */
                                t[1] = tc(X1 + 1, y1, X2);
                                s[1] = sc(X1 + 1, y1, X2);
                                /* Point 2 */
                                t[2] = tc(X1, y1, X2 + 1);
                                s[2] = sc(X1, y1, X2 + 1);
                                /* Point 3 */
                                t[3] = tc(X1 + 1, y1, X2 + 1);
                                s[3] = sc(X1 + 1, y1, X2 + 1);

                                d01 = rX2 * sinqy1*df;
                                d02 = h;
                                d13 = h;
                                d23 = rXp1 * sinqy1*df;

                                try = fdsph2d(t, s, d01, d02, d13, d23);
                                if (try < guess) {
                                    fhead = (guess - try) / (d01 * slow0[index]);
                                    guess = try;
                                }
                            }
                        }
                        /* Parallel - East Top */
                        if (X2 > 0 && X1 < nx - 1 && time0[index + 1] < 1.e9 && time0[index - nxy + 1] < 1.e9
                                && time0[index - nxy] < 1.e9) {
                            s[0] = sc(X1, y1, X2);
                            /* Point 1 */
                            t[1] = tc(X1 + 1, y1, X2);
                            s[1] = sc(X1 + 1, y1, X2);
                            /* Point 2 */
                            t[2] = tc(X1, y1, X2 - 1);
                            s[2] = sc(X1, y1, X2 - 1);
                            /* Point 3 */
                            t[3] = tc(X1 + 1, y1, X2 - 1);
                            s[3] = sc(X1 + 1, y1, X2 - 1);

                            d01 = rX2 * sinqy1*df;
                            d02 = h;
                            d13 = h;
                            d23 = rXm1 * sinqy1*df;

                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) {
                                fhead = (guess - try) / (d01 * slow0[index]);
                                guess = try;
                            }
                        }
                        /* Parallel - West Bottom */
                        /*if(time0[index-1] < 1.e9 && time0[index+nxy-1] < 1.e9
                           && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) { */
                        if (X2 < nz - 1 && X1 > 0) {
                            if (time0[index - 1] < 1.e9 && time0[index + nxy - 1] < 1.e9
                                    && time0[index + nxy] < 1.e9) {
                                s[0] = sc(X1, y1, X2);
                                /* Point 1 */
                                t[1] = tc(X1 - 1, y1, X2);
                                s[1] = sc(X1 - 1, y1, X2);
                                /* Point 2 */
                                t[2] = tc(X1, y1, X2 + 1);
                                s[2] = sc(X1, y1, X2 + 1);
                                /* Point 3 */
                                t[3] = tc(X1 - 1, y1, X2 + 1);
                                s[3] = sc(X1 - 1, y1, X2 + 1);
                                d01 = rX2 * sinqy1*df;
                                d02 = h;
                                d13 = h;
                                d23 = rXp1 * sinqy1*df;
                                try = fdsph2d(t, s, d01, d02, d13, d23);
                                if (try < guess) {
                                    fhead = (guess - try) / (d01 * slow0[index]);
                                    guess = try;
                                }
                            }
                        }
                        /* Parallel - West Top */
                        if (X2 > 0 && X1 > 0 && time0[index - 1] < 1.e9 && time0[index - nxy - 1] < 1.e9
                                && time0[index - nxy] < 1.e9) {
                            s[0] = sc(X1, y1, X2);
                            /* Point 1 */
                            t[1] = tc(X1 - 1, y1, X2);
                            s[1] = sc(X1 - 1, y1, X2);
                            /* Point 2 */
                            t[2] = tc(X1, y1, X2 - 1);
                            s[2] = sc(X1, y1, X2 - 1);
                            /* Point 3 */
                            t[3] = tc(X1 - 1, y1, X2 - 1);
                            s[3] = sc(X1 - 1, y1, X2 - 1);

                            d01 = rX2 * sinqy1*df;
                            d02 = h;
                            d13 = h;
                            d23 = rXm1 * sinqy1*df;

                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) {
                                fhead = (guess - try) / (d01 * slow0[index]);
                                guess = try;
                            }
                        }

                        /*** New Face Trials ***/
                        if (guess > 1.0e9) {
                            if (X1 > x1 + 1 && X1 < x2 - 1 && X2 > z1 + 1 && X2 < z2 - 1) {
                                /* Point 0 */
                                t[0] = tc(X1, y1 + 1, X2 + 1);
                                s[0] = sc(X1, y1 + 1, X2 + 1);
                                /* Point 1 */
                                t[1] = tc(X1 + 1, y1 + 1, X2);
                                s[1] = sc(X1 + 1, y1 + 1, X2);
                                /* Point 2 */
                                t[2] = tc(X1, y1 + 1, X2);
                                s[2] = sc(X1, y1 + 1, X2);
                                /* Point 3 */
                                t[3] = tc(X1 - 1, y1 + 1, X2);
                                s[3] = sc(X1 - 1, y1 + 1, X2);
                                /* Point 4 */
                                t[4] = tc(X1, y1 + 1, X2 - 1);
                                s[4] = sc(X1, y1 + 1, X2 - 1);
                                /* Point 5 */
                                s[5] = sc(X1, y1, X2);
                                d02 = h;
                                d12 = rX2 * sin(qy1p1) * df;
                                d25 = rX2*dq;
                                try = fdsphnf(t, s, d02, d12, d25);
                                if (try < guess) guess = try;
                            }
                        }

                        /*** Edge Trials ***/
                        /* Change in Latitude */
                        try = tc(X1, y1 + 1, X2) + .5 * (sc(X1, y1 + 1, X2) + sc(X1, y1, X2)) * rX2*dq;
                        if (try < guess) guess = try;
                        /* Change in East Longitude */
                        if (X1 < nx - 1 && time0[index + 1] < 1.e9) {
                            try = tc(X1 + 1, y1, X2) + .5 * (sc(X1, y1, X2) + sc(X1 + 1, y1, X2)) * rX2 * df*sinqy1;
                            if (try < guess) {
                                fhead = (guess - try) / (rX2 * df * sinqy1 * slow0[index]);
                                guess = try;
                            }
                        }
                        /* Change in West Longitude */
                        if (X1 > 0 && time0[index - 1] < 1.e9) {
                            try = tc(X1 - 1, y1, X2) + .5 * (sc(X1, y1, X2) + sc(X1 - 1, y1, X2)) * rX2 * df*sinqy1;
                            if (try < guess) {
                                fhead = (guess - try) / (rX2 * df * sinqy1 * slow0[index]);
                                guess = try;
                            }
                        }
                        /* Change in depth (increasing) */
                        /*if ( time0[index+nxy]<1.e9 && X2<nz-1 )  { */
                        if (X2 < nz - 1) {
                            if (time0[index + nxy] < 1.e9) {
                                try = tc(X1, y1, X2 + 1) + .5 * (sc(X1, y1, X2) + sc(X1, y1, X2 + 1)) * h;
                                if (try < guess) {
                                    fhead = (guess - try) / (h * slow0[index]);
                                    guess = try;
                                }
                            }
                        }
                        /* Change in depth (decreasing) */
                        if (X2 > 0 && time0[index - nxy] < 1.e9) {
                            try = tc(X1, y1, X2 - 1) + .5 * (sc(X1, y1, X2) + sc(X1, y1, X2 - 1)) * h;
                            if (try < guess) {
                                fhead = (guess - try) / (h * slow0[index]);
                                guess = try;
                            }
                        }
                        if (guess < time0[index]) {
                            time0[index] = guess;
                            if (fhead > headtest) headw[3]++;
                        }
                    }
                    if (y1 == 0) dy1 = 0;
                    y1--;
                }
                /* End y1 if condition */
            }
            /* End growth loop */


            /* SOUTH (BACK) SIDE */
            for (igrow = 1; igrow <= jplus; igrow++) {
                if (dy2) {
                    ii = 0;
                    for (k = z1 + 1; k <= z2 - 1; k++) {
                        for (i = x1 + 1; i <= x2 - 1; i++) {
                            sort[ii].time = tc(i, y2 - 1, k);
                            sort[ii].i1 = i;
                            sort[ii].i2 = k;
                            ii++;
                        }
                    }
                    qsort((char *) sort, ii, sizeof (struct sorted), compar);
                    for (i = 0; i < ii; i++) {
                        X1 = sort[i].i1;
                        X2 = sort[i].i2;
                        index = X2 * nxy + y2 * nx + X1;
                        fhead = 0.;
                        guess = time0[index];

                        rX2 = rc(X2);
                        rXp1 = rc(X2 + 1);
                        rXm1 = rc(X2 - 1);
                        qy2 = qc(y2);
                        qy2m1 = qc(y2 - 1);
                        fx1 = fc(X1);
                        fx1p1 = fc(X1 + 1);
                        fx1m1 = fc(X1 - 1);

                        sinqy2 = sin(qy2);

                        /*** 3D Trials ***/
                        /* Top West corner (point 5) */
                        /*if(time0[index+1] < 1.e9 && time0[index+nxy+1] < 1.e9
                            && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<nx-1 ) { */
                        if (X2 < nz - 1 && X1 < nx - 1) {
                            if (time0[index + 1] < 1.e9 && time0[index + nxy + 1] < 1.e9
                                    && time0[index + nxy] < 1.e9) {
                                /* Point 0 */
                                t[0] = tc(X1 + 1, y2, X2 + 1);
                                s[0] = sc(X1 + 1, y2, X2 + 1);
                                r[0] = rXp1;
                                q[0] = qy2;
                                f[0] = fx1p1;
                                g[0] = -1;
                                n[0] = 1;
                                m[0] = 1;
                                /* Point 1 */
                                t[1] = tc(X1, y2, X2 + 1);
                                s[1] = sc(X1, y2, X2 + 1);
                                r[1] = rXp1;
                                q[1] = qy2;
                                f[1] = fx1;
                                g[1] = -1;
                                n[1] = 1;
                                m[1] = -1;
                                /* Point 2 */
                                t[2] = tc(X1, y2 - 1, X2 + 1);
                                s[2] = sc(X1, y2 - 1, X2 + 1);
                                r[2] = rXp1;
                                q[2] = qy2m1;
                                f[2] = fx1;
                                g[2] = -1;
                                n[2] = -1;
                                m[2] = -1;
                                /* Point 3 */
                                t[3] = tc(X1 + 1, y2 - 1, X2 + 1);
                                s[3] = sc(X1 + 1, y2 - 1, X2 + 1);
                                r[3] = rXp1;
                                q[3] = qy2m1;
                                f[3] = fx1p1;
                                g[3] = -1;
                                n[3] = -1;
                                m[3] = 1;
                                /* Point 4 */
                                t[4] = tc(X1 + 1, y2, X2);
                                s[4] = sc(X1 + 1, y2, X2);
                                r[4] = rX2;
                                q[4] = qy2;
                                f[4] = fx1p1;
                                g[4] = 1;
                                n[4] = 1;
                                m[4] = 1;
                                /* Point 6 */
                                t[5] = tc(X1, y2 - 1, X2);
                                s[5] = sc(X1, y2 - 1, X2);
                                r[5] = rX2;
                                q[5] = qy2m1;
                                f[5] = fx1;
                                g[5] = 1;
                                n[5] = -1;
                                m[5] = -1;
                                /* Point 7 */
                                t[6] = tc(X1 + 1, y2 - 1, X2);
                                s[6] = sc(X1 + 1, y2 - 1, X2);
                                r[6] = rX2;
                                q[6] = qy2m1;
                                f[6] = fx1p1;
                                g[6] = 1;
                                n[6] = -1;
                                m[6] = 1;
                                /* Point 5 */
                                s[7] = sc(X1, y2, X2);
                                r[7] = rX2;
                                q[7] = qy2;
                                f[7] = fx1;
                                g[7] = 1;
                                n[7] = 1;
                                m[7] = -1;

                                try = fdsph3d(t, s, r, q, f, g, n, m, h, dq, df);
                                if (try < guess) guess = try;
                            }
                        }

                        /* Top East corner (point 4) */
                        /*if(time0[index-1] < 1.e9 && time0[index+nxy-1] < 1.e9
                           && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) { */
                        if (X2 < nz - 1 && X1 > 0) {
                            if (time0[index - 1] < 1.e9 && time0[index + nxy - 1] < 1.e9
                                    && time0[index + nxy] < 1.e9) {
                                /* Point 3 */
                                t[0] = tc(X1, y2 - 1, X2 + 1);
                                s[0] = sc(X1, y2 - 1, X2 + 1);
                                r[0] = rXp1;
                                q[0] = qy2m1;
                                f[0] = fx1;
                                g[0] = -1;
                                n[0] = -1;
                                m[0] = 1;
                                /* Point 2 */
                                t[1] = tc(X1 - 1, y2 - 1, X2 + 1);
                                s[1] = sc(X1 - 1, y2 - 1, X2 + 1);
                                r[1] = rXp1;
                                q[1] = qy2m1;
                                f[1] = fx1m1;
                                g[1] = -1;
                                n[1] = -1;
                                m[1] = -1;
                                /* Point 1 */
                                t[2] = tc(X1 - 1, y2, X2 + 1);
                                s[2] = sc(X1 - 1, y2, X2 + 1);
                                r[2] = rXp1;
                                q[2] = qy2;
                                f[2] = fx1m1;
                                g[2] = -1;
                                n[2] = 1;
                                m[2] = -1;
                                /* Point 0 */
                                t[3] = tc(X1, y2, X2 + 1);
                                s[3] = sc(X1, y2, X2 + 1);
                                r[3] = rXp1;
                                q[3] = qy2;
                                f[3] = fx1;
                                g[3] = -1;
                                n[3] = 1;
                                m[3] = 1;
                                /* Point 7 */
                                t[4] = tc(X1, y2 - 1, X2);
                                s[4] = sc(X1, y2 - 1, X2);
                                r[4] = rX2;
                                q[4] = qy2m1;
                                f[4] = fx1;
                                g[4] = 1;
                                n[4] = -1;
                                m[4] = 1;
                                /* Point 6 */
                                t[5] = tc(X1 - 1, y2 - 1, X2);
                                s[5] = sc(X1 - 1, y2 - 1, X2);
                                r[5] = rX2;
                                q[5] = qy2m1;
                                f[5] = fx1m1;
                                g[5] = 1;
                                n[5] = -1;
                                m[5] = -1;
                                /* Point 5 */
                                t[6] = tc(X1 - 1, y2, X2);
                                s[6] = sc(X1 - 1, y2, X2);
                                r[6] = rX2;
                                q[6] = qy2;
                                f[6] = fx1m1;
                                g[6] = 1;
                                n[6] = 1;
                                m[6] = -1;
                                /* Point 4 */
                                s[7] = sc(X1, y2, X2);
                                r[7] = rX2;
                                q[7] = qy2;
                                f[7] = fx1;
                                g[7] = 1;
                                n[7] = 1;
                                m[7] = 1;

                                try = fdsph3d(t, s, r, q, f, g, n, m, h, dq, df);
                                if (try < guess) guess = try;
                            }
                        }

                        /*  Bottom West corner (point 1) */
                        if (X2 > 0 && X1 < nx - 1 && time0[index + 1] < 1.e9 && time0[index - nxy + 1] < 1.e9
                                && time0[index - nxy] < 1.e9) {
                            /* Point 3 */
                            t[0] = tc(X1 + 1, y2 - 1, X2);
                            s[0] = sc(X1 + 1, y2 - 1, X2);
                            r[0] = rX2;
                            q[0] = qy2m1;
                            f[0] = fx1p1;
                            g[0] = -1;
                            n[0] = -1;
                            m[0] = 1;
                            /* Point 2 */
                            t[1] = tc(X1, y2 - 1, X2);
                            s[1] = sc(X1, y2 - 1, X2);
                            r[1] = rX2;
                            q[1] = qy2m1;
                            f[1] = fx1;
                            g[1] = -1;
                            n[1] = -1;
                            m[1] = -1;
                            /* Point 0 */
                            t[2] = tc(X1 + 1, y2, X2);
                            s[2] = sc(X1 + 1, y2, X2);
                            r[2] = rX2;
                            q[2] = qy2;
                            f[2] = fx1p1;
                            g[2] = -1;
                            n[2] = 1;
                            m[2] = 1;
                            /* Point 7 */
                            t[3] = tc(X1 + 1, y2 - 1, X2 - 1);
                            s[3] = sc(X1 + 1, y2 - 1, X2 - 1);
                            r[3] = rXm1;
                            q[3] = qy2m1;
                            f[3] = fx1p1;
                            g[3] = 1;
                            n[3] = -1;
                            m[3] = 1;
                            /* Point 6 */
                            t[4] = tc(X1, y2 - 1, X2 - 1);
                            s[4] = sc(X1, y2 - 1, X2 - 1);
                            r[4] = rXm1;
                            q[4] = qy2m1;
                            f[4] = fx1;
                            g[4] = 1;
                            n[4] = -1;
                            m[4] = -1;
                            /* Point 5 */
                            t[5] = tc(X1, y2, X2 - 1);
                            s[5] = sc(X1, y2, X2 - 1);
                            r[5] = rXm1;
                            q[5] = qy2;
                            f[5] = fx1;
                            g[5] = 1;
                            n[5] = 1;
                            m[5] = -1;
                            /* Point 4 */
                            t[6] = tc(X1 + 1, y2, X2 - 1);
                            s[6] = sc(X1 + 1, y2, X2 - 1);
                            r[6] = rXm1;
                            q[6] = qy2;
                            f[6] = fx1p1;
                            g[6] = 1;
                            n[6] = 1;
                            m[6] = 1;
                            /* Point 1 */
                            s[7] = sc(X1, y2, X2);
                            r[7] = rX2;
                            q[7] = qy2;
                            f[7] = fx1;
                            g[7] = -1;
                            n[7] = 1;
                            m[7] = -1;

                            try = fdsph3d(t, s, r, q, f, g, n, m, h, dq, df);
                            if (try < guess) guess = try;
                        }

                        /* Bottom East corner (point 0) */
                        if (X2 > 0 && X1 > 0 && time0[index - 1] < 1.e9 && time0[index - nxy - 1] < 1.e9
                                && time0[index - nxy] < 1.e9) {
                            /* Point 3 */
                            t[0] = tc(X1, y2 - 1, X2);
                            s[0] = sc(X1, y2 - 1, X2);
                            r[0] = rX2;
                            q[0] = qy2m1;
                            f[0] = fx1;
                            g[0] = -1;
                            n[0] = -1;
                            m[0] = 1;
                            /* Point 2 */
                            t[1] = tc(X1 - 1, y2 - 1, X2);
                            s[1] = sc(X1 - 1, y2 - 1, X2);
                            r[1] = rX2;
                            q[1] = qy2m1;
                            f[1] = fx1m1;
                            g[1] = -1;
                            n[1] = -1;
                            m[1] = -1;
                            /* Point 1 */
                            t[2] = tc(X1 - 1, y2, X2);
                            s[2] = sc(X1 - 1, y2, X2);
                            r[2] = rX2;
                            q[2] = qy2;
                            f[2] = fx1m1;
                            g[2] = -1;
                            n[2] = 1;
                            m[2] = -1;
                            /* Point 7 */
                            t[3] = tc(X1, y2 - 1, X2 - 1);
                            s[3] = sc(X1, y2 - 1, X2 - 1);
                            r[3] = rXm1;
                            q[3] = qy2m1;
                            f[3] = fx1;
                            g[3] = 1;
                            n[3] = -1;
                            m[3] = 1;
                            /* Point 6 */
                            t[4] = tc(X1 - 1, y2 - 1, X2 - 1);
                            s[4] = sc(X1 - 1, y2 - 1, X2 - 1);
                            r[4] = rXm1;
                            q[4] = qy2m1;
                            f[4] = fx1m1;
                            g[4] = 1;
                            n[4] = -1;
                            m[4] = -1;
                            /* Point 5 */
                            t[5] = tc(X1 - 1, y2, X2 - 1);
                            s[5] = sc(X1 - 1, y2, X2 - 1);
                            r[5] = rXm1;
                            q[5] = qy2;
                            f[5] = fx1m1;
                            g[5] = 1;
                            n[5] = 1;
                            m[5] = -1;
                            /* Point 4 */
                            t[6] = tc(X1, y2, X2 - 1);
                            s[6] = sc(X1, y2, X2 - 1);
                            r[6] = rXm1;
                            q[6] = qy2;
                            f[6] = fx1;
                            g[6] = 1;
                            n[6] = 1;
                            m[6] = 1;
                            /* Point 0 */
                            s[7] = sc(X1, y2, X2);
                            r[7] = rX2;
                            q[7] = qy2;
                            f[7] = fx1;
                            g[7] = -1;
                            n[7] = 1;
                            m[7] = 1;
                            try = fdsph3d(t, s, r, q, f, g, n, m, h, dq, df);
                            if (try < guess) guess = try;
                        }

                        /*** New Edge Trials ***/
                        if (guess > 1.0e9) {
                            /* Top edge */
                            /*if(time0[index+nxy] < 1.e9 && X2<nz-1 && X1>x1+1 && X1<x2-1 )  { */
                            if (X2 < nz - 1 && X1 > x1 + 1 && X1 < x2 - 1) {
                                if (time0[index + nxy] < 1.e9) {
                                    /* Point 0 */
                                    t[0] = tc(X1 + 1, y2 - 1, X2 + 1);
                                    s[0] = sc(X1 + 1, y2 - 1, X2 + 1);
                                    /* Point 1 */
                                    t[1] = tc(X1, y2 - 1, X2 + 1);
                                    s[1] = sc(X1, y2 - 1, X2 + 1);
                                    /* Point 2 */
                                    t[2] = tc(X1, y2 - 1, X2);
                                    s[2] = sc(X1, y2 - 1, X2);
                                    /* Point 3 */
                                    t[3] = tc(X1 - 1, y2 - 1, X2 + 1);
                                    s[3] = sc(X1 - 1, y2 - 1, X2 + 1);
                                    /* Point 4 */
                                    t[4] = tc(X1, y2, X2 + 1);
                                    s[4] = sc(X1, y2, X2 + 1);
                                    /* Point 5 */
                                    s[5] = sc(X1, y2, X2);

                                    d01 = rXp1 * sin(qy2m1) * df;
                                    d12 = h;
                                    d14 = rXp1*dq;
                                    d25 = rX2*dq;
                                    d45 = h;

                                    try = fdsphne(t, s, d01, d12, d14, d25, d45);
                                    if (try < guess) guess = try;
                                }
                            }
                            /* Bottom Edge */
                            if (X2 > 0 && X1 > x1 + 1 && X1 < x2 - 1 && time0[index - nxy] < 1.e9) {
                                /* Point 0 */
                                t[0] = tc(X1 - 1, y2 - 1, X2 - 1);
                                s[0] = sc(X1 - 1, y2 - 1, X2 - 1);
                                /* Point 1 */
                                t[1] = tc(X1, y2 - 1, X2 - 1);
                                s[1] = sc(X1, y2 - 1, X2 - 1);
                                /* Point 2 */
                                t[2] = tc(X1, y2 - 1, X2);
                                s[2] = sc(X1, y2 - 1, X2);
                                /* Point 3 */
                                t[3] = tc(X1 + 1, y2 - 1, X2 - 1);
                                s[3] = sc(X1 + 1, y2 - 1, X2 - 1);
                                /* Point 4 */
                                t[4] = tc(X1, y2, X2 - 1);
                                s[4] = sc(X1, y2, X2 - 1);
                                /* Point 5 */
                                s[5] = sc(X1, y2, X2);

                                d01 = rXm1 * sin(qy2m1) * df;
                                d12 = h;
                                d14 = rXm1*dq;
                                d25 = rX2*dq;
                                d45 = h;

                                try = fdsphne(t, s, d01, d12, d14, d25, d45);
                                if (try < guess) guess = try;
                            }
                            /* Eastern Edge */
                            if (X1 > 0 && X2 > z1 + 1 && X2 < z2 - 1 && time0[index - 1] < 1.e9) {
                                /* Point 0 */
                                t[0] = tc(X1 - 1, y2 - 1, X2 + 1);
                                s[0] = sc(X1 - 1, y2 - 1, X2 + 1);
                                /* Point 1 */
                                t[1] = tc(X1 - 1, y2 - 1, X2);
                                s[1] = sc(X1 - 1, y2 - 1, X2);
                                /* Point 2 */
                                t[2] = tc(X1, y2 - 1, X2);
                                s[2] = sc(X1, y2 - 1, X2);
                                /* Point 3 */
                                t[3] = tc(X1 - 1, y2 - 1, X2 - 1);
                                s[3] = sc(X1 - 1, y2 - 1, X2 - 1);
                                /* Point 4 */
                                t[4] = tc(X1 - 1, y2, X2);
                                s[4] = sc(X1 - 1, y2, X2);
                                /* Point 5 */
                                s[5] = sc(X1, y2, X2);

                                d01 = h;
                                d12 = rX2 * sin(qy2m1) * df;
                                d14 = rX2*dq;
                                d25 = rX2*dq;
                                d45 = rX2 * sin(qy2) * df;

                                try = fdsphne(t, s, d01, d12, d14, d25, d45);
                                if (try < guess) guess = try;
                            }
                            /* Western Edge */
                            if (X1 < nx - 1 && X2 > z1 + 1 && X2 < z2 - 1 && time0[index + 1] < 1.e9) {
                                /* Point 0 */
                                t[0] = tc(X1 + 1, y2 - 1, X2 - 1);
                                s[0] = sc(X1 + 1, y2 - 1, X2 - 1);
                                /* Point 1 */
                                t[1] = tc(X1 + 1, y2 - 1, X2);
                                s[1] = sc(X1 + 1, y2 - 1, X2);
                                /* Point 2 */
                                t[2] = tc(X1, y2 - 1, X2);
                                s[2] = sc(X1, y2 - 1, X2);
                                /* Point 3 */
                                t[3] = tc(X1 + 1, y2 - 1, X2 + 1);
                                s[3] = sc(X1 + 1, y2 - 1, X2 + 1);
                                /* Point 4 */
                                t[4] = tc(X1 + 1, y2, X2);
                                s[4] = sc(X1 + 1, y2, X2);
                                /* Point 5 */
                                s[5] = sc(X1, y2, X2);

                                d01 = h;
                                d12 = rX2 * sin(qy2m1) * df;
                                d14 = rX2*dq;
                                d25 = rX2*dq;
                                d45 = rX2 * sin(qy2) * df;

                                try = fdsphne(t, s, d01, d12, d14, d25, d45);
                                if (try < guess) guess = try;
                            }
                        }

                        /*** Two-D Trials ***/
                        /* South-East */
                        if (X1 < nx - 1 && time0[index + 1] < 1.e9) {
                            s[0] = sc(X1, y2, X2);
                            /* Point 1 */
                            t[1] = tc(X1 + 1, y2, X2);
                            s[1] = sc(X1 + 1, y2, X2);
                            /* Point 2 */
                            t[2] = tc(X1, y2 - 1, X2);
                            s[2] = sc(X1, y2 - 1, X2);
                            /* Point 3 */
                            t[3] = tc(X1 + 1, y2 - 1, X2);
                            s[3] = sc(X1 + 1, y2 - 1, X2);

                            d01 = rX2 * sinqy2*df;
                            d02 = rX2*dq;
                            d13 = rX2*dq;
                            d23 = rX2 * sin(qy2m1) * df;

                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) guess = try;
                        }
                        /* South-West */
                        if (X1 > 0 && time0[index - 1] < 1.e9) {
                            s[0] = sc(X1, y2, X2);
                            /* Point 1 */
                            t[1] = tc(X1 - 1, y2, X2);
                            s[1] = sc(X1 - 1, y2, X2);
                            /* Point 2 */
                            t[2] = tc(X1, y2 - 1, X2);
                            s[2] = sc(X1, y2 - 1, X2);
                            /* Point 3 */
                            t[3] = tc(X1 - 1, y2 - 1, X2);
                            s[3] = sc(X1 - 1, y2 - 1, X2);

                            d01 = rX2 * sinqy2*df;
                            d02 = rX2*dq;
                            d13 = rX2*dq;
                            d23 = rX2 * sin(qy2m1) * df;

                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) guess = try;
                        }
                        /* South-Top */
                        if (X2 > 0 && time0[index - nxy] < 1.e9) {
                            s[0] = sc(X1, y2, X2);
                            /* Point 1 */
                            t[1] = tc(X1, y2, X2 - 1);
                            s[1] = sc(X1, y2, X2 - 1);
                            /* Point 2 */
                            t[2] = tc(X1, y2 - 1, X2);
                            s[2] = sc(X1, y2 - 1, X2);
                            /* Point 3 */
                            t[3] = tc(X1, y2 - 1, X2 - 1);
                            s[3] = sc(X1, y2 - 1, X2 - 1);

                            d01 = h;
                            d02 = rX2*dq;
                            d13 = rXm1*dq;
                            d23 = h;

                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) guess = try;
                        }
                        /* South-Bottom */
                        /*if(time0[index+nxy] < 1.e9 && X2<nz-1 )  { */
                        if (X2 < nz - 1) {
                            if (time0[index + nxy] < 1.e9) {
                                s[0] = sc(X1, y2, X2);
                                /* Point 1 */
                                t[1] = tc(X1, y2, X2 + 1);
                                s[1] = sc(X1, y2, X2 + 1);
                                /* Point 2 */
                                t[2] = tc(X1, y2 - 1, X2);
                                s[2] = sc(X1, y2 - 1, X2);
                                /* Point 3 */
                                t[3] = tc(X1, y2 - 1, X2 + 1);
                                s[3] = sc(X1, y2 - 1, X2 + 1);

                                d01 = h;
                                d02 = rX2*dq;
                                d13 = rXp1*dq;
                                d23 = h;

                                try = fdsph2d(t, s, d01, d02, d13, d23);
                                if (try < guess) guess = try;
                            }
                        }
                        /* Parallel - East Bottom */
                        /*if(time0[index+1] < 1.e9 && time0[index+nxy+1] < 1.e9
                           && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<nx-1 ) { */
                        if (X2 < nz - 1 && X1 < nx - 1) {
                            if (time0[index + 1] < 1.e9 && time0[index + nxy + 1] < 1.e9
                                    && time0[index + nxy] < 1.e9) {
                                s[0] = sc(X1, y2, X2);
                                /* Point 1 */
                                t[1] = tc(X1 + 1, y2, X2);
                                s[1] = sc(X1 + 1, y2, X2);
                                /* Point 2 */
                                t[2] = tc(X1, y2, X2 + 1);
                                s[2] = sc(X1, y2, X2 + 1);
                                /* Point 3 */
                                t[3] = tc(X1 + 1, y2, X2 + 1);
                                s[3] = sc(X1 + 1, y2, X2 + 1);
                                d01 = rX2 * sinqy2*df;
                                d02 = h;
                                d13 = h;
                                d23 = rXp1 * sinqy2*df;
                                try = fdsph2d(t, s, d01, d02, d13, d23);
                                if (try < guess) {
                                    fhead = (guess - try) / (d01 * slow0[index]);
                                    guess = try;
                                }
                            }
                        }
                        /* Parallel - East Top */
                        if (X2 > 0 && X1 < nx - 1 && time0[index + 1] < 1.e9 && time0[index - nxy + 1] < 1.e9
                                && time0[index - nxy] < 1.e9) {
                            s[0] = sc(X1, y2, X2);
                            /* Point 1 */
                            t[1] = tc(X1 + 1, y2, X2);
                            s[1] = sc(X1 + 1, y2, X2);
                            /* Point 2 */
                            t[2] = tc(X1, y2, X2 - 1);
                            s[2] = sc(X1, y2, X2 - 1);
                            /* Point 3 */
                            t[3] = tc(X1 + 1, y2, X2 - 1);
                            s[3] = sc(X1 + 1, y2, X2 - 1);
                            d01 = rX2 * sinqy2*df;
                            d02 = h;
                            d13 = h;
                            d23 = rXm1 * sinqy2*df;
                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) {
                                fhead = (guess - try) / (d01 * slow0[index]);
                                guess = try;
                            }
                        }
                        /* Parallel - West Bottom */
                        /*if(time0[index-1] < 1.e9 && time0[index+nxy-1] < 1.e9
                           && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) { */
                        if (X2 < nz - 1 && X1 > 0) {
                            if (time0[index - 1] < 1.e9 && time0[index + nxy - 1] < 1.e9
                                    && time0[index + nxy] < 1.e9) {
                                s[0] = sc(X1, y2, X2);
                                /* Point 1 */
                                t[1] = tc(X1 - 1, y2, X2);
                                s[1] = sc(X1 - 1, y2, X2);
                                /* Point 2 */
                                t[2] = tc(X1, y2, X2 + 1);
                                s[2] = sc(X1, y2, X2 + 1);
                                /* Point 3 */
                                t[3] = tc(X1 - 1, y2, X2 + 1);
                                s[3] = sc(X1 - 1, y2, X2 + 1);
                                d01 = rX2 * sinqy2*df;
                                d02 = h;
                                d13 = h;
                                d23 = rXp1 * sinqy2*df;
                                try = fdsph2d(t, s, d01, d02, d13, d23);
                                if (try < guess) {
                                    fhead = (guess - try) / (d01 * slow0[index]);
                                    guess = try;
                                }
                            }
                        }
                        /* Parallel - West Top */
                        if (X2 > 0 && X1 > 0 && time0[index - 1] < 1.e9 && time0[index - nxy - 1] < 1.e9
                                && time0[index - nxy] < 1.e9) {
                            s[0] = sc(X1, y2, X2);
                            /* Point 1 */
                            t[1] = tc(X1 - 1, y2, X2);
                            s[1] = sc(X1 - 1, y2, X2);
                            /* Point 2 */
                            t[2] = tc(X1, y2, X2 - 1);
                            s[2] = sc(X1, y2, X2 - 1);
                            /* Point 3 */
                            t[3] = tc(X1 - 1, y2, X2 - 1);
                            s[3] = sc(X1 - 1, y2, X2 - 1);
                            d01 = rX2 * sinqy2*df;
                            d02 = h;
                            d13 = h;
                            d23 = rXm1 * sinqy2*df;
                            try = fdsph2d(t, s, d01, d02, d13, d23);
                            if (try < guess) {
                                fhead = (guess - try) / (d01 * slow0[index]);
                                guess = try;
                            }
                        }

                        /*** New Face Trials ***/
                        if (guess > 1.0e9) {
                            if (X1 > x1 + 1 && X1 < x2 - 1 && X2 > z1 + 1 && X2 < z2 - 1) {
                                /* Point 0 */
                                t[0] = tc(X1, y2 - 1, X2 + 1);
                                s[0] = sc(X1, y2 - 1, X2 + 1);
                                /* Point 1 */
                                t[1] = tc(X1 - 1, y2 - 1, X2);
                                s[1] = sc(X1 - 1, y2 - 1, X2);
                                /* Point 2 */
                                t[2] = tc(X1, y2 - 1, X2);
                                s[2] = sc(X1, y2 - 1, X2);
                                /* Point 3 */
                                t[3] = tc(X1 + 1, y2 - 1, X2);
                                s[3] = sc(X1 + 1, y2 - 1, X2);
                                /* Point 4 */
                                t[4] = tc(X1, y2 - 1, X2 - 1);
                                s[4] = sc(X1, y2 - 1, X2 - 1);
                                /* Point 5 */
                                s[5] = sc(X1, y2, X2);
                                d02 = h;
                                d12 = rX2 * sin(qy2m1) * df;
                                d25 = rX2*dq;
                                try = fdsphnf(t, s, d02, d12, d25);
                                if (try < guess) guess = try;
                            }
                        }

                        /*** Edge Trials ***/
                        /* Change in Latitude */
                        try = tc(X1, y2 - 1, X2) + .5 * (sc(X1, y2 - 1, X2) + sc(X1, y2, X2)) * rX2*dq;
                        if (try < guess) guess = try;
                        /* Change in East Longitude */
                        if (X1 < nx - 1 && time0[index + 1] < 1.e9) {
                            try = tc(X1 + 1, y2, X2) + .5 * (sc(X1, y2, X2) + sc(X1 + 1, y2, X2)) * rX2 * df*sinqy2;
                            if (try < guess) {
                                fhead = (guess - try) / (rX2 * df * sinqy2 * slow0[index]);
                                guess = try;
                            }
                        }
                        /* Change in West Longitude */
                        if (X1 > 0 && time0[index - 1] < 1.e9) {
                            try = tc(X1 - 1, y2, X2) + .5 * (sc(X1, y2, X2) + sc(X1 - 1, y2, X2)) * rX2 * df*sinqy2;
                            if (try < guess) {
                                fhead = (guess - try) / (rX2 * df * sinqy2 * slow0[index]);
                                guess = try;
                            }
                        }
                        /* Change in depth (increasing) */
                        /*if ( time0[index+nxy]<1.e9 && X2<nz-1 )  { */
                        if (X2 < nz - 1) {
                            if (time0[index + nxy] < 1.e9) {
                                try = tc(X1, y2, X2 + 1) + .5 * (sc(X1, y2, X2) + sc(X1, y2, X2 + 1)) * h;
                                if (try < guess) {
                                    fhead = (guess - try) / (h * slow0[index]);
                                    guess = try;
                                }
                            }
                        }
                        /* Change in depth (decreasing) */
                        if (X2 > 0 && time0[index - nxy] < 1.e9) {
                            try = tc(X1, y2, X2 - 1) + .5 * (sc(X1, y2, X2) + sc(X1, y2, X2 - 1)) * h;
                            if (try < guess) {
                                fhead = (guess - try) / (h * slow0[index]);
                                guess = try;
                            }
                        }
                        if (guess < time0[index]) {
                            time0[index] = guess;
                            if (fhead > headtest) headw[4]++;
                        }
                    }
                    if (y2 == ny - 1) dy2 = 0;
                    y2++;
                }
                /* End y2 if condition */
            }
            /* End growth loop */


            /* UPDATE RADIUS */
            radius++;
            if (radius % 10 == 0) fprintf(stderr, "Completed radius = %d\n", radius);
            if (radius == maxrad) rad0 = 0;

        } /* END BIG LOOP */


        /* TEST IF REVERSE PROPAGATION IS NEEDED */

        if (headw[1] == 0 && headw[2] == 0 && headw[3] == 0 && headw[4] == 0
                && headw[5] == 0 && headw[6] == 0)
            reverse = 0;
        else {
            head = 0;
            if (headw[1] > 0) {
                fprintf(stderr, "Head waves found on left: %d\n", headw[1]);
                if (headw[1] > head) {
                    head = headw[1];
                    srcwall = 1;
                }
            }
            if (headw[2] > 0) {
                fprintf(stderr, "Head waves found on right: %d\n", headw[2]);
                if (headw[2] > head) {
                    head = headw[2];
                    srcwall = 2;
                }
            }
            if (headw[3] > 0) {
                fprintf(stderr, "Head waves found on front: %d\n", headw[3]);
                if (headw[3] > head) {
                    head = headw[3];
                    srcwall = 3;
                }
            }
            if (headw[4] > 0) {
                fprintf(stderr, "Head waves found on back: %d\n", headw[4]);
                if (headw[4] > head) {
                    head = headw[4];
                    srcwall = 4;
                }
            }
            if (headw[5] > 0) {
                fprintf(stderr, "Head waves found on top: %d\n", headw[5]);
                if (headw[5] > head) {
                    head = headw[5];
                    srcwall = 5;
                }
            }
            if (headw[6] > 0) {
                fprintf(stderr, "Head waves found on bottom: %d\n", headw[6]);
                if (headw[6] > head) {
                    head = headw[6];
                    srcwall = 6;
                }
            }
            if (headpref > 0 && headw[headpref] > 0) {
                fprintf(stderr, "Preference to restart on wall opposite source\n");
                srcwall = headpref;
            }
            /* SET LOCATIONS OF SIDES OF THE SPHERICAL ELEMENT SO THAT THE ELEMENT IS A FACE */
            dx1 = 1;
            dx2 = 1;
            dy1 = 1;
            dy2 = 1;
            dz1 = 1;
            dz2 = 1;
            rad0 = 1;
            radius = 1;
            if (srcwall == 1) {
                x2 = 1;
                fprintf(stderr, "RESTART at left side of model\n");
            } else {
                x2 = nx;
                dx2 = 0;
            }
            if (srcwall == 2) {
                x1 = nx - 2;
                fprintf(stderr, "RESTART at right side of model\n");
            } else {
                x1 = -1;
                dx1 = 0;
            }
            if (srcwall == 3) {
                y2 = 1;
                fprintf(stderr, "RESTART at front side of model\n");
            } else {
                y2 = ny;
                dy2 = 0;
            }
            if (srcwall == 4) {
                y1 = ny - 2;
                fprintf(stderr, "RESTART at back side of model\n");
            } else {
                y1 = -1;
                dy1 = 0;
            }
            if (srcwall == 5) {
                z2 = 1;
                fprintf(stderr, "RESTART at top side of model\n");
            } else {
                z2 = nz;
                dz2 = 0;
            }
            if (srcwall == 6) {
                z1 = nz - 2;
                fprintf(stderr, "RESTART at bottom side of model\n");
            } else {
                z1 = -1;
                dz1 = 0;
            }
            if (reverse == 0)
                fprintf(stderr, "WARNING:  RESTART CANCELLED by choice of input parameter 'reverse'\n");
        }
        reverse--;

    } /* END BIGGER LOOP */

    /* Debugging step - compare with constant velocity model */
    /*
            ii = 0;
            for (k = 0; k < nz; k++) {
              rz1 =   rc(k);
              for (j = 0; j < ny; j++) {
                qX2 =   qc(j);
                for (i = 0; i < nx; i++) {
                  fx1 =   fc(i);
                  x2p = rz1*cos(qX2);
                  y2p = rz1*sin(qX2)*cos(fx1);
                  z2p = rz1*sin(qX2)*sin(fx1);
                  delt =  DIST(x1p,y1p,z1p,x2p,y2p,z2p);
                  timec = s000*delt;
                  time0[ii] = fabs(100*(time0[ii] - timec)/timec);
                  ii  = ii + 1;
                }
              }
            }
     */
    /* End of debugging step */

    /* parse the time file if required */
    nxp = nx;
    nyp = ny;
    nzp = nz;
    if (parse > 1) {
        nzp = 0;
        ii = 0;
        for (k = 0; k < nz; k += parse) {
            nzp++;
            rz1 = rc(k);
            for (j = 0; j < ny; j += parse) {
                qX2 = qc(j);
                for (i = 0; i < nx; i += parse) {
                    time0[ii] = time0[i + j * nx + k * nxy];
                    /* Testing */
                    /*
                                  fx1 =   fc(i);
                                  x2p = rz1*cos(qX2);
                                  y2p = rz1*sin(qX2)*cos(fx1);
                                  z2p = rz1*sin(qX2)*sin(fx1);
                                  delt =  DIST(x1p,y1p,z1p,x2p,y2p,z2p);
                                  timec = s000*delt;
                                  if (k==nz-1) {
                                     fprintf(stderr," k, j, i, ii, x, y, z, tc, t0, dt = %d %d %d %d %g %g %g %g %g %g\n",
                                      k, j, i, ii, fc(i), qc(j), rc(k), timec, time0[ii], fabs(100*(time0[ii] - timec)/timec));
                                    }
                     */
                    /* End Testing*/
                    ii++;
                }
            }
        }
        nxp = 0;
        nyp = 0;
        for (j = 0; j < ny; j += parse) nyp++;
        for (i = 0; i < nx; i += parse) nxp++;
        nxyz = ii;
    }

    /* OUTPUT COMPLETED WAVEFRONT */
    strncpy(headout.header, headin.header, 120);
    headout.clat = headin.clat;
    headout.clon = headin.clon;
    headout.cz = headin.cz;
    headout.az = headin.az;
    headout.x0 = headin.x0;
    headout.y0 = headin.y0;
    headout.z0 = headin.z0;
    headout.dx = headin.dx*parse;
    headout.dy = headin.dy*parse;
    headout.dz = headin.dz*parse;
    headout.nx = nxp;
    headout.ny = nyp;
    headout.nz = nzp;
    /* Reassign the quantity */
    strcpy(&headout.header[12], "P");
    strcpy(&headout.header[13], "T");
    strcpy(&headout.header[14], "I");
    strcpy(&headout.header[15], "M");

    for (i = 16; i < 120; i++) {
        strcpy(&headout.header[i], &headin.header[i]);
    }
    // Adding a null probably not necessary here
    // strcpy (&headout.header[120], "\0");
    // fprintf(stderr, " Header  =  %s \n", &headout.header[0]);
    headout.fxs = fxss;
    headout.fys = fyss;
    headout.fzs = fzss;
    fprintf(stderr, " Float Header fxs, fys, fzs = %f %f %f \n", headout.fxs, headout.fys, headout.fzs);

    /* Swap back to proper byte order if required */
    if (swab == 2 || swab == 3) {
        fprintf(stdout, "Swapping bytes on output header and times ..\n");
        FIX_FLOAT(headout.az);
        FIX_INT(headout.nx);
        FIX_INT(headout.ny);
        FIX_INT(headout.nz);
        FixDouble(&headout.clat);
        FixDouble(&headout.clon);
        FixDouble(&headout.cz);
        FixDouble(&headout.x0);
        FixDouble(&headout.y0);
        FixDouble(&headout.z0);
        FixDouble(&headout.dx);
        FixDouble(&headout.dy);
        FixDouble(&headout.dz);
        FixDouble(&headout.fxs);
        FixDouble(&headout.fys);
        FixDouble(&headout.fzs);
        for (i = 0; i < nxyz; i++) FIX_FLOAT(time0[i]);
        fprintf(stdout, ".. Done\n");
    }
    fprintf(stdout, "fxss =  %g\n", headout.fxs);
    fprintf(stdout, "fyss =  %g\n", headout.fys);
    fprintf(stdout, "fzss =  %g\n", headout.fzs);

    // 20190405 AJL    write(tfint, &headout, 232);
    // 20190405 AJL    write(tfint, time0, nxyz * 4);
    fprintf(stderr, "wavefront done \n\n");
    writeNLLtimeGrid(&headout, timefile, time0);
    /* CONVERT TO VELOCITY */
    for (i = 0; i < nxyz; i++) slow0[i] = 1 / slow0[i];
    writeNLLmodelGrid(&headout, timefile, slow0);

    endpar(); // 20190405 AJL
}

#define PNAME  "sphfd_SWR_NLL"

/** function to write NonLinLoc format time grid buffer and header to disk
 *
 * 20190404 AJL - added
 *
 */

int writeNLLtimeGrid(struct vhead* headout, char *timefile, float *time0) {

    strcpy(prog_name, PNAME);

    char fn_time_output[MAXLINE_LONG];
    sprintf(fn_time_output, "%s", timefile);
    sprintf(MsgStr, "Creating time grid files: %s.time.*", fn_time_output);
    nll_putmsg(0, MsgStr);


    GridDesc time_grid;
    SourceDesc source;

    // set NLL constants
    SetConstants();

    // virtual NLL control file settings
    if (get_transform(0, "GLOBAL") < 0) {
        nll_puterr("ERROR: setting GLOBAL transformation parameters.");
    }


    // initialize 3D grid

    time_grid.type = GRID_TIME;
    strcpy(time_grid.chr_type, "TIME");
    strcpy(time_grid.float_type, "FLOAT");

    time_grid.numx = headout->nx;
    time_grid.numy = headout->ny;
    time_grid.numz = headout->nz;

    double x0, y0, z0;
    getpar("x0", "F", &x0);
    getpar("y0", "F", &y0);
    getpar("z0", "F", &z0);

    double dz, dq, df;
    mstpar("h", "F", &dz);
    getpar("dq", "F", &dq);
    getpar("df", "F", &df);
    if (dq <= FLT_MIN) {
        dq = dz * KM2DEG;
    }
    if (df <= FLT_MIN) {
        df = (dz / cos(y0 * DE2RA)) * KM2DEG;
    }
    time_grid.dx = df;
    time_grid.dy = dq;
    time_grid.dz = dz;

    time_grid.origx = x0;
    time_grid.origy = y0 - (double) (time_grid.numy - 1) * time_grid.dy;
    time_grid.origz = z0;

    time_grid.flagGridCascading = IS_NOT_CASCADING;
    time_grid.iSwapBytes = 0;

    display_grid_param(&time_grid);

    // must get source coordinates from original input because they were converted to the geocentric frame used internally
    strcpy(source.label, strrchr(timefile, '.') + 1);
    double fxs, fys, fzs;
    mstpar("fxs", "F", &fxs);
    mstpar("fys", "F", &fys);
    mstpar("fzs", "F", &fzs);
    source.dlat = fys;
    source.dlong = fxs;
    source.depth = fzs;
    source.is_coord_latlon = 1;
    source.is_coord_xyz = 0;
    ConvertASourceLocation(0, &source, 1, 0);

    // allocate time grid
    time_grid.buffer = AllocateGrid(&time_grid);
    if (time_grid.buffer == NULL) {
        nll_puterr(
                "ERROR: allocating memory for 3D slowness grid buffer.");
        exit(EXIT_ERROR_MEMORY);
    }
    // create array access pointers
    time_grid.array = CreateGridArray(&time_grid);
    if (time_grid.array == NULL) {
        nll_puterr(
                "ERROR: creating array for accessing 3D vel/slowness grid buffer.");
        exit(EXIT_ERROR_MEMORY);
    }

    // covert sphfd time grid to NLL grid
    int nx = headout->nx;
    int nxy = nx * headout->ny;

    int nwritten = 0;
    for (int ix = 0; ix < time_grid.numx; ix++) {
        printf("ix = %d/%d\r", ix, time_grid.numx);
        for (int iy = 0; iy < time_grid.numy; iy++) {
            int iy_sphfd = time_grid.numy - iy - 1; // sphfd has elements increasing in latitude from north to south, NLL is south to north.
            for (int iz = 0; iz < time_grid.numz; iz++) {
                ((GRID_FLOAT_TYPE ***) time_grid.array)[ix][iy][iz] = tc(ix, iy_sphfd, iz);
                nwritten++;
            }
        }
    }
    printf("\n");
    printf("%d/%d cell values written to NLL grid\n", nwritten, headout->nx * headout->ny * headout->nz);

    // save grid to disk

    if (WriteGrid3dBuf(&time_grid, &source, fn_time_output, "time") < 0) {
        nll_puterr("ERROR: writing slowness grid to disk.");
    }

    // write pseudo control file for plotting
    char fn_hdr_output[MAXLINE_LONG];
    sprintf(fn_hdr_output, "%s.time.in", timefile);
    sprintf(MsgStr, "Creating plotting pseudo control file: %s.time.*", fn_hdr_output);
    nll_putmsg(0, MsgStr);
    FILE *fp_hdr_output;
    if ((fp_hdr_output = fopen(fn_hdr_output, "w")) == NULL) {
        nll_putmsg2(1, "INFO: cannot open plotting pseudo control file", fn_hdr_output);
    }


    double xsource_igrid = (source.dlong - time_grid.origx) / time_grid.dx;
    double ysource_igrid = (source.dlat - time_grid.origy) / time_grid.dy;
    double zsource_igrid = (source.depth - time_grid.origz) / time_grid.dz;
    sprintf(MsgStr,
            "Source:  GridLoc: ix=%lf iy=%lf iz=%lf",
            xsource_igrid, ysource_igrid, zsource_igrid);
    nll_putmsg(0, MsgStr);

    fprintf(fp_hdr_output, "CONTROL 1 54321\n");
    fprintf(fp_hdr_output, "TRANSFORM  GLOBAL\n");
    fprintf(fp_hdr_output, "MAPTRANS AZIMUTHAL_EQUIDIST WGS-84 %f %f 0\n", time_grid.origy, time_grid.origx);
    fprintf(fp_hdr_output, "MAPGRID  %d %d %d   0 0 %f  %f %f %f TIME FLOAT\n",
            time_grid.numx, time_grid.numy, time_grid.numz,
            time_grid.origz,
            time_grid.dz, time_grid.dz, time_grid.dz);

    fclose(fp_hdr_output);

}

/** function to write NonLinLoc format model grid buffer and header to disk
 *
 * 20190410 AJL - added
 *
 */

int writeNLLmodelGrid(struct vhead* headout, char *modelfile, float *slow0) {

    strcpy(prog_name, PNAME);

    char fn_model_output[MAXLINE_LONG];
    sprintf(fn_model_output, "%s", modelfile);
    sprintf(MsgStr, "Creating model grid files: %s.mod.*", fn_model_output);
    nll_putmsg(0, MsgStr);


    GridDesc model_grid;

    // set NLL constants
    SetConstants();

    // virtual NLL control file settings
    if (get_transform(0, "GLOBAL") < 0) {
        nll_puterr("ERROR: setting GLOBAL transformation parameters.");
    }


    // initialize 3D grid

    model_grid.type = GRID_VELOCITY;
    strcpy(model_grid.chr_type, "VELOCITY");
    strcpy(model_grid.float_type, "FLOAT");

    model_grid.numx = headout->nx;
    model_grid.numy = headout->ny;
    model_grid.numz = headout->nz;

    double x0, y0, z0;
    getpar("x0", "F", &x0);
    getpar("y0", "F", &y0);
    getpar("z0", "F", &z0);

    double dz, dq, df;
    mstpar("h", "F", &dz);
    getpar("dq", "F", &dq);
    getpar("df", "F", &df);
    if (dq <= FLT_MIN) {
        dq = dz * KM2DEG;
    }
    if (df <= FLT_MIN) {
        df = (dz / cos(y0 * DE2RA)) * KM2DEG;
    }
    model_grid.dx = df;
    model_grid.dy = dq;
    model_grid.dz = dz;

    model_grid.origx = x0;
    model_grid.origy = y0 - (double) (model_grid.numy - 1) * model_grid.dy;
    model_grid.origz = z0;

    model_grid.flagGridCascading = IS_NOT_CASCADING;
    model_grid.iSwapBytes = 0;

    display_grid_param(&model_grid);

    // allocate model grid
    model_grid.buffer = AllocateGrid(&model_grid);
    if (model_grid.buffer == NULL) {
        nll_puterr(
                "ERROR: allocating memory for 3D model grid buffer.");
        exit(EXIT_ERROR_MEMORY);
    }
    // create array access pointers
    model_grid.array = CreateGridArray(&model_grid);
    if (model_grid.array == NULL) {
        nll_puterr(
                "ERROR: creating array for accessing 3D model grid buffer.");
        exit(EXIT_ERROR_MEMORY);
    }

    // covert sphfd model grid to NLL grid
    int nx = headout->nx;
    int nxy = nx * headout->ny;

    int nwritten = 0;
    for (int ix = 0; ix < model_grid.numx; ix++) {
        printf("ix = %d/%d\r", ix, model_grid.numx);
        for (int iy = 0; iy < model_grid.numy; iy++) {
            int iy_sphfd = model_grid.numy - iy - 1; // sphfd has elements increasing in latitude from north to south, NLL is south to north.
            for (int iz = 0; iz < model_grid.numz; iz++) {
                ((GRID_FLOAT_TYPE ***) model_grid.array)[ix][iy][iz] = sc(ix, iy_sphfd, iz);
                nwritten++;
            }
        }
    }
    printf("\n");
    printf("%d/%d cell values written to NLL grid\n", nwritten, headout->nx * headout->ny * headout->nz);

    // save grid to disk

    if (WriteGrid3dBuf(&model_grid, NULL, fn_model_output, "mod") < 0) {
        nll_puterr("ERROR: writing model grid to disk.");
    }

    // write pseudo control file for plotting
    char fn_hdr_output[MAXLINE_LONG];
    sprintf(fn_hdr_output, "%s.mod.in", modelfile);
    sprintf(MsgStr, "Creating plotting pseudo control file: %s", fn_hdr_output);
    nll_putmsg(0, MsgStr);
    FILE *fp_hdr_output;
    if ((fp_hdr_output = fopen(fn_hdr_output, "w")) == NULL) {
        nll_putmsg2(1, "INFO: cannot open plotting pseudo control file", fn_hdr_output);
    }

    fprintf(fp_hdr_output, "CONTROL 1 54321\n");
    fprintf(fp_hdr_output, "TRANSFORM  GLOBAL\n");
    fprintf(fp_hdr_output, "MAPTRANS AZIMUTHAL_EQUIDIST %f %f 0\n", model_grid.origy, model_grid.origx);
    fprintf(fp_hdr_output, "MAPGRID  %d %d %d   0 0 %f  %f %f %f VELOCITY FLOAT\n",
            model_grid.numx, model_grid.numy, model_grid.numz,
            model_grid.origz,
            model_grid.dz, model_grid.dz, model_grid.dz);
    fclose(fp_hdr_output);


    // get velocity profile under source
    char fn_vprof_output[MAXLINE_LONG];
    sprintf(fn_vprof_output, "%s.mod.profile", modelfile);
    sprintf(MsgStr, "Creating model profile file: %s", fn_vprof_output);
    nll_putmsg(0, MsgStr);
    FILE *fp_vprof_output;
    if ((fp_vprof_output = fopen(fn_vprof_output, "w")) == NULL) {
        nll_putmsg2(1, "INFO: cannot open model profile file", fn_vprof_output);
    }
        fprintf(fp_vprof_output, "depth velocity\n");
    double fxs, fys, fzs;
    mstpar("fxs", "F", &fxs);
    mstpar("fys", "F", &fys);
    mstpar("fzs", "F", &fzs);
    GRID_FLOAT_TYPE value;
    double zloc;
    for (int iz = 0; iz < model_grid.numz; iz++) {
        zloc = model_grid.origz + (double) (iz) * model_grid.dz;
        value = ReadAbsInterpGrid3d(NULL, &model_grid, fxs, fys, zloc, 0);
        fprintf(fp_vprof_output, "%f %f\n", zloc, value);
    }
    fclose(fp_vprof_output);

}
