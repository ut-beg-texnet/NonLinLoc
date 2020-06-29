/*----------------------------------------------------------------------------*/
/*  TIME_3D: FINITE DIFFERENCE COMPUTATION OF FIRST ARRIVAL TIMES IN 3D.      */
/*----------------------------------------------------------------------------*/
/*  P.PODVIN, Geophysique, Ecole des Mines de Paris, Fontainebleau.           */
/*  E-MAIL: podvin@geophy.ensmp.fr          Tel: 33-(1) 64 69 49 25.          */
/*  April 1991, last revision: 7 May 1993                                     */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* Copyright (C) 1993 Pascal Podvin   <pascal@geophy.ensmp.fr>                */
/*                                                                            */
/* 1997-1999; incorporated in NonLinLoc software by Anthony Lomax             */
/*                                                                            */
/* (Please contact P. Podvin if you intend to use this code for use other     */
/*    than with the NonLinLoc software package)                               */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*  This is a plain K&R C (not ANSI-C) code.                                  */
/*                                                                            */
/*  USAGE: (from a program written in C or in FORTRAN)                        */
/*       (int)i_status=time_3d(HS,T,NX,NY,NZ,XS,YS,ZS,HS_EPS_INIT,MSG)        */
/*                                                                            */
/*  ARGUMENTS: (C description. All FORTRAN arguments are pointers)            */
/*    (int)      NX,NY,NZ : Dimensions of the timefield. Cells are cubic.     */
/*                          Real model dimensions are NX-1,NY-1,NZ-1.         */
/*    (GRID_FLOAT_TYPE *)      HS,T : 1-D arrays. Length: NX*NY*NZ                      */
/*                          HS[] is interpreted as slowness*mesh parameter,   */
/*                          so that associated physical dimension is Time.    */
/*                          Both arrays must be arranged as a succession of   */
/*                          planes x (resp. z)=Ct, each plane consisting of   */
/*                          a suite of y=Ct vectors when this program is      */
/*                          called from a program written in C (resp.FORTRAN).*/
/*                          hs and t are 3D pointers allocated in order       */
/*                          to address points and cells straightforwardly.    */
/*                          Implicitly staggered grids are used: t[][][]      */
/*                          represent arrival times at grid-points, whereas   */
/*                          hs[][][] describe (constant) slownesses in cells. */
/*                          Cells such that x=NX-1 or y=NY-1 or z=NZ-1 are    */
/*                          thus dummy cells located OUT of the model.        */
/*                          Their associated slowness is treated as INFINITY  */
/*                          (absorbing boundary conditions).                  */
/*                          Negative hs values are unphysical and forbidden.  */
/*    (GRID_FLOAT_TYPE)    XS,YS,ZS : real point source coordinates, i.e. distances     */
/*                          along axes between source point and corner 0,0,0  */
/*                          of the timefield, measured in h units.            */
/*                          e.g. 0.0<=XS<=(GRID_FLOAT_TYPE)NX                           */
/*                          If source coordinates are illicit (out of model   */
/*                          bounds), the program uses T as a pre-initialized  */
/*                          timefield (multiple source, e.g. exploding        */
/*                          reflector). At least one element of T must then   */
/*                          be finite.                                        */
/*    (GRID_FLOAT_TYPE) HS_EPS_INIT : fraction (typically 1.0E-3) defining the toler-   */
/*                          ated model inhomogeneity for exact initialization.*/
/*                          A tolerance larger than 0.01 will potentially     */
/*                          create errors larger than those involved by the   */
/*                          F.D. scheme without any exact initialization.     */
/*                          (this arg. is used only in the point source case) */
/*    (int)           MSG : Message flag (0:silent,1:few messages,2:verbose)  */
/*                          A negative value inhibits "clever" initialization.*/
/*                                                                            */
/*  VALUE: time_3d() returns a nonzero integer if an error occurred.          */
/*      In such case, an explicit message is printed on 'stderr'.             */
/*                                                                            */
/*  FORTRAN INTERFACE: see function time_3d_(). Standard FORTRAN conventions  */
/*      are adopted: indexes start at 1 (but this is irrelevant here!),       */
/*      and arrays are transposed to fit FORTRAN-style memory mapping.        */
/*      Please note that time_3d is an INTEGER FUNCTION which name DOES       */
/*      NOT BEGIN WITH [I,J,K,L,M,N] ! You must declare it explicitly !       */
/*      e.g.: MAIN program                                                    */
/*            parameter (MAXNODES=...)                                        */
/*            real*4 hs(MAXNODES),t(MAXNODES)                                 */
/*            integer*4 time_3d                                               */
/*            nx=... ny=... nz=...                        ! #grid-points      */
/*            xs=...(0<=xs<=nx-1) ys=... zs=...           ! source point      */
/*            eps_hs=0.001                                ! tolerance         */
/*            call build_model(hs,nx,ny,nz)                                   */
/*            itime_3d_status=time_3d(hs,t,nx,ny,nz,xs,ys,zs,eps_hs,msg)      */
/*            if(itime_3d_status.ne.0) ----> an error occurred !              */
/*            ...                                                             */
/*            stop                                                            */
/*            end                                                             */
/*                                                                            */
/*            SUBROUTINE build_model(hs,nx,ny,nz)                             */
/*            real*4 hs(nx,ny,nz)                                             */
/*            nested loops: k=...(1<=k<=nz-1) j=... i=... ! cell coordinates  */
/*                  hs(i,j,k)=...                                             */
/*            end loops                                                       */
/*            return                                                          */
/*            end                                                             */
/*                                                                            */
/*      With SunOS, the FORTRAN compiler appends the final '_' to FORTRAN     */
/*      function name 'time_3d'. Full name 'time_3d_' may be necessary on     */
/*      other systems. Remember that NO standard C/Fortran interface exists ! */
/*                                                                            */
/*  PREPROCESSOR OPTIONS:       (/bin/cpp options applicable at compile time) */
/*      -DNO_IEEE_PROTOCOL : for systems where IEEE standards are not defined */
/*      -DINIT_MIN=value   :         these two options control initialization */
/*      -DINIT_RECURS_LIMIT=value     (see comments after #ifndef directives) */
/*                                                                            */
/*  BUGS :                                                                    */
/*       All known bugs are ascribed to bad argument list  or C/Fortran       */
/*       interface problems . Most frequent sources of confusing errors :     */
/*       -incorrect argument types (NX MUST be an integer, XS a real number)  */
/*       -incorrect initialization of array T when using a multiple source    */
/*       (e.g. T being initialized with zeroes,if source coordinates are out  */
/*        of model boundaries, then array T will be left unchanged !)         */
/*                                                                            */
/*  REFERENCE: Podvin & Lecomte, Geophys.J.Intern. 105, 271-284, 1991.        */
/*----------------------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>

// 20100614 AJL -
// the following sets globally the float size (float or double) for NLL grids - modify with care!
// normal NonLinLoc usage currenlty require float type
// takw-off angle support requires float type
#ifdef GRID_FLOAT_TYPE_DOUBLE
#define GRID_FLOAT_TYPE double
#else
#define GRID_FLOAT_TYPE float
#endif

#ifndef M_SQRT2
#define M_SQRT2     1.41421356237309504880
#endif
#ifndef M_SQRT3
#define M_SQRT3     1.732050807568877076
#endif


#ifdef NO_IEEE_PROTOCOL
#ifndef INFINITY
#define INFINITY    1.0e+19
    /* INFINITY should be < sqrt(MAX_FLOAT_VALUE), machine-dependent */
#endif
#define ISINF(x)    ((x)>1.0e+18)
#define NINT(x)     (int)floor((x)+0.5)
    /* NINT should exactly be "Nearest INTeger" */
#else
#ifdef sun
/* Note : do NOT use option -Xs with Sun cc as this option UNDEFINES "sun" ! */
#include <sunmath.h>
#endif
#ifndef INFINITY
#define INFINITY    (GRID_FLOAT_TYPE)infinity()
#endif
#define ISINF(x)    isinf(x)
#define NINT(x)     nint(x)
#endif /* NO_IEEE_PROTOCOL */

#define min(x,y) (((x)<(y)) ? (x):(y))
#define max(x,y) (((x)>(y)) ? (x):(y))
#define min3(x,y,z) (min(x,min(y,z)))
#define min4(x,y,z,t) (min(x,min(y,min(z,t))))

# ifdef __APPLE__
#include <stdlib.h>
# else
#include <malloc.h>
#endif

/*-------------------------------------Static functions-----------------------*/

static int
    pre_init(void ),
    init_point(void ),
    recursive_init(void ),
    propagate_point(int ),
    x_side(int,int,int,int,int,int ),
    y_side(int,int,int,int,int,int ),
    z_side(int,int,int,int,int,int ),
    scan_x_ff(int,int,int,int,int,int,int ),
    scan_x_fb(int,int,int,int,int,int,int ),
    scan_x_bf(int,int,int,int,int,int,int ),
    scan_x_bb(int,int,int,int,int,int,int ),
    scan_y_ff(int,int,int,int,int,int,int ),
    scan_y_fb(int,int,int,int,int,int,int ),
    scan_y_bf(int,int,int,int,int,int,int ),
    scan_y_bb(int,int,int,int,int,int,int ),
    scan_z_ff(int,int,int,int,int,int,int ),
    scan_z_fb(int,int,int,int,int,int,int ),
    scan_z_bf(int,int,int,int,int,int,int ),
    scan_z_bb(int,int,int,int,int,int,int );
    /* the only fully commented "side" functions are x_side() ans scan_x_ff() */

static void
    error(int ),
    init_nearest(void ),
    init_cell(GRID_FLOAT_TYPE,GRID_FLOAT_TYPE,GRID_FLOAT_TYPE,int,int,int ),
    free_ptrs(int );

static GRID_FLOAT_TYPE
    exact_delay(GRID_FLOAT_TYPE,GRID_FLOAT_TYPE,GRID_FLOAT_TYPE,int,int,int );

static int
    t_1d(int,int,int,GRID_FLOAT_TYPE,GRID_FLOAT_TYPE,GRID_FLOAT_TYPE,GRID_FLOAT_TYPE,GRID_FLOAT_TYPE ),
    t_2d(int,int,int,GRID_FLOAT_TYPE,GRID_FLOAT_TYPE,GRID_FLOAT_TYPE,GRID_FLOAT_TYPE ),
    diff_2d(int,int,int,GRID_FLOAT_TYPE,GRID_FLOAT_TYPE,GRID_FLOAT_TYPE ),
    t_3d_(int,int,int,GRID_FLOAT_TYPE,GRID_FLOAT_TYPE,GRID_FLOAT_TYPE,GRID_FLOAT_TYPE,GRID_FLOAT_TYPE,int ),
    t_3d_part1(int,int,int,GRID_FLOAT_TYPE,GRID_FLOAT_TYPE,GRID_FLOAT_TYPE,GRID_FLOAT_TYPE ),
    point_diff(int,int,int,GRID_FLOAT_TYPE,GRID_FLOAT_TYPE ),
    edge_diff(int,int,int,GRID_FLOAT_TYPE,GRID_FLOAT_TYPE,GRID_FLOAT_TYPE );

/*-------------------------------------Static variables-----------------------*/

/* MODEL */

static int
    nmesh_x,nmesh_y,nmesh_z;          /* Model dimensions (cells) */
static GRID_FLOAT_TYPE
    ***hs,*hs_buf,                    /* 1D and 3D arrays */
    *hs_keep=(GRID_FLOAT_TYPE *)NULL;           /* to save boundary values */

/* TIMEFIELD */

static int
    nx,ny,nz;                         /* Timefield dimensions (nodes) */
static GRID_FLOAT_TYPE
    ***t,*t_buf;                      /* 1D and 3D arrays */

/* SOURCE */

static GRID_FLOAT_TYPE
    fxs,fys,fzs;                      /* Point source coordinates */
static int
    xs,ys,zs;                         /* Nearest node */
static int
    mult=0;                           /* Flag used for multiple source */

/* PARAMETERS */

#ifndef INIT_MIN
#define INIT_MIN             7
#endif /* INIT_MIN */
#define N_INIT_X    (4*INIT_MIN+3)
#define N_INIT      N_INIT_X*N_INIT_X*N_INIT_X
                                      /* This conventional value defines  */
                                      /* the maximum size of the box that */
                                      /* will be initialized recursively. */
                                      /* (ADJUSTABLE at compile time).    */
                                      /* Default value is reasonable.     */
                                      /* Cost is already high: 3584 to    */
                                      /* 26416 more points according to   */
                                      /* source position.                 */

#ifndef INIT_RECURS_LIMIT
#define INIT_RECURS_LIMIT    1
#endif /* INIT_RECURS_LIMIT */
                                      /* This parameter defines the maximal */
                                      /* level of recursivity during init.  */
                                      /* (ADJUSTABLE at compile time).      */
                                      /* Value zero would practically rest- */
                                      /* rain initialization to the source  */
                                      /* point in heterogeneous models.     */
                                      /* Value 2 is advisable only when     */
                                      /* VERY severe heterogeneities are    */
                                      /* located close to the source point. */

static int
    messages,                         /* message flag (0:silent)              */
    source_at_node=0,                 /* are source coordinate int's ? (0/1)  */
    no_init=0,                        /* 1: inhibition of "clever" init.      */
    init_stage=0,                     /* level of recursivity during init.    */
    current_side_limit,               /* actual boundary of computations      */
    X0,X1,Y0,Y1,Z0,Z1,                /* inclusive boundaries of timed region */
    sum_updated,                      /* total count of adopted FD stencils   */
    reverse_order,                    /* level of recursivity in FD scheme    */
    *longflags,                       /* local headwave flags.                */
    flag_fb,x_start_fb,y_start_fb,z_start_fb,
    flag_bf,x_start_bf,y_start_bf,z_start_bf,
    flag_ff,x_start_ff,y_start_ff,z_start_ff,
    flag_bb,x_start_bb,y_start_bb,z_start_bb;
                                      /* control current side scanning.       */

static GRID_FLOAT_TYPE
    hs_eps_init;                      /* tolerance on homogeneity
                                       (fraction of slowness at source point) */

#define SMALLTALK        messages
#define VERBOSE          messages==2
/*------------------------------------------------Error flags---------------*/

#define NO_ERROR         0
#define ERR_MULT         (-1)
#define ERR_INF          (-2)
#define ERR_RECURS       (-3)
#define ERR_MALLOC       (-4)
#define ERR_HS_EPS       (-5)
#define ERR_NONPHYSICAL  (-6)

static char *err_msg[]=
            {
                "\ntime_3d: Computations terminated normally.\n",
                "\ntime_3d: Multiple source but no source at finite time.\n",
                "\ntime_3d: Source point is in a zero velocity zone.\n",
                "\ntime_3d: Fatal error during recursive init.\n",
                "\ntime_3d: Memory Allocation failed.\n",
                "\ntime_3d: Init: Illegal tolerance on inhomogeneity.\n",
                "\ntime_3d: Illegal negative slowness value.\n"
            };

/*-------------------------------------------------Error()------------------*/

static void
error(flag)

int flag;

{
    fflush(stdout);
    if(messages || flag) fprintf(stderr,"%s",err_msg[-flag]);
}

/*-------------------------------------------------Time_3d()----------------*/

// 20100317 ALomax - converted to ANSI C call, previously arguments were not passed correctly in Mac OS X (i686-apple-darwin10-gcc-4.2.1)
int time_3d(GRID_FLOAT_TYPE *HS, GRID_FLOAT_TYPE *T, int NX, int NY, int NZ, GRID_FLOAT_TYPE XS, GRID_FLOAT_TYPE YS, GRID_FLOAT_TYPE ZS, GRID_FLOAT_TYPE HS_EPS_INIT, int MSG)

/*GRID_FLOAT_TYPE *HS,*T;
GRID_FLOAT_TYPE HS_EPS_INIT,XS,YS,ZS;
int   NX,NY,NZ,MSG;*/

/* This function merely does nothing else than copying its arguments  */
/* to internal (static) variables. This allows you to alter the user  */
/* interface (e.g., pass on vector N[3] instead of NX,NY,NZ, or use   */
/* members of a struct or whatever else fits your needs) very easily. */
/* Arrays are passed as 1D vectors in order to make life easier with  */
/* Fortran calling programs...                                        */

{
    int    signal;

    //fprintf(stdout, "DEBUG: NX %d  NY %d  NZ %d  XS %.3f  YS %.3f  ZS %.3f  HS_EPS_INIT %.2e\n", NX, NY, NZ, XS, YS, ZS, HS_EPS_INIT);

/* copy args (with a few preliminary tests) to internal variables */
    hs_buf=HS;
    t_buf=T;
    nx=NX;
    ny=NY;
    nz=NZ;
    fxs=XS;
    fys=YS;
    fzs=ZS;
    hs_eps_init=HS_EPS_INIT;
    if(hs_eps_init<0.0 || hs_eps_init>1.0) {
        error(ERR_HS_EPS);
        return(ERR_HS_EPS);
    }
    if(MSG<0){
        no_init=1;
        MSG= -MSG;
    }/* this trick inhibits search for homogeneous region around the source */
    messages=MSG;

    //fprintf(stdout, "DEBUG: nx %d  ny %d  nz %d  fxs %.3f  fys %.3f  fzs %.3f  hs_eps_init %.2e\n", nx, ny, nz, fxs, fys, fzs, hs_eps_init);

    /* compute */
    if((signal=pre_init())==NO_ERROR){
        signal=propagate_point(init_point());
        free_ptrs(nx);
    }
    if(init_stage==0 || signal!=NO_ERROR) error(signal);
    return(signal);
}

/*---------------------------- Time_3d_(): FORTRAN INTERFACE ---------------*/

int
time_3d_(HS,T,NZ,NY,NX,ZS,YS,XS,HS_EPS_INIT,MSG)

/* All FORTRAN arguments are pointers. Dimensions X and Z are */
/* swapped in order to fit FORTRAN mapping conventions.       */

GRID_FLOAT_TYPE *HS,*T,*HS_EPS_INIT,*XS,*YS,*ZS;
int   *NX,*NY,*NZ,*MSG;

{
    int     signal;

    hs_buf=HS;
    t_buf=T;
    nx= *NX;
    ny= *NY;
    nz= *NZ;
    fxs= *XS;
    fys= *YS;
    fzs= *ZS;
    hs_eps_init= *HS_EPS_INIT;
    if(hs_eps_init<0.0 || hs_eps_init>1.0) {
        error(ERR_HS_EPS);
        return(ERR_HS_EPS);
    }
    if(*MSG<0){
        no_init=1;
        messages= -(*MSG);
    }
    else messages= *MSG;

    if((signal=pre_init())==NO_ERROR){
        signal=propagate_point(init_point());
        free_ptrs(nx);
    }
    if(init_stage==0 || signal!=NO_ERROR) error(signal);
    return(signal);
}

/*------------------------------------------------Pre_init()----------------*/

static int
pre_init()

{
    int
        x,y,z,
        np,nt,
        n0,n1;
    GRID_FLOAT_TYPE
        *pf;

    nmesh_x=nx-1;
    nmesh_y=ny-1;
    nmesh_z=nz-1;
    np=ny*nz;
    nt=nx*np;
    n1=max(nx,ny);
    n0=max(nx,nz);
    if(n1==n0) n0=max(ny,nz);
    n1*=n0;

/* allocate pointers */
    if(!(hs=(GRID_FLOAT_TYPE ***)malloc((unsigned)nx*sizeof(GRID_FLOAT_TYPE **))))
        return(ERR_MALLOC);
    if(!(t =(GRID_FLOAT_TYPE ***)malloc((unsigned)nx*sizeof(GRID_FLOAT_TYPE **)))){
        free((char *)hs);
        return(ERR_MALLOC);
    }
    if(!(longflags =(int *)malloc((unsigned)n1*sizeof(int)))){
        free((char *)t);
        free((char *)hs);
        return(ERR_MALLOC);
    }/* size of the largest side of the model */
    for(x=0;x<nx;x++)
        if(   !(hs[x]=(GRID_FLOAT_TYPE **)malloc((unsigned)ny*sizeof(GRID_FLOAT_TYPE *)))
           || !(t[x] =(GRID_FLOAT_TYPE **)malloc((unsigned)ny*sizeof(GRID_FLOAT_TYPE *)))
            ){
                free_ptrs(x);
                return(ERR_MALLOC);
            }
    for(x=0;x<nx;x++)
        for(y=0;y<ny;y++){
            hs[x][y]=hs_buf+x*np+y*nz;
            t[x][y]=t_buf+x*np+y*nz;
        }

/* stop here if recursive call */
    if(init_stage) return(NO_ERROR);

/* initialize all times as INFINITY if licit point source */
    if(fxs>=0.0 && fxs<=nx-1 && fys>=0 && fys<=ny-1 && fzs>=0 && fzs<=nz-1)
        for(x=0,pf=t_buf;x<nt;x++) *pf++=INFINITY;

/* assign INFINITY to hs in dummy meshes (x=nmesh_x|y=nmesh_y|z=nmesh_z) */
/* and keep masked values in hs_keep[].                                  */
    x=((nx+1)*(ny+1)+(nx+1)*nz+nz*ny)*sizeof(GRID_FLOAT_TYPE);
    if(!(hs_keep=(GRID_FLOAT_TYPE *)malloc((unsigned)x))) {
        free_ptrs(nx);
        return(ERR_MALLOC);
    }
    pf=hs_keep;
    for(x=0;x<nx;x++){
        for(y=0;y<ny;y++) {
            *pf++=hs[x][y][nmesh_z];
            hs[x][y][nmesh_z]=INFINITY;
        }
        for(z=0;z<nmesh_z;z++) {
            *pf++=hs[x][nmesh_y][z];
            hs[x][nmesh_y][z]=INFINITY;
        }
    }
    for(y=0;y<nmesh_y;y++)
        for(z=0;z<nmesh_z;z++) {
            *pf++=hs[nmesh_x][y][z];
            hs[nmesh_x][y][z]=INFINITY;
        }

/* test for negative slowness value */
    for(x=0,pf=hs_buf;x<nx*ny*nz;x++,pf++)
        if(*pf<0.0){
            free_ptrs(nx);
            return(ERR_NONPHYSICAL);
        }/* a negative value would provoke an infinitely recursive call */
         /* and act as a "black hole" driving all times to -INFINITY !! */

    return(NO_ERROR);
}

/*------------------------------------------------Init_point()--------------*/

static int
init_point()

{
    int
        signal=NO_ERROR,
        x,y,z,
        test,
        t_X0,t_X1,t_Y0,t_Y1,t_Z0,t_Z1;
    GRID_FLOAT_TYPE
        min_t,
        hs0,
        allowed_delta_hs,
        dist;

    // 20100204 AJL Satriano Bug Fix.
    hs0 = 0.0;

/* test relevance of source position or locate minimum time source point */
    if(fxs<0.0 || fxs>nx-1 || fys<0 || fys>ny-1 || fzs<0 || fzs>nz-1){
        for(x=0,min_t=INFINITY;x<nx;x++)
            for(y=0;y<ny;y++)
                for(z=0;z<nz;z++)
                    if(t[x][y][z]<min_t){
                        min_t=t[x][y][z];
                        xs=x;
                        ys=y;
                        zs=z;
                    }
        if(ISINF(min_t)) return(ERR_MULT);
        source_at_node=1;
        mult=1;
        if(SMALLTALK)
            printf("\nMultiple source starting at node [%d,%d,%d] at time %g.",
                xs,ys,zs,min_t);
    }

    else {
/* locate node closest to source */
        xs=NINT(fxs);
        ys=NINT(fys);
        zs=NINT(fzs);
        if(xs==fxs && ys==fys && zs==fzs){
            source_at_node=1;
            if(SMALLTALK) printf("\nSource located exactly at node [%d,%d,%d].",
                    xs,ys,zs);
        }
        mult=0;
    }

/* test relevance of slowness at the vicinity of the source */
/* (not tested in the case of multiple source)              */
    if(source_at_node){
        hs0=hs[xs][ys][zs];
        if(ISINF(hs0) && zs) hs0=hs[xs][ys][zs-1];
        if(ISINF(hs0) && ys) hs0=hs[xs][ys-1][zs];
        if(ISINF(hs0) && xs) hs0=hs[xs-1][ys][zs];
        if(ISINF(hs0) && zs && ys) hs0=hs[xs][ys-1][zs-1];
        if(ISINF(hs0) && zs && xs) hs0=hs[xs-1][ys][zs-1];
        if(ISINF(hs0) && ys && xs) hs0=hs[xs-1][ys-1][zs];
        if(ISINF(hs0) && zs && ys && xs) hs0=hs[xs-1][ys-1][zs-1];
    }
    else if(!mult){
        x=(fxs<xs) ? xs-1:xs;
        y=(fys<ys) ? ys-1:ys;
        z=(fzs<zs) ? zs-1:zs;
        hs0=hs[x][y][z];
        if(ISINF(hs0) && fxs==xs && xs){
            hs0=hs[x-1][y][z];
            if(ISINF(hs0) && fys==ys && ys) hs0=hs[x-1][y-1][z];
            if(ISINF(hs0) && fzs==zs && zs) hs0=hs[x-1][y][z-1];
        }
        if(ISINF(hs0) && fys==ys && ys){
            hs0=hs[x][y-1][z];
            if(ISINF(hs0) && fzs==zs && zs) hs0=hs[x][y-1][z-1];
        }
        if(ISINF(hs0) && fzs==zs && zs) hs0=hs[x][y][z-1];
    }
    if(ISINF(hs0)) return(ERR_INF);

/* if source is multiple, do not initialize at all the timefield */
    if(mult){
        X0=X1=xs;
        Y0=Y1=ys;
        Z0=Z1=zs;
        return(NO_ERROR);
    }/* this case has higher priority than the no_init directive */

/* if asked for, use minimal initialization : t=0.0 at the source point */
    if(no_init){
        X0=X1=xs;
        Y0=Y1=ys;
        Z0=Z1=zs;
        init_nearest();
        return(NO_ERROR);
    }

/* initialize inclusive boundaries of explored region */
    X0=max(xs-1,0);
    Y0=max(ys-1,0);
    Z0=max(zs-1,0);
    X1=min(xs+1,nmesh_x-1);
    Y1=min(ys+1,nmesh_y-1);
    Z1=min(zs+1,nmesh_z-1);

/* search largest parallelepipedic homogeneous box centered on the source */
    t_X0=t_X1=t_Y0=t_Y1=t_Z0=t_Z1=0;
    /* these flags will signal that a heterogeneity has been reached */
    allowed_delta_hs=hs0*hs_eps_init;
    /* defines tolerated inhomogeneity for exact initialization */
    do{
        test=0;
        if(X0 && !t_X0){
            test++;
            x=X0;
            for(y=Y0;y<=Y1 && y<nmesh_y && !t_X0;y++)
                for(z=Z0;z<=Z1 && z<nmesh_z && !t_X0;z++)
                    if(fabs(hs[x][y][z]-hs0)>allowed_delta_hs) t_X0=1;
            if(!t_X0) X0--;
        }
        if(Y0 && !t_Y0){
            test++;
            y=Y0;
            for(x=X0;x<=X1 && x<nmesh_x && !t_Y0;x++)
                for(z=Z0;z<=Z1 && z<nmesh_z && !t_Y0;z++)
                    if(fabs(hs[x][y][z]-hs0)>allowed_delta_hs) t_Y0=1;
            if(!t_Y0) Y0--;
        }
        if(Z0 && !t_Z0){
            test++;
            z=Z0;
            for(x=X0;x<=X1 && x<nmesh_x && !t_Z0;x++)
                for(y=Y0;y<=Y1 && y<nmesh_y && !t_Z0;y++)
                    if(fabs(hs[x][y][z]-hs0)>allowed_delta_hs) t_Z0=1;
            if(!t_Z0) Z0--;
        }
        if(X1<nmesh_x && !t_X1){
            test++;
            X1++;
            x=X1;
            for(y=Y0;y<=Y1 && y<nmesh_y && !t_X1;y++)
                for(z=Z0;z<=Z1 && z<nmesh_z && !t_X1;z++)
                    if(fabs(hs[x][y][z]-hs0)>allowed_delta_hs) t_X1=1;
        }
        if(Y1<nmesh_y && !t_Y1){
            test++;
            Y1++;
            y=Y1;
            for(x=X0;x<=X1 && x<nmesh_x && !t_Y1;x++)
                for(z=Z0;z<=Z1 && z<nmesh_z && !t_Y1;z++)
                    if(fabs(hs[x][y][z]-hs0)>allowed_delta_hs) t_Y1=1;
        }
        if(Z1<nmesh_z && !t_Z1){
            test++;
            Z1++;
            z=Z1;
            for(x=X0;x<=X1 && x<nmesh_x && !t_Z1;x++)
                for(y=Y0;y<=Y1 && y<nmesh_y && !t_Z1;y++)
                    if(fabs(hs[x][y][z]-hs0)>allowed_delta_hs) t_Z1=1;
        }
    } while(test);

    if(X0) X0++;
    if(Y0) Y0++;
    if(Z0) Z0++;
    if(X1<nmesh_x) X1--;
    if(Y1<nmesh_y) Y1--;
    if(Z1<nmesh_z) Z1--;
    /* limits are decremented so that interfaces where heterogeneities     */
    /* were detected are dealt with by finite differences (cf. headwaves). */
    /* (but this is not necessary when located at model boundaries !)      */

    if( init_stage>=INIT_RECURS_LIMIT ||
        (   (X0==0 || (xs-X0)>=INIT_MIN) &&
            (Y0==0 || (ys-Y0)>=INIT_MIN) &&
            (Z0==0 || (zs-Z0)>=INIT_MIN) &&
            (X1==nmesh_x || (X1-xs)>=INIT_MIN) &&
            (Y1==nmesh_y || (Y1-ys)>=INIT_MIN) &&
            (Z1==nmesh_z || (Z1-zs)>=INIT_MIN)     )   ) {
        if((X1-X0+1)*(Y1-Y0+1)*(Z1-Z0+1)==1)
            init_nearest();
        else for(x=X0;x<=X1;x++)
            for(y=Y0;y<=Y1;y++)
                for(z=Z0;z<=Z1;z++){
                    dist=(x-fxs)*(x-fxs)+(y-fys)*(y-fys)+(z-fzs)*(z-fzs);
                    t[x][y][z]=hs0*sqrt(dist);
                }
        if(SMALLTALK)
            printf("\nHomogeneous region: x[%d->%d] y[%d->%d] z[%d->%d]\n",
                    X0,X1,Y0,Y1,Z0,Z1);
    } /* if smallest distance from source to boundaries of the homogeneous */
      /* box exceeds conventional limit INIT_MIN, OR if no further recursi-*/
      /* vity is allowed, then exact arrivals are computed in this region. */

    else {
        if((signal=recursive_init())!=NO_ERROR) return(signal);
        X0=max(xs-INIT_MIN,0);
        Y0=max(ys-INIT_MIN,0);
        Z0=max(zs-INIT_MIN,0);
        X1=min(xs+INIT_MIN,nmesh_x);
        Y1=min(ys+INIT_MIN,nmesh_y);
        Z1=min(zs+INIT_MIN,nmesh_z);
    } /* otherwise, time_3d() is used recursively   */
      /* on a re-discretized (2*INIT_MIN+1)^3 cube. */

    return(signal);

}

/*------------------------------------------------Init_nearest()------------*/

static void
init_nearest()

/* initialize the 8|12|18 nearest neighbour nodes of the source    */
/* according to source position (inside a mesh or at a boundary).  */
/* WARNING: errors are maximal when the source is located close to */
/* a grid-point. Best configurations are close to the centre of a  */
/* mesh face, or of a mesh. Errors increase (anisotropically) when */
/* the source gets close to a grid-point. Better use the grid-     */
/* point itself as the source in such case...                      */
{
    int   x,y,z;
    GRID_FLOAT_TYPE distx,disty,distz;

    if(source_at_node){
        t[xs][ys][zs]=0.0;
        return;
    }
    x=(fxs<xs) ? xs-1:xs;
    y=(fys<ys) ? ys-1:ys;
    z=(fzs<zs) ? zs-1:zs;
    /* x,y,z : coordinates of current cell */
    distx=fabs(fxs-x);
    disty=fabs(fys-y);
    distz=fabs(fzs-z);
    /* dist* : distances from source to node minx,miny,minz of current cell */

    init_cell(distx,disty,distz,x,y,z);
    /* this is enough if the source is strictly located */
    /* within the current cell (init: 8 neighbours).    */

    if(fxs==xs){
        if(fys==ys){
            if(x) init_cell(1.,0.,distz,x-1,y,z);
            if(y) init_cell(0.,1.,distz,x,y-1,z);
            if(x && y) init_cell(1.,1.,distz,x-1,y-1,z);
        }/* source located on cell edge parallel to z (18 neighbours) */
        else
        if(fzs==zs){
            if(x) init_cell(1.,disty,0.,x-1,y,z);
            if(z) init_cell(0.,disty,1.,x,y,z-1);
            if(z && x) init_cell(1.,disty,1.,x-1,y,z-1);
        }/* source located on cell edge parallel to y (18 neighbours) */
        else{
            if(x) init_cell(1.,disty,distz,x-1,y,z);
        }/* source located on cell face perpendicular to x (12 neighbours) */
    }
    else
    if(fys==ys){
        if(fzs==zs){
            if(y) init_cell(distx,1.,0.,x,y-1,z);
            if(z) init_cell(distz,0.,1.,x,y,z-1);
            if(y && z) init_cell(distx,1.,1.,x,y-1,z-1);
        }/* source located on cell edge parallel to x (18 neighbours) */
        else {
            if (y) init_cell(distx, 1., distz, x, y - 1, z);
        }/* source located on cell face perpendicular to y (12 neighbours) */
    } else
        if (fzs == zs) {
        if (y) init_cell(distx, disty, 1., x, y, z - 1);
    }/* source located on cell face perpendicular to z (12 neighbours) */

}

/*------------------------------------------------Init_cell()---------------*/

static void
init_cell(GRID_FLOAT_TYPE vx, GRID_FLOAT_TYPE vy, GRID_FLOAT_TYPE vz, int xl, int yl, int zl)

/* compute delays between floating source and nodes of current cell     */
/* xl,yl,zl are current cell coordinates,                               */
/* vx,vy,vz are distances from source to node xl,yl,zl (0<=vx<=1.0,...) */

{
    GRID_FLOAT_TYPE est;
    est = exact_delay(vx, vy, vz, xl, yl, zl);
    if (est < t[xl][yl][zl]) t[xl][yl][zl] = est;
    est = exact_delay(1.0 - vx, vy, vz, xl, yl, zl);
    if (est < t[xl + 1][yl][zl]) t[xl + 1][yl][zl] = est;
    est = exact_delay(vx, 1.0 - vy, vz, xl, yl, zl);
    if (est < t[xl][yl + 1][zl]) t[xl][yl + 1][zl] = est;
    est = exact_delay(vx, vy, 1.0 - vz, xl, yl, zl);
    if (est < t[xl][yl][zl + 1]) t[xl][yl][zl + 1] = est;
    est = exact_delay(1.0 - vx, 1.0 - vy, vz, xl, yl, zl);
    if (est < t[xl + 1][yl + 1][zl]) t[xl + 1][yl + 1][zl] = est;
    est = exact_delay(1.0 - vx, vy, 1.0 - vz, xl, yl, zl);
    if (est < t[xl + 1][yl][zl + 1]) t[xl + 1][yl][zl + 1] = est;
    est = exact_delay(vx, 1.0 - vy, 1.0 - vz, xl, yl, zl);
    if (est < t[xl][yl + 1][zl + 1]) t[xl][yl + 1][zl + 1] = est;
    est = exact_delay(1.0 - vx, 1.0 - vy, 1.0 - vz, xl, yl, zl);
    if (est < t[xl + 1][yl + 1][zl + 1]) t[xl + 1][yl + 1][zl + 1] = est;
}

/*------------------------------------------------Recursive_init()----------*/

static int
recursive_init() {
    int
    signal,
            nx_, ny_, nz_,
            xs_, ys_, zs_,
            X0_, X1_, Y0_, Y1_, Z0_, Z1_,
            n, d,
            i, ii, ihs, i0,
            j, jj, jhs, j0,
            k, kk, khs, k0;
    GRID_FLOAT_TYPE
    *hs_buf_, *t_buf_,
            fxs_, fys_, fzs_,
            HS[N_INIT], T[N_INIT];

    /* increment count of recursivity level */
    init_stage++;
    if (SMALLTALK)
        printf("\nRecursive initialization: level %d", init_stage);

    /* free locally allocated pointers (GRID_FLOAT_TYPE ***) */
    free_ptrs(nx);

    /* save static parameters at this stage */
    nx_ = nx;
    ny_ = ny;
    nz_ = nz;
    hs_buf_ = hs_buf;
    t_buf_ = t_buf;
    xs_ = xs;
    ys_ = ys;
    zs_ = zs;
    fxs_ = fxs;
    fys_ = fys;
    fzs_ = fzs;
    X0_ = X0;
    X1_ = X1;
    Y0_ = Y0;
    Y1_ = Y1;
    Z0_ = Z0;
    Z1_ = Z1;

    /* build the re-discretized local model and the associated source position */
    for (i = 0; i < N_INIT; i++) HS[i] = T[i] = INFINITY;
    nx = ny = nz = N_INIT_X;
    xs = ys = zs = 2 * INIT_MIN + 1;
    i0 = j0 = k0 = 1;
    ihs = xs_ - INIT_MIN - 1;
    if ((d = INIT_MIN - xs_) >= 0) {
        ihs += d + 1;
        d = 1 + 2 * d;
        nx -= d;
        xs -= d;
        i0 = 0;
    }
    if ((d = xs_ + INIT_MIN - nx_ + 1) >= 0) nx -= 1 + 2 * d;
    jhs = ys_ - INIT_MIN - 1;
    if ((d = INIT_MIN - ys_) >= 0) {
        jhs += d + 1;
        d = 1 + 2 * d;
        ny -= d;
        ys -= d;
        j0 = 0;
    }
    if ((d = ys_ + INIT_MIN - ny_ + 1) >= 0) ny -= 1 + 2 * d;
    khs = zs_ - INIT_MIN - 1;
    if ((d = INIT_MIN - zs_) >= 0) {
        khs += d + 1;
        d = 1 + 2 * d;
        nz -= d;
        zs -= d;
        k0 = 0;
    }
    if ((d = zs_ + INIT_MIN - nz_ + 1) >= 0) nz -= 1 + 2 * d;
    for (i = ihs, n = ii = 0; ii < nx; ii++) {
        for (j = jhs, jj = 0; jj < ny; jj++) {
            for (k = khs, kk = 0; kk < nz; kk++, n++) {
                HS[n] = 0.5 * hs_buf_[i * ny_ * nz_ + j * nz_ + k];
                if (kk % 2 != k0) k++;
            }
            if (jj % 2 != j0) j++;
        }
        if (ii % 2 != i0) i++;
    }/* No smoothing is associated with this re-discretization */
    fxs = xs + 2.0 * (fxs_ - xs_);
    fys = ys + 2.0 * (fys_ - ys_);
    fzs = zs + 2.0 * (fzs_ - zs_);

    if (VERBOSE)
        printf("\nRediscretized timefield dimensions: %d %d %d", nx, ny, nz);

    /* recursively compute times on this rediscretized model */
    signal = time_3d(HS, T, nx, ny, nz, fxs, fys, fzs, hs_eps_init, messages);

    /* assign relevant times to parent timefield */
    if (signal == NO_ERROR) {
        for (i = ihs + i0, ii = i0; ii < nx; ii += 2, i++)
            for (j = jhs + j0, jj = j0; jj < ny; jj += 2, j++)
                for (k = khs + k0, kk = k0; kk < nz; kk += 2, k++)
                    t_buf_[i * ny_ * nz_ + j * nz_ + k] = T[ii * ny * nz + jj * nz + kk];
    }

    /* retrieve initial static parameters */
    nx = nx_;
    ny = ny_;
    nz = nz_;
    hs_buf = hs_buf_;
    t_buf = t_buf_;
    xs = xs_;
    ys = ys_;
    zs = zs_;
    fxs = fxs_;
    fys = fys_;
    fzs = fzs_;
    X0 = X0_;
    X1 = X1_;
    Y0 = Y0_;
    Y1 = Y1_;
    Z0 = Z0_;
    Z1 = Z1_;

    /* reallocate pointers (but do not re-initialize!) */
    signal = pre_init();

    /* decrement count of recursivity level */
    init_stage--;

    return (signal);

}

/*------------------------------------------------Propagate_point()---------*/

static int
propagate_point(start)

int start;

{
    int
    msg, test;

    if (start != NO_ERROR) return (start); /* Initialization failed */

    sum_updated = 0;

    /* Make recursive_init silent */
    if (SMALLTALK) printf("\nStarting F.D. computation...");
    msg = messages;
    if (init_stage) messages = 0;

    /* Increment boundaries of timed zone as long as necessary... */
    /* (Outwards propagation is adopted as an initial guess).     */
    do {
        test = 0;

        if (X0 > 0) {
            X0--;
            if (VERBOSE) printf("\nx_side %d->%d: ", X0 + 1, X0);
            x_side(Y0, Y1, Z0, Z1, X0, -1);
            test++;
        }

        if (Y0 > 0) {
            Y0--;
            if (VERBOSE) printf("\ny_side %d->%d: ", Y0 + 1, Y0);
            y_side(X0, X1, Z0, Z1, Y0, -1);
            test++;
        }

        if (Z0 > 0) {
            Z0--;
            if (VERBOSE) printf("\nz_side %d->%d: ", Z0 + 1, Z0);
            z_side(X0, X1, Y0, Y1, Z0, -1);
            test++;
        }

        if (X1 < nmesh_x) {
            X1++;
            if (VERBOSE) printf("\nx_side %d->%d: ", X1 - 1, X1);
            x_side(Y0, Y1, Z0, Z1, X1, 1);
            test++;
        }

        if (Y1 < nmesh_y) {
            Y1++;
            if (VERBOSE) printf("\ny_side %d->%d: ", Y1 - 1, Y1);
            y_side(X0, X1, Z0, Z1, Y1, 1);
            test++;
        }

        if (Z1 < nmesh_z) {
            Z1++;
            if (VERBOSE) printf("\nz_side %d->%d: ", Z1 - 1, Z1);
            z_side(X0, X1, Y0, Y1, Z1, 1);
            test++;
        }

    } while (test);

    messages = msg;

    return (NO_ERROR);

}

/*---------------------------------------------- Free_ptrs()------------------*/

static void
free_ptrs(max_x)

int max_x;

{
    int x, y, z;
    GRID_FLOAT_TYPE *pf;

    /* if relevant, retrieve INFINITY-masked hs values at model boundaries */
    if (init_stage == 0 && hs_keep) {
        pf = hs_keep;
        for (x = 0; x < nx; x++) {
            for (y = 0; y < ny; y++) hs[x][y][nmesh_z] = *pf++;
            for (z = 0; z < nmesh_z; z++) hs[x][nmesh_y][z] = *pf++;
        }
        for (y = 0; y < nmesh_y; y++)
            for (z = 0; z < nmesh_z; z++) hs[nmesh_x][y][z] = *pf++;
        free((char *) hs_keep);
    }

    /* free pointers */
    for (x = 0; x < max_x; x++) {
        free((char *) hs[x]);
        free((char *) t[x]);
    }
    free((char *) hs);
    free((char *) t);
    free((char *) longflags);

}
/****end mail1/3****/
/*--------------------LOCAL 3-D STENCILS (FUNCTIONS AND MACROS)---------------*/
/* Counts of stencils refer to the local inwards isotropic formulation of the */
/* local finite difference computation function (170 stencils for a given     */
/* point, taken as the center of a 2*2*2 cube, with 26 timed neighbours).     */
/* See Podvin and Lecomte, 1991.                                              */
/* In this implementation, this function is never fully computed, because we  */
/* sequentially select only a limited number of relevant directions of propa- */
/* gation, according to the recent history of the signal.                     */
/* As a consequence, only 28 stencils are generally tested at each point.     */
/* The corresponding count is indicated between brackets.                     */
/*----------------------------------------------------------------------------*/
/* Code was restructured in order to minimize the number of actually computed */
/* sqrt()s. Some redundancy was introduced in tests in order to cope properly */
/* with problems arising from numerical errors (e.g. sqrt(x*x) may be lower   */
/* than x, a fact leading to infinite recursions in some awkward cases).      */
/*----------------------------------------------------------------------------*/

/*----------------------------------------- exact_delay() ------------------- */

static GRID_FLOAT_TYPE
exact_delay(GRID_FLOAT_TYPE vx, GRID_FLOAT_TYPE vy, GRID_FLOAT_TYPE vz, int xm, int ym, int zm)

{
    GRID_FLOAT_TYPE estimate;

    if (xm < 0 || xm >= nmesh_x || ym < 0 || ym >= nmesh_y || zm < 0 || zm >= nmesh_z)
        return (INFINITY);
    estimate = (vx * vx + vy * vy + vz * vz) * hs[xm][ym][zm] * hs[xm][ym][zm];
    return (sqrt(estimate));
}

/*----------------------------------------- 1-D transmission : 6 stencils [3] */

/*------------------------------------- (Direct arrival from first neighbour) */
static int
t_1d(int x, int y, int z, GRID_FLOAT_TYPE t0, GRID_FLOAT_TYPE hs0, GRID_FLOAT_TYPE hs1, GRID_FLOAT_TYPE hs2, GRID_FLOAT_TYPE hs3)
{
    GRID_FLOAT_TYPE estimate;
    estimate = t0 + min4(hs0, hs1, hs2, hs3);
    if (estimate < t[x][y][z]) {
        t[x][y][z] = estimate;
        return (1);
    }
    return (0);
}

/*----------------------------------------- 2-D diffraction : 12 stencils [3] */

/*------------------------------------ (Direct arrival from second neighbour) */
static int
diff_2d(int x, int y, int z, GRID_FLOAT_TYPE t0, GRID_FLOAT_TYPE hs0, GRID_FLOAT_TYPE hs1)
{
    GRID_FLOAT_TYPE estimate;
    estimate = t0 + M_SQRT2 * min(hs0, hs1);
    if (estimate < t[x][y][z]) {
        t[x][y][z] = estimate;
        return (1);
    }
    return (0);
}

/*------------------------------------ 3-D point diffraction : 8 stencils [1] */

/*------------------------------------- (Direct arrival from third neighbour) */
static int
point_diff(int x, int y, int z, GRID_FLOAT_TYPE t0, GRID_FLOAT_TYPE hs0)
{
    GRID_FLOAT_TYPE estimate;
    estimate = t0 + hs0*M_SQRT3;
    if (estimate < t[x][y][z]) {
        t[x][y][z] = estimate;
        return (1);
    }
    return (0);
}

/*---------------------------------------- 2-D transmission : 24 stencils [6] */

/*----------------------------------------- (Arrival from coplanar mesh edge) */
static int
t_2d(int x, int y, int z, GRID_FLOAT_TYPE t0, GRID_FLOAT_TYPE t1, GRID_FLOAT_TYPE hs0, GRID_FLOAT_TYPE hs1)
{
    GRID_FLOAT_TYPE estimate, dt, hsm, test2, u2;
    dt = t1 - t0;
    test2 = t[x][y][z] - t1;
    if (dt < 0.0 || test2 < 0.0) return (0);
    test2 *= test2;
    hsm = min(hs0, hs1);
    u2 = hsm * hsm - dt*dt;
    if (dt <= hsm / M_SQRT2 && u2 <= test2) {
        estimate = t1 + sqrt(u2);
        if (estimate < t[x][y][z]) {
            t[x][y][z] = estimate;
            return (1);
        }
    }
    return (0);
}

/*------------------------------------ 3-D edge diffraction : 24 stencils [3] */

/*------------------------------------- (Arrival from non-coplanar mesh edge) */
static int
edge_diff(int x, int y, int z, GRID_FLOAT_TYPE t0, GRID_FLOAT_TYPE t1, GRID_FLOAT_TYPE hs0)
{
    GRID_FLOAT_TYPE estimate, u2, test2, dt;
    dt = t1 - t0;
    test2 = t[x][y][z] - t1;
    if (dt < 0.0 || test2 < 0.0) return (0);
    test2 *= test2;
    u2 = hs0 * hs0 - dt*dt;
    if (dt <= hs0 / M_SQRT3 && 2.0 * u2 <= test2) {
        estimate = t1 + M_SQRT2 * sqrt(u2);
        if (estimate < t[x][y][z]) {
            t[x][y][z] = estimate;
            return (1);
        }
    }
    return (0);
}

/*--------------------------------------- 3-D transmission : 96 stencils [12] */
/*------------------------------------- (Arrival from non-coplanar interface) */
/* 4 stencils per function call or 1+3 using two function calls. */

#define t_3d(x,y,z,a,b,c,d,e)        t_3d_(x,y,z,a,b,c,d,e,0)
#define t_3d_part2(x,y,z,a,b,c,d,e)  t_3d_(x,y,z,a,b,c,d,e,1)

static int
t_3d_(int x, int y, int z, GRID_FLOAT_TYPE t0, GRID_FLOAT_TYPE tl, GRID_FLOAT_TYPE tr, GRID_FLOAT_TYPE td, GRID_FLOAT_TYPE hs0, int redundant)
/* The current point is in diagonal position with respect to t0     */
/* and it is a first neighbour of td. tl,tr are second neighbours.  */
/* One of these estimators is redundant during first step of *_side */
/* functions. See t_3d_part1() which also computes it.              */
/* This function is always called through macros t_3d or t_3d_part2 */
{
    GRID_FLOAT_TYPE test2, r2, s2, t2, u2, dta, dtb, dta2, dtb2, estimate;
    int action;
    action = 0;
    hs0 *= hs0;

    dta = tl - t0;
    dtb = tr - t0;
    dta2 = dta*dta;
    dtb2 = dtb*dtb;
    if (dta >= 0.0 && dtb >= 0.0 && dta2 + dtb2 + dta * dtb >= 0.5 * hs0
            && 2.0 * dta2 + dtb2 <= hs0 && 2.0 * dtb2 + dta2 <= hs0) {
        test2 = t[x][y][z] - tr - tl + t0;
        if (test2 >= 0.0) {
            test2 *= test2;
            r2 = hs0 - dta2 - dtb2;
            if (r2 < test2) {
                estimate = tr + tl - t0 + sqrt(r2);
                if (estimate < t[x][y][z]) {
                    t[x][y][z] = estimate;
                    action++;
                }
            }
        }
    }

    test2 = t[x][y][z] - td;
    if (test2 < 0.0) return (action);
    test2 *= test2;
    s2 = t2 = u2 = INFINITY;

    dtb = td - tl;
    dtb2 = dtb*dtb;
    if (dta >= 0.0 && dtb >= dta && 2.0 * dtb2 + dta2 <= hs0)
        s2 = hs0 - dta2 - dtb2;

    dta = td - tr;
    dta2 = dta*dta;
    if (!redundant && dta >= 0.0 && dtb >= 0.0 && dta2 + dtb2 + dta * dtb <= 0.5 * hs0)
        t2 = hs0 - dta2 - dtb2;

    dtb = tr - t0;
    dtb2 = dtb*dtb;
    if (dtb >= 0.0 && dta >= dtb && 2.0 * dta2 + dtb2 <= hs0)
        u2 = hs0 - dta2 - dtb2;

    u2 = min3(s2, t2, u2);
    if (u2 < test2) {
        estimate = td + sqrt(u2);
        if (estimate < t[x][y][z]) {
            t[x][y][z] = estimate;
            action++;
        }
    }
    return (action);
}

/* 3-D transmission: partial stencil, introduced because initial scheme */
/* failed to fullfill the exhaustivity condition requested by Fermat's  */

/* principle. (See "a-causal" step in *_side() functions; 18/07/91)     */
static int
t_3d_part1(int x, int y, int z, GRID_FLOAT_TYPE t0, GRID_FLOAT_TYPE tl, GRID_FLOAT_TYPE tr, GRID_FLOAT_TYPE hs0)
/* The current point is a first neighbour of t0; tl,tr are two other */
/* first neighbours of t0. Transmission through 0-l-r is tested.     */
{
    GRID_FLOAT_TYPE dtl, dtr, s2, u2, estimate, test2;
    dtl = t0 - tl;
    dtr = t0 - tr;
    test2 = t[x][y][z] - t0;
    if (test2 < 0.0 || dtl < 0.0 || dtr < 0.0) return (0);
    test2 *= test2;
    hs0 *= hs0;
    s2 = dtl * dtl + dtr*dtr;
    if (s2 + dtl * dtr > 0.5 * hs0) return (0);
    /* illumination condition */
    u2 = hs0 - s2;
    if (u2 < test2) {
        estimate = t0 + sqrt(u2);
        if (estimate < t[x][y][z]) {
            t[x][y][z] = estimate;
            return (1);
        }
    }
    return (0);
}

/*----------------------------------------------X_SIDE()--------------------*/

static int
x_side(y_begin, y_end, z_begin, z_end, x, future)

int y_begin, y_end, z_begin, z_end, x, future;

/* Propagates computations from side x-future to current side x */
/* between *_begin and *_end coordinates. Returns a nonzero     */
/* integer if something actually happened (a time was lowered). */
/* Extensions _bb, _fb etc... define simple orientation rules:  */
/* _bf means backwards along y axis and forwards along z axis.  */
/* So-called "longitudinal" headwaves refer to first arrivals   */
/* due to a headwave propagating along the current side.        */

{
    int
    updated, /* counts adopted FD stencils */
            longhead, /* counts "longitudinal" headwaves */
            x0, /* past side coordinate */
            x_s, /* current meshes coordinate */
            y, z, /* current point coordinate */
            sign_ff, sign_bf, sign_bb, sign_fb, /* sign flags for time differences */
            past, /* opposite to future ! */
            test;
    GRID_FLOAT_TYPE
    hs_ff, hs_bf, hs_bb, hs_fb; /* local slownesses */

    if (reverse_order == 0)
        current_side_limit = x + future;
    updated = 0;
    x0 = x - future;
    if (future == 1) x_s = x0;
    else x_s = x;

    flag_fb = flag_bf = flag_ff = flag_bb = 0;
    y_start_ff = y_start_fb = y_end;
    y_start_bf = y_start_bb = y_begin;
    z_start_ff = z_start_bf = z_end;
    z_start_fb = z_start_bb = z_begin;

    /* First,  Compute all stencils using only nodes of side x0.   */
    /* As only times on side x will be changed, these stencils     */
    /* are computed initially, once for all, in any order (no      */
    /* causality problem involved).                                */
    /* During this first pass, future directions of propagation    */
    /* are diagnosed, according to the time pattern on side x0.    */
    /* Borders of zones to be scanned are computed.                */
    /* This part may be seen as the explicit part of the FD scheme.*/

    for (y = y_begin; y <= y_end; y++) {
        for (z = z_begin; z <= z_end; z++) {

            hs_ff = hs[x_s][y][z];
            if (y > 0) hs_bf = hs[x_s][y - 1][z];
            else hs_bf = INFINITY;
            if (z > 0 && y > 0) hs_bb = hs[x_s][y - 1][z - 1];
            else hs_bb = INFINITY;
            if (z > 0) hs_fb = hs[x_s][y][z - 1];
            else hs_fb = INFINITY;
            sign_fb = sign_bf = sign_ff = sign_bb = 0;

            /* illuminate first neighbours */
            /* 1 1D transmission and 4 partial 3D transmission */
            updated += t_1d(x, y, z, t[x0][y][z], hs_ff, hs_bf, hs_bb, hs_fb);
            if (y < y_end && z < z_end)
                updated += t_3d_part1(x, y, z,
                    t[x0][y][z], t[x0][y + 1][z], t[x0][y][z + 1], hs_ff);
            if (y > y_begin && z < z_end)
                updated += t_3d_part1(x, y, z,
                    t[x0][y][z], t[x0][y - 1][z], t[x0][y][z + 1], hs_bf);
            if (y > y_begin && z > z_begin)
                updated += t_3d_part1(x, y, z,
                    t[x0][y][z], t[x0][y - 1][z], t[x0][y][z - 1], hs_bb);
            if (y < y_end && z > z_begin)
                updated += t_3d_part1(x, y, z,
                    t[x0][y][z], t[x0][y + 1][z], t[x0][y][z - 1], hs_fb);

            /* illuminate second neighbours (if necessary)    */
            /* 4 2D diffraction and 4 2D transmission */
            if (y < y_end && t[x0][y][z] <= t[x0][y + 1][z]) {
                sign_fb++;
                sign_ff++;
                if (y < y_start_ff) y_start_ff = y;
                if (y < y_start_fb) y_start_fb = y;
                updated += diff_2d(x, y + 1, z, t[x0][y][z], hs_ff, hs_fb);
                updated += t_2d(x, y + 1, z, t[x0][y][z], t[x0][y + 1][z], hs_ff, hs_fb);
            }
            if (y > y_begin && t[x0][y][z] <= t[x0][y - 1][z]) {
                sign_bb++;
                sign_bf++;
                if (y > y_start_bf) y_start_bf = y;
                if (y > y_start_bb) y_start_bb = y;
                updated += diff_2d(x, y - 1, z, t[x0][y][z], hs_bf, hs_bb);
                updated += t_2d(x, y - 1, z, t[x0][y][z], t[x0][y - 1][z], hs_bf, hs_bb);
            }
            if (z < z_end && t[x0][y][z] <= t[x0][y][z + 1]) {
                sign_bf++;
                sign_ff++;
                if (z < z_start_ff) z_start_ff = z;
                if (z < z_start_bf) z_start_bf = z;
                updated += diff_2d(x, y, z + 1, t[x0][y][z], hs_ff, hs_bf);
                updated += t_2d(x, y, z + 1, t[x0][y][z], t[x0][y][z + 1], hs_ff, hs_bf);
            }
            if (z > z_begin && t[x0][y][z] <= t[x0][y][z - 1]) {
                sign_bb++;
                sign_fb++;
                if (z > z_start_fb) z_start_fb = z;
                if (z > z_start_bb) z_start_bb = z;
                updated += diff_2d(x, y, z - 1, t[x0][y][z], hs_bb, hs_fb);
                updated += t_2d(x, y, z - 1, t[x0][y][z], t[x0][y][z - 1], hs_bb, hs_fb);
            }

            /* illuminate third neighbours (if necessary) */
            /* 4 3D point diffraction, 8 3D edge diffraction and 12 3D transmission */
            if (sign_ff == 2) {
                flag_ff = 1;
                updated += point_diff(x, y + 1, z + 1, t[x0][y][z], hs_ff);
                updated += edge_diff(x, y + 1, z + 1, t[x0][y][z], t[x0][y + 1][z], hs_ff);
                updated += edge_diff(x, y + 1, z + 1, t[x0][y][z], t[x0][y][z + 1], hs_ff);
                updated += t_3d_part2(x, y + 1, z + 1, t[x0][y][z],
                        t[x0][y + 1][z], t[x0][y][z + 1], t[x0][y + 1][z + 1], hs_ff);
            }
            if (sign_bf == 2) {
                flag_bf = 1;
                updated += point_diff(x, y - 1, z + 1, t[x0][y][z], hs_bf);
                updated += edge_diff(x, y - 1, z + 1, t[x0][y][z], t[x0][y - 1][z], hs_bf);
                updated += edge_diff(x, y - 1, z + 1, t[x0][y][z], t[x0][y][z + 1], hs_bf);
                updated += t_3d_part2(x, y - 1, z + 1, t[x0][y][z],
                        t[x0][y - 1][z], t[x0][y][z + 1], t[x0][y - 1][z + 1], hs_bf);
            }
            if (sign_bb == 2) {
                flag_bb = 1;
                updated += point_diff(x, y - 1, z - 1, t[x0][y][z], hs_bb);
                updated += edge_diff(x, y - 1, z - 1, t[x0][y][z], t[x0][y - 1][z], hs_bb);
                updated += edge_diff(x, y - 1, z - 1, t[x0][y][z], t[x0][y][z - 1], hs_bb);
                updated += t_3d_part2(x, y - 1, z - 1, t[x0][y][z],
                        t[x0][y - 1][z], t[x0][y][z - 1], t[x0][y - 1][z - 1], hs_bb);
            }
            if (sign_fb == 2) {
                flag_fb = 1;
                updated += point_diff(x, y + 1, z - 1, t[x0][y][z], hs_fb);
                updated += edge_diff(x, y + 1, z - 1, t[x0][y][z], t[x0][y + 1][z], hs_fb);
                updated += edge_diff(x, y + 1, z - 1, t[x0][y][z], t[x0][y][z - 1], hs_fb);
                updated += t_3d_part2(x, y + 1, z - 1, t[x0][y][z],
                        t[x0][y + 1][z], t[x0][y][z - 1], t[x0][y + 1][z - 1], hs_fb);
            }
        }
    }

    /* Now, all remaining stencils depend on nodes located on the current  */
    /* side. They must be propagated causally, and occurrences of critical */
    /* conditions must be diagnosed, in order to propagate associated head */
    /* waves exhaustively. Four independent scanning directions are succes-*/
    /* sively explored, and headwave flags are used to generate the supple-*/
    /* mentary scans requested by exhaustivity ("Reverse" Propagations).   */
    /* flag_* flag  is non-zero while the corresponding direction remains  */
    /* to be examined (or reexamined). Its examination may detect critical */
    /* conditions which in their turn will make another scan necessary.    */
    /* Initialization of this process was achieved during first step.      */
    /* This second step may be seen as the implicit part of the FD scheme. */

    /* initialize local headwave flags */
    for (y = 0; y < ny * nz; y++) longflags[y] = 0;

    /* enforce all scans if current side has a null surface */
    /* (This may only be encountered near the source point) */
    if (y_begin == y_end || z_begin == z_end) {
        flag_ff = flag_fb = flag_bf = flag_bb = 1;
        y_start_ff = y_start_fb = y_begin;
        y_start_bf = y_start_bb = y_end;
        z_start_ff = z_start_bf = z_begin;
        z_start_fb = z_start_bb = z_end;
    }

    /* Reexamine each direction, while necessary */
    do {
        test = 0;
        if (flag_ff) {
            test++;
            if (VERBOSE) printf("ff ");
            updated += scan_x_ff(y_start_ff, y_end, z_start_ff, z_end, x0, x, x_s);
        }
        if (flag_fb) {
            test++;
            if (VERBOSE) printf("fb ");
            updated += scan_x_fb(y_start_fb, y_end, z_begin, z_start_fb, x0, x, x_s);
        }
        if (flag_bb) {
            test++;
            if (VERBOSE) printf("bb ");
            updated += scan_x_bb(y_begin, y_start_bb, z_begin, z_start_bb, x0, x, x_s);
        }
        if (flag_bf) {
            test++;
            if (VERBOSE) printf("bf ");
            updated += scan_x_bf(y_begin, y_start_bf, z_start_bf, z_end, x0, x, x_s);
        }
    } while (test);

    sum_updated += updated;

    /* At this stage, all points of the current side have been timed.     */
    /* Now, Reverse propagation must be invoked if a headwave propagating */
    /* along the current side was generated, because a new branch of the  */
    /* timefield (conical wave) propagates towards the already timed box. */

    for (y = longhead = 0; y < ny * nz; y++) longhead += longflags[y];

    if (longhead) {

        reverse_order++;

        if (VERBOSE) printf("\nReverse#%d from x_side %d", reverse_order, x);
        past = -future;
        for (x = x0; x != current_side_limit; x += past) {
            if (x < 0 || x >= nx) break;
            if (VERBOSE) printf("\nupdate side x=%d: ", x);
            if (x_side(y_begin, y_end, z_begin, z_end, x, past) == 0) break;
            if (VERBOSE) printf("x=%d <R#%d>updated.", x, reverse_order);
        }
        if (VERBOSE) printf("\nEnd Reverse#%d\n", reverse_order);

        reverse_order--;

    }

    return (updated);

}

/*--------------------------------------X_SIDE() : SCAN_X_EE()--------------*/

static int
scan_x_ff(y_start, y_end, z_start, z_end, x0, x, x_s)
int y_start, y_end, z_start, z_end, x0, x, x_s;

/* scan x_side by increasing y and z ("ff"=forwards, forwards)      */
/* propagating causal stencils with provisional a-priori that some  */
/* significant wavefronts propagate in this direction in the zone   */
/* defined by y_start,y_end,z_start,z_end.                          */
/* Critical conditions on any interface are detected so that other  */
/* relevant directions of propagation due to headwave generation    */
/* may be exhaustively taken into account at a later stage.         */

{
    int
    updated = 0,
            alert0, alert1,
            y, z,
            x_sf;
    GRID_FLOAT_TYPE
    hs_bf, hs_bb, hs_fb,
            hs_ube, hs_ubb, hs_ueb;

    x_sf = x_s + x - x0;

    /* We first propagate headwaves along the two borders of the current zone */
    /* These headwaves are usually relevant only when a local minimum valley  */
    /* is present along these borders on the preceding x_side (x=x0).         */
    /* This is analogous with timing local minima in the 2-D implementation.  */

    /* interface waves along y_side: 1 1D transmission and 1 2D transmission */
    hs_bb = hs_ubb = hs_ube = INFINITY;
    for (y = y_start, z = z_start; y < y_end; y++) {
        hs_bf = hs[x_s][y][z];
        if (z) hs_bb = hs[x_s][y][z - 1];
        if (x_sf >= 0 && x_sf < nmesh_x) {
            hs_ube = hs[x_sf][y][z];
            if (z) hs_ubb = hs[x_sf][y][z - 1];
        }
        alert1 = t_1d(x, y + 1, z, t[x][y][z], hs_bb, hs_bf, hs_ubb, hs_ube);
        alert0 = t_2d(x, y + 1, z, t[x0][y][z], t[x][y][z], hs_bb, hs_bf);
        updated += alert0 + alert1;
        if (alert1) longflags[y * nz + nz + z] = 1;
        if (alert0) longflags[y * nz + nz + z] = 0;
    }

    /* interface waves along z_side: 1 1D transmission and 1 2D transmission */
    hs_bb = hs_ubb = hs_ueb = INFINITY;
    for (y = y_start, z = z_start; z < z_end; z++) {
        hs_fb = hs[x_s][y][z];
        if (y) hs_bb = hs[x_s][y - 1][z];
        if (x_sf >= 0 && x_sf < nmesh_x) {
            hs_ueb = hs[x_sf][y][z];
            if (y) hs_ubb = hs[x_sf][y - 1][z];
        }
        alert1 = t_1d(x, y, z + 1, t[x][y][z], hs_bb, hs_fb, hs_ubb, hs_ueb);
        alert0 = t_2d(x, y, z + 1, t[x0][y][z], t[x][y][z], hs_bb, hs_fb);
        updated += alert0 + alert1;
        if (alert1) longflags[y * nz + z + 1] = 1;
        if (alert0) longflags[y * nz + z + 1] = 0;
    }

    /* We now propagate bulk and head waves into the region of interest. */

    for (y = y_start; y < y_end; y++) {
        for (z = z_start; z < z_end; z++) {

            hs_bb = hs[x_s][y][z];
            hs_bf = hs[x_s][y][z + 1];
            hs_fb = hs[x_s][y + 1][z];
            if (x_sf >= 0 && x_sf < nmesh_x) {
                hs_ubb = hs[x_sf][y][z];
                hs_ube = hs[x_sf][y][z + 1];
                hs_ueb = hs[x_sf][y + 1][z];
            } else hs_ubb = hs_ube = hs_ueb = INFINITY;

            /* bulk waves: 1 3D edge diffraction and 2 (*4) 3D transmission */
            alert0 = edge_diff(x, y + 1, z + 1, t[x0][y][z], t[x][y][z], hs_bb)
                    + t_3d(x, y + 1, z + 1, t[x0][y][z],
                    t[x0][y + 1][z], t[x][y][z], t[x][y + 1][z], hs_bb)
                    + t_3d(x, y + 1, z + 1, t[x0][y][z],
                    t[x0][y][z + 1], t[x][y][z], t[x][y][z + 1], hs_bb);
            if (alert0) {
                updated += alert0;
                longflags[y * nz + nz + z + 1] = 0;
            }

            /* interface waves along y_side: 1 1D transmission and 1 2D transmission */
            alert1 = t_1d(x, y + 1, z + 1, t[x][y][z + 1], hs_bb, hs_bf, hs_ubb, hs_ube);
            alert0 = t_2d(x, y + 1, z + 1, t[x0][y][z + 1], t[x][y][z + 1], hs_bb, hs_bf);
            if (alert0 + alert1) {
                updated += alert0 + alert1;
                flag_fb++; /* this scan must be (re-)examined */
                if (y_start_fb > y) y_start_fb = y;
                if (z_start_fb < z + 1) z_start_fb = z + 1;
                if (alert1) longflags[y * nz + nz + z + 1] = 1;
                else longflags[y * nz + nz + z + 1] = 0;
            }

            /* interface waves along z_side: 1 1D transmission and 1 2D transmission */
            alert1 = t_1d(x, y + 1, z + 1, t[x][y + 1][z], hs_bb, hs_fb, hs_ubb, hs_ueb);
            alert0 = t_2d(x, y + 1, z + 1, t[x0][y + 1][z], t[x][y + 1][z], hs_bb, hs_fb);
            if (alert0 + alert1) {
                updated += alert0 + alert1;
                flag_bf++; /* this scan must be (re-)examined */
                if (y_start_bf < y + 1) y_start_bf = y + 1;
                if (z_start_bf > z) z_start_bf = z;
                if (alert1) longflags[y * nz + nz + z + 1] = 1;
                else longflags[y * nz + nz + z + 1] = 0;
            }

            /* interface waves along x_side : 2 2D transmission and 1 2D diffraction */
            alert1 = diff_2d(x, y + 1, z + 1, t[x][y][z], hs_bb, hs_ubb)
                    + t_2d(x, y + 1, z + 1, t[x][y][z], t[x][y + 1][z], hs_bb, hs_ubb)
                    + t_2d(x, y + 1, z + 1, t[x][y][z], t[x][y][z + 1], hs_bb, hs_ubb);
            if (alert1) {
                updated += alert1;
                longflags[y * nz + nz + z + 1] = 1;
            }
        }
    }

    flag_ff = 0;
    y_start_ff = y_end;
    z_start_ff = z_end;
    /* this direction has been examined: unset corresponding flag */

    return (updated);
}

/*--------------------------------------X_SIDE() : SCAN_X_BE()--------------*/

static int
scan_x_bf(y_begin, y_start, z_start, z_end, x0, x, x_s)
int y_start, y_begin, z_start, z_end, x0, x, x_s;

{
    int
    updated = 0,
            alert0, alert1,
            y, z,
            x_sf;
    GRID_FLOAT_TYPE
    hs_ff, hs_bb, hs_fb,
            hs_uee, hs_ubb, hs_ueb;

    x_sf = x_s + x - x0;

    hs_fb = hs_uee = hs_ueb = INFINITY;
    for (y = y_start, z = z_start; y > y_begin; y--) {
        hs_ff = hs[x_s][y - 1][z];
        if (z) hs_fb = hs[x_s][y - 1][z - 1];
        if (x_sf >= 0 && x_sf < nmesh_x) {
            hs_uee = hs[x_sf][y - 1][z];
            if (z) hs_ueb = hs[x_sf][y - 1][z - 1];
        }
        alert1 = t_1d(x, y - 1, z, t[x][y][z], hs_fb, hs_ff, hs_ueb, hs_uee);
        alert0 = t_2d(x, y - 1, z, t[x0][y][z], t[x][y][z], hs_fb, hs_ff);
        updated += alert0 + alert1;
        if (alert1) longflags[y * nz - nz + z] = 1;
        if (alert0) longflags[y * nz - nz + z] = 0;
    }

    hs_bb = hs_ubb = hs_ueb = INFINITY;
    for (y = y_start, z = z_start; z < z_end; z++) {
        if (y) hs_bb = hs[x_s][y - 1][z];
        hs_fb = hs[x_s][y][z];
        if (x_sf >= 0 && x_sf < nmesh_x) {
            if (y) hs_ubb = hs[x_sf][y - 1][z];
            hs_ueb = hs[x_sf][y][z];
        }
        alert1 = t_1d(x, y, z + 1, t[x][y][z], hs_bb, hs_fb, hs_ubb, hs_ueb);
        alert0 = t_2d(x, y, z + 1, t[x0][y][z], t[x][y][z], hs_bb, hs_fb);
        updated += alert0 + alert1;
        if (alert1) longflags[y * nz + z + 1] = 1;
        if (alert0) longflags[y * nz + z + 1] = 0;
    }

    for (y = y_start; y > y_begin; y--) {
        for (z = z_start; z < z_end; z++) {

            hs_ff = hs[x_s][y - 1][z + 1];
            if (y > 1) hs_bb = hs[x_s][y - 2][z];
            else hs_bb = INFINITY;
            hs_fb = hs[x_s][y - 1][z];
            if (x_sf >= 0 && x_sf < nmesh_x) {
                hs_uee = hs[x_sf][y - 1][z + 1];
                if (y > 1) hs_ubb = hs[x_sf][y - 2][z];
                else hs_ubb = INFINITY;
                hs_ueb = hs[x_sf][y - 1][z];
            } else hs_ubb = hs_uee = hs_ueb = INFINITY;

            /* bulk waves: 1 3D edge diffraction and 2 (*4) 3D transmission */
            alert0 = edge_diff(x, y - 1, z + 1, t[x0][y][z], t[x][y][z], hs_fb)
                    + t_3d(x, y - 1, z + 1, t[x0][y][z],
                    t[x0][y - 1][z], t[x][y][z], t[x][y - 1][z], hs_fb)
                    + t_3d(x, y - 1, z + 1, t[x0][y][z],
                    t[x0][y][z + 1], t[x][y][z], t[x][y][z + 1], hs_fb);
            if (alert0) {
                updated += alert0;
                longflags[y * nz - nz + z + 1] = 0;
            }

            /* interface waves along y_side: 1 1D transmission and 1 2D transmission */
            alert1 = t_1d(x, y - 1, z + 1, t[x][y][z + 1], hs_ff, hs_fb, hs_uee, hs_ueb);
            alert0 = t_2d(x, y - 1, z + 1, t[x0][y][z + 1], t[x][y][z + 1], hs_ff, hs_fb);
            if (alert0 + alert1) {
                // 20100204 AJL Satriano Bug Fix.
                updated += alert0 + alert1;
                //updated+alert0+alert1;
                flag_bb++; /* this scan must be (re-)examined */
                if (y_start_bb < y) y_start_bb = y;
                if (z_start_bb < z + 1) z_start_bb = z + 1;
                if (alert1) longflags[y * nz - nz + z + 1] = 1;
                else longflags[y * nz - nz + z + 1] = 0;
            }

            /* interface waves along z_side: 1 1D transmission and 1 2D transmission */
            alert1 = t_1d(x, y - 1, z + 1, t[x][y - 1][z], hs_bb, hs_fb, hs_ubb, hs_ueb);
            alert0 = t_2d(x, y - 1, z + 1, t[x0][y - 1][z], t[x][y - 1][z], hs_bb, hs_fb);
            if (alert0 + alert1) {
                // 20100204 AJL Satriano Bug Fix.
                updated += alert0 + alert1;
                //updated+alert0+alert1;
                flag_ff++; /* this scan must be (re-)examined */
                if (y_start_ff > y - 1) y_start_ff = y - 1;
                if (z_start_ff > z) z_start_ff = z;
                if (alert1) longflags[y * nz - nz + z + 1] = 1;
                else longflags[y * nz - nz + z + 1] = 0;
            }

            /* interface waves along x_side : 2 2D transmission and 1 2D diffraction */
            alert1 = diff_2d(x, y - 1, z + 1, t[x][y][z], hs_fb, hs_ueb)
                    + t_2d(x, y - 1, z + 1, t[x][y][z], t[x][y - 1][z], hs_fb, hs_ueb)
                    + t_2d(x, y - 1, z + 1, t[x][y][z], t[x][y][z + 1], hs_fb, hs_ueb);
            if (alert1) {
                updated += alert1;
                longflags[y * nz - nz + z + 1] = 1;
            }
        }
    }

    flag_bf = 0;
    y_start_bf = y_begin;
    z_start_bf = z_end;
    /* this direction has been examined: unset corresponding flag */

    return (updated);
}

/*--------------------------------------X_SIDE() : SCAN_X_BB()--------------*/

static int
scan_x_bb(y_begin, y_start, z_begin, z_start, x0, x, x_s)
int y_start, y_begin, z_begin, z_start, x0, x, x_s;

{
    int
    updated = 0,
            alert0, alert1,
            y, z,
            x_sf;
    GRID_FLOAT_TYPE
    hs_ff, hs_bf, hs_fb,
            hs_uee, hs_ube, hs_ueb;

    x_sf = x_s + x - x0;

    hs_ff = hs_uee = hs_ueb = INFINITY;
    for (y = y_start, z = z_start; y > y_begin; y--) {
        if (z) hs_ff = hs[x_s][y - 1][z - 1];
        hs_fb = hs[x_s][y - 1][z];
        if (x_sf >= 0 && x_sf < nmesh_x) {
            if (z) hs_uee = hs[x_sf][y - 1][z - 1];
            hs_ueb = hs[x_sf][y - 1][z];
        }
        alert1 = t_1d(x, y - 1, z, t[x][y][z], hs_fb, hs_ff, hs_ueb, hs_uee);
        alert0 = t_2d(x, y - 1, z, t[x0][y][z], t[x][y][z], hs_fb, hs_ff);
        updated += alert0 + alert1;
        if (alert1) longflags[y * nz - nz + z] = 1;
        if (alert0) longflags[y * nz - nz + z] = 0;
    }

    hs_ff = hs_uee = hs_ube = INFINITY;
    for (y = y_start, z = z_start; z > z_begin; z--) {
        if (y) hs_ff = hs[x_s][y - 1][z - 1];
        hs_bf = hs[x_s][y][z - 1];
        if (x_sf >= 0 && x_sf < nmesh_x) {
            if (y) hs_uee = hs[x_sf][y - 1][z - 1];
            hs_ube = hs[x_sf][y][z - 1];
        }
        alert1 = t_1d(x, y, z - 1, t[x][y][z], hs_bf, hs_ff, hs_ube, hs_uee);
        alert0 = t_2d(x, y, z - 1, t[x0][y][z], t[x][y][z], hs_bf, hs_ff);
        updated += alert0 + alert1;
        if (alert1) longflags[y * nz + z - 1] = 1;
        if (alert0) longflags[y * nz + z - 1] = 0;
    }

    for (y = y_start; y > y_begin; y--) {
        for (z = z_start; z > z_begin; z--) {

            hs_ff = hs[x_s][y - 1][z - 1];
            if (y > 1) hs_bf = hs[x_s][y - 2][z - 1];
            else hs_bf = INFINITY;
            if (z > 1) hs_fb = hs[x_s][y - 1][z - 2];
            else hs_fb = INFINITY;
            if (x_sf >= 0 && x_sf < nmesh_x) {
                hs_uee = hs[x_sf][y - 1][z - 1];
                if (y > 1) hs_ube = hs[x_sf][y - 2][z - 1];
                else hs_ube = INFINITY;
                if (z > 1) hs_ueb = hs[x_sf][y - 1][z - 2];
                else hs_ueb = INFINITY;
            } else hs_uee = hs_ube = hs_ueb = INFINITY;

            /* bulk waves: 1 3D edge diffraction and 2 (*4) 3D transmission */
            alert0 = edge_diff(x, y - 1, z - 1, t[x0][y][z], t[x][y][z], hs_ff)
                    + t_3d(x, y - 1, z - 1, t[x0][y][z],
                    t[x0][y - 1][z], t[x][y][z], t[x][y - 1][z], hs_ff)
                    + t_3d(x, y - 1, z - 1, t[x0][y][z],
                    t[x0][y][z - 1], t[x][y][z], t[x][y][z - 1], hs_ff);
            if (alert0) {
                updated += alert0;
                longflags[y * nz - nz + z - 1] = 0;
            }

            /* interface waves along y_side: 1 1D transmission and 1 2D transmission */
            alert1 = t_1d(x, y - 1, z - 1, t[x][y][z - 1], hs_ff, hs_fb, hs_uee, hs_ueb);
            alert0 = t_2d(x, y - 1, z - 1, t[x0][y][z - 1], t[x][y][z - 1], hs_ff, hs_fb);
            if (alert0 + alert1) {
                updated += alert0 + alert1;
                flag_bf++; /* this scan must be (re-)examined */
                if (y_start_bf < y) y_start_bf = y;
                if (z_start_bf > z - 1) z_start_bf = z - 1;
                if (alert1) longflags[y * nz - nz + z - 1] = 1;
                else longflags[y * nz - nz + z - 1] = 0;
            }

            /* interface waves along z_side: 1 1D transmission and 1 2D transmission */
            alert1 = t_1d(x, y - 1, z - 1, t[x][y - 1][z], hs_ff, hs_bf, hs_uee, hs_ube);
            alert0 = t_2d(x, y - 1, z - 1, t[x0][y - 1][z], t[x][y - 1][z], hs_ff, hs_bf);
            if (alert0 + alert1) {
                updated += alert0 + alert1;
                flag_fb++; /* this scan must be (re-)examined */
                if (y_start_fb > y - 1) y_start_fb = y - 1;
                if (z_start_fb < z) z_start_fb = z;
                if (alert1) longflags[y * nz - nz + z - 1] = 1;
                else longflags[y * nz - nz + z - 1] = 0;
            }

            /* interface waves along x_side : 2 2D transmission and 1 2D diffraction */
            alert1 = diff_2d(x, y - 1, z - 1, t[x][y][z], hs_ff, hs_uee)
                    + t_2d(x, y - 1, z - 1, t[x][y][z], t[x][y - 1][z], hs_ff, hs_uee)
                    + t_2d(x, y - 1, z - 1, t[x][y][z], t[x][y][z - 1], hs_ff, hs_uee);
            if (alert1) {
                updated += alert1;
                longflags[y * nz - nz + z - 1] = 1;
            }
        }
    }

    flag_bb = 0;
    y_start_bb = y_begin;
    z_start_bb = z_begin;
    /* this direction has been examined: unset corresponding flag */

    return (updated);
}

/*--------------------------------------X_SIDE() : SCAN_X_EB()--------------*/

static int
scan_x_fb(y_start, y_end, z_begin, z_start, x0, x, x_s)
int y_start, y_end, z_begin, z_start, x0, x, x_s;

{
    int
    updated = 0,
            alert0, alert1,
            y, z,
            x_sf;
    GRID_FLOAT_TYPE
    hs_ff, hs_bb, hs_bf,
            hs_uee, hs_ubb, hs_ube;

    x_sf = x_s + x - x0;

    hs_bf = hs_ubb = hs_ube = INFINITY;
    for (y = y_start, z = z_start; y < y_end; y++) {
        if (z) hs_bf = hs[x_s][y][z - 1];
        hs_bb = hs[x_s][y][z];
        if (x_sf >= 0 && x_sf < nmesh_x) {
            if (z) hs_ube = hs[x_sf][y][z - 1];
            hs_ubb = hs[x_sf][y][z];
        }
        alert1 = t_1d(x, y + 1, z, t[x][y][z], hs_bf, hs_bb, hs_ube, hs_ubb);
        alert0 = t_2d(x, y + 1, z, t[x0][y][z], t[x][y][z], hs_bf, hs_bb);
        updated += alert0 + alert1;
        if (alert1) longflags[y * nz + nz + z] = 1;
        if (alert0) longflags[y * nz + nz + z] = 0;
    }

    hs_ff = hs_uee = hs_ube = INFINITY;
    for (y = y_start, z = z_start; z > z_begin; z--) {
        if (y) hs_ff = hs[x_s][y - 1][z - 1];
        hs_bf = hs[x_s][y][z - 1];
        if (x_sf >= 0 && x_sf < nmesh_x) {
            if (y) hs_uee = hs[x_sf][y - 1][z - 1];
            hs_ube = hs[x_sf][y][z - 1];
        }
        alert1 = t_1d(x, y, z - 1, t[x][y][z], hs_bf, hs_ff, hs_ube, hs_uee);
        alert0 = t_2d(x, y, z - 1, t[x0][y][z], t[x][y][z], hs_bf, hs_ff);
        updated += alert0 + alert1;
        if (alert1) longflags[y * nz + z - 1] = 1;
        if (alert0) longflags[y * nz + z - 1] = 0;
    }

    for (y = y_start; y < y_end; y++) {
        for (z = z_start; z > z_begin; z--) {

            hs_bf = hs[x_s][y][z - 1];
            if (z > 1) hs_bb = hs[x_s][y][z - 2];
            else hs_bb = INFINITY;
            hs_ff = hs[x_s][y + 1][z - 1];
            if (x_sf >= 0 && x_sf < nmesh_x) {
                hs_ube = hs[x_sf][y][z - 1];
                if (z > 1) hs_ubb = hs[x_sf][y][z - 2];
                else hs_ubb = INFINITY;
                hs_uee = hs[x_sf][y + 1][z - 1];
            } else hs_ubb = hs_ube = hs_uee = INFINITY;

            /* bulk waves: 1 3D edge diffraction and 2 (*4) 3D transmission */
            alert0 = edge_diff(x, y + 1, z - 1, t[x0][y][z], t[x][y][z], hs_bf)
                    + t_3d(x, y + 1, z - 1, t[x0][y][z],
                    t[x0][y + 1][z], t[x][y][z], t[x][y + 1][z], hs_bf)
                    + t_3d(x, y + 1, z - 1, t[x0][y][z],
                    t[x0][y][z - 1], t[x][y][z], t[x][y][z - 1], hs_bf);
            if (alert0) {
                updated += alert0;
                longflags[y * nz + nz + z - 1] = 0;
            }

            /* interface waves along y_side: 1 1D transmission and 1 2D transmission */
            alert1 = t_1d(x, y + 1, z - 1, t[x][y][z - 1], hs_bb, hs_bf, hs_ubb, hs_ube);
            alert0 = t_2d(x, y + 1, z - 1, t[x0][y][z - 1], t[x][y][z - 1], hs_bb, hs_bf);
            if (alert0 + alert1) {
                updated += alert0 + alert1;
                flag_ff++; /* this scan must be (re-)examined */
                if (y_start_ff > y) y_start_ff = y;
                if (z_start_ff > z - 1) z_start_ff = z - 1;
                if (alert1) longflags[y * nz + nz + z - 1] = 1;
                else longflags[y * nz + nz + z - 1] = 0;
            }

            /* interface waves along z_side: 1 1D transmission and 1 2D transmission */
            alert1 = t_1d(x, y + 1, z - 1, t[x][y + 1][z], hs_ff, hs_bf, hs_uee, hs_ube);
            alert0 = t_2d(x, y + 1, z - 1, t[x0][y + 1][z], t[x][y + 1][z], hs_ff, hs_bf);
            if (alert0 + alert1) {
                updated += alert0 + alert1;
                flag_bb++; /* this scan must be (re-)examined */
                if (y_start_bb < y + 1) y_start_bb = y + 1;
                if (z_start_bb < z) z_start_bb = z;
                if (alert1) longflags[y * nz + nz + z - 1] = 1;
                else longflags[y * nz + nz + z - 1] = 0;
            }

            /* interface waves along x_side : 2 2D transmission and 1 2D diffraction */
            alert1 = diff_2d(x, y + 1, z - 1, t[x][y][z], hs_bf, hs_ube)
                    + t_2d(x, y + 1, z - 1, t[x][y][z], t[x][y + 1][z], hs_bf, hs_ube)
                    + t_2d(x, y + 1, z - 1, t[x][y][z], t[x][y][z - 1], hs_bf, hs_ube);
            if (alert1) {
                updated++;
                longflags[y * nz + nz + z - 1] = 1;
            }
        }
    }

    flag_fb = 0;
    y_start_fb = y_end;
    z_start_fb = z_begin;
    /* this direction has been examined: unset corresponding flag */

    return (updated);
}

/****end mail2/3***-------------------------END_X_SIDE()--------------------*/

/*----------------------------------------------Y_SIDE()--------------------*/

static int
y_side(x_begin, x_end, z_begin, z_end, y, future)

int x_begin, x_end, z_begin, z_end, y, future;

/* Propagates computations from side y-future to side y        */
/* between *_begin and *_end coordinates. Returns a nonzero    */
/* integer if something actually happened (a time was lowered).*/
/* Extensions _bb, _fb etc... define simple orientation rules: */
/* _bf means backwards along x axis and forwards along z axis. */
/* See complete comments in function x_side().                 */

{
    int
    updated,
            longhead,
            y0,
            y_s,
            x, z,
            sign_ff, sign_bf, sign_bb, sign_fb,
            test,
            past;
    GRID_FLOAT_TYPE
    hs_ff, hs_bf, hs_bb, hs_fb;

    if (reverse_order == 0)
        current_side_limit = y + future;
    updated = 0;
    y0 = y - future;
    if (future == 1) y_s = y0;
    else y_s = y;

    flag_fb = flag_bf = flag_ff = flag_bb = 0;
    x_start_ff = x_start_fb = x_end;
    x_start_bf = x_start_bb = x_begin;
    z_start_ff = z_start_bf = z_end;
    z_start_fb = z_start_bb = z_begin;

    /* First Step: "a-causal" stencils */

    for (x = x_begin; x <= x_end; x++) {
        for (z = z_begin; z <= z_end; z++) {
            hs_ff = hs[x][y_s][z];
            if (x > 0) hs_bf = hs[x - 1][y_s][z];
            else hs_bf = INFINITY;
            if (x > 0 && z > 0) hs_bb = hs[x - 1][y_s][z - 1];
            else hs_bb = INFINITY;
            if (z > 0) hs_fb = hs[x][y_s][z - 1];
            else hs_fb = INFINITY;

            sign_ff = sign_bf = sign_bb = sign_fb = 0;

            /* illuminate first neighbours */
            /* 1 1D transmission and 4 partial 3D transmission */
            updated += t_1d(x, y, z, t[x][y0][z], hs_ff, hs_bf, hs_bb, hs_fb);
            if (x < x_end && z < z_end)
                updated += t_3d_part1(x, y, z,
                    t[x][y0][z], t[x + 1][y0][z], t[x][y0][z + 1], hs_ff);
            if (x > x_begin && z < z_end)
                updated += t_3d_part1(x, y, z,
                    t[x][y0][z], t[x - 1][y0][z], t[x][y0][z + 1], hs_bf);
            if (x > x_begin && z > z_begin)
                updated += t_3d_part1(x, y, z,
                    t[x][y0][z], t[x - 1][y0][z], t[x][y0][z - 1], hs_bb);
            if (x < x_end && z > z_begin)
                updated += t_3d_part1(x, y, z,
                    t[x][y0][z], t[x + 1][y0][z], t[x][y0][z - 1], hs_fb);

            /* illuminate second neighbours */
            /* 4 2D diffraction and 4 2D transmission */
            if (x < x_end && t[x][y0][z] <= t[x + 1][y0][z]) {
                sign_fb++;
                sign_ff++;
                if (x < x_start_ff) x_start_ff = x;
                if (x < x_start_fb) x_start_fb = x;
                updated += diff_2d(x + 1, y, z, t[x][y0][z], hs_ff, hs_fb);
                updated += t_2d(x + 1, y, z, t[x][y0][z], t[x + 1][y0][z], hs_ff, hs_fb);
            }
            if (x > x_begin && t[x][y0][z] <= t[x - 1][y0][z]) {
                sign_bb++;
                sign_bf++;
                if (x > x_start_bf) x_start_bf = x;
                if (x > x_start_bb) x_start_bb = x;
                updated += diff_2d(x - 1, y, z, t[x][y0][z], hs_bf, hs_bb);
                updated += t_2d(x - 1, y, z, t[x][y0][z], t[x - 1][y0][z], hs_bf, hs_bb);
            }
            if (z < z_end && t[x][y0][z] <= t[x][y0][z + 1]) {
                sign_bf++;
                sign_ff++;
                if (z < z_start_ff) z_start_ff = z;
                if (z < z_start_bf) z_start_bf = z;
                updated += diff_2d(x, y, z + 1, t[x][y0][z], hs_ff, hs_bf);
                updated += t_2d(x, y, z + 1, t[x][y0][z], t[x][y0][z + 1], hs_ff, hs_bf);
            }
            if (z > z_begin && t[x][y0][z] <= t[x][y0][z - 1]) {
                sign_bb++;
                sign_fb++;
                if (z > z_start_fb) z_start_fb = z;
                if (z > z_start_bb) z_start_bb = z;
                updated += diff_2d(x, y, z - 1, t[x][y0][z], hs_bb, hs_fb);
                updated += t_2d(x, y, z - 1, t[x][y0][z], t[x][y0][z - 1], hs_bb, hs_fb);
            }

            /* illuminate third neighbours */
            /* 4 3D point diffraction, 8 3D edge diffraction and 12 3D transmission */
            if (sign_ff == 2) {
                flag_ff = 1;
                updated += point_diff(x + 1, y, z + 1, t[x][y0][z], hs_ff)
                        + edge_diff(x + 1, y, z + 1, t[x][y0][z], t[x + 1][y0][z], hs_ff)
                        + edge_diff(x + 1, y, z + 1, t[x][y0][z], t[x][y0][z + 1], hs_ff)
                        + t_3d_part2(x + 1, y, z + 1, t[x][y0][z], t[x + 1][y0][z],
                        t[x][y0][z + 1], t[x + 1][y0][z + 1], hs_ff);
            }
            if (sign_bf == 2) {
                flag_bf = 1;
                updated += point_diff(x - 1, y, z + 1, t[x][y0][z], hs_bf)
                        + edge_diff(x - 1, y, z + 1, t[x][y0][z], t[x - 1][y0][z], hs_bf)
                        + edge_diff(x - 1, y, z + 1, t[x][y0][z], t[x][y0][z + 1], hs_bf)
                        + t_3d_part2(x - 1, y, z + 1, t[x][y0][z], t[x - 1][y0][z],
                        t[x][y0][z + 1], t[x - 1][y0][z + 1], hs_bf);
            }
            if (sign_bb == 2) {
                flag_bb = 1;
                updated += point_diff(x - 1, y, z - 1, t[x][y0][z], hs_bb)
                        + edge_diff(x - 1, y, z - 1, t[x][y0][z], t[x - 1][y0][z], hs_bb)
                        + edge_diff(x - 1, y, z - 1, t[x][y0][z], t[x][y0][z - 1], hs_bb)
                        + t_3d_part2(x - 1, y, z - 1, t[x][y0][z], t[x - 1][y0][z],
                        t[x][y0][z - 1], t[x - 1][y0][z - 1], hs_bb);
            }
            if (sign_fb == 2) {
                flag_fb = 1;
                updated += point_diff(x + 1, y, z - 1, t[x][y0][z], hs_fb)
                        + edge_diff(x + 1, y, z - 1, t[x][y0][z], t[x + 1][y0][z], hs_fb)
                        + edge_diff(x + 1, y, z - 1, t[x][y0][z], t[x][y0][z - 1], hs_fb)
                        + t_3d_part2(x + 1, y, z - 1, t[x][y0][z], t[x + 1][y0][z],
                        t[x][y0][z - 1], t[x + 1][y0][z - 1], hs_fb);
            }
        }
    }

    /* Second Step: causal propagation */

    for (x = 0; x < nx * nz; x++) longflags[x] = 0;

    if (x_begin == x_end || z_begin == z_end) {
        flag_ff = flag_fb = flag_bf = flag_bb = 1;
        x_start_ff = x_start_fb = x_begin;
        x_start_bf = x_start_bb = x_end;
        z_start_ff = z_start_bf = z_begin;
        z_start_fb = z_start_bb = z_end;
    }

    do {
        test = 0;
        if (flag_ff) {
            test++;
            if (VERBOSE) printf("ff ");
            updated += scan_y_ff(x_start_ff, x_end, z_start_ff, z_end, y0, y, y_s);
        }
        if (flag_fb) {
            test++;
            if (VERBOSE) printf("fb ");
            updated += scan_y_fb(x_start_fb, x_end, z_begin, z_start_fb, y0, y, y_s);
        }
        if (flag_bb) {
            test++;
            if (VERBOSE) printf("bb ");
            updated += scan_y_bb(x_begin, x_start_bb, z_begin, z_start_bb, y0, y, y_s);
        }
        if (flag_bf) {
            test++;
            if (VERBOSE) printf("bf ");
            updated += scan_y_bf(x_begin, x_start_bf, z_start_bf, z_end, y0, y, y_s);
        }
    } while (test);

    sum_updated += updated;

    /* Third Step: Reverse propagation, if necessary */

    for (x = longhead = 0; x < nx * nz; x++) longhead += longflags[x];

    if (longhead) {

        reverse_order++;

        if (VERBOSE) printf("\nReverse#%d from y_side %d", reverse_order, y);
        past = -future;
        for (y = y0; y != current_side_limit; y += past) {
            if (y < 0 || y >= ny) break;
            if (VERBOSE) printf("\nupdate side y=%d: ", y);
            if (y_side(x_begin, x_end, z_begin, z_end, y, past) == 0) break;
            if (VERBOSE) printf("y=%d <R#%d>updated.", y, reverse_order);
        }
        if (VERBOSE) printf("\nEnd Reverse#%d\n", reverse_order);

        reverse_order--;

    }

    return (updated);

}

/*--------------------------------------Y_SIDE() : SCAN_Y_EE()--------------*/

static int
scan_y_ff(x_start, x_end, z_start, z_end, y0, y, y_s)
int x_start, x_end, z_start, z_end, y0, y, y_s;

{
    int
    updated = 0,
            alert0, alert1,
            x, z,
            y_sf;
    GRID_FLOAT_TYPE
    hs_bf, hs_bb, hs_fb,
            hs_ube, hs_ubb, hs_ueb;

    y_sf = y_s + y - y0;

    hs_bb = hs_ubb = hs_ube = INFINITY;
    for (x = x_start, z = z_start; x < x_end; x++) {
        hs_bf = hs[x][y_s][z];
        if (z) hs_bb = hs[x][y_s][z - 1];
        if (y_sf >= 0 && y_sf < nmesh_y) {
            hs_ube = hs[x][y_sf][z];
            if (z) hs_ubb = hs[x][y_sf][z - 1];
        }
        alert1 = t_1d(x + 1, y, z, t[x][y][z], hs_bb, hs_bf, hs_ubb, hs_ube);
        alert0 = t_2d(x + 1, y, z, t[x][y0][z], t[x][y][z], hs_bb, hs_bf);
        updated += alert0 + alert1;
        if (alert1) longflags[x * nz + nz + z] = 1;
        if (alert0) longflags[x * nz + nz + z] = 0;
    }

    hs_bb = hs_ubb = hs_ueb = INFINITY;
    for (x = x_start, z = z_start; z < z_end; z++) {
        hs_fb = hs[x][y_s][z];
        if (x) hs_bb = hs[x - 1][y_s][z];
        if (y_sf >= 0 && y_sf < nmesh_y) {
            hs_ueb = hs[x][y_sf][z];
            if (x) hs_ubb = hs[x - 1][y_sf][z];
        }
        alert1 = t_1d(x, y, z + 1, t[x][y][z], hs_bb, hs_fb, hs_ubb, hs_ueb);
        alert0 = t_2d(x, y, z + 1, t[x][y0][z], t[x][y][z], hs_bb, hs_fb);
        updated += alert0 + alert1;
        if (alert1) longflags[x * nz + z + 1] = 1;
        if (alert0) longflags[x * nz + z + 1] = 0;
    }

    for (x = x_start; x < x_end; x++) {
        for (z = z_start; z < z_end; z++) {
            hs_bb = hs[x][y_s][z];
            hs_bf = hs[x][y_s][z + 1];
            hs_fb = hs[x + 1][y_s][z];
            if (y_sf >= 0 && y_sf < nmesh_y) {
                hs_ubb = hs[x][y_sf][z];
                hs_ube = hs[x][y_sf][z + 1];
                hs_ueb = hs[x + 1][y_sf][z];
            } else hs_ubb = hs_ube = hs_ueb = INFINITY;

            /* bulk waves: 1 3D edge diffraction and 2 (*4) 3D transmission */
            alert0 = edge_diff(x + 1, y, z + 1, t[x][y0][z], t[x][y][z], hs_bb)
                    + t_3d(x + 1, y, z + 1, t[x][y0][z],
                    t[x + 1][y0][z], t[x][y][z], t[x + 1][y][z], hs_bb)
                    + t_3d(x + 1, y, z + 1, t[x][y0][z],
                    t[x][y0][z + 1], t[x][y][z], t[x][y][z + 1], hs_bb);
            if (alert0) {
                updated++;
                longflags[x * nz + nz + z + 1] = 0;
            }

            /* interface waves along x_side: 1 1D transmission and 1 2D transmission */
            alert1 = t_1d(x + 1, y, z + 1, t[x][y][z + 1], hs_bb, hs_bf, hs_ubb, hs_ube);
            alert0 = t_2d(x + 1, y, z + 1, t[x][y0][z + 1], t[x][y][z + 1], hs_bb, hs_bf);
            if (alert0 + alert1) {
                updated += alert0 + alert1;
                flag_fb++; /* this scan must be (re-)examined */
                if (x_start_fb > x) x_start_fb = x;
                if (z_start_fb < z + 1) z_start_fb = z + 1;
                if (alert1) longflags[x * nz + nz + z + 1] = 1;
                else longflags[x * nz + nz + z + 1] = 0;
            }

            /* interface waves along z_side: 1 1D transmission and 1 2D transmission */
            alert1 = t_1d(x + 1, y, z + 1, t[x + 1][y][z], hs_bb, hs_fb, hs_ubb, hs_ueb);
            alert0 = t_2d(x + 1, y, z + 1, t[x + 1][y0][z], t[x + 1][y][z], hs_bb, hs_fb);
            if (alert0 + alert1) {
                updated += alert0 + alert1;
                flag_bf++; /* this scan must be (re-)examined */
                if (x_start_bf < x + 1) x_start_bf = x + 1;
                if (z_start_bf > z) z_start_bf = z;
                if (alert1) longflags[x * nz + nz + z + 1] = 1;
                else longflags[x * nz + nz + z + 1] = 0;
            }

            /* interface waves along y_side : 2 2D transmission and 1 2D diffraction */
            alert1 = diff_2d(x + 1, y, z + 1, t[x][y][z], hs_bb, hs_ubb)
                    + t_2d(x + 1, y, z + 1, t[x][y][z], t[x + 1][y][z], hs_bb, hs_ubb)
                    + t_2d(x + 1, y, z + 1, t[x][y][z], t[x][y][z + 1], hs_bb, hs_ubb);
            if (alert1) {
                updated += alert1;
                longflags[x * nz + nz + z + 1] = 1;
            }
        }
    }

    flag_ff = 0;
    x_start_ff = x_end;
    z_start_ff = z_end;

    return (updated);
}

/*--------------------------------------Y_SIDE() : SCAN_Y_BE()--------------*/

static int
scan_y_bf(x_begin, x_start, z_start, z_end, y0, y, y_s)
int x_start, x_begin, z_start, z_end, y0, y, y_s;

{
    int
    updated = 0,
            alert0, alert1,
            x, z,
            y_sf;
    GRID_FLOAT_TYPE
    hs_ff, hs_bb, hs_fb,
            hs_uee, hs_ubb, hs_ueb;

    y_sf = y_s + y - y0;

    hs_fb = hs_uee = hs_ueb = INFINITY;
    for (x = x_start, z = z_start; x > x_begin; x--) {
        hs_ff = hs[x - 1][y_s][z];
        if (z) hs_fb = hs[x - 1][y_s][z - 1];
        if (y_sf >= 0 && y_sf < nmesh_y) {
            hs_uee = hs[x - 1][y_sf][z];
            if (z) hs_ueb = hs[x - 1][y_sf][z - 1];
        }
        alert1 = t_1d(x - 1, y, z, t[x][y][z], hs_fb, hs_ff, hs_ueb, hs_uee);
        alert0 = t_2d(x - 1, y, z, t[x][y0][z], t[x][y][z], hs_fb, hs_ff);
        updated += alert0 + alert1;
        if (alert1) longflags[x * nz - nz + z] = 1;
        if (alert0) longflags[x * nz - nz + z] = 0;
    }

    hs_bb = hs_ubb = hs_ueb = INFINITY;
    for (x = x_start, z = z_start; z < z_end; z++) {
        if (x) hs_bb = hs[x - 1][y_s][z];
        hs_fb = hs[x][y_s][z];
        if (y_sf >= 0 && y_sf < nmesh_y) {
            if (x) hs_ubb = hs[x - 1][y_sf][z];
            hs_ueb = hs[x][y_sf][z];
        }
        alert1 = t_1d(x, y, z + 1, t[x][y][z], hs_bb, hs_fb, hs_ubb, hs_ueb);
        alert0 = t_2d(x, y, z + 1, t[x][y0][z], t[x][y][z], hs_bb, hs_fb);
        updated += alert0 + alert1;
        if (alert1) longflags[x * nz + z + 1] = 1;
        if (alert0) longflags[x * nz + z + 1] = 0;
    }

    for (x = x_start; x > x_begin; x--) {
        for (z = z_start; z < z_end; z++) {

            hs_ff = hs[x - 1][y_s][z + 1];
            if (x > 1) hs_bb = hs[x - 2][y_s][z];
            else hs_bb = INFINITY;
            hs_fb = hs[x - 1][y_s][z];
            if (y_sf >= 0 && y_sf < nmesh_y) {
                hs_uee = hs[x - 1][y_sf][z + 1];
                if (x > 1) hs_ubb = hs[x - 2][y_sf][z];
                else hs_ubb = INFINITY;
                hs_ueb = hs[x - 1][y_sf][z];
            } else hs_ubb = hs_uee = hs_ueb = INFINITY;

            /* bulk waves: 1 3D edge diffraction and 2 (*4) 3D transmission */
            alert0 = edge_diff(x - 1, y, z + 1, t[x][y0][z], t[x][y][z], hs_fb)
                    + t_3d(x - 1, y, z + 1, t[x][y0][z],
                    t[x - 1][y0][z], t[x][y][z], t[x - 1][y][z], hs_fb)
                    + t_3d(x - 1, y, z + 1, t[x][y0][z],
                    t[x][y0][z + 1], t[x][y][z], t[x][y][z + 1], hs_fb);
            if (alert0) {
                updated += alert0;
                longflags[x * nz - nz + z + 1] = 0;
            }

            /* interface waves along x_side: 1 1D transmission and 1 2D transmission */
            alert1 = t_1d(x - 1, y, z + 1, t[x][y][z + 1], hs_ff, hs_fb, hs_uee, hs_ueb);
            alert0 = t_2d(x - 1, y, z + 1, t[x][y0][z + 1], t[x][y][z + 1], hs_ff, hs_fb);
            if (alert0 + alert1) {
                updated += alert0 + alert1;
                flag_bb++; /* this scan must be (re-)examined */
                if (x_start_bb < x) x_start_bb = x;
                if (z_start_bb < z + 1) z_start_bb = z + 1;
                if (alert1) longflags[x * nz - nz + z + 1] = 1;
                else longflags[x * nz - nz + z + 1] = 0;
            }

            /* interface waves along z_side: 1 1D transmission and 1 2D transmission */
            alert1 = t_1d(x - 1, y, z + 1, t[x - 1][y][z], hs_bb, hs_fb, hs_ubb, hs_ueb);
            alert0 = t_2d(x - 1, y, z + 1, t[x - 1][y0][z], t[x - 1][y][z], hs_bb, hs_fb);
            if (alert0 + alert1) {
                updated += alert0 + alert1;
                flag_ff++; /* this scan must be (re-)examined */
                if (x_start_ff > x - 1) x_start_ff = x - 1;
                if (z_start_ff > z) z_start_ff = z;
                if (alert1) longflags[x * nz - nz + z + 1] = 1;
                else longflags[x * nz - nz + z + 1] = 0;
            }

            /* interface waves along y_side : 2 2D transmission and 1 2D diffraction */
            alert1 = diff_2d(x - 1, y, z + 1, t[x][y][z], hs_fb, hs_ueb)
                    + t_2d(x - 1, y, z + 1, t[x][y][z], t[x - 1][y][z], hs_fb, hs_ueb)
                    + t_2d(x - 1, y, z + 1, t[x][y][z], t[x][y][z + 1], hs_fb, hs_ueb);
            if (alert1) {
                updated += alert1;
                longflags[x * nz - nz + z + 1] = 1;
            }
        }
    }

    flag_bf = 0;
    x_start_bf = x_begin;
    z_start_bf = z_end;

    return (updated);
}

/*--------------------------------------Y_SIDE() : SCAN_Y_BB()--------------*/

static int
scan_y_bb(x_begin, x_start, z_begin, z_start, y0, y, y_s)
int x_start, x_begin, z_begin, z_start, y0, y, y_s;

{
    int
    updated = 0,
            alert0, alert1,
            x, z,
            y_sf;
    GRID_FLOAT_TYPE
    hs_ff, hs_bf, hs_fb,
            hs_uee, hs_ube, hs_ueb;

    y_sf = y_s + y - y0;

    hs_ff = hs_uee = hs_ueb = INFINITY;
    for (x = x_start, z = z_start; x > x_begin; x--) {
        if (z) hs_ff = hs[x - 1][y_s][z - 1];
        hs_fb = hs[x - 1][y_s][z];
        if (y_sf >= 0 && y_sf < nmesh_y) {
            if (z) hs_uee = hs[x - 1][y_sf][z - 1];
            hs_ueb = hs[x - 1][y_sf][z];
        }
        alert1 = t_1d(x - 1, y, z, t[x][y][z], hs_fb, hs_ff, hs_ueb, hs_uee);
        alert0 = t_2d(x - 1, y, z, t[x][y0][z], t[x][y][z], hs_fb, hs_ff);
        updated += alert0 + alert1;
        if (alert1) longflags[x * nz - nz + z] = 1;
        if (alert0) longflags[x * nz - nz + z] = 0;
    }

    hs_ff = hs_uee = hs_ube = INFINITY;
    for (x = x_start, z = z_start; z > z_begin; z--) {
        if (x) hs_ff = hs[x - 1][y_s][z - 1];
        hs_bf = hs[x][y_s][z - 1];
        if (y_sf >= 0 && y_sf < nmesh_y) {
            if (x) hs_uee = hs[x - 1][y_sf][z - 1];
            hs_ube = hs[x][y_sf][z - 1];
        }
        alert1 = t_1d(x, y, z - 1, t[x][y][z], hs_bf, hs_ff, hs_ube, hs_uee);
        alert0 = t_2d(x, y, z - 1, t[x][y0][z], t[x][y][z], hs_bf, hs_ff);
        updated += alert0 + alert1;
        if (alert1) longflags[x * nz + z - 1] = 1;
        if (alert0) longflags[x * nz + z - 1] = 0;
    }

    for (x = x_start; x > x_begin; x--) {
        for (z = z_start; z > z_begin; z--) {

            hs_ff = hs[x - 1][y_s][z - 1];
            if (x > 1) hs_bf = hs[x - 2][y_s][z - 1];
            else hs_bf = INFINITY;
            if (z > 1) hs_fb = hs[x - 1][y_s][z - 2];
            else hs_fb = INFINITY;
            if (y_sf >= 0 && y_sf < nmesh_y) {
                hs_uee = hs[x - 1][y_sf][z - 1];
                if (x > 1) hs_ube = hs[x - 2][y_sf][z - 1];
                else hs_ube = INFINITY;
                if (z > 1) hs_ueb = hs[x - 1][y_sf][z - 2];
                else hs_ueb = INFINITY;
            } else hs_uee = hs_ube = hs_ueb = INFINITY;

            /* bulk waves: 1 3D edge diffraction and 2 (*4) 3D transmission */
            alert0 = edge_diff(x - 1, y, z - 1, t[x][y0][z], t[x][y][z], hs_ff)
                    + t_3d(x - 1, y, z - 1, t[x][y0][z],
                    t[x - 1][y0][z], t[x][y][z], t[x - 1][y][z], hs_ff)
                    + t_3d(x - 1, y, z - 1, t[x][y0][z],
                    t[x][y0][z - 1], t[x][y][z], t[x][y][z - 1], hs_ff);
            if (alert0) {
                updated += alert0;
                longflags[x * nz - nz + z - 1] = 0;
            }

            /* interface waves along x_side: 1 1D transmission and 1 2D transmission */
            alert1 = t_1d(x - 1, y, z - 1, t[x][y][z - 1], hs_ff, hs_fb, hs_uee, hs_ueb);
            alert0 = t_2d(x - 1, y, z - 1, t[x][y0][z - 1], t[x][y][z - 1], hs_ff, hs_fb);
            if (alert0 + alert1) {
                updated += alert0 + alert1;
                flag_bf++; /* this scan must be (re-)examined */
                if (x_start_bf < x) x_start_bf = x;
                if (z_start_bf > z - 1) z_start_bf = z - 1;
                if (alert1) longflags[x * nz - nz + z - 1] = 1;
                else longflags[x * nz - nz + z - 1] = 0;
            }

            /* interface waves along z_side: 1 1D transmission and 1 2D transmission */
            alert1 = t_1d(x - 1, y, z - 1, t[x - 1][y][z], hs_ff, hs_bf, hs_uee, hs_ube);
            alert0 = t_2d(x - 1, y, z - 1, t[x - 1][y0][z], t[x - 1][y][z], hs_ff, hs_bf);
            if (alert0 + alert1) {
                updated += alert0 + alert1;
                flag_fb++; /* this scan must be (re-)examined */
                if (x_start_fb > x - 1) x_start_fb = x - 1;
                if (z_start_fb < z) z_start_fb = z;
                if (alert1) longflags[x * nz - nz + z - 1] = 1;
                else longflags[x * nz - nz + z - 1] = 0;
            }

            /* interface waves along y_side : 2 2D transmission and 1 2D diffraction */
            alert1 = diff_2d(x - 1, y, z - 1, t[x][y][z], hs_ff, hs_uee)
                    + t_2d(x - 1, y, z - 1, t[x][y][z], t[x - 1][y][z], hs_ff, hs_uee)
                    + t_2d(x - 1, y, z - 1, t[x][y][z], t[x][y][z - 1], hs_ff, hs_uee);
            if (alert1) {
                updated += alert1;
                longflags[x * nz - nz + z - 1] = 1;
            }
        }
    }

    flag_bb = 0;
    x_start_bb = x_begin;
    z_start_bb = z_begin;

    return (updated);
}

/*--------------------------------------Y_SIDE() : SCAN_Y_EB()--------------*/

static int
scan_y_fb(x_start, x_end, z_begin, z_start, y0, y, y_s)
int x_start, x_end, z_begin, z_start, y0, y, y_s;

{
    int
    updated = 0,
            alert0, alert1,
            x, z,
            y_sf;
    GRID_FLOAT_TYPE
    hs_ff, hs_bb, hs_bf,
            hs_uee, hs_ubb, hs_ube;

    y_sf = y_s + y - y0;

    hs_bf = hs_ubb = hs_ube = INFINITY;
    for (x = x_start, z = z_start; x < x_end; x++) {
        if (z) hs_bf = hs[x][y_s][z - 1];
        hs_bb = hs[x][y_s][z];
        if (y_sf >= 0 && y_sf < nmesh_y) {
            if (z) hs_ube = hs[x][y_sf][z - 1];
            hs_ubb = hs[x][y_sf][z];
        }
        alert1 = t_1d(x + 1, y, z, t[x][y][z], hs_bf, hs_bb, hs_ube, hs_ubb);
        alert0 = t_2d(x + 1, y, z, t[x][y0][z], t[x][y][z], hs_bf, hs_bb);
        updated += alert0 + alert1;
        if (alert1) longflags[x * nz + nz + z] = 1;
        if (alert0) longflags[x * nz + nz + z] = 0;
    }

    hs_ff = hs_uee = hs_ube = INFINITY;
    for (x = x_start, z = z_start; z > z_begin; z--) {
        if (x) hs_ff = hs[x - 1][y_s][z - 1];
        hs_bf = hs[x][y_s][z - 1];
        if (y_sf >= 0 && y_sf < nmesh_y) {
            if (x) hs_uee = hs[x - 1][y_sf][z - 1];
            hs_ube = hs[x][y_sf][z - 1];
        }
        alert1 = t_1d(x, y, z - 1, t[x][y][z], hs_bf, hs_ff, hs_ube, hs_uee);
        alert0 = t_2d(x, y, z - 1, t[x][y0][z], t[x][y][z], hs_bf, hs_ff);
        updated += alert0 + alert1;
        if (alert1) longflags[x * nz + z - 1] = 1;
        if (alert0) longflags[x * nz + z - 1] = 0;
    }

    for (x = x_start; x < x_end; x++) {
        for (z = z_start; z > z_begin; z--) {

            hs_bf = hs[x][y_s][z - 1];
            if (z > 1) hs_bb = hs[x][y_s][z - 2];
            else hs_bb = INFINITY;
            hs_ff = hs[x + 1][y_s][z - 1];
            if (y_sf >= 0 && y_sf < nmesh_y) {
                hs_ube = hs[x][y_sf][z - 1];
                if (z > 1) hs_ubb = hs[x][y_sf][z - 2];
                else hs_ubb = INFINITY;
                hs_uee = hs[x + 1][y_sf][z - 1];
            } else hs_ubb = hs_ube = hs_uee = INFINITY;

            /* bulk waves: 1 3D edge diffraction and 2 (*4) 3D transmission */
            alert0 = edge_diff(x + 1, y, z - 1, t[x][y0][z], t[x][y][z], hs_bf)
                    + t_3d(x + 1, y, z - 1, t[x][y0][z],
                    t[x + 1][y0][z], t[x][y][z], t[x + 1][y][z], hs_bf)
                    + t_3d(x + 1, y, z - 1, t[x][y0][z],
                    t[x][y0][z - 1], t[x][y][z], t[x][y][z - 1], hs_bf);
            if (alert0) {
                updated += alert0;
                longflags[x * nz + nz + z - 1] = 0;
            }

            /* interface waves along x_side: 1 1D transmission and 1 2D transmission */
            alert1 = t_1d(x + 1, y, z - 1, t[x][y][z - 1], hs_bb, hs_bf, hs_ubb, hs_ube);
            alert0 = t_2d(x + 1, y, z - 1, t[x][y0][z - 1], t[x][y][z - 1], hs_bb, hs_bf);
            if (alert0 + alert1) {
                updated += alert0 + alert1;
                flag_ff++; /* this scan must be (re-)examined */
                if (x_start_ff > x) x_start_ff = x;
                if (z_start_ff > z - 1) z_start_ff = z - 1;
                if (alert1) longflags[x * nz + nz + z - 1] = 1;
                else longflags[x * nz + nz + z - 1] = 0;
            }

            /* interface waves along z_side: 1 1D transmission and 1 2D transmission */
            alert1 = t_1d(x + 1, y, z - 1, t[x + 1][y][z], hs_ff, hs_bf, hs_uee, hs_ube);
            alert0 = t_2d(x + 1, y, z - 1, t[x + 1][y0][z], t[x + 1][y][z], hs_ff, hs_bf);
            if (alert0 + alert1) {
                updated += alert0 + alert1;
                flag_bb++; /* this scan must be (re-)examined */
                if (x_start_bb < x + 1) x_start_bb = x + 1;
                if (z_start_bb < z) z_start_bb = z;
                if (alert1) longflags[x * nz + nz + z - 1] = 1;
                else longflags[x * nz + nz + z - 1] = 0;
            }

            /* interface waves along y_side : 2 2D transmission and 1 2D diffraction */
            alert1 = diff_2d(x + 1, y, z - 1, t[x][y][z], hs_bf, hs_ube)
                    + t_2d(x + 1, y, z - 1, t[x][y][z], t[x + 1][y][z], hs_bf, hs_ube)
                    + t_2d(x + 1, y, z - 1, t[x][y][z], t[x][y][z - 1], hs_bf, hs_ube);
            if (alert1) {
                updated += alert1;
                longflags[x * nz + nz + z - 1] = 1;
            }
        }
    }

    flag_fb = 0;
    x_start_fb = x_end;
    z_start_fb = z_begin;
    /* this scan has been examined */

    return (updated);
}

/*--------------------------------------------END_Y_SIDE()--------------------*/

/*----------------------------------------------Z_SIDE()--------------------*/

static int
z_side(x_begin, x_end, y_begin, y_end, z, future)

int x_begin, x_end, y_begin, y_end, z, future;

/* Propagates computations from side z-future to side z.       */
/* between *_begin and *_end coordinates. Returns a nonzero    */
/* integer if something actually happened (a time was lowered).*/
/* Extensions _bb, _fb etc... define simple orientation rules: */
/* _bf means backwards along x axis and forwards along y axis. */
/* See complete comments in function x_side().                 */

{
    int
    updated,
            longhead,
            z0,
            z_s,
            x, y,
            sign_ff, sign_bf, sign_bb, sign_fb,
            test,
            past;
    GRID_FLOAT_TYPE
    hs_ff, hs_bf, hs_bb, hs_fb;

    if (reverse_order == 0)
        current_side_limit = z + future;
    updated = 0;
    z0 = z - future;
    if (future == 1) z_s = z0;
    else z_s = z;

    flag_fb = flag_bf = flag_ff = flag_bb = 0;
    x_start_ff = x_start_fb = x_end;
    x_start_bf = x_start_bb = x_begin;
    y_start_ff = y_start_bf = y_end;
    y_start_fb = y_start_bb = y_begin;

    /* First Step: "a-causal" stencils */

    for (x = x_begin; x <= x_end; x++) {
        for (y = y_begin; y <= y_end; y++) {

            hs_ff = hs[x][y][z_s];
            if (x > 0) hs_bf = hs[x - 1][y][z_s];
            else hs_bf = INFINITY;
            if (x > 0 && y > 0) hs_bb = hs[x - 1][y - 1][z_s];
            else hs_bb = INFINITY;
            if (y > 0) hs_fb = hs[x][y - 1][z_s];
            else hs_fb = INFINITY;
            sign_ff = sign_bf = sign_bb = sign_fb = 0;

            /* illuminate first neighbours */
            /* 1 1D transmission and 4 partial 3D transmission */
            updated += t_1d(x, y, z, t[x][y][z0], hs_ff, hs_bf, hs_bb, hs_fb);
            if (x < x_end && y < y_end)
                updated += t_3d_part1(x, y, z, t[x][y][z0],
                    t[x + 1][y][z0], t[x][y + 1][z0], hs_ff);
            if (x > x_begin && y < y_end)
                updated += t_3d_part1(x, y, z, t[x][y][z0],
                    t[x - 1][y][z0], t[x][y + 1][z0], hs_bf);
            if (x > x_begin && y > y_begin)
                updated += t_3d_part1(x, y, z, t[x][y][z0],
                    t[x - 1][y][z0], t[x][y - 1][z0], hs_bb);
            if (x < x_end && y > y_begin)
                updated += t_3d_part1(x, y, z, t[x][y][z0],
                    t[x + 1][y][z0], t[x][y - 1][z0], hs_fb);

            /* illuminate second neighbours */
            /* 4 2D diffraction and 4 2D transmission */
            if (x < x_end && t[x][y][z0] <= t[x + 1][y][z0]) {
                sign_fb++;
                sign_ff++;
                if (x < x_start_ff) x_start_ff = x;
                if (x < x_start_fb) x_start_fb = x;
                updated += diff_2d(x + 1, y, z, t[x][y][z0], hs_ff, hs_fb);
                updated += t_2d(x + 1, y, z, t[x][y][z0], t[x + 1][y][z0], hs_ff, hs_fb);
            }
            if (x > x_begin && t[x][y][z0] <= t[x - 1][y][z0]) {
                sign_bb++;
                sign_bf++;
                if (x > x_start_bf) x_start_bf = x;
                if (x > x_start_bb) x_start_bb = x;
                updated += diff_2d(x - 1, y, z, t[x][y][z0], hs_bf, hs_bb);
                updated += t_2d(x - 1, y, z, t[x][y][z0], t[x - 1][y][z0], hs_bf, hs_bb);
            }
            if (y < y_end && t[x][y][z0] <= t[x][y + 1][z0]) {
                sign_bf++;
                sign_ff++;
                if (y < y_start_ff) y_start_ff = y;
                if (y < y_start_bf) y_start_bf = y;
                updated += diff_2d(x, y + 1, z, t[x][y][z0], hs_ff, hs_bf);
                updated += t_2d(x, y + 1, z, t[x][y][z0], t[x][y + 1][z0], hs_ff, hs_bf);
            }
            if (y > y_begin && t[x][y][z0] <= t[x][y - 1][z0]) {
                sign_bb++;
                sign_fb++;
                if (y > y_start_fb) y_start_fb = y;
                if (y > y_start_bb) y_start_bb = y;
                updated += diff_2d(x, y - 1, z, t[x][y][z0], hs_bb, hs_fb);
                updated += t_2d(x, y - 1, z, t[x][y][z0], t[x][y - 1][z0], hs_bb, hs_fb);
            }

            /* illuminate third neighbours */
            /* 4 3D point diffraction, 8 3D edge diffraction and 12 3D transmission */
            if (sign_ff == 2) {
                flag_ff = 1;
                updated += point_diff(x + 1, y + 1, z, t[x][y][z0], hs_ff)
                        + edge_diff(x + 1, y + 1, z, t[x][y][z0], t[x + 1][y][z0], hs_ff)
                        + edge_diff(x + 1, y + 1, z, t[x][y][z0], t[x][y + 1][z0], hs_ff)
                        + t_3d_part2(x + 1, y + 1, z, t[x][y][z0], t[x + 1][y][z0],
                        t[x][y + 1][z0], t[x + 1][y + 1][z0], hs_ff);
            }
            if (sign_bf == 2) {
                flag_bf = 1;
                updated += point_diff(x - 1, y + 1, z, t[x][y][z0], hs_bf)
                        + edge_diff(x - 1, y + 1, z, t[x][y][z0], t[x - 1][y][z0], hs_bf)
                        + edge_diff(x - 1, y + 1, z, t[x][y][z0], t[x][y + 1][z0], hs_bf)
                        + t_3d_part2(x - 1, y + 1, z, t[x][y][z0], t[x - 1][y][z0],
                        t[x][y + 1][z0], t[x - 1][y + 1][z0], hs_bf);
            }
            if (sign_bb == 2) {
                flag_bb = 1;
                updated += point_diff(x - 1, y - 1, z, t[x][y][z0], hs_bb)
                        + edge_diff(x - 1, y - 1, z, t[x][y][z0], t[x - 1][y][z0], hs_bb)
                        + edge_diff(x - 1, y - 1, z, t[x][y][z0], t[x][y - 1][z0], hs_bb)
                        + t_3d_part2(x - 1, y - 1, z, t[x][y][z0], t[x - 1][y][z0],
                        t[x][y - 1][z0], t[x - 1][y - 1][z0], hs_bb);
            }
            if (sign_fb == 2) {
                flag_fb = 1;
                updated += point_diff(x + 1, y - 1, z, t[x][y][z0], hs_fb)
                        + edge_diff(x + 1, y - 1, z, t[x][y][z0], t[x + 1][y][z0], hs_fb)
                        + edge_diff(x + 1, y - 1, z, t[x][y][z0], t[x][y - 1][z0], hs_fb)
                        + t_3d_part2(x + 1, y - 1, z, t[x][y][z0], t[x + 1][y][z0],
                        t[x][y - 1][z0], t[x + 1][y - 1][z0], hs_fb);
            }
        }
    }

    /* Second Step: causal propagation */

    for (x = 0; x < nx * ny; x++) longflags[x] = 0;

    if (x_begin == x_end || y_begin == y_end) {
        flag_ff = flag_fb = flag_bf = flag_bb = 1;
        x_start_ff = x_start_fb = x_begin;
        x_start_bf = x_start_bb = x_end;
        y_start_ff = y_start_bf = y_begin;
        y_start_fb = y_start_bb = y_end;
    }

    do {
        test = 0;
        if (flag_ff) {
            test++;
            if (VERBOSE) printf("ff ");
            updated += scan_z_ff(x_start_ff, x_end, y_start_ff, y_end, z0, z, z_s);
        }
        if (flag_fb) {
            test++;
            if (VERBOSE) printf("fb ");
            updated += scan_z_fb(x_start_fb, x_end, y_begin, y_start_fb, z0, z, z_s);
        }
        if (flag_bb) {
            test++;
            if (VERBOSE) printf("bb ");
            updated += scan_z_bb(x_begin, x_start_bb, y_begin, y_start_bb, z0, z, z_s);
        }
        if (flag_bf) {
            test++;
            if (VERBOSE) printf("bf ");
            updated += scan_z_bf(x_begin, x_start_bf, y_start_bf, y_end, z0, z, z_s);
        }
    } while (test);

    sum_updated += updated;

    /* Third Step: Reverse Propagation if necessary */

    for (x = longhead = 0; x < nx * ny; x++) longhead += longflags[x];

    if (longhead) {

        reverse_order++;

        if (VERBOSE) printf("\nReverse#%d from z_side %d", reverse_order, z);
        past = -future;
        for (z = z0; z != current_side_limit; z += past) {
            if (z < 0 || z >= nz) break;
            if (VERBOSE) printf("\nupdate side z=%d: ", z);
            if (z_side(x_begin, x_end, y_begin, y_end, z, past) == 0) break;
            if (VERBOSE) printf("z=%d <R#%d>updated.", z, reverse_order);
        }
        if (VERBOSE) printf("\nEnd Reverse#%d\n", reverse_order);

        reverse_order--;

    }

    return (updated);

}

/*--------------------------------------Z_SIDE() : SCAN_Z_EE()--------------*/

static int
scan_z_ff(x_start, x_end, y_start, y_end, z0, z, z_s)
int x_start, x_end, y_start, y_end, z0, z, z_s;

{
    int
    updated = 0,
            alert0, alert1,
            x, y,
            z_sf;
    GRID_FLOAT_TYPE
    hs_bf, hs_bb, hs_fb, hs_ube, hs_ubb, hs_ueb;

    z_sf = z_s + z - z0;

    hs_bb = hs_ubb = hs_ube = INFINITY;
    for (x = x_start, y = y_start; x < x_end; x++) {
        hs_bf = hs[x][y][z_s];
        if (y) hs_bb = hs[x][y - 1][z_s];
        if (z_sf >= 0 && z_sf < nmesh_z) {
            hs_ube = hs[x][y][z_sf];
            if (y) hs_ubb = hs[x][y - 1][z_sf];
        }
        alert1 = t_1d(x + 1, y, z, t[x][y][z], hs_bb, hs_bf, hs_ubb, hs_ube);
        alert0 = t_2d(x + 1, y, z, t[x][y][z0], t[x][y][z], hs_bb, hs_bf);
        if (alert1) longflags[x * ny + ny + y] = 1;
        if (alert0) longflags[x * ny + ny + y] = 0;
    }

    hs_bb = hs_ubb = hs_ueb = INFINITY;
    for (x = x_start, y = y_start; y < y_end; y++) {
        hs_fb = hs[x][y][z_s];
        if (x) hs_bb = hs[x - 1][y][z_s];
        if (z_sf >= 0 && z_sf < nmesh_z) {
            hs_ueb = hs[x][y][z_sf];
            if (x) hs_ubb = hs[x - 1][y][z_sf];
        }
        alert1 = t_1d(x, y + 1, z, t[x][y][z], hs_bb, hs_fb, hs_ubb, hs_ueb);
        alert0 = t_2d(x, y + 1, z, t[x][y][z0], t[x][y][z], hs_bb, hs_fb);
        if (alert1) longflags[x * ny + y + 1] = 1;
        if (alert0) longflags[x * ny + y + 1] = 0;
    }

    for (x = x_start; x < x_end; x++) {
        for (y = y_start; y < y_end; y++) {

            hs_bb = hs[x][y][z_s];
            hs_bf = hs[x][y + 1][z_s];
            hs_fb = hs[x + 1][y][z_s];
            if (z_sf >= 0 && z_sf < nmesh_z) {
                hs_ubb = hs[x][y][z_sf];
                hs_ube = hs[x][y + 1][z_sf];
                hs_ueb = hs[x + 1][y][z_sf];
            } else hs_ubb = hs_ube = hs_ueb = INFINITY;

            /* bulk waves: 1 3D edge diffraction and 2 (*4) 3D transmission */
            alert0 = edge_diff(x + 1, y + 1, z, t[x][y][z0], t[x][y][z], hs_bb)
                    + t_3d(x + 1, y + 1, z, t[x][y][z0],
                    t[x + 1][y][z0], t[x][y][z], t[x + 1][y][z], hs_bb)
                    + t_3d(x + 1, y + 1, z, t[x][y][z0],
                    t[x][y + 1][z0], t[x][y][z], t[x][y + 1][z], hs_bb);
            if (alert0) {
                updated += alert0;
                longflags[x * ny + ny + y + 1] = 0;
            }

            /* interface waves along x_side: 1 1D transmission and 1 2D transmission */
            alert1 = t_1d(x + 1, y + 1, z, t[x][y + 1][z], hs_bb, hs_bf, hs_ubb, hs_ube);
            alert0 = t_2d(x + 1, y + 1, z, t[x][y + 1][z0], t[x][y + 1][z], hs_bb, hs_bf);
            if (alert0 + alert1) {
                updated += alert0 + alert1;
                flag_fb++; /* this scan must be (re-)examined */
                if (x_start_fb > x) x_start_fb = x;
                if (y_start_fb < y + 1) y_start_fb = y + 1;
                if (alert1) longflags[x * ny + ny + y + 1] = 1;
                else longflags[x * ny + ny + y + 1] = 0;
            }

            /* interface waves along y_side: 1 1D transmission and 1 2D transmission */
            alert1 = t_1d(x + 1, y + 1, z, t[x + 1][y][z], hs_bb, hs_fb, hs_ubb, hs_ueb);
            alert0 = t_2d(x + 1, y + 1, z, t[x + 1][y][z0], t[x + 1][y][z], hs_bb, hs_fb);
            if (alert0 + alert1) {
                updated += alert0 + alert1;
                flag_bf++; /* this scan must be (re-)examined */
                if (x_start_bf < x + 1) x_start_bf = x + 1;
                if (y_start_bf > y) y_start_bf = y;
                if (alert1) longflags[x * ny + ny + y + 1] = 1;
                else longflags[x * ny + ny + y + 1] = 0;
            }

            /* interface waves along z_side : 2 2D transmission and 1 2D diffraction */
            alert1 = diff_2d(x + 1, y + 1, z, t[x][y][z], hs_bb, hs_ubb)
                    + t_2d(x + 1, y + 1, z, t[x][y][z], t[x + 1][y][z], hs_bb, hs_ubb)
                    + t_2d(x + 1, y + 1, z, t[x][y][z], t[x][y + 1][z], hs_bb, hs_ubb);
            if (alert1) {
                updated += alert1;
                longflags[x * ny + ny + y + 1] = 1;
            }
        }
    }

    flag_ff = 0;
    x_start_ff = x_end;
    y_start_ff = y_end;

    return (updated);
}

/*--------------------------------------Z_SIDE() : SCAN_Z_BE()--------------*/

static int
scan_z_bf(x_begin, x_start, y_start, y_end, z0, z, z_s)
int x_start, x_begin, y_start, y_end, z0, z, z_s;

{
    int
    updated = 0,
            alert0, alert1,
            x, y,
            z_sf;
    GRID_FLOAT_TYPE
    hs_ff, hs_bb, hs_fb, hs_uee, hs_ubb, hs_ueb;

    z_sf = z_s + z - z0;

    hs_fb = hs_uee = hs_ueb = INFINITY;
    for (x = x_start, y = y_start; x > x_begin; x--) {
        hs_ff = hs[x - 1][y][z_s];
        if (y) hs_fb = hs[x - 1][y - 1][z_s];
        if (z_sf >= 0 && z_sf < nmesh_z) {
            hs_uee = hs[x - 1][y][z_sf];
            if (y) hs_ueb = hs[x - 1][y - 1][z_sf];
        }
        alert1 = t_1d(x - 1, y, z, t[x][y][z], hs_fb, hs_ff, hs_ueb, hs_uee);
        alert0 = t_2d(x - 1, y, z, t[x][y][z0], t[x][y][z], hs_fb, hs_ff);
        if (alert1) longflags[x * ny - ny + y] = 1;
        if (alert0) longflags[x * ny - ny + y] = 0;
    }

    hs_bb = hs_ubb = hs_ueb = INFINITY;
    for (x = x_start, y = y_start; y < y_end; y++) {
        if (x) hs_bb = hs[x - 1][y][z_s];
        hs_fb = hs[x][y][z_s];
        if (z_sf >= 0 && z_sf < nmesh_z) {
            if (x) hs_ubb = hs[x - 1][y][z_sf];
            hs_ueb = hs[x][y][z_sf];
        }
        alert1 = t_1d(x, y + 1, z, t[x][y][z], hs_bb, hs_fb, hs_ubb, hs_ueb);
        alert0 = t_2d(x, y + 1, z, t[x][y][z0], t[x][y][z], hs_bb, hs_fb);
        if (alert1) longflags[x * ny + y + 1] = 1;
        if (alert0) longflags[x * ny + y + 1] = 0;
    }

    for (x = x_start; x > x_begin; x--) {
        for (y = y_start; y < y_end; y++) {

            hs_ff = hs[x - 1][y + 1][z_s];
            if (x > 1) hs_bb = hs[x - 2][y][z_s];
            else hs_bb = INFINITY;
            hs_fb = hs[x - 1][y][z_s];
            if (z_sf >= 0 && z_sf < nmesh_z) {
                hs_uee = hs[x - 1][y + 1][z_sf];
                if (x > 1) hs_ubb = hs[x - 2][y][z_sf];
                else hs_ubb = INFINITY;
                hs_ueb = hs[x - 1][y][z_sf];
            } else hs_ubb = hs_uee = hs_ueb = INFINITY;

            /* bulk waves: 1 3D edge diffraction and 2 (*4) 3D transmission */
            alert0 = edge_diff(x - 1, y + 1, z, t[x][y][z0], t[x][y][z], hs_fb)
                    + t_3d(x - 1, y + 1, z, t[x][y][z0],
                    t[x - 1][y][z0], t[x][y][z], t[x - 1][y][z], hs_fb)
                    + t_3d(x - 1, y + 1, z, t[x][y][z0],
                    t[x][y + 1][z0], t[x][y][z], t[x][y + 1][z], hs_fb);
            if (alert0) {
                updated += alert0;
                longflags[x * ny - ny + y + 1] = 0;
            }

            /* interface waves along x_side: 1 1D transmission and 1 2D transmission */
            alert1 = t_1d(x - 1, y + 1, z, t[x][y + 1][z], hs_ff, hs_fb, hs_uee, hs_ueb);
            alert0 = t_2d(x - 1, y + 1, z, t[x][y + 1][z0], t[x][y + 1][z], hs_ff, hs_fb);
            if (alert0 + alert1) {
                updated += alert0 + alert1;
                flag_bb++; /* this scan must be (re-)examined */
                if (x_start_bb < x) x_start_bb = x;
                if (y_start_bb < y + 1) y_start_bb = y + 1;
                if (alert1) longflags[x * ny - ny + y + 1] = 1;
                else longflags[x * ny - ny + y + 1] = 0;
            }

            /* interface waves along y_side: 1 1D transmission and 1 2D transmission */
            alert1 = t_1d(x - 1, y + 1, z, t[x - 1][y][z], hs_bb, hs_fb, hs_ubb, hs_ueb);
            alert0 = t_2d(x - 1, y + 1, z, t[x - 1][y][z0], t[x - 1][y][z], hs_bb, hs_fb);
            if (alert0 + alert1) {
                updated += alert0 + alert1;
                flag_ff++; /* this scan must be (re-)examined */
                if (x_start_ff > x - 1) x_start_ff = x - 1;
                if (y_start_ff > y) y_start_ff = y;
                if (alert1) longflags[x * ny - ny + y + 1] = 1;
                else longflags[x * ny - ny + y + 1] = 0;
            }

            /* interface waves along z_side : 2 2D transmission and 1 2D diffraction */
            alert1 = diff_2d(x - 1, y + 1, z, t[x][y][z], hs_fb, hs_ueb)
                    + t_2d(x - 1, y + 1, z, t[x][y][z], t[x - 1][y][z], hs_fb, hs_ueb)
                    + t_2d(x - 1, y + 1, z, t[x][y][z], t[x][y + 1][z], hs_fb, hs_ueb);
            if (alert1) {
                updated += alert1;
                longflags[x * ny - ny + y + 1] = 1;
            }
        }
    }

    flag_bf = 0;
    x_start_bf = x_begin;
    y_start_bf = y_end;

    return (updated);
}

/*--------------------------------------Z_SIDE() : SCAN_Z_EB()--------------*/

static int
scan_z_bb(x_begin, x_start, y_begin, y_start, z0, z, z_s)
int x_start, x_begin, y_begin, y_start, z0, z, z_s;

{
    int
    updated = 0,
            alert0, alert1,
            x, y,
            z_sf;
    GRID_FLOAT_TYPE
    hs_ff, hs_bf, hs_fb, hs_uee, hs_ube, hs_ueb;

    z_sf = z_s + z - z0;

    hs_ff = hs_uee = hs_ueb = INFINITY;
    for (x = x_start, y = y_start; x > x_begin; x--) {
        if (y) hs_ff = hs[x - 1][y - 1][z_s];
        hs_fb = hs[x - 1][y][z_s];
        if (z_sf >= 0 && z_sf < nmesh_z) {
            if (y) hs_uee = hs[x - 1][y - 1][z_sf];
            hs_ueb = hs[x - 1][y][z_sf];
        }
        alert1 = t_1d(x - 1, y, z, t[x][y][z], hs_fb, hs_ff, hs_ueb, hs_uee);
        alert0 = t_2d(x - 1, y, z, t[x][y][z0], t[x][y][z], hs_fb, hs_ff);
        if (alert1) longflags[x * ny - ny + y] = 1;
        if (alert0) longflags[x * ny - ny + y] = 0;
    }

    hs_ff = hs_uee = hs_ube = INFINITY;
    for (x = x_start, y = y_start; y > y_begin; y--) {
        if (x) hs_ff = hs[x - 1][y - 1][z_s];
        hs_bf = hs[x][y - 1][z_s];
        if (z_sf >= 0 && z_sf < nmesh_z) {
            if (x) hs_uee = hs[x - 1][y - 1][z_sf];
            hs_ube = hs[x][y - 1][z_sf];
        }
        alert1 = t_1d(x, y - 1, z, t[x][y][z], hs_bf, hs_ff, hs_ube, hs_uee);
        alert0 = t_2d(x, y - 1, z, t[x][y][z0], t[x][y][z], hs_bf, hs_ff);
        if (alert1) longflags[x * ny + y - 1] = 1;
        if (alert0) longflags[x * ny + y - 1] = 0;
    }

    for (x = x_start; x > x_begin; x--) {
        for (y = y_start; y > y_begin; y--) {

            hs_ff = hs[x - 1][y - 1][z_s];
            if (x > 1) hs_bf = hs[x - 2][y - 1][z_s];
            else hs_bf = INFINITY;
            if (y > 1) hs_fb = hs[x - 1][y - 2][z_s];
            else hs_fb = INFINITY;
            if (z_sf >= 0 && z_sf < nmesh_z) {
                hs_uee = hs[x - 1][y - 1][z_sf];
                if (x > 1) hs_ube = hs[x - 2][y - 1][z_sf];
                else hs_ube = INFINITY;
                if (y > 1) hs_ueb = hs[x - 1][y - 2][z_sf];
                else hs_ueb = INFINITY;
            } else hs_uee = hs_ube = hs_ueb = INFINITY;

            /* bulk waves: 1 3D edge diffraction and 2 (*4) 3D transmission */
            alert0 = edge_diff(x - 1, y - 1, z, t[x][y][z0], t[x][y][z], hs_ff)
                    + t_3d(x - 1, y - 1, z, t[x][y][z0],
                    t[x - 1][y][z0], t[x][y][z], t[x - 1][y][z], hs_ff)
                    + t_3d(x - 1, y - 1, z, t[x][y][z0],
                    t[x][y - 1][z0], t[x][y][z], t[x][y - 1][z], hs_ff);
            if (alert0) {
                updated += alert0;
                longflags[x * ny - ny + y - 1] = 0;
            }

            /* interface waves along x_side: 1 1D transmission and 1 2D transmission */
            alert1 = t_1d(x - 1, y - 1, z, t[x][y - 1][z], hs_ff, hs_fb, hs_uee, hs_ueb);
            alert0 = t_2d(x - 1, y - 1, z, t[x][y - 1][z0], t[x][y - 1][z], hs_ff, hs_fb);
            if (alert0 + alert1) {
                updated += alert0 + alert1;
                flag_bf++; /* this scan must be (re-)examined */
                if (x_start_bf < x) x_start_bf = x;
                if (y_start_bf > y - 1) y_start_bf = y - 1;
                if (alert1) longflags[x * ny - ny + y - 1] = 1;
                else longflags[x * ny - ny + y - 1] = 0;
            }

            /* interface waves along y_side: 1 1D transmission and 1 2D transmission */
            alert1 = t_1d(x - 1, y - 1, z, t[x - 1][y][z], hs_ff, hs_bf, hs_uee, hs_ube);
            alert0 = t_2d(x - 1, y - 1, z, t[x - 1][y][z0], t[x - 1][y][z], hs_ff, hs_bf);
            if (alert0 + alert1) {
                updated += alert0 + alert1;
                flag_fb++; /* this scan must be (re-)examined */
                if (x_start_fb > x - 1) x_start_fb = x - 1;
                if (y_start_fb < y) y_start_fb = y;
                if (alert1) longflags[x * ny - ny + y - 1] = 1;
                else longflags[x * ny - ny + y - 1] = 0;
            }

            /* interface waves along z_side : 2 2D transmission and 1 2D diffraction */
            alert1 = diff_2d(x - 1, y - 1, z, t[x][y][z], hs_ff, hs_uee)
                    + t_2d(x - 1, y - 1, z, t[x][y][z], t[x - 1][y][z], hs_ff, hs_uee)
                    + t_2d(x - 1, y - 1, z, t[x][y][z], t[x][y - 1][z], hs_ff, hs_uee);
            if (alert1) {
                updated += alert1;
                longflags[x * ny - ny + y - 1] = 1;
            }
        }
    }

    flag_bb = 0;
    x_start_bb = x_begin;
    y_start_bb = y_begin;

    return (updated);
}

/*--------------------------------------Z_SIDE() : SCAN_Z_EB()--------------*/

static int
scan_z_fb(x_start, x_end, y_begin, y_start, z0, z, z_s)
int x_start, x_end, y_begin, y_start, z0, z, z_s;

{
    int
    updated = 0,
            alert0, alert1,
            x, y,
            z_sf;
    GRID_FLOAT_TYPE
    hs_ff, hs_bb, hs_bf, hs_uee, hs_ubb, hs_ube;

    z_sf = z_s + z - z0;

    hs_bf = hs_ubb = hs_ube = INFINITY;
    for (x = x_start, y = y_start; x < x_end; x++) {
        if (y) hs_bf = hs[x][y - 1][z_s];
        hs_bb = hs[x][y][z_s];
        if (z_sf >= 0 && z_sf < nmesh_z) {
            if (y) hs_ube = hs[x][y - 1][z_sf];
            hs_ubb = hs[x][y][z_sf];
        }
        alert1 = t_1d(x + 1, y, z, t[x][y][z], hs_bf, hs_bb, hs_ube, hs_ubb);
        alert0 = t_2d(x + 1, y, z, t[x][y][z0], t[x][y][z], hs_bf, hs_bb);
        if (alert1) longflags[x * ny + ny + y] = 1;
        if (alert0) longflags[x * ny + ny + y] = 0;
    }

    hs_ff = hs_uee = hs_ube = INFINITY;
    for (x = x_start, y = y_start; y > y_begin; y--) {
        if (x) hs_ff = hs[x - 1][y - 1][z_s];
        hs_bf = hs[x][y - 1][z_s];
        if (z_sf >= 0 && z_sf < nmesh_z) {
            if (x) hs_uee = hs[x - 1][y - 1][z_sf];
            hs_ube = hs[x][y - 1][z_sf];
        }
        alert1 = t_1d(x, y - 1, z, t[x][y][z], hs_bf, hs_ff, hs_ube, hs_uee);
        alert0 = t_2d(x, y - 1, z, t[x][y][z0], t[x][y][z], hs_bf, hs_ff);
        if (alert1) longflags[x * ny + y - 1] = 1;
        if (alert0) longflags[x * ny + y - 1] = 0;
    }

    for (x = x_start; x < x_end; x++) {
        for (y = y_start; y > y_begin; y--) {

            hs_bf = hs[x][y - 1][z_s];
            if (y > 1) hs_bb = hs[x][y - 2][z_s];
            else hs_bb = INFINITY;
            hs_ff = hs[x + 1][y - 1][z_s];
            if (z_sf >= 0 && z_sf < nmesh_z) {
                hs_ube = hs[x][y - 1][z_sf];
                if (y > 1) hs_ubb = hs[x][y - 2][z_sf];
                else hs_ubb = INFINITY;
                hs_uee = hs[x + 1][y - 1][z_sf];
            } else hs_ubb = hs_ube = hs_uee = INFINITY;

            /* bulk waves: 1 3D edge diffraction and 2 (*4) 3D transmission */
            alert0 = edge_diff(x + 1, y - 1, z, t[x][y][z0], t[x][y][z], hs_bf)
                    + t_3d(x + 1, y - 1, z, t[x][y][z0],
                    t[x + 1][y][z0], t[x][y][z], t[x + 1][y][z], hs_bf)
                    + t_3d(x + 1, y - 1, z, t[x][y][z0],
                    t[x][y - 1][z0], t[x][y][z], t[x][y - 1][z], hs_bf);
            if (alert0) {
                updated += alert0;
                longflags[x * ny + ny + y - 1] = 0;
            }

            /* interface waves along x_side: 1 1D transmission and 1 2D transmission */
            alert1 = t_1d(x + 1, y - 1, z, t[x][y - 1][z], hs_bb, hs_bf, hs_ubb, hs_ube);
            alert0 = t_2d(x + 1, y - 1, z, t[x][y - 1][z0], t[x][y - 1][z], hs_bb, hs_bf);
            if (alert0 + alert1) {
                updated += alert0 + alert1;
                flag_ff++; /* this scan must be (re-)examined */
                if (x_start_ff > x) x_start_ff = x;
                if (y_start_ff > y - 1) y_start_ff = y - 1;
                if (alert1) longflags[x * ny + ny + y - 1] = 1;
                else longflags[x * ny + ny + y - 1] = 0;
            }

            /* interface waves along y_side: 1 1D transmission and 1 2D transmission */
            alert1 = t_1d(x + 1, y - 1, z, t[x + 1][y][z], hs_ff, hs_bf, hs_uee, hs_ube);
            alert0 = t_2d(x + 1, y - 1, z, t[x + 1][y][z0], t[x + 1][y][z], hs_ff, hs_bf);
            if (alert0 + alert1) {
                updated += alert0 + alert1;
                flag_bb++; /* this scan must be (re-)examined */
                if (x_start_bb < x + 1) x_start_bb = x + 1;
                if (y_start_bb < y) y_start_bb = y;
                if (alert1) longflags[x * ny + ny + y - 1] = 1;
                else longflags[x * ny + ny + y - 1] = 0;
            }

            /* interface waves along z_side : 2 2D transmission and 1 2D diffraction */
            alert1 = diff_2d(x + 1, y - 1, z, t[x][y][z], hs_bf, hs_ube)
                    + t_2d(x + 1, y - 1, z, t[x][y][z], t[x + 1][y][z], hs_bf, hs_ube)
                    + t_2d(x + 1, y - 1, z, t[x][y][z], t[x][y - 1][z], hs_bf, hs_ube);
            if (alert1) {
                updated += alert1;
                longflags[x * ny + ny + y - 1] = 1;
            }
        }
    }

    flag_fb = 0;
    x_start_fb = x_end;
    y_start_fb = y_begin;

    return (updated);
}
/*--------------------------------------------END_Z_SIDE()--------------------*/
/*********************************************END_TIME_3D**********************/

