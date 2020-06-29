#ifdef MAIN
#	define GLOBAL
#else
#	define GLOBAL extern
#endif

#define MAXPHASE	300
#include "gsetype.h"
//#include "group_velocity.h"

/* for iaspei travel times: */
typedef struct{
   float save_depth;
   char  save_phase[10];
} s_dep;

/* for iaspei data input */      
typedef struct{
   int nasgr;
   int file_data;
} rec_siz;
 
typedef struct {
   int 		nlay;
   int		nmoho;
   float	vp[10];
   float	vp2[10];
   float	vs[10];
   float	vs2[10];
   float	z[10];
   float	h[10];
   float	f[10];
   char		refrakp[10][8];
   char		reflecp[10][8];
   char		refraks[10][8];
   char		reflecs[10][8];
} VZmodel;

typedef struct {	/* station list structure */
   char name[6];
   float lat;
   float lon;
   float xlon;			/* swiss coordinates */
   float xlat;			/* swiss coordinates */
   float elevation;
   double dircos[3];
   float tcor[3];
   int 	model;
   VZmodel lmod;
   float mohodepth;
} StnRec;

typedef struct {	/* location structure */
   int	o_time_mn;
   float o_time_sc;
   float o_lat;
   float o_lon;
   float o_depth;
   float o_mag;
   char  o_mag_typ[3];
   char  o_event_type[5];
   int	p_time_mn;
   char  agency[7];
   int	 gap;
   int   ev_nr;
   int 	nstat;
   int	nphas;
   int	reg_nr;
   float rms;
   float err_0;
   float err_x;
   float err_y;
   float err_z;
   float err_mag;
   int 	model;
   int	stlist;
   float dmin;
   float xnear;
   float xfar;
   float pos;
   int   s_weight;
   int o_fix;
   char type;
} LRec;

#define PS_STR_LEN 9
typedef struct {	/* phase data structure */
   char name[6];
   char phase[PS_STR_LEN];
   char phr[PS_STR_LEN];
   char phi[PS_STR_LEN];
   char onset[3];
   int  use;
   float p_time_sc;
   float corr;
   float period;
   float amplitude;
   float weight;
   char  rec_sys[10];
   int	 clipped;
   float mag;
   char  magrmk;
   float delta;
   float azim;
   float angle;
   StnRec *tab;
   int   stntab_index;
   VZmodel lmod;
   float elevcor;
   float res;
   int	 nrph;
   char  amp_chan[5];      /* added SH 25.05.2006 */
} PRec;

typedef struct {		/* grid point structure */
   float lat;
   float lon;
   float depth;
   float xlon;			/* swiss coordinates */
   float xlat;			/* swiss coordinates */
   int nearest;			/* number of closest station to this gridpoint */
   double dircos[3];
   float rms;
   int o_time_mn;
   float o_time_sc;
   int  do_it;
   float azimuth;
   float delta;
   float t_time_sc[MAXPHASE]; 	/* travel times from the stations */
   float tcor[MAXPHASE]; 	/* azimuthal travel time correction to stations*/
   float resid[MAXPHASE];	/* travel time residuals */
   float dweight[MAXPHASE];	/* phase weights for distance*/
   float rweight[MAXPHASE];	/* phase weights for residuals*/
   int   do_out;                /* calculate outlyers: 1 yes; 0 no */
   int   outwt[MAXPHASE];	/* outlying phases =0; otherwise= 1*/
   float outvar;		/* best variance after outlayers */
   float outskew;		/* best skew     after outlayers */
   float outcurt_limit;		/* sqrt(96/n) */
   float outcurt;		/* best curtosis after outlayers */
   float outave;		/* average after outlayers */
   int   noutlay;		/* number of outlayers */
} GridStr;
   
typedef struct {
   int to;
   int do_it;
} FromTo;

typedef struct {
   FromTo gm[27];
} GridMove;

GLOBAL s_dep iaspei_parameters;
GLOBAL rec_siz iaspei_input;

GLOBAL FILE *fp_in;	/* Pointer to input  file   */
GLOBAL FILE *fp_out;	/* Pointer to output file   */
GLOBAL FILE *fp_err;	/* Pointer to error-field file   */
GLOBAL GridMove move[27];
// AJL 20040527
//GLOBAL LRec location[5];
GLOBAL PRec phases[MAXPHASE];
GLOBAL StnRec stn[1000];

/** h-r maurer **/
GLOBAL cal_struct *calinfo;
GLOBAL station_struct *stntab;
GLOBAL int nsttab, ncalib;

GLOBAL float WA_gain;
GLOBAL float twopi;

GLOBAL GridStr grid[27], *g[27], *gs[27];
GLOBAL VZmodel Model[2];
GLOBAL int nlocal_mod;
GLOBAL int check_out, KU_max, KU_modulo;
GLOBAL int nloc, nphas, nstat, early_stat, niter, maxiter,node_comp;
GLOBAL int travel_time_corr, no_corr, final, auto_mode;
// AJL 20040518
GLOBAL char author[50];
//GLOBAL char author[4];
//
GLOBAL char linp[100], fixms_string[10];
GLOBAL char *station_list_file,*local_traveltime_file,*tele_traveltime_file;
GLOBAL char sed_path[256], calib_file[256];
GLOBAL float dlat,dlon,ddepth,dlat_limit,dz_limit;
GLOBAL float lat_ch,lon_ch;
GLOBAL float ta[60], sina[60], pi2;
GLOBAL float warr[MAXPHASE];
GLOBAL int   idx[MAXPHASE];
GLOBAL char pname[60][8], event_string[25];
GLOBAL int nph, do_distance_weight, do_residual_weight, start_dw, start_rw;
GLOBAL float ieq_weight;
GLOBAL int minstat, minphas, terminate;
GLOBAL int fix_depth_search_tele;
GLOBAL int debug_output;
GLOBAL int nr_error_field;
GLOBAL float deltah_km_error_field,deltaz_km_error_field,pirad;
GLOBAL float fixml, fixmb, fixms;
GLOBAL int  do_fixml, do_fixmb, do_fixms;
GLOBAL int  center_node;

#ifdef MAIN
   int	fix_0[]  = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,
                     18,19,20,21,22,23,24,25,26};
   int  fix_1[]  = { 9,10,11,12,13,14,15,16,17};
   int  fix_2[]  = { 4,13,22 };
   int  fix_3[]  = { 13 };
   
   float dlatf[] = { 1.0, 1.0, 1.0,
                     0.0, 0.0, 0.0,
		    -1.0,-1.0,-1.0,
		    
                     1.0, 1.0, 1.0,
                     0.0, 0.0, 0.0,
		    -1.0,-1.0,-1.0,
		    
                     1.0, 1.0, 1.0,
                     0.0, 0.0, 0.0,
		    -1.0,-1.0,-1.0
                    };   
   float dlonf[] = {-1.0, 0.0, 1.0,
                    -1.0, 0.0, 1.0,
		    -1.0, 0.0, 1.0,
		    
		    -1.0, 0.0, 1.0,
                    -1.0, 0.0, 1.0,
		    -1.0, 0.0, 1.0,
		    
		    -1.0, 0.0, 1.0,
                    -1.0, 0.0, 1.0,
		    -1.0, 0.0, 1.0 
                    };
   float depthf[] ={-1.0,-1.0,-1.0,
		    -1.0,-1.0,-1.0,
		    -1.0,-1.0,-1.0,
		     
		     0.0, 0.0, 0.0,
		     0.0, 0.0, 0.0,
		     0.0, 0.0, 0.0,
		     
		     1.0, 1.0, 1.0,
		     1.0, 1.0, 1.0,
		     1.0, 1.0, 1.0 
                    };
#else
   GLOBAL int   fix_0[], fix_1[], fix_2[], fix_3[];
   GLOBAL float dlatf[], dlonf[], depthf[];
   GLOBAL char *VERSION;
#endif


