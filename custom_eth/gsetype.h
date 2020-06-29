// AJL 20040527
#include <stddef.h>
#include <stdio.h>

/*------------------------------------------------------------------------*/
/* Typedefinitions                                                        */
/*------------------------------------------------------------------------*/

typedef char str9[9];
typedef char str7[7];
typedef char string50[128];


typedef struct {
	short int   year;
	short int   month;
	short int   jday;
	short int   day;
	short int   hour;
	short int   min;
	long  int   jmin;
	short int   sec;
	short int   msec;
	float fsec;
}  GSEtime;

typedef struct {
	char wavetype;
	int statnr;
	int statnr_sed;
	int ncnr;
	GSEtime phtime;
	char pol[4];
	char pid;
	char pidstring[8];
	short int ampl;
	float per;
	short int weight;
} phi;

typedef struct {
        int used;
	char name[7];
	char chan[4];
	float lat;
	float lon;
	float elev;
	char seismotype[32];
	char instrtype[64];
	char auxid[8];
	int jminb;
	int jmine;
	char resptype[4];
	int npoles;
	int nzeros;
	float *realpole;
	float *imagpole;
	float *realzero;
	float *imagzero;
	float scalefactor;
	float seismoperiod;
	float seismodamp;
	float seismocoil;
	float volts;
	float gain;
	float resolution;
	float cpv;
	float sra;
	char units;
	float calper;
	float widg;
} cal_struct;

typedef struct {
	GSEtime  tfsamp;
	long     nosamp;
	char     name[7];
	char     chid[9];
	char     chan[4];
	float    sra;
	char     systype[7];
	char     form[5];
	short    dflag;
	float    gain;
	short    units;
	float    cperiod;
	float    lat;
	float    lon;
	float    alt;
	float    depsens;
	long     checksum;
	short    nblock;
	short    dummyval;
	long     cmpbytes;
	short    stat_nr;
	size_t   no_sec;
	int      use_flag;
	int      *trace_i4;
	char     *cbuf;
	cal_struct *cal_data;
} gseh;

typedef struct {
	char type[4];
	char stlist;
	char model;
	float xnear;
	float xfar;
	float trlat;
	float trlon;
	float trdepth;
	float pos;
	char  stopweight;
	char  locstring[82];
	float epilat;
	float epilong;
	float epidepth;
	GSEtime  epitime;
	float magn;
	char  magntype[3];
	short int   map;
	short int   xkm;
	short int   ykm;
	float xerr;
	float yerr;
	float zerr;
	short int flen;
	float raypar;
	float errpar;
	float dazim;
	short int azim;
	float delta;
	float rms;
	short int   gap;
	short int   dmin;
	short int   nrread;
	short int   nrofreadphases;
	phi   *phaselist;
	int   event_nr;
	str7  *new_stat;
} hyda;


typedef struct {
	double  azimuth;
	double  dist;
	double  rsecp;
	double  rsecs;
	long    firstsamp;
	long    lastsamp;
	long    tr_len;
} fsp;

typedef struct {
	double shotlat;
	double shotlon;
	double profaz;
	double  vp;
	double  vs;
	char   pname[19];
	GSEtime    shottime;
	double firstsec;
	double lastsec;
} fsm;

typedef struct {
	short nst;
	short nwf;
	short firstwf;
	short lastwf;
	str9  *wfindex;
	gseh  *gseinfo;
	hyda  hydainfo;
	fsp   *fsparam;
	fsm   fsmain;
	int   no_cal;
	short  **bads;
	int   no_use_counter;
} masterheader;


typedef struct {
   char name[7];
   char type[5];
   float lat;
   float lon;
   float elev;
   long on_date;
   long off_date;
   char country[4];
   char network[13];
   char full_name[17];
   int v_model;
   float MOHO_depth;
   char bull_ind[5];
   float az_corr[3];
   char reserved[17];
   int nchan;
   cal_struct *chan;
} station_struct;

typedef struct {
   char  auxid[8];
   int   npoles;
   int   nzeros;
   float scalefactor;
   float *realpole;
   float *imagpole;
   float *realzero;
   float *imagzero;
} SedAuxid;
/*------------------------------------------------------------------------*/
/* End ofTypedefinitions                                                  */
/*------------------------------------------------------------------------*/


station_struct *GetGseStationList();
void WriteBinStationList();


#ifdef Motif
#define malloc XtMalloc
#define realloc XtRealloc
#define free XtFree
#endif

