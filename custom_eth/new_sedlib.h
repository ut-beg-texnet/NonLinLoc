/*************************
   AnyPicks.h 
 *************************/
#include <dirent.h>
#include <libgen.h>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>

#ifdef MAIN
#       define GLOBAL
#else
#       define GLOBAL extern
#endif

GLOBAL char sedlibPath[256];

typedef struct {
        char  stname[8];    
        int   minute;		/* 4	4 */
   	float time;		/* 4	8 */
	float tmin;
	float tmax;
   	char  phase[8];		/* 8	16*/
	char  quali[4];		/* 2	18*/
	int   weight;		/* 2	20*/
	int   dmy;		/* 4    24*/
	int   amplitude;	/*20	44*/
	char  ampchan[8];
	char  amptype[8];
	float periode;
	float amptime;
	int   ampuse;
	char  recsys[8];
	float clipp_threshold;
	float dev_fac;
} NewPick;

GLOBAL NewPick *nmxTrigger, *autoPicks, *manualPicks;
GLOBAL int    NnmxTrigger, NautoPicks, NmanualPicks;
GLOBAL int    *sort_idx;
GLOBAL float  *sort_time;

typedef struct {
      int 	n_picks;	
      int  	st_list;	
      int	model;		
      float	lat;		
      float	lon;		
      float	depth;		
      float 	xnear;		
      float	xfar;		
      float	pos;		
      float	dray;		
      float	daz;		
      char 	type[8];	
      int	fix;
      int	yr;
      int	mo;
      int	dy;
      int	hr;
      int	mn;
      float	sec;		
      int	swt;		
      char      AmpFilt;	
      char	NoMag;		
      int	dumy2;
} NewLocPar;

GLOBAL NewLocPar locNmxT, locAutoP, locManualP;
GLOBAL char author[50];


/***************************
 Calibration
 ***************************/
typedef struct{
   char		units[2];
   float	srate;
   float	counts_volt;
} Digitizer;
typedef struct{
   char		units[2];
   float	factor;
   float 	*poles;
   float	*zeros;
   int		npoles;
   int		nzeros;
} PAZ_str;
typedef struct{
   char		units[2];
   float 	factor;
   float 	*coeff;	
   int		ncoeff;
   int		decim;
} FIR_str;
typedef struct {
   Digitizer	*digit;
   PAZ_str	*paz;
   FIR_str	*fir;
} Stage_str;
typedef struct {
   char 	name[8];
   char		comp[8];
   char		auxid[8];
   char         recsys[8];
   float	nmct;
   float	calper;
   float	srate;
   int  	time1;
   int		time2;
   int		resolution;
   float	volts;
   float	sconst;
   float  	periode;
   float	damping;
   float	gain;
   Stage_str	stage[40];
   int 		Nstage;
} CalEntry;

GLOBAL CalEntry *calib;
GLOBAL int ncalib;



/***********************
 file descritors 
 ***********************/

#define MAX_STAT	256

#ifndef True
#define True	1
#define False	0
#endif

#define TRIGGER		0
#define APICK		1
#define MPICK		2
#define AFILTER		3
#define MFILTER		4
#define ALOC		5
#define APDE		6
#define MLOC		7
#define MPDE		8
#define MAGNITUDE	9
#define TIMESTATUS	10
#define SIGNALSTATUS	11
#define WAVEFORM	12
#define NOEXT		13

typedef struct {
   char extension[20];
   int  limit;
} ExtLimit;

GLOBAL ExtLimit fext[13];
GLOBAL int nfext;

typedef int BoolSED;
typedef char c8[8];
typedef char c4[4];
typedef char FSize[256];

GLOBAL char workingfile[256], cleanName[256], extension[256];
GLOBAL BoolSED is_compressed;
GLOBAL int maxseismo;
GLOBAL c8 channel_to_name[MAX_STAT];
GLOBAL int channel_to_station[MAX_STAT], channel_to_compnr[MAX_STAT];
GLOBAL BoolSED IsGSEFormat, UpdatePermission;
GLOBAL int NrEvent;

/**************************
 Filter
 **************************/

typedef struct {
   	char	typ[3];
	char	zerophase;
   	float	low;
	float	high;
	int	order;
} NewFilter;

   
GLOBAL  NewFilter   *autofil, *manufil;
GLOBAL  int   Nautofil, Nmanufil;



/*==================================
  locations:
  =================================*/

typedef struct {
   char locsum[81];
   char region[81];
} NewLoc;
typedef struct {
	float fix_ml;
	float fix_mb;
	float fix_ms;
	char  fix_ms_str[8];
} FixMag;
GLOBAL NewLoc autoLoc, manuLoc;
GLOBAL FixMag mag_fixed;

/***********************************************************
 * new signal header structure for seismic event files:
 ***********************************************************/

typedef struct {
   char		s_name[8];		/*	station name		*/
   char 	s_comp[8];		/*	station component	*/
   char 	s_typ[8];		/*	seismometer type	*/
   char 	s_recsys[8];		/*	seismometer type	*/
   char 	s_compress[8];		/*	data compression	*/
   float 	s_lat;			/*	station latitude	*/
   float	s_lon;			/*	station longitude	*/
   float	s_elev;			/*	station elevation	*/
   float	s_burrial;		/*	station burrial depth	*/
   long   	s_start_min;		/*	start time minutes	*/
   float      	s_start_sec;		/*	start time seconds	*/
   long   	s_end_min;		/*	end time minutes	*/
   float      	s_end_sec;		/*	end time seconds	*/
   float      	s_freq;			/*	sampling frequency	*/
   long   	s_nsampl;		/*	number of samples	*/
   int   	s_resol;		/*	resolution		*/
   int     	s_position;		/*	data position 		*/
   int  	s_nbytes;		/*	bytes per data word	*/
   int          s_nblocks;		/*	number of blocks (512 b)*/
   int     	s_status;		/*	status bit		*/
   int     	s_ndiff;		/*	differences taken	*/
   long   	s_n_compress;		/*	nr of character compressed*/
   int  	s_time_base;		/* 	time base channel reference */
   int     	s_gps_sync;		/*	time syncronization (gps status) */
   float	s_volts;		/*	volts at max dig value	*/
   float	s_gain;			/*	voltage gain		*/
   float	s_damp;			/* 	damping			*/
   float	s_const;		/* 	seimometer constant	*/
   float	s_period;		/*	seismometer eigenperiod	*/
   long   	s_checksm;		/*	signal checksum		*/
   float  	s_delay;		/*	time delay of signal (sec) */
   int          s_equip;                /*      equipment id            */
   int          s_z_sync;               /*      time sync 0=ok; 1= no time sync */
   int          s_z_typ;                /*      time base equippment (not used)*/ 
   float	s_saturation;		/*	saturation based on resolution */  
   CalEntry     *calib;			/*      Calibration section */
   int		*status_out;		/*	sample number of signal-status out */
   int		*status_back;		/*      sample number when signal-status back */
   int          Nstatus_out;
   float	*timing_out;		/*	seconds after s_start_min/s_start_sec that time-sync is lost */
   float	*timing_back;		/*	seconds after s_start_min/s_start_sec that time-sync is back */   
   int          Ntiming_out;
} SignalHead;

/*=========================
  New Eventfile-signals
 ==========================*/
typedef struct {
        char     name[8];             /*station name */       
	SignalHead shead[32];    /*signal header */
	char     *wave[32];
	c4	 chanid[32];
	int	 nsig;
	int	 available[32];	/*selected channels */
	int	 trace_nr[32];	/*trace number in input file */
	int	 comp_tab[32];  /*component tabel entry */
	int	 selected;
} Stat_Str;

GLOBAL Stat_Str station[MAX_STAT];
GLOBAL int nstat;

GLOBAL int ncomp_table;
GLOBAL c4 comp_table[32];

GLOBAL int time_minimum, time_maximum;
GLOBAL float sec_minimum, sec_maximum;
GLOBAL int DoStatusCheck;
GLOBAL float s_freq_max;
GLOBAL int MaxSample;

typedef struct {
   char		name[8];
   char		comp[4];
   float	lat;
   float	lon;
   float	elev;
   char		ondate[20];
   char		ofdate[20];
   char		nation[4];
   int		v_model;
   float	MOHO_depth;
   char		bulletin;
   float	coeff[3];
} GseStatList;

GLOBAL GseStatList *gse_stn;
GLOBAL int  gse_nst;

typedef struct {
   char recsys[8];
   int	resol;
   float volts;
   float clip;
} RecSysDefinition;

GLOBAL RecSysDefinition *recSys;
GLOBAL int NrecSys;

GLOBAL char line[256], str1[256], str2[256], str3[256], str4[256];
/*=========================
  KP- Filter.h 
 ==========================*/

typedef struct {
   	char	typ[3];
	char	zerophase;
   	float	low;
	float	high;
	int	order;
} FilterStruct;

typedef struct{
      char      f_typ[4];
      float     f_corner;
} one_filter;
typedef struct{
      one_filter        auto_f[32];
      FilterStruct      manu_f[16];
} filt_rec;
   
GLOBAL filt_rec   frec;

/*==================================
  KP- location records:
  =================================*/
typedef struct {
	char auto_loc[80];
	char auto_reg[80];
	short auto_fill[48];
	char manu_loc[80];
	char manu_reg[80];
	float fix_ml;
	float fix_mb;
	float fix_ms;
	char  fix_ms_str[4];
	short manu_fill[40];
} kpLocRec;

GLOBAL kpLocRec location;

/*=====================================================
 * KP master structure 
 ======================================================*/

typedef unsigned short R256[256];	
typedef struct {
	long 	         starttime;		 /*   0 */
	float            startsec;		 /*   2 */
	long 	         endtime;		 /*   4 */
	float            endsec;		 /*   6 */
	unsigned short  lastrec;		/*   7 */
	char	         timecode[2];		 /*   8 */
	unsigned short 	nseismo;		/*   9 */
	unsigned short 	whereseismo;		/*  10 */
	unsigned short 	station_rec;		/*  11 */
	unsigned short 	trigger_rec;		/*  12 */
	unsigned short 	filter_rec;		/*  13 */
	unsigned short 	a_pick_rec;		/*  14 */
	unsigned short 	m_pick_rec;		/*  15 */
	unsigned short 	loc_rec;		/*  16 */
	unsigned short  pick_type;		/*  17 */
	unsigned short  ev_nr;			/*  18 */
	unsigned short  n_pick_rec;		/*  19 */
	unsigned short  nmx_tr_rec;		/*  20 */
	unsigned short 	m_dummy[233];		/*  21 */
	unsigned short 	next_master;
} Master;

GLOBAL R256 *mrec;
GLOBAL int Nmaster;
/*=====================================================
   KP new pick records:
   definition of the manual pick record for the seismic
   data files:
   each record consists of:
    	10 records of     48 bytes length
	 1 dummy array of 30   "      "
	 1 pointer to next record (2 bytes)
   the first record of the first block contains
   the location parameters according to the
   structure 'loc_par' below.
   all the following records contain the manual picks
   according to the structure 'one_mpick'.
   
   typedef struct {
   	int   amp;		* 4	4 *
        char  achan[4];		* 4	8 *
   	float period;		* 4	12*
	float time;		* 4	16*
	short use		* 2     18*
	short dmy		* 2     20*
	} Amplitude;
   typedef struct {
        int   minute;		* 4	4 *
   	float time;		* 4	8 *
   	char  phase[8];		* 8	16*
	char  quali[2];		* 2	18*
	short weight;		* 2	20*
	int   dmy;		* 4     24*
	Amplitude amplitude;	*20	44*
	} Pick;

=====================================================*/
   typedef struct {
        int   minute;		/* 4	4 */
   	float time;		/* 4	8 */
   	char  phase[8];		/* 8	16*/
	char  quali[2];		/* 2	18*/
	short weight;		/* 2	20*/
	int   dmy;		/* 4    24*/
   	int   amplitude;	/* 4	4 */
        char  ampchan[4];	/* 4	8 */
   	float period;		/* 4	12*/
	float amptime;		/* 4	16*/
	short ampuse;		/* 2    18*/
	char  dmy2[2];		/* 4    20*/
   } Pick;
	
   /* one pick record: */
   typedef struct {
      short 	chan;		/* 2	2  */
      short	dummy;		/* 2	4  (for alignment reasons) */
      Pick	m_phase;	/* 44   48 */
   } one_mpick;
   
   /* location parameters: */
   typedef struct{
      short 	n_picks;	/* 2 	2  */
      char  	st_list;	/* 1 	3  */
      char	model;		/* 1 	4  */
      float	lat;		/* 4 	8 */
      float	lon;		/* 4 	12 */
      float	depth;		/* 4 	16 */
      float 	xnear;		/* 4	20 */
      float	xfar;		/* 4 	24 */
      float	pos;		/* 4 	28 */
      float	dray;		/* 4 	32 */
      float	daz;		/* 4 	36 */
      char 	type[4];	/* 4 	40 */	
      char	fix;		/* 1 	41 */
      char	swt;		/* 1 	42 */
      char      AmpFilt;	/* 1 	43 */
      char	NoMag;		/* 1	44 */
      int	dumy2;		/* 4	48 */
   } loc_par;
   
   /* make the union:
      the first pick record will have the location parameters
      as the first 48 bytes, all other pick blocks 
      will contain only pick structures: */
   union par_pick{
      loc_par	par;
      one_mpick pick;
   };
   
   /* pick record in data file: */
   typedef struct{
      union par_pick  p[10];            /* 10*48 480 */
      short	      dumy3[15];	/* 30	 510 */      
      unsigned short p_next;		/* 2	 512 */
   } mpick_rec;
/*==========================================================
 * KP-signal header structure for seismic event files:
 ===========================================================*/

typedef struct {
   char		s_name[6];		/*	station name		*/
   char 	s_comp[2];		/*	station component	*/
   char 	s_typ[6];		/*	seismometer type	*/
   short   	s_1;			/*	dummy 1			*/
   float 	s_lat;			/*	station latitude	*/
   float	s_lon;			/*	station longitude	*/
   float	s_elev;			/*	station elevation	*/
   float	s_burrial;		/*	station burrial depth	*/
   long   	s_start_min;		/*	start time minutes	*/
   float      	s_start_sec;		/*	start time seconds	*/
   long   	s_end_min;		/*	end time minutes	*/
   float      	s_end_sec;		/*	end time seconds	*/
   float      	s_freq;			/*	sampling frequency	*/
   long   	s_nsampl;		/*	number of samples	*/
   short   	s_resol;		/*	resolution		*/
   short   	s_position;		/*	data position 		*/
   short   	s_nbytes;		/*	bytes per data word	*/
   unsigned short s_nblocks;		/*	number of blocks (512 b)*/
   char 	s_compress[4];		/*	data compression	*/
   short   	s_status;		/*	status bit		*/
   short   	s_ndiff;		/*	differences taken	*/
   long   	s_n_compress;		/*	nr of character compressed*/
   short	s_time_base;		/* 	time base channel reference */
   short   	s_gps_sync;		/*	time syncronization (gps status) */
   float	s_volts;		/*	volts at max dig value	*/
   float	s_gain;			/*	voltage gain		*/
   float	s_damp;			/* 	damping			*/
   float	s_const;		/* 	seimometer constant	*/
   float	s_period;		/*	seismometer eigenperiod	*/
   long   	s_4;			/*	dummy 1			*/
   long   	s_checksm;		/*	signal checksum		*/
   float  	s_delay;		/*	time delay of signal (sec) */
   short   	s_wh_poles;		/*	where is pole represent.*/
   short        s_equip;                /*      equipment id            */
   short        s_z_sync;               /*      time sync 0=ok; 1= no time sync */
   short        s_z_typ;                /*      time base equippment (not used)*/   
   short   	s_nzero;		/*	number of zeros		*/
   short   	s_npoles;		/*	number of poles		*/
   float	s_factor;		/*	factor			*/
   float	s_zero_pole[95];	/*	zero and poles array	*/	
   short   	s_7;			/*	dummy 1			*/
   unsigned short s_next_rec;		/*	next record flag	*/
} kpSignal;
/*=====================================================
 * KP-station list structure for sed seismic data files
 =====================================================*/
typedef struct {
   	char stname[6];
   	char stcomp[2];
} one_stn;
typedef struct {
   	one_stn stat_name[63];
   	short	dummy[3];
	unsigned short	next_stat;
} one_stn_rec;

/*=====================================================
 * general runstring structure
 =====================================================*/

typedef char OptString[256];
typedef struct {
   char def[6];
   char set[300];
   } Opts;
GLOBAL Opts opt[50];
GLOBAL int nopt;

/* SH 07212004 function declarations in new_sedlib.c */
int    NewOpenFile();
BoolSED GetPermission();
int    ReadExtensions();
int    AddComp2Table();
int    ReadFullFile();
void   CeckCompressed();
void   LastExt();
void   PointerInitialize();
int    AddCalibToSignals();
int    AddRecSysToPicks();
int    AddRecSysToThesePicks();
int    AddRecSys();
int    TimeWindowOverLap();
int    DefaultLocPar();
int    DecodeGSEfile();
int    ReadMagnitudefromFile();
int    ReadTimeStatusFile();
int    ReadSignalStatusFile();
int    ReadPicksfromFile();
int    ReadFiltersfromFile();
int    ReadLocationfromFile();
int    ReadCalibrationfromFile();
int    ReadWIDfromFile();
int    ReadWIDoneChannel();
int    checkCompTable();
int    checkComponent();
int    checkStationinList();
int    ReadCALfromFile();
char  *Comp6toBuf();
// AJL 20040518
//void   Decomp6toArray();
void Decomp6toArray(long lb, long *lout, long *iout, char *cbuf);
// AJL
void   MakeDiff();
void   RemoveDiff();
int    WritePicksToFile();
char  *FixMagnitude();
int    MagnitudeToFile();
int    UpdateLocation();
int    UpdateLocationOutput();
int    SaveLocKPFile();
int    WriteLocToFile();
int    WriteFilterToFile();
int    WriteSignalStatusToFile();
int    ConvertToFloat();
int    ConvertToInt();
float  AsignalC();
int    sortPicks();
void   bubble_sort();
void   bubble_sortI();
int    DecodeKPfile();
void   DecodeKPTrigger();
int    ReadOldPicks();
int    ReadNewPicks();
void   DecodeKPaPicks();
void   DecodeKPmPicks();
void   DecodeKPamLoc();
void   DecodeKPFilter();
void   DecodeKPsignalHead();
void   TryRecordingSystem();
int    SignalHeadFromKP();
int    SaveKPFilter();
int    WritePicksToKPFile();
void   cparm();
int    ccheck_opt();
char  *cget_opt();
void   cset_opt();
int    juliam();
int    juliam4();
void   datum();
void   datum4();
void   build_name_chanid();
int    csplitstring();
int    remove_leading();
char  *nfs_mount();
void   trimlencc();
void   trimlen0c();
void   trimlenc();
void   trim();
void   expandString();
FSize  *FindFiles();
void   LimitFiles();
FILE   *OpenLimiter();
void   GetStationList();
void   GetCal4StnTime();

