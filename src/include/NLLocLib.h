/*
 * Copyright (C) 1999-2010 Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
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


/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
161 All?e du Micocoulier, 06370 Mouans-Sartoux, France
tel: +33(0)493752502  e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


/*  NLLocLib.h

        include file for program NLLoc, NLDiffLoc

 */


/* the following should be included before NLLocLib.h in all sources that include NLLocLib.h
#include "GridLib.h"
#include "ran1.h"
#include "velmod.h"
#include "GridMemLib.h"
#include "calc_crust_corr.h"
#include "phaseloclist.h"
#include "otime_limit.h"
 * */


/* defines */




/*------------------------------------------------------------*/
/* structures */

/* gaussian errors location parameters */

/*	see (TV82, eq. 10-14; MEN92, eq. 22) */
typedef struct {
    double SigmaT; /* theoretical error coeff for travel times */
    double CorrLen; /* model corellation length */
    MatrixDouble EDTMtrx; /* EDT covariance (row=col) or correlation (row!=col) matrix */
    MatrixDouble WtMtrx; /* weight matrix */
    double WtMtrxSum; /* sum of elements of weight matrix */
    long double meanObs; /* weighted mean of obs arrival times */
    double meanPred; /* weighted mean of predicted travel times */
    double arrivalWeightMax; /* maximum station (row) weight */
}
GaussLocParams;

/* gaussian error travel time parameters */
typedef struct {
    double SigmaTfraction; /* fraction of travel times to use as travel times error */
    double SigmaTmin; /* minimum travel times error */
    double SigmaTmax; /* maximum travel times error */
}
Gauss2LocParams;

/* station/inst/component parameters */
typedef struct {
    char label[ARRIVAL_LABEL_LEN]; /* label (i.e. station name) */
    char inst[INST_LABEL_LEN]; /* instrument */
    char comp[COMP_LABEL_LEN]; /* component */
    double amp_fact_ml_hb; /* amplitude scale factor */
    double sta_corr_ml_hb; /* station correction */
    double sta_corr_md_fmag; /* station correction MD_FMAG */
}
CompDesc;

/* station label alias */
typedef struct {
    char name[ARRIVAL_LABEL_LEN]; /* original label (i.e. station name) */
    char alias[ARRIVAL_LABEL_LEN]; /* alias, i.e. new label */
    int byr, bmo, bday; /* begin year, month, day of validity */
    int eyr, emo, eday; /* end year, month, day of validity */
}
AliasDesc;

/* excluded phase desc */
typedef struct {
    char label[ARRIVAL_LABEL_LEN]; /* label (i.e. station name) */
    char phase[ARRIVAL_LABEL_LEN]; /* phase */
}
ExcludeDesc;

/* time delays */
typedef struct {
    char label[ARRIVAL_LABEL_LEN];
    char phase[ARRIVAL_LABEL_LEN];
    int n_residuals;
    double delay; /* time delay (sec) */
    double std_dev; /* std dev of time delay (sec) */
}
TimeDelayDesc;



/* Event time information extracted from filename or hypocenter line in obs file */

/*  SH 03/05/2004  added event number ev_nr which is used by SED_LOC format */
typedef struct {
    int ev_nr;
    int year, month, day;
    int hour, min;
    double sec;
}
EventTimeExtract;



/* Octtree */

#define OCTREE_UNDEF_VALUE -VERY_SMALL_DOUBLE

typedef struct {
    int init_num_cells_x, init_num_cells_y, init_num_cells_z;
    // num nodes for each side of initial Tree3D
    double min_node_size; // size of smallest side of smallest allowed node
    int max_num_nodes; // maximum number of nodes to evaluate
    int num_scatter; // number of scatter points to output
    int use_stations_density; // if 1, weight oct node order in result tree by station density in node
    int stop_on_min_node_size; // if 1, stop search when first min_node_size reached, if 0 stop subdividing a given cell when min_node_size reached
    double mean_cell_velocity; // mean cell velocity (>0 = Increases misfit sigma in proportion to travel time across diagonal of cell.
}
OcttreeParams;

/* Metropolis */
typedef struct {
    double x, y, z;
    long double dt;
    double dx;
    double likelihood;
    double velocity;
    double initial_temperature;
}
WalkParams;

/* Magnitude */
typedef struct {
    int type; /* magnitude calculation method */

    /* ML - Hutton & Boore, BSSA, v77, n6, Dec 1987 */
    double amp_fact_ml_hb; /* amplitude scale factor */
    double hb_n; /* "n" in H&B eq. (2) */
    double hb_K; /* "K" in H&B eq. (2) */
    double hb_Ro; /* reference distance "100" in H&B eq. (2) */
    double hb_Mo; /* reference magnitude "3.0" in H&B eq. (2) */

    /* coda duration (FMAG) - HYPOELLIPSE users manual chap 4;
            Lee et al., 1972; Lahr et al., 1975; Bakun and Lindh, 1977 */
    double fmag_c1, fmag_c2, fmag_c3, fmag_c4, fmag_c5;
}
MagDesc;

/* Search PDF grid (PRIOR or POSTERIOR */
// 20190510 AJL - added

typedef struct {
    int gridType; // grid type, e.g. GRID, OCT_TREE
    double default_value; // prior value to use where prior undefined
    char grid_file_path[FILENAME_MAX];
    // GRID
    GridDesc grid; // prior grid containing relative prior values over superset or subset of location search volume
    // OCT-TREE
    int nGrids;
    Tree3D **tree3D;    // array of nGrids tree3D oct-trees
    double *coherence;  // array of nGrids coherences for each corresponding tree3D
    double coherence_min;  // minimum coherence to use grid
    double max_total_other_weight;  // maximum total of coherence weight for other events; if exceeded, other event coherences normalized to sum to this value
    double *weight;  // array of nGrids weights (function of coherence and coherence_min) for each corresponding tree3D
    // arrivals storage
    ArrivalDesc** first_motion_arrivals; // arrivals with first-motion readings
    int *nfirst_motion_arrivals; // number of arrivals with first-motion readings
}
SearchPdfGridDesc;



/*------------------------------------------------------------*/
/* globals  */

/* output file path  */
EXTERN_TXT char f_outpath[FILENAME_MAX];

/* Gaussian error parameters */
EXTERN_TXT GaussLocParams Gauss;
EXTERN_TXT Gauss2LocParams Gauss2;
EXTERN_TXT int iUseGauss2;

/* Scatter parameters */
EXTERN_TXT ScatterParams Scatter;


/* events */
EXTERN_TXT int NumEvents;
EXTERN_TXT int NumEventsLocated;
EXTERN_TXT int NumLocationsCompleted;

// 20200107 AJL  #define MAX_NUM_OBS_FILES 10000
#define MAX_NUM_OBS_FILES 20000  // 20200107 AJL
EXTERN_TXT int NumObsFiles;

/* number of arrivals read from obs file */
EXTERN_TXT int NumArrivalsRead;

/* number of arrivals used for location */
EXTERN_TXT int NumArrivalsLocation;

/* observations filenames */
EXTERN_TXT char fn_loc_obs[MAX_NUM_OBS_FILES][FILENAME_MAX];
/* filetype */
EXTERN_TXT char ftype_obs[MAXLINE];

/* filenames */
EXTERN_TXT char fn_loc_grids[FILENAME_MAX], fn_path_output[FILENAME_MAX];
EXTERN_TXT int iSwapBytesOnInput;

// model files
EXTERN_TXT FILE *fp_model_grid_P;
EXTERN_TXT FILE *fp_model_hdr_P;
EXTERN_TXT GridDesc model_grid_P;
EXTERN_TXT FILE *fp_model_grid_S;
EXTERN_TXT FILE *fp_model_hdr_S;
EXTERN_TXT GridDesc model_grid_S;

/* location search type (grid, simulated annealing, Metropolis, etc) */
#define SEARCH_GRID   0
#define SEARCH_MET   1
#define SEARCH_OCTTREE  2
EXTERN_TXT int SearchType;

#define PDF_GRID_UNDEF   0
#define PDF_GRID_GRID   1
#define PDF_GRID_OCT_TREE   2
#define PDF_GRID_PRIOR   0
#define PDF_GRID_POSTERIOR   1
#define MAX_NUM_PDF_GRID_FILES 5000
// location search prior  // 20190510 AJL - added
EXTERN_TXT SearchPdfGridDesc SearchPrior;
EXTERN_TXT int iUseSearchPrior;
EXTERN_TXT SearchPdfGridDesc SearchPosterior;
EXTERN_TXT int iUseSearchPosterior;

/* location method (misfit, etc) */
#define METH_UNDEF    0
#define METH_GAU_ANALYTIC   1
#define METH_GAU_TEST    2
#define METH_EDT    3
#define METH_EDT_BOX    4
#define METH_ML_OT    5
#define METH_OT_STACK    6
#define METH_L1_NORM    7         // 20140515 AJL - added for NLDiffLoc
EXTERN_TXT int LocMethod;
EXTERN_TXT int EDT_use_otime_weight;
EXTERN_TXT int EDT_otime_weight_active;
EXTERN_TXT double DistStaGridMin;
EXTERN_TXT double DistStaGridMax;
EXTERN_TXT int MinNumArrLoc;
EXTERN_TXT int MaxNumArrLoc;
EXTERN_TXT int MinNumSArrLoc;
EXTERN_TXT double VpVsRatio;

/* location signature */
EXTERN_TXT char LocSignature[MAXLINE_LONG];

/* location grids */
#define MAX_NUM_LOCATION_GRIDS 10
EXTERN_TXT GridDesc LocGrid[MAX_NUM_LOCATION_GRIDS];
EXTERN_TXT int NumLocGrids;
EXTERN_TXT int LocGridSave[MAX_NUM_LOCATION_GRIDS]; /* !should be in GridDesc */
//EXTERN_TXT int Num3DGridReadToMemory, MaxNum3DGridMemory;

/* related hypocenter file pointers */
FILE *pSumFileHypNLLoc[MAX_NUM_LOCATION_GRIDS];
FILE *pSumFileHypo71[MAX_NUM_LOCATION_GRIDS];
FILE *pSumFileHypoEll[MAX_NUM_LOCATION_GRIDS];
FILE *pSumFileHypoInv[MAX_NUM_LOCATION_GRIDS];
FILE *pSumFileHypoInvY2K[MAX_NUM_LOCATION_GRIDS];
FILE *pSumFileAlberto4[MAX_NUM_LOCATION_GRIDS];
FILE *pSumFileFmamp[MAX_NUM_LOCATION_GRIDS];


/* related flags */
EXTERN_TXT int iWriteHypHeader[MAX_NUM_LOCATION_GRIDS];

/* format specific event data */
EXTERN_TXT char HypoInverseArchiveSumHdr[MAXLINE_LONG];

/* hypocenter filetype saving flags */
/* SH 02/26/2004  added iSaveSnapSum for output to be read
    by SNAP */
EXTERN_TXT int iSaveNLLocEvent, iSaveNLLocSum, iSaveNLLocOctree,
iSaveHypo71Event, iSaveHypo71Sum,
iSaveHypoEllEvent, iSaveHypoEllSum,
iSaveHypoInvSum, iSaveHypoInvY2KArc, iSaveAlberto4Sum, iSaveFmamp,
iSaveSnapSum, iCalcSedOrigin, iSaveDecSec, iSaveNone;
// 20170811 AJL - added to allow saving of expectation hypocenter results instead of maximum likelihood
EXTERN_TXT int iSaveNLLocExpectation;


// Arrival prior weighting flag (NLL_FORMAT_VER_2)
EXTERN_TXT int iUseArrivalPriorWeights;

/* station distance weighting flag. */
EXTERN_TXT int iSetStationDistributionWeights;
EXTERN_TXT double stationDistributionWeightCutoff;

/* station density weghting */
EXTERN_TXT double AveInterStationDistance;
EXTERN_TXT int NumForceOctTreeStaDenWt;

EXTERN_TXT int iRejectDuplicateArrivals;

/* Event Information extracted from phase file */
EXTERN_TXT EventTimeExtract EventTime;
EXTERN_TXT long int EventID;

/* magnitude calculation */
#define MAG_UNDEF   0
#define MAG_ML_HB   1
#define MAG_MD_FMAG   2
EXTERN_TXT int NumMagnitudeMethods;
#define MAX_NUM_MAG_METHODS   2
EXTERN_TXT MagDesc Magnitude[MAX_NUM_MAG_METHODS];

/* station/inst/component parameters */
#define MAX_NUM_COMP_DESC 1000
EXTERN_TXT CompDesc Component[MAX_NUM_COMP_DESC];
EXTERN_TXT int NumCompDesc;

/* arrival label alias */
#define MAX_NUM_LOC_ALIAS 1000
#define MAX_NUM_LOC_ALIAS_CHECKS 2*MAX_NUM_LOC_ALIAS
EXTERN_TXT AliasDesc LocAlias[MAX_NUM_LOC_ALIAS];
EXTERN_TXT int NumLocAlias;

/* exclude arrivals */
#define MAX_NUM_LOC_EXCLUDE 1000
EXTERN_TXT ExcludeDesc LocExclude[MAX_NUM_LOC_EXCLUDE];
EXTERN_TXT int NumLocExclude;

/* include arrivals */
#define MAX_NUM_LOC_INCLUDE 1000
EXTERN_TXT ExcludeDesc LocInclude[MAX_NUM_LOC_INCLUDE];
EXTERN_TXT int NumLocInclude;

/* station delays */
#define WRITE_RESIDUALS 0
#define WRITE_RES_DELAYS 1
#define WRITE_PDF_RESIDUALS 2
#define WRITE_PDF_DELAYS 3
#define MAX_NUM_STA_DELAYS 10000
EXTERN_TXT TimeDelayDesc TimeDelay[MAX_NUM_STA_DELAYS];
EXTERN_TXT int NumTimeDelays;

EXTERN_TXT char TimeDelaySurfacePhase[MAX_SURFACES][PHASE_LABEL_LEN];
EXTERN_TXT double TimeDelaySurfaceMultiplier[MAX_SURFACES];
EXTERN_TXT int NumTimeDelaySurface;

/* crustal and elev corrections */
EXTERN_TXT int ApplyElevCorrFlag;
EXTERN_TXT double ElevCorrVelP;
EXTERN_TXT double ElevCorrVelS;
EXTERN_TXT int ApplyCrustElevCorrFlag;
EXTERN_TXT double MinDistCrustElevCorr;

/* topo surface */
EXTERN_TXT struct surface *topo_surface;
EXTERN_TXT int topo_surface_index; // topo surface index is velmod.h.MAX_SURFACES-1 so as not to interferce with any TimeDelaySurfaces read in



/* station list */
int NumStationPhases;
SourceDesc StationPhaseList[X_MAX_NUM_ARRIVALS];

/* fixed origin time parameters */
EXTERN_TXT int FixOriginTimeFlag;

/* Metropolis */
EXTERN_TXT WalkParams Metrop; /* walk parameters */
EXTERN_TXT int MetNumSamples; /* number of samples to evaluate */
EXTERN_TXT int MetLearn; /* learning length in number of samples for
					calculation of sample statistics */
EXTERN_TXT int MetEquil; /* number of samples to equil before using */
EXTERN_TXT int MetStartSave; /* number of sample to begin saving */
EXTERN_TXT int MetSkip; /* number of samples to wait between saves */
EXTERN_TXT double MetStepInit; /* initial step size (km) (< 0.0 for auto) */
EXTERN_TXT double MetStepMin; /* minimum step size (km) */
EXTERN_TXT double MetStepMax; /* maximum step size (km) (NLDiffLoc) */
EXTERN_TXT double MetStepFact; /* step size factor */
EXTERN_TXT double MetProbMin; /* minimum likelihood necessary after learn */
EXTERN_TXT double MetVelocity; /* velocity for conversion of distance to time */
EXTERN_TXT double MetInititalTemperature; /* initial temperature */
EXTERN_TXT int MetUse; /* number of samples to use
					= MetNumSamples - MetEquil */


/* Octtree */
EXTERN_TXT OcttreeParams octtreeParams; /* Octtree parameters */
EXTERN_TXT Tree3D* octTree; /* Octtree */
EXTERN_TXT ResultTreeNode* resultTreeRoot; /* Octtree likelihood*volume results tree root node */
//EXTERN_TXT ResultTreeNode* resultTreeLikelihoodRoot;	/* Octtree likelihood results tree root node */


/* take-off angles */
EXTERN_TXT int angleMode; /* angle mode - ANGLE_MODE_NO, ANGLE_MODE_YES */
EXTERN_TXT int iAngleQualityMin; /* minimum quality for angles to be used */
#define ANGLE_MODE_NO 0
#define ANGLE_MODE_YES 1
#define ANGLE_MODE_UNDEF -1


/* otime list */
EXTERN_TXT OtimeLimit** OtimeLimitList;
EXTERN_TXT int NumOtimeLimit;



/*------------------------------------------------------------*/
/** hashtable routines for accumulating station statistics */

/* from Kernigham and Ritchie, C prog lang, 2nd ed, 1988, sec 6.6 */

struct staStatNode { /* station statistics node for linked list */

    struct staStatNode *next; /* next entry in list */
    char label[ARRIVAL_LABEL_LEN]; /* arrival label (station name) */
    char phase[ARRIVAL_LABEL_LEN]; /* arrival phase id */
    int flag_ignore; /* ignore flag  = 1 if phase not used for misfit calc */
    /* standard (max like point) residuals */
    double residual_min; /* minimum residual */
    double residual_max; /* maximum residual */
    double residual_sum; /* accumulated residuals */
    double residual_square_sum; /* accumulated square of residuals */
    double weight_sum; /* accumulated weights */
    int num_residuals; /* number of residuals accumulated */
    /* statistical (PDF weighted) residuals */
    double pdf_residual_sum; /* accumulated pdf residuals */
    double pdf_residual_square_sum; /* accumulated square of pdf residuals */
    int num_pdf_residuals; /* number of residuals accumulated */
    double delay; /* input time delay */

};
typedef struct staStatNode StaStatNode;

#define HASHSIZE 46
/* table of pointers to list of StaStatNode */
EXTERN_TXT StaStatNode *hashtab[MAX_NUM_LOCATION_GRIDS][HASHSIZE];
/* maxumum residual values to include in statistics */
EXTERN_TXT int NRdgs_Min;
EXTERN_TXT double RMS_Max, Gap_Max;
EXTERN_TXT double P_ResidualMax;
EXTERN_TXT double S_ResidualMax;
EXTERN_TXT double Ell_Len3_Max;
EXTERN_TXT double Hypo_Depth_Min;
EXTERN_TXT double Hypo_Depth_Max;
EXTERN_TXT double Hypo_Dist_Max;

/* hashtable function declarations */
StaStatNode *InstallStaStatInTable(int, char*, char*, int, double,
        double, double, double, double);
int FreeStaStatTable(int ntable);
int WriteStaStatTable(int, FILE *, double, int, double,
        double, double, double, double, double, double, int);
void UpdateStaStat(int, ArrivalDesc *, int, double, double, double, double);

/** end of hashtable routines */
/*------------------------------------------------------------*/




/*------------------------------------------------------------*/
/* function declarations */

int NLLoc(char *pid_main, char *fn_control_main, char **param_line_array, int n_param_lines, char **obs_line_array, int n_obs_lines,
        int return_locations, int return_oct_tree_grid, int return_scatter_sample, LocNode **ploc_list_head);

int Locate(int ngrid, char* fn_loc_obs, char* fn_root_out, int numArrivalsReject, int return_locations, int return_oct_tree_grid, int return_scatter_sample, LocNode **ploc_list_head);

int checkObs(ArrivalDesc *arrival, int nobs);
int ExtractFilenameInfo(char*, char*);
int ReadNLLoc_Input(FILE* fp_input, char **param_line_array, int n_param_lines);
int GetNLLoc_Grid(char*);
int GetNLLoc_HypOutTypes(char*);
int GetStaWeight(char* line1);
int GetNLLoc_Gaussian(char*);
int GetNLLoc_Gaussian2(char*);
int GetNLLoc_PhaseStats(char*);
int GetNLLoc_Angles(char*);
int GetNLLoc_Magnitude(char*);
int GetNLLoc_Files(char*);
int GetNLLoc_Method(char*);
int GetNLLoc_SearchType(char*);
int GetNLLoc_PdfGrid(char*, int);
int GetNLLoc_FixOriginTime(char*);
int GetObservations(FILE*, char*, char*, ArrivalDesc*, int*, int*, int*, int, HypoDesc*, int*, int*, int);
int GetNextObs(HypoDesc* phypo, FILE*, ArrivalDesc *, char*, int);
int IsGoodDate(int, int, int);
int ReadArrivalSheets(int, ArrivalDesc*, double);
int IsSameArrival(ArrivalDesc *, int, int, char *);
int IsDuplicateArrival(ArrivalDesc *, int, int, int);
int FindDuplicateTimeGrid(ArrivalDesc *arrival, int num_arrivals, int ntest);

int WriteHypo71(FILE *, HypoDesc*, ArrivalDesc*, int, char*, int, int);
int WriteHypoEll(FILE *, HypoDesc*, ArrivalDesc*, int, char*, int, int);
int WriteHypoInverseArchive(FILE *fpio, HypoDesc *phypo, ArrivalDesc *parrivals, int narrivals,
        char *filename, int writeY2000, int write_arrivals, double arrivalWeightMax);
int WriteHypoAlberto4(FILE *, HypoDesc*, ArrivalDesc*, int, char*);
int WriteHypoFmamp(FILE *fpio, HypoDesc* phypo, ArrivalDesc* parrivals, int narrivals, char* filename, int write_header);
int WriteHypoFmampSearchPosterior(SearchPdfGridDesc *searchPdfGrid, FILE *fpio, HypoDesc* phypo, char* filename, int write_header);
int OpenSummaryFiles(char *, char *);
int CloseSummaryFiles();

int GetCompDesc(char*);
int GetLocAlias(char*);
int GetLocExclude(char* line1);
int GetLocInclude(char* line1);
int GetTimeDelays(char*);
int GetTimeDelaySurface(char*);

int GetTopoSurface(char*);

int LocGridSearch(int, int, int, ArrivalDesc*, GridDesc*,
        GaussLocParams*, HypoDesc*);
int LocMetropolis(int, int, int, ArrivalDesc *,
        GridDesc*, GaussLocParams*, HypoDesc*, WalkParams*, float*);
int SaveBestLocation(OctNode* poct_node, int num_arr_total, int num_arr_loc, ArrivalDesc *arrival,
        GridDesc* ptgrid, GaussLocParams* gauss_par, HypoDesc* phypo,
        double misfit_max, int iGridType, int ignore_pred_travel_time_best,
        double cell_diagonal_time_var_best, double cell_diagonal_best, double cell_volume_best);
int ConstWeightMatrix(int, ArrivalDesc*, GaussLocParams*);
int CleanWeightMatrix();
void CalcCenteredTimesObs(int, ArrivalDesc*, GaussLocParams*, HypoDesc*);
void CalcCenteredTimesPred(int, ArrivalDesc*, GaussLocParams*);
double CalcSolutionQuality(double hypo_x, double hypo_y, double hypo_z, OctNode* poct_node, int num_arrivals, ArrivalDesc *arrival, GaussLocParams* gauss_par, int itype,
        double* pmisfit, double* potime, double* potime_var, double cell_diagonal_time_var, double cell_diagonal, double cell_volume, double* effective_cell_size, double *pot_variance_factor, double *prior);
double CalcSolutionQuality_GAU_ANALYTIC(int, ArrivalDesc*, GaussLocParams*, int, double*, double*);
double CalcSolutionQuality_GAU_TEST(int, ArrivalDesc*, GaussLocParams*, int, double*, double*);
double CalcSolutionQuality_L1_NORM(int num_arrivals, ArrivalDesc *arrival,
        GaussLocParams* gauss_par, int itype, double* pmisfit, double* potime);
double CalcSolutionQuality_EDT(int num_arrivals, ArrivalDesc *arrival, GaussLocParams* gauss_par, int itype, double* pmisfit, double* potime, double* potime_var, double cell_diagonal_time_var, int method_box);
double CalcSolutionQuality_OT_STACK(OctNode* poct_node, int num_arrivals, ArrivalDesc *arrival,
        GaussLocParams* gauss_par, int itype, double* pmisfit, double* potime, double* potime_var,
        double cell_half_diagonal_time_range, double cell_diagonal, double cell_volume, double* effective_cell_size, double *pot_variance_factor);
double CalcSolutionQuality_ML_OT(int num_arrivals, ArrivalDesc *arrival, GaussLocParams* gauss_par, int itype, double* pmisfit, double* potime, double* potime_var, double cell_diagonal_time_var, int method_box);
double calc_maximum_likelihood_ot_sort(
        OctNode* poct_node, int num_arrivals, ArrivalDesc *arrival,
        double cell_half_diagonal_time_range, double cell_diagonal, double cell_volume, double *pot_var, int icalc_otime,
        double *plog_prob_max, double *pot_stack_weight, double* effective_cell_size, double *pot_variance_factor);
double calc_maximum_likelihood_ot(double *ot_ml_arrival, double *ot_ml_arrival_edt_sum, int num_arrivals, ArrivalDesc *arrival, MatrixDouble edtmtx, double *pot_ml_var, int iwrite_errors,
        double *pprob_max);
double calc_likelihood_ot(double *ot_ml_arrival, double *ot_ml_arrival_edt_sum, int num_arrivals, ArrivalDesc *arrival, MatrixDouble edtmtx, double time);
double calc_variance_ot(double *ot_ml_arrival, double *ot_ml_arrival_edt_sum, int num_arrivals, ArrivalDesc *arrival, MatrixDouble edtmtx, double expectation_time);
long double CalcMaxLikeOriginTime(int, ArrivalDesc*, GaussLocParams*);
int NormalizeWeights(int num_arrivals, ArrivalDesc *arrival);
void UpdateProbabilisticResiduals(int, ArrivalDesc *, double);
int CalcConfidenceIntrvl(GridDesc*, HypoDesc*, char*);
int HomogDateTime(ArrivalDesc*, int, HypoDesc*);
int CheckAbsoluteTiming(ArrivalDesc *arrival, int num_arrivals);
int hypotime2hrminsec(long double phypo_time, int *phypo_hour, int *phypo_min, double *phypo_sec);
int StdDateTime(ArrivalDesc*, int, HypoDesc*);
int SetOutName(ArrivalDesc *arrival, char* out_file_root, char* out_file, char lastfile[FILENAME_MAX], int isec, int *pncount);
int SaveLocation(HypoDesc* hypo, int ngrid, char* fnobs, char *fnout, int numArrivalsReject,
        char* typename, int isave_phases, GaussLocParams* gauss_par);
void InitializeArrivalFields(ArrivalDesc *);
int isExcluded(char *label, char *phase);
int EvaluateArrivalAlias(ArrivalDesc *);
int ApplyTimeDelays(ArrivalDesc *);
double ApplySurfaceTimeDelay(int nsurface, ArrivalDesc *arrival);
double CalcSimpleElevCorr(ArrivalDesc *arrival, int narr, double pvel, double svel);

double CalcAzimuthGap(ArrivalDesc *arrival, int num_arrivals, double *pgap_secondary);
double CalcArrivalDistances(ArrivalDesc *arrival, int num_arrivals, double *pmaximumDistance, double *pmedianDistance, int usedStationCount);
int CalcArrivalCounts(ArrivalDesc *arrival, int num_arrivals, int num_arrivals_read, int* passociatedPhaseCount, int* passociatedStationCount, int* pusedStationCount, int* pdepthPhaseCount);
void InitializeMetropolisWalk(GridDesc*, ArrivalDesc*, int,
        WalkParams*, int, double);
int GetNextMetropolisSample(WalkParams*, double, double,
        double, double, double, double, double*, double*, double*);
int MetropolisTest(double, double);

double CalculateVpVsEstimate(HypoDesc* phypo, ArrivalDesc* parrivals, int narrivals);

int CalculateMagnitude(HypoDesc*, ArrivalDesc*, int, CompDesc*, int, MagDesc*);
int findStaInstComp(ArrivalDesc* parr, CompDesc* pcomp, int nCompDesc);
double Calc_ML_HuttonBoore(double amplitude, double dist, double depth, double sta_corr,
        double hb_n, double hb_K, double hb_Ro, double hb_Mo);
double Calc_MD_FMAG(double coda_dur, double dist, double depth, double sta_corr,
        double fmag_c1, double fmag_c2, double fmag_c3, double fmag_c4, double fmag_c5);

int setStationDistributionWeights(SourceDesc *stations, int numStations, ArrivalDesc *arrival, int nArrivals);

int getTravelTimes(ArrivalDesc *arrival, int num_arr_loc, double xval, double yval, double zval);
double applyCrustElevCorrection(ArrivalDesc* parrival, double xval, double yval, double zval);
int isAboveTopo(double xval, double yval, double zval);

Tree3D* InitializeOcttree(GridDesc* ptgrid, OcttreeParams* pParams);
int LocOctree(int ngrid, int num_arr_total, int num_arr_loc,
        ArrivalDesc *arrival,
        GridDesc* ptgrid, GaussLocParams* gauss_par, HypoDesc* phypo,
        OcttreeParams* pParams, Tree3D* pOctTree, float* fdata,
        double *poct_node_value_max, double *poct_tree_integral);
long double LocOctree_core(int ngrid, double xval, double yval, double zval,
        int num_arr_loc, ArrivalDesc *arrival,
        OctNode* poct_node,
        int icalc_cell_diagonal_time_var, double *volume_min,
        double *diagonal, double *cell_diagonal_time_var,
        OcttreeParams* pParams, GaussLocParams* gauss_par, int iGridType,
        double *misfit, double logWtMtrxSum);
double getOctTreeStationDensityWeight(OctNode* poct_node, SourceDesc *stations, int numStations, GridDesc *pgrid, int iOctLevelMax);
int GenEventScatterOcttree(OcttreeParams* pParams, double oct_node_value_max, float* fscatterdata, double integral, HypoDesc* Hypocenter);

int GetElevCorr(char* line1);

/*------------------------------------------------------------*/




/*------------------------------------------------------------*/
/** custom build header items */

//#ifdef CUSTOM_ETH
/* PID of SNAP
   this needs to be given as second argument
   and is needed to construct output file name for SNAP
   SH  02/27/2004 */
//EXTERN_TXT char snap_pid[255];
/** function declarations */
/* added this function for output in SED SNAP format SH 02/27/2004 */
//int WriteSnapSum (FILE *fpio, HypoDesc* phypo, ArrivalDesc* parrivals);
//#endif

/*------------------------------------------------------------*/




/* */
/*------------------------------------------------------------/ */


