/*------------------------------------------------------------*/
/** ETH build header items */

/* PID of SNAP
   this needs to be given as second argument
   and is needed to construct output file name for SNAP
   SH  02/27/2004 */
EXTERN_TXT char snap_pid[255];
/* EXTERN_TXT char snap_param_file[1024];  SH 02AUG2004 not needed any more */
//EXTERN_TXT char snap_file_PATH1[1024];
#define snap_file_PATH1 "/usr/local/lib/fereg96.poly"
//EXTERN_TXT char snap_file_PATH2[1024];
#define snap_file_PATH2 "/usr/local/lib/ch_reg.poly"
//EXTERN_TXT char snap_file_PATH3[1024];
#define snap_file_PATH3 "/usr/local/lib/regions.bin"
//EXTERN_TXT char snap_file_SED_CAL[1024];
#define snap_file_SED_CAL "/usr/local/lib/SED.CAL2"


EXTERN_TXT char last_snap_filename[FILENAME_MAX];
EXTERN_TXT int n_snap_files;

/** function declarations */
/* added this function for output in SED SNAP format SH 02/27/2004 */
int WriteSnapSum (FILE *fpio, HypoDesc *phypo, ArrivalDesc *parrivals, int narrivals);
//int CalculateMagnitudeSNAP(char o_event_type, HypoDesc* phypo, ArrivalDesc* parrivals, int narrivals, char* snap_file_SED_CAL);
int CalculateMagnitudeSNAP(char o_event_type, HypoDesc* phypo, ArrivalDesc* parrivals, int narrivals);
int ArrivalsToPhases(ArrivalDesc* parrivals, int narrivals);
void magnitude(int jmin, float* mag, float* magerr, char* cmaxmag, int* num_mag, char o_event_type, double g13_depth);


/*------------------------------------------------------------*/

