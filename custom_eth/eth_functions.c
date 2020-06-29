/*
    29JUL2004 SH   did some bug fixing; added use of new variables ArrivalDesc->clipped and 
                   HypoDesc->mag_err
				   
        30JUL2004 SH   did changes to cope with teleseismic events: print out in WriteSnapSum is
                                   different, o_event_type is now defined using map_trans_type
				   
        02AUG2004 SH   discard the use of the SNAP parameter file introduced by AJL; calibration and region files are
                       now defined in eth_functions.h
				   
        25MAY2006 SH   changes to cope cope with modifications in magnitude.c by Manfred Baer
	
        15JUN2006 SH   bug fix in routine WriteSnapSum
	
        18JUL2006 SH   bug fix in routine WriteSnapSum:
                         undefined first motion is printed out as " ", but left as "?" for NLLoc output
 */


#define EXTERN_MODE 1
#include "../GridLib.h"
#include "../ran1.h"
#include "../velmod.h"
#include "../GridMemLib.h"
#include "../calc_crust_corr.h"
#include "../phaseloclist.h"
#include "../otime_limit.h"
#include "../NLLocLib.h"

#define MAIN
#include "eth_functions.h"
#include "get_region_name_nr.h"
#include "grid_search.h"
#include "new_sedlib.h"
//#include "gsetype.h"
//#include "gselib.h"

/* SH 02/26/2004
   new function to write hypocenter summary to file (SNAP format)  */

int WriteSnapSum(FILE *fpio, HypoDesc* phypo, ArrivalDesc* parrivals, int narrivals) {

    //	 FILE *fp_param;
    ArrivalDesc* parr;
    int ifile = 0;
    //int jfile = 0;
    int i, narr;
    int reg_nr;
    int eth_use_loc;
    float swissx, swissy;

    char region_name[81];
    char magrmk, qual, firstmot[2];
    //double errx,erry,errh,errz;
    double errx, erry, errz;
    double diff;
    char new_snap_filename[FILENAME_MAX];
    char snap_filename[FILENAME_MAX];

    /* SH 30JUL2004
        o_event_type can be L for local events and T for teleseismic
             char o_event_type = 'L';  */
    char o_event_type;
    char ns, ew;


    // set globals
    twopi = (float) 2.0 * PI;
    WA_gain = (float) 2800.0; // ???

    /* SH 02AUG2004 use of snap_parameter file has been discarded; files names and locations
                    are now defined in eth_function.h */

    /* open snap_param file */
    /*	   if ((fp_param = fopen(snap_param_file, "r")) == NULL) {
                 puterr2("ERROR: opening SNAP parameter intput file", snap_param_file);
                     return(-1);
               } else {
                 jfile = 1;
                 NumFilesOpen++;
               } */
    /* read parameters */
    /*	if (fscanf(fp_param, "%s", snap_file_PATH1) < 1)
                 puterr2("ERROR: reading snap_file_PATH1 from SNAP parameter intput file", snap_param_file);
            if (fscanf(fp_param, "%s", snap_file_PATH2) < 1)
                 puterr2("ERROR: reading snap_file_PATH2 from SNAP parameter intput file", snap_param_file);
            if (fscanf(fp_param, "%s", snap_file_PATH3) < 1)
                 puterr2("ERROR: reading snap_file_PATH3 from SNAP parameter intput file", snap_param_file);
            if (fscanf(fp_param, "%s", snap_file_SED_CAL) < 1)
                 puterr2("ERROR: reading snap_file_SED_CAL from SNAP parameter intput file", snap_param_file); */
    /* close file */
    /*	if (jfile) {
                    fclose(fp_param);
                    NumFilesOpen--;
            } */



    /* open output file */

    if (fpio == NULL) {
        sprintf(new_snap_filename, "%shyprint%s", fn_path_output, snap_pid);
        printf("DEBUG: new_snap_filename %s\n", new_snap_filename);
        // check if this filename already used in this run
        if (strcmp(new_snap_filename, last_snap_filename) == 0) {
            sprintf(snap_filename, "%s.%d", new_snap_filename, ++n_snap_files);
        } else {
            n_snap_files = 0;
            strcpy(snap_filename, new_snap_filename);
        }
        if ((fpio = fopen(snap_filename, "w")) == NULL) {
            nll_puterr("ERROR: opening SNAP summary output file.");
            /*			if (jfile) {
                                    fclose(fp_param);
                                    NumFilesOpen--;
                                    } */
            return (-1);
        } else {
            ifile = 1;
            NumFilesOpen++;
        }
    }
    strcpy(last_snap_filename, new_snap_filename);

    /* SH 30JUL2004
       get type of event  */
    if (strcmp(map_trans_type[0], "GLOBAL") == 0) {
        o_event_type = 'T'; /* teleseismic events */
    } else {
        o_event_type = 'L'; /* local events */
    }


    /* compute region name and number, and Swiss coordinates of event
         take M. Baer routine get_region_name_nr.c  */

    i = get_region_name_nr(phypo->dlat, phypo->dlong, -1, &reg_nr, region_name, &swissx, &swissy);

    /* AJL 20040513 compute magnitude using M. Baer routines */
    //CalculateMagnitudeSNAP(o_event_type, phypo, parrivals, narrivals, snap_file_SED_CAL);
    CalculateMagnitudeSNAP(o_event_type, phypo, parrivals, narrivals);

    /* compute error in epicenter and depth based on covariance matrix, i.e.
       projection of 3D error ellipsoid as defined by covariance matrix
       onto single confidence intervals in x,y,z
       (see NR Chap. 15.6, eq. 15.6.4*/

    if (o_event_type == 'L') { /* local events in km */
        errx = sqrt(phypo->cov.xx);
        erry = sqrt(phypo->cov.yy);
        errz = sqrt(phypo->cov.zz);
    } else { /* teleseismic in deg */
        errx = sqrt(phypo->cov.xx) / 111.17;
        erry = sqrt(phypo->cov.yy) / 111.17;
        errz = sqrt(phypo->cov.zz) / 111.17;
    }

    /* determine quality factor
           A: RMS < 0.5 s; diff < 0.5 km; errh < 2.0 km & errz < 2.0 km
           B: RMS < 0.5 s; diff < 0.5 km; errh >= 2.0 km & errz >= 2.0 km
           C: RMS < 0.5 s; diff >= 0.5 km;
           C: RMS >= 0.5 s   */

    /* difference between maximum likelihood and expectation hypocenter locations */
    diff = sqrt((phypo->expect.x - phypo->x)*(phypo->expect.x - phypo->x) +
            (phypo->expect.y - phypo->y)*(phypo->expect.y - phypo->y) +
            (phypo->expect.z - phypo->z)*(phypo->expect.z - phypo->z));
    /* SH 29JUL2004 added */
    printf("\nWriteSnapSum: diff maximum likelihood and expectation hypo: %f\n", diff);

    qual = '-';
    if (phypo->rms >= 0.5) {
        qual = 'D';
    } else if (diff > 0.5) {
        qual = 'C';
    } else if ((errx > 2.0 || erry > 2.0) && errz > 2.0) {
        qual = 'B';
    } else {
        qual = 'A';
    }

    /* now start with output */

    fprintf(fpio, "  E V E N T   N R . %8d                  MODEL: 1     STATION LIST: 1\n\n",
            EventTime.ev_nr);

    /* SH 07/26/2004
       for SNAP this should read Grid Search Solution, not Global-Search Soluion */
    fprintf(fpio, " Grid Search Solution: NonLinLoc (Ver. %s %s, (c) Anthony Lomax - anthony@alomax.net) \n\n",
            PVER, PDATE);

    if (strcmp(phypo->locStat, "LOCATED") == 0) { /* successful location */

        /* SH 30JUL2004
            print out depending on event type */

        /* info on hypocenter location  */
        switch (o_event_type) {
            case 'L': /* local events */
                fprintf(fpio, "   DATE  ORIGIN TIME    LAT      LON    Z   AZ  D-MIN  Ml   RMS GAP  NO NI Q\n");
                fprintf(fpio, "%5d%02d%02d %02d:%02d:%04.1f%7.3fN%8.3fE%5.1f%4d%7.1f %3.1f %5.2f%4d%4d 99 %c\n",
                        phypo->year,
                        phypo->month,
                        phypo->day,
                        phypo->hour,
                        phypo->min,
                        phypo->sec,
                        phypo->dlat,
                        phypo->dlong,
                        phypo->depth,
                        (int) parrivals->azim,
                        phypo->dist,
                        phypo->amp_mag > 0.0 ? phypo->amp_mag : 0.0,
                        phypo->rms,
                        (int) (0.5 + phypo->gap),
                        phypo->nreadings,
                        qual
                        );
                fprintf(fpio, "   ACCURACY: (+/-)  %7.3f%9.3f%6.2f           %4.1f\n\n",
                        erry, errx, errz, phypo->mag_err);
                break;
            case 'T': /* teleseismic events */
                if (phypo->dlat < 0) {
                    ns = 'S'; /* south */
                } else {
                    ns = 'N'; /* north */
                }
                if (phypo->dlong < 0) {
                    ew = 'W'; /* west */
                } else {
                    ew = 'E'; /* east */
                }
                fprintf(fpio, "   DATE  ORIGIN TIME  LAT    LON    Z   AZ  DELTA  Mb   RMS  NO NI Q\n");
                fprintf(fpio, "%5d%02d%02d %02d:%02d:%04.1f%5.1f%c%6.1f%c%4d%5d%7.1f%4.1f %5.2f%4d 99 %c\n",
                        phypo->year,
                        phypo->month,
                        phypo->day,
                        phypo->hour,
                        phypo->min,
                        phypo->sec,
                        phypo->dlat >= 0.0 ? phypo->dlat : -phypo->dlat,
                        ns,
                        phypo->dlong >= 0.0 ? phypo->dlong : -phypo->dlong,
                        ew,
                        (int) phypo->depth,
                        (int) parrivals->azim,
                        phypo->dist,
                        phypo->amp_mag > 0.0 ? phypo->amp_mag : 0.0,
                        phypo->rms,
                        phypo->nreadings,
                        qual);
                fprintf(fpio, "   ACCURACY: (+/-)  %6.3f%7.3f%5.2f           %4.1f\n\n",
                        erry, errx, errz, phypo->mag_err);
                break;
        } /* switch */

        /* SH 30JUL2004
           added */
        if (reg_nr >= 1000) {
            /* Swiss region and name, and Swiss coordinates */
            fprintf(fpio, " L+T NR:%4d  %-41sKOORDINATEN:%4d/%3d KM\n\n\n",
                    reg_nr, region_name, (int) (swissx + 0.5), (int) (swissy + 0.5));
        } else {
            /* global events */
            fprintf(fpio, " F-E NR:%4d  %s\n\n",
                    reg_nr, region_name);
        }

        /* phase list w/ residuals etc*/
        /* SH 30JUL2004 added */
        if (o_event_type == 'L') {
            fprintf(fpio, "STA    PHASE  RMK W MN P-SEC P-CAL P-RES   DT   WT     AMP    PER MAG    DELTA  AZIM ANGLE\n");
        } else { /* teleseismic events */
            fprintf(fpio, "STA    PHASE  RMK W MN P-SEC P-CAL P-RES   DT   WT     AMP    PER MAG    DELTA\n");
        }
        // AJL 20040527
        //	   for (narr = 0; narr < phypo->nreadings; narr++) {
        for (narr = 0; narr < narrivals; narr++) {
            parr = parrivals + narr;
            /* determmine wether arrival time pick and amplitude pick have been used or not */
            if (parr->error == 9999.9) {
                eth_use_loc = 0;
            } else {
                eth_use_loc = 1;
            }
            /* SH 07212004  magrmk is now based on value in parr->clipped */
            if (parr->clipped < 0) {
                magrmk = '<';
            } else if (parr->clipped > 0) {
                magrmk = '?';
            } else if (parr->phase[0] == 'S') {
                magrmk = '-';
            } else {
                magrmk = ' ';
            }
            /* SH 15JUN2006
                change '?' to ' ' for print out  */
            if (strcmp(parr->first_mot, "?") == 0) {
                strcpy(firstmot, " ");
            } else {
                strcpy(firstmot, parr->first_mot);
            }
            /* arrival time and resiudal */
            /* SH 30JUL2004
                change to print calc. arrival time relative to origin time, not calc. travel time */
            fprintf(fpio, "%-6s %-8s%1s%1s %1d %2d %5.2f%6.2f%6.2f 9.99%5.2f",
                    parr->label, parr->phase, parr->onset, firstmot, eth_use_loc, parr->min,
                    parr->sec, parr->sec - parr->residual, parr->residual, parr->weight);
            /* amplitude and magnitude */
            /* SH 15JUN2006
                 deal w/ amplitude values that are too small to compute magnitude */
            if ((int) parr->amplitude <= 0) {
                parr->amplitude = 0.0;
                parr->period = 0.0;
                fprintf(fpio, "%8d %c%5.2f     ",
                        (int) parr->amplitude, parr->inst[0], parr->period);
            } else if (parr->amp_mag <= -9.8) {
                fprintf(fpio, "%8d %c%5.2f    %c",
                        (int) parr->amplitude, parr->inst[0], parr->period, magrmk);
            } else {
                fprintf(fpio, "%8d %c%5.2f%4.1f%c",
                        (int) parr->amplitude, parr->inst[0], parr->period, parr->amp_mag, magrmk);
            }
            /* distance, azimuth, and take off angle */
            /* SH 30JUL2004 added */
            if (o_event_type == 'L') {
                fprintf(fpio, "%8.1f %5.1f %5.1f\n",
                        parr->dist, parr->azim, parr->ray_dip);
            } else { /* teleseismic events */
                fprintf(fpio, "%8.3f\n",
                        parr->dist);
            }
        } /* for (narr... */
    } else { /* location was not successful.... */
        fprintf(fpio, " ERROR\n");
        fprintf(fpio, "    %s\n", phypo->locStatComm);
    }
    /* close files */
    if (ifile) {
        fclose(fpio);
        NumFilesOpen--;
    }

    return (0);

} /* end routine WriteSnapSum */

/*** function to convert array of NLL ArrivalDesc to SNAP phases */

int ArrivalsToPhases(ArrivalDesc* parrivals, int narrivals) {
    /*
    #define PS_STR_LEN 9
    typedef struct {	// phase data structure
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
    } PRec;
     */

    ArrivalDesc* parr;
    int narr;
    PRec* pphs;

    VZmodel dummy_VZmodel;

    nphas = 0;
    for (narr = 0; narr < narrivals && narr < MAXPHASE; narr++) {

        parr = parrivals + narr;
        pphs = phases + narr;

        strncpy(pphs->name, parr->label, PS_STR_LEN < ARRIVAL_LABEL_LEN ? PS_STR_LEN : ARRIVAL_LABEL_LEN);
        strncpy(pphs->phase, parr->phase, PS_STR_LEN < PHASE_LABEL_LEN ? PS_STR_LEN : PHASE_LABEL_LEN);
        strncpy(pphs->phase, parr->phase, PS_STR_LEN < PHASE_LABEL_LEN ? PS_STR_LEN : PHASE_LABEL_LEN);
        strcpy(pphs->phr, "\0");
        strcpy(pphs->phi, "\0");
        strncpy(pphs->onset, parr->onset, 3 < 2 ? 3 : 2);
        pphs->use = 0;
        pphs->p_time_sc = (float) parr->sec;
        pphs->corr = (float) parr->delay; // sign +/- ???
        pphs->period = (float) parr->period;
        pphs->amplitude = (float) parr->amplitude;
        pphs->weight = (float) parr->weight;
        strncpy(pphs->rec_sys, parr->inst, 10 < INST_LABEL_LEN ? 10 : INST_LABEL_LEN);
        pphs->clipped = parr->clipped; /* SH 07212004 added "clipped" */
        pphs->mag = (float) parr->amp_mag;
        pphs->magrmk = '?'; /* SH 07212004 dummy value */
        pphs->delta = (float) parr->dist; // deg
        if (GeometryMode != MODE_GLOBAL)
            pphs->delta /= 111.17; // deg
        pphs->azim = (float) parr->azim;
        pphs->angle = (float) parr->ray_azim; // ???
        pphs->tab = NULL; // ???
        pphs->stntab_index = 0; // ???
        pphs->lmod = dummy_VZmodel; // ???
        pphs->elevcor = (float) 0.0;
        pphs->res = (float) parr->residual; // ???
        pphs->nrph = 0; // ???
        strcpy(pphs->amp_chan, parr->comp); /* SH 05252006 added component for amplitude */

        nphas++;

    }

    return (0);

}


/*** function to calculate magintude */
/* SH 07232004 
   added double magerr -> error in magnitude estimate (?= standard deviation?) */

//int CalculateMagnitudeSNAP(char o_event_type, HypoDesc* phypo, ArrivalDesc* parrivals, int narrivals, char* snap_file_SED_CAL)

int CalculateMagnitudeSNAP(char o_event_type, HypoDesc* phypo, ArrivalDesc* parrivals, int narrivals) {

    // NLL
    int narr;
    ArrivalDesc* parr;

    // ETH
    PRec* pphs;
    // calib
    int i;

    // magnitude
    int o_time_mn, num_mag;
    float mag, magerr;
    char cmaxmag[1];

    //convert array of NLL ArrivalDesc to SNAP phases (declared in grid_search.h)
    ArrivalsToPhases(parrivals, narrivals);


    // get calibration from file
    i = ReadCALfromFile(snap_file_SED_CAL, "", 0);
    if (i > 0) {
        nll_puterr2("ERROR: opening SED calibration file", snap_file_SED_CAL);
        return (1);
    }



    // create and set variables and structures needed by ETH routines
    o_time_mn = juliam4(&phypo->year, &phypo->month, &phypo->day, &phypo->hour, &phypo->min);


    // call ETH magnitude routine here
    magnitude(o_time_mn, &mag, &magerr, cmaxmag, &num_mag, o_event_type, phypo->depth);
    //printf("call mag() %d %f %f %c %d %c %lf\n", o_time_mn, mag, magerr, cmaxmag, num_mag, o_event_type, phypo->depth);


    // load magnitude results

    if (num_mag > 0) {

        // load results to NLL hypo structure
        phypo->amp_mag = mag;
        phypo->num_amp_mag = num_mag;
        /* SH 23072004
           added magerr to hypo structure */
        phypo->mag_err = magerr;

        // load results to NLL arrivals structure
        for (narr = 0; narr < narrivals && narr < MAXPHASE; narr++) {
            parr = parrivals + narr;
            pphs = phases + narr;
            parr->amp_mag = pphs->mag;
            ///??? = pphs->magrmk;
        }
    }

    return (0);

}






