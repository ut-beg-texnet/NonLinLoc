/*
 * Copyright (C) 1999-2015 Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
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


/*   NLLocOutput.c

        Hypocenter output-format writers for NonLinLoc: WriteHypoAlberto4,
        WriteHypoFmamp, WriteHypoEll, WriteHypo71 and WriteHypoInverseArchive.

        Extracted verbatim from NLLocLib.c during the modernization effort
        (Phase 3: modularization). Pure code movement, no logic change; all
        function prototypes remain in NLLocLib.h, so callers are unaffected.
 */


/* NLLocLib.h is not self-contained; these must be included before it,
   in the same order as NLLocLib.c */
#include "GridLib.h"
#include "ran1/ran1.h"
#include "velmod.h"
#include "GridMemLib.h"
#include "calc_crust_corr.h"
#include "phaseloclist.h"
#include "otime_limit.h"
#include "NLLocLib.h"
#include "json_io.h"

#ifdef CUSTOM_ETH
#include "custom_eth/eth_functions.h"
#endif


/** function to write hypocenter and arrivals to file (Alberto 4 SIMULPS) */

int WriteHypoAlberto4(FILE *fpio, HypoDesc* phypo, ArrivalDesc* parrivals, int narrivals, char* filename) {

    int ifile = 0, narr;
    char fname[FILENAME_MAX];
    double mag;
    ArrivalDesc* parr;
    int nlat, nlon;


    /* set hypocenter parameters */
    if (phypo->amp_mag != MAGNITUDE_NULL)
        mag = phypo->amp_mag;
    else if (phypo->dur_mag != MAGNITUDE_NULL)
        mag = phypo->dur_mag;
    else
        mag = 0.0;


    /* write hypocenter to file */

    if (fpio == NULL) {
        snprintf(fname, sizeof(fname), "%s.loc.sim", filename);
        if ((fpio = fopen(fname, "w")) == NULL) {
            nll_puterr("ERROR: opening Alberto 4 hypocenter output file.");
            return (-1);
        } else {
            NumFilesOpen++;
        }
        ifile = 1;
    }


    /* write hypocenter parameters */

    //87 1 1  1 2  0.00 35N14.26 120W44.00   5.15   0.00
    nlat = (int) fabs(phypo->dlat);
    nlon = (int) fabs(phypo->dlong);
    fprintf(fpio,
            "%2.2d%2.2d%2.2d %2.2d%2.2d%6.2f %2.2d%c%5.2f %3.3d%c%5.2f %6.2f %6.2f",
            phypo->year % 100, phypo->month, phypo->day,
            phypo->hour, phypo->min, phypo->sec,
            nlat, phypo->dlat > 0.0 ? 'N' : 'S',
            60.0 * (fabs(phypo->dlat) - (double) nlat),
            nlon, phypo->dlong > 0.0 ? 'E' : 'W',
            60.0 * (fabs(phypo->dlong) - (double) nlon),
            phypo->depth, mag
            );

    // write arrivals

    // 20130228 AJL - bug fix
    // for (narr = 0; narr < phypo->nreadings; narr++) {
    for (narr = 0; narr < narrivals; narr++) {
        if (narr % 5 == 0)
            fprintf(fpio, "\n");
        parr = parrivals + narr;
        fprintf(fpio, "%4s%1s%1s%2.2d%7.4f",
                parr->label,
                strcmp(parr->onset, ARRIVAL_NULL_STR) == 0 ? "i" : parr->onset,
                parr->phase,
                parr->min, parr->sec
                );
    }

    fprintf(fpio, "\n");

    if (ifile) {

        fclose(fpio);
        NumFilesOpen--;
    }

    return (0);

}

/** function to write hypocenter summary and arrivals to file (fmamp format)
 *
 * 20160920 AJL - created
 */
int WriteHypoFmamp(FILE *fpio, HypoDesc* phypo, ArrivalDesc* parrivals, int narrivals, char* filename, int write_header) {


    // write hypocenter to file
    int ifile = 0;
    char fname[FILENAME_MAX];
    if (fpio == NULL) {
        snprintf(fname, sizeof(fname), "%s.loc.fmamp", filename);
        if ((fpio = fopen(fname, "w")) == NULL) {
            nll_puterr("ERROR: opening hypocenter output file.");
            return (-1);
        } else {
            NumFilesOpen++;
        }
        ifile = 1;
    }


    // write hypocenter parameters

    /* fmamp Hypocenter
        typedef struct {
           long unique_id;
           int year, month, day, hour, min; // origin time
           double dec_sec;
           double rms;
           double lat;
           double lon;
           double errh;
           double depth;
           double errz;
           int nassoc_P; // number of picks that contributed with weight > 0 to location
           double dist_min; // minimum distance of associated phase counted as nassoc_P
           double dist_max; // maximum distance of associated phase counted as nassoc_P
           double gap_primary; // maximum azimuth gap
           double gap_secondary; // secondary azimuth gap - largest azumth gap filled by a single station
           double ampAttenPower; // attenuation decay power from linear regression fit of P amplitudes to distance
           double magnitude; // 1/n_stations * distance of sum of vectors from epicenter to stations
           char mag_type[16];
           // picks
           Pick pick_list[MAX_NUM_PICKS];
           int pick_list_size;
        } Hypocenter;
     */
    if (write_header) {
        fprintf(fpio,
                "event_unique_id year month day hour min dec_sec rms lat lon errh depth errz "
                "nassoc_P dist_min dist_max gap_primary gap_secondary ampAttenPower magnitude mag_type\n");
        fprintf(fpio,
                "event_unique_id station location channel network "
                "phase "
                "year month day hour min dec_sec " // pick time
                "pick_error pick_error_type "
                "residual "
                "fmpolarity fmquality fmtype "
                "amplitude "
                "take_off_angle_az " // degrees CW from N
                "take_off_angle_inc " // degrees (0/down->180/up)
                "epicentral_distance " // degrees
                "epicentral_azimuth " // degrees CW from N
                "\n");
    }
    fprintf(fpio, "\n"); // event separator
    char event_unique_id[64];
    snprintf(event_unique_id, sizeof(event_unique_id), "%4.4d%2.2d%2.2d%2.2d%2.2d%5.5d", phypo->year, phypo->month, phypo->day, phypo->hour, phypo->min, (int) (phypo->sec * 1000.0)); // unique_id
    fprintf(fpio, "%s ", event_unique_id); // unique_id
    fprintf(fpio, "%4.4d %2.2d %2.2d %2.2d %2.2d %8.4f %f ", phypo->year, phypo->month, phypo->day, phypo->hour, phypo->min, phypo->sec, phypo->rms); // year month day hour min dec_sec rms
    double errz = -1.0;
    if (phypo->cov.zz > FLT_MIN) {
        errz = sqrt(phypo->cov.zz);
    }
    fprintf(fpio, "%f %f %f %f %f ", phypo->dlat, phypo->dlong, phypo->ellipse.len2, phypo->depth, errz); // lat lon errh depth errz
    double dist_max = -1.0;
    fprintf(fpio, "%d %f %f %f %f ", phypo->associatedPhaseCount, phypo->dist, dist_max, phypo->gap, phypo->gap_secondary); // nassoc_P dist_min dist_max gap_primary gap_secondary
    double atten = -999.0;
    double mag = 0.0;
    char mag_type[] = "NA";
    if (phypo->amp_mag != MAGNITUDE_NULL) {
        mag = phypo->amp_mag;
        snprintf(mag_type, sizeof(mag_type), "%s", "ML");
    } else if (phypo->dur_mag != MAGNITUDE_NULL) {
        mag = phypo->dur_mag;
        snprintf(mag_type, sizeof(mag_type), "%s", "MD");
    }
    fprintf(fpio, "%f %f %s ", atten, mag, mag_type); // ampAttenPower magnitude mag_type
    fprintf(fpio, "\n"); // event separator


    // write arrival parameters

    int write_arrivals = 1;
    if (write_arrivals) {

        /* fmamp Pick
            typedef struct {
                long event_unique_id;
                char station[16];
                char location[16];
                char channel[16];
                char network[16];
                char phase[16];
                int year, month, day, hour, min; // pick time
                double dec_sec;
                char pick_error_type[32];
                double pick_error;
                int fmpolarity; // broad-band polarity measure
                double fmquality; // broad-band polarity measure
                char fmtype[32];
                double amplitude; // amplitude
                double take_off_angle_inc; // degrees (0/down->180/up)
                double take_off_angle_az; // degrees CW from N
                double epicentral_distance; // degrees
                double epicentral_azimuth; // degrees CW from N
                double residual;
                double weight; // calculated first motion obs quality weight
                double take_off_angle_distrib_weight; // calculated take-off angle distribution weight on sphere
            } Pick;
         */

        char loc[] = "--";
        char fmtype[] = "F";
        double amplitude = -1.0;
        ArrivalDesc* parr;
        for (int narr = 0; narr < narrivals; narr++) {
            parr = parrivals + narr;
            if (parr->ray_qual < iAngleQualityMin || parr->first_mot_quality < FLT_MIN) {
                continue;
            }
            fprintf(fpio, "%s ", event_unique_id); // unique_id
            fprintf(fpio, "%s %s %s %s%s ", parr->label, loc, parr->network, parr->inst, parr->comp);
            fprintf(fpio, "%s ", parr->phase);
            fprintf(fpio, "%4.2d %2.2d %2.2d %2.2d %2.2d %8.4f ", parr->year, parr->month, parr->day, parr->hour, parr->min, parr->sec); // year month day hour min dec_sec
            fprintf(fpio, "%f %s ", parr->error, parr->error_type);
            fprintf(fpio, "%f ", parr->residual);
            fprintf(fpio, "%s %f %s ", parr->first_mot, parr->first_mot_quality, fmtype);
            fprintf(fpio, "%f ", amplitude);
            fprintf(fpio, "%f %f ", rect2latlonAngle(0, parr->ray_azim), parr->ray_dip);
            fprintf(fpio, "%f %f ", parr->dist, rect2latlonAngle(0, parr->azim));
            fprintf(fpio, "\n");
        }

    }


    if (ifile) {

        fclose(fpio);
        NumFilesOpen--;
    }

    return (0);

}

/** function to write hypocenter summary to file (quasi HypoEllipse format) */

int WriteHypoEll(FILE *fpio, HypoDesc* phypo, ArrivalDesc* parrivals, int narrivals,
        char* filename, int write_header, int write_arrivals) {

    int ifile = 0, narr;
    char fname[FILENAME_MAX];
    double mag;
    ArrivalDesc* parr;
    double tpobs, resid;


    /* set hypocenter parameters */
    if (phypo->amp_mag != MAGNITUDE_NULL)
        mag = phypo->amp_mag;
    else if (phypo->dur_mag != MAGNITUDE_NULL)
        mag = phypo->dur_mag;
    else
        mag = 0.0;


    /* write hypocenter to file */

    if (fpio == NULL) {
        snprintf(fname, sizeof(fname), "%s.loc.hypo_ell", filename);
        if ((fpio = fopen(fname, "w")) == NULL) {
            nll_puterr("ERROR: opening hypocenter output file.");
            return (-1);
        } else {
            NumFilesOpen++;
        }
        ifile = 1;
    }


    /* write hypocenter parameters */

    if (write_header) {
        fprintf(fpio,
                "DATE     ORIGIN     LAT         LONG         DEPTH   ");
        fprintf(fpio,
                "MAG  NO  GAP D1     RMS   ");
        fprintf(fpio,
                "AZ1  DIP1 SE1    AZ2  DIP2 SE2    SE3    \n");
        /*"ERH  ERZ Q SQD  ADJ IN NR  AVR  AAR NM AVXM SDXM NF AVFM SDFM I\n");*/
    }
    fprintf(fpio,
            "%4.4d%2.2d%2.2d %2.2d%2.2d %5.2lf %3d %1c %5.2lf %4d %1c %5.2lf %7.3lf ",
            phypo->year, phypo->month, phypo->day,
            phypo->hour, phypo->min, phypo->sec,
            (int) fabs(phypo->dlat), (phypo->dlat >= 0.0 ? 'N' : 'S'),
            (fabs(phypo->dlat) - (int) fabs(phypo->dlat)) * 60.0,
            (int) fabs(phypo->dlong), (phypo->dlong >= 0.0 ? 'E' : 'W'),
            (fabs(phypo->dlong) - (int) fabs(phypo->dlong)) * 60.0,
            phypo->depth);
    fprintf(fpio, "%4.2lf %3d %3d %6.2lf %5.2lf ",
            mag, phypo->nreadings, (int) (0.5 + phypo->gap), phypo->dist, phypo->rms);
    fprintf(fpio, "%4d %4d %6.2lf %4d %4d %6.2lf %6.2lf ",
            (int) (0.5 + phypo->ellipsoid.az1),
            (int) (0.5 + phypo->ellipsoid.dip1),
            phypo->ellipsoid.len1,
            (int) (0.5 + phypo->ellipsoid.az2),
            (int) (0.5 + phypo->ellipsoid.dip2),
            phypo->ellipsoid.len2,
            phypo->ellipsoid.len3);
    fprintf(fpio, "\n");



    if (write_arrivals) {

        fprintf(fpio, "\n");


        fprintf(fpio,
                "  STN  DIST AZM AIN PRMK HRMN P-SEC TPOBS TPCAL DLY/H1 P-RES P-WT AMX PRX CALX K XMAG RMK FMP FMAG\n");
        /*"  STN  DIST AZM AIN PRMK HRMN P-SEC TPOBS TPCAL DLY/H1 P-RES P-WT AMX PRX CALX K XMAG RMK FMP FMAG SRMK S-SEC TSOBS S-RES  S-WT    DT\n"); */

        // 20130228 AJL - bug fix
        // for (narr = 0; narr < phypo->nreadings; narr++) {
        for (narr = 0; narr < narrivals; narr++) {
            parr = parrivals + narr;
            tpobs = parr->obs_travel_time > -9.99 ? parr->obs_travel_time : 0.0;
            resid = parr->residual > -99.99 ? parr->residual : -99.99;
            fprintf(fpio,
                    "%5s %5.1lf %3d %3d %2s%1s%1d %2.2d%2.2d %5.2lf %5.2lf %5.2lf       %-6.2lf %5.2lf\n",
                    parr->label, parr->dist,
                    (int) (0.5 + rect2latlonAngle(0, parr->ray_azim)),
                    (int) (0.5 + parr->ray_dip),
                    parr->phase, parr->first_mot, parr->quality,
                    parr->hour, parr->min, parr->sec,
                    tpobs, parr->pred_travel_time,
                    resid, parr->weight
                    );
        }

    }


    if (ifile) {

        fclose(fpio);
        NumFilesOpen--;
    }

    return (0);

}

/** function to write hypocenter summary to file (HYPO71 format) */

int WriteHypo71(FILE *fpio, HypoDesc* phypo, ArrivalDesc* parrivals, int narrivals,
        char* filename, int write_header, int write_arrivals) {

    int ifile = 0, narr;
    char fname[FILENAME_MAX];
    ArrivalDesc* parr;
    double tpobs, resid;
    double mag, erh, erz;
    char qualS, qualD, qual;
    int pha_qual;
    double xmag, fmag;

    /* set hypocenter parameters */
    if (phypo->amp_mag != MAGNITUDE_NULL)
        mag = phypo->amp_mag;
    else if (phypo->dur_mag != MAGNITUDE_NULL)
        mag = phypo->dur_mag;
    else
        mag = 0.0;

    /* write hypocenter to file */

    if (fpio == NULL) {
        snprintf(fname, sizeof(fname), "%s.loc.h71", filename);
        if ((fpio = fopen(fname, "w")) == NULL) {
            nll_puterr("ERROR: opening hypocenter output file.");
            return (-1);
        } else {
            NumFilesOpen++;
        }
        ifile = 1;
    }


    /* write hypocenter parameters */

    if (write_header) {
        fprintf(fpio,
                "  DATE    ORIGIN    LAT      LONG      DEPTH    ");
        fprintf(fpio,
                "MAG NO DM GAP M  RMS  ERH  ERZ Q SQD  ADJ IN NR  AVR  AAR NM AVXM SDXM NF AVFM SDFM I\n");
    }

    fprintf(fpio,
            " %2.2d%2.2d%2.2d %2.2d%2.2d %5.2lf%3d %5.2lf%4d %5.2lf %6.2lf",
            phypo->year % 100, phypo->month, phypo->day,
            phypo->hour, phypo->min, phypo->sec,
            (int) phypo->dlat, /*(phypo->dlat >= 0.0 ? 'N' : 'S'),*/
            (phypo->dlat - (int) phypo->dlat) * 60.0,
            (int) phypo->dlong, /*(phypo->dlong >= 0.0 ? 'E' : 'W'),*/
            (phypo->dlong - (int) phypo->dlong) * 60.0,
            phypo->depth);

    fprintf(fpio, " %6.2lf%3d%3d %3d 0%5.2lf",
            mag, phypo->nreadings, (int) (0.5 + phypo->dist),
            (int) (0.5 + phypo->gap), phypo->rms);

    // 20100204 AJL Satriano Bug Fix
    //erh = sqrt(phypo->cov.xx * phypo->cov.xx
    //        + phypo->cov.yy * phypo->cov.yy);
    //erz = sqrt(phypo->cov.zz * phypo->cov.zz);
    erh = sqrt(phypo->cov.xx + phypo->cov.yy);
    erz = sqrt(phypo->cov.zz);
    // End - 20100204
    fprintf(fpio, "%5.1lf%5.1lf", erh, erz);

    /* ABCD quality levels (from HYPO71 open file report 75-311, p. 27) */
    if (phypo->rms < 0.15 && erh <= 1.0 && erz <= 2.0)
        qualS = 'A';
    else if (phypo->rms < 0.30 && erh <= 2.5 && erz <= 5.0)
        qualS = 'B';
    else if (phypo->rms < 0.50 && erh <= 5.0)
        qualS = 'C';
    else
        qualS = 'D';
    if (phypo->nreadings >= 6 && phypo->gap <= 90
            && (phypo->dist <= phypo->depth || phypo->dist <= 5.0))
        qualD = 'A';
    else if (phypo->nreadings >= 6 && phypo->gap <= 135
            && (phypo->dist <= 2.0 * phypo->depth
            || phypo->dist <= 10.0))
        qualD = 'B';
    else if (phypo->nreadings >= 6 && phypo->gap <= 180
            && phypo->dist <= 50.0)
        qualD = 'C';
    else
        qualD = 'D';
    /*if (abs(qualS - qualD) == 1)
    qual = qualS < qualD ? qualS : qualD;
    else*/
    qual = 1 + (qualS + qualD) / 2;
    fprintf(fpio, " %1c %1c %1c", qual, qualS, qualD);

    /* dummy values for remaining fields */
    fprintf(fpio,
            " %4.2lf %2d %2d-%4.2lf %4.2lf %2d %4.1lf %4.1lf %2d %4.1lf %4.1lf%2d\n",
            0.0, 0, 0, 0.0, 0.0, 0, 0.0, 0.0, 0, 0.0, 0.0, 0);



    if (write_arrivals) {

        fprintf(fpio, "\n");


        fprintf(fpio,
                "  STN  DIST AZM AIN PRMK HRMN P-SEC TPOBS TPCAL DLY/H1 P-RES P-WT AMX PRX CALX K XMAG RMK FMP FMAG SRMK S-SEC TSOBS S-RES  S-WT    DT\n");

        // 20130228 AJL - bug fix
        // for (narr = 0; narr < phypo->nreadings; narr++) {
        for (narr = 0; narr < narrivals; narr++) {
            parr = parrivals + narr;
            pha_qual = (parr->quality >= 0 && parr->quality <= 4) ?
                    parr->quality : Err2Qual(parr);
            if (pha_qual < 0)
                pha_qual = 4;
            tpobs = parr->obs_travel_time > -9.99 ?
                    parr->obs_travel_time : 0.0;
            resid = parr->residual > -99.99 ? parr->residual : -99.99;
            fprintf(fpio,
                    /*" %-4s %5.1lf %3d %3d %2s%1s%1d %2.2d%2.2d %5.2lf %5.2lf %5.2lf  0.00 %6.3lf %5.2lf\n", */
                    " %-4s %5.1lf %3d %3d %2s%1s%1d %2.2d%2.2d %5.2lf %5.2lf %5.2lf  0.00 %-6.2lf %4.2lf",
                    parr->label, parr->dist,
                    (int) (0.5 + rect2latlonAngle(0, parr->ray_azim)),
                    (int) (0.5 + parr->ray_dip),
                    parr->phase, parr->first_mot, pha_qual,
                    parr->hour, parr->min, parr->sec,
                    tpobs, parr->pred_travel_time,
                    resid, parr->weight
                    );

            /* set magnitudes */
            xmag = (parr->amp_mag != MAGNITUDE_NULL) ? parr->amp_mag : 0.0;
            fmag = (parr->dur_mag != MAGNITUDE_NULL) ? parr->dur_mag : 0.0;
            fprintf(fpio,
                    " 0.0 0.0 0.00 0 %3.2lf 000 00.0 %3.2lf ??4 00.00 00.00 00.00   0.0      \n",
                    xmag, fmag
                    );
        }

    }


    if (ifile) {

        fclose(fpio);
        NumFilesOpen--;
    }

    return (0);

}

/** function to write hypocenter summary to file (HYPOINVERSE archive format) */

int WriteHypoInverseArchive(FILE *fpio, HypoDesc *phypo, ArrivalDesc *parrivals, int narrivals,
        char *filename, int writeY2000, int write_arrivals, double arrivalWeightMax) {

    int ifile = 0, narr;
    char fname[FILENAME_MAX];
    ArrivalDesc *parr, *psarr;
    double rms, resid, amplitude;
    double dtemp;
    char chrtmp[MAXSTRING];

    char first_mot;

    int haveSumHdr = 0;

    /* hypocenter fields */
    char loc_remark[] = "NLL", aux_remark[] = " ", aux_remark_prog[] = " ";
    int num_S_wt = 0, num_P_fmot = 0;
    double amp_mag = 0.0, dur_mag = 0.0;
    double amp_mag_wt = 0.0, dur_mag_wt = 0.0;
    double err_horiz = 0.0, err_vert = 0.0;
    int az_err_prin = 0.0, az_err_inter = 0.0;
    int dip_err_prin = 0.0, dip_err_inter = 0.0;

    /* arrival fields */
    char sta_remark[] = " ";
    int amp_mag_wt_code = 0, dur_mag_wt_code = 0;


    /* if Y2000 and saved HypoInverseArchiveSumHdr, extract hypo fields */
    /* see ftp://ehzftp.wr.usgs.gov/klein/hyp2000/docs/hyp2000-1.0.htm#_Toc7234741 */
    if (writeY2000 && strlen(HypoInverseArchiveSumHdr) > 0) {
        haveSumHdr = 1;
        //printf("HypoInverseArchiveSumHdr:\n|%s|\n", HypoInverseArchiveSumHdr);
        // 81 1 A1 Auxiliary remark from analyst (i.e. Q for quarry).
        // 20240129 AJL - bug fix //strncpy(aux_remark, HypoInverseArchiveSumHdr + 80, 1);
        snprintf(aux_remark, 2, "%1s", HypoInverseArchiveSumHdr + 80);
    }

    /* if Y2000, intialize fields */
    if (writeY2000) {
        for (narr = 0; narr < narrivals; narr++) {
            parr = parrivals + narr;
            // num_S_wt - 83 3 I3 Number of S times with weights greater than 0.1.
            if (IsPhaseID(parr->phase, "S") && parr->weight > 0.1)
                num_S_wt++;
            // 94 3 I3 Number of P first motions.
            if (IsPhaseID(parr->phase, "P")
                    && parr->first_mot[0] != ' '
                    && parr->first_mot[0] != ARRIVAL_NULL_CHR)
                num_P_fmot++;
        }
    }


    /* set hypocenter parameters */
    //printf("Ma %4.2f  Md %4.2f  MAGNITUDE_NULL %4.2f\n", phypo->amp_mag, phypo->dur_mag, MAGNITUDE_NULL);
    if (fabs(phypo->amp_mag - MAGNITUDE_NULL) > 0.01) {
        amp_mag = phypo->amp_mag;
        // 97 4 F4.1 Total of amplitude mag weights ~number of readings.*
        amp_mag_wt = (double) phypo->num_amp_mag;
    } else if (haveSumHdr) {
        // 37 3 F3.2 Amplitude magnitude.
        // JMS 20210304 - bug fix: added checks for no mag
        strncpy(chrtmp, HypoInverseArchiveSumHdr + 36, 3);
        chrtmp[3] = '\0';
        if (sscanf(chrtmp, "%3lf", &amp_mag) != 0) {
            amp_mag /= 100.0;
        } else {
            amp_mag = 0.0;
        }
    } else {
        amp_mag = 0.0;
    }
    if (fabs(phypo->dur_mag - MAGNITUDE_NULL) > 0.01) {
        dur_mag = phypo->dur_mag;
        // 101 4 F4.1 Total of duration mag weights ~number of readings. *
        dur_mag_wt = (double) phypo->num_dur_mag;
    } else if (haveSumHdr) {
        // 71 3 F3.2 Coda duration magnitude.
        // JMS 20210304 - bug fix: added checks for no mag
        strncpy(chrtmp, HypoInverseArchiveSumHdr + 70, 3);
        chrtmp[3] = '\0';
        if (sscanf(chrtmp, "%3lf", &dur_mag) != 0) {
            dur_mag /= 100.0;
        } else {
            dur_mag = 0.0;
        }
    } else {
        dur_mag = 0.0;
    }
    //printf("-> amp_mag %4.2f  dur_mag %4.2f\n", amp_mag, dur_mag);



    /* write hypocenter to file */

    if (fpio == NULL) {
        snprintf(fname, sizeof(fname), "%s.loc.hypo_inv", filename);
        if ((fpio = fopen(fname, "w")) == NULL) {
            nll_puterr("ERROR: opening hypocenter output file.");
            return (-1);
        } else {
            NumFilesOpen++;
        }
        ifile = 1;
    }


    /* write hypocenter parameters */

    // 1 4 I4 Year. *
    if (writeY2000)
        fprintf(fpio, "%4.4d", phypo->year);
    else
        fprintf(fpio, "%2.2d", phypo->year % 100);
    // 5 8 4I2 Month, day, hour and minute.
    // 13 4 F4.2 origin time seconds.
    // 17 2 F2.0 Latitude (deg). First character must not be blank.
    // 19 1 A1 S for south, blank otherwise.
    // 20 4 F4.2 Latitude (min).
    // 24 3 F3.0 Longitude (deg).
    // 27 1 A1 E for east, blank otherwise.
    // 28 4 F4.2 Longitude (min).
    // 32 5 F5.2 Depth (km).
    fprintf(fpio,
            "%2.2d%2.2d%2.2d%2.2d%4.0lf%2.2d%1c%4.0lf%3.3d%1c%4.0lf%5.0lf",
            phypo->month, phypo->day,
            phypo->hour, phypo->min, 100.0 * phypo->sec,
            (int) fabs(phypo->dlat), (phypo->dlat >= 0.0 ? ' ' : 'S'),
            (fabs(phypo->dlat) - (int) fabs(phypo->dlat)) * 6000.0,
            (int) fabs(phypo->dlong), (phypo->dlong >= 0.0 ? 'E' : ' '),
            (fabs(phypo->dlong) - (int) fabs(phypo->dlong)) * 6000.0,
            100.0 * phypo->depth);
    // 37 3 F3.2 Amplitude magnitude. *
    // AJL 20080304 - bug fix: added checks for MAGNITUDE_NULL
    // JMS 20210304 - improvement: seems unecessary now, check performed earlier at mag reading
    // JMS 20210304 - bug fix: if no magnitude, print spaces
    if (amp_mag == 0.0)
        fprintf(fpio, "%s", writeY2000 ? "   " : "  ");
    else {
        if (writeY2000)
            fprintf(fpio, "%3.0lf", 100.0 * amp_mag);
        else
            fprintf(fpio, "%2.0lf", 10.0 * amp_mag);
    }
    // 40 3 I3 Number of P & S times with final weights greater than 0.1.
    // 43 3 I3 Maximum azimuthal gap, degrees.
    // 46 3 F3.0 Distance to nearest station (km).
    // 49 4 F4.2 RMS travel time residual.
    rms = phypo->rms < 99.99 ? phypo->rms : 99.99;
    fprintf(fpio, "%3.3d%3.3d%3.0lf%4.0lf",
            phypo->nreadings, (int) (0.5 + phypo->gap), phypo->dist, 100.0 * rms);
    // 53 3 F3.0 Azimuth of largest principal error (deg E of N).
    // 56 2 F2.0 Dip of largest principal error (deg).
    // TODO this is smallest princiapl error !!
    az_err_prin = (int) (0.5 + phypo->ellipsoid.az1);
    dip_err_prin = (int) (0.5 + phypo->ellipsoid.dip1);
    if (dip_err_prin < 0) {
        dip_err_prin += 90;
        az_err_prin += 180;
        if (az_err_prin > 360)
            az_err_prin -= 360;
    }
    az_err_inter = (int) (0.5 + phypo->ellipsoid.az2);
    dip_err_inter = (int) (0.5 + phypo->ellipsoid.dip2);
    if (dip_err_inter < 0) {
        dip_err_inter += 90;
        az_err_inter += 180;
        if (az_err_inter > 360)
            az_err_inter -= 360;
    }
    // 58 4 F4.2 Size of largest principal error (km).
    // 62 3 F3.0 Azimuth of intermediate principal error.
    // 65 2 F2.0 Dip of intermediate principal error.
    // 67 4 F4.2 Size of intermediate principal error (km).
    /*fprintf(fpio, "%3.3d%2.2d%4.0lf%3.3d%2.2d%4.0lf",
    az_err_prin, dip_err_prin,
    100.0 * phypo->ellipsoid.len1,
    az_err_inter, dip_err_inter,
    100.0 * phypo->ellipsoid.len2);
     */
    // AJL 20060829
    // AJL 20061124
    fprintf(fpio, "000000000000000000");

    // 71 3 F3.2 Coda duration magnitude. *
    // JMS 20210304 - bug fix: if no magnitude, print spaces
    if (dur_mag == 0.0)
        fprintf(fpio, "%s", writeY2000 ? "   " : "  ");
    else {
        if (writeY2000)
            fprintf(fpio, "%3.0lf", 100.0 * dur_mag);
        else
            fprintf(fpio, "%2.0lf", 10.0 * dur_mag);
    }
    // 74 3 A3 Event location remark (region), derived from location.
    // 77 4 F4.2 Size of smallest principal error (km).
    // 81 1 A1 Auxiliary remark from analyst (i.e. Q for quarry).
    // 82 1 A1 Auxiliary remark from program (i.e. “-“ for depth fixed, etc.).
    // 83 3 I3 Number of S times with weights greater than 0.1.
    // 86 4 F4.2 Horizontal error (km).
    // 90 4 F4.2 Vertical error (km).
    if (writeY2000) {
        // TODO this is largest principal error !!
        dtemp = 100.0 * phypo->ellipsoid.len3;
        fprintf(fpio, "%3.3s%4.0lf", loc_remark, dtemp < 9999.0 ? dtemp : 9999.0);
        fprintf(fpio, "%1.1s%1.1s%3.3d%4.0lf%4.0lf",
                aux_remark, aux_remark_prog, num_S_wt, err_horiz, err_vert);
    } else {
        dtemp = 100.0 * phypo->ellipsoid.len3;
        fprintf(fpio, "%3.3s%4.0lf", loc_remark, dtemp < 9999.0 ? dtemp : 9999.0);
        fprintf(fpio, "%1.1s%1.1s%2.2d%4.0lf%4.0lf",
                aux_remark, aux_remark_prog, num_S_wt, err_horiz, err_vert);
    }


    if (writeY2000) {
        // 94 3 I3 Number of P first motions. *
        // 97 4 F4.1 Total of amplitude mag weights ~number of readings.*
        // 101 4 F4.1 Total of duration mag weights ~number of readings. *
        fprintf(fpio, "%3.3d%4.0lf%4.0lf", num_P_fmot, 10.0 * amp_mag_wt, 10.0 * dur_mag_wt);
        // TODO - not iniitalized from NLL results
        // 105 3 F3.2 Median-absolute-difference of amplitude magnitudes.
        // 108 3 F3.2 Median-absolute-difference of duration magnitudes.
        // 111 3 A3 3-letter code of crust and delay model.
        fprintf(fpio, "  0  0   ");
        // 114 1 A1 Authority” code, i.e. what network furnished the information.  Hypoinverse passes this code through.
        // 115 1 A1 Most common P & S data source code. (See table 1 below).
        // 116 1 A1 Most common duration data source code. (See cols. 68-69)
        // 117 1 A1 Most common amplitude data source code.
        if (haveSumHdr)
            fprintf(fpio, "%4.4s", HypoInverseArchiveSumHdr + 113);
        else
            fprintf(fpio, "   ");
        // TODO - not iniitalized from NLL results
        // 118 1 A1 Primary coda duration magnitude type code
        fprintf(fpio, "D");
        // 119 3 I3 Number of valid P & S readings (assigned weight > 0)
        fprintf(fpio, "%3.3d", phypo->nreadings);
        // TODO - not iniitalized from NLL results
        // 122 1 A1 Primary amplitude magnitude type code
        fprintf(fpio, "X");
        // 123 1    A1  "External" magnitude label or type code. Typically "L" (=ML) This information is not computed by Hypoinverse, but passed along, as computed by UCB.
        // 124 3 F3.2 "External" magnitude.
        // 127 3 F3.1 Total of "external" magnitude weights (~ number of readings).
        // 130 1 A1  Alternate amplitude magnitude label or type code.
        // 131 3 F3.2 Alternate amplitude magnitude.
        // 134 3 F3.1 Total of the alternate amplitude mag weights ~no. of readings.
        // 137 10 I10  Event identification number
        if (haveSumHdr)
            fprintf(fpio, "%24.24s", HypoInverseArchiveSumHdr + 122);
        else
            fprintf(fpio, "%24.24s", "");
        // 147 1 A1 Preferred magnitude label code chosen from those available.
        // 148 3 F3.2 Preferred magnitude, chosen by the Hypoinverse PRE command.
        // 151 4 F4.1 Total of the preferred mag weights (~ number of readings).
        if (phypo->num_amp_mag > 0 && phypo->num_amp_mag >= phypo->num_dur_mag) {
            fprintf(fpio, "X%3.0lf%4.0lf", 100.0 * amp_mag, 10.0 * amp_mag_wt);
        } else if (phypo->num_dur_mag > 0 && phypo->num_dur_mag >= phypo->num_amp_mag) {
            fprintf(fpio, "D%3.0lf%4.0lf", 100.0 * dur_mag, 10.0 * dur_mag_wt);
        } else {
            if (haveSumHdr)
                fprintf(fpio, "%8.8s", HypoInverseArchiveSumHdr + 146);
            else
                fprintf(fpio, "%8.8s", "");
        }
        // 155 1 A1 Alternate coda duration magnitude label or type code.
        // 156 3 F3.2 Alternate coda duration magnitude.
        // 159 4 F4.1 Total of the alternate coda duration magnitude weights. *
        // 163 1 A1 Version” of the information, i.e. the stage of processing.  This can either be passed through, or assigned by Hypoinverse with the LAB command.
        // 164 1 A1  Version of last human review. Hypoinverse passes this through.
        if (haveSumHdr)
            fprintf(fpio, "%10.10s", HypoInverseArchiveSumHdr + 154);
        else
            fprintf(fpio, "%10.10s", "");
        // 164 is the last filled column
    } else {
        fprintf(fpio, "%2.2d", num_P_fmot);
        fprintf(fpio, "                   D   X");
    }

    fprintf(fpio, "\n");



    if (write_arrivals) {

        // following format commens are for Y2000 (station) archive format

        // AJL 21JUN2000 bug fix (not all phases are written, looses fm's from phases not used for misfit!
        //	for (narr = 0; narr < phypo->nreadings; narr++) {
        for (narr = 0; narr < narrivals; narr++) {
            parr = parrivals + narr;

            /* skip non-P phases */
            if (!IsPhaseID(parr->phase, "P"))
                continue;

            /* check for following S arrival at same station */
            psarr = NULL;
            // AJL 21JUN2000
            //		if (narr + 1 < phypo->nreadings) {
            if (narr + 1 < narrivals) {
                psarr = parr + 1;
                if (!IsPhaseID(psarr->phase, "S")
                        || strcmp(parr->label, psarr->label) != 0)
                    psarr = NULL;
            }

            // 20060619 AJL
            /* no first motion for phases with take-off angle quality < min qual */
            if (angleMode == ANGLE_MODE_YES && parr->ray_qual < iAngleQualityMin)
                first_mot = ' ';
            else if (strpbrk(parr->first_mot, "cCuU+"))
                first_mot = 'U';
            else if (strpbrk(parr->first_mot, "dD-"))
                first_mot = 'D';
            else
                first_mot = ' ';

            if (writeY2000) {
                // 1 5 A5 5-letter station site code, left justified. *
                // 6 2 A2 2-letter seismic network code. *
                // 8 1 1X Blank *
                // 9 1 A1 One letter station component code.
                // 10 3 A3 3-letter station component code. *
                // 13 1 1X Blank *
                /*  BUG FIX - 20081029 signaled by Roman Racine (racine@sed.ethz.ch)
                fprintf(fpio, "%-5.5s%-2.2s%1.1s%1.1s%-3.3s%1.1s",
                        parr->label,
                        parr->network != ARRIVAL_NULL_STR ? parr->network : "",
                        " ",
                        parr->comp != ARRIVAL_NULL_STR ? parr->comp : "",
                        parr->inst != ARRIVAL_NULL_STR ? parr->inst : "",
                        " ");
                 */
                fprintf(fpio, "%-5.5s%-2.2s%1.1s%1.1s%-3.3s%1.1s",
                        parr->label,
                        strcmp(parr->network, ARRIVAL_NULL_STR) ? parr->network : "",
                        " ",
                        strcmp(parr->comp, ARRIVAL_NULL_STR) ? parr->comp : "",
                        strcmp(parr->inst, ARRIVAL_NULL_STR) ? parr->inst : "",
                        " ");
                // END - BUG FIX - 20081029
                //
                // 14 2 A2 P remark such as "IP".
                // try and reconstruct P remark
                // INGV - "Pn"
                // NCAL - "IP"
                if (strcmp(parr->onset, ARRIVAL_NULL_STR) == 0) { // NCAL form
                    fprintf(fpio, "%1.1s%1.1s", "", parr->phase);
                } else if (strcmp(parr->onset, ARRIVAL_NULL_STR) == 0
                        || strpbrk(parr->onset, " iIeE")) { // NCAL form
                    fprintf(fpio, "%1.1s%1.1s", parr->onset, parr->phase);
                } else { // INGV form
                    fprintf(fpio, "%-2.2s", parr->phase);
                }
                // 16 1 A1 P first motion.
                // 17 1 I1 Assigned P weight code.
                fprintf(fpio, "%1.1s%1d", parr->first_mot, parr->quality);
            } else {
                fprintf(fpio, "%4.4s%1.1s%1.1s%1c%1d%1.1s",
                        parr->label, " ", parr->phase, first_mot,
                        parr->quality > 9 ? 9 : parr->quality, " ");
            }
            //
            // 18 4 I4 Year. *
            // 22 8 4I2 Month, day, hour and minute.
            // 30 5 F5.2 Second of P arrival.
            if (writeY2000) {
                fprintf(fpio,
                        "%4.4d%2.2d%2.2d%2.2d%2.2d%5.0lf",
                        parr->year, parr->month, parr->day,
                        parr->hour, parr->min, 100.0 * parr->sec);
            } else {
                fprintf(fpio,
                        "%2.2d%2.2d%2.2d%2.2d%2.2d%5.0lf",
                        parr->year % 100, parr->month, parr->day,
                        parr->hour, parr->min, 100.0 * parr->sec);
            }
            //
            // 35 4 F4.2 P travel time residual.
            // 39 3 F3.2 P weight actually used.
            resid = parr->residual > -9.99 ? parr->residual : 0.0;
            resid = resid < 99.99 ? resid : 99.99;
            fprintf(fpio, "%4.0lf%3.0lf",
                    100.0 * resid, 100.0 * parr->weight);
            //
            // 42 5 F5.2 Second of S arrival.
            // 47 2 A2 S remark such as "ES".
            // 49 1 1X Blank
            // 50 1 I1 Assigned S weight code.
            // 51 4 F4.2 S travel time residual.
            //
            // 55 7 F7.2 Amplitude (Normally peak-to-peak). *
            // 62 2 I2 Amp units code. 0=PP mm, 1=0 to peak mm (UCB), 2=digital counts. *
            // 64 3 F3.2 S weight actually used.
            // 67 4 F4.2 P delay time.
            // 71 4 F4.2 S delay time.
            if (psarr != NULL) {
                resid = psarr->residual > -9.99 ? psarr->residual : 0.0;
                fprintf(fpio,
                        "%5.0lf%1s%1c %1d%4.0lf",
                        100.0 * psarr->sec, " ", psarr->phase[0],
                        psarr->quality, 100.0 * resid);
                if (writeY2000) {
                    fprintf(fpio, "%7.0lf%2d", 100.0 * parr->amplitude,
                            9); // TODO add true amp units code
                } else {
                    amplitude = parr->amplitude != AMPLITUDE_NULL
                            ? 100.0 * parr->amplitude : 0.0;
                    amplitude = amplitude < 99.9 ? amplitude : 99.9;
                    fprintf(fpio, "%3.0lf", amplitude);
                }
                fprintf(fpio,
                        "%3.0lf%4.0lf%4.0lf",
                        100.0 * psarr->weight,
                        100.0 * parr->delay, 100.0 * psarr->delay);
            } else {
                if (writeY2000) {
                    fprintf(fpio,
                            "    0   0   0%7.0lf%2d  0%4.0lf%4.0lf",
                            parr->amplitude != AMPLITUDE_NULL
                            ? 100.0 * parr->amplitude : 0.0,
                            9, 100.0 * parr->delay, 0.0);
                    // TODO add true amp units code
                } else {
                    amplitude = parr->amplitude != AMPLITUDE_NULL
                            ? 100.0 * parr->amplitude : 0.0;
                    amplitude = amplitude < 99.9 ? amplitude : 99.9;
                    fprintf(fpio,
                            "    0   0   0%3.0lf   %4.0lf%4.0lf", amplitude,
                            100.0 * parr->delay, 0.0);
                }
            }
            //
            // 75 4 F4.1 Epicentral distance (km).
            // 79 3 F3.0 Emergence angle at source.
            // 82 1 I1 Amplitude magnitude weight code.
            // 83 1 I1 Duration magnitude weight code.
            // 84 3 F3.2 Period at which the amplitude was measured for this station.
            // 87 1 A1 1-letter station remark.
            // 88 4 F4.0 Coda duration in seconds.
            // 92 3 F3.0 Azimuth to station in degrees E of N.
            fprintf(fpio,
                    "%4.0lf%3d%1d%1d%3.0lf%1s%4.0lf%3d",
                    GeometryMode == MODE_GLOBAL ? parr->dist * KM2DEG : 10.0 * parr->dist,
                    (int) (0.5 + parr->ray_dip),
                    amp_mag_wt_code, dur_mag_wt_code,
                    parr->period != PERIOD_NULL ? 100.0 * parr->period : 0.0,
                    sta_remark, parr->coda_dur,
                    (int) (0.5 + rect2latlonAngle(0, parr->ray_azim)));
            //
            // 95 3 F3.2 Duration magnitude for this station. *
            // 98 3 F3.2 Amplitude magnitude for this station. *
            if (writeY2000) {
                fprintf(fpio, "%3.0lf%3.0lf",
                        100.0 * parr->dur_mag,
                        fabs(parr->amp_mag - MAGNITUDE_NULL) > 0.01 ?
                        100.0 * parr->amp_mag : 0.0);
            } else {
                fprintf(fpio, "%2.0lf%2.0lf",
                        fabs(parr->dur_mag - MAGNITUDE_NULL) > 0.01 ?
                        10.0 * parr->dur_mag : 0.0,
                        fabs(parr->amp_mag - MAGNITUDE_NULL) > 0.01 ?
                        10.0 * parr->amp_mag : 0.0);
            }
            //
            // 101 4 F4.3 Importance of P arrival.
            // 105 4 F4.3 Importance of S arrival.
            fprintf(fpio,
                    "%4.0lf%4.0lf",
                    1000.0 * parr->weight / arrivalWeightMax,
                    psarr != NULL ? 1000.0 * psarr->weight / arrivalWeightMax : 0.0);
            //
            // 109 1 A1 Data source code.
            // 110 1 A1 Label code for duration magnitude from FC1 or FC2 command.
            // 111 1 A1 Label code for amplitude magnitude from XC1 or XC2 command.
            // 112 2 A2 2-letter station location code (component extension).
            // 113 is the last filled column.
            fprintf(fpio, "     "); // TODO - add data source?

            fprintf(fpio, "\n");
        }

    }


    fprintf(fpio, "\n");
    if (ifile) {

        fclose(fpio);
        NumFilesOpen--;
    }

    return (0);

}
