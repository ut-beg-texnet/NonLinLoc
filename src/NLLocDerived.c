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


/*   NLLocDerived.c

        Derived-quantity calculations for a NonLinLoc hypocenter: azimuth gap,
        arrival distances, arrival counts, Vp/Vs estimate, and magnitude
        (CalculateMagnitude plus the findStaInstComp / Calc_ML_HuttonBoore /
        Calc_MD_FMAG helpers).

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


/** function to determine primary and secondary azimuth gap as seen from epicenter (deg.
 *
 *     primary gap - largest azimuth separation between stations
 *     secondary gap - largest azimuth separation between stations filled by a single station
 *
 */

double* azimuths = NULL;

double CalcAzimuthGap(ArrivalDesc *arrival, int num_arrivals, double *pgap_secondary) {

    if (azimuths == NULL) {
        if ((azimuths = (double *) malloc(MAX_NUM_ARRIVALS * sizeof (double))) == NULL) {
            nll_puterr("ERROR: allocating memory for azimuths in CalcAzimuthGap().");
            return (-1);
        }
    }

    //printf("DEBUG: CalcAzimuthGap =========\n");

    /* load azimuths to working array */
    int naz = 0;
    for (int narr = 0; narr < num_arrivals; narr++) {
        // AJL 20091208 Zero weight phase modification
        //if (!(arrival + narr)->flag_ignore)
        if (!(arrival + narr)->flag_ignore && (arrival + narr)->weight > VERY_SMALL_DOUBLE) { // 20100521 AJL
            //if ((arrival + narr)->weight > 0.001)
            azimuths[naz++] = (arrival + narr)->azim;
            //printf("DEBUG: Azimuth added: %s %f\n", (arrival + narr)->label, (arrival + narr)->azim);
        }
    }

    /* sort */
    SortDoubles(azimuths, naz);

    double az_last, az_last2, az;
    double gap_primary, gap_primary_max = -1.0;
    double gap_secondary, gap_secondary_max = -1.0;

    // find largest gap and secondary gap
    az_last2 = azimuths[naz - 2] - 360.0;
    az_last = azimuths[naz - 1] - 360.0;
    // 20230712 AJL - bug fix, skip identical azimuths (assume from same station)
    int ndx = 3;
    while (az_last2 == az_last && ndx <= naz) {
        az_last2 = azimuths[naz - ndx] - 360.0;
        ndx++;
    }
    for (int narr = 0; narr < naz; narr++) {
        // 20230712 AJL - bug fix, skip identical azimuths (assume from same station)
        if (azimuths[narr] == az_last) {
            continue;
        }
        az = azimuths[narr];
        gap_primary = az - az_last;
        if (gap_primary > gap_primary_max)
            gap_primary_max = gap_primary;
        gap_secondary = az - az_last2;
        if (gap_secondary > gap_secondary_max)
            gap_secondary_max = gap_secondary;
        //printf("DEBUG: Azimuths tested: %f %f %f\n", az_last2, az_last, az);
        az_last2 = az_last;
        az_last = az;
    }

    *pgap_secondary = gap_secondary_max;
    //printf("DEBUG: Azimuth gaps: gap %f gap2 %f\n", gap_primary_max, gap_secondary_max);

    return (gap_primary_max);

    /*
    az_last = azimuth[0];
    for (narr = 1; narr < naz; narr++) {
        az = azimuth[narr];
        gap = az - az_last;
        if (gap > gap_max)
            gap_max = gap;
        az_last = az;
    }
    az = azimuth[0] + 360.0;
    gap = az - az_last;
    if (gap > gap_max)
        gap_max = gap;

    return (gap_max);
     * */

}

/** function to determine distance of closest non-ignored arrival
 *
 *  IMPORTANT! - Assumes arrivals sorted by distance, SortArrivalsDist()
 *
    // QML fields added for compatibility with QuakeML OriginQuality attributes (AJL 201005)
    double minimumDistance; // QML - Epicentral distance of station closest to the epicenter. Unit: km
    double maximumDistance; // QML - Epicentral distance of station farthest from the epicenter. Unit: km
    double medianDistance; // QML - Median epicentral distance of used stations. Unit: km
 */

double CalcArrivalDistances(ArrivalDesc *arrival, int num_arrivals, double *pmaximumDistance, double *pmedianDistance, int usedStationCount) {

    double minimumDistance = VERY_LARGE_DOUBLE;
    *pmaximumDistance = -1.0;
    *pmedianDistance = VERY_LARGE_DOUBLE;

    int stationCount = 0;
    char label_last[ARRIVAL_LABEL_LEN] = "!!!!!!";

    double dist;
    int narr;
    for (narr = 0; narr < num_arrivals; narr++) {
        // AJL 20091208 Zero weight phase modification
        //if (!(arrival + narr)->flag_ignore) {
        if (!(arrival + narr)->flag_ignore && (arrival + narr)->weight > VERY_SMALL_DOUBLE) { // 20100521 AJL
            dist = (arrival + narr)->dist;
            if (dist < minimumDistance)
                minimumDistance = dist;
            if (dist > *pmaximumDistance)
                *pmaximumDistance = dist;
            if (strcmp((arrival + narr)->label, label_last)) {
                stationCount++;
                if (usedStationCount % 2 == 1) { // case of usedStationCount odd
                    if (stationCount == 1 + usedStationCount / 2)
                        *pmedianDistance = dist;
                } else { // case of usedStationCount even
                    if (stationCount == usedStationCount / 2)
                        *pmedianDistance = dist;
                    else

                        if (stationCount == 1 + usedStationCount / 2)
                        *pmedianDistance = (*pmedianDistance + dist) / 2.0;
                }
            }
            snprintf(label_last, sizeof(label_last), "%s", (arrival + narr)->label);
        }
    }

    return (minimumDistance);

}

/** function to determine counts of arrivals used for location
 *
 *  IMPORTANT! - Assumes arrivals sorted by distance, SortArrivalsDist()
 *
    // QML fields added for compatibility with QuakeML OriginQuality attributes (AJL 201005)
    int associatedPhaseCount; // QML - Number of associated phases, regardless of their use for origin computation.
    // [->nreadings] int usedPhaseCount; // QML - Number of defining phases, i. e., phase observations that were actually used for computing
    // the origin. Note that there may be more than one defining phase per station.
    int associatedStationCount; // QML - Number of stations at which the event was observed.
    int usedStationCount; // QML - Number of stations from which data was used for origin computation.
    int depthPhaseCount; // QML - Number of depth phases (typically pP, sometimes sP) used in depth computation.
 */

int CalcArrivalCounts(ArrivalDesc *arrival, int num_arrivals, int num_arrivals_read,

        int* passociatedPhaseCount, // QML - Number of associated phases, regardless of their use for origin computation.
        // [->nreadings] int usedPhaseCount; // QML - Number of defining phases, i. e., phase observations that were actually used for computing
        // the origin. Note that there may be more than one defining phase per station.
        int* passociatedStationCount, // QML - Number of stations at which the event was observed.
        int* pusedStationCount, // QML - Number of stations from which data was used for origin computation.
        int* pdepthPhaseCount // QML - Number of depth phases (typically pP, sometimes sP) used in depth computation.
        ) {

    *passociatedPhaseCount = 0;
    *passociatedStationCount = 0;
    *pusedStationCount = 0;

    *pdepthPhaseCount = -1; // not supported by NLL
    *passociatedStationCount = -1; // not supported by NLL
    *passociatedPhaseCount = num_arrivals_read;

    int usedPhaseCount = 0;

    char label_last[ARRIVAL_LABEL_LEN] = "!!!!!!";
    char label_last_used[ARRIVAL_LABEL_LEN] = "!!!!!!";

    int ignored;
    int narr;
    for (narr = 0; narr < num_arrivals; narr++) {
        // AJL 20091208 Zero weight phase modification
        ignored = (arrival + narr)->flag_ignore;
        if (!ignored) {
            usedPhaseCount++;
            if (strcmp((arrival + narr)->label, label_last_used)) {
                (

                        *pusedStationCount)++;
            }
            snprintf(label_last_used, sizeof(label_last_used), "%s", (arrival + narr)->label);
        }
        //if (strcmp((arrival + narr)->label, label_last)) {
        //    (*passociatedStationCount)++;
        //}
        snprintf(label_last, sizeof(label_last), "%s", (arrival + narr)->label);
    }

    return (usedPhaseCount);

}



/** function to estimate Vp/Vs ratio */

/*  follows chapter 5 of Lahr, J.C., 1989. HYPOELLIPSE/Version 2.0: A computer program for
determining local earthquake hypocentral parameters, magnitude and first motion pattern, U.S.
Geological Survey Open-File Report 89-116, 92p.
 */

double CalculateVpVsEstimate(HypoDesc* phypo, ArrivalDesc* arrival, int narrivals) {

    int narr, npair;
    int k;
    double p_time[MAX_NUM_ARRIVALS], p_time_wt;
    double p_time_cent, p_error[MAX_NUM_ARRIVALS];
    double s_time[MAX_NUM_ARRIVALS], s_time_wt;
    double s_time_cent, s_error[MAX_NUM_ARRIVALS];
    double weight[MAX_NUM_ARRIVALS], weight_sum;
    double B, Btest, Bmin, dB, Tmin, T, temp;
    double tsp, tsp_min_max_diff;
    double tsp_min = VERY_LARGE_DOUBLE;
    double tsp_max = -VERY_LARGE_DOUBLE;

    //printf("CALCULATING VpVs\n");

    /* load P S pairs to working arrays (assumes arrivals sorted by distance) */

    npair = 0;
    for (narr = 1; narr < narrivals; narr++) {
        if (strcmp(arrival[narr - 1].label, arrival[narr].label) == 0
                && IsPhaseID(arrival[narr - 1].phase, "P")
                && IsPhaseID(arrival[narr].phase, "S")) {
            // DELAY_CORR
            // removing delay so add (Tobs = Tcorr + (O-C))
            p_time[npair] = arrival[narr - 1].obs_time + arrival[narr - 1].delay;
            p_error[npair] = arrival[narr - 1].error;
            // DELAY_CORR
            // removing delay so add (Tobs = Tcorr + (O-C))
            s_time[npair] = arrival[narr].obs_time + arrival[narr].delay;
            s_error[npair] = arrival[narr].error;
            //printf("PAIR %s %s (%lf +/- %lf) + %s %s (%lf +/- %lf)\n", arrival[narr - 1].label, arrival[narr - 1].phase, p_time[npair], p_error[npair], arrival[narr].label, arrival[narr].phase, s_time[npair], s_error[npair]);

            tsp = s_time[npair] - p_time[npair];
            tsp_min = tsp < tsp_min ? tsp : tsp_min;
            tsp_max = tsp > tsp_max ? tsp : tsp_max;

            npair++;
        }
    }

    tsp_min_max_diff = tsp_max - tsp_min;
    phypo->tsp_min_max_diff = tsp_min_max_diff;

    /* not enough pairs found */
    if (npair < 2) {
        phypo->VpVs = -1.0;
        phypo->nVpVs = npair;
        return (-1.0);
    }


    /* search for optimal VpVs */

    B = Bmin = 2.0;
    Tmin = VERY_LARGE_DOUBLE;
    for (dB = 0.5; dB > 0.00001; dB *= 0.4) {

        for (k = -3; k < 4; k++) {

            Btest = B + (double) k * dB;

            // form weights
            p_time_wt = 0.0;
            s_time_wt = 0.0;
            weight_sum = 0.0;
            // 20261026 Bug fix (Julián Olivar & Claude): misses first PS pair //for (narr = 1; narr < npair; narr++) {
            for (narr = 0; narr < npair; narr++) {
                weight[narr] = 1.0 /
                        (s_error[narr] * s_error[narr]
                        + Btest * p_error[narr] * p_error[narr]);
                p_time_wt += weight[narr] * p_time[narr];
                s_time_wt += weight[narr] * s_time[narr];
                weight_sum += weight[narr];
            }

            // form centered times and accumulate T
            T = 0.0;
            // 20261026 Bug fix (Julián Olivar & Claude): misses first PS pair //for (narr = 1; narr < npair; narr++) {
            for (narr = 0; narr < npair; narr++) {
                p_time_cent = p_time[narr] - p_time_wt / weight_sum;
                s_time_cent = s_time[narr] - s_time_wt / weight_sum;
                temp = s_time_cent - Btest * p_time_cent;
                T += weight[narr] * temp * temp;
            }

            if (T < Tmin) {
                Tmin = T;
                Bmin = Btest;
                //printf("  NEW MIN VpVs = %lf dB = %lf T = %le  k = %d\n", Btest, dB, T, k);
            }

        }

        B = Bmin;
    }


    //printf("  OPTIMAL VpVs = %.3lf last dB = %lf npair = %d\n", B, dB / 0.4, npair);

    phypo->VpVs = B;
    phypo->nVpVs = npair;

    return (B);

}

/** function to calculate magintude */

int CalculateMagnitude(HypoDesc* phypo, ArrivalDesc* parrivals,
        int narrivals, CompDesc* pcomp, int nCompDesc, MagDesc * pmagnitude) {

    int nmag, narr, nComp;
    double amp_mag_sum, dur_mag_sum;
    double amp_fact_ml_hb, sta_corr;
    ArrivalDesc* parr;



    /* calculate magnitude for specified calculation type */

    if (pmagnitude == MAG_UNDEF) {
        /* no magnitude calculation type given */

        return (0);

    } else if (pmagnitude->type == MAG_ML_HB) {
        /* ML - Hutton & Boore, BSSA, v77, n6, Dec 1987 */

        /* write message */
        if (message_flag >= 3) {
            snprintf(MsgStr, sizeof(MsgStr), "\nComponent results for: ML - Hutton & Boore, BSSA, v77, n6, Dec 1987:");
            nll_putmsg(3, MsgStr);
        }

        nmag = 0;
        amp_mag_sum = 0.0;
        for (narr = 0; narr < narrivals; narr++) {

            /* check for valid amplitude */
            if ((parr = parrivals + narr)->amplitude > 0.0 && parr->amplitude != AMPLITUDE_NULL) {

                nComp = findStaInstComp(parr, pcomp, nCompDesc);
                if (nComp >= 0) {
                    /* sta/inst/comp found */
                    amp_fact_ml_hb = (pcomp + nComp)->amp_fact_ml_hb;
                    sta_corr = (pcomp + nComp)->sta_corr_ml_hb;
                } else {
                    amp_fact_ml_hb = 1.0;
                    sta_corr = 0.0;
                }

                /* calc ML */
                parr->amp_mag = Calc_ML_HuttonBoore(
                        parr->amplitude * pmagnitude->amp_fact_ml_hb * amp_fact_ml_hb,
                        parr->dist, phypo->depth, sta_corr,
                        pmagnitude->hb_n, pmagnitude->hb_K,
                        pmagnitude->hb_Ro, pmagnitude->hb_Mo);

                /* write message */
                if (message_flag >= 3) {
                    snprintf(MsgStr, sizeof(MsgStr), "%s %s %s amp %.2e f %.2e f_sta %.2e dist %.2f depth %.2f sta_corr %.4f hb_n %.2f hb_K %.5f mag %.2f",
                            parr->label, parr->inst, parr->comp, parr->amplitude, pmagnitude->amp_fact_ml_hb, amp_fact_ml_hb, parr->dist,
                            phypo->depth, sta_corr, pmagnitude->hb_n, pmagnitude->hb_K, parr->amp_mag);
                    nll_putmsg(3, MsgStr);
                }


                /* update event magnitude */
                amp_mag_sum += parr->amp_mag;
                nmag++;
            }
        }
        if (nmag > 0) {
            phypo->amp_mag = amp_mag_sum / (double) nmag;
            phypo->num_amp_mag = nmag;
        }
    } else if (pmagnitude->type == MAG_MD_FMAG) {
        /* coda duration (FMAG) - HYPOELLIPSE users manual chap 4;
        Lee et al., 1972; Lahr et al., 1975; Bakun and Lindh, 1977 */

        nmag = 0;
        dur_mag_sum = 0.0;
        for (narr = 0; narr < narrivals; narr++) {

            /* check for valid amplitude */
            if ((parr = parrivals + narr)->coda_dur > 0.0 && parr->coda_dur != CODA_DUR_NULL) {

                nComp = findStaInstComp(parr, pcomp, nCompDesc);
                if (nComp >= 0) {
                    /* sta/inst/comp found */
                    sta_corr = (pcomp + nComp)->sta_corr_md_fmag;
                } else {
                    sta_corr = 1.0;
                }

                /* calc ML */
                parr->dur_mag = Calc_MD_FMAG(
                        parr->coda_dur, parr->dist, phypo->depth, sta_corr,
                        pmagnitude->fmag_c1, pmagnitude->fmag_c2, pmagnitude->fmag_c3,
                        pmagnitude->fmag_c4, pmagnitude->fmag_c5);

                /*printf("%s %s %s coda_dur %lf dist %lf depth %lf sta_corr %lf c1 %lf c2 %lf c3 %lf c4 %lf c5 %lf mag %lf\n",
                                                parr->label, parr->inst, parr->comp, parr->coda_dur, parr->dist,
                                                phypo->depth, sta_corr, pmagnitude->fmag_c1, pmagnitude->fmag_c2, pmagnitude->fmag_c3,
                                                pmagnitude->fmag_c4, pmagnitude->fmag_c5, parr->dur_mag);
                 */

                /* update event magnitude */
                dur_mag_sum += parr->dur_mag;
                nmag++;
            }
        }
        if (nmag > 0) {

            phypo->dur_mag = dur_mag_sum / (double) nmag;
            phypo->num_dur_mag = nmag;
        }
    }

    return (0);
}

/** function to find component parameters */

int findStaInstComp(ArrivalDesc* parr, CompDesc* pcomp, int nCompDesc) {

    int nComp;
    char *pchr, test_label[ARRIVAL_LABEL_LEN];

    snprintf(test_label, sizeof(test_label), "%s", parr->time_grid_label);

    for (nComp = 0; nComp < nCompDesc; nComp++) {
        snprintf(test_label, sizeof(test_label), "%s", parr->time_grid_label);
        if ((pchr = strrchr(test_label, '_')) != NULL)
            *pchr = '\0';
        //printf("comp %s  arr_test %s %s %s\n", (pcomp + nComp)->label, parr->label, parr->time_grid_label, test_label);
        if ((pcomp + nComp)->label[0] != ARRIVAL_NULL_CHR &&
                (pcomp + nComp)->label[0] != '*' &&
                strcmp((pcomp + nComp)->label, test_label) != 0)
            continue;
        if ((pcomp + nComp)->inst[0] != ARRIVAL_NULL_CHR &&
                (pcomp + nComp)->inst[0] != '*' &&
                strcmp((pcomp + nComp)->inst, parr->inst) != 0)
            continue;
        if ((pcomp + nComp)->comp[0] != ARRIVAL_NULL_CHR &&
                (pcomp + nComp)->comp[0] != '*' &&
                strcmp((pcomp + nComp)->comp, parr->comp) != 0)

            continue;
        // found
        return (nComp);
    }

    // not found
    return (-1);
}


/** function to calculate local magnitude ML */

/* ML - Hutton & Boore, BSSA, v77, n6, Dec 1987 */

double Calc_ML_HuttonBoore(double amplitude, double dist, double depth, double sta_corr,
        double hb_n, double hb_K, double hb_Ro, double hb_Mo) {

    double hyp_dist, magnitude;

    hyp_dist = sqrt(dist * dist + depth * depth);

    if (hyp_dist < SMALL_DOUBLE)
        return (MAGNITUDE_NULL);

    magnitude = log10(amplitude)
            + hb_n * log10(hyp_dist / hb_Ro) + hb_K * (hyp_dist - hb_Ro)
            + hb_Mo + sta_corr;

    return (magnitude);

}



/** function to calculate coda duration magnitude MD */

/* coda duration (FMAG) - HYPOELLIPSE users manual chap 4;
        Lee et al., 1972; Lahr et al., 1975; Bakun and Lindh, 1977 */

double Calc_MD_FMAG(double coda_dur, double dist, double depth, double sta_corr,
        double fmag_c1, double fmag_c2, double fmag_c3, double fmag_c4, double fmag_c5) {

    double magnitude;

    if (coda_dur * sta_corr < SMALL_DOUBLE)
        return (MAGNITUDE_NULL);

    magnitude = fmag_c1
            + fmag_c2 * log10(coda_dur * sta_corr)
            + fmag_c3 * dist
            + fmag_c4 * depth
            + fmag_c5 * pow(log10(coda_dur * sta_corr), 2);

    return (magnitude);

}
