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


/*   NLLocQuality.c

        Solution-quality / likelihood evaluation for NonLinLoc: the
        CalcSolutionQuality* family (EDT, OT_STACK, ML_OT, L1_NORM, GAU_TEST,
        GAU_ANALYTIC), the maximum-likelihood / likelihood / variance origin-
        time helpers, weight normalization, CalcMaxLikeOriginTime,
        UpdateProbabilisticResiduals, CalcConfidenceIntrvl, and the file-local
        getLogPdfValue (guarded by STACK_POSTERIOR / TEST_*_POSTERIOR).

        Extracted verbatim from NLLocLib.c during the modernization effort
        (Phase 3: modularization). Pure code movement, no logic change; the
        CalcSolutionQuality* prototypes remain in NLLocLib.h, so callers are
        unaffected. getLogPdfValue and the STACK_POSTERIOR define are local to
        this block and move with it, so the same variant compiles as before.
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

/* moved with the getLogPdfValue/CalcSolutionQuality block from NLLocLib.c */
#define EDT_OT_WT_FLOOR log(0.00001)


#define STACK_POSTERIOR
// various test to satisfy Reviewer #1
//#define TEST_TARGET_PRODUCT_POSTERIOR
//#define TEST_FULL_PRODUCT_POSTERIOR

#ifdef STACK_POSTERIOR

// stack over all other events and target event with coherence weights

double getLogPdfValue(SearchPdfGridDesc *searchPdfGrid, double hypo_x, double hypo_y, double hypo_z) {

    double log_pdf_value = 0.0;

    if (searchPdfGrid->gridType == PDF_GRID_GRID) {
        double pdf_value = (double) ReadAbsInterpGrid3d(NULL, &searchPdfGrid->grid, hypo_x, hypo_y, hypo_z, 1);
        if (pdf_value < searchPdfGrid->default_value) {
            pdf_value = searchPdfGrid->default_value;
        }
        if (pdf_value > FLT_MIN) {
            log_pdf_value = log(pdf_value);
        }
    } else if (searchPdfGrid->gridType == PDF_GRID_OCT_TREE) {
        double pdf_value = 0.0;
        double weight_sum = 0.0;
        OctNode* node;
        Vect3D coords;
        coords.x = hypo_x;
        coords.y = hypo_y;
        coords.z = hypo_z;
        for (int ngrid = 0; ngrid < searchPdfGrid->nGrids; ngrid++) {
            //printf("DEBUG: ngrid : %d\n", ngrid);
            if (searchPdfGrid->coherence[ngrid] > searchPdfGrid->coherence_min) {
                node = getLeafNodeContaining(searchPdfGrid->tree3D[ngrid], coords);
                if (node != NULL) { // 20200528 AJL - bug fix.
                    double value = (double) node->value;
                    //                    if (value > maxvalue) {
                    //                        printf("DEBUG: node %ld  value %le  weight %lf  coords.x y z: %f %f %f\n", (long) node, value, searchPdfGrid->weight[ngrid], coords.x, coords.y, coords.z);
                    //                        maxvalue = value;
                    //                    }
                    if (value < searchPdfGrid->default_value) {
                        value = searchPdfGrid->default_value;
                    }
                    pdf_value += value * searchPdfGrid->weight[ngrid];
                    weight_sum += searchPdfGrid->weight[ngrid];
                    //pdf_value += value * searchPdfGrid->coherence[ngrid];
                    // pdf_value += value * searchPdfGrid->coherence[ngrid] * searchPdfGrid->coherence[ngrid];   // square of coherence
                }
            }
        }
        if (pdf_value > FLT_MIN) {

            log_pdf_value = log(pdf_value);
            // 20200617 AJL - raise pdf value to power of weight_sum (like raising EDT to power of N -> concentrates pdf)
            log_pdf_value *= weight_sum;
        }

    }

    return (log_pdf_value);

}

#else
#ifdef TEST_TARGET_PRODUCT_POSTERIOR

// weighted stack of other event pdf's multiplied by target event pdf
// gives much worse results relative to STACK_POSTERIOR
//    /Users/anthony/work_temp/nlloc_tmp/Hukkakero_2007/20210809B_FP_3D_PRODUCT_TEST/ak135/pdf_prior_posterior/Hukkakero/RUN1000/001_2.0-25Hz_dmsl4.0_cm0.45__coherence_0.45_TARGET_PRODUCT_ONLY/Hukkakero_2007.sum_ALL.grid0.loc.hyp

// if use "// EXP raise pdf value to power of weight_sum", gives results only slightly more scattered than STACK_POSTERIOR
//    /Users/anthony/work_temp/nlloc_tmp/Hukkakero_2007/20210809B_FP_3D_PRODUCT_TEST/ak135/pdf_prior_posterior/Hukkakero/RUN1000/001_2.0-25Hz_dmsl4.0_cm0.45__coherence_0.45_TARGET_PRODUCT_ONLY_EXP/Hukkakero_2007.sum_ALL.grid0.loc.hyp

double getLogPdfValue(SearchPdfGridDesc *searchPdfGrid, double hypo_x, double hypo_y, double hypo_z) {

    double log_pdf_value = 0.0;

    if (searchPdfGrid->gridType == PDF_GRID_GRID) {
        double pdf_value = (double) ReadAbsInterpGrid3d(NULL, &searchPdfGrid->grid, hypo_x, hypo_y, hypo_z, 1);
        if (pdf_value < searchPdfGrid->default_value) {
            pdf_value = searchPdfGrid->default_value;
        }
        if (pdf_value > FLT_MIN) {
            log_pdf_value = log(pdf_value);
        }
    } else if (searchPdfGrid->gridType == PDF_GRID_OCT_TREE) {
        double pdf_value = 0.0;
        double weight_sum = 0.0;
        OctNode* node;
        Vect3D coords;
        coords.x = hypo_x;
        coords.y = hypo_y;
        coords.z = hypo_z;
        double target_pdf_value = 0.0;
        for (int ngrid = 0; ngrid < searchPdfGrid->nGrids; ngrid++) {
            //printf("DEBUG: ngrid : %d\n", ngrid);
            if (searchPdfGrid->coherence[ngrid] > searchPdfGrid->coherence_min) {
                node = getLeafNodeContaining(searchPdfGrid->tree3D[ngrid], coords);
                if (node != NULL) { // 20200528 AJL - bug fix.
                    double value = (double) node->value;
                    // stack pdf over all but target event
                    if (ngrid == 0) { // target event
                        target_pdf_value = value;
                        continue;
                    }
                    //                    if (value > maxvalue) {
                    //                        printf("DEBUG: node %ld  value %le  weight %lf  coords.x y z: %f %f %f\n", (long) node, value, searchPdfGrid->weight[ngrid], coords.x, coords.y, coords.z);
                    //                        maxvalue = value;
                    //                    }
                    if (value < searchPdfGrid->default_value) {
                        value = searchPdfGrid->default_value;
                    }
                    pdf_value += value * searchPdfGrid->weight[ngrid];
                    weight_sum += searchPdfGrid->weight[ngrid];
                }
            }
        }

        if (pdf_value > FLT_MIN) {

            log_pdf_value = log(pdf_value);
            // EXP raise pdf value to power of weight_sum
            //log_pdf_value *= weight_sum;

            // multiply by target event pdf
            if (target_pdf_value > FLT_MIN) {
                log_pdf_value += log(target_pdf_value);
            }

        }

    }

    return (log_pdf_value);

}

#else
#ifdef TEST_FULL_PRODUCT_POSTERIOR

// product of all events with no weight
// gives worse results than to FULL_PRODUCT_POSTERIOR
//    /Users/anthony/work_temp/nlloc_tmp/Hukkakero_2007/20210809B_FP_3D_PRODUCT_TEST/ak135/pdf_prior_posterior/Hukkakero/RUN1000/001_2.0-25Hz_dmsl4.0_cm0.45__coherence_0.45_FULL_PRODUCT_NO_WT/Hukkakero_2007.sum_ALL.grid0.loc.hyp

// product of all events with exp weight
// gives moderately worse results than to FULL_PRODUCT_POSTERIOR, better than FULL_PRODUCT_NO_WT
//    /Users/anthony/work_temp/nlloc_tmp/Hukkakero_2007/20210809B_FP_3D_PRODUCT_TEST/ak135/pdf_prior_posterior/Hukkakero/RUN1000/001_2.0-25Hz_dmsl4.0_cm0.45__coherence_0.45_FULL_PRODUCT_EXP_WT/Hukkakero_2007.sum_ALL.grid0.loc.hyp

double getLogPdfValue(SearchPdfGridDesc *searchPdfGrid, double hypo_x, double hypo_y, double hypo_z) {

    double log_pdf_value = 0.0;

    if (searchPdfGrid->gridType == PDF_GRID_GRID) {
        double pdf_value = (double) ReadAbsInterpGrid3d(NULL, &searchPdfGrid->grid, hypo_x, hypo_y, hypo_z, 1);
        if (pdf_value < searchPdfGrid->default_value) {
            pdf_value = searchPdfGrid->default_value;
        }
        if (pdf_value > FLT_MIN) {
            log_pdf_value = log(pdf_value);
        }
    } else if (searchPdfGrid->gridType == PDF_GRID_OCT_TREE) {
        double pdf_value = 0.0;
        double weight_sum = 0.0;
        OctNode* node;
        Vect3D coords;
        coords.x = hypo_x;
        coords.y = hypo_y;
        coords.z = hypo_z;
        double target_pdf_value = 0.0;
        for (int ngrid = 0; ngrid < searchPdfGrid->nGrids; ngrid++) {
            //printf("DEBUG: ngrid : %d\n", ngrid);
            if (searchPdfGrid->coherence[ngrid] > searchPdfGrid->coherence_min) {
                node = getLeafNodeContaining(searchPdfGrid->tree3D[ngrid], coords);
                if (node != NULL) { // 20200528 AJL - bug fix.
                    double value = (double) node->value;
                    //                    if (value > maxvalue) {
                    //                        printf("DEBUG: node %ld  value %le  weight %lf  coords.x y z: %f %f %f\n", (long) node, value, searchPdfGrid->weight[ngrid], coords.x, coords.y, coords.z);
                    //                        maxvalue = value;
                    //                    }
                    if (value < searchPdfGrid->default_value) {
                        value = searchPdfGrid->default_value;
                    }
                    if (value > FLT_MIN) {

                        pdf_value += log(value) * searchPdfGrid->weight[ngrid];
                        // NO_WT  pdf_value += log(value);
                    }
                    //weight_sum += searchPdfGrid->weight[ngrid];
                }
            }
        }

        log_pdf_value = pdf_value;

    }

    return (log_pdf_value);

}

#endif
#endif
#endif

/** function to calculate probability density */

double CalcSolutionQuality(double hypo_x, double hypo_y, double hypo_z, OctNode* poct_node, int num_arrivals, ArrivalDesc *arrival,
        GaussLocParams* gauss_par, int itype, double* pmisfit, double* potime, double* potime_var,
        double cell_half_diagonal_time_range, double cell_diagonal, double cell_volume,
        double* peffective_cell_size, double *pot_variance_factor, double *log_prior) {

    *log_prior = 0.0;
    if (iUseSearchPrior) {
        SearchPdfGridDesc *searchPdfGrid = &SearchPrior;
        *log_prior = getLogPdfValue(searchPdfGrid, hypo_x, hypo_y, hypo_z);
    }
    if (potime == NULL) { // do not need hypo stats (e.g. for SaveBestLocation()), just return posterior
        if (iUseSearchPosterior) {
            SearchPdfGridDesc *searchPdfGrid = &SearchPosterior;
            double log_posterior = getLogPdfValue(searchPdfGrid, hypo_x, hypo_y, hypo_z);
            return (log_posterior);
        }
    }
    //printf("DEBUG: hypo_x %f, hypo_y %f, hypo_z %f\n, ", hypo_x, hypo_y, hypo_z);
    double value;
    if (LocMethod == METH_GAU_ANALYTIC) {
        value = CalcSolutionQuality_GAU_ANALYTIC(num_arrivals, arrival, gauss_par, itype, pmisfit, potime);
    } else if (LocMethod == METH_GAU_TEST) {
        value = CalcSolutionQuality_GAU_TEST(num_arrivals, arrival, gauss_par, itype, pmisfit, potime);
    } else if (LocMethod == METH_L1_NORM) {
        value = CalcSolutionQuality_L1_NORM(num_arrivals, arrival, gauss_par, itype, pmisfit, potime);
    } else if (LocMethod == METH_OT_STACK) {
        value = CalcSolutionQuality_OT_STACK(poct_node, num_arrivals, arrival,
                gauss_par, itype, pmisfit, potime, potime_var, cell_half_diagonal_time_range, cell_diagonal, cell_volume, peffective_cell_size, pot_variance_factor);
        return (value);
    } else if (LocMethod == METH_ML_OT) {
        value = CalcSolutionQuality_ML_OT(num_arrivals, arrival,
                gauss_par, itype, pmisfit, potime, potime_var, cell_half_diagonal_time_range, 0);
    } else if (LocMethod == METH_EDT) {
        value = CalcSolutionQuality_EDT(num_arrivals, arrival,
                gauss_par, itype, pmisfit, potime, potime_var, cell_half_diagonal_time_range, 0);
    } else if (LocMethod == METH_EDT_BOX) {
        value = CalcSolutionQuality_EDT(num_arrivals, arrival,
                gauss_par, itype, pmisfit, potime, potime_var, cell_half_diagonal_time_range, 1);
    } else {
        return (-1.0);
    }

    if (potime != NULL) { // need hypo stats (e.g. for SaveBestLocation()), also return posterior
        if (iUseSearchPosterior) {
            SearchPdfGridDesc *searchPdfGrid = &SearchPosterior;
            double log_posterior = getLogPdfValue(searchPdfGrid, hypo_x, hypo_y, hypo_z);

            return (log_posterior);
        }
    }

    return (value);

}





/** function to calculate probability density */

/*	EDT - sum of probabilities of difference of obs - difference of travel times
                for all pairs of obs
 */

double CalcSolutionQuality_EDT(int num_arrivals, ArrivalDesc *arrival,
        GaussLocParams* gauss_par, int itype, double* pmisfit, double* potime,
        double* potime_var, double cell_half_diagonal_time_range, int method_box) {

    double cell_diagonal_time_var = cell_half_diagonal_time_range * cell_half_diagonal_time_range;

    int nrow, ncol;

    long double edt_sum, prob, edtSumMisfit;

    double edt_misfit, edt_weight;
    double weight, weight2;
    double ln_prob_density, rms_misfit;

    MatrixDouble edtmtx;
    double sigma2_row;
    double obs_minus_pred;

    int no_abs_time_row;

    // EDT_OT_WT
    int num_otime_error;
    long double ot_row, ot_prob, ot_row_2, ot_error_2;
    long double ot_sum, ot_weight, ot_var, ot_2_sum, ot_var_weight;

    // EDT_OT_WT_ML
    double ot_ml = 0.0, ot_ml_var;

    // OT additions 20071220
    double ot_prob_max = 0.0;

    // Gauss2
    double tt_error;


    // method_box
    //double error_row;
    double amp_row, unc_limit;

    // search pdf different from true
    int iuse_cell_diagonal_time_var;
    //double sigma2_row_search = 0.0;
    //double weight_search = 0.0, weight2_search, edt_weight_search;
    //long double prob_search = 0.0L, edt_sum_search = 0.0L, edtSumMisfit_search = 0.0L;


    // initialize otime flags
    int icalc_otime = 0;
    int icalc_otime_default = 0;
    int icalc_otime_force_ml = 0;
    if (potime != NULL) {
        icalc_otime = 1;
        icalc_otime_default = 1; // final OT is maximum likelihood OT
        if (0) { // Non-standard, use only for testing
            icalc_otime_default = 0;
            icalc_otime_force_ml = 1; // final OT is ot_sum / ot_weight for EDT_OT_WT or EDT
        }
    }

    edtmtx = gauss_par->EDTMtrx;

    // check if use_cell_diagonal_time_var
    iuse_cell_diagonal_time_var = 0;
    if (cell_diagonal_time_var > 0.0)
        iuse_cell_diagonal_time_var = 1;

    // check size of EDT_OT_WT_ML static arrays
    if ((EDT_use_otime_weight == 2 || icalc_otime_default)) {
        if (isize_ot_ml_array < num_arrivals) {
            isize_ot_ml_array = num_arrivals;
            free(ot_ml_arrival);
            if ((ot_ml_arrival = (double *) calloc(isize_ot_ml_array, sizeof (double))) == NULL)
                nll_puterr("ERROR: allocating double storage array for EDT_OT_WT_ML ot_ml_arrival.");
            free(ot_ml_arrival_edt_sum);
            if ((ot_ml_arrival_edt_sum = (double *) calloc(isize_ot_ml_array, sizeof (double))) == NULL)
                nll_puterr("ERROR: allocating double storage array for EDT_OT_WT_ML ot_ml_arrival_edt_sum.");
        }
        for (nrow = 0; nrow < num_arrivals; nrow++)
            ot_ml_arrival_edt_sum[nrow] = 0.0;
    }

    if (icalc_otime) {
        for (nrow = 0; nrow < num_arrivals; nrow++)
            arrival[nrow].weight = 0.0;
    }


    /* calculate weighted mean of predicted travel times  */
    /*		(TV82, eq. A-38) */
    CalcCenteredTimesPred(num_arrivals, arrival, gauss_par); // not used for EDT


    /* calculate EDT prop sum */

    edt_sum = 0.0L;
    edt_weight = 0.0;
    ot_prob = 0.0L;
    ot_row = 0.0L;
    ot_row_2 = 0.0L;
    ot_sum = 0.0L;
    ot_2_sum = 0.0L;
    ot_var = 0.0L;
    ot_var_weight = 0.0L;
    ot_weight = 0.0L;
    ot_error_2 = 0.0L;
    num_otime_error = 0;
    //#define TEST_COUNT_ONLY_USED_ARRIVALS
#ifdef TEST_COUNT_ONLY_USED_ARRIVALS
    int num_arrivals_used = 0;
#endif
    for (nrow = 0; nrow < num_arrivals; nrow++) {

        //printf("DEBUG: arrival[%d].pred_travel_time %f\n", nrow, arrival[nrow].pred_travel_time);

        // AJL 20041115 bug fix!
        if (arrival[nrow].pred_travel_time <= 0.0) {
            // iniitalize EDT_OT_WT_ML values
            if (EDT_use_otime_weight == 2 || icalc_otime_default) {
                ot_ml_arrival_edt_sum[nrow] = -1.0;
            }
            continue; // ignore obs without predicted times
        }
        // END

#ifdef TEST_COUNT_ONLY_USED_ARRIVALS
        num_arrivals_used++;
#endif
        // set error
        //printf("iUseGauss2 %d\n", iUseGauss2);
        if (iUseGauss2) {
            tt_error = arrival[nrow].pred_travel_time * Gauss2.SigmaTfraction;
            if (tt_error < Gauss2.SigmaTmin)
                tt_error = Gauss2.SigmaTmin;
            if (tt_error > Gauss2.SigmaTmax)
                tt_error = Gauss2.SigmaTmax;
            if (icalc_otime) {
                arrival[nrow].tt_error = tt_error;
                //printf("DEBUG: arrival[nrow].pred_travel_time %f\t  tt_error %f, arrival[nrow].error %f\n", arrival[nrow].pred_travel_time, tt_error, arrival[nrow].error);
            }
            tt_error *= tt_error;
            edtmtx[nrow][nrow] = arrival[nrow].error * arrival[nrow].error + tt_error;
            sigma2_row = edtmtx[nrow][nrow];
        }
        sigma2_row = edtmtx[nrow][nrow];

        if (iuse_cell_diagonal_time_var)
            sigma2_row += cell_diagonal_time_var;
        /*if (iuse_cell_diagonal_time_var) {
        sigma2_row_search = edtmtx[nrow][nrow] + cell_diagonal_time_var;
}*/
        //error_row = arrival[nrow].error;
        amp_row = arrival[nrow].amplitude;
        obs_minus_pred = arrival[nrow].obs_centered - arrival[nrow].pred_centered;
        no_abs_time_row = !arrival[nrow].abs_time;
        // SP-TEST PATCH: picks without absolute time give no origin-time estimate
        if (no_abs_time_row) {
            ;
        } else if (EDT_use_otime_weight == 2 || icalc_otime_default) { // EDT_OT_WT_ML or otime
            ot_ml_arrival[nrow] = arrival[nrow].obs_time - (long double) arrival[nrow].pred_travel_time;
            //ot_ml_arrival_edt_sum[nrow] = 0.0;
            ot_error_2 += sigma2_row;
            num_otime_error++;
        } else if (EDT_use_otime_weight == 1 || icalc_otime_force_ml) { // EDT_OT_WT or EDT
            ot_prob = 0.0;
            ot_row = arrival[nrow].obs_time - (long double) arrival[nrow].pred_travel_time;
            ot_row_2 = ot_row * ot_row;
            ot_error_2 += sigma2_row;
            num_otime_error++;
        }
        for (ncol = nrow + 1; ncol < num_arrivals; ncol++) {
            // AJL 20041115 bug fix!
            if (arrival[ncol].pred_travel_time <= 0.0)
                continue; // ignore obs without predicted times
            // END
            // check absolute timing
            // SP-TEST PATCH: symmetric - skip mixed pairs whichever pick comes first
            if (no_abs_time_row || !arrival[ncol].abs_time) {
                if (!no_abs_time_row || arrival[ncol].abs_time) // mixed pair: one has absolute time, the other not
                    continue;
                if (strcmp(arrival[nrow].label, arrival[ncol].label) != 0
                        || strcmp(arrival[nrow].inst, arrival[ncol].inst) != 0)
                    continue; // not same sta/inst
            }
            // calculate EDT misfit:  (obs1 - obs2) - (pred1 - pred2)
            edt_misfit = (double) (obs_minus_pred + arrival[ncol].pred_centered - arrival[ncol].obs_centered);
            // set error
            if (iUseGauss2) {
                tt_error = arrival[ncol].pred_travel_time * Gauss2.SigmaTfraction;
                if (tt_error < Gauss2.SigmaTmin)
                    tt_error = Gauss2.SigmaTmin;
                if (tt_error > Gauss2.SigmaTmax)
                    tt_error = Gauss2.SigmaTmax;
                tt_error *= tt_error;
                edtmtx[ncol][ncol] = arrival[ncol].error * arrival[ncol].error + tt_error;
            }
            // calculate probability
            if (method_box) {
                unc_limit = amp_row + arrival[ncol].amplitude; // sum of mean pick unc for each box
                //unc_limit = error_row + arrival[ncol].error;	// sum of box widths
                prob = fabs(edt_misfit) <= unc_limit ? 1.0 : 0.0;
                weight = amp_row * arrival[ncol].amplitude; // product of mean pick unc for each box
                weight *= (1.0 - edtmtx[nrow][ncol]); // correlation coeff
            } else {
                if (iuse_cell_diagonal_time_var)
                    weight2 = 1.0 / (sigma2_row + edtmtx[ncol][ncol] + cell_diagonal_time_var); // sum of errors**2
                else
                    weight2 = 1.0 / (sigma2_row + edtmtx[ncol][ncol]); // sum of errors**2
                prob = exp(-0.5 * edt_misfit * edt_misfit * weight2);
                weight = sqrt(weight2); // errors factor
                weight *= (1.0 - edtmtx[nrow][ncol]); // correlation coeff
                /*if (iuse_cell_diagonal_time_var) {	// duplicate above 4 lines
                weight2_search = 1.0 / (sigma2_row_search + edtmtx[ncol][ncol] + cell_diagonal_time_var);	// sum of errors**2
                prob_search = exp(-0.5 * edt_misfit * edt_misfit * weight2_search);
                weight_search = sqrt(weight2_search);		// errors factor
                weight_search *= (1.0 - edtmtx[nrow][ncol]);		// correlation coeff
        }*/
            }
            // 20130627 AJL - change weighing from sum to product
            //if (iSetStationDistributionWeights)
            //    weight *= (arrival[nrow].station_weight + arrival[ncol].station_weight) / 2.0;
            if (iSetStationDistributionWeights)
                weight *= sqrt(arrival[nrow].station_weight * arrival[ncol].station_weight);
            // 20130627 AJL - add prior weighting as product
            if (iUseArrivalPriorWeights && arrival[nrow].apriori_weight >= -VERY_SMALL_DOUBLE && arrival[ncol].apriori_weight >= -VERY_SMALL_DOUBLE)
                weight *= sqrt(arrival[nrow].apriori_weight * arrival[ncol].apriori_weight);
            prob *= weight;
            edt_sum += prob;
            edt_weight += weight;
            /*if (iuse_cell_diagonal_time_var) {	// duplicate above 5 lines
            if (iSetStationDistributionWeights)
            weight_search *= (arrival[nrow].station_weight + arrival[ncol].station_weight) / 2.0;
            prob_search *= weight_search;
            edt_sum_search += prob_search;
            edt_weight_search += weight_search;
    }*/
            // accumulate EDT weights
            if (icalc_otime) {
                //arrival[ncol].weight += weight;
                //arrival[nrow].weight += weight;
                arrival[ncol].weight += prob;
                arrival[nrow].weight += prob;
            }
            // otime
            if (EDT_use_otime_weight == 2 || icalc_otime_default) { // EDT_OT_WT_ML or otime
                /*if (iuse_cell_diagonal_time_var)	// ???? TEST
                ot_ml_arrival_edt_sum[nrow] += prob_search;
                else*/
                ot_ml_arrival_edt_sum[ncol] += prob;
                // AJL 20070326 bug fix!
                ot_ml_arrival_edt_sum[nrow] += prob;
            } else if (EDT_use_otime_weight == 1 || icalc_otime_force_ml) { // EDT_OT_WT or EDT
                ot_prob += prob;
            }
        }
        if (EDT_use_otime_weight == 1 && !no_abs_time_row) { // EDT_OT_WT or EDT  (SP-TEST PATCH: abs rows only)
            ot_sum += ot_prob * ot_row;
            ot_2_sum += ot_prob * ot_row_2;
            ot_weight += ot_prob;
        }
    }

    // SP-TEST PATCH: picks without absolute time are excluded from the origin-time likelihood
    if (EDT_use_otime_weight == 2 || icalc_otime_default) {
        for (nrow = 0; nrow < num_arrivals; nrow++)
            if (!arrival[nrow].abs_time)
                ot_ml_arrival_edt_sum[nrow] = -1.0;
    }

    // OT_WT methods
    if (EDT_use_otime_weight == 2 || icalc_otime_default) { // EDT_OT_WT_ML
        // EDT_OT_WT_ML method
        ot_ml = calc_maximum_likelihood_ot(ot_ml_arrival, ot_ml_arrival_edt_sum, num_arrivals, arrival, edtmtx, &ot_ml_var, icalc_otime, &ot_prob_max);
        if (icalc_otime && potime_var != NULL) {
            snprintf(MsgStr, sizeof(MsgStr), "INFO: EDT_otime_weight: ot_ml_std %lf\n", sqrt(ot_ml_var));
            nll_putmsg(2, MsgStr);
            *potime_var = ot_ml_var;
        }
        ot_var_weight = -ot_ml_var / (ot_error_2 / (long double) num_otime_error);
        if (ot_var_weight > EDT_OT_WT_FLOOR) {
            if (!EDT_otime_weight_active) {
                EDT_otime_weight_active = 1;
                snprintf(MsgStr, sizeof(MsgStr), "INFO: EDT_otime_weight activated, OT_WT exceeds EDT_OT_WT_FLOOR.");
                nll_putmsg(2, MsgStr);
            }
        } else {
            ot_var_weight = EDT_OT_WT_FLOOR;
        }
    } else if ((EDT_use_otime_weight == 1 || icalc_otime_force_ml) && num_otime_error > 0) { // EDT_OT_WT
        // EDT_OT_WT method
        ot_var = ot_2_sum / ot_weight - (ot_sum / ot_weight) * (ot_sum / ot_weight);
        if (icalc_otime && potime_var != NULL) {
            printf("ot_ml_std %lf\n", sqrt(ot_var));
            *potime_var = ot_var;
        }
        ot_var_weight = -ot_var / (ot_error_2 / (long double) num_otime_error);
        if (ot_var_weight > EDT_OT_WT_FLOOR) {
            if (!EDT_otime_weight_active) {
                EDT_otime_weight_active = 1;
                snprintf(MsgStr, sizeof(MsgStr), "INFO: EDT_otime_weight activated, OT_WT exceeds EDT_OT_WT_FLOOR.");
                nll_putmsg(2, MsgStr);
            }
        } else {
            ot_var_weight = EDT_OT_WT_FLOOR;
        }
        //edt_sum *= ot_var_weight;
        //edt_weight *= ot_var_weight;
    }

    // avoid numerical problems
    if (edt_sum < SMALL_DOUBLE)
        edt_sum = SMALL_DOUBLE;
    /*if (edt_sum_search < SMALL_DOUBLE)
    edt_sum_search = SMALL_DOUBLE;*/
    if (num_arrivals == 1) {
        edt_weight = 1.0; // 20190826 AJL - bug fix when locating with only one arrival (e.g. LOCPOSTERIOR)
    }
    // create EDT misfit
    if (edt_weight > 0.0 && num_arrivals > 0) {
        // AJL 20051115 bug fix!
        //edtSumMisfit = -log(edt_sum / edt_weight);
        edtSumMisfit = -log(edt_sum); // non-normalized sum of edt probs - favors solutions with more readings
        // END
        //edtSumMisfit *= 2.0;
        //edtSumMisfit *= num_arrivals;	// best, give power of N
        //edtSumMisfit *= sqrt(num_arrivals);	// works much better for INGV Italy region scale
        /*if (iuse_cell_diagonal_time_var) { 	// duplicate above lines
        edtSumMisfit_search = -log(edt_sum_search);	// non-normalized sum of edt probs - favors solutions with more readings
}*/
        // otime
        if (icalc_otime) {
            if (EDT_use_otime_weight == 2 || icalc_otime_default) { // EDT_OT_WT_ML or otime
                *potime = ot_ml;
            } else { // EDT_OT_WT or EDT
                *potime = ot_sum / ot_weight;
            }
            // normalize weights
            //if (iUseGauss2) {
            NormalizeWeights(num_arrivals, arrival);
            //}
        }
    } else {
        edtSumMisfit = VERY_LARGE_DOUBLE;
        //edtSumMisfit = edtSumMisfit_search = VERY_LARGE_DOUBLE;
        if (icalc_otime)
            *potime = VERY_LARGE_DOUBLE;
    }
    //printf("edt_sum %lf  edtSumMisfit %lf\n", edt_sum, edtSumMisfit);


    // return misfit or ln(prob density)
    if (itype == GRID_MISFIT) {
        /* convert misfit to rms misfit */
        rms_misfit = sqrt(edtSumMisfit + edt_weight); // edt_weight term similar to divide by num readings in rms
        *pmisfit = rms_misfit;
        return (rms_misfit);
    } else if (itype == GRID_PROB_DENSITY) {
#ifdef TEST_COUNT_ONLY_USED_ARRIVALS
        ln_prob_density = -edtSumMisfit * (long double) num_arrivals_used; // best, give power of N
        if (icalc_otime) {
            printf("DEBUG ln_prob_density: num_arrivals %d  num_arrivals_used %d\n", num_arrivals, num_arrivals_used);
        }
#else
        ln_prob_density = -edtSumMisfit * (long double) num_arrivals; // best, give power of N
#endif
        //if (EDT_otime_weight_active)
        ln_prob_density += ot_var_weight;
        // 20200619 AJL - Test weighting by edt_weight (~number of readings) to avoid discontinuities in pdf when travel-times for a phase become unavailable
        //ln_prob_density -= edt_weight;
        // 20200619 AJL - END Test
        /*if (ln_prob_density > 1.0) {
            printf("DEBUG: ln_prob_density %lf = -edtSumMisfit %lf * num_arrivals %d + ot_var_weight %lf\n", ln_prob_density, -(double) edtSumMisfit, num_arrivals, (double) ot_var_weight);
        }*/
        //ln_prob_density = -0.5 * edtSumMisfit;
        rms_misfit = sqrt(edtSumMisfit + edt_weight); // edt_weight term similar to divide by num readings in rms
        *pmisfit = rms_misfit;
        return (ln_prob_density);
    } else {

        return (-1.0);
    }



}







/** function to calculate probability density */

/*	OT_STACK - maximum of stack of otime estimates
 */

double CalcSolutionQuality_OT_STACK(OctNode* poct_node, int num_arrivals, ArrivalDesc *arrival,
        GaussLocParams* gauss_par, int itype, double* pmisfit, double* potime, double* potime_var,
        double cell_half_diagonal_time_range, double cell_diagonal, double cell_volume, double* peffective_cell_size, double *pot_variance_factor) {

    // initialize
    int icalc_otime = 0;
    if (potime != NULL)
        icalc_otime = 1;

    // ot calculation
    double ot_var;
    double prob_max;
    double ot_stack_weight;
    double ot_ml = calc_maximum_likelihood_ot_sort(poct_node, num_arrivals, arrival,
            cell_half_diagonal_time_range, cell_diagonal, cell_volume, &ot_var, icalc_otime, &prob_max, &ot_stack_weight, peffective_cell_size, pot_variance_factor);
    if (icalc_otime && potime_var != NULL) {
        printf("ot_ml_sort_std %lf\n", sqrt(ot_var));
        printf("ot_ml_sort_ot_prob_max %lf\n", prob_max);
        printf("cell_half_diagonal_time_range %lf\n", cell_half_diagonal_time_range);
        *potime_var = ot_var;
    }

    // create misfit
    if (ot_stack_weight > 0.0) {
        // otime
        if (icalc_otime) {
            *potime = ot_ml;
        }
    } else {
        if (icalc_otime)
            *potime = VERY_LARGE_DOUBLE;
    }

    // return misfit or ln(prob density)
    if (itype == GRID_MISFIT) {
        double rms_misfit = sqrt(ot_var);
        *pmisfit = rms_misfit;
        return (rms_misfit);
    } else if (itype == GRID_PROB_DENSITY) {
        //double ln_prob_density = log_ot_prob_max - log(cell_volume); // weight by cell voume
        //double ln_prob_density = log_ot_prob_max - log(cell_diagonal); // weight by cell linear dimension
        if (icalc_otime && potime_var != NULL) {
            printf(">>> prob_max %le   ", prob_max);
            printf(">>> sqrt(ot_var) %lf   ", sqrt(ot_var));
            printf(">>> cell_diagonal %le   ", cell_diagonal);
            printf(">>> cell_volume %le\n", cell_volume);
        }
        double rms_misfit = sqrt(ot_var);
        *pmisfit = rms_misfit;
        return (prob_max);
    } else {

        return (-1.0);
    }



}

/** function to calculate origin time based on sort search for maximum likelihood peak of a set of ot estimates */

double calc_maximum_likelihood_ot_sort(
        OctNode* poct_node, int num_arrivals, ArrivalDesc *arrival,
        double cell_half_diagonal_time_range, double cell_diagonal, double cell_volume, double *pot_var, int icalc_otime,
        double *pprob_max, double *pot_stack_weight, double* peffective_cell_size, double *pot_variance_factor) {

    /*
    // set parent value
    double parent_value = 0.0;
    if (poct_node->parent != NULL) {
        if (poct_node->parent->pdata == NULL) { // should not get here
            fprintf(stderr, "Error: parent OctTree node exists but has no pdata value!\n");
        } else {
            parent_value = *((double *) poct_node->parent->pdata);
        }
    }

    // initalize data
    if (poct_node->pdata == NULL)
        poct_node->pdata = (void *) malloc(sizeof (double));
    if (poct_node->pdata != NULL) {
     *((double *) poct_node->pdata) = 0.0;
    } else {
        fprintf(stderr, "Error: allocating double storage for OctTree node pdata.\n");
    }
     */

    // set ot limits for used arrivals
    double smallest_pick_error = VERY_LARGE_DOUBLE;
    double arr_weight_sum = 0.0;
    //double arr_weight_mean = 0.0;
    //double time_error_mean = 0.0;
    double half_diagonal_time_range = 0.0;

    int narrr_used = 0;
    int narr;
    ArrivalDesc * parr;

    for (narr = 0; narr < num_arrivals; narr++) {
        parr = arrival + narr;
        // skip ignored arrivals
        if (parr->pred_travel_time <= 0.0 || !parr->abs_time)
            continue;
        narrr_used++;
        // set travel time error
        double tt_error;
        if (iUseGauss2) {
            tt_error = parr->pred_travel_time * Gauss2.SigmaTfraction;
            if (tt_error < Gauss2.SigmaTmin)
                tt_error = Gauss2.SigmaTmin;
            if (tt_error > Gauss2.SigmaTmax)
                tt_error = Gauss2.SigmaTmax;
            if (icalc_otime)
                parr->tt_error = tt_error;
            //printf("arrival[nrow].pred_travel_time %f\t  tt_error %f\n", arrival[nrow].pred_travel_time, tt_error);
        } else {
            tt_error = parr->tt_error;
        }
        // set pick eror
        double pick_error = parr->error;
        // set half diagonal time range
        half_diagonal_time_range = cell_half_diagonal_time_range;
        if (parr->slowness > 0.0) {
            half_diagonal_time_range = 0.5 * cell_diagonal * parr->slowness;
            /*static int icount = 0;
            if (narr == 5 && icount++ % 1000 == 0) {
                //printf("depth=%f  slowness_P=%f  slowness_S=%f\n", zval, slowness_P, slowness_S);
                printf("half_diagonal_time_range=%f  cell_diagonal=%f  1/parr->slowness=%f\n", half_diagonal_time_range, cell_diagonal, 1.0 / parr->slowness);
                //printf("model_grid_P.numx=%d  model_grid_P.dy=%f  yval_grid=%f\n", model_grid_P.numx, model_grid_P.dy, yval_grid);
            }*/

        }
        // set ot limits time range
        double time_range = half_diagonal_time_range + tt_error + pick_error;
        // set arrival weights
        //double weight = 1.0 / time_range;
        //double weight = 1.0 / (2.0 * (tt_error + pick_error));
        //double tt_error_ref = iUseGauss2 ? Gauss2.SigmaTmin : Gauss.SigmaT;
        //weight = (octtreeParams.min_node_size / cell_diagonal) * (tt_error_ref / (tt_error + pick_error));  // !!! TEST
        //double weight = tt_error_ref / (tt_error + pick_error); // !!! TEST
        double weight = half_diagonal_time_range / time_range; // !!! TEST - independent of tt_error_ref, but changes as a function of cell size and slowness
        weight = 1.0;
        if (iSetStationDistributionWeights) {
            //if (parr->station_weight < 0.5 || parr->station_weight > 2.0)
            //    printf("NOTE!!!!!!!!!!! parr->station_weight < 0.5 || parr->station_weight > 2.0 %s %f\n", parr->label, parr->station_weight);
            weight *= parr->station_weight;
        }
        // 20130627 AJL - add prior weighting
        if (iUseArrivalPriorWeights && parr->apriori_weight >= -VERY_SMALL_DOUBLE)
            weight *= parr->apriori_weight;
        parr->weight = weight;
        /*{
            static int icount = 0;
            if (narr == 5 && icount++ % 1000 == 0)
                printf("half_diagonal_time_range=%f  tt_error=%f  pick_error=%f  parr->weight=%f\n", half_diagonal_time_range, tt_error, pick_error, weight);
        }*/
        arr_weight_sum += weight;
        // set OtimeLimits
        double ot_arr = parr->obs_time - (long double) parr->pred_travel_time;
        double dist_range = 0.0;
        if (parr->slowness > 0.0)
            dist_range = 2.0 * time_range / parr->slowness;
        OtimeLimit * otimeLimitMin = new_OtimeLimit(narr, ot_arr - time_range, ot_arr, 1, dist_range, 2.0 * time_range);
        OtimeLimit * otimeLimitMax = new_OtimeLimit(narr, ot_arr + time_range, ot_arr, -1, dist_range, 2.0 * time_range);
        //otimeLimitMin->pair = otimeLimitMax;
        //otimeLimitMax->pair = otimeLimitMin;
        addOtimeLimitToList(otimeLimitMin, &OtimeLimitList, &NumOtimeLimit);
        addOtimeLimitToList(otimeLimitMax, &OtimeLimitList, &NumOtimeLimit);
        // set smallest pick error
        if (pick_error < smallest_pick_error)
            smallest_pick_error = pick_error;
        // accumulate total time error
        //time_error_mean += tt_error + pick_error;
        //time_error_mean += time_range;
    }
    //time_error_mean /= (double) narrr_used;

    /*{
        static int icount = 0;
        if (icount++ % 1000 == 0)
            printf("arr_weight_sum=%f  narrr_used=%d\n", arr_weight_sum, narrr_used);
    }*/
    // normalize weights
    for (narr = 0; narr < num_arrivals; narr++) {
        parr = arrival + narr;
        // skip ignored arrivals
        if (parr->pred_travel_time <= 0.0 || !parr->abs_time)
            continue;
        parr->weight = (double) narrr_used * parr->weight / arr_weight_sum;
    }
    arr_weight_sum = (double) narrr_used;
    //arr_weight_mean = arr_weight_sum / (double) narrr_used;

    // find max of otime limit histogram
    int narrival_stack = 0;
    int best_nstation = 0;
    double ot_sum = 0.0;
    double ot_sum_sqr = 0.0;
    double time_range_sum_sqr = 0.0;
    double weight_sum = 0.0;
    double best_weight_sum = 0.0;
    double best_ot_mean = 0.0;
    double best_prob = -LARGE_DOUBLE;
    double best_ot_variance = -1.0;
    double dist_range_sum = 0.0;
    double best_dist_range_sum = 0.0;
    double best_ot_variance_factor = 0.0;
    int i;
    OtimeLimit* otimeLimit;
    int data_id;
    double otime, weight;
    for (i = 0; i < NumOtimeLimit; i++) { // parse otime limits in time order
        otimeLimit = *(OtimeLimitList + i);
        data_id = otimeLimit->data_id;
        otime = otimeLimit->otime;
        weight = arrival[data_id].weight;
        if (otimeLimit->polarity > 0) { // enter otime limit for this datum
            ot_sum += weight * otime;
            ot_sum_sqr += weight * otime * otime;
            weight_sum += weight;
            dist_range_sum += weight * otimeLimit->dist_range;
            time_range_sum_sqr += weight * otimeLimit->time_range * otimeLimit->time_range;
            narrival_stack++;
        } else { // leave otime limit for this datum
            ot_sum -= weight * otime;
            ot_sum_sqr -= weight * otime * otime;
            weight_sum -= weight;
            dist_range_sum -= weight * otimeLimit->dist_range;
            time_range_sum_sqr -= weight * otimeLimit->time_range * otimeLimit->time_range;
            narrival_stack--;
            // check if not enough data remaining to get more stations than best
            // i.e. nstation + max_new_possible < best_nstation
            //if (nstation + (NumOtimeLimit - i - 1 - nstation) / 2 < best_nstation)
            //    break;
        }
        // need at least 2 stations for location
        // WORK_POINT
        double min_weight_sum_assoc = 2.0;
        min_weight_sum_assoc += 0.01; // prevent divide by small values below
        //if (narrival_stack > 1 && weight_sum > best_weight_sum && weight_sum > arr_weight_mean) {
        if (narrival_stack > 1 && weight_sum > min_weight_sum_assoc) {
            //if (narrival_stack > 1 && weight_sum >= min_weight_sum_assoc) { // 20101217  this is used in warning monitor for speed and efficiency in convergence, with the risk of less thorough search

            double ot_mean = ot_sum / weight_sum;
            //double ot_variance = (ot_sum_sqr - weight_sum * ot_mean * ot_mean) / weight_sum;
            double ot_variance = (ot_sum_sqr - weight_sum * ot_mean * ot_mean) / (weight_sum - min_weight_sum_assoc + 1.0);
            double time_range_variance = time_range_sum_sqr / (weight_sum - 2.0); // 20101224 AJL - changed  - 1.0  to  - 2.0
            double ot_variance_factor = exp(-ot_variance / time_range_variance);
            double adjusted_weight_sum = weight_sum - 1.0;
            double prob = ot_variance_factor * adjusted_weight_sum;

            double best_effective_cell_size = dist_range_sum / weight_sum;
            double effective_cell_volume = pow(best_effective_cell_size, 3);
            if (effective_cell_volume < cell_volume) // 20101217
                effective_cell_volume = cell_volume;

            prob -= log(effective_cell_volume); // division by effective_cell_volume converts prob to prob density

            // need at least 2 stations for location
            //if (nassociated_P_work >= best_nassociated_P_work && nassociated_P_work > 1 && weight_sum > 1.0) {
            //    if (nassociated_P_work > best_nassociated_P_work || weight_sum > weight_sum) { // use weight_sum to select between solutions with same num P
            //if (nassociated_P_work > 1 && weight_sum > weight_sum && weight_sum > 1.0) { // same as NLL OT_STACK
            if (prob > best_prob) { // same as NLL OT_STACK

                // save raw value
                if (poct_node->pdata != NULL) {
                    *((double *) poct_node->pdata) = weight_sum - 1.0;
                }
                best_ot_mean = ot_mean;
                best_ot_variance = ot_variance;
                best_prob = prob;
                best_ot_variance_factor = ot_variance_factor;
                best_nstation = narrival_stack;
                best_weight_sum = weight_sum;
                best_dist_range_sum = dist_range_sum;
                //best_weight_sum -= log((best_ot_variance + time_error_mean) / time_error_mean);
                /*static int icount = 0;
                if (icount++ % 100 == 0) {
                    printf("\nbest_log_prob_max=%lg,  ot_variance %f\n", best_log_prob_max, best_ot_variance);
                }
                 */
            }
        }
    }

    /*static int icount = 0;
    if (icount++ % 1000 == 0) {
        printf("best_ot_mean=%lg,  best_ot_variance %f  best_weight_sum %f  best_log_prob_max %f\n", best_ot_mean, best_ot_variance, best_weight_sum, best_log_prob);
    }*/



    free_OtimeLimitList(&OtimeLimitList, &NumOtimeLimit);



    // set returned reference values
    *pprob_max = best_prob;
    *pot_var = best_ot_variance;
    *pot_stack_weight = best_weight_sum;
    *peffective_cell_size = best_dist_range_sum / best_weight_sum;
    *pot_variance_factor = best_ot_variance_factor;

    /*DEBUG*/
    if (icalc_otime) {
        printf("=================\nNumOtimeLimit %d  ", i);
        printf("cell_half_diagonal_time_range=%e  ", cell_half_diagonal_time_range);
        printf("half_diagonal_time_range=%e  ", half_diagonal_time_range);
        printf("best_nstation=%d  ", best_nstation);
        printf("best_weight_sum=%f  ", best_weight_sum);
        printf("ot_mean=%f  ", best_ot_mean);
        printf("best_log_prob_max=%f  ", best_prob);
        printf("best_ot_variance=%f  ", best_ot_variance);
        printf("effective_cell_size=%f  ", *peffective_cell_size);
        printf("\n");
    }

    if (icalc_otime && best_nstation < 2)
        nll_puterr("ERROR: calc_maximum_likelihood_ot_stack: best_nstation < 2.");

    return (best_ot_mean);


}





/** function to calculate probability density */
/*	ML_OT - maximum of sum of probabilities in otime domain
 */

#define ML_OT_WT_FLOOR log(0.00001)

double CalcSolutionQuality_ML_OT(int num_arrivals, ArrivalDesc *arrival,
        GaussLocParams* gauss_par, int itype, double* pmisfit, double* potime,
        double* potime_var, double cell_half_diagonal_time_range, int method_box) {

    double cell_diagonal_time_var = cell_half_diagonal_time_range * cell_half_diagonal_time_range;

    int nrow;

    long double edt_sum, prob, edtSumMisfit;

    double edt_weight;
    double weight, weight2;
    double ln_prob_density, rms_misfit;

    MatrixDouble edtmtx;
    double sigma2_row;
    //double obs_minus_pred;

    int icalc_otime = 0;

    // EDT_OT_WT
    int num_otime_error;
    long double ot_row, ot_prob, ot_row_2, ot_error_2;
    long double ot_sum, ot_weight, ot_var, ot_2_sum, ot_var_weight;

    // EDT_OT_WT_ML
    double ot_ml = 0.0, ot_ml_var;

    // OT additions 20071220
    double ot_prob_max = 0.0;

    // Gauss2
    double tt_error;


    // search pdf different from true
    int iuse_cell_diagonal_time_var;
    //double sigma2_row_search = 0.0;
    //double weight_search = 0.0, weight2_search, edt_weight_search;
    //long double prob_search = 0.0L, edt_sum_search = 0.0L, edtSumMisfit_search = 0.0L;


    // initialize
    if (potime != NULL)
        icalc_otime = 1;
    edtmtx = gauss_par->EDTMtrx;

    // check if use_cell_diagonal_time_var
    iuse_cell_diagonal_time_var = 0;
    if (cell_diagonal_time_var > 0.0)
        iuse_cell_diagonal_time_var = 1;

    // check size of EDT_OT_WT_ML static arrays
    if ((EDT_use_otime_weight == 2 || icalc_otime)) {
        if (isize_ot_ml_array < num_arrivals) {
            isize_ot_ml_array = num_arrivals;
            free(ot_ml_arrival);
            if ((ot_ml_arrival = (double *) calloc(isize_ot_ml_array, sizeof (double))) == NULL)
                nll_puterr("ERROR: allocating double storage array for EDT_OT_WT_ML ot_ml_arrival.");
            free(ot_ml_arrival_edt_sum);
            if ((ot_ml_arrival_edt_sum = (double *) calloc(isize_ot_ml_array, sizeof (double))) == NULL)
                nll_puterr("ERROR: allocating double storage array for EDT_OT_WT_ML ot_ml_arrival_edt_sum.");
        }
        for (nrow = 0; nrow < num_arrivals; nrow++)
            ot_ml_arrival_edt_sum[nrow] = 0.0;
    }

    if (icalc_otime) {
        for (nrow = 0; nrow < num_arrivals; nrow++)
            arrival[nrow].weight = 0.0;
    }


    /* calculate weighted mean of predicted travel times  */
    /*		(TV82, eq. A-38) */
    CalcCenteredTimesPred(num_arrivals, arrival, gauss_par); // not used for EDT


    /* calculate EDT prop sum */

    edt_sum = 0.0L;
    edt_weight = 0.0;
    ot_prob = 0.0L;
    ot_row = 0.0L;
    ot_row_2 = 0.0L;
    ot_sum = 0.0L;
    ot_2_sum = 0.0L;
    ot_var = 0.0L;
    ot_var_weight = 0.0L;
    ot_weight = 0.0L;
    ot_error_2 = 0.0L;
    num_otime_error = 0;
    for (nrow = 0; nrow < num_arrivals; nrow++) {

        // AJL 20041115 bug fix!
        if (arrival[nrow].pred_travel_time <= 0.0
                || !arrival[nrow].abs_time) // AJL 20071220
        {
            // iniitalize EDT_OT_WT_ML values
            if (EDT_use_otime_weight == 2 || icalc_otime) {
                ot_ml_arrival_edt_sum[nrow] = -1.0;
            }
            continue; // ignore obs without predicted times
        }
        // END

        // set error
        // printf("iUseGauss2 %d\n", iUseGauss2);
        if (iUseGauss2) {
            tt_error = arrival[nrow].pred_travel_time * Gauss2.SigmaTfraction;
            if (tt_error < Gauss2.SigmaTmin)
                tt_error = Gauss2.SigmaTmin;
            if (tt_error > Gauss2.SigmaTmax)
                tt_error = Gauss2.SigmaTmax;
            if (icalc_otime)
                arrival[nrow].tt_error = tt_error;
            //printf("arrival[nrow].pred_travel_time %f\t  tt_error %f\n", arrival[nrow].pred_travel_time, tt_error);
            tt_error *= tt_error;
            edtmtx[nrow][nrow] = arrival[nrow].error * arrival[nrow].error + tt_error;
        }
        sigma2_row = edtmtx[nrow][nrow];

        if (iuse_cell_diagonal_time_var)
            sigma2_row += cell_diagonal_time_var;
        /*if (iuse_cell_diagonal_time_var) {
        sigma2_row_search = edtmtx[nrow][nrow] + cell_diagonal_time_var;
}*/
        //error_row = arrival[nrow].error;
        //obs_minus_pred = arrival[nrow].obs_centered - arrival[nrow].pred_centered;
        if (EDT_use_otime_weight == 2 || icalc_otime) { // EDT_OT_WT_ML or otime
            ot_ml_arrival[nrow] = arrival[nrow].obs_time - (long double) arrival[nrow].pred_travel_time;
            //ot_ml_arrival_edt_sum[nrow] = 0.0;
            ot_error_2 += sigma2_row;
            num_otime_error++;
        } else if (EDT_use_otime_weight == 1) { // EDT_OT_WT or EDT
            ot_prob = 0.0;
            ot_row = arrival[nrow].obs_time - (long double) arrival[nrow].pred_travel_time;
            ot_row_2 = ot_row * ot_row;
            ot_error_2 += sigma2_row;
            num_otime_error++;
        }

        weight2 = 1.0 / sigma2_row; // sum of errors**2
        weight = sqrt(weight2); // errors factor

        if (iSetStationDistributionWeights)
            weight *= arrival[nrow].station_weight;
        // 20130627 AJL - add prior weighting
        if (iUseArrivalPriorWeights && arrival[nrow].apriori_weight >= -VERY_SMALL_DOUBLE)
            weight *= arrival[nrow].apriori_weight;

        prob = weight;
        edt_sum += prob;
        edt_weight += weight;

        // accumulate EDT weights
        if (icalc_otime) {
            arrival[nrow].weight += prob;
        }
        // otime
        if (EDT_use_otime_weight == 2 || icalc_otime) { // EDT_OT_WT_ML or otime
            //ot_ml_arrival_edt_sum[nrow] += prob;
            ot_ml_arrival_edt_sum[nrow] = 1.0; // AJL 20071220 not needed?  redundant with normalisation in
        } else if (EDT_use_otime_weight == 1) { // EDT_OT_WT or EDT
            ot_prob += prob;
        }

        if (EDT_use_otime_weight == 1) { // EDT_OT_WT or EDT
            ot_sum += ot_prob * ot_row;
            ot_2_sum += ot_prob * ot_row_2;
            ot_weight += ot_prob;
        }
    }

    // OT_WT methods
    if (EDT_use_otime_weight == 2 || icalc_otime) { // EDT_OT_WT_ML
        // EDT_OT_WT_ML method
        ot_ml = calc_maximum_likelihood_ot(ot_ml_arrival, ot_ml_arrival_edt_sum, num_arrivals, arrival, edtmtx, &ot_ml_var, icalc_otime, &ot_prob_max);
        if (icalc_otime && potime_var != NULL) {
            printf("ot_ml_std %lf\n", sqrt(ot_ml_var));
            *potime_var = ot_ml_var;
        }
        ot_var_weight = -ot_ml_var / (ot_error_2 / (long double) num_otime_error);
        ot_var_weight = log(ot_prob_max);
        if (1 || ot_var_weight > ML_OT_WT_FLOOR) {
            if (!EDT_otime_weight_active) {
                EDT_otime_weight_active = 1;
                snprintf(MsgStr, sizeof(MsgStr), "EDT_otime_weight activated, OT_WT exceeds ML_OT_WT_FLOOR.");
                nll_putmsg(2, MsgStr);
            }
        } else {
            ot_var_weight = ML_OT_WT_FLOOR;
        }
    } else if (EDT_use_otime_weight == 1 && num_otime_error > 0) { // EDT_OT_WT
        // EDT_OT_WT method
        ot_var = ot_2_sum / ot_weight - (ot_sum / ot_weight) * (ot_sum / ot_weight);
        if (icalc_otime && potime_var != NULL) {
            printf("ot_ml_std %lf\n", sqrt(ot_var));
            *potime_var = ot_var;
        }
        ot_var_weight = -ot_var / (ot_error_2 / (long double) num_otime_error);
        if (1 || ot_var_weight > ML_OT_WT_FLOOR) {
            if (!EDT_otime_weight_active) {
                EDT_otime_weight_active = 1;
                snprintf(MsgStr, sizeof(MsgStr), "EDT_otime_weight activated, OT_WT exceeds ML_OT_WT_FLOOR.");
                nll_putmsg(2, MsgStr);
            }
        } else {
            ot_var_weight = ML_OT_WT_FLOOR;
        }
        //edt_sum *= ot_var_weight;
        //edt_weight *= ot_var_weight;
    }

    // avoid numerical problems
    if (edt_sum < SMALL_DOUBLE)
        edt_sum = SMALL_DOUBLE;
    /*if (edt_sum_search < SMALL_DOUBLE)
    edt_sum_search = SMALL_DOUBLE;*/
    // create EDT misfit
    if (edt_weight > 0.0 && num_arrivals > 0) {
        // AJL 20051115 bug fix!
        //edtSumMisfit = -log(edt_sum / edt_weight);
        edtSumMisfit = -log(edt_sum); // non-normalized sum of edt probs - favors solutions with more readings
        // END
        //edtSumMisfit *= 2.0;
        //edtSumMisfit *= num_arrivals;	// best, give power of N
        //edtSumMisfit *= sqrt(num_arrivals);	// works much better for INGV Italy region scale
        /*if (iuse_cell_diagonal_time_var) { 	// duplicate above lines
        edtSumMisfit_search = -log(edt_sum_search);	// non-normalized sum of edt probs - favors solutions with more readings
}*/
        // otime
        if (icalc_otime) {
            if (EDT_use_otime_weight == 2 || icalc_otime) { // EDT_OT_WT_ML or otime
                *potime = ot_ml;
            } else { // EDT_OT_WT or EDT
                *potime = ot_sum / ot_weight;
            }
            // normalize weights
            //if (iUseGauss2) {
            NormalizeWeights(num_arrivals, arrival);
            //}
        }
    } else {
        edtSumMisfit = VERY_LARGE_DOUBLE;
        //edtSumMisfit = edtSumMisfit_search = VERY_LARGE_DOUBLE;
        if (icalc_otime)
            *potime = VERY_LARGE_DOUBLE;
    }
    //printf("edt_sum %lf  edtSumMisfit %lf\n", edt_sum, edtSumMisfit);


    // return misfit or ln(prob density)
    if (itype == GRID_MISFIT) {
        /* convert misfit to rms misfit */
        rms_misfit = sqrt(edtSumMisfit + edt_weight); // edt_weight term similar to divide by num readings in rms
        *pmisfit = rms_misfit;
        return (rms_misfit);
    } else if (itype == GRID_PROB_DENSITY) {
        ln_prob_density = -edtSumMisfit * num_arrivals; // best, give power of N
        //if (EDT_otime_weight_active)
        //ln_prob_density += ot_var_weight;
        ln_prob_density = ot_var_weight * num_arrivals;
        //ln_prob_density = -0.5 * edtSumMisfit;
        rms_misfit = sqrt(edtSumMisfit + edt_weight); // edt_weight term similar to divide by num readings in rms
        *pmisfit = rms_misfit;
        return (ln_prob_density);
    } else {

        return (-1.0);
    }



}

/** function to calculate origin time based on grid search for maximum likelihood peak of a set of ot estimates */

double calc_maximum_likelihood_ot(double *pot_ml_arrival, double *pot_ml_arrival_edt_sum,
        int num_arrivals, ArrivalDesc *arrival, MatrixDouble edtmtx, double *pot_ml_var, int iwrite_errors,
        double *pprob_max) {

    int narr;
    double arr_prob_max, prob;
    double ot_arr, ot_arr_prob_max = 0.0;
    double range, step, time, tlimit, prob_max, ot_max_like;
    double edt_matrix_rms_error, edt_matrix_diagonal_sum;


    // coarse search: find ot estimate with largest likelihood over set of ot's for each arrival
    arr_prob_max = -1.0;
    edt_matrix_diagonal_sum = 0.0;
    for (narr = 0; narr < num_arrivals; narr++) {
        if (pot_ml_arrival_edt_sum[narr] < 0.0) // skip ignored arrivals
            continue;
        ot_arr = pot_ml_arrival[narr];
        prob = calc_likelihood_ot(ot_ml_arrival, pot_ml_arrival_edt_sum, num_arrivals, arrival, edtmtx, ot_arr);
        //printf("DEBUG: calc_likelihood_ot arr: narr %d, ot_ml_arrival %f, prob %f, edtmtx[narr][narr] %f ", narr, ot_arr, prob, edtmtx[narr][narr]);
        // store if largest cumulative ot likelihood
        if (prob > arr_prob_max) {
            arr_prob_max = prob;
            ot_arr_prob_max = ot_arr;
            //printf("---------> max");
        }
        //printf("\n");
        edt_matrix_diagonal_sum += edtmtx[narr][narr];
    }
    if (iwrite_errors && arr_prob_max < 0.0)
        nll_puterr("ERROR: calc_maximum_likelihood_ot: failed to find arr_prob_max.");

    edt_matrix_rms_error = sqrt(edt_matrix_diagonal_sum / num_arrivals);


    // search numerically for precise ot likelihood peak
    range = 3.0 * edt_matrix_rms_error;
    step = edt_matrix_rms_error / 100.0;
    prob_max = arr_prob_max;
    ot_max_like = ot_arr_prob_max;
    // search increasing time
    time = ot_arr_prob_max;
    tlimit = ot_arr_prob_max + range;
    while ((time = time + step) < tlimit) {
        prob = calc_likelihood_ot(ot_ml_arrival, pot_ml_arrival_edt_sum, num_arrivals, arrival, edtmtx, time);
        //printf("DEBUG: calc_likelihood_ot: inc time: time %f, prob %f\n", time, prob);
        if (prob < prob_max)
            break;
        prob_max = prob;
        ot_max_like = time;
        //printf("DEBUG: calc_likelihood_ot 1: ot_max_like %f\n", ot_max_like);
    }
    if (iwrite_errors && time >= tlimit) {
        snprintf(MsgStr, sizeof(MsgStr), "ot_arr_prob_max: %f, range %f, tlimit %f", ot_arr_prob_max, range, tlimit);
        nll_puterr2("ERROR: calc_maximum_likelihood_ot: reached end of increasing-time search limit:", MsgStr);
    }

    // search decreasing time
    time = ot_arr_prob_max;
    tlimit = ot_arr_prob_max - range;
    while ((time = time - step) > tlimit) {
        prob = calc_likelihood_ot(ot_ml_arrival, pot_ml_arrival_edt_sum, num_arrivals, arrival, edtmtx, time);
        //printf("DEBUG: calc_likelihood_ot: dec time: time %f, prob %f\n", time, prob);
        if (prob < prob_max)
            break;
        prob_max = prob;
        ot_max_like = time;
        //printf("DEBUG: calc_likelihood_ot 2: ot_max_like %f\n", ot_max_like);
    }
    if (iwrite_errors && time <= tlimit) {
        snprintf(MsgStr, sizeof(MsgStr), "ot_arr_prob_max: %f, range %f, tlimit %f", ot_arr_prob_max, range, tlimit);
        nll_puterr2("ERROR: calc_maximum_likelihood_ot: reached end of decreasing-time search limit:", MsgStr);
    }


    // set prob_max - added AJL 20071220
    *pprob_max = prob_max;

    // get variance
    *pot_ml_var = calc_variance_ot(ot_ml_arrival, pot_ml_arrival_edt_sum, num_arrivals, arrival, edtmtx, ot_max_like);

    //printf("DEBUG: calc_likelihood_ot 3: ot_max_like %f\n", ot_max_like);

    return (ot_max_like);


}

/** function to calculate origin time likelihood at at given time */

double calc_likelihood_ot(double *pot_ml_arrival, double *pot_ml_arrival_edt_sum, int num_arrivals, ArrivalDesc *arrival, MatrixDouble edtmtx, double time) {

    int narr1;
    ArrivalDesc *parr1;
    double prob, prob_arr1, sigma2, temp;

    prob = 0.0;
    for (narr1 = 0; narr1 < num_arrivals; narr1++) {
        if (pot_ml_arrival_edt_sum[narr1] < 0.0) // skip ignored arrivals
            continue;
        parr1 = arrival + narr1;
        sigma2 = edtmtx[narr1][narr1];
        temp = pot_ml_arrival[narr1] - time;
        //printf("DEBUG: calc_likelihood_ot: narr %d, temp %f = pot_ml_arrival[narr1] %f - time %f", narr1, temp, pot_ml_arrival[narr1], time);
        // 20121005 AJL - bug fix, check for valid temp
        if (temp > -1.0e8 && temp < 1.0e8) {
            // normal dist (normalized relative) of arrival error around ot est for arrival
            prob_arr1 = exp(-0.5 * temp * temp / sigma2) / sqrt(sigma2);
            //printf(", -> A prob_arr1 %f", prob_arr1);
            // weight by edt_prob for arrival
            if (num_arrivals > 1) { // no edt prob if < 2 arrivals  20190826 AJL - bug fix when locating with only one arrival (e.g. LOCPOSTERIOR)
                prob_arr1 *= pot_ml_arrival_edt_sum[narr1];
            }
            //printf(", -> B prob_arr1 %f", prob_arr1);
            // weight by station distribution
            // 20161013 AJL - bug fix: removed: pot_ml_arrival_edt_sum already includes station dist weight!
            //if (iSetStationDistributionWeights)
            //    prob_arr1 *= parr1->station_weight;
            //printf(", swt %f -> C prob_arr1 %f", parr1->station_weight, prob_arr1);
            // 20130627 AJL - add prior weighting
            if (iUseArrivalPriorWeights && parr1->apriori_weight >= -VERY_SMALL_DOUBLE)
                prob_arr1 *= parr1->apriori_weight;
        } else {

            prob_arr1 = 0.0;
        }
        //printf(", -> D prob_arr1 %f\n", prob_arr1);

        // cumulate sum of ot likelihood
        prob += prob_arr1;
    }

    return (prob);


}

/** function to calculate origin time variance for at given time */

double calc_variance_ot(double *pot_ml_arrival, double *pot_ml_arrival_edt_sum, int num_arrivals, ArrivalDesc *arrival, MatrixDouble edtmtx, double expectation_time) {

    int narr1;
    double variance_sum, temp, sigma2, weight, weight_sum;


    variance_sum = 0.0;
    weight_sum = 0.0;
    for (narr1 = 0; narr1 < num_arrivals; narr1++) {
        if (pot_ml_arrival_edt_sum[narr1] < 0.0) // skip ignored arrivals
            continue;
        temp = pot_ml_arrival[narr1] - expectation_time;
        // normal dist (normalized relative) of arrival error around ot est for arrival
        sigma2 = edtmtx[narr1][narr1];
        weight = 1.0 / sqrt(sigma2);
        // weight by edt_prob for arrival
        if (num_arrivals > 1) { // no edt prob if < 2 arrivals  20190826 AJL - bug fix when locating with only one arrival (e.g. LOCPOSTERIOR)
            weight *= pot_ml_arrival_edt_sum[narr1];
        }
        // weight by station distribution

        // 20161013 AJL - bug fix: removed: pot_ml_arrival_edt_sum already includes station dist weight!
        //if (iSetStationDistributionWeights)
        //weight *= arrival[narr1].station_weight;
        // 20130627 AJL - add prior weighting

        if (iUseArrivalPriorWeights && arrival[narr1].apriori_weight >= -VERY_SMALL_DOUBLE)
            weight *= arrival[narr1].apriori_weight;
        variance_sum += weight * temp * temp;
        weight_sum += weight;
    }

    return (variance_sum / weight_sum);


}

/** function to normalize arrival weights */

int NormalizeWeights(int num_arrivals, ArrivalDesc * arrival) {
    int narr;
    //double sigmaT2, corr_len2;
    //int corr_len_nonzero = 1;
    //double dx, dy, dz, dist2;
    double weight_sum;
    //double sta_wt;
    //MatrixDouble null_mtrx = NULL;
    //SourceDesc *sta1, *sta2;
    //double arrivalWeightMax = -1.0;

    //static MatrixDouble wt_matrix, edt_matrix = NULL;
    //static int last_matrix_alloc_size = -1;


    // get row weights & sum of weights

    weight_sum = 0.0;
    for (narr = 0; narr < num_arrivals; narr++) {
        weight_sum += arrival[narr].weight;
        //printf("row %d col %d: wt_tx(r,c) %f  arr(row)_wt %f   wt_sum %lf\n", nrow, ncol, wt_matrix[nrow][ncol], arrival[nrow].weight, weight_sum);
    }
    for (narr = 0; narr < num_arrivals; narr++) {
        arrival[narr].weight = (double) num_arrivals * arrival[narr].weight / weight_sum;
        //printf("observation weight: %s %s %s weight: %lf\n", arrival[nrow].label, arrival[nrow].inst, arrival[nrow].comp, arrival[nrow].weight);
    }
    if (message_flag >= 4) {

        snprintf(MsgStr, sizeof(MsgStr), "EDT Posterior Weight Matrix sum: %f", weight_sum);
        nll_putmsg(4, MsgStr);
    }


    return (0);

}




/** function to calculate probability density using L1 norm */

// 20150323 AJL - added

double CalcSolutionQuality_L1_NORM(int num_arrivals, ArrivalDesc *arrival,
        GaussLocParams* gauss_par, int itype, double* pmisfit, double* potime) {

    int nrow, ncol, narr;
    int narr_misfit;

    double misfit = 0.0;
    double ln_prob_density, l1_misfit;

    double arr_row_res;
    MatrixDouble wtmtx;
    double *wtmtxrow;

    wtmtx = gauss_par->WtMtrx;


    /* calculate weighted mean of predicted travel times  */
    /*		(TV82, eq. A-38) */

    CalcCenteredTimesPred(num_arrivals, arrival, gauss_par);


    /* calculate residuals (TV82, eqs. 10-12, 10-13; MEN92, eq. 15) */

    for (narr = 0; narr < num_arrivals; narr++) {
        // AJL 20041115 bug fix!
        if (arrival[narr].pred_travel_time <= 0.0)
            arrival[narr].cent_resid = 0.0; // ignore obs without predicted times
        else
            arrival[narr].cent_resid = arrival[narr].obs_centered - arrival[narr].pred_centered;
    }

    narr_misfit = 0;
    for (nrow = 0; nrow < num_arrivals; nrow++) {
        // AJL 20041115 bug fix!
        if (arrival[nrow].pred_travel_time <= 0.0) {
            //printf("IGNORE: %s %s\n", arrival[nrow].label, arrival[nrow].phase);
            continue; // ignore obs without predicted times
        }
        // END
        if (!arrival[nrow].abs_time)
            continue;
        narr_misfit++;
        wtmtxrow = wtmtx[nrow];
        arr_row_res = arrival[nrow].cent_resid;
        for (ncol = 0; ncol <= nrow; ncol++) {
            // AJL 20041115 bug fix!
            if (arrival[ncol].pred_travel_time <= 0.0)
                continue; // ignore obs without predicted times
            // END
            if (!arrival[ncol].abs_time)
                continue;
            if (ncol != nrow) {
                // misfit += 2.0 * (double) *(wtmtxrow + ncol) * arr_row_res * arrival[ncol].cent_resid;    // L2  METH_GAU_ANALYTIC
                misfit += 2.0 * (double) *(wtmtxrow + ncol) * sqrt(fabs(arr_row_res * arrival[ncol].cent_resid)); // L1  METH_L1_NORM  20150324 AJL - TODO: should be sqrt(PRODUCT) or (SUM)/2 ???
            } else {
                //misfit += (double) *(wtmtxrow + ncol) * arr_row_res * arrival[ncol].cent_resid;    // L2  METH_GAU_ANALYTIC
                misfit += (double) *(wtmtxrow + ncol) * fabs(arr_row_res); // L1  METH_L1_NORM
            }
        }
    }


    if (potime != NULL)
        *potime = CalcMaxLikeOriginTime(num_arrivals, arrival, gauss_par);

    /* return misfit or ln(prob density) */

    if (itype == GRID_MISFIT) {
        /* convert misfit to rms misfit */
        if (narr_misfit > 0) { // 20120724 AJL - bug fix!
            //rms_misfit = sqrt(misfit / narr_misfit);    // L2  METH_GAU_ANALYTIC
            l1_misfit = misfit / narr_misfit; // L1  METH_L1_NORM
        } else {
            l1_misfit = VERY_LARGE_DOUBLE;
        }
        *pmisfit = l1_misfit;
        return (l1_misfit);
    } else if (itype == GRID_PROB_DENSITY) {
        if (narr_misfit > 0) { // 20120724 AJL - bug fix!
            //ln_prob_density = -0.5 * misfit;    // L2  METH_GAU_ANALYTIC
            ln_prob_density = -1.0 * misfit; // L1  METH_L1_NORM
            //rms_misfit = sqrt(misfit / narr_misfit);    // L2  METH_GAU_ANALYTIC
            l1_misfit = misfit / narr_misfit; // L1  METH_L1_NORM
        } else {
            ln_prob_density = -VERY_LARGE_DOUBLE;
            l1_misfit = VERY_LARGE_DOUBLE;
        }
        *pmisfit = l1_misfit;
        return (ln_prob_density);
    } else {

        return (-1.0);
    }



}




/** function to calculate probability density */

/*	sum of individual L2 residual probablities */

double CalcSolutionQuality_GAU_TEST(int num_arrivals, ArrivalDesc *arrival,
        GaussLocParams* gauss_par, int itype, double* pmisfit, double* potime) {

    int nrow, ncol, narr;

    double misfit;
    double ln_prob_density, rms_misfit;

    double arr_row_res;
    MatrixDouble wtmtx;
    double *wtmtxrow;

    double pred_min = VERY_LARGE_DOUBLE, pred_max = -VERY_LARGE_DOUBLE;
    double prob_max = -VERY_LARGE_DOUBLE;
    double tshift, tshift_at_prob_max = 0.0;
    double tstep, tstart, tstop;
    double misfit_min = VERY_LARGE_DOUBLE;
    double misfit_tmp, prob;


    wtmtx = gauss_par->WtMtrx;


    /* calculate weighted mean of predicted travel times  */
    /*		(TV82, eq. A-38) */

    CalcCenteredTimesPred(num_arrivals, arrival, gauss_par);


    /* calculate residuals (TV82, eqs. 10-12, 10-13; MEN92, eq. 15) */

    for (narr = 0; narr < num_arrivals; narr++) {
        // AJL 20041115 bug fix!
        if (arrival[narr].pred_travel_time <= 0.0)
            arrival[narr].cent_resid = 0.0; // ignore obs without predicted times
        else {
            arrival[narr].cent_resid = arrival[narr].obs_centered - arrival[narr].pred_centered;
            if (arrival[narr].pred_centered < pred_min)
                pred_min = arrival[narr].pred_centered;
            if (arrival[narr].pred_centered > pred_max)
                pred_max = arrival[narr].pred_centered;
        }
    }
    //printf("pred_min %lf  pred_max %lf\n", pred_min, pred_max);


    // find maximum prob time shift

    tstep = (pred_max - pred_min) / 10.0;
    tstart = pred_min;
    tstop = pred_max;
    while (tstep > (pred_max - pred_min) / 1000000.0) {
        for (tshift = tstart; tshift <= tstop; tshift += tstep) {
            misfit = 0.0;
            prob = 0.0;
            for (nrow = 0; nrow < num_arrivals; nrow++) {
                // AJL 20041115 bug fix!
                if (arrival[nrow].pred_travel_time <= 0.0)
                    continue; // ignore obs without predicted times
                // END
                wtmtxrow = wtmtx[nrow];
                arr_row_res = arrival[nrow].cent_resid + tshift;
                for (ncol = 0; ncol <= nrow; ncol++) {
                    // AJL 20041115 bug fix!
                    if (arrival[ncol].pred_travel_time <= 0.0)
                        continue; // ignore obs without predicted times
                    // END
                    if (ncol != nrow) {
                        misfit_tmp = (double) *(wtmtxrow + ncol) *
                                arr_row_res * (arrival[ncol].cent_resid + tshift);
                        //prob += 2.0 * exp(-misfit_tmp);
                        //misfit += 2.0 * misfit_tmp;
                    } else {
                        misfit_tmp = (double) *(wtmtxrow + ncol) *
                                arr_row_res * (arrival[ncol].cent_resid + tshift);
                        prob += exp(-0.5 * misfit_tmp);
                        misfit += misfit_tmp;
                    }
                }
            }
            prob /= num_arrivals;
            if (prob > prob_max) {
                //if (misfit < misfit_min) {
                misfit_min = misfit;
                //misfit_min = -log(prob / num_arrivals);
                prob_max = prob;
                tshift_at_prob_max = tshift;
            }
        }
        tstart = tshift_at_prob_max - tstep;
        tstop = tshift_at_prob_max + tstep;
        tstep /= 10.0;
    }
    misfit = misfit_min;
    //printf("tshift_at_prob_max %lf  prob_max %lf\n", tshift_at_prob_max, prob_max);

    // get origin time
    if (potime != NULL)
        *potime = CalcMaxLikeOriginTime(num_arrivals, arrival, gauss_par);

    /* return misfit or ln(prob density) */

    if (itype == GRID_MISFIT) {
        /* convert misfit to rms misfit */
        rms_misfit = sqrt(misfit / num_arrivals);
        *pmisfit = rms_misfit;
        return (rms_misfit);
    } else if (itype == GRID_PROB_DENSITY) {
        //ln_prob_density = -0.5 * misfit;
        ln_prob_density = log(prob_max) * num_arrivals * num_arrivals;
        //rms_misfit = sqrt(misfit / num_arrivals);
        rms_misfit = sqrt(misfit);
        *pmisfit = rms_misfit;
        return (ln_prob_density);
    } else {

        return (-1.0);
    }



}



/** function to calculate probability density */

/*		(MEN92, eq. 14) */

double CalcSolutionQuality_GAU_ANALYTIC(int num_arrivals, ArrivalDesc *arrival,
        GaussLocParams* gauss_par, int itype, double* pmisfit, double* potime) {

    int nrow, ncol, narr;
    int narr_misfit;

    double misfit = 0.0;
    double ln_prob_density, rms_misfit;

    double arr_row_res;
    MatrixDouble wtmtx;
    double *wtmtxrow;

    wtmtx = gauss_par->WtMtrx;


    /* calculate weighted mean of predicted travel times  */
    /*		(TV82, eq. A-38) */

    CalcCenteredTimesPred(num_arrivals, arrival, gauss_par);


    /* calculate residuals (TV82, eqs. 10-12, 10-13; MEN92, eq. 15) */

    for (narr = 0; narr < num_arrivals; narr++) {
        // AJL 20041115 bug fix!
        if (arrival[narr].pred_travel_time <= 0.0)
            arrival[narr].cent_resid = 0.0; // ignore obs without predicted times
        else
            arrival[narr].cent_resid = arrival[narr].obs_centered - arrival[narr].pred_centered;
    }

    narr_misfit = 0;
    for (nrow = 0; nrow < num_arrivals; nrow++) {
        // AJL 20041115 bug fix!
        if (arrival[nrow].pred_travel_time <= 0.0) {
            //printf("IGNORE: %s %s\n", arrival[nrow].label, arrival[nrow].phase);
            continue; // ignore obs without predicted times
        }
        // END
        if (!arrival[nrow].abs_time)
            continue;
        narr_misfit++;
        wtmtxrow = wtmtx[nrow];
        arr_row_res = arrival[nrow].cent_resid;
        for (ncol = 0; ncol <= nrow; ncol++) {
            // AJL 20041115 bug fix!
            if (arrival[ncol].pred_travel_time <= 0.0)
                continue; // ignore obs without predicted times
            // END
            if (!arrival[ncol].abs_time)
                continue;
            if (ncol != nrow)
                misfit += 2.0 * (double) *(wtmtxrow + ncol) *
                arr_row_res * arrival[ncol].cent_resid;
            else
                misfit += (double) *(wtmtxrow + ncol) *
                arr_row_res * arrival[ncol].cent_resid;
        }
    }


    if (potime != NULL)
        *potime = CalcMaxLikeOriginTime(num_arrivals, arrival, gauss_par);

    /* return misfit or ln(prob density) */

    if (itype == GRID_MISFIT) {
        /* convert misfit to rms misfit */
        if (narr_misfit > 0) // 20120724 AJL - bug fix!
            rms_misfit = sqrt(misfit / narr_misfit);
        else
            rms_misfit = VERY_LARGE_DOUBLE;
        *pmisfit = rms_misfit;
        return (rms_misfit);
    } else if (itype == GRID_PROB_DENSITY) {
        if (narr_misfit > 0) { // 20120724 AJL - bug fix!
            ln_prob_density = -0.5 * misfit;
            rms_misfit = sqrt(misfit / narr_misfit);
        } else {
            ln_prob_density = -VERY_LARGE_DOUBLE;
            rms_misfit = VERY_LARGE_DOUBLE;
        }
        *pmisfit = rms_misfit;
        return (ln_prob_density);
    } else {

        return (-1.0);
    }



}



/** function to calculate maximum likelihood estimate for origin time */

/*		(MEN92, eq. 19) */

long double CalcMaxLikeOriginTime(int num_arrivals, ArrivalDesc *arrival, GaussLocParams * gauss_par) {

    //MatrixDouble wtmtx;

    //wtmtx = gauss_par->WtMtrx;

    /* NOTE: the following (MEN92, eq. 19) is unnecessary

                                                                            for (nrow = 0; nrow < num_arrivals; nrow++)
                                                                            for (ncol = 0; ncol <= nrow; ncol++) {
                                                                            if (ncol != nrow)
                                                                            time_est += 2.0L * (long double) wtmtx[nrow][ncol]
     * (arrival[nrow].obs_time -
                                                                            (long double) arrival[nrow].pred_travel_time);
                                                                            else
                                                                            time_est += (long double) wtmtx[nrow][ncol] *
                                                                            (arrival[nrow].obs_time -
                                                                            (long double) arrival[nrow].pred_travel_time);
                                                                    }
                                                                            return(time_est / (long double) gauss_par->WtMtrxSum);
     */

    return (gauss_par->meanObs - (long double) gauss_par->meanPred);

}

/** function to update arrival probabilistic residuals */

void UpdateProbabilisticResiduals(int num_arrivals, ArrivalDesc *arrival, double prob) {

    int narr;


    /* update probabilistc residuals */

    for (narr = 0; narr < num_arrivals; narr++) {

        arrival[narr].pdf_residual_sum += prob *
                arrival[narr].cent_resid;
        arrival[narr].pdf_weight_sum += prob;
    }

}





/** function to calculate confidence intervals and save to file */
/*		(MEN92, eq. 25ff) */

#define N_STEPS_SRCH 101
#define N_STEPS_CONF 11

int CalcConfidenceIntrvl(GridDesc* ptgrid, HypoDesc* phypo, char* filename) {
    FILE *fpio;
    char fname[FILENAME_MAX];

    int ix, iy, iz;
    int iconf, isrch;
    double srch_level, srch_incr, conf_level, conf_incr, prob_den;
    double srch_sum[N_STEPS_SRCH];
    double contour[N_STEPS_CONF];
    double sum_volume;


    /* write message */
    if (message_flag >= 3) {
        nll_putmsg(3, "");
        nll_putmsg(3, "Calculating confidence intervals over grid...");
    }


    for (isrch = 0; isrch < N_STEPS_SRCH; isrch++)
        srch_sum[isrch] = 0.0;


    /* accumulate approx integral of probability density in search bins */
    /*  and normalize sum of bin values * dx*dy*dz */

    sum_volume = ptgrid->sum * ptgrid->dx * ptgrid->dy * ptgrid->dz;
    phypo->probmax /= sum_volume;
    srch_incr = phypo->probmax / (N_STEPS_SRCH - 1);
    for (ix = 0; ix < ptgrid->numx; ix++) {
        for (iy = 0; iy < ptgrid->numy; iy++) {
            for (iz = 0; iz < ptgrid->numz; iz++) {
                ((GRID_FLOAT_TYPE ***) ptgrid->array)[ix][iy][iz] = (GRID_FLOAT_TYPE) (
                        exp((double) ((GRID_FLOAT_TYPE ***) ptgrid->array)[ix][iy][iz])
                        / sum_volume);
                prob_den = ((GRID_FLOAT_TYPE ***) ptgrid->array)[ix][iy][iz];
                srch_level = 0.0;
                for (isrch = 0; isrch < N_STEPS_SRCH; isrch++) {
                    if (prob_den >= srch_level)
                        srch_sum[isrch] += prob_den;
                    srch_level += srch_incr;
                }

            }
        }
    }
    ptgrid->sum = 1.0;

    /* normalize by 100% confidence level sum */

    for (isrch = 1; isrch < N_STEPS_SRCH; isrch++)
        srch_sum[isrch] /= srch_sum[0];
    srch_sum[0] = 1.0;


    /* open confidence interval file */

    snprintf(fname, sizeof(fname), "%s.loc.conf", filename);
    if ((fpio = fopen(fname, "w")) == NULL) {
        nll_puterr("ERROR: opening confidence interval output file.");
        return (-1);
    } else {
        NumFilesOpen++;
    }

    /* find confidence levels and write to file */

    conf_incr = 1.0 / (N_STEPS_CONF - 1);
    conf_level = 1.0;
    iconf = N_STEPS_CONF - 1;
    for (isrch = 0; isrch < N_STEPS_SRCH; isrch++) {
        if (srch_sum[isrch] <= conf_level) {
            contour[iconf] = (double) isrch * srch_incr;
            fprintf(fpio, "%lf C %.2lf\n",
                    contour[iconf], conf_level);
            if (--iconf < 0)
                break;
            conf_level -= conf_incr;
        }
    }


    fclose(fpio);
    NumFilesOpen--;

    return (0);

}
