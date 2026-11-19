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


/*   NLLocStaStat.c

        Station-statistics hash-table routines for NonLinLoc.

        Extracted verbatim from NLLocLib.c during the modernization effort
        (Phase 3: modularization). Pure code movement, no logic change; the
        struct, HASHSIZE, extern hashtab and the function prototypes remain
        in NLLocLib.h, so callers are unaffected.
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


/* table of pointers to list of StaStatNode (declared extern in NLLocLib.h) */
StaStatNode *hashtab[MAX_NUM_LOCATION_GRIDS][HASHSIZE];

/*------------------------------------------------------------/ */
/** hashtable routines for accumulating station statistics */

/* from Kernigham and Ritchie, C prog lang, 2nd ed, 1988, sec 6.6 */

/** funtion to form ordered hash value from firt char of label */

static unsigned hash(char* label, char* phase) {

    unsigned hashval;

    if (isdigit(label[0]))
        hashval = label[0] - '0';
    else if (isalpha(label[0]))
        hashval = 10 + toupper(label[0]) - 'A';
    else
        hashval = 36 + label[0] % 10;

    hashval = hashval % HASHSIZE;

    return hashval;
}

/** function to lookup labelphase in hashtable */

static StaStatNode * lookup(int ntable, char* label, char* phase) {

    StaStatNode *np;

    for (np = hashtab[ntable][hash(label, phase)]; np != NULL; np = np->next)
        if (strcmp(label, np->label) == 0
                && strcmp(phase, np->phase) == 0)
            return (np); /* found */

    return (NULL); /* not found */

}

/** function to install or add (labelphase, residual, weight) in hashtable */

StaStatNode * InstallStaStatInTable(int ntable, char* label, char* phase, int flag_ignore,
        double residual, double weight,
        double pdf_residual_sum, double pdf_weight_sum, double delay) {
    int icomp;
    StaStatNode *np, *npcheck, *nplast;
    unsigned hashval;

    if ((np = lookup(ntable, label, phase)) == NULL) {
        /* not found, create new StaStatNode */
        if ((np = (StaStatNode *) malloc(sizeof (StaStatNode))) == NULL)
            return (NULL);
        snprintf(np->label, sizeof(np->label), "%s", label);
        snprintf(np->phase, sizeof(np->phase), "%s", phase);
        np->flag_ignore = flag_ignore;
        np->residual_min = residual;
        np->residual_max = residual;
        np->residual_sum = residual * weight;
        np->residual_square_sum = residual * residual * weight;
        np->weight_sum = weight;
        np->num_residuals = 1;
        np->next = NULL; // 20101124 racine@sed.ethz.ch - added: Otherwise the next pointer will be left uninitialized under some circumstances
        // which leads to a crash when FreeStaStatTable() is called and the pointer is non-zero.
        // After this fix, it seems that this solves the memory problems.
        if (pdf_weight_sum > VERY_SMALL_DOUBLE) {
            np->pdf_residual_sum =
                    pdf_residual_sum / pdf_weight_sum;
            np->pdf_residual_square_sum =
                    pdf_residual_sum * pdf_residual_sum / (pdf_weight_sum * pdf_weight_sum);
            np->num_pdf_residuals = 1;
        } else {
            np->num_pdf_residuals = 0;
        }
        np->delay = delay;

        /* put in table in alphabetical order */
        hashval = hash(label, phase);
        nplast = NULL;
        npcheck = hashtab[ntable][hashval];
        while (npcheck != NULL &&
                (icomp = strcmp(npcheck->label, label)) <= 0) {
            if (icomp == 0 && strcmp(npcheck->phase, phase) >= 0)
                break;
            nplast = npcheck;
            npcheck = npcheck->next;
        }
        np->next = npcheck;
        if (nplast != NULL)
            nplast->next = np;
        else
            hashtab[ntable][hashval] = np;
    } else {
        /* already there */
        if (residual < np->residual_min)
            np->residual_min = residual;
        if (residual > np->residual_max)
            np->residual_max = residual;
        np->residual_sum += residual * weight;
        np->residual_square_sum += residual * residual * weight;
        np->weight_sum += weight;
        np->num_residuals++;
        if (pdf_weight_sum > VERY_SMALL_DOUBLE) {

            np->pdf_residual_sum +=
                    pdf_residual_sum / pdf_weight_sum;
            np->pdf_residual_square_sum +=
                    pdf_residual_sum * pdf_residual_sum / (pdf_weight_sum * pdf_weight_sum);
            np->num_pdf_residuals++;
        }
    }

    return (np);

}

/** function to free hashtable */

int FreeStaStatTable(int ntable) {
    int nnodes;
    unsigned hashval;
    StaStatNode *np, *np_curr;


    nnodes = 0;
    for (hashval = 0; hashval < HASHSIZE; hashval++) {
        for (np = hashtab[ntable][hashval]; np != NULL;) {

            np_curr = np;
            np = np->next;
            free(np_curr);
            np_curr = NULL;
            nnodes++;
        }
        hashtab[ntable][hashval] = NULL; // 20101129 AJL Bug fix.
    }

    return (nnodes);

}

/** function to output hashtable values */

int WriteStaStatTable(int ntable, FILE *fpio,
        double rms_max, int nRdgs_min, double gap_max,
        double p_residual_max, double s_residual_max,
        double ell_len3_max, double hypo_depth_min, double hypo_depth_max,
        double hypo_dist_max, int imode) {
    int nnodes;
    unsigned hashval;
    char frmt1[MAXLINE], frmt2[MAXLINE];
    double res_temp, res_std_temp;
    StaStatNode *np;

    /* 20160919 AJL  snprintf(frmt1, sizeof(frmt1), "LOCDELAY  %%-%ds %%-%ds %%-8d %%-12lf %%-12lf\n",
            ARRIVAL_LABEL_LEN, ARRIVAL_LABEL_LEN);
    snprintf(frmt2, sizeof(frmt2), "LOCDELAY  %%-%ds %%-%ds %%-8d %%-12lf %%-12lf %%-12lf %%-12lf %%d\n",
            ARRIVAL_LABEL_LEN, ARRIVAL_LABEL_LEN);*/
    snprintf(frmt1, sizeof(frmt1), "LOCDELAY  %%-s %%-s %%-8d %%-12lf %%-12lf\n");
    snprintf(frmt2, sizeof(frmt2), "LOCDELAY  %%-s %%-s %%-8d %%-12lf %%-12lf %%-12lf %%-12lf %%d\n");

    if (imode == WRITE_RESIDUALS) {
        fprintf(fpio,
                "\n#Average Phase Residuals (CalcResidual)  RMS_Max: %lf  NRdgs_Min: %d  Gap_Max: %lf  P_Res_Max: %lf  S_Res_Max: %lf  Ell_Len3_Max: %lf  Hypo_Depth_Min: %lf  Hypo_Depth_Max: %lf  Hypo_Dist_Max: %lf\n",
                rms_max, nRdgs_min, gap_max,
                p_residual_max, s_residual_max, ell_len3_max,
                hypo_depth_min, hypo_depth_max, hypo_dist_max);
        fprintf(fpio,
                "#         ID      Phase   Nres      AveRes       StdDev       ResMin       ResMax     ignored\n");
    } else if (imode == WRITE_RES_DELAYS) {
        fprintf(fpio,
                "\n#Total Phase Corrections (CalcResidual + InputDelay)  RMS_Max: %lf  NRdgs_Min: %d  Gap_Max: %lf  P_Res_Max: %lf  S_Res_Max: %lf  Ell_Len3_Max: %lf  Hypo_Depth_Min: %lf  Hypo_Depth_Max: %lf  Hypo_Dist_Max: %lf\n",
                rms_max, nRdgs_min, gap_max,
                p_residual_max, s_residual_max, ell_len3_max,
                hypo_depth_min, hypo_depth_max, hypo_dist_max);
        fprintf(fpio,
                "#         ID      Phase   Nres      TotCorr      StdDev\n");
    } else if (imode == WRITE_PDF_RESIDUALS) {
        fprintf(fpio,
                "\n#Average Phase Residuals PDF (CalcPDFResidual)  RMS_Max: %lf  NRdgs_Min: %d  Gap_Max: %lf  P_Res_Max: %lf  S_Res_Max: %lf  Ell_Len3_Max: %lf  Hypo_Depth_Min: %lf  Hypo_Depth_Max: %lf  Hypo_Dist_Max: %lf\n",
                rms_max, nRdgs_min, gap_max,
                p_residual_max, s_residual_max, ell_len3_max,
                hypo_depth_min, hypo_depth_max, hypo_dist_max);
        fprintf(fpio,
                "#         ID      Phase   Nres      AveRes       StdDev       ResMin       ResMax     ignored\n");
    } else if (imode == WRITE_PDF_DELAYS) {
        fprintf(fpio,
                "\n#Total Phase Corrections PDF (CalcPDFResidual + InputDelay)  RMS_Max: %lf  NRdgs_Min: %d  Gap_Max: %lf  P_Res_Max: %lf  S_Res_Max: %lf  Ell_Len3_Max: %lf  Hypo_Depth_Min: %lf  Hypo_Depth_Max: %lf  Hypo_Dist_Max: %lf\n",
                rms_max, nRdgs_min, gap_max,
                p_residual_max, s_residual_max, ell_len3_max,
                hypo_depth_min, hypo_depth_max, hypo_dist_max);
        fprintf(fpio,
                "#         ID      Phase   Nres      TotCorr      StdDev\n");
    }

    nnodes = 0;
    for (hashval = 0; hashval < HASHSIZE; hashval++) {
        for (np = hashtab[ntable][hashval]; np != NULL; np = np->next) {
            if (imode == WRITE_RESIDUALS || imode == WRITE_RES_DELAYS) {
                res_temp = np->residual_sum / np->weight_sum;
                res_std_temp = np->residual_square_sum / np->weight_sum - res_temp * res_temp;
                if (np->num_residuals > 1)
                    res_std_temp = sqrt(np->residual_square_sum / np->weight_sum - res_temp * res_temp);
                else
                    res_std_temp = -1.0;
                if (imode == WRITE_RESIDUALS) {
                    fprintf(fpio, frmt2, np->label, np->phase,
                            np->num_residuals, res_temp, res_std_temp,
                            np->residual_min, np->residual_max, np->flag_ignore);
                } else if (imode == WRITE_RES_DELAYS) {
                    fprintf(fpio, frmt1, np->label, np->phase,
                            np->num_residuals, res_temp + np->delay, res_std_temp);
                }
                //printf("LOCDELAY  %s %s %d %f = %f + %f/%f\n", np->label, np->phase, np->num_residuals, res_temp,
                //np->delay, np->residual_sum, np->weight_sum);
            } else if (imode == WRITE_PDF_RESIDUALS || imode == WRITE_PDF_DELAYS) {
                if (np->num_pdf_residuals > 0) {
                    res_temp = np->pdf_residual_sum / (double) np->num_pdf_residuals;
                } else {
                    res_temp = 0.0;
                }
                if (np->num_pdf_residuals > 1)
                    res_std_temp = sqrt(np->pdf_residual_square_sum
                        / (double) (np->num_pdf_residuals - 1)
                        - res_temp * res_temp);
                else
                    res_std_temp = -1.0;
                if (imode == WRITE_PDF_RESIDUALS) {
                    fprintf(fpio, frmt2, np->label, np->phase,
                            np->num_pdf_residuals, res_temp, res_std_temp,
                            np->residual_min, np->residual_max, np->flag_ignore);
                } else if (imode == WRITE_PDF_DELAYS) {

                    fprintf(fpio, frmt1, np->label, np->phase,
                            np->num_pdf_residuals, res_temp + np->delay, res_std_temp);
                }
            }
            nnodes++;
        }
    }


    return (nnodes);

}

/** function to update hashtable */

void UpdateStaStat(int ntable, ArrivalDesc *arrival, int num_arrivals,
        double p_residual_max, double s_residual_max, double hypo_dist_max, double weight) {
    int narr;
    StaStatNode *np;

    // 20181005 AJL - moved to function argument
    //double weight = 1.0;

    for (narr = 0; narr < num_arrivals; narr++) {
        if (
                (
                (IsPhaseID((arrival + narr)->phase, "P")
                && fabs((arrival + narr)->residual) <= p_residual_max)
                ||
                (IsPhaseID((arrival + narr)->phase, "S")
                && fabs((arrival + narr)->residual) <= s_residual_max)
                )
                &&
                (arrival + narr)->dist <= hypo_dist_max
                )

            if ((np = InstallStaStatInTable(ntable,
                    (arrival + narr)->label,
                    (arrival + narr)->phase,
                    (arrival + narr)->flag_ignore,
                    (arrival + narr)->residual,
                    weight,
                    (arrival + narr)->pdf_residual_sum,
                    (arrival + narr)->pdf_weight_sum,
                    (arrival + narr)->delay)) == NULL)
                nll_puterr("ERROR: cannot put arrival statistics in table");
    }

}


/** end of hashtable routines */
/*------------------------------------------------------------/ */
