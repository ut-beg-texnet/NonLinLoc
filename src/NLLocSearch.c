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


/*   NLLocSearch.c

        Grid-search and Metropolis (simulated-annealing) search methods for
        NonLinLoc: LocGridSearch, LocMetropolis, GetNextMetropolisSample,
        MetropolisTest and SaveBestLocation.

        Extracted verbatim from NLLocLib.c during the modernization effort
        (Phase 3: modularization). Pure code movement, no logic change; the
        public prototypes remain in NLLocLib.h, so callers are unaffected.
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


/** function to perform grid search location */

int LocGridSearch(int ngrid, int num_arr_total, int num_arr_loc,
        ArrivalDesc *arrival,
        GridDesc* ptgrid, GaussLocParams* gauss_par, HypoDesc * phypo) {

    int istat;
    int ix = -1, iy = -1, iz = -1, narr;
    int iGridType;
    int nReject, numGridReject = 0, numStaReject = 0;
    double xval, yval, zval;
    /*double travel_time;*/
    double value;
    double misfit;
    double misfit_min = VERY_LARGE_DOUBLE, misfit_max = -VERY_LARGE_DOUBLE;
    double dlike;



    /* get solution quality at each grid point */

    if (message_flag >= 4) {
        nll_putmsg(4, "");
        nll_putmsg(4, "Calculating solution over grid...");
    }

    iGridType = ptgrid->type;

    xval = ptgrid->origx;

    /* loop over grid points */

    for (ix = 0; ix < ptgrid->numx; ix++) {

        /* read y-z sheets for arrival travel-times (3D grids) */
        if ((istat = ReadArrivalSheets(num_arr_loc, arrival, xval)) < 0)
            nll_puterr("ERROR: reading arrival travel time sheets.");

        yval = ptgrid->origy;
        for (iy = 0; iy < ptgrid->numy; iy++) {
            zval = ptgrid->origz;
            for (iz = 0; iz < ptgrid->numz; iz++) {


                // get travel times for observed arrivals

                if (isAboveTopo(xval, yval, zval)) {

                    misfit = -1.0;
                    value = 0.0;
                    if (iGridType == GRID_MISFIT)
                        value = -1.0;
                    else if (iGridType == GRID_PROB_DENSITY)
                        value = -LARGE_FLOAT;
                    ((GRID_FLOAT_TYPE ***) ptgrid->array)[ix][iy][iz] = value;

                } else {

                    nReject = getTravelTimes(arrival, num_arr_loc, xval, yval, zval);

                    if (nReject) {

                        numGridReject++;
                        numStaReject += nReject;
                        misfit = -1.0;
                        value = 0.0;
                        if (iGridType == GRID_MISFIT)
                            value = -1.0;
                        else if (iGridType == GRID_PROB_DENSITY)
                            value = -LARGE_FLOAT;
                        ((GRID_FLOAT_TYPE ***) ptgrid->array)[ix][iy][iz] = value;

                    } else {

                        /* calc misfit or prob density */

                        double log_prior;
                        value = CalcSolutionQuality(xval, yval, zval, NULL, num_arr_loc,
                                arrival, gauss_par,
                                iGridType, &misfit, NULL, NULL, 0.0, 0.0, 0.0, NULL, NULL, &log_prior);
                        if (iGridType == GRID_MISFIT) {
                            ptgrid->sum += value;
                        } else if (iGridType == GRID_PROB_DENSITY) {
                            value += log_prior; // 20190513 AJL
                            dlike = exp(value);
                            ptgrid->sum += dlike;
                            /* update  probabilistic residuals */
                            UpdateProbabilisticResiduals(num_arr_loc, arrival, dlike);
                        }
                        ((GRID_FLOAT_TYPE ***) ptgrid->array)[ix][iy][iz] = value;

                        /* check for minimum misfit */
                        if (misfit < misfit_min) {
                            misfit_min = misfit;
                            phypo->misfit = misfit;
                            phypo->ix = ix;
                            phypo->iy = iy;
                            phypo->iz = iz;
                            phypo->x = xval;
                            phypo->y = yval;
                            phypo->z = zval;
                            for (narr = 0; narr < num_arr_loc; narr++)
                                arrival[narr].pred_travel_time_best =
                                    arrival[narr].pred_travel_time;
                        }
                        if (misfit > misfit_max)
                            misfit_max = misfit;

                    }
                }

                zval += ptgrid->dz;
            }
            yval += ptgrid->dy;
        }
        xval += ptgrid->dx;
    }


    /* give warning if grid points rejected */

    if (numGridReject > 0) {
        snprintf(MsgStr, sizeof(MsgStr), "WARNING: %d grid locations rejected; travel times for an average of %.2lf arrival observations were not valid.",
                numGridReject, (double) numStaReject / numGridReject);
        nll_putmsg(1, MsgStr);
    }


    /* construct search information string */
    snprintf(phypo->searchInfo, sizeof(phypo->searchInfo), "GRID nPts %d%c", ix * iy *iz, '\0');
    /* write message */
    /*nll_putmsg(2, phypo->searchInfo);*/


    /* re-calculate solution and arrival statistics for best location */

    double cell_diagonal_time_var_best = 0.0; // TODO: add to Grid Search ?
    double cell_diagonal_best = 0.0; // TODO: add to Grid Search ?
    double cell_volume_best = 0.0; // TODO: add to Grid Search ?
    SaveBestLocation(NULL, num_arr_total, num_arr_loc, arrival, ptgrid, gauss_par, phypo, misfit_max,
            iGridType, 0, cell_diagonal_time_var_best, cell_diagonal_best, cell_volume_best);

    return (0);

}



/** function to perform Metropolis location */

#define MAX_NUM_MET_TRIES 1000

int LocMetropolis(int ngrid, int num_arr_total, int num_arr_loc,
        ArrivalDesc *arrival,
        GridDesc* ptgrid, GaussLocParams* gauss_par, HypoDesc* phypo,
        WalkParams* pMetrop, float* fdata) {

    int istat;
    int ntry, nSamples, nSampStat, narr, ipos;
    long int ngenerated;
    int maxNumTries;
    int writeMessage = 0;
    int iGridType;
    int nReject, numClipped = 0, numGridReject = 0, numStaReject = 0;
    int iAbort = 0, iReject = 0;
    int iBoundary = 0;
    int iAccept, numAcceptDeepMinima = 0;
    double xval, yval, zval;
    double currentMetStepFact;

    double value, dlike, dlike_max = -VERY_LARGE_DOUBLE;

    double misfit;
    double misfit_min = VERY_LARGE_DOUBLE, misfit_max = -VERY_LARGE_DOUBLE;

    double xmin, xmax, ymin, ymax, zmin, zmax;
    double dx_init, dx_test;

    int nScatterSaved;

    double xmean_sum = 0.0, ymean_sum = 0.0, zmean_sum = 0.0;
    double xvar_sum = 0.0, yvar_sum = 0.0, zvar_sum = 0.0;
    double xvar = 0.0, yvar = 0.0, zvar = 0.0;
    double dsamp = 0.0, dsamp2;



    /* get solution quality at each sample on random walk */

    if (message_flag >= 4) {
        nll_putmsg(4, "");
        nll_putmsg(4, "Calculating solution along Metropolis walk...");
    }

    iGridType = GRID_PROB_DENSITY;

    /* set walk limits equal to grid limits */
    xmin = ptgrid->origx;
    xmax = xmin + (double) (ptgrid->numx - 1) * ptgrid->dx;
    ymin = ptgrid->origy;
    ymax = ymin + (double) (ptgrid->numy - 1) * ptgrid->dy;
    zmin = ptgrid->origz;
    zmax = zmin + (double) (ptgrid->numz - 1) * ptgrid->dz;

    /* save intiial values */
    currentMetStepFact = MetStepFact;
    dx_init = pMetrop->dx;


    /* loop over walk samples */

    nSamples = 0;
    nSampStat = 0;
    nScatterSaved = 0;
    ipos = 0;
    ntry = 0;
    ngenerated = 0;
    maxNumTries = MAX_NUM_MET_TRIES;
    while (nSamples < MetNumSamples
            && (nSamples <= MetLearn || ntry < maxNumTries)) {

        ntry++;
        ngenerated++;
        istat = GetNextMetropolisSample(pMetrop,
                xmin, xmax, ymin, ymax,
                zmin, zmax, &xval, &yval, &zval);
        if (nSamples > MetEquil && istat > 0)
            numClipped += istat;

        /* get travel times for observed arrivals */

        if (isAboveTopo(xval, yval, zval)) {

            misfit = -1.0;
            dlike = 0.0;

        } else {

            nReject = getTravelTimes(arrival, num_arr_loc, xval, yval, zval);

            if (nReject) {
                numGridReject++;
                numStaReject += nReject;
                misfit = -1.0;
                dlike = 0.0;
            } else {

                /* calc misfit or prob density */
                double log_prior;
                value = CalcSolutionQuality(xval, yval, zval, NULL, num_arr_loc, arrival, gauss_par,
                        iGridType, &misfit, NULL, NULL, 0.0, 0.0, 0.0, NULL, NULL, &log_prior);
                value += log_prior; // 20190513 AJL
                dlike = gauss_par->WtMtrxSum * exp(value);

                /* apply Metropolis test */
                iAccept = MetropolisTest(pMetrop->likelihood, dlike);

                /* if not accepted, but at maxNumTries... */
                if (!iAccept && ntry == maxNumTries) {
                    /* if not learning, accept anyway since
                    may be stuck in a deep minima */
                    if (nSamples >= MetLearn && numAcceptDeepMinima++ < 5) {
                        iAccept = 1;
                        //printf("Max Num Tries: accept deep minima\n");

                        /* try reducing step size */
                        currentMetStepFact /= 2.0;
                        ntry = 0;
                        //printf("            +: step ch: was %lf\n", pMetrop->dx);

                        /*if (pMetrop->dx > MetStepMin) {
                        pMetrop->dx /= 2.0;
                        //printf("            +: step ch: %lf -> %lf\n", 2.0 * pMetrop->dx, pMetrop->dx);
                        ntry = 0;
                                                }*/

                        /* if learning, try reducing step size */
                    } else if (nSamples < MetLearn && pMetrop->dx > MetStepMin) {
                        pMetrop->dx /= 2.0;
                        //printf("Max Num Tries: step ch: %lf -> %lf\n", 2.0 * pMetrop->dx, pMetrop->dx);
                        ntry = 0;
                    }
                }


                if (iAccept) {

                    ntry = 0;
                    nSamples++;

                    /* check for minimum misfit */
                    if (misfit < misfit_min) {
                        misfit_min = misfit;
                        dlike_max = dlike;
                        phypo->misfit = misfit;
                        phypo->x = xval;
                        phypo->y = yval;
                        phypo->z = zval;
                        for (narr = 0; narr < num_arr_loc; narr++)
                            arrival[narr].pred_travel_time_best =
                                arrival[narr].pred_travel_time;
                    }
                    if (misfit > misfit_max)
                        misfit_max = misfit;

                    /* update sample location */
                    pMetrop->x = xval;
                    pMetrop->y = yval;
                    pMetrop->z = zval;
                    pMetrop->likelihood = dlike;

                    /* if learning, update sample statistics */
                    if (nSamples > MetLearn / 2 && nSamples <= MetLearn + MetEquil) {

                        xmean_sum += xval;
                        ymean_sum += yval;
                        zmean_sum += zval;
                        xvar_sum += xval * xval;
                        yvar_sum += yval * yval;
                        zvar_sum += zval * zval;
                        nSampStat++;
                    }

                    /* if equilibrating, update Met step */
                    if (nSamples > MetLearn
                            && nSamples <= MetLearn + MetEquil) {

                        /* update Met step */
                        dsamp = (double) nSampStat;
                        dsamp2 = dsamp * dsamp;
                        xvar = xvar_sum / dsamp -
                                xmean_sum * xmean_sum / dsamp2;
                        yvar = yvar_sum / dsamp -
                                ymean_sum * ymean_sum / dsamp2;
                        zvar = zvar_sum / dsamp -
                                zmean_sum * zmean_sum / dsamp2;
                        dx_test = currentMetStepFact * pow(
                                sqrt(xvar) * sqrt(yvar) * sqrt(zvar)
                                / (double) MetUse, 1.0 / 3.0);
                        /*/ (double) (MetUse / MetSkip),*/
                        //if (pMetrop->dx != dx_test) printf("equil step ch: %lf -> %lf\n", pMetrop->dx, dx_test);

                        if (dx_test > MetStepMin)
                            pMetrop->dx = dx_test;
                        else
                            pMetrop->dx = MetStepMin;
                    }

                    /* if saving samples */
                    if (nSamples > MetStartSave
                            && nSamples % MetSkip == 0) {

                        /* save sample to scatter file */
                        fdata[ipos++] = xval;
                        fdata[ipos++] = yval;
                        fdata[ipos++] = zval;
                        fdata[ipos++] = dlike;

                        /* update  probabilitic residuals */
                        if (1)
                            UpdateProbabilisticResiduals(
                                num_arr_loc, arrival, 1.0);


                        nScatterSaved++;
                    }

                    if (nSamples % 1000 == 1
                            || nSamples == MetLearn / 2)
                        writeMessage = 1;

                }


                if (writeMessage || ntry == maxNumTries - 1) {
                    if (message_flag >= 4) {
                        snprintf(MsgStr, sizeof(MsgStr),
                                "Metropolis: n %d x %.2lf y %.2lf z %.2lf  xm %.2lf ym %.2lf zm %.2lf  xdv %.2lf ydv %.2lf zdv %.2lf  dx %.2lf  li %.2le", nSamples, pMetrop->x, pMetrop->y, pMetrop->z, xmean_sum / dsamp, ymean_sum / dsamp, zmean_sum / dsamp, sqrt(xvar), sqrt(yvar), sqrt(zvar), pMetrop->dx, pMetrop->likelihood);
                        nll_putmsg(4, MsgStr);
                    }
                    writeMessage = 0;
                }

            }
        }


        /* check abort search conditions */

        /* failure to accept sample after maxNumTries */
        if (nSamples > MetLearn && ntry >= maxNumTries) {
            snprintf(MsgStr, sizeof(MsgStr),
                    "ERROR: failed to accept new Metropolis sample after %d tries, aborting location.", ntry);
            nll_puterr(MsgStr);
            snprintf(phypo->locStatComm, sizeof (phypo->locStatComm), "%s", MsgStr);
            iAbort = 1;
            break;
        }

        /* maximum likelihood too low after learning stage */
        if (nSamples == MetLearn && dlike_max < MetProbMin) {
            snprintf(MsgStr, sizeof(MsgStr),
                    "ERROR: after learning stage (%d samples), best probability = %.2le is less than ProbMin = %.2le, aborting location.",
                    MetLearn, dlike_max, MetProbMin);
            nll_puterr(MsgStr);
            snprintf(phypo->locStatComm, sizeof (phypo->locStatComm), "%s", MsgStr);
            iAbort = 1;
            break;
        }

    }


    /* give warning if sample points clipped */

    if (numClipped > 0) {
        snprintf(MsgStr, sizeof(MsgStr), "WARNING: %d Metropolis samples clipped at search grid boundary.",
                numClipped);
        nll_putmsg(1, MsgStr);
    }


    /* give warning if grid points rejected */

    if (numGridReject > 0) {
        snprintf(MsgStr, sizeof(MsgStr), "WARNING: %d Metropolis samples rejected; travel times for an average of %.2lf arrival observations were not valid.",
                numGridReject, (double) numStaReject / numGridReject);
        nll_putmsg(1, MsgStr);
    }


    /* check reject location conditions */

    /* maximum like hypo on edge of grid */
    if ((iBoundary = isOnGridBoundary(phypo->x, phypo->y, phypo->z,
            ptgrid, pMetrop->dx, pMetrop->dx, 0))) {
        snprintf(MsgStr, sizeof(MsgStr), "WARNING: max prob location on grid boundary %d, rejecting location.", iBoundary);
        nll_putmsg(1, MsgStr);
        snprintf(phypo->locStatComm, sizeof (phypo->locStatComm), "%s", MsgStr);
        iReject = 1;
    }

    /* construct search information string */
    snprintf(phypo->searchInfo, sizeof(phypo->searchInfo),
            "METROPOLIS nSamp %ld nAcc %d nSave %d nClip %d Dstep0 %lf Dstep %lf%c",
            ngenerated, nSamples, nScatterSaved, numClipped, dx_init, pMetrop->dx, '\0');
    /* write message */
    nll_putmsg(2, phypo->searchInfo);


    /* check for termination */
    if (iAbort) {
        snprintf(Hypocenter.locStat, sizeof(Hypocenter.locStat), "ABORTED");
    } else if (iReject) {
        snprintf(Hypocenter.locStat, sizeof(Hypocenter.locStat), "REJECTED");
    }


    /* re-calculate solution and arrival statistics for best location */

    double cell_diagonal_time_var_best = 0.0; // TODO: add to Metropolis Search ?
    double cell_diagonal_best = 0.0; // TODO: add to Metropolis Search ?
    double cell_volume_best = 0.0; // TODO: add to Metropolis Search ?
    SaveBestLocation(NULL, num_arr_total, num_arr_loc, arrival, ptgrid,
            gauss_par, phypo, misfit_max, iGridType, 0, cell_diagonal_time_var_best, cell_diagonal_best, cell_volume_best);

    return (nScatterSaved);

}




/** function to create next metropolis sample */

/* move sample random distance and direction */

int GetNextMetropolisSample(WalkParams* pMetrop, double xmin, double xmax,
        double ymin, double ymax, double zmin, double zmax,
        double* pxval, double* pyval, double* pzval) {

    int iClip = 0;
    double valx, valy, valz, valsum, norm;
    double x, y, z;


    /* get unit vector in random direction */

    do {
        valx = get_rand_double(-1.0, 1.0);
        valy = get_rand_double(-1.0, 1.0);
        valz = get_rand_double(-1.0, 1.0);
        valsum = valx * valx + valy * valy + valz * valz;
    } while (valsum < SMALL_DOUBLE);

    norm = pMetrop->dx / sqrt(valsum);

    /* add step to last sample location */

    x = pMetrop->x + norm * valx;
    y = pMetrop->y + norm * valy;
    z = pMetrop->z + norm * valz;


    /* crude clip against grid boundary */
    /* clip needed because travel time lookup requires that
    location is within initial search grid */
    if (x < xmin) {
        x = xmin;
        iClip = 1;
    } else if (x > xmax) {
        x = xmax;
        iClip = 1;
    }
    if (y < ymin) {
        y = ymin;
        iClip = 1;
    } else if (y > ymax) {
        y = ymax;
        iClip = 1;
    }
    if (z < zmin) {
        z = zmin;
        iClip = 1;
    } else if (z > zmax) {
        z = zmax;
        iClip = 1;
    }


    /* update sample location */

    *pxval = x;
    *pyval = y;
    *pzval = z;

    return (iClip);

}

/** function to test new metropolis string */

int MetropolisTest(double likelihood_last, double likelihood_new) {

    double prob;

    /* compare with last sample using Mosegaard & Tarantola eq (17) */

    if (likelihood_new >= likelihood_last)
        return (1);
    else if ((prob = get_rand_double(0.0, 1.0)) < likelihood_new / likelihood_last)
        return (1);

    else
        return (0);

}


/** function to re-calculate solution and arrival statistics for best location */

/* some quantities are calculated only for arrivals used in location
                (num_arr_loc) others for all arrivals (num_arr_total) */

int SaveBestLocation(OctNode* poct_node, int num_arr_total, int num_arr_loc, ArrivalDesc *arrival,
        GridDesc* ptgrid, GaussLocParams* gauss_par, HypoDesc* phypo,
        double misfit_max, int iGridType, int ignore_pred_travel_time_best,
        double cell_diagonal_time_var_best, double cell_diagonal_best, double cell_volume_best) {

    int istat, narr, n_compan;
    char filename[2 * FILENAME_MAX];

    SourceDesc station;

    //printf("SaveBestLocation num_arr_total %d num_arr_loc %d\n", num_arr_total, num_arr_loc);

    // force longitude within -180->180 for global mode  // 20160922 AJL - added
    if (GeometryMode == MODE_GLOBAL) {
        if (phypo->x < -180.0)
            phypo->x += 360.0;
        else if (phypo->x > 180.0)
            phypo->x -= 360.0;
    }


    // 20101005 AJL - added calculation of mean slowness
    double slowness_P = -1.0;
    double slowness_S = -1.0;
    if (LocMethod == METH_OT_STACK) {
        double yval_grid;
        if (fp_model_grid_P != NULL) {
            if (model_grid_P.numx > 2) {
                // 3D grid
                slowness_P = (double) ReadAbsInterpGrid3d(fp_model_grid_P, &model_grid_P, phypo->x, phypo->y, phypo->z, 1);
            } else {
                // 2D grid (1D model)
                yval_grid = model_grid_P.dy; // aribitrary, small y grid value
                slowness_P = ReadAbsInterpGrid2d(fp_model_grid_P, &model_grid_P, yval_grid, phypo->z);
                if (GeometryMode != MODE_GLOBAL)
                    slowness_P /= model_grid_P.dy; // value in model file is slowness * ds
            }
        }
        if (fp_model_grid_S != NULL) {
            if (model_grid_S.numx > 2) {
                // 3D grid
                slowness_S = (double) ReadAbsInterpGrid3d(fp_model_grid_S, &model_grid_S, phypo->x, phypo->y, phypo->z, 1);
            } else {
                // 2D grid (1D model)
                yval_grid = model_grid_S.dy; // aribitrary, small y grid value
                slowness_S = ReadAbsInterpGrid2d(fp_model_grid_S, &model_grid_S, yval_grid, phypo->z);
                if (GeometryMode != MODE_GLOBAL)
                    slowness_S /= model_grid_P.dy; // value in model file is slowness * ds
            }
        }
        if (slowness_P <= SMALL_FLOAT)
            slowness_P = -1.0;
        if (slowness_S < 0.0 && VpVsRatio > 0.0)
            slowness_S = slowness_P * VpVsRatio;
        if (slowness_S <= SMALL_FLOAT)
            slowness_S = -1.0;
    }

    /* loop over observed arrivals */
    for (narr = 0; narr < num_arr_total; narr++) {

        arrival[narr].dist = GetEpiDist(&(arrival[narr].station), phypo->x, phypo->y);
        // 20060619 AJL - dist changed to always km, output converted to degrees for GLOBAL in
        //if (GeometryMode == MODE_GLOBAL)
        //	arrival[narr].dist *= KM2DEG;
        arrival[narr].azim = GetEpiAzim(&(arrival[narr].station), phypo->x, phypo->y);

        /* get best travel time */

        arrival[narr].pred_travel_time = 0.0;
        n_compan = arrival[narr].n_companion;
        //iopened = 0;
        /* check for stored best travel time */
        if (!ignore_pred_travel_time_best && arrival[narr].pred_travel_time_best > 0.0) {
            /* load stored best travel time */
            arrival[narr].pred_travel_time = arrival[narr].pred_travel_time_best;
            /* check for companion travel time */
        } else if (!ignore_pred_travel_time_best && n_compan >= 0 && arrival[n_compan].pred_travel_time_best > 0.0) {
            /* load companion stored best travel time */
            arrival[narr].pred_travel_time =
                    arrival[n_compan].pred_travel_time_best;
            arrival[narr].pred_travel_time *= arrival[narr].tfact;
        } else {
            // temporarily open time grid file and read time for ignored arrivals
            // save station information (will be overwritten in OpenGrid3dFile()
            station = arrival[narr].station;
            snprintf(filename, sizeof(filename), "%s.time", arrival[narr].fileroot);
            // 20250215 AJL - Bug fix, check if need to open companion time grid file
            if (n_compan >= 0)
                snprintf(filename, sizeof (filename), "%s.time", arrival[n_compan].fileroot);
            else
                snprintf(filename, sizeof (filename), "%s.time", arrival[narr].fileroot);
            //  20260106 AJL - Bug Fix: check if time grid already open; prevent orphan open files!
            if (arrival[narr].fpgrid == NULL) {
                if ((istat = OpenGrid3dFile(filename,
                        &(arrival[narr].fpgrid),
                        &(arrival[narr].fphdr),
                        &(arrival[narr].gdesc), "time",
                        &(arrival[narr].station),
                        iSwapBytesOnInput)) < 0)
                    continue;
            }
            arrival[narr].station = station;
            //iopened = 1;
            /* check grid type, read travel time */
            if (arrival[narr].gdesc.type == GRID_TIME) {
                /* 3D grid */
                if (arrival[narr].fpgrid != NULL)
                    arrival[narr].pred_travel_time = (double) ReadAbsInterpGrid3d(
                        arrival[narr].fpgrid, &(arrival[narr].gdesc),
                        phypo->x, phypo->y, phypo->z, 1);
                if (arrival[narr].pred_travel_time < -LARGE_DOUBLE)
                    arrival[narr].pred_travel_time = 0.0;
            } else {
                /* 2D grid (1D model) */
                // AJL 20060602 dist stored as km, KM2DEG added
                if (arrival[narr].fpgrid != NULL)
                    arrival[narr].pred_travel_time =
                        ReadAbsInterpGrid2d(
                        arrival[narr].fpgrid, &(arrival[narr].gdesc),
                        GeometryMode == MODE_GLOBAL ? arrival[narr].dist * KM2DEG : arrival[narr].dist,
                        phypo->z);
                if (arrival[narr].pred_travel_time < -LARGE_DOUBLE)
                    arrival[narr].pred_travel_time = 0.0;
            }
            arrival[narr].pred_travel_time *= arrival[narr].tfact;
            // apply crustal correction
            if (ApplyCrustElevCorrFlag && GeometryMode == MODE_GLOBAL
                    && arrival[narr].pred_travel_time > 0.0) {
                if (arrival[narr].dist > MinDistCrustElevCorr)
                    arrival[narr].pred_travel_time +=
                        applyCrustElevCorrection(arrival + narr, phypo->x, phypo->y, phypo->z);
            } else if (ApplyElevCorrFlag) {
                if (arrival[narr].pred_travel_time > 0.0) // ignore arrivals with no pred tt
                    arrival[narr].pred_travel_time += arrival[narr].elev_corr;
            }
            CloseGrid3dFile(&(arrival[narr].gdesc), &(arrival[narr].fpgrid), &(arrival[narr].fphdr));

        }

        /* read angles */
        /* angle grid file name */
        if (n_compan >= 0)
            snprintf(filename, sizeof (filename), "%s.angle", arrival[n_compan].fileroot);
        else
            snprintf(filename, sizeof (filename), "%s.angle", arrival[narr].fileroot);
        if (angleMode == ANGLE_MODE_YES) {
            if (arrival[narr].gdesc.type == GRID_TIME) {
                /* 3D grid */
                ReadTakeOffAnglesFile(filename,
                        phypo->x, phypo->y, phypo->z,
                        &(arrival[narr].ray_azim),
                        &(arrival[narr].ray_dip),
                        &(arrival[narr].ray_qual), -1.0, iSwapBytesOnInput);
            } else {
                /* 2D grid (1D model) */
                // AJL 20060828 dist stored as km, KM2DEG added
                ReadTakeOffAnglesFile(filename,
                        0.0,
                        GeometryMode == MODE_GLOBAL ? arrival[narr].dist * KM2DEG : arrival[narr].dist,
                        phypo->z,
                        &(arrival[narr].ray_azim),
                        &(arrival[narr].ray_dip),
                        &(arrival[narr].ray_qual), arrival[narr].azim, iSwapBytesOnInput);
            }
        }

        //		/* close time grid file for ignored arrivals */
        //		if (iopened)
        //			CloseGrid3dFile(&(arrival[narr].fpgrid),
        //				&(arrival[narr].fphdr));

        // set slowness
        if (arrival[narr].isS)
            arrival[narr].slowness = slowness_S;
        else
            arrival[narr].slowness = slowness_P;

    }

    /* calc misfit or prob density */
    double value, misfit, otime, otime_var, effective_cell_size, ot_variance_factor;
    otime_var = -1.0;
    double log_prior;
    value = CalcSolutionQuality(phypo->x, phypo->y, phypo->z, poct_node, num_arr_loc, arrival, gauss_par, iGridType, &misfit, &otime, &otime_var,
            cell_diagonal_time_var_best, cell_diagonal_best, cell_volume_best, &effective_cell_size, &ot_variance_factor, &log_prior);
    value += log_prior; // 20190513 AJL

    // set rms if otime variance is available
    if (otime_var > 0.0
            && !FixOriginTimeFlag // 20201201 AJL - bug fix.
            )
        phypo->rms = sqrt(otime_var);
    else
        phypo->rms = -1.0;

    /* set origin time */
    if (!FixOriginTimeFlag)
        phypo->time = otime;
    if (iGridType == GRID_PROB_DENSITY) {
        phypo->probmax = expl(value); // 20130314 C Satriano, AJL - changed to long double
    }

    /* set misc hypo fields */
    phypo->grid_misfit_max = misfit_max;
    istat = rect2latlon(0, phypo->x, phypo->y, &(phypo->dlat), &(phypo->dlong));
    phypo->depth = phypo->z;
    phypo->nreadings = num_arr_loc;

    return (0);

}
