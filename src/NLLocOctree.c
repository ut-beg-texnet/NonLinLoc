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


/*   NLLocOctree.c

        Oct-tree search routines for NonLinLoc: InitializeOcttree, LocOctree,
        LocOctree_core, the getOctTreeStationDensityWeight family and
        GenEventScatterOcttree.

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


/*------------------------------------------------------------/ */
/** Octtree search routines */

/** function to initialize Octtree search */

Tree3D * InitializeOcttree(GridDesc* ptgrid, OcttreeParams * pParams) {

    double dx, dy, dz;
    Tree3D* newTree;
    void *pdata = NULL;
    double integral;

    // set up oct-tree x, y, z grid
    dx = ptgrid->dx * (double) (ptgrid->numx - 1) / (double) (pParams->init_num_cells_x);
    dy = ptgrid->dy * (double) (ptgrid->numy - 1) / (double) (pParams->init_num_cells_y);
    //dz = ptgrid->dz * (double) (ptgrid->numz - 1) / (double) (pParams->init_num_cells_z);
    // 20101004 AJL - Fixed bug: Oct-tree search was not initialized to bottom of LOCSEARCH grid, bottom LOCSEARCH grid layer was lost.
    //dz = ptgrid->dz * (double) ptgrid->numz / (double) (pParams->init_num_cells_z);
    // 20200302 AJL - Un-Fixed bug: Seems that Oct-tree search was initialized to bottom of LOCSEARCH grid.
    dz = ptgrid->dz * (double) (ptgrid->numz - 1) / (double) (pParams->init_num_cells_z);

    integral = 0.0;
    // 20160927 AJL   if (LocMethod == METH_OT_STACK && GeometryMode == MODE_GLOBAL) {
    //if (GeometryMode == MODE_GLOBAL && !iSaveNLLocOctree) { // 20160927 AJL - bug fix?  octtree.c->writeTree3D() (called if iSaveNLLocOctree) does not yet support spherical trees
    if (GeometryMode == MODE_GLOBAL) { // 20190508 AJL - bug fix?  octtree.c->writeTree3D() now handles spherical trees issues
        newTree = newTree3D_spherical(ptgrid->type, pParams->init_num_cells_x,
                pParams->init_num_cells_y, pParams->init_num_cells_z,
                ptgrid->origx, ptgrid->origy, ptgrid->origz,
                dx, dy, dz, OCTREE_UNDEF_VALUE, integral, pdata);
    } else {

        newTree = newTree3D(ptgrid->type, pParams->init_num_cells_x,
                pParams->init_num_cells_y, pParams->init_num_cells_z,
                ptgrid->origx, ptgrid->origy, ptgrid->origz,
                dx, dy, dz, OCTREE_UNDEF_VALUE, integral, pdata);
    }

    return (newTree);
}

/** function to perform Octree location */

int LocOctree(int ngrid, int num_arr_total, int num_arr_loc,
        ArrivalDesc *arrival,
        GridDesc* ptgrid, GaussLocParams* gauss_par, HypoDesc* phypo,
        OcttreeParams* pParams, Tree3D* pOctTree, float* fdata,
        double *poct_node_value_max, double *poct_tree_integral) {

    int nSamples, narr, ipos;
    int nInitial;
    int iGridType;
    //int nReject;
    int iReject = 0;
    int iBoundary = 0;
    double xval, yval, zval;

    long double value, dlike, value_max = (long double) -VERY_LARGE_DOUBLE;
    double misfit = -1.0;
    //double misfit_min = VERY_LARGE_DOUBLE;
    double misfit_max = -VERY_LARGE_DOUBLE;
    double hypo_dx = -1.0, hypo_dz = -1.0;
    double cell_diagonal_time_var_best = 0.0;
    double cell_diagonal_best = 0.0;
    double cell_volume_best = 0.0;
    OctNode* poct_node_best = NULL;

    int nScatterSaved;

    int ix, iy, iz;
    double logWtMtrxSum;
    //double volume, log_value_volume;
    int icalc_cell_diagonal_time_var = 0;
    double diagonal, cell_half_diagonal_time_range, volume_min;
    //double dsx, dsy, dsz;
    //double dsx_global, dsy_global, depth_corr;
    ResultTreeNode* presult_node;
    OctNode* poct_node;
    OctNode* pparent_oct_node;

    int n_neigh;
    OctNode* neighbor_node;
    Vect3D coords;
    double min_node_size_x;
    double min_node_size_y;
    double min_node_size_z;
    double smallest_node_size_x = -1.0;
    double smallest_node_size_y = -1.0;
    double smallest_node_size_z = -1.0;

    //double stationDensityWeight = 0.0;



    // reset EDT_otime_weight_active flag
    EDT_otime_weight_active = 0;

    // cell diagonal variance defaults to null - i.e. no effect
    diagonal = cell_half_diagonal_time_range = 0.0;
    icalc_cell_diagonal_time_var = pParams->mean_cell_velocity > 0.0;
    volume_min = VERY_LARGE_DOUBLE;


    iGridType = GRID_PROB_DENSITY;

    if (message_flag >= 4) {
        nll_putmsg(4, "");
        nll_putmsg(4, "Calculating solution in Octree...");
    }

    if (LocMethod == METH_OT_STACK) {
        logWtMtrxSum = 0.0;
    } else {
        logWtMtrxSum = log(gauss_par->WtMtrxSum);
    }

    // set min node size (in km)
    min_node_size_x = min_node_size_y = min_node_size_z = pParams->min_node_size;
    // following neglects convergence of longitude towards the poles and convergence with depth
    if (GeometryMode == MODE_GLOBAL)
        min_node_size_x = min_node_size_y = pParams->min_node_size * KM2DEG;

    /* first get solutions at each cell in Tree3D */

    nSamples = 0;
    resultTreeRoot = NULL;
    for (ix = 0; ix < pOctTree->numx; ix++) {
        for (iy = 0; iy < pOctTree->numy; iy++) {
            for (iz = 0; iz < pOctTree->numz; iz++) {
                poct_node = pOctTree->nodeArray[ix][iy][iz];
                if (poct_node == NULL) // case of Tree3D_spherical
                    continue;
                // $$$ NOTE: this block must be identical to block $$$ below
                xval = poct_node->center.x;
                yval = poct_node->center.y;
                zval = poct_node->center.z;
                value = LocOctree_core(ngrid, xval, yval, zval, num_arr_loc, arrival, poct_node,
                        icalc_cell_diagonal_time_var, &volume_min, &diagonal,
                        &cell_half_diagonal_time_range, pParams, gauss_par, iGridType, &misfit, logWtMtrxSum);
                nSamples++;
                // END - this block must be identical to block $$$ below

                // save node size
                smallest_node_size_x = poct_node->ds.x;
                smallest_node_size_y = poct_node->ds.y;
                smallest_node_size_z = poct_node->ds.z;

                if (message_flag >= 1 && nSamples % 5000 == 0) {
                    fprintf(stdout,
                            "OctTree num samples = %d / %d\r", nSamples, pParams->max_num_nodes);
                    fflush(stdout);
                }

            }
        }
    }
    nInitial = nSamples;


    /* loop over oct-tree nodes */

    nScatterSaved = 0;
    ipos = 0;

    while (nSamples < pParams->max_num_nodes) {

        if (pParams->stop_on_min_node_size)
            presult_node = getHighestLeafValue(resultTreeRoot);
        else
            presult_node = getHighestLeafValueMinSize(resultTreeRoot,
                min_node_size_x, min_node_size_y, min_node_size_z);
        // check if null node
        if (presult_node == NULL) {
            if (message_flag >= 1)
                fprintf(stdout, "\nINFO: No more nodes larger than min_node_size, terminating Octree search.");
            break;
        }


        //if (nSamples % 100 == 0)
        //fprintf(stderr, "%d getHighestLeafValue %lf\n", nSamples, presult_node->value);

        if (presult_node == NULL)
            fprintf(stderr, "\npresult_node == NULL!!\n");

        pparent_oct_node = presult_node->pnode;

        // subdivide all HighestLeafValue neighbors

        int n_neigh_max = 7;
        if (LocMethod == METH_OT_STACK) // this is in warning monitor for speed and efficiency in convergence, with the risk of less thorough search
            n_neigh_max = 1;
        for (n_neigh = 0; n_neigh < n_neigh_max; n_neigh++) {

            if (n_neigh == 0) {
                neighbor_node = pparent_oct_node;
            } else {
                coords.x = pparent_oct_node->center.x;
                coords.y = pparent_oct_node->center.y;
                coords.z = pparent_oct_node->center.z;
                if (n_neigh == 1) {
                    coords.x = pparent_oct_node->center.x
                            + (pparent_oct_node->ds.x + smallest_node_size_x) / 2.0;
                } else if (n_neigh == 2) {
                    coords.x = pparent_oct_node->center.x
                            - (pparent_oct_node->ds.x + smallest_node_size_x) / 2.0;
                } else if (n_neigh == 3) {
                    coords.y = pparent_oct_node->center.y
                            + (pparent_oct_node->ds.y + smallest_node_size_y) / 2.0;
                } else if (n_neigh == 4) {
                    coords.y = pparent_oct_node->center.y
                            - (pparent_oct_node->ds.y + smallest_node_size_y) / 2.0;
                } else if (n_neigh == 5) {
                    coords.z = pparent_oct_node->center.z
                            + (pparent_oct_node->ds.z + smallest_node_size_z) / 2.0;
                } else if (n_neigh == 6) {
                    coords.z = pparent_oct_node->center.z
                            - (pparent_oct_node->ds.z + smallest_node_size_z) / 2.0;
                }
                // check for longitude wrap-around
                if (GeometryMode == MODE_GLOBAL) {
                    if (coords.x > 180.0)
                        coords.x -= 360.0;
                    if (coords.x < -180.0)
                        coords.x += 360.0;
                }

                // find neighbor node
                neighbor_node = getLeafNodeContaining(pOctTree, coords);
                // outside of octTree volume
                if (neighbor_node == NULL)
                    continue;
                // already subdivided
                if (neighbor_node->ds.z < 0.99 * pparent_oct_node->ds.z)
                    continue;
            }


            // subdivide node and evaluate solution at each child
            subdivide(neighbor_node, OCTREE_UNDEF_VALUE, NULL);

            for (ix = 0; ix < 2; ix++) {
                for (iy = 0; iy < 2; iy++) {
                    for (iz = 0; iz < 2; iz++) {

                        poct_node = neighbor_node->child[ix][iy][iz];

                        //if (poct_node->ds.x < pParams->min_node_size || poct_node->ds.y < pParams->min_node_size || poct_node->ds.z < pParams->min_node_size)
                        //fprintf(stderr, "\nnode size too small!! %lf %lf %lf\n", poct_node->ds.x, poct_node->ds.y, poct_node->ds.z);

                        // save node size if smallest so far
                        if (poct_node->ds.x < smallest_node_size_x)
                            smallest_node_size_x = poct_node->ds.x;
                        if (poct_node->ds.y < smallest_node_size_y)
                            smallest_node_size_y = poct_node->ds.y;
                        if (poct_node->ds.z < smallest_node_size_z)
                            smallest_node_size_z = poct_node->ds.z;

                        // $$$ NOTE: this block must be identical to block $$$ above
                        xval = poct_node->center.x;
                        yval = poct_node->center.y;
                        zval = poct_node->center.z;
                        value = LocOctree_core(ngrid, xval, yval, zval, num_arr_loc, arrival, poct_node,
                                icalc_cell_diagonal_time_var, &volume_min, &diagonal,
                                &cell_half_diagonal_time_range, pParams, gauss_par, iGridType, &misfit, logWtMtrxSum);
                        nSamples++;
                        // END - this block must be identical to block $$$ above

                        if (message_flag >= 1 && nSamples % 5000 == 0) {
                            fprintf(stdout,
                                    "OctTree num samples = %d / %d\r", nSamples, pParams->max_num_nodes);
                            fflush(stdout);
                        }

                        // check value
                        /*if (value < -LARGE_FLOAT) {
                            snprintf(MsgStr, sizeof(MsgStr), "ERROR: log(prob_density) at (%lf,%lf,%lf) is too small %lg.", xval, yval, zval, (double) value);
                            nll_puterr(MsgStr);
                        }*/
                        /*if (isnan(value)) {
                            snprintf(MsgStr, sizeof(MsgStr), "WARNNG: log(prob_density) at (%lf,%lf,%lf) is NaN (%lg), reset to %g.", xval, yval, zval, (double) value, -VERY_LARGE_DOUBLE);
                            nll_puterr(MsgStr);
                            value = -VERY_LARGE_DOUBLE;
                        }*/

                        /* check for maximum likelihood */
                        //printf("value=%lg, value_max=%lg, diagonal=%f\r", (double) value, (double) value_max, diagonal);
                        if (value >= value_max) {
                            //printf(">>>>>>>>>>>>>>>>> value=%lg > value_max=%lg!!, diagonal=%f, xyz= %f %g %g\n", (double) value, (double) value_max, diagonal, xval, yval, zval);
                            value_max = value;
                            //misfit_min = misfit;
                            phypo->misfit = misfit;
                            phypo->x = xval;
                            phypo->y = yval;
                            phypo->z = zval;
                            hypo_dx = poct_node->ds.x;
                            hypo_dz = poct_node->ds.z;
                            for (narr = 0; narr < num_arr_loc; narr++)
                                arrival[narr].pred_travel_time_best = arrival[narr].pred_travel_time;
                            poct_node_best = poct_node;
                            *poct_node_value_max = poct_node->value;
                            cell_diagonal_time_var_best = cell_half_diagonal_time_range * cell_half_diagonal_time_range;
                            cell_diagonal_best = diagonal;
                            cell_volume_best = volume_min;
                        }
                        if (misfit > 0.0 && misfit > misfit_max) // misfit < 0 for topo masking
                            misfit_max = misfit;


                        /* set to TRUE to save all samples, REMEMBER to set OCT num_scatter high enough in control file */
                        if (0) {
                            /* save sample to scatter file */
                            fdata[ipos++] = xval;
                            fdata[ipos++] = yval;
                            fdata[ipos++] = zval;
                            dlike = (long double) gauss_par->WtMtrxSum * (long double) exp(value);
                            fdata[ipos++] = dlike;

                            /* update  probabilitic residuals */
                            if (1)
                                UpdateProbabilisticResiduals(num_arr_loc, arrival, 1.0);

                            nScatterSaved++;
                        }


                    } // end triple loop over node children
                }
            }

        } // end loop over HighestLeafValue neighbors

        // check if minimum node size reached
        if (pParams->stop_on_min_node_size && (smallest_node_size_x < min_node_size_x
                || smallest_node_size_y < min_node_size_y
                || smallest_node_size_z < min_node_size_z)) {
            if (message_flag >= 1)
                fprintf(stdout, "\nINFO: Min node size reached, terminating Octree search.");
            break;
        }

    } // end while (nSamples < pParams->max_num_nodes)

    if (message_flag >= 1)
        fprintf(stdout, "\n");




    // for OT_STACK best location must be leaf node
    /*
    if (LocMethod == METH_OT_STACK) {
        ResultTreeNode* presultTreeNode = getHighestLeafValue(resultTreeRoot);
        poct_node = presultTreeNode->pnode;
        value_max = poct_node->value;
        //misfit_min = misfit;
        phypo->misfit = -1.0;
        phypo->x = poct_node->center.x;
        phypo->y = poct_node->center.y;
        phypo->z = poct_node->center.z;
        hypo_dx = poct_node->ds.x;
        hypo_dz = poct_node->ds.z;
        for (narr = 0; narr < num_arr_loc; narr++)
            arrival[narr].pred_travel_time_best = -1.0;
     *poct_node_value_max = poct_node->value;
        //cell_diagonal_time_var_best = cell_half_diagonal_time_range * cell_half_diagonal_time_range;
        //cell_diagonal_best = diagonal;
        cell_volume_best = presultTreeNode->value;
    }
     */



    /* check reject location conditions */

    /* maximum like hypo on edge of grid */

    if ((iBoundary = isOnGridBoundary(phypo->x, phypo->y, phypo->z, ptgrid, hypo_dx, hypo_dz, 0))) {
        snprintf(MsgStr, sizeof(MsgStr),
                "WARNING: max prob location on grid boundary %d, rejecting location.", iBoundary);
        nll_putmsg(1, MsgStr);
        snprintf(phypo->locStatComm, sizeof (phypo->locStatComm), "%s", MsgStr);
        iReject = 1;
    }

    // determine integral of all oct-tree leaf node pdf values
    *poct_tree_integral = integrateResultTree(resultTreeRoot, VALUE_IS_LOG_PROB_DENSITY_IN_NODE, 0.0, *poct_node_value_max);
    snprintf(MsgStr, sizeof(MsgStr), "Octree oct_node_value_max= %le oct_tree_integral= %le", *poct_node_value_max, *poct_tree_integral);
    nll_putmsg(1, MsgStr);



    // construct search information string
    if (GeometryMode == MODE_GLOBAL) {
        smallest_node_size_x *= DEG2KM;
        smallest_node_size_y *= DEG2KM;
    }
    snprintf(phypo->searchInfo, sizeof(phypo->searchInfo), "OCTREE nInitial %d nEvaluated %d smallestNodeSide %lf/%lf/%lf oct_tree_integral %le",
            nInitial, nSamples, smallest_node_size_x, smallest_node_size_y, smallest_node_size_z, *poct_tree_integral);
    // set values in hypo
    phypo->oct_tree_integral = *poct_tree_integral;
    /* write message */
    nll_putmsg(2, phypo->searchInfo);


    /* check for termination */
    if (iReject) {
        snprintf(Hypocenter.locStat, sizeof(Hypocenter.locStat), "REJECTED");
    }


    /* re-calculate solution and arrival statistics for best location */

    //int XX_last = NumAllocations;
    SaveBestLocation(poct_node_best, num_arr_total, num_arr_loc, arrival, ptgrid,
            gauss_par, phypo, misfit_max, iGridType, 0, cell_diagonal_time_var_best, cell_diagonal_best, cell_volume_best);
    //printf("XXX: SaveBestLocation: NumAllocations %d->%d\n", XX_last, NumAllocations);

    return (nScatterSaved);

}

/** function to perform Octree core solution evaluation */

long double LocOctree_core(int ngrid, double xval, double yval, double zval,
        int num_arr_loc, ArrivalDesc *arrival,
        OctNode* poct_node,
        int icalc_cell_diagonal_time_var, double *volume_min,
        double *pdiagonal, double *cell_half_diagonal_time_range,
        OcttreeParams* pParams, GaussLocParams* gauss_par, int iGridType,
        double *misfit, double logWtMtrxSum) {

    long double value;

    int iAboveTopo;
    int nReject;
    double volume, log_value_volume;
    double dsx, dsy, dsz;
    double dsx_global, dsy_global, depth_corr;


    /* get travel times for observed arrivals */
    iAboveTopo = isAboveTopo(xval, yval, zval);
    if (!iAboveTopo) {
        nReject = getTravelTimes(arrival, num_arr_loc, xval, yval, zval);
        if (message_flag > 3 && nReject && GeometryMode != MODE_GLOBAL) {
            snprintf(MsgStr, sizeof(MsgStr),
                    "WARNING: oct-tree sample at (%lf,%lf,%lf) is outside of %d travel time grids.",
                    xval, yval, zval, nReject);
            nll_putmsg(4, MsgStr);
        }
    }

    // calculate cell volume and diagonal
    dsx = poct_node->ds.x;
    dsy = poct_node->ds.y;
    dsz = poct_node->ds.z;
    if (GeometryMode == MODE_GLOBAL) {
        //depth_corr = (ERAD - poct_node->center.z) / ERAD;   // NOTE: ERAD is max Earth radius, could use average radius AVG_ERAD, difference should be very minor
        depth_corr = (AVG_ERAD - poct_node->center.z) / AVG_ERAD; // 20151106 AJL - changed to AVG_ERAD, difference should be very minor
        dsx_global = dsx * DEG2KM * cos(DE2RA * poct_node->center.y) * depth_corr;
        dsy_global = dsy * DEG2KM * depth_corr;
        volume = dsx_global * dsy_global * dsz;
        if (icalc_cell_diagonal_time_var || LocMethod == METH_OT_STACK) {
            //if (volume < *volume_min)
            *volume_min = volume;
            *pdiagonal = pow(volume, 0.33333333);
            //*diagonal = sqrt(dsx_global * dsx_global + dsy_global * dsy_global + dsz * dsz);
            //*diagonal = dsx_global < dsy_global ? dsx_global : dsy_global;
            //*diagonal = dsz < *diagonal ? dsz : *diagonal;
            //*diagonal = dsx_global > dsy_global ? dsx_global : dsy_global;
            //*diagonal = dsz > *diagonal ? dsz : *diagonal;
        }
    } else {
        volume = dsx * dsy * dsz;
        if (icalc_cell_diagonal_time_var || LocMethod == METH_OT_STACK) {
            //if (volume < *volume_min)
            *volume_min = volume;
            *pdiagonal = pow(volume, 0.33333333);
            //*diagonal = sqrt(dsx * dsx + dsy * dsy + dsz * dsz);
        }
    }
    if (icalc_cell_diagonal_time_var) {
        // 20101005 AJL
        *cell_half_diagonal_time_range = 0.5 * *pdiagonal / pParams->mean_cell_velocity;
        //printf("*cell_diagonal_time_var %lf\n", *cell_diagonal_time_var);
    }
    // calc misfit and prob density
    //EDT_use_otime_weight = 0;
    double effective_cell_size = -1.0;
    double ot_variance_factor = 0.0;
    if (!iAboveTopo) { // not above topo
        double log_prior;
        value = CalcSolutionQuality(xval, yval, zval, poct_node, num_arr_loc, arrival, gauss_par, iGridType, misfit, NULL, NULL,
                *cell_half_diagonal_time_range, *pdiagonal, volume, &effective_cell_size, &ot_variance_factor, &log_prior);
        /*        if (LocMethod == METH_OT_STACK) {
                    if (poct_node->parent != NULL) {
                        //value += (poct_node->parent->value - logWtMtrxSum);
                        //value -= log(8); // volume parent / volume
                    }
                }
         */
        value += log_prior; // 20190513 AJL
    } else {
        value = -VERY_LARGE_DOUBLE;
        *misfit = -VERY_LARGE_DOUBLE;
    }
    double logStationDensityWeight = 0.0;
    if (pParams->use_stations_density > 0) {
        logStationDensityWeight =
                getOctTreeStationDensityWeight(poct_node, StationPhaseList, NumStationPhases, LocGrid + ngrid, pParams->use_stations_density);
        // 20101022 AJL - limit small values of log st wt, needed to stabilize OT_STACK method
        if (logStationDensityWeight < -10.0)
            logStationDensityWeight = -10.0;
    }
    // calculate cell volume * prob density and put cell in resultTree
    log_value_volume = logWtMtrxSum + value + log(volume);
    poct_node->value = logWtMtrxSum + value;
    //EDT_use_otime_weight = 1;
    // AJL 20060531 bug fix ? - changed *= to +=
    log_value_volume += logStationDensityWeight;
    poct_node->value += logStationDensityWeight;

    resultTreeRoot = addResult(resultTreeRoot, log_value_volume, volume, poct_node);

    /*static int icount_value = 0;
                                                                                                                                                                    if (icount_value < 10 && poct_node->value < -1.0e50) {
                                                                                                                                                                        printf("poct_node->value < -1.0e50 !!! poct_node->value %lg  logStationDensityWeight %lg  logWtMtrxSum %lg  log(volume) %lg\n",
                                                                                                                                                                                poct_node->value, logStationDensityWeight, logWtMtrxSum, log(volume));
                                                                                                                                                                        icount_value++;
                                                                                                                                                                    }*/

    return (value);

}

/** function to calculate (logarithmic) station density weight value for an oct tree node */

double getOctTreeStationDensityWeight_OLD1(OctNode* poct_node, SourceDesc *stations, int numStations, GridDesc * ptgrid) {

    int n, n_inside, numStations_this_event;
    double staDensityWeight;
    SourceDesc *station;

    // check if parent node contains no stations

    if (poct_node->parent != NULL
            && poct_node->parent->pdata != NULL && *((int *) poct_node->parent->pdata) <= 1)
        return (1.0);

    // count number of stations inside this node
    numStations_this_event = 0;
    n_inside = 0;
    for (n = 0; n < numStations; n++) {
        station = stations + n;
        // check if station has not-ignored reading for this event
        if (station->ignored)
            continue;
        numStations_this_event++;
        // check if station location is unknown
        if (station->x <= -LARGE_DOUBLE)
            continue;
        // check if station is above top of grid
        if (station->z < ptgrid->origz) {
            if (extendedNodeContains(poct_node, station->x, station->y, ptgrid->origz, 0))
                n_inside++;
        } else {
            if (extendedNodeContains(poct_node, station->x, station->y, station->z, 0))
                n_inside++;
        }
        //printf("n %d  node %f %f %f  stations %f %f %f  ptgrid->origz %f\n", n, poct_node->center.x, poct_node->center.y, poct_node->center.z, station->x, station->y, station->z, ptgrid->origz);
    }

    if (poct_node->pdata == NULL)
        poct_node->pdata = (void *) malloc(sizeof (int));
    if (poct_node->pdata != NULL)
        *((int *) poct_node->pdata) = n_inside;
    else
        nll_puterr("ERROR: allocating int storage for OctTree Station Density Weight count.");

    staDensityWeight = log((double) (n_inside + 1));

    //if (n_inside > 0) printf("OctTreeStationDensityWeight: node %f %f %f  n_inside %d  numStations_this_event %d\n", poct_node->center.x, poct_node->center.y, poct_node->center.z, n_inside, numStations_this_event);

    return (staDensityWeight);

}

/** function to calculate (logarithmic) station density weight value for an oct tree node */

double getOctTreeStationDensityWeight_OLD2(OctNode* poct_node, SourceDesc *stations, int numStations, GridDesc * ptgrid) {

    int n, n_inside, numStations_this_event;
    double staDensityWeight;
    SourceDesc *station;


    n_inside = 0;

    // check if parent node contains stations data
    if (poct_node->parent != NULL) {
        // this must be a child node, parent should have station count data
        if (poct_node->parent->pdata == NULL) { // should not get here
            nll_puterr("ERROR: parent node exists but has no OctTree Station Density Weight count!");
        } else {
            n_inside = *((int *) poct_node->parent->pdata);
        }
    } else {
        // this must be a root node
        // count number of stations inside this node
        numStations_this_event = 0;
        n_inside = 0;
        for (n = 0; n < numStations; n++) {
            station = stations + n;
            // check if station has not-ignored reading for this event
            if (station->ignored)
                continue;
            numStations_this_event++;
            // check if station location is unknown
            if (station->x <= -LARGE_DOUBLE)
                continue;
            // check if station is above top of grid
            if (station->z < ptgrid->origz) {
                if (extendedNodeContains(poct_node, station->x, station->y, ptgrid->origz, 0)) {
                    n_inside++;
                    //printf("n %d  node %f %f %f  stations %f %f %f  ptgrid->origz %f\n", n, poct_node->center.x, poct_node->center.y, poct_node->center.z, station->x, station->y, station->z, ptgrid->origz);
                }
            } else {
                if (extendedNodeContains(poct_node, station->x, station->y, station->z, 0)) {
                    n_inside++;
                    //printf("n %d  node %f %f %f  stations %f %f %f  ptgrid->origz %f\n", n, poct_node->center.x, poct_node->center.y, poct_node->center.z, station->x, station->y, station->z, ptgrid->origz);
                }
            }
        }
        //if (n_inside > 0) printf("OctTreeStationDensityWeight: node %f %f %f  n_inside %d  numStations_this_event %d\n", poct_node->center.x, poct_node->center.y, poct_node->center.z, n_inside, numStations_this_event);
    }

    if (poct_node->pdata == NULL)
        poct_node->pdata = (void *) malloc(sizeof (int));
    if (poct_node->pdata != NULL)
        *((int *) poct_node->pdata) = n_inside;
    else
        nll_puterr("ERROR: allocating int storage for OctTree Station Density Weight count.");

    staDensityWeight = log((double) (n_inside + 1));

    // multiply by arbitrary constant... (seems to work!)
    //staDensityWeight *= 10.0;
    staDensityWeight *= 2.0;

    return (staDensityWeight);

}

/** function to calculate (logarithmic) station density weight value for an oct tree node */

double getOctTreeStationDensityWeight(OctNode* poct_node, SourceDesc *stations, int numStations, GridDesc *ptgrid, int iOctLevelMax) {

    int node_level, n, numStations_this_event;
    double log_station_density_weight, return_weight;
    SourceDesc *station;
    double epi_dist, depth_diff, hypo_dist, hypo_dist_min, cut_off_dist;
    double x_node_cent, y_node_cent, z_node_cent, mean_node_horiz_ds;
    OctNode* pnode;

    static double mean_root_node_horiz_ds = -VERY_LARGE_DOUBLE;
    // !!! shoud be initialized for each event????



    if (mean_root_node_horiz_ds == -VERY_LARGE_DOUBLE) {
        pnode = poct_node;
        while (pnode->parent != NULL)
            pnode = pnode->parent;
        mean_root_node_horiz_ds = (pnode->ds.x + pnode->ds.y);
        if (GeometryMode == MODE_GLOBAL)
            mean_root_node_horiz_ds *= DEG2KM;
        snprintf(MsgStr, sizeof(MsgStr), "Station Density Weight:  Mean Root Node Horiz dS: %lf", mean_root_node_horiz_ds);
        nll_putmsg(1, MsgStr);
    }
    if (mean_root_node_horiz_ds < SMALL_DOUBLE) { // should not get here
        nll_puterr("ERROR: cannot apply OctTree Station Density Weight: Mean Root Node Horiz dS is zero!");
    }



    log_station_density_weight = 0.0;
    return_weight = log_station_density_weight;

    // get level of this node in tree
    node_level = 0;
    pnode = poct_node;
    while (pnode->parent != NULL) {
        pnode = pnode->parent;
        node_level++;
    }


    // node level dependent calculation
    if (node_level >= iOctLevelMax) {
        // node above max level, get data from parent
        if (poct_node->parent->pdata == NULL) { // should not get here
            nll_puterr("ERROR: parent node exists but has no OctTree Station Density Weight value!");
        } else {
            log_station_density_weight = *((double *) poct_node->parent->pdata);
        }
        return_weight = log_station_density_weight;
        //TESTreturn_weight = 0.0;
    } else {
        // node below max level
        // calcualte average station distance from center of this node
        x_node_cent = poct_node->center.x;
        y_node_cent = poct_node->center.y;
        z_node_cent = poct_node->center.z;
        numStations_this_event = 0;
        hypo_dist_min = VERY_LARGE_DOUBLE;
        for (n = 0; n < numStations; n++) {
            station = stations + n;
            // check if station has ignored reading for this event
            if (station->ignored)
                continue;
            numStations_this_event++;
            // check if station location is unknown
            if (station->x <= -LARGE_DOUBLE)
                continue;
            // get distance from station to center of node
            epi_dist = GetEpiDist(station, x_node_cent, y_node_cent);
            depth_diff = z_node_cent - station->z;
            hypo_dist = sqrt(epi_dist * epi_dist + depth_diff * depth_diff);
            hypo_dist_min = hypo_dist < hypo_dist_min ? hypo_dist : hypo_dist_min;
            //printf("n %d  node %f %f %f  stations %f %f %f  ptgrid->origz %f\n", n, poct_node->center.x, poct_node->center.y, poct_node->center.z, station->x, station->y, station->z, ptgrid->origz);
        }
        //if (ndist == 0) {	// should not get here
        //	nll_puterr("ERROR: no stations found for OctTree Station Density Weight calculation!");
        //}
        if (hypo_dist_min > VERY_SMALL_DOUBLE) {
            mean_node_horiz_ds = (poct_node->ds.x + poct_node->ds.y);
            if (GeometryMode == MODE_GLOBAL)
                mean_node_horiz_ds *= DEG2KM;
            cut_off_dist = mean_node_horiz_ds > AveInterStationDistance
                    ? mean_node_horiz_ds : AveInterStationDistance;
            /*
            if (hypo_dist_min <= cut_off_dist)
            log_station_density_weight = 1.0;
            else
            log_station_density_weight = cut_off_dist / hypo_dist_min;
            //log_station_density_weight *= mean_root_node_horiz_ds / AveInterStationDistance;
            log_station_density_weight *= 10.0;
            //	if (log_station_density_weight > 10.0)
            //		log_station_density_weight = 10.0;
             */
            log_station_density_weight = hypo_dist_min / cut_off_dist;
            log_station_density_weight = -log_station_density_weight * log_station_density_weight;
            //log_station_density_weight *= mean_root_node_horiz_ds / AveInterStationDistance;
            //log_station_density_weight *= 10.0;
            //	if (log_station_density_weight > 10.0)
            //		log_station_density_weight = 10.0;
            //if (log_station_density_weight >= -0.5) printf("OctTreeStationDensityWeight: node %f %f %f  sta_den_wt %f  numStations_this_event %d  AveInterStationDistance %f  hypo_dist_min %f  mean_node_horiz_ds %f  node_level %d\n", poct_node->center.x, poct_node->center.y, poct_node->center.z, exp(log_station_density_weight), numStations_this_event, AveInterStationDistance, hypo_dist_min, mean_node_horiz_ds, node_level);
            // force cell division
            if (node_level < iOctLevelMax && hypo_dist_min < cut_off_dist) {
                return_weight = (double) (iOctLevelMax - node_level);
                return_weight = return_weight * return_weight;
                NumForceOctTreeStaDenWt++;
            } else {
                return_weight = log_station_density_weight;
            }
        }
    }

    if (poct_node->pdata == NULL)
        poct_node->pdata = (void *) malloc(sizeof (double));
    if (poct_node->pdata != NULL)
        *((double *) poct_node->pdata) = log_station_density_weight;
    else
        nll_puterr("ERROR: allocating int storage for OctTree Station Density Weight count.");

    // multiply by arbitrary constant... (seems to work!)
    //staDensityWeight *= 10.0;
    //log_station_density_weight *= 10.0;

    return (return_weight);

}

/** function to generate sample (scatter) of OctTree results */

int GenEventScatterOcttree(OcttreeParams* pParams, double oct_node_value_max, float* fscatterdata, double integral, HypoDesc * phypo) {

    int tot_npoints;
    int fdata_index;
    double oct_tree_scatter_volume;
    char scatter_volume_text[32];


    oct_tree_scatter_volume = 0.0;

    /* return if no scatter samples requested */
    if (pParams->num_scatter < 1)
        return (0);

    // return if integral is nan    // 20201022 AJL - bug fix
    if (isnan(integral)) {
        nll_puterr("ERROR: Generating event scatter: oct_tree_integral is nan.");
        return (0);
    }

    /* write message */
    if (message_flag >= 3) {
        nll_putmsg(3, "");
        nll_putmsg(3, "Generating event scatter file...");
    }


    /* generate scatter points at uniformly-randomly chosen locations in each leaf node */

    tot_npoints = 0;
    fdata_index = 0;
    tot_npoints = getScatterSampleResultTree(resultTreeRoot, VALUE_IS_LOG_PROB_DENSITY_IN_NODE, pParams->num_scatter, integral,
            fscatterdata, tot_npoints, &fdata_index, oct_node_value_max, &oct_tree_scatter_volume);

    /* write message */
    if (message_flag >= 3) {
        snprintf(MsgStr, sizeof(MsgStr), "  %d points generated, %d points requested, oct_tree_scatter_volume= %le",
                tot_npoints, pParams->num_scatter, oct_tree_scatter_volume);
        nll_putmsg(3, MsgStr);
    }

    // update hypocenter searchInfo
    snprintf(scatter_volume_text, sizeof(scatter_volume_text), " scatter_volume %le", oct_tree_scatter_volume);
    strncat(phypo->searchInfo, scatter_volume_text, sizeof(phypo->searchInfo) - strlen(phypo->searchInfo) - 1);
    // set values in hypo
    phypo->oct_tree_scatter_volume = oct_tree_scatter_volume;

    return (tot_npoints);

}



/** end of Octree search routines */
/*------------------------------------------------------------/ */
