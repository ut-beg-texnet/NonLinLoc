/***************************************************************************
 * polarisation.c:
 *
 * Polarization analysis library.
 *
 * Copyright (C) 2016 Anthony Lomax <anthony@alomax.net www.alomax.net>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "alomax_matrix.h"
#include "polarization.h"

#ifndef RAD2DEG
#define RAD2DEG (180.0 / M_PI)
#endif

/**
 * Method to apply polarization analysis to single set of arrays using real Covariance analysis assumes arrays are in ZYX order
 *
 * follows Covariance algorithms of :
 *
 * Hendrick, N. and S. Hearn (1999), Polarisation Analysis: What is it? Why do you need it? How do you do it?, Explor. Geophys., 30, 177-190.
 *
 * Nguyen, D. T., R. J. Brown, D. C. Lawton, (1989), Polarization filter for multi-component seismic data, in CREWES Research Report 1989, chap.
 * 7, 93-101, (http://www.crewes.org/Reports/1989/1989-07.pdf).
 *
 *
 */
PolarizationAnalysis *polarization_covariance(double *data_z, double *data_x, double *data_y, int istart, int nsamp, int ipolarity_z, double y_azimuth, int verbose) {

    // check some things

    // initialize return value
    PolarizationAnalysis *polar_analysis = malloc(sizeof (PolarizationAnalysis));
    if (polar_analysis == NULL) {
        return (NULL);
    }
    polar_analysis->azimuth = -1.0;
    polar_analysis->dip = -999.0;
    polar_analysis->degree_of_linearity = -1.0;
    polar_analysis->degree_of_planarity = -1.0;
    polar_analysis->mean_vect_amp = -1.0;

    // declare matrix
    MatrixDouble A_matrix = matrix_double(3, 3);

    double mean_vect_amp = 0.0;
    int namp = 0;
    double polarity_z = ipolarity_z > 0.0 ? 1.0 : -1.0;
    for (int i = istart; i < istart + nsamp; i++) { // loop over samples

        // accumulate covariance matrix
        double u = data_x[i];
        double v = data_y[i];
        double w = data_z[i] * polarity_z;
        A_matrix[0][0] += u * u;
        A_matrix[0][1] += u * v;
        A_matrix[0][2] += u * w;
        A_matrix[1][0] += v * u;
        A_matrix[1][1] += v * v;
        A_matrix[1][2] += v * w;
        A_matrix[2][0] += w * u;
        A_matrix[2][1] += w * v;
        A_matrix[2][2] += w * w;
        // accumulate mean_vect_amp
        mean_vect_amp += u * u + v * v + w * w;
        namp ++;

    }
    if (namp > 0 && mean_vect_amp > FLT_MIN) {
        polar_analysis->mean_vect_amp = sqrt(mean_vect_amp) / (double) namp ;
    }

    // find eigenvalues and eigenvectors
    VectorDouble S_vector = vector_double(3);
    MatrixDouble V_matrix = matrix_double(3, 3);

    // S_vector - vector of eigenvalues of size isize, in ascending order
    // V_matrix - orthogonal matrix of right singular vectors of size isize x isize
    real_symmetric_eigen_helper(A_matrix, 3, S_vector, V_matrix);

    if (verbose) {
        display_matrix_double("doubleSquareMatrix: ", A_matrix, 3, 3);
        display_vector_double("eigenvalues", S_vector, 3);
        display_vector_double("eigenvector[0] (min)", V_matrix[0], 3);
        display_vector_double("eigenvector[1] (int)", V_matrix[1], 3);
        display_vector_double("eigenvector[2] (max)", V_matrix[2], 3);
    }
    // set eigenvector corresponding to largest eigenvalue
    int indexEigenValueMax = 2;
    double eigenValueMax = S_vector[indexEigenValueMax];
    // check for degenerate solution
    if (eigenValueMax < FLT_MIN) {
        goto cleanup;
    }
    VectorDouble eigenvectorMax = V_matrix[indexEigenValueMax];
    // set eigenvector corresponding to intermediate eigenvalue
    int indexEigenValueInter = 1;
    double eigenValueInter = S_vector[indexEigenValueInter];
    // set eigenvector corresponding to smallest eigenvalue
    int indexEigenValueMin = 0;
    double eigenValueMin = S_vector[indexEigenValueMin];

    // check for effect of float to real roundoff
    double test = 1.0e-10 * eigenValueMax;
    for (int n = 0; n < 3; n++) {
        if (fabs(S_vector[n]) < test) {
            S_vector[n] = 0.0;
        } else if (S_vector[n] < 0.0) {
            fprintf(stdout, "ERROR: polarization_covariance: negative eigenvalue(s): %lf %lf %lf\n", S_vector[0], S_vector[1], S_vector[2]);
            display_matrix_double("doubleSquareMatrix: ", A_matrix, 3, 3);
            display_vector_double("eigenvalues", S_vector, 3);
            display_vector_double("eigenvector[0]", V_matrix[0], 3);
            display_vector_double("eigenvector[1]", V_matrix[1], 3);
            display_vector_double("eigenvector[2]", V_matrix[2], 3);
            continue;
        }
    }
    //fprintf(stdout, "eigenvalues: %lf %lf %lf\n", S_vector[0], S_vector[1], S_vector[2]);

    // normalize
    //double norm = eigenvectorMax.norm();
    //if (norm > Double.MIN_VALUE)
    //    eigenvectorMax = (ComplexVector) eigenvectorMax.scalarDivide(norm);
    // get components
    double rx0 = eigenvectorMax[0];
    double ry0 = eigenvectorMax[1];
    double rz0 = eigenvectorMax[2];
    // make dip angle always positive -> up, always in direction away from source for rays arriving upwards from below surface
    if (rz0 < 0.0) {
        rx0 *= -1.0;
        ry0 *= -1.0;
        rz0 *= -1.0;
    }

    // calculate results values
    // azimuth (0 -> 360 E of N)
    double azimuth = 90.0 - RAD2DEG * atan2(ry0, rx0);
    //printf("DEBUG: azimuth %f", azimuth);
    azimuth += y_azimuth;
    //printf(", y_azimuth-> %f", azimuth);
    int icount = 0;
    while (azimuth < 0.0 && icount++ < 3) {
        azimuth += 360.0;
    }
    icount = 0;
    while (azimuth > 360.0 && icount++ < 3) {
        azimuth -= 360.0;
    }
    //printf(", +/-360-> %f\n", azimuth);
    // dip (-90 - down -> 90 up, following Vidale, 1986)
    double dip;
    double horizMag = sqrt(rx0 * rx0 + ry0 * ry0);
    if (horizMag > FLT_MIN) {
        dip = (float) (RAD2DEG * atan(rz0 / horizMag));
    } else {
        dip = 90.0; // no horizontal magnitude, use up for dip  // 20160308 AJL - added to support 1- and 2- component analysis
    }
    polar_analysis->azimuth = azimuth;
    polar_analysis->dip = dip;
    // linear polarization - degree of linearity
    if (eigenValueMax > FLT_MIN) {
        polar_analysis->degree_of_linearity = 1.0 - eigenValueInter / eigenValueMax;
    }
    // strength of polarization - degree of polarization
    /*if (eigenValueMax > FLT_MIN) {
        polar_analysis->degree_of_planarity = 1.0 - (eigenValueInter + eigenValueMin) / eigenValueMax;
    }*/
    // degree of planarity (Vidale, 1986)
    if (eigenValueInter > FLT_MIN) {
        polar_analysis->degree_of_planarity = 1.0 - eigenValueMin / eigenValueInter;
    }

cleanup:

    free_matrix_double(A_matrix, 3, 3);
    free_vector_double(S_vector);
    free_matrix_double(V_matrix, 3, 3);

    return (polar_analysis);

}
