/***************************************************************************
 * polarisation.h:
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

/*
 * File:   polarization.h
 * Author: anthony
 *
 * Created on August 8, 2016, 9:33 AM
 */

#ifndef POLARIZATION_H
#define POLARIZATION_H

#ifdef __cplusplus
extern "C" {
#endif


// convenience structure to hold PolarizationAnalysis parameters

typedef struct {
    double azimuth; // azimuth (0 -> 360deg E of N)
    double dip; // dip (-90deg - down -> 90deg up, following Vidale, 1986)
    double degree_of_linearity; // linear polarization - degree of linearity
    double degree_of_planarity; // degree of planarity (Vidale, 1986)
    double mean_vect_amp;        // maximum zxy vector amplitude in analysis window
}
PolarizationAnalysis;

PolarizationAnalysis *polarization_covariance(double *data_z, double *data_x, double *data_y, int istart, int nsamp, int ipolarity_z, double y_azimuth, int verbose);


#ifdef __cplusplus
}
#endif

#endif /* POLARIZATION_H */

