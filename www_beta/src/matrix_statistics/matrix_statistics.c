/***************************************************************************
 * matrix_statistics.c:
 *
 * TODO: add doc
 *
 * Written by Anthony Lomax
 *   ALomax Scientific www.alomax.net
 *
 * modified: 2010.12.16
 *           2014.10.30
 ***************************************************************************/




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include "../vector/vector.h"
#include "../geometry/geometry.h"
#include "../alomax_matrix/alomax_matrix.h"
#include "matrix_statistics.h"

#define RA2DE 57.2957795129
#define DE2RA 0.01745329252
#define KM2DEG (90.0/10000.0)
#define DEG2KM (10000.0/90.0)

#ifndef SMALL_DOUBLE
#define SMALL_DOUBLE 1.0e-20
#endif
#ifndef LARGE_DOUBLE
#define LARGE_DOUBLE 1.0e20
#endif

static char error_message[4096];

/** function to print error and return last error message */
char *get_matrix_statistics_error_mesage() {
    return (error_message);
}

/** function to calculate the expectation (mean)  of a set of samples */

Vect3D CalcExpectationSamples(float* fdata, int nSamples) {

    int nsamp, ipos;

    float x, y, z, prob;
    Vect3D expect = {0.0, 0.0, 0.0};


    ipos = 0;
    for (nsamp = 0; nsamp < nSamples; nsamp++) {
        x = fdata[ipos++];
        y = fdata[ipos++];
        z = fdata[ipos++];
        prob = fdata[ipos++];
        expect.x += (double) x;
        expect.y += (double) y;
        expect.z += (double) z;
    }

    expect.x /= (double) nSamples;
    expect.y /= (double) nSamples;
    expect.z /= (double) nSamples;

    return (expect);
}

/** function to calculate the expectation (mean)  of a set of samples */

Vect3D CalcExpectationSamplesWeighted(float* fdata, int nSamples) {

    int nsamp, ipos;

    float x, y, z;
    Vect3D expect = {0.0, 0.0, 0.0};

    double weight;
    double weight_sum = 0.0;

    ipos = 0;
    for (nsamp = 0; nsamp < nSamples; nsamp++) {
        x = fdata[ipos++];
        y = fdata[ipos++];
        z = fdata[ipos++];
        weight = fdata[ipos++];
        expect.x += (double) x * weight;
        expect.y += (double) y * weight;
        expect.z += (double) z * weight;
        weight_sum += weight;
    }

    expect.x /= weight_sum;
    expect.y /= weight_sum;
    expect.z /= weight_sum;

    return (expect);
}

/** function to calculate the expectation (mean) of a set of samples (lon,lat,depth,weight)
 *
 * global case - checks for wrap around in longitude (x) using specified xReference as correct longitude zone
 * TODO: uses rectangular lat/lon geometry, does not try and correct for change in longitude distance with latitude.
 * TODO: does not try and correct for problems in latitude near poles.
 *
 */

Vect3D CalcExpectationSamplesGlobal(float* fdata, int nSamples, double xReference) {

    int nsamp, ipos;

    double x, y, z;
    Vect3D expect = {0.0, 0.0, 0.0};

    ipos = 0;
    for (nsamp = 0; nsamp < nSamples; nsamp++) {
        x = fdata[ipos++];
        if (x - xReference > 180.0)
            x -= 360.0;
        else if (x - xReference < -180.0)
            x += 360.0;
        y = fdata[ipos++];
        z = fdata[ipos++];
        ipos++; // fdata value is in 4th position
        expect.x += x;
        expect.y += y;
        expect.z += z;
    }

    expect.x /= (double) nSamples;
    expect.y /= (double) nSamples;
    expect.z /= (double) nSamples;

    return (expect);
}

/** function to calculate the weighted expectation (mean) of a set of samples (lon,lat,depth,weight)
 *
 * global case - checks for wrap around in longitude (x) using specified xReference as correct longitude zone
 * TODO: does not try and correct for problems in latitude near poles.
 *
 */

Vect3D CalcExpectationSamplesGlobalWeighted(float* fdata, int nSamples, double xReference) {

    int nsamp, ipos;

    double x, y, z;
    Vect3D expect = {0.0, 0.0, 0.0};

    double weight;
    double weight_sum = 0.0;

    ipos = 0;
    for (nsamp = 0; nsamp < nSamples; nsamp++) {
        x = fdata[ipos++];
        if (x - xReference > 180.0)
            x -= 360.0;
        else if (x - xReference < -180.0)
            x += 360.0;
        y = fdata[ipos++];
        z = fdata[ipos++];
        weight = fdata[ipos++];
        expect.x += x * weight;
        expect.y += y * weight;
        expect.z += z * weight;
        weight_sum += weight;
    }

    expect.x /= weight_sum;
    expect.y /= weight_sum;
    expect.z /= weight_sum;

    return (expect);
}

/** function to calculate the covariance of a set of samples assumed to be distributed following a target PDF
 *
 * 20141030 AJL - Bug fix: new version which subtracts the expectation from each data value before summing,
 *      instead of correcting for expectation after summing and dividing by nSamples.
 *      Should prevent precision errors when expectation is far from coordinates origin.
 */

Mtrx3D CalcCovarianceSamplesRect(float* fdata, int nSamples, Vect3D* pexpect) {

    int nsamp, ipos;

    float x, y, z, prob;
    Mtrx3D cov = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};


    /* calculate covariance following eq. (6-12), T & V, 1982 */

    ipos = 0;
    for (nsamp = 0; nsamp < nSamples; nsamp++) {
        x = fdata[ipos++] - pexpect->x;
        y = fdata[ipos++] - pexpect->y;
        z = fdata[ipos++] - pexpect->z;
        prob = fdata[ipos++]; // do not use prob since samples follow target PDF

        cov.xx += (double) (x * x);
        cov.xy += (double) (x * y);
        cov.xz += (double) (x * z);

        cov.yy += (double) (y * y);
        cov.yz += (double) (y * z);

        cov.zz += (double) (z * z);

    }

    cov.xx = cov.xx / (double) nSamples;
    cov.xy = cov.xy / (double) nSamples;
    cov.xz = cov.xz / (double) nSamples;

    cov.yx = cov.xy;
    cov.yy = cov.yy / (double) nSamples;
    cov.yz = cov.yz / (double) nSamples;

    cov.zx = cov.xz;
    cov.zy = cov.yz;
    cov.zz = cov.zz / (double) nSamples;


    return (cov);
}

/** function to calculate the covariance of a set of samples
 *
 * !!! DO NOT USE - subject to precision errors!
 */

Mtrx3D CalcCovarianceSamplesRect_OLD(float* fdata, int nSamples, Vect3D* pexpect) {

    int nsamp, ipos;

    float x, y, z, prob;
    Mtrx3D cov = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};


    /* calculate covariance following eq. (6-12), T & V, 1982 */

    ipos = 0;
    for (nsamp = 0; nsamp < nSamples; nsamp++) {
        x = fdata[ipos++];
        y = fdata[ipos++];
        z = fdata[ipos++];
        prob = fdata[ipos++];

        cov.xx += (double) (x * x);
        cov.xy += (double) (x * y);
        cov.xz += (double) (x * z);

        cov.yy += (double) (y * y);
        cov.yz += (double) (y * z);

        cov.zz += (double) (z * z);

    }

    cov.xx = cov.xx / (double) nSamples - pexpect->x * pexpect->x;
    cov.xy = cov.xy / (double) nSamples - pexpect->x * pexpect->y;
    cov.xz = cov.xz / (double) nSamples - pexpect->x * pexpect->z;

    cov.yx = cov.xy;
    cov.yy = cov.yy / (double) nSamples - pexpect->y * pexpect->y;
    cov.yz = cov.yz / (double) nSamples - pexpect->y * pexpect->z;

    cov.zx = cov.xz;
    cov.zy = cov.yz;
    cov.zz = cov.zz / (double) nSamples - pexpect->z * pexpect->z;


    return (cov);
}

/** simple great-circle distance and azimuth calculation on a sphere
 *
 * returns distance and azimuth in degrees for great-circle from latA/lonA to latB/lonB
 *
 */

double GCDistanceAzimuth__(double latA, double lonA, double latB, double lonB, double *pazimuth) {

    lonA *= DE2RA;
    latA *= DE2RA;
    lonB *= DE2RA;
    latB *= DE2RA;

    // distance
    double dist = sin(latA) * sin(latB) + cos(latA) * cos(latB) * cos(lonA - lonB);
    dist = acos(dist);

    // 20141106 AJL - added following check, to prevent div by 0 of sin(dist))
    if (dist < FLT_MIN) {
        *pazimuth = 0.0;
        return (dist * RA2DE);
    }

    // azimuth
    double cosAzimuth =
            (cos(latA) * sin(latB)
            - sin(latA) * cos(latB)
            * cos((lonB - lonA)))
            / sin(dist);
    double sinAzimuth =
            cos(latB) * sin((lonB - lonA)) / sin(dist);
    double az = atan2(sinAzimuth, cosAzimuth) * RA2DE;

    if (isnan(az) && fabs(lonB - lonA) < 0.000001) {
        if (latA > latB) {
            az = 180.0;
        } else {
            az = 0.0;
        }
    }

    if (az < 0.0) {
        az += 360.0;
    }

    *pazimuth = az;
    return (dist * RA2DE);

}




/** function to calculate the covariance of a set of samples in long(deg)/lat(deg)/depth(km) coordinates
 * samples assumed to be distributed following a target PDF
 *
 * 20141107 AJL - New version which calculates x and y coords using distance and azimuth of each sample
 *      from expectation x,y.  Centers x,y and avoids problems near poles.
 *
 */
//Mtrx3D CalcCovarianceSamplesGlobal_GCD(float* fdata, int nSamples, Vect3D* pexpect) {

Mtrx3D CalcCovarianceSamplesGlobal(float* fdata, int nSamples, Vect3D* pexpect) {

    int nsamp, ipos;

    double lat, lon, x, y, z, prob;
    Mtrx3D cov = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    double distance, azimuth;
    double xReference = pexpect->x;

    // calculate covariance following eq. (6-12), T & V, 1982

    ipos = 0;
    for (nsamp = 0; nsamp < nSamples; nsamp++) {

        lon = fdata[ipos++];
        if (lon - xReference > 180.0)
            lon -= 360.0;
        else if (lon - xReference < -180.0)
            lon += 360.0;
        lat = fdata[ipos++];

        distance = GCDistanceAzimuth__(pexpect->y, pexpect->x, lat, lon, &azimuth);
        distance *= DEG2KM;

        x = distance * sin(azimuth * DE2RA); // azimuth is deg CW from North
        y = distance * cos(azimuth * DE2RA);
        z = (fdata[ipos++] - pexpect->z);
        prob = fdata[ipos++]; // do not use prob since samples follow target PDF

        cov.xx += (double) (x * x);
        cov.xy += (double) (x * y);
        cov.xz += (double) (x * z);

        cov.yy += (double) (y * y);
        cov.yz += (double) (y * z);

        cov.zz += (double) (z * z);

    }

    cov.xx = cov.xx / (double) nSamples;
    cov.xy = cov.xy / (double) nSamples;
    cov.xz = cov.xz / (double) nSamples;

    cov.yx = cov.xy;
    cov.yy = cov.yy / (double) nSamples;
    cov.yz = cov.yz / (double) nSamples;

    cov.zx = cov.xz;
    cov.zy = cov.yz;
    cov.zz = cov.zz / (double) nSamples;


    return (cov);
}

/** function to calculate the covariance of a set of samples in long(deg)/lat(deg)/depth(km) coordinates
 * samples assumed to be distributed following a target PDF
 *
 * 20141030 AJL - Bug fix: new version which subtracts the expectation from each data value before summing,
 *      instead of correcting for expectation after summing and dividing by nSamples.
 *      Should prevent precision errors when expectation is far from coordinates origin.
 */

Mtrx3D CalcCovarianceSamplesGlobal_NEW(float* fdata, int nSamples, Vect3D* pexpect) {
    //Mtrx3D CalcCovarianceSamplesGlobal(float* fdata, int nSamples, Vect3D* pexpect) {

    int nsamp, ipos;

    float x, y, z, prob;
    Mtrx3D cov = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    double cos_lat = cos(pexpect->y * DE2RA);
    double xReference = pexpect->x;

    /* calculate covariance following eq. (6-12), T & V, 1982 */

    ipos = 0;
    for (nsamp = 0; nsamp < nSamples; nsamp++) {
        x = fdata[ipos++];
        if (x - xReference > 180.0)
            x -= 360.0;
        else if (x - xReference < -180.0)
            x += 360.0;
        x = (x - pexpect->x) * DEG2KM * cos_lat;
        y = (fdata[ipos++] - pexpect->y) * DEG2KM;
        z = (fdata[ipos++] - pexpect->z);
        prob = fdata[ipos++]; // do not use prob since samples follow target PDF

        cov.xx += (double) (x * x);
        cov.xy += (double) (x * y);
        cov.xz += (double) (x * z);

        cov.yy += (double) (y * y);
        cov.yz += (double) (y * z);

        cov.zz += (double) (z * z);

    }

    cov.xx = cov.xx / (double) nSamples;
    cov.xy = cov.xy / (double) nSamples;
    cov.xz = cov.xz / (double) nSamples;

    cov.yx = cov.xy;
    cov.yy = cov.yy / (double) nSamples;
    cov.yz = cov.yz / (double) nSamples;

    cov.zx = cov.xz;
    cov.zy = cov.yz;
    cov.zz = cov.zz / (double) nSamples;


    return (cov);
}

/** function to calculate the covariance of a set of samples in long(deg)/lat(deg)/depth(km) coordinates
 *
 * !!! DO NOT USE - subject to precision errors!
 */

Mtrx3D CalcCovarianceSamplesGlobal_OLD(float* fdata, int nSamples, Vect3D* pexpect) {
    //Mtrx3D CalcCovarianceSamplesGlobal(float* fdata, int nSamples, Vect3D* pexpect) {

    int nsamp, ipos;

    float x, y, z, prob;
    Mtrx3D cov = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    double cos_lat = cos(pexpect->y * DE2RA);
    double xReference = pexpect->x;

    /* calculate covariance following eq. (6-12), T & V, 1982 */

    ipos = 0;
    for (nsamp = 0; nsamp < nSamples; nsamp++) {
        x = fdata[ipos++];
        if (x - xReference > 180.0)
            x -= 360.0;
        else if (x - xReference < -180.0)
            x += 360.0;
        x = x * DEG2KM * cos_lat;
        y = fdata[ipos++] * DEG2KM;
        z = fdata[ipos++];
        prob = fdata[ipos++];

        cov.xx += (double) (x * x);
        cov.xy += (double) (x * y);
        cov.xz += (double) (x * z);

        cov.yy += (double) (y * y);
        cov.yz += (double) (y * z);

        cov.zz += (double) (z * z);

    }

    cov.xx = cov.xx / (double) nSamples - pexpect->x * pexpect->x * DEG2KM * cos_lat * DEG2KM * cos_lat;
    cov.xy = cov.xy / (double) nSamples - pexpect->x * pexpect->y * DEG2KM * cos_lat * DEG2KM;
    cov.xz = cov.xz / (double) nSamples - pexpect->x * pexpect->z * DEG2KM * cos_lat;

    cov.yx = cov.xy;
    cov.yy = cov.yy / (double) nSamples - pexpect->y * pexpect->y * DEG2KM * DEG2KM;
    cov.yz = cov.yz / (double) nSamples - pexpect->y * pexpect->z * DEG2KM;

    cov.zx = cov.xz;
    cov.zy = cov.yz;
    cov.zz = cov.zz / (double) nSamples - pexpect->z * pexpect->z;


    return (cov);
}

/** function to calculate the covariance of a set of samples in long(deg)/lat(deg)/depth(km) coordinates
 *
 * 20141030 AJL - Bug fix: new version which subtracts the expectation from each data value before summing,
 *      instead of correcting for expectation after summing and dividing by nSamples.
 *      Should prevent precision errors when expectation is far from coordinates origin.
 */

Mtrx3D CalcCovarianceSamplesGlobalWeighted(float* fdata, int nSamples, Vect3D* pexpect) {

    int nsamp, ipos;

    double x, y, z;
    Mtrx3D cov = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    double weight;
    double weight_sum = 0.0;

    double cos_lat = cos(pexpect->y * DE2RA);
    double xReference = pexpect->x;

    /* calculate covariance following eq. (6-12), T & V, 1982 */

    ipos = 0;
    for (nsamp = 0; nsamp < nSamples; nsamp++) {
        x = fdata[ipos++];
        if (x - xReference > 180.0)
            x -= 360.0;
        else if (x - xReference < -180.0)
            x += 360.0;
        x = (x - pexpect->x) * DEG2KM * cos_lat;
        y = (fdata[ipos++] - pexpect->y) * DEG2KM;
        z = (fdata[ipos++] - pexpect->z);
        weight = fdata[ipos++];

        cov.xx += (x * x) * weight;
        cov.xy += (x * y) * weight;
        cov.xz += (x * z) * weight;

        cov.yy += (y * y) * weight;
        cov.yz += (y * z) * weight;

        cov.zz += (z * z) * weight;

        weight_sum += weight;

    }

    cov.xx = cov.xx / weight_sum;
    cov.xy = cov.xy / weight_sum;
    cov.xz = cov.xz / weight_sum;

    cov.yx = cov.xy;
    cov.yy = cov.yy / weight_sum;
    cov.yz = cov.yz / weight_sum;

    cov.zx = cov.xz;
    cov.zy = cov.yz;
    cov.zz = cov.zz / weight_sum;


    return (cov);
}

/** function to calculate the covariance of a set of samples in long(deg)/lat(deg)/depth(km) coordinates
 *
 * !!! DO NOT USE - subject to precision errors!
 */

Mtrx3D CalcCovarianceSamplesGlobalWeighted_OLD(float* fdata, int nSamples, Vect3D* pexpect) {

    int nsamp, ipos;

    double x, y, z;
    Mtrx3D cov = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    double weight;
    double weight_sum = 0.0;

    double cos_lat = cos(pexpect->y * DE2RA);
    double xReference = pexpect->x;

    /* calculate covariance following eq. (6-12), T & V, 1982 */

    ipos = 0;
    for (nsamp = 0; nsamp < nSamples; nsamp++) {
        x = fdata[ipos++];
        if (x - xReference > 180.0)
            x -= 360.0;
        else if (x - xReference < -180.0)
            x += 360.0;
        x = x * DEG2KM * cos_lat;
        y = fdata[ipos++] * DEG2KM;
        z = fdata[ipos++];
        weight = fdata[ipos++];

        cov.xx += (x * x) * weight;
        cov.xy += (x * y) * weight;
        cov.xz += (x * z) * weight;

        cov.yy += (y * y) * weight;
        cov.yz += (y * z) * weight;

        cov.zz += (z * z) * weight;

        weight_sum += weight;

    }

    cov.xx = cov.xx / weight_sum - pexpect->x * pexpect->x * DEG2KM * cos_lat * DEG2KM * cos_lat;
    cov.xy = cov.xy / weight_sum - pexpect->x * pexpect->y * DEG2KM * cos_lat * DEG2KM;
    cov.xz = cov.xz / weight_sum - pexpect->x * pexpect->z * DEG2KM * cos_lat;

    cov.yx = cov.xy;
    cov.yy = cov.yy / weight_sum - pexpect->y * pexpect->y * DEG2KM * DEG2KM;
    cov.yz = cov.yz / weight_sum - pexpect->y * pexpect->z * DEG2KM;

    cov.zx = cov.xz;
    cov.zy = cov.yz;
    cov.zz = cov.zz / weight_sum - pexpect->z * pexpect->z;


    return (cov);
}



/** function to calculate confidence ellipsoid from covariance matrix */

/* 	finds confidence ellipsoid from SVD of Cov mtrx.  See
                Num Rec, 2nd ed, secs 2.6 & 15.6

                del_chi_2 is delta Chi-square (see Num Rec, 2nd ed, fig 15.6.5)
 */


Ellipse2D CalcHorizontalErrorEllipse(Mtrx3D *pcov, double del_chi_2) {

    int ndx, iSwitched;
    MatrixDouble A_matrix, V_matrix;
    VectorDouble W_vector;
    double wtemp, vtemp;
    Ellipse2D ell;

    int ierr = 0;


    /* allocate A mtrx */
    A_matrix = matrix_double(2, 2);

    /* load A matrix in NumRec format */
    A_matrix[0][0] = pcov->xx;
    A_matrix[0][1] = A_matrix[1][0] = pcov->xy;
    A_matrix[1][1] = pcov->yy;


    /* allocate V mtrx and W vector */
    V_matrix = matrix_double(2, 2);
    W_vector = vector_double(2);

    /* do SVD */
    //if ((istat = nll_svdcmp0(A_matrix, 2, 2, W_vector, V_matrix)) < 0) {
    svd_helper(A_matrix, 2, 2, W_vector, V_matrix);
    if (W_vector[0] < SMALL_DOUBLE || W_vector[1] < SMALL_DOUBLE) {
        fprintf(stderr, "ERROR: invalid SVD singular value for confidence ellipsoids.");
        ierr = 1;
    } else {

        /* sort by singular values W */
        iSwitched = 1;
        while (iSwitched) {
            iSwitched = 0;
            for (ndx = 0; ndx < 1; ndx++) {
                if (W_vector[ndx] > W_vector[ndx + 1]) {
                    wtemp = W_vector[ndx];
                    W_vector[ndx] = W_vector[ndx + 1];
                    W_vector[ndx + 1] = wtemp;
                    vtemp = V_matrix[0][ndx];
                    V_matrix[0][ndx] = V_matrix[0][ndx + 1];
                    V_matrix[0][ndx + 1] = vtemp;
                    vtemp = V_matrix[1][ndx];
                    V_matrix[1][ndx] = V_matrix[1][ndx + 1];
                    V_matrix[1][ndx + 1] = vtemp;
                    iSwitched = 1;
                }
            }
        }


        /* calculate ellipsoid axes */
        /* length: w in Num Rec, 2nd ed, fig 15.6.5 must be replaced
                by 1/sqrt(w) since we are using SVD of Cov mtrx and not
                SVD of A mtrx (compare eqns 2.6.1  & 15.6.10) */

        ell.az1 = atan2(V_matrix[0][0], V_matrix[1][0]) * RA2DE;
        if (ell.az1 < 0.0)
            ell.az1 += 360.0;
        else if (ell.az1 >= 360.0)
            ell.az1 -= 360.0;
        if (ell.az1 >= 180.0) // force in range [0, 180)
            ell.az1 -= 180.0;
        ell.len1 = sqrt(del_chi_2) / sqrt(1.0 / W_vector[0]);
        ell.len2 = sqrt(del_chi_2) / sqrt(1.0 / W_vector[1]);

    }

    free_matrix_double(A_matrix, 2, 2);
    free_matrix_double(V_matrix, 2, 2);
    free_vector_double(W_vector);

    if (ierr) {
        Ellipse2D EllipseNULL = {-1.0, -1.0, -1.0};
        return (EllipseNULL);
    }

    return (ell);

}

/** function to calculate confidence ellipsoid from covariance matrix */

/* 	finds confidence ellipsoid from SVD of Cov mtrx.  See
                Num Rec, 2nd ed, secs 2.6 & 15.6

                del_chi_2 is delta Chi-square (see Num Rec, 2nd ed, fig 15.6.5)
 */


Ellipsoid3D CalcErrorEllipsoid(Mtrx3D *pcov, double del_chi_2) {
    int ndx, iSwitched;
    MatrixDouble A_matrix, V_matrix;
    VectorDouble W_vector;
    double wtemp, vtemp;
    Ellipsoid3D ell;

    int ierr = 0;


    /* allocate A mtrx */
    A_matrix = matrix_double(3, 3);

    /* load A matrix in NumRec format */
    A_matrix[0][0] = pcov->xx;
    A_matrix[0][1] = A_matrix[1][0] = pcov->xy;
    A_matrix[0][2] = A_matrix[2][0] = pcov->xz;
    A_matrix[1][1] = pcov->yy;
    A_matrix[1][2] = A_matrix[2][1] = pcov->yz;
    A_matrix[2][2] = pcov->zz;


    /* allocate V mtrx and W vector */
    V_matrix = matrix_double(3, 3);
    W_vector = vector_double(3);

    /* do SVD */
    //if ((istat = nll_svdcmp0(A_matrix, 3, 3, W_vector, V_matrix)) < 0) {
    svd_helper(A_matrix, 3, 3, W_vector, V_matrix);
    if (W_vector[0] < SMALL_DOUBLE || W_vector[1] < SMALL_DOUBLE || W_vector[2] < SMALL_DOUBLE) {
        fprintf(stderr, "ERROR: invalid SVD singular value for confidence ellipsoids.");
        ierr = 1;
    } else {

        /* sort by singular values W */
        iSwitched = 1;
        while (iSwitched) {
            iSwitched = 0;
            for (ndx = 0; ndx < 2; ndx++) {
                if (W_vector[ndx] > W_vector[ndx + 1]) {
                    wtemp = W_vector[ndx];
                    W_vector[ndx] = W_vector[ndx + 1];
                    W_vector[ndx + 1] = wtemp;
                    vtemp = V_matrix[0][ndx];
                    V_matrix[0][ndx] = V_matrix[0][ndx + 1];
                    V_matrix[0][ndx + 1] = vtemp;
                    vtemp = V_matrix[1][ndx];
                    V_matrix[1][ndx] = V_matrix[1][ndx + 1];
                    V_matrix[1][ndx + 1] = vtemp;
                    vtemp = V_matrix[2][ndx];
                    V_matrix[2][ndx] = V_matrix[2][ndx + 1];
                    V_matrix[2][ndx + 1] = vtemp;
                    iSwitched = 1;
                }
            }
        }


        /* calculate ellipsoid axes */
        /* length: w in Num Rec, 2nd ed, fig 15.6.5 must be replaced
                by 1/sqrt(w) since we are using SVD of Cov mtrx and not
                SVD of A mtrx (compare eqns 2.6.1  & 15.6.10) */

        ell.az1 = atan2(V_matrix[0][0], V_matrix[1][0]) * RA2DE;
        if (ell.az1 < 0.0)
            ell.az1 += 360.0;
        ell.dip1 = asin(V_matrix[2][0]) * RA2DE;
        ell.len1 = sqrt(del_chi_2) / sqrt(1.0 / W_vector[0]);
        ell.az2 = atan2(V_matrix[0][1], V_matrix[1][1]) * RA2DE;
        if (ell.az2 < 0.0)
            ell.az2 += 360.0;
        ell.dip2 = asin(V_matrix[2][1]) * RA2DE;
        ell.len2 = sqrt(del_chi_2) / sqrt(1.0 / W_vector[1]);
        ell.len3 = sqrt(del_chi_2) / sqrt(1.0 / W_vector[2]);

        // 20150601 AJL - semi-major axis az and dip added to support conversion to QuakeML Tait-Bryan representation
        ell.az3 = atan2(V_matrix[0][2], V_matrix[1][2]) * RA2DE;
        ell.dip3 = asin(V_matrix[2][2]) * RA2DE;


    }

    free_matrix_double(A_matrix, 3, 3);
    free_matrix_double(V_matrix, 3, 3);
    free_vector_double(W_vector);

    if (ierr) {
        Ellipsoid3D EllipsoidNULL = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
        return (EllipsoidNULL);
    }

    return (ell);

}


/** method to convert ellipsoid axes in az/dip/se to vector axes */

/* converted from method init in Java class Ellipsoid3D (04DEC1998) */

void ellipsiod2Axes(Ellipsoid3D *pellipsoid, Vect3D *paxis1, Vect3D *paxis2, Vect3D *paxis3) {

    double az1, az2, dip1, dip2;
    double cosd1, cosd2;


    /* strike angles converted to positive CCW from East = 0 */
    az1 = 90.0 - pellipsoid->az1;
    az2 = 90.0 - pellipsoid->az2;
    /* dip angles increasing downwards from horiz = 0, convert to upwards */
    dip1 = -pellipsoid->dip1;
    dip2 = -pellipsoid->dip2;

    /* get 3D vector axes */

    cosd1 = cos(DE2RA * dip1);
    paxis1->x = cos(DE2RA * az1) * cosd1;
    paxis1->y = sin(DE2RA * az1) * cosd1;
    paxis1->z = -sin(DE2RA * dip1);

    cosd2 = cos(DE2RA * dip2);
    paxis2->x = cos(DE2RA * az2) * cosd2;
    paxis2->y = sin(DE2RA * az2) * cosd2;
    paxis2->z = -sin(DE2RA * dip2);

    cross_product_3d(
            paxis1->x, paxis1->y, paxis1->z,
            paxis2->x, paxis2->y, paxis2->z,
            &paxis3->x, &paxis3->y, &paxis3->z);

    paxis1->x *= pellipsoid->len1;
    paxis1->y *= pellipsoid->len1;
    paxis1->z *= pellipsoid->len1;
    paxis2->x *= pellipsoid->len2;
    paxis2->y *= pellipsoid->len2;
    paxis2->z *= pellipsoid->len2;
    paxis3->x *= pellipsoid->len3;
    paxis3->y *= pellipsoid->len3;
    paxis3->z *= pellipsoid->len3;

}

/** method to convert ellipsoid to an XML (pseudo-QuakeML) ConfidenceEllipsoid
 *
 *  !!! Very incomplete.  Only converts NLL Ellipsoid axes parameters minor(3)/intermediate(3)/major(1) to
 *      major(3)/intermediate(3)/minor(1) ordering

 */

void nllEllipsiod2XMLConfidenceEllipsoid(Ellipsoid3D *pellipsoid,
        double* psemiMajorAxisLength, double* pmajorAxisPlunge, double* pmajorAxisAzimuth,
        double* psemiIntermediateAxisLength, double* pintermediateAxisPlunge, double* pintermediateAxisAzimuth,
        double* psemiMinorAxisLength) {

    Vect3D axis1;
    Vect3D axis2;
    Vect3D axis3;
    ellipsiod2Axes(pellipsoid, &axis1, &axis2, &axis3);

    *psemiMajorAxisLength = pellipsoid->len3;
    *psemiIntermediateAxisLength = pellipsoid->len2;
    *psemiMinorAxisLength = pellipsoid->len1;

    double plunge = axis3.z >= 0.0 ? 90.0 : -90.0;
    double hypot = sqrt(axis3.x * axis3.x + axis3.y * axis3.y);
    if (hypot > FLT_MIN) {
        plunge = RA2DE * atan(axis3.z / hypot);
    }
    double azim = RA2DE * atan2(axis3.x, axis3.y);
    if (azim < 0.0)
        azim += 360.0;
    if (plunge < 0.0) {
        plunge *= -1.0;
        azim -= 180.0;
        if (azim < 0.0)
            azim += 360.0;
    }
    *pmajorAxisPlunge = plunge;
    *pmajorAxisAzimuth = azim;


    plunge = axis2.z >= 0.0 ? 90.0 : -90.0;
    hypot = sqrt(axis2.x * axis2.x + axis2.y * axis2.y);
    if (hypot > FLT_MIN) {
        plunge = RA2DE * atan(axis2.z / hypot);
    }
    azim = RA2DE * atan2(axis2.x, axis2.y);
    if (azim < 0.0)
        azim += 360.0;
    if (plunge < 0.0) {
        plunge *= -1.0;
        azim -= 180.0;
        if (azim < 0.0)
            azim += 360.0;
    }
    *pintermediateAxisPlunge = plunge;
    *pintermediateAxisAzimuth = azim;

}

/** method to perform dot product on two 3x3 matrices
 *
 * WARING: no check if input matrices are 3x3!
 */

int matrix_dot_3_3(MatrixDouble first, MatrixDouble second, MatrixDouble mtx_dot) {

    int c, d, k;
    double sum;
    for (c = 0; c < 3; c++) {
        for (d = 0; d < 3; d++) {
            sum = 0.0;
            for (k = 0; k < 3; k++) {
                sum = sum + first[c][k] * second[k][d];
            }
            mtx_dot[c][d] = sum;
        }
    }

    return (0);

}

/** method to convert NonLinLoc/Hypoellipse ellipsoid to an into the Tait-Bryan representation of QuakeML ConfidenceEllipsoid
 *
 * sets QuakeML ConfidenceEllipsoid fields
 * returns 0 on success
 *          -1 on error
 */

int nllEllipsiod2QMLConfidenceEllipsoid(Ellipsoid3D *pellipsoid,
        double* psemiMajorAxisLength,
        double* psemiMinorAxisLength,
        double* psemiIntermediateAxisLength,
        double* pmajorAxisAzimuth,
        double* pmajorAxisPlunge,
        double* pmajorAxisRotation) {

    /*
       // 3D ellipsoid
       typedef struct {
           double az1, dip1, len1;   // semi-minor axis km
           double az2, dip2, len2;   // semi-intermediate axis km
           double len3;   // semi-major axis km
           double az3, dip3;   // 20150601 AJL - semi-major axis az and dip added to support conversion to QuakeML Tait-Bryan representation
       } Ellipsoid3D;
     */
    /* adapted from: https://github.com/usgs/libcomcat libcomcat.ellipse.py
     * by A Lomax / ISTI 20150601

                def vec2tait(ellipsoid):
                 """Convert earthquake origin error ellipsoid 3x3 matrix into Tait-Bryan representation.

                 Input argument:
                 ellipsoid - 3x3 Numpy array containing the elements describing earthquake origin error ellipsoid.
                 [SemiMajorAxisAzimuth SemiMajorAxisPlunge SemiMajorAxisLength;
                  SemiMinorAxisAzimuth SemiMinorAxisPlunge SemiMinorAxisLength;
                  SemiIntermediateAxisAzimuth SemiIntermediateAxisPlunge SemiIntermediateAxisLength]
                 (distance values in kilometers, angles in degrees)

                 Returns: 6 element tuple:
                 semiMajorAxisLength (km)
                 semiMinorAxisLength (km)
                 semiIntermediateAxisLength (km)
                 majorAxisAzimuth (degrees)
                 majorAxisPlunge (degrees)
                 majorAxisRotation (degrees)

                 #      Tait-Bryan angles use in Section 3.3.9 ConfidenceEllipsoid
                 #      QuakeML-BED-20130214a.pdf
                 #
                 #      The image showing two views of the rotation in the QuakeML document is wrong
                 #      and disagrees with the text, which says
                 #      X is major axis, Y is the minor axis and thus Z is the intermediate axis
                 #      The rotations are as follow:
                 #      1. about z-axis by angle PSI (heading) to give (x', y', z)
                 #      2, about y' with angle with angle PHI (elevation) to give (x'', y', z'')
                 #      3. about x'' with angle THETA (bank) to give (x'', y''', z''')
                 #
                 #      This sequence is known as z-y'-x'' intrinsic rotation (http://en.wikipedia.org/wiki/Euler_angles)
                 #      The figure in the QuakeML document is the z-x'-y'' rotation. Note the order
                 #
                 #      azimuth is measured positive in direction from x to y
                 #      plunge is measured positive such that the x' moves in the positive z-direction
                 #      rotation is measured positive such that the y'' moves in the positive z-direction
                 #
                 #
                 # azimuth is heading, measure in degrees from north
                 # plunge is elevation
                 # rotation is roll

                 # Author: R.B.Herrmann (rbh@eas.slu.edu)
                 # Created: 23 November 2014
                 # Adapted-By: rbh
                 # Translated to Python by Mike Hearne (mhearne@usgs.gov)

     */

    double PHI, PSI;
    double PHImaj, PSImaj;
    double PHIint, PSIint;
    double PHImin, PSImin;


    //ellipsearray = ellipsoid.flatten()

    // get indices of the minor, intermediate, major axes
    //smallest,intermediate,largest = ellipsearray[2:9:3].argsort()

    // we have two of the angles already
    // convert from km to meters for the lengths
    //
    //k=3*(largest)
    *psemiMajorAxisLength = pellipsoid->len3;
    *psemiMinorAxisLength = pellipsoid->len1;
    *psemiIntermediateAxisLength = pellipsoid->len2;
    *pmajorAxisPlunge = pellipsoid->dip3; // PHI
    *pmajorAxisAzimuth = pellipsoid->az3; // PSI
    if (*pmajorAxisAzimuth >= 360.0) {
        *pmajorAxisAzimuth -= 360.0;
    } else if (*pmajorAxisAzimuth < 0.0) {
        *pmajorAxisAzimuth += 360.0;
    }

    PHI = *pmajorAxisPlunge;
    PSI = *pmajorAxisAzimuth;

    //
    //    to get the THETA angle, we need to do some transformations
    //
    double RPHI[3][3] = {
        {cos(DE2RA * PHI), 0, sin(DE2RA * PHI)},
        {0, 1, 0},
        {sin(DE2RA * PHI), 0, cos(DE2RA * PHI)}
    };
    double RPSI[3][3] = {
        {cos(DE2RA * PSI), sin(DE2RA * PSI), 0},
        {-sin(DE2RA * PSI), cos(DE2RA * PSI), 0},
        {0, 0, 1}
    };

    //
    //    reconstruct the T matrix
    //
    // major axis
    //k=3*(largest)
    PHImaj = pellipsoid->dip3;
    PSImaj = pellipsoid->az3;
    MatrixDouble T = matrix_double(3, 3);
    T[0][0] = cos(DE2RA * PSImaj) * cos(DE2RA * PHImaj);
    T[0][1] = sin(DE2RA * PSImaj) * cos(DE2RA * PHImaj);
    T[0][2] = sin(DE2RA * PHImaj);

    // minor axis
    //k=3*(smallest)
    PHImin = pellipsoid->dip1;
    PSImin = pellipsoid->az1;
    T[1][0] = cos(DE2RA * PSImin) * cos(DE2RA * PHImin);
    T[1][1] = sin(DE2RA * PSImin) * cos(DE2RA * PHImin);
    T[1][2] = sin(DE2RA * PHImin);

    // minor axis
    //k=3*(intermediate)
    PHIint = pellipsoid->dip2;
    PSIint = pellipsoid->az2;
    T[2][0] = cos(DE2RA * PSIint) * cos(DE2RA * PHIint);
    T[2][1] = sin(DE2RA * PSIint) * cos(DE2RA * PHIint);
    T[2][2] = sin(DE2RA * PHIint);

    // set an invert RPSI
    int nrow, ncol;
    MatrixDouble inv_RPSI = matrix_double(3, 3);
    for (nrow = 0; nrow < 3; nrow++) {
        for (ncol = 0; ncol < 3; ncol++) {
            inv_RPSI[nrow][ncol] = RPSI[nrow][ncol];
        }
    }
    // invert
    if (matrix_double_inverse(inv_RPSI, 3, 3) < 0) {
        snprintf(error_message, sizeof (error_message), "ERROR: in matrix_double_check_diagonal_non_zero_inverse()");
        return (-1);
    }
    // set an invert RPHI
    MatrixDouble inv_RPHI = matrix_double(3, 3);
    for (nrow = 0; nrow < 3; nrow++) {
        for (ncol = 0; ncol < 3; ncol++) {
            inv_RPHI[nrow][ncol] = RPHI[nrow][ncol];
        }
    }
    // invert
    if (matrix_double_inverse(inv_RPHI, 3, 3) < 0) {
        snprintf(error_message, sizeof (error_message), "ERROR: in matrix_double_check_diagonal_non_zero_inverse()");
        return (-1);
    }
    // dot
    MatrixDouble dot_inv_RPSI_inv_RPHI = matrix_double(3, 3);
    if (matrix_dot_3_3(inv_RPSI, inv_RPHI, dot_inv_RPSI_inv_RPHI)) {
        snprintf(error_message, sizeof (error_message), "ERROR: in matrix_double_check_diagonal_non_zero_inverse()");
    }
    // dot
    MatrixDouble dot_T__inv_RPSI_inv_RPHI = matrix_double(3, 3);
    if (matrix_dot_3_3(T, dot_inv_RPSI_inv_RPHI, dot_T__inv_RPSI_inv_RPHI)) {
        snprintf(error_message, sizeof (error_message), "ERROR: in matrix_double_check_diagonal_non_zero_inverse()");
        return (-1);
    }
    //double G[3][3][3] = np.dot(T, np.dot(inv_RPSI, inv_RPHI));
    double THETA = atan2(dot_T__inv_RPSI_inv_RPHI[1][2], dot_T__inv_RPSI_inv_RPHI[1][1]) * RA2DE;
    if (THETA >= 360.0) {
        THETA -= 360.0;
    } else if (THETA < 0.0) {
        THETA += 360.0;
    }
    *pmajorAxisRotation = THETA;

    free_matrix_double(T, 3, 3);
    free_matrix_double(inv_RPSI, 3, 3);
    free_matrix_double(inv_RPHI, 3, 3);
    free_matrix_double(dot_inv_RPSI_inv_RPHI, 3, 3);
    free_matrix_double(dot_T__inv_RPSI_inv_RPHI, 3, 3);

    return (0);

}


