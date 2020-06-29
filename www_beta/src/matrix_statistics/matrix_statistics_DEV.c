/***************************************************************************
 * matrix_statistics.c:
 *
 * TODO: add doc
 *
 * Written by Anthony Lomax
 *   ALomax Scientific www.alomax.net
 *
 * modified: 2010.12.16
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

/** function to calculate the covariance of a set of samples */

Mtrx3D CalcCovarianceSamplesRect(float* fdata, int nSamples, Vect3D* pexpect) {

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

/** function to calculate the covariance of a set of samples with longitude (x) wrap around of 360 deg
 *  if isGeographic=1 use long(deg)/lat(deg)/depth(km) coordinates
 *  if isGeographic=1 use x=0,360 deg wrap coordinate and y,z = linear coordinates
 */

Mtrx3D CalcCovarianceSamplesWrappedX(float* fdata, int nSamples, Vect3D* pexpect, int isGeographic) {

    int nsamp, ipos;

    float x, y, z, prob;
    Mtrx3D cov = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    double cos_lat = 1.0;
    double WORK_DEG2KM = 1.0;
    if (isGeographic) {
        cos_lat = cos(pexpect->y * DE2RA);
        WORK_DEG2KM = DEG2KM;
    }
    double xReference = pexpect->x;

    /* calculate covariance following eq. (6-12), T & V, 1982 */

    ipos = 0;
    for (nsamp = 0; nsamp < nSamples; nsamp++) {
        x = fdata[ipos++];
        if (x - xReference > 180.0)
            x -= 360.0;
        else if (x - xReference < -180.0)
            x += 360.0;
        x = x * WORK_DEG2KM * cos_lat;
        y = fdata[ipos++] * WORK_DEG2KM;
        z = fdata[ipos++];
        prob = fdata[ipos++];

        cov.xx += (double) (x * x);
        cov.xy += (double) (x * y);
        cov.xz += (double) (x * z);

        cov.yy += (double) (y * y);
        cov.yz += (double) (y * z);

        cov.zz += (double) (z * z);

    }

    cov.xx = cov.xx / (double) nSamples - pexpect->x * pexpect->x * WORK_DEG2KM * cos_lat * WORK_DEG2KM * cos_lat;
    cov.xy = cov.xy / (double) nSamples - pexpect->x * pexpect->y * WORK_DEG2KM * cos_lat * WORK_DEG2KM;
    cov.xz = cov.xz / (double) nSamples - pexpect->x * pexpect->z * WORK_DEG2KM * cos_lat;

    cov.yx = cov.xy;
    cov.yy = cov.yy / (double) nSamples - pexpect->y * pexpect->y * WORK_DEG2KM * WORK_DEG2KM;
    cov.yz = cov.yz / (double) nSamples - pexpect->y * pexpect->z * WORK_DEG2KM;

    cov.zx = cov.xz;
    cov.zy = cov.yz;
    cov.zz = cov.zz / (double) nSamples - pexpect->z * pexpect->z;


    return (cov);
}

/** function to calculate the covariance of a set of samples in long(deg)/lat(deg)/depth(km) coordinates */

Mtrx3D CalcCovarianceSamplesGlobal(float* fdata, int nSamples, Vect3D* pexpect) {
    return(CalcCovarianceSamplesWrappedX(fdata, nSamples, pexpect, 1));
}


/** function to calculate the covariance of a set of samples in long(deg)/lat(deg)/depth(km) coordinates */

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


    /* strike angles positive CCW from East = 0 */
    az1 = 90.0 - pellipsoid->az1;
    az2 = 90.0 - pellipsoid->az2;
    /* dip angles increasing downwards from horiz = 0 */
    dip1 = -pellipsoid->dip1;
    dip2 = -pellipsoid->dip2;

    /* get 3D vector axes */

    cosd1 = cos(DE2RA * pellipsoid->dip1);
    paxis1->x = cos(DE2RA * az1) * cosd1;
    paxis1->y = sin(DE2RA * az1) * cosd1;
    paxis1->z = -sin(DE2RA * dip1);

    cosd2 = cos(DE2RA * pellipsoid->dip2);
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

/** method to convert ellipsoid to an XML (pseudo-QuakeML) ConfidenceEllipsoid */

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

    double plunge = 0.0;
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


    plunge = 0.0;
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

