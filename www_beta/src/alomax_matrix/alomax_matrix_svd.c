/*
 * Converted from Java to C by Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
 *
 *
 */

#include <stdlib.h>
//#include <stdio.h>
//#include <string.h>
#include <math.h>
//#include <float.h>
#include "alomax_matrix.h"
#include "alomax_matrix_svd.h"


#define Math_min(A,B) ((B)<(A)?(B):(A))
#define Math_max(A,B) ((B)>(A)?(B):(A))

/*
Copyright � 1999 CERN - European Organization for Nuclear Research.
Permission to use, copy, modify, distribute and sell this software and its documentation for any purpose
is hereby granted without fee, provided that the above copyright notice appear in all copies and
that both that copyright notice and this permission notice appear in supporting documentation.
CERN makes no representations about the suitability of this software for any purpose.
It is provided "as is" without expressed or implied warranty.
 */

/**
 * Returns sqrt(a^2 + b^2) without under/overflow.
 */
static double Algebra_hypot(double a, double b) {
    double r;
    if (fabs(a) > fabs(b)) {
        r = b / a;
        r = fabs(a) * sqrt(1 + r * r);
    } else if (b != 0) {
        r = a / b;
        r = fabs(b) * sqrt(1 + r * r);
    } else {
        r = 0.0;
    }
    return r;
}



/**
For an <tt>m x n</tt> matrix <tt>A</tt> with <tt>m >= n</tt>, the singular value decomposition is
an <tt>m x n</tt> orthogonal matrix <tt>U</tt>, an <tt>n x n</tt> diagonal matrix <tt>S</tt>, and
an <tt>n x n</tt> orthogonal matrix <tt>V</tt> so that <tt>A = U*S*V'</tt>.
<P>
The singular values, <tt>sigma[k] = S[k][k]</tt>, are ordered so that
<tt>sigma[0] >= sigma[1] >= ... >= sigma[n-1]</tt>.
<P>
The singular value decomposition always exists, so the constructor will
never fail.  The matrix condition number and the effective numerical
rank can be computed from this decomposition.
 */

/** Arrays for internal storage of U and V.
@serial internal storage of U.
@serial internal storage of V.
 */
static MatrixDouble U_matrix = NULL;
static MatrixDouble V_matrix = NULL;
static MatrixDouble S_matrix = NULL;

/** Array for internal storage of singular values.
@serial internal storage of singular values.
 */
static VectorDouble singular_values = NULL;

/** Row and column dimensions.
@serial row dimension.
@serial column dimension.
 */
static int num_rows, num_columns;

/**
Constructs and returns a new singular value decomposition object;
The decomposed matrices can be retrieved via instance methods of the returned decomposition object.
@param A    A rectangular matrix.
@return     A decomposition object to access <tt>U</tt>, <tt>S</tt> and <tt>V</tt>.
@throws IllegalArgumentException if <tt>A.rows() < A.columns()</tt>.
 */
void SingularValueDecomposition(MatrixDouble A_matrix_orig, int nrows, int ncolumns) {

    int i, j, k;
    //Property.DEFAULT.checkRectangular(Arg);


    // Derived from LINPACK code.
    // Initialize.
    num_rows = nrows;
    num_columns = ncolumns;

    // make local copy of original A matrix
    MatrixDouble A_matrix = matrix_double(num_rows, num_columns);
    for (i = 0; i < num_rows; i++) {
        for (j = 0; j < num_columns; j++) {
            A_matrix[i][j] = A_matrix_orig[i][j];
        }
    }

    clean_SingularValueDecomposition();

    int nu = Math_min(num_rows, num_columns);
    singular_values = vector_double(Math_min(num_rows + 1, num_columns));
    U_matrix = matrix_double(num_rows, nu);
    V_matrix = matrix_double(num_columns, num_columns);
    double *e = calloc(num_columns, sizeof (double));
    double *work = calloc(num_rows, sizeof (double));

    int wantu = 1;
    int wantv = 1;

    // Reduce A to bidiagonal form, storing the diagonal elements
    // in s and the super-diagonal elements in e.

    int nct = Math_min(num_rows - 1, num_columns);
    int nrt = Math_max(0, Math_min(num_columns - 2, num_rows));
    for (k = 0; k < Math_max(nct, nrt); k++) {
        if (k < nct) {

            // Compute the transformation for the k-th column and
            // place the k-th diagonal in s[k].
            // Compute 2-norm of k-th column without under/overflow.
            singular_values[k] = 0;
            for (i = k; i < num_rows; i++) {
                singular_values[k] = Algebra_hypot(singular_values[k], A_matrix[i][k]);
            }
            if (singular_values[k] != 0.0) {
                if (A_matrix[k][k] < 0.0) {
                    singular_values[k] = -singular_values[k];
                }
                for (i = k; i < num_rows; i++) {
                    A_matrix[i][k] /= singular_values[k];
                }
                A_matrix[k][k] += 1.0;
            }
            singular_values[k] = -singular_values[k];
        }
        for (j = k + 1; j < num_columns; j++) {
            if ((k < nct) & (singular_values[k] != 0.0)) {

                // Apply the transformation.

                double t = 0;
                for (i = k; i < num_rows; i++) {
                    t += A_matrix[i][k] * A_matrix[i][j];
                }
                t = -t / A_matrix[k][k];
                for (i = k; i < num_rows; i++) {
                    A_matrix[i][j] += t * A_matrix[i][k];
                }
            }

            // Place the k-th row of A into e for the
            // subsequent calculation of the row transformation.

            e[j] = A_matrix[k][j];
        }
        if (wantu & (k < nct)) {

            // Place the transformation in U for subsequent back
            // multiplication.

            for (i = k; i < num_rows; i++) {
                U_matrix[i][k] = A_matrix[i][k];
            }
        }
        if (k < nrt) {

            // Compute the k-th row transformation and place the
            // k-th super-diagonal in e[k].
            // Compute 2-norm without under/overflow.
            e[k] = 0;
            for (i = k + 1; i < num_columns; i++) {
                e[k] = Algebra_hypot(e[k], e[i]);
            }
            if (e[k] != 0.0) {
                if (e[k + 1] < 0.0) {
                    e[k] = -e[k];
                }
                for (i = k + 1; i < num_columns; i++) {
                    e[i] /= e[k];
                }
                e[k + 1] += 1.0;
            }
            e[k] = -e[k];
            if ((k + 1 < num_rows) & (e[k] != 0.0)) {

                // Apply the transformation.

                for (i = k + 1; i < num_rows; i++) {
                    work[i] = 0.0;
                }
                for (j = k + 1; j < num_columns; j++) {
                    for (i = k + 1; i < num_rows; i++) {
                        work[i] += e[j] * A_matrix[i][j];
                    }
                }
                for (j = k + 1; j < num_columns; j++) {
                    double t = -e[j] / e[k + 1];
                    for (i = k + 1; i < num_rows; i++) {
                        A_matrix[i][j] += t * work[i];
                    }
                }
            }
            if (wantv) {

                // Place the transformation in V for subsequent
                // back multiplication.

                for (i = k + 1; i < num_columns; i++) {
                    V_matrix[i][k] = e[i];
                }
            }
        }
    }

    // Set up the final bidiagonal matrix or order p.

    int p = Math_min(num_columns, num_rows + 1);
    if (nct < num_columns) {
        singular_values[nct] = A_matrix[nct][nct];
    }
    if (num_rows < p) {
        singular_values[p - 1] = 0.0;
    }
    if (nrt + 1 < p) {
        e[nrt] = A_matrix[nrt][p - 1];
    }
    e[p - 1] = 0.0;

    // If required, generate U.

    if (wantu) {
        for (j = nct; j < nu; j++) {
            for (i = 0; i < num_rows; i++) {
                U_matrix[i][j] = 0.0;
            }
            U_matrix[j][j] = 1.0;
        }
        for (k = nct - 1; k >= 0; k--) {
            if (singular_values[k] != 0.0) {
                for (j = k + 1; j < nu; j++) {
                    double t = 0;
                    for (i = k; i < num_rows; i++) {
                        t += U_matrix[i][k] * U_matrix[i][j];
                    }
                    t = -t / U_matrix[k][k];
                    for (i = k; i < num_rows; i++) {
                        U_matrix[i][j] += t * U_matrix[i][k];
                    }
                }
                for (i = k; i < num_rows; i++) {
                    U_matrix[i][k] = -U_matrix[i][k];
                }
                U_matrix[k][k] = 1.0 + U_matrix[k][k];
                for (i = 0; i < k - 1; i++) {
                    U_matrix[i][k] = 0.0;
                }
            } else {
                for (i = 0; i < num_rows; i++) {
                    U_matrix[i][k] = 0.0;
                }
                U_matrix[k][k] = 1.0;
            }
        }
    }

    // If required, generate V.

    if (wantv) {
        for (k = num_columns - 1; k >= 0; k--) {
            if ((k < nrt) & (e[k] != 0.0)) {
                for (j = k + 1; j < nu; j++) {
                    double t = 0;
                    for (i = k + 1; i < num_columns; i++) {
                        t += V_matrix[i][k] * V_matrix[i][j];
                    }
                    t = -t / V_matrix[k + 1][k];
                    for (i = k + 1; i < num_columns; i++) {
                        V_matrix[i][j] += t * V_matrix[i][k];
                    }
                }
            }
            for (i = 0; i < num_columns; i++) {
                V_matrix[i][k] = 0.0;
            }
            V_matrix[k][k] = 1.0;
        }
    }

    // Main iteration loop for the singular values.

    int pp = p - 1;
    int iter = 0;
    double eps = pow(2.0, -52.0);
    while (p > 0) {
        int k, kase;

        // Here is where a test for too many iterations would go.

        // This section of the program inspects for
        // negligible elements in the s and e arrays.  On
        // completion the variables kase and k are set as follows.

        // kase = 1     if s(p) and e[k-1] are negligible and k<p
        // kase = 2     if s(k) is negligible and k<p
        // kase = 3     if e[k-1] is negligible, k<p, and
        //              s(k), ..., s(p) are not negligible (qr step).
        // kase = 4     if e(p-1) is negligible (convergence).

        for (k = p - 2; k >= -1; k--) {
            if (k == -1) {
                break;
            }
            if (fabs(e[k]) <= eps * (fabs(singular_values[k]) + fabs(singular_values[k + 1]))) {
                e[k] = 0.0;
                break;
            }
        }
        if (k == p - 2) {
            kase = 4;
        } else {
            int ks;
            for (ks = p - 1; ks >= k; ks--) {
                if (ks == k) {
                    break;
                }
                double t = (ks != p ? fabs(e[ks]) : 0.) +
                        (ks != k + 1 ? fabs(e[ks - 1]) : 0.);
                if (fabs(singular_values[ks]) <= eps * t) {
                    singular_values[ks] = 0.0;
                    break;
                }
            }
            if (ks == k) {
                kase = 3;
            } else if (ks == p - 1) {
                kase = 1;
            } else {
                kase = 2;
                k = ks;
            }
        }
        k++;

        // Perform the task indicated by kase.

        switch (kase) {

                // Deflate negligible s(p).

            case 1:
            {
                double f = e[p - 2];
                e[p - 2] = 0.0;
                for (j = p - 2; j >= k; j--) {
                    double t = Algebra_hypot(singular_values[j], f);
                    double cs = singular_values[j] / t;
                    double sn = f / t;
                    singular_values[j] = t;
                    if (j != k) {
                        f = -sn * e[j - 1];
                        e[j - 1] = cs * e[j - 1];
                    }
                    if (wantv) {
                        for (i = 0; i < num_columns; i++) {
                            t = cs * V_matrix[i][j] + sn * V_matrix[i][p - 1];
                            V_matrix[i][p - 1] = -sn * V_matrix[i][j] + cs * V_matrix[i][p - 1];
                            V_matrix[i][j] = t;
                        }
                    }
                }
            }
                break;

                // Split at negligible s(k).

            case 2:
            {
                double f = e[k - 1];
                e[k - 1] = 0.0;
                for (j = k; j < p; j++) {
                    double t = Algebra_hypot(singular_values[j], f);
                    double cs = singular_values[j] / t;
                    double sn = f / t;
                    singular_values[j] = t;
                    f = -sn * e[j];
                    e[j] = cs * e[j];
                    if (wantu) {
                        for (i = 0; i < num_rows; i++) {
                            t = cs * U_matrix[i][j] + sn * U_matrix[i][k - 1];
                            U_matrix[i][k - 1] = -sn * U_matrix[i][j] + cs * U_matrix[i][k - 1];
                            U_matrix[i][j] = t;
                        }
                    }
                }
            }
                break;

                // Perform one qr step.

            case 3:
            {

                // Calculate the shift.

                double scale = Math_max(Math_max(Math_max(Math_max(
                        fabs(singular_values[p - 1]), fabs(singular_values[p - 2])), fabs(e[p - 2])),
                        fabs(singular_values[k])), fabs(e[k]));
                double sp = singular_values[p - 1] / scale;
                double spm1 = singular_values[p - 2] / scale;
                double epm1 = e[p - 2] / scale;
                double sk = singular_values[k] / scale;
                double ek = e[k] / scale;
                double b = ((spm1 + sp)*(spm1 - sp) + epm1 * epm1) / 2.0;
                double c = (sp * epm1)*(sp * epm1);
                double shift = 0.0;
                if ((b != 0.0) | (c != 0.0)) {
                    shift = sqrt(b * b + c);
                    if (b < 0.0) {
                        shift = -shift;
                    }
                    shift = c / (b + shift);
                }
                double f = (sk + sp)*(sk - sp) + shift;
                double g = sk*ek;

                // Chase zeros.

                for (j = k; j < p - 1; j++) {
                    double t = Algebra_hypot(f, g);
                    double cs = f / t;
                    double sn = g / t;
                    if (j != k) {
                        e[j - 1] = t;
                    }
                    f = cs * singular_values[j] + sn * e[j];
                    e[j] = cs * e[j] - sn * singular_values[j];
                    g = sn * singular_values[j + 1];
                    singular_values[j + 1] = cs * singular_values[j + 1];
                    if (wantv) {
                        for (i = 0; i < num_columns; i++) {
                            t = cs * V_matrix[i][j] + sn * V_matrix[i][j + 1];
                            V_matrix[i][j + 1] = -sn * V_matrix[i][j] + cs * V_matrix[i][j + 1];
                            V_matrix[i][j] = t;
                        }
                    }
                    t = Algebra_hypot(f, g);
                    cs = f / t;
                    sn = g / t;
                    singular_values[j] = t;
                    f = cs * e[j] + sn * singular_values[j + 1];
                    singular_values[j + 1] = -sn * e[j] + cs * singular_values[j + 1];
                    g = sn * e[j + 1];
                    e[j + 1] = cs * e[j + 1];
                    if (wantu && (j < num_rows - 1)) {
                        for (i = 0; i < num_rows; i++) {
                            t = cs * U_matrix[i][j] + sn * U_matrix[i][j + 1];
                            U_matrix[i][j + 1] = -sn * U_matrix[i][j] + cs * U_matrix[i][j + 1];
                            U_matrix[i][j] = t;
                        }
                    }
                }
                e[p - 2] = f;
                iter = iter + 1;
            }
                break;

                // Convergence.

            case 4:
            {

                // Make the singular values positive.

                if (singular_values[k] <= 0.0) {
                    singular_values[k] = (singular_values[k] < 0.0 ? -singular_values[k] : 0.0);
                    if (wantv) {
                        for (i = 0; i <= pp; i++) {
                            V_matrix[i][k] = -V_matrix[i][k];
                        }
                    }
                }

                // Order the singular values.

                while (k < pp) {
                    if (singular_values[k] >= singular_values[k + 1]) {
                        break;
                    }
                    double t = singular_values[k];
                    singular_values[k] = singular_values[k + 1];
                    singular_values[k + 1] = t;
                    if (wantv && (k < num_columns - 1)) {
                        for (i = 0; i < num_columns; i++) {
                            t = V_matrix[i][k + 1];
                            V_matrix[i][k + 1] = V_matrix[i][k];
                            V_matrix[i][k] = t;
                        }
                    }
                    if (wantu && (k < num_rows - 1)) {
                        for (i = 0; i < num_rows; i++) {
                            t = U_matrix[i][k + 1];
                            U_matrix[i][k + 1] = U_matrix[i][k];
                            U_matrix[i][k] = t;
                        }
                    }
                    k++;
                }
                iter = 0;
                p--;
            }
                break;
        }
    }


    // clean up
    free(e);
    free(work);
    free_matrix_double(A_matrix, num_rows, nu);

}

/**
 * Clean up - free allocated matrices and vectors
 */
void clean_SingularValueDecomposition() {
    if (singular_values != NULL) {
        free_vector_double(singular_values);
        singular_values = NULL;
    }
    if (U_matrix != NULL) {
        int nu = Math_min(num_rows, num_columns);
        free_matrix_double(U_matrix, num_rows, nu);
        U_matrix = NULL;
    }
    if (V_matrix != NULL) {
        free_matrix_double(V_matrix, num_columns, num_columns);
        V_matrix = NULL;
    }
}

/**
Returns the two norm condition number, which is <tt>max(S) / min(S)</tt>.
 */
double svd_cond() {
    return singular_values[0] / singular_values[Math_min(num_rows, num_columns) - 1];
}

/**
Returns the diagonal matrix of singular values.
@return     S
 */
MatrixDouble svd_getS() {

    int i, j;

    if (S_matrix != NULL)
        free_matrix_double(S_matrix, num_columns, num_columns);
    S_matrix = matrix_double(num_columns, num_columns);
    for (i = 0; i < num_columns; i++) {
        for (j = 0; j < num_columns; j++) {
            S_matrix[i][j] = 0.0;
        }
        S_matrix[i][i] = singular_values[i];
    }
    return (S_matrix);
}

/**
Returns the diagonal of <tt>S</tt>, which is a one-dimensional array of singular values
@return     diagonal of <tt>S</tt>.
 */
VectorDouble svd_getSingularValues() {
    return singular_values;
}

/**
Returns the left singular vectors <tt>U</tt>.
@return     <tt>U</tt>
 */
MatrixDouble svd_getU() {
    //return new DoubleMatrix2D(U,m,Math_min(m+1,n));
    //return DoubleFactory2D.dense.make(U).viewPart(0, 0, num_rows, Math_min(num_rows + 1, num_columns));
    return (U_matrix); // TODO: should be partial matrix, see lines above
}

/**
Returns the right singular vectors <tt>V</tt>.
@return     <tt>V</tt>
 */
MatrixDouble svd_getV() {
    return (V_matrix);
}

/**
Returns the two norm, which is <tt>max(S)</tt>.
 */
double svd_norm2() {
    return singular_values[0];
}

/**
Returns the effective numerical matrix rank, which is the number of nonnegligible singular values.
 */
int svd_rank() {
    double eps = pow(2.0, -52.0);
    double tol = Math_max(num_rows, num_columns) * singular_values[0] * eps;
    int r = 0;
    int i;
    for (i = 0; i < Math_min(num_rows + 1, num_columns); i++) {
        if (singular_values[i] > tol) {
            r++;
        }
    }
    return r;
}




