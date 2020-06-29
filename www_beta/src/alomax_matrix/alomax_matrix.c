/*
 * Written by Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
 *
 * Released into the Public Domain
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "alomax_matrix.h"
#include "alomax_matrix_svd.h"

//#define EPSILON_BIG  (1.0e-5)
//#define EPSILON_BIG  (1.0e-3)
//#define EPSILON_BIG  (1.0e-6)   // 20100601 AJL (pred analy)
#define EPSILON_BIG  (1.0e-3)   // 20100617 AJL (pred analy)
#define EPSILON  FLT_MIN

static char error_message[4096];

/** function to print error and return last error message */
char *get_matrix_error_mesage() {
    return (error_message);
}

MatrixDouble matrix_double(int nrow, int ncol) {
    MatrixDouble mtx;

    mtx = (MatrixDouble) calloc(nrow, sizeof (double*));
    if (mtx == NULL) {
        snprintf(error_message, sizeof (error_message), "ERROR: in matrix_double(): allocating rows.");
        return (NULL);
    }

    int n, m;
    for (n = 0; n < nrow; n++) {
        mtx[n] = (double *) calloc(ncol, sizeof (double));
        if (mtx[n] == NULL) {
            snprintf(error_message, sizeof (error_message), "ERROR: in matrix_double():  allocating columns.");
            free_matrix_double(mtx, nrow, ncol);
            return (NULL);
        }
        for (m = 0; m < ncol; m++) {
            mtx[n][m] = 0.0;
        }
    }
    return (mtx);

}

void free_matrix_double(MatrixDouble mtx, int nrow, int ncol) {

    if ((mtx) == NULL) return;

    int n;
    for (n = nrow - 1; n >= 0; n--) {
        if (mtx[n] != NULL)
            free(mtx[n]);
    }
    free((mtx));

}

/** function  to display double matrix */

void display_matrix_double(char* name, MatrixDouble matrix, int num_rows, int num_cols) {
    int nrow, ncol;

    fprintf(stdout, "\n%s Matrix: %d rows X %d columns\n",
            name, num_rows, num_cols);
    for (nrow = 0; nrow < num_rows; nrow++) {
        for (ncol = 0; ncol < num_cols; ncol++) {
            if (ncol == nrow)
                fprintf(stdout, "\\ ");
            fprintf(stdout, "%g ", matrix[nrow][ncol]);
            if (ncol == nrow)
                fprintf(stdout, "\\ ");
        }
        fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n");

}

VectorDouble vector_double(int nsize) {
    VectorDouble vect;

    vect = (VectorDouble) calloc(nsize, sizeof (double));
    if (vect == NULL) {
        snprintf(error_message, sizeof (error_message), "ERROR: in vector_double(): allocating elements.");
        return (NULL);
    }

    return (vect);

}

void free_vector_double(VectorDouble vect) {

    if ((vect) == NULL) return;

    free((vect));

}

/** function  to display double matrix */

void display_vector_double(char* name, VectorDouble vect, int nsize) {
    int n;

    fprintf(stdout, "\n%s Vector: %d elements\n", name, nsize);
    for (n = 0; n < nsize; n++) {
        fprintf(stdout, "%g ", vect[n]);
    }
    fprintf(stdout, "\n");

}

/*
 * Simple Gauss-Jordan elimination
 * Modified from Python written by Jarno Elonen <elonen@iki.fi>, april 2005, released into the Public Domain
 *
 * The following ultra-compact Python function performs in-place Gaussian elimination for given matrix,
 * putting it into the Reduced Row Echelon Form. It can be used to solve linear equation systems or to invert a matrix.
 *
 * Returns 0 if successful, -1 if 'm' is singular.
 */

int gauss_jordan(MatrixDouble dmtx, int num_rows, int num_cols) {

    //(h, w) = (len(m), len(m[0]))

    double c, temp;
    int ix, iy, iy2;
    int i;

    int maxrow;
    /*{
        fprintf(stderr, "DEBUG:  original state of matrix:\n");
        int k, l;
        for (k = 0; k < num_rows; k++) {
            for (l = 0; l < num_cols; l++)
                fprintf(stderr, "\t  %g", dmtx[k][l]);
            fprintf(stderr, "\n");
        }
    }*/

    for (iy = 0; iy < num_rows; iy++) {
        maxrow = iy;
        for (iy2 = iy + 1; iy2 < num_rows; iy2++) { // Find max pivot
            if (fabs(dmtx[iy2][iy]) > fabs(dmtx[maxrow][iy]))
                maxrow = iy2;
        }
        //(m[iy], m[maxrow]) = (m[maxrow], m[iy]);
        for (i = 0; i < num_cols; i++) {
            temp = dmtx[maxrow][i];
            dmtx[maxrow][i] = dmtx[iy][i];
            dmtx[iy][i] = temp;
        }
        if (fabs(dmtx[iy][iy]) <= EPSILON) { // Singular?
            snprintf(error_message, sizeof (error_message), "ERROR: in gauss_jordan(): singular matrix: element %d %d with value %f.", iy, iy, dmtx[iy][iy]);
            /*{
                fprintf(stderr, "DEBUG:  current state of matrix:\n");
                int k, l;
                for (k = 0; k < num_rows; k++) {
                    for (l = 0; l < num_cols; l++)
                        fprintf(stderr, "\t  %g", dmtx[k][l]);
                    fprintf(stderr, "\n");
                }
            }*/
            return (-1);
        }
        for (iy2 = iy + 1; iy2 < num_rows; iy2++) { // Eliminate column y
            c = dmtx[iy2][iy] / dmtx[iy][iy];
            for (ix = iy; ix < num_cols; ix++) {
                dmtx[iy2][ix] -= dmtx[iy][ix] * c;
            }
        }
    }
    for (iy = num_rows - 1; iy >= 0; iy--) { // Backsubstitute
        c = dmtx[iy][iy];
        for (iy2 = 0; iy2 < iy; iy2++) {
            for (ix = num_cols - 1; ix >= iy; ix--) {
                dmtx[iy2][ix] -= dmtx[iy][ix] * dmtx[iy2][iy] / c;
            }
        }
        dmtx[iy][iy] /= c;
        for (ix = num_rows; ix < num_cols; ix++) { // Normalize row y
            dmtx[iy][ix] /= c;
        }
    }

    return (0);

}

/*
 *   return the inv of the matrix M
 */

int matrix_double_inverse(MatrixDouble dmtx, int num_rows, int num_cols) {

    //clone the matrix and append the identity matrix
    MatrixDouble augmented_mtx;
    if ((augmented_mtx = matrix_double(num_rows, 2 * num_cols)) < 0) {
        snprintf(error_message, sizeof (error_message), "ERROR: in matrix_double_inverse(): allocating matrix: augmented_mtx.");
        return (-1);
    }

    // [int(i==j) for j in range_M] is nothing but the i(th row of the identity matrix
    int nrow, ncol;
    for (nrow = 0; nrow < num_rows; nrow++) {
        for (ncol = 0; ncol < num_cols; ncol++) {
            augmented_mtx[nrow][ncol] = dmtx[nrow][ncol];
        }
        for (; ncol < 2 * num_cols; ncol++) {
            augmented_mtx[nrow][ncol] = nrow == ncol - num_cols ? 1.0 : 0.0;
        }
    }
    // apply Gauss-Jordan elimination
    int istat = gauss_jordan(augmented_mtx, num_rows, 2 * num_cols);
    // extract the appended inverse matrix
    for (nrow = 0; nrow < num_rows; nrow++) {
        for (ncol = 0; ncol < num_cols; ncol++) {
            dmtx[nrow][ncol] = augmented_mtx[nrow][num_cols + ncol];
        }
    }

    free_matrix_double(augmented_mtx, num_rows, 2 * num_cols);

    if (istat < 0) {
        char error_message2[4096];
        sprintf(error_message2, "ERROR: in matrix_double_inverse(): %s", error_message);
        strcpy(error_message, error_message2);
        return (-1);
    }


    return (istat);

}

/*
 *   return the inv of the square matrix M for the rows/columns which have non-zero diagonal
 */

int matrix_double_check_diagonal_non_zero_inverse(MatrixDouble mtx_original, int i_original_size, int verify_inverse, int verbose) {

    int i;

    // estimate tolerance
    double tolerance = 0.0;
    for (i = 0; i < i_original_size; i++) {
        if (fabs(mtx_original[i][i]) > tolerance)
            tolerance = fabs(mtx_original[i][i]);
    }
    //tolerance /= 10000.0;
    tolerance /= 1.0e30;


    // check for zeros on diagonal
    int inon_zero[i_original_size];
    int i_checked_size = 0;
    for (i = 0; i < i_original_size; i++) {
        if (fabs(mtx_original[i][i]) < tolerance) {
            inon_zero[i] = 0;
        } else {
            inon_zero[i] = 1;
            i_checked_size++;
        }
    }

    if (i_checked_size < 1) {
        snprintf(error_message, sizeof (error_message), "ERROR: in matrix_double_check_diagonal_non_zero_inverse(): no non-zero diagonal elements.");
        return (-1);
    }

    // initialized checked matrix
    MatrixDouble mtx_checked = mtx_original;
    // load original matrix to checked matrix
    if (i_checked_size != i_original_size) {
        if ((mtx_checked = matrix_double(i_checked_size, i_checked_size)) < 0) {
            snprintf(error_message, sizeof (error_message), "ERROR: in matrix_double_check_diagonal_non_zero_inverse(): allocating matrix: mtx_checked.");
            return (-1);
        }
        int nrow_original, ncol_original;
        int nrow_checked = 0;
        for (nrow_original = 0; nrow_original < i_original_size; nrow_original++) {
            if (inon_zero[nrow_original]) {
                int ncol_checked = 0;
                for (ncol_original = 0; ncol_original < i_original_size; ncol_original++) {
                    if (inon_zero[ncol_original]) {
                        mtx_checked[nrow_checked][ncol_checked] = mtx_original[nrow_original][ncol_original];
                        ncol_checked++;
                    }
                }
                nrow_checked++;
            }
        }
    }
    /*{
        fprintf(stderr, "DEBUG:  mtx_checked matrix:\n");
        int k, l;
        for (k = 0; k < i_checked_size; k++) {
            for (l = 0; l < i_checked_size; l++)
                fprintf(stderr, "\t  %g", mtx_checked[k][l]);
            fprintf(stderr, "\n");
        }
    }*/

    // save matrix for later display and testing of results
    MatrixDouble mtx_checked_before_inverse = NULL;
    if (verify_inverse) {
        mtx_checked_before_inverse = matrix_double(i_checked_size, i_checked_size);
        int l, k;
        for (l = 0; l < i_checked_size; l++) {
            for (k = 0; k < i_checked_size; k++) mtx_checked_before_inverse[k][l] = mtx_checked[k][l];
        }
    }
    /*{
        if (i_checked_size == 2)
            fprintf(stderr, "DEBUG:  mtx_checked det: %g\n", mtx_checked[0][0] * mtx_checked[1][1] - mtx_checked[0][1] * mtx_checked[1][0]);
    }*/

    // invert
    if (matrix_double_inverse(mtx_checked, i_checked_size, i_checked_size) < 0) {
        char error_message2[4096];
        sprintf(error_message2, "ERROR: in matrix_double_check_diagonal_non_zero_inverse(): %s", error_message);
        strcpy(error_message, error_message2);
        if (verify_inverse)
            free_matrix_double(mtx_checked_before_inverse, i_checked_size, i_checked_size);
        return (-1);
    }

    // check inverse
    int ierror = 0;
    if (verify_inverse) {
        if (square_inverse_not_ok(mtx_checked, mtx_checked_before_inverse, i_checked_size, verbose)) {
            snprintf(error_message, sizeof (error_message), "ERROR: in matrix_double_check_diagonal_non_zero_inverse(): square_inverse_not_ok.");
            ierror = -1;
        }
        free_matrix_double(mtx_checked_before_inverse, i_checked_size, i_checked_size);
    }

    // load inverted checked matrix to original matrix
    if (i_checked_size != i_original_size) {
        int nrow_original, ncol_original;
        int nrow_checked = 0;
        for (nrow_original = 0; nrow_original < i_original_size; nrow_original++) {
            if (inon_zero[nrow_original]) {
                int ncol_checked = 0;
                for (ncol_original = 0; ncol_original < i_original_size; ncol_original++) {
                    if (inon_zero[ncol_original]) {
                        mtx_original[nrow_original][ncol_original] = mtx_checked[nrow_checked][ncol_checked];
                        ncol_checked++;
                    }
                }
                nrow_checked++;
            }
        }
        free_matrix_double(mtx_checked, i_checked_size, i_checked_size);
    }

    return (ierror);

}

/*
 * check if a specified inverse matrix is the inverse of a specified original matrix
 *   return 0 if inverse is ok
 *   return 1 if inverse in not ok
 *   return -1 on error
 */

int square_inverse_not_ok(MatrixDouble inverse_mtrx, MatrixDouble original_mtx, int nsize, int verbose) {

    MatrixDouble test_matrix;
    if ((test_matrix = matrix_double(nsize, nsize)) == NULL) {
        if (verbose)
            fprintf(stderr, "ERROR: %s\n", get_matrix_error_mesage());
        return (-1);
    }

    // check inverse
    int mtx_error = 0;
    int k, l, j;
    for (k = 0; k < nsize; k++) {
        for (l = 0; l < nsize; l++) {
            test_matrix[k][l] = 0.0;
            for (j = 0; j < nsize; j++) {
                test_matrix[k][l] += (inverse_mtrx[k][j] * original_mtx[j][l]);
            }
            if ((k == l && fabs(test_matrix[k][l]) < (1.0 - EPSILON_BIG))
                    || (k != l && fabs(test_matrix[k][l]) > EPSILON_BIG))
                mtx_error = 1;
        }
    }
    if (mtx_error) {
        if (verbose) {
            fprintf(stderr, "ERROR: inverse_mtrx times original_mtx not identity matrix:\n");
            for (k = 0; k < nsize; k++) {
                for (l = 0; l < nsize; l++)
                    fprintf(stderr, "\t  %g", test_matrix[k][l]);
                fprintf(stderr, "\n");
            }
        }
        free_matrix_double(test_matrix, nsize, nsize);
        return (1);
    }

    free_matrix_double(test_matrix, nsize, nsize);
    return (0);

}

#undef EPSILON
#undef EPSILON_BIG

/**
 *  Helper function to apply singular value decomposition
 *
 * input:
 *    a_matrix - num_rows x num_cols matrix to which to apply svs
 *
 * output:
 *   s_vector - vector of singular values of size num_cols)
 *   v_matrix - orthogonal matrix of right singular vectors of size num_columns x num_columns
 *
 */

void svd_helper(MatrixDouble A_matrix, int num_rows, int num_cols, VectorDouble S_vector, MatrixDouble V_matrix) {

    SingularValueDecomposition(A_matrix, num_rows, num_cols);
    VectorDouble s_vector = svd_getSingularValues();
    MatrixDouble v_matrix = svd_getV();

    int i, j;
    for (i = 0; i < num_rows; i++) {
        for (j = 0; j < num_cols; j++) {
            V_matrix[i][j] = v_matrix[i][j];
        }
    }
    for (j = 0; j < num_cols; j++) {
        S_vector[j] = s_vector[j];
    }

    clean_SingularValueDecomposition();

}