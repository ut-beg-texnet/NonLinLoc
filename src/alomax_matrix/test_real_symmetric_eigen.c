/*
 * File:   test_real_symmetric_eigen.c
 * Author: anthony
 *
 * Created on February 17, 2010, 7:03 AM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "alomax_matrix.h"
#include "eigv.h"

void usage() {
    fprintf(stderr, "USAGE: test_real_symmetric_eigen a11 a21 a22 a31 a32 a33\n");
    fprintf(stderr, "Compare to: http://www.akiti.ca/Eig3RSSolv.html\n");
}

/*
 *
 */
int main(int argc, char** argv) {

    if (argc < 7) {
        usage();
        exit(2);
    }

    int n = 0;
    double a11 = atof(argv[++n]);
    double a21 = atof(argv[++n]);
    double a22 = atof(argv[++n]);
    double a31 = atof(argv[++n]);
    double a32 = atof(argv[++n]);
    double a33 = atof(argv[++n]);

    printf("The input [A] Matrix components follow:\n");
    printf("a11 %f\n", a11);
    printf("a21 %f a22 %f\n", a21, a22);
    printf("a31 %f a32 %f a33 %f\n", a31, a32, a33);

    MatrixDouble A_matrix = matrix_double(3, 3);
    A_matrix[0][0] = a11;
    A_matrix[0][1] = a21;
    A_matrix[0][2] = a31;
    A_matrix[1][0] = a21;
    A_matrix[1][1] = a22;
    A_matrix[1][2] = a32;
    A_matrix[2][0] = a31;
    A_matrix[2][1] = a32;
    A_matrix[2][2] = a33;

    // find eigenvalues and eigenvectors
    VectorDouble S_vector = vector_double(3);
    MatrixDouble V_matrix = matrix_double(3, 3);

    // S_vector - vector of eigenvalues of size isize, in ascending order
    // V_matrix - orthogonal matrix of right singular vectors of size isize x isize
    real_symmetric_eigen_helper(A_matrix, 3, S_vector, V_matrix);


    printf("The eigenvalues follow:\n");
    printf("%f\n", S_vector[0]);
     printf("%f\n", S_vector[1]);
    printf("%f\n", S_vector[2]);

    return (EXIT_SUCCESS);
}

