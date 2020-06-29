/*
 * Written by Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
 *
 */

void clean_SingularValueDecomposition();
void SingularValueDecomposition(MatrixDouble A_matrix, int nrows, int ncolumns);
double svd_cond();
MatrixDouble svd_getS();
VectorDouble svd_getSingularValues();
MatrixDouble svd_getU();
MatrixDouble svd_getV();
double svd_norm2();
int svd_rank();