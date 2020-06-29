/*
 * File:   test_ellipsoid.c
 * Author: anthony
 *
 * Created on 2014.11.06
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../vector/vector.h"
#include "../geometry/geometry.h"
#include "../alomax_matrix/alomax_matrix.h"
#include "matrix_statistics.h"

// EllAz1  288.007 Dip1  84.5079 Len1  57.5551 Az2  273.534 Dip2  -5.3188 Len2  84.7635 Len3  1.984275e+02

void usage() {
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "test_ellipsoid EllAz1  288.007 Dip1  84.5079 Len1  57.5551 Az2  273.534 Dip2  -5.3188 Len2  84.7635 Len3  1.984275e+02\n");
    fprintf(stderr, "test_ellipsoid CovXX 2063.45 XY 583.753 XZ 85.5223 YY 11110.7 YZ -248.964 ZZ 953.632\n");
    fprintf(stderr, "test_ellipsoid EllAz1  288.007 Dip1  84.5079 Len1  57.5551 Az2  273.534 Dip2  -5.3188 Len2  84.7635 Len3  1.984275e+02\n");
    fprintf(stderr, "test_ellipsoid\n");
}

int assert_almost_equal(double value, double target, int ndecimal, char* cname) {

    if (fabs(value - target) / fabs((value + target) / 2.0)> pow(10.0, -ndecimal)) {
        fprintf(stderr, "Failed: %s: %f != %f\n", cname, value, target);
        return (-1);
    }

    return (1);
}

/*
 *
 */
int main(int argc, char** argv) {

    Ellipsoid3D ellipsoid;

    int do_test = 0;
    if (argc < 2) {
        do_test = 1;
    } else if (strcmp(argv[1], "EllAz1") == 0) {
        if (argc < 15) {
            usage();
            exit(2);
        }
        int n = 2;
        ellipsoid.az1 = atof(argv[n]);
        n += 2;
        ellipsoid.dip1 = atof(argv[n]);
        n += 2;
        ellipsoid.len1 = atof(argv[n]);
        n += 2;
        ellipsoid.az2 = atof(argv[n]);
        n += 2;
        ellipsoid.dip2 = atof(argv[n]);
        n += 2;
        ellipsoid.len2 = atof(argv[n]);
        n += 2;
        ellipsoid.az3 = 0.0;
        ellipsoid.dip3 = 0.0;
        ellipsoid.len3 = atof(argv[n]);
    } else if (strcmp(argv[1], "CovXX") == 0) {
        if (argc < 13) {
            usage();
            exit(2);
        }
        Mtrx3D cov;
        int n = 2;
        cov.xx = atof(argv[n]);
        n += 2;
        cov.xy = cov.yx = atof(argv[n]);
        n += 2;
        cov.xz = cov.zx = atof(argv[n]);
        n += 2;
        cov.yy = atof(argv[n]);
        n += 2;
        cov.yz = cov.zy = atof(argv[n]);
        n += 2;
        cov.zz = atof(argv[n]);
        ellipsoid = CalcErrorEllipsoid(&cov, DELTA_CHI_SQR_68_3);
    }

    double semiMajorAxisLength, majorAxisPlunge, majorAxisAzimuth,
            semiIntermediateAxisLength,
            semiMinorAxisLength, majorAxisRotation;

    if (!do_test) {
        fprintf(stdout, "\nNLL Confidence Ellipsoid:\n");
        fprintf(stdout, "   ellipsoid.az1= %.3g\n", ellipsoid.az1);
        fprintf(stdout, "   ellipsoid.dip1= %.3g\n", ellipsoid.dip1);
        fprintf(stdout, "   ellipsoid.len1= %.3g\n", ellipsoid.len1);
        fprintf(stdout, "   ellipsoid.az2= %.3g\n", ellipsoid.az2);
        fprintf(stdout, "   ellipsoid.dip2= %.3g\n", ellipsoid.dip2);
        fprintf(stdout, "   ellipsoid.len2= %.3g\n", ellipsoid.len2);
        fprintf(stdout, "   ellipsoid.az3= %.3g\n", ellipsoid.az3);
        fprintf(stdout, "   ellipsoid.dip3= %.3g\n", ellipsoid.dip3);
        fprintf(stdout, "   ellipsoid.len3= %.3g\n", ellipsoid.len3);


        /*
        double intermediateAxisPlunge, intermediateAxisAzimuth;
        nllEllipsiod2XMLConfidenceEllipsoid(&ellipsoid,
                &semiMajorAxisLength, &majorAxisPlunge, &majorAxisAzimuth,
                &semiIntermediateAxisLength, &intermediateAxisPlunge, &intermediateAxisAzimuth,
                &semiMinorAxisLength);
        fprintf(stdout, "\nXMLConfidenceEllipsoid:\n");
        fprintf(stdout, "   semiMajorAxisLength= %.3g\n", semiMajorAxisLength);
        fprintf(stdout, "   majorAxisPlunge= %.3g\n", majorAxisPlunge);
        fprintf(stdout, "   majorAxisAzimuth= %.3g\n", majorAxisAzimuth);
        fprintf(stdout, "   semiIntermediateAxisLength= %.3g\n", semiIntermediateAxisLength);
        fprintf(stdout, "   intermediateAxisPlunge= %.3g\n", intermediateAxisPlunge);
        fprintf(stdout, "   intermediateAxisAzimuth= %.3g\n", intermediateAxisAzimuth);
        fprintf(stdout, "   semiMinorAxisLength= %.3g\n", semiMinorAxisLength);
         */

        nllEllipsiod2QMLConfidenceEllipsoid(&ellipsoid,
                &semiMajorAxisLength,
                &semiMinorAxisLength,
                &semiIntermediateAxisLength,
                &majorAxisAzimuth,
                &majorAxisPlunge,
                &majorAxisRotation);
        fprintf(stdout, "\nQMLConfidenceEllipsoid:\n");
        fprintf(stdout, "   semiMajorAxisLength= %.3g\n", semiMajorAxisLength);
        fprintf(stdout, "   semiMinorAxisLength= %.3g\n", semiMinorAxisLength);
        fprintf(stdout, "   semiIntermediateAxisLength= %.3g\n", semiIntermediateAxisLength);
        fprintf(stdout, "   majorAxisPlunge= %.3g\n", majorAxisPlunge);
        fprintf(stdout, "   majorAxisAzimuth= %.3g\n", majorAxisAzimuth);
        fprintf(stdout, "   majorAxisRotation= %.3g\n", majorAxisRotation);

    } else { // do_test (from libcomcat_ellipse.py)
        //CONSTANTS FOR TESTING
        double ALEN = 16000.0;
        double BLEN = 6100.0;
        double CLEN = 10900.0;
        double APLUNGE = 6.0;
        double AAZIMUTH = 139.0;
        double AROT = 88.9075;
        ellipsoid.az3 = 139.0;
        ellipsoid.dip3 = 6.0;
        ellipsoid.len3 = 16000 / 1000.;
        ellipsoid.az1 = 310.0;
        ellipsoid.dip1 = 83.0;
        ellipsoid.len1 = 6100 / 1000.;
        ellipsoid.az2 = 49.0;
        ellipsoid.dip2 = 1.0;
        ellipsoid.len2 = 10900 / 1000.;

        nllEllipsiod2QMLConfidenceEllipsoid(&ellipsoid,
                &semiMajorAxisLength,
                &semiMinorAxisLength,
                &semiIntermediateAxisLength,
                &majorAxisAzimuth,
                &majorAxisPlunge,
                &majorAxisRotation);

        assert_almost_equal(semiMajorAxisLength, ALEN / 1000.0, 3, "semiMajorAxisLength");
        assert_almost_equal(semiMinorAxisLength, BLEN / 1000.0, 3, "semiMinorAxisLength");
        assert_almost_equal(semiIntermediateAxisLength, CLEN / 1000.0, 3, "semiIntermediateAxisLength");
        assert_almost_equal(majorAxisPlunge, APLUNGE, 3, "majorAxisPlunge");
        assert_almost_equal(majorAxisAzimuth, AAZIMUTH, 3, "majorAxisAzimuth");
        assert_almost_equal(majorAxisRotation, AROT, 3, "majorAxisRotation");
    }

    return (EXIT_SUCCESS);
}

