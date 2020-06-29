/***************************************************************************
 * matrix_statistics.h:
 *
 * TODO: add doc
 *
 * Written by Anthony Lomax
 *   ALomax Scientific www.alomax.net
 *
 * modified: 2010.12.16
 ***************************************************************************/

#define CONFIDENCE_LEVEL 68.0   // 68% confidence level used throughout

/** The relevant scaling factors for 68% are DELTA_CHI_SQR_68_2 and DELTA_CHI_SQR_68_3 in matrix_statistics/matrix_statistics.h
 * These are explained in Presss et al, Num Rec, 2nd ed, sec 15.6,
 * Chapter 3. Error Estimates in:  https://pubs.usgs.gov/of/1999/ofr-99-0023/HYPOELLIPSE_Full_Manual.pdf,
 * and also tabulated at https://en.wikipedia.org/wiki/Chi-squared_distribution
 * I think the values would be DELTA_CHI_SQR_68_3 = 6.25, and DELTA_CHI_SQR_68_2 = 4.61 for 90% confidence.  TODO: confirm this!
*/

// 2D ellipse

typedef struct {
    double az1, len1;   // semi-minor axis km
    double len2;   // semi-major axis km
} Ellipse2D;
#define DELTA_CHI_SQR_68_2 2.30    // value for 68% conf (see Num Rec, 2nd ed, sec 15.6, table)


// 3D ellipsoid

typedef struct {
    double az1, dip1, len1;   // semi-minor axis km
    double az2, dip2, len2;   // semi-intermediate axis km
    double len3;   // semi-major axis km
    double az3, dip3;   // 20150601 AJL - semi-major axis az and dip added to support conversion to QuakeML Tait-Bryan representation
} Ellipsoid3D;
#define DELTA_CHI_SQR_68_3 3.53    // value for 68% conf (see Num Rec, 2nd ed, sec 15.6, table)



char *get_matrix_statistics_error_mesage();
Vect3D CalcExpectationSamples(float*, int);
Vect3D CalcExpectationSamplesWeighted(float* fdata, int nSamples);
Vect3D CalcExpectationSamplesGlobal(float* fdata, int nSamples, double xReference);
Vect3D CalcExpectationSamplesGlobalWeighted(float* fdata, int nSamples, double xReference);
Mtrx3D CalcCovarianceSamplesRect(float* fdata, int nSamples, Vect3D* pexpect);
Mtrx3D CalcCovarianceSamplesGlobal(float* fdata, int nSamples, Vect3D* pexpect);
Mtrx3D CalcCovarianceSamplesGlobalWeighted(float* fdata, int nSamples, Vect3D* pexpect);
Ellipsoid3D CalcErrorEllipsoid(Mtrx3D *, double);
Ellipse2D CalcHorizontalErrorEllipse(Mtrx3D *pcov, double del_chi_2);
void ellipsiod2Axes(Ellipsoid3D *, Vect3D *, Vect3D *, Vect3D *);
void nllEllipsiod2XMLConfidenceEllipsoid(Ellipsoid3D *pellipsoid,
        double* psemiMajorAxisLength, double* pmajorAxisPlunge, double* pmajorAxisAzimuth,
        double* psemiIntermediateAxisLength, double* pintermediateAxisPlunge, double* pintermediateAxisAzimuth,
        double* psemiMinorAxisLength);
int nllEllipsiod2QMLConfidenceEllipsoid(Ellipsoid3D *pellipsoid,
        double* psemiMajorAxisLength,
        double* psemiMinorAxisLength,
        double* psemiIntermediateAxisLength,
        double* pmajorAxisAzimuth,
        double* pmajorAxisPlunge,
        double* pmajorAxisRotation);
