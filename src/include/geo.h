//	geo.h
//
// Adapted from:
//    Geographic Distance and Azimuth Calculations
//    by Andy McGovern
//    http://www.codeguru.com/algorithms/GeoCalc.html

//#ifndef GEO_CALCULATIONS_H_
//#define GEO_CALCULATIONS_H_

  //
  // some geo constants
  //

// 20171122 AJL  #define PI 3.14159265359
#define PI M_PI
// 20171122 AJL  #define TWOPI 6.28318530718
#define TWOPI (2.0*M_PI)
// 20171122 AJL  #define DE2RA 0.01745329252
// 20171122 AJL  #define RA2DE 57.2957795129
#define DE2RA (M_PI/180.0)
#define RA2DE (180.0/M_PI)
#define ERAD 6378.135       // WGS-72
#define ERADM 6378135.0     // WGS-72
// 20171122 AJL  #define AVG_ERAD 6371.0
#define AVG_ERAD 6371.0087714    // GMT "Sphere"
#define FLATTENING (1.0/298.26)	// Earth flattening (WGS-72)
#define EPS 0.000000000005
#define KM2MI 0.621371
// 20151106 AJL - changed km/deg scaling to be based on sphere with radius 6371, average Earth radius.
//#define KM2DEG (90.0/10000.0)
//#define DEG2KM (10000.0/90.0)
#define KM2DEG (180.0/(PI*AVG_ERAD))
#define DEG2KM (PI*AVG_ERAD/180.0)

double GCDistance(double lat1, double lon1, double lat2, double lon2);

double ApproxDistance(double lat1, double lon1, double lat2, double lon2);

double EllipsoidDistance(double lat1, double lon1, double lat2, double lon2);

double GCAzimuth(double lat1, double lon1, double lat2, double lon2);

double EllipsoidAzimuth(double lat1, double lon1, double lat2, double lon2);

double ApproxAzimuth(double lat1, double lon1, double lat2, double lon2);

void PointAtGCDistanceAzimuth(double lat1, double lon1, double dist, double az, double* lat2, double* lon2);
