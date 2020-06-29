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

#define PI 3.14159265359
#define TWOPI 6.28318530718
#define DE2RA 0.01745329252
#define RA2DE 57.2957795129
#define ERAD 6378.135
#define ERADM 6378135.0
#define AVG_ERAD 6371.0
#define FLATTENING 1.0/298.26	// Earth flattening (WGS '72)
#define EPS 0.000000000005
#define KM2MI 0.621371
#define KM2DEG (90.0/10000.0)
#define DEG2KM (10000.0/90.0)

double GCDistance(double lat1, double lon1, double lat2, double lon2);

double ApproxDistance(double lat1, double lon1, double lat2, double lon2);

double EllipsoidDistance(double lat1, double lon1, double lat2, double lon2);

double GCAzimuth(double lat1, double lon1, double lat2, double lon2);

double EllipsoidAzimuth(double lat1, double lon1, double lat2, double lon2);

double ApproxAzimuth(double lat1, double lon1, double lat2, double lon2);

