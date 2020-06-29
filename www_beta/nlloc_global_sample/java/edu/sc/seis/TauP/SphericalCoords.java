/*
 * The TauP Toolkit: Flexible Seismic Travel-Time and Raypath Utilities.
 * Copyright (C) 1998-2000 University of South Carolina
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place - Suite 330, Boston, MA 02111-1307, USA.
 * 
 * The current version can be found at <A
 * HREF="www.seis.sc.edu">http://www.seis.sc.edu</A>
 * 
 * Bug reports and comments should be directed to H. Philip Crotwell,
 * crotwell@seis.sc.edu or Tom Owens, owens@seis.sc.edu
 * 
 */
package edu.sc.seis.TauP;

/**
 * Utility class for spherical coordinate (lat-lon) transformations. Given lat,
 * lon, lat, lon you can find the distance or azimuth and given lat, lon,
 * distance, azimuth you can find the lat lon of the resultant point. Just uses
 * spherical relations, no ellpticity correction is applied.
 * 
 * See Appendix A of "Seismology and Plate Tectonics" by David Gubbins Cambridge
 * University Press, 1990
 * 
 * and Chapter 3 of "Plate Tectonics: How it Works" by Allan Cox and Robert
 * Brian Hart Blackwell Scientific Publications, 1986
 * 
 * @author H. Philip Crotwell
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * 
 * 
 * 
 */
public class SphericalCoords {

    protected static final double dtor = Math.PI / 180.0;

    protected static final double rtod = 180.0 / Math.PI;

    /** Calculates angular distance between two lat lon pairs. */
    public static double distance(double latA,
                                  double lonA,
                                  double latB,
                                  double lonB) {
        return rtod
                * Math.acos(Math.sin(latA * dtor) * Math.sin(latB * dtor)
                        + Math.cos(latA * dtor) * Math.cos(latB * dtor)
                        * Math.cos((lonB - lonA) * dtor));
    }

    /** Calculates azimuth between two lat lon pairs. */
    public static double azimuth(double latA,
                                 double lonA,
                                 double latB,
                                 double lonB) {
        double cosAzimuth = (Math.cos(latA * dtor) * Math.sin(latB * dtor) - Math.sin(latA
                * dtor)
                * Math.cos(latB * dtor) * Math.cos((lonB - lonA) * dtor))
                / Math.sin(distance(latA, lonA, latB, lonB) * dtor);
        double sinAzimuth = Math.cos(latB * dtor)
                * Math.sin((lonB - lonA) * dtor)
                / Math.sin(distance(latA, lonA, latB, lonB) * dtor);
        return rtod * Math.atan2(sinAzimuth, cosAzimuth);
    }

    /**
     * Find the rotation pole required to rotate the first lat lon pair to the
     * second. Just does a cross product.
     * 
     * @return a 3 element double array with the X, Y and Z components of the
     *          pole.
     */
    public static double[] rotationPole(double latA,
                                        double lonA,
                                        double latB,
                                        double lonB) {
        double[] pointA = new double[3];
        double[] pointB = new double[3];
        double[] pole = new double[3];
        double dToR = Math.PI / 180.0;
        pointA[0] = Math.cos(latA * dToR) * Math.cos(lonA * dToR);
        pointA[1] = Math.cos(latA * dToR) * Math.sin(lonA * dToR);
        pointA[2] = Math.sin(latA * dToR);
        pointB[0] = Math.cos(latB * dToR) * Math.cos(lonB * dToR);
        pointB[1] = Math.cos(latB * dToR) * Math.sin(lonB * dToR);
        pointB[2] = Math.sin(latB * dToR);
        pole[0] = pointA[1] * pointB[2] - pointA[2] * pointB[1];
        pole[1] = pointA[2] * pointB[0] - pointA[0] * pointB[2];
        pole[2] = pointA[0] * pointB[1] - pointA[1] * pointB[0];
        return pole;
    }

    /**
     * rotates a point about a pole by an angle.
     * 
     * @param pole
     *            is a 3 element double array with X, Y and Z components of the
     *            pole.
     * @return [lat, lon] in array.
     */
    public static double[] rotate(double latA,
                                  double lonA,
                                  double[] pole,
                                  double angleDeg) {
        double[][] R = new double[3][3]; /* rotation matrix. */
        double[] point = new double[3];
        double[] newPoint = new double[3];
        double rToDeg = 180.0 / Math.PI;
        double angle = angleDeg / rToDeg;
        R[0][0] = pole[0] * pole[0] * (1 - Math.cos(angle)) + Math.cos(angle);
        R[0][1] = pole[0] * pole[1] * (1 - Math.cos(angle)) - pole[2]
                * Math.sin(angle);
        R[0][2] = pole[0] * pole[2] * (1 - Math.cos(angle)) + pole[1]
                * Math.sin(angle);
        R[1][0] = pole[1] * pole[0] * (1 - Math.cos(angle)) + pole[2]
                * Math.sin(angle);
        R[1][1] = pole[1] * pole[1] * (1 - Math.cos(angle)) + Math.cos(angle);
        R[1][2] = pole[1] * pole[2] * (1 - Math.cos(angle)) - pole[0]
                * Math.sin(angle);
        R[2][0] = pole[2] * pole[0] * (1 - Math.cos(angle)) - pole[1]
                * Math.sin(angle);
        R[2][1] = pole[2] * pole[1] * (1 - Math.cos(angle)) + pole[0]
                * Math.sin(angle);
        R[2][2] = pole[2] * pole[2] * (1 - Math.cos(angle)) + Math.cos(angle);
        point[0] = Math.cos(latA / rToDeg) * Math.cos(lonA / rToDeg);
        point[1] = Math.cos(latA / rToDeg) * Math.sin(lonA / rToDeg);
        point[2] = Math.sin(latA / rToDeg);
        newPoint[0] = R[0][0] * point[0] + R[0][1] * point[1] + R[0][2]
                * point[2];
        newPoint[1] = R[1][0] * point[0] + R[1][1] * point[1] + R[1][2]
                * point[2];
        newPoint[2] = R[2][0] * point[0] + R[2][1] * point[1] + R[2][2]
                * point[2];
        double newLat = Math.asin(newPoint[2]) * 180.0 / Math.PI;
        double newLon = Math.atan2(newPoint[1], newPoint[0]) * 180.0 / Math.PI;
        newPoint = new double[2];
        newPoint[0] = newLat;
        newPoint[1] = newLon;
        return newPoint;
    }

    /**
     * Calculates the latitude of a point a given distance along a given azimuth
     * from a starting lat lon.
     */
    public static double latFor(double latA,
                                double lonA,
                                double distance,
                                double azimuth) {
        return rtod
                * Math.asin(Math.cos(azimuth * dtor)
                        * Math.sin(distance * dtor) * Math.cos(latA * dtor)
                        + Math.cos(distance * dtor) * Math.sin(latA * dtor));
    }

    /**
     * Calculates the longitude of a point a given distance along a given
     * azimuth from a starting lat lon.
     */
    public static double lonFor(double latA,
                                double lonA,
                                double distance,
                                double azimuth) {
        double tempLat = latFor(latA, lonA, distance, azimuth);
        double sinLon = Math.sin(azimuth * dtor) * Math.sin(distance * dtor)
                / Math.cos(tempLat * dtor);
        double cosLon = (Math.cos(distance * dtor) - Math.sin(latA * dtor)
                * Math.sin(tempLat * dtor))
                / (Math.cos(latA * dtor) * Math.cos(tempLat * dtor));
        // make sure answer (lon) is between -180 and 180
        double lon = lonA + rtod * Math.atan2(sinLon, cosLon);
        if(lon <= -180.0)
            lon += 360.0;
        if(lon > 180.0)
            lon -= 360.0;
        return lon;
    }

    public static void main(String args[]) {
        System.out.println(distance(0, 0, 0, 45) + "  " + azimuth(0, 0, 0, 45)
                + "   " + azimuth(0, 45, 0, 0));
        System.out.println(latFor(0, 0, 45, 90) + "   " + lonFor(0, 0, 45, 90));
        System.out.println("(35,42,36,43)  " + distance(35, 42, 36, 43) + "  "
                + azimuth(35, 42, 36, 43) + "   " + azimuth(36, 43, 35, 42));
        System.out.println(latFor(35, 42, distance(35, 42, 36, 43), azimuth(35,
                                                                            42,
                                                                            36,
                                                                            43))
                + "   "
                + lonFor(35, 42, distance(35, 42, 36, 43), azimuth(35,
                                                                   42,
                                                                   36,
                                                                   43)));
    }
}
