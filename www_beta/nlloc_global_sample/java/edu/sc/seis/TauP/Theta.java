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
/**
 * Theta.java
 * 
 * 
 * Created: Mon Feb 15 13:48:06 1999
 * 
 * @author Philip Crotwell
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * 
 * 
 * 
 */
package edu.sc.seis.TauP;

public class Theta {

    protected double radians;

    protected double[] thetaAtX;

    protected double[] rayParams;

    public Theta(SeismicPhase phase, double radians) {
        this.radians = radians;
        int minRayParamIndex = phase.getMinRayParamIndex();
        int maxRayParamIndex = phase.getMaxRayParamIndex();
        rayParams = phase.getRayParams();
        thetaAtX = phase.getTau();
        TauModel tMod = phase.getTauModel();
        // change tau to theta, theta = tau + p*x0
        for(int i = 0; i < thetaAtX.length; i++) {
            thetaAtX[i] += rayParams[i] * radians;
        }
    }

    /**
     * Get the value of radians.
     * 
     * @return Value of radians.
     */
    public double getRadians() {
        return radians;
    }

    public double getMaxRayParam() {
        return rayParams[0];
    }

    public double getStepRayParam(double rayParam, double timeStep) {
        double thetaStart = getTheta(rayParam);
        // loop until we find a theta s.t. abs(thetaStart-theta) == timeStep
        // or we fall off the end of the array, ie ArrayIndexOutOfBounds
        boolean found = false;
        int i = getThetaIndex(rayParam);
        while(Math.abs(thetaAtX[i + 1] - thetaStart) <= timeStep) {
            i++;
        }
        double newTheta;
        if(thetaStart < thetaAtX[i + 1]) {
            newTheta = thetaStart + timeStep;
        } else {
            newTheta = thetaStart - timeStep;
        }
        // return the interpolated ray parameter.
        return linInterp(thetaAtX[i],
                         thetaAtX[i + 1],
                         rayParams[i],
                         rayParams[i + 1],
                         newTheta);
    }

    public double getTheta(double rayParam) {
        if(rayParam > rayParams[0]
                || rayParam < rayParams[rayParams.length - 1]) {
            throw new ArrayIndexOutOfBoundsException(rayParam
                    + " not in range " + rayParams[0] + " to "
                    + rayParams[rayParams.length - 1]);
        }
        int currentNum = getThetaIndex(rayParam);
        // find theta at given rayParam
        double thetaStart = linInterp(rayParams[currentNum],
                                      rayParams[currentNum + 1],
                                      thetaAtX[currentNum],
                                      thetaAtX[currentNum + 1],
                                      rayParam);
        return thetaStart;
    }

    protected int getThetaIndex(double rayParam) {
        // find index containing rayParam
        int tooSmallNum = 0;
        int tooLargeNum = rayParams.length - 1;
        int currentNum = 0;
        boolean found = false;
        while(!found) {
            currentNum = (int)Math.floor((tooSmallNum + tooLargeNum) / 2.0f);
            if(rayParams[currentNum] >= rayParam
                    && rayParams[currentNum + 1] < rayParam) {
                found = true;
            } else if(rayParams[currentNum] > rayParam) {
                tooSmallNum = currentNum;
            } else if(rayParams[currentNum] <= rayParam) {
                tooLargeNum = currentNum;
            }
        }
        return currentNum;
    }

    public static double linInterp(double xa,
                                   double xb,
                                   double ya,
                                   double yb,
                                   double x) {
        return (yb - ya) * (x - xa) / (xb - xa) + ya;
    }
} // Theta
