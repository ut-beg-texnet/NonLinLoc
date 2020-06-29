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

import java.io.Serializable;

/**
 * Class to hold a single slowness layer sample.
 * 
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * 
 * 
 * 
 * @author H. Philip Crotwell
 */
public class SlownessLayer implements Serializable {

    /** top slowness, top depth, bottom slowness, bottom depth */
    public SlownessLayer(double topP,
                         double topDepth,
                         double botP,
                         double botDepth) {
        Assert.isFalse(topDepth > botDepth, "topDepth > botDepth: " + topDepth
                + " > " + botDepth);
        Assert.isFalse(topDepth < 0.0 || Double.isNaN(topDepth)
                               || Double.isInfinite(topDepth),
                       "topDepth is not a number or is negative: " + topDepth);
        Assert.isFalse(botDepth < 0.0 || Double.isNaN(botDepth)
                               || Double.isInfinite(botDepth),
                       "botDepth is not a number or is negative: " + botDepth);
        this.topP = topP;
        this.topDepth = topDepth;
        this.botP = botP;
        this.botDepth = botDepth;
    }

    /**
     * Compute the slowness layer from a velocity layer.
     */
    public SlownessLayer(VelocityLayer vLayer,
                         boolean spherical,
                         double radiusOfEarth,
                         boolean isPWave) {
        Assert.isFalse(vLayer.getTopDepth() > vLayer.getBotDepth(),
                       "vLayer.topDepth > vLayer.botDepth :"
                               + vLayer.getTopDepth() + " "
                               + vLayer.getBotDepth());
        topDepth = vLayer.getTopDepth();
        botDepth = vLayer.getBotDepth();
        char waveType;
        if(isPWave) {
            waveType = 'P';
        } else {
            waveType = 'S';
        }
        try {
            if(spherical) {
                topP = (radiusOfEarth - getTopDepth())
                        / vLayer.evaluateAtTop(waveType);
                botP = (radiusOfEarth - getBotDepth())
                        / vLayer.evaluateAtBottom(waveType);
            } else {
                topP = 1.0 / vLayer.evaluateAtTop(waveType);
                botP = 1.0 / vLayer.evaluateAtBottom(waveType);
            }
            Assert.isFalse(Double.isNaN(getTopP()) || Double.isNaN(getBotP()),
                           "Slowness sample is NaN: topP=" + getTopP()
                                   + " botP=" + getBotP()+" depth "+topDepth+" to "+botDepth);
        } catch(NoSuchMatPropException e) {
            // Can't happen
            throw new RuntimeException(e);
        }
    }

    /**
     * Compute the slowness layer from a velocity layer. Since radiusOfEarth is
     * given we assume a spherical model.
     */
    public SlownessLayer(VelocityLayer vLayer,
                         boolean isPWave,
                         double radiusOfEarth) {
        this(vLayer, true, radiusOfEarth, isPWave);
    }

    /**
     * Compute the slowness layer from a velocity layer. Since radiusOfEarth is
     * not given we assume a flat model.
     */
    public SlownessLayer(VelocityLayer vLayer, boolean isPWave) {
        this(vLayer, false, 0.0, isPWave);
    }

    public double getTopP() {
        return topP;
    }

    public double getBotP() {
        return botP;
    }

    public double getTopDepth() {
        return topDepth;
    }

    public double getBotDepth() {
        return botDepth;
    }

    /** Is the layer a zero thickness layer, ie a total reflection? */
    public boolean isZeroThickness() {
        if(getTopDepth() == getBotDepth()) {
            return true;
        } else {
            return false;
        }
    }
    
    public boolean containsDepth(double depth) {
        return depth >= getTopDepth() && depth <= getBotDepth();
    }

    /**
     * Finds the slowness at the given depth. radiusOfEarth is needed as a
     * slowness layer doesn't have access to the slowness model. Note that this
     * method assumes a Bullen type of slowness interpolation, ie p(r) = a*r^b.
     * This will produce results consistent with a tau model that uses this
     * interpolant, but it may differ slightly from going directly to the
     * velocity model. Also, if the tau model is generated using another
     * interpolant, linear for instance, then the result may not be consistent
     * with the tau model.
     */
    public double evaluateAt_bullen(double depth, double radiusOfEarth)
            throws SlownessModelException {
        Assert.isFalse(getBotDepth() > radiusOfEarth,
                       "SlownessLayer.evaluateAt_bullen:"
                               + " radiusOfEarth="
                               + radiusOfEarth
                               + " is smaller than the maximum depth of this layer."
                               + " topDepth=" + getTopDepth() + " botDepth="
                               + getBotDepth());
        Assert.isFalse((getTopDepth() - depth) * (depth - getBotDepth()) < 0.0,
                       "SlownessLayer.evaluateAt_bullen:" + " depth=" + depth
                               + " is not contained within this layer."
                               + " topDepth=" + getTopDepth() + " botDepth="
                               + getBotDepth());
        if(depth == getTopDepth()) {
            return getTopP();
        } else if(depth == getBotDepth()) {
            return getBotP();
        } else {
            double B = Math.log(getTopP() / getBotP())
                    / Math.log((radiusOfEarth - getTopDepth())
                            / (radiusOfEarth - getBotDepth()));
            double ADenominator = Math.pow((radiusOfEarth - getTopDepth()), B);
            double A = getTopP() / ADenominator;
            double answer = A * Math.pow((radiusOfEarth - depth), B);
            if(answer < 0.0 || Double.isNaN(answer)
                    || Double.isInfinite(answer)) {
                // numerical instability in power law calculation???
                // try a linear interpolation if the layer is small ( <2 km)
                // or if denominator of A is infinity as we probably overflowed
                // the double in that case.
                if((getBotDepth() - getTopDepth()) < 2.0
                        || Double.isInfinite(ADenominator) || getBotP() == 0.0) {
                    double linear = (getBotP() - getTopP())
                            / (getBotDepth() - getTopDepth())
                            * (depth - getTopDepth()) + getTopP();
                    if(linear < 0.0 || Double.isNaN(linear)
                            || Double.isInfinite(linear)) {} else {
                        return linear;
                    }
                }
                throw new SlownessModelException("calculated slowness at depth="
                        + depth
                        + " is not a number or is negative: "
                        + answer
                        + "\n"
                        + this.toString()
                        + "\n A="
                        + A
                        + "   B="
                        + B
                        + " rad-dep="
                        + (radiusOfEarth - depth)
                        + "  a demon="
                        + ADenominator);
            }
            return answer;
        }
    }

    /**
     * Calculates the time and distance (in radians) increments accumulated by a
     * ray of spherical ray parameter p when passing through this layer. Note
     * that this gives 1/2 of the true range and time increments since there
     * will be both an up going and a downgoing path. Here we use the
     * Mohorovicic or Bullen law p=A*r^B
     * 
     * @exception SlownessModelException
     *                occurs if the calculated distance or time increments are
     *                negative or NaN, this indicates a bug in the code (and
     *                hopefully will never happen).
     */
    public TimeDist bullenRadialSlowness(double p, double radiusOfEarth)
            throws SlownessModelException {
        // To hold the return values.
        TimeDist timedist = new TimeDist(p);
        if(getTopDepth() == getBotDepth()) {
            timedist.distRadian = 0.0;
            timedist.time = 0.0;
            return timedist;
        }
        // only do bullen radial slowness if the layer is not too thin
        // here we use 1 micron = .000000001
        // just return 0 in this case
        if(getBotDepth() - getTopDepth() < .000000001) {
            return timedist;
        }
        double B = Math.log(getTopP() / getBotP())
                / Math.log((radiusOfEarth - getTopDepth())
                        / (radiusOfEarth - getBotDepth()));
        double sqrtTopTopMpp = Math.sqrt(getTopP() * getTopP() - p * p);
        double sqrtBotBotMpp = Math.sqrt(getBotP() * getBotP() - p * p);
        timedist.distRadian = (Math.atan2(p, sqrtBotBotMpp) - Math.atan2(p,
                                                                   sqrtTopTopMpp))
                / B;
        timedist.time = (sqrtTopTopMpp - sqrtBotBotMpp) / B;
        if(timedist.distRadian < 0.0 || timedist.time < 0.0
                || Double.isNaN(timedist.time) || Double.isNaN(timedist.distRadian)) {
            throw new SlownessModelException("timedist <0.0 or NaN: "
                    + "\n RayParam= " + p + "\n topDepth = " + getTopDepth()
                    + "\n botDepth = " + getBotDepth() + "\n dist="
                    + timedist.distRadian + "\n time=" + timedist.time + "\n topP = "
                    + getTopP() + "\n botP = " + getBotP() + "\n B = " + B
                    + " " + toString());
        }
        return timedist;
    }

    /**
     * Finds the depth for a ray parameter within this layer. Uses a Bullen
     * interpolant, Ar^B. Special case for botP == 0 or botDepth ==
     * radiusOfEarth as these cause div by 0, use linear interpolation in this
     * case.
     */
    public double bullenDepthFor(double rayParam, double radiusOfEarth)
            throws SlownessModelException {
        if((getTopP() - rayParam) * (rayParam - getBotP()) >= 0) {
            double tempDepth;
            // easy case for 0 thickness layer
            if(getTopDepth() == getBotDepth()) {
                return getBotDepth();
            }
            if (getTopP() == rayParam) {return getTopDepth();}
            if (getBotP() == rayParam) {return getBotDepth();}
            if(getBotP() != 0.0 && getBotDepth() != radiusOfEarth) {
                double B = Math.log(getTopP() / getBotP())
                        / Math.log((radiusOfEarth - getTopDepth())
                                / (radiusOfEarth - getBotDepth()));
                double A = getTopP()
                        / Math.pow((radiusOfEarth - getTopDepth()), B);
                tempDepth = radiusOfEarth
                        - Math.exp(1.0 / B * Math.log(rayParam / A));
                /*
                 * tempDepth = radiusOfEarth - Math.pow(rayParam/A, 1.0/B);
                 */
                // check for slightly outside layer due to rounding or numerical instability
                if (tempDepth < getTopDepth() && tempDepth > getTopDepth()-0.000001) {tempDepth = getTopDepth();}
                if (tempDepth > getBotDepth() && tempDepth < getBotDepth()+0.000001) {tempDepth = getBotDepth();}
                if(tempDepth < 0.0 || Double.isNaN(tempDepth)
                        || Double.isInfinite(tempDepth)
                        || tempDepth < getTopDepth()
                        || tempDepth > getBotDepth()) {
                    // numerical instability in power law calculation???
                    // try a linear interpolation if the layer is small ( <5
                    // km).
                    if((getBotDepth() - getTopDepth()) < 5.0) {
                        double linear = (getBotDepth() - getTopDepth())
                                / (getBotP() - getTopP())
                                * (rayParam - getTopP()) + getTopDepth();
                        if(linear < 0.0 || Double.isNaN(linear)
                                || Double.isInfinite(linear)) {} else {
                            return linear;
                        }
                    }
                    throw new SlownessModelException("claculated depth is outside layer, not a number or is negative: "
                            + tempDepth
                            + "\n"
                            + this
                            + "\n"
                            + A
                            + "  "
                            + B
                            + "\n" + rayParam);
                }
                // check for tempDepth just above top depth
                if(tempDepth < getTopDepth()
                        && (getTopDepth() - tempDepth) < 1e-10) {
                    return getTopDepth();
                }
                // check for tempDepth just below bottom depth
                if(tempDepth > getBotDepth()
                        && (tempDepth - getBotDepth()) < 1e-10) {
                    return getBotDepth();
                }
                return tempDepth;
            } else {
                // a special case for the center of the earth, since ar^b
                // might blow up at r=0
                if(getTopP() != getBotP()) {
                    return getBotDepth() + (rayParam - getBotP())
                            * (getTopDepth() - getBotDepth())
                            / (getTopP() - getBotP());
                } else {
                    // weird case, return botDepth???
                    return getBotDepth();
                }
            }
        } else {
            throw new SlownessModelException("Ray parameter = " + rayParam
                    + " is not contained within this slowness layer. topP="
                    + getTopP() + " botP=" + getBotP());
        }
    }

    /** returns a String description of this SlownessLayer. */
    public String toString() {
        // String desc = "top p "+ (float)topP +", topDepth " + (float)topDepth
        // +", bot p "+ (float)botP +", botDepth " + (float)botDepth;
        String desc = "top p " + getTopP() + ", topDepth " + getTopDepth()
                + ", bot p " + getBotP() + ", botDepth " + getBotDepth();
        return desc;
    }

    public boolean validate() throws SlownessModelException {
        if(Double.isNaN(getTopP()) || Double.isNaN(getTopDepth())
                || Double.isNaN(getBotP()) || Double.isNaN(getBotDepth())) {
            throw new SlownessModelException("Slowness layer has NaN values."
                    + "\n " + this);
        }
        if(getTopP() < 0.0 || getBotP() < 0.0) {
            throw new SlownessModelException("Slowness layer has negative slownesses. \n "
                    + this);
        }
        if(getTopDepth() > getBotDepth()) {
            throw new SlownessModelException("Slowness layer has negative thickness. \n"
                    + this);
        }
        return true;
    }

    /** Slowness at the top of the layer. */
    private double topP;

    /** Slowness at the bottom of the layer. */
    private double botP;

    /** Depth at the top of the layer. */
    private double topDepth;

    /** Depth at the bottom of the layer. */
    private double botDepth;
}
