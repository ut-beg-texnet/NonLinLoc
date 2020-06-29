/*
 * The TauP Toolkit: Flexible Seismic Travel-Time and Raypath Utilities.
 * Copyright (C) 1998-2000 University of South Carolina This program is free
 * software; you can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version. This program
 * is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU General Public License for more details. You
 * should have received a copy of the GNU General Public License along with this
 * program; if not, write to the Free Software Foundation, Inc., 59 Temple Place -
 * Suite 330, Boston, MA 02111-1307, USA. The current version can be found at <A
 * HREF="www.seis.sc.edu">http://www.seis.sc.edu </A> Bug reports and comments
 * should be directed to H. Philip Crotwell, crotwell@seis.sc.edu or Tom Owens,
 * owens@seis.sc.edu
 */
package edu.sc.seis.TauP;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.Serializable;

/**
 * provides storage and methods for distance, time and tau increments for a
 * branch. A branch is a group of layers bounded by discontinuities or reversals
 * in slowness gradient.
 * 
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * @author H. Philip Crotwell
 */
public class TauBranch implements Serializable, Cloneable {

    /** Turns on debugging output. */
    transient public boolean DEBUG = false;

    /** The type of wave for this branch, P or S. */
    protected boolean isPWave;

    /** The minimum depth of this layer. */
    private double topDepth;

    /** The maximum depth of this layer. */
    private double botDepth;

    /**
     * The maximum ray parameter that can penetrate into this branch. Time,
     * distance, and tau increments are undefined (0.0) for ray parameters
     * larger than this.
     */
    private double maxRayParam;

    /**
     * The minimum ray parameter that is turned, but not reflected, in this
     * branch. If no rays turn, then it is set equal to maxRayParam.
     */
    private double minTurnRayParam;

    /**
     * The minimum ray parameter that is turned or critically reflected in this
     * branch. If no rays turn, then it is set equal to maxRayParam.
     */
    private double minRayParam;

    /**
     * Holds distance (radians) evaluated at the ith ray parameter for this
     * branch.
     */
    protected double[] dist;

    /**
     * Holds time (seconds) evaluated at the ith ray parameter for this branch.
     */
    protected double[] time;

    /** Holds tau evaluated at the ith ray parameter for this branch. */
    protected double[] tau;

    // Constructors --------------------------------------------------------
    public TauBranch(boolean isPWave,
                     double topDepth,
                     double botDepth,
                     double maxRayParam,
                     double minTurnRayParam,
                     double minRayParam,
                     double[] dist,
                     double[] time,
                     double[] tau) {
        this.isPWave = isPWave;
        this.topDepth = topDepth;
        this.botDepth = botDepth;
        this.maxRayParam = maxRayParam;
        this.minTurnRayParam = minTurnRayParam;
        this.minRayParam = minRayParam;
        this.dist = dist;
        this.time = time;
        this.tau = tau;
    }
    
    public TauBranch(double topDepth, double botDepth, boolean isPWave) {
        this.topDepth = topDepth;
        this.botDepth = botDepth;
        this.isPWave = isPWave;
    }

    // Methods -------------------------------------------------------------
    // Accessor methods
    /** @return the minimum (top) depth of this layer. */
    public double getTopDepth() {
        return topDepth;
    }

    /** @return the maximum (bottom) depth of this layer. */
    public double getBotDepth() {
        return botDepth;
    }

    /**
     * @return the maximum ray parameter that can penetrate into this branch.
     *          Time, distance, and tau increments are undefined, set to 0.0,
     *          for ray parameters larger than this.
     */
    public double getMaxRayParam() {
        return maxRayParam;
    }

    /**
     * @return the minimum ray parameter that is turned, but not reflected, in
     *          this branch.
     */
    public double getMinTurnRayParam() {
        return minTurnRayParam;
    }

    /**
     * @return the minimum ray parameter that is turned or reflected in this
     *          branch.
     */
    public double getMinRayParam() {
        return minRayParam;
    }

    /**
     * @return an array, cloned, containing distance (radians) evaluated at the
     *          i_th ray parameter for this branch.
     */
    public double[] getDist() {
        return (double[])dist.clone();
    }

    /**
     * @return the distance (radians) evaluated at the i_th ray parameter for
     *          this branch.
     */
    public double getDist(int i) {
        return dist[i];
    }

    /**
     * @return an array, cloned, containing time (seconds) evaluated at the
     *          i_th ray parameter for this branch.
     */
    public double[] getTime() {
        return (double[])time.clone();
    }

    /**
     * @return the time (seconds) evaluated at the i_th ray parameter for this
     *          branch.
     */
    public double getTime(int i) {
        return time[i];
    }

    /**
     * @return an array, cloned, containing tau (seconds) evaluated at the i_th
     *          ray parameter for this branch.
     */
    public double[] getTau() {
        return (double[])tau.clone();
    }

    /**
     * @return tau (seconds) evaluated at the i_th ray parameter for this
     *          branch.
     */
    public double getTau(int i) {
        return tau[i];
    }

    // normal methods
    /**
     * Calculates tau for this branch, between slowness layers topLayerNum and
     * botLayerNum, inclusive.
     * 
     * @exception NoSuchLayerException
     *                if a needed slowness layer cannot be found.
     * @exception SlownessModelException
     *                if there is a problem with the slowness model
     * @exception TauModelException
     *                if the slownessmodel and taumodel are not compatible
     */
    public void createBranch(SlownessModel sMod,
                             double minPSoFar,
                             double rayParams[]) throws NoSuchLayerException,
            SlownessModelException, TauModelException {
        TimeDist timeDist;
        double p;
        int topLayerNum = sMod.layerNumberBelow(getTopDepth(), isPWave);
        int botLayerNum = sMod.layerNumberAbove(getBotDepth(), isPWave);
        SlownessLayer topSLayer = sMod.getSlownessLayer(topLayerNum,
                                                             isPWave);
        SlownessLayer botSLayer = sMod.getSlownessLayer(botLayerNum,
                                                             isPWave);
        if(topSLayer.getTopDepth() != getTopDepth()
                || botSLayer.getBotDepth() != getBotDepth()) {
            if(topSLayer.getTopDepth() != getTopDepth()
                    && Math.abs(topSLayer.getTopDepth() - getTopDepth()) < 0.000001) {
                // really close, so just move top
                System.err.println("Changing topDepth: " + "\ntopDepth: "
                        + getTopDepth() + " " + topSLayer.getTopDepth()
                        + "\nbotDepth: " + getBotDepth() + " "
                        + botSLayer.getBotDepth());
                topDepth = topSLayer.getTopDepth();
            } else if(botSLayer.getBotDepth() != getBotDepth()
                    && Math.abs(botSLayer.getBotDepth() - getBotDepth()) < 0.000001) {
                // really close, so just move bottom
                System.err.println("Changing botDepth: " + "\ntopDepth: "
                        + getTopDepth() + " " + topSLayer.getTopDepth()
                        + "\nbotDepth: " + getBotDepth() + " "
                        + botSLayer.getBotDepth());
                botDepth = botSLayer.getBotDepth();
            } else {
                // bad match, throw exception
                throw new TauModelException("createBranch: TauBranch not compatible with slowness sampling for "
                        + "SP".charAt(isPWave ? 1 : 0)
                        + ":"
                        + "\ntopDepth: "
                        + getTopDepth()
                        + " "
                        + topSLayer.getTopDepth()
                        + "\nbotDepth: "
                        + getBotDepth()
                        + " "
                        + botSLayer.getBotDepth());
            }
        }
        /*
         * Here we set minTurnRayParam to be the ray parameter that turns within
         * the layer, not including total reflections off of the bottom.
         * maxRayParam is the largest ray parameter that can penetrate this
         * branch. minRayParam is the minimum ray parameter that turns or is
         * totally reflected in this branch.
         */
        maxRayParam = minPSoFar;
        minTurnRayParam = sMod.getMinTurnRayParam(getBotDepth(), isPWave);
        minRayParam = sMod.getMinRayParam(getBotDepth(), isPWave);
        tau = new double[rayParams.length];
        dist = new double[rayParams.length];
        time = new double[rayParams.length];
        for(int rayNum = 0; rayNum < rayParams.length; rayNum++) {
            p = rayParams[rayNum];
            timeDist = calcTimeDist(sMod, topLayerNum, botLayerNum, p);
            dist[rayNum] = timeDist.distRadian;
            time[rayNum] = timeDist.time;
            tau[rayNum] = time[rayNum] - p * dist[rayNum];
            if(DEBUG && (rayNum % ((int)rayParams.length / 10) == 0)) {
                System.out.print(rayNum + ", ");
            }
        }
    }

    /**
     * calculates the time and distance increments for the given ray parameter.
     * The topDepth and botDepth must be correct as they determine the bounds on
     * the integration/summing.
     * 
     * @exception NoSuchLayerException
     *                if topLayerNum or botLayerNum are not in the slowness
     *                model.
     * @exception SlownessModelException
     *                if the ray with ray parameter p turns within a layer
     *                instead of at the bottom.
     */
    public TimeDist calcTimeDist(SlownessModel sMod,
                                 int topLayerNum,
                                 int botLayerNum,
                                 double p) throws NoSuchLayerException,
            SlownessModelException {
        int layerNum;
        TimeDist timeDist = new TimeDist(p);
        SlownessLayer layer;
        if(p <= getMaxRayParam()) {
            layerNum = topLayerNum;
            layer = sMod.getSlownessLayer(layerNum, isPWave);
            while(layerNum <= botLayerNum && p <= layer.getTopP()
                    && p <= layer.getBotP()) {
                timeDist.add(sMod.layerTimeDist(p, layerNum, isPWave));
                layerNum++;
                if(layerNum <= botLayerNum) {
                    layer = sMod.getSlownessLayer(layerNum, isPWave);
                }
            }
            if((layer.getTopP() - p) * (p - layer.getBotP()) > 0) {
                throw new SlownessModelException("Ray turns in the middle of this"
                        + " layer. layerNum = "
                        + layerNum
                        + " sphericalRayParam " + p + " layer =" + layer);
            }
        }
        return timeDist;
    }

    /**
     * Inserts the distance, time, and tau increment for the slowness sample
     * given to the branch. This is used for making the depth correction to a
     * tau model for a non-surface source.
     * 
     * @throws TauModelException
     *             if the tau branch is not compatable with the slowness
     *             sampling
     * @see edu.sc.seis.TauP.TauModel.depthCorrect(double)
     */
    protected void insert(double rayParam, SlownessModel sMod, int index)
            throws NoSuchLayerException, SlownessModelException,
            TauModelException {
        int topLayerNum = sMod.layerNumberBelow(getTopDepth(), isPWave);
        int botLayerNum = sMod.layerNumberAbove(getBotDepth(), isPWave);
        SlownessLayer topSLayer = sMod.getSlownessLayer(topLayerNum,
                                                             isPWave);
        SlownessLayer botSLayer = sMod.getSlownessLayer(botLayerNum,
                                                             isPWave);
        if(topSLayer.getTopDepth() != getTopDepth()
                || botSLayer.getBotDepth() != getBotDepth()) {
            throw new TauModelException("insert: TauBranch depths not compatible with slowness sampling:"
                    + "\ntopDepth: "
                    + getTopDepth()
                    + " "
                    + topSLayer.getTopDepth()
                    + "\nbotDepth: "
                    + getBotDepth()
                    + " " + botSLayer.getBotDepth());
        }
        TimeDist td = new TimeDist(rayParam, 0.0, 0.0);
        TimeDist temptd;
        if(topSLayer.getBotP() >= rayParam && topSLayer.getTopP() >= rayParam) {
            for(int i = topLayerNum; i <= botLayerNum; i++) {
                if(sMod.getSlownessLayer(i, isPWave).getBotP() < rayParam) {
                    // so we don't sum below the turning depth
                    break;
                } else {
                    temptd = sMod.layerTimeDist(rayParam, i, isPWave);
                    td.distRadian += temptd.distRadian;
                    td.time += temptd.time;
                }
            }
        }
        shiftBranch(index);
        dist[index] = td.distRadian;
        time[index] = td.time;
        tau[index] = td.time - rayParam * td.distRadian;
    }

    /**
     * generates a new tau branch by "subtracting" the given tau branch from
     * this tau branch. The given tau branch is assumed to by the upper part of
     * this branch.
     * 
     * @argument indexP specifies where a new ray coresponding to a P wave
     *           sample has been added, it is -1 if no ray parameter has been
     *           added to topBranch.
     * @arguement indexS is similar to indexP except for a S wave sample. Note
     *            that although the ray parameters for indexP and indexS were
     *            for the P and S waves that turned at the source depth, both
     *            ray parameters need to be added to both P and S branches.
     */
    protected TauBranch difference(TauBranch topBranch,
                                   int indexP,
                                   int indexS,
                                   SlownessModel sMod,
                                   double minPSoFar,
                                   double rayParams[])
            throws NoSuchLayerException, SlownessModelException,
            TauModelException {
        if(topBranch.getTopDepth() != getTopDepth()
                || topBranch.getBotDepth() > getBotDepth()) {
            if(topBranch.getTopDepth() != getTopDepth()
                    && Math.abs(topBranch.getTopDepth() - getTopDepth()) < 0.000001) {
                // really close, so just move top
                topDepth = topBranch.getTopDepth();
            } else {
                // bad match, throw exception
                throw new TauModelException("difference: TauBranch not compatible with slowness sampling:"
                        + "\ntopDepth: "
                        + getTopDepth()
                        + " "
                        + topBranch.getTopDepth()
                        + "\nbotDepth: "
                        + getBotDepth() + " " + topBranch.getBotDepth());
            }
        }
        if(topBranch.isPWave != isPWave) {
            throw new TauModelException("Can't difference branches: "
                    + "topBranch.topDepth=" + topBranch.getTopDepth()
                    + " topDepth=" + getTopDepth() + " topBranch.botDepth="
                    + topBranch.getBotDepth() + " botDepth=" + getBotDepth()
                    + " waveTypes:" + topBranch.isPWave + " " + isPWave);
        }
        // find the top and bottom slowness layers of bottom half
        int topLayerNum = sMod.layerNumberBelow(topBranch.getBotDepth(),
                                                isPWave);
        int botLayerNum = sMod.layerNumberBelow(getBotDepth(), isPWave); // branch include zero thickness layers at bottom
        SlownessLayer topSLayer = sMod.getSlownessLayer(topLayerNum, isPWave);
        SlownessLayer botSLayer = sMod.getSlownessLayer(botLayerNum, isPWave);
        if (botSLayer.getTopDepth() == getBotDepth() && botSLayer.getBotDepth() > getBotDepth()) {
            // gone one too far
            botLayerNum--;
            botSLayer = sMod.getSlownessLayer(botLayerNum, isPWave);
        }
        if(topSLayer.getTopDepth() != topBranch.getBotDepth()
                || botSLayer.getBotDepth() != getBotDepth()) {
            throw new TauModelException("difference: TauBranch not compatible with slowness sampling:"
                    + "\ntopDepth: "
                    + topBranch.getBotDepth()
                    + " "
                    + topSLayer.getTopDepth()
                    + "\nbotDepth: "
                    + getBotDepth()
                    + " "
                    + botSLayer.getBotDepth()
                    + "\n"
                    + topSLayer
                    + "\n"
                    + botSLayer);
        }
        // make sure indexP and indexS really correspond to
        // new ray parameters at the top of this branch
        SlownessLayer sLayer = sMod.getSlownessLayer(sMod.layerNumberBelow(topBranch.getBotDepth(),
                                                                           true),
                                                     true);
        if(indexP >= 0 && sLayer.getTopP() != rayParams[indexP]) {
            throw new TauModelException("difference: P wave index doesn't match top layer.\n "
                    + sMod.layerNumberBelow(topBranch.getBotDepth(), true)
                    + "\n"
                    + rayParams[indexP - 1]
                    + "\n"
                    + rayParams[indexP]
                    + "\n"
                    + rayParams[indexP + 1]
                    + "\nP="
                    + sLayer
                    + "\nS="
                    + sMod.getSlownessLayer(sMod.layerNumberBelow(topBranch.getBotDepth(),
                                                                  false),
                                            false));
        }
        sLayer = sMod.getSlownessLayer(sMod.layerNumberBelow(topBranch.getBotDepth(),
                                                             false),
                                       false);
        if(indexS >= 0 && sLayer.getTopP() != rayParams[indexS]) {
            throw new TauModelException("difference: S wave index doesn't match top layer. "
                    + sMod.layerNumberBelow(topBranch.getBotDepth(), false)
                    + " "
                    + rayParams[indexS - 1]
                    + " "
                    + rayParams[indexS]
                    + " " + rayParams[indexS + 1] + "\n" + sLayer);
        }
        sLayer = null;
        // construct the new TauBranch, going from the bottom of
        // the top half to the bottom of the whole branch
        TauBranch botBranch = new TauBranch(topBranch.getBotDepth(),
                                            getBotDepth(),
                                            isPWave);
        botBranch.maxRayParam = topBranch.getMinRayParam();
        botBranch.minTurnRayParam = getMinTurnRayParam();
        botBranch.minRayParam = getMinRayParam();
        double PRayParam = -1.0, SRayParam = -1.0;
        TimeDist timeDistP = new TimeDist();
        TimeDist timeDistS = new TimeDist();
        int arrayLength = dist.length;
        if(indexP != -1) {
            arrayLength++;
            PRayParam = rayParams[indexP];
            timeDistP = botBranch.calcTimeDist(sMod,
                                               topLayerNum,
                                               botLayerNum,
                                               PRayParam);
        }
        if(indexS != -1 && indexS != indexP) {
            arrayLength++;
            SRayParam = rayParams[indexS];
            timeDistS = botBranch.calcTimeDist(sMod,
                                               topLayerNum,
                                               botLayerNum,
                                               SRayParam);
        } else {
            // in case indexS==indexP then we only need one
            indexS = -1;
        }
        // allocate enough space
        botBranch.dist = new double[arrayLength];
        botBranch.time = new double[arrayLength];
        botBranch.tau = new double[arrayLength];
        if(indexP == -1) {
            // then both are -1 so no new ray parameters added
            for(int i = 0; i < dist.length; i++) {
                botBranch.dist[i] = dist[i] - topBranch.dist[i];
                botBranch.time[i] = time[i] - topBranch.time[i];
                botBranch.tau[i] = tau[i] - topBranch.tau[i];
            }
        } else {
            if(indexS == -1) {
                // only indexP != 0
                for(int i = 0; i < indexP; i++) {
                    botBranch.dist[i] = dist[i] - topBranch.dist[i];
                    botBranch.time[i] = time[i] - topBranch.time[i];
                    botBranch.tau[i] = tau[i] - topBranch.tau[i];
                }
                botBranch.dist[indexP] = timeDistP.distRadian;
                botBranch.time[indexP] = timeDistP.time;
                botBranch.tau[indexP] = timeDistP.time - PRayParam
                        * timeDistP.distRadian;
                for(int i = indexP; i < dist.length; i++) {
                    botBranch.dist[i + 1] = dist[i] - topBranch.dist[i + 1];
                    botBranch.time[i + 1] = time[i] - topBranch.time[i + 1];
                    botBranch.tau[i + 1] = tau[i] - topBranch.tau[i + 1];
                }
            } else {
                // both indexP and indexS != -1 so we have two new samples
                for(int i = 0; i < indexS; i++) {
                    botBranch.dist[i] = dist[i] - topBranch.dist[i];
                    botBranch.time[i] = time[i] - topBranch.time[i];
                    botBranch.tau[i] = tau[i] - topBranch.tau[i];
                }
                botBranch.dist[indexS] = timeDistS.distRadian;
                botBranch.time[indexS] = timeDistS.time;
                botBranch.tau[indexS] = timeDistS.time - SRayParam
                        * timeDistS.distRadian;
                for(int i = indexS; i < indexP; i++) {
                    botBranch.dist[i + 1] = dist[i] - topBranch.dist[i + 1];
                    botBranch.time[i + 1] = time[i] - topBranch.time[i + 1];
                    botBranch.tau[i + 1] = tau[i] - topBranch.tau[i + 1];
                }
                // put into indexP+1 as we have already shifted by 1
                // due to indexS
                botBranch.dist[indexP + 1] = timeDistP.distRadian;
                botBranch.time[indexP + 1] = timeDistP.time;
                botBranch.tau[indexP + 1] = timeDistP.time - PRayParam
                        * timeDistP.distRadian;
                for(int i = indexP; i < dist.length; i++) {
                    botBranch.dist[i + 2] = dist[i] - topBranch.dist[i + 2];
                    botBranch.time[i + 2] = time[i] - topBranch.time[i + 2];
                    botBranch.tau[i + 2] = tau[i] - topBranch.tau[i + 2];
                }
            }
        }
        return botBranch;
    }

    public void shiftBranch(int index) {
        double[] newDist = new double[dist.length + 1];
        System.arraycopy(dist, 0, newDist, 0, index);
        newDist[index] = 0.0;
        System.arraycopy(dist, index, newDist, index + 1, dist.length - index);
        dist = newDist;
        double[] newTime = new double[time.length + 1];
        System.arraycopy(time, 0, newTime, 0, index);
        newTime[index] = 0.0;
        System.arraycopy(time, index, newTime, index + 1, time.length - index);
        time = newTime;
        double[] newTau = new double[tau.length + 1];
        System.arraycopy(tau, 0, newTau, 0, index);
        newTau[index] = 0.0;
        System.arraycopy(tau, index, newTau, index + 1, tau.length - index);
        tau = newTau;
    }

    public TimeDist[] path(double rayParam,
                           boolean downgoing,
                           SlownessModel sMod) throws SlownessModelException {
        if(rayParam > getMaxRayParam()) {
            return null;
        }
        Assert.isTrue(rayParam >= 0.0, "ray parameter must not be negative.");
        int topLayerNum;
        int botLayerNum;
        try {
            topLayerNum = sMod.layerNumberBelow(getTopDepth(), isPWave);
            botLayerNum = sMod.layerNumberAbove(getBotDepth(), isPWave);
        } catch(NoSuchLayerException e) {
            throw new SlownessModelException("Caught NoSuchLayerException. This likely means the"
                    + "SlownessModel and TauModel are out of sync. "
                    + e.getMessage());
        }
        TimeDist[] thePath = new TimeDist[botLayerNum - topLayerNum + 1];
        int sLayerNum;
        int pathIndex = 0;
        double turnDepth;
        SlownessLayer sLayer, turnSLayer;
        /** check to make sure layers and branches are compatable. */
        sLayer = sMod.getSlownessLayer(topLayerNum, isPWave);
        if(sLayer.getTopDepth() != getTopDepth()) {
            throw new SlownessModelException("Branch and Slowness model are not compatible! "
                    + sLayer.getTopDepth()
                    + " != "
                    + getTopDepth()
                    + "=topDepth");
        }
        sLayer = sMod.getSlownessLayer(botLayerNum, isPWave);
        if(sLayer.getBotDepth() != getBotDepth()) {
            throw new SlownessModelException("Branch and Slowness model are not compatible! "
                    + sLayer.getBotDepth()
                    + " != "
                    + getBotDepth()
                    + "=botDepth");
        }
        if(downgoing) {
            sLayerNum = topLayerNum;
            sLayer = sMod.getSlownessLayer(sLayerNum, isPWave);
            while(sLayer.getBotP() >= rayParam && sLayerNum <= botLayerNum) {
                if(!sLayer.isZeroThickness()) {
                    thePath[pathIndex] = sMod.layerTimeDist(rayParam,
                                                            sLayerNum,
                                                            isPWave);
                    thePath[pathIndex].depth = sLayer.getBotDepth();
                    pathIndex++;
                }
                sLayerNum++;
                if(sLayerNum <= botLayerNum) {
                    sLayer = sMod.getSlownessLayer(sLayerNum, isPWave);
                }
            }
            if(sLayerNum <= botLayerNum && !sLayer.isZeroThickness()) {
                turnDepth = sLayer.bullenDepthFor(rayParam,
                                                  sMod.getRadiusOfEarth());
                turnSLayer = new SlownessLayer(sLayer.getTopP(),
                                               sLayer.getTopDepth(),
                                               rayParam,
                                               turnDepth);
                thePath[pathIndex] = turnSLayer.bullenRadialSlowness(rayParam,
                                                                     sMod.getRadiusOfEarth());
                thePath[pathIndex].depth = turnSLayer.getBotDepth();
                pathIndex++;
            }
        } else {
            // upgoing
            sLayerNum = botLayerNum;
            sLayer = sMod.getSlownessLayer(sLayerNum, isPWave);
            while(sLayer.getTopP() <= rayParam || sLayer.isZeroThickness()) {
                sLayerNum--;
                sLayer = sMod.getSlownessLayer(sLayerNum, isPWave);
            }
            if(sLayer.getBotP() < rayParam) {
                turnDepth = sLayer.bullenDepthFor(rayParam,
                                                  sMod.getRadiusOfEarth());
                turnSLayer = new SlownessLayer(sLayer.getTopP(),
                                               sLayer.getTopDepth(),
                                               rayParam,
                                               turnDepth);
                thePath[pathIndex] = turnSLayer.bullenRadialSlowness(rayParam,
                                                                     sMod.getRadiusOfEarth());
                thePath[pathIndex].depth = turnSLayer.getTopDepth();
                pathIndex++;
                sLayerNum--;
                if(sLayerNum >= topLayerNum) {
                    sLayer = sMod.getSlownessLayer(sLayerNum, isPWave);
                }
            }
            while(sLayerNum >= topLayerNum) {
                if(!sLayer.isZeroThickness()) {
                    thePath[pathIndex] = sMod.layerTimeDist(rayParam,
                                                            sLayerNum,
                                                            isPWave);
                    thePath[pathIndex].depth = sLayer.getTopDepth();
                    pathIndex++;
                }
                sLayerNum--;
                if(sLayerNum >= topLayerNum) {
                    sLayer = sMod.getSlownessLayer(sLayerNum, isPWave);
                }
            }
        }
        TimeDist[] tempPath = new TimeDist[pathIndex];
        System.arraycopy(thePath, 0, tempPath, 0, pathIndex);
        return tempPath;
    }

    public void writeToStream(DataOutputStream dos) throws IOException {
        dos.writeInt(getClass().getName().length());
        dos.writeBytes(getClass().getName());
        dos.writeDouble(getTopDepth());
        dos.writeDouble(getBotDepth());
        dos.writeDouble(getMaxRayParam());
        dos.writeDouble(getMinRayParam());
        dos.writeDouble(getMinTurnRayParam());
        dos.writeInt(dist.length);
        for(int i = 0; i < dist.length; i++) {
            dos.writeDouble(dist[i]);
        }
        dos.writeInt(time.length);
        for(int i = 0; i < time.length; i++) {
            dos.writeDouble(time[i]);
        }
        dos.writeInt(tau.length);
        for(int i = 0; i < tau.length; i++) {
            dos.writeDouble(tau[i]);
        }
    }

    public static TauBranch readFromStream(DataInputStream dis)
            throws IOException, ClassNotFoundException, InstantiationException,
            IllegalAccessException {
        int length;
        byte[] classString = new byte[dis.readInt()];
        dis.read(classString);
        Class tBranchClass = Class.forName(new String(classString));
        TauBranch tBranch = (TauBranch)tBranchClass.newInstance();
        tBranch.topDepth = dis.readDouble();
        tBranch.botDepth = dis.readDouble();
        tBranch.maxRayParam = dis.readDouble();
        tBranch.minRayParam = dis.readDouble();
        tBranch.minTurnRayParam = dis.readDouble();
        length = dis.readInt();
        tBranch.dist = new double[length];
        for(int i = 0; i < tBranch.dist.length; i++) {
            tBranch.dist[i] = dis.readDouble();
        }
        length = dis.readInt();
        tBranch.time = new double[length];
        for(int i = 0; i < tBranch.time.length; i++) {
            tBranch.time[i] = dis.readDouble();
        }
        length = dis.readInt();
        tBranch.tau = new double[length];
        for(int i = 0; i < tBranch.tau.length; i++) {
            tBranch.tau[i] = dis.readDouble();
        }
        return tBranch;
    }

    /**
     * Returns a clone of this TauBranch object. Note that super.clone() handles
     * all normal variables while the arrays need to be cloned separately to
     * generate a new array as opposed to a new reference to the old array.
     * 
     * @see Cloneable
     */
    public TauBranch clone() {
//        double[] newDist = new double[dist.length];
//        System.arraycopy(dist, 0, newDist, 0, dist.length);
//        double[] newTime = new double[time.length];
//        System.arraycopy(time, 0, newTime, 0, dist.length);
//        double[] newTau = new double[tau.length];
//        System.arraycopy(tau, 0, newTau, 0, dist.length);
        return new TauBranch(isPWave,
                             topDepth,
                             botDepth,
                             maxRayParam,
                             minTurnRayParam,
                             minRayParam,
                             dist, time, tau);
    }

    public String toString() {
        String desc = "Tau Branch\n";
        desc += " topDepth = " + getTopDepth() + "\n";
        desc += " botDepth = " + getBotDepth() + "\n";
        desc += " maxRayParam=" + getMaxRayParam() + " minTurnRayParam="
                + getMinTurnRayParam();
        desc += " minRayParam=" + getMinRayParam() + "\n";
        /*
         * for (int i=0;i <tau.length;i++) { desc += "i = "+i+" time="+time[i]+"
         * dist="+dist[i]+" tau="+tau[i]+"\n"; }
         */
        return desc;
    }
}
