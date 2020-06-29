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
import java.util.ArrayList;
import java.util.List;

/**
 * This class provides storage and methods for generating slowness-depth pairs.
 * 
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * 
 * 
 * 
 * @author H. Philip Crotwell
 * 
 */
public abstract class SlownessModel implements Serializable {

    public SlownessModel(VelocityModel vMod,
                         double minDeltaP,
                         double maxDeltaP,
                         double maxDepthInterval,
                         double maxRangeInterval,
                         double maxInterpError,
                         boolean allowInnerCoreS,
                         double slownessTolerance) throws NoSuchMatPropException, NoSuchLayerException, SlownessModelException {
        this.vMod = vMod;
        this.minDeltaP = minDeltaP;
        this.maxDeltaP = maxDeltaP;
        this.maxDepthInterval = maxDepthInterval;
        this.maxRangeInterval = maxRangeInterval;
        this.maxInterpError = maxInterpError;
        this.allowInnerCoreS = allowInnerCoreS;
        this.slownessTolerance = slownessTolerance;
        createSample();
    }
    public SlownessModel(double radiusOfEarth,
                         VelocityModel vMod,
                         List<CriticalDepth> criticalDepth,
                         List<DepthRange> highSlownessLayerDepthsP,
                         List<DepthRange> highSlownessLayerDepthsS,
                         List<DepthRange> fluidLayerDepths,
                         List<SlownessLayer> pLayers,
                         List<SlownessLayer> sLayers,
                         double minDeltaP,
                         double maxDeltaP,
                         double maxDepthInterval,
                         double maxRangeInterval,
                         double maxInterpError,
                         boolean allowInnerCoreS,
                         double slownessTolerance) {
        this.radiusOfEarth = radiusOfEarth;
        this.vMod = vMod;
        this.criticalDepths = criticalDepth;
        this.highSlownessLayerDepthsP = highSlownessLayerDepthsP;
        this.highSlownessLayerDepthsS = highSlownessLayerDepthsS;
        this.fluidLayerDepths = fluidLayerDepths;
        this.PLayers = pLayers;
        this.SLayers = sLayers;
        this.minDeltaP = minDeltaP;
        this.maxDeltaP = maxDeltaP;
        this.maxDepthInterval = maxDepthInterval;
        this.maxRangeInterval = maxRangeInterval;
        this.maxInterpError = maxInterpError;
        this.allowInnerCoreS = allowInnerCoreS;
        this.slownessTolerance = slownessTolerance;
    }

    /** True to enable debugging output. */
    transient static public boolean DEBUG = false;

    /** True to enable verbose output. */
    transient static public boolean verbose = false;

    /** Radius of the Earth in km, usually input from the velocity model. */
    protected double radiusOfEarth = 6371.0;

    /**
     * Velocity Model used to get slowness model. Usually set in
     * createSlowness().
     */
    protected VelocityModel vMod;

    /**
     * Stores the layer number for layers in the velocity model with a critical
     * point at their top. These form the "branches" of slowness sampling.
     * 
     * @see edu.sc.seis.TauP.CriticalDepth
     */
    protected List<CriticalDepth> criticalDepths = new ArrayList<CriticalDepth>();

    /**
     * Stores depth ranges that contains a high slowness zone for P. Stored as
     * DepthRange objects, containing the top depth and bottom depth.
     * 
     * @see DepthRange
     */
    protected List<DepthRange> highSlownessLayerDepthsP = new ArrayList<DepthRange>();

    /**
     * Stores depth ranges that contains a high slowness zone for S. Stored as
     * DepthRange objects, containing the top depth and bottom depth.
     * 
     * @see DepthRange
     */
    protected List<DepthRange> highSlownessLayerDepthsS = new ArrayList<DepthRange>();

    /**
     * Stores depth ranges that are fluid, ie S velocity is zero. Stored as
     * DepthRange objects, containing the top depth and bottom depth.
     * 
     * @see DepthRange
     */
    protected List<DepthRange> fluidLayerDepths = new ArrayList<DepthRange>();

    /** Initial length of the slowness vectors. */
    protected static int vectorLength = 256;

    /**
     * Stores the final slowness-depth layers for P waves. Stored as
     * SlownessLayer objects.
     * 
     * @see edu.sc.seis.TauP.SlownessLayer
     */
    protected List<SlownessLayer> PLayers = new ArrayList<SlownessLayer>(vectorLength);

    /**
     * Stores the final slowness-depth layers for S waves. Stored as
     * SlownessLayer objects. Note that SLayers and PLayers share the same
     * SlownessLayer object within fluid layers, so changes made to one will
     * affect the other.
     * 
     * @see edu.sc.seis.TauP.SlownessLayer
     */
    protected List<SlownessLayer> SLayers = new ArrayList<SlownessLayer>(vectorLength);

    /**
     * Minimum difference between successive slowness samples. The default is
     * 0.1 (km-sec/km or sec/rad for spherical, sec/km for flat models). This
     * keeps the sampling from becoming too fine. For example, a strong negative
     * S velocity gradient just above the CMB will cause the totally reflected
     * ScS too have an extremely large range of distances, over a very small
     * range of ray parameters. The distance check would otherwise force a very
     * fine sampling of this region. However since in this case time and
     * distance are likely to be very close to being linearly related, this sort
     * of sampling is overkill. So we ignore the distance check if the ray
     * parameter becomes smaller than minDeltaP.
     */
    protected double minDeltaP = 0.1;

    /**
     * Maximum difference between successive slowness samples. The default is
     * 11.0 (km-sec/km or sec/rad for spherical, sec/km for flat models). See
     * Buland and Chapman p1292
     */
    protected double maxDeltaP = 11.0;

    /**
     * Maximum difference between successive depth samples, default is 115 km.
     * See Buland and Chapman p1292
     */
    protected double maxDepthInterval = 115.0;

    /**
     * Maximum difference between successive ranges, in radians. The default is
     * 200 km / radiusOfEarth. See Buland and Chapman p1292.
     * 
     * @see radiusOfEarth
     */
    protected double maxRangeInterval = 200.0 / radiusOfEarth;

    protected double maxInterpError = .5;

    /**
     * Should we allow J phases, S waves in the inner core? If true, then the
     * slowness sampling for S will use the S velocity structure for the inner
     * core. If false, then we will use the P velocity structure for both the
     * inner and outer core for S waves as well as P waves. Disallowing inner
     * core S phases reduces the number of slowness samples significantly due to
     * the large geometrical spreading of S waves in the inner core. The default
     * is false.
     * 
     * @see minInnerCoreDepth
     */
    protected boolean allowInnerCoreS = true;

    public static final double DEFAULT_SLOWNESS_TOLERANCE = 1e-16;
    
    /**
     * Tolerance for slownesses. If two slownesses are closer that this value,
     * then we consider them to be identical. Basically this just provides some
     * protection against numerical "chatter".
     */
    protected double slownessTolerance = DEFAULT_SLOWNESS_TOLERANCE;

    /**
     * Just useful for calling methods that need to know whether to use P or S
     * waves.
     */
    public static final boolean PWAVE = true;

    /**
     * Just useful for calling methods that need to know whether to use P or S
     * waves.
     */
    public static final boolean SWAVE = false;

    // METHODS ----------------------------------------------------------------
    // Accessor methods
    public void setRadiusOfEarth(double radiusOfEarth) {
        this.radiusOfEarth = radiusOfEarth;
    }

    public void setMinDeltaP(double minDeltaP) {
        this.minDeltaP = minDeltaP;
    }

    public void setMaxDeltaP(double maxDeltaP) {
        this.maxDeltaP = maxDeltaP;
    }

    public void setMaxDepthInterval(double maxDepthInterval) {
        this.maxDepthInterval = maxDepthInterval;
    }

    /**
     * sets the maximum range interval for surface focus turning waves between
     * slowness samples, input in degrees.
     */
    public void setMaxRangeInterval(double maxRangeInterval) {
        this.maxRangeInterval = maxRangeInterval * Math.PI / 180.0;
    }

    /**
     * sets the maximum value of the estimated error due to linear
     * interpolation. Care should be taken not to set this too small as a very
     * large number of samples may be required. Note also that this is only an
     * estimate of the error, and thus the bound is by no means assured.
     */
    public void setMaxInterpError(double maxInterpError) {
        this.maxInterpError = maxInterpError;
    }

    public void setAllowInnerCoreS(boolean allowInnerCoreS) {
        this.allowInnerCoreS = allowInnerCoreS;
    }

    public void setSlownessTolerance(double slownessTolerance) {
        this.slownessTolerance = slownessTolerance;
    }

    // get accessor methods
    public VelocityModel getVelocityModel() {
        return vMod;
    }

    public final double getRadiusOfEarth() {
        return radiusOfEarth;
    }

    public final double getMinDeltaP() {
        return minDeltaP;
    }

    public final double getMaxDeltaP() {
        return maxDeltaP;
    }

    public final double getMaxDepthInterval() {
        return maxDepthInterval;
    }

    /**
     * @return the maximum range interval for surface focus turning waves
     *          between slowness samples output in degrees.
     */
    public final double getMaxRangeInterval() {
        return 180.0 * maxRangeInterval / Math.PI;
    }

    /**
     * gets the maximum value of the estimated error due to linear
     * interpolation. Care should be taken not to set this too small as a very
     * large number of samples may be required. Note also that this is only an
     * estimate of the error, and thus the bound is by no means assured.
     */
    public final double getMaxInterpError() {
        return maxInterpError;
    }

    public final boolean isAllowInnerCoreS() {
        return allowInnerCoreS;
    }

    public final double getSlownessTolerance() {
        return slownessTolerance;
    }

    public final int getNumCriticalDepths() {
        return criticalDepths.size();
    }

    public final CriticalDepth getCriticalDepth(int i) {
        return criticalDepths.get(i);
    }

    public final int getNumLayers(boolean isPWave) {
        if(isPWave) {
            return PLayers.size();
        } else {
            return SLayers.size();
        }
    }

    /**
     * @return the minimum ray parameter that turns, but is not reflected, at
     *          or above the given depth. Normally this is the slowness sample
     *          at the given depth, but if the depth is within a high slowness
     *          zone, then it may be smaller.
     */
    public double getMinTurnRayParam(double depth, boolean isPWave)
            throws NoSuchLayerException, SlownessModelException {
        double minPSoFar = Double.MAX_VALUE;
        SlownessLayer sLayer;
        List<SlownessLayer> layers;
        if(isPWave) {
            layers = PLayers;
        } else {
            layers = SLayers;
        }
        if(depthInHighSlowness(depth, Double.MAX_VALUE, isPWave)) {
            for(int i = 0; i < layers.size(); i++) {
                sLayer = getSlownessLayer(i, isPWave);
                if(sLayer.getBotDepth() == depth) {
                    minPSoFar = Math.min(minPSoFar, sLayer.getBotP());
                    return minPSoFar;
                } else if(sLayer.getBotDepth() > depth) {
                    minPSoFar = Math.min(minPSoFar,
                                         sLayer.evaluateAt_bullen(depth,
                                                                  getRadiusOfEarth()));
                    return minPSoFar;
                } else {
                    minPSoFar = Math.min(minPSoFar, sLayer.getBotP());
                }
            }
        } else {
            sLayer = getSlownessLayer(layerNumberAbove(depth, isPWave), isPWave);
            if(depth == sLayer.getBotDepth()) {
                minPSoFar = sLayer.getBotP();
            } else {
                minPSoFar = sLayer.evaluateAt_bullen(depth, getRadiusOfEarth());
            }
        }
        return minPSoFar;
    }

    /**
     * @return the minimum ray parameter that turns or is reflected at or above
     *          the given depth. Normally this is the slowness sample at the
     *          given depth, but if the depth is within a high slowness zone,
     *          then it may be smaller. Also, at first order discontinuities,
     *          there may be many slowness samples at the same depth.
     */
    public double getMinRayParam(double depth, boolean isPWave)
            throws NoSuchLayerException, SlownessModelException {
        double minPSoFar = getMinTurnRayParam(depth, isPWave);
        int i = layerNumberAbove(depth, isPWave);
        int j = layerNumberBelow(depth, isPWave);
        SlownessLayer sLayerAbove = getSlownessLayer(i, isPWave);
        SlownessLayer sLayerBelow = getSlownessLayer(j, isPWave);
        if(sLayerAbove.getBotDepth() == depth) {
            minPSoFar = Math.min(Math.min(minPSoFar, sLayerAbove.getBotP()),
                                 sLayerBelow.getTopP());
        }
        return minPSoFar;
    }

    /**
     * @return the DepthRange objects for all high slowness zones within the
     *          slowness model.
     */
    public DepthRange[] getHighSlowness(boolean isPWave) {
        List<DepthRange> highSlownessLayerDepths;
        if(isPWave) {
            highSlownessLayerDepths = highSlownessLayerDepthsP;
        } else {
            highSlownessLayerDepths = highSlownessLayerDepthsS;
        }
        DepthRange[] hsz = new DepthRange[highSlownessLayerDepths.size()];
        for(int i = 0; i < highSlownessLayerDepths.size(); i++) {
            hsz[i] = (DepthRange)highSlownessLayerDepths.get(i).clone();
        }
        return hsz;
    }

    /**
     * Returns the SlownessLayer of the requested waveType. This is NOT a clone
     * and any changes will possibly corrupt the SlownessModel.
     */
    protected SlownessLayer getSlownessLayer(int layerNum, boolean isPWave) {
        if(isPWave) {
            return PLayers.get(layerNum);
        } else {
            return SLayers.get(layerNum);
        }
    }

    protected List<SlownessLayer> getAllSlownessLayers(boolean isPWave) {
        if(isPWave) {
            return PLayers;
        } else {
            return SLayers;
        }
    }

    // Abstract methods
    public abstract double toSlowness(double velocity, double depth)
            throws SlownessModelException;

    public abstract double toVelocity(double slowness, double depth)
            throws SlownessModelException;

    public abstract TimeDist layerTimeDist(double rayParam,
                                           int layerNum,
                                           boolean isPWave)
            throws SlownessModelException;

    public abstract SlownessLayer toSlownessLayer(VelocityLayer vLayer,
                                                  boolean isPWave)
            throws SlownessModelException;

    public abstract double interpolate(double p,
                                       double topVelocity,
                                       double topDepth,
                                       double slope)
            throws SlownessModelException;

    // Defined methods
    /**
     * generate approximate distance, in radians, for a ray from a surface
     * source that turns at the bottom of the given slowness layer.
     * 
     * @exception NoSuchLayerException
     *                occurs if no layer in the velocity model contains the
     *                given depth.
     * @exception SlownessModelException
     *                occurs if getNumLayers() == 0 as we cannot compute a
     *                distance without a layer.
     */
    public TimeDist approxDistance(int slownessTurnLayer,
                                   double p,
                                   boolean isPWave)
            throws NoSuchLayerException, SlownessModelException {
        /*
         * First, if slowness contains less than slownessTurnLayer elements then
         * we can't calculate a distance, otherwise we must signal an exception.
         */
        if(slownessTurnLayer >= getNumLayers(isPWave)) {
            throw new SlownessModelException("Can't calculate a distance when "
                    + "slownessTurnLayer >= getNumLayers(" + isPWave + ")\n"
                    + " slownessTurnLayer=" + slownessTurnLayer
                    + " getNumLayers()=" + getNumLayers(isPWave));
        }
        if(p < 0.0) {
            throw new SlownessModelException("approxDistance: Ray parameter is negative!!!"
                    + p + " slownessTurnLayer=" + slownessTurnLayer);
        }
        /*
         * OK, now we are able to do the calculations for the approximate
         * distance, hopefully without errors.
         */
        TimeDist td = new TimeDist(p);
        for(int layerNum = 0; layerNum <= slownessTurnLayer; layerNum++) {
            td.add(layerTimeDist(p, layerNum, isPWave));
        }
        /*
         * Return 2.0*distance and time because there is a downgoing as well as
         * up going leg, which are equal because this is for a surface source.
         */
        td.distRadian *= 2.0;
        td.time *= 2.0;
        return td;
    }

    /**
     * Determines if the given depth and corresponding slowness is contained
     * within a high slowness zone. Whether the high slowness zone includes its
     * upper boundary and its lower boundaries depends upon the ray parameter.
     * The slowness at the depth is needed because if depth happens to
     * correspond to a discontinuity that marks the bottom of the high slowness
     * zone but the ray is actually a total reflection then it is not part of
     * the high slowness zone. Calls depthInHighSlowness(double, double,
     * DepthRange).
     * 
     * @see depthInHighSlowness.
     */
    public boolean depthInHighSlowness(double depth,
                                       double rayParam,
                                       boolean isPWave) {
        DepthRange highSZoneDepth = new DepthRange();
        return depthInHighSlowness(depth, rayParam, highSZoneDepth, isPWave);
    }

    /**
     * Determines if the given depth and corresponding slowness is contained
     * within a high slowness zone. Whether the high slowness zone includes its
     * upper boundary and its lower boundaries depends upon the ray parameter.
     * The slowness at the depth is needed because if depth happens to
     * correspond to a discontinuity that marks the bottom of the high slowness
     * zone but the ray is actually a total reflection then it is not part of
     * the high slowness zone. The ray parameter that delimits the zone, ie it
     * can turn at the top and the bottom, is in the zone at the top, but out of
     * the zone at the bottom.
     */
    public boolean depthInHighSlowness(double depth,
                                       double rayParam,
                                       DepthRange highSZoneDepth,
                                       boolean isPWave) {
        DepthRange tempRange;
        List<DepthRange> highSlownessLayerDepths;
        if(isPWave) {
            highSlownessLayerDepths = highSlownessLayerDepthsP;
        } else {
            highSlownessLayerDepths = highSlownessLayerDepthsS;
        }
        for(int i = 0; i < highSlownessLayerDepths.size(); i++) {
            tempRange = highSlownessLayerDepths.get(i);
            if(tempRange.topDepth <= depth && depth <= tempRange.botDepth) {
                highSZoneDepth.topDepth = tempRange.topDepth;
                highSZoneDepth.botDepth = tempRange.botDepth;
                highSZoneDepth.rayParam = tempRange.rayParam;
                if(rayParam > tempRange.rayParam
                        || (rayParam == tempRange.rayParam && depth == tempRange.topDepth)) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * Determines if the given depth is contained within a fluid zone. The fluid
     * zone includes its upper boundary but not its lower boundary. Calls
     * depthInFluid(double, DepthRange).
     * 
     * @see depthInFluid(double, DepthRange).
     */
    public boolean depthInFluid(double depth) {
        DepthRange fluidZoneDepth = new DepthRange();
        return depthInFluid(depth, fluidZoneDepth);
    }

    /**
     * Determines if the given depth is contained within a fluid zone. The fluid
     * zone includes its upper boundary but not its lower boundary. The top and
     * bottom of the fluid zone are returned in DepthRange.
     */
    public boolean depthInFluid(double depth, DepthRange fluidZoneDepth) {
        DepthRange tempRange;
        for(int i = 0; i < fluidLayerDepths.size(); i++) {
            tempRange = fluidLayerDepths.get(i);
            if(tempRange.topDepth <= depth && depth < tempRange.botDepth) {
                fluidZoneDepth.topDepth = tempRange.topDepth;
                fluidZoneDepth.botDepth = tempRange.botDepth;
                return true;
            }
        }
        return false;
    }

    /*
     * Splits a slowness layer into two slowness layers. returns a
     * SplitLayerInfo object with neededSplit=true if a layer was actually
     * split, false otherwise, movedSample=true if a layer was very close, and
     * so moving the layers depth is better than making a very thin layer,
     * rayParam= the new ray parameter, if the layer was split. The
     * interpolation for splitting a layer is a Bullen p=Ar^B and so does not
     * directly use information from the VelocityModel.
     */
    public SplitLayerInfo splitLayer(double depth, boolean isPWave)
            throws SlownessModelException, NoSuchLayerException {
        int layerNum = layerNumberAbove(depth, isPWave);
        SlownessLayer sLayer = getSlownessLayer(layerNum, isPWave);
        if(sLayer.getTopDepth() == depth || sLayer.getBotDepth() == depth) {
            /*
             * depth is already on a slowness layer boundary so we don't need to
             * split any slowness layers.
             */
            return new SplitLayerInfo(this, false, false, 0.0);
        } else if(Math.abs(sLayer.getTopDepth() - depth) < 0.000001) {
            /*
             * check for very thin layers, just move the layer to hit the
             * boundary
             */
            List<SlownessLayer> allLayers = getAllSlownessLayers(isPWave);
            List<SlownessLayer> outLayers = new ArrayList<SlownessLayer>(allLayers.size());
            outLayers.addAll(allLayers);
            outLayers.set(layerNum, new SlownessLayer(sLayer.getTopP(),
                                                      depth,
                                                      sLayer.getBotP(),
                                                      sLayer.getBotDepth()));
            sLayer = getSlownessLayer(layerNum - 1, isPWave);
            outLayers.set(layerNum - 1, new SlownessLayer(sLayer.getTopP(),
                                                          sLayer.getTopDepth(),
                                                          sLayer.getBotP(),
                                                          depth));
            List<SlownessLayer> outPLayers;
            List<SlownessLayer> outSLayers;
            if(isPWave) {
                outPLayers = outLayers;
                outSLayers = SLayers;
            } else {
                outPLayers = PLayers;
                outSLayers = outLayers;
            }
            SlownessModel out = new SphericalSModel(radiusOfEarth,
                                                    vMod,
                                                    criticalDepths,
                                                    highSlownessLayerDepthsP,
                                                    highSlownessLayerDepthsS,
                                                    fluidLayerDepths,
                                                    outPLayers,
                                                    outSLayers,
                                                    minDeltaP,
                                                    maxDeltaP,
                                                    maxDepthInterval,
                                                    maxRangeInterval,
                                                    maxInterpError,
                                                    allowInnerCoreS,
                                                    slownessTolerance);
            return new SplitLayerInfo(out, false, true, sLayer.getBotP());
        } else if(Math.abs(depth - sLayer.getBotDepth()) < 0.000001) {
            /*
             * check for very thin layers, just move the layer to hit the
             * boundary
             */
            List<SlownessLayer> allLayers = getAllSlownessLayers(isPWave);
            List<SlownessLayer> outLayers = new ArrayList<SlownessLayer>(allLayers.size());
            outLayers.addAll(allLayers);
            outLayers.set(layerNum, new SlownessLayer(sLayer.getTopP(),
                                                      sLayer.getTopDepth(),
                                                      sLayer.getBotP(),
                                                      depth));
            sLayer = getSlownessLayer(layerNum + 1, isPWave);
            outLayers.set(layerNum + 1, new SlownessLayer(sLayer.getTopP(),
                                                          depth,
                                                          sLayer.getBotP(),
                                                          sLayer.getBotDepth()));
            List<SlownessLayer> outPLayers, outSLayers;
            if(isPWave) {
                outPLayers = outLayers;
                outSLayers = SLayers;
            } else {
                outPLayers = PLayers;
                outSLayers = outLayers;
            }
            SlownessModel out = new SphericalSModel(radiusOfEarth,
                                                    vMod,
                                                    criticalDepths,
                                                    highSlownessLayerDepthsP,
                                                    highSlownessLayerDepthsS,
                                                    fluidLayerDepths,
                                                    outPLayers,
                                                    outSLayers,
                                                    minDeltaP,
                                                    maxDeltaP,
                                                    maxDepthInterval,
                                                    maxRangeInterval,
                                                    maxInterpError,
                                                    allowInnerCoreS,
                                                    slownessTolerance);
            return new SplitLayerInfo(out, false, true, sLayer.getBotP());
        } else {
            double p = sLayer.evaluateAt_bullen(depth, radiusOfEarth);
            SlownessLayer topLayer, botLayer;
            topLayer = new SlownessLayer(sLayer.getTopP(),
                                         sLayer.getTopDepth(),
                                         p,
                                         depth);
            botLayer = new SlownessLayer(p,
                                         depth,
                                         sLayer.getBotP(),
                                         sLayer.getBotDepth());
            
            List<SlownessLayer> allLayers = getAllSlownessLayers(isPWave);
            List<SlownessLayer> outLayers = new ArrayList<SlownessLayer>(allLayers.size());
            outLayers.addAll(allLayers);
            outLayers.remove(layerNum);
            outLayers.add(layerNum, botLayer);
            outLayers.add(layerNum, topLayer);
            
            List<SlownessLayer> outPLayers, outSLayers;
            // fix critical layers since we have added a slowness layer
            List<CriticalDepth> outCriticalDepths = new ArrayList<CriticalDepth>();
            outCriticalDepths.addAll(criticalDepths);
            fixCriticalDepths(outCriticalDepths, layerNum, isPWave);
            if(isPWave) {
                outPLayers = outLayers;
                outSLayers = fixOtherLayers(SLayers,
                                            p,
                                            sLayer,
                                            topLayer,
                                            botLayer,
                                            outCriticalDepths,
                                            ! isPWave);
            } else {
                outPLayers = fixOtherLayers(PLayers,
                                            p,
                                            sLayer,
                                            topLayer,
                                            botLayer,
                                            outCriticalDepths,
                                            ! isPWave);
                outSLayers = outLayers;
            }
            SlownessModel out = new SphericalSModel(radiusOfEarth,
                                                    vMod,
                                                    outCriticalDepths,
                                                    highSlownessLayerDepthsP,
                                                    highSlownessLayerDepthsS,
                                                    fluidLayerDepths,
                                                    outPLayers,
                                                    outSLayers,
                                                    minDeltaP,
                                                    maxDeltaP,
                                                    maxDepthInterval,
                                                    maxRangeInterval,
                                                    maxInterpError,
                                                    allowInnerCoreS,
                                                    slownessTolerance);
            return new SplitLayerInfo(out, true, false, p);
        }
    }

    private void fixCriticalDepths(List<CriticalDepth> criticalDepths, int layerNum, boolean isPWave) {
        for(int i = 0; i < criticalDepths.size(); i++) {
            CriticalDepth cd = criticalDepths.get(i);
            if(cd.getLayerNum(isPWave) > layerNum) {
                if(isPWave) {
                    criticalDepths.set(i,
                                          new CriticalDepth(cd.getDepth(),
                                                            cd.getVelLayerNum(),
                                                            cd.getPLayerNum() + 1,
                                                            cd.getSLayerNum()));
                } else {
                    criticalDepths.set(i,
                                          new CriticalDepth(cd.getDepth(),
                                                            cd.getVelLayerNum(),
                                                            cd.getPLayerNum(),
                                                            cd.getSLayerNum() + 1));
                }
            }
        }
    }
    private List<SlownessLayer> fixOtherLayers(List<SlownessLayer> otherLayers,
                                  double p,
                                  SlownessLayer changedLayer,
                                  SlownessLayer newTopLayer,
                                  SlownessLayer newBotLayer,
                                  List<CriticalDepth> criticalDepths,
                                  boolean isPWave) throws SlownessModelException {
        List<SlownessLayer> out = new ArrayList<SlownessLayer>();
        out.addAll(otherLayers);
        int otherIndex = otherLayers.indexOf(changedLayer);
        // now make sure we keep the sampling consistant
        // if in a fluid, then both wavetypes will share a single
        // slowness layer object. Otherwise indexOf returns -1
        if(otherIndex != -1) {
            out.remove(otherIndex);
            out.add(otherIndex, newBotLayer);
            out.add(otherIndex, newTopLayer);
        }
        for(int otherLayerNum = 0; otherLayerNum < out.size(); otherLayerNum++) {
            SlownessLayer sLayer = out.get(otherLayerNum);
            if((sLayer.getTopP() - p) * (p - sLayer.getBotP()) > 0.0) {
                // found a slowness layer with the other wave type that
                // contains the new slowness sample
                SlownessLayer topLayer, botLayer;
                topLayer = new SlownessLayer(sLayer.getTopP(),
                                             sLayer.getTopDepth(),
                                             p,
                                             sLayer.bullenDepthFor(p, radiusOfEarth));
                botLayer = new SlownessLayer(p,
                                             topLayer.getBotDepth(),
                                             sLayer.getBotP(),
                                             sLayer.getBotDepth());
                out.remove(otherLayerNum);
                out.add(otherLayerNum, botLayer);
                out.add(otherLayerNum, topLayer);
                
                // fix critical layers since we have added a slowness layer
                fixCriticalDepths(criticalDepths, otherLayerNum, ! isPWave);
                otherLayerNum++; // skip next layer as it was just added
            }
        }
        return out;
    }

    /**
     * Finds all critical points within a velocity model. Critical points are
     * first order discontinuities in velocity/slowness, local extrema in
     * slowness. A high slowness zone is a low velocity zone, but it is possible
     * to have a slight low velocity zone within a spherical earth that is not a
     * high slowness zone and thus does not exhibit any of the pathological
     * behavior of a low velocity zone.
     * 
     * @exception NoSuchMatPropException
     *                occurs if wavetype is not recognized.
     * @exception SlownessModelException
     *                occurs if validate() returns false, this indicates a bug
     *                in the code.
     */
    protected void findCriticalPoints() throws SlownessModelException {
        double minPSoFar = Double.MAX_VALUE;
        double minSSoFar = Double.MAX_VALUE;
        DepthRange highSlownessZoneP = new DepthRange();
        DepthRange highSlownessZoneS = new DepthRange();
        boolean inHighSlownessZoneP = false;
        boolean inHighSlownessZoneS = false;
        DepthRange fluidZone = new DepthRange();
        boolean inFluidZone = false;
        boolean belowOuterCore = false; /*
                                         * are we in the inner core, see
                                         * allowInnerCoreS.
                                         */
        VelocityLayer prevVLayer, currVLayer;
        SlownessLayer prevSLayer, currSLayer, prevPLayer, currPLayer;
        // First remove any critical points previously stored
        highSlownessLayerDepthsP.clear();
        highSlownessLayerDepthsS.clear();
        criticalDepths.clear();
        fluidLayerDepths.clear();
        // Initialize the current velocity layer
        // to be zero thickness layer with values at the surface
        currVLayer = vMod.getVelocityLayer(0);
        currVLayer = new VelocityLayer(0,
                                       currVLayer.getTopDepth(),
                                       currVLayer.getTopDepth(),
                                       currVLayer.getTopPVelocity(),
                                       currVLayer.getTopPVelocity(),
                                       currVLayer.getTopSVelocity(),
                                       currVLayer.getTopSVelocity(),
                                       currVLayer.getTopDensity(),
                                       currVLayer.getTopDensity(),
                                       currVLayer.getTopQp(),
                                       currVLayer.getTopQp(),
                                       currVLayer.getTopQs(),
                                       currVLayer.getTopQs());
        currSLayer = toSlownessLayer(currVLayer, SWAVE);
        currPLayer = toSlownessLayer(currVLayer, PWAVE);
        // We know that the top is always a critical slowness so add 0
        criticalDepths.add(new CriticalDepth(0.0, 0, 0, 0));
        /* Check to see if we start in a fluid zone. */
        if(!inFluidZone && currVLayer.getTopSVelocity() == 0.0) {
            inFluidZone = true;
            fluidZone = new DepthRange();
            fluidZone.topDepth = currVLayer.getTopDepth();
            currSLayer = currPLayer;
        }
        if(minSSoFar > currSLayer.getTopP()) {
            minSSoFar = currSLayer.getTopP();
        }
        if(minPSoFar > currPLayer.getTopP()) {
            minPSoFar = currPLayer.getTopP();
        }
        for(int layerNum = 0; layerNum < vMod.getNumLayers(); layerNum++) {
            prevVLayer = currVLayer;
            prevSLayer = currSLayer;
            prevPLayer = currPLayer;
            currVLayer = vMod.getVelocityLayerClone(layerNum);
            /*
             * If we are not already in a fluid check to see if we have just
             * entered a fluid zone.
             */
            if(!inFluidZone && currVLayer.getTopSVelocity() == 0.0) {
                inFluidZone = true;
                fluidZone = new DepthRange();
                fluidZone.topDepth = currVLayer.getTopDepth();
            }
            /*
             * If we are already in a fluid check to see if we have just exited
             * it.
             */
            if(inFluidZone && currVLayer.getTopSVelocity() != 0.0) {
                if(prevVLayer.getBotDepth() > vMod.getIocbDepth()) {
                    belowOuterCore = true;
                }
                inFluidZone = false;
                fluidZone.botDepth = prevVLayer.getBotDepth();
                fluidLayerDepths.add(fluidZone);
            }

            currPLayer = toSlownessLayer(currVLayer, PWAVE);
            
            /*
             * If we are in a fluid zone ( S velocity = 0.0 ) or if we are below
             * the outer core and allowInnerCoreS=false then use the P velocity
             * structure to look for critical points.
             */
            if(inFluidZone || (belowOuterCore && !allowInnerCoreS)) {
                currSLayer = currPLayer;
            } else {
                currSLayer = toSlownessLayer(currVLayer, SWAVE);
            }
            if(prevSLayer.getBotP() != currSLayer.getTopP()
                    || prevPLayer.getBotP() != currPLayer.getTopP()) {
                // first order discontinuity
                criticalDepths.add(new CriticalDepth(currSLayer.getTopDepth(),
                                                                 layerNum,
                                                                 -1,
                                                                 -1));
                if(DEBUG) {
                    System.out.println("first order discontinuity, depth="
                            + currSLayer.getTopDepth());
                    System.out.println(prevSLayer + "\n" + currSLayer);
                    System.out.println(prevPLayer + "\n" + currPLayer);
                }
                if(inHighSlownessZoneS && (currSLayer.getTopP() < minSSoFar)) {
                    // top of current layer is the bottom of a high slowness
                    // zone.
                    if(DEBUG) {
                        System.out.println("top of current layer is the bottom"
                                + " of a high slowness zone.");
                    }
                    highSlownessZoneS.botDepth = currSLayer.getTopDepth();
                    highSlownessLayerDepthsS.add(highSlownessZoneS);
                    inHighSlownessZoneS = false;
                }
                if(inHighSlownessZoneP && (currPLayer.getTopP() < minPSoFar)) {
                    // top of current layer is the bottom of a high slowness
                    // zone.
                    if(DEBUG) {
                        System.out.println("top of current layer is the bottom"
                                + " of a high slowness zone.");
                    }
                    highSlownessZoneP.botDepth = currSLayer.getTopDepth();
                    highSlownessLayerDepthsP.add(highSlownessZoneP);
                    inHighSlownessZoneP = false;
                }
                /*
                 * Update minPSoFar and minSSoFar as all total reflections off
                 * of the top of the discontinuity are ok even though below the
                 * discontinuity could be the start of a high slowness zone.
                 */
                if(minPSoFar > currPLayer.getTopP()) {
                    minPSoFar = currPLayer.getTopP();
                }
                if(minSSoFar > currSLayer.getTopP()) {
                    minSSoFar = currSLayer.getTopP();
                }
                if(!inHighSlownessZoneS
                        && (prevSLayer.getBotP() < currSLayer.getTopP() || currSLayer.getTopP() < currSLayer.getBotP())) {
                    // start of a high slowness zone
                    if(DEBUG) {
                        System.out.println("Found S high slowness at first order "
                                + "discontinuity, layer = " + layerNum);
                    }
                    inHighSlownessZoneS = true;
                    highSlownessZoneS = new DepthRange();
                    highSlownessZoneS.topDepth = currSLayer.getTopDepth();
                    highSlownessZoneS.rayParam = minSSoFar;
                }
                if(!inHighSlownessZoneP
                        && (prevPLayer.getBotP() < currPLayer.getTopP() || currPLayer.getTopP() < currPLayer.getBotP())) {
                    // start of a high slowness zone
                    if(DEBUG) {
                        System.out.println("Found P high slowness at first order "
                                + "discontinuity, layer = " + layerNum);
                    }
                    inHighSlownessZoneP = true;
                    highSlownessZoneP = new DepthRange();
                    highSlownessZoneP.topDepth = currPLayer.getTopDepth();
                    highSlownessZoneP.rayParam = minPSoFar;
                }
            } else {
                if((prevSLayer.getTopP() - prevSLayer.getBotP())
                        * (prevSLayer.getBotP() - currSLayer.getBotP()) < 0.0
                        || (prevPLayer.getTopP() - prevPLayer.getBotP())
                                * (prevPLayer.getBotP() - currPLayer.getBotP()) < 0.0) {
                    // local slowness extrema
                    criticalDepths.add(new CriticalDepth(currSLayer.getTopDepth(),
                                                                     layerNum,
                                                                     -1,
                                                                     -1));
                    if(DEBUG) {
                        System.out.println("local slowness extrema, depth="
                                + currSLayer.getTopDepth());
                    }
                    if(!inHighSlownessZoneP
                            && (currPLayer.getTopP() < currPLayer.getBotP())) {
                        // start of a high slowness zone
                        if(DEBUG) {
                            System.out.println("start of a P high slowness zone,"
                                    + " local slowness extrema, minPSoFar="
                                    + minPSoFar);
                        }
                        inHighSlownessZoneP = true;
                        highSlownessZoneP = new DepthRange();
                        highSlownessZoneP.topDepth = currPLayer.getTopDepth();
                        highSlownessZoneP.rayParam = minPSoFar;
                    }
                    if(!inHighSlownessZoneS
                            && (currSLayer.getTopP() < currSLayer.getBotP())) {
                        // start of a high slowness zone
                        if(DEBUG) {
                            System.out.println("start of a S high slowness zone,"
                                    + " local slowness extrema, minSSoFar="
                                    + minSSoFar);
                        }
                        inHighSlownessZoneS = true;
                        highSlownessZoneS = new DepthRange();
                        highSlownessZoneS.topDepth = currSLayer.getTopDepth();
                        highSlownessZoneS.rayParam = minSSoFar;
                    }
                }
            }
            if(inHighSlownessZoneP && (currPLayer.getBotP() < minPSoFar)) {
                // layer contains the bottom of a high slowness zone.
                if(DEBUG) {
                    System.out.println("layer contains the bottom of a P "
                            + "high slowness zone. minPSoFar=" + minPSoFar
                            + " " + currPLayer);
                }
                highSlownessZoneP.botDepth = findDepth(minPSoFar,
                                                       layerNum,
                                                       layerNum,
                                                       PWAVE);
                highSlownessLayerDepthsP.add(highSlownessZoneP);
                inHighSlownessZoneP = false;
            }
            if(inHighSlownessZoneS && (currSLayer.getBotP() < minSSoFar)) {
                // layer contains the bottom of a high slowness zone.
                if(DEBUG) {
                    System.out.println("layer contains the bottom of a S "
                            + "high slowness zone. minSSoFar=" + minSSoFar
                            + " " + currSLayer);
                }
                // in fluid layers we want to check PWAVE structure when looking for S wave critical points
                highSlownessZoneS.botDepth = findDepth(minSSoFar,
                                                       layerNum,
                                                       layerNum,
                                                       (currSLayer == currPLayer)?PWAVE:SWAVE);
                highSlownessLayerDepthsS.add(highSlownessZoneS);
                inHighSlownessZoneS = false;
            }
            if(minPSoFar > currPLayer.getBotP()) {
                minPSoFar = currPLayer.getBotP();
            }
            if(minPSoFar > currPLayer.getTopP()) {
                minPSoFar = currPLayer.getTopP();
            }
            if(minSSoFar > currSLayer.getBotP()) {
                minSSoFar = currSLayer.getBotP();
            }
            if(minSSoFar > currSLayer.getTopP()) {
                minSSoFar = currSLayer.getTopP();
            }
            if(DEBUG && inHighSlownessZoneS) {
                System.out.println("In S high slowness zone, layerNum = "
                        + layerNum + " minSSoFar=" + minSSoFar);
            }
            if(DEBUG && inHighSlownessZoneP) {
                System.out.println("In P high slowness zone, layerNum = "
                        + layerNum + " minPSoFar=" + minPSoFar);
            }
        }
        // We know that the bottommost depth is always a critical slowness,
        // so we add vMod.getNumLayers()
        criticalDepths.add(new CriticalDepth(getRadiusOfEarth(),
                                                         vMod.getNumLayers(),
                                                         -1,
                                                         -1));
        // Check if the bottommost depth is contained within a high slowness
        // zone, might happen in a flat non-whole-earth model
        if(inHighSlownessZoneS) {
            highSlownessZoneS.botDepth = currVLayer.getBotDepth();
            highSlownessLayerDepthsS.add(highSlownessZoneS);
        }
        if(inHighSlownessZoneP) {
            highSlownessZoneP.botDepth = currVLayer.getBotDepth();
            highSlownessLayerDepthsP.add(highSlownessZoneP);
        }
        /*
         * Check if the bottommost depth is contained within a fluid zone, this
         * would be the case if we have a non whole earth model with the bottom
         * in the outer core or if allowInnerCoreS == false and we want to use
         * the P velocity structure in the inner core.
         */
        if(inFluidZone) {
            fluidZone.botDepth = currVLayer.getBotDepth();
            fluidLayerDepths.add(fluidZone);
        }
        if(DEBUG && criticalDepths.size() != 0) {
            int botCriticalLayerNum, topCriticalLayerNum;
            String desc = "**** Critical Velocity Layers ************************\n";
            botCriticalLayerNum = criticalDepths.get(0).getVelLayerNum() - 1;
            for(int criticalNum = 1; criticalNum < criticalDepths.size(); criticalNum++) {
                topCriticalLayerNum = botCriticalLayerNum + 1;
                botCriticalLayerNum = criticalDepths.get(criticalNum).getVelLayerNum() - 1;
                desc += " " + topCriticalLayerNum + "," + botCriticalLayerNum;
            }
            System.out.println(desc);
        }
        if(DEBUG && highSlownessLayerDepthsP.size() != 0) {
            for(int layerNum = 0; layerNum < highSlownessLayerDepthsP.size(); layerNum++) {
                System.out.println(highSlownessLayerDepthsP.get(layerNum));
            }
        }
        if(DEBUG && highSlownessLayerDepthsS.size() != 0) {
            for(int layerNum = 0; layerNum < highSlownessLayerDepthsS.size(); layerNum++) {
                System.out.println(highSlownessLayerDepthsS.get(layerNum));
            }
        }
        if(!validate()) {
            throw new SlownessModelException("Validation Failed!");
        }
    }

    /**
     * Finds a depth corresponding to a slowness over the whole VelocityModel.
     * Calls findDepth(double, int, int, char).
     */
    public double findDepth(double rayParam, boolean isPWave)
            throws SlownessModelException {
        return findDepth(rayParam, 0, vMod.getNumLayers() - 1, isPWave);
    }

    /**
     * Finds a depth corresponding to a slowness between two given depths in the
     * Velocity Model. Calls findDepth(double, int, int, char).
     */
    public double findDepth(double rayParam,
                            double topDepth,
                            double botDepth,
                            boolean isPWave) throws SlownessModelException {
        try {
            int topLayerNum = vMod.layerNumberBelow(topDepth);
            if(vMod.getVelocityLayer(topLayerNum).getBotDepth() == topDepth) {
                topLayerNum++;
            }
            int botLayerNum = vMod.layerNumberAbove(botDepth);
            return findDepth(rayParam, topLayerNum, botLayerNum, isPWave);
        } catch(NoSuchLayerException e) {
            throw new SlownessModelException(e.getMessage());
        }
    }

    /**
     * Finds a depth corresponding to a slowness between two given velocity
     * layers, including the top and the bottom. We also check to see if the
     * slowness is less than the bottom slowness of these layers but greater
     * than the top slowness of the next deeper layer. This corresponds to a
     * total reflection. In this case a check needs to be made to see if this is
     * an S wave reflecting off of a fluid layer, use P velocity below in this
     * case. We assume that slowness is monotonic within these layers and
     * therefore there is only one depth with the given slowness. This means we
     * return the first depth that we find.
     * 
     * @exception SlownessModelException
     *                occurs if topCriticalLayer > botCriticalLayer because
     *                there are no layers to search, or if there is an increase
     *                in slowness, ie a negative velocity gradient, that just
     *                balances the decrease in slowness due to the spherical
     *                earth, or if the ray parameter p is not contained within
     *                the specified layer range.
     */
    public double findDepth(double p,
                            int topCriticalLayer,
                            int botCriticalLayer,
                            boolean isPWave) throws SlownessModelException {
        VelocityLayer velLayer = null;
        double topP = Double.MAX_VALUE, botP = Double.MAX_VALUE;
        double topVelocity, botVelocity;
        double depth;
        double slope;
        char waveType;
        if(isPWave) {
            waveType = 'P';
        } else {
            waveType = 'S';
        }
        try {
            if(topCriticalLayer > botCriticalLayer) {
                throw new SlownessModelException("findDepth: no layers to search!: "
                        + "topCriticalLayer = "
                        + topCriticalLayer
                        + "botCriticalLayer = " + botCriticalLayer);
            }
            for(int layerNum = topCriticalLayer; layerNum <= botCriticalLayer; layerNum++) {
                velLayer = (VelocityLayer)vMod.getVelocityLayer(layerNum);
                topVelocity = velLayer.evaluateAtTop(waveType);
                botVelocity = velLayer.evaluateAtBottom(waveType);
                topP = toSlowness(topVelocity, velLayer.getTopDepth());
                botP = toSlowness(botVelocity, velLayer.getBotDepth());
                /*
                 * check to see if we are within chatter level of the top or
                 * bottom and if so then return that depth.
                 */
                if(Math.abs(topP - p) < slownessTolerance) {
                    return velLayer.getTopDepth();
                }
                if(Math.abs(p - botP) < slownessTolerance) {
                    return velLayer.getBotDepth();
                }
                if((topP - p) * (p - botP) >= 0.0) { // found the layer
                    // containing p
                    /*
                     * We interpolate assuming that velocity is linear within
                     * this interval. So slope is the slope for velocity versus
                     * depth.
                     */
                    slope = (botVelocity - topVelocity)
                            / (velLayer.getBotDepth() - velLayer.getTopDepth());
                    depth = interpolate(p,
                                        topVelocity,
                                        velLayer.getTopDepth(),
                                        slope);
                    return depth;
                } else if(layerNum == topCriticalLayer
                        && Math.abs(p - topP) < slownessTolerance) {
                    /*
                     * Check to see if p is just outside the topmost layer. If
                     * so than return the top depth.
                     */
                    return velLayer.getTopDepth();
                }
                /*
                 * Is p a total reflection? botP is the slowness at the bottom
                 * of the last velocity layer from the previous loop, set topP
                 * to be the slowness at the top of the next layer.
                 */
                if(layerNum < vMod.getNumLayers() - 1) {
                    velLayer = (VelocityLayer)vMod.getVelocityLayer(layerNum + 1);
                    topVelocity = velLayer.evaluateAtTop(waveType);
                    if(!isPWave && depthInFluid(velLayer.getTopDepth())) {
                        /*
                         * Special case for S waves above a fluid. If top next
                         * layer is in a fluid then we should set topVelocity to
                         * be the P velocity at the top of the layer.
                         */
                        topVelocity = velLayer.evaluateAtTop('P');
                    }
                    topP = toSlowness(topVelocity, velLayer.getTopDepth());
                    if(botP >= p && p >= topP) {
                        return velLayer.getTopDepth();
                    }
                }
            }
            if(Math.abs(p - botP) < slownessTolerance) {
                /*
                 * Check to see if p is just outside the bottommost layer. If so
                 * than return the bottom depth.
                 */
                System.out.println(" p is just outside the bottommost layer."
                        + " This probably shouldn't be allowed to happen!\n");
                return velLayer.getBotDepth();
            }
        } catch(NoSuchMatPropException e) {
            // can't happen...
            e.printStackTrace();
        }
        throw new SlownessModelException("slowness p=" + p
                + " is not contained within the specified layers." + "\np=" + p
                + " topCriticalLayer=" + topCriticalLayer
                + " botCriticalLayer=" + botCriticalLayer + " isPWave="
                + isPWave + " topP=" + topP + " botP=" + botP);
    }

    /**
     * This method takes a velocity model and creates a vector containing
     * slowness-depth layers that, hopefully, adequately sample both slowness
     * and depth so that the travel time as a function of distance can be
     * reconstructed from the theta function. It catches NoSuchLayerException
     * which might be generated in the velocity model. This shouldn't happen
     * though.
     * 
     * @see VelocityModel
     * @exception SlownessModelException
     *                occurs if the validation on the velocity model fails, or
     *                if the velocity model has no layers.
     * @exception NoSuchMatPropException
     *                occurs if wavetype is not recognized.
     */
    private void createSample()
            throws SlownessModelException, NoSuchMatPropException,
            NoSuchLayerException {
        // First check to make sure velocity model is ok.
        if(vMod.validate() == false) {
            throw new SlownessModelException("Error in velocity model!");
        }
        if(vMod.getNumLayers() == 0) {
            throw new SlownessModelException("velModel.getNumLayers()==0");
        }
        if (vMod.getVelocityLayer(0).getTopSVelocity() == 0) {
            throw new SlownessModelException("Unable to handle zero S velocity layers at surface. This should be fixed at some point, but is a limitation of TauP at this point.");
        }
        if(DEBUG) {
            System.out.println("start createSample");
        }
        setRadiusOfEarth(vMod.getRadiusOfEarth());
        if(DEBUG) {
            System.out.println("findCriticalPoints");
        }
        findCriticalPoints();
        if(DEBUG) {
            System.out.println("coarseSample");
        }
        coarseSample();
        boolean isOK = false;
        if(DEBUG) {
            isOK = validate();
            System.out.println("rayParamCheck");
        }
        rayParamIncCheck();
        if(DEBUG) {
            isOK &= validate();
            System.out.println("depthIncCheck");
        }
        depthIncCheck();
        if(DEBUG) {
            isOK &= validate();
            System.out.println("distanceCheck");
        }
        distanceCheck();
        if(DEBUG) {
            isOK &= validate();
            System.out.println("fixCriticalPoints");
        }
        fixCriticalPoints();
        if(DEBUG) {
            System.out.println("done createSample");
        }
    }

    /**
     * Creates a coarse slowness sampling of the velocity model (vMod). The
     * resultant slowness layers will satisfy the maximum depth increments as
     * well as sampling each point specified within the VelocityModel. The P and
     * S sampling will also be compatible.
     */
    protected void coarseSample() throws SlownessModelException,
            NoSuchLayerException, NoSuchMatPropException {
        VelocityLayer prevVLayer;
        VelocityLayer origVLayer;
        VelocityLayer currVLayer;
        SlownessLayer currPLayer, currSLayer;
        PLayers.clear();
        SLayers.clear();
        // to initialize prevVLayer
        origVLayer = vMod.getVelocityLayer(0);
        origVLayer = new VelocityLayer(0,
                                       origVLayer.getTopDepth(),
                                       origVLayer.getTopDepth(),
                                       origVLayer.getTopPVelocity(),
                                       origVLayer.getTopPVelocity(),
                                       origVLayer.getTopSVelocity(),
                                       origVLayer.getTopSVelocity(),
                                       origVLayer.getTopDensity(),
                                       origVLayer.getTopDensity(),
                                       origVLayer.getTopQp(),
                                       origVLayer.getTopQp(),
                                       origVLayer.getTopQs(),
                                       origVLayer.getTopQs());
        for(int layerNum = 0; layerNum < vMod.getNumLayers(); layerNum++) {
            prevVLayer = origVLayer;
            origVLayer = vMod.getVelocityLayer(layerNum);
            /*
             * Check for first order discontinuity. However, we only
             * consider S discontinuities in the inner core if
             * allowInnerCoreS is true.
             */
            if(prevVLayer.getBotPVelocity() != origVLayer.getTopPVelocity()
                    || (prevVLayer.getBotSVelocity() != origVLayer.getTopSVelocity() 
                            && (allowInnerCoreS || origVLayer.getTopDepth() < vMod.getIocbDepth()))) {
                double topSVel, botSVel;
                /*
                 * if we are going from a fluid to a solid or solid to
                 * fluid, ex core mantle or outer core to inner core then we
                 * need to use the P velocity for determining the S
                 * discontinuity.
                 */
                if(prevVLayer.getBotSVelocity() == 0.0) {
                    topSVel = prevVLayer.getBotPVelocity();
                } else {
                    topSVel = prevVLayer.getBotSVelocity();
                }
                if(origVLayer.getTopSVelocity() == 0.0) {
                    botSVel = origVLayer.getTopPVelocity();
                } else {
                    botSVel = origVLayer.getTopSVelocity();
                }
                currVLayer = new VelocityLayer(layerNum,
                                               prevVLayer.getBotDepth(),
                                               prevVLayer.getBotDepth(),
                                               prevVLayer.getBotPVelocity(),
                                               origVLayer.getTopPVelocity(),
                                               topSVel,
                                               botSVel);
                /*
                 * Add the zero thickness, but with nonzero slowness step,
                 * layer corresponding to the discontinuity.
                 */
                currPLayer = toSlownessLayer(currVLayer, PWAVE);
                PLayers.add(currPLayer);
                if((prevVLayer.getBotSVelocity() == 0.0 && origVLayer.getTopSVelocity() == 0.0)
                        || (!allowInnerCoreS && currVLayer.getTopDepth() >= vMod.getIocbDepth())) {
                    currSLayer = currPLayer;
                } else {
                    currSLayer = toSlownessLayer(currVLayer, SWAVE);
                }
                SLayers.add(currSLayer);
            }
            currPLayer = toSlownessLayer(origVLayer, PWAVE);
            PLayers.add(currPLayer);
            if(depthInFluid(origVLayer.getTopDepth())
                    || (!allowInnerCoreS && origVLayer.getTopDepth() >= vMod.getIocbDepth())) {
                currSLayer = currPLayer;
            } else {
                currSLayer = toSlownessLayer(origVLayer, SWAVE);
            }
            SLayers.add(currSLayer);
        }
        // make sure that all high slowness layers are sampled exactly
        // at their bottom
        int highZoneNum, SLayerNum;
        SlownessLayer highSLayer;
        DepthRange highZone;
        for(highZoneNum = 0; highZoneNum < highSlownessLayerDepthsS.size(); highZoneNum++) {
            highZone = highSlownessLayerDepthsS.get(highZoneNum);
            SLayerNum = layerNumberAbove(highZone.botDepth, SWAVE);
            highSLayer = getSlownessLayer(SLayerNum, SWAVE);
            while(highSLayer.getTopDepth() == highSLayer.getBotDepth()
                    && (highSLayer.getTopP() - highZone.rayParam)
                            * (highZone.rayParam - highSLayer.getBotP()) < 0) {
                SLayerNum++;
                highSLayer = getSlownessLayer(SLayerNum, SWAVE);
            }
            if(highZone.rayParam != highSLayer.getBotP()) {
                addSlowness(highZone.rayParam, SWAVE);
            }
        }
        for(highZoneNum = 0; highZoneNum < highSlownessLayerDepthsP.size(); highZoneNum++) {
            highZone = highSlownessLayerDepthsP.get(highZoneNum);
            SLayerNum = layerNumberAbove(highZone.botDepth, PWAVE);
            highSLayer = getSlownessLayer(SLayerNum, PWAVE);
            while(highSLayer.getTopDepth() == highSLayer.getBotDepth()
                    && (highSLayer.getTopP() - highZone.rayParam)
                            * (highZone.rayParam - highSLayer.getBotP()) < 0) {
                SLayerNum++;
                highSLayer = getSlownessLayer(SLayerNum, PWAVE);
            }
            if(highZone.rayParam != highSLayer.getBotP()) {
                addSlowness(highZone.rayParam, PWAVE);
            }
        }
        // make sure P and S sampling are consistant
        double botP = -1;
        double topP = -1;
        for(int j = 0; j < PLayers.size(); j++) {
            topP = PLayers.get(j).getTopP();
            if(topP != botP) {
                addSlowness(topP, SWAVE);
            }
            botP = PLayers.get(j).getBotP();
            addSlowness(botP, SWAVE);
        }
        botP = -1;
        for(int j = 0; j < SLayers.size(); j++) {
            topP = SLayers.get(j).getTopP();
            if(topP != botP) {
                addSlowness(topP, PWAVE);
            }
            botP = SLayers.get(j).getBotP();
            addSlowness(botP, PWAVE);
        }
    }

    /**
     * Checks to make sure that no slowness layer spans more than maxDeltaP.
     */
    protected void rayParamIncCheck() throws SlownessModelException,
            NoSuchLayerException {
        SlownessLayer sLayer;
        double numNewP;
        double deltaP;
        for(int j = 0; j < SLayers.size(); j++) {
            sLayer = SLayers.get(j);
            if(Math.abs(sLayer.getTopP() - sLayer.getBotP()) > maxDeltaP) {
                numNewP = Math.ceil(Math.abs(sLayer.getTopP()
                        - sLayer.getBotP())
                        / maxDeltaP);
                deltaP = (sLayer.getTopP() - sLayer.getBotP()) / numNewP;
                for(int rayNum = 1; rayNum < numNewP; rayNum++) {
                    addSlowness(sLayer.getTopP() + rayNum * deltaP, PWAVE);
                    addSlowness(sLayer.getTopP() + rayNum * deltaP, SWAVE);
                }
            }
        }
        for(int j = 0; j < PLayers.size(); j++) {
            sLayer = PLayers.get(j);
            if(Math.abs(sLayer.getTopP() - sLayer.getBotP()) > maxDeltaP) {
                numNewP = Math.ceil(Math.abs(sLayer.getTopP()
                        - sLayer.getBotP())
                        / maxDeltaP);
                deltaP = (sLayer.getTopP() - sLayer.getBotP()) / numNewP;
                for(int rayNum = 1; rayNum < numNewP; rayNum++) {
                    addSlowness(sLayer.getTopP() + rayNum * deltaP, PWAVE);
                    addSlowness(sLayer.getTopP() + rayNum * deltaP, SWAVE);
                }
            }
        }
    }

    /**
     * Checks to make sure no slowness layer spans more than maxDepthInterval.
     */
    protected void depthIncCheck() throws SlownessModelException,
            NoSuchLayerException {
        SlownessLayer sLayer;
        int numNewDepths;
        double deltaDepth;
        double velocity;
        double p;
        try {
            for(int j = 0; j < SLayers.size(); j++) {
                sLayer = SLayers.get(j);
                if((sLayer.getBotDepth() - sLayer.getTopDepth()) > maxDepthInterval) {
                    numNewDepths = (int)Math.ceil((sLayer.getBotDepth() - sLayer.getTopDepth())
                            / maxDepthInterval);
                    deltaDepth = (sLayer.getBotDepth() - sLayer.getTopDepth())
                            / numNewDepths;
                    for(int depthNum = 1; depthNum < numNewDepths; depthNum++) {
                        velocity = vMod.evaluateAbove(sLayer.getTopDepth()
                                + depthNum * deltaDepth, 'S');
                        if(velocity == 0.0
                                || (!allowInnerCoreS && sLayer.getTopDepth()
                                        + depthNum * deltaDepth >= vMod.getIocbDepth())) {
                            velocity = vMod.evaluateAbove(sLayer.getTopDepth()
                                    + depthNum * deltaDepth, 'P');
                        }
                        p = toSlowness(velocity, sLayer.getTopDepth()
                                + depthNum * deltaDepth);
                        addSlowness(p, PWAVE);
                        addSlowness(p, SWAVE);
                    }
                }
            }
            for(int j = 0; j < PLayers.size(); j++) {
                sLayer = PLayers.get(j);
                if((sLayer.getBotDepth() - sLayer.getTopDepth()) > maxDepthInterval) {
                    numNewDepths = (int)Math.ceil((sLayer.getBotDepth() - sLayer.getTopDepth())
                            / maxDepthInterval);
                    deltaDepth = (sLayer.getBotDepth() - sLayer.getTopDepth())
                            / numNewDepths;
                    for(int depthNum = 1; depthNum < numNewDepths; depthNum++) {
                        p = toSlowness(vMod.evaluateAbove(sLayer.getTopDepth()
                                               + depthNum * deltaDepth, 'P'),
                                       sLayer.getTopDepth() + depthNum
                                               * deltaDepth);
                        addSlowness(p, PWAVE);
                        addSlowness(p, SWAVE);
                    }
                }
            }
        } catch(NoSuchMatPropException e) {
            throw new RuntimeException("can't happen", e);
        }
    }

    /**
     * Checks to make sure no slowness layer spans more than maxRangeInterval
     * and that the (estimated) error due to linear interpolation is less than
     * maxInterpError.
     */
    protected void distanceCheck() throws SlownessModelException,
            NoSuchMatPropException, NoSuchLayerException {
        SlownessLayer sLayer, prevSLayer;
        int j;
        TimeDist prevTD;
        TimeDist currTD;
        TimeDist prevPrevTD;
        boolean isCurrOK;
        boolean isPrevOK;
        boolean currWaveType; // TRUE=P and FALSE=S
        /* do SWAVE and then PWAVE, waveN is ONLY used on the next 2 lines */
        for(int waveN = 0; waveN < 2; waveN++) {
            currWaveType = waveN == 0 ? SWAVE : PWAVE;
            prevPrevTD = null;
            prevTD = null;
            currTD = null;
            isCurrOK = false;
            isPrevOK = false;
            j = 0;
            sLayer = getSlownessLayer(0, currWaveType); // preset sLayer so
            // prevSLayer is ok
            while(j < getNumLayers(currWaveType)) {
                prevSLayer = sLayer;
                sLayer = getSlownessLayer(j, currWaveType);
                if(!depthInHighSlowness(sLayer.getBotDepth(),
                                        sLayer.getBotP(),
                                        currWaveType)
                        && !depthInHighSlowness(sLayer.getTopDepth(),
                                                sLayer.getTopP(),
                                                currWaveType)) {
                    // Don't calculate prevTD if we can avoid it
                    if(isCurrOK) {
                        if(isPrevOK) {
                            prevPrevTD = prevTD;
                        } else {
                            prevPrevTD = null;
                        }
                        prevTD = currTD;
                        isPrevOK = true;
                    } else {
                        prevTD = approxDistance(j - 1,
                                                sLayer.getTopP(),
                                                currWaveType);
                        isPrevOK = true;
                    }
                    currTD = approxDistance(j, sLayer.getBotP(), currWaveType);
                    isCurrOK = true;
                    // check for too great of distance jump
                    if(Math.abs(prevTD.distRadian - currTD.distRadian) > maxRangeInterval
                            && Math.abs(sLayer.getTopP() - sLayer.getBotP()) > 2 * minDeltaP) {
                        if (DEBUG) {
                            System.out.println(" " + j+"Dist jump too great: "+Math.abs(prevTD.distRadian - currTD.distRadian)+" > "+maxRangeInterval
                                               +"  adding slowness: "+(sLayer.getTopP() + sLayer.getBotP()) / 2.0);
                        }
                        addSlowness((sLayer.getTopP() + sLayer.getBotP()) / 2.0,
                                    PWAVE);
                        addSlowness((sLayer.getTopP() + sLayer.getBotP()) / 2.0,
                                    SWAVE);
                        currTD = prevTD;
                        prevTD = prevPrevTD;
                    } else {
                        // make guess as to error estimate due to linear
                        // interpolation
                        // if it is not ok, then we split both the previous and
                        // current
                        // slowness layers, this has the nice, if unintended,
                        // consequense
                        // of adding extra samples in the neighborhood of poorly
                        // sampled
                        // caustics
                        double splitRayParam = (sLayer.getTopP()+sLayer.getBotP())/2;
                        TimeDist allButLayer = approxDistance(j-1, splitRayParam, currWaveType);
                        SlownessLayer splitLayer = new SlownessLayer(sLayer.getTopP(), sLayer.getTopDepth(), splitRayParam, sLayer.bullenDepthFor(splitRayParam, getRadiusOfEarth()));
                        TimeDist justLayer = splitLayer.bullenRadialSlowness(splitRayParam, getRadiusOfEarth());
                        TimeDist splitTD = new TimeDist(splitRayParam, allButLayer.time+2*justLayer.time, allButLayer.distRadian+2*justLayer.distRadian);
                        //                        if(Math.abs(prevTD.time
//                                    - ((currTD.time - prevPrevTD.time)
//                                            * (prevTD.dist - prevPrevTD.dist)
//                                            / (currTD.dist - prevPrevTD.dist) + prevPrevTD.time)) > maxInterpError) {
                        if(Math.abs(currTD.time
                                        - ((splitTD.time - prevTD.time)
                                                * (currTD.distRadian - prevTD.distRadian)
                                                / (splitTD.distRadian - prevTD.distRadian) + prevTD.time)) > maxInterpError) {

                            if(DEBUG ) {
                                System.out.print(" " + j+" add slowness "+Math.abs(currTD.time
                                                                                   - ((splitTD.time - prevTD.time)
                                                                                           * (currTD.distRadian - prevTD.distRadian)
                                                                                           / (splitTD.distRadian - prevTD.distRadian) + prevTD.time))+" > "+maxInterpError);
                            }
                            addSlowness((prevSLayer.getTopP() + prevSLayer.getBotP()) / 2.0,
                                        PWAVE);
                            addSlowness((prevSLayer.getTopP() + prevSLayer.getBotP()) / 2.0,
                                        SWAVE);
                            addSlowness((sLayer.getTopP() + sLayer.getBotP()) / 2.0,
                                        PWAVE);
                            addSlowness((sLayer.getTopP() + sLayer.getBotP()) / 2.0,
                                        SWAVE);
                            currTD = prevPrevTD;
                            isPrevOK = false;
                            if (j>0) {
                                // back up one step unless we are at beginning, then stay put
                                j--;
                                sLayer = getSlownessLayer(((j - 1 >= 0) ? j - 1 : 0),
                                                          currWaveType);
                                // ^^^ make sure j != 0
                                // this sLayer will become prevSLayer in next loop
                            } else {
                                isPrevOK = false;
                                isCurrOK = false;
                            }
                        } else {
                            j++;
                            if(DEBUG && (j % 10 == 0)) {
                                System.out.println(j);
                            }
                        }
                    }
                } else {
                    prevPrevTD = null;
                    prevTD = null;
                    currTD = null;
                    isCurrOK = false;
                    isPrevOK = false;
                    j++;
                    if(DEBUG && (j % 100 == 0)) {
                        System.out.print(" " + j);
                    }
                }
            }
            if(DEBUG) {
                System.out.println("\nNumber of " + (currWaveType ? 'P' : 'S')
                        + " slowness layers: " + j);
            }
        }
    }

    /**
     * Adds the given ray parameter, p, to the slowness sampling for the given
     * waveType. It splits slowness layers as needed and keeps P and S sampling
     * consistant within fluid layers. Note, this makes use of the velocity
     * model, so all interpolation is linear in velocity, not in slowness!
     * 
     */
    protected void addSlowness(double p, boolean isPWave)
            throws SlownessModelException, NoSuchLayerException {
        List<SlownessLayer> layers, otherLayers;
        SlownessLayer sLayer, topLayer, botLayer;
        double slope;
        double topVelocity, botVelocity;
        int otherIndex;
        if(isPWave) {
            layers = PLayers;
            otherLayers = SLayers;
        } else {
            layers = SLayers;
            otherLayers = PLayers;
        }
        for(int i = 0; i < layers.size(); i++) {
            sLayer = layers.get(i);
            try {
                if(sLayer.getTopDepth() != sLayer.getBotDepth()) {
                    topVelocity = vMod.evaluateBelow(sLayer.getTopDepth(),
                                                     (isPWave ? 'P' : 'S'));
                    botVelocity = vMod.evaluateAbove(sLayer.getBotDepth(),
                                                     (isPWave ? 'P' : 'S'));
                } else {
                    // if depths are same we really only need topVelocity,
                    // and just to verify that we are not in a fluid.
                    topVelocity = vMod.evaluateAbove(sLayer.getBotDepth(),
                                                     (isPWave ? 'P' : 'S'));
                    botVelocity = vMod.evaluateBelow(sLayer.getTopDepth(),
                                                     (isPWave ? 'P' : 'S'));
                }
            } catch(NoSuchMatPropException e) {
                // Can't happen but...
                throw new SlownessModelException("Caught NoSuchMatPropException: "
                        + e.getMessage());
            }
            // We don't need to check for S waves in a fluid or
            // in inner core if allowInnerCoreS==false.
            if(!isPWave) {
                if(!allowInnerCoreS
                        && sLayer.getBotDepth() > vMod.getIocbDepth()) {
                    break;
                } else if(topVelocity == 0.0) {
                    continue;
                }
            }
            if((sLayer.getTopP() - p) * (p - sLayer.getBotP()) > 0) {
                double botDepth = sLayer.getBotDepth();
                if(sLayer.getBotDepth() != sLayer.getTopDepth()) {
                    /*
                     * not a zero thickness layer, so calculate the depth for
                     * the ray parameter.
                     */
                    slope = (botVelocity - topVelocity)
                            / (sLayer.getBotDepth() - sLayer.getTopDepth());
                    botDepth = interpolate(p,
                                                     topVelocity,
                                                     sLayer.getTopDepth(),
                                                     slope);
                }
                botLayer = new SlownessLayer(p, botDepth, sLayer.getBotP(), sLayer.getBotDepth() );
                topLayer = new SlownessLayer(sLayer.getTopP(), sLayer.getTopDepth(), p, botDepth);
                layers.remove(i);
                layers.add(i, botLayer);
                layers.add(i, topLayer);
                otherIndex = otherLayers.indexOf(sLayer);
                if(otherIndex != -1) {
                    otherLayers.remove(otherIndex);
                    otherLayers.add(otherIndex, botLayer);
                    otherLayers.add(otherIndex, topLayer);
                }
            }
        }
    }

    /**
     * Resets the slowness layers that correspond to critical points.
     */
    protected void fixCriticalPoints() throws NoSuchLayerException {
        CriticalDepth cd;
        SlownessLayer sLayer;
        for(int i = 0; i < criticalDepths.size(); i++) {
            cd = criticalDepths.get(i);
            cd.setPLayerNum(layerNumberBelow(cd.getDepth(), PWAVE));
            sLayer = getSlownessLayer(cd.getPLayerNum(), PWAVE);
            if(cd.getPLayerNum() == PLayers.size() - 1
                    && sLayer.getBotDepth() == cd.getDepth()) {
                cd.setPLayerNum(cd.getPLayerNum() + 1); // want the last
                // critical point to be
                // the bottom of the
                // last layer
            }
            cd.setSLayerNum(layerNumberBelow(cd.getDepth(), SWAVE));
            sLayer = getSlownessLayer(cd.getSLayerNum(), SWAVE);
            if(cd.getSLayerNum() == SLayers.size() - 1
                    && sLayer.getBotDepth() == cd.getDepth()) {
                cd.setSLayerNum(cd.getSLayerNum() + 1); // want the last
                // critical point to be
                // the bottom of the
                // last layer
            }
        }
    }

    /** finds a layer that contains the depth. This may not be unique in the case of a depth on
     * a boundary in the velocity model due to zero thickness layers. If the uppermost or
     * lowermost layer containing the depth is needed, use layerNumberAbove() or layerNumberBelow().
     */
    public int layerNumForDepth(double depth, boolean isPWave)  throws NoSuchLayerException {
        SlownessLayer tempLayer;
        List<SlownessLayer> layers;
        if(isPWave) {
            layers = PLayers;
        } else {
            layers = SLayers;
        }
        // check to make sure depth is within the range available
        if(depth < layers.get(0).getTopDepth()
                || layers.get(layers.size() - 1).getBotDepth() < depth) {
            throw new NoSuchLayerException(depth);
        }
        int tooSmallNum = 0;
        int tooLargeNum = layers.size() - 1;
        int currentNum = 0;
        while(true) {
            if (tooLargeNum - tooSmallNum < 3) {
                // end of Newton, just check
                for (currentNum = tooSmallNum; currentNum <= tooLargeNum; currentNum++) {
                    tempLayer = getSlownessLayer(currentNum, isPWave);
                    if (tempLayer.containsDepth(depth)) {
                        return currentNum;
                    }
                }
            } else {
                currentNum = Math.round((tooSmallNum + tooLargeNum) / 2.0f);
            }
            tempLayer = getSlownessLayer(currentNum, isPWave);
            if(tempLayer.getTopDepth() > depth) {
                tooLargeNum = currentNum - 1;
            } else if(tempLayer.getBotDepth() < depth) {
                tooSmallNum = currentNum + 1;
            } else {
                return currentNum;
            }
            if (tooSmallNum > tooLargeNum) {
                throw new RuntimeException("tooSmallNum ("+tooSmallNum+") >= tooLargeNum ("+tooLargeNum+")");
            }
        }
    }
    
    /**
     * Finds the index of the slowness layer that contains the given depth Note
     * that if the depth is a layer boundary, it returns the shallower of the
     * two or possibly more (since total reflections are zero thickness layers)
     * layers.
     * 
     * @return the layer number.
     * @exception NoSuchLayerException
     *                occurs if no layer in the slowness model contains the
     *                given depth.
     */
    public int layerNumberAbove(double depth, boolean isPWave)
            throws NoSuchLayerException {
        int foundLayerNum = layerNumForDepth(depth, isPWave);
        SlownessLayer tempLayer = getSlownessLayer(foundLayerNum, isPWave);
        while(tempLayer.getTopDepth() == depth && foundLayerNum > 0) {
            foundLayerNum--;
            tempLayer = getSlownessLayer(foundLayerNum, isPWave);
        }
        return foundLayerNum;
    }

    /**
     * Finds the index of the slowness layer that contains the given depth Note
     * that if the depth is a layer boundary, it returns the deeper of the two
     * or possibly more (since total reflections are zero thickness layers)
     * layers.
     * 
     * @return the layer number.
     * @exception NoSuchLayerException
     *                occurs if no layer in the slowness model contains the
     *                given depth.
     */
    public int layerNumberBelow(double depth, boolean isPWave)
            throws NoSuchLayerException {
        int foundLayerNum = layerNumForDepth(depth, isPWave);
        List<SlownessLayer> layers;
        if(isPWave) {
            layers = PLayers;
        } else {
            layers = SLayers;
        }
        SlownessLayer tempLayer = getSlownessLayer(foundLayerNum, isPWave);
        while(tempLayer.getBotDepth() == depth && foundLayerNum < layers.size()-1) {
            foundLayerNum++;
            tempLayer = getSlownessLayer(foundLayerNum, isPWave);
        }
        return foundLayerNum;
    }

    /**
     * Performs consistency check on the slowness model.
     * 
     * @return true if successful, throws SlownessModelException otherwise.
     * @exception SlownessModelException
     *                if any check fails
     */
    public boolean validate() throws SlownessModelException {
        boolean isOK = true;
        double prevDepth;
        DepthRange highSZoneDepth, fluidZone;
        /* is radiusOfEarth positive? */
        if(radiusOfEarth <= 0.0) {
            throw new SlownessModelException("Radius of earth is not positive. radiusOfEarth = "
                    + radiusOfEarth);
        }
        /* is maxDepthInterval positive? */
        if(maxDepthInterval <= 0.0) {
            throw new SlownessModelException("maxDepthInterval is not positive. maxDepthInterval = "
                    + maxDepthInterval);
        }
        /* Check for inconsistencies in high slowness zones. */
        List<DepthRange> highSlownessLayerDepths = highSlownessLayerDepthsP;
        boolean isPWave = PWAVE;
        for(int j = 0; j < 2; j++, isPWave = SWAVE) {
            if(isPWave) {
                highSlownessLayerDepths = highSlownessLayerDepthsP;
            } else {
                highSlownessLayerDepths = highSlownessLayerDepthsS;
            }
            prevDepth = -1 * Double.MAX_VALUE;
            for(int i = 0; i < highSlownessLayerDepths.size(); i++) {
                highSZoneDepth = highSlownessLayerDepths.get(i);
                if(highSZoneDepth.topDepth >= highSZoneDepth.botDepth) {
                    throw new SlownessModelException("High slowness zone has zero or negative thickness. Num "
                            + i
                            + " isPWave="
                            + isPWave
                            + " top depth "
                            + highSZoneDepth.topDepth
                            + " bottom depth "
                            + highSZoneDepth.botDepth);
                }
                if(highSZoneDepth.topDepth <= prevDepth) {
                    throw new SlownessModelException("High slowness zone overlaps previous zone. Num "
                            + i
                            + " isPWave="
                            + isPWave
                            + " top depth "
                            + highSZoneDepth.topDepth
                            + " bottom depth "
                            + highSZoneDepth.botDepth);
                }
                prevDepth = highSZoneDepth.botDepth;
            }
        }
        /* Check for inconsistencies in fluid zones. */
        prevDepth = -1 * Double.MAX_VALUE;
        for(int i = 0; i < fluidLayerDepths.size(); i++) {
            fluidZone = fluidLayerDepths.get(i);
            if(fluidZone.topDepth >= fluidZone.botDepth) {
                throw new SlownessModelException("Fluid zone has zero or negative thickness. Num "
                        + i
                        + " top depth "
                        + fluidZone.topDepth
                        + " bottom depth " + fluidZone.botDepth);
            }
            if(fluidZone.topDepth <= prevDepth) {
                throw new SlownessModelException("Fluid zone overlaps previous zone. Num "
                        + i
                        + " top depth "
                        + fluidZone.topDepth
                        + " bottom depth " + fluidZone.botDepth);
            }
            prevDepth = fluidZone.botDepth;
        }
        /* Check for inconsistencies in slowness layers. */
        isPWave = PWAVE;
        double prevBotP;
        for(int j = 0; j < 2; j++, isPWave = SWAVE) {
            prevDepth = 0.0;
            SlownessLayer sLayer = null;
            if(getNumLayers(isPWave) > 0) {
                sLayer = getSlownessLayer(0, isPWave);
                prevBotP = sLayer.getTopP();
            } else {
                prevBotP = -1;
            }
            SlownessLayer prevSLayer = null;
            for(int i = 0; i < getNumLayers(isPWave); i++) {
                sLayer = getSlownessLayer(i, isPWave);
                isOK &= sLayer.validate();
                if(sLayer.getTopDepth() > prevDepth) {
                    throw new SlownessModelException("Gap of "
                            + (sLayer.getTopDepth() - prevDepth)
                            + " between slowness layers. Num " + i
                            + " isPWave=" + isPWave + "  top "
                            + prevSLayer+ " bottom "
                            + sLayer);
                }
                if(sLayer.getTopDepth() < prevDepth) {
                    throw new SlownessModelException("Slowness layer overlaps previous layer by "
                            + (prevDepth - sLayer.getTopDepth())
                            + ". Num "
                            + i
                            + " isPWave="
                            + isPWave
                            + " top depth "
                            + sLayer.getTopDepth()
                            + " bottom depth "
                            + sLayer.getBotDepth());
                }
                if(sLayer.getTopP() != prevBotP) {
                    throw new SlownessModelException("Slowness layer gap/overlaps previous layer in slowness "
                            + ". Num "
                            + i
                            + " isPWave="
                            + isPWave
                            + " prevBotP= " + prevBotP + " prevSLayer= " + prevSLayer+ " sLayer= " + sLayer);
                }
                if(Double.isNaN(sLayer.getTopDepth())) {
                    throw new SlownessModelException("Top depth is NaN, layerNum="
                            + i + " waveType=" + (isPWave ? 'P' : 'S'));
                }
                if(Double.isNaN(sLayer.getBotDepth())) {
                    throw new SlownessModelException("Top depth is NaN, layerNum="
                            + i + " waveType=" + (isPWave ? 'P' : 'S'));
                }
                prevSLayer = sLayer;
                prevBotP = sLayer.getBotP();
                prevDepth = sLayer.getBotDepth();
            }
        }
        /* Everything checks out OK so return true. */
        return isOK;
    }

    public String toString() {
        int topCriticalLayerNum;
        int botCriticalLayerNum;
        String desc = "";
        desc = "radiusOfEarth=" + radiusOfEarth + "\n maxDeltaP=" + maxDeltaP
                + "\n minDeltaP=" + minDeltaP + "\n maxDepthInterval="
                + maxDepthInterval + "\n maxRangeInterval=" + maxRangeInterval
                + "\n allowInnerCoreS=" + allowInnerCoreS
                + "\n slownessTolerance=" + slownessTolerance
                + "\n getNumLayers('P')=" + getNumLayers(PWAVE)
                + "\n getNumLayers('S')=" + getNumLayers(SWAVE)
                + "\n fluidLayerDepths.size()=" + fluidLayerDepths.size()
                + "\n highSlownessLayerDepthsP.size()="
                + highSlownessLayerDepthsP.size()
                + "\n highSlownessLayerDepthsS.size()="
                + highSlownessLayerDepthsS.size()
                + "\n criticalDepths.size()=" + criticalDepths.size()
                + "\n";
        if(criticalDepths.size() != 0) {
            desc += ("**** Critical Depth Layers ************************\n");
            botCriticalLayerNum = criticalDepths.get(0).getVelLayerNum() - 1;
            for(int criticalNum = 1; criticalNum < criticalDepths.size(); criticalNum++) {
                topCriticalLayerNum = botCriticalLayerNum + 1;
                botCriticalLayerNum = criticalDepths.get(criticalNum).getVelLayerNum() - 1;
                desc += " " + topCriticalLayerNum + "," + botCriticalLayerNum;
            }
        }
        desc += "\n";
        if(fluidLayerDepths.size() != 0) {
            desc += "\n**** Fluid Layer Depths ************************\n";
            for(int i = 0; i < fluidLayerDepths.size(); i++) {
                desc += fluidLayerDepths.get(i).topDepth
                        + ","
                        + fluidLayerDepths.get(i).botDepth
                        + " ";
            }
        }
        desc += "\n";
        if(highSlownessLayerDepthsP.size() != 0) {
            desc += "\n**** P High Slowness Layer Depths ****************\n";
            for(int i = 0; i < highSlownessLayerDepthsP.size(); i++) {
                desc += highSlownessLayerDepthsP.get(i).topDepth
                        + ","
                        + highSlownessLayerDepthsP.get(i).botDepth
                        + " ";
            }
        }
        desc += "\n";
        if(highSlownessLayerDepthsS.size() != 0) {
            desc += "\n**** S High Slowness Layer Depths ****************\n";
            for(int i = 0; i < highSlownessLayerDepthsS.size(); i++) {
                desc += highSlownessLayerDepthsS.get(i).topDepth
                        + ","
                        + highSlownessLayerDepthsS.get(i).botDepth
                        + " ";
            }
        }
        desc += "\n";
        desc += "\n**** P Layers ****************\n";
        for (SlownessLayer l : PLayers) {
            desc+=l.toString()+"\n";
        }
        return desc;
    }
}
