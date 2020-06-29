/*
 * <pre> The TauP Toolkit: Flexible Seismic Travel-Time and Raypath Utilities.
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
 * owens@seis.sc.edu </pre>
 */
package edu.sc.seis.TauP;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OptionalDataException;
import java.io.Serializable;
import java.io.StreamCorruptedException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * Stores and transforms seismic phase names to and from their corresponding
 * sequence of branches.
 * 
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * @author H. Philip Crotwell
 * 
 * Modified to add "expert" mode wherein paths may start in the core. Principal
 * use is to calculate leg contributions for scattered phases. Nomenclature: "K" -
 * downgoing wave from source in core; "k" - upgoing wave from source in core.
 * 
 * G. Helffrich/U. Bristol 24 Feb. 2007
 */
public class SeismicPhase implements Serializable, Cloneable {

    /** Enables debugging output. */
    public transient boolean DEBUG = TauP_Time.DEBUG;

    /** Enables verbose output. */
    public transient boolean verbose = false;

    /** Enables phases originating in core. */
    public static transient boolean expert = TauP_Time.expert;

    /** TauModel to generate phase for. */
    protected TauModel tMod;

    /**
     * Used by addToBranch when the path turns within a segment. We assume that
     * no ray will turn downward so turning implies turning from downward to
     * upward, ie U.
     */
    public static final int TURN = 0;

    /**
     * Used by addToBranch when the path reflects off the top of the end of a
     * segment, ie ^.
     */
    public static final int REFLECTTOP = 1;

    /**
     * Used by addToBranch when the path reflects off the bottom of the end of a
     * segment, ie v.
     */
    public static final int REFLECTBOT = 2;

    /**
     * Used by addToBranch when the path transmits up through the end of a
     * segment.
     */
    public static final int TRANSUP = 3;

    /**
     * Used by addToBranch when the path transmits down through the end of a
     * segment.
     */
    public static final int TRANSDOWN = 4;

    /**
     * The maximum degrees that a Pn or Sn can refract along the moho. Note this
     * is not the total distance, only the segment along the moho. The default
     * is 20 degrees.
     */
    protected static double maxRefraction = 20;

    /**
     * The maximum degrees that a Pdiff or Sdiff can diffract along the CMB.
     * Note this is not the total distance, only the segment along the CMB. The
     * default is 60 degrees.
     */
    protected static double maxDiffraction = 60;

    /**
     * The source depth within the TauModel that was used to generate this
     * phase.
     */
    protected double sourceDepth;

    /**
     * Array of distances corresponding to the ray parameters stored in
     * rayParams.
     */
    protected double[] dist = new double[0];

    /**
     * Array of times corresponding to the ray parameters stored in rayParams.
     */
    protected double[] time = new double[0];

    /** Array of possible ray parameters for this phase. */
    protected double[] rayParams = new double[0];

    /** Minimum ray parameter that exists for this phase. */
    protected double minRayParam;

    /** Maximum ray parameter that exists for this phase. */
    protected double maxRayParam;

    /**
     * Index within TauModel.rayParams that corresponds to maxRayParam. Note
     * that maxRayParamIndex < minRayParamIndex as ray parameter decreases with
     * increasing index.
     */
    protected int maxRayParamIndex = -1;

    /**
     * Index within TauModel.rayParams that corresponds to minRayParam. Note
     * that maxRayParamIndex < minRayParamIndex as ray parameter decreases with
     * increasing index.
     */
    protected int minRayParamIndex = -1;

    /** The minimum distance that this phase can be theoretically observed. */
    protected double minDistance = 0.0;

    /** The maximum distance that this phase can be theoretically observed. */
    protected double maxDistance = Double.MAX_VALUE;

    /**
     * Array of branch numbers for the given phase. Note that this depends upon
     * both the earth model and the source depth.
     */
    protected ArrayList<Integer> branchSeq = new ArrayList<Integer>();

    /** The phase name, ie PKiKP. */
    protected String name;

    /**
     * name with depths corrected to be actuall discontinuities in the model.
     */
    protected String puristName;

    /** ArrayList containing Strings for each leg. */
    protected ArrayList<String> legs = new ArrayList<String>();

    /**
     * temporary branch number so we know where to start add to the branch
     * sequence. Used in addToBranch() and parseName().
     */
    protected transient int currBranch;

    /**
     * temporary end action so we know what we did at the end of the last
     * section of the branch sequence. Used in addToBranch() and parseName().
     */
    protected transient int endAction;

    /**
     * records the end action for the current leg. Will be one of
     * SeismicPhase.TURN, SeismicPhase.TRANSDOWN, SeismicPhase.TRANSUP,
     * SeismicPhase.REFLECTBOT, or SeismicPhase.REFLECTTOP. This allows a check
     * to make sure the path is correct. Used in addToBranch() and parseName().
     */
    protected ArrayList<Integer> legAction = new ArrayList<Integer>();

    /**
     * true if the current leg of the phase is down going. This allows a check
     * to make sure the path is correct. Used in addToBranch() and parseName().
     */
    protected ArrayList<Boolean> downGoing = new ArrayList<Boolean>();

    /**
     * ArrayList of wave types corresponding to each leg of the phase.
     * 
     * @see legs
     */
    protected ArrayList<Boolean> waveType = new ArrayList<Boolean>();

    public static final boolean PWAVE = true;

    public static final boolean SWAVE = false;

    public SeismicPhase(String name, String modelName, double depth) throws TauModelException {
        this(name, TauModelLoader.load(modelName).depthCorrect(depth));
    }
    /**
     * @param phaseName
     *            String containing a name of the phase.
     * @param tMod
     *            Tau model to be used to construct the phase. This should be corrected for the source
     *            depth.
     * @throws TauModelException 
     */
    public SeismicPhase(String name, TauModel tMod) throws TauModelException {
        this.name = name;
        this.sourceDepth = tMod.getSourceDepth();
        this.tMod = tMod;
        legs = legPuller(name);
        createPuristName(tMod);
        parseName(tMod);
        sumBranches(tMod);
    }
    
    public Arrival getEarliestArrival(double degrees) {
        double soonest = 999999999.0;
        Arrival soonestArrival = null;
        List<Arrival> arrivals = calcTime(degrees);
        for (Arrival a : arrivals) {
            if (a.getTime() < soonest) {
                soonestArrival = a;
                soonest = a.getTime();
            }
        }
        return soonestArrival;
    }

    public TauModel getTauModel() {
        return tMod;
    }

    public double getMinDistanceDeg() {
        return getMinDistance() * 180.0 / Math.PI;
    }
    
    public double getMinDistance() {
        return minDistance;
    }

    public double getMaxDistanceDeg() {
        return getMaxDistance() * 180.0 / Math.PI;
    }
    
    public double getMaxDistance() {
        return maxDistance;
    }

    public double getMaxRayParam() {
        return maxRayParam;
    }

    public double getMinRayParam() {
        return minRayParam;
    }

    public int getMaxRayParamIndex() {
        return maxRayParamIndex;
    }

    public int getMinRayParamIndex() {
        return minRayParamIndex;
    }

    public static double getMaxRefraction() {
        return maxRefraction;
    }

    public static void setMaxRefraction(double max) {
        maxRefraction = max;
    }

    public static double getMaxDiffraction() {
        return maxDiffraction;
    }

    public static void setMaxDiffraction(double max) {
        maxDiffraction = max;
    }

    public String getName() {
        return name;
    }

    public String getPuristName() {
        return name;
    }

    public List<String> getLegs() {
        return Collections.unmodifiableList(legs);
    }

    public double[] getRayParams() {
        return (double[])rayParams.clone();
    }

    public double[] getDist() {
        return (double[])dist.clone();
    }

    public double[] getTime() {
        return (double[])time.clone();
    }

    public double[] getTau() {
        double[] tau = new double[dist.length];
        for(int i = 0; i < dist.length; i++) {
            tau[i] = time[i] - rayParams[i] * dist[i];
        }
        return tau;
    }

    /**
     * Direction of the leg between pierce point i and i+1, true is downgoing,
     * false if upgoing
     */
    public boolean[] getDownGoing() {
        Boolean[] b = (Boolean[])downGoing.toArray(new Boolean[0]);
        boolean[] out = new boolean[b.length];
        for(int i = 0; i < b.length; i++) {
            out[i] = b[i].booleanValue();
        }
        return out;
    }

    /**
     * Wave type of the leg between pierce point i and i+1, true is P, false if
     * S
     */
    public boolean[] getWaveType() {
        Boolean[] b = (Boolean[])waveType.toArray(new Boolean[0]);
        boolean[] out = new boolean[b.length];
        for(int i = 0; i < b.length; i++) {
            out[i] = b[i].booleanValue();
        }
        return out;
    }

    /**
     * Leg type i layer interaction, one of TURN, REFLECTTOP, REFLECTBOT,
     * TRANSUP, TRANSDOWN
     */
    public int[] getLegAction() {
        Integer[] b = (Integer[])legAction.toArray(new Integer[0]);
        int[] out = new int[b.length];
        for(int i = 0; i < b.length; i++) {
            out[i] = b[i].intValue();
        }
        return out;
    }
    
    public boolean hasArrivals() {
        return dist != null && dist.length != 0;
    }

    // Normal methods

    /** calculates arrival times for this phase, sorted by time. 
     * @throws NoSuchMatPropException 
     * @throws NoSuchLayerException */
    public List<Arrival> calcTime(double deg) {
        double tempDeg = deg;
        if(tempDeg < 0.0) {
            tempDeg *= -1.0;
        } // make sure deg is positive
        while(tempDeg > 360.0) {
            tempDeg -= 360.0;
        } // make sure it is less than 360
        if(tempDeg > 180.0) {
            tempDeg = 360.0 - tempDeg;
        } // make sure less than or equal to 180
        // now we have 0.0 <= deg <= 180
        double radDist = tempDeg * Math.PI / 180.0;
        List<Arrival> arrivals = new ArrayList<Arrival>();
        /*
         * Search all distances 2n*PI+radDist and 2(n+1)*PI-radDist that are
         * less than the maximum distance for this phase. This insures that we
         * get the time for phases that accumulate more than 180 degrees of
         * distance, for instance PKKKKP might wrap all of the way around. A
         * special case exists at 180, so we skip the second case if
         * tempDeg==180.
         */
        int n = 0;
        double searchDist;
        while(n * 2.0 * Math.PI + radDist <= maxDistance) {
            /*
             * Look for arrivals that are radDist + 2nPi, ie rays that have done
             * more than n laps.
             */
            searchDist = n * 2.0 * Math.PI + radDist;
            for(int rayNum = 0; rayNum < (dist.length - 1); rayNum++) {
                if(searchDist == dist[rayNum + 1]
                        && rayNum + 1 != dist.length - 1) {
                    /* So we don't get 2 arrivals for the same ray. */
                    continue;
                } else if((dist[rayNum] - searchDist)
                        * (searchDist - dist[rayNum + 1]) >= 0.0) {
                    /* look for distances that bracket the search distance */
                    if((rayParams[rayNum] == rayParams[rayNum + 1])
                            && rayParams.length > 2) {
                        /*
                         * Here we have a shadow zone, so it is not really an
                         * arrival.
                         */
                        continue;
                    }
                    if(DEBUG) {
                        System.err.println("SeismicPhase " + name
                                + ", found arrival:\n" + "dist "
                                + (float)(180 / Math.PI * dist[rayNum]) + " "
                                + (float)(180 / Math.PI * searchDist) + " "
                                + (float)(180 / Math.PI * dist[rayNum + 1]));
                    }
                    arrivals.add(linearInterpArrival(searchDist, rayNum, name, puristName, sourceDepth));
                }
            }
            /*
             * Look for arrivals that are 2(n+1)Pi-radDist, ie rays that have
             * done more than one half lap plus some number of whole laps.
             */
            searchDist = (n + 1) * 2.0 * Math.PI - radDist;
            if(tempDeg != 180) {
                for(int rayNum = 0; rayNum < (dist.length - 1); rayNum++) {
                    if(searchDist == dist[rayNum + 1]
                            && rayNum + 1 != dist.length - 1) {
                        /* So we don't get 2 arrivals for the same ray. */
                        continue;
                    } else if((dist[rayNum] - searchDist)
                            * (searchDist - dist[rayNum + 1]) >= 0.0) {
                        if((rayParams[rayNum] == rayParams[rayNum + 1])
                                && rayParams.length > 2) {
                            /*
                             * Here we have a shadow zone, so it is not really
                             * an arrival.
                             */
                            continue;
                        }
                        if(DEBUG) {
                            System.err.println("SeismicPhase " + name
                                    + ", found arrival:\n" + "dist "
                                    + (float)(180 / Math.PI * dist[rayNum])
                                    + " " + (float)(180 / Math.PI * searchDist)
                                    + " "
                                    + (float)(180 / Math.PI * dist[rayNum + 1]));
                        }
                        arrivals.add(linearInterpArrival(searchDist, rayNum, name, puristName, sourceDepth));
                        
                    }
                }
            }
            n++;
        }
        Collections.sort(arrivals, new Comparator<Arrival>() {
            public int compare(Arrival o1, Arrival o2) {
                return Double.compare(o1.getTime(), o2.getTime());
            }});
        return arrivals;
    }
    
    private Arrival linearInterpArrival(double searchDist,
                                        int rayNum,
                                        String name,
                                        String puristName,
                                        double sourceDepth) {
        double arrivalTime = (searchDist - dist[rayNum])
                / (dist[rayNum + 1] - dist[rayNum])
                * (time[rayNum + 1] - time[rayNum]) + time[rayNum];
        double arrivalRayParam = (searchDist - dist[rayNum + 1])
                * (rayParams[rayNum] - rayParams[rayNum + 1])
                / (dist[rayNum] - dist[rayNum + 1])
                + rayParams[rayNum + 1];
        double takeoffVelocity;
        double takeoffAngle = -1;
        double incidentAngle = -1;
        if (name.endsWith("kmps")) {
            takeoffAngle = 0;
            incidentAngle = 0;
        } else {
            VelocityModel vMod = getTauModel().getVelocityModel();
            try {
                if (getDownGoing()[0]) {
                    takeoffVelocity = vMod.evaluateBelow(sourceDepth, name.charAt(0));
                } else { 
                    //fake neg velocity so angle is neg in case of upgoing
                    takeoffVelocity = -1*vMod.evaluateAbove(sourceDepth, name.charAt(0));
                }
                takeoffAngle = 180/Math.PI*Math.asin(takeoffVelocity*arrivalRayParam/(getTauModel().getRadiusOfEarth()-sourceDepth));
                char lastLeg = getLegs().get(getLegs().size()-2).charAt(0); // last item is "E", assume first char is P or S
                incidentAngle = 180/Math.PI*Math.asin(vMod.evaluateBelow(0, lastLeg)*arrivalRayParam/getTauModel().getRadiusOfEarth());
            } catch(NoSuchLayerException e) {
                throw new RuntimeException("Should not happen", e);
            } catch(NoSuchMatPropException e) {
                throw new RuntimeException("Should not happen", e);
            }
        }
        return new Arrival(this,
                                 arrivalTime,
                                 searchDist,
                                 arrivalRayParam,
                                 rayNum,
                                 name,
                                 puristName,
                                 sourceDepth,
                                 takeoffAngle,
                                 incidentAngle);
    }

    /**
     * changes maxRayParam and minRayParam whenever there is a phase conversion.
     * For instance, SKP needs to change the maxRayParam because there are SKS
     * ray parameters that cannot propagate from the cmb into the mantle as a p
     * wave.
     */
    protected void phaseConversion(TauModel tMod,
                                   int fromBranch,
                                   int endAction,
                                   boolean isPtoS) throws TauModelException {
        if(endAction == TURN) {
            // can't phase convert for just a turn point
            throw new TauModelException("Illegal endAction: endAction="
                    + endAction
                    + "\nphase conversion are not allowed at turn points.");
        } else if(endAction == REFLECTTOP) {
            maxRayParam = Math.min(maxRayParam, tMod.getTauBranch(fromBranch,
                                                                  isPtoS)
                    .getMaxRayParam());
            maxRayParam = Math.min(maxRayParam, tMod.getTauBranch(fromBranch,
                                                                  !isPtoS)
                    .getMaxRayParam());
        } else if(endAction == REFLECTBOT) {
            maxRayParam = Math.min(maxRayParam, tMod.getTauBranch(fromBranch,
                                                                  isPtoS)
                    .getMinTurnRayParam());
            maxRayParam = Math.min(maxRayParam, tMod.getTauBranch(fromBranch,
                                                                  !isPtoS)
                    .getMinTurnRayParam());
        } else if(endAction == TRANSUP) {
            maxRayParam = Math.min(maxRayParam, tMod.getTauBranch(fromBranch,
                                                                  isPtoS)
                    .getMaxRayParam());
            maxRayParam = Math.min(maxRayParam,
                                   tMod.getTauBranch(fromBranch - 1, !isPtoS)
                                           .getMinTurnRayParam());
        } else if(endAction == TRANSDOWN) {
            maxRayParam = Math.min(maxRayParam, tMod.getTauBranch(fromBranch,
                                                                  isPtoS)
                    .getMinRayParam());
            maxRayParam = Math.min(maxRayParam,
                                   tMod.getTauBranch(fromBranch + 1, !isPtoS)
                                           .getMaxRayParam());
        } else {
            throw new TauModelException("Illegal endAction: endAction="
                    + endAction);
        }
    }

    /*
     * Adds the branch numbers from startBranch to endBranch, inclusive, to
     * branchSeq, in order. Also, currBranch is set correctly based on the value
     * of endAction. endAction can be one of TRANSUP, TRANSDOWN, REFLECTTOP,
     * REFLECTBOT, or TURN.
     */
    protected void addToBranch(TauModel tMod,
                               int startBranch,
                               int endBranch,
                               boolean isPWave,
                               int endAction) throws TauModelException {
        int endOffset;
        boolean isDownGoing;
        this.endAction = endAction;
        if(DEBUG) {
            System.out.print("start=" + startBranch + " end=" + endBranch
                    + " endOffset=");
            if(endAction == TURN) {
                System.out.println("TURN");
            } else if(endAction == REFLECTTOP) {
                System.out.println("REFLECTTOP");
            } else if(endAction == REFLECTBOT) {
                System.out.println("REFLECTBOT");
            } else if(endAction == TRANSUP) {
                System.out.println("TRANSUP");
            } else if(endAction == TRANSDOWN) {
                System.out.println("TRANSDOWN");
            } else {
                System.out.println(endAction);
            }
        }
        if(endAction == TURN) {
            endOffset = 0;
            isDownGoing = true;
            minRayParam = Math.max(minRayParam, tMod.getTauBranch(endBranch,
                                                                  isPWave)
                    .getMinTurnRayParam());
        } else if(endAction == REFLECTTOP) {
            endOffset = 0;
            isDownGoing = false;
            maxRayParam = Math.min(maxRayParam, tMod.getTauBranch(endBranch,
                                                                  isPWave)
                    .getMaxRayParam());
        } else if(endAction == REFLECTBOT) {
            endOffset = 0;
            isDownGoing = true;
            maxRayParam = Math.min(maxRayParam, tMod.getTauBranch(endBranch,
                                                                  isPWave)
                    .getMinTurnRayParam());
        } else if(endAction == TRANSUP) {
            endOffset = -1;
            isDownGoing = false;
            maxRayParam = Math.min(maxRayParam, tMod.getTauBranch(endBranch,
                                                                  isPWave)
                    .getMaxRayParam());
        } else if(endAction == TRANSDOWN) {
            endOffset = 1;
            isDownGoing = true;
            maxRayParam = Math.min(maxRayParam, tMod.getTauBranch(endBranch,
                                                                  isPWave)
                    .getMinRayParam());
        } else {
            throw new TauModelException("Illegal endAction: endAction="
                    + endAction);
        }
        if(isDownGoing) {
            /* Must be downgoing, so use i++. */
            for(int i = startBranch; i <= endBranch; i++) {
                branchSeq.add(new Integer(i));
                downGoing.add(new Boolean(isDownGoing));
                waveType.add(new Boolean(isPWave));
                legAction.add(new Integer(endAction));
            }
            if(DEBUG) {
                for(int i = startBranch; i <= endBranch; i++) {
                    System.out.println("i=" + i + " isDownGoing=" + isDownGoing
                            + " isPWave=" + isPWave + " startBranch="
                            + startBranch + " endBranch=" + endBranch + " "
                            + endAction);
                }
            }
        } else {
            /* Must be up going so use i--. */
            for(int i = startBranch; i >= endBranch; i--) {
                branchSeq.add(new Integer(i));
                downGoing.add(new Boolean(isDownGoing));
                waveType.add(new Boolean(isPWave));
                legAction.add(new Integer(endAction));
            }
            if(DEBUG) {
                for(int i = startBranch; i >= endBranch; i--) {
                    System.out.println("i=" + i + " isDownGoing=" + isDownGoing
                            + " isPWave=" + isPWave + " startBranch="
                            + startBranch + " endBranch=" + endBranch + " "
                            + endAction);
                }
            }
        }
        currBranch = endBranch + endOffset;
    }

    /**
     * Finds the closest discontinuity to the given depth that can have
     * refletions and phase transformations.
     * 
     * @return the branch number with the closest top depth.
     */
    public int closestBranchToDepth(TauModel tMod, String depthString) {
        if(depthString.equals("m")) {
            return tMod.getMohoBranch();
        } else if(depthString.equals("c")) {
            return tMod.getCmbBranch();
        } else if(depthString.equals("i")) {
            return tMod.getIocbBranch();
        }
        // nonstandard boundary, given by a number, so we must look for it
        int disconBranch = -1;
        double disconMax = Double.MAX_VALUE;
        double disconDepth = (Double.valueOf(depthString)).doubleValue();
        TauBranch tBranch;
        for(int i = 0; i < tMod.getNumBranches(); i++) {
            tBranch = tMod.getTauBranch(i, PWAVE);
            if(Math.abs(disconDepth - tBranch.getTopDepth()) < disconMax
                    && !tMod.isNoDisconDepth(tBranch.getTopDepth())) {
                disconBranch = i;
                disconMax = Math.abs(disconDepth - tBranch.getTopDepth());
            }
        }
        return disconBranch;
    }

    /**
     * Constructs a branch sequence from the given phase name and tau model.
     */
    protected void parseName(TauModel tMod) throws TauModelException {
        String prevLeg;
        String currLeg = (String)legs.get(0);
        String nextLeg = currLeg;
        branchSeq.clear();
        boolean isPWave = PWAVE;
        boolean isPWavePrev = isPWave;
        /*
         * Deal with surface wave velocities first, since they are a special
         * case.
         */
        if(legs.size() == 2 && currLeg.endsWith("kmps")) {
            return;
        }
        /* Make a check for J legs if the model doesn not allow J */
        if(name.indexOf('J') != -1
                && !tMod.getSlownessModel().isAllowInnerCoreS()) {
            throw new TauModelException("'J' phases were not created for this model: "
                    + name);
        }
        /* set currWave to be the wave type for this leg, 'P' or 'S'. */
        if(currLeg.equals("p") || currLeg.startsWith("P")
                || currLeg.equals("K") || currLeg.equals("k")) {
            isPWave = PWAVE;
            isPWavePrev = isPWave;
        } else if(currLeg.equals("s") || currLeg.startsWith("S")
                || currLeg.equals("J")) {
            isPWave = SWAVE;
            isPWavePrev = isPWave;
        }
        /*
         * First, decide whether the ray is up going or downgoing from the
         * source. If it is up going then the first branch number would be
         * tMod.getSourceBranch()-1 and downgoing would be
         * tMod.getSourceBranch().
         */
        if(currLeg.startsWith("s") || currLeg.startsWith("S")) {
            // Exclude S sources in fluids
            double sdep = tMod.getSourceDepth();
            if(sdep > tMod.getCmbDepth() && sdep < tMod.getIocbDepth()) {
                maxRayParam = minRayParam = -1;
                return;
            }
        }
        if(currLeg.startsWith("P")
                || currLeg.startsWith("S")
                || (expert && (currLeg.startsWith("K") || currLeg.startsWith("I")))) {
            // Downgoing from source
            currBranch = tMod.getSourceBranch();
            endAction = REFLECTBOT; // treat initial downgoing as if it were a
            // underside reflection
        } else if(currLeg.equals("p") || currLeg.equals("s")
                || (expert && currLeg.startsWith("k"))) {
            // Up going from source
            endAction = REFLECTTOP; // treat initial upgoing as if it were a
            // topside reflection
            if(tMod.getSourceBranch() != 0) {
                currBranch = tMod.getSourceBranch() - 1;
            } else {
                /*
                 * p and s for zero source depth are only at zero distance and
                 * then can be called P or S.
                 */
                maxRayParam = -1;
                minRayParam = -1;
                return;
            }
        } else {
            throw new TauModelException("First phase not recognized: "
                    + currLeg
                    + " must be one of P, Pg, Pn, Pdiff, p or the S equivalents");
        }
        /*
         * Set maxRayParam to be a horizontal ray leaving the source and set
         * minRayParam to be a vertical (p=0) ray.
         */
        if(tMod.getSourceBranch() != 0) {
            maxRayParam = Math.max(tMod.getTauBranch(tMod.getSourceBranch() - 1,
                                                     isPWave)
                                           .getMinTurnRayParam(),
                                   tMod.getTauBranch(tMod.getSourceBranch(),
                                                     isPWave).getMaxRayParam());
        } else {
            maxRayParam = tMod.getTauBranch(tMod.getSourceBranch(), isPWave)
                    .getMaxRayParam();
        }
        minRayParam = 0.0;
        int disconBranch = 0;
        double nextLegDepth = 0.0;
        boolean isLegDepth, isNextLegDepth = false;
        endAction = TRANSDOWN;
        /*
         * Now loop over all of the phase legs and construct the proper branch
         * sequence.
         */
        currLeg = "START"; // So the prevLeg isn't wrong on the first pass
        for(int legNum = 0; legNum < legs.size() - 1; legNum++) {
            prevLeg = currLeg;
            currLeg = nextLeg;
            nextLeg = (String)legs.get(legNum + 1);
            if(DEBUG) {
                System.out.println(legNum + "  " + prevLeg + "  " + currLeg
                        + "  " + nextLeg);
            }
            isLegDepth = isNextLegDepth;
            // find out if the next leg represents a phase conversion depth
            try {
                nextLegDepth = (new Double(nextLeg)).doubleValue();
                isNextLegDepth = true;
            } catch(NumberFormatException e) {
                nextLegDepth = -1;
                isNextLegDepth = false;
            }
            /* set currWave to be the wave type for this leg, 'P' or 'S'. */
            isPWavePrev = isPWave;
            if(currLeg.equals("p") || currLeg.startsWith("P")
                    || currLeg.equals("k") || currLeg.equals("I")) {
                isPWave = PWAVE;
            } else if(currLeg.equals("s") || currLeg.startsWith("S")
                    || currLeg.equals("J")) {
                isPWave = SWAVE;
            } else if(currLeg.equals("K")) {
                /*
                 * here we want to use whatever isPWave was on the last leg so
                 * do nothing. This makes sure we us the correct maxRayParam
                 * from the correct TauBranch within the outer core. In other
                 * words K has a high slowness zone if it entered the outer core
                 * as a mantle P wave, but doesn't if it entered as a mantle S
                 * wave. It shouldn't matter for inner core to outer core type
                 * legs.
                 */
            }
            // check to see if there has been a phase conversion
            if(branchSeq.size() > 0 && isPWavePrev != isPWave) {
                phaseConversion(tMod,
                                ((Integer)branchSeq.get(branchSeq.size() - 1)).intValue(),
                                endAction,
                                isPWavePrev);
            }
            /* Deal with p and s case first. */
            if(currLeg.equals("p") || currLeg.equals("s")
                    || currLeg.equals("k")) {
                if(nextLeg.startsWith("v")) {
                    throw new TauModelException("p and s must always be up going "
                            + " and cannot come immediately before a top-side reflection."
                            + " currLeg=" + currLeg + " nextLeg=" + nextLeg);
                } else if(nextLeg.startsWith("^")) {
                    disconBranch = closestBranchToDepth(tMod,
                                                        nextLeg.substring(1));
                    if(currBranch >= disconBranch) {
                        addToBranch(tMod,
                                    currBranch,
                                    disconBranch,
                                    isPWave,
                                    REFLECTTOP);
                    } else {
                        throw new TauModelException("Phase not recognized: "
                                + currLeg + " followed by " + nextLeg
                                + " when currBranch=" + currBranch
                                + " > disconBranch=" + disconBranch);
                    }
                } else if(nextLeg.equals("m")
                        && currBranch >= tMod.getMohoBranch()) {
                    addToBranch(tMod,
                                currBranch,
                                tMod.getMohoBranch(),
                                isPWave,
                                TRANSUP);
                } else if(nextLeg.startsWith("P") || nextLeg.startsWith("S")
                        || nextLeg.equals("K") || nextLeg.equals("END")) {
                    disconBranch = nextLeg.equals("K") ? tMod.getCmbBranch()
                            : 0;
                    addToBranch(tMod,
                                currBranch,
                                disconBranch,
                                isPWave,
                                (currLeg.equals("k") & !nextLeg.equals("K") ? TRANSUP
                                        : REFLECTTOP));
                } else if(isNextLegDepth) {
                    disconBranch = closestBranchToDepth(tMod, nextLeg);
                    addToBranch(tMod,
                                currBranch,
                                disconBranch,
                                isPWave,
                                TRANSUP);
                } else {
                    throw new TauModelException("Phase not recognized: "
                            + currLeg + " followed by " + nextLeg);
                }
                /* Now deal with P and S case. */
            } else if(currLeg.equals("P") || currLeg.equals("S")) {
                if(nextLeg.equals("P") || nextLeg.equals("S")
                        || nextLeg.equals("Pn") || nextLeg.equals("Sn")
                        || nextLeg.equals("END")) {
                    if(endAction == TRANSDOWN || endAction == REFLECTTOP) {
                        // downgoing, so must first turn in mantle
                        addToBranch(tMod,
                                    currBranch,
                                    tMod.getCmbBranch() - 1,
                                    isPWave,
                                    TURN);
                    }
                    addToBranch(tMod, currBranch, 0, isPWave, REFLECTTOP);
                } else if(nextLeg.startsWith("v")) {
                    disconBranch = closestBranchToDepth(tMod,
                                                        nextLeg.substring(1));
                    if(currBranch <= disconBranch - 1) {
                        addToBranch(tMod,
                                    currBranch,
                                    disconBranch - 1,
                                    isPWave,
                                    REFLECTBOT);
                    } else {
                        throw new TauModelException("Phase not recognized: "
                                + currLeg + " followed by " + nextLeg
                                + " when currBranch=" + currBranch
                                + " < disconBranch=" + disconBranch);
                    }
                } else if(nextLeg.startsWith("^")) {
                    disconBranch = closestBranchToDepth(tMod,
                                                        nextLeg.substring(1));
                    if(prevLeg.equals("K")) {
                        addToBranch(tMod,
                                    currBranch,
                                    disconBranch,
                                    isPWave,
                                    REFLECTTOP);
                    } else if(prevLeg.startsWith("^") || prevLeg.equals("P")
                            || prevLeg.equals("S") || prevLeg.equals("p")
                            || prevLeg.equals("s") || prevLeg.equals("START")) {
                        addToBranch(tMod,
                                    currBranch,
                                    tMod.getCmbBranch() - 1,
                                    isPWave,
                                    TURN);
                        addToBranch(tMod,
                                    currBranch,
                                    disconBranch,
                                    isPWave,
                                    REFLECTTOP);
                    } else if((prevLeg.startsWith("v") && disconBranch < closestBranchToDepth(tMod,
                                                                                              prevLeg.substring(1)))
                            || (prevLeg.equals("m") && disconBranch < tMod.getMohoBranch())
                            || (prevLeg.equals("c") && disconBranch < tMod.getCmbBranch())) {
                        addToBranch(tMod,
                                    currBranch,
                                    disconBranch,
                                    isPWave,
                                    REFLECTTOP);
                    } else {
                        throw new TauModelException("Phase not recognized: "
                                + currLeg + " followed by " + nextLeg
                                + " when currBranch=" + currBranch
                                + " > disconBranch=" + disconBranch);
                    }
                } else if(nextLeg.equals("c")) {
                    addToBranch(tMod,
                                currBranch,
                                tMod.getCmbBranch() - 1,
                                isPWave,
                                REFLECTBOT);
                } else if(nextLeg.equals("K")) {
                    addToBranch(tMod,
                                currBranch,
                                tMod.getCmbBranch() - 1,
                                isPWave,
                                TRANSDOWN);
                } else if(nextLeg.equals("m")
                        || (isNextLegDepth && nextLegDepth < tMod.getCmbDepth())) {
                    // treat the moho in the same wasy as 410 type
                    // discontinuities
                    disconBranch = closestBranchToDepth(tMod, nextLeg);
                    if(DEBUG) {
                        System.out.println("DisconBranch=" + disconBranch
                                + " for " + nextLeg);
                        System.out.println(tMod.getTauBranch(disconBranch,
                                                             isPWave)
                                .getTopDepth());
                    }
                    if(endAction == TURN || endAction == REFLECTBOT
                            || endAction == TRANSUP) {
                        // upgoing section
                        if(disconBranch > currBranch) {
                            // check for discontinuity below the current branch
                            // when the ray should be upgoing
                            throw new TauModelException("Phase not recognized: "
                                    + currLeg
                                    + " followed by "
                                    + nextLeg
                                    + " when currBranch="
                                    + currBranch
                                    + " > disconBranch=" + disconBranch);
                        }
                        addToBranch(tMod,
                                    currBranch,
                                    disconBranch,
                                    isPWave,
                                    TRANSUP);
                    } else {
                        // downgoing section, must look at the leg after the
                        // next
                        // leg to determine whether to convert on the downgoing
                        // or
                        // upgoing part of the path
                        String nextNextLeg = (String)legs.get(legNum + 2);
                        if(nextNextLeg.equals("p") || nextNextLeg.equals("s")) {
                            // convert on upgoing section
                            addToBranch(tMod,
                                        currBranch,
                                        tMod.getCmbBranch() - 1,
                                        isPWave,
                                        TURN);
                            addToBranch(tMod,
                                        currBranch,
                                        disconBranch,
                                        isPWave,
                                        TRANSUP);
                        } else if(nextNextLeg.equals("P")
                                || nextNextLeg.equals("S")) {
                            if(disconBranch > currBranch) {
                                // discon is below current loc
                                addToBranch(tMod,
                                            currBranch,
                                            disconBranch - 1,
                                            isPWave,
                                            TRANSDOWN);
                            } else {
                                // discon is above current loc, but we have a
                                // downgoing ray, so this is an illegal ray for
                                // this source depth
                                maxRayParam = -1;
                                if(DEBUG) {
                                    System.out.println("Cannot phase convert on the "
                                            + "downgoing side if the discontinuity is above "
                                            + "the phase leg starting point, "
                                            + currLeg
                                            + " "
                                            + nextLeg
                                            + " "
                                            + nextNextLeg
                                            + ", so this phase, "
                                            + getName()
                                            + " is illegal for this sourceDepth.");
                                }
                                return;
                            }
                        } else {
                            throw new TauModelException("Phase not recognized: "
                                    + currLeg
                                    + " followed by "
                                    + nextLeg
                                    + " followed by " + nextNextLeg);
                        }
                    }
                } else {
                    throw new TauModelException("Phase not recognized: "
                            + currLeg + " followed by " + nextLeg);
                }
            } else if(currLeg.startsWith("P") || currLeg.startsWith("S")) {
                if(currLeg.equals("Pdiff") || currLeg.equals("Sdiff")) {
                    /*
                     * in the diffracted case we trick addToBranch into thinking
                     * we are turning, but then make the maxRayParam equal to
                     * minRayParam, which is the deepest turning ray.
                     */
                    if(maxRayParam >= tMod.getTauBranch(tMod.getCmbBranch() - 1,
                                                        isPWave)
                            .getMinTurnRayParam()
                            && minRayParam <= tMod.getTauBranch(tMod.getCmbBranch() - 1,
                                                                isPWave)
                                    .getMinTurnRayParam()) {
                        addToBranch(tMod,
                                    currBranch,
                                    tMod.getCmbBranch() - 1,
                                    isPWave,
                                    TURN);
                        maxRayParam = minRayParam;
                        if(nextLeg.equals("END") || nextLeg.startsWith("P")
                                || nextLeg.startsWith("S")) {
                            addToBranch(tMod,
                                        currBranch,
                                        0,
                                        isPWave,
                                        REFLECTTOP);
                        }
                    } else {
                        // can't have head wave as ray param is not within range
                        maxRayParam = -1;
                        if(DEBUG) {
                            System.out.println("Cannot have the head wave "
                                    + currLeg + " within phase " + name
                                    + " for this sourceDepth and/or path.");
                        }
                        return;
                    }
                } else if(currLeg.equals("Pg") || currLeg.equals("Sg")
                        || currLeg.equals("Pn") || currLeg.equals("Sn")) {
                    if(currBranch >= tMod.getMohoBranch()) {
                        /*
                         * Pg, Pn, Sg and Sn must be above the moho and so is
                         * not valid for rays coming upwards from below,
                         * possibly due to the source depth. Setting maxRayParam =
                         * -1 effectively disallows this phase.
                         */
                        maxRayParam = -1;
                        if(DEBUG) {
                            System.out.println("(currBranch >= tMod.getMohoBranch() "
                                    + currBranch
                                    + " "
                                    + tMod.getMohoBranch()
                                    + " so there cannot be a "
                                    + currLeg
                                    + " phase for this sourceDepth and/or path.");
                        }
                        return;
                    }
                    if(currLeg.equals("Pg") || currLeg.equals("Sg")) {
                        addToBranch(tMod,
                                    currBranch,
                                    tMod.getMohoBranch() - 1,
                                    isPWave,
                                    TURN);
                        addToBranch(tMod, currBranch, 0, isPWave, REFLECTTOP);
                    } else if(currLeg.equals("Pn") || currLeg.equals("Sn")) {
                        /*
                         * in the refracted case we trick addToBranch into
                         * thinking we are turning below the moho, but then make
                         * the minRayParam equal to maxRayParam, which is the
                         * head wave ray.
                         */
                        if(maxRayParam >= tMod.getTauBranch(tMod.getMohoBranch(),
                                                            isPWave)
                                .getMaxRayParam()
                                && minRayParam <= tMod.getTauBranch(tMod.getMohoBranch(),
                                                                    isPWave)
                                        .getMaxRayParam()) {
                            addToBranch(tMod,
                                        currBranch,
                                        tMod.getMohoBranch(),
                                        isPWave,
                                        TURN);
                            addToBranch(tMod,
                                        currBranch,
                                        tMod.getMohoBranch(),
                                        isPWave,
                                        TRANSUP);
                            minRayParam = maxRayParam;
                            if(nextLeg.equals("END") || nextLeg.startsWith("P")
                                    || nextLeg.startsWith("S")) {
                                addToBranch(tMod,
                                            currBranch,
                                            0,
                                            isPWave,
                                            REFLECTTOP);
                            }
                        } else {
                            // can't have head wave as ray param is not within
                            // range
                            maxRayParam = -1;
                            if(DEBUG) {
                                System.out.println("Cannot have the head wave "
                                        + currLeg + " within phase " + name
                                        + " for this sourceDepth and/or path.");
                            }
                            return;
                        }
                    }
                } else {
                    throw new TauModelException("Phase not recognized: "
                            + currLeg + " followed by " + nextLeg);
                }
            } else if(currLeg.equals("K")) {
                /* Now deal with K. */
                if(nextLeg.equals("P") || nextLeg.equals("S")) {
                    if(prevLeg.equals("P") || prevLeg.equals("S")
                            || prevLeg.equals("K") || prevLeg.equals("k")
                            || prevLeg.equals("START")) {
                        addToBranch(tMod,
                                    currBranch,
                                    tMod.getIocbBranch() - 1,
                                    isPWave,
                                    TURN);
                    }
                    addToBranch(tMod,
                                currBranch,
                                tMod.getCmbBranch(),
                                isPWave,
                                TRANSUP);
                } else if(nextLeg.equals("K")) {
                    if(prevLeg.equals("P") || prevLeg.equals("S")
                            || prevLeg.equals("K")) {
                        addToBranch(tMod,
                                    currBranch,
                                    tMod.getIocbBranch() - 1,
                                    isPWave,
                                    TURN);
                    }
                    addToBranch(tMod,
                                currBranch,
                                tMod.getCmbBranch(),
                                isPWave,
                                REFLECTTOP);
                } else if(nextLeg.equals("I") || nextLeg.equals("J")) {
                    addToBranch(tMod,
                                currBranch,
                                tMod.getIocbBranch() - 1,
                                isPWave,
                                TRANSDOWN);
                } else if(nextLeg.equals("i")) {
                    addToBranch(tMod,
                                currBranch,
                                tMod.getIocbBranch() - 1,
                                isPWave,
                                REFLECTBOT);
                } else {
                    throw new TauModelException("Phase not recognized: "
                            + currLeg + " followed by " + nextLeg);
                }
            } else if(currLeg.equals("I") || currLeg.equals("J")) {
                /* And now consider inner core, I and J. */
                addToBranch(tMod,
                            currBranch,
                            tMod.getNumBranches() - 1,
                            isPWave,
                            TURN);
                if(nextLeg.equals("I") || nextLeg.equals("J")) {
                    addToBranch(tMod,
                                currBranch,
                                tMod.getIocbBranch(),
                                isPWave,
                                REFLECTTOP);
                } else if(nextLeg.equals("K")) {
                    addToBranch(tMod,
                                currBranch,
                                tMod.getIocbBranch(),
                                isPWave,
                                TRANSUP);
                }
            } else if(currLeg.equals("m")) {} else if(currLeg.equals("c")) {} else if(currLeg.equals("i")) {} else if(currLeg.startsWith("^")) {} else if(currLeg.startsWith("v")) {} else if(isLegDepth) {} else {
                throw new TauModelException("Phase not recognized: " + currLeg
                        + " followed by " + nextLeg);
            }
        }
    }

    /**
     * Tokenizes a phase name into legs, ie PcS becomes 'P'+'c'+'S' while p^410P
     * would become 'p'+'^410'+'P'. Once a phase name has been broken into
     * tokens we can begin to construct the sequence of branches to which it
     * corresponds. Only minor error checking is done at this point, for
     * instance PIP generates an exception but ^410 doesn't. It also appends
     * "END" as the last leg.
     * 
     * @throws TauModelException
     *             if the phase name cannot be tokenized.
     */
    protected static ArrayList<String> legPuller(String name) throws TauModelException {
        int offset = 0;
        ArrayList<String> legs = new ArrayList<String>();
        /* Special case for surface wave velocity. */
        if(name.endsWith("kmps")) {
            try {
                legs.add(name);
            } catch(NumberFormatException e) {
                throw new TauModelException("Invalid phase name:\n" + name);
            }
        } else
            while(offset < name.length()) {
                if(name.charAt(offset) == 'K' || name.charAt(offset) == 'I'
                        || name.charAt(offset) == 'k'
                        || name.charAt(offset) == 'J'
                        || name.charAt(offset) == 'p'
                        || name.charAt(offset) == 's'
                        || name.charAt(offset) == 'm'
                        || name.charAt(offset) == 'c'
                        || name.charAt(offset) == 'i') {
                    // Do the easy ones, ie K,k,I,J,p,s,m,c,i
                    legs.add(name.substring(offset, offset + 1));
                    offset = offset + 1;
                } else if(name.charAt(offset) == 'P'
                        || name.charAt(offset) == 'S') {
                    /*
                     * Now it gets complicated, first see if the next char is
                     * part of a different leg or we are at the end.
                     */
                    if(offset + 1 == name.length()
                            || name.charAt(offset + 1) == 'P'
                            || name.charAt(offset + 1) == 'S'
                            || name.charAt(offset + 1) == 'K'
                            || name.charAt(offset + 1) == 'm'
                            || name.charAt(offset + 1) == 'c'
                            || name.charAt(offset + 1) == '^'
                            || name.charAt(offset + 1) == 'v'
                            || Character.isDigit(name.charAt(offset + 1))) {
                        legs.add(name.substring(offset, offset + 1));
                        offset++;
                    } else if(name.charAt(offset + 1) == 'p'
                            || name.charAt(offset + 1) == 's') {
                        throw new TauModelException("Invalid phase name:\n"
                                + name.charAt(offset)
                                + " cannot be followed by "
                                + name.charAt(offset + 1) + " in " + name);
                    } else if(name.charAt(offset + 1) == 'g'
                            || name.charAt(offset + 1) == 'b'
                            || name.charAt(offset + 1) == 'n') {
                        /* The leg is not described by one letter, check for 2. */
                        legs.add(name.substring(offset, offset + 2));
                        offset = offset + 2;
                    } else if(name.length() >= offset + 5
                            && (name.substring(offset, offset + 5)
                                    .equals("Sdiff") || name.substring(offset,
                                                                       offset + 5)
                                    .equals("Pdiff"))) {
                        legs.add(name.substring(offset, offset + 5));
                        offset = offset + 5;
                    } else {
                        throw new TauModelException("Invalid phase name:\n"
                                + name.substring(offset) + " in " + name);
                    }
                } else if(name.charAt(offset) == '^'
                        || name.charAt(offset) == 'v') {
                    /*
                     * Top side or bottom side reflections, check for standard
                     * boundaries and then check for numerical ones.
                     */
                    if(name.charAt(offset + 1) == 'm'
                            || name.charAt(offset + 1) == 'c'
                            || name.charAt(offset + 1) == 'i') {
                        legs.add(name.substring(offset, offset + 2));
                        offset = offset + 2;
                    } else if(Character.isDigit(name.charAt(offset + 1))
                            || name.charAt(offset + 1) == '.') {
                        String numString = name.substring(offset, offset + 1);
                        offset++;
                        while(Character.isDigit(name.charAt(offset))
                                || name.charAt(offset) == '.') {
                            numString += name.substring(offset, offset + 1);
                            offset++;
                        }
                        try {
                            legs.add(numString);
                        } catch(NumberFormatException e) {
                            throw new TauModelException("Invalid phase name: "
                                    + numString + "\n" + e.getMessage()
                                    + " in " + name);
                        }
                    } else {
                        throw new TauModelException("Invalid phase name:\n"
                                + name.substring(offset) + " in " + name);
                    }
                } else if(Character.isDigit(name.charAt(offset))
                        || name.charAt(offset) == '.') {
                    String numString = name.substring(offset, offset + 1);
                    offset++;
                    while(Character.isDigit(name.charAt(offset))
                            || name.charAt(offset) == '.') {
                        numString += name.substring(offset, offset + 1);
                        offset++;
                    }
                    try {
                        legs.add(numString);
                    } catch(NumberFormatException e) {
                        throw new TauModelException("Invalid phase name: "
                                + numString + "\n" + e.getMessage() + " in "
                                + name);
                    }
                } else {
                    throw new TauModelException("Invalid phase name:\n"
                            + name.substring(offset) + " in " + name);
                }
            }
        legs.add(new String("END"));
        String validationMsg = phaseValidate(legs);
        if(validationMsg != null) {
            throw new TauModelException("Phase failed validation: " + name
                    + "  " + validationMsg);
        }
        return legs;
    }

    protected void createPuristName(TauModel tMod) {
        String currLeg = (String)legs.get(0);
        /*
         * Deal with surface wave velocities first, since they are a special
         * case.
         */
        if(legs.size() == 2 && currLeg.endsWith("kmps")) {
            puristName = name;
            return;
        }
        puristName = "";
        double legDepth;
        int intLegDepth;
        int disconBranch;
        // only loop to size()-1 as last leg is always "END"
        for(int legNum = 0; legNum < legs.size() - 1; legNum++) {
            currLeg = (String)legs.get(legNum);
            // find out if the next leg represents a
            // phase conversion or reflection depth
            if(currLeg.startsWith("v") || currLeg.startsWith("^")) {
                disconBranch = closestBranchToDepth(tMod, currLeg.substring(1));
                legDepth = tMod.getTauBranch(disconBranch, true).getTopDepth();
                puristName += currLeg.substring(0, 1);
                if(legDepth == Math.rint(legDepth)) {
                    intLegDepth = (int)legDepth;
                    puristName += intLegDepth;
                } else {
                    puristName += legDepth;
                }
            } else {
                try {
                    legDepth = (new Double(currLeg)).doubleValue();
                    // only get this far if the currLeg is a number,
                    // otherwise exception
                    disconBranch = closestBranchToDepth(tMod, currLeg);
                    legDepth = tMod.getTauBranch(disconBranch, true)
                            .getTopDepth();
                    if(legDepth == Math.rint(legDepth)) {
                        intLegDepth = (int)legDepth;
                        puristName += intLegDepth;
                    } else {
                        puristName += legDepth;
                    }
                } catch(NumberFormatException e) {
                    puristName += currLeg;
                }
            }
        }
    }

    /**
     * Sums the appropriate branches for this phase.
     * 
     * @throws TauModelException
     *             if the topDepth of the high slowness zone is not contained
     *             within the TauModel. This should never happen and would
     *             indicate an invalid TauModel.
     */
    protected void sumBranches(TauModel tMod) throws TauModelException {
        if(maxRayParam < 0.0 || minRayParam > maxRayParam) {
            /* Phase has no arrivals, possibly due to source depth. */
            rayParams = new double[0];
            minRayParam = -1;
            maxRayParam = -1;
            dist = new double[0];
            time = new double[0];
            maxDistance = -1;
            return;
        }
        /* Special case for surface waves. */
        if(name.endsWith("kmps")) {
            dist = new double[2];
            time = new double[2];
            rayParams = new double[2];
            dist[0] = 0.0;
            time[0] = 0.0;
            rayParams[0] = tMod.radiusOfEarth
                    / Double.valueOf(name.substring(0, name.length() - 4))
                            .doubleValue();
            dist[1] = 2 * Math.PI;
            time[1] = 2
                    * Math.PI
                    * tMod.radiusOfEarth
                    / Double.valueOf(name.substring(0, name.length() - 4))
                            .doubleValue();
            rayParams[1] = rayParams[0];
            minDistance = 0.0;
            maxDistance = 2 * Math.PI;
            downGoing.add(true);
            return;
        }
        /*
         * Find the ray parameter index that corresponds to the minRayParam and
         * maxRayParam.
         */
        for(int i = 0; i < tMod.rayParams.length; i++) {
            if(tMod.rayParams[i] >= minRayParam) {
                minRayParamIndex = i;
            }
            if(tMod.rayParams[i] >= maxRayParam) {
                maxRayParamIndex = i;
            }
        }
        if(maxRayParamIndex == 0
                && minRayParamIndex == tMod.rayParams.length - 1) {
            // all ray parameters are valid so just copy
            rayParams = new double[tMod.rayParams.length];
            System.arraycopy(tMod.rayParams,
                             0,
                             rayParams,
                             0,
                             tMod.rayParams.length);
        } else if(maxRayParamIndex == minRayParamIndex) {
            if(name.indexOf("Sdiff") != -1 || name.indexOf("Pdiff") != -1) {
                rayParams = new double[2];
                rayParams[0] = minRayParam;
                rayParams[1] = minRayParam;
            } else if(name.indexOf("Pn") != -1 || name.indexOf("Sn") != -1) {
                rayParams = new double[2];
                rayParams[0] = minRayParam;
                rayParams[1] = minRayParam;
            } else if(name.endsWith("kmps")) {
                rayParams = new double[2];
                rayParams[0] = 0;
                rayParams[1] = maxRayParam;
            } else {
                rayParams = new double[2];
                rayParams[0] = minRayParam;
                rayParams[1] = minRayParam;
            }
        } else {
            if(DEBUG) {
                System.out.println("maxRayParamIndex=" + maxRayParamIndex
                        + " minRayParamIndex=" + minRayParamIndex
                        + " tMod.rayParams.length=" + tMod.rayParams.length
                        + " tMod.rayParams[0]=" + tMod.rayParams[0]
                        + " maxRayParam=" + maxRayParam);
            }
            // only a subset of ray parameters are valid so only use those
            rayParams = new double[minRayParamIndex - maxRayParamIndex + 1];
            System.arraycopy(tMod.rayParams,
                             maxRayParamIndex,
                             rayParams,
                             0,
                             minRayParamIndex - maxRayParamIndex + 1);
        }
        dist = new double[rayParams.length];
        time = new double[rayParams.length];
        /* initialize the counter for each branch to 0. 0 is P and 1 is S. */
        int[][] timesBranches = new int[2][tMod.getNumBranches()];
        for(int i = 0; i < timesBranches[0].length; i++) {
            timesBranches[0][i] = 0;
            timesBranches[1][i] = 0;
        }
        /* Count how many times each branch appears in the path. */
        for(int i = 0; i < branchSeq.size(); i++) {
            if(((Boolean)waveType.get(i)).booleanValue()) {
                timesBranches[0][((Integer)branchSeq.get(i)).intValue()]++;
            } else {
                timesBranches[1][((Integer)branchSeq.get(i)).intValue()]++;
            }
        }
        /* Sum the branches with the appropriate multiplier. */
        for(int j = 0; j < tMod.getNumBranches(); j++) {
            if(timesBranches[0][j] != 0) {
                for(int i = maxRayParamIndex; i < minRayParamIndex + 1; i++) {
                    dist[i - maxRayParamIndex] += timesBranches[0][j]
                            * tMod.getTauBranch(j, PWAVE).getDist(i);
                    time[i - maxRayParamIndex] += timesBranches[0][j]
                            * tMod.getTauBranch(j, PWAVE).time[i];
                }
            }
            if(timesBranches[1][j] != 0) {
                for(int i = maxRayParamIndex; i < minRayParamIndex + 1; i++) {
                    dist[i - maxRayParamIndex] += timesBranches[1][j]
                            * tMod.getTauBranch(j, SWAVE).getDist(i);
                    time[i - maxRayParamIndex] += timesBranches[1][j]
                            * tMod.getTauBranch(j, SWAVE).time[i];
                }
            }
        }
        if(name.indexOf("Sdiff") != -1 || name.indexOf("Pdiff") != -1) {
            if(tMod.getSlownessModel()
                    .depthInHighSlowness(tMod.cmbDepth - 1e-10,
                                         minRayParam,
                                         (name.charAt(0) == 'P'))) {
                /*
                 * No diffraction if there is a high slowness zone at the CMB.
                 */
                minRayParam = -1;
                maxRayParam = -1;
                maxDistance = -1;
                dist = new double[0];
                time = new double[0];
                rayParams = new double[0];
                return;
            } else {
                dist[1] = dist[0] + getMaxDiffraction() * Math.PI / 180.0;
                time[1] = time[0] + getMaxDiffraction() * Math.PI / 180.0
                        * minRayParam;
            }
        } else if(name.indexOf("Pn") != -1 || name.indexOf("Sn") != -1) {
            dist[1] = dist[0] + maxRefraction * Math.PI / 180.0;
            time[1] = time[0] + maxRefraction * Math.PI / 180.0 * minRayParam;
        } else if(maxRayParamIndex == minRayParamIndex) {
            dist[1] = dist[0];
            time[1] = time[0];
        }
        minDistance = Double.MAX_VALUE;
        maxDistance = 0.0;
        for(int j = 0; j < dist.length; j++) {
            if(dist[j] < minDistance) {
                minDistance = dist[j];
            }
            if(dist[j] > maxDistance) {
                maxDistance = dist[j];
            }
        }
        /*
         * Now check to see if our ray parameter range includes any ray
         * parameters that are associated with high slowness zones. If so, then
         * we will need to insert a "shadow zone" into our time and distance
         * arrays. It is represented by a repeated ray parameter.
         */
        DepthRange[] hsz;
        int hSZIndex;
        int indexOffset;
        boolean foundOverlap = false;
        boolean isPWave;
        int branchNum;
        int dummy;
        for(dummy = 0, isPWave = true; dummy < 2; dummy++, isPWave = false) {
            hsz = tMod.getSlownessModel().getHighSlowness(isPWave);
            hSZIndex = 0;
            indexOffset = 0;
            for(int i = 0; i < hsz.length; i++) {
                if(maxRayParam > hsz[i].rayParam
                        && hsz[i].rayParam > minRayParam) {
                    /*
                     * There is a high slowness zone within our ray parameter
                     * range so we might need to add a shadow zone. We need to
                     * check to see if this wave type, P or S, is part of the
                     * phase at this depth/ray parameter.
                     */
                    branchNum = tMod.findBranch(hsz[i].topDepth);
                    foundOverlap = false;
                    for(int legNum = 0; legNum < branchSeq.size(); legNum++) {
                        // check for downgoing legs that cross the high slowness
                        // zone
                        // with the same wave type
                        if(((Integer)branchSeq.get(legNum)).intValue() == branchNum
                                && ((Boolean)waveType.get(legNum)).booleanValue() == isPWave
                                && ((Boolean)downGoing.get(legNum)).booleanValue() == true
                                && ((Integer)branchSeq.get(legNum - 1)).intValue() == branchNum - 1
                                && ((Boolean)waveType.get(legNum - 1)).booleanValue() == isPWave
                                && ((Boolean)downGoing.get(legNum - 1)).booleanValue() == true) {
                            foundOverlap = true;
                            break;
                        }
                    }
                    if(foundOverlap) {
                        double[] newdist = new double[dist.length + 1];
                        double[] newtime = new double[time.length + 1];
                        double[] newrayParams = new double[rayParams.length + 1];
                        for(int j = 0; j < rayParams.length; j++) {
                            if(rayParams[j] == hsz[i].rayParam) {
                                hSZIndex = j;
                                break;
                            }
                        }
                        System.arraycopy(dist, 0, newdist, 0, hSZIndex);
                        System.arraycopy(time, 0, newtime, 0, hSZIndex);
                        System.arraycopy(rayParams,
                                         0,
                                         newrayParams,
                                         0,
                                         hSZIndex);
                        newrayParams[hSZIndex] = hsz[i].rayParam;
                        /* Sum the branches with the appropriate multiplier. */
                        newdist[hSZIndex] = 0.0;
                        newtime[hSZIndex] = 0.0;
                        for(int j = 0; j < tMod.getNumBranches(); j++) {
                            if(timesBranches[0][j] != 0
                                    && tMod.getTauBranch(j, PWAVE)
                                            .getTopDepth() < hsz[i].topDepth) {
                                newdist[hSZIndex] += timesBranches[0][j]
                                        * tMod.getTauBranch(j, PWAVE).dist[maxRayParamIndex
                                                + hSZIndex - indexOffset];
                                newtime[hSZIndex] += timesBranches[0][j]
                                        * tMod.getTauBranch(j, PWAVE).time[maxRayParamIndex
                                                + hSZIndex - indexOffset];
                            }
                            if(timesBranches[1][j] != 0
                                    && tMod.getTauBranch(j, SWAVE)
                                            .getTopDepth() < hsz[i].topDepth) {
                                newdist[hSZIndex] += timesBranches[1][j]
                                        * tMod.getTauBranch(j, SWAVE).dist[maxRayParamIndex
                                                + hSZIndex - indexOffset];
                                newtime[hSZIndex] += timesBranches[1][j]
                                        * tMod.getTauBranch(j, SWAVE).time[maxRayParamIndex
                                                + hSZIndex - indexOffset];
                            }
                        }
                        System.arraycopy(dist,
                                         hSZIndex,
                                         newdist,
                                         hSZIndex + 1,
                                         dist.length - hSZIndex);
                        System.arraycopy(time,
                                         hSZIndex,
                                         newtime,
                                         hSZIndex + 1,
                                         time.length - hSZIndex);
                        System.arraycopy(rayParams,
                                         hSZIndex,
                                         newrayParams,
                                         hSZIndex + 1,
                                         rayParams.length - hSZIndex);
                        indexOffset++;
                        dist = newdist;
                        time = newtime;
                        rayParams = newrayParams;
                    }
                }
            }
        }
    }

    /**
     * Calculates the "pierce points" for the arrivals stored in arrivals. The
     * pierce points are stored within each arrival object.
     */
    public List<Arrival> calcPierce(double deg) throws TauModelException {
        List<Arrival> arrivals = calcTime(deg);
        for (Arrival a : arrivals) {
            calcPierce(a);
        }
        return arrivals;
    }
    
    /** Calculates the pierce points for a particular arrival. The returned arrival is the same
     * as the input arguement but now has the pierce points filled in.
     * @param currArrival
     * @return same arrival with pierce points
     */
    public Arrival calcPierce(Arrival currArrival) {
        double branchDist = 0.0;
        double branchTime = 0.0;
        double prevBranchTime = 0.0;
        List<TimeDist> pierce = new ArrayList<TimeDist>();
        /*
         * Find the ray parameter index that corresponds to the arrival ray
         * parameter in the TauModel, ie it is between rayNum and rayNum+1,
         * We know that it must be <tMod.rayParams.length-1 since the last
         * ray parameter sample is 0, at least in a spherical model...
         */
        int rayNum = 0;
        for(int i = 0; i < tMod.rayParams.length - 1; i++) {
            if(tMod.rayParams[i] >= currArrival.getRayParam()) {
                rayNum = i;
            } else {
                break;
            }
        }
        // here we use ray parameter and dist info stored within the
        // SeismicPhase so we can use currArrival.rayParamIndex, which
        // may not correspond to rayNum (for tMod.rayParams).
        double rayParamA = rayParams[currArrival.getRayParamIndex()];
        double rayParamB = rayParams[currArrival.getRayParamIndex() + 1];
        double distA = dist[currArrival.getRayParamIndex()];
        double distB = dist[currArrival.getRayParamIndex() + 1];
        double distRatio = (currArrival.getDist() - distA) / (distB - distA);
        double distRayParam = distRatio * (rayParamB - rayParamA) + rayParamA;
        /* First pierce point is always 0 distance at the source depth. */
        pierce.add(new TimeDist(distRayParam,
                                             0.0,
                                             0.0,
                                             tMod.getSourceDepth()));
        /*
         * Loop from 0 but already done 0, so the pierce point when the ray
         * leaves branch i is stored in i+1. Use linear interpolation
         * between rays that we know.
         */
        for(int i = 0; i < branchSeq.size(); i++) {
            int branchNum = ((Integer)branchSeq.get(i)).intValue();
            boolean isPWave = ((Boolean)waveType.get(i)).booleanValue();
            if(DEBUG) {
                System.out.println(i + " branchNum =" + branchNum
                        + " downGoing=" + (Boolean)downGoing.get(i)
                        + "  isPWave=" + isPWave);
            }
            /*
             * Save the turning depths for the ray parameter for both P and
             * S waves. This way we get the depth correct for any rays that
             * turn within a layer. We have to do this on a per branch basis
             * because of converted phases, e.g. SKS.
             */
            double turnDepth;
            try {
                if(distRayParam > tMod.getTauBranch(branchNum, isPWave)
                        .getMaxRayParam()) {
                    turnDepth = tMod.getTauBranch(branchNum, isPWave)
                            .getTopDepth();
                } else if(distRayParam <= tMod.getTauBranch(branchNum,
                                                            isPWave)
                        .getMinRayParam()) {
                    turnDepth = tMod.getTauBranch(branchNum, isPWave)
                            .getBotDepth();
                } else {
                    if(isPWave
                            || tMod.getSlownessModel()
                                    .depthInFluid((tMod.getTauBranch(branchNum,
                                                                     isPWave)
                                            .getTopDepth() + tMod.getTauBranch(branchNum,
                                                                               isPWave)
                                            .getBotDepth()) / 2.0)) {
                        turnDepth = tMod.getSlownessModel()
                                .findDepth(distRayParam,
                                           tMod.getTauBranch(branchNum,
                                                             isPWave)
                                                   .getTopDepth(),
                                           tMod.getTauBranch(branchNum,
                                                             isPWave)
                                                   .getBotDepth(),
                                           PWAVE);
                    } else {
                        turnDepth = tMod.getSlownessModel()
                                .findDepth(distRayParam,
                                           tMod.getTauBranch(branchNum,
                                                             isPWave)
                                                   .getTopDepth(),
                                           tMod.getTauBranch(branchNum,
                                                             isPWave)
                                                   .getBotDepth(),
                                           isPWave);
                    }
                }
            } catch(SlownessModelException e) {
                // shouldn't happen but...
                throw new RuntimeException("SeismicPhase.calcPierce: Caught SlownessModelException. "
                        , e);
            }
            double timeA, timeB;
            if(name.indexOf("Pdiff") != -1 || name.indexOf("Pn") != -1
                    || name.indexOf("Sdiff") != -1
                    || name.indexOf("Sn") != -1) {
                /* head waves and diffracted waves are a special case. */
                distA = tMod.getTauBranch(branchNum, isPWave)
                        .getDist(rayNum);
                timeA = tMod.getTauBranch(branchNum, isPWave).time[rayNum];
                distB = tMod.getTauBranch(branchNum, isPWave)
                        .getDist(rayNum);
                timeB = tMod.getTauBranch(branchNum, isPWave).time[rayNum];
            } else {
                distA = tMod.getTauBranch(branchNum, isPWave)
                        .getDist(rayNum);
                timeA = tMod.getTauBranch(branchNum, isPWave).time[rayNum];
                distB = tMod.getTauBranch(branchNum, isPWave)
                        .getDist(rayNum + 1);
                timeB = tMod.getTauBranch(branchNum, isPWave).time[rayNum + 1];
            }
            branchDist += distRatio * (distB - distA) + distA;
            prevBranchTime = branchTime;
            branchTime += distRatio * (timeB - timeA) + timeA;
            double branchDepth;
            if(((Boolean)downGoing.get(i)).booleanValue()) {
                branchDepth = Math.min(tMod.getTauBranch(branchNum, isPWave)
                                               .getBotDepth(),
                                       turnDepth);
            } else {
                branchDepth = Math.min(tMod.getTauBranch(branchNum, isPWave)
                                               .getTopDepth(),
                                       turnDepth);
            }
            // make sure ray actually propagates in this branch, leave
            // a little room for numerical "chatter"
            if(Math.abs(prevBranchTime - branchTime) > 1e-10) {
                pierce.add(new TimeDist(distRayParam,
                                                              branchTime,
                                                              branchDist,
                                                              branchDepth));
                if(DEBUG) {
                    System.out.println(" branchTime=" + branchTime
                            + " branchDist=" + branchDist + " branchDepth="
                            + branchDepth);
                    System.out.println("incrementTime = "
                            + (distRatio * (timeB - timeA)) + " timeB="
                            + timeB + " timeA=" + timeA);
                }
            }
        }
        if(name.indexOf("Pdiff") != -1 || name.indexOf("Pn") != -1
                || name.indexOf("Sdiff") != -1 || name.indexOf("Sn") != -1) {
            pierce = handleHeadOrDiffractedWave(currArrival, pierce);
        } else if(name.indexOf("kmps") != -1) {
            pierce.add(new TimeDist(distRayParam,
                                                 currArrival.getTime(),
                                                 currArrival.getDist(),
                                                 0));
        }
        currArrival.pierce = pierce.toArray(new TimeDist[0]);
        return currArrival;
    }

    /**
     * Here we worry about the special case for head and diffracted
     * waves. It is assumed that a phase can be a diffracted wave or a
     * head wave, but not both. Nor can it be a head wave or diffracted
     * wave for both P and S.
     */
    List<TimeDist> handleHeadOrDiffractedWave(Arrival currArrival, List<TimeDist> orig) {
        String[] phaseSegments = new String[] {"Pn", "Sn", "Pdiff", "Sdiff"};
        String phaseSeg = "";
        for (int i = 0; i < phaseSegments.length; i++) {
            if (name.indexOf(phaseSegments[i]) != -1) {
                phaseSeg = phaseSegments[i];
                break;
            }
        }
        if (phaseSeg.equals("")) {throw new RuntimeException("no head/diff segment in "+name); }
        double headDepth;
        if (phaseSeg.equals("Pn") || phaseSeg.equals("Sn")) {
            headDepth = tMod.getMohoDepth();
        } else {
            headDepth = tMod.getCmbDepth();
        }
        int i = 0;
        int numFound = 0;
        int indexInString = -1;
        // can't have both Pxxx and Sxxx in a head wave phase, so one of these
        // should do nothing
        while((indexInString = name.indexOf(phaseSeg, indexInString + 1)) != -1) {
            numFound++;
        }
        double refractDist = currArrival.getDist() - dist[0];
        double refractTime = refractDist*currArrival.getRayParam();
        List<TimeDist> out = new ArrayList<TimeDist>();
        int j = 0;
        for (TimeDist td : orig) {
            // this is a little weird as we are not checking where we are in the phase name, but simply
            // if the depth matches. This likely works in most cases, but may not for head/diffracted
            // waves that undergo a phase change, if that type of phase can even exist
            out.add(new TimeDist(td.p, td.time+j * refractTime / numFound, td.distRadian + j * refractDist / numFound, td.depth));
            if (td.depth == headDepth) {
                j++;
                out.add(new TimeDist(td.p, td.time+j * refractTime / numFound, td.distRadian + j * refractDist / numFound, td.depth));
            }
        }
        return out;
    }

    /** calculates the paths this phase takes through the earth model. */
    public List<Arrival> calcPath(double deg) {
        
        List<Arrival> arrivals = calcTime(deg);
        for (Arrival a : arrivals) {
            calcPath(a);
        }
        return arrivals;
    }
    
    public Arrival calcPath(Arrival currArrival) {
        ArrayList<TimeDist[]> pathList = new ArrayList<TimeDist[]>();
        /*
         * Find the ray parameter index that corresponds to the arrival ray
         * parameter in the TauModel, ie it is between rayNum and rayNum+1.
         */
        TimeDist[] tempTimeDist = new TimeDist[1];
        tempTimeDist[0] = new TimeDist(currArrival.getRayParam(),
                                       0.0,
                                       0.0,
                                       tMod.getSourceDepth());
        pathList.add(tempTimeDist);
        TimeDist prevTimeDist = tempTimeDist[0];
            for(int i = 0; i < branchSeq.size(); i++) {
                int branchNum = ((Integer)branchSeq.get(i)).intValue();
                boolean isPWave = ((Boolean)waveType.get(i)).booleanValue();
                if(DEBUG) {
                    System.out.println("i=" + i + " branchNum=" + branchNum
                            + " isPWave=" + isPWave + " downgoing="
                            + ((Boolean)downGoing.get(i)).booleanValue());
                }
                try {
                    tempTimeDist = tMod.getTauBranch(branchNum, isPWave)
                            .path(currArrival.getRayParam(),
                                  ((Boolean)downGoing.get(i)).booleanValue(),
                                  tMod.getSlownessModel());
                } catch(SlownessModelException e) {
                    // shouldn't happen but...
                    throw new RuntimeException("SeismicPhase.calcPath: Caught SlownessModelException. "
                            , e);
                }
                if(tempTimeDist != null) {
                    pathList.add(tempTimeDist);
                    for (int j = 0; j < tempTimeDist.length; j++) {
                        if (tempTimeDist[j].getDistDeg() < 0) {
                            throw new RuntimeException("Path is backtracking, no possible: "+j+" ("+tempTimeDist[j]+")");
                        }
                        prevTimeDist = tempTimeDist[j];
                    }
                }
                /*
                 * Here we worry about the special case for head and
                 * diffracted waves.
                 */
                if(branchNum == tMod.cmbBranch - 1 && i < branchSeq.size()-1
                        && ((Integer)branchSeq.get(i + 1)).intValue() == tMod.cmbBranch - 1
                        && (name.indexOf("Pdiff") != -1 || name.indexOf("Sdiff") != -1)) {
                    TimeDist[] diffTD = new TimeDist[1];
                    diffTD[0] = new TimeDist(currArrival.getRayParam(),
                                             (currArrival.getDist() - dist[0])
                                                     * currArrival.getRayParam(),
                                             currArrival.getDist()
                                                     - dist[0],
                                             tMod.cmbDepth);
                    pathList.add(diffTD);
                } else if(branchNum == tMod.mohoBranch - 1 && i < branchSeq.size()-1
                        && ((Integer)branchSeq.get(i + 1)).intValue() == tMod.mohoBranch - 1
                        && (name.indexOf("Pn") != -1 || name.indexOf("Sn") != -1)) {
                    int numFound = 0;
                    int indexInString = -1;
                    // can't have both Pn and Sn in a phase, so one of these
                    // should do nothing
                    while((indexInString = name.indexOf("Pn",
                                                        indexInString + 1)) != -1) {
                        numFound++;
                    }
                    while((indexInString = name.indexOf("Sn",
                                                        indexInString + 1)) != -1) {
                        numFound++;
                    }
                    TimeDist[] headTD = new TimeDist[1];
                    headTD[0] = new TimeDist(currArrival.getRayParam(),
                                             (currArrival.getDist() - dist[0])
                                                     / numFound
                                                     * currArrival.getRayParam(),
                                             (currArrival.getDist() - dist[0])
                                                     / numFound,
                                             tMod.mohoDepth);
                    pathList.add(headTD);
                } 
            }
            if (name.indexOf("kmps") != -1) {
                // kmps phases have no branches, so need to end them at the arrival distance
                TimeDist[] headTD = new TimeDist[1];
                headTD[0] = new TimeDist(currArrival.getRayParam(),
                                         currArrival.getDist()
                                         * currArrival.getRayParam(),
                                         currArrival.getDist(),
                                         0);
                pathList.add(headTD);
            }
            int arraySize = 0;
            for(int i = 0; i < pathList.size(); i++) {
                arraySize += ((TimeDist[])pathList.get(i)).length;
            }
            currArrival.path = new TimeDist[arraySize];
            TimeDist cummulative = new TimeDist(currArrival.getRayParam(),
                                                0.0,
                                                0.0,
                                                currArrival.getSourceDepth());
            TimeDist[] branchPath;
            int numAdded = 0;
            for(int i = 0; i < pathList.size(); i++) {
                branchPath = (TimeDist[])pathList.get(i);
                for(int j = 0; j < branchPath.length; j++) {
                    cummulative.add(branchPath[j]);
                    cummulative.depth = branchPath[j].depth;
                    currArrival.path[numAdded] = (TimeDist)cummulative.clone();
                    if (numAdded > 0 && currArrival.path[numAdded].getDistRadian() < currArrival.path[numAdded-1].getDistRadian()) {
                        throw new RuntimeException("Backtracking ray, not possible: "+numAdded+" "+currArrival.path[numAdded-1]+") > ("+currArrival.path[numAdded]+")");
                    }
                    numAdded++;
                }
            }
        return currArrival;
    }

    private String validationFailMessage = "";

    public String getValidationFailMessage() {
        return validationFailMessage;
    }

    /**
     * Performs consistency checks on the previously tokenized phase name stored
     * in legs. Returns null if all is ok, a message if there is a problem.
     */
    public static String phaseValidate(ArrayList<String> legs) {
        String currToken = (String)legs.get(0);
        String prevToken;
        boolean prevIsReflect = false;
        /* Special cases for diffracted waves. */
        if(legs.size() == 2
                && (currToken.equals("Pdiff") || currToken.equals("Sdiff") || currToken.endsWith("kmps"))
                && ((String)legs.get(1)).equals("END")) {
            return null;
        }
        /* Check first leg. */
        if(!(currToken.equals("Pg") || currToken.equals("Pb")
                || currToken.equals("Pn") || currToken.equals("Pdiff")
                || currToken.equals("Sg") || currToken.equals("Sb")
                || currToken.equals("Sn") || currToken.equals("Sdiff")
                || currToken.equals("P") || currToken.equals("S")
                || currToken.equals("p") || currToken.equals("s") || (expert && (currToken.equals("K")
                || currToken.equals("k") || currToken.equals("I"))))) {
            String validationFailMessage = "First leg ("
                    + currToken
                    + ") must be one of Pg, Pb, Pn, Pdiff, Sg, Sb, Sn, Sdiff, P, S, p, s";
            if(expert) {
                validationFailMessage += ", K, k, I";
            }
            return validationFailMessage;
        }
        for(int i = 1; i < legs.size(); i++) {
            prevToken = currToken;
            currToken = (String)legs.get(i);
            /* Check for 2 reflections with no leg between them. */
            if(currToken.startsWith("^") || currToken.startsWith("v")
                    || currToken.equals("m") || currToken.equals("c")
                    || currToken.equals("i")) {
                if(prevIsReflect) {
                    return "Two reflections with no leg in between: "
                            + prevToken + ", " + currToken;
                } else {
                    prevIsReflect = true;
                }
            } else {
                prevIsReflect = false;
            }
            /* Check for "END" before the end. */
            if(prevToken.equals("END")) {
                return "Legs ended but more tokens exist: " + currToken;
            }
            /* Check for P or S next to I or J */
            if((currToken.startsWith("P") || currToken.startsWith("S")
                    || currToken.startsWith("p") || currToken.startsWith("s")
                    || currToken.equals("m") || currToken.equals("c"))
                    && (prevToken.equals("I") || prevToken.equals("J") || prevToken.equals("i"))) {
                return "Cannot have I,J,i followed by P,S,p,s,m,c: "
                        + prevToken + ", " + currToken;
            }
            if((prevToken.startsWith("P") || prevToken.startsWith("S")
                    || prevToken.startsWith("p") || prevToken.startsWith("s")
                    || prevToken.equals("m") || prevToken.equals("c"))
                    && (currToken.equals("I") || currToken.equals("J") || currToken.equals("i"))) {
                return "Cannot have P,S,p,s,m,c followed by I,J,i: "
                        + prevToken + ", " + currToken;
            }
            /* Check for m next to K. */
            if(prevToken.equals("m") && currToken.equals("K")) {
                return "Cannot have m followed by K";
            }
            if(currToken.equals("m") && prevToken.equals("K")) {
                return "Cannot have K followed by m";
            }
        }
        /* Make sure legs end in "END". */
        if(!currToken.equals("END")) {
            return "Last token must be END";
        }
        return null;
    }
    
    public static Arrival getEarliestArrival(List<SeismicPhase> phases, double degrees) {
        Arrival minArrival = null;
        for (SeismicPhase seismicPhase : phases) {
            seismicPhase.calcTime(degrees);
            Arrival currArrival = seismicPhase.getEarliestArrival(degrees);
            if (currArrival != null && ( minArrival == null || minArrival.getTime() > currArrival.getTime())) {
                minArrival = currArrival;
            }
        }
        return minArrival;
    }

    public String toString() {
        String desc = name + ": ";
        for(int i = 0; i < legs.size(); i++) {
            desc += legs.get(i) + " ";
        }
        desc += "\n";
        for(int i = 0; i < branchSeq.size(); i++) {
            desc += (Integer)branchSeq.get(i) + " ";
        }
        desc += "\n";
        desc += "minRayParam=" + minRayParam + " maxRayParam=" + maxRayParam;
        desc += "\n";
        desc += "minDistance=" + (minDistance * 180.0 / Math.PI)
                + " maxDistance=" + (maxDistance * 180.0 / Math.PI);
        return desc;
    }

    public void dump() {
        for(int j = 0; j < dist.length; j++) {
            System.out.println(j + "  " + dist[j] + "  " + rayParams[j]);
        }
    }

    public static void main(String args[]) {
        TauModel tMod;
        TauModel tModDepth;
        try {
            if(args.length < 3) {
                System.out.println("Usage: SeismicPhase modelfile depth phasename [phasename ...]");
            }
            tMod = TauModel.readModel(args[0]);
            tModDepth = tMod.depthCorrect(Double.valueOf(args[1]).doubleValue());
            for(int i = 2; i < args.length; i++) {
                System.out.println("-----");
                SeismicPhase sp = new SeismicPhase(args[i], tModDepth);
                System.out.println(sp);
                sp.dump();
            }
            System.out.println("-----");
        } catch(FileNotFoundException e) {
            System.out.println(e.getMessage());
        } catch(OptionalDataException e) {
            System.out.println(e.getMessage());
        } catch(StreamCorruptedException e) {
            System.out.println(e.getMessage());
        } catch(IOException e) {
            System.out.println(e.getMessage());
        } catch(ClassNotFoundException e) {
            System.out.println(e.getMessage());
        } catch(TauModelException e) {
            System.out.println(e.getMessage());
            e.printStackTrace();
        }
    }
}
