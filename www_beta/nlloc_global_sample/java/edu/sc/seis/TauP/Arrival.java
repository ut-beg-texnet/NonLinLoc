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
 * convenience class for storing the parameters associated with a phase arrival.
 *
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 *
 *
 *
 * @author H. Philip Crotwell
 *
 */
public class Arrival {

    public Arrival(SeismicPhase phase,
                   double time,
                   double dist,
                   double rayParam,
                   int rayParamIndex,
                   String name,
                   String puristName,
                   double sourceDepth,
                   double takeoffAngle,
                   double incidentAngle) {
        this.phase = phase;
        this.time = time;
        this.dist = dist;
        this.rayParam = rayParam;
        this.rayParamIndex = rayParamIndex;
        this.name = name;
        this.puristName = puristName;
        this.sourceDepth = sourceDepth;
        this.takeoffAngle = takeoffAngle;
        this.incidentAngle = incidentAngle;
    }

    /** phase that generated this arrival. */
    protected SeismicPhase phase;

    /** travel time in seconds */
    protected double time;

    /** angular distance (great circle) in radians */
    protected double dist;

    /** ray parameter in seconds per radians. */
    protected double rayParam;

    protected int rayParamIndex;

    /** phase name */
    protected String name;

    /** phase name changed for true depths */
    protected String puristName;

    /** source depth in kilometers */
    protected double sourceDepth;

    /** pierce and path points */
    protected TimeDist[] pierce, path;

    protected double incidentAngle;

    protected double takeoffAngle;

    // get set methods
    /** @return the phase used to calculate this arrival. */
    public SeismicPhase getPhase() {
        return phase;
    }

    /** @return travel time in seconds */
    public double getTime() {
        return time;
    }

    /** returns travel distance in radians */
    public double getDist() {
        return dist;
    }

    /**
     * returns travel distance in degrees.
     */
    public double getDistDeg() {
        return RtoD * getDist();
    }

    /**
     * returns distance in radians and in the range 0-PI. Note this may not be
     * the actual distance traveled.
     */
    public double getModuloDist() {
        double moduloDist = getDist() % TWOPI;
        if(moduloDist > Math.PI) {
            moduloDist = TWOPI - moduloDist;
        }
        return moduloDist;
    }

    /**
     * returns distance in degrees and in the range 0-180. Note this may not be
     * the actual distance traveled.
     */
    public double getModuloDistDeg() {
        double moduloDist = (RtoD * getDist()) % 360;
        if(moduloDist > 180) {
            moduloDist = 360 - moduloDist;
        }
        return moduloDist;
    }

    /** returns ray parameter in seconds per radian */
    public double getRayParam() {
        return rayParam;
    }

    /** returns ray parameter in seconds per deg */
    public double getRayParamDeg() {
        return getRayParam()/RtoD;
    }

    public double getIncidentAngle() {
        return incidentAngle;
    }

    public double getTakeoffAngle() {
        return takeoffAngle;
    }

    public int getRayParamIndex() {
        return rayParamIndex;
    }

    /** returns phase name */
    public String getName() {
        return name;
    }

    /**
     * returns purist's version of name. Depths are changed to reflect the true
     * depth of the interface.
     */
    public String getPuristName() {
        return puristName;
    }

    /** returns source depth in kilometers */
    public double getSourceDepth() {
        return sourceDepth;
    }

    /** returns pierce points as TimeDist objects. */
    public TimeDist[] getPierce() {
        if (pierce == null) {
            this.pierce = getPhase().calcPierce(this).getPierce();
        }
        return pierce;
    }

    /** returns pierce points as TimeDist objects. */
    public TimeDist[] getPath() {
        if (path == null) {
            this.path = getPhase().calcPath(this).getPath();
        }
        return path;
    }

    public String toString() {
        String desc =  Outputs.formatDistance(getModuloDistDeg()) + Outputs.formatDepth(getSourceDepth()) + "   " + getName()
                + "  " + Outputs.formatTime(getTime()) + "  " + Outputs.formatRayParam(Math.PI / 180.0 * getRayParam())
                + "  " + Outputs.formatDistance(getTakeoffAngle()) + " " + Outputs.formatDistance(getIncidentAngle())
                + " " + Outputs.formatDistance(getDistDeg());
        if (getName().equals(getPuristName())) {
            desc += "   = ";
        } else {
            desc += "   * ";
        }
        desc += getPuristName();
        return desc;
    }

    public int getNumPiercePoints() {
        if(pierce != null) {
            return pierce.length;
        } else {
            return 0;
        }
    }

    public int getNumPathPoints() {
        if(path != null) {
            return path.length;
        } else {
            return 0;
        }
    }

    public TimeDist getPiercePoint(int i) {
        // don't check for i> length since we want an ArrayOutOfBounds anyway
        return pierce[i];
    }

    /**
     * finds the first pierce point at the given depth.
     *
     * @throws ArrayIndexOutOfBoundsException
     *             if depth is not found
     */
    public TimeDist getFirstPiercePoint(double depth) {
        for(int i = 0; i < pierce.length; i++) {
            if(pierce[i].depth == depth) {
                return pierce[i];
            }
        }
        throw new ArrayIndexOutOfBoundsException("No Pierce point found for depth "
                + depth);
    }

    /**
     * finds the last pierce point at the given depth.
     *
     * @throws ArrayIndexOutOfBoundsException
     *             if depth is not found
     */
    public TimeDist getLastPiercePoint(double depth) {
        TimeDist piercepoint = null;
        for(int i = 0; i < pierce.length; i++) {
            if(pierce[i].depth == depth) {
                piercepoint = pierce[i];
            }
        }
        if(piercepoint == null) {
            throw new ArrayIndexOutOfBoundsException("No Pierce point found for depth "
                    + depth);
        }
        return piercepoint;
    }

    public TimeDist getPathPoint(int i) {
        // don't check for i> length since we want an ArrayOutOfBounds anyway
        return path[i];
    }

    protected static final double TWOPI = 2.0 * Math.PI;

    protected static final double DtoR = Math.PI / 180.0;

    protected static final double RtoD = 180.0 / Math.PI;
}
