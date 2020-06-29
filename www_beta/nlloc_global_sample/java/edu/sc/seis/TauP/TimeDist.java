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
 * Holds the ray parameter, time and distance increments, and optionally a
 * depth, for a ray passing through some layer.
 * 
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * 
 * 
 * 
 * @author H. Philip Crotwell
 */
public class TimeDist implements Cloneable {

    public double p;

    public double depth;

    public double time;

    public double distRadian;

    public TimeDist() {
        this.p = 0;
        this.depth = 0;
        this.time = 0;
        this.distRadian = 0;
    }

    public TimeDist(double p) {
        this.p = p;
        this.depth = 0;
        this.time = 0;
        this.distRadian = 0;
    }

    public TimeDist(double p, double time, double dist) {
        this.p = p;
        this.depth = 0;
        this.time = time;
        this.distRadian = dist;
    }

    public TimeDist(double p, double time, double dist, double depth) {
        this.p = p;
        this.depth = depth;
        this.time = time;
        this.distRadian = dist;
    }

    public void add(TimeDist td) {
        this.time += td.time;
        this.distRadian += td.distRadian;
    }

    public String toString() {
        String desc = "p= " + p + " time=" + time + " dist=" +getDistDeg()+"("+distRadian+" rad) depth="
                + depth;
        return desc;
    }

    
    public double getP() {
        return p;
    }

    
    public double getDepth() {
        return depth;
    }

    
    public double getTime() {
        return time;
    }


    public double getDistRadian() {
        return distRadian;
    }
    
    public double getDistDeg() {
        return Arrival.RtoD * getDistRadian();
    }

    public Object clone() {
        try {
            return super.clone();
        } catch(CloneNotSupportedException e) {
            // Can't happen, but...
            throw new InternalError(e.toString());
        }
    }
}
