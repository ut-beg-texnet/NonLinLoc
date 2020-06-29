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
 * Convenience class for storing a depth range. It has a top and a bottom and
 * can have an associated ray parameter.
 * 
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * 
 * 
 * 
 * @author H. Philip Crotwell
 */
public class DepthRange implements Serializable, Cloneable {

    /** Top of the depth range. */
    public double topDepth;

    /** Bottom of the depth range. */
    public double botDepth;

    /**
     * rayParam associated with the depth range. If this were a high slowness
     * depth range, then rayParam would be the largest ray parameter that would
     * penetrate the depth range.
     */
    public double rayParam = -1;

    public DepthRange() {}

    public DepthRange(double topDepth, double botDepth) {
        this.topDepth = topDepth;
        this.botDepth = botDepth;
    }

    public DepthRange(double topDepth, double botDepth, double rayParam) {
        this.topDepth = topDepth;
        this.botDepth = botDepth;
        this.rayParam = rayParam;
    }

    public Object clone() {
        DepthRange newObject;
        try {
            newObject = (DepthRange)super.clone();
            return newObject;
        } catch(CloneNotSupportedException e) {
            // Can't happen, but...
            System.err.println("Caught CloneNotSupportedException: "
                    + e.getMessage());
            throw new InternalError(e.toString());
        }
    }

    public String toString() {
        if(rayParam != -1) {
            return "topDepth=" + topDepth + " botDepth=" + botDepth
                    + " rayParam=" + rayParam;
        } else {
            return "topDepth=" + topDepth + " botDepth=" + botDepth;
        }
    }
}
