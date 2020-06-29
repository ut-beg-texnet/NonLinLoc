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
 * Convenience class that allows a sac header variable to be associated with a
 * seismic phase name.
 * 
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * 
 * 
 * 
 * @author H. Philip Crotwell
 * 
 */
public class PhaseName implements Serializable {

    /** name of the phase, ie PKIKP */
    public String name;

    /** sac t header to be associated with the phase */
    public int sacTNum = -1;

    public PhaseName(String name) throws TauModelException {
        this.name = name;
        // check name is valid
        SeismicPhase.legPuller(name);
    }

    public PhaseName(String name, int sacTNum) throws TauModelException {
        this(name);
        this.sacTNum = sacTNum;
    }

    public boolean equals(PhaseName obj) {
        if(obj.name.equals(this.name) && obj.sacTNum == this.sacTNum) {
            return true;
        } else {
            return false;
        }
    }

    public String getName() {
        return name;
    }

    public String toString() {
        return name;
    }
}
