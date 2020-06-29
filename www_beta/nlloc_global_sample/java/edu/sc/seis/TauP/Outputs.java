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

import java.util.Properties;

/**
 * Outputs.java contains formating, similar to printf, routines for the output
 * types in the TauP package.
 * 
 * 
 * Created: Tue Sep 21 11:45:35 1999
 * 
 * @author Philip Crotwell
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * 
 * 
 * 
 */
public class Outputs {
    
    public static void configure(Properties props) {
        String formString;
        formString = "%8." + props.getProperty("taup.depth.precision", "1")
                + "f";
        depthFormat = new Format(formString);
        formString = "%8." + props.getProperty("taup.distance.precision", "2")
                + "f";
        distanceFormat = new Format(formString);
        formString = "%8." + props.getProperty("taup.time.precision", "2")
                + "f";
        timeFormat = new Format(formString);
        formString = "%8." + props.getProperty("taup.rayparam.precision", "3")
                + "f";
        rayParamFormat = new Format(formString);
        formString = "%8." + props.getProperty("taup.latlon.precision", "2")
                + "f";
        latLonFormat = new Format(formString);
    }

    public static String formatDepth(double depth) {
        return depthFormat.form(depth);
    }

    public static String formatDistance(double distance) {
        return distanceFormat.form(distance);
    }

    public static String formatTime(double time) {
        return timeFormat.form(time);
    }

    public static String formatRayParam(double rayParam) {
        return rayParamFormat.form(rayParam);
    }

    public static String formatLatLon(double latlon) {
        return latLonFormat.form(latlon);
    }

    protected static Format depthFormat = new Format("%8.1f");

    protected static Format distanceFormat = new Format("%8.2f");;

    protected static Format timeFormat = new Format("%8.2f");;

    protected static Format rayParamFormat = new Format("%8.3f");;

    protected static Format latLonFormat = new Format("%8.2f");;
} // Outputs
