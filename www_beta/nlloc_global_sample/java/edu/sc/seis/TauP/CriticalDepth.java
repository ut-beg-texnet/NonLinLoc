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
 * Utility class to keep track of criticalpoints (discontinuities or reversals
 * in slowness gradient) within slowness and velocity models.
 * 
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * 
 * 
 * 
 * @author H. Philip Crotwell
 */
public class CriticalDepth implements Cloneable, Serializable {

    /** depth in kilometers at which there is a critical point. */
    private double depth;

    /** layer number within the velocity model with this depth at its top. */
    private int velLayerNum;

    /**
     * slowness layer for P waves with this depth at its top. This can be
     * PLayers.size() for the last critical layer.
     */
    private int PLayerNum;

    /**
     * slowness layer for S waves with this depth at its top. This can be
     * SLayers.size() for the last critical layer.
     */
    private int SLayerNum;

    // Constructors
    public CriticalDepth() {}

    public CriticalDepth(double depth,
                         int velLayerNum,
                         int PLayerNum,
                         int SLayerNum) {
        this.setDepth(depth);
        this.velLayerNum = velLayerNum;
        this.PLayerNum = PLayerNum;
        this.SLayerNum = SLayerNum;
    }

    // Accessor methods
    public void setVelLayerNum(int layerNum) {
        velLayerNum = layerNum;
    }

    public void setPLayerNum(int layerNum) {
        PLayerNum = layerNum;
    }

    public void setSLayerNum(int layerNum) {
        SLayerNum = layerNum;
    }

    public int getVelLayerNum() {
        return velLayerNum;
    }

    public int getPLayerNum() {
        return PLayerNum;
    }

    public int getSLayerNum() {
        return SLayerNum;
    }

    /**
     * sets slowness layer for waveType waves with this depth at its top.
     */
    public void setLayerNum(int layerNum, boolean isPWave) {
        if(isPWave) {
            PLayerNum = layerNum;
        } else {
            SLayerNum = layerNum;
        }
    }

    /**
     * @return slowness layer for waveType waves with this depth at its top.
     */
    public int getLayerNum(boolean isPWave) {
        if(isPWave) {
            return PLayerNum;
        } else {
            return SLayerNum;
        }
    }

    public Object clone() {
        CriticalDepth newObject;
        try {
            newObject = (CriticalDepth)super.clone();
            return newObject;
        } catch(CloneNotSupportedException e) {
            // Can't happen, but...
            System.err.println("Caught CloneNotSupportedException: "
                    + e.getMessage());
            throw new InternalError(e.toString());
        }
    }

    public void setDepth(double depth) {
        this.depth = depth;
    }

    public double getDepth() {
        return depth;
    }
}
