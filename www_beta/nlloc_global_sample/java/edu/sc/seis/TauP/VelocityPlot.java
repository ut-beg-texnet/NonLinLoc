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

import java.awt.Container;
import java.util.Vector;

/**
 * Velocity versus depth plot.
 * 
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * 
 * 
 * 
 * @author H. Philip Crotwell
 */
public class VelocityPlot extends XYPlot {

    public VelocityPlot(Container parent) {
        super(parent);
        title = "Velocity";
        xLabel = "velocity (km/sec)";
        yLabel = "depth (km)";
        xTickWidth = 1;
        yTickWidth = 500;
    }

    public VelocityPlot(Container parent, int width, int height) {
        super(parent, width, height);
        title = "Velocity";
        xLabel = "velocity (km/sec)";
        yLabel = "depth (km)";
        xTickWidth = 1;
        yTickWidth = 500;
    }

    public void plot(VelocityModel velModel, char waveTypeA, char waveTypeB)
            throws NoSuchMatPropException {
        plot(velModel, waveTypeA);
        xData = new double[2 * velModel.getNumLayers()];
        yData = new double[2 * velModel.getNumLayers()];
        int j = 0;
        for(int i = 0; i < velModel.getNumLayers(); i++) {
            yData[j] = velModel.radiusOfEarth - velModel.depthAtTop(i);
            xData[j] = velModel.evaluateAtTop(i, waveTypeB);
            if(DEBUG)
                System.out.println("x " + xData[j] + " y " + yData[j]);
            j++;
            yData[j] = velModel.radiusOfEarth - velModel.depthAtBottom(i)
                    - minY;
            xData[j] = velModel.evaluateAtBottom(i, waveTypeB);
            if(DEBUG)
                System.out.println("x " + xData[j] + " y " + yData[j]);
            j++;
        }
        xSegments.addElement(xData);
        ySegments.addElement(yData);
        xData = null;
        yData = null;
    }

    public void plot(VelocityModel velModel, char waveType)
            throws NoSuchMatPropException {
        minX = 0.0;
        zoomMinX = minX;
        maxX = 15.0;
        zoomMaxX = maxX;
        minY = 0.0;
        zoomMinY = minY;
        maxY = velModel.radiusOfEarth;
        zoomMaxY = maxY;
        xSegments = new Vector();
        ySegments = new Vector();
        xData = new double[2 * velModel.getNumLayers()];
        yData = new double[2 * velModel.getNumLayers()];
        int j = 0;
        for(int i = 0; i < velModel.getNumLayers(); i++) {
            yData[j] = velModel.radiusOfEarth - velModel.depthAtTop(i);
            xData[j] = velModel.evaluateAtTop(i, waveType);
            if(DEBUG)
                System.out.println(xData[j] + " " + yData[j]);
            j++;
            yData[j] = velModel.radiusOfEarth - velModel.depthAtBottom(i)
                    - minY;
            xData[j] = velModel.evaluateAtBottom(i, waveType);
            if(DEBUG)
                System.out.println(xData[j] + " " + yData[j]);
            j++;
        }
        if(DEBUG)
            System.out.println("> ");
        xSegments.addElement(xData);
        ySegments.addElement(yData);
        xData = null;
        yData = null;
        repaint();
    }
}
