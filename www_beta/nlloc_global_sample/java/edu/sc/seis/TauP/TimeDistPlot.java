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
import java.awt.Graphics;
import java.util.Vector;

/**
 * Time versus Distance plot.
 * 
 * 
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * 
 * 
 * 
 * @author H. Philip Crotwell
 * 
 */
public class TimeDistPlot extends XYPlot {

    public TimeDistPlot(Container parent) {
        super(parent);
        title = "Time/Distance";
        xLabel = "distance (deg)";
        yLabel = "time (sec)";
        xTickWidth = 10;
        yTickWidth = 200;
        minX = 0.0;
        maxX = 270.0;
        minY = 0.0;
        maxY = 2800.0;
    }

    public TimeDistPlot(Container parent, int width, int height) {
        super(parent, width, height);
        title = "Time/Distance";
        xLabel = "distance (deg)";
        yLabel = "time (sec)";
        xTickWidth = 10;
        yTickWidth = 200;
    }

    public void plot(TauModel tModel, boolean isPWave) {
        int waveNum;
        if(!isPWave) {
            waveNum = 1;
            minX = 0.0;
            maxX = 270.0;
            minY = 0.0;
            maxY = 2800.0;
        } else {
            waveNum = 0;
            minX = 0.0;
            maxX = 200.0;
            minY = 0.0;
            maxY = 1300.0;
        }
        zoomMinX = minX;
        zoomMaxX = maxX;
        zoomMinY = minY;
        zoomMaxY = maxY;
        int jj = 0;
        double x, y;
        double[] tempXData, tempYData;
        xSegments = new Vector();
        ySegments = new Vector();
        for(int i = 0; i < tModel.getNumBranches(); i++) {
            jj = 0;
            tempXData = new double[tModel.rayParams.length];
            tempYData = new double[tModel.rayParams.length];
            for(int j = 0; j < tModel.rayParams.length; j++) {
                y = 0;
                x = 0;
                for(int k = 0; k <= i; k++) {
                    x += 2 * tModel.tauBranches[waveNum][k].getDist(j);
                    y += 2 * tModel.tauBranches[waveNum][k].time[j];
                }
                x *= 180 / Math.PI;
                if(tModel.tauBranches[waveNum][i].time[j] != 0
                        || (i == 0 && x == 0)) {
                    // System.out.println("branch "+i+" dist="+x+" time="+y+
                    // " p="+tModel.rayParams[j]);
                    tempXData[jj] = x;
                    tempYData[jj] = y;
                    if(DEBUG)
                        System.out.println(x + " " + y);
                    jj++;
                }
            }
            if(DEBUG)
                System.out.println("> ");
            xData = new double[jj];
            System.arraycopy(tempXData, 0, xData, 0, jj);
            yData = new double[jj];
            System.arraycopy(tempYData, 0, yData, 0, jj);
            xSegments.addElement(xData);
            ySegments.addElement(yData);
        }
        repaint();
    }

    public void paint(Graphics g) {
        super.paint(g);
    }
}
