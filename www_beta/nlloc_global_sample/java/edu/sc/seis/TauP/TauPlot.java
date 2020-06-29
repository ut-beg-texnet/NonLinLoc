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
 * generic xy plot. Probably should be subclass to get better behavior.
 * 
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * 
 * 
 * 
 * @author H. Philip Crotwell
 * 
 */
public class TauPlot extends XYPlot {

    public TauPlot(Container parent) {
        super(parent);
        title = "Tau";
        xLabel = "p (sec/km or km-sec/km)";
        yLabel = "tau (sec)";
    }

    public TauPlot(Container parent, int width, int height) {
        super(parent, width, height);
        title = "Tau";
        xLabel = "p (sec/km or km-sec/km)";
        yLabel = "depth (sec)";
    }

    public void plot(TauModel tModel, boolean isPWave) {
        int waveNum;
        if(!isPWave) {
            waveNum = 1;
            minX = 0.0;
            zoomMinX = minX;
            maxX = 2500.0;
            zoomMaxX = maxX;
            minY = 0.0;
            zoomMinY = minY;
            maxY = 2500.0;
            zoomMaxY = maxY;
        } else {
            waveNum = 0;
            minX = 0.0;
            zoomMinX = minX;
            maxX = 1300.0;
            zoomMaxX = maxX;
            minY = 0.0;
            zoomMinY = minY;
            maxY = 1300.0;
            zoomMaxY = maxY;
        }
        int jj = 0;
        double x, y;
        double[] tempXData, tempYData;
        xSegments = new Vector();
        ySegments = new Vector();
        for(int i = 0; i < tModel.tauBranches[0].length; i++) {
            jj = 0;
            tempXData = new double[tModel.rayParams.length];
            tempYData = new double[tModel.rayParams.length];
            for(int j = 0; j < tModel.rayParams.length; j++) {
                x = tModel.rayParams[j];
                y = 0;
                for(int k = 0; k <= i; k++) {
                    y += 2 * tModel.tauBranches[waveNum][k].tau[j];
                }
                if((y != 0 || x == 0)
                        && tModel.tauBranches[waveNum][i].tau[j] != 0.0) {
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
