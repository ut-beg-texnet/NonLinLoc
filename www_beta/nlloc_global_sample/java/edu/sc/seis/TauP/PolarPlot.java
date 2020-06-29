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

import java.awt.Canvas;
import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Insets;
import java.util.Vector;

/**
 * Simple polar plot widget.
 * 
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * 
 * 
 * 
 * @author H. Philip Crotwell
 */
public class PolarPlot extends Canvas {

    Container pappy;

    Dimension minSize;

    int width; // circles, so height=width

    protected Vector segments = new Vector();

    double[] circleRadius;

    double outsideRadius = 6371.0;

    int centerX, centerY;

    public static final short FULL = 0;

    public static final short HALF = 1;

    public static final short QUARTER = 2;

    short displayMode;

    public PolarPlot(Container parent) {
        pappy = parent;
        width = 250;
        minSize = new Dimension(width, width);
        setDisplayMode(FULL);
    }

    public PolarPlot(Container parent, int width) {
        pappy = parent;
        this.width = width;
        minSize = new Dimension(width, width);
        setDisplayMode(FULL);
    }

    public void setCircles(double[] circles) {
        circleRadius = new double[circles.length];
        for(int i = 0; i < circles.length; i++) {
            circleRadius[i] = circles[i];
        }
    }

    public void setOutsideRadius(double r) {
        outsideRadius = r;
    }

    public void appendSegment(TimeDist[] td) {
        segments.addElement(td);
    }

    public void clearSegments() {
        segments.removeAllElements();
    }

    public void setDisplayMode(short mode) {
        displayMode = mode;
        switch(displayMode){
            case FULL:
                centerX = width / 2;
                centerY = width / 2;
                break;
            case HALF:
                centerX = width / 2;
                centerY = width;
                break;
            case QUARTER:
                centerX = 0;
                centerY = width;
                break;
        }
    }

    public Dimension preferredSize() {
        return minimumSize();
    }

    public synchronized Dimension minimumSize() {
        return minSize;
    }

    public Insets insets() {
        return new Insets(5, 5, 5, 5);
    }

    public void paint(Graphics g) {
        FontMetrics fontMetrics = g.getFontMetrics();
        g.drawRect(0, 0, width - 1, width - 1);
        plotCircles(g);
        plotData(g);
    }

    protected void plotCircles(Graphics g) {
        int radius;
        if(circleRadius != null) {
            switch(displayMode){
                case FULL:
                    for(int i = 0; i < circleRadius.length; i++) {
                        radius = (int)Math.rint(circleRadius[i] / outsideRadius
                                * width / 2);
                        g.drawOval(width / 2 - radius,
                                   width / 2 - radius,
                                   2 * radius,
                                   2 * radius);
                    }
                    break;
                case HALF:
                    for(int i = 0; i < circleRadius.length; i++) {
                        radius = (int)Math.rint(circleRadius[i] / outsideRadius
                                * width / 2);
                        g.drawArc(width / 2 - radius,
                                  width - radius,
                                  2 * radius,
                                  2 * radius,
                                  0,
                                  180);
                    }
                    break;
                case QUARTER:
                    for(int i = 0; i < circleRadius.length; i++) {
                        radius = (int)Math.rint(circleRadius[i] / outsideRadius
                                * width);
                        g.drawArc(-1 * radius,
                                  width - radius,
                                  2 * radius,
                                  2 * radius,
                                  0,
                                  90);
                    }
                    break;
            }
        }
    }

    protected void plotData(Graphics g) {
        int numSegments = segments.size();
        double scale = centerY / outsideRadius;
        for(int segNum = 0; segNum < numSegments; segNum++) {
            g.setColor(colorForSegment(segNum));
            TimeDist[] data = (TimeDist[])segments.elementAt(segNum);
            for(int i = 0; i < data.length - 1; i++) {
                g.drawLine(centerX
                                   - (int)Math.rint((outsideRadius - data[i].depth)
                                           * Math.cos(data[i].distRadian + Math.PI
                                                   / 2) * scale),
                           centerY
                                   - (int)Math.rint((outsideRadius - data[i].depth)
                                           * Math.sin(data[i].distRadian + Math.PI
                                                   / 2) * scale),
                           centerX
                                   - (int)Math.rint((outsideRadius - data[i + 1].depth)
                                           * Math.cos(data[i + 1].distRadian
                                                   + Math.PI / 2) * scale),
                           centerY
                                   - (int)Math.rint((outsideRadius - data[i + 1].depth)
                                           * Math.sin(data[i + 1].distRadian
                                                   + Math.PI / 2) * scale));
            }
        }
        g.setColor(Color.black);
    }

    public Color colorForSegment(int segNum) {
        if(segments.size() == 1)
            return Color.black;
        int colorInc = (int)Math.floor(255 / (segments.size() - 1));
        return new Color(segNum * colorInc,
                         255 - segNum * colorInc,
                         (segNum % 3) * 127);
    }
}
