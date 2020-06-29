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
import java.awt.Event;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Insets;
import java.util.Vector;

/**
 * Simple y versus x plot widget.
 * 
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * 
 * 
 * 
 * @author H. Philip Crotwell
 */
public class XYPlot extends Canvas {

    Container pappy;

    Dimension minSize;

    int width, height;

    Vector xSegments, ySegments;

    double[] xData, yData;

    int xOffset = 40;

    int yOffset = 40;

    int xTickWidth = 100;

    int yTickWidth = 100;

    double minX = 0.0;

    double zoomMinX = minX;

    double maxX = 10.0;

    double zoomMaxX = maxX;

    double minY = 0.0;

    double zoomMinY = minY;

    double maxY = 10.0;

    double zoomMaxY = maxY;

    private double xScale;

    private double yScale;

    String title = "Title";

    String xLabel = "X Label";

    String yLabel = "Y Label";

    double mouseDownX, mouseDownY;

    boolean DEBUG = false;

    public XYPlot(Container parent) {
        pappy = parent;
        width = 600;
        height = 600;
        minSize = new Dimension(width, height);
    }

    public XYPlot(Container parent, int width, int height) {
        pappy = parent;
        this.width = width;
        this.height = height;
        minSize = new Dimension(width, height);
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
        xScale = (width - 2.0 * xOffset) / (zoomMaxX - zoomMinX);
        yScale = (height - 2.0 * yOffset) / (zoomMaxY - zoomMinY);
        if(DEBUG)
            System.out.println("xOffset " + xOffset + " yOffset " + yOffset
                    + " minX " + minX + " maxX " + maxX + " minY " + minY
                    + " maxY " + maxY + " xScale " + xScale + " yScale "
                    + yScale);
        g.drawRect(0, 0, width - 1, height - 1);
        g.drawLine(xOffset, yOffset, xOffset, height - yOffset);
        g.drawLine(xOffset, height - yOffset, width - xOffset, height - yOffset);
        g.drawString(title,
                     (width - fontMetrics.stringWidth(title)) / 2,
                     yOffset - fontMetrics.getDescent()
                             - fontMetrics.getLeading());
        g.drawString(xLabel,
                     (width - fontMetrics.stringWidth(xLabel)) / 2,
                     height - fontMetrics.getDescent()
                             - fontMetrics.getLeading());
        for(int i = 0; i * xTickWidth >= zoomMinX && i * xTickWidth <= zoomMaxX; i++) {
            g.drawLine(xOffset + (int)(xScale * (i * xTickWidth - zoomMinX)),
                       height - yOffset,
                       xOffset + (int)(xScale * (i * xTickWidth - zoomMinX)),
                       height - (int)(3.0 / 4 * yOffset));
            g.drawString("" + (i * xTickWidth), xOffset
                    + (int)(xScale * (i * xTickWidth - zoomMinX)), height
                    - fontMetrics.getAscent() - fontMetrics.getLeading());
        }
        for(int i = 0; i * yTickWidth >= zoomMinY && i * yTickWidth <= zoomMaxY; i++) {
            g.drawLine(xOffset,
                       height - yOffset
                               - (int)(yScale * (i * yTickWidth - zoomMinY)),
                       (int)(3.0 / 4 * xOffset),
                       height - yOffset
                               - (int)(yScale * (i * yTickWidth - zoomMinY)));
            g.drawString("" + (i * yTickWidth), 0, height - yOffset
                    - (int)(yScale * (i * yTickWidth - zoomMinY)));
        }
        plotData(g);
    }

    public void plotData(Graphics g) {
        if(xSegments != null && ySegments != null) {
            if(xSegments.size() != ySegments.size()) {
                System.out.println("xSegments.size() != ySegments.size()");
            }
            for(int segNum = 0; segNum < xSegments.size(); segNum++) {
                g.setColor(colorForSegment(segNum));
                xData = (double[])xSegments.elementAt(segNum);
                yData = (double[])ySegments.elementAt(segNum);
                if(xData != null) {
                    if(xData.length == yData.length && xData.length > 0) {
                        for(int i = 0; i < xData.length - 1; i++)
                            g.drawLine(xOffset
                                               + (int)(xScale * (xData[i] - zoomMinX)),
                                       height
                                               - yOffset
                                               - (int)(yScale * (yData[i] - zoomMinY)),
                                       xOffset
                                               + (int)(xScale * (xData[i + 1] - zoomMinX)),
                                       height
                                               - yOffset
                                               - (int)(yScale * (yData[i + 1] - zoomMinY)));
                    }
                } else {
                    System.out.println("null data");
                }
            }
        }
        g.setColor(Color.black);
    }

    public Color colorForSegment(int segNum) {
        if(xSegments.size() == 1)
            return Color.black;
        int colorInc = (int)Math.floor(255 / (xSegments.size() - 1));
        return new Color(segNum * colorInc,
                         255 - segNum * colorInc,
                         (segNum % 3) * 127);
    }

    public void plot() {
        xData = new double[2];
        yData = new double[2];
        xData[0] = minX;
        yData[0] = minY;
        xData[1] = maxX;
        yData[1] = maxY;
        xSegments.addElement(xData);
        ySegments.addElement(yData);
        xData = null;
        yData = null;
    }

    public boolean mouseDown(Event evt, int x, int y) {
        mouseDownX = x;
        mouseDownY = y;
        return true;
    }

    public boolean mouseUp(Event evt, int mouseUpX, int mouseUpY) {
        /* If the difference is too small then it must have been a mistake. */
        if(Math.abs(mouseDownX - mouseUpX) + Math.abs(mouseDownY - mouseUpY) > 5) {
            double tempDown, tempUp;
            tempDown = (mouseDownX - xOffset) / xScale + zoomMinX;
            tempUp = (mouseUpX - xOffset) / xScale + zoomMinX;
            if(tempDown > tempUp) {
                zoomMinX = tempUp;
                zoomMaxX = tempDown;
            } else {
                zoomMinX = tempDown;
                zoomMaxX = tempUp;
            }
            tempDown = (mouseDownY - yOffset) / yScale + zoomMinY;
            tempUp = (mouseUpY - yOffset) / yScale + zoomMinY;
            if(tempDown > tempUp) {
                zoomMinY = tempUp;
                zoomMaxY = tempDown;
            } else {
                zoomMinY = tempDown;
                zoomMaxY = tempUp;
            }
            repaint();
        }
        return true;
    }
}
