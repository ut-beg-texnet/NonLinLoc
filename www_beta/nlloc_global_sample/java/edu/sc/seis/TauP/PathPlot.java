package edu.sc.seis.TauP;

import java.awt.Color;
import java.awt.Graphics;

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
/**
 * PathPlot.java
 * 
 * 
 * Created: Fri May 7 15:45:43 1999
 * 
 * @author Philip Crotwell
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * 
 * 
 * 
 */
public class PathPlot extends ArrivalPlot {

    public PathPlot() {
        super();
        // setDebugGraphicsOptions(DebugGraphics.LOG_OPTION);
        setOpaque(true);
        // setBackground(java.awt.Color.white);
    }

    public void paintBackground(Graphics g) {
        int xOffset = getSize().width / 2;
        int yOffset = getSize().height / 2;
        int pixelRad = Math.min(xOffset, yOffset);
        double roe = 6371;
        Color origColor = g.getColor();
        Color aColor, bColor, fillColor;
        aColor = new Color(220, 220, 220);
        bColor = Color.white;
        fillColor = aColor;
        boolean whichColor = true;
        int disconRad;
        if(tMod != null) {
            roe = tMod.getRadiusOfEarth();
            disconRad = pixelRad;
            fillColor = whichColor ? aColor : bColor;
            g.setColor(fillColor);
            whichColor = !whichColor;
            g.fillOval(xOffset - disconRad,
                       yOffset - disconRad,
                       2 * disconRad,
                       2 * disconRad);
            for(int i = 0; i < tMod.tauBranches[0].length - 1; i++) {
                if(tMod.tauBranches[0][i].getBotDepth() != tMod.tauBranches[1][i].getBotDepth()) {
                    disconRad = (int)Math.round((tMod.getRadiusOfEarth() - tMod.tauBranches[0][i].getBotDepth())
                            / roe * pixelRad);
                    fillColor = whichColor ? aColor : bColor;
                    g.setColor(fillColor);
                    whichColor = !whichColor;
                    g.fillOval(xOffset - disconRad,
                               yOffset - disconRad,
                               2 * disconRad,
                               2 * disconRad);
                } else {
                    disconRad = (int)Math.round((tMod.getRadiusOfEarth() - tMod.tauBranches[1][i].getBotDepth())
                            / roe * pixelRad);
                    fillColor = whichColor ? aColor : bColor;
                    g.setColor(fillColor);
                    whichColor = !whichColor;
                    g.fillOval(xOffset - disconRad,
                               yOffset - disconRad,
                               2 * disconRad,
                               2 * disconRad);
                }
            }
        }
        g.setColor(origColor);
    }

    public void paintArrivals(Graphics g) {
        for(int i = 0; i < arrivals.size(); i++) {
            paintPaths(g, i);
        }
    }

    public void paintForeground(Graphics g) {
        if(selectedIndex >= 0 && selectedIndex < arrivals.size()) {
            Color orig = g.getColor();
            g.setColor(Color.red);
            paintPaths(g, selectedIndex);
            g.setColor(orig);
        }
    }

    protected void paintPaths(Graphics g, int i) {
        Arrival a;
        int[] x, y;
        int xOffset = getSize().width / 2;
        int yOffset = getSize().height / 2;
        int pixelRad = Math.min(xOffset, yOffset);
        double roe = 6371;
        a = (Arrival)arrivals.elementAt(i);
        x = new int[a.getNumPathPoints()];
        y = new int[a.getNumPathPoints()];
        if((a.getDist() * 180 / Math.PI) % 360 > 180) {
            // long way around
            for(int j = 0; j < x.length; j++) {
                x[j] = xOffset
                        + (int)Math.rint(Math.sin(-1 * a.getPathPoint(j).distRadian)
                                * (roe - a.getPathPoint(j).depth) / roe
                                * pixelRad);
                y[j] = yOffset
                        - (int)Math.rint(Math.cos(-1 * a.getPathPoint(j).distRadian)
                                * (roe - a.getPathPoint(j).depth) / roe
                                * pixelRad);
            }
        } else {
            for(int j = 0; j < x.length; j++) {
                x[j] = xOffset
                        + (int)Math.rint(Math.sin(a.getPathPoint(j).distRadian)
                                * (roe - a.getPathPoint(j).depth) / roe
                                * pixelRad);
                y[j] = yOffset
                        - (int)Math.rint(Math.cos(a.getPathPoint(j).distRadian)
                                * (roe - a.getPathPoint(j).depth) / roe
                                * pixelRad);
            }
        }
        g.drawPolyline(x, y, x.length);
    }
} // PathPlot
