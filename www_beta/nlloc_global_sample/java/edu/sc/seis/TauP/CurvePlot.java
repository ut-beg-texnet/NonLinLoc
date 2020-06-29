package edu.sc.seis.TauP;

import java.awt.Color;
import java.awt.Graphics;

/**
 * CurvePlot.java
 * 
 * 
 * Created: Thu Jun 22 13:19:37 2000
 * 
 * @author Philip Crotwell
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * 
 * 
 * 
 */
public class CurvePlot extends ArrivalPlot {

    public CurvePlot() {
        super();
        setOpaque(true);
    }

    public void paintBackground(Graphics g) {}

    public void paintArrivals(Graphics g) {
        for(int i = 0; i < arrivals.size(); i++) {
            paintCurve(g, i);
        }
    }

    public void paintForeground(Graphics g) {
        Color orig = g.getColor();
        g.setColor(Color.red);
        paintCurve(g, selectedIndex);
        g.setColor(orig);
    }

    protected void paintCurve(Graphics g, int i) {
        Arrival a;
        int[] x, y;
        int xOffset = getSize().width / 2;
        int yOffset = getSize().height / 2;
        int pixelRad = Math.min(xOffset, yOffset);
        double roe = 6371;
        a = (Arrival)arrivals.elementAt(i);
        x = new int[a.getNumPathPoints()];
        y = new int[a.getNumPathPoints()];
        for(int j = 0; j < x.length; j++) {
            x[j] = xOffset
                    + (int)Math.rint(Math.sin(a.getPathPoint(j).distRadian)
                            * (roe - a.getPathPoint(j).depth) / roe * pixelRad);
            y[j] = yOffset
                    - (int)Math.rint(Math.cos(a.getPathPoint(j).distRadian)
                            * (roe - a.getPathPoint(j).depth) / roe * pixelRad);
        }
        g.drawPolyline(x, y, x.length);
    }
} // CurvePlot
