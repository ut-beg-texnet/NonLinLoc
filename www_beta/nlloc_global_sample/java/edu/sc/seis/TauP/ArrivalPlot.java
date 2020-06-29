package edu.sc.seis.TauP;

import java.awt.Graphics;
import java.util.Vector;

import javax.swing.JPanel;

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
 * ArrivalPlot.java
 * 
 * 
 * Created: Thu Jun 22 13:04:16 EDT 2000
 * 
 * @author Philip Crotwell
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * 
 * 
 * 
 */
public abstract class ArrivalPlot extends JPanel {

    public void setSelectedIndex(int index) {
        this.selectedIndex = index;
        repaint();
    }

    public void addElement(Arrival a) {
        arrivals.addElement(a);
    }

    public void removeAllElements() {
        arrivals.removeAllElements();
        repaint();
    }

    public void setTauModel(TauModel tMod) {
        this.tMod = tMod;
    }

    public abstract void paintBackground(Graphics g);

    public abstract void paintArrivals(Graphics g);

    public abstract void paintForeground(Graphics g);

    public void paintComponent(Graphics g) {
        super.paintComponent(g);
        paintBackground(g);
        paintArrivals(g);
        paintForeground(g);
    }

    protected Vector arrivals = new Vector();

    protected TauModel tMod;

    protected int selectedIndex = -1;
} // ArrivalPlot
