package edu.sc.seis.TauP;

import java.util.ArrayList;
import java.util.List;

import javax.swing.table.AbstractTableModel;

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
 * PierceTableModel
 * 
 * 
 * Created: Thu May 6 14:42:17 1999
 * 
 * @author Philip Crotwell
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * 
 * 
 * 
 */
public class PierceTableModel extends AbstractTableModel {

    public PierceTableModel() {}

    public int getRowCount() {
        if(arrivals.size() != 0) {
            return arrivals.get(selectedIndex).getNumPiercePoints();
        } else {
            return 0;
        }
    }

    public int getColumnCount() {
        return 2;
    }

    public Object getValueAt(int row, int col) {
        switch(col){
            case 0:
                double dist = (arrivals.get(selectedIndex).getPiercePoint(row).getDistDeg());
                if(arrivals.get(selectedIndex).getDistDeg() % 360 > 180
                        && dist != 0.0) {
                    dist *= -1.0;
                }
                return Outputs.formatDistance(dist);
            case 1:
                return Outputs.formatDepth(arrivals.get(selectedIndex).getPiercePoint(row).depth);
            default:
                return "";
        }
    }

    public String getColumnName(int col) {
        switch(col){
            case 0:
                return "Dist";
            case 1:
                return "Depth";
            default:
                return "";
        }
    }

    public void setArrivals(List<Arrival> arrivals) {
        this.arrivals = arrivals;
        setSelectedIndex(0);
        fireTableDataChanged();
    }

    public void setSelectedIndex(int index) {
        this.selectedIndex = index;
        fireTableDataChanged();
    }

    private int selectedIndex = 0;

    private List<Arrival> arrivals = new ArrayList<Arrival>();
    
} // PierceTableModel
