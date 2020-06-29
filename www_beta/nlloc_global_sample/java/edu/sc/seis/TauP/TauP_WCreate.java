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

import java.applet.Applet;
import java.awt.Button;
import java.awt.CardLayout;
import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.Frame;
import java.awt.Panel;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.io.IOException;

/**
 * GUI for TauP_Create, a re-implementation of the seismic travel time
 * calculation method described in "The Computation of Seismic Travel Times" by
 * Buland and Chapman, BSSA vol. 73, No. 5, October 1983, pp 1271-1302
 * 
 * TauP_WCreate can run as an applet or as a java application, except for the
 * security restrictions on file access.
 * 
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * 
 * 
 * 
 * @author H. Philip Crotwell
 */
public class TauP_WCreate extends Applet implements ActionListener {

    VelocityPlot velocityPlot;

    SlownessPlot slownessPlot;

    TauPlot tauPPlot;

    DistPlot distPlot;

    TimePlot timePlot;

    TimeDistPlot timeDistPlot;

    String modelFilename = "iasp91.tvel";

    String directory = ".";

    TauP_Create tauPCreate = new TauP_Create();

    Panel cardPanel = new Panel();

    CardLayout cards = new CardLayout();

    Button nextButton;

    Button previousButton;

    Button quitButton;

    public void init() {
        setBackground(Color.white);
        setLayout(new FlowLayout());
        cardPanel.setLayout(cards);
        add(cardPanel);
        previousButton = new Button("Previous");
        add(previousButton);
        previousButton.addActionListener(this);
        nextButton = new Button("Next");
        add(nextButton);
        nextButton.addActionListener(this);
        quitButton = new Button("Quit");
        add(quitButton);
        quitButton.addActionListener(this);
        velocityPlot = new VelocityPlot(this);
        velocityPlot.DEBUG = false;
        cardPanel.add("Velocity", velocityPlot);
        slownessPlot = new SlownessPlot(this);
        slownessPlot.DEBUG = false;
        cardPanel.add("Slowness", slownessPlot);
        tauPPlot = new TauPlot(this);
        tauPPlot.DEBUG = false;
        cardPanel.add("Tau", tauPPlot);
        distPlot = new DistPlot(this);
        distPlot.DEBUG = false;
        cardPanel.add("Dist", distPlot);
        timePlot = new TimePlot(this);
        timePlot.DEBUG = false;
        cardPanel.add("Time", timePlot);
        timeDistPlot = new TimeDistPlot(this);
        timeDistPlot.DEBUG = false;
        cardPanel.add("TimeDist", timeDistPlot);
        validate();
    }

    public void start() {
        tauPCreate.modelFilename = modelFilename;
        tauPCreate.directory = directory;
        try {
            VelocityModel vMod = tauPCreate.loadVMod();
            TauModel tMod = tauPCreate.createTauModel(vMod);
            velocityPlot.plot(vMod, 'P', 'S');
            slownessPlot.plot(tMod.getSlownessModel(), true);
            tauPPlot.plot(tMod, true);
            timePlot.plot(tMod, true);
            distPlot.plot(tMod, true);
            timeDistPlot.plot(tMod, true);
        } catch(NoSuchMatPropException e) {
            throw new RuntimeException("Can't happen", e);
        } catch(VelocityModelException e) {
            throw new RuntimeException("Can't happen", e);
        } catch(IOException e) {
            throw new RuntimeException("Can't happen", e);
        } catch(SlownessModelException e) {
            throw new RuntimeException("Can't happen", e);
        } catch(TauModelException e) {
            throw new RuntimeException("Can't happen", e);
        }
    }

    public void actionPerformed(ActionEvent action) {
        String arg = action.getActionCommand();
        if("Previous".equals(arg)) {
            cards.previous(cardPanel);
        } else if("Next".equals(arg)) {
            cards.next(cardPanel);
        } else if("Quit".equals(arg)) {
            System.exit(0);
        }
    }

    public static void main(String[] args) {
        if(args.length != 1) {
            System.out.println("Usage java TauP_WCreate velocitymodel.tvel [wavetype]");
            System.exit(1);
        }
        MainFrame f = new MainFrame("TauP_WCreate");
        TauP_WCreate tauPWCreate = new TauP_WCreate();
        int j = args[0].lastIndexOf(System.getProperty("file.separator"));
        if(j == -1) {
            tauPWCreate.modelFilename = args[0];
            tauPWCreate.directory = ".";
        } else {
            tauPWCreate.modelFilename = args[0].substring(j + 1);
            tauPWCreate.directory = args[0].substring(0, j);
        }
        tauPWCreate.init();
        f.add("Center", tauPWCreate);
        f.pack();
        f.show();
        tauPWCreate.start();
    }
}

class MainFrame extends Frame implements WindowListener {

    MainFrame(String title) {
        super(title);
        addWindowListener(this);
    }

    public void windowClosed(WindowEvent event) {}

    public void windowDeiconified(WindowEvent event) {}

    public void windowIconified(WindowEvent event) {}

    public void windowActivated(WindowEvent event) {}

    public void windowDeactivated(WindowEvent event) {}

    public void windowOpened(WindowEvent event) {}

    public void windowClosing(WindowEvent event) {
        System.exit(0);
    }
}
