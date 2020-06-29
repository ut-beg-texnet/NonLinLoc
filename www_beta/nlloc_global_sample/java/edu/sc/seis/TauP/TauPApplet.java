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
import java.awt.BorderLayout;
import java.awt.Button;
import java.awt.CardLayout;
import java.awt.Choice;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.Label;
import java.awt.Panel;
import java.awt.TextArea;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.IOException;
import java.io.InputStream;
import java.io.InvalidClassException;
import java.io.PrintWriter;
import java.io.StreamTokenizer;
import java.io.StringWriter;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.net.URL;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

/**
 * Simple applet to run TauP tools.
 * 
 * 
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * 
 * 
 * 
 * @author H. Philip Crotwell
 */
public class TauPApplet extends Applet implements ActionListener, ItemListener {

    Choice getModel;

    Choice toolChoice;

    TauP_Time tool;

    TauModel tMod;

    Panel inputPanel;

    TextField modelField, depthField, distanceField, phasesField;

    TextArea textArea;

    PolarPlot plotArea;

    CardLayout outputCards;

    Panel outputPanel;

    Button calculate;

    Button textOrPlot;

    Button clearTextArea;

    String newline;

    String modelName, phases = "";

    double distance, depth;
    
    protected void resetInputPanel() {
        String modelName = getModel.getSelectedItem();
        inputPanel.removeAll();
        // Add Components to the Applet.
        Label l;
        l = new Label("Model Name", Label.CENTER);
        inputPanel.add(l);
        l = new Label("Source Depth", Label.CENTER);
        inputPanel.add(l);
        l = new Label("Distance", Label.CENTER);
        inputPanel.add(l);
        l = new Label("Phase List", Label.CENTER);
        inputPanel.add(l);
        modelField = new TextField(modelName, 5);
        modelField.setEditable(false);
        inputPanel.add(modelField);
        depthField = new TextField("0.0", 5);
        inputPanel.add(depthField);
        distanceField = new TextField("10.0", 5);
        inputPanel.add(distanceField);
        phasesField = new TextField("P,S,PKP,SKS", 5);
        inputPanel.add(phasesField);
        inputPanel.validate();
    }

    public void init() {
        GridBagLayout gridbag = new GridBagLayout();
        setLayout(gridbag);
        GridBagConstraints constraints = new GridBagConstraints();
        constraints.gridwidth = GridBagConstraints.REMAINDER;
        constraints.weightx = 1.0;
        constraints.insets = new Insets(4, 4, 5, 5);
        Panel choicePanel = new Panel();
        toolChoice = new Choice();
        toolChoice.addItem("TauP_Time");
        toolChoice.addItem("TauP_Pierce");
        toolChoice.addItem("TauP_Path");
        toolChoice.addItem("TauP_Curve");
        toolChoice.select(0);
        toolChoice.addItemListener(this);
        choicePanel.add(new Label("Tool"));
        choicePanel.add(toolChoice);
        getModel = new Choice();
        choicePanel.add(new Label("Choose model."));
        choicePanel.add(getModel);
        getModel.addItemListener(this);
        gridbag.setConstraints(choicePanel, constraints);
        add(choicePanel);
        inputPanel = new Panel();
        inputPanel.setLayout(new GridLayout(0, 4));
        constraints.weighty = 0.0;
        gridbag.setConstraints(inputPanel, constraints);
        add(inputPanel);
        resetInputPanel();
        calculate = new Button("Calculate");
        calculate.addActionListener(this);
        clearTextArea = new Button("Clear");
        clearTextArea.addActionListener(this);
        textOrPlot = new Button("Plot");
        textOrPlot.addActionListener(this);
        constraints.weighty = 0.0;
        constraints.gridwidth = 1;
        add(clearTextArea);
        add(textOrPlot);
        constraints.gridwidth = GridBagConstraints.REMAINDER;
        gridbag.setConstraints(calculate, constraints);
        add(calculate);
        outputPanel = new Panel();
        outputCards = new CardLayout();
        outputPanel.setLayout(outputCards);
        textArea = new TextArea();
        textArea.setEditable(false);
        outputPanel.add(textArea, "Text");
        outputCards.first(outputPanel);
        Panel plotPanel = new Panel();
        plotPanel.setLayout(new BorderLayout());
        Panel buttonPanel = new Panel();
        buttonPanel.setLayout(new GridLayout(0, 1));
        Button fullButton = new Button("Full");
        fullButton.addActionListener(this);
        buttonPanel.add(fullButton);
        Button halfButton = new Button("Half");
        halfButton.addActionListener(this);
        buttonPanel.add(halfButton);
        Button quarterButton = new Button("Quarter");
        quarterButton.addActionListener(this);
        buttonPanel.add(quarterButton);
        plotPanel.add(buttonPanel, BorderLayout.WEST);
        plotArea = new PolarPlot(plotPanel, 250);
        plotPanel.add(plotArea);
        outputPanel.add(plotPanel, "Plot");
        constraints.weighty = 1.0;
        constraints.fill = GridBagConstraints.BOTH;
        gridbag.setConstraints(outputPanel, constraints);
        add(outputPanel);
        findTauModel();
        if(getModel.getItemCount() > 0) {
            getModel.select(0);
            getModel.select("prem");
        }
        validate();
        newline = System.getProperty("line.separator");
    }

    public void loadTauModel(String modelName) throws IOException,
            InvalidClassException {
        try {
            InputStream modelStream = getClass().getResourceAsStream("/StdModels/"
                    + modelName + ".taup");
            if(modelStream != null) {
                textArea.append("loading " + modelName + "...");
                tMod = TauModel.readModelFromStream(modelStream);
                this.modelName = modelName;
                modelField.setText(modelName);
                textArea.append("Got it.\n");
            } else {
                textArea.append("Couldn't find model, InputStream is null.\n");
            }
        } catch(ClassNotFoundException exptn) {
            System.out.println("itemStateChanged: caught ClassNotFoundException:"
                    + exptn.getMessage());
        }
    }

    public void findTauModel() {
        try {
            StreamTokenizer manifest;
            InputStream manifestStream;
            if(getParameter("ARCHIVE") != null) {
                URL jarURL = new URL(getCodeBase() + getParameter("ARCHIVE"));
                ZipInputStream jarStream = new ZipInputStream(jarURL.openStream());
                ZipEntry jarEntry = jarStream.getNextEntry();
                while(jarEntry != null) {
                    if(jarEntry.getName().startsWith("StdModels")
                            && jarEntry.getName().endsWith(".taup")) {
                        getModel.addItem(jarEntry.getName()
                                .substring(jarEntry.getName().lastIndexOf('/') + 1,
                                           jarEntry.getName().length() - 5));
                    }
                    textArea.append(jarEntry.getName() + "\n");
                    System.out.println("jar entry name is "
                            + jarEntry.getName());
                    jarEntry = jarStream.getNextEntry();
                }
                jarStream.close();
                textArea.append("Got models from jar.\n");
            } else {
                // no manifest, so not a jar? try getting StdModels directory
                manifestStream = getClass().getResourceAsStream("/StdModels/");
                if(manifestStream != null) {
                    manifest = new StreamTokenizer(manifestStream);
                    manifest.ordinaryChars('.', '.');
                    manifest.wordChars('.', '.');
                    manifest.ordinaryChars('0', '9');
                    manifest.wordChars('0', '9');
                    while(manifest.nextToken() != manifest.TT_EOF) {
                        if(manifest.ttype == StreamTokenizer.TT_WORD
                                && manifest.sval.endsWith(".taup")) {
                            getModel.addItem(manifest.sval.substring(0,
                                                                     manifest.sval.length() - 5));
                        }
                    }
                    textArea.append("Got models from local directory.\n");
                } else {
                    textArea.append("Error: Couldn't find models.\n");
                }
                manifestStream.close();
            }
        } catch(IOException exptn) {
            System.out.println("findTauModel: caught IOException:"
                    + exptn.getMessage() + "\ngetCodeBase=" + getCodeBase()
                    + "\ngetParameter(archive)=" + getParameter("archive"));
        }
    }

    public void loadTool(String toolName, TauModel newTauModel) {
        try {
            textArea.append("tool change to " + toolChoice.getSelectedItem());
            // create an instance of the tool class
            Class toolClass = Class.forName("edu.sc.seis.TauP." + toolName);
            if(tool == null || !tool.getClass().equals(toolClass)) {
                // invoke the constructor with single arguement (TauModel)
                Class[] argClasses = {Class.forName("edu.sc.seis.TauP.TauModel")};
                Object[] argObjects = {newTauModel};
                Constructor toolConstructor = toolClass.getConstructor(argClasses);
                TauP_Time tempTool = (TauP_Time)toolConstructor.newInstance(argObjects);
                tool = tempTool;
                tMod = newTauModel;
                textArea.append(" successful.");
            } else {
                textArea.append(" already done.");
            }
        } catch(ClassNotFoundException ex) {
            textArea.append(" failed. ClassNotFoundException:\n"
                    + ex.getMessage());
            return;
        } catch(InstantiationException ex) {
            textArea.append(" failed. InstantiationException:\n"
                    + ex.getMessage());
            return;
        } catch(IllegalAccessException ex) {
            textArea.append(" failed. IllegalAccessException:\n"
                    + ex.getMessage());
            return;
        } catch(InvocationTargetException ex) {
            textArea.append(" failed. InvocationTargetException:\n"
                    + ex.getTargetException().getMessage());
            ex.printStackTrace();
            ex.getTargetException().printStackTrace();
            return;
        } catch(NoSuchMethodException ex) {
            textArea.append(" failed. NoSuchMethodException:\n"
                    + ex.getMessage());
            return;
        } finally {
            textArea.append("\n");
        }
    }

    public void makePlotActive() {
        outputCards.show(outputPanel, "Plot");
        textOrPlot.setLabel("Text");
    }

    public void makeTextActive() {
        outputCards.show(outputPanel, "Text");
        textOrPlot.setLabel("Plot");
    }

    public void itemStateChanged(ItemEvent event) {
        try {
            if(event.getItemSelectable() == getModel) {
                loadTauModel(getModel.getSelectedItem());
                if(tool == null) {
                    loadTool(toolChoice.getSelectedItem(), tMod);
                } else {
                    tool.setTauModel(tMod);
                }
                phases = phasesField.getText();
                tool.parsePhaseList(phases);
                tool.depthCorrect(depth);
            } else if(event.getItemSelectable() == toolChoice) {
                if(tMod == null) {
                    loadTauModel(getModel.getSelectedItem());
                }
                loadTool(toolChoice.getSelectedItem(), tMod);
                phases = phasesField.getText();
                tool.parsePhaseList(phases);
                tool.depthCorrect(depth);
            }
        } catch(TauModelException e) {
            textArea.append("\n\ndepthCorrection failed, TauModelException: "
                    + e.getMessage());
        } catch(IOException e) {
            textArea.append("\n\ndepthCorrection failed, IOException: "
                    + e.getMessage());
        }
    }

    public void actionPerformed(ActionEvent e) {
        String target = e.getActionCommand();
        if(target == "Calculate") {
            try {
                // make sure we have loaded the model
                if(tMod == null) {
                    loadTauModel(getModel.getSelectedItem());
                    if(tool != null) {
                        tool.setTauModel(tMod);
                        phases = phasesField.getText();
                        tool.parsePhaseList(phases);
                        tool.recalcPhases();
                    }
                }
                if(tool == null) {
                    loadTool(toolChoice.getSelectedItem(), tMod);
                    phases = phasesField.getText();
                    tool.parsePhaseList(phases);
                    tool.depthCorrect(depth);
                }
                // update the phase list if it has changed
                if(phases != phasesField.getText()) {
                    phases = phasesField.getText();
                    tool.clearPhaseNames();
                    tool.parsePhaseList(phases);
                    tool.recalcPhases();
                }
                // update depth. Note this only recomputes the TauModel if the
                // depth has changed
                depth = Double.valueOf(depthField.getText()).doubleValue();
                tool.depthCorrect(depth);
                // calculate the times for each phase at the given distance
                distance = Double.valueOf(distanceField.getText())
                        .doubleValue();
                tool.calculate(distance);
                if(textOrPlot.getLabel().equals("Plot")) {
                    // print it out
                    StringWriter result = new StringWriter();
                    tool.printResult(new PrintWriter(result));
                    textArea.append(result.toString());
                    textArea.append("Done!\n");
                } else {
                    if(tool.getClass()
                            .getName()
                            .equals("edu.sc.seis.TauP.TauP_Time")) {
                        // plotting doesn't make sense for times, so just print
                        outputCards.show(outputPanel, "Text");
                        textOrPlot.setLabel("Plot");
                        validate();
                        StringWriter result = new StringWriter();
                        tool.printResult(new PrintWriter(result));
                        textArea.append(result.toString());
                        textArea.append("Done!\n");
                    } else if(tool.getClass()
                            .getName()
                            .equals("edu.sc.seis.TauP.TauP_Pierce")) {
                        double[] disconDepths = tool.getDisconDepths();
                        double[] disconRadius = new double[disconDepths.length];
                        for(int i = 0; i < disconDepths.length; i++) {
                            disconRadius[i] = tool.getTauModel()
                                    .getRadiusOfEarth()
                                    - disconDepths[i];
                        }
                        plotArea.setCircles(disconRadius);
                        List<Arrival> arrivals = tool.getArrivals();
                        plotArea.clearSegments();
                        for (Arrival arrival : arrivals) {
                            if(arrival.pierce != null) {
                                plotArea.appendSegment(arrival.pierce);
                            }
                        }
                        plotArea.repaint();
                    } else if(tool.getClass()
                            .getName()
                            .equals("edu.sc.seis.TauP.TauP_Path")) {
                        double[] disconDepths = tool.getDisconDepths();
                        double[] disconRadius = new double[disconDepths.length];
                        for(int i = 0; i < disconDepths.length; i++) {
                            disconRadius[i] = tool.getTauModel()
                                    .getRadiusOfEarth()
                                    - disconDepths[i];
                        }
                        plotArea.setCircles(disconRadius);
                        List<Arrival> arrivals = tool.getArrivals();
                        plotArea.clearSegments();
                        for (Arrival arrival : arrivals) {
                            if(arrival.path != null) {
                                plotArea.appendSegment(arrival.path);
                            }
                        }
                        plotArea.repaint();
                    } else if(tool.getClass()
                            .getName()
                            .equals("edu.sc.seis.TauP.TauP_Curve")) {}
                }
            } catch(TauModelException exptn) {
                System.out.println("actionPerformed: caught TauModelException:"
                        + exptn.getMessage());
                exptn.printStackTrace();
            } catch(IOException exptn) {
                System.out.println("actionPerformed: caught IOException:"
                        + exptn.getMessage());
                exptn.printStackTrace();
            }
        } else if(target == "depthField") {
            try {
                // update depth. Note this only recomputes the TauModel if the
                // depth has changed
                depth = Double.valueOf(depthField.getText()).doubleValue();
                tool.depthCorrect(depth);
            } catch(TauModelException exptn) {
                System.out.println("actionPerformed: caught TauModelException:"
                        + exptn.getMessage());
                exptn.printStackTrace();
            }
        } else if(target == "Clear") {
            textArea.setText("");
        } else if(target == "Plot") {
            makePlotActive();
        } else if(target == "Text") {
            makeTextActive();
        } else if(target == "Full") {
            plotArea.setDisplayMode(PolarPlot.FULL);
        } else if(target == "Half") {
            plotArea.setDisplayMode(PolarPlot.HALF);
        } else if(target == "Quarter") {
            plotArea.setDisplayMode(PolarPlot.QUARTER);
        }
    }
}
