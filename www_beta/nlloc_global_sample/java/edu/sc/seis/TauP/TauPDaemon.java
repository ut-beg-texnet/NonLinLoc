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

import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.InvalidClassException;
import java.io.OptionalDataException;
import java.io.StreamCorruptedException;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.net.ServerSocket;
import java.net.Socket;
import java.util.StringTokenizer;

/**
 * Daemon for travel time calculations. Listens to port 6371 by default and
 * returns travel time calculations in response to client requests.
 * 
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * 
 * 
 * 
 * @author H. Philip Crotwell
 * 
 */
public class TauPDaemon extends Thread {

    Socket theConnection;

    protected static TauModel[] tauModels = new TauModel[10];

    public TauPDaemon(Socket s) {
        theConnection = s;
    }

    /**
     * loads a TauModel, after checking a cache to see if it is already loaded.
     * 
     * @throws FileNotFoundException
     *             if file not found
     * @InvalidClassException if class invalid
     * @throws OptionalDataException
     *             if...
     * @throws TauModelException 
     * @StreamCorruptedException if...
     * @ClassNotFoundException if...
     * @IOException if...
     */
    static TauModel getTauModel(String modelName) throws FileNotFoundException,
            InvalidClassException, OptionalDataException,
            StreamCorruptedException, ClassNotFoundException, IOException, TauModelException {
        for(int i = 0; i < tauModels.length; i++) {
            if(tauModels[i] != null
                    && tauModels[i].getModelName().equals(modelName)) {
                TauModel tMod = tauModels[i];
                for(int j = i; j > 0; j--) {
                    tauModels[j] = tauModels[j - 1];
                }
                tauModels[0] = tMod;
                return tMod;
            }
        }
        // not already loaded, so see if we can find it
        TauModel tMod = TauModelLoader.load(modelName,
                                            System.getProperty("taup.model.path"));
        // add it to front, so older models fall off the end
        for(int j = tauModels.length - 1; j > 0; j--) {
            tauModels[j] = tauModels[j - 1];
        }
        tauModels[0] = tMod;
        return tMod;
    }

    public static void main(String[] args) {
        int thePort;
        ServerSocket server;
        // set port to listen to
        try {
            thePort = Integer.parseInt(args[0]);
            if(thePort < 0 || thePort > 65535)
                thePort = 6371;
        } catch(Exception e) {
            thePort = 6371;
        }
        // load some more common earth models
        String[] modelnames = {"ak135", "prem", "iasp91"};
        TauModel tempTMod;
        for(int i = 0; i < modelnames.length; i++) {
            try {
                tempTMod = getTauModel(modelnames[i]);
            } catch(Exception e) {
                System.err.println("Couldn't load tau model. " + e.getMessage());
            }
        }
        try {
            server = new ServerSocket(thePort);
            System.out.println("Accepting connections on port "
                    + server.getLocalPort());
            while(true) {
                TauPDaemon taupd = new TauPDaemon(server.accept());
                taupd.start();
            }
        } catch(IOException e) {
            System.err.println("Server aborted prematurely. " + e.getMessage());
        }
    }

    public void run() {
        try {
            DataOutputStream out = new DataOutputStream(theConnection.getOutputStream());
            BufferedReader in = new BufferedReader(new InputStreamReader(theConnection.getInputStream()));
            out.writeBytes("Welcome to TauP Travel Time Daemon, version 0.92\n");
            String modelName = "iasp91", phaseString = "P,S,PKP,SKS,PKIKP,SKIKS";
            double distance, depth = 0.0;
            TauModel tMod = null;
            String toolString = "TauP_Time";
            TauP_Time tool = null;
            String line, word;
            StringTokenizer tokenIn;
            while((line = in.readLine()) != null) {
                tokenIn = new StringTokenizer(line);
                word = tokenIn.nextToken();
                if(word.equals("MODEL")) {
                    try {
                        String tempModelName = tokenIn.nextToken();
                        tMod = getTauModel(tempModelName);
                        modelName = tempModelName;
                    } catch(ClassNotFoundException e) {
                        System.err.println("Couldn't load tau model. "
                                + e.getMessage());
                    }
                } else if(word.equals("TOOL")) {
                    toolString = tokenIn.nextToken();
                } else if(word.equals("PHASES")) {
                    phaseString = tokenIn.nextToken();
                } else if(word.equals("DEPTH")) {
                    depth = Double.valueOf(tokenIn.nextToken()).doubleValue();
                } else if(word.equals("DISTANCE")) {
                    distance = Double.valueOf(tokenIn.nextToken())
                            .doubleValue();
                    if(tMod == null) {
                        try {
                            tMod = getTauModel(modelName);
                        } catch(ClassNotFoundException e) {
                            System.err.println("Couldn't load tau model. "
                                    + e.getMessage());
                        }
                    }
                    if(tool == null
                            || !tool.getClass().getName().endsWith(toolString)) {
                        tool = loadTool(toolString, tMod);
                        tool.depthCorrect(depth);
                    }
                    if(!tMod.equals(tool.getTauModel())) {
                        tool.setTauModel(tMod);
                        System.out.println("Changing tmods");
                    }
                    if(phaseString != tool.getPhaseNameString()) {
                        tool.parsePhaseList(phaseString);
                        tool.depthCorrect(depth);
                        phaseString = tool.getPhaseNameString();
                    }
                    tool.calculate(distance);
                    tool.printResult(tool.getWriter());
                }
            }
        } catch(IOException e) {
            System.err.println(e.getMessage());
        } catch(TauModelException e) {
            System.err.println(e.getMessage());
        }
    }

    /**
     * attempts to create a new tool based on the toolName.
     * 
     * @returns a subclass of TauP_Time, or null if unsucessful.
     */
    TauP_Time loadTool(String toolName, TauModel tMod) {
        try {
            // create an instance of the tool class
            Class toolClass = Class.forName("edu.sc.seis.TauP." + toolName);
            // invoke the constructor with single arguement (TauModel)
            Class[] argClasses = {Class.forName("edu.sc.seis.TauP.TauModel")};
            Object[] argObjects = {tMod};
            Constructor toolConstructor = toolClass.getConstructor(argClasses);
            TauP_Time tool = (TauP_Time)toolConstructor.newInstance(argObjects);
            return tool;
        } catch(ClassNotFoundException ex) {} catch(InstantiationException ex) {} catch(IllegalAccessException ex) {} catch(InvocationTargetException ex) {} catch(NoSuchMethodException ex) {}
        return null;
    }
}
