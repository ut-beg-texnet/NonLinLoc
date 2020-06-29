/*
 * This file is part of the Anthony Lomax Java Library.
 *
 * Copyright (C) 2003 Anthony Lomax <anthony@alomax.net>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

/*
 The TauP Toolkit: Flexible Seismic Travel-Time and Raypath Utilities.
 Copyright (C) 1998-2000 University of South Carolina

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

 The current version can be found at
 <A HREF="www.seis.sc.edu">http://www.seis.sc.edu</A>

 Bug reports and comments should be directed to
 H. Philip Crotwell, crotwell@seis.sc.edu or
 Tom Owens, owens@seis.sc.edu

 */
/*
 * TauP_Table_NLL.java
 *
 * Created on 01 December 2003, 11:08
 */
package edu.sc.seis.TauP;

import java.io.*;
import java.util.*;

/**
 *
 * @author Anthony Lomax
 */
public class TauP_Table_NLL extends TauP_Table {

    public static final int NONLINLOC = 1002;
    public static final int NONLINLOC_ANGLES = 1003;
    public static final int EARLYEST = 1004;
    public static final int EARLYEST_ANGLES = 1005;
    protected String nllGridDesc = "181,0.0,180.0,61,0.0,600.0";
    protected int numDistance;
    protected double distance0;
    protected double deltaDistance;
    protected int numDepth;
    protected double depth0;
    protected double deltaDepth;
    protected String stationName = "DEFAULT";
    protected double stationLong = 0.0;
    protected double stationLat = 0.0;
    protected double stationDepth = 0.0;

    /**
     * Creates a new instance of TauP_Table_NLL
     */
    public TauP_Table_NLL() throws TauModelException {
        super();
    }

    public void init() throws IOException {
        super.init();
    }

    public void start() throws TauModelException, TauPException, IOException {
        switch (outputType) {
            case TauP_Table.GENERIC:
                genericTable(getWriter());
                break;
            case TauP_Table.LOCSAT:
                locsatTable(getWriter());
                break;
            case TauP_Table_NLL.NONLINLOC:
                nonLinLocTable();
                break;
            case TauP_Table_NLL.NONLINLOC_ANGLES:
                nonLinLocTable();
                break;
            case TauP_Table_NLL.EARLYEST:
                earlyEstTable();
                break;
            case TauP_Table_NLL.EARLYEST_ANGLES:
                earlyEstTable();
                break;
            default:
                throw new TauPException(
                        "TauP_Table_NLL: undefined state for output type: " + outputType);
        }
    }

    protected void fillTables(float[][] timeArray, short[][][] angleArray, float[][] slownessArray, char materialProperty,
            boolean fillSlownessArray, int numDistanceSlowness) throws TauPException, TauModelException {

        // loop over depth (z)
        for (int depthNum = 0; depthNum < numDepth; depthNum++) {
            double depth = depth0 + deltaDepth * (double) depthNum;
            depthCorrect(depth0 + deltaDepth * (double) depthNum);
            double velocity = 0.0;
            try {
                int layerNumber = tMod.getVelocityModel().layerNumberBelow(depth);
                VelocityLayer velocityLayer = tMod.getVelocityModel().getVelocityLayer(layerNumber);
                velocity = velocityLayer.evaluateAt(depth, materialProperty);
            } catch (Exception e) {
                throw (new TauPException(e.getMessage()));
            }
            // loop over distance (x)
            for (int distNum = 0; distNum < numDistance; distNum++) {
                //
                if (fillSlownessArray && distNum < numDistanceSlowness) {
                    slownessArray[depthNum][distNum] = (float) (1.0 / velocity);
                }
                //
                calculate(distance0 + deltaDistance * (double) distNum);
                //Arrival[] arrivals = getArrivals();
                List<Arrival> arrivals = getArrivals();  // 20130510 AJL - changed for compatibility with TauP-2.1.1
                double time = Double.MAX_VALUE;
                Arrival minArrival = null;
                for (int aNum = 0; aNum < getNumArrivals(); aNum++) {
                    Arrival currArrival = arrivals.get(aNum);
                    if (currArrival.time < time) {
                        time = currArrival.time;
                        minArrival = currArrival;
                    }
                }
                if (outputType == TauP_Table_NLL.NONLINLOC_ANGLES || outputType == TauP_Table_NLL.EARLYEST_ANGLES) {
                    /* 1983__Buland_Chapman__The_Computation_of_Seismic_Travel_Times__BSSA
                     *  rayParam = r * sin(i) / velocity
                     *  thus,
                     *  i = asin(rayParam * velocity / r)
                     */
                    // If there is no arrival, set the angle to -999 - C.Satriano 2011/11/11
                    if (minArrival == null) {
                        if (outputType == TauP_Table_NLL.NONLINLOC_ANGLES) {
                            setTakeOffAngles(angleArray[depthNum][distNum], 1.0, -999, 9);
                        } else {
                            angleArray[depthNum][distNum][0] = (short) -1;
                            angleArray[depthNum][distNum][1] = (short) -999;
                        }
                        continue;
                    }
                    // 20130514 AJL - replaced dip calculation with TauP takeoffAngle
                    double dip = minArrival.getTakeoffAngle();
                    if (dip < 0.0) {
                        dip += 180.0;   // TauP applies "fake neg velocity so angle is neg in case of upgoing"
                    }                    /*
                     double dip = (180.0 / Math.PI) * Math.asin(minArrival.rayParam * velocity / tMod.getSlownessModel().getRadiusOfEarth());
                     // check down- or up-going
                     boolean downgoing = true;
                     int branch = minArrival.phase.tMod.getSourceBranch();
                     while (true) {
                     try {
                     downgoing = ((Boolean) (minArrival.phase.downGoing.get(branch))).booleanValue();
                     break;
                     } catch (Exception e) {
                     if (branch > 0) {
                     branch--;
                     continue;
                     }
                     System.out.println("downgoing error: " + e);
                     break;
                     }
                     }
                     if (!downgoing) {
                     dip = 180.0 - dip;
                     }
                     */

                    if (outputType == TauP_Table_NLL.NONLINLOC_ANGLES) {
                        setTakeOffAngles(angleArray[depthNum][distNum], 1.0, dip, 9);
                    } else {
                        angleArray[depthNum][distNum][0] = (short) 1;
                        angleArray[depthNum][distNum][1] = (short) dip;
                    }
                    //System.out.println("dist " + (distance0 + deltaDistance * (double) distNum) + " d " + depth + " v " + velocity
                    //+ " rp " + minArrival.rayParam + " dip " + dip + " downgoing " + downgoing);
                } else {
                    if (time == Double.MAX_VALUE) {
                        time = -1.0;
                    }
                    timeArray[depthNum][distNum] = (float) time;
                }
            }
        }
    }

    /**
     * generates and save to disk time, angles and model tables for NLL
     *
     */
    protected void nonLinLocTable() throws TauPException, TauModelException, IOException {

        parseGridHeader();

        DataOutputStream bufOut;
        DataOutputStream slownessBufOut = null;

        int numDistanceSlowness = 3;

        char materialProperty = (((PhaseName) phaseNames.get(0)).getName()).charAt(0);

        String timeFileRoot = "";
        String type = "ERROR";
        if (outputType == TauP_Table_NLL.NONLINLOC) {
            timeFileRoot = outFile + "." + phaseNames.get(0) + "." + stationName;
            timeFileRoot += ".time";
            type = "TIME2D";
        } else if (outputType == TauP_Table_NLL.NONLINLOC_ANGLES) {
            timeFileRoot = outFile + "." + phaseNames.get(0) + "." + stationName;
            timeFileRoot += ".angle";
            type = "ANGLE2D";
        }
        // open output file
        if (outFile != null && outFile.length() != 0) {
            bufOut = new DataOutputStream(new BufferedOutputStream(
                    new FileOutputStream(timeFileRoot + ".buf")));
            PrintStream hdrOut = new PrintStream(new FileOutputStream(timeFileRoot + ".hdr"));
            writeGridHeader(hdrOut, numDistance, type);
            hdrOut.close();
            // slowness model file (used by NLLoc LOCMETH OT_STACK, etc.)
            if ((outputType == TauP_Table_NLL.NONLINLOC) && (materialProperty == 'p' || materialProperty == 'P')) {
                String slownessFileRoot = outFile + ".P.mod";
                slownessBufOut = new DataOutputStream(new BufferedOutputStream(
                        new FileOutputStream(slownessFileRoot + ".buf")));
                hdrOut = new PrintStream(new FileOutputStream(slownessFileRoot + ".hdr"));
                writeGridHeader(hdrOut, numDistanceSlowness, "SLOWNESS");
                hdrOut.close();
            } else if ((outputType == TauP_Table_NLL.NONLINLOC) && (materialProperty == 's' || materialProperty == 'S')) {
                String slownessFileRoot = outFile + ".S.mod";
                slownessBufOut = new DataOutputStream(new BufferedOutputStream(
                        new FileOutputStream(slownessFileRoot + ".buf")));
                hdrOut = new PrintStream(new FileOutputStream(slownessFileRoot + ".hdr"));
                writeGridHeader(hdrOut, numDistanceSlowness, "SLOWNESS");
                hdrOut.close();
            }

        } else {
            bufOut = new DataOutputStream(System.out);
            writeGridHeader(System.out, numDistance, type);
        }

        float[][] timeArray = null;
        short[][][] angleArray = null;
        if (outputType == TauP_Table_NLL.NONLINLOC_ANGLES) {
            angleArray = new short[numDepth][numDistance][2];
        } else {
            timeArray = new float[numDepth][numDistance];
        }
        float[][] slownessArray = null;
        if (slownessBufOut != null) {
            slownessArray = new float[numDepth][numDistanceSlowness];
        }

        // fill requires table arrays
        fillTables(timeArray, angleArray, slownessArray, materialProperty, slownessBufOut != null, numDistanceSlowness);

        // write to grid file
        if (outputType == TauP_Table_NLL.NONLINLOC) {
            // loop over distance (y)
            for (int distNum = 0; distNum < numDistance; distNum++) {
                // loop over depth (z)
                for (int depthNum = 0; depthNum < numDepth; depthNum++) {
                    bufOut.writeFloat(timeArray[depthNum][distNum]);
                }
            }
            if (slownessBufOut != null) {
                // loop over distance (y)
                for (int distNum = 0; distNum < numDistanceSlowness; distNum++) {
                    // loop over depth (z)
                    for (int depthNum = 0; depthNum < numDepth; depthNum++) {
                        slownessBufOut.writeFloat(slownessArray[depthNum][distNum]);
                    }
                }
            }
        } else if (outputType == TauP_Table_NLL.NONLINLOC_ANGLES) {
            // loop over distance (y)
            for (int distNum = 0; distNum < numDistance; distNum++) {
                // loop over depth (z)
                for (int depthNum = 0; depthNum < numDepth; depthNum++) {
                    bufOut.writeShort(angleArray[depthNum][distNum][0]);
                    bufOut.writeShort(angleArray[depthNum][distNum][1]);
                    //System.out.println("i0 " + angleArray[depthNum][distNum][0] + " i1 " + angleArray[depthNum][distNum][1]);
                }
            }
        }

        bufOut.close();
        if (slownessBufOut != null) {
            slownessBufOut.close();
        }
    }

    /**
     * generates and save to disk time, angles and model tables for NLL
     *
     */
    protected void earlyEstTable() throws TauPException, TauModelException, IOException {

        parseGridHeader();

        PrintStream tableCheaderOut;
        PrintStream slownessBufOut = null;

        int numDistanceSlowness = 3;

        char materialProperty = (phaseNames.get(0).getName()).charAt(0);

        String tableCheaderFileName = "";
        String type = "ERROR";
        if (outputType == TauP_Table_NLL.EARLYEST) {
            // <outFile>_times_phases.h
            tableCheaderFileName = outFile + "_times_phases.h_part_" + phaseNames.get(0);
            type = "TIME2D";
        } else if (outputType == TauP_Table_NLL.EARLYEST_ANGLES) {
            // <outFile>_times_phases.h
            tableCheaderFileName = outFile + "_toang_phases.h_part_" + phaseNames.get(0);
            type = "ANGLE2D";
        }
        // open output file
        if (outFile != null && outFile.length() != 0) {
            tableCheaderOut = new PrintStream(new FileOutputStream(tableCheaderFileName));
            /*
             // slowness model file (used by NLLoc LOCMETH OT_STACK, etc.)
             if ((outputType == TauP_Table_NLL.NONLINLOC) && (materialProperty == 'p' || materialProperty == 'P')) {
             String slownessFileRoot = outFile + ".P.mod";
             slownessBufOut = new DataOutputStream(new BufferedOutputStream(
             new FileOutputStream(slownessFileRoot + ".buf")));
             hdrOut = new PrintStream(new FileOutputStream(slownessFileRoot + ".hdr"));
             writeGridHeader(hdrOut, numDistanceSlowness, "SLOWNESS");
             hdrOut.close();
             } else if ((outputType == TauP_Table_NLL.NONLINLOC) && (materialProperty == 's' || materialProperty == 'S')) {
             String slownessFileRoot = outFile + ".S.mod";
             slownessBufOut = new DataOutputStream(new BufferedOutputStream(
             new FileOutputStream(slownessFileRoot + ".buf")));
             hdrOut = new PrintStream(new FileOutputStream(slownessFileRoot + ".hdr"));
             writeGridHeader(hdrOut, numDistanceSlowness, "SLOWNESS");
             hdrOut.close();
             }
             */
        } else {
            tableCheaderOut = new PrintStream((System.out));
        }

        float[][] timeArray = null;
        short[][][] angleArray = null;
        if (outputType == TauP_Table_NLL.EARLYEST_ANGLES) {
            angleArray = new short[numDepth][numDistance][2];
        } else {
            timeArray = new float[numDepth][numDistance];
        }
        float[][] slownessArray = null;
        /*
         if (slownessBufOut != null) {
         slownessArray = new float[numDepth][numDistanceSlowness];
         }*/

        // fill requires table arrays
        fillTables(timeArray, angleArray, slownessArray, materialProperty, slownessBufOut != null, numDistanceSlowness);

        // write to grid file
        java.text.DecimalFormat formatter = new java.text.DecimalFormat(".##");
        tableCheaderOut.println("// Created: " + (new Date()) + " by: " + this.getClass().getName());
        if (outputType == TauP_Table_NLL.EARLYEST) {
            tableCheaderOut.println("// " + "times " + phaseNames.get(0));
            tableCheaderOut.println("{");
            // loop over distance (y)
            for (int distNum = 0; distNum < numDistance; distNum++) {
                tableCheaderOut.print("{");
                tableCheaderOut.print(formatter.format(distance0 + deltaDistance * (double) distNum));
                tableCheaderOut.print("," + 1);
                // loop over depth (z)
                for (int depthNum = 0; depthNum < numDepth; depthNum++) {
                    tableCheaderOut.print(",");
                    tableCheaderOut.print(formatter.format(timeArray[depthNum][distNum]));
                }
                tableCheaderOut.println("},");
            }
            tableCheaderOut.println("},");
            /*
             if (slownessBufOut != null) {
             // loop over distance (y)
             for (int distNum = 0; distNum < numDistanceSlowness; distNum++) {
             // loop over depth (z)
             for (int depthNum = 0; depthNum < numDepth; depthNum++) {
             slownessBufOut.writeFloat(slownessArray[depthNum][distNum]);
             }
             }
             }*/
        } else if (outputType == TauP_Table_NLL.EARLYEST_ANGLES) {
            tableCheaderOut.println("// " + "toang " + phaseNames.get(0));
            tableCheaderOut.println("{");
            // loop over distance (y)
            for (int distNum = 0; distNum < numDistance; distNum++) {
                tableCheaderOut.print("{");
                tableCheaderOut.print(formatter.format(distance0 + deltaDistance * (double) distNum));
                tableCheaderOut.print("," + 1);
                // loop over depth (z)
                for (int depthNum = 0; depthNum < numDepth; depthNum++) {
                    tableCheaderOut.print(",");
                    tableCheaderOut.print(angleArray[depthNum][distNum][1]);
                }
                tableCheaderOut.println("},");
            }
            tableCheaderOut.println("},");
        }

        tableCheaderOut.close();
        /*
         if (slownessBufOut != null) {
         slownessBufOut.close();
         }*/

        // write model depth_Vp_Vs_rho header file
        /*
         //ak135 - depth, Vp, Vs, rho from ak135.tvel in iaspei-tau ttimes
         // from 0 to 800km depth only

         #define NUM_TVEL_DEPTH 24

         static double depth_Vp_Vs_rho[NUM_TVEL_DEPTH][4] = {
         { 0.00, 5.80, 3.46, 2.72 },
         ...
         { 809.50, 11.14, 6.24, 4.46 },
         };
         */
        tableCheaderFileName = outFile + "_tvel.h";
        PrintStream modelCheaderOut = new PrintStream(new FileOutputStream(tableCheaderFileName));
        modelCheaderOut.println("// Created: " + (new Date()) + " by: " + this.getClass().getName());
        modelCheaderOut.println("// " + tMod.getVelocityModel().getModelName());
        int nlayers = tMod.getVelocityModel().getNumLayers();
        modelCheaderOut.println("#define NUM_TVEL_DEPTH " + nlayers);
        modelCheaderOut.println("static double depth_Vp_Vs_rho[NUM_TVEL_DEPTH][4] = {");
        for (int n = 0; n < nlayers; n++) {
            VelocityLayer velocityLayer = tMod.getVelocityModel().getVelocityLayer(n);
            modelCheaderOut.println("{ " + velocityLayer.getTopDepth() + ", " + velocityLayer.getTopPVelocity() + ", "
                    + velocityLayer.getTopSVelocity() + ", " + velocityLayer.getTopDensity() + " },");
            // use only top of layer - Early-est assumes constant velocity layers!
            //modelCheaderOut.println("{ " + velocityLayer.getBotDepth() + ", " + velocityLayer.getBotPVelocity() + ", "
            //        + velocityLayer.getBotSVelocity() + ", " + velocityLayer.getBotDensity() + " },");
        }
        modelCheaderOut.println("};");
        modelCheaderOut.close();
    }

    /**
     * * function to set angle values in take-off angles union
     */
    public void setTakeOffAngles(short[] ival, double azim, double dip, int iqual) {

        int iazim = (int) Math.round(10.0 * azim);
        int idip = (int) Math.round(10.0 * dip);

        // remove sign bits
        iazim = (iazim << 1) >>> 1;
        idip = (idip << 1) >>> 1;
        iqual = (iqual << 1) >>> 1;

        ival[0] = (short) iazim;
        ival[1] = (short) (idip * 16 + iqual);

        //System.out.println("iazim " + iazim + " idip " + idip + " iqual " + iqual + " i0 " + ival[0] + " i1 " + ival[1]);
    }

    /**
     * writes NLL Grid Header file contents to hdrOUT
     */
    public void writeGridHeader(PrintStream hdrOut, int numDistanceHdr, String type) {

        // write NLL grid header lines
        hdrOut.print(1);
        hdrOut.print(" ");
        hdrOut.print(numDistanceHdr);
        hdrOut.print(" ");
        hdrOut.print(numDepth);
        hdrOut.print("  ");

        hdrOut.print(0.0);
        hdrOut.print(" ");
        hdrOut.print(distance0);
        hdrOut.print(" ");
        hdrOut.print(depth0);
        hdrOut.print("  ");

        hdrOut.print(0.0);
        hdrOut.print(" ");
        hdrOut.print(deltaDistance);
        hdrOut.print(" ");
        hdrOut.print(deltaDepth);
        hdrOut.print("  ");

        hdrOut.print(type);
        hdrOut.print("  ");

        hdrOut.print("FLOAT");

        hdrOut.println();

        hdrOut.println(stationName + " " + stationLong + " " + stationLat + " " + stationDepth);

    }

    /**
     * parse NLL Grid Header string
     */
    public void parseGridHeader() throws TauPException {

        double distance1;
        double depth1;

        // parse grid description
        StringTokenizer stkzr = new StringTokenizer(nllGridDesc, ",");
        try {
            // "181,0.0,180.0,61,0.0,600.0"
            numDistance = Integer.parseInt(stkzr.nextToken());
            distance0 = Double.parseDouble(stkzr.nextToken());
            distance1 = Double.parseDouble(stkzr.nextToken());
            numDepth = Integer.parseInt(stkzr.nextToken());
            depth0 = Double.parseDouble(stkzr.nextToken());
            depth1 = Double.parseDouble(stkzr.nextToken());
        } catch (Exception e) {
            throw new TauPException(
                    "TauP_Table_NLL: error parsing grid description string: " + nllGridDesc);
        }
        deltaDistance = (distance1 - distance0) / (double) (numDistance - 1);
        deltaDepth = (depth1 - depth0) / (double) (numDepth - 1);
    }

    public void printUsage() {
        printStdUsageHead();
        System.out.println(
                "-ph phase list     -- comma separated phase list\n"
                + "-pf phasefile      -- file containing phases\n\n"
                + "-mod[el] modelname -- use velocity model \"modelname\" for calculations\n"
                + "                      Default is iasp91.\n\n");

        System.out.println(
                "-header filename   -- reads depth and distance spacing data\n"
                + "                      from a LOCSAT style file.");
        System.out.println(
                "-generic           -- outputs a \"generic\" ascii table\n");
        System.out.println(
                "-locsat            -- outputs a \"locsat\" style ascii table\n");
        System.out.println(
                "-nll griddesc      -- outputs a \"NonLinLoc\" 3D Grid Data buffer file of travel-times for all phases\n");
        System.out.println(
                "-nll_angles griddesc      -- outputs a \"NonLinLoc\" 3D Grid Data buffer file of take-off angles for all phases\n");
        System.out.println(
                "-earlyest griddesc      -- outputs an \"Early-Est\" C header file of travel-times for all phases\n");
        System.out.println(
                "-earlyest_angles griddesc      -- outputs an \"Early-Est\" C header file of travel-times for all phases\n");
        printStdUsageTail();
    }

    public String[] parseCmdLineArgs(String[] args) throws IOException {
        int i = 0;
        String[] leftOverArgs;
        int numNoComprendoArgs = 0;
        File tempFile;

        leftOverArgs = super.parseCmdLineArgs(args);
        String[] noComprendoArgs = new String[leftOverArgs.length];

        while (i < leftOverArgs.length) {
            if (leftOverArgs[i].equals("-nll")) {
                outputType = NONLINLOC;
                nllGridDesc = args[i + 1];
                i++;
            } else if (leftOverArgs[i].equals("-nll_angles")) {
                outputType = NONLINLOC_ANGLES;
                nllGridDesc = args[i + 1];
                i++;
            } else if (leftOverArgs[i].equals("-earlyest")) {
                outputType = EARLYEST;
                nllGridDesc = args[i + 1];
                i++;
            } else if (leftOverArgs[i].equals("-earlyest_angles")) {
                outputType = EARLYEST_ANGLES;
                nllGridDesc = args[i + 1];
                i++;
            } else {
                noComprendoArgs[numNoComprendoArgs++] = leftOverArgs[i];
            }
            i++;
        }

        if (numNoComprendoArgs > 0) {
            String[] temp = new String[numNoComprendoArgs];
            System.arraycopy(noComprendoArgs, 0, temp, 0, numNoComprendoArgs);
            return temp;
        } else {
            return new String[0];
        }
    }

    public static void main(String[] args) {
        TauP_Table_NLL me;
        try {
            me = new TauP_Table_NLL();
            String[] noComprendoArgs = me.parseCmdLineArgs(args);
            if (noComprendoArgs.length > 0) {
                for (int i = 0; i < noComprendoArgs.length; i++) {
                    if (noComprendoArgs[i].equals("-help")
                            || noComprendoArgs[i].equals("-version")) {
                        System.exit(0);
                    }
                }
                System.out.println("I don't understand the following arguments, continuing:");
                for (int i = 0; i < noComprendoArgs.length; i++) {
                    System.out.print(noComprendoArgs[i] + " ");
                    if (noComprendoArgs[i].equals("-help")
                            || noComprendoArgs[i].equals("-version")) {
                        System.out.println();
                        System.exit(0);
                    }
                }
                System.out.println();
                noComprendoArgs = null;
            }

            me.init();
            me.start();
        } catch (Exception e) {
            System.err.println("Caught Exception: " + e);
            e.printStackTrace();
        }
    }
}
