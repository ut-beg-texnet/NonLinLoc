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
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StreamTokenizer;
import java.util.List;

/**
 * Creates a table of travel times for a phase. Only uses the first arrival at
 * any distance.
 * 
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * 
 * 
 * 
 * @author H. Philip Crotwell
 * 
 */
public class TauP_Table extends TauP_Time {

    public static final int GENERIC = 0;

    public static final int LOCSAT = 1;

    protected int outputType = TauP_Table.GENERIC;

    protected String headerFile;

    protected double[] depths = {0.00f,
                                 1.00f,
                                 2.00f,
                                 3.00f,
                                 4.00f,
                                 5.00f,
                                 10.00f,
                                 15.00f,
                                 19.00f,
                                 21.00f,
                                 25.00f,
                                 30.00f,
                                 33.00f,
                                 35.00f,
                                 40.00f,
                                 45.00f,
                                 49.00f,
                                 51.00f,
                                 55.00f,
                                 60.00f,
                                 70.00f,
                                 80.00f,
                                 90.00f,
                                 100.00f,
                                 120.00f,
                                 140.00f,
                                 160.00f,
                                 180.00f,
                                 200.00f,
                                 220.00f,
                                 240.00f,
                                 260.00f,
                                 280.00f,
                                 300.00f,
                                 350.00f,
                                 400.00f,
                                 450.00f,
                                 500.00f,
                                 550.00f,
                                 600.00f};

    protected double[] distances = {0.00f,
                                    0.10f,
                                    0.20f,
                                    0.30f,
                                    0.40f,
                                    0.50f,
                                    0.60f,
                                    0.70f,
                                    0.80f,
                                    0.90f,
                                    1.00f,
                                    1.10f,
                                    1.20f,
                                    1.30f,
                                    1.40f,
                                    1.50f,
                                    1.60f,
                                    1.70f,
                                    1.80f,
                                    1.90f,
                                    2.00f,
                                    2.10f,
                                    2.20f,
                                    2.30f,
                                    2.40f,
                                    2.50f,
                                    2.60f,
                                    2.70f,
                                    2.80f,
                                    2.90f,
                                    3.00f,
                                    3.10f,
                                    3.20f,
                                    3.30f,
                                    3.40f,
                                    3.50f,
                                    3.60f,
                                    3.70f,
                                    3.80f,
                                    3.90f,
                                    4.00f,
                                    4.10f,
                                    4.20f,
                                    4.30f,
                                    4.40f,
                                    4.50f,
                                    4.60f,
                                    4.70f,
                                    4.80f,
                                    4.90f,
                                    5.00f,
                                    5.10f,
                                    5.20f,
                                    5.30f,
                                    5.40f,
                                    5.50f,
                                    5.60f,
                                    5.70f,
                                    5.80f,
                                    5.90f,
                                    6.00f,
                                    6.10f,
                                    6.20f,
                                    6.30f,
                                    6.40f,
                                    6.50f,
                                    6.60f,
                                    6.70f,
                                    6.80f,
                                    6.90f,
                                    7.00f,
                                    7.10f,
                                    7.20f,
                                    7.30f,
                                    7.40f,
                                    7.50f,
                                    7.60f,
                                    7.70f,
                                    7.80f,
                                    7.90f,
                                    8.00f,
                                    8.10f,
                                    8.20f,
                                    8.30f,
                                    8.40f,
                                    8.50f,
                                    8.60f,
                                    8.70f,
                                    8.80f,
                                    8.90f,
                                    9.00f,
                                    9.10f,
                                    9.20f,
                                    9.30f,
                                    9.40f,
                                    9.50f,
                                    9.60f,
                                    9.70f,
                                    9.80f,
                                    9.90f,
                                    10.00f,
                                    11.00f,
                                    12.00f,
                                    13.00f,
                                    14.00f,
                                    15.00f,
                                    16.00f,
                                    17.00f,
                                    18.00f,
                                    19.00f,
                                    20.00f,
                                    21.00f,
                                    22.00f,
                                    23.00f,
                                    24.00f,
                                    25.00f,
                                    26.00f,
                                    27.00f,
                                    28.00f,
                                    29.00f,
                                    30.00f,
                                    31.00f,
                                    32.00f,
                                    33.00f,
                                    34.00f,
                                    35.00f,
                                    36.00f,
                                    37.00f,
                                    38.00f,
                                    39.00f,
                                    40.00f,
                                    41.00f,
                                    42.00f,
                                    43.00f,
                                    44.00f,
                                    45.00f,
                                    46.00f,
                                    47.00f,
                                    48.00f,
                                    49.00f,
                                    50.00f,
                                    51.00f,
                                    52.00f,
                                    53.00f,
                                    54.00f,
                                    55.00f,
                                    56.00f,
                                    57.00f,
                                    58.00f,
                                    59.00f,
                                    60.00f,
                                    61.00f,
                                    62.00f,
                                    63.00f,
                                    64.00f,
                                    65.00f,
                                    66.00f,
                                    67.00f,
                                    68.00f,
                                    69.00f,
                                    70.00f,
                                    71.00f,
                                    72.00f,
                                    73.00f,
                                    74.00f,
                                    75.00f,
                                    76.00f,
                                    77.00f,
                                    78.00f,
                                    79.00f,
                                    80.00f,
                                    81.00f,
                                    82.00f,
                                    83.00f,
                                    84.00f,
                                    85.00f,
                                    86.00f,
                                    87.00f,
                                    88.00f,
                                    89.00f,
                                    90.00f,
                                    91.00f,
                                    92.00f,
                                    93.00f,
                                    94.00f,
                                    95.00f,
                                    96.00f,
                                    97.00f,
                                    98.00f,
                                    99.00f,
                                    100.00f,
                                    105.00f,
                                    110.00f,
                                    115.00f,
                                    120.00f,
                                    125.00f,
                                    130.00f,
                                    135.00f,
                                    140.00f,
                                    145.00f,
                                    150.00f,
                                    155.00f,
                                    160.00f,
                                    165.00f,
                                    170.00f,
                                    175.00f,
                                    180.00f};

    public TauP_Table() throws TauModelException {
        super();
    }

    public void init() throws IOException {
        super.init();
        if(headerFile != null) {
            try {
                StreamTokenizer head = new StreamTokenizer(new BufferedReader(new FileReader(headerFile)));
                head.commentChar('#');
                head.nextToken();
                if(head.ttype == StreamTokenizer.TT_WORD
                        && head.sval.equals("n")) {
                    head.nextToken();
                } else {
                    Alert.warning("First character of header file is not 'n'",
                                  "'" + head.ttype + "'  " + head.sval);
                }
                if(head.ttype != StreamTokenizer.TT_NUMBER) {
                    if(head.ttype == StreamTokenizer.TT_WORD) {
                        Alert.error("Expected a number of depth samples, but got ",
                                    head.sval);
                    } else {
                        Alert.error("Expected a number of depth samples, but got ",
                                    "'" + head.ttype + "'");
                    }
                    System.exit(1);
                }
                if(head.nval != Math.rint(head.nval)) {
                    Alert.error("Expected a number of depth samples, but got ",
                                head.sval);
                }
                depths = new double[(int)head.nval];
                for(int i = 0; i < depths.length; i++) {
                    head.nextToken();
                    if(head.ttype != StreamTokenizer.TT_NUMBER) {
                        if(head.ttype == StreamTokenizer.TT_WORD) {
                            Alert.error("Expected a number of depth samples, but got ",
                                        head.sval);
                        } else {
                            Alert.error("Expected a number of depth samples, but got ",
                                        "'" + head.ttype + "'");
                        }
                        System.exit(1);
                    }
                    depths[i] = head.nval;
                }
                head.nextToken();
                if(head.ttype != StreamTokenizer.TT_NUMBER) {
                    if(head.ttype == StreamTokenizer.TT_WORD) {
                        Alert.error("Expected a number of distance samples, but got ",
                                    head.sval);
                    } else {
                        Alert.error("Expected a number of distance samples, but got ",
                                    "'" + head.ttype + "'");
                    }
                    System.exit(1);
                }
                if(head.nval != Math.rint(head.nval)) {
                    Alert.error("Expected a number of distance samples, but got ",
                                head.sval);
                }
                distances = new double[(int)head.nval];
                for(int i = 0; i < distances.length; i++) {
                    head.nextToken();
                    if(head.ttype != StreamTokenizer.TT_NUMBER) {
                        if(head.ttype == StreamTokenizer.TT_WORD) {
                            Alert.error("Expected a distance sample, but got ",
                                        head.sval);
                        } else {
                            Alert.error("Expected a distance sample, but got ",
                                        "'" + head.ttype + "'");
                        }
                        System.exit(1);
                    }
                    distances[i] = head.nval;
                }
            } catch(FileNotFoundException e) {
                Alert.error("Couldn't find file.", e.getMessage());
                System.exit(1);
            } catch(IOException e) {
                Alert.error("Caught IOException.", e.getMessage());
                System.exit(1);
            }
        }
    }

    public void start() throws TauModelException, TauPException, IOException {
        switch(outputType){
            case TauP_Table.GENERIC:
                genericTable(getWriter());
                break;
            case TauP_Table.LOCSAT:
                locsatTable(getWriter());
                break;
            default:
                throw new TauPException("TauP_Table: undefined state for output type: "
                        + outputType);
        }
    }

    protected void genericTable(PrintWriter out) throws TauModelException,
            IOException {
        for(int depthNum = 0; depthNum < depths.length; depthNum++) {
            depthCorrect(depths[depthNum]);
            for(int distNum = 0; distNum < distances.length; distNum++) {
                calculate(distances[distNum]);
                List<Arrival> arrivals = getArrivals();
                for (Arrival currArrival : arrivals) {
                    double moduloDist = currArrival.getModuloDistDeg();
                    out.print(modelName + " "
                                   + outForms.formatDistance(moduloDist) + " "
                                   + outForms.formatDepth(depth) + " ");
                    out.print(currArrival.getName());
                    out.print("  "
                                   + outForms.formatTime(currArrival.getTime())
                                   + "  ");
                    out.print(outForms.formatRayParam(Math.PI / 180.0
                                                           * currArrival.getRayParam())
                                                           + "   ");
                    out.print(outForms.formatDistance(currArrival.getDistDeg()));
                    out.print("  " + currArrival.getPuristName());
                    out.println();
                }
            }
        }
        out.flush();
    }

    protected void locsatTable(PrintWriter out) throws TauModelException,
            IOException {
        Format float15_4 = new Format("%15.4f");
        Format float7_2 = new Format("%7.2f");
        Format decimal7 = new Format("%-7d");
        double maxDiff = Double.valueOf(toolProps.getProperty("taup.table.locsat.maxdiff",
                                                              "105.0"))
                .doubleValue();
        out.print("n # " + getPhaseNameString()
                + " travel-time tables for " + modelName
                + " structure. (From TauP_Table)\n");
        out.print(decimal7.form(depths.length)
                + "# number of depth samples\n");
        for(int depthNum = 0; depthNum < depths.length; depthNum++) {
            out.print(float7_2.form(depths[depthNum]));
            if(depthNum % 10 == 9) {
                out.println();
            }
        }
        if((depths.length - 1) % 10 != 9) {
            out.println();
        }
        out.println(decimal7.form(distances.length)
                + "# number of distances");
        for(int distNum = 0; distNum < distances.length; distNum++) {
            out.print(float7_2.form(distances[distNum]));
            if(distNum % 10 == 9) {
                out.println();
            }
        }
        if((distances.length - 1) % 10 != 9) {
            out.println();
        }
        for(int depthNum = 0; depthNum < depths.length; depthNum++) {
            depthCorrect(depths[depthNum]);
            out.println("#  Travel time for z =    " + depths[depthNum]);
            for(int distNum = 0; distNum < distances.length; distNum++) {
                calculate(distances[distNum]);
                List<Arrival>arrivals = getArrivals();
                String outString = float15_4.form(-1.0) + "    none\n";
                for (Arrival arrival : arrivals) {
                    if(distances[distNum] > maxDiff
                            && (arrival.getName().endsWith("diff"))) {
                        continue;
                    } else {
                        outString = float15_4.form(arrival.getTime())+ "    " + arrival.getName() + "\n";
                        break;
                    }
                }
                out.print(outString);
            }
        }
        out.flush();
    }

    public void printUsage() {
        printStdUsageHead();
        System.out.println("-ph phase list     -- comma separated phase list\n"
                + "-pf phasefile      -- file containing phases\n\n"
                + "-mod[el] modelname -- use velocity model \"modelname\" for calculations\n"
                + "                      Default is iasp91.\n\n");
        System.out.println("-header filename   -- reads depth and distance spacing data\n"
                + "                      from a LOCSAT style file.");
        System.out.println("-generic           -- outputs a \"generic\" ascii table\n");
        System.out.println("-locsat            -- outputs a \"locsat\" style ascii table\n");
        printStdUsageTail();
    }

    public String[] parseCmdLineArgs(String[] args) throws IOException {
        int i = 0;
        String[] leftOverArgs;
        int numNoComprendoArgs = 0;
        leftOverArgs = super.parseCmdLineArgs(args);
        String[] noComprendoArgs = new String[leftOverArgs.length];
        while(i < leftOverArgs.length) {
            if(dashEquals("header", leftOverArgs[i])
                    && i < leftOverArgs.length - 1) {
                headerFile = leftOverArgs[i + 1];
                i++;
            } else if(dashEquals("locsat", leftOverArgs[i])) {
                outputType = LOCSAT;
            } else if(dashEquals("generic", leftOverArgs[i])) {
                outputType = GENERIC;
            } else if(dashEquals("help", leftOverArgs[i])) {
                noComprendoArgs[numNoComprendoArgs++] = leftOverArgs[i];
            } else {
                noComprendoArgs[numNoComprendoArgs++] = leftOverArgs[i];
            }
            i++;
        }
        if(numNoComprendoArgs > 0) {
            String[] temp = new String[numNoComprendoArgs];
            System.arraycopy(noComprendoArgs, 0, temp, 0, numNoComprendoArgs);
            return temp;
        } else {
            return new String[0];
        }
    }

    public static void main(String[] args) {
        TauP_Table me;
        try {
            me = new TauP_Table();
            String[] noComprendoArgs = me.parseCmdLineArgs(args);
            printNoComprendoArgs(noComprendoArgs);
            me.init();
            me.start();
        } catch(TauModelException e) {
            System.err.println("Caught TauModelException: " + e.getMessage());
            e.printStackTrace();
        } catch(TauPException e) {
            System.err.println("Caught TauModelException: " + e.getMessage());
            e.printStackTrace();
        } catch(IOException e) {
            System.err.println("Caught IOException: " + e.getMessage());
            e.printStackTrace();
        }
    }
}
