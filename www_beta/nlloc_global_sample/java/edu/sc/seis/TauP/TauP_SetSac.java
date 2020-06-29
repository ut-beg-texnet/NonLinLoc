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

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OptionalDataException;
import java.io.StreamCorruptedException;
import java.util.ArrayList;
import java.util.List;

import edu.sc.seis.seisFile.sac.SacConstants;
import edu.sc.seis.seisFile.sac.SacHeader;
import edu.sc.seis.seisFile.sac.SacTimeSeries;

/**
 * Calculate times for phases and set sac headers based on gcarc or dist or
 * station lat and lon and event lat and lon.
 * 
 * Note that triplicated phases will cause problems, as there is only one spot
 * to put a time. An improved method would allow a phase to have several t#'s
 * associated with it, so that all arrivals could be marked. Currently however,
 * only the first arrival for a phase name is used.
 * 
 * Warning: I assume the evdp header has depth in meters unless the -evdpkm flag
 * is set, in which case I assume kilometers. This may be a problem for users
 * that improperly use kilometers for the depth units. Due to much abuse of the
 * sac depth header units I output a warning message if the depth appears to be
 * in kilometers, ie it is < 1000. This can be safely ignored if the event
 * really is less than 1000 meters deep.
 * 
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * 
 * 
 * 
 * @author H. Philip Crotwell
 * 
 */
public class TauP_SetSac extends TauP_Time {

    protected List<String> sacFileNames = new ArrayList<String>();

    protected boolean evdpkm = false;

    public static final int A_HEADER = 10;

    public boolean getEvdpkm() {
        return evdpkm;
    }

    public void setEvdpkm(boolean evdpkm) {
        this.evdpkm = evdpkm;
    }

    public void setSacFileNames(String[] sacFileNames) {
        this.sacFileNames = new ArrayList<String>();
        for(int i = 0; i < sacFileNames.length; i++) {
            this.sacFileNames.add(sacFileNames[i]);
        }
    }

    protected TauP_SetSac() {
        super();
    }

    public TauP_SetSac(TauModel tMod) throws TauModelException {
        super(tMod);
    }

    public TauP_SetSac(String modelName) throws TauModelException {
        super(modelName);
    }

    protected void setSacVarNums() {
        boolean[] headersUsed = new boolean[10];
        for(int i = 0; i < headersUsed.length; i++) {
            headersUsed[i] = false;
        }
        for(int i = 0; i < phaseNames.size() && i < 10; i++) {
            if(((PhaseName)phaseNames.get(i)).sacTNum != -1) {
                headersUsed[((PhaseName)phaseNames.get(i)).sacTNum] = true;
            }
        }
        int j;
        for(int i = 0; i < phaseNames.size(); i++) {
            if(((PhaseName)phaseNames.get(i)).sacTNum == -1) {
                // find a j that hasn't been used
                for(j = 0; j < headersUsed.length && headersUsed[j]; j++) {}
                if(j < 10) {
                    ((PhaseName)phaseNames.get(i)).sacTNum = j;
                    headersUsed[j] = true;
                }
            }
        }
    }

    public void calculate(double degrees) {
        recalcPhases();
        calcTime(degrees);
    }

    public void init() throws IOException {
        super.init();
        setSacVarNums();
    }

    public void start() throws IOException, TauModelException {
        for (String filename : sacFileNames) {
            if(verbose) {
                System.out.println(filename);
            }
            processSacFile(new File(filename));
        }
    }
    
    public void processSacFile(File f) throws FileNotFoundException, IOException, TauModelException {
        
        if (f.isDirectory()) {
            File[] subfiles = f.listFiles();
            for (int j = 0; j < subfiles.length; j++) {
                if (subfiles[j].getName().startsWith(".")) {
                    continue;
                }
                processSacFile(subfiles[j]);
            }
            return;
        }
        // regular file, hopefully
        SacTimeSeries sacFile = new SacTimeSeries(f);
        SacHeader header = sacFile.getHeader();
        if(SacConstants.isUndef(header.getEvdp())) {
            System.out.println("Depth not set in "
                    + f.getName() + ", skipping");
            return;
        }
        if(SacConstants.isUndef(header.getO())) {
            System.out.println("O marker not set in "
                    + f + ", skipping");
            return;
        }
        double deg;
        if(! SacConstants.isUndef(header.getGcarc())) {
            if(verbose) {
                System.out.println("Using gcarc: " + header.getGcarc());
            }
            deg = header.getGcarc();
        } else if(! SacConstants.isUndef(header.getDist())) {
            if(verbose) {
                System.out.println("Using dist: " + header.getDist());
            }
            deg = header.getDist() / 6371.0 * 180.0 / Math.PI;
        } else if( ! SacConstants.isUndef(sacFile.getHeader().getStla()) && ! SacConstants.isUndef(sacFile.getHeader().getStlo())
                && ! SacConstants.isUndef(sacFile.getHeader().getEvla()) && ! SacConstants.isUndef(sacFile.getHeader().getEvlo())) {
            if(verbose) {
                System.out.println("Using stla,stlo, evla,evlo to calculate");
            }
            Alert.warning("Warning: Sac header gcarc is not set,",
                          "using lat and lons to calculate distance.");
            Alert.warning("No ellipticity correction will be applied.",
                          "This may introduce errors. Please see the manual.");
            deg = SphericalCoords.distance(header.getStla(),
                                           header.getStlo(),
                                           header.getEvla(),
                                           header.getEvlo());
        } else {
            /* can't get a distance, skipping */
            Alert.warning("Can't get a distance, all distance fields are undef.",
                          "skipping " + f);
            return;
        }
        if(!((evdpkm && depth == header.getEvdp()) || (!evdpkm && depth == 1000 * header.getEvdp()))) {
            if(!evdpkm && header.getEvdp() != 0 && header.getEvdp() < 1000.0) {
                Alert.warning("Sac header evdp is < 1000 in "
                                      + f,
                              "If the depth is in kilometers instead of meters "
                                      + "(default), you should use the -evdpkm flag");
            }
            if(evdpkm) {
                depthCorrect(header.getEvdp());
            } else {
                depthCorrect(header.getEvdp() / 1000.0);
            }
        }
        if(verbose) {
            System.out.println(f
                    + " searching for " + getPhaseNameString());
        }
        calculate(deg);
        // calcTime(deg);
        if(verbose) {
            System.out.println(f + " "
                    + arrivals.size() + " arrivals found.");
        }
        for(int arrivalNum = arrivals.size() - 1; arrivalNum >= 0; arrivalNum--) {
            int phaseNum = -1;
            for(int j = phaseNames.size() - 1; j >= 0; j--) {
                if(getArrival(arrivalNum).getName()
                        .equals(((PhaseName)phaseNames.get(j)).name)) {
                    phaseNum = j;
                    break;
                }
            }
            if(phaseNum != -1) {
                if(verbose) {
                    System.out.println(f
                            + " phase found "
                            + getArrival(arrivalNum).getName()
                            + " -> t"
                            + ((PhaseName)phaseNames.get(phaseNum)).sacTNum
                            + ", travel time="
                            + (float)getArrival(arrivalNum).getTime());
                }
                setSacTHeader(sacFile,
                              ((PhaseName)phaseNames.get(phaseNum)).sacTNum,
                              getArrival(arrivalNum));
            }
        }
        sacFile.write(f);
    }

    public static void setSacTHeader(SacTimeSeries sacFile,
                                     int headerNum,
                                     Arrival arrival) {
        float arrivalTime = sacFile.getHeader().getO() + (float)arrival.getTime();
        if(headerNum == A_HEADER) {
                // there is no t10, so use that for the A header
            sacFile.getHeader().setA( arrivalTime);
            sacFile.getHeader().setKa( arrival.getName());
                // no place to put the ray param
        } else {
            sacFile.getHeader().setTHeader(headerNum, arrivalTime, arrival.getName());
            sacFile.getHeader().setUserHeader(headerNum, (float)arrival.getRayParam());
        }
    }

    public void printStdUsage() {
        String className = this.getClass().getName();
        className = className.substring(className.lastIndexOf('.') + 1,
                                        className.length());
        System.out.println("Usage: " + className.toLowerCase() + " [arguments]");
        System.out.println("  or, for purists, java "
                + this.getClass().getName() + " [arguments]");
        System.out.println("\nArguments are:");
        System.out.println("-ph phase list     -- comma separated phase list,\n"
                + "                      use phase-# to specify the sac header,\n"
                + "                      for example, ScS-8 puts ScS in t8\n"
                + "-pf phasefile      -- file containing phases\n\n"
                + "-mod[el] modelname -- use velocity model \"modelname\" for calculations\n"
                + "                      Default is iasp91.\n\n");
    }

    public void printStdUsageTail() {
        System.out.println("\n-debug             -- enable debugging output\n"
                + "-verbose           -- enable verbose output\n"
                + "-version           -- print the version\n"
                + "-help              -- print this out, but you already know that!\n");
    }

    public void printUsage() {
        printStdUsage();
        System.out.println("-evdpkm            -- sac depth header is in km, default is meters\n");
        printStdUsageTail();
        System.out.println("sacfilename [sacfilename ...]");
        System.out.println("\nEx: taup_setsac "
                + "-mod S_prem -ph S-8,ScS-9 wmq.r wmq.t wmq.z");
        System.out.println("puts the first S arrival in T8 and ScS in T9");
    }

    public String[] parseCmdLineArgs(String[] args) throws IOException {
        int i = 0;
        String[] leftOverArgs;
        int numNoComprendoArgs = 0;
        File tempFile;
        leftOverArgs = super.parseCmdLineArgs(args);
        String[] noComprendoArgs = new String[leftOverArgs.length];
        while(i < leftOverArgs.length) {
            if(dashEquals("evdpkm", leftOverArgs[i])) {
                evdpkm = true;
            } else if(dashEquals("help", leftOverArgs[i])) {
                noComprendoArgs[numNoComprendoArgs++] = leftOverArgs[i];
            } else {
                tempFile = new File(leftOverArgs[i]);
                if(tempFile.exists() && tempFile.isFile() && tempFile.canRead()) {
                    sacFileNames.add(leftOverArgs[i]);
                } else {
                    if(! tempFile.exists()) {
                        System.err.println(leftOverArgs[i]+" does not exist. "+tempFile.getAbsolutePath() );
                    } else if( ! tempFile.isFile() ) {
                        System.err.println(leftOverArgs[i]+" is not a file.");
                    } else if( ! tempFile.canRead()) {
                        System.err.println(leftOverArgs[i]+" is not readable.");
                    }
                    noComprendoArgs[numNoComprendoArgs++] = leftOverArgs[i];
                }
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

    /**
     * Allows TauP_SetSac to run as an application. Creates an instance of
     * TauP_SetSac. .
     */
    public static void main(String[] args) throws FileNotFoundException,
            IOException, StreamCorruptedException, ClassNotFoundException,
            OptionalDataException {
        TauP_SetSac tauPSetSac = new TauP_SetSac();
        if(args.length == 0) {
            tauPSetSac.printUsage();
            System.exit(1);
        } else
            try {
                String[] noComprendoArgs = tauPSetSac.parseCmdLineArgs(args);
                printNoComprendoArgs(noComprendoArgs);
                if(TauP_Time.DEBUG) {
                    System.out.println("Done reading " + tauPSetSac.modelName);
                }
                tauPSetSac.init();
                tauPSetSac.start();
            } catch(TauModelException e) {
                System.out.println("Caught TauModelException: "
                        + e.getMessage());
                e.printStackTrace();
            }
    }
}
