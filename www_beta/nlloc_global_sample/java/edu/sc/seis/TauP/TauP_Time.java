/*
 * The TauP Toolkit: Flexible Seismic Travel-Time and Raypath Utilities.
 * Copyright (C) 1998-2000 University of South Carolina This program is free
 * software; you can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version. This program
 * is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU General Public License for more details. You
 * should have received a copy of the GNU General Public License along with this
 * program; if not, write to the Free Software Foundation, Inc., 59 Temple Place -
 * Suite 330, Boston, MA 02111-1307, USA. The current version can be found at <A
 * HREF="www.seis.sc.edu">http://www.seis.sc.edu </A> Bug reports and comments
 * should be directed to H. Philip Crotwell, crotwell@seis.sc.edu or Tom Owens,
 * owens@seis.sc.edu
 */
package edu.sc.seis.TauP;

import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.InvalidClassException;
import java.io.OptionalDataException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.StreamCorruptedException;
import java.io.StreamTokenizer;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Properties;

/**
 * Calculate travel times for different branches using linear interpolation
 * between known slowness samples.
 * 
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * @author H. Philip Crotwell
 */
public class TauP_Time {

    /** Turns on debugging output. */
    public static boolean DEBUG = false;

    /** Turns on verbose output. */
    public boolean verbose = false;

    /** Turns on expert mode. */
    public static boolean expert = false;

    protected String modelName = "iasp91";

    /**
     * Tau model calculated previously.
     * 
     * @see TauModel
     */
    protected TauModel tMod;

    /**
     * TauModel derived from tMod by correcting it for a non-surface source.
     */
    protected transient TauModel tModDepth;

    /**
     * vector to hold the SeismicPhases for the phases named in phaseNames.
     */
    protected List<SeismicPhase> phases = new ArrayList<SeismicPhase>(10);

    /** names of phases to be used, ie PKIKP. */
    protected List<PhaseName> phaseNames = new ArrayList<PhaseName>();

    protected double depth = 0.0;

    protected double degrees = Double.MAX_VALUE;

    protected double azimuth = Double.MAX_VALUE;

    protected double backAzimuth = Double.MAX_VALUE;

    protected double stationLat = Double.MAX_VALUE;

    protected double stationLon = Double.MAX_VALUE;

    protected double eventLat = Double.MAX_VALUE;

    protected double eventLon = Double.MAX_VALUE;

    protected List<Arrival> arrivals = new ArrayList<Arrival>();

    protected boolean GUI = false;

    protected boolean onlyPrintRayP = false;

    protected boolean onlyPrintTime = false;
    
    protected String relativePhaseName = "";
    
    protected Arrival relativeArrival;

    protected String outFile = "";

    protected PrintWriter writer;

    protected Properties toolProps;

    protected Outputs outForms;

    /* Constructors */
    protected TauP_Time() {
        try {
            toolProps = PropertyLoader.load();
        } catch(Exception e) {
            Alert.warning("Unable to load properties, using defaults.",
                          e.getMessage());
            toolProps = new Properties();
        }
        Outputs.configure(toolProps);
    }

    public TauP_Time(TauModel tMod) throws TauModelException {
        this();
        this.tMod = tMod;
        this.tModDepth = tMod;
        modelName = tMod.getModelName();
    }

    /**
     * creates a TauP_Time object with the tau model specified by modelName
     * already loaded.
     * 
     * @throws TauModelException
     *             if the file can't be found or is corrupted in some way.
     */
    public TauP_Time(String modelName) throws TauModelException {
        this();
        try {
            loadTauModel(modelName);
        } catch(FileNotFoundException e) {
            throw new TauModelException("FileNotFoundException:"
                    + e.getMessage(), e);
        } catch(InvalidClassException e) {
            throw new TauModelException("InvalidClassException:"
                    + e.getMessage(), e);
        } catch(StreamCorruptedException e) {
            throw new TauModelException("StreamCorruptedException:"
                    + e.getMessage(), e);
        } catch(OptionalDataException e) {
            throw new TauModelException("OptionalDataException:"
                    + e.getMessage(), e);
        } catch(IOException e) {
            throw new TauModelException("IOException:" + e.getMessage(), e);
        }
    }

    /* Get/Set methods */
    public String[] getPhaseNames() {
        String[] phases = new String[phaseNames.size()];
        for(int i = 0; i < phaseNames.size(); i++) {
            phases[i] = phaseNames.get(i).getName();
        }
        return phases;
    }

    public String getPhaseNameString() {
        // in case of empty phase list
        if(getNumPhases() == 0)
            return "";
        String phases = phaseNames.get(0).getName();
        for(int i = 1; i < getNumPhases(); i++) {
            phases += "," + phaseNames.get(i).getName();
        }
        return phases;
    }

    public void setPhaseNames(String[] phaseNames) throws TauModelException {
        this.phaseNames.clear();
        for(int i = 0; i < phaseNames.length; i++) {
            appendPhaseName(phaseNames[i]);
        }
    }

    public void setPhaseNames(PhaseName[] phaseNames) {
        this.phaseNames.clear();
        for(int i = 0; i < phaseNames.length; i++) {
            this.phaseNames.add(phaseNames[i]);
        }
    }

    public static List<String> getPhaseNames(String phaseName) {
        List<String> names = new ArrayList<String>();
        if(phaseName.equalsIgnoreCase("ttp")
                || phaseName.equalsIgnoreCase("tts")
                || phaseName.equalsIgnoreCase("ttbasic")
                || phaseName.equalsIgnoreCase("tts+")
                || phaseName.equalsIgnoreCase("ttp+")
                || phaseName.equalsIgnoreCase("ttall")) {
            if(phaseName.equalsIgnoreCase("ttp")
                    || phaseName.equalsIgnoreCase("ttp+")
                    || phaseName.equalsIgnoreCase("ttbasic")
                    || phaseName.equalsIgnoreCase("ttall")) {
                names.add("p");
                names.add("P");
                names.add("Pn");
                names.add("Pdiff");
                names.add("PKP");
                names.add("PKiKP");
                names.add("PKIKP");
            }
            if(phaseName.equalsIgnoreCase("tts")
                    || phaseName.equalsIgnoreCase("tts+")
                    || phaseName.equalsIgnoreCase("ttbasic")
                    || phaseName.equalsIgnoreCase("ttall")) {
                names.add("s");
                names.add("S");
                names.add("Sn");
                names.add("Sdiff");
                names.add("SKS");
                names.add("SKIKS");
            }
            if(phaseName.equalsIgnoreCase("ttp+")
                    || phaseName.equalsIgnoreCase("ttbasic")
                    || phaseName.equalsIgnoreCase("ttall")) {
                names.add("PcP");
                names.add("pP");
                names.add("pPdiff");
                names.add("pPKP");
                names.add("pPKIKP");
                names.add("pPKiKP");
                names.add("sP");
                names.add("sPdiff");
                names.add("sPKP");
                names.add("sPKIKP");
                names.add("sPKiKP");
            }
            if(phaseName.equalsIgnoreCase("tts+")
                    || phaseName.equalsIgnoreCase("ttbasic")
                    || phaseName.equalsIgnoreCase("ttall")) {
                names.add("sS");
                names.add("sSdiff");
                names.add("sSKS");
                names.add("sSKIKS");
                names.add("ScS");
                names.add("pS");
                names.add("pSdiff");
                names.add("pSKS");
                names.add("pSKIKS");
            }
            if(phaseName.equalsIgnoreCase("ttbasic")
                    || phaseName.equalsIgnoreCase("ttall")) {
                names.add("ScP");
                names.add("SKP");
                names.add("SKIKP");
                names.add("PKKP");
                names.add("PKIKKIKP");
                names.add("SKKP");
                names.add("SKIKKIKP");
                names.add("PP");
                names.add("PKPPKP");
                names.add("PKIKPPKIKP");
            }
            if(phaseName.equalsIgnoreCase("ttall")) {
                names.add("SKiKP");
                names.add("PP");
                names.add("ScS");
                names.add("PcS");
                names.add("PKS");
                names.add("PKIKS");
                names.add("PKKS");
                names.add("PKIKKIKS");
                names.add("SKKS");
                names.add("SKIKKIKS");
                names.add("SKSSKS");
                names.add("SKIKSSKIKS");
                names.add("SS");
                names.add("SP");
                names.add("PS");
            }
        } else {
            names.add(phaseName);
        }
        return names;
    }

    public synchronized void appendPhaseName(String phaseName)
            throws TauModelException {
        Iterator<String> it = getPhaseNames(phaseName).iterator();
        while(it.hasNext()) {
            appendPhaseName(new PhaseName(it.next()));
        }
    }

    public synchronized void appendPhaseName(PhaseName phaseName) {
        boolean unique = true;
        if(phaseName.name == null || phaseName.name.length() == 0) {
            // make sure not null string
            return;
        }
        for(int i = 0; i < phaseNames.size(); i++) {
            if(phaseNames.get(i).equals(phaseName)) {
                unique = false;
                return;
            }
        }
        if(unique) {
            this.phaseNames.add(phaseName);
        }
    }

    public int getNumPhases() {
        return phaseNames.size();
    }

    public void clearPhaseNames() {
        phases.clear();
        phaseNames.clear();
    }

    public double getSourceDepth() {
        return Double.valueOf(toolProps.getProperty("taup.source.depth", "0.0"))
                .doubleValue();
    }

    public void setSourceDepth(double depth) {
        this.depth = depth;
        toolProps.put("taup.source.depth", Double.toString(depth));
    }

    public String getTauModelName() {
        return modelName;
    }

    public TauModel getTauModel() {
        return tMod;
    }

    public void setTauModel(TauModel tMod) {
        this.tMod = tMod;
        this.tModDepth = tMod;
        modelName = tMod.getModelName();
        toolProps.put("taup.model.name", modelName);
    }

    public void loadTauModel(String modelName) throws FileNotFoundException,
            InvalidClassException, IOException, StreamCorruptedException,
            OptionalDataException, TauModelException {
        this.modelName = modelName;
        readTauModel();
        this.modelName = tMod.getModelName();
    }

    public double[] getDisconDepths() {
        return tMod.getVelocityModel().getDisconDepths();
    }

    public void clearArrivals() {
        arrivals = new ArrayList<Arrival>();
    }

    public int getNumArrivals() {
        return arrivals.size();
    }

    public Arrival getArrival(int i) {
        return (Arrival)arrivals.get(i);
    }

    public List<Arrival> getArrivals() {
        return Collections.unmodifiableList(arrivals);
    }
    
    public List<SeismicPhase> getSeismicPhases() {
        return Collections.unmodifiableList(phases);
    }

    /* Normal methods */
    /**
     * Reads the velocity model, slowness model, and tau model from a file saved
     * using Java's Serializable interface. Performs a depth correction if the
     * current depth is not 0.0
     */
    protected void readTauModel() throws TauModelException {
            TauModel tModLoad = TauModelLoader.load(modelName,
                                                    toolProps.getProperty("taup.model.path"),
                                                    verbose);
            if(tModLoad != null) {
                tMod = tModLoad;
                tModDepth = tMod;
                this.modelName = tMod.getModelName();
            } else {
                throw new TauModelException("Unable to load "+modelName);
            }
    }

    /**
     * Reads in list of phase names from a text file. So long as each phase name
     * is separated by some whitespace, " " or newline or tab, it should read
     * them fine. Also, comments are allowed, either # or // are comments to the
     * end of the line while c style slash-star make a block a comment.
     */
    protected void readPhaseFile(String filename) throws IOException {
        FileReader fileIn = new FileReader(filename);
        StreamTokenizer tokenIn = new StreamTokenizer(fileIn);
        tokenIn.commentChar('#'); // '#' means ignore to end of line
        tokenIn.slashStarComments(true); // '/*...*/' means a comment
        tokenIn.slashSlashComments(true); // '//' means ignore to end of line
        tokenIn.wordChars('^', '^');
        tokenIn.wordChars('0', '9');
        tokenIn.wordChars('.', '.');
        tokenIn.wordChars('[', '[');
        tokenIn.wordChars(']', ']');
        while(tokenIn.nextToken() != StreamTokenizer.TT_EOF) {
            if(tokenIn.sval != null) {
                parsePhaseList(tokenIn.sval);
            } else {
                if(DEBUG) {
                    Alert.info("Token.sval was null! nval=" + tokenIn.nval);
                }
            }
        }
    }

    /**
     * parses a comma separated list of phase names and adds them to the
     * phaseNames vector. Each phase can have an optional argument after a dash.
     * This would be used for specifying which sac header to put the time in, or
     * for other unforeseen uses. This may be called multiple times to append
     * more phases. For example: P-0,PcP-1,ScP-4,Sn,SS,S^410S would, assuming no
     * previous phases have been added, put P in T0, PcP in T1, ScP in T5, Sn in
     * T2, SS in T3, and S^410S in T5.
     */
    public void parsePhaseList(String phaseList) {
        String phaseEntry = "";
        phaseList = phaseList.trim();
        phaseList = phaseList.replace(' ', ',');
        // remove any empty phases, ie two commas next to each other
        // should be replaced with one comma
        phaseList = phaseList.replaceAll(",,+", ",");
        // remove comma at beginning
        if(phaseList.startsWith(",")) {
            if(phaseList.length() > 1) {
                phaseList = phaseList.substring(1);
            } else {
                // phaseList is just a single comma, no phases, so just return
                return;
            }
        }
        // and comma at end
        if(phaseList.charAt(phaseList.length() - 1) == ',') {
            // we know that the length is > 1 as if not then we would have
            // returned from the previous if
            phaseList = phaseList.substring(0, phaseList.length() - 1);
        }
        String[] namesInList = phaseList.split(",");
        for(int i = 0; i < namesInList.length; i++) {
            String[] phaseAndHeader = namesInList[i].split("-");
            try {
                if(phaseAndHeader.length == 1) {
                    /* no optional dash argument, so just add the name. */
                    appendPhaseName(phaseAndHeader[0]);
                } else {
                    if(phaseAndHeader[1].length() == 1
                            && Character.isDigit(phaseAndHeader[1].charAt(0))) {
                        /*
                         * There is an optional argument, so store it and the
                         * phase name.
                         */
                        appendPhaseName(new PhaseName(phaseAndHeader[0],
                                                      Integer.valueOf(phaseAndHeader[1])
                                                              .intValue()));
                    } else if(phaseAndHeader[1].length() == 1
                            && phaseEntry.charAt(phaseEntry.length() - 1) == 'a') {
                        /*
                         * There is an optional argument, use 10 for sac A, so
                         * store it and the phase name.
                         */
                        appendPhaseName(new PhaseName(phaseAndHeader[0],
                                                      TauP_SetSac.A_HEADER));
                    } else {
                        Alert.warning("Problem with phase=" + phaseEntry,
                                      "Skipping this phase.");
                    }
                }
            } catch(TauModelException e) {
                Alert.warning("Problem with phase=" + phaseEntry + " "
                        + e.getMessage(), "Skipping this phase: ");
            }
        }
    }

    /**
     * Parses a comma separated list of distances and returns them in an array.
     */
    public List<Double> parseDegreeList(String degList) {
        String[] split = degList.split(",");
        List<Double> degreesFound = new ArrayList<Double>(split.length);
        for (int i = 0; i < split.length; i++) {
            degreesFound.add(Double.parseDouble(split[i]));
        }
        return degreesFound;
    }
    
    public static boolean dashEquals(String argName, String arg) {
        return arg.equalsIgnoreCase("-"+argName) || arg.equalsIgnoreCase("--"+argName);
    }

    /*
     * parses the standard command line args for the taup package. Other tools
     * that subclass this class will likely override this.
     */
    protected String[] parseCmdLineArgs(String[] args) throws IOException {
        int i = 0;
        String[] noComprendoArgs = new String[args.length];
        int numNoComprendoArgs = 0;
        boolean cmdLineArgPhase = false;
        boolean cmdLineArgPhaseFile = false;
        while(i < args.length) {
            if(dashEquals("help", args[i])) {
                printUsage();
                noComprendoArgs[numNoComprendoArgs++] = args[i];
            } else if(dashEquals("version", args[i])) {
                Alert.info(BuildVersion.getDetailedVersion());
                noComprendoArgs[numNoComprendoArgs++] = args[i];
            } else if(dashEquals("verbose", args[i])) {
                verbose = true;
            } else if(dashEquals("expert", args[i])) {
                expert = true;
            } else if(dashEquals("debug", args[i])) {
                verbose = true;
                DEBUG = true;
            } else if(dashEquals("gui", args[i])) {
                GUI = true;
            } else if(dashEquals("rayp", args[i])) {
                onlyPrintRayP = true;
                onlyPrintTime = false;
            } else if(dashEquals("time", args[i])) {
                onlyPrintTime = true;
                onlyPrintRayP = false;
            } else if(i < args.length - 1) {
                if(dashEquals("mod", args[i]) || dashEquals("model", args[i])) {
                    toolProps.put("taup.model.name", args[i + 1]);
                    i++;
                } else if(args[i].equalsIgnoreCase("-h")) {
                    toolProps.put("taup.source.depth", args[i + 1]);
                    i++;
                } else if(dashEquals("deg", args[i])) {
                    degrees = Double.valueOf(args[i + 1]).doubleValue();
                    i++;
                } else if(dashEquals("km", args[i])) {
                    degrees = Double.valueOf(args[i + 1]).doubleValue() / 6371
                            * 180.0 / Math.PI;
                    i++;
                } else if(dashEquals("az", args[i])) {
                    azimuth = Double.valueOf(args[i + 1]).doubleValue();
                    i++;
                } else if(dashEquals("baz", args[i])) {
                    backAzimuth = Double.valueOf(args[i + 1]).doubleValue();
                    i++;
                } else if(args[i].equalsIgnoreCase("-o")) {
                    outFile = args[i + 1];
                    i++;
                } else if(dashEquals("rel", args[i])) {
                    relativePhaseName = args[i + 1];
                    i++;
                } else if(dashEquals("ph", args[i])) {
                    if(cmdLineArgPhase) {
                        // previous cmd line -ph so append
                        toolProps.put("taup.phase.list",
                                      toolProps.getProperty("taup.phase.list",
                                                            "")
                                              + "," + args[i + 1]);
                    } else {
                        // no previous cmd line -ph so replace defaults
                        toolProps.put("taup.phase.list", args[i + 1]);
                    }
                    cmdLineArgPhase = true;
                    i++;
                } else if(dashEquals("pf", args[i])) {
                    cmdLineArgPhaseFile = true;
                    toolProps.put("taup.phase.file", args[i + 1]);
                    i++;
                } else if(i < args.length - 2) {
                    if(args[i].equalsIgnoreCase("-sta")
                            || args[i].equalsIgnoreCase("-station")) {
                        stationLat = Double.valueOf(args[i + 1]).doubleValue();
                        stationLon = Double.valueOf(args[i + 2]).doubleValue();
                        i += 2;
                    } else if(args[i].equalsIgnoreCase("-evt")
                            || args[i].equalsIgnoreCase("-event")) {
                        eventLat = Double.valueOf(args[i + 1]).doubleValue();
                        eventLon = Double.valueOf(args[i + 2]).doubleValue();
                        i += 2;
                    } else {
                        /*
                         * I don't know how to interpret this argument, so pass
                         * it back
                         */
                        noComprendoArgs[numNoComprendoArgs++] = args[i];
                    }
                } else {
                    /*
                     * I don't know how to interpret this argument, so pass it
                     * back
                     */
                    noComprendoArgs[numNoComprendoArgs++] = args[i];
                }
            } else {
                /* I don't know how to interpret this argument, so pass it back */
                noComprendoArgs[numNoComprendoArgs++] = args[i];
            }
            i++;
        }
        // check to see if there were phases or a phase file as an argument.
        // if so then dump the defaults
        if(cmdLineArgPhaseFile || cmdLineArgPhase) {
            if(cmdLineArgPhaseFile && !cmdLineArgPhase) {
                toolProps.remove("taup.phase.list");
            }
            if(!cmdLineArgPhaseFile && cmdLineArgPhase) {
                toolProps.remove("taup.phase.file");
            }
        }
        if(numNoComprendoArgs > 0) {
            String[] temp = new String[numNoComprendoArgs];
            System.arraycopy(noComprendoArgs, 0, temp, 0, numNoComprendoArgs);
            return temp;
        } else {
            return new String[0];
        }
    }

    public synchronized void sortArrivals() {
        Collections.sort(arrivals, new Comparator<Arrival>() {
            public int compare(Arrival o1, Arrival o2) {
                return Double.compare(o1.getTime(), o2.getTime());
            }
        });
    }

    public void calculate(double degrees) throws TauModelException {
        depthCorrect(getSourceDepth());
        recalcPhases();
        calcTime(degrees);
        if (relativePhaseName != "") {
            List<SeismicPhase> relPhases = new ArrayList<SeismicPhase>();
            List<String> splitNames = getPhaseNames(relativePhaseName);
            for (String sName : splitNames) {
                SeismicPhase relPhase = new SeismicPhase(sName, tModDepth);
                relPhases.add(relPhase);
            }
            relativeArrival = SeismicPhase.getEarliestArrival(relPhases, degrees);
        }
    }

    public void calcTime(double degrees) {
        this.degrees = degrees;
        SeismicPhase phase;
        clearArrivals();
        for(int phaseNum = 0; phaseNum < phases.size(); phaseNum++) {
            phase = phases.get(phaseNum);
            List<Arrival> phaseArrivals = phase.calcTime(degrees);
            for (Arrival arrival : phaseArrivals) {
                arrivals.add(arrival);
            }
        }
        sortArrivals();
    }

    /**
     * corrects the TauModel for the given source depth. It only performs the
     * correction of the model is not already corrected to that depth.
     */
    public void depthCorrect(double depth) throws TauModelException {
        if(tModDepth == null || tModDepth.getSourceDepth() != depth) {
            tModDepth = tMod.depthCorrect(depth);
            clearArrivals();
            recalcPhases();
        }
        setSourceDepth(depth);
    }

    /**
     * reclaulates the given phases using a possibly new or changed tau model.
     * This should not need to be called by outside classes as it is called by
     * depthCorrect, and calculate.
     */
    public synchronized void recalcPhases() {
        SeismicPhase seismicPhase;
        List<SeismicPhase> newPhases = new ArrayList<SeismicPhase>(phases.size());
        boolean alreadyAdded;
        String tempPhaseName;
        for(int phaseNameNum = 0; phaseNameNum < phaseNames.size(); phaseNameNum++) {
            tempPhaseName = phaseNames.get(phaseNameNum).getName();
            alreadyAdded = false;
            for(int phaseNum = 0; phaseNum < phases.size(); phaseNum++) {
                seismicPhase = phases.get(phaseNum);
                if(seismicPhase.name.equals(tempPhaseName)) {
                    phases.remove(phaseNum);
                    if(seismicPhase.sourceDepth == depth
                            && seismicPhase.tMod.equals(tModDepth)) {
                        // ok so copy to newPhases
                        newPhases.add(seismicPhase);
                        alreadyAdded = true;
                        if(verbose) {
                            Alert.info(seismicPhase.toString());
                        }
                        break;
                    }
                }
            }
            if(!alreadyAdded) {
                // didn't find it precomputed, so recalculate
                try {
                    seismicPhase = new SeismicPhase(tempPhaseName, tModDepth);
                    newPhases.add(seismicPhase);
                    if(verbose) {
                        Alert.info(seismicPhase.toString());
                    }
                } catch(TauModelException e) {
                    Alert.warning("Error with phase=" + tempPhaseName,
                                  e.getMessage() + "\nSkipping this phase");
                } finally {
                    if(verbose) {
                        Alert.info("-----------------");
                    }
                }
            }
        }
        phases = newPhases;
    }

    public void printResult(PrintWriter out) throws IOException {
        Arrival currArrival;
        int maxNameLength = 5;
        int maxPuristNameLength = 5;
        for(int j = 0; j < arrivals.size(); j++) {
            if(((Arrival)arrivals.get(j)).getName().length() > maxNameLength) {
                maxNameLength = ((Arrival)arrivals.get(j)).getName()
                        .length();
            }
            if(((Arrival)arrivals.get(j)).getPuristName().length() > maxPuristNameLength) {
                maxPuristNameLength = ((Arrival)arrivals.get(j)).getPuristName()
                        .length();
            }
        }
        Format phaseFormat = new Format("%-" + maxNameLength + "s");
        Format phasePuristFormat = new Format("%-" + maxPuristNameLength + "s");
        if(!(onlyPrintRayP || onlyPrintTime)) {
            out.println("\nModel: " + modelName);
            String lineOne = "Distance   Depth   " + phaseFormat.form("Phase")
                    + "   Travel    Ray Param  Takeoff  Incident  Purist    Purist";
            String lineTwo = "  (deg)     (km)   " + phaseFormat.form("Name ")
                    + "   Time (s)  p (s/deg)   (deg)    (deg)   Distance   Name ";
            if (relativePhaseName != "") {
                lineOne += " Relative to";
                for (int s=0; s<(11-relativePhaseName.length())/2;s++) {
                    lineTwo += " ";
                }
                lineTwo += "  "+phaseFormat.form(relativePhaseName);
            }
            out.println(lineOne);
            out.println(lineTwo);
            for(int i = 0; i < lineOne.length(); i++) {
                out.write("-");
            }
            out.write("\n");
            for(int j = 0; j < arrivals.size(); j++) {
                currArrival = (Arrival)arrivals.get(j);
                out.print(Outputs.formatDistance(currArrival.getModuloDistDeg())
                        + Outputs.formatDepth(depth) + "   ");
                out.print(phaseFormat.form(currArrival.getName()));
                out.print("  "
                        + Outputs.formatTime(currArrival.getTime())
                        + "  "
                        + Outputs.formatRayParam(Math.PI / 180.0
                                * currArrival.getRayParam()) + "  ");
                out.print(Outputs.formatDistance(currArrival.getTakeoffAngle())+" ");
                out.print(Outputs.formatDistance(currArrival.getIncidentAngle())+" ");
                out.print(outForms.formatDistance(currArrival.getDistDeg()));
                if(currArrival.getName().equals(currArrival.getPuristName())) {
                    out.print("   = ");
                } else {
                    out.print("   * ");
                }
                out.print(phasePuristFormat.form(currArrival.getPuristName()));
                if (relativePhaseName != "") {
                    out.print(outForms.formatTime(currArrival.getTime() - relativeArrival.getTime()));
                }
                out.println();
            }
        } else if(onlyPrintTime) {
            for(int j = 0; j < arrivals.size(); j++) {
                currArrival = (Arrival)arrivals.get(j);
                out.print(String.valueOf((float)(currArrival.getTime())) + " ");
            }
            out.println();
        } else if(onlyPrintRayP) {
            for(int j = 0; j < arrivals.size(); j++) {
                currArrival = (Arrival)arrivals.get(j);
                out.write(String.valueOf((float)(Math.PI / 180.0 * currArrival.getRayParam()))
                        + " ");
            }
            out.println();
        }
        out.println();
    }

    /**
     * preforms intialization of the tool. Properties are queried for the the
     * default model to load, source depth to use, phases to use, etc. Note that
     * because of the IO inherent in these operations, this method is not
     * appropriate for Applets. Applets should load TauModels themselves and use
     * the setTauModel(TauModel) method.
     */
    public void init() throws IOException {
        if(phaseNames.size() == 0) {
            if(toolProps.containsKey("taup.phase.file")) {
                if(toolProps.containsKey("taup.phase.list")) {
                    parsePhaseList(toolProps.getProperty("taup.phase.list"));
                }
                try {
                    readPhaseFile(toolProps.getProperty("taup.phase.file"));
                } catch(IOException e) {
                    Alert.warning("Caught IOException while attempting to reading phase file "
                                          + toolProps.getProperty("taup.phase.file"),
                                  e.getMessage());
                    if(phaseNames.size() <= 0) {
                        parsePhaseList(toolProps.getProperty("taup.phase.list",
                                                             "p,s,P,S,Pn,Sn,PcP,ScS,Pdiff,Sdiff,PKP,SKS,PKiKP,SKiKS,PKIKP,SKIKS"));
                    }
                }
            } else {
                parsePhaseList(toolProps.getProperty("taup.phase.list",
                                                     "p,s,P,S,Pn,Sn,PcP,ScS,Pdiff,Sdiff,PKP,SKS,PKiKP,SKiKS,PKIKP,SKIKS"));
            }
        }
        depth = Double.valueOf(toolProps.getProperty("taup.source.depth", "0.0"))
                .doubleValue();
        if(tMod == null
                || tMod.getVelocityModel().getModelName() != toolProps.getProperty("taup.model.name",
                                                                                   "iasp91")) {
            modelName = toolProps.getProperty("taup.model.name", "iasp91");
            try {
                readTauModel();
            } catch(TauModelException ee) {
                if (ee.getCause() instanceof InvalidClassException) {
                    Alert.error("Model file "
                                + modelName
                                + " is not compatible with the current version.",
                        "Recreate using taup_create.");
                } else {
                    Alert.error("Caught TauModelException", ee.getMessage());
                }
                throw new RuntimeException(ee);
            }
        }
    }
    
    public PrintWriter getWriter() throws IOException {
        if (writer == null) {
            if(outFile != null && outFile.length() != 0 && !outFile.equals("stdout")) {
                writer = new PrintWriter(new BufferedWriter(new FileWriter(outFile)));
            } else {
                writer = new PrintWriter(new OutputStreamWriter(System.out));
            }
            printScriptBeginning(writer);
        }
        return writer;
    }
    
    /** a noop that allows overriding classes to print things
     * before results are calculated. For example to set up GMT commands before drawing paths.
     * @param out
     * @throws IOException
     */
    public void printScriptBeginning(PrintWriter out)  throws IOException {}

    public void printHelp() {
        Alert.info("Enter:\nh for new depth\nr to recalculate\n"
                + "p to append phases, \nc to clear phases\n"
                + "l to list phases\n"
                + "s for new station lat lon\ne for new event lat lon\n"
                + "a for new azimuth\nb for new back azimuth\n"
                + "m for new model or \nq to quit.\n");
    }

    void doLotsODepths() {
        for(int ii = 0; ii < 1000; ii++) {
            try {
                depthCorrect(Math.random() * 200);
            } catch(TauModelException e) {
                throw new RuntimeException(e);
            }
        }
    }

    public void start() throws IOException, TauModelException, TauPException {
        if((degrees != Double.MAX_VALUE || (stationLat != Double.MAX_VALUE
                && stationLon != Double.MAX_VALUE
                && eventLat != Double.MAX_VALUE && eventLon != Double.MAX_VALUE))) {
            /* enough info given on cmd line, so just do one calc. */
            if(degrees == Double.MAX_VALUE) {
                degrees = SphericalCoords.distance(stationLat,
                                                   stationLon,
                                                   eventLat,
                                                   eventLon);
                azimuth = SphericalCoords.azimuth(eventLat,
                                                  eventLon,
                                                  stationLat,
                                                  stationLon);
                backAzimuth = SphericalCoords.azimuth(stationLat,
                                                      stationLon,
                                                      eventLat,
                                                      eventLon);
            }
            depthCorrect(depth);
            calculate(degrees);
            printResult(getWriter());
        } else {
            /* interactive mode... */
            long prevTime = 0;
            long currTime;
            char readMode = 'd';
            double tempDepth = depth;
            depthCorrect(depth);
            StreamTokenizer tokenIn = new StreamTokenizer(new InputStreamReader(System.in));
            tokenIn.parseNumbers();
            tokenIn.wordChars(',', ',');
            tokenIn.wordChars('_', '_');
            tokenIn.wordChars('^', '^');
            tokenIn.ordinaryChar('/');
            tokenIn.wordChars('/', '/');
            tokenIn.commentChar('#');
            printHelp();
            do {
                switch(readMode){
                    case 'x':
                        doLotsODepths();
                    case 'h':
                        // new source depth
                        System.out.print("Enter Depth: ");
                        tokenIn.nextToken();
                        tempDepth = tokenIn.nval;
                        if(tempDepth < 0.0
                                || tempDepth > tMod.getRadiusOfEarth()) {
                            Alert.warning("Depth must be >= 0.0 and <= tMod.getRadiusOfEarth().",
                                          "depth = " + tempDepth
                                                  + " getRadiusOfEarth= "
                                                  + tMod.getRadiusOfEarth());
                            continue;
                        }
                        prevTime = System.currentTimeMillis();
                        depthCorrect(tempDepth);
                        currTime = System.currentTimeMillis();
                        if(verbose) {
                            Alert.info("depthCorrect time="
                                    + (currTime - prevTime));
                        }
                        readMode = 'd';
                        break;
                    case 'd':
                        // new distance or option
                        System.out.print("Enter Distance or Option [hrpclseabmq]: ");
                        tokenIn.nextToken();
                        if(tokenIn.ttype == StreamTokenizer.TT_NUMBER) {
                            degrees = tokenIn.nval;
                            if(DEBUG) {
                                Alert.info("degrees=" + degrees);
                            }
                            calculate(degrees);
                            printResult(getWriter());
                        } else {
                            if(tokenIn.ttype == StreamTokenizer.TT_EOF
                                    || (tokenIn.ttype == StreamTokenizer.TT_WORD && (tokenIn.sval.equalsIgnoreCase("q")
                                            || tokenIn.sval.equalsIgnoreCase("quit")
                                            || tokenIn.sval.equalsIgnoreCase("exit") || tokenIn.sval.equalsIgnoreCase("bye")))) {
                                readMode = 'q';
                            } else if(tokenIn.ttype == StreamTokenizer.TT_WORD) {
                                if(tokenIn.sval.equalsIgnoreCase("l")) {
                                    readMode = 'l';
                                } else if(tokenIn.sval.equalsIgnoreCase("c")) {
                                    readMode = 'c';
                                } else if(tokenIn.sval.equalsIgnoreCase("s")) {
                                    readMode = 's';
                                } else if(tokenIn.sval.equalsIgnoreCase("e")) {
                                    readMode = 'e';
                                } else if(tokenIn.sval.equalsIgnoreCase("a")) {
                                    readMode = 'a';
                                } else if(tokenIn.sval.equalsIgnoreCase("b")) {
                                    readMode = 'b';
                                } else if(tokenIn.sval.equalsIgnoreCase("r")) {
                                    readMode = 'r';
                                } else if(tokenIn.sval.equalsIgnoreCase("p")) {
                                    readMode = 'p';
                                } else if(tokenIn.sval.equalsIgnoreCase("m")) {
                                    readMode = 'm';
                                } else if(tokenIn.sval.equalsIgnoreCase("h")) {
                                    readMode = 'h';
                                } else if(tokenIn.sval.equalsIgnoreCase("x")) {
                                    readMode = 'x';
                                } else if(tokenIn.sval.equalsIgnoreCase("?")) {
                                    printHelp();
                                } else {
                                    Alert.warning("I don't understand this option",
                                                  tokenIn.sval);
                                    printHelp();
                                }
                            } else {
                                printHelp();
                            }
                        }
                        break;
                    case 'r':
                        // recalulate
                        if(degrees != Double.MAX_VALUE) {
                            calculate(degrees);
                            printResult(getWriter());
                        }
                        readMode = 'd';
                        break;
                    case 'p':
                        // append phases
                        System.out.print("Enter phases (ie P,p,PcP,S): ");
                        tokenIn.ordinaryChars('0', '9');
                        tokenIn.ordinaryChar('.');
                        tokenIn.ordinaryChar('-');
                        tokenIn.wordChars('0', '9');
                        tokenIn.wordChars('.', '.');
                        tokenIn.wordChars('-', '-');
                        tokenIn.ordinaryChar(' ');
                        tokenIn.wordChars(' ', ' ');
                        tokenIn.nextToken();
                        if(tokenIn.ttype == StreamTokenizer.TT_WORD) {
                            parsePhaseList(tokenIn.sval);
                            recalcPhases();
                        } else {
                            Alert.warning("Input phases not recognized.",
                                          "Please retry.");
                        }
                        tokenIn.parseNumbers();
                        tokenIn.ordinaryChar(' ');
                        tokenIn.whitespaceChars(' ', ' ');
                        readMode = 'd';
                        break;
                    case 'l':
                        // list phases
                        int numPhases = phaseNames.size();
                        String output = numPhases + " phases.";
                        Alert.info(output);
                        output = "";
                        for(int i = 0; i < numPhases; i++) {
                            output += phaseNames.get(i).getName();
                            if(i < numPhases - 1) {
                                output += ",";
                            }
                        }
                        Alert.info(output);
                        readMode = 'd';
                        break;
                    case 'c':
                        // clear phases and then enter new phases
                        clearPhaseNames();
                        readMode = 'p';
                        break;
                    case 'a':
                        // event to station azimuth
                        System.out.print("Enter azimuth: ");
                        tokenIn.nextToken();
                        if(tokenIn.ttype == StreamTokenizer.TT_NUMBER) {
                            azimuth = tokenIn.nval;
                            stationLat = Double.MAX_VALUE;
                            stationLon = Double.MAX_VALUE;
                            if(DEBUG) {
                                Alert.info("azimuth=" + azimuth);
                            }
                        } else {
                            Alert.warning("Expected a number.", "got "
                                    + tokenIn + " instead.");
                            printHelp();
                            break;
                        }
                        if(eventLat == Double.MAX_VALUE
                                || eventLon == Double.MAX_VALUE) {
                            readMode = 'e';
                        } else if(degrees == Double.MAX_VALUE) {
                            readMode = 'd';
                        } else {
                            calculate(degrees);
                            printResult(getWriter());
                        }
                        readMode = 'd';
                        break;
                    case 'b':
                        // event to station back azimuth (ie station to event
                        // azimuth)
                        System.out.print("Enter back azimuth: ");
                        tokenIn.nextToken();
                        if(tokenIn.ttype == StreamTokenizer.TT_NUMBER) {
                            backAzimuth = tokenIn.nval;
                            eventLat = Double.MAX_VALUE;
                            eventLon = Double.MAX_VALUE;
                            if(DEBUG) {
                                Alert.info("backAzimuth=" + backAzimuth);
                            }
                        } else {
                            Alert.warning("Expected a number.", "got "
                                    + tokenIn + " instead");
                            printHelp();
                            break;
                        }
                        if(stationLat == Double.MAX_VALUE
                                || stationLon == Double.MAX_VALUE) {
                            readMode = 's';
                        } else if(degrees == Double.MAX_VALUE) {
                            readMode = 'd';
                        } else {
                            calculate(degrees);
                            printResult(getWriter());
                        }
                        readMode = 'd';
                        break;
                    case 'e':
                        // event lat and lon
                        System.out.print("Enter event lat and lon: ");
                        tokenIn.nextToken();
                        if(tokenIn.ttype == StreamTokenizer.TT_NUMBER) {
                            eventLat = tokenIn.nval;
                            if(DEBUG) {
                                Alert.info("eventLat=" + eventLat);
                            }
                            tokenIn.nextToken();
                            if(tokenIn.ttype == StreamTokenizer.TT_NUMBER) {
                                eventLon = tokenIn.nval;
                                if(DEBUG) {
                                    Alert.info("eventLon=" + eventLon);
                                }
                            } else {
                                printHelp();
                            }
                        } else {
                            printHelp();
                        }
                        if(stationLat != Double.MAX_VALUE
                                && stationLon != Double.MAX_VALUE) {
                            degrees = SphericalCoords.distance(stationLat,
                                                               stationLon,
                                                               eventLat,
                                                               eventLon);
                            azimuth = SphericalCoords.azimuth(eventLat,
                                                              eventLon,
                                                              stationLat,
                                                              stationLon);
                            backAzimuth = SphericalCoords.azimuth(stationLat,
                                                                  stationLon,
                                                                  eventLat,
                                                                  eventLon);
                            calculate(degrees);
                            printResult(getWriter());
                        }
                        readMode = 'd';
                        break;
                    case 's':
                        // station lat and lon
                        System.out.print("Enter station lat and lon: ");
                        tokenIn.nextToken();
                        if(tokenIn.ttype == StreamTokenizer.TT_NUMBER) {
                            stationLat = tokenIn.nval;
                            if(DEBUG) {
                                Alert.info("stationLat=" + stationLat);
                            }
                            tokenIn.nextToken();
                            if(tokenIn.ttype == StreamTokenizer.TT_NUMBER) {
                                stationLon = tokenIn.nval;
                                if(DEBUG) {
                                    Alert.info("stationLon=" + stationLon);
                                }
                            } else {
                                printHelp();
                                break;
                            }
                        } else {
                            printHelp();
                            break;
                        }
                        if(eventLat != Double.MAX_VALUE
                                && eventLon != Double.MAX_VALUE) {
                            degrees = SphericalCoords.distance(stationLat,
                                                               stationLon,
                                                               eventLat,
                                                               eventLon);
                            azimuth = SphericalCoords.azimuth(eventLat,
                                                              eventLon,
                                                              stationLat,
                                                              stationLon);
                            backAzimuth = SphericalCoords.azimuth(stationLat,
                                                                  stationLon,
                                                                  eventLat,
                                                                  eventLon);
                            calculate(degrees);
                            printResult(getWriter());
                        }
                        readMode = 'd';
                        break;
                    case 'm':
                        // change model
                        tokenIn.ordinaryChars('0', '9');
                        tokenIn.wordChars('0', '9');
                        tokenIn.ordinaryChars('.', '.');
                        tokenIn.wordChars('.', '.');
                        tokenIn.ordinaryChars('-', '-');
                        tokenIn.wordChars('-', '-');
                        String oldModelName = modelName;
                        TauModel oldTMod = tMod;
                        TauModel oldTModDepth = tModDepth;
                        System.out.print("Enter model name: ");
                        tokenIn.nextToken();
                        if(tokenIn.ttype == StreamTokenizer.TT_WORD) {
                            modelName = tokenIn.sval;
                        }
                        tokenIn.ordinaryChars('0', '9');
                        tokenIn.ordinaryChars('.', '.');
                        tokenIn.ordinaryChars('-', '-');
                        tokenIn.parseNumbers();
                        if(!modelName.equals(oldModelName)) {
                            try {
                                readTauModel();
                                depthCorrect(depth);
                            } catch(TauModelException e) {
                                if (e.getCause() instanceof InvalidClassException) {
                                    Alert.warning("Model file "
                                                  + modelName
                                                  + " is not compatible with the current version.",
                                          "Recreate using taup_create. Still using model "
                                                  + oldModelName + ".");
                                } else {
                                    Alert.warning("I can't load model file "
                                                  + modelName, "Still using model "
                                                  + oldModelName + ".");
                                }
                                modelName = oldModelName;
                                tMod = oldTMod;
                                tModDepth = oldTModDepth;
                            }
                        }
                        readMode = 'd';
                        break;
                    case 'q':
                        return;
                }
            } while(tokenIn.ttype == StreamTokenizer.TT_NUMBER
                    || tokenIn.ttype != StreamTokenizer.TT_WORD
                    || (tokenIn.ttype == StreamTokenizer.TT_WORD && !tokenIn.sval.equalsIgnoreCase("q")));
        }
    }

    public void destroy() throws IOException {
        if(writer != null) {
            writer.close();
            writer = null;
        }
    }

    public void printStdUsageHead() {
        printStdUsageHead(this.getClass());
    }
    
    public static void printStdUsageHead(Class toolClass) {
        String className = toolClass.getName();
        className = className.substring(className.lastIndexOf('.') + 1,
                                        className.length());
        Alert.info("Usage: " + className.toLowerCase() + " [arguments]");
        Alert.info("  or, for purists, java " + toolClass.getName()
                + " [arguments]");
        Alert.info("\nArguments are:");
    }

    /** Prints the command line arguments common to all TauP tools. */
    public void printStdUsage() {
        printStdUsageHead();
        Alert.info("-ph phase list     -- comma separated phase list\n"
                + "-pf phasefile      -- file containing phases\n\n"
                + "-mod[el] modelname -- use velocity model \"modelname\" for calculations\n"
                + "                      Default is iasp91.\n\n"
                + "-h depth           -- source depth in km\n\n"
                + "Distance is given by:\n\n"
                + "-deg degrees       -- distance in degrees,\n"
                + "-km kilometers     -- distance in kilometers,\n"
                + "                      assumes radius of earth is 6371km,\n\n"
                + "or by giving the station and event latitude and lonitude,\n"
                + "                      assumes a spherical earth,\n\n"
                + "-sta[tion] lat lon -- sets the station latitude and longitude\n"
                + "-evt       lat lon -- sets the event latitude and longitude\n\n");
    }

    public void printStdUsageTail() {
        Alert.info("\n-o [stdout|outfile]         -- output is redirected to stdout or to the \"outfile\" file\n"
                + "--debug             -- enable debugging output\n"
                + "--verbose           -- enable verbose output\n"
                + "--version           -- print the version\n"
                + "--help              -- print this out, but you already know that!\n");
    }

    public void printUsage() {
        printStdUsage();
        Alert.info("-rayp              -- only output the ray parameter\n"
                 + "-time              -- only output travel time\n"
                 + "-rel phasename     -- also output relative travel time");
        printStdUsageTail();
    }

    public static void printNoComprendoArgs(String[] noComprendoArgs) {
        if(noComprendoArgs.length > 0) {
            for(int i = 0; i < noComprendoArgs.length; i++) {
                if(noComprendoArgs[i].equals("-help") || noComprendoArgs[i].equals("--help")
                        || noComprendoArgs[i].equals("-version")
                        || noComprendoArgs[i].equals("--version")) {
                    System.exit(0);
                }
            }
            String outStringA = "I don't understand the following arguments, continuing:";
            String outStringB = "";
            for(int i = 0; i < noComprendoArgs.length; i++) {
                outStringB += noComprendoArgs[i] + " ";
            }
            Alert.warning(outStringA, outStringB);
            noComprendoArgs = null;
        }
    }
    /**
     * Allows TauP_Time to run as an application. Creates an instance of
     * TauP_Time. .
     */
    public static void main(String[] args) throws FileNotFoundException,
            IOException, StreamCorruptedException, ClassNotFoundException,
            OptionalDataException {
        // BasicConfigurator.configure();
        try {
            long prevTime = 0;
            long currTime;
            prevTime = System.currentTimeMillis();
            TauP_Time tauPTime = new TauP_Time();
            String[] noComprendoArgs = tauPTime.parseCmdLineArgs(args);
            printNoComprendoArgs(noComprendoArgs);
            currTime = System.currentTimeMillis();
            prevTime = System.currentTimeMillis();
            tauPTime.init();
            currTime = System.currentTimeMillis();
            if(TauP_Time.DEBUG) {
                Alert.info("taup model read time=" + (currTime - prevTime));
            }
            tauPTime.start();
            tauPTime.destroy();
        } catch(TauModelException e) {
            Alert.error("Caught TauModelException", e.getMessage());
            e.printStackTrace();
        } catch(TauPException e) {
            Alert.error("Caught TauPException", e.getMessage());
            e.printStackTrace();
        }
    }

}
