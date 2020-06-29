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

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OptionalDataException;
import java.io.StreamCorruptedException;
import java.io.StreamTokenizer;
import java.util.List;

/**
 * Allows peeking into the taumodel, slowness model and velocity model
 * previously create by TauP_Create. Mainly used for debugging purposes and is
 * probably not very useful to end users.
 * 
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * 
 * 
 * 
 * @author H. Philip Crotwell
 * 
 */
public class TauP_Peek {

    TauModel tMod;

    public void readTauModel(String filename) throws FileNotFoundException,
            IOException, StreamCorruptedException, ClassNotFoundException,
            OptionalDataException {
        tMod = TauModel.readModel(filename);
    }

    public static void main(String[] args) throws FileNotFoundException,
            IOException, StreamCorruptedException, ClassNotFoundException,
            OptionalDataException {
        String modelFilename;
        if(args.length == 2
                && (args[0].equals("-mod") || args[0].equals("-model"))) {
            modelFilename = args[1];
        } else {
            modelFilename = "iasp91.taup";
        }
        boolean DEBUG = false;
        TauModel tModDepth;
        SeismicPhase sp;
        double depth = 0.0;
        try {
            TauP_Peek peek = new TauP_Peek();
            peek.readTauModel(modelFilename);
            tModDepth = peek.tMod;
            // just so sp is initialized for compiler
            sp = new SeismicPhase("S", tModDepth);
            StreamTokenizer tokenIn = new StreamTokenizer(new InputStreamReader(System.in));
            tokenIn.parseNumbers();
            System.out.print("seismic (p)hase or (t)au model or (s)lowness model "
                    + "or (v)elocity model or tau (b)ranch? ");
            char modelType = 't';
            int branch, rayNum;
            int layer;
            tokenIn.nextToken();
            while(tokenIn.nval != -1 && !tokenIn.sval.equalsIgnoreCase("q")) {
                if(tokenIn.sval.equalsIgnoreCase("t")) {
                    modelType = 't';
                } else if(tokenIn.sval.equalsIgnoreCase("b")) {
                    modelType = 'b';
                } else if(tokenIn.sval.equalsIgnoreCase("s")) {
                    modelType = 's';
                } else if(tokenIn.sval.equalsIgnoreCase("v")) {
                    modelType = 'v';
                } else if(tokenIn.sval.equalsIgnoreCase("p")) {
                    modelType = 'p';
                } else if(tokenIn.sval.equalsIgnoreCase("r")) {
                    modelType = 'r';
                } else {
                    System.out.println("Unrecognized model type: " + modelType);
                    System.out.println("Using (t)au model.");
                    modelType = 't';
                }
                switch(modelType){
                    case 'r':
                    case 'R':
                        System.out.println("Enter source depth");
                        tokenIn.nextToken();
                        depth = (double)tokenIn.nval;
                        tModDepth = peek.tMod.depthCorrect(depth);
                        double[] rayParams = tModDepth.getRayParams();
                        for(int i = 0; i < rayParams.length; i++) {
                            System.out.print(rayParams[i] + "  ");
                            if(i % 5 == 0) {
                                System.out.println();
                            }
                        }
                        modelType = 'T';
                        break;
                    case 't':
                    case 'T':
                        System.out.println("spherical=" + peek.tMod.isSpherical()
                                + " sourceDepth=" + peek.tMod.getSourceDepth()
                                + " radiusOfEarth="
                                + peek.tMod.getRadiusOfEarth() + " DEBUG="
                                + peek.tMod.DEBUG + " rayParams.length="
                                + peek.tMod.rayParams.length
                                + " tauBranches[0].length="
                                + peek.tMod.tauBranches[0].length
                                + " tauBranches[1].length="
                                + peek.tMod.tauBranches[1].length);
                        for(int i = 0; i < peek.tMod.getNumBranches(); i++) {
                            System.out.println("peek.tMod.tauBranches[0][" + i
                                    + "].dist.length="
                                    + peek.tMod.tauBranches[0][i].dist.length
                                    + " peek.tMod.tauBranches[1][" + i
                                    + "].dist.length="
                                    + peek.tMod.tauBranches[1][i].dist.length);
                        }
                        System.out.println("Enter source depth");
                        tokenIn.nextToken();
                        depth = (double)tokenIn.nval;
                        tModDepth = peek.tMod.depthCorrect(depth);
                        System.out.println("Enter branch rayNum");
                        break;
                    case 's':
                    case 'S':
                        System.out.println(peek.tMod.getSlownessModel());
                        System.out.println("Enter slowness layer");
                        break;
                    case 'v':
                    case 'V':
                        System.out.println(peek.tMod.getVelocityModel());
                        System.out.println("Enter velocity layer");
                        break;
                    case 'p':
                    case 'P':
                        System.out.println("Enter depth");
                        tokenIn.nextToken();
                        depth = tokenIn.nval;
                        tModDepth = peek.tMod.depthCorrect(depth);
                        System.out.println("Enter phase name");
                        tokenIn.nextToken();
                        sp = new SeismicPhase(tokenIn.sval, tModDepth);
                        System.out.println("Enter degrees");
                        break;
                    case 'b':
                    case 'B':
                        System.out.println("Enter source depth");
                        tokenIn.nextToken();
                        depth = (double)tokenIn.nval;
                        tModDepth = peek.tMod.depthCorrect(depth);
                        System.out.println("Enter Branch");
                        break;
                    default:
                        System.out.println("Unrecognized model type: "
                                + modelType);
                        System.out.println("Using (t)au model.");
                        modelType = 't';
                        break;
                }
                tokenIn.nextToken();
                while(tokenIn.ttype == StreamTokenizer.TT_NUMBER) {
                    switch(modelType){
                        case 't':
                        case 'T':
                            branch = (int)tokenIn.nval;
                            tokenIn.nextToken();
                            rayNum = (int)tokenIn.nval;
                            if(branch >= tModDepth.getNumBranches()
                                    || rayNum >= tModDepth.rayParams.length) {
                                System.out.println("Out of bounds!");
                            } else {
                                System.out.println("ray parameter="
                                        + tModDepth.rayParams[rayNum]
                                        + " distance="
                                        + tModDepth.tauBranches[0][branch].getDist(rayNum)
                                        + " time="
                                        + tModDepth.tauBranches[0][branch].time[rayNum]
                                        + " tau="
                                        + tModDepth.tauBranches[0][branch].tau[rayNum]);
                                System.out.println("ray parameter="
                                        + tModDepth.rayParams[rayNum]
                                        + " distance="
                                        + tModDepth.tauBranches[1][branch].getDist(rayNum)
                                        + " time="
                                        + tModDepth.tauBranches[1][branch].time[rayNum]
                                        + " tau="
                                        + tModDepth.tauBranches[1][branch].tau[rayNum]);
                            }
                            System.out.println("Enter branch rayNum");
                            break;
                        case 's':
                        case 'S':
                            layer = (int)tokenIn.nval;
                            if(layer >= peek.tMod.getSlownessModel().getNumLayers(true)) {
                                System.out.println("P wave Out of bounds!");
                            } else {
                                System.out.println(peek.tMod.getSlownessModel().getSlownessLayer(layer,
                                                                                   true));
                            }
                            if(layer >= peek.tMod.getSlownessModel().getNumLayers(false)) {
                                System.out.println("S wave Out of bounds!");
                            } else {
                                System.out.println(peek.tMod.getSlownessModel().getSlownessLayer(layer,
                                                                                   false));
                            }
                            System.out.println("Enter slowness layer");
                            break;
                        case 'b':
                        case 'B':
                            branch = (int)tokenIn.nval;
                            if(branch >= tModDepth.getNumBranches()) {
                                System.out.println("Out of bounds!");
                            } else {
                                System.out.println(tModDepth.tauBranches[0][branch]);
                                System.out.println(tModDepth.tauBranches[1][branch]);
                            }
                            System.out.println("Enter Branch");
                            break;
                        case 'v':
                        case 'V':
                            layer = (int)tokenIn.nval;
                            if(layer >= peek.tMod.getSlownessModel().vMod.getNumLayers()) {
                                System.out.println("Out of bounds!");
                            } else {
                                System.out.println(peek.tMod.getSlownessModel().vMod.getVelocityLayer(layer));
                            }
                            System.out.println("Enter velocity layer");
                            break;
                        case 'p':
                        case 'P':
                            List<Arrival> arr = sp.calcPierce(tokenIn.nval);
                            System.out.println("MaxRayParamIndex="
                                    + sp.getMaxRayParamIndex()
                                    + " MinRayParamIndex="
                                    + sp.getMinRayParamIndex());
                            for (Arrival arrival : arr) {
                                System.out.println(arrival);
                            }
                            System.out.println("Enter degrees");
                    }
                    tokenIn.nextToken();
                }
                System.out.print("(t)au model or (s)lowness model or (v)elocity model? ");
                tokenIn.nextToken();
            }
        } catch(TauModelException e) {
            System.out.println("Caught TauModelException " + e.getMessage());
            e.printStackTrace();
        } catch(IOException e) {
            System.out.println("Tried to read!\n Caught IOException "
                    + e.getMessage());
        } finally {
            System.out.println("Done!\n");
        }
    }
}
