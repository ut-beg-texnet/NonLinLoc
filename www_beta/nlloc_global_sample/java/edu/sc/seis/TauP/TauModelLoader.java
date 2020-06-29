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
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.InvalidClassException;
import java.io.OptionalDataException;
import java.io.Reader;
import java.io.StreamCorruptedException;
import java.lang.ref.SoftReference;
import java.util.HashMap;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

/**
 * static class that loads a tau model, after searching for it. It can be extended to change the search mechanism.
 *
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 *
 *
 *
 * @author H. Philip Crotwell
 *
 */
public class TauModelLoader {

    protected static String packageName = "/edu/sc/seis/TauP/StdModels";

    public static TauModel load(String modelName) throws TauModelException {
        return load(modelName, System.getProperty("taup.model.path"));
    }

    /**
     * Reads the velocity model, slowness model, and tau model from a file saved using Java's Serializable interface.
     */
    public static TauModel load(String modelName, String searchPath)
            throws TauModelException {
        return load(modelName, searchPath, false);
    }

    public static TauModel load(String modelName,
            String searchPath,
            boolean verbose) throws TauModelException {
        TauModel out = loadFromCache(modelName);
        if (out == null) {
            out = internalLoad(modelName, searchPath, verbose);
            tModCache.put(modelName, new SoftReference<TauModel>(out));
        }
        return out;
    }

    public static TauModel internalLoad(String modelName,
            String searchPath,
            boolean verbose) throws TauModelException {
        try {
            String filename;
            /* Append ".taup" to modelname if it isn't already there. */
            if (modelName.endsWith(".taup")) {
                filename = modelName;
            } else {
                filename = modelName + ".taup";
            }
            /* 20140514 AJL - avoid pre-mature exception in applet
             String classPath = System.getProperty("java.class.path");
             20140514 AJL */
            String taupPath = searchPath;
            int offset = 0;
            int pathSepIndex;
            String pathEntry;
            File jarFile;
            File modelFile;
            /* First we try to find the model in the distributed taup.jar file. */
            try {
                Class c = null;
                c = Class.forName("edu.sc.seis.TauP.TauModelLoader");
                InputStream in = c.getResourceAsStream(packageName + "/" + filename);
                if (in != null) {
                    return TauModel.readModelFromStream(in);
                }
            } catch (Exception ex) {
                // couldn't get as a resource, so keep going
                if (verbose) {
                    System.out.println("couldn't load as resource: " + filename
                            + "\n message: " + ex.getMessage());
                }
            }
            /* couldn't find as a resource, try in classpath. */
            // 20140514 AJL - avoid pre-mature exception in applet
            String classPath = System.getProperty("java.class.path");
            // 20140514 AJL
            while (offset < classPath.length()) {
                pathSepIndex = classPath.indexOf(File.pathSeparatorChar, offset);
                if (pathSepIndex != -1) {
                    pathEntry = classPath.substring(offset, pathSepIndex);
                    offset = pathSepIndex + 1;
                } else {
                    pathEntry = classPath.substring(offset);
                    offset = classPath.length();
                }
                jarFile = new File(pathEntry);
                if (jarFile.exists() && jarFile.isFile()
                        && jarFile.getName().equals("taup.jar")
                        && jarFile.canRead()) {
                    ZipFile zippy = new ZipFile(jarFile);
                    ZipEntry zipEntry = zippy.getEntry("StdModels/" + filename);
                    if (zipEntry != null) {
                        return TauModel.readModelFromStream(zippy.getInputStream(zipEntry));
                    }
                }
            }
            /*
             * It isn't in the taup.jar so we try to find it within the paths
             * specified in the taup.model.path property.
             */
            offset = 0;
            if (taupPath != null) {
                while (offset < taupPath.length()) {
                    pathSepIndex = taupPath.indexOf(File.pathSeparatorChar, offset);
                    if (pathSepIndex != -1) {
                        pathEntry = taupPath.substring(offset, pathSepIndex);
                        offset = pathSepIndex + 1;
                    } else {
                        pathEntry = taupPath.substring(offset);
                        offset = taupPath.length();
                    }
                    /* Check each jar file. */
                    if (pathEntry.endsWith(".jar") || pathEntry.endsWith(".zip")) {
                        jarFile = new File(pathEntry);
                        if (jarFile.exists() && jarFile.isFile()
                                && jarFile.canRead()) {
                            ZipFile zippy = new ZipFile(jarFile);
                            ZipEntry zipEntry = zippy.getEntry("Models/" + filename);
                            if (zipEntry != null) {
                                return TauModel.readModelFromStream(zippy.getInputStream(zipEntry));
                            }
                        }
                    } else {
                        /* Check for regular files. */
                        modelFile = new File(pathEntry + "/" + filename);
                        if (modelFile.exists() && modelFile.isFile()
                                && modelFile.canRead()) {
                            return TauModel.readModel(modelFile.getCanonicalPath());
                        }
                    }
                }
            }
            /*
             * Couldn't find it in the taup.model.path either, look in the current
             * directory.
             */
            modelFile = new File(filename);
            if (modelFile.exists() && modelFile.isFile() && modelFile.canRead()) {
                return TauModel.readModel(modelFile.getCanonicalPath());
            }
            // try to load velocity model of same name and do a create
            try {
                VelocityModel vmod = loadVelocityModel(modelName);
                if (vmod != null) {
                    TauP_Create tauPCreate = new TauP_Create();
                    return tauPCreate.createTauModel(vmod);
                }
            } catch (Exception e) {
                throw new TauModelException("Can't find any saved models for "
                        + modelName + " and creation from velocity model failed.", e);
            }

            throw new TauModelException("Can't find any saved models for "
                    + modelName);

        } catch (InvalidClassException e) {
            throw new TauModelException("Unable to load '" + modelName + "'", e);
        } catch (StreamCorruptedException e) {
            throw new TauModelException("Unable to load '" + modelName + "'", e);
        } catch (OptionalDataException e) {
            throw new TauModelException("Unable to load '" + modelName + "'", e);
        } catch (FileNotFoundException e) {
            throw new TauModelException("Unable to load '" + modelName + "'", e);
        } catch (ClassNotFoundException e) {
            throw new TauModelException("Unable to load '" + modelName + "'", e);
        } catch (IOException e) {
            throw new TauModelException("Unable to load '" + modelName + "'", e);
        }
    }

    static VelocityModel loadVelocityModel(String modelName) throws IOException, VelocityModelException {
        /* First we try to find the model in the distributed taup.jar file. */
        try {
            Class c = Class.forName("edu.sc.seis.TauP.TauModelLoader");
            String filename = modelName + ".nd";
            InputStream in = c.getResourceAsStream(packageName + "/" + filename);
            if (in != null) {
                Reader inReader = new InputStreamReader(in);
                VelocityModel vmod = VelocityModel.readNDFile(inReader, modelName);
                inReader.close();
                return vmod;
            } else {
                // try tvel
                filename = modelName + ".tvel";
                in = c.getResourceAsStream(packageName + "/" + filename);
                if (in != null) {
                    Reader inReader = new InputStreamReader(in);
                    VelocityModel vmod = VelocityModel.readTVelFile(inReader, modelName);
                    inReader.close();
                    return vmod;
                }
            }
        } catch (Exception ex) {
            // couldn't get as a resource, so keep going
        }
        // try a .tvel or .nd file in current directory, or no suffix
        String[] types = new String[]{"", ".nd", ".tvel"};
        for (int i = 0; i < types.length; i++) {
            String vmodFile = modelName + types[i];
            File modelFile = new File(vmodFile);
            if (modelFile.exists() && modelFile.isFile() && modelFile.canRead()) {
                VelocityModel vMod = VelocityModel.readVelocityFile(modelFile.getPath(), types[i]);
                return vMod;
            }
        }
        return null;
    }

    protected static TauModel loadFromCache(String modelName) {
        SoftReference<TauModel> sr = tModCache.get(modelName);
        if (sr != null) {
            TauModel out = sr.get();
            if (out == null) {
                tModCache.remove(modelName);
            }
            return out;
        }
        return null;
    }

    static HashMap<String, SoftReference<TauModel>> tModCache = new HashMap<String, SoftReference<TauModel>>();
}
