/*
 * <pre> The TauP Toolkit: Flexible Seismic Travel-Time and Raypath Utilities.
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
 * owens@seis.sc.edu </pre>
 */
/**
 * ReflTransCoefficient.java Reflection and transmission coefficients for body
 * waves. Methods for calculating coefficients for each of the possible
 * interactions are provided. Calculations are done using the
 * com.visualnumerics.javagrande.Complex class from VisualNumerics. It is
 * further assume that the incoming ray is coming from the "top" for solid-solid
 * interactions and from the bottom for free surface interactions. If the ray is
 * actually coming from the bottom, the flip the velocities. The convention for
 * free surface and solid solid is a little strange, but can be thought of as
 * the top velocities correspond to the layer that they ray starts in.
 * 
 * @see "Aki and Richards page 144-151"
 * @see "Lay and Wallace page 98 "
 * @see <A HREF="http://math.nist.gov/javanumerics/">Java Numerics </A> Created:
 *      Wed Feb 17 12:25:27 1999
 * @author Philip Crotwell
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 */
package edu.sc.seis.TauP;

import java.io.Serializable;

public class ReflTransCoefficient implements Serializable {

    // IMPORTANT!!!!
    // Where ever "CX" appears in this class, it is used as a shorthand for
    // the Complex class, so CX.times() is the same as Complex.times, but
    // the code is, IMHO, less cluttered.
    /** just to avoid having Complex all over the place. */
    private static final Complex CX = new Complex();

    protected double topVp;

    protected double topVs;

    protected double topDensity;

    protected double botVp;

    protected double botVs;

    protected double botDensity;

    // "flat earth" ray parameter
    protected double rp;

    // temp variables to make calculations less ugly
    // first 3 lines follow both Aki and Richards and Lay and Wallace
    protected double a, b, c, d;

    protected Complex det, E, F, G, H;

    /** used only in free surface calculations */
    protected Complex fsA;

    /**
     * delta for SH-SH equations
     */
    protected Complex shDelta;

    // store the vertical slownesses for both the top and bottom halfspaces
    // for both P and S waves
    protected Complex topVertSlownessP, topVertSlownessS;

    protected Complex botVertSlownessP, botVertSlownessS;

    // we need the squared terms so often that it is worthwhile to store them
    protected double sqBotVs; // botVs squared

    protected double sqTopVs; // topVs squared

    protected double sqBotVp; // botVp squared

    protected double sqTopVp; // topVp squared

    protected double sqRP; // rp squared

    // remember last calculated ray param and wave type to avoid repeating
    protected double lastRayParam = -1.0;

    protected boolean lastInIsPWave = true;

    // CM
    protected boolean firstTime = true;

    // CM
    public ReflTransCoefficient(double topVp,
                                double topVs,
                                double topDensity,
                                double botVp,
                                double botVs,
                                double botDensity) {
        this.topVp = topVp;
        this.topVs = topVs;
        this.topDensity = topDensity;
        this.botVp = botVp;
        this.botVs = botVs;
        this.botDensity = botDensity;
    }

    /**
     * Flips the sense of the layers, useful when you have a ray going through
     * the same layer in the opposite direction.
     */
    public ReflTransCoefficient flip() {
        return new ReflTransCoefficient(botVp,
                                        botVs,
                                        botDensity,
                                        topVp,
                                        topVs,
                                        topDensity);
    }

    protected void calcTempVars(double rayParam, boolean inIsPWave) {
        if(rayParam < 0) {
            throw new IllegalArgumentException("rayParam cannot be negative");
        }
        this.rp = rayParam; // ray parameter
        // CM
        // if (rayParam != lastRayParam && inIsPWave == lastInIsPWave ) {
        // if ( (rayParam != lastRayParam || inIsPWave != lastInIsPWave ) ||
        // firstTime ) {
        if(rayParam != lastRayParam || inIsPWave != lastInIsPWave) {
            lastRayParam = -1.0; // in case of failure in method
            // CM
            firstTime = false;
            sqBotVs = botVs * botVs; // botVs squared
            sqTopVs = topVs * topVs; // topVs squared
            sqBotVp = botVp * botVp; // botVp squared
            sqTopVp = topVp * topVp; // topVp squared
            sqRP = rp * rp; // rp squared
            topVertSlownessP = Complex.sqrt(new Complex(1.0 / sqTopVp - sqRP));
            topVertSlownessS = Complex.sqrt(new Complex(1.0 / sqTopVs - sqRP));
            botVertSlownessP = Complex.sqrt(new Complex(1.0 / sqBotVp - sqRP));
            botVertSlownessS = Complex.sqrt(new Complex(1.0 / sqBotVs - sqRP));
            a = botDensity * (1.0 - 2 * sqBotVs * sqRP) - topDensity
                    * (1.0 - 2 * sqTopVs * sqRP);
            b = botDensity * (1.0 - 2 * sqBotVs * sqRP) + 2 * topDensity
                    * sqTopVs * sqRP;
            c = topDensity * (1.0 - 2 * sqTopVs * sqRP) + 2 * botDensity
                    * sqBotVs * sqRP;
            d = 2 * (botDensity * sqBotVs - topDensity * sqTopVs);
            // math with complex objects is hard to read, so we give
            // the formulas as comments
            // E = b * topVertSlownessP + c * botVertSlownessP
            // CM E = b * cos(i1)/alpha1 + c * cos(i2)/alpha2
            E = CX.plus(topVertSlownessP.times(b), botVertSlownessP.times(c));
            // F = b * topVertSlownessS + c * botVertSlownessS
            F = CX.plus(topVertSlownessS.times(b), botVertSlownessS.times(c));
            // G = a - d * topVertSlownessP * botVertSlownessS
            G = CX.minus(new Complex(a),
                         CX.times(d, CX.times(topVertSlownessP,
                                              botVertSlownessS)));
            // H = a - d * botVertSlownessP * topVertSlownessS
            H = CX.minus(new Complex(a),
                         CX.times(d, CX.times(botVertSlownessP,
                                              topVertSlownessS)));
            // det = E * F + G * H * sqRP
            det = CX.plus(CX.times(E, F), CX.times(G, H).times(sqRP));
            // free surface denominator
            // fsA = ((1/sqBotVs) - 2 * sqRP)^2 +
            // 4 * sqRP * botVertSlownessP * botVertSlownessS
            fsA = CX.plus(new Complex(((1 / sqTopVs) - 2 * sqRP)
                    * ((1 / sqTopVs) - 2 * sqRP)), CX.times(topVertSlownessP,
                                                            topVertSlownessS)
                    .times(4 * sqRP));
            // SH delta
            shDelta = CX.plus(CX.plus(topDensity * topVs * topVs,
                                      topVertSlownessS), CX.plus(botDensity
                    * botVs * botVs, botVertSlownessS));
            lastRayParam = rayParam;
            lastInIsPWave = inIsPWave;
        }
    }

    // FREE SURFACE
    /**
     * Calculates incident P wave to reflected P wave complex coefficient at
     * free surface. Only topVp, topVs, and topDensity are used, the bottom
     * values are ignored. This is a little strange as free surface rays are
     * always upgoing, but it mantains consistency with the solid-solid
     * calculations.
     * <P>= (-1*((1/sqTopVs) - 2 * sqRP)^2 +<BR>
     * 4 * sqRP * topVertSlownessP * topVertSlownessS) / A
     */
    public Complex getComplexFreePtoPRefl(double rayParam) {
        calcTempVars(rayParam, true);
        Complex numerator = CX.plus(-1.0 * ((1 / sqTopVs) - 2 * sqRP)
                * ((1 / sqTopVs) - 2 * sqRP), CX.times(topVertSlownessP,
                                                       topVertSlownessS)
                .times(4 * sqRP));
        return CX.over(numerator, fsA);
    }

    /**
     * Calculates incident P wave to reflected P wave coefficient at free
     * surface. This just returns the real part, assuming the imaginary part is
     * zero.
     * 
     * @see #getComplexFreePtoPRefl(double)
     */
    public double getFreePtoPRefl(double rayParam) {
        return CX.real(getComplexFreePtoPRefl(rayParam));
    }

    /**
     * Calculates incident P wave to reflected SV wave complex coefficient at
     * free surface. = (4 * (topVp/topVs) * rp * topVertSlownessP * ((1/sqTopVs) -
     * 2 * sqRP)) / fsA
     */
    public Complex getComplexFreePtoSVRefl(double rayParam) {
        calcTempVars(rayParam, true);
        double realNumerator = 4 * (topVp / topVs) * rp
                * ((1 / sqTopVs) - 2 * sqRP);
        return CX.over(CX.times(realNumerator, topVertSlownessP), fsA);
    }

    /**
     * Calculates incident P wave to reflected SV wave coefficient at free
     * surface.
     * 
     * @see #getComplexFreePtoSVRefl(double)
     */
    public double getFreePtoSVRefl(double rayParam) {
        return CX.real(getComplexFreePtoSVRefl(rayParam));
    }

    /**
     * Calculates incident SV wave to reflected P wave complex coefficient at
     * free surface.
     * <P>= (4 * (topVs/topVp) * rp * topVertSlownessS *<BR>
     * ((1/sqTopVs) - 2 * sqRP)) / fsA
     */
    public Complex getComplexFreeSVtoPRefl(double rayParam) {
        calcTempVars(rayParam, false);
        double realNumerator = 4 * (topVs / topVp) * rp
                * ((1 / sqTopVs) - 2 * sqRP);
        return CX.over(CX.times(realNumerator, topVertSlownessS), fsA);
    }

    /**
     * Calculates incident SV wave to reflected P wave coefficient at free
     * surface.
     * 
     * @see #getComplexFreeSVtoPRefl(double)
     */
    public double getFreeSVtoPRefl(double rayParam) {
        return CX.real(getComplexFreeSVtoPRefl(rayParam));
    }

    /**
     * Calculates incident SV wave to reflected SV wave complex coefficient at
     * free surface.
     * <P>= (-1 * ((1/sqTopVs) - 2 * sqRP)^2 +<BR>
     * 4 * sqRP * topVertSlownessP * topVertSlownessS) / fsA
     */
    public Complex getComplexFreeSVtoSVRefl(double rayParam) {
        calcTempVars(rayParam, false);
        // Aki and Richards don't have -1
        double realNumerator = ((1 / sqTopVs) - 2 * sqRP)
                * ((1 / sqTopVs) - 2 * sqRP);
        Complex numerator = CX.plus(realNumerator,
                                    CX.times(4 * sqRP,
                                             CX.times(topVertSlownessP,
                                                      topVertSlownessS)));
        return CX.over(numerator, fsA);
    }

    /**
     * Calculates incident SV wave to reflected SV wave coefficient at free
     * surface.
     */
    public double getFreeSVtoSVRefl(double rayParam) {
        return CX.real(getComplexFreeSVtoSVRefl(rayParam));
    }

    /**
     * Calculates incident SH wave to reflected SH wave complex coefficient at
     * free surface. Free surface SH is always 1.
     */
    public Complex getComplexFreeSHtoSHRefl(double rayParam) {
        return new Complex(1);
    }

    /**
     * Calculates incident SH wave to reflected SH wave coefficient at free
     * surface. Free surface SH is always 1.
     */
    public double getFreeSHtoSHRefl(double rayParam) {
        return 1;
    }

    // Solid-Solid interface
    /**
     * Calculates incident P wave to reflected P wave Complex coefficient.
     * <P>= ((b*topVertSlownessP - c*botVertSlownessP)*F -<BR>
     * (a + d*topVertSlownessP * botVertSlownessS)*H*sqRP) / det
     */
    public Complex getComplexPtoPRefl(double rayParam) {
        calcTempVars(rayParam, true);
        Complex FTerm = CX.times(CX.minus(CX.times(b, topVertSlownessP),
                                          CX.times(c, botVertSlownessP)), F);
        Complex HTerm = CX.times(CX.plus(a,
                                         CX.times(d, CX.times(topVertSlownessP,
                                                              botVertSlownessS))),
                                 CX.times(H, sqRP));
        return CX.over(CX.minus(FTerm, HTerm), det);
    }

    /**
     * Calculates incident P wave to reflected P wave coefficient.
     */
    public double getPtoPRefl(double rayParam) {
        // return CX.abs(getComplexPtoPRefl(rayParam));
        return CX.real(getComplexPtoPRefl(rayParam));
    }

    /**
     * Calculates incident P wave to reflected SV wave Complex coefficient.
     * <P>= -2 * topVertSlownessP *<BR>
     * (a * b + c * d *botVertSlownessP *botVertSlownessS) *<BR>
     * rp * (topVp/topVs)) / det
     */
    public Complex getComplexPtoSVRefl(double rayParam) {
        calcTempVars(rayParam, true);
        double realNumerator = -2 * rp * (topVp / topVs);
        Complex middleTerm = CX.plus(a * b,
                                     CX.times(c * d, CX.times(botVertSlownessP,
                                                              botVertSlownessS)));
        Complex numerator = CX.times(CX.times(realNumerator, topVertSlownessP),
                                     middleTerm);
        return CX.over(numerator, det);
    }

    /**
     * Calculates incident P wave to reflected SV wave coefficient.
     */
    public double getPtoSVRefl(double rayParam) {
        return CX.real(getComplexPtoSVRefl(rayParam));
    }

    /**
     * Calculates incident P wave to transmitted P wave Complex coefficient.
     * <P>= ( 2 * topDensity * topVertSlownessP * F *<BR>
     * (topVp / botVp)) / det
     */
    public Complex getComplexPtoPTrans(double rayParam) {
        calcTempVars(rayParam, true);
        double realNumerator = 2 * topDensity * (topVp / botVp);
        return CX.over(CX.times(realNumerator, CX.times(topVertSlownessP, F)),
                       det);
    }

    /**
     * Calculates incident P wave to transmitted P wave coefficient.
     */
    public double getPtoPTrans(double rayParam) {
        // return CX.abs(getComplexPtoPTrans(rayParam));
        return CX.real(getComplexPtoPTrans(rayParam));
    }

    /**
     * Calculates incident P wave to transmitted SV wave Complex coefficient.
     * <P>= (2 * topDensity * topVertSlownessP * H * rp * (topVp / botVs)) /
     * <BR>
     * det
     */
    public Complex getComplexPtoSVTrans(double rayParam) {
        calcTempVars(rayParam, true);
        double realNumerator = 2 * topDensity * rp * (topVp / botVs);
        Complex numerator = CX.times(realNumerator, CX.times(topVertSlownessP,
                                                             H));
        return CX.over(numerator, det);
    }

    /**
     * Calculates incident P wave to transmitted SV wave coefficient.
     */
    public double getPtoSVTrans(double rayParam) {
        // return CX.abs(getComplexPtoSVTrans(rayParam));
        return CX.real(getComplexPtoSVTrans(rayParam));
    }

    /**
     * Calculates incident SV wave to reflected P wave Complex coefficient.
     * <P>= (-2 * topVertSlownessS *<BR>
     * (a * b + c * d * botVertSlownessP * botVertSlownessS) *<BR>
     * rp * (topVs / topVp)) /<BR>
     * det
     */
    public Complex getComplexSVtoPRefl(double rayParam) {
        calcTempVars(rayParam, false);
        double realNumerator = -2 * rp * (topVs / topVp);
        // double realNumerator = -2 * rp * (topVs / topVp);
        Complex middleTerm = CX.plus(a * b,
                                     CX.times(c * d, CX.times(botVertSlownessP,
                                                              botVertSlownessS)));
        Complex numerator = CX.times(realNumerator, CX.times(topVertSlownessS,
                                                             middleTerm));
        return CX.over(numerator, det);
    }

    /**
     * Calculates incident SV wave to reflected P wave coefficient.
     */
    public double getSVtoPRefl(double rayParam) {
        // return CX.abs(getComplexSVtoPRefl(rayParam));
        return CX.real(getComplexSVtoPRefl(rayParam));
    }

    /**
     * Calculates incident SV wave to reflected SV wave Complex coefficient.
     * <P>= -1 * ((b * topVertSlownessS - c * botVertSlownessS) * E -<BR>
     * (a + b * botVertSlownessP * topVertSlownessS) * G * sqRP) /<BR>
     * det
     */
    public Complex getComplexSVtoSVRefl(double rayParam) {
        calcTempVars(rayParam, false);
        Complex adNumerator = CX.times(CX.plus(a,
                                               CX.times(d,
                                                        CX.times(botVertSlownessP,
                                                                 topVertSlownessS))),
                                       CX.times(G, sqRP));
        Complex bcNumerator = CX.times(CX.minus(CX.times(b, topVertSlownessS),
                                                CX.times(c, botVertSlownessS)),
                                       E);
        return CX.over(CX.minus(adNumerator, bcNumerator), det);
    }

    /**
     * Calculates incident SV wave to reflected SV wave coefficient.
     */
    public double getSVtoSVRefl(double rayParam) {
        // return CX.abs(getComplexSVtoSVRefl(rayParam));
        return CX.real(getComplexSVtoSVRefl(rayParam));
    }

    /**
     * Calculates incident SV wave to transmitted P wave Complex coefficient.
     * <P>= -2 * topDensity * topVertSlownessS * G * rp * (topVs / botVp) /
     * <BR>
     * det
     */
    public Complex getComplexSVtoPTrans(double rayParam) {
        calcTempVars(rayParam, false);
        double realNumerator = -2 * topDensity * rp * (topVs / botVp);
        Complex numerator = CX.times(realNumerator, CX.times(topVertSlownessS,
                                                             G));
        return CX.over(numerator, det);
    }

    /**
     * Calculates incident SV wave to transmitted P wave coefficient.
     */
    public double getSVtoPTrans(double rayParam) {
        // return CX.abs(getComplexSVtoPTrans(rayParam));
        return CX.real(getComplexSVtoPTrans(rayParam));
    }

    /**
     * Calculates incident SV wave to transmitted SV wave Complex coefficient.
     * <P>= 2 * topDensity * topVertSlownessS * E * (topVs / botVs) /<BR>
     * det
     */
    public Complex getComplexSVtoSVTrans(double rayParam) {
        calcTempVars(rayParam, false);
        double realNumerator = 2 * topDensity * rp * (topVs / botVs);
        Complex numerator = CX.times(realNumerator, CX.times(topVertSlownessS,
                                                             E));
        return CX.over(numerator, det);
    }

    /**
     * Calculates incident SV wave to transmitted SV wave coefficient.
     */
    public double getSVtoSVTrans(double rayParam) {
        // return CX.abs(getComplexSVtoSVTrans(rayParam));
        return CX.real(getComplexSVtoSVTrans(rayParam));
    }

    // SH waves
    /**
     * Calculates incident SH wave to reflected SH wave Complex coefficient.
     * <P>
     * mu = Vs * Vs * density
     * <P>= (topMu * topVertSlownessS - botMu * botVertSlownessS) /<BR>
     * (topMu * topVertSlownessS + botMu * botVertSlownessS)
     */
    public Complex getComplexSHtoSHRefl(double rayParam) {
        calcTempVars(rayParam, false);
        double topMu = topVs * topVs * topDensity;
        double botMu = botVs * botVs * botDensity;
        Complex topTerm = CX.times(topMu, topVertSlownessS);
        Complex botTerm = CX.times(botMu, botVertSlownessS);
        return CX.over(CX.minus(topTerm, botTerm), CX.plus(topTerm, botTerm));
    }

    /**
     * Calculates incident SH wave to reflected SH wave coefficient.
     */
    public double getSHtoSHRefl(double rayParam) {
        // return CX.abs(getComplexSHtoSHRefl(rayParam));
        return CX.real(getComplexSHtoSHRefl(rayParam));
    }

    /**
     * Calculates incident SH wave to transmitted SH wave Complex coefficient.
     * <P>
     * mu = Vs * Vs * density
     * <P>= 2 * topMu * topVertSlownessS /<BR>
     * (topMu * topVertSlownessS + botMu * botVertSlownessS)
     */
    public Complex getComplexSHtoSHTrans(double rayParam) {
        calcTempVars(rayParam, false);
        double topMu = topVs * topVs * topDensity;
        double botMu = botVs * botVs * botDensity;
        Complex topTerm = CX.times(topMu, topVertSlownessS);
        Complex botTerm = CX.times(botMu, botVertSlownessS);
        return CX.over(CX.times(topTerm, 2), CX.plus(topTerm, botTerm));
    }

    /**
     * Calculates incident SH wave to transmitted SH wave coefficient.
     */
    public double getSHtoSHTrans(double rayParam) {
        // return CX.abs(getComplexSHtoSHTrans(rayParam));
        return CX.real(getComplexSHtoSHTrans(rayParam));
    }

    public static void main(String[] args) {
        double topVp = 4.98;
        double topVs = 2.9;
        double topDensity = 2.667;
        double botVp = 8.0;
        double botVs = 4.6;
        double botDensity = 3.38;
        double depth;
        double radiusOfEarth;
        double DtoR = Math.PI / 180.0;
        double RtoD = 180.0 / Math.PI;
        double[] RPP = new double[91];
        double[] RPS = new double[91];
        double[] RSP = new double[91];
        double[] RSS = new double[91];
        double[] TPP = new double[91];
        double[] TPS = new double[91];
        double[] TSP = new double[91];
        double[] TSS = new double[91];
        ReflTransCoefficient coeff = new ReflTransCoefficient(topVp,
                                                              topVs,
                                                              topDensity,
                                                              botVp,
                                                              botVs,
                                                              botDensity);
        double rayParam;
        for(int i = 0; i <= 90; i++) {
            rayParam = Math.sin(DtoR * i) / topVp;
            RPP[i] = coeff.getPtoPRefl(rayParam);
            RPS[i] = coeff.getPtoSVRefl(rayParam);
            TPP[i] = coeff.getPtoPTrans(rayParam);
            TPS[i] = coeff.getPtoSVTrans(rayParam);
            rayParam = Math.sin(DtoR * i) / topVs;
            RSP[i] = coeff.getSVtoPRefl(rayParam);
            RSS[i] = coeff.getSVtoSVRefl(rayParam);
            TSP[i] = coeff.getSVtoPTrans(rayParam);
            TSS[i] = coeff.getSVtoSVTrans(rayParam);
        }
        try {
            java.io.Writer out = new java.io.BufferedWriter(new java.io.FileWriter("refltrans.gmt"));
            out.write("#!/bin/sh\n\n");
            out.write("/bin/rm -f refltrans.ps\n\n");
            out.write("psbasemap -K -P -JX6 -R0/90/0/2  -B10/1 > refltrans.ps\n");
            // out.write("psxy -K -O -JX -R -M -W1/255/0/0 >> refltrans.ps
            // <<END\n");
            // for (int i=0; i<=90; i++) {
            // out.write(i+" "+RPP[i]+"\n");
            // }
            // out.write("END\n");
            // out.write("psxy -K -O -JX -R -M -W1/0/255/0 >> refltrans.ps
            // <<END\n");
            // for (int i=0; i<=90; i++) {
            // out.write(i+" "+RPS[i]+"\n");
            // }
            // out.write("END\n");
            // out.write("psxy -K -O -JX -R -M -W1/0/0/255 >> refltrans.ps
            // <<END\n");
            // for (int i=0; i<=90; i++) {
            // out.write(i+" "+TPP[i]+"\n");
            // }
            // out.write("END\n");
            // out.write("psxy -K -O -JX -R -M -W1/255/255/0 >> refltrans.ps
            // <<END\n");
            // for (int i=0; i<=90; i++) {
            // out.write(i+" "+TPS[i]+"\n");
            // }
            // out.write("END\n");
            out.write("psxy -K -O -JX -R -M -W1/255/0/0 >> refltrans.ps <<END\n");
            for(int i = 0; i <= 90; i++) {
                out.write(i + " " + RSP[i] + "\n");
            }
            out.write("END\n");
            out.write("psxy -K -O -JX -R -M -W1/0/255/0 >> refltrans.ps <<END\n");
            for(int i = 0; i <= 90; i++) {
                out.write(i + " " + RSS[i] + "\n");
            }
            out.write("END\n");
            out.write("psxy -K -O -JX -R -M -W1/0/0/255 >> refltrans.ps <<END\n");
            for(int i = 0; i <= 90; i++) {
                out.write(i + " " + TSP[i] + "\n");
            }
            out.write("END\n");
            out.write("psxy -O -JX -R -M -W1/255/0/255 >> refltrans.ps <<END\n");
            for(int i = 0; i <= 90; i++) {
                out.write(i + " " + TSS[i] + "\n");
            }
            out.write("END\n");
            out.close();
        } catch(java.io.IOException e) {
            System.err.println(e);
        }
    }

    public Complex getBotVertSlownessP(double rayParam) {
        calcTempVars(rayParam, true);
        return botVertSlownessP;
    }

    public Complex getBotVertSlownessS(double rayParam) {
        calcTempVars(rayParam, false);
        return botVertSlownessS;
    }

    public Complex getTopVertSlownessP(double rayParam) {
        calcTempVars(rayParam, true);
        return topVertSlownessP;
    }

    public Complex getTopVertSlownessS(double rayParam) {
        calcTempVars(rayParam, false);
        return topVertSlownessS;
    }
} // ReflTransCoefficient
