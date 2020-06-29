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
/**
 * package for storage and manipulation of seismic earth models.
 * 
 */
package edu.sc.seis.TauP;

import java.io.Serializable;

/**
 * The VelocityModelLayer class stores and manipulates a singly layer. An entire
 * velocity model is implemented as an Vector of layers.
 * 
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * 
 * 
 * 
 * @author H. Philip Crotwell
 */
public class VelocityLayer implements Cloneable, Serializable {

    private int myLayerNumber;

    private double topDepth;

    private double botDepth;

    private double topPVelocity;

    private double botPVelocity;

    private double topSVelocity;

    private double botSVelocity;

    private double topDensity = 2.6;

    private double botDensity = 2.6;

    private double topQp = 1000;

    private double botQp = 1000;

    private double topQs = 2000;

    private double botQs = 2000;

    public VelocityLayer(int myLayerNumber,
                         double topDepth,
                         double botDepth,
                         double topPVelocity,
                         double botPVelocity,
                         double topSVelocity,
                         double botSVelocity) {
        this(myLayerNumber,
             topDepth,
             botDepth,
             topPVelocity,
             botPVelocity,
             topSVelocity,
             botSVelocity,
             2.6,
             2.6);
    }

    public VelocityLayer(int myLayerNumber,
                         double topDepth,
                         double botDepth,
                         double topPVelocity,
                         double botPVelocity,
                         double topSVelocity,
                         double botSVelocity,
                         double topDensity,
                         double bottomDensity) {
        this(myLayerNumber,
             topDepth,
             botDepth,
             topPVelocity,
             botPVelocity,
             topSVelocity,
             botSVelocity,
             topDensity,
             bottomDensity,
             1000,
             1000,
             2000,
             2000);
    }

    public VelocityLayer(int myLayerNumber,
                         double topDepth,
                         double botDepth,
                         double topPVelocity,
                         double botPVelocity,
                         double topSVelocity,
                         double botSVelocity,
                         double topDensity,
                         double botDensity,
                         double topQp,
                         double botQp,
                         double topQs,
                         double botQs) {
        if(topPVelocity <= 0) {
            throw new IllegalArgumentException("topPVelocity must be positive: "
                    + topPVelocity);
        }
        if(botPVelocity <= 0) {
            throw new IllegalArgumentException("botPVelocity must be positive: "
                    + botPVelocity);
        }
        if(topSVelocity < 0) {
            throw new IllegalArgumentException("topSVelocity must be nonnegative: "
                    + topSVelocity);
        }
        if(botSVelocity < 0) {
            throw new IllegalArgumentException("botSVelocity must be nonnegative: "
                    + botSVelocity);
        }
        this.myLayerNumber = myLayerNumber;
        this.topDepth = topDepth;
        this.botDepth = botDepth;
        this.topPVelocity = topPVelocity;
        this.botPVelocity = botPVelocity;
        this.topSVelocity = topSVelocity;
        this.botSVelocity = botSVelocity;
        this.topDensity = topDensity;
        this.botDensity = botDensity;
        this.topQp = topQp;
        this.botQp = botQp;
        this.topQs = topQs;
        this.botQs = botQs;
    }

    public Object clone() {
        try {
            VelocityLayer newObject = (VelocityLayer)super.clone();
            return newObject;
        } catch(CloneNotSupportedException e) {
            // Cannot happen, we support clone
            // and our parent is Object, which supports clone.
            throw new InternalError(e.toString());
        }
    }

    public double evaluateAtBottom(char materialProperty)
            throws NoSuchMatPropException {
        double answer;
        switch(materialProperty){
            case 'P':
            case 'p':
                answer = getBotPVelocity();
                break;
            case 's':
            case 'S':
                answer = getBotSVelocity();
                break;
            case 'r':
            case 'R':
            case 'D':
            case 'd':
                answer = getBotDensity();
                break;
            default:
                throw new NoSuchMatPropException(materialProperty);
        }
        return answer;
    }

    public double evaluateAtTop(char materialProperty)
            throws NoSuchMatPropException {
        double answer;
        switch(materialProperty){
            case 'P':
            case 'p':
                answer = getTopPVelocity();
                break;
            case 's':
            case 'S':
                answer = getTopSVelocity();
                break;
            case 'r':
            case 'R':
            case 'D':
            case 'd':
                answer = getTopDensity();
                break;
            default:
                throw new NoSuchMatPropException(materialProperty);
        }
        return answer;
    }

    public double evaluateAt(double depth, char materialProperty)
            throws NoSuchMatPropException {
        double slope, answer;
        switch(materialProperty){
            case 'P':
            case 'p':
                slope = (getBotPVelocity() - getTopPVelocity())
                        / (getBotDepth() - getTopDepth());
                answer = slope * (depth - getTopDepth()) + getTopPVelocity();
                break;
            case 's':
            case 'S':
                slope = (getBotSVelocity() - getTopSVelocity())
                        / (getBotDepth() - getTopDepth());
                answer = slope * (depth - getTopDepth()) + getTopSVelocity();
                break;
            case 'r':
            case 'R':
            case 'D':
            case 'd':
                slope = (getBotDensity() - getTopDensity())
                        / (getBotDepth() - getTopDepth());
                answer = slope * (depth - getTopDepth()) + getTopDensity();
                break;
            default:
                System.out.println("I don't understand this material property: "
                        + materialProperty + "\nUse one of P p S s R r D d");
                throw new NoSuchMatPropException(materialProperty);
        }
        return answer;
    }

    public String toString() {
        String description;
        description = myLayerNumber + " " + getTopDepth() + " " + getBotDepth();
        description += " P " + getTopPVelocity() + " " + getBotPVelocity();
        description += " S " + getTopSVelocity() + " " + getBotSVelocity();
        description += " Density " + getTopDensity() + " " + getBotDensity();
        return description;
    }
    
    public int getLayerNum() {
        return myLayerNumber;
    }

    public void setTopDepth(double topDepth) {
        this.topDepth = topDepth;
    }

    public double getTopDepth() {
        return topDepth;
    }

    public void setBotDepth(double botDepth) {
        this.botDepth = botDepth;
    }

    public double getBotDepth() {
        return botDepth;
    }

    public void setTopPVelocity(double topPVelocity) {
        this.topPVelocity = topPVelocity;
    }

    public double getTopPVelocity() {
        return topPVelocity;
    }

    public void setBotPVelocity(double botPVelocity) {
        this.botPVelocity = botPVelocity;
    }

    public double getBotPVelocity() {
        return botPVelocity;
    }

    public void setTopSVelocity(double topSVelocity) {
        this.topSVelocity = topSVelocity;
    }

    public double getTopSVelocity() {
        return topSVelocity;
    }

    public void setBotSVelocity(double botSVelocity) {
        this.botSVelocity = botSVelocity;
    }

    public double getBotSVelocity() {
        return botSVelocity;
    }

    public void setTopDensity(double topDensity) {
        this.topDensity = topDensity;
    }

    public double getTopDensity() {
        return topDensity;
    }

    public void setBotDensity(double botDensity) {
        this.botDensity = botDensity;
    }

    public double getBotDensity() {
        return botDensity;
    }

    public void setTopQp(double topQp) {
        this.topQp = topQp;
    }

    public double getTopQp() {
        return topQp;
    }

    public void setBotQp(double botQp) {
        this.botQp = botQp;
    }

    public double getBotQp() {
        return botQp;
    }

    public void setTopQs(double topQs) {
        this.topQs = topQs;
    }

    public double getTopQs() {
        return topQs;
    }

    public void setBotQs(double botQs) {
        this.botQs = botQs;
    }

    public double getBotQs() {
        return botQs;
    }

    public double getThickness() {
        return getBotDepth()-getTopDepth();
    }
}
