#!/usr/bin/env python

import numpy as np
from scipy import interpolate
from scipy import linalg
import warnings

#radians to degrees conversion
RAD2DEG = 180./3.1415927;

#CONSTANTS FOR TESTING
ALEN = 16000.0
BLEN = 6100.0
CLEN = 10900.0
APLUNGE = 6.0
AAZIMUTH = 139.0
AROT = 88.9075

PRAX = np.zeros((3,3))
PRAX[0,0] = 139.0
PRAX[0,1] = 6.0
PRAX[0,2] = 16000 / 1000.
PRAX[1,0] = 310.0
PRAX[1,1] = 83.0
PRAX[1,2] = 6100 / 1000.
PRAX[2,0] = 49.0
PRAX[2,1] = 1.0
PRAX[2,2] = 10900 / 1000.

NDEF=247  # number of defining phases
STDERR=0.76
ISFIX=0

SFCMAJOR = 13.664099354621868
SFCMINOR = 9.3490986970307421
SFCAZIMUTH = 318.96374831744004

SFCMAJOR2 = 13.66054224573919
SFCMINOR2 = 9.3488176089599762
SFCAZIMUTH2 = 319.06909376893924

def vec2tait(ellipsoid):
    """Convert earthquake origin error ellipsoid 3x3 matrix into Tait-Bryan representation.

    Input argument:
    ellipsoid - 3x3 Numpy array containing the elements describing earthquake origin error ellipsoid.
    [SemiMajorAxisAzimuth SemiMajorAxisPlunge SemiMajorAxisLength;
     SemiMinorAxisAzimuth SemiMinorAxisPlunge SemiMinorAxisLength;
     SemiIntermediateAxisAzimuth SemiIntermediateAxisPlunge SemiIntermediateAxisLength]
    (distance values in kilometers, angles in degrees)
    
    Returns: 6 element tuple:
    semiMajorAxisLength (km)
    semiMinorAxisLength (km)
    semiIntermediateAxisLength (km)
    majorAxisAzimuth (degrees)
    majorAxisPlunge (degrees)
    majorAxisRotation (degrees)
    """
    ellipsearray = ellipsoid.flatten()

    #get indices of the minor, intermediate, major axes
    smallest,intermediate,largest = ellipsearray[2:9:3].argsort()

    # we have two of the angles already
    # convert from km to meters for the lengths
    #
    k=3*(largest)
    semiMajorAxisLength = ellipsearray[2+k]
    semiMinorAxisLength = ellipsearray[2 + 3*(smallest)]
    semiIntermediateAxisLength = ellipsearray[2 + 3*(intermediate)]
    majorAxisPlunge  = ellipsearray[1 + k]       # PHI
    majorAxisAzimuth = ellipsearray[0 + k]       # PSI
    PHI = majorAxisPlunge 
    PSI = majorAxisAzimuth 

    #####
    #    to get the THETA angle, we need to do some transformations
    #####
    RPHI=np.array([
        [np.cos(np.radians(PHI)), 0,  np.sin(np.radians(PHI))],
        [0, 1,  0],
        [-np.sin(np.radians(PHI)), 0, np.cos(np.radians(PHI))]])
    RPSI = np.array([
        [np.cos(np.radians(PSI)), np.sin(np.radians(PSI)), 0],
        [-np.sin(np.radians(PSI)),np.cos(np.radians(PSI)),0],
        [0 , 0, 1]])

    #####
    #    reconstruct the T matrix
    #####
    # major axis
    k=3*(largest)
    PHImaj=ellipsearray[1 + k]
    PSImaj=ellipsearray[0 + k]
    T = np.zeros((3,3))
    T[0,0] = np.cos(np.radians(PSImaj))*np.cos(np.radians(PHImaj))
    T[0,1] = np.sin(np.radians(PSImaj))*np.cos(np.radians(PHImaj))
    T[0,2] = np.sin(np.radians(PHImaj))

    # minor axis
    k=3*(smallest)
    PHImin=ellipsearray[1 + k]
    PSImin=ellipsearray[0 + k]
    T[1,0] = np.cos(np.radians(PSImin))*np.cos(np.radians(PHImin))
    T[1,1] = np.sin(np.radians(PSImin))*np.cos(np.radians(PHImin))
    T[1,2] = np.sin(np.radians(PHImin))

    # minor axis
    k=3*(intermediate)
    PHIint=ellipsearray[1 + k]
    PSIint=ellipsearray[0 + k]
    T[2,0] = np.cos(np.radians(PSIint))*np.cos(np.radians(PHIint))
    T[2,1] = np.sin(np.radians(PSIint))*np.cos(np.radians(PHIint))
    T[2,2] = np.sin(np.radians(PHIint))

    G=np.dot(T,np.dot(np.linalg.inv(RPSI),np.linalg.inv(RPHI)))
    THETA = np.arctan2(G[1,2], G[1,1])*RAD2DEG % 360 
    majorAxisRotation = THETA 
    
    return (semiMajorAxisLength,
            semiMinorAxisLength,
            semiIntermediateAxisLength,
            majorAxisAzimuth,
            majorAxisPlunge,
            majorAxisRotation)

def tait2vec(alen,blen,clen,azimuth,plunge,rotation):
    """Convert Tait-Bryan representation of earthquake origin error ellipsoid into 3x3 matrix representation.

    Input arguments:
    alen - semiMajorAxisLength (km)
    blen - semiMinorAxisLength (km)
    clen - semiIntermediateAxisLength (km)
    azimuth - majorAxisPlunge (degrees)
    plunge - majorAxisAzimuth (degrees)
    majorAxisRotation (degrees)
    
    Returns: 6 element tuple:
    ellipsoid - 3x3 Numpy array containing the elements describing earthquake origin error ellipsoid.
    [SemiMajorAxisAzimuth SemiMajorAxisPlunge SemiMajorAxisLength;
     SemiMinorAxisAzimuth SemiMinorAxisPlunge SemiMinorAxisLength;
     SemiIntermediateAxisAzimuth SemiIntermediateAxisPlunge SemiIntermediateAxisLength]
    (distance values in kilometers, angles in degrees)
    
    """
    
    # usage: T = tait2vec(PSI, PHI, THETA)
    #
    # return the 3x3 transformation matrix corresponding to the
    #      Tait-Bryan angles use in Section 3.3.9 ConfidenceEllipsoid
    #      QuakeML-BED-20130214a.pdf
    #     
    #      The image showing two views of the rotation in the QuakeML document is wrong
    #      and disagrees with the text, which says
    #      X is major axis, Y is the minor axis and thus Z is the intermediate axis
    #      The rotations are as follow:
    #      1. about z-axis by angle PSI (heading) to give (x', y', z)
    #      2, about y' with angle with angle PHI (elevation) to give (x'', y', z'')
    #      3. about x'' with angle THETA (bank) to give (x'', y''', z''')
    #
    #      This sequence is known as z-y'-x'' intrinsic rotation (http://en.wikipedia.org/wiki/Euler_angles)
    #      The figure in the QuakeML document is the z-x'-y'' rotation. Note the order
    #
    #      azimuth is measured positive in direction from x to y
    #      plunge is measured positive such that the x' moves in the positive z-direction
    #      rotation is measured positive such that the y'' moves in the positive z-direction
    #      
    #
    # azimuth is heading, measure in degrees from north
    # plunge is elevation
    # rotation is roll

    # Author: R.B.Herrmann (rbh@eas.slu.edu)
    # Created: 23 November 2014
    # Adapted-By: rbh
    # Translated to Python by Mike Hearne (mhearne@usgs.gov)
    alen *= 1000.0
    blen *= 1000.0
    clen *= 1000.0
    
    c1=np.cos(np.radians(azimuth))
    s1=np.sin(np.radians(azimuth))
    c2=np.cos(np.radians(plunge))
    s2=np.sin(np.radians(plunge))
    c3=np.cos(np.radians(rotation))
    s3=np.sin(np.radians(rotation))
    ##X1 Y2 Z3 rotation Table Tait-Bryan angles http://en.wikipedia.org/wiki/Euler_angles
    r11=c2*c1;
    r12=c2*s1;
    r13=s2;
    r21=-c3*s1-s3*s2*c1;
    r22=c3*c1-s3*s2*s1;
    r23=s3*c2;
    r31=s3*s1 - c3*s2*c1;
    r32=-s3*c1 - c3*s2*s1;
    r33=c3*c2;
    T= np.array([[ r11, r12, r13], [r21, r22, r23] , [r31, r32, r33 ]])

    # major axis
    prax = np.zeros((3,3))
    prax[0,0] = np.arctan2( T[0,1], T[0,0]) * RAD2DEG % 360
    prax[0,1] = np.degrees(np.arcsin(T[0,2]))
    prax[0,2] = alen / 1000.0 ;

    #minor axis
    prax[1,0] = np.arctan2(T[1,1], T[1,0]) * RAD2DEG % 360;
    prax[1,1] = np.degrees(np.arcsin(T[1,2]))
    prax[1,2] = blen/ 1000.0 ;

    #intermediate axis
    prax[2,0] = np.arctan2( T[2,1], T[2,0]) * RAD2DEG % 360
    prax[2,1] = np.degrees(np.arcsin(T[2,2]))
    prax[2,2] = clen / 1000.0 ;
    # arrange do that all dip angles are positive
    # since dip is measured from horizontal, changing dip
    # means chaninging the sign and then adding 180 degrees to the
    # azimuth
    for k in range(0,3):
       if ( prax[k,2] < 0 ):
            prax[k,1] = prax[k,0] +180 % 360
            prax[k,2] = - prax[k,1] ;
       
    return prax

def tait2surface(alen,blen,clen,aazimuth,aplunge,arot,ndef,stderr,isfix):
    """Project Tait-Bryan representation of error ellipse onto surface.
    Input parameters:
    alen - semiMajorAxisLength (km)
    blen - semiMinorAxisLength (km)
    clen - semiIntermediateAxisLength (km)
    azimuth - majorAxisPlunge (degrees)
    plunge - majorAxisAzimuth (degrees)
    arot - majorAxisRotation (degrees)
    ndef - Number of defining phases
    stderr - stderr of fit (??)
    isfix - True if depth was fixed in inversion, False if not
    
    Returns tuple of three values:
    Surface error ellipse major axis length
    Surface error ellipse minor axis length
    Surface error ellipse major axis azimuth
    """
    prax = tait2vec(alen,blen,clen,aazimuth,aplunge,arot)
    major,minor,azimuth = vec2surface(prax,ndef,stderr,isfix)
    return (major,minor,azimuth)
    
def vec2surface(prax,ndef,stderr,isfix):
    """Project 3x3 matrix representation of error ellipse onto surface.
    Input parameters:
    ellipsoid - 3x3 Numpy array containing the elements describing earthquake origin error ellipsoid.
    [SemiMajorAxisAzimuth SemiMajorAxisPlunge SemiMajorAxisLength;
     SemiMinorAxisAzimuth SemiMinorAxisPlunge SemiMinorAxisLength;
     SemiIntermediateAxisAzimuth SemiIntermediateAxisPlunge SemiIntermediateAxisLength]
     (distance values in kilometers, angles in degrees)
    ndef - Number of defining phases
    stderr - stderr of fit (??)
    isfix - True if depth was fixed in inversion, False if not
    
    Returns tuple of three values:
    Surface error ellipse major axis length
    Surface error ellipse minor axis length
    Surface error ellipse major axis azimuth
    """
    nFree = 8
    m = 3
    m1 = 2
    m2 = 4    
    xval = np.array([1,2,3,4,5,10,15,20,25,30,120,1000])
    f902 = np.array([49.50,9.00,5.46,4.32,3.78,2.92,2.70,2.59,2.53,2.49,2.35,2.30])
    f903 = np.array([53.59,9.16,5.39,4.19,3.62,2.73,2.49,2.38,2.32,2.20,2.13,2.08])
    func1 = interpolate.PchipInterpolator(xval,f903,extrapolate=True)
    fac3d = func1(nFree + ndef-4)
    func2 = interpolate.PchipInterpolator(xval,f902,extrapolate=True)
    fac2d = func2(nFree + ndef-4)
    s2 = ( nFree + ( ndef-m2)*stderr*stderr)/(nFree + ndef - m2)
    corr=(m1 *s2*fac2d)/(m *s2*fac3d)
    L = np.array([[np.power(prax[0,2],2), 0, 0 ],
                  [0, np.power(prax[1,2],2), 0 ],
                  [0, 0, np.power(prax[2,2],2) ]])

    X = np.zeros((3,3))
    for i in range(0,3):
        X[i,:] = np.array([np.cos(np.radians(prax[i,0]))*np.cos(np.radians(prax[i,1])) , 
                            np.sin(np.radians(prax[i,0]))*np.cos(np.radians(prax[i,1])), 
                            np.sin(np.radians(prax[i,1]))])
    C = np.dot(np.dot(X.T,L),X)
    C2x2 = np.array([[C[0,0],C[0,1]],
                     [C[1,0],C[1,1]]])
    [tmpEVAL,EVECT] = linalg.eig(C2x2)
    EVAL = np.zeros((2,2))
    EVAL[0,0] = np.real(tmpEVAL[0])
    EVAL[1,1] = np.real(tmpEVAL[1])
 
    if EVAL[0,0] > EVAL[1,1]:
        MajL = np.sqrt(np.abs(corr*EVAL[0,0]))
        warnings.simplefilter("error", RuntimeWarning)
        MinL = np.sqrt(np.abs(corr*EVAL[1,1]))
        Ang = np.arctan2(EVECT[1,0], EVECT[0,0])*RAD2DEG % 360
    else:
        MajL = np.sqrt(np.abs(corr*EVAL[1,1]))
        MinL = np.sqrt(np.abs(corr*EVAL[0,0]))
        Ang = np.arctan2(EVECT[1,1], EVECT[0,1])*RAD2DEG % 360

    return (MajL,MinL,Ang)

def test_vec2tait():
    alen,blen,clen,aazimuth,aplunge,arot = vec2tait(PRAX)

    np.testing.assert_almost_equal(alen,ALEN/1000.0,decimal=3)
    np.testing.assert_almost_equal(blen,BLEN/1000.0,decimal=3)
    np.testing.assert_almost_equal(clen,CLEN/1000.0,decimal=3)
    np.testing.assert_almost_equal(aplunge,APLUNGE,decimal=3)
    np.testing.assert_almost_equal(aazimuth,AAZIMUTH,decimal=3)
    np.testing.assert_almost_equal(arot,AROT,decimal=3)

def test_tait2vec():
    praxout = tait2vec(ALEN/1000.0,BLEN/1000.0,CLEN/1000.0,AAZIMUTH,APLUNGE,AROT)
    m,n = praxout.shape
    for i in range(0,m):
        for j in range(0,n):
            np.testing.assert_approx_equal(praxout[i][j],PRAX[i][j],significant=2)

def test_vec2surface():
    major,minor,azimuth = vec2surface(PRAX,NDEF,STDERR,ISFIX)
    np.testing.assert_almost_equal(major,SFCMAJOR,decimal=3)
    np.testing.assert_almost_equal(minor,SFCMINOR,decimal=3)
    np.testing.assert_almost_equal(azimuth,SFCAZIMUTH,decimal=3)
            
def test_tait2surface():
    major,minor,azimuth = tait2surface(ALEN/1000.0,BLEN/1000.0,CLEN/1000.0,AAZIMUTH,APLUNGE,AROT,NDEF,STDERR,ISFIX)
    np.testing.assert_almost_equal(major,SFCMAJOR2,decimal=3)
    np.testing.assert_almost_equal(minor,SFCMINOR2,decimal=3)
    np.testing.assert_almost_equal(azimuth,SFCAZIMUTH2,decimal=3)
    
    
if __name__ == '__main__':
    test_vec2tait()   
    test_tait2vec()
    test_vec2surface()
    test_tait2surface()
    alen, blen, clen, aazimuth, aplunge,arot = vec2tait(PRAX)

    print 'Major Axis Length: %f' % alen
    print 'Minor Axis Length: %f' % blen
    print 'Inter Axis Length: %f' % clen
    print 'Major Axis Plunge: %f' % aplunge
    print 'Major Axis Azimuth: %f' % aazimuth
    print 'Major Axis Rotation: %f' % arot

    np.set_printoptions(precision=1)
    prax2 = tait2vec(ALEN/1000.0,BLEN/1000.0,CLEN/1000.0,AAZIMUTH, APLUNGE, AROT)
    print 'INPUT ARRAY:'
    print PRAX
    print
    print 'OUTPUT ARRAY:'
    print prax2

    
