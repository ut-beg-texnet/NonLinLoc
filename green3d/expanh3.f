      subroutine expanh3 (x1,x2,dx1,dx2,x)

c-------------------------------------------------------------------------
c-     Programed by P.S.Lucio and G.C.Lambare in Dec/94                  -
c-     Last modification by P.S.Lucio in Oct/08/95                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c-------------------------------------------------------------------------

C*       THE CUBIC HERMITE INTERPOLATION BETWEEN RAYS X1 AND X2      *

      implicit none
      real x1, x2, dx1, dx2, x

c*********************************************************************
c   "u" must be between x1 and x2, in this way we take "u" equal the 
c   mid-point and "t", for the reparameterization, between 0 and 1.
c*********************************************************************

c      u = (x1+x2)/2.
c      t = (u-x1)/(x2-x1)
c      h0 = (1-t)**3+3*t*(1-t)**2
c      h1 = t*(1-t)**2
c      h2 = -t**2*(1-t)
c      h3 = 3*t**2*(1-t)+t**3
c      x = x1*h0 + dx1*h1 + dx2*h2 + x2*h3

c   In fact we take "t" equal 0.5, then the scheme is the following:

c        x = .125 * (4*(x1+x2) + (dx1-dx2))

      x = .5 * (x1+x2) + .125 * (dx1-dx2)

      return
      end
