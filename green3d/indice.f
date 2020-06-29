      subroutine indice(x1,x2,x3,x4,
     &              nxr,nyr,nzr,
     &              icmin,icmax,jcmin,jcmax,kcmin,kcmax)
      implicit none
      integer icmin,icmax,jcmin,jcmax,kcmin,kcmax
      integer nxr,nyr,nzr
      real    xcmin,xcmax,ycmin,ycmax,zcmin,zcmax

c-------------------------------------------------------------------------
c-     Programed by P.S.Lucio and G.C.Lambare in Dec/94                  -
c-     Last modification by P.S.Lucio in Oct/20/95                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c-------------------------------------------------------------------------

c************************************************************************
c*            DETERMINATION OF A LOCAL NEIGHBOROOD OF THE CELL          *
c************************************************************************
c--
      real x1(3), x2(3), x3(3), x4(3)
c--
      xcmin=amin1(x1(1),x2(1),x3(1),x4(1))
      xcmax=amax1(x1(1),x2(1),x3(1),x4(1))
      ycmin=amin1(x1(2),x2(2),x3(2),x4(2))
      ycmax=amax1(x1(2),x2(2),x3(2),x4(2))
      zcmin=amin1(x1(3),x2(3),x3(3),x4(3))
      zcmax=amax1(x1(3),x2(3),x3(3),x4(3))

      icmin=max(  1,nint(-.5+xcmin+2))
      icmax=min(nxr,nint(-.5+xcmax+1))
      jcmin=max(  1,nint(-.5+ycmin+2))
      jcmax=min(nyr,nint(-.5+ycmax+1))
      kcmin=max(  1,nint(-.5+zcmin+2))
      kcmax=min(nzr,nint(-.5+zcmax+1))

      return
      end
