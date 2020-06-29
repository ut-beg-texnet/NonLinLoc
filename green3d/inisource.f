      subroutine inisource (ip,xsrc,usrc,
     &        vm,xm,ym,zm,nxm,nym,nzm,dxm,dym,dzm,
     &        fi1min,fi2min,fi1max,fi2max,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)

c-------------------------------------------------------------------------
c-     Programed by P.S.Lucio and G.C.Lambare in Nov/94                  -
c-     Last modification by P.S.Lucio in Oct/20/95                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c-------------------------------------------------------------------------

C        **********************************************************
C        *   PROGRAM OF INITIALISATION OF THE FIRST WAVEFRONT     *
C        **********************************************************

c*******************************************************************
c INPUT  : ip - the wavefront
c
c OUTPUT : nrai - number of the "rays" at the wavefront
c          nseg - number of the segments at the wavefront
c          ncel - number of the cells at the wavefront
c          iseg - points of each segment
c          icel - points of each cell
c*******************************************************************
      implicit none
      integer ip
 
      integer npn,nsn,ncn,nfront,nd,nk,npar   
      real ys(npar,npn,nfront)
      integer iys(nk,npn,nfront)
      double precision dys(nd,npn,nfront)
      integer  iseg(nsn,6,nfront),icel(ncn,6,nfront)
      integer nrai(nfront), ncel(nfront), nseg(nfront)
      integer irai(3), iraf(3), irpop(3), irais(2,3)
      integer isatr(nfront), isats(nfront), isatc(nfront)

c--
      integer nxm,nym,nzm
      real xm(2),ym(2),zm(2)
      real vm(nzm,nxm,nym)
      real dxm,dym,dzm
c--
      real xsrc(3),usrc
      real fi1min,fi2min,fi1max,fi2max,pi,pi2
      integer i,j,ipoint
C***********************************************************************
C* The choice of a initial wavefront is not straightforward, it should *
C* be a sphere in the 3D slowness space (it's not possible to define a *
C*  regular parameterization of this wavefront). As our applications   *
C*   concern seismic reflexion data we've decided to take simply two   *
C*    triangles for the initial wavefront defined by a 2D aperture.    *
C***********************************************************************
c **********************************************************************
c *                                                                    *
c *                    1           (4)           4                     *
c *                     #------------------------#                     *
c *                     |  \__                   |                     *
c *                     |     \____  (3)    II   |                     *
c *                  (1)|          \_____        |(5)                  *
c *               |     |                \____   |                     *
c *               |     |   I                 \_ |                     *
c *               |     #_______________________\#                     *
c *               |    2         (2)             3                     *
c *               |                                                    *
c *               v fi2             -------> fi1                       *
c *                                                                    *
c *                                                                    *
c **********************************************************************
 
      pi=3.14159
      pi2=pi/2.
      nrai(ip)=4

      dys(1,1,ip)=pi2+fi1min
      dys(2,1,ip)=pi2+fi2min
      dys(1,2,ip)=pi2+fi1min
      dys(2,2,ip)=pi2+fi2max
      dys(1,3,ip)=pi2+fi1max
      dys(2,3,ip)=pi2+fi2max
      dys(1,4,ip)=pi2+fi1max
      dys(2,4,ip)=pi2+fi2min

      do i=1,2
        do j=1,4
          iys(i,j,ip)=0
        enddo
      enddo
c--
      do ipoint=1,nrai(ip)
         call newpoint (ip,xsrc,ipoint,usrc,
     &             vm,xm,ym,zm,nxm,nym,nzm,dxm,dym,dzm,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)
      enddo
c--
      ncel(ip)=2
      nseg(ip)=5
c--
      icel(1,1,ip)=1
      icel(1,2,ip)=2
      icel(1,3,ip)=3
      icel(1,4,ip)=1
      icel(1,5,ip)=2
      icel(1,6,ip)=3

      icel(2,1,ip)=1
      icel(2,2,ip)=4
      icel(2,3,ip)=3
      icel(2,4,ip)=3
      icel(2,5,ip)=4
      icel(2,6,ip)=5
c--
      iseg(1,1,ip)=1
      iseg(1,2,ip)=2
      iseg(1,3,ip)=0
      iseg(1,4,ip)=0
      iseg(1,5,ip)=0
      iseg(1,6,ip)=0

      iseg(2,1,ip)=2 
      iseg(2,2,ip)=3 
      iseg(2,3,ip)=0 
      iseg(2,4,ip)=0 
      iseg(2,5,ip)=0 
      iseg(2,6,ip)=0 

      iseg(3,1,ip)=1 
      iseg(3,2,ip)=3 
      iseg(3,3,ip)=0 
      iseg(3,4,ip)=0 
      iseg(3,5,ip)=0 
      iseg(3,6,ip)=0 

      iseg(4,1,ip)=1 
      iseg(4,2,ip)=4 
      iseg(4,3,ip)=0 
      iseg(4,4,ip)=0 
      iseg(4,5,ip)=0 
      iseg(4,6,ip)=0 

      iseg(5,1,ip)=4 
      iseg(5,2,ip)=3 
      iseg(5,3,ip)=0 
      iseg(5,4,ip)=0 
      iseg(5,5,ip)=0 
      iseg(5,6,ip)=0 

      return
      end       
