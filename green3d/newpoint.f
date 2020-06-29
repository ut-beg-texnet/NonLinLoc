      subroutine newpoint (ip,xsrc,ipoint,usrc,
     &       vm,xm,ym,zm,nxm,nym,nzm,dxm,dym,dzm,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)
c-------------------------------------------------------------------------
c-     Programed by P.S.Lucio and G.C.Lambare in Aug/94                  -
c-     Last modification by P.S.lucio in Oct/20/95                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c-------------------------------------------------------------------------

C      **********************************************************
C       CALCULUS OF A POINT PARAMETERS ON THE INITIAL WAVEFRONT
C      **********************************************************
      implicit none
 
      integer npn,nsn,ncn,nfront,nd,nk,npar   
      real ys(npar,npn,nfront)
      integer iys(nk,npn,nfront)
      double precision dys(nd,npn,nfront)
      integer  iseg(nsn,6,nfront),icel(ncn,6,nfront)
      integer nrai(nfront), ncel(nfront), nseg(nfront)
      integer irai(3), iraf(3), irpop(3), irais(2,3)
      integer isatr(nfront), isats(nfront), isatc(nfront)

      integer ip,ipoint

      real xm(2),ym(2),zm(2)
      integer nxm,nym,nzm
      real dxm,dym,dzm
      real vm(nzm,nxm,nym)
      real xsrc(3),usrc
      real v,dvx,dvy,dvz
      real d2vxx,d2vxy,d2vxz,d2vyy,d2vyz,d2vzz
      integer ier,j

C**********************************************************************
C                    CALCULUS OF THE SQUARED SLOWNESS                 *
C**********************************************************************

      ier=1
      call lenteur2(xsrc(1),xsrc(2),xsrc(3),v,dvx,dvy,dvz,
     &              d2vxx,d2vxy,d2vxz,d2vyy,d2vyz,d2vzz,ier,
     &              vm,xm,ym,zm,nxm,nym,nzm,dxm,dym,dzm)
      if (ier.ne.0) then
        write(*,*) 'source en dehors du milieu'
        write(*,*)xsrc(1),xsrc(2),xsrc(3)
        stop
      endif
      usrc=sqrt(v)

      ys(1,ipoint,ip)=xsrc(1)
      ys(2,ipoint,ip)=xsrc(2)
      ys(3,ipoint,ip)=xsrc(3)
      do j=7,12
         ys(j,ipoint,ip)=0.
      enddo
      ys(4,ipoint,ip)=cos(dys(1,ipoint,ip))*usrc
      ys(5,ipoint,ip)=sin(dys(1,ipoint,ip))*cos(dys(2,ipoint,ip))*usrc
      ys(6,ipoint,ip)=sin(dys(1,ipoint,ip))*sin(dys(2,ipoint,ip))*usrc
      ys(13,ipoint,ip)=-sin(dys(1,ipoint,ip))*usrc
      ys(14,ipoint,ip)=cos(dys(1,ipoint,ip))*cos(dys(2,ipoint,ip))*usrc
      ys(15,ipoint,ip)=cos(dys(1,ipoint,ip))*sin(dys(2,ipoint,ip))*usrc
      ys(16,ipoint,ip)=0.
      ys(17,ipoint,ip)=-sin(dys(1,ipoint,ip))*sin(dys(2,ipoint,ip))*usrc
      ys(18,ipoint,ip)=sin(dys(1,ipoint,ip))*cos(dys(2,ipoint,ip))*usrc
      ys(19,ipoint,ip)=0.
      ys(20,ipoint,ip)=0.
      do j=1,6
         ys(20+j,ipoint,ip)=ys(j,ipoint,ip)
      enddo
c--
      ys(27,ipoint,ip)=dvx*.5
      ys(28,ipoint,ip)=dvy*.5
      ys(29,ipoint,ip)=dvz*.5
      ys(30,ipoint,ip)=v
      ys(31,ipoint,ip)=d2vxx*.5
      ys(32,ipoint,ip)=d2vxy*.5
      ys(33,ipoint,ip)=d2vxz*.5
      ys(34,ipoint,ip)=d2vyy*.5
      ys(35,ipoint,ip)=d2vyz*.5
      ys(36,ipoint,ip)=d2vzz*.5
c--
c     --- Initialisation of Jacobian for 3D preserved amplitude migration
      ys(37,ipoint,ip)=1.
      ys(38,ipoint,ip)=0.
      ys(39,ipoint,ip)=0.
      ys(40,ipoint,ip)=0.
      ys(41,ipoint,ip)=1.
      ys(42,ipoint,ip)=0.
      ys(43,ipoint,ip)=.25*dvx/v*ys(4,ipoint,ip)
      ys(44,ipoint,ip)=.25*dvx/v*ys(5,ipoint,ip)
      ys(45,ipoint,ip)=.25*dvx/v*ys(6,ipoint,ip)
      ys(46,ipoint,ip)=.25*dvy/v*ys(4,ipoint,ip)
      ys(47,ipoint,ip)=.25*dvy/v*ys(5,ipoint,ip)
      ys(48,ipoint,ip)=.25*dvy/v*ys(6,ipoint,ip)
c     ------------------------------------------

      usrc=usrc/(4.*3.14159)**2
      return
      end       
