      subroutine update (yf,ier,vm,xm,ym,zm,nxm,nym,nzm,dxm,dym,dzm)

c-------------------------------------------------------------------------
c-     programed by p.s.lucio and g.c.lambare in aug/94                  -
c-     last modification by p.s.lucio in dec/29/94                       -
c-------------------------------------------------------------------------

C      **********************************************************
C       CALCULUS OF A POINT UPDATED PARAMETERS ON THE WAVEFRONT
C      **********************************************************
c--
      implicit none
      integer nxm,nym,nzm
      real dxm,dym,dzm
      real vm(nzm,nxm,nym)
      real xm(2),ym(2),zm(2)
      real yf(36)
      integer ier
      real v,dvx,dvy,dvz
      real d2vxx,d2vxy,d2vxz,d2vyy,d2vyz,d2vzz

c*********************************************************************
c*                        interpolation values                       *
c*********************************************************************

      call lenteur2(yf(1),yf(2),yf(3),v,dvx,dvy,dvz,
     &              d2vxx,d2vxy,d2vxz,d2vyy,d2vyz,d2vzz,ier,
     &              vm,xm,ym,zm,nxm,nym,nzm,dxm,dym,dzm)
      if (ier.ne.0) return

      yf(27)=dvx*.5
      yf(28)=dvy*.5
      yf(29)=dvz*.5
      yf(30)=v
      yf(31)=d2vxx*.5
      yf(32)=d2vxy*.5
      yf(33)=d2vxz*.5
      yf(34)=d2vyy*.5
      yf(35)=d2vyz*.5
      yf(36)=d2vzz*.5
      return
      end       
