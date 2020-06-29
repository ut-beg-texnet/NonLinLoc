      subroutine runge4 (ip,ipoint,ind,yi,iyi,yf,iyf,ier,
     &                   vm,xm,ym,zm,nxm,nym,nzm,dxm,dym,dzm,
     &                   dtemps,xmin,xmax,ymin,ymax,zmin,zmax,zrmin,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)

c-------------------------------------------------------------------------
c-     Programed by P.S.Lucio and G.C.Lambare in May/94                  -
c-     Last modification by P.S.Lucio in Oct/20/95                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c-------------------------------------------------------------------------

C*********************************************************************
C*       PROGRAM FOR THE CALCULUS OF THE RAYS ON A 3D FIELD BY       *
C*                       A 4TH ORDER RUNGE-KUTTA                     *
C*  (INTEGRATION BY CONSTANT STEP FOR THE PROPAGATOR INTEGRATION)    *
C*********************************************************************
c--
      implicit none

      integer npn,nsn,ncn,nfront,nd,nk,npar
      real ys(npar,npn,nfront)
      integer iys(nk,npn,nfront)
      double precision dys(nd,npn,nfront)
      integer  iseg(nsn,6,nfront),icel(ncn,6,nfront)
      integer nrai(nfront), ncel(nfront), nseg(nfront)
      integer irai(3), iraf(3), irpop(3), irais(2,3)
      integer isatr(nfront), isats(nfront), isatc(nfront)


      integer iyi,iyf
      real hmahi,hmahf
      real yi(48), yf(48)
      real y0(48),y1(48),y2(48),y3(48),y4(48)
      real deyi(48), dey0(48),dey3(48),dey4(48)
c--
      integer nxm,nym,nzm
      real dxm,dym,dzm
      real xm(2),ym(2),zm(2)
      real vm(nzm,nxm,nym)
c--
      integer i,ier
      real  dtemps,xmin,xmax,ymin,ymax,zmin,zmax,zrmin
      real temps,dttau,dta2,dta3,dta6
      integer ind,ip,ipoint


      call deter (yi(4),yi(7),yi(10),hmahi)
      call derivebis (ip,ipoint,deyi,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)

      temps=(ind-1)*dtemps-yi(19)
      dttau=temps/deyi(19)

      dta2=dttau/2.
      dta3=dttau/3.
      dta6=dttau/6.

c***********************************************************************
c*                     FIRST STEP FOR RUNGE KUTTA                      *
c***********************************************************************
 
      do 1 i=1,20
         y0(i)=yi(i)+deyi(i)*dta2
1     continue
      do 11 i=37,48
         y0(i)=yi(i)+deyi(i)*dta2
11     continue

      ier=0
      call derive (y0,dey0,ier,vm,xm,ym,zm,nxm,nym,nzm,dxm,dym,dzm)
      if (ier.ne.0) return

c***********************************************************************
c*                     SECOND STEP FOR RUNGE KUTTA                     *
c***********************************************************************

      do 2 i=1,20
         y1(i)=yi(i)+dey0(i)*dta3+deyi(i)*dta6
         y3(i)=yi(i)+dey0(i)*dta2
2     continue 
      do 22 i=37,48
         y1(i)=yi(i)+dey0(i)*dta3+deyi(i)*dta6
         y3(i)=yi(i)+dey0(i)*dta2
22    continue

      call derive (y3,dey3,ier,vm,xm,ym,zm,nxm,nym,nzm,dxm,dym,dzm)
      if (ier.ne.0) return

c***********************************************************************
c*                     THIRD STEP FOR RUNGE KUTTA                      *
c***********************************************************************
 
      do 3 i=1,20
         y4(i)=yi(i)+dey3(i)*dttau
         y2(i)=y1(i)+dey3(i)*dta3
3     continue
      do 33 i=37,48
         y4(i)=yi(i)+dey3(i)*dttau
         y2(i)=y1(i)+dey3(i)*dta3
33    continue

      call derive (y4,dey4,ier,vm,xm,ym,zm,nxm,nym,nzm,dxm,dym,dzm)
      if (ier.ne.0) return
c***********************************************************************
c*                     FOURTH STEP FOR RUNGE KUTTA                     *
c***********************************************************************
 
      do 4 i=1,20
         yf(i)=y2(i)+dey4(i)*dta6
4     continue
      do 44	 i=37,48
         yf(i)=y2(i)+dey4(i)*dta6
44    continue
c*********************************************************************
c*                        updating the values                        *
c*********************************************************************
      do i=1,6
        yf(20+i)=yf(i)
      enddo

      call deter (yf(4),yf(7),yf(10),hmahf)
      iyf = iyi
      if ((hmahi.ne.0.).and.((hmahf*hmahi).le.0)) then
        iyf = iyf + 1
      endif

      call update (yf,ier,vm,xm,ym,zm,nxm,nym,nzm,dxm,dym,dzm)
      if (ier.ne.0) return

      if ( (yf(1)-xmin)*(yf(1)-xmax).gt.0.0 .or.
     &     (yf(2)-ymin)*(yf(2)-ymax).gt.0.0 .or.
     &     (yf(3)-zmin)*(yf(3)-zmax).gt.0.0 ) ier=99
      if   ((yf(3).lt.zrmin).and.(yf(6).le.0.0)) ier=99

      return
      end                                    
