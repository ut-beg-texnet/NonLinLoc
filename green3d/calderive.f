       subroutine calderive(ip,ipoint,y,dy,dfi,dto,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)

c-------------------------------------------------------------------------
c-     Programed by P.S.Lucio and G.C.Lambare in Jun/95                  -
c-     Last modification by P.S.Lucio in Jun/20/95                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c-------------------------------------------------------------------------

C*************************************************************************
C This is a SP for the calculus of the derivatives that will be employed *
C      in the SR for the Cubic Hermite Interpolation between rays.       *
C*************************************************************************

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
       real dto
       real y(20),dey(20),dy(20)
       double precision dfi(2)
       integer i
       real pscal
       
c      **** Derivatives of x and p *************

       call derivebis (ip,ipoint,dey,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)

       do 1 i=1,3
          dy(i)=y(i+6)*dfi(1)+y(i+9)*dfi(2)+y(i+3)*dto
          dy(i+3)=y(i+12)*dfi(1)+y(i+15)*dfi(2)+dey(i+3)*dto
1      continue

c      **** Derivatives of time ***********

       dy(19)=pscal(y(4),y( 7))*dfi(1)
     &       +pscal(y(4),y(10))*dfi(2)+dey(19)*dto

       return
       end

c *************************************************************
c *************************************************************

       function pscal(a,b)
       implicit none

       real a(3),b(3),pscal

       pscal=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)

       return
       end
