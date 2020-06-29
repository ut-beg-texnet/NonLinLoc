      subroutine cellule 
     & (ip1,ip2,indix,iechec,
     &  vm,xm,ym,zm,nxm,nym,nzm,dxm,dym,dzm,dtemps,
     &      xmin,xmax,ymin,ymax,zmin,zmax,zrmin,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)

c-------------------------------------------------------------------------
c-     Programed by P.S.Lucio and G.C.Lambare in Nov/94                  -
c-     last modification by p.s.lucio in jun/19/95                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c-------------------------------------------------------------------------

c     ***************************************************************
c     *      PROGRAM FOR THE PROPAGATION OF A ELEMENTARY CELL       *
c     ***************************************************************

      implicit none
       

      integer npn,nsn,ncn,nfront,nd,nk,npar
      real ys(npar,npn,nfront)
      integer iys(nk,npn,nfront)
      double precision dys(nd,npn,nfront)
      integer  iseg(nsn,6,nfront),icel(ncn,6,nfront)
      integer nrai(nfront), ncel(nfront), nseg(nfront)
      integer irai(3), iraf(3), irpop(3), irais(2,3)
      integer isatr(nfront), isats(nfront), isatc(nfront)

c--
      integer ip1,ip2,indix,iechec
      integer nxm,nym,nzm
      real vm(nzm,nxm,nym)
      real xm(2),ym(2),zm(2)
      real dxm,dym,dzm
      real xmin,xmax,ymin,ymax,zmin,zmax,zrmin
      real dtemps

C      **************************************************************
C      *    iYs[irai](1) : the KMAH index                           *
C      *    iYs[irai](2) : the indicator of the ray propagation     *
C      *                   -1 echec for the ray [irai]              *
C      *                    number of the ray [irai] propageted     *
C      *                    0 ray to be propageted                  *
C      *    iYs[irai](i) : output                                   *
C      **************************************************************

      iechec=-1
c--
      if (iys(2,irai(1),ip1).eq.-1) return
      if (iys(2,irai(2),ip1).eq.-1) return
      if (iys(2,irai(3),ip1).eq.-1) return
c--
      iechec=0
      if (iys(2,irai(1),ip1).eq.0) then
        call runge4
     &       (ip1,irai(1),indix,ys(1,irai(1),ip1),iys(1,irai(1),ip1),
     &                  ys(1,nrai(ip2)+1,ip2),iys(1,nrai(ip2)+1,ip2),
     &                  iechec,vm,xm,ym,zm,nxm,nym,nzm,dxm,dym,dzm,
     &                  dtemps,xmin,xmax,ymin,ymax,zmin,zmax,zrmin,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)
        call propage (ip1,ip2,iechec,irai(1),
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)
      endif
      if (iechec.ne.0) return

      if (iys(2,irai(2),ip1).eq.0) then
        call runge4
     &       (ip1,irai(2),indix,ys(1,irai(2),ip1),iys(1,irai(2),ip1),
     &                  ys(1,nrai(ip2)+1,ip2),iys(1,nrai(ip2)+1,ip2),
     &                  iechec,vm,xm,ym,zm,nxm,nym,nzm,dxm,dym,dzm,
     &                  dtemps,xmin,xmax,ymin,ymax,zmin,zmax,zrmin,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)
        call propage (ip1,ip2,iechec,irai(2),
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)
      endif
      if (iechec.ne.0) return

      if (iys(2,irai(3),ip1).eq.0) then
        call runge4
     &       (ip1,irai(3),indix,ys(1,irai(3),ip1),iys(1,irai(3),ip1),
     &                  ys(1,nrai(ip2)+1,ip2),iys(1,nrai(ip2)+1,ip2),
     &                  iechec,vm,xm,ym,zm,nxm,nym,nzm,dxm,dym,dzm,
     &                  dtemps,xmin,xmax,ymin,ymax,zmin,zmax,zrmin,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)
        call propage (ip1,ip2,iechec,irai(3),
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)
      endif
      if (iechec.ne.0) return

      return
      end
