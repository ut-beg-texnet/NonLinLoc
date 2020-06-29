c****************************************************************************
c****************************************************************************

      subroutine main(nxm,nym,nzm,nxr,nyr,nzr,nir,npr,
     & vm,xmap,imap,ntab,tab)
c--
c-------------------------------------------------------------------------
c-     Programed by P.S.Lucio and G.C.Lambare in Apr/95                  -
c-     Last modification by P.S.Lucio in Jul/10/95                       -
c-------------------------------------------------------------------------
c--
      implicit none
      integer ntab
      real tab(ntab)
      integer nxm,nym,nzm
      integer nxr,nyr,nzr,nir,npr
      real xm(2),ym(2),zm(2)
      real vm(nzm,nxm,nym)
      real dxm,dym,dzm

      real xsrc(3)

      real xmap(npr,nzr,nxr,nyr,nir)
      integer imap(nzr,nxr,nyr)
      real targ_orient(3,7)
      real fi1min,fi2min,fi1max,fi2max,dxmin2,dpmin2,dtemps

c     --- Enter the velocity filed description ---------------------
      call modelvit(vm,xm,ym,zm,dxm,dym,dzm,nxm,nym,nzm)

c     --- Enter the reservoir description --------------------------
      call modelrese(targ_orient)

c     --- Enter the ray tracing parameters -------------------------
      call modelray(xsrc,
     &      fi1min,fi2min,fi1max,fi2max,dxmin2,dpmin2,dtemps)
      
      call green3a(nxm,nym,nzm,dxm,dym,dzm,xm,ym,zm,vm,
     &     nxr,nyr,nzr,nir,npr,targ_orient,xmap,imap,
     &     xsrc,fi1min,fi2min,fi1max,fi2max,dxmin2,dpmin2,dtemps,
     &     ntab,tab)

      write(*,*)'fi1min,fi2min,fi1max,fi2max (FIN)'
      write(*,*) fi1min,fi2min,fi1max,fi2max 

      call sortie(nxr,nyr,nzr,nir,npr,xmap,imap)

1     continue

      return
      end
