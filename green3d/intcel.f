      subroutine intcel (ip1,ip2,nxr,nyr,nzr,nir,npr,
     &           usrc,targ_orient,xmap,imap,
     &           fi1min,fi1max,fi2min,fi2max,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)
c-------------------------------------------------------------------------
c-     Programed by P.S.Lucio and G.C.Lambare in Dec/94                  -
c-     Last modification by P.S.Lucio in Jun/21/95                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c-------------------------------------------------------------------------
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

      integer ip1,ip2,nxr,nyr,nzr,nir,npr
      real usrc,fi1min,fi1max,fi2min,fi2max

c     ***************************************************************
c     *    PROGRAM TO PREPARE THE INTERPOLATION INSIDE A CELL       *
c     ***************************************************************
      real targ_orient(3,7)
      real xmap(npr,nzr,nxr,nyr,nir)
      integer imap(nzr,nxr,nyr)
c--

      real yt1(36,4),yt2(36,4),yt3(36,4)
      integer iyt1(1,4),iyt2(1,4),iyt3(1,4)
c--
      call critere (ip1,ip2,yt1,yt2,yt3,iyt1,iyt2,iyt3,targ_orient,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)
c--
c     ----  ---> Interpolation inside the first tetrahedrum -----

      call intetra (yt1,iyt1,nxr,nyr,nzr,nir,npr,usrc,xmap,imap,
     &           fi1min,fi1max,fi2min,fi2max)

c     ----  ---> Interpolation inside the second tetrahedrum ----

      call intetra (yt2,iyt2,nxr,nyr,nzr,nir,npr,usrc,xmap,imap,
     &           fi1min,fi1max,fi2min,fi2max)

c     ----  ---> Interpolation inside the third tetrahedrum ---

      call intetra (yt3,iyt3,nxr,nyr,nzr,nir,npr,usrc,xmap,imap,
     &           fi1min,fi1max,fi2min,fi2max)
c--

      return
      end
