       subroutine testseg (ip,nrai1,nrai2,indix,itest,dxmin2,dpmin2,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)

c-------------------------------------------------------------------------
c-     Programed by P.S.Lucio and G.C.Lambare in Oct/94                  -
c-     Last modification by P.S.Lucio in Oct/10/95                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c-------------------------------------------------------------------------

       implicit none
       integer ip,nrai1,nrai2,indix,itest
       real dxmin2,dpmin2
c--
      integer npn,nsn,ncn,nfront,nd,nk,npar
      real ys(npar,npn,nfront)
      integer iys(nk,npn,nfront)
      double precision dys(nd,npn,nfront)
      integer  iseg(nsn,6,nfront),icel(ncn,6,nfront)
      integer nrai(nfront), ncel(nfront), nseg(nfront)
      integer irai(3), iraf(3), irpop(3), irais(2,3)
      integer isatr(nfront), isats(nfront), isatc(nfront)
c--

       real yi1(6),yi2(6),dey(20)
       double precision fi(2)
       integer i,ecartx2,ecartp2,dto

c *******************************************************************
c *    INPUT :     ip    : number of the front                      *
c *                nrai1 : number of the first ray                  *
c *                nrai2 : number of the second ray                 *
c *                nrai3 : number of the third ray                  *
c *    OUTPUT :    itest :  -1 we subdivise,  0 the segment is "ok" *
c *******************************************************************

c      --- Calculus of the misfit in tau, phi1 and phi2 between the ---
c                      --- adjacents rays ---

       do i=1,2
          fi(i)=dys(i,nrai2,ip)-dys(i,nrai1,ip)
       enddo

       dto=ys(20,nrai2,ip)-ys(20,nrai1,ip)

       call derivebis (ip,nrai1,dey,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)

       do i=1,3
          yi1(i)=ys(i,nrai1,ip)+ys(i+6 ,nrai1,ip)*fi(1)
     &                         +ys(i+9 ,nrai1,ip)*fi(2)
     &                         +dey(i)*dto
          yi1(i+3)=ys(i+3,nrai1,ip)+ys(i+12 ,nrai1,ip)*fi(1)
     &                             +ys(i+15 ,nrai1,ip)*fi(2)
     &                             +dey(i+3)*dto
       enddo
c--
       do i=1,3
          yi2(i)=ys(i,nrai2,ip)
          yi2(i+3)=ys(i+3,nrai2,ip)
       enddo

c--
       ecartx2=0.
       ecartp2=0.
c--
       do i=1,3
          ecartx2=ecartx2+(yi1(i  )-yi2(i  ))**2
          ecartp2=ecartp2+(yi1(i+3)-yi2(i+3))**2
       enddo
c--
       if ((ecartx2.gt.dxmin2).or.(ecartp2.gt.dpmin2)) then
          itest=-1
       endif

       return
       end
