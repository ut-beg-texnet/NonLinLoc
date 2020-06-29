       subroutine addrai (ip,irai1,irai2,indix,xsrc,
     &                    vm,xm,ym,zm,nxm,nym,nzm,dxm,dym,dzm,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)
 

c-------------------------------------------------------------------------
c-     Programed by P.S.Lucio and G.C.Lambare in Dec/94                  -
c-     Last modification by P.S.Lucio in Oct/20/95                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c-------------------------------------------------------------------------

c    ***************************************************************
c       an additional "ray" is inserted on the previous wavefront
c    ***************************************************************

c*******************************************************************
c INPUT  : ip - the wavefront                                      *
c          irai - the points of the triangular cell                *
c          indix - the wavefront indicator                         *
c                                                                  *
c OUTPUT : nrai(ip) - the midpoint of the segment                  *
c          isatr(ip) - index of the wavefront saturation           *
c*******************************************************************

       implicit none
       integer ip,irai1,irai2,indix

      integer npn,nsn,ncn,nfront,nd,nk,npar
      real ys(npar,npn,nfront)
      integer iys(nk,npn,nfront)
      double precision dys(nd,npn,nfront)
      integer  iseg(nsn,6,nfront),icel(ncn,6,nfront)
      integer nrai(nfront), ncel(nfront), nseg(nfront)
      integer irai(3), iraf(3), irpop(3), irais(2,3)
      integer isatr(nfront), isats(nfront), isatc(nfront)
c-

       integer nxm,nym,nzm
       real vm(nzm,nxm,nym)
       real xm(2),ym(2),zm(2)
       real dxm,dym,dzm
c--
       real xsrc(3),usrc
c--
       double precision dfi(2), dphi(2)
       real dys1(48),dys2(48)
       integer kmahi,i
       integer isigni,ier,isigni1,isigni2,ind

       real hmahi1,hmahi2,dto,hmahi

c--       
c **********************************************************************
       call raisat (ip,indix,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)


       if (isatr(ip).ne.0) return
c **********************************************************************

       nrai(ip)=nrai(ip)+1

c      ### calculus of the derivatives of fi1 fi2 tau under the segment ### 

       dto=ys(20,irai2,ip)-ys(20,irai1,ip)   

c      ---     interpolation of tau     ---

       ys(20,nrai(ip),ip)=ys(20,irai1,ip)+dto*.5

c     --- medians angles ------------------------------------
       do i=1,2
          dfi(i)=dys(i,irai2,ip)-dys(i,irai1,ip)
          dys(i,nrai(ip),ip)=dys(i,irai1,ip)+dfi(i)*.5
          dphi(i) = dys(i,nrai(ip),ip)
       enddo

       call deter (ys(4,irai1,ip),ys(7,irai1,ip),ys(10,irai1,ip),hmahi1)
       isigni1=1 
       if (hmahi1.lt.0) isigni1=-1

       call deter (ys(4,irai2,ip),ys(7,irai2,ip),ys(10,irai2,ip),hmahi2)
       isigni2=1
       if (hmahi2.lt.0) isigni2=-1

c     ######### ON THE PREVIOUS WAVEFRONT ####################
c     ######### SOURCE POINTS ################################

       if (indix.eq.2) then
          call newpoint(ip,xsrc,nrai(ip),usrc,
     &                 vm,xm,ym,zm,nxm,nym,nzm,dxm,dym,dzm,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)
       else

c     ######### FRONT SUCCESSIF  ############################
 
       call calderive (ip,irai1,ys(1,irai1,ip),dys1,dfi,dto,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)

       call calderive (ip,irai2,ys(1,irai2,ip),dys2,dfi,dto,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)

       do 2 ind=1,19
          if ((ind.ge.7).and.(ind.le.18)) goto 2

c      --- propagation ----------------------------------------

          call expanh3 (ys(ind,irai1,ip),ys(ind,irai2,ip),dys1(ind),
     &                  dys2(ind),ys(ind,nrai(ip),ip))
2      continue

c      --- interpolation of the paraxials ----------------------
c      --- linear interpolation (weigth = .5) ------------------

       do 3 ind=7,18
          ys(ind,nrai(ip),ip)=
     &          (ys(ind,irai2,ip)+ys(ind,irai1,ip))*.5
          ys(ind+30,nrai(ip),ip)=
     &          (ys(ind+30,irai2,ip)+ys(ind+30,irai1,ip))*.5
3      continue

c      --- interpolation ---------------------------------------

       do 4 ind=1,6
          ys(20+ind,nrai(ip),ip)=(ys(ind,irai2,ip)+ys(ind,irai1,ip))*.5
4      continue

       iys(2,nrai(ip),ip)=0
       endif

c      --- the kmah index -------------------------------------- 

       call deter (ys( 4,nrai(ip),ip),ys(7,nrai(ip),ip),
     &             ys(10,nrai(ip),ip),hmahi)

       isigni=1  
       if (hmahi.lt.0) isigni=-1
       if (isigni.eq.isigni1) then
         kmahi = iys(1,irai1,ip)
       elseif (isigni.eq.isigni2) then
         kmahi = iys(1,irai2,ip)
       else
         kmahi = max(iys(1,irai1,ip),iys(1,irai2,ip)) + 1
       endif

       iys(1,nrai(ip),ip)=kmahi
c--       
c      --- updating -------------------------------------------

       ier=0
       call update 
     & (ys(1,nrai(ip),ip),ier,vm,xm,ym,zm,nxm,nym,nzm,dxm,dym,dzm)
       if (ier.ne.0) return
c--
       return
       end
