       subroutine segsat (ip,isegnew,indix,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)

c-------------------------------------------------------------------------
c-     Programed by P.S.Lucio and G.C.Lambare in Dec/94                  -
c-     Last modification by P.S.Lucio in Jun/12/95                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c-------------------------------------------------------------------------

c*************************************************************************
c*              addition of a segment on the second wavefront            *
c*************************************************************************

c****************************************************** 
c INPUT  : ip - the wavefront 
c
c OUTPUT : isegnew - the new segment on the front
c          isats(ip) - index of wavefront saturation
c******************************************************
       implicit none
       integer ip,isegnew,indix
 
      integer npn,nsn,ncn,nfront,nd,nk,npar 
      real ys(npar,npn,nfront) 
      integer iys(nk,npn,nfront) 
      double precision dys(nd,npn,nfront) 
      integer  iseg(nsn,6,nfront),icel(ncn,6,nfront) 
      integer nrai(nfront), ncel(nfront), nseg(nfront) 
      integer irai(3), iraf(3), irpop(3), irais(2,3) 
      integer isatr(nfront), isats(nfront), isatc(nfront)

c --  the segment satisfy the criterion then this segment is included  --
c --                        to the wavefront                           --

       if (isegnew.gt.nsn) then
         isats(ip)=1
         write(*,*) 'THE WAVEFRONT',indix-1,' IS SATURETED'
         write(*,*) 'NSEG=',nseg(ip)
         write(*,*) 'WE WILL NOT SUBDIVISE ANYMORE'
         return
       else
         isats(ip)=0
       endif
       
       return
       end

