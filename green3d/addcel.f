       subroutine addcel (ip,irai1,irai2,irai3,iseg1,iseg2,iseg3,
     &                    icelnew,indix,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)

c-------------------------------------------------------------------------
c-     Programed by P.S.Lucio and G.C.Lambare in Dec/94                  -
c-     Last modification by P.S.Lucio in Jun/20/95                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c-------------------------------------------------------------------------

c*************************************************************************
c*            addition of a triangle on the courent wavefront            *
c*************************************************************************

c*************************************************************************
c INPUT  : ip - the wavefront                                            *
c          irai - the points of the triangular cell                      *
c          iseg - the segments of this cell                              *
c                                                                        *
c OUTPUT : icelnew - the new triangular cell                             *
c          isatc(ip) - index of wavefront saturation                     *
c*************************************************************************

       implicit none
       integer ip,irai1,irai2,irai3,iseg1,iseg2,iseg3,icelnew,indix

      integer npn,nsn,ncn,nfront,nd,nk,npar
      real ys(npar,npn,nfront)
      integer iys(nk,npn,nfront)
      double precision dys(nd,npn,nfront)
      integer  iseg(nsn,6,nfront),icel(ncn,6,nfront)
      integer nrai(nfront), ncel(nfront), nseg(nfront)
      integer irai(3), iraf(3), irpop(3), irais(2,3)
      integer isatr(nfront), isats(nfront), isatc(nfront)
c--


c*************************************************************************
       call celsat (ip,icelnew,indix,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)

       if (isatc(ip).ne.0) return
c*************************************************************************

c --  the segments satisfy the criterion then the triangle is included   --
c --                          to the wavefront                           --

       icel(icelnew,1,ip)=irai1
       icel(icelnew,2,ip)=irai2
       icel(icelnew,3,ip)=irai3
       icel(icelnew,4,ip)=iseg1
       icel(icelnew,5,ip)=iseg2 
       icel(icelnew,6,ip)=iseg3
       return
       end

