       subroutine addseg (ip,irai1,irai2,isegnew,indix,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)

c-------------------------------------------------------------------------
c-     Programed by P.S.Lucio and G.C.Lambare in Dec/94                  -
c-     last modification by P.S.Lucio in jun/20/95                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c-------------------------------------------------------------------------

c*************************************************************************
c*              addition of a segment on the second wavefront            *
c*************************************************************************

c*************************************************************************
c INPUT  : ip - the wavefront                                            *
c          irai - the points of the triangular cell                      *
c                                                                        *
c OUTPUT : isegnew - the new segment on the front                        *
c          isats(ip) - index of wavefront saturation                     *
c*************************************************************************

       implicit none
       integer ip,irai1,irai2,isegnew,indix

      integer npn,nsn,ncn,nfront,nd,nk,npar
      real ys(npar,npn,nfront)
      integer iys(nk,npn,nfront)
      double precision dys(nd,npn,nfront)
      integer  iseg(nsn,6,nfront),icel(ncn,6,nfront)
      integer nrai(nfront), ncel(nfront), nseg(nfront)
      integer irai(3), iraf(3), irpop(3), irais(2,3)
      integer isatr(nfront), isats(nfront), isatc(nfront)
c-

c --  the segment satisfy the criterion then this segment is included  --
c --                        to the wavefront                           --

c*************************************************************************
       call segsat (ip,isegnew,indix,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)
       if (isats(ip).ne.0) return
c*************************************************************************
c--
       iseg(isegnew,1,ip)=irai1
       iseg(isegnew,2,ip)=irai2
       iseg(isegnew,3,ip)=0
       iseg(isegnew,4,ip)=0
       iseg(isegnew,5,ip)=0
       iseg(isegnew,6,ip)=0
c--
       return
       end

