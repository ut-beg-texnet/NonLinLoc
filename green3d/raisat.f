       subroutine raisat (ip,indix,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)

c-------------------------------------------------------------------------
c-     Programed by P.S.Lucio and G.C.Lambare in Dec/94                  -
c-     Last modification by P.S.Lucio in Jun/12/95                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c-------------------------------------------------------------------------

c    ***************************************************************
c       an additional "ray" is inserted on the previous wavefront
c    ***************************************************************

c   ****************************************************************
c   The segment is associated to the two extremities of rays (X1,X2)
c   ****************************************************************

c*******************************************************************
c INPUT  : ip - the wavefront
c          indix - the wavefront indicator
c
c OUTPUT : nrai(ip) - the midpoint of the segment
c          isatr(ip) - index of the wavefront saturation
c******************************************************************

      implicit none
      integer ip,indix
 
      integer npn,nsn,ncn,nfront,nd,nk,npar 
      real ys(npar,npn,nfront) 
      integer iys(nk,npn,nfront) 
      double precision dys(nd,npn,nfront) 
      integer  iseg(nsn,6,nfront),icel(ncn,6,nfront) 
      integer nrai(nfront), ncel(nfront), nseg(nfront) 
      integer irai(3), iraf(3), irpop(3), irais(2,3) 
      integer isatr(nfront), isats(nfront), isatc(nfront) 
c--

      if (nrai(ip).ge.npn) then
        isatr(ip)=1
        write(*,*) 'THE WAVEFRONT',indix-1,' IS SATURETED'
        write(*,*) 'NRAY=',nrai(ip)
        write(*,*) 'WE WILL NOT SUBDIVISE ANYMORE'
      else
        isatr(ip)=0
      endif

      return
      end
