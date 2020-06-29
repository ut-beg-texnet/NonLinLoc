       subroutine stock (ip1,ip2,ncrai,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)
c-------------------------------------------------------------------------
c-     programed by P.S.Lucio and G.C.Lambare in Oct/94                  -
c-     Last modification by P.S.Lucio in Jan/10/95                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c-------------------------------------------------------------------------

C ******************************************************************
C * INPUT  ncrai : number of the cell                              *
C *        ip    : index of wavefront                              *
C *        isatc(ip) : index of the wavefront saturation (cells)   *
C *                                                                *
C * OUTPUT The endpoints of the cell                               *
C ******************************************************************

C ******************************************************************
C *                 INITIALISATION OF A NEW FRONT                  *
c ******************************************************************
C     **********************************************************
C     *            STOCKING OF THE STUMPS OF THE RAYS          *
C     **********************************************************
      implicit none
      integer ip1,ip2,ncrai,i
      integer npn,nsn,ncn,nfront,nd,nk,npar
      real ys(npar,npn,nfront)
      integer iys(nk,npn,nfront)
      double precision dys(nd,npn,nfront)
      integer  iseg(nsn,6,nfront),icel(ncn,6,nfront)
      integer nrai(nfront), ncel(nfront), nseg(nfront)
      integer irai(3), iraf(3), irpop(3), irais(2,3)
      integer isatr(nfront), isats(nfront), isatc(nfront)
c--
      if (iys(2,irai(1),ip1).gt.0.and.iys(2,irai(2),ip1).gt.0
     &    .and.iys(2,irai(3),ip1).gt.0) then

c     ---- add of a cell in the front ----

       ncel(ip2)=ncel(ip2)+1

       do i=1,3
          icel(ncel(ip2),i,ip2)=iys(2,irai(i),ip1)
          iraf(i)=icel(ncel(ip2),i,ip2)
       enddo 

       do i=1,3
          icel(ncel(ip2),i+3,ip2)=iseg(icel(ncrai,i+3,ip1),6,ip1)
       enddo 
      endif

      return
      end       
