       subroutine infronti (ip,ncrai,ini,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)

c-------------------------------------------------------------------------
c-     Programed by P.S.Lucio and G.C.Lambare in Oct/94                  -
c-     Last modification by P.S.Lucio in Dec/12/94                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c-------------------------------------------------------------------------

C ******************************************************************
C *                  INITIALISATION OF A WAVEFRONT                 *
C ******************************************************************
       implicit none
       integer ip,ncrai,ini

      integer npn,nsn,ncn,nfront,nd,nk,npar
      real ys(npar,npn,nfront)
      integer iys(nk,npn,nfront)
      double precision dys(nd,npn,nfront)
      integer  iseg(nsn,6,nfront),icel(ncn,6,nfront)
      integer nrai(nfront), ncel(nfront), nseg(nfront)
      integer irai(3), iraf(3), irpop(3), irais(2,3)
      integer isatr(nfront), isats(nfront), isatc(nfront)

       integer i,j,k,isomme

C ******************************************************************
C * INPUT  ncrai :                                                 *
C *        ip    : index of wavefront                              *
C *        isatr(ip) - index of the wavefront saturation           *
C *        isatr(ip) - index of the wavefront saturation           *
C *        isatr(ip) - index of the wavefront saturation           *
C *                                                                *
C * OUTPUT ini:  indicator of the end of the wavefront             *
C *            0 :  there is still to be calculated                *
C *            1 :  that's all!                                    *
C *            2 :  saturated wavefront!                           *
C ******************************************************************

C********************************************************************
C*    irai: number of the point of the cell                         *
C     irpop: number of opposites points to each one of the segments *
C*    irais: number of the two "rays" of each one segment           *
C********************************************************************

      ini=0
      isomme=0
      ncrai=ncrai+1

      if (ncrai.gt.ncel(ip)) then 
        ini=1 
        return 
      endif

      do j=1,3
         irai(j)=icel(ncrai,j,ip)
         isomme=isomme+irai(j)
      enddo

      do i=1,2
        do j=1,3
          irais(i,j)=iseg(icel(ncrai,j+3,ip),i,ip)       
        enddo
      enddo

      do k=1,3
        irpop(k)=isomme-(irais(1,k)+irais(2,k))
      enddo

      return
      end
