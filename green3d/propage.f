      subroutine propage (ip1,ip2,ier,iirai,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)

c-------------------------------------------------------------------------
c-     programed by P.S.Lucio and G.C.Lambare in Oct/94                  -
c-     Last modification by P.S.lucio in Jun/19/95                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c-------------------------------------------------------------------------

c      ************************************************************
c      *           PROGRAM TO VERIFIY THE RAY PROPAGATION         *
c      ************************************************************
      implicit none
      integer ip1,ip2,ier,iirai

      integer npn,nsn,ncn,nfront,nd,nk,npar
      real ys(npar,npn,nfront)
      integer iys(nk,npn,nfront)
      double precision dys(nd,npn,nfront)
      integer  iseg(nsn,6,nfront),icel(ncn,6,nfront)
      integer nrai(nfront), ncel(nfront), nseg(nfront)
      integer irai(3), iraf(3), irpop(3), irais(2,3)
      integer isatr(nfront), isats(nfront), isatc(nfront)
c--


      if (ier.ne.0) then 
        iys(2,iirai,ip1)=-1
      else

c      -----  a ray is included to the wavefront  -----

        nrai(ip2)=nrai(ip2)+1
        iys(2,nrai(ip2),ip2)=0
        iys(2,iirai,ip1)=nrai(ip2)
        dys(1,nrai(ip2),ip2)=dys(1,iirai,ip1)
        dys(2,nrai(ip2),ip2)=dys(2,iirai,ip1)
      endif

      return
      end
