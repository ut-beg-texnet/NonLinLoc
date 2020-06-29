      subroutine cleanfront
     & (npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)
      implicit none
      integer i,j,k

      integer npn,nsn,ncn,nfront,nd,nk,npar   
      real ys(npar,npn,nfront)
      integer iys(nk,npn,nfront)
      double precision dys(nd,npn,nfront) 
      integer  iseg(nsn,6,nfront),icel(ncn,6,nfront) 
      integer nrai(nfront), ncel(nfront), nseg(nfront) 
      integer irai(3), iraf(3), irpop(3), irais(2,3) 
      integer isatr(nfront), isats(nfront), isatc(nfront)
 
      do k=1,nfront
        do j=1,npn
          do i=1,npar
             ys(i,j,k) = 0
          enddo
          do i=1,nd
             dys(i,j,k) = 0
          enddo
          do i=1,nk
             iys(i,j,k) = 0
          enddo
        enddo
      enddo

      return
      end
