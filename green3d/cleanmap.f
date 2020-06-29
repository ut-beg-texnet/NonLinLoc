      subroutine cleanmap(nxr,nyr,nzr,nir,npr,xmap,imap) 
      implicit none
      integer nxr,nyr,nzr,nir,npr
      integer i,j,k,l,m
      real xmap(npr,nzr,nxr,nyr,nir)
      integer imap(nzr,nxr,nyr)
      write(*,*)'initialisation of the maps in cleanmap.f'
      write(*,*)'present version'
      write(*,*) 'time to 1.e10'
      write(*,*) 'all the others to 0'
 
      do l=1,nir
         do j=1,nyr
            do i=1,nxr
               do k=1,nzr
                  xmap(1,k,i,j,l)=1.E10
                  xmap(2,k,i,j,l)=0.
                  xmap(3,k,i,j,l)=0.
                  do m=4,npr
                     xmap(m,k,i,j,l)=0.
                  enddo
               enddo
            enddo
         enddo
      enddo
      do j=1,nyr
         do i=1,nxr
            do k=1,nzr
               imap(k,i,j)=0
          enddo  
        enddo 
      enddo 


      return
      end
