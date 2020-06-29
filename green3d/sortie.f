      subroutine sortie(nxr,nyr,nzr,nir,npr,xmap,imap)
c-------------------------------------------------------------------------
c-     Programed by P.S.Lucio and G.C.Lambare in Dec/94                  -
c-     Last modification by P.S.Lucio in Jun/22/95                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c-------------------------------------------------------------------------

c ********************************************************************
c *              PROGRAM OF TRI AND OUTPUT OF THE RESULTS            *
c ********************************************************************
c--
      implicit none
      integer npr,nzr,nxr,nyr,nir
      real  xmap(npr,nzr,nxr,nyr,nir)
      integer imap(nzr,nxr,nyr)
      integer i,j,k,l,in,ix,ipara

      character*80 ficmap
5     format(a80)
 
      write(*,*)'name of the map of the number of arrivals'
         read(*,5)ficmap
         call centra(ficmap,in,ix)
         write(*,*)ficmap(in:ix)
         open(7,file=ficmap(in:ix),access='direct',form='unformatted',
     &      recl=4*nzr*nxr*nyr)            
      write(7,rec=1) (((float(imap(k,i,j)),k=1,nzr)
     &             ,i=1,nxr),j=1,nyr)
      close(7)

      write(*,*)'names of the  ',npr,' results file'
      do 1 ipara=1,npr
         read(*,5)ficmap
         call centra(ficmap,in,ix)
         write(*,*)ficmap(in:ix)
         open(7,file=ficmap(in:ix),access='direct',form='unformatted',
     &      recl=4*nir*nzr*nxr*nyr)            
         write(7,rec=1) ((((xmap(ipara,k,i,j,l),k=1,nzr)
     &               ,i=1,nxr),j=1,nyr),l=1,nir)
         close(7)
1     continue

      return
      end
