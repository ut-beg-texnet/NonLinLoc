      subroutine modelvit(vm,xm,ym,zm,dxm,dym,dzm,nxm,nym,nzm)

c-------------------------------------------------------------------------
c-     Programed by P.S.Lucio and G.C.Lambare in Oct/94                  -
c-     Last modification by P.S.Lucio in Dec/29/94                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c-------------------------------------------------------------------------
      implicit none 
      integer nxm,nym,nzm
      real vm(nzm,nxm,nym)
      real xm(2),ym(2),zm(2)
      character *80 nomfic
      integer i,j,k
      real xmin,ymin,zmin
      real dxm,dym,dzm


C******************************************************************
C*          READING THE MEDIUM (DATA FILE OF THE BORDS)           *
C******************************************************************

      print*,'origine et pas  du modele de vitesse xmin ddx ?'
      read(5,*)xmin, dxm
      print*,xmin, dxm
      print*,'origine et pas  du modele de vitesse ymin dym ?'
      read(5,*)ymin, dym
      print*,ymin, dym
      print*,'origine et pas  du modele de vitesse zmin dzm ?'
      read(5,*)zmin, dzm
      print*,zmin, dzm

c ****************************************************************
c           READING THE DATA FILE OF THE 3D VELOCITY MODEL
c ****************************************************************
c--
      xm(1)=xmin+dxm
      xm(2)=xmin+dxm*(nxm-2)
      ym(1)=ymin+dym
      ym(2)=ymin+dym*(nym-2)
      zm(1)=zmin+dzm
      zm(2)=zmin+dzm*(nzm-2)

Cv******************************************************************
C*          INITIALISATION OF THE VELOCITY IN THE MEDIUM          *
C******************************************************************
 
      print*,'nom fic (weight of square slowness Bspline'
      read(*,'(a80)')nomfic
      print*,nomfic
      open(7,file=nomfic,access='direct'
     &       ,form='unformatted',recl=4*nxm*nzm)

C******************************************************************
C*          INITIALISATION OF THE VELOCITY IN THE MEDIUM          *
C******************************************************************
         do 4 j=1,nym
           read(7,rec=j) ((vm(k,i,j),k=1,nzm),i=1,nxm)
           do 4 i=1,nxm
              do 4 k=1,nzm
                 vm(k,i,j)=1./vm(k,i,j)**2
4        continue
      close(7)

      return
      end
