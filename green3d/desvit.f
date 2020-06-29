      subroutine desvit (nxm,nym,nzm,nxr,nyr,nzr,
     &       vm,xm,ym,zm,dxm,dym,dzm,vit,targ_orient)

c-------------------------------------------------------------------------
c-     Programed by P.S.Lucio and G.C.Lambare in Jan/95                  -
c-     Last modification by P.S.Lucio in Jun/26/95                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c-------------------------------------------------------------------------

c ********************************************************************
c *         PROGRAM OF THE CALCULUS OF THE SQUARED SLOWNESS          *
c *       INTERPOLATION OF THE VELOCITY MODEL IN THE RESERVOIR       *
c *               AND DESIGN OF THE VELOCITY MODEL                   *
c ********************************************************************
c--      
      implicit none
      integer nxm,nym,nzm
      integer nxr,nyr,nzr
      real dxm,dym,dzm
      real vm(nzm,nxm,nym)
      real xm(2),ym(2),zm(2)
      real vit(nzr)
      real targ_orient(3,7)
      character *80 nomvit 
      integer i,j,k,ier
      real xr,yr,zr,v,dvx,dvy,dvz
      real d2vxx,d2vxy,d2vxz,d2vyy,d2vyz,d2vzz

            
      print*,'name target velocity file'
      read(*,'(a80)')nomvit
      call centra(nomvit,i,j)
      write(*,*)nomvit(i:j)
      open(18,file=nomvit(i:j),access='direct',form=
     &'unformatted',recl=4*nzr)
 

      do 1 i=1,nxr
        do 1 j=1,nyr
          do 2 k=1,nzr
        xr=targ_orient(1,1)+targ_orient(1,2)*(i-1)
     &                     +targ_orient(1,3)*(j-1)
     &                     +targ_orient(1,4)*(k-1)
           yr=targ_orient(2,1)+targ_orient(2,2)*(i-1)
     &                        +targ_orient(2,3)*(j-1)
     &                        +targ_orient(2,4)*(k-1)
             zr=targ_orient(3,1)+targ_orient(3,2)*(i-1)
     &                          +targ_orient(3,3)*(j-1)
     &                          +targ_orient(3,4)*(k-1)
            call lenteur2(xr,yr,zr,v,dvx,dvy,dvz,
     &            d2vxx,d2vxy,d2vxz,d2vyy,d2vyz,d2vzz,ier,
     &            vm,xm,ym,zm,nxm,nym,nzm,dxm,dym,dzm)
            if (ier.ne.0) v=1.E20
            vit(k)=1./sqrt(v)
2         continue
          write(18,rec=(j-1)*nxr+i) (vit(k),k=1,nzr)
1     continue

      close(18)

      return
      end
