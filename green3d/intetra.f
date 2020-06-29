      subroutine intetra (yt,iyt,nxr,nyr,nzr,nir,npr,usrc,xmap,imap,
     &             fi1min,fi1max,fi2min,fi2max) 

c-------------------------------------------------------------------------
c-     programed by P.S.Lucio and G.C.Lambare in Dec/94                  -
c-     Last modification by P.S.Lucio in Jun/20/95                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c-------------------------------------------------------------------------
      implicit none
      real pi2,usrc
      parameter (pi2=3.14159*.5)
c--
      integer npr,nzr,nxr,nyr,nir
      real xmap(npr,nzr,nxr,nyr,nir) 
      integer imap(nzr,nxr,nyr) 

      real xrr(3)
      real yt(36,4),ymap(25,4)
      integer iyt(1,4),iymap(1,4)
      real vol,vol1,vol2,vol3,vol4

      real fi1min,fi1max,fi2min,fi2max
      integer icmin,icmax,jcmin,jcmax,kcmin,kcmax 
      integer ic,jc,kc
c--
c *************************************************************
c *             INTERPOLATION INSIDE A TETRAHEDRUM            *
c *************************************************************

      call volume (yt(1,1),yt(1,2),yt(1,3),yt(1,4),vol)
      if (abs(vol).le.1.e-25) return
      vol=1./vol

      call indice(yt(1,1),yt(1,2),yt(1,3),yt(1,4),
     &            nxr,nyr,nzr,
     &            icmin,icmax,jcmin,jcmax,kcmin,kcmax)

c************************************************************************
c     THE CRITERIUM TO VERIFY THE INTERIORISATION IS BASED ON THE       *
c                      BARYCENTRICS COORDINATES                         *
c************************************************************************
 
      do 99 ic=icmin,icmax
         xrr(1)=ic-1
         do 99 jc=jcmin,jcmax
            xrr(2)=jc-1
            do 99 kc=kcmin,kcmax
               xrr(3)=kc-1

c************************************************************************
c          CALCULUS OF THE ALGEBRICS VOLUMES OF THE TETRAHEDRAS         *
c************************************************************************
 
               call volume (yt(1,1),yt(1,2),yt(1,3),xrr,vol1)
               vol1=vol1*vol
               if (vol1.lt.0.) goto 99

               call volume (yt(1,2),yt(1,3),xrr,yt(1,4),vol2)
               vol2=vol2*vol
               if (vol2.lt.0.) goto 99

               call volume (yt(1,3),xrr,yt(1,4),yt(1,1),vol3)
               vol3=vol3*vol
               if (vol3.lt.0.) goto 99

               call volume (xrr,yt(1,1),yt(1,2),yt(1,4),vol4)
               vol4=vol4*vol
               if (vol4.lt.0.) goto 99

c************************************************************************
c               THE INTERPOLATION INSIDE A TETRAHEDRUM                  *
c************************************************************************
C_AJL      fi1min=amin0(yt(35,1),yt(35,2),yt(35,3),yt(35,4),fi1min+pi2)-pi2
      fi1min=amin1(yt(35,1),yt(35,2),yt(35,3),yt(35,4),fi1min+pi2)-pi2
C_AJL      fi1max=amax0(yt(35,1),yt(35,2),yt(35,3),yt(35,4),fi1max+pi2)-pi2
      fi1max=amax1(yt(35,1),yt(35,2),yt(35,3),yt(35,4),fi1max+pi2)-pi2
C_AJL      fi2min=amin0(yt(36,1),yt(36,2),yt(36,3),yt(36,4),fi2min+pi2)-pi2
      fi2min=amin1(yt(36,1),yt(36,2),yt(36,3),yt(36,4),fi2min+pi2)-pi2
C_AJL      fi2max=amax0(yt(36,1),yt(36,2),yt(36,3),yt(36,4),fi2max+pi2)-pi2
      fi2max=amax1(yt(36,1),yt(36,2),yt(36,3),yt(36,4),fi2max+pi2)-pi2
               call calcray(yt(1,1),iyt(1,1),ymap(1,1),iymap(1,1))
               call calcray(yt(1,2),iyt(1,2),ymap(1,2),iymap(1,2))
               call calcray(yt(1,3),iyt(1,3),ymap(1,3),iymap(1,3))
               call calcray(yt(1,4),iyt(1,4),ymap(1,4),iymap(1,4))
               call ajuste (ic,jc,kc,vol1,vol2,vol3,vol4,ymap,iymap,
     &                      nxr,nyr,nzr,nir,npr,usrc,xmap,imap)
99    continue

      return
      end
