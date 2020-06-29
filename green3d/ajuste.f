      subroutine ajuste (ic,jc,kc,vol1,vol2,vol3,vol4,ymap,iymap,
     &                   nxr,nyr,nzr,nir,npr,usrc,xmap,imap)

c-------------------------------------------------------------------------
c-     Programed by P.S.Lucio and G.C.Lambare in Jan/95                  -
c-     Last modification by P.S.Lucio in Jun/21/95                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c-------------------------------------------------------------------------
c  **********************************************************************
c  *    THIS INTERPOLATION IS BASED ON THE BARYCENTRICS COORDINATES     *
c  *           (Curves and Surfaces for CAGD: Gerald FARIN)             *
c  **********************************************************************
c**************************************************************************
c INPUT  : ic, jc, kc - index for the position on the reservoir           *
c                     -  the variables for interpolation                  *
c                     -  the weigths                                      *
c OUTPUT : the variables in for the maps                                  *
c**************************************************************************
      implicit none
c--
      integer npr,nzr,nxr,nyr,nir
      real  xmap(npr,nzr,nxr,nyr,nir),ymap(25,4)
      integer  imap(nzr,nxr,nyr),iymap(1,4)
      real tempsprov,amplprov,vol1,vol2,vol3,vol4
      integer i,ic,jc,kc,icur
      real usrc
c--

c     --- Temps ----------------------------------------------
      tempsprov=vol1*ymap(1,4)+vol2*ymap(1,1)
     &         +vol3*ymap(1,2)+vol4*ymap(1,3)

      imap(kc,ic,jc)=imap(kc,ic,jc)+1
      icur=imap(kc,ic,jc)

      if (icur.gt.nir) then
         icur=nir
         if (tempsprov.gt.xmap(1,kc,ic,jc,icur) ) return
      endif

      xmap(1,kc,ic,jc,icur)=tempsprov
      amplprov=(vol1*ymap(2,4)+vol2*ymap(2,1)+
     &             vol3*ymap(2,2)+vol4*ymap(2,3)) 
c     --- KMAH ---------------------------------------------------
      xmap(3,kc,ic,jc,icur)=iymap(1,1)
      do i=2,4
      if ((amplprov*ymap(2,i)).gt.0.) xmap(3,kc,ic,jc,icur)=iymap(1,i)
      enddo

c     --- Amplitude ----------------------------------------------
      if (abs(amplprov).lt.1.e-25) amplprov=1.e-25
      xmap(2,kc,ic,jc,icur)=sqrt(usrc/abs(amplprov))


c     --- Ray parameters -----------------------------------------
      do i=3,npr-1
         xmap(i+1,kc,ic,jc,icur)=vol1*ymap(i,4)+vol2*ymap(i,1)
     &                         +vol3*ymap(i,2)+vol4*ymap(i,3)
      enddo

c ***********************************************************
c *  xmap( 1) : temps
c *  xmap( 2) : ampl
c *  xmap( 3) : KMAH
c *  xmap( 4) : px
c *  xmap( 5) : py
c *  xmap( 6) : pz
c *  xmap( 7) : phi1
c *  xmap( 8) : phi2
c *  xmap( 9) : dx/dfi1
c *  xmap( 0) : dy/dfi1
c *  xmap(11) : dz/dfi1
c *  xmap(12) : dx/dfi2
c *  xmap(13) : dy/dfi2
c *  xmap(14) : dz/dfi2
c *  xmap(15) : dpx/dfi1
c *  xmap(16) : dpy/dfi1
c *  xmap(17) : dpz/dfi1
c *  xmap(18) : dpx/dfi2
c *  xmap(19) : dpy/dfi2
c *  xmap(20) : dpz/dfi2
c *  xmap(21) : dpx/dxs
c *  xmap(22) : dpy/dxs
c *  xmap(23) : dpz/dxs
c *  xmap(24) : dpx/dys
c *  xmap(25) : dpy/dys
c *  xmap(26) : dpz/dys
c *******************************************************************
      return 
      end
