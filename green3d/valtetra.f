      subroutine valtetra (ir,ipoint,yt,iyt,targ_orient,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg, 
     & irai,iraf,irpop,irais,isatr,isats,isatc)
c-------------------------------------------------------------------------
c-     Programed by P.S.Lucio and G.C.Lambare in Jan/95                  -
c-     Last modification by P.S.Lucio in Jun/21/95                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c-------------------------------------------------------------------------

c **************************************
c *  yt( 1) : X in target coordinate
c *  yt( 2) : Y in target coordinate
c *  yt( 3) : Z in target coordinate
c *  yt( 4) : px
c *  yt( 5) : py
c *  yt( 6) : pz
c *  yt( 7) : temps
c *  yt( 8) : dx/dfi1
c *  yt( 9) : dy/dfi1
c *  yt(10) : dz/dfi1
c *  yt(11) : dx/dfi2
c *  yt(12) : dy/dfi2
c *  yt(13) : dz/dfi2
c *  yt(14) : dpx/dfi1
c *  yt(15) : dpy/dfi1
c *  yt(16) : dpz/dfi1
c *  yt(17) : dpx/dfi2
c *  yt(18) : dpy/dfi2
c *  yt(19) : dpz/dfi2
c *  yt(20) : dx/dxs
c *  yt(21) : dy/dxs
c *  yt(22) : dz/dxs
c *  yt(23) : dx/dys
c *  yt(24) : dy/dys
c *  yt(25) : dz/dys
c *  yt(26) : dpx/dxs
c *  yt(27) : dpy/dxs
c *  yt(28) : dpz/dxs
c *  yt(29) : dpx/dys
c *  yt(30) : dpy/dys
c *  yt(31) : dpz/dys
c *  yt(32) : dpx/dtau
c *  yt(33) : dpy/dtau
c *  yt(34) : dpz/dtau
c *  yt(35) : phi1
c *  yt(36) : phi2
c * iyt(1)  : KMAH
c***************************************************************************
c*          ATTRIBUTION OF VALUES AT THE VERTEX OF THE TETRAHEDRUM         *
c***************************************************************************
      implicit none
      integer ir,ipoint,i
      integer npn,nsn,ncn,nfront,nd,nk,npar 
      real ys(npar,npn,nfront) 
      integer iys(nk,npn,nfront) 
      double precision dys(nd,npn,nfront) 
      integer  iseg(nsn,6,nfront),icel(ncn,6,nfront) 
      integer nrai(nfront), ncel(nfront), nseg(nfront) 
      integer irai(3), iraf(3), irpop(3), irais(2,3) 
      integer isatr(nfront), isats(nfront), isatc(nfront)
c--
      real yt(36)
      integer iyt(1)
      real targ_orient(3,7)
      real xi0,yi0,zi0,xi1,yi1,zi1,dxi,dyi,dzi

c*************************************************************************

     
c     --- POSITION IN TARGET COORDINATES -------------
      xi0=targ_orient(1,1) 
      yi0=targ_orient(2,1) 
      zi0=targ_orient(3,1) 
      xi1=ys(20+1,ipoint,ir)-xi0
      yi1=ys(20+2,ipoint,ir)-yi0
      zi1=ys(20+3,ipoint,ir)-zi0
      do i=1,3
        dxi=targ_orient(i,1+4)
        dyi=targ_orient(i,2+4)
        dzi=targ_orient(i,3+4)
        yt(i)= dxi*xi1+dyi*yi1+dzi*zi1
      enddo
c     --- SLOWNESS -----------------------------------
      do i=1,3
        yt(i+3)=ys(i+23,ipoint,ir)
      enddo
c     --- TRAVEL TIME --------------------------------
        yt(7)=ys(19,ipoint,ir)
c     --- dX/dfi1 ------------------------------------
      do i=1,3
        yt(i+7)=ys(i+6,ipoint,ir)
      enddo
c     --- dX/dfi2 ------------------------------------
      do i=1,3
        yt(i+10)=ys(i+9,ipoint,ir)
      enddo
c     --- dP/dfi1 ------------------------------------
      do i=1,3
        yt(i+13)=ys(i+12,ipoint,ir)
      enddo
c     --- dP/dfi2 ------------------------------------
      do i=1,3
        yt(i+16)=ys(i+15,ipoint,ir)
      enddo
c     --- dP/dtau -----------------------------------
      do i=1,3
        yt(i+31)=ys(i+26,ipoint,ir)
      enddo
c     --- dX/dxs ------------------------------------
      do i=1,3
        yt(i+19)=ys(i+36,ipoint,ir)
      enddo
c     --- dX/dys ------------------------------------
      do i=1,3
        yt(i+22)=ys(i+39,ipoint,ir)
      enddo
c     --- dP/dxs ------------------------------------
      do i=1,3
        yt(i+25)=ys(i+42,ipoint,ir)
      enddo   
c     --- dP/dys ------------------------------------
      do i=1,3
        yt(i+28)=ys(i+45,ipoint,ir)
      enddo
c     --- PHI1 --------------------------------------
        yt(35)=dys(1,ipoint,ir)
        yt(36)=dys(2,ipoint,ir)
c     --- KMAH --------------------------------------
       iyt(1)=iys(1,ipoint,ir)
        
      return
      end
