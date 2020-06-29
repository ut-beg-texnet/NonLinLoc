       subroutine calcray(yt,iyt,ymap,iymap)
       real yt(36),ymap(25)
       integer iyt(1),iymap(1)
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
c*****************
c *  ymap( 1) : temps
c *  ymap( 2) : 1/ampl
c *  ymap( 3) : px
c *  ymap( 4) : py
c *  ymap( 5) : pz
c *  ymap( 6) : phi1
c *  ymap( 7) : phi2
c *  ymap( 8) : dx/dfi1
c *  ymap( 9) : dy/dfi1
c *  ymap(10) : dz/dfi1
c *  ymap(11) : dx/dfi2
c *  ymap(12) : dy/dfi2
c *  ymap(13) : dz/dfi2
c *  ymap(14) : dpx/dfi1
c *  ymap(15) : dpy/dfi1
c *  ymap(16) : dpz/dfi1
c *  ymap(17) : dpx/dfi2
c *  ymap(18) : dpy/dfi2
c *  ymap(19) : dpz/dfi2
c *  ymap(20) : dpx/dxs
c *  ymap(21) : dpy/dxs
c *  ymap(22) : dpz/dxs
c *  ymap(23) : dpx/dys
c *  ymap(24) : dpy/dys
c *  ymap(25) : dpz/dys
c *  iymap(1) : KMAH
c*****************
       ymap(1)=yt(7)
c           --- calculus of the geometrical spreading ---
       call deter (yt(4),yt(8),yt(11),pond)
       ymap(2)=pond/abs(sin(yt(35)))
       ymap(3)=yt(4)
       ymap(4)=yt(5)
       ymap(5)=yt(6)
       ymap(6)=yt(35)
       ymap(7)=yt(36)
       do 1 i=8,19
          ymap(i)=yt(i)
1      continue
       call extirp(ymap(20),yt)
       iymap(1)=iyt(1)
       return
       end



