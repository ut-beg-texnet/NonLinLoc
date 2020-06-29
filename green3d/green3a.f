c****************************************************************************
c****************************************************************************

      subroutine green3a(nxm,nym,nzm,dxm,dym,dzm,xm,ym,zm,vm,
     &     nxr,nyr,nzr,nir,npr,targ_orient,xmap,imap,
     &     xsrc,fi1min,fi2min,fi1max,fi2max,dxmin2,dpmin2,dtemps,
     &     ntab,tab)
c ************************************************************************
c      INPUT
c                             VELOCITY MODEL
c      nxm,nym,nzm : size of the velocity model defined by weight of
c                    cubic cardinal Bsplines in square slowness
c      dxm,dym,dzm : step of the Bplines knots           in m
c      xm(2) = (xmin+dxm,xmax-dxm)           in m
c      ym(2) = (ymin+dym,ymax-dym)           in m
c      zm(2) = (zmin+dzm,zmax-dzm)           in m
c      vm(nzm,nxm,nym) : weights of the Bsplines (Attention to the order
c                        of dimensions) in (s/m)**2
c                        Attention the interpolation is not valid in the
c                        first layer of the cube !
c
c                             TARGET
c      nxr,nyr,nzr : size of the target
c      nir         : maximum number of stored arrival
c                   (by default it always keep the first arrival traveltime
c                    see in ajuste.f for any modification of this default)
c      npr         : number of parameters for which to compute maps
c                    it should be less than 27 and greater or equal to 3
c                    see in ajuste for any modification
c
c                             RAY TRACING
c      xsrc(3)         : position of the source (in m)
c      fi1min,fi2min,fi1max,fi2max : initial angular aperture (in radian)
c                              see (Lucio et al, 1996)
c      dxmin2,dpmin2   = dxmax**2,dpmax**2  see (Lucio et al, 1996)
c                        it gives some idea of the precision of the ray
c                        field sampling (in m**2 and (s/m)**2)
c                        example (5**2 , (5E-E06)**2)
c                        increase this number in order to go faster
c      dtemps          : travel time step (in s) 
c      ntab            : size of the memory allocation 
c      tab(ntab)       : memory allocation
c                        a minimum value is required (see green3a.f)
c                        if your computation overpass the storage
c                        allocation the computation continues but 
c                        you get holes in your maps
c                        (the code advise you in such a case)
c
c
c      OUTPUT
c      xmap(npr,nzr,nxr,nyr,nir) : maps 
c          Attention to the order of the indices
c                                  the first indice denotes the parameter
c                       xmap( 1) : temps
c                       xmap( 2) : ampl see (lucio et al, 1996)
c                       xmap( 3) : KMAH
c                       xmap( 4) : px
c                       xmap( 5) : py
c                       xmap( 6) : pz
c                       xmap( 7) : phi1 see (lucio et al, 1996)
c                       xmap( 8) : phi2
c                       xmap( 9) : dx/dfi1
c                       xmap( 0) : dy/dfi1
c                       xmap(11) : dz/dfi1
c                       xmap(12) : dx/dfi2
c                       xmap(13) : dy/dfi2
c                       xmap(14) : dz/dfi2
c                       xmap(15) : dpx/dfi1
c                       xmap(16) : dpy/dfi1
c                       xmap(17) : dpz/dfi1
c                       xmap(18) : dpx/dfi2
c                       xmap(19) : dpy/dfi2
c                       xmap(20) : dpz/dfi2
c                       xmap(21) : dpx/dxs
c                       xmap(22) : dpy/dxs
c                       xmap(23) : dpz/dxs
c                       xmap(24) : dpx/dys
c                       xmap(25) : dpy/dys
c                       xmap(26) : dpz/dys
c        distances are in m, time in s, angles in radians, ...
c        Attention the arrivals are not ordered. The only thing than you
c        can expect is the code fills the maps by incrementing the last
c        indice.
c        defaults values are fixed in cleanmap.f (see for modifications)
c        When nir is not sufficient for recovering all the arrival the first
c        arrival is always preserved (see ajuste.f for modification)
c        imap(nzr,nxr,nyr)      : total number of arrivals by points
c ************************************************************************
c--
c-------------------------------------------------------------------------
c-     Programed by P.S.Lucio and G.C.Lambare in Apr/95                  -
c-     Last modification by P.S.Lucio in Jul/10/95                       -
c-------------------------------------------------------------------------
c--
      implicit none
      integer nfront,nd,nk,npar
      parameter (nfront=2,nd=2,nk=2,npar=48)
      integer nrai(nfront), ncel(nfront), nseg(nfront)
      integer irai(3), iraf(3), irpop(3), irais(2,3)
      integer isatr(nfront), isats(nfront), isatc(nfront)
      integer nxm,nym,nzm
      integer nxr,nyr,nzr,nir,npr
      real xm(2),ym(2),zm(2)
      real vm(nzm,nxm,nym)
      real dxm,dym,dzm

      real xsrc(3)

      real xmap(npr,nzr,nxr,nyr,nir)
      integer imap(nzr,nxr,nyr)
      real targ_orient(3,7)
      real fi1min,fi2min,fi1max,fi2max,dxmin2,dpmin2,dtemps

      integer npn,ncn,nsn
      integer n1,n2,n3,n4,n5,ntot,ntab
      real tab(ntab)
      
      npn=ntab/(nfront*(npar+nk+2*nd+30))
      ncn=2*npn
      nsn=3*npn
      write(*,*)'greena, maximum size for isochron mesh'
      write(*,*)'rays :',npn,' triangles :',ncn,' segments :',nsn
      if (npn.lt.10) then
         write(*,*)'NOT ENOUGH MEMORY !'
         write(*,*)'increase ntab in tabu.f'
         stop
      endif

      n1=1
      n2=n1+npar*npn*nfront
      n3=n2+nk  *npn*nfront
      n4=n3+2*nd*npn*nfront
      n5=n4+nsn*6*nfront
      ntot=n5+ncn*6*nfront-n1
      
c--
      print*,'we use an extra memory of:',ntot
      print*,'for storing the wavefront'
      if (ntot.gt.ntab)then
         print*, 'error in greena.f'
         stop
      endif
c--


      call green3(nxm,nym,nzm,dxm,dym,dzm,xm,ym,zm,vm,
     &     nxr,nyr,nzr,nir,npr,targ_orient,xmap,imap,
     &     xsrc,fi1min,fi2min,fi1max,fi2max,dxmin2,dpmin2,dtemps,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & tab(n1),tab(n2),tab(n3),tab(n4),tab(n5),nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)


      return
      end
