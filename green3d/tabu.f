      program tabu
c--
c        #########################################################
c        This program was developed by P.S. LUCIO and G.C.Lambare
c        in the Centre de Geophysique at Ecole des Mines de Paris
c                and A. Hanyga at University of Bergen
c        This is a tool for the calculus of the 3D Greens function
c        #########################################################
c--
      implicit none
      integer nx,ny,nz
      integer nxr,nyr,nzr,nir,npr
      integer n1,n2,n3,ntot,ntab,ntabnew
C_AJL      parameter (ntab=2000000)
      parameter (ntab=3000000)
      real      tab(ntab)
c******************************************************************
c ...... reading the number of points in x, in y and in z ................

      print*,'dimension of velocity model nx,ny,nz ???'
      read(5,*)nx,ny,nz
      print*,nx,ny,nz
      print*,"dimension of reservoir nxr,nyr,nzr,nir ???"
      read(5,*)nxr,nyr,nzr,nir
      print*,nxr,nyr,nzr,nir
      print*,"number of map parameters npr (<12)"
      read(5,*)npr
      print*,npr

 
      n1=1
      n2=n1+nx*ny*nz
      n3 =n2+nxr*nyr*nzr*nir*npr
      ntot =n3+nxr*nyr*nzr-n1
      ntabnew=ntab-ntot
c--
      print*,'we need a memory of:',ntot
      if (ntot.gt.ntab)then
         print*, 'not enough place change ntab in tabu.f'
         stop
      endif
c--
      call main(nx,ny,nz,nxr,nyr,nzr,nir,npr,
     &          tab(n1),tab(n2),tab(n3),ntabnew,tab(ntot+1))
      stop
      end

