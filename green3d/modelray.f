       subroutine modelray(xsrc,
     & fi1min,fi2min,fi1max,fi2max,dxmin2,dpmin2,dtemps)
      implicit none
      real xsrc(3)
      real fi1min,fi2min,fi1max,fi2max,dxmin2,dpmin2,dtemps
      real dxmin,dpmin
c   ***************************************************************
c   *             INPUT THE PARAMETERS FOR THE 3D MAPS            *
c   ***************************************************************
      print*,'source position xs, ys, zs???'
      read(*,*)xsrc(1),xsrc(2),xsrc(3)
      print*,  xsrc(1),xsrc(2),xsrc(3)

      print*,'angle / vertical direction x fi1min,fi1max ??'
      read(*,*)fi1min,fi1max
      print*,fi1min,fi1max
      print*,'angle / vertical direction y fi2min,fi2max ??'
      read(*,*)fi2min,fi2max
      print*,fi2min,fi2max
      print*,'precision en x dxmin ?'
      read(*,*)dxmin
      dxmin2=dxmin**2
      print*,dxmin,'m'
      print*,'precision en p dpmin ?'
      read(*,*)dpmin
      dpmin2=dpmin**2
      print*,dpmin,'s/m'
      print*,'pas du schemas de Runge ?'
      read(*,*)dtemps
      print*,dtemps

      return
      end
