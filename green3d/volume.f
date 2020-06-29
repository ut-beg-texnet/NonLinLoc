      subroutine volume (x1,x2,x3,x4,vol)
      implicit none
      real px,py,pz,qx,qy,qz,rx,ry,rz,wx,wy,wz,vol

c-------------------------------------------------------------------------
c-     Programed by P.S.Lucio and G.C.lambare in Dec/94                  -
c-     Last modification by P.S.Lucio in Jan/05/95                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c-------------------------------------------------------------------------

      real x1(3), x2(3), x3(3), x4(3)

c*******************************************************************
c*      CALCULUS OF THE ALGEBRICS VOLUMES OF THE TETRAHEDRUM       *
c*******************************************************************

      px=x2(1)-x1(1)
      py=x2(2)-x1(2)
      pz=x2(3)-x1(3)

      qx=x2(1)-x3(1)
      qy=x2(2)-x3(2)
      qz=x2(3)-x3(3)

      rx=py*qz-pz*qy
      ry=pz*qx-px*qz
      rz=px*qy-py*qx

      wx=x4(1)-x1(1)
      wy=x4(2)-x1(2)
      wz=x4(3)-x1(3)

      vol=wx*rx+wy*ry+wz*rz

      return
      end 
