      subroutine deter (x1,x2,x3,det)

c-------------------------------------------------------------------------
c-     Programed by P.S.Lucio and G.C.Lambare in Jan/95                  -
c-     Last modification by P.S.Lucio in Jan/10/95                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c-------------------------------------------------------------------------

      implicit none
      real x1(3), x2(3), x3(3),det,d1,d2,d3,d4,d5,d6

c*******************************************************************
c*           CALCULUS OF THE DETERMINANT OF A 3X3 MATRIX           *
c*******************************************************************

      d1=x1(1)*x2(2)*x3(3)
      d2=x1(2)*x2(3)*x3(1)
      d3=x1(3)*x2(1)*x3(2) 
      d4=x1(3)*x2(2)*x3(1)
      d5=x1(2)*x2(1)*x3(3) 
      d6=x1(1)*x2(3)*x3(2)

      det=d1+d2+d3-(d4+d5+d6)

      return
      end 
