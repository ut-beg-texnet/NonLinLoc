      subroutine derive (yv,dy,ier,vm,xm,ym,zm,nxm,nym,nzm,dxm,dym,dzm)

c-------------------------------------------------------------------------
c-     programed by P.S.Lucio and G.C.Lambare in Sep/94                  -
c-     last modification by p.s.lucio in sep/08/94                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c-------------------------------------------------------------------------

C    *****************************************************************
C    *             CALCULUS OF DERIVATIVES FOR RUNGE KUTTA           *
C    *            AND FOR THE TEST OF THE SIZE OF A SEGMENT          *
C    *****************************************************************

C*****************************************************************
C*      THE VARIABLES Y ET DY ARE THE RAY PARAMETERS WHERE:      *
C*                                                               *
C*    Y(1),  Y(2),  Y(3) : X,    Y,    Z                         *
C*    Y(4),  Y(5),  Y(6) : PX,   PY,   PZ                        *
C*    Y(7),  Y(8),  Y(9) : DX1,  DY1,  DZ1                       *
C*    Y(10), Y(11), Y(12): DX2,  DY2,  DZ2                       *
C*    Y(13), Y(14), Y(15): DPX1, DPY1, DPZ1                      *
C*    Y(16), Y(17), Y(18): DPX2, DPY2, DPZ2                      *
C*    Y(19)              : TIME                                  *
C*    Y(20)              : TAU                                   *
C*    Y(21), Y(22), Y(23) : X,    Y,    Z       }                *
C*    Y(24), Y(25), Y(26) : PX,   PY,   PZ      } interpolation  *
C*    Y(27), Y(28), Y(29) : DVX/2, DVY/2, DVZ/2                  *
C*    Y(30)               : V                                    *
C*    Y(31), Y(32), Y(33) :                     } the cross      *
C*    Y(34), Y(35), Y(36) :                     } derivatives    *
C   *    Y(37), Y(38), Y(39) : DXxs,  DYxs,  DPZxs                  *
C   *    Y(40), Y(41), Y(42) : DXys,  DYys,  DPZys                  *
C   *    Y(43), Y(44), Y(45) : DPXxs, DPYxs, DPZxs                  *
C   *    Y(46), Y(47), Y(48) : DPXys, DPYys, DPZys                  *
C*****************************************************************
C*            (THE "D" ARE THEIR DERIVATIVES IN TAU)             *
C*****************************************************************
      implicit none
      integer nxm,nym,nzm
      real vm(nzm,nxm,nym)
      real dxm,dym,dzm
      real xm(2),ym(2),zm(2)
      real yv(48), dy(48)
      integer ier
      real v,dvx,dvy,dvz
      real d2vxx,d2vxy,d2vxz,d2vyy,d2vyz,d2vzz
c--
      ier=0
      call lenteur2(yv(1),yv(2),yv(3),v,dvx,dvy,dvz,
     &              d2vxx,d2vxy,d2vxz,d2vyy,d2vyz,d2vzz,ier,
     &              vm,xm,ym,zm,nxm,nym,nzm,dxm,dym,dzm) 
      if (ier.ne.0) return
c--
      dy(1)=yv(4)
      dy(2)=yv(5)
      dy(3)=yv(6)
      dy(4)=dvx*.5
      dy(5)=dvy*.5
      dy(6)=dvz*.5
      dy(7)=yv(13)
      dy(8)=yv(14)
      dy(9)=yv(15)
      dy(10)=yv(16)
      dy(11)=yv(17)
      dy(12)=yv(18)
      dy(13)=(d2vxx*yv(7)+d2vxy*yv(8)+d2vxz*yv(9))*.5
      dy(14)=(d2vxy*yv(7)+d2vyy*yv(8)+d2vyz*yv(9))*.5
      dy(15)=(d2vxz*yv(7)+d2vyz*yv(8)+d2vzz*yv(9))*.5
      dy(16)=(d2vxx*yv(10)+d2vxy*yv(11)+d2vxz*yv(12))*.5
      dy(17)=(d2vxy*yv(10)+d2vyy*yv(11)+d2vyz*yv(12))*.5
      dy(18)=(d2vxz*yv(10)+d2vyz*yv(11)+d2vzz*yv(12))*.5
      dy(19)=v
      dy(20)=1.
      dy(37)=yv(43)
      dy(38)=yv(44)
      dy(39)=yv(45)
      dy(40)=yv(46)
      dy(41)=yv(47)
      dy(42)=yv(48)
      dy(43)=(d2vxx*yv(37)+d2vxy*yv(38)+d2vxz*yv(39))*.5
      dy(44)=(d2vxy*yv(37)+d2vyy*yv(38)+d2vyz*yv(39))*.5
      dy(45)=(d2vxz*yv(37)+d2vyz*yv(38)+d2vzz*yv(39))*.5
      dy(46)=(d2vxx*yv(40)+d2vxy*yv(41)+d2vxz*yv(42))*.5
      dy(47)=(d2vxy*yv(40)+d2vyy*yv(41)+d2vyz*yv(42))*.5
      dy(48)=(d2vxz*yv(40)+d2vyz*yv(41)+d2vzz*yv(42))*.5
    
c--
      return
      end                                    
