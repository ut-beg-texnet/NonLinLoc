      subroutine derivebis (ip,ipoint,dey,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)

c---------------------------------------------------------------------
c-     programed by p.s.lucio and g.c.lambare in sep/94               
c-     last modification by p.s.lucio in sep/08/94                   
c---------------------------------------------------------------------

C    *****************************************************************
C    *             CALCULUS OF DERIVATIVES FOR RUNGE KUTTA           *
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
C*****************************************************************
C*            (THE "D" ARE THEIR DERIVATIVES IN TAU)             *
C*****************************************************************
      implicit none
      integer ip,ipoint

      integer npn,nsn,ncn,nfront,nd,nk,npar
      real ys(npar,npn,nfront)
      integer iys(nk,npn,nfront)
      double precision dys(nd,npn,nfront)
      integer  iseg(nsn,6,nfront),icel(ncn,6,nfront)
      integer nrai(nfront), ncel(nfront), nseg(nfront)
      integer irai(3), iraf(3), irpop(3), irais(2,3)
      integer isatr(nfront), isats(nfront), isatc(nfront)


      real dx,dy,dz,d2vxx,d2vxy,d2vxz,d2vyy,d2vyz,d2vzz
      real dey(48)
c--
      dey(1)=ys(4,ipoint,ip)
      dey(2)=ys(5,ipoint,ip)
      dey(3)=ys(6,ipoint,ip)
      dey(4)=ys(27,ipoint,ip)
      dey(5)=ys(28,ipoint,ip)
      dey(6)=ys(29,ipoint,ip)
      dey(7)=ys(13,ipoint,ip)
      dey(8)=ys(14,ipoint,ip)
      dey(9)=ys(15,ipoint,ip)
      dey(10)=ys(16,ipoint,ip)
      dey(11)=ys(17,ipoint,ip)
      dey(12)=ys(18,ipoint,ip)
         dx=ys(7,ipoint,ip)
         dy=ys(8,ipoint,ip)
         dz=ys(9,ipoint,ip)
         d2vxx=ys(31,ipoint,ip)
         d2vxy=ys(32,ipoint,ip)
         d2vxz=ys(33,ipoint,ip)
         d2vyy=ys(34,ipoint,ip)
         d2vyz=ys(35,ipoint,ip)
         d2vzz=ys(36,ipoint,ip)
      dey(13)=d2vxx*dx+d2vxy*dy+d2vxz*dz
      dey(14)=d2vxy*dx+d2vyy*dy+d2vyz*dz
      dey(15)=d2vxz*dx+d2vyz*dy+d2vzz*dz
         dx=ys(10,ipoint,ip)
         dy=ys(11,ipoint,ip)
         dz=ys(12,ipoint,ip)
      dey(16)=d2vxx*dx+d2vxy*dy+d2vxz*dz
      dey(17)=d2vxy*dx+d2vyy*dy+d2vyz*dz
      dey(18)=d2vxz*dx+d2vyz*dy+d2vzz*dz
      dey(19)=ys(30,ipoint,ip)
      dey(20)=1.

      dey(37)=ys(43,ipoint,ip)
      dey(38)=ys(44,ipoint,ip)
      dey(39)=ys(45,ipoint,ip)
      dey(40)=ys(46,ipoint,ip)
      dey(41)=ys(47,ipoint,ip)
      dey(42)=ys(48,ipoint,ip)
         dx=ys(37,ipoint,ip)
         dy=ys(38,ipoint,ip)
         dz=ys(39,ipoint,ip)
      dey(43)=d2vxx*dx+d2vxy*dy+d2vxz*dz
      dey(44)=d2vxy*dx+d2vyy*dy+d2vyz*dz
      dey(45)=d2vxz*dx+d2vyz*dy+d2vzz*dz
         dx=ys(40,ipoint,ip) 
         dy=ys(41,ipoint,ip) 
         dz=ys(42,ipoint,ip) 
      dey(46)=d2vxx*dx+d2vxy*dy+d2vxz*dz
      dey(47)=d2vxy*dx+d2vyy*dy+d2vyz*dz
      dey(48)=d2vxz*dx+d2vyz*dy+d2vzz*dz

c--
      return
      end                                  
