      subroutine lenteur2(x,y,z,vi,dvx,dvy,dvz,
     &                    d2vxx,d2vxy,d2vxz,d2vyy,d2vyz,d2vzz,ier,
     &                    vm,xm,ym,zm,nxm,nym,nzm,dxm,dym,dzm)

c-------------------------------------------------------------------------
c-     Programed by P.S.Lucio and G.C.Lambare in Jan/94                  -
c-     Based on the "de Boor" algorithm - Carl de Boor                   -
c-                                        A pratical Guide to Splines    -
c-     last modification by P.S.Lucio in Dec/12/94                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c-------------------------------------------------------------------------
c--
      implicit none
      integer nxm,nym,nzm
      real dxm,dym,dzm
      real vm(nzm,nxm,nym)
      real xm(2),ym(2),zm(2)
      real x,y,z
      real vi,dvz,dvx,dvy,d2vzz,d2vxx,d2vyy,d2vxz,d2vxy,d2vyz
      integer ier
c--
	  call bspl3d2(z,x,y,nzm,nxm,nym,vm,zm,xm,ym,dzm,dxm,dym,
     & vi,dvz,dvx,dvy,d2vzz,d2vxx,d2vyy,d2vxz,d2vxy,d2vyz,ier)

      return
      end
c--
C-------------------------------------------------------------------------------
C               3-D CARDINAL CUBIC B-SPLINE INTERPOLATOR
C-------------------------------------------------------------------------------

      subroutine bspl3d2(x,y,z,nxm,nym,nzm,vm,xm,ym,zm,dx,dy,dz,
     &w,dxw,dyw,dzw,d2xxw,d2yyw,d2zzw,d2xyw,d2yzw,d2zxw,ier)

C----------------------------------------------------------------------
C This is a different interpolator than 'FNOYAU+FINTERP'
C because it manages the borders of the interpolated region exactly 
C in the same way as what is currently used in 1D and 2D.   
C----------------------------------------------------------------------
C The algorithm is entirely formulated as embedded 1D Bspline interpo-
C lations (see subroutines lspl1d*()).
C----------------------------------------------------------------------
C Cost:                   multadd, add, mult
C lspl1d2   -> 21 times :     105, 189,   42
C lspl1d1   ->  6 times :      30,  42,    6
C lspl1d0   ->  7 times :      21,  49,    7
C autres... ->          :       0,   6,   24
C Total :                     156, 286,   79  ====>  442 add, 235 mult
C Old (FINTERP) :             412,  12,   15  ====>  424 add, 427 mult
C-----------------------------------------------------------------------

	implicit none
	integer nxm,nym,nzm,ier
	real*4 x,y,z,vm(nxm,nym,nzm),xm(2),ym(2),zm(2),dx,dy,dz
        real*4 w,dxw,dyw,dzw,d2xxw,d2yyw,d2zzw,d2xyw,d2yzw,d2zxw

C---------------------------------------------------------------------
C x must be in [xm(1),xm(2)], otherwise, nothing is done and 'ier'
C is set to 99, else 'ier' is set to 0.
C vm(*) contains spline weights at nodes. Its contains
C (xm(2)-xm(1))/dx+3 elements; vm(1) corresponds to position xm(1)-dx.
C Same applies for each dimension.
C---------------------------------------------------------------------

C intermediate nodes after x-interpolation
	real*4 w11,w12,w13,w14,w21,w22,w23,w24
        real*4 w31,w32,w33,w34,w41,w42,w43,w44
	real*4 dxw11,dxw12,dxw13,dxw14,dxw21,dxw22,dxw23,dxw24
        real*4 dxw31,dxw32,dxw33,dxw34,dxw41,dxw42,dxw43,dxw44
	real*4 d2xxw11,d2xxw12,d2xxw13,d2xxw14
        real*4 d2xxw21,d2xxw22,d2xxw23,d2xxw24
        real*4 d2xxw31,d2xxw32,d2xxw33,d2xxw34
        real*4 d2xxw41,d2xxw42,d2xxw43,d2xxw44
C intermediate nodes after x- then y-interpolation
	real*4 w1,w2,w3,w4
        real*4 dyw1,dyw2,dyw3,dyw4
	real*4 dxw1,dxw2,dxw3,dxw4
        real*4 d2yyw1,d2yyw2,d2yyw3,d2yyw4
        real*4 d2xyw1,d2xyw2,d2xyw3,d2xyw4
        real*4 d2xxw1,d2xxw2,d2xxw3,d2xxw4

        real*4 ux,ux_2,ux_3,uy,uy_2,uy_3,uz,uz_2,uz_3
	integer ix,iy,iz,ix1,ix2,ix3,ix4,ly,lz

C check bounds, if error, do nothing, set nonzero `ier' flag
	if((x-xm(1))*(x-xm(2)).gt.0.0
     &.or.(y-ym(1))*(y-ym(2)).gt.0.0
     &.or.(z-zm(1))*(z-zm(2)).gt.0.0) then
	    ier=99
	    return
	endif

C x->ix,ux; y->iy,uy; z->iz,uz = local positioning
C e.g., x=xm(1)+ix*dx+ux, ux in [0,1[  [except for x=xm(2)=xm(1)+ix*dx+1.0
C i.e, ux=1.0 in order to avoid reference to undefined elements of vw(,,)].
C -> used spline-weight x coordinates will be vm(ix+1,,)->vm(ix+4,,)
C -> wheras ux will be the argument used in base-polynomials P1->P4
	ux=(x-xm(1))/dx
	ix=ux
	ux=ux-ix
	if(x.eq.xm(2)) then
	    ux=1.0
	    ix=ix-1
	endif
	uy=(y-ym(1))/dy
	iy=uy
	uy=uy-iy
	if(y.eq.ym(2)) then
	    uy=1.0
	    iy=iy-1
	endif
	uz=(z-zm(1))/dz
	iz=uz
	uz=uz-iz
	if(z.eq.zm(2)) then
	    uz=1.0
	    iz=iz-1
	endif

	ux_2=0.5*ux
	uy_2=0.5*uy
	uz_2=0.5*uz
	ux_3=ux/3.0
	uy_3=uy/3.0
	uz_3=uz/3.0

C First Step : interpolating in the x direction
C One stage = 16 1D full interpolations (i.e., at orders 0,1,2)

	ix1=ix+1
	ix2=ix+2
	ix3=ix+3
	ix4=ix+4
C y=1,2,3,4 ; z=1
        lz=iz+1
        ly=iy+1
	call lspl1d2(vm(ix1,ly,lz),vm(ix2,ly,lz),vm(ix3,ly,lz),vm(ix4,ly,lz),
     +              ux,ux_2,ux_3,w11,dxw11,d2xxw11)
        ly=ly+1
	call lspl1d2(vm(ix1,ly,lz),vm(ix2,ly,lz),vm(ix3,ly,lz),vm(ix4,ly,lz),
     +              ux,ux_2,ux_3,w21,dxw21,d2xxw21)
        ly=ly+1
	call lspl1d2(vm(ix1,ly,lz),vm(ix2,ly,lz),vm(ix3,ly,lz),vm(ix4,ly,lz),
     +              ux,ux_2,ux_3,w31,dxw31,d2xxw31)
        ly=ly+1
	call lspl1d2(vm(ix1,ly,lz),vm(ix2,ly,lz),vm(ix3,ly,lz),vm(ix4,ly,lz),
     +              ux,ux_2,ux_3,w41,dxw41,d2xxw41)
C y=1,2,3,4 ; z=2
        lz=lz+1
        ly=iy+1
	call lspl1d2(vm(ix1,ly,lz),vm(ix2,ly,lz),vm(ix3,ly,lz),vm(ix4,ly,lz),
     +              ux,ux_2,ux_3,w12,dxw12,d2xxw12)
        ly=ly+1
	call lspl1d2(vm(ix1,ly,lz),vm(ix2,ly,lz),vm(ix3,ly,lz),vm(ix4,ly,lz),
     +              ux,ux_2,ux_3,w22,dxw22,d2xxw22)
        ly=ly+1
	call lspl1d2(vm(ix1,ly,lz),vm(ix2,ly,lz),vm(ix3,ly,lz),vm(ix4,ly,lz),
     +              ux,ux_2,ux_3,w32,dxw32,d2xxw32)
        ly=ly+1
	call lspl1d2(vm(ix1,ly,lz),vm(ix2,ly,lz),vm(ix3,ly,lz),vm(ix4,ly,lz),
     +              ux,ux_2,ux_3,w42,dxw42,d2xxw42)
C y=1,2,3,4 ; z=3
        lz=lz+1
        ly=iy+1
	call lspl1d2(vm(ix1,ly,lz),vm(ix2,ly,lz),vm(ix3,ly,lz),vm(ix4,ly,lz),
     +              ux,ux_2,ux_3,w13,dxw13,d2xxw13)
        ly=ly+1
	call lspl1d2(vm(ix1,ly,lz),vm(ix2,ly,lz),vm(ix3,ly,lz),vm(ix4,ly,lz),
     +              ux,ux_2,ux_3,w23,dxw23,d2xxw23)
        ly=ly+1
	call lspl1d2(vm(ix1,ly,lz),vm(ix2,ly,lz),vm(ix3,ly,lz),vm(ix4,ly,lz),
     +              ux,ux_2,ux_3,w33,dxw33,d2xxw33)
        ly=ly+1
	call lspl1d2(vm(ix1,ly,lz),vm(ix2,ly,lz),vm(ix3,ly,lz),vm(ix4,ly,lz),
     +              ux,ux_2,ux_3,w43,dxw43,d2xxw43)
C y=1,2,3,4 ; z=4
        lz=lz+1
        ly=iy+1
	call lspl1d2(vm(ix1,ly,lz),vm(ix2,ly,lz),vm(ix3,ly,lz),vm(ix4,ly,lz),
     +              ux,ux_2,ux_3,w14,dxw14,d2xxw14)
        ly=ly+1
	call lspl1d2(vm(ix1,ly,lz),vm(ix2,ly,lz),vm(ix3,ly,lz),vm(ix4,ly,lz),
     +              ux,ux_2,ux_3,w24,dxw24,d2xxw24)
        ly=ly+1
	call lspl1d2(vm(ix1,ly,lz),vm(ix2,ly,lz),vm(ix3,ly,lz),vm(ix4,ly,lz),
     +              ux,ux_2,ux_3,w34,dxw34,d2xxw34)
        ly=ly+1
	call lspl1d2(vm(ix1,ly,lz),vm(ix2,ly,lz),vm(ix3,ly,lz),vm(ix4,ly,lz),
     +              ux,ux_2,ux_3,w44,dxw44,d2xxw44)

C Second step : interpolating in the y direction, 3 stages
C First stage : 4 full interpolations
C z=1,2,3,4
	call lspl1d2(w11,w21,w31,w41,uy,uy_2,uy_3,w1,dyw1,d2yyw1)
	call lspl1d2(w12,w22,w32,w42,uy,uy_2,uy_3,w2,dyw2,d2yyw2)
	call lspl1d2(w13,w23,w33,w43,uy,uy_2,uy_3,w3,dyw3,d2yyw3)
	call lspl1d2(w14,w24,w34,w44,uy,uy_2,uy_3,w4,dyw4,d2yyw4)
C Second stage : 4 interpolations at orders (0,1)
C z=1,2,3,4
	call lspl1d1(dxw11,dxw21,dxw31,dxw41,uy,uy_2,uy_3,dxw1,d2xyw1)
	call lspl1d1(dxw12,dxw22,dxw32,dxw42,uy,uy_2,uy_3,dxw2,d2xyw2)
	call lspl1d1(dxw13,dxw23,dxw33,dxw43,uy,uy_2,uy_3,dxw3,d2xyw3)
	call lspl1d1(dxw14,dxw24,dxw34,dxw44,uy,uy_2,uy_3,dxw4,d2xyw4)
C Third stage : 4 interpolations at order 0
C z=1,2,3,4
	call lspl1d0(d2xxw11,d2xxw21,d2xxw31,d2xxw41,uy,uy_2,uy_3,d2xxw1)
	call lspl1d0(d2xxw12,d2xxw22,d2xxw32,d2xxw42,uy,uy_2,uy_3,d2xxw2)
	call lspl1d0(d2xxw13,d2xxw23,d2xxw33,d2xxw43,uy,uy_2,uy_3,d2xxw3)
	call lspl1d0(d2xxw14,d2xxw24,d2xxw34,d2xxw44,uy,uy_2,uy_3,d2xxw4)

C Third and last step : interpolating in the z direction, 3 stages
C First stage : 1 full interpolation -> w, dzw, d2zzw
	call lspl1d2(w1,w2,w3,w4,uz,uz_2,uz_3,w,dzw,d2zzw)
C Second stage : 2 interpolations at orders (0,1) -> dyw, d2yzw, dxw, d2zxw
	call lspl1d1(dyw1,dyw2,dyw3,dyw4,uz,uz_2,uz_3,dyw,d2yzw)
	call lspl1d1(dxw1,dxw2,dxw3,dxw4,uz,uz_2,uz_3,dxw,d2zxw)
C Third stage : 3 interpolations at order 0 -> d2yyw, d2xyw, d2xxw
	call lspl1d0(d2yyw1,d2yyw2,d2yyw3,d2yyw4,uz,uz_2,uz_3,d2yyw)
	call lspl1d0(d2xyw1,d2xyw2,d2xyw3,d2xyw4,uz,uz_2,uz_3,d2xyw)
	call lspl1d0(d2xxw1,d2xxw2,d2xxw3,d2xxw4,uz,uz_2,uz_3,d2xxw)

C normalize
	dxw=dxw/dx
	dyw=dyw/dy
	dzw=dzw/dz
	d2xxw=d2xxw/(dx*dx)
	d2yyw=d2yyw/(dy*dy)
	d2zzw=d2zzw/(dz*dz)
	d2xyw=d2xyw/(dx*dy)
	d2yzw=d2yzw/(dy*dz)
	d2zxw=d2zxw/(dz*dx)

	ier=0
	return
	end

C-----------------------------------------------------------
C the following routines are 1D basic Bspline interpolators  
C (purely local) computing orders 0,1,2 or 0,1 or 0 only.
C-----------------------------------------------------------

	subroutine lspl1d2(w1,w2,w3,w4,u,u_2,u_3,S,dS,d2S)
C Full 1D interpolator (f,gradf,gradgradf)
C Cost : 5 multadd, 9 add, 2 mult
        implicit none
	real*4 w1,w2,w3,w4,u,u_2,u_3,S,dS,d2S
	real*4 d32,f0,f1,f2,f3,tmp0,tmp1
	f1=0.5*(w3-w1)
	d32=w3-w2
	f2=w2-w1-d32
	f3=f2-d32-w3+w4
	f0=f2/6.0-w2
	S=((f3*u_3-f2)*u_2+f1)*u-f0
	tmp0=f3*u_2
	tmp1=tmp0-f2
	dS=tmp1*u+f1
	d2S=tmp0+tmp1
	return
	end

	subroutine lspl1d1(w1,w2,w3,w4,u,u_2,u_3,S,dS)
C 1D interpolator (f,gradf)
C Cost : 5 multadd, 7 add, 1 mult
        implicit none
	real*4 w1,w2,w3,w4,u,u_2,u_3,S,dS
	real*4 d32,f0,f1,f2,f3
	f1=0.5*(w3-w1)
	d32=w3-w2
	f2=w2-w1-d32
	f3=f2-d32-w3+w4
	f0=f2/6.0-w2
	S=((f3*u_3-f2)*u_2+f1)*u-f0
	dS=(f3*u_2-f2)*u+f1
	return
	end

	subroutine lspl1d0(w1,w2,w3,w4,u,u_2,u_3,S)
C 1D interpolator (f)
C Cost : 3 multadd, 7 add, 1 mult
        implicit none
	real*4 w1,w2,w3,w4,u,u_2,u_3,S
	real*4 d32,f0,f1,f2,f3
	f1=0.5*(w3-w1)
	d32=w3-w2
	f2=w2-w1-d32
	f3=f2-d32-w3+w4
	f0=f2/6.0-w2
	S=((f3*u_3-f2)*u_2+f1)*u-f0
	return
	end
