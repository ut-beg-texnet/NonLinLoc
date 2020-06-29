         subroutine limit(targ_orient,nxr,nyr,nzr,xsrc,
     &      xmin,xmax,ymin,ymax,zmin,zmax,zrmin)
         implicit none
         real xmarg,ymarg,zmarg
         parameter (xmarg=500,ymarg=500,zmarg=500)
         real targ_orient(3,7),xsrc(3)
         real x1(3),x2(3),x3(3),x4(3),x5(3),x6(3),x7(3),x8(3)
         integer nxr,nyr,nzr
         real xmin,xmax,ymin,ymax,zmin,zmax,zrmin

         call position(targ_orient,  1,  1,  1,x1)
         call position(targ_orient,  1,  1,nzr,x2)
         call position(targ_orient,  1,nyr,nzr,x3)
         call position(targ_orient,  1,nyr,  1,x4)
         call position(targ_orient,nxr,  1,  1,x5)
         call position(targ_orient,nxr,  1,nzr,x6)
         call position(targ_orient,nxr,nyr,nzr,x7)
         call position(targ_orient,nxr,nyr,  1,x8)

C_AJL         xmin=amin0(xsrc(1),x1(1),x2(1),x3(1),x4(1),
C_AJL     &                      x5(1),x6(1),x7(1),x8(1))-xmarg
         xmin=amin1(xsrc(1),x1(1),x2(1),x3(1),x4(1),
     &                      x5(1),x6(1),x7(1),x8(1))-xmarg
C_AJL         xmax=amax0(xsrc(1),x1(1),x2(1),x3(1),x4(1),
C_AJL     &                      x5(1),x6(1),x7(1),x8(1))+xmarg
         xmax=amax1(xsrc(1),x1(1),x2(1),x3(1),x4(1),
     &                      x5(1),x6(1),x7(1),x8(1))+xmarg

C_AJL         ymin=amin0(xsrc(2),x1(2),x2(2),x3(2),x4(2),
C_AJL     &                      x5(2),x6(2),x7(2),x8(2))-ymarg
         ymin=amin1(xsrc(2),x1(2),x2(2),x3(2),x4(2),
     &                      x5(2),x6(2),x7(2),x8(2))-ymarg
C_AJL         ymax=amax0(xsrc(2),x1(2),x2(2),x3(2),x4(2),
C_AJL     &                      x5(2),x6(2),x7(2),x8(2))+ymarg
         ymax=amax1(xsrc(2),x1(2),x2(2),x3(2),x4(2),
     &                      x5(2),x6(2),x7(2),x8(2))+ymarg

C_AJL         zmin=amin0(xsrc(3),x1(3),x2(3),x3(3),x4(3),
C_AJL     &                      x5(3),x6(3),x7(3),x8(3))-zmarg
         zmin=amin1(xsrc(3),x1(3),x2(3),x3(3),x4(3),
     &                      x5(3),x6(3),x7(3),x8(3))-zmarg
C_AJL         zmax=amax0(xsrc(3),x1(3),x2(3),x3(3),x4(3),
C_AJL     &                      x5(3),x6(3),x7(3),x8(3))+zmarg
         zmax=amax1(xsrc(3),x1(3),x2(3),x3(3),x4(3),
     &                      x5(3),x6(3),x7(3),x8(3))+zmarg
         write(*,*)'effective velocity field (cf limit.f)'
         write(*,*)'XMIN,XMAX',xmin,xmax
         write(*,*)'YMIN,YMAX',ymin,ymax
         write(*,*)'ZMIN,ZMAX',zmin,zmax

C_AJL         zrmin=amin0(x1(3),x2(3),x3(3),x4(3),
C_AJL     &              x5(3),x6(3),x7(3),x8(3))
         zrmin=amin1(x1(3),x2(3),x3(3),x4(3),
     &              x5(3),x6(3),x7(3),x8(3))

         return
         end


         subroutine position(targ_orient,i,j,k,x)
         implicit none
         integer i,j,k,l
         real targ_orient(3,7),x(3)

         do l=1,3
         x(l)=targ_orient(l,1)+(i-1)*targ_orient(l,2)
     &                        +(j-1)*targ_orient(l,3)
     &                        +(k-1)*targ_orient(l,4)
         enddo

         return
         end

         
