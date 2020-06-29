        subroutine repos(xm,ym,zm,targ_orient,xsrc,xref,yref,zref)
        implicit none
        real xref,yref,zref
        real xm(2),ym(2),zm(2)
        real targ_orient(3,7)
        real xsrc(3)
        integer i

        do i=1,2
           xm(i)=xm(i)-xref
           ym(i)=ym(i)-yref
           zm(i)=zm(i)-zref
        enddo

        targ_orient(1,1)=targ_orient(1,1)-xref 
        targ_orient(2,1)=targ_orient(2,1)-yref 
        targ_orient(3,1)=targ_orient(3,1)-zref 

        xsrc(1)=xsrc(1)-xref
        xsrc(2)=xsrc(2)-yref
        xsrc(3)=xsrc(3)-zref

        return
        end

        subroutine depos(xm,ym,zm,targ_orient,xsrc,xref,yref,zref) 
        implicit none
        real xref,yref,zref
        real xm(2),ym(2),zm(2)
        real targ_orient(3,7)
        real xsrc(3)
        integer i

        do i=1,2
           xm(i)=xm(i)+xref
           ym(i)=ym(i)+yref
           zm(i)=zm(i)+zref
        enddo
 
        targ_orient(1,1)=targ_orient(1,1)+xref
        targ_orient(2,1)=targ_orient(2,1)+yref
        targ_orient(3,1)=targ_orient(3,1)+zref
 
        xsrc(1)=xsrc(1)+xref
        xsrc(2)=xsrc(2)+yref
        xsrc(3)=xsrc(3)+zref
 
        return
        end
