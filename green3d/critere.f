      subroutine critere 
     & (ip1,ip2,yt1,yt2,yt3,iyt1,iyt2,iyt3,targ_orient,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)
c-------------------------------------------------------------------------
c-     Programed by P.S.Lucio and G.C.Lambare in jan/95                  -
c-     Last modification by P.S.Lucio in Jun/21/95                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c--------------------------------------------------------------------------

c*************************************************************************
c* CRITERION FOR THE CONSTRUCTION OF TETRAHEDRA FROM A TRIANGULAR PRISME *
c*************************************************************************
      implicit none
      integer ip1,ip2
      integer isumrai,isumraf
      integer i,j
      integer irai1,irai2,irai3
      integer iraf1,iraf2,iraf3
      integer iraixx,irafxx

      integer npn,nsn,ncn,nfront,nd,nk,npar 
      real ys(npar,npn,nfront) 
      integer iys(nk,npn,nfront) 
      double precision dys(nd,npn,nfront) 
      integer  iseg(nsn,6,nfront),icel(ncn,6,nfront) 
      integer nrai(nfront), ncel(nfront), nseg(nfront) 
      integer irai(3), iraf(3), irpop(3), irais(2,3) 
      integer isatr(nfront), isats(nfront), isatc(nfront) 
c--
      real yt1(36,4),yt2(36,4),yt3(36,4)
      integer iyt1(1,4),iyt2(1,4),iyt3(1,4)
c--
      real targ_orient(3,7)
c--
      double precision minphi1, maxphi1

c--
      isumrai=irai(1)+irai(2)+irai(3)
      isumraf=iraf(1)+iraf(2)+iraf(3)

c--
      minphi1=dmin1(dys(1,irai(1),ip1),dys(1,irai(2),ip1),
     &              dys(1,irai(3),ip1))
      maxphi1=dmax1(dys(1,irai(1),ip1),dys(1,irai(2),ip1),
     &              dys(1,irai(3),ip1))

c*************************************************************************
c    PROCEDURE FOR NEVER FALLS AT THE STRUCTURES WITCH THE CUTTING OUT   *
c                 OF THE PRISME WILL BE NOT A TETRAHEDRON                *
c                     (TO AVOID THE BORDS PROBLEMES).                    *
c*************************************************************************
 
      do i=1,3
        if (minphi1.eq.dys(1,irai(i),ip1)) then
          irai1=irai(i)
          iraf1=iraf(i)
          goto 10
        endif
      enddo
10    continue
      do j=1,3
        if (maxphi1.eq.dys(1,irai(j),ip1)) then 
          irai3=irai(j)
          iraf3=iraf(j)
          irai2=isumrai-(irai1+irai3)
          iraf2=isumraf-(iraf1+iraf3)
          goto 20
        endif
      enddo
20    continue
      if (dys(1,irai1,ip1).eq.dys(1,irai2,ip1)) then
         if (dys(2,irai1,ip1).gt.dys(2,irai2,ip1)) then
            iraixx=irai1
            irai1=irai2
            irai2=iraixx
            irafxx=iraf1
            iraf1=iraf2
            iraf2=irafxx
          endif
      endif
      if (dys(1,irai2,ip1).eq.dys(1,irai3,ip1)) then
         if (dys(2,irai2,ip1).gt.dys(2,irai3,ip1)) then
            iraixx=irai2
            irai2=irai3
            irai3=iraixx
            irafxx=iraf2
            iraf2=iraf3
            iraf3=irafxx
          endif
      endif

c*************************************************************************
 
c     ----    ---> points of the first tetrahedrum  -----
 
      call valtetra(ip1,irai1,yt1(1,1),iyt1(1,1),targ_orient,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg, 
     & irai,iraf,irpop,irais,isatr,isats,isatc)
      call valtetra(ip1,irai2,yt1(1,2),iyt1(1,2),targ_orient,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg, 
     & irai,iraf,irpop,irais,isatr,isats,isatc)
      call valtetra(ip1,irai3,yt1(1,3),iyt1(1,3),targ_orient,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg, 
     & irai,iraf,irpop,irais,isatr,isats,isatc)
      call valtetra(ip2,iraf3,yt1(1,4),iyt1(1,4),targ_orient,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg, 
     & irai,iraf,irpop,irais,isatr,isats,isatc)

c     ----    ---> points of the second tetrahedrum   -----
 
      call valtetra(ip1,irai1,yt2(1,1),iyt2(1,1),targ_orient,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg, 
     & irai,iraf,irpop,irais,isatr,isats,isatc)
      call valtetra(ip1,irai2,yt2(1,2),iyt2(1,2),targ_orient,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg, 
     & irai,iraf,irpop,irais,isatr,isats,isatc)
      call valtetra(ip2,iraf2,yt2(1,3),iyt2(1,3),targ_orient,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg, 
     & irai,iraf,irpop,irais,isatr,isats,isatc)
      call valtetra(ip2,iraf3,yt2(1,4),iyt2(1,4),targ_orient,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg, 
     & irai,iraf,irpop,irais,isatr,isats,isatc)

c     ----    ---> points of the third tetrahedrum    -----

      call valtetra(ip1,irai1,yt3(1,1),iyt3(1,1),targ_orient,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg, 
     & irai,iraf,irpop,irais,isatr,isats,isatc)
      call valtetra(ip2,iraf1,yt3(1,2),iyt3(1,2),targ_orient,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg, 
     & irai,iraf,irpop,irais,isatr,isats,isatc)
      call valtetra(ip2,iraf2,yt3(1,3),iyt3(1,3),targ_orient,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg, 
     & irai,iraf,irpop,irais,isatr,isats,isatc)
      call valtetra(ip2,iraf3,yt3(1,4),iyt3(1,4),targ_orient,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg, 
     & irai,iraf,irpop,irais,isatr,isats,isatc)
c--

      return 
      end
