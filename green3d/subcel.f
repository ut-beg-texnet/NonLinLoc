       subroutine subcel (ncrai,ip1,ip2,isubcel,indix,xsrc,
     &                    vm,xm,ym,zm,nxm,nym,nzm,dxm,dym,dzm,
     &                    dxmin2,dpmin2,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)

c-------------------------------------------------------------------------
c-     Programed by P.S.Lucio and G.C.Lambare in Sep/94                  -
c-     Last modification by P.S.Lucio in Oct/10/95                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c-------------------------------------------------------------------------

c*************************************************************************
c*             THIS SR IS TO VERIFY THE SIZE OF A SEGMENT                *
c*                 (MISFIT IN POSITION AND SLOWNESS)                     *
c*************************************************************************
       implicit none

      integer npn,nsn,ncn,nfront,nd,nk,npar
      real ys(npar,npn,nfront)
      integer iys(nk,npn,nfront)
      double precision dys(nd,npn,nfront)
      integer  iseg(nsn,6,nfront),icel(ncn,6,nfront)
      integer nrai(nfront), ncel(nfront), nseg(nfront)
      integer irai(3), iraf(3), irpop(3), irais(2,3)
      integer isatr(nfront), isats(nfront), isatc(nfront)
c--
       integer ncrai,ip1,ip2,isubcel,indix
       integer nxm,nym,nzm
       real dxm,dym,dzm
       real vm(nzm,nxm,nym)
       real xm(2),ym(2),zm(2)
c--  
       real xsrc(3)
       real dxmin2,dpmin2

       integer i,j,irainew,isegnew1,isegnew2,isegnew3
       integer isegment,indtest,icelnew1,isegpop1,isegpop2
       integer icelnew2,itest,isegnew

c********************************************************************
c INPUT  : icrai -
c          ip1 - index of the previous wavefront
c          ip2 - index of the courent wavefront
c          indix - the wavefront indicator
c
c OUTPUT : isubcel - index of subdivision of the cell
c          isatr(ip) - index of wavefront saturation on rays
c          isats(ip) - index of wavefront saturation on segments
c          isatc(ip) - index of wavefront saturation on cells
c********************************************************************
 
C********************************************************************
C*     iseg(nsn,1): indicator about the subdivision (first point)   *
C*     iseg(nsn,2): indicator about the subdivision (second point)  *
C*     iseg(nsn,3): indicator about the subdivision of the segment  *
C*                  definies by two rays                            *
C*                 -1           : success on test                   *
C*                  0           : it's not have been tested         *
C*                  n           : number of the midpoint            *
C*     iseg(nsn,4): number of the first segment included            *
C*     iseg(nsn,5): number of the second segment included           *
C*     iseg(nsn,6): number of the segment on the propaged front     *
C********************************************************************
 
c--- we include the rays, segments and cells to the wavefront ---
 
       do 10 i=1,3
          isegment=icel(ncrai,i+3,ip1)
          indtest=iseg(isegment,3,ip1)

*********************************************************************
          if (indtest.eq.-1) then
*********************************************************************

c            ----  the segment have been tested with success ----

             goto 10

*********************************************************************
          elseif (indtest.eq.0) then
*********************************************************************

c       -- test about the saturation of the previous wavefront --
 
               if (isatr(ip1).ne.0.or.isats(ip1).ne.0.or.
     &             isatc(ip1).ne.0) goto 100

c       -- the segment have not been still tested (we must test it) --
 
             itest=0

             call testseg(ip2,iys(2,irais(1,i),ip1),
     &                    iys(2,irais(2,i),ip1),indix,itest,
     &                    dxmin2,dpmin2,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)
 
             if (itest.eq.0) then
               call testseg(ip2,iys(2,irais(2,i),ip1),
     &                      iys(2,irais(1,i),ip1),indix,itest,
     &                      dxmin2,dpmin2,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)
 
             endif

             if (itest.eq.-1) then

c  ----  the test is negatif, we add a point and three segments ----
 
               call addrai (ip1,irais(1,i),irais(2,i),indix,xsrc,
     &                      vm,xm,ym,zm,nxm,nym,nzm,dxm,dym,dzm,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)
 

               if (isatr(ip1).ne.0) goto 100
c--
               irainew=nrai(ip1)
               isegnew1=nseg(ip1)+1
               isegnew2=nseg(ip1)+2
               isegnew3=nseg(ip1)+3
               nseg(ip1)=nseg(ip1)+3
               iseg(isegment,3,ip1)=irainew
               iseg(isegment,4,ip1)=isegnew1
               iseg(isegment,5,ip1)=isegnew2

               call addseg (ip1,irais(1,i),irainew,isegnew1,indix,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)


               if (isats(ip1).ne.0) goto 100

               call addseg (ip1,irais(2,i),irainew,isegnew2,indix,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)


               if (isats(ip1).ne.0) goto 100

               call addseg (ip1,irpop(i)  ,irainew,isegnew3,indix,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)


               if (isats(ip1).ne.0) goto 100

c             --- we add two cells ---
 
               do j=1,3
                  if (j.ne.i) then
                    if (irais(1,j).eq.irais(1,i).or.
     &                  irais(2,j).eq.irais(1,i)) then
                      isegpop1=icel(ncrai,j+3,ip1)
                    elseif (irais(1,j).eq.irais(2,i).or.
     &                      irais(2,j).eq.irais(2,i)) then
                      isegpop2=icel(ncrai,j+3,ip1)
                    endif
                  endif
               enddo
               icelnew1=ncrai
               ncrai=ncrai-1
               ncel(ip1)=ncel(ip1)+1
               icelnew2=ncel(ip1)

               call addcel (ip1,irpop(i),irainew,irais(1,i),isegpop1,
     &                      isegnew1,isegnew3,icelnew1,indix,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)


               if (isatc(ip1).ne.0) goto 100

               call addcel (ip1,irpop(i),irainew,irais(2,i),isegpop2,
     &                      isegnew3,isegnew2,icelnew2,indix,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)


               if (isatc(ip1).ne.0) goto 100

               isubcel=-1

               return

c                      ----  the test is positif ----
 
             elseif (itest.eq.0) then 

100          continue

c             ---- we add a segment on the courent front ----

                nseg(ip2)=nseg(ip2)+1

                call addseg (ip2,iys(2,irais(1,i),ip1),
     &                       iys(2,irais(2,i),ip1),nseg(ip2),indix,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)


                if (isats(ip2).ne.0) return

                iseg(isegment,3,ip1)=-1
                iseg(isegment,6,ip1)=nseg(ip2)
                goto 10
             endif

*********************************************************************
           elseif (indtest.gt.0) then 
*********************************************************************

c            ----  the segment have been divided  ----
 
             irainew=indtest

c   -we include a segment and two cells and we cut out the courent segment-
 
             isegnew=isegment
             isegnew1=iseg(isegment,4,ip1)
             isegnew2=iseg(isegment,5,ip1)

             call addseg (ip1,irpop(i),irainew,isegnew,indix,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)


             if (isats(ip1).ne.0) goto 5

             do j=1,3
                if (j.ne.i) then
                  if (irais(1,j).eq.irais(1,i).or.
     &                irais(2,j).eq.irais(1,i)) then
                      isegpop1=icel(ncrai,j+3,ip1)
                  elseif (irais(1,j).eq.irais(2,i).or. 
     &                    irais(2,j).eq.irais(2,i)) then
                          isegpop2=icel(ncrai,j+3,ip1)
                  endif
                endif
             enddo
             icelnew1=ncrai
             ncrai=ncrai-1
             ncel(ip1)=ncel(ip1)+1
             icelnew2=ncel(ip1)

             call addcel (ip1,irpop(i),irainew,irais(1,i),isegpop1,
     &                    isegnew1,isegnew,icelnew1,indix,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)


             if (isatc(ip1).ne.0) goto 5

             call addcel (ip1,irpop(i),irainew,irais(2,i),isegpop2, 
     &                    isegnew,isegnew2,icelnew2,indix,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)


             if (isatc(ip1).ne.0) goto 5

5            continue

             isubcel=-1

             return
           endif
10     continue

       isubcel=0

       return
       end
