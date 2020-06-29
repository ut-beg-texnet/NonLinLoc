      subroutine green3(nxm,nym,nzm,dxm,dym,dzm,xm,ym,zm,vm,
     &          nxr,nyr,nzr,nir,npr,targ_orient,xmap,imap,
     &          xsrc,fi1min,fi2min,fi1max,fi2max,dxmin2,dpmin2,dtemps,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)
c-------------------------------------------------------------------------
c-     Programed by p.S.Lucio and G.C.Lambare in Apr/94                  -
c-     Last modification by P.S.Lucio in Jun/20/95                       -
c-     GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP)    -
c-------------------------------------------------------------------------
       implicit none

c       #########################################################
c       This program was developed by P.S. LUCIO and G.C.Lambare
c       in the Centre de Geophysique at Ecole des Mines de Paris
c               and A. Hanyga at University of Bergen
c       This is a tool for the calculus of the 3D Greens function
c       #########################################################

C   *****************************************************************
C   *      LES VARIABLES Y ET DY SONT DES PARAMETRES DE RAI OU:     *
C   *                                                               *
C   *    Y(1),  Y(2),  Y(3) : X,    Y,    Z        }                *
C   *    Y(4),  Y(5),  Y(6) : PX,   PY,   PZ       } propagation    *
C   *    Y(7),  Y(8),  Y(9) : DX1,  DY1,  DPZ1                      *
C   *    Y(10), Y(11), Y(12): DX2,  DY2,  DPZ2                      *
C   *    Y(13), Y(14), Y(15): DPX1, DPY1, DPZ1                      *
C   *    Y(16), Y(17), Y(18): DPX2, DPY2, DPZ2                      *
C   *    Y(19)              : TEMPS                                 *
C   *    Y(20)              : TAU                                   *
C   *    Y(21), Y(22), Y(23) : X,    Y,    Z       }                *
C   *    Y(24), Y(25), Y(26) : PX,   PY,   PZ      } interpolation  *
C   *    Y(27), Y(28), Y(29) : DVX/2, DVY/2, DVZ/2                  * 
C   *    Y(30)               : V                                    *
C   *    Y(31), Y(32), Y(33) :                     } the cross      *
C   *    Y(34), Y(35), Y(36) :                     } derivatives    *
C   *    Y(37), Y(38), Y(39) : DXxs,  DYxs,  DPZxs                  *
C   *    Y(40), Y(41), Y(42) : DXys,  DYys,  DPZys                  *
C   *    Y(43), Y(44), Y(45) : DPXxs, DPYxs, DPZxs                  *
C   *    Y(46), Y(47), Y(48) : DPXys, DPYys, DPZys                  *
C   *****************************************************************
C      **************************************************************
C      *   relation number of the final ray and initial ray         *
C      **************************************************************
C      *    iYs(1) : the KMAH index                                 *
C      *    iYs(2) : indicator of the ray propagation               *
C      *                   -1 echec at the ray [irai]               *
C      *                    number of the propaged ray              *
C      *                    0 it will be propage                    *
C      **************************************************************
C      **************************************************************
C*                     ABOUT THE SEGMENTS                           
C********************************************************************
C*     iseg(nsn,1): number of the first point                       *
C*     iseg(nsn,2): number of the second point                      *
C*     iseg(nsn,3): indicator about the subdivision of the segment  *
C*                  definies by two rays                            *
C*                 -1           : success on test                   *
C*                  0           : it's not have been tested         *
C*                  n           : number of the midpoint            *
C*     iseg(nsn,4): number of the first segment included            *
C*     iseg(nsn,5): number of the second segment included           *
C*     iseg(nsn,6): number of the segment on the propaged fron      *
C********************************************************************
C      **************************************************************
C*                      ABOUT THE CELLS                            
C********************************************************************
C      icel(ncn,1): number of the first point                       *
C      icel(ncn,2): number of the second point                      *
C      icel(ncn,3): number of the third point                       *
C             NUMBER OF THE THREE SEGMENTS OF THE CELL              *
C      icel(ncn,4): number of the first segment of the cell         *
C      icel(ncn,5): number of the second segment of thwe cell       *
C      icel(ncn,6): number of the third segment of the cell         *
C********************************************************************
C********************************************************************
C*    nrai: number of "rays" at the front                           *
C*    nseg: number of segments at the front                         *
C*    ncel: number of cells at the front                            *
C********************************************************************
C********************************************************************
C*    irai: number of the point of the cell                         *
C     irpop: number of opposed points to each one segment           *
C*    irais: number of the two "rays" of each one of the three      *
C*           segments                                               *
C********************************************************************
C********************************************************************
C*    isatr(ip) : index of wavefront saturation (rays)              *
C*    isats(ip) : index of wavefront saturation (segments)          *
C*    isatc(ip) : index of wavefront saturation (cells)             *
C********************************************************************
      integer npn,nsn,ncn,nfront,nd,nk,npar   
      real ys(npar,npn,nfront)
      integer iys(nk,npn,nfront)
      double precision dys(nd,npn,nfront)
      integer  iseg(nsn,6,nfront),icel(ncn,6,nfront)
      integer nrai(nfront), ncel(nfront), nseg(nfront)
      integer irai(3), iraf(3), irpop(3), irais(2,3)
      integer isatr(nfront), isats(nfront), isatc(nfront)
c-- 
c     --- Source point ------------------------------------------
      real xsrc(3),usrc

c     --- Velocity field ----------------------------------------
      integer nxm,nym,nzm
      real dxm,dym,dzm
      real vm(nzm,nxm,nym)
      real xm(2),ym(2),zm(2)

c     --- Target ------------------------------------------------
      integer nxr,nyr,nzr,nir,npr
      real xmap(npr,nzr,nxr,nyr,nir)
      integer imap(nzr,nxr,nyr)
      real targ_orient(3,7)

c     --- Ray tracing -------------------------------------------
      real fi1min,fi2min,fi1max,fi2max,dxmin2,dpmin2,dtemps
      real xref,yref,zref
      real xmin,xmax,ymin,ymax,zmin,zmax,zrmin
      integer ip1,ip2,ips,indix,ncrai
      integer isubcel,iechec,ini

      xref=xm(1)
      yref=ym(1)
      zref=zm(1)
      call repos(xm,ym,zm,targ_orient,xsrc,xref,yref,zref)

C******************************************************************
C*    nrai: number of "rays" at the front                         *
C*    nseg: number of segments at the front                       *
C*    ncel: number of cells at the front                          *
C******************************************************************

      if (npr.gt.26) then
         write(*,*)'with this version npr has to be <26'
         stop
      endif

      call limit(targ_orient,nxr,nyr,nzr,xsrc,
     &      xmin,xmax,ymin,ymax,zmin,zmax,zrmin)
C******************************************************************
C*             INITIALIZATION OF WAVEFRONT SAMPLING
C******************************************************************
      call cleanfront
     & (npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)
C******************************************************************
C*             INITIALIZATION OF TARGET 
C******************************************************************
      call cleanmap(nxr,nyr,nzr,nir,npr,xmap,imap)

C******************************************************************
c     --- Index of the front (1) or (2) ---

      ip1=1
      ip2=2

C**********************************************************************
C*              INITIALISATION OF THE RAY AT THE SOURCE               *
C**********************************************************************
      call inisource (ip2,xsrc,usrc,
     &      vm,xm,ym,zm,nxm,nym,nzm,dxm,dym,dzm,
     &                fi1min,fi2min,fi1max,fi2max,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)
         fi1min=1000
         fi1max=-1000
         fi2min=1000
         fi2max=-1000

C        ********************************************************
C        *          INITIALISATION OF A BEAM OF RAYS            *
C        *                   AT THE WAVEFRONT                   *
C        ********************************************************
 
c     --- Number of the front ---

      indix=1
200   continue 

C        ********************************************************
C        *          COMMUTATION OF THE TWO WAVEFRONTS           *
C        ********************************************************
 
      ips=ip1
      ip1=ip2
      ip2=ips

      if (mod(indix,1).eq.0)
     &   write(*,*)
     &                'Front',indix,' : ',ncel(ip1),' cells',
     &                                 nseg(ip1),' segments',
     &                                 nrai(ip1),' rays'

      if (ncel(ip1).eq.0) goto 400

      indix=indix+1

C        ********************************************************
C        *    CALCULUS OF THE REFERENCE RAY AND THE PARAXIALS   *
C        *                   (THE TRACING)                      *
C        *          INITIALISATION OF A BEAM OF RAYS            *
C        ********************************************************
 
C*              NONE RAY ARE STILL PROPAGETED AT THE NEW FRONT
      nrai(ip2)=0
      isatr(ip2)=0
      isatr(ip1)=0
C*              NONE SEGMENT ARE STILL DEFINED AT THE NEW FRONT
      nseg(ip2)=0
      isats(ip2)=0
      isats(ip1)=0
C*              NONE CELL ARE STILL DEFINED AT THE NEW FRONT
 
      ncel(ip2)=0
      isatc(ip2)=0
      isatc(ip1)=0
C*              WE LOOPS AT THE CELLS ON THE PREVIOUS FRONT
      ncrai=0

500   call infronti (ip1,ncrai,ini,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)

      if (ini.eq.1) goto 200
      if (ini.eq.2) then
        write(*,*) 'WAVEFRONT SATURED: SEE THE DIMENSION OF TABLES'
        goto 200
      endif

C*                        PROPAGATION OF A CELL

      iechec=0

      call cellule 
     &  (ip1,ip2,indix,iechec,
     &   vm,xm,ym,zm,nxm,nym,nzm,dxm,dym,dzm,dtemps,
     &      xmin,xmax,ymin,ymax,zmin,zmax,zrmin,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)
      if (iechec.ne.0) goto 500

c     ********************************************************
c     *               TEST THE SIZE OF A CELL                *
c     ********************************************************
  
      isubcel=0
      call subcel 
     & (ncrai,ip1,ip2,isubcel,indix,
     & xsrc,vm,xm,ym,zm,nxm,nym,nzm,dxm,dym,dzm,
     & dxmin2,dpmin2,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)

c     ********************************************************
c     *                   CELL TOO LARGE !!!                 *
c     ********************************************************
 
      if (isubcel.ne.0) goto 500

c     ********************************************************
c     *          THE COURENT WAVEFRONT IS SATURED!!          *
c     ********************************************************
 
      if (isatr(ip2).ne.0.or.
     &    isats(ip2).ne.0.or.
     &    isatc(ip2).ne.0) then
        write(*,*) 'END: WAVEFRONT',indix,'IS SATURED'
        goto 200
      endif

c     ********************************************************
c     *                     STOCK THE RAYS                  *
c     ********************************************************
 
      call stock (ip1,ip2,ncrai,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)

c     ********************************************************
c     *                     INTERPOLATION                    *
c     ********************************************************
 
      call intcel 
     & (ip1,ip2,nxr,nyr,nzr,nir,npr,usrc,targ_orient,xmap,imap,
     &           fi1min,fi1max,fi2min,fi2max,
     & npn,nsn,ncn,nfront,nd,nk,npar,
     & ys,iys,dys,iseg,icel,nrai,ncel,nseg,
     & irai,iraf,irpop,irais,isatr,isats,isatc)
      goto 500

400   continue
      call depos(xm,ym,zm,targ_orient,xsrc,xref,yref,zref)

      return
      end
