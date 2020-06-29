      subroutine extirp(b1,yt)
      implicit none
      real a1(3,3),a2(3,3),b1(3,2),b2(3,2)
      real a1inv(3,3),b3(3,2)
      real yt(36)
      integer i,j,k

c     ----- dX/dphi1 -----------
      a1(1,1)=yt(8)
      a1(2,1)=yt(9)
      a1(3,1)=yt(10)
c     ----- dX/dphi2 -----------
      a1(1,2)=yt(11)
      a1(2,2)=yt(12)
      a1(3,2)=yt(13)
c     ----- dX/dtau -----------
      a1(1,3)=yt(4)
      a1(2,3)=yt(5)
      a1(3,3)=yt(6)

c     ----- dP/dphi1 -----------
      a2(1,1)=yt(14)
      a2(2,1)=yt(15)
      a2(3,1)=yt(16)
c     ----- dP/dphi2 -----------
      a2(1,2)=yt(17)
      a2(2,2)=yt(18)
      a2(3,2)=yt(19)
c     ----- dP/dtau -----------
      a2(1,3)=yt(32)
      a2(2,3)=yt(33)
      a2(3,3)=yt(34)

c     ----- dX/dxs -----------  
      b1(1,1)=yt(20)  
      b1(2,1)=yt(21)  
      b1(3,1)=yt(22)  
c     ----- dX/dpys -----------  
      b1(1,2)=yt(23)  
      b1(2,2)=yt(24)  
      b1(3,2)=yt(25)  
 
c     ----- dP/dxs ----------- 
      b2(1,1)=yt(26) 
      b2(2,1)=yt(27) 
      b2(3,1)=yt(28) 
c     ----- dP/dpys ----------- 
      b2(1,2)=yt(29) 
      b2(2,2)=yt(30) 
      b2(3,2)=yt(31) 

      call inv3x3(a1,a1inv)

      do 1 i=1,3
         do 1 j=1,2
            b3(i,j)=0.
            do 1 k=1,3
               b3(i,j)=b3(i,j)+a1inv(i,k)*b1(k,j)
1     continue

      do 2 i=1,3 
         do 2 j=1,2 
            b1(i,j)=0. 
            do 2 k=1,3 
               b1(i,j)=b1(i,j)+a2(i,k)*b3(k,j) 
2     continue

      do 3 i=1,3  
         do 3 j=1,2  
            b1(i,j)=-b1(i,j)+b2(i,j)
3     continue

      return
      end
