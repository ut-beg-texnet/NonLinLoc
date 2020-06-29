c****************************************************************************
c****************************************************************************
      subroutine modelrese(targ_orient)
c   ***************************************************************
c   *             INPUT OF THE 3D RESERVOIR PARAMETERS            *
c   ***************************************************************
      implicit none
      real targ_orient(3,7)
 
      print*,'origine point of target x,y,z ???'
      read(*,*) targ_orient(1,1),targ_orient(2,1),targ_orient(3,1)
      write(*,*) targ_orient(1,1),targ_orient(2,1),targ_orient(3,1)
      print*,'vector for first increment in target ???'
      read(*,*) targ_orient(1,2),targ_orient(2,2),targ_orient(3,2)
      write(*,*) targ_orient(1,2),targ_orient(2,2),targ_orient(3,2)
      print*,'vector for second increment in target ???'
      read(*,*) targ_orient(1,3),targ_orient(2,3),targ_orient(3,3)
      write(*,*) targ_orient(1,3),targ_orient(2,3),targ_orient(3,3)
      print*,'vector for third increment in target ???'
      read(*,*) targ_orient(1,4),targ_orient(2,4),targ_orient(3,4)
      write(*,*) targ_orient(1,4),targ_orient(2,4),targ_orient(3,4)

      call inv3x3(targ_orient(1,2),targ_orient(1,5))

      return
      end
