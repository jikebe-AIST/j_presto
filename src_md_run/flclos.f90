
      subroutine flclos(iounit,idacc,ier)

!*******************************************************************
!
!     THIS SUBROUTINE IS FOR CLOSE FILE
!                                                                    
!*******************************************************************

      implicit none

      ! Logical unit number
        integer(4),intent(in):: iounit
      ! ID number for access (10: KEEP, 11: DELETE)
        integer(4),intent(in):: idacc
      ! Condition code (0: NO ERROR, 1: PARAMETER ERROR, 2: CLOSE ERROR)
        integer(4),intent(out):: ier

!*********************************

      ier = 0
!     1) CHECK PARAMETER 
      if ( (iounit.lt.0 .or. iounit.ge.100) .or.                       &
     &     (idacc.ne.10 .and. idacc.ne.11) ) then
        ier = -1 ; return
      endif
                          
!     2) CLOSE FILE
      if ( idacc .eq. 10 ) then
        close(unit=iounit,status="KEEP",iostat=ier)
      elseif ( idacc .eq. 11 ) then
        close(unit=iounit,status="DELETE",iostat=ier)
      endif

      if ( ier .ge. 0 ) then
        ier = 0
      else
        ier = -2
      endif

!******************************

      return
      end subroutine flclos
