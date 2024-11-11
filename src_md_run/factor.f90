
      function factor(n)
 
!*******************************************************************
!
!     THIS FUNCTION IS FOR CALCULATION FACTORIAL
!
!*******************************************************************

      implicit none

      integer(4),intent(in):: n
      integer(4):: factor,i

!**************************************

      factor = 1
      do i = 1,n
        factor = factor * i
      enddo

!**************************************

      return
      end function factor
