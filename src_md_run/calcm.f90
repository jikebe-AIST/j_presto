
      subroutine calcm(numatm,weight,center,ier)
 
!*******************************************************************
!
!     CALCULATE CENTER OF MASS
!
!*******************************************************************

      use COMCMMC, only: cord

      implicit none

      integer(4),intent(in):: numatm
      real(8),intent(in):: weight(numatm)
      real(8),intent(out):: center(3)
      integer(4),intent(out):: ier

      integer(4):: i
      real(8):: totwet
 
!**********************************************************

      ier = 0 ; center(1:3) = 0.d0
      totwet = sum(weight(1:numatm))
      do i = 1,numatm
        center(1:3) = center(1:3) + weight(i)*cord(1:3,i)
      enddo

      if ( totwet .gt. 0.d0 ) then
        totwet = 1.d0 / totwet
        center(1:3) = center(1:3) * totwet
      else
        ier = -1
      endif

!*********************************

      return
      end subroutine calcm
