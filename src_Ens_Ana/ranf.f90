
      function ranf(iseed)

!*************************************************************
!
!     Random Number Generator (Linear Congruential Method)
!
!*************************************************************

      implicit none

      integer(4),parameter:: l = 32771
      integer(4),parameter:: ic = 1234567891
      real(8),parameter:: rm = 2147483648.0d0

      real(8):: ranf

      integer(4):: iseed,k
      real(8):: w

!***************************

      w = iseed
      w = w * l + ic

      if ( w .ge. rm ) then
        k = w / rm
        w = w - k * rm
      endif
      iseed = w
      ranf = dble(iseed) / rm

!**************************

      return
      end function ranf
