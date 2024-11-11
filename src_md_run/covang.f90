
      subroutine covang(dihed,iopt)
 
!*******************************************************************
!
!     THIS SUBROUTINE IS FOR CONVERSION ANGLE
!
!*******************************************************************

      use PHYCNS, only: pi

      implicit none

      ! Angle (radian)
        real(8),intent(inout):: dihed
      ! Option (1: (0,pi2), 2: (-pi,pi))
        integer(4),intent(in):: iopt

      integer(4):: jrot
      real(8),save:: invpi,invpi0_5,pi2
      logical(1),save:: fstflg = .true.
 
!*****************************************************

      ! first call treatment
      if ( fstflg ) then
        invpi = 1.d0 / pi
        invpi0_5 = 0.5d0 * invpi
        pi2 = 2.d0 * pi
        fstflg = .false.
      endif

!     * RANGE OF DIHEDRAL ANGLE IS FROM 0 TO 2PI *
      if ( iopt .eq. 1 ) then
        if ( dihed .ge. 0.d0 ) then
          jrot = int(dihed * invpi0_5)
        else
          jrot = int(dihed * invpi0_5) + 1
        endif
        dihed = dihed - (jrot*pi2)

!     * RANGE OF DIHEDRAL ANGLE IS FROM -PI TO PI *
      elseif ( iopt .eq. 2 ) then
        if ( abs(dihed) .ge. pi ) then
          if ( dihed .eq. -pi ) then
            dihed = pi
          else
            if ( dihed .ne. pi ) then
              jrot = int(dihed*invpi)
              if ( mod(jrot,2) .eq. 0 ) then
                dihed = dihed - (pi*jrot)
              else
                if (dihed .ge. 0.d0 ) then
                  dihed = dihed - (pi*(jrot+1))
                else
                  dihed = dihed - (pi*(jrot-1))
                endif
              endif
            endif
          endif
        endif
      endif

!**************************************

      return
      end subroutine covang
