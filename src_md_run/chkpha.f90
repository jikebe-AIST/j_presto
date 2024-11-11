
      subroutine chkpha(numtor,atmlis,phase,eps,iprint,ier)

!*******************************************************************
!
!     CHECK PHASE OF TORSIONAL PARAMETER
!
!*******************************************************************

      use PHYCNS, only: pi

      ! Number of (improper) torsion
        integer(4),intent(in):: numtor
      ! Atom list of torsion
        integer(4),intent(in):: atmlis(numtor,4)
      ! Phase of torsion (radian)
        real(8),intent(inout):: phase(numtor)
      ! Small value ( if difference bet. phase & 0.d0 or pi is less
      !               equal eps, then reset phase to 0.0 or pi)
        real(8),intent(in):: eps
      ! Logical unit number for output log
        integer(4),intent(in):: iprint
      ! Condition code (0: NO ERROR, -1: PHASE ERROR)
        integer(4),intent(out):: ier 

      integer(4):: i
 
!******************************************

      ier = 0
      ! PHASE CONVERSION ( 0 <= phase <= 2*pi)
      do i = 1,numtor
        call covang(phase(i),1)
      enddo

      ! CHECK AND CONVERT PHASE
      do i = 1,numtor
        if ( abs(phase(i)) .le. eps ) then
          phase(i) = 0.d0
        elseif ( abs(phase(i)-pi) .le. eps ) then
          phase(i) = pi
        endif
      enddo

!***********************************

      return
      end subroutine chkpha
