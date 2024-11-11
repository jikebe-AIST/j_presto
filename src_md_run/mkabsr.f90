
      subroutine mkabsr
 
!*******************************************************************
!
!     MAKE ABSOLUTE RESIDUE NUMBER LIST OF EACH ATOM
!
!*******************************************************************
 
      use COMBAS ; use COMCMM

      implicit none

      integer(4):: ichnold,numadd,iatm,ichn
 
!************************************************

      allocate(absres(ixnatm))
      ichnold = ixachn(1) ; numadd = 0
      do iatm = 1,ixnatm
        ichn = ixachn(iatm)
        if ( ichn .ne. ichnold ) then
          numadd = numadd + ixares(iatm-1)
          ichnold = ichn
        endif
        absres(iatm) = ixares(iatm) + numadd
      enddo

!*************************************

      return
      end subroutine mkabsr
