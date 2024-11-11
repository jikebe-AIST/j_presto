
      subroutine energy(eneval)

!***********************************************************************
!
!     CALCULATION OF ENERGY
!
!***********************************************************************
 
      use COMBAS ; use COMERG 

      implicit none

      real(8),intent(inout):: eneval(maxene)
 
!******************************************
 
!     <<<  INITIALIZATION  >>>
      eneval(1:maxene) = 0.d0
 
!     <<<  BASIC ENERGY  >>>
!          BOND  ANGLE  TORSION  IMPROPER  1-4 VAN DER WAALS
!          1-4 ELECTROSTATIC
!          1-5 VAN DER WAALS  1-5 ELECTROSTATIC
!          HYDROGEN BONDED
      call energv(eneval)
 
!     <<<  CONSTRAINED ENERGY  >>>
!          POSITION RESTRAINTS
!          DISTANCE RESTRAINTS
!          TORSION RESTRAINTS
!          SOFT REPULSION
!          CAP RESTRAINTS
      call eneres(eneval)
 
!     <<<  SUM-UP ENERGY  >>>
      eneval(1) = sum(eneval(2:maxene))

!************************************

      return
      end subroutine energy
