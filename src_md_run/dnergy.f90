
      subroutine dnergy(eneval,enevcc)

!***********************************************************************
!
!     CALCULATION OF ENERGY AND GRADIENT
!
!***********************************************************************
 
      use COMBAS ; use COMERG
      use COMCMM, only: nfrg
      use COMCMMC,only: grad
      !$ use omp_lib

      implicit none

      ! Energy value (KCAL/MOL)
!         ENEVAL(1)  ; TOTAL ENERGY
!         ENEVAL(2)  ; BOND ENERGY
!         ENEVAL(3)  ; ANGLE ENERGY
!         ENEVAL(4)  ; TORSIONAL ENERGY
!         ENEVAL(5)  ; IMPROPER ENERGY
!         ENEVAL(6)  ; 1-4 VAN DER WAALS ENERGY
!         ENEVAL(7)  ; 1-4 ELECTROSTATIC ENERGY
!         ENEVAL(8)  ; 1-5 VAN DER WAALS ENERGY
!         ENEVAL(9)  ; 1-5 ELECTROSTATIC ENERGY
!         ENEVAL(10) ; HYDROGEN BONDED ENERGY
!         ENEVAL(11) ; 1-5 VAN DER WAALS ENERGY (NO CUT-OFF)
!         ENEVAL(12) ; 1-5 ELECTROSTATIC ENERGY (NO CUT-OFF)
!         ENEVAL(13) ; HYDROGEN BONDED ENERGY   (NO CUT-OFF)
!         ENEVAL(14) ; POSITION RESTRAINTS ENERGY
!         ENEVAL(15) ; DISTANCE RESTRAINTS ENERGY
!         ENEVAL(16) ; TORSIONAL RESTRAINTS ENERGY
!         ENEVAL(17) ; SOFT REPULSION ENERGY
!         ENEVAL(18) ; CAP CONSTARINTS ENERGY
        real(8),intent(out):: eneval(maxene)
        real(8),intent(out):: enevcc(maxene,nfrg)

      integer(4):: i,j

!****************************************************************

      ! <<<  INITIALIZATION  >>>
      enevcc(1:maxene,1:nfrg) = 0.d0
      !$OMP parallel default (none)                                  & !
      !$OMP private(i,j)                                             & !
      !$OMP shared(nfrg,ixnatm,grad)
      !$OMP do schedule (static) collapse(2)
      do j = 1,nfrg
      do i = 1,ixnatm
        grad(:,i,j) = 0.d0
      enddo
      enddo
      !$OMP end do
      !$OMP end parallel

!     dneres <<<  RESTRAINED ENERGY  >>>
!          POSITION RESTRAINTS
!          DISTANCE RESTRAINTS
!          TORSION RESTRAINTS
!          SOFT REPULSION
!          CAP RESTRAINTS
!     dnergv <<<  BASIC ENERGY  >>>
!          BOND  ANGLE  TORSION  IMPROPER  1-4 VAN DER WAALS
!          1-4 ELECTROSTATIC
!          1-5 VAN DER WAALS  1-5 ELECTROSTATIC
!          HYDROGEN BONDED
      select case (ixfbou)
        case(1)
          call dneres_PB(enevcc)
          call dnergv_PB(enevcc)
        case default
          call dneres_CB(enevcc)
          call dnergv_CB(enevcc)
      end select

!     <<<  SUM-UP ENERGY  >>>
      forall ( i=1:nfrg ) enevcc(1,i) = sum(enevcc(2:maxene,i))
      forall ( i=1:maxene ) eneval(i) = sum(enevcc(i,1:nfrg))

!****************************************

      return
      end subroutine dnergy
