
      subroutine shake_PB(r1,r2,iprint,ier) ! #SLC2 #PB
      subroutine shake_CB(r1,r2,iprint,ier) ! #SLC2 #CB
 
!***********************************************************************
!
!     Constraint molecular dynamics for bonds involving
!       hydrogen atoms using simultaneous linear equations
!
!***********************************************************************

      use COMBAS ; use COMMIS
!!      use CALC_TIME

      implicit none

      ! Coordinates at T+dT which are to be modified by SHAKE ALGORITHM
        real(8),intent(inout):: r1(3,ixnatm)
      ! Coordinates at T (reference coordinates)
        real(8),intent(in):: r2(3,ixnatm)
      ! Logical unit number for output log
        integer(4),intent(in):: iprint
      ! Condition code 
      !  ( 0:NO ERROR, 1:NOT CONVERGED, -2:MATRIX IS SIGULAR)
        integer(4),intent(inout):: ier
 
!****************************************

!!      call system_clock(tim1)
      ier = 0
      ! Two particles
      if ( iugsk2 .gt. 0 ) then
        call shake2_PB(r1,r2,iprint,ier) ! #SLC2 #PB
        call shake2_CB(r1,r2,iprint,ier) ! #SLC2 #CB
        if ( ier .ne. 0 ) return
      endif

      ! Three particles
      if ( iugsk3-iugsk2 .gt. 0 ) then
        call shake3_PB(r1,r2,iprint,ier) ! #SLC2 #PB
        call shake3_CB(r1,r2,iprint,ier) ! #SLC2 #CB
        if ( ier .ne. 0 ) return
      endif
 
      ! Four particles
      if ( iugshk-iugsk3 .gt. 0 ) then
        call shake4_PB(r1,r2,iprint,ier) ! #SLC2 #PB
        call shake4_CB(r1,r2,iprint,ier) ! #SLC2 #CB
        if ( ier .ne. 0 ) return
      endif

!!      call system_clock(tim2)
!!      tmptime = tmptime + tim2 - tim1

!***************************

      return
      end subroutine shake_PB ! #SLC2 #PB
      end subroutine shake_CB ! #SLC2 #CB


!======================================================================= 


      subroutine shake2_PB(r1,r2,iprint,ier) ! #SLC2 #PB
      subroutine shake2_CB(r1,r2,iprint,ier) ! #SLC2 #CB

!****************************************************************
!
!     Constant molecular dynamics for one bond ( two atoms )
!       1) SHAK2E : Exact solution is given
!                   This subroutine can be used only when all
!                   constraints are independent each other
!       2) SHAK2I : Iteration method is used
!                   This subroutine can be used for all types of
!                   constraints
!
!****************************************************************

      use COMBAS ; use COMERG

      implicit none

      ! Coordinates at T+dT which are to be modified by SHAKE ALGORITHM
        real(8),intent(inout):: r1(3,ixnatm)
      ! Coordinates at T (reference coordinates)
        real(8),intent(in):: r2(3,ixnatm)
      ! Logical unit number for output
        integer(4),intent(in):: iprint
      ! Condition code ( 0 : NO ERROR, 1 : NOT CONVERGED)
        integer(4),intent(inout):: ier

!****************************************************

      if ( iyfshk .eq. 1 ) then
        call shak2e_PB(r1,r2)            ! #SLC2 #PB 
        call shak2e_CB(r1,r2)            ! #SLC2 #CB 
      elseif ( iyfshk .eq. 2 ) then
        call shak2i_PB(r1,r2,iprint,ier) ! #SLC2 #PB
        call shak2i_CB(r1,r2,iprint,ier) ! #SLC2 #CB
      endif

!************************

      return
      end subroutine shake2_PB ! #SLC2 #PB
      end subroutine shake2_CB ! #SLC2 #CB


!======================================================================= 


      subroutine shake3_PB(r1,r2,iprint,ier) ! #SLC2 #PB
      subroutine shake3_CB(r1,r2,iprint,ier) ! #SLC2 #CB

!***********************************************************************
!
!     Constraint molecular dynamics for three bonds (three atoms)
!     Newton method is applied
!     Exact inverse matrix (3x3) is given
!
!***********************************************************************

      use COMBAS ; use COMMIS
!      use CALC_TIME
      !$ use omp_lib

      implicit none

      ! Coordinates at T+dT which are to be modified by SHAKE ALGORITHM
        real(8),intent(inout):: r1(3,ixnatm)
      ! Coordinates at T (reference coordinates)
        real(8),intent(in):: r2(3,ixnatm)
      ! Logical unit number for output
        integer(4),intent(in):: iprint
      ! Condition code ( 0: NO ERROR, 1: NOT CONVERGED, -2: Singular
      !                  maxtix)
        integer(4),intent(inout):: ier

      integer(4):: numchk,numerr,ierr,iloop,igrp
      ! Flag for conversion 
      ! (iok(igrp) = .true. , then igrp-th constraint is converged)
        logical(1):: iok(iugsk2+1:iugsk3)
      real(8),save,allocatable:: invmas(:,:),coeg(:,:,:)
      ! Distance vector of r1 & r2
        real(8):: d1(3,3),d2(3,3)
      ! Gradient matrix of SHAKE constraint
        real(8):: gradf(3,3)
      real(8):: deta,detb(3,3)
      ! Coefficient of shake constraint
      !   G(1-3): difference of coefficient
      !   G(4-6): coefficient
        real(8):: G(6)
      ! Right hand vector
        real(8):: B(3)
      ! List of constraint number whose gradient matrix is singular
        integer(4):: ierlis(iugsk3-iugsk2)
      integer(4):: Tiuashk(3)
      real(8):: Tcoeg(3,3),Tfudshk(3),TTfudshk(3)
      
!*************************************

      ! If first loop ...
      if ( .not. allocated(invmas) ) then
        !! Make invmas & coeg
        allocate(invmas(3,iugsk2+1:iugsk3),coeg(3,3,iugsk2+1:iugsk3))

        do igrp = iugsk2+1,iugsk3
          invmas(1:3,igrp) = ifxmass(iuashk(1:3,igrp))
          coeg(1,1,igrp) = invmas(1,igrp) + invmas(2,igrp)
          coeg(1,2,igrp) = - invmas(2,igrp)
          coeg(1,3,igrp) = - invmas(1,igrp)
          coeg(2,1,igrp) = - invmas(2,igrp)
          coeg(2,2,igrp) = invmas(2,igrp) + invmas(3,igrp)
          coeg(2,3,igrp) = - invmas(3,igrp)
          coeg(3,1,igrp) = - invmas(1,igrp)
          coeg(3,2,igrp) = - invmas(3,igrp)
          coeg(3,3,igrp) = invmas(3,igrp) + invmas(1,igrp)
        enddo
      endif

!******************

!      call system_clock(tim1)
      numerr = 0
      !$OMP parallel default(none)                                   & !
      !$OMP private(iloop,igrp,Tiuashk,Tcoeg,Tfudshk,TTfudshk,gradf, & !
      !$OMP         B,deta,detb,d1,d2,G)                             & !
      !$OMP shared(iuslop,iugsk2,iugsk3,iuashk,r1,r2,fxcell,invcel,  & !
      !$OMP        iok,coeg,fudshk,fustol,numerr,ierlis,invmas)
      !$OMP do
      do igrp = iugsk2+1,iugsk3
        Tcoeg(:,:) = coeg(:,:,igrp) ; Tiuashk(1:3) = iuashk(1:3,igrp)
        Tfudshk(1:3) = fudshk(1:3,igrp)
        TTfudshk(:) = fustol*Tfudshk(:)
        G(:) = 0.d0 ; iok(igrp) = .false.
        ! Calculate distance vector
        d1(1,1:3) = r1(1:3,Tiuashk(1))-r1(1:3,Tiuashk(2))
        d1(2,1:3) = r1(1:3,Tiuashk(2))-r1(1:3,Tiuashk(3))
        d1(3,1:3) = r1(1:3,Tiuashk(3))-r1(1:3,Tiuashk(1))
        d2(1,1:3) = r2(1:3,Tiuashk(1))-r2(1:3,Tiuashk(2))
        d2(2,1:3) = r2(1:3,Tiuashk(2))-r2(1:3,Tiuashk(3))
        d2(3,1:3) = r2(1:3,Tiuashk(3))-r2(1:3,Tiuashk(1))
        d1(1,1:3) = d1(1,1:3) -                             & ! #SLC2 #PB
                    fxcell(1:3)*nint(d1(1,1:3)*invcel(1:3))   ! #SLC2 #PB
        d1(2,1:3) = d1(2,1:3) -                             & ! #SLC2 #PB
                    fxcell(1:3)*nint(d1(2,1:3)*invcel(1:3))   ! #SLC2 #PB
        d1(3,1:3) = d1(3,1:3) -                             & ! #SLC2 #PB
                    fxcell(1:3)*nint(d1(3,1:3)*invcel(1:3))   ! #SLC2 #PB
        d2(1,1:3) = d2(1,1:3) -                             & ! #SLC2 #PB
                    fxcell(1:3)*nint(d2(1,1:3)*invcel(1:3))   ! #SLC2 #PB
        d2(2,1:3) = d2(2,1:3) -                             & ! #SLC2 #PB
                    fxcell(1:3)*nint(d2(2,1:3)*invcel(1:3))   ! #SLC2 #PB
        d2(3,1:3) = d2(3,1:3) -                             & ! #SLC2 #PB
                    fxcell(1:3)*nint(d2(3,1:3)*invcel(1:3))   ! #SLC2 #PB

        ! Calculation of G-value
        do iloop = 1,iuslop
          !! 1) Calculate difference of distance**2
          B(1:3) = Tfudshk(1:3) - ( d1(1:3,1)*d1(1:3,1) +              &
                     d1(1:3,2)*d1(1:3,2) + d1(1:3,3)*d1(1:3,3) )
 
          !! 2) Check convergence
          if ( abs(B(1)).lt.TTfudshk(1) .and. abs(B(2)).lt.TTfudshk(2) &
                 .and. abs(B(3)).lt.TTfudshk(3) ) then
            iok(igrp) = .true.
            ! Modify r1 coordinate
            r1(1:3,Tiuashk(1)) = r1(1:3,Tiuashk(1)) + invmas(1,igrp) * &
                                   ( G(4)*d2(1,1:3)-G(6)*d2(3,1:3))
            r1(1:3,Tiuashk(2)) = r1(1:3,Tiuashk(2)) + invmas(2,igrp) * &
                                   (-G(4)*d2(1,1:3)+G(5)*d2(2,1:3))
            r1(1:3,Tiuashk(3)) = r1(1:3,Tiuashk(3)) + invmas(3,igrp) * &
                                   (-G(5)*d2(2,1:3)+G(6)*d2(3,1:3))
            exit
          else
            !! 3) Gradient matrix gradf(icon,ig) = df(icon)/dg(ig)
            gradf(1:3,1) = 2.d0*Tcoeg(1:3,1) * ( d1(1:3,1)*d2(1,1) +   &
                             d1(1:3,2)*d2(1,2) + d1(1:3,3)*d2(1,3))
            gradf(1:3,2) = 2.d0*Tcoeg(1:3,2) * ( d1(1:3,1)*d2(2,1) +   &
                             d1(1:3,2)*d2(2,2) + d1(1:3,3)*d2(2,3))
            gradf(1:3,3) = 2.d0*Tcoeg(1:3,3) * ( d1(1:3,1)*d2(3,1) +   &
                             d1(1:3,2)*d2(3,2) + d1(1:3,3)*d2(3,3))

            !! 4) Solve linear equation ( gradf * G = B )
            deta = ( gradf(1,1)*gradf(2,2)*gradf(3,3) )                &
                 + ( gradf(1,2)*gradf(2,3)*gradf(3,1) )                &
                 + ( gradf(1,3)*gradf(2,1)*gradf(3,2) )                &
                 - ( gradf(1,1)*gradf(2,3)*gradf(3,2) )                &
                 - ( gradf(1,2)*gradf(2,1)*gradf(3,3) )                &
                 - ( gradf(1,3)*gradf(2,2)*gradf(3,1) )
            if ( deta .eq. 0.d0 ) then
              !$OMP critical
              numerr = numerr + 1 ; ierlis(numerr) = igrp
              !$OMP end critical
            else
              deta = 1.d0 / deta
            endif

            detb(1,1) = gradf(2,2)*gradf(3,3) - gradf(2,3)*gradf(3,2)
            detb(1,2) = gradf(3,2)*gradf(1,3) - gradf(3,3)*gradf(1,2)
            detb(1,3) = gradf(1,2)*gradf(2,3) - gradf(1,3)*gradf(2,2)
            detb(2,1) = gradf(2,3)*gradf(3,1) - gradf(2,1)*gradf(3,3)
            detb(2,2) = gradf(3,3)*gradf(1,1) - gradf(3,1)*gradf(1,3)
            detb(2,3) = gradf(1,3)*gradf(2,1) - gradf(1,1)*gradf(2,3)
            detb(3,1) = gradf(2,1)*gradf(3,2) - gradf(2,2)*gradf(3,1)
            detb(3,2) = gradf(3,1)*gradf(1,2) - gradf(3,2)*gradf(1,1)
            detb(3,3) = gradf(1,1)*gradf(2,2) - gradf(1,2)*gradf(2,1)
            G(1:3) = deta * ( detb(1:3,1)*B(1) + detb(1:3,2)*B(2) +    &
                              detb(1:3,3)*B(3) )
            G(4:6) = G(4:6) + G(1:3)

            !! 5) Update distance vector d1
            d1(1,1:3) = d1(1,1:3) + Tcoeg(1,1)*d2(1,1:3)*G(1) +        &
                  Tcoeg(1,2)*d2(2,1:3)*G(2) + Tcoeg(1,3)*d2(3,1:3)*G(3)
            d1(2,1:3) = d1(2,1:3) + Tcoeg(2,1)*d2(1,1:3)*G(1) +        &
                  Tcoeg(2,2)*d2(2,1:3)*G(2) + Tcoeg(2,3)*d2(3,1:3)*G(3)
            d1(3,1:3) = d1(3,1:3) + Tcoeg(3,1)*d2(1,1:3)*G(1) +        &
                  Tcoeg(3,2)*d2(2,1:3)*G(2) + Tcoeg(3,3)*d2(3,1:3)*G(3)
          endif
        enddo

      enddo
      !$OMP end do
      !$OMP end parallel

      if ( numerr .ne. 0 ) then
        ier = -2
        write(iprint,*)"ERROR> SHAKE (V3.0) "
        write(iprint,*)"     ",numerr," GRADIENT MATRIX IS SINGULAR "
        write(iprint,*)"     ",iloop," TH LOOP "
!        do ierr = 1,numerr
!          igrp = ierlis(ierr)
          igrp = ierlis(1)
!          write(iprint,*)"  ",igrp," -TH MATRIX IS SINGULAR "
          write(iprint,*)"     ATOM ",iuashk(1,igrp)," - ",          &
                         iuashk(2,igrp)," - ",iuashk(3,igrp)
!        enddo
        write(iprint,*)" "
        return
      endif

      if ( .not. all(iok(iugsk2+1:iugsk3)) ) then
        numchk = count(iok(iugsk2+1:iugsk3))
        ier = 1
        write(iprint,*)" "
        write(iprint,*)"WARNING> SHAKE (V3.0) "
        write(iprint,*)"     ",iugsk3-iugsk2-numchk," CONSTRAINTS "
        write(iprint,*)"       ARE NOT CONVERGED"
        do igrp = iugsk2+1,iugsk3
          if ( .not. iok(igrp) ) then
            write(iprint,*)"  ",igrp," -TH CONSTRAINT IS NOT CONVERGED "
            write(iprint,*)"     ATOM ",iuashk(1,igrp)," - ",          &
                           iuashk(2,igrp)," - ",iuashk(3,igrp)
            exit
          endif
        enddo
        write(iprint,*)" "
        return
      endif
!      call system_clock(tim2)
!      tmptime = tmptime + tim2 - tim1

!*******************************

      return 
      end subroutine shake3_PB ! #SLC2 #PB
      end subroutine shake3_CB ! #SLC2 #CB


!======================================================================= 


      subroutine shake4_PB(r1,r2,iprint,ier) ! #SLC2 #PB
      subroutine shake4_CB(r1,r2,iprint,ier) ! #SLC2 #CB

!***********************************************************************
!
!     Constraint molecular dynamics for six bonds ( four atoms )
!       Newton method (iteration) is applied
!
!***********************************************************************

      use COMBAS ; use COMMIS ; use PHYCNS
      !$ use omp_lib

      implicit none

      real(8),intent(inout):: r1(3,ixnatm)
      real(8),intent(in):: r2(3,ixnatm)
      integer(4),intent(in):: iprint
      integer(4),intent(inout):: ier

      ! Atom numebr list in each SHAKE constraint
        integer(4):: atmlis(mxashk)
      ! Atom pair list in each SHAKE constraint ( atom number are input)
        integer(4):: parlis(maxequ,2)
      ! Inverse of mass & Coefficient matrix for G
        real(8),save,allocatable:: invmas(:,:),coeg(:,:,:)
        real(8):: Tcoeg(maxequ,maxequ),Tfudshk(maxequ),TTfudshk(maxequ)
      ! Distance vector of r1 & r2, & work distance vector of r1
        real(8):: wkd1(3,maxequ),d1(3,maxequ),d2(3,maxequ)
      ! Coefficient of SHAKE constraint ( & for work)
        real(8):: G(maxequ)
      ! Gradient matrix of SHAKE cons.
        real(8):: gradf(maxequ,maxequ)
      ! Right vector for calc. of G
        real(8):: B(maxequ)
      logical(1):: flg,iok(iugsk3+1:iugshk)

      integer(4):: igrp
      data parlis/ 1,2,3,4,3,2,                                        &
                   2,3,1,1,4,4 /
      ! If diagonal element is less than epsl,
      !    then this matrix is singular (SIMULU)
        real(8):: epsl = eps*eps*eps
      ! Index vector for partial pivoting (SIMULU)
        integer(4):: ip(maxequ)
      integer(4):: ierlis(iugshk-iugsk3)

      integer(4):: numatm,numcon,numchk,iloop,icon,iatm,jatm,ig,       &
                   numerr,ierr
 
!***************************************

      ! If first loop ...
      if ( .not. allocated(invmas) ) then
        !! Make invmas & coeg
        allocate(invmas(mxashk,iugsk3+1:iugshk),                       &
                 coeg(maxequ,maxequ,iugsk3+1:iugshk))
        ! 1) Calculation of coefficient of G
        do igrp = iugsk3+1,iugshk
          atmlis(1:4) = iuashk(1:4,igrp)
          invmas(1:4,igrp) = 1.d0 / fxmass(atmlis(1:4))
 
          coeg(1,1,igrp) =   invmas(1,igrp) + invmas(2,igrp)
          coeg(1,2,igrp) = - invmas(2,igrp)
          coeg(1,3,igrp) = - invmas(1,igrp)
          coeg(1,4,igrp) = - invmas(1,igrp)
          coeg(1,5,igrp) =   0.d0
          coeg(1,6,igrp) = - invmas(2,igrp)
          coeg(2,1,igrp) = - invmas(2,igrp)
          coeg(2,2,igrp) =   invmas(2,igrp) + invmas(3,igrp)
          coeg(2,3,igrp) = - invmas(3,igrp)
          coeg(2,4,igrp) =   0.d0
          coeg(2,5,igrp) = - invmas(3,igrp)
          coeg(2,6,igrp) =   invmas(2,igrp)
          coeg(3,1,igrp) = - invmas(1,igrp)
          coeg(3,2,igrp) = - invmas(3,igrp)
          coeg(3,3,igrp) =   invmas(3,igrp) + invmas(1,igrp)
          coeg(3,4,igrp) =   invmas(1,igrp)
          coeg(3,5,igrp) =   invmas(3,igrp)
          coeg(3,6,igrp) =   0.d0
          coeg(4,1,igrp) = - invmas(1,igrp)
          coeg(4,2,igrp) =   0.d0
          coeg(4,3,igrp) =   invmas(1,igrp)
          coeg(4,4,igrp) =   invmas(4,igrp) + invmas(1,igrp)
          coeg(4,5,igrp) = - invmas(4,igrp)
          coeg(4,6,igrp) = - invmas(4,igrp)
          coeg(5,1,igrp) =   0.d0
          coeg(5,2,igrp) = - invmas(3,igrp)
          coeg(5,3,igrp) =   invmas(3,igrp)
          coeg(5,4,igrp) = - invmas(4,igrp)
          coeg(5,5,igrp) =   invmas(3,igrp) + invmas(4,igrp)
          coeg(5,6,igrp) =   invmas(4,igrp)
          coeg(6,1,igrp) = - invmas(2,igrp)
          coeg(6,2,igrp) =   invmas(2,igrp)
          coeg(6,3,igrp) =   0.d0
          coeg(6,4,igrp) = - invmas(4,igrp)
          coeg(6,5,igrp) =   invmas(4,igrp)
          coeg(6,6,igrp) =   invmas(2,igrp) + invmas(4,igrp)
        enddo
      endif

      numerr = 0
      !$OMP parallel default (none)                                  & !
      !$OMP private(igrp,atmlis,numatm,numcon,Tcoeg,Tfudshk,TTfudshk,& !
      !$OMP         icon,iatm,jatm,d1,d2,iloop,G,ig,gradf,B,flg,ip,  & !
      !$OMP         wkd1,ierr)                                       & !
      !$OMP shared(iugsk3,iugshk,iuashk,iuhshk,fustol,parlis,r1,r2,  & !
      !$OMP        fxcell,invcel,iuslop,coeg,fudshk,iok,invmas,epsl, & !
      !$OMP        numerr,ierlis)
      !$OMP do
      do igrp = iugsk3+1,iugshk
        ! 2) Calculation of G
        atmlis(1:4) = iuashk(1:4,igrp)
        numatm = iuhshk(igrp) + 1
        numcon = (numatm*(numatm-1)) / 2
        Tcoeg(:,:) = coeg(:,:,igrp)
        Tfudshk(1:numcon) = fudshk(1:numcon,igrp)
        TTfudshk(1:numcon) = fustol*Tfudshk(1:numcon)
        G(1:numcon) = 0.d0 ; iok(igrp) = .false.

        !! Distance vector
        do icon = 1,numcon
          iatm = atmlis(parlis(icon,1)) ; jatm = atmlis(parlis(icon,2))
          d1(1:3,icon) = r1(1:3,iatm) - r1(1:3,jatm)
          d2(1:3,icon) = r2(1:3,iatm) - r2(1:3,jatm)
          d1(1:3,icon) = d1(1:3,icon) -                           & ! #SLC2 #PB
                         fxcell(1:3)*nint(d1(1:3,icon)*invcel(1:3)) ! #SLC2 #PB
          d2(1:3,icon) = d2(1:3,icon) -                           & ! #SLC2 #PB
                         fxcell(1:3)*nint(d2(1:3,icon)*invcel(1:3)) ! #SLC2 #PB
        enddo
        wkd1(1:3,1:numcon) = d1(1:3,1:numcon)

!**************************

        ! Iteration of newton method start
        do iloop = 1,iuslop
          !! Check convergence
          B(1:numcon) = Tfudshk(1:numcon) -                            &
                      ( wkd1(1,1:numcon)*wkd1(1,1:numcon) +            &
                        wkd1(2,1:numcon)*wkd1(2,1:numcon) +            &
                        wkd1(3,1:numcon)*wkd1(3,1:numcon) )
          flg = .true.
          do icon = 1,numcon
            if ( abs(B(icon)) .ge. TTfudshk(icon) ) flg = .false.
          enddo
          if ( flg ) then
            iok(igrp) = .true.
            ! 3) Modify r1 coordinates
            r1(1:3,atmlis(1)) = r1(1:3,atmlis(1)) + invmas(1,igrp) *   &
               (G(1)*d2(1:3,1) - G(3)*d2(1:3,3) - G(4)*d2(1:3,4))
            r1(1:3,atmlis(2)) = r1(1:3,atmlis(2)) + invmas(2,igrp) *   &
               (-G(1)*d2(1:3,1) + G(2)*d2(1:3,2) + G(6)*d2(1:3,6))
            r1(1:3,atmlis(3)) = r1(1:3,atmlis(3)) + invmas(3,igrp) *   &
               (-G(2)*d2(1:3,2) + G(3)*d2(1:3,3) + G(5)*d2(1:3,5))
            r1(1:3,atmlis(4)) = r1(1:3,atmlis(4)) + invmas(4,igrp) *   &
               (  G(4)*d2(1:3,4) - G(5)*d2(1:3,5) - G(6)*d2(1:3,6))
            exit
          else
            !! Gradient matrix df(icon)/dg(ig) 
            do ig = 1,numcon
            do icon = 1,numcon
              gradf(icon,ig) = 2.d0*Tcoeg(icon,ig) *                   &
                (wkd1(1,icon)*d2(1,ig) + wkd1(2,icon)*d2(2,ig) +       &
                 wkd1(3,icon)*d2(3,ig))
            enddo
            enddo
            !! Right hand vector
            do ig = 1,numcon
            do icon = 1,numcon
              B(icon) = B(icon) + gradf(icon,ig)*G(ig)
            enddo
            enddo
            !! Linear equation solver
            call simulu(maxequ,numcon,gradf,B,G,epsl,1,ip,ierr)
            if ( ierr .ne. 0 ) then
              !$OMP critical
              numerr = numerr + 1 ; ierlis(numerr) = igrp
              !$OMP end critical
            endif

            !! distance d1
            wkd1(1:3,1:numcon) = d1(1:3,1:numcon)
            do ig = 1,numcon
            do icon = 1,numcon
              wkd1(1:3,icon) = wkd1(1:3,icon) +                        &
                             Tcoeg(icon,ig)*d2(1:3,ig)*G(ig)
            enddo
            enddo
          endif
        enddo
      enddo
      !$OMP end do
      !$OMP end parallel
      if ( numerr .ne. 0 ) then
        write(iprint,*)'ERROR> SHAKE (V3.0) '
        write(iprint,*)"     ",numerr," GRADIENT MATRIX IS SINGULAR "
        write(iprint,*)"     ",iloop," TH LOOP "
!        do ierr = 1,numerr
!          igrp = ierlis(ierr)
          igrp = ierlis(1)
!          write(iprint,*)"  ",igrp," -TH MATRIX IS SINGULAR "
          write(iprint,*)"     ATOM ",iuashk(1,igrp)," - ",          &
            iuashk(2,igrp)," - ",iuashk(3,igrp)," - ",iuashk(4,igrp)
!        enddo
        write(iprint,*)' '
        ier = -2 ; return
       endif
       if ( .not. all(iok(iugsk3+1:iugshk)) ) then
         ier = 1 ; numchk = count(iok(iugsk3+1:iugshk))
         write(iprint,*)" "
         write(iprint,*)"WARNIG> SHAKE (V3.0) "
         write(iprint,*)"     ",iugshk-iugsk3-numchk," CONSTRAINTS "
         write(iprint,*)"       ARE NOT CONVERGED"
         do igrp = iugsk3+1,iugshk
           if ( .not. iok(igrp) ) then
             write(iprint,*)"   ",igrp," -TH CONSTRAINT "
             write(iprint,*)"   ATOM ",iuashk(1:4,igrp)
             write(iprint,*)" "
             exit
           endif
         enddo
         write(iprint,*)" "
         return
       endif

!**********************************

      return
      end subroutine shake4_PB ! #SLC2 #PB
      end subroutine shake4_CB ! #SLC2 #CB
 

!=======================================================================


      subroutine shak2i_PB(r1,r2,iprint,ier) ! #SLC2 #PB
      subroutine shak2i_CB(r1,r2,iprint,ier) ! #SLC2 #CB

!******************************************************
!
!     Constraint molecular dynamics for bonds
!     Any bonds can be constrainted
!
!******************************************************

      use COMBAS ; use COMMIS

      implicit none

      ! Coordinates at T+dT which are to be modified by SHAKE ALGORITHM
        real(8),intent(inout):: r1(3,ixnatm)
      ! Coordinates at T (reference coordinates)
        real(8),intent(in):: r2(3,ixnatm)
      ! Logical unit number for output
        integer(4),intent(in):: iprint
      ! Condition code (0: NO ERROR, 1: NOT CONVERGED)
        integer(4),intent(inout):: ier

      integer(4):: igrp,iloop,numchk
      ! Flag for conversion ( .true. : converged)
        logical(1):: iok(iugsk2)
      real(8):: difd2,inpro
      real(8):: d1(3)
      ! Distance vector of r2
        real(8):: d2(3,3,iugsk2)
      real(8),save,allocatable:: invmas(:,:)
      ! Coefficient of SHAKE constraint
        real(8):: G(iugsk2)

!********************************************

      ! If first loop ...
      if ( .not. allocated(invmas) ) then
        !! Make invmas
        allocate(invmas(3,iugsk2))
        do igrp = 1,iugsk2
          invmas(1:2,igrp) = 1.d0 / fxmass(iuashk(1:2,igrp))
          invmas(3,igrp) = 0.5d0 / ( invmas(1,igrp) + invmas(2,igrp) ) 
        enddo
      endif

      do igrp = 1,iugsk2
        d2(1,1:3,igrp) = r2(1:3,iuashk(1,igrp)) - r2(1:3,iuashk(2,igrp))
        d2(1,1:3,igrp) = d2(1,1:3,igrp) -                           & ! #SLC2 #PB
                         fxcell(1:3)*nint(d2(1,1:3,igrp)*invcel(1:3)) ! #SLC2 #PB
      enddo

      ! Iteration part
      do iloop = 1,iuslop

        !! 1) Iteration over constraints
        do igrp = 1,iugsk2
          !!! 1-1) Calculate distance vector
          d1(1:3) = r1(1:3,iuashk(1,igrp)) - r1(1:3,iuashk(2,igrp))
          d1(1:3) = d1(1:3) - fxcell(1:3)*nint(d1(1:3)*invcel(1:3)) ! #SLC2 #PB
          difd2 = d1(1)*d1(1) + d1(2)*d1(2) + d1(3)*d1(3)

          !!! 1-2) Calculate G-value
          difd2 = fudshk(1,igrp) - difd2
          inpro = d1(1)*d2(1,1,igrp) + d1(2)*d2(1,2,igrp) +            &
                  d1(3)*d2(1,3,igrp)
          G(igrp) = ( difd2*invmas(3,igrp)) / inpro

          !!! 1-3) Update r1 coordinates
          r1(1:3,iuashk(1,igrp)) = r1(1:3,iuashk(1,igrp)) +            &
                            (invmas(1,igrp)*G(igrp)*d2(1,1:3,igrp))
          r1(1:3,iuashk(2,igrp)) = r1(1:3,iuashk(2,igrp)) -            &
                            (invmas(2,igrp)*G(igrp)*d2(1,1:3,igrp))
        enddo

        ! Check convergence
        do igrp = 1,iugsk2
          if ( iok(igrp) ) cycle
          d1(1:3) = r1(1:3,iuashk(1,igrp)) - r1(1:3,iuashk(2,igrp))
          d1(1:3) = d1(1:3) - fxcell(1:3)*nint(d1(1:3)*invcel(1:3))
          difd2 = d1(1)*d1(1) + d1(2)*d1(2) + d1(3)*d1(3)
          difd2 = fudshk(1,igrp) - difd2
          difd2 = abs( difd2 / fudshk(1,igrp) )
          if ( difd2 .lt. fustol ) then
            iok(igrp) = .true.
          else
            iok(igrp) = .false.
          endif
        enddo

        numchk = count(iok)
        if ( numchk .eq. iugsk2 ) return

      enddo 

      ier = 1
      write(iprint,*)" "
      write(iprint,*)"WARNIG> SHAKE (V3.0) "
      write(iprint,*)"   ",(iugsk2-numchk)," CONSTRAINTS "
      write(iprint,*)"     ARE NOT CONVERGED "
      do igrp = 1,iugsk2
        if ( .not. iok(igrp) ) then
          write(iprint,*)"  ",igrp," -TH CONSTRAINT IS NOT CONVERGED "
          write(iprint,*)"     ATOM ",iuashk(1,igrp)," - ",            &
                         iuashk(2,igrp)
        endif
      enddo
      write(iprint,*)" "

!*********************************
 
      return
      end subroutine shak2i_PB ! #SLC2 #PB
      end subroutine shak2i_CB ! #SLC2 #CB


!=======================================================================


      subroutine shak2e_PB(r1,r2) ! #SLC2 #PB
      subroutine shak2e_CB(r1,r2) ! #SLC2 #CB

!*********************************************************
!
!     Constant molecular dynamics for one bond (two atoms)
!     Exact solution is given
!
!*********************************************************

      use COMBAS ; use COMMIS
      !$ use omp_lib

      implicit none

      real(8),intent(inout):: r1(3,ixnatm)
      real(8),intent(in):: r2(3,ixnatm) 

      integer(4):: igrp,iatm,jatm
      real(8):: d1(3),d2(3),d1d1,d1d2,d2d2,A,B,C,G
      real(8),save,allocatable:: invmas(:,:)

!*************************************

      ! If first loop ...
      if ( .not. allocated(invmas) ) then
        !! Make invmas
        allocate(invmas(4,iugsk2))
        do igrp = 1,iugsk2
          invmas(1:2,igrp) = 1.d0 / fxmass(iuashk(1:2,igrp))
          invmas(3,igrp) = invmas(1,igrp) + invmas(2,igrp)
          invmas(4,igrp) = invmas(3,igrp) * invmas(3,igrp)
        enddo
      endif

      !$OMP parallel default (none)                                  & !
      !$OMP private(igrp,iatm,jatm,d1,d2,d1d1,d1d2,d2d2,A,B,C,G)     & !
      !$OMP shared(iugsk2,iuashk,r1,r2,fxcell,invcel,invmas,fudshk)
      !$OMP do
      do igrp = 1,iugsk2
        ! 1) Calculate distance vector
        iatm = iuashk(1,igrp) ; jatm = iuashk(2,igrp)
        d1(1:3) = r1(1:3,jatm) - r1(1:3,iatm)
        d2(1:3) = r2(1:3,jatm) - r2(1:3,iatm)
        d1(1:3) = d1(1:3) - fxcell(1:3)*nint(d1(1:3)*invcel(1:3)) ! #SLC2 #PB
        d2(1:3) = d2(1:3) - fxcell(1:3)*nint(d2(1:3)*invcel(1:3)) ! #SLC2 #PB
        d1d1 = d1(1)*d1(1) + d1(2)*d1(2) + d1(3)*d1(3)
        d1d2 = d1(1)*d2(1) + d1(2)*d2(2) + d1(3)*d2(3)
        d2d2 = d2(1)*d2(1) + d2(2)*d2(2) + d2(3)*d2(3)

        ! 2) Solve quadratic equation
        !    A*G**2 + 2*B*G + C = 0
        A = invmas(4,igrp) * d2d2
        B = invmas(3,igrp) * d1d2
        C = d1d1 - fudshk(1,igrp)
        G = ( B - sqrt(B*B - C*A) ) / A

        ! 3) Modify r1 coordinate
        r1(1:3,iatm) = r1(1:3,iatm) + (invmas(1,igrp) * G * d2(1:3))
        r1(1:3,jatm) = r1(1:3,jatm) - (invmas(2,igrp) * G * d2(1:3))
      enddo
      !$OMP end do
      !$OMP end parallel

!************************

      return
      end subroutine shak2e_PB ! #SLC2 #PB
      end subroutine shak2e_CB ! #SLC2 #CB
