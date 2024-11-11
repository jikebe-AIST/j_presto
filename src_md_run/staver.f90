
      subroutine staver(iprint,ier,convcv,deltt,numfre,epot,idfstp,    &
                        idfvsl,temps,taut,lamda)
 
!*******************************************************************
!
!     INITIAL STEP OF MOLECULAR DYNAMICS VERLET METHOD (LEAP FROG)
!
!*******************************************************************
 
      use COMBAS ; use COMERG ; use PHYCNS ; use COMCMMC ; use COMCMM

      implicit none

      integer(4),intent(in):: iprint,numfre,idfstp,idfvsl
      integer(4),intent(out):: ier
      real(8),intent(in):: convcv,deltt,temps,taut
      real(8),intent(inout):: lamda,epot(maxene)

      real(8):: cordm1(3,ixnatm),tempc,ek,epotcc(maxene,nfrg),rtmp
      integer(4):: i,j

!      logical(4),save:: first_touch = .true.
      logical(4),save:: first_touch = .false.
      integer(4):: jjj(1)
      real(8),allocatable:: tg(:,:),Temp_epot(:),Tepot(:)

!**********************************************************
 
      ier = 0
      cordm1(1:3,1:ixnatm) = cord(1:3,1:ixnatm)
 
      ! <<<  CALCULATION OF TEMERATURE AT TIME ZERO  >>>
      if ( idfvsl .eq. 1 ) then
        call caltmp(iynvar,iylvar,vel,fxmass(1:ixnatm),numfre,tempc,ek,&
                    lambda_v,lambda_m)
      endif

      ! <<<  CALCULATION OF NO-CONSTRAINED VELOCITY AT TIME ZERO  >>>
      !      EPOT IS POTENTIAL ENERGY AT TIME T
      call dnergy(epot,epotcc)
      if ( cluster_method .eq. "LMD" ) then
        rtmp = lambda*lambda
        select case (nfrg)
        case (1)
          !$OMP parallel default(none)                               & !
          !$OMP private(i)                                           & !
          !$OMP shared(grad,rtmp,ixnatm)
          !$OMP do schedule (static)
          do i = 1,ixnatm
            grad(:,i,1) = grad(:,i,1)*rtmp
          enddo
          !$OMP end do
          !$OMP end parallel
        case (3)
          !$OMP parallel default(none)                               & !
          !$OMP private(i)                                           & !
          !$OMP shared(grad,rtmp,lambda,ixnatm)
          !$OMP do schedule (static)
          do i = 1,ixnatm
            grad(:,i,3) = grad(:,i,1)*rtmp +                           &
                          grad(:,i,2)*lambda + grad(:,i,3)
          enddo
          !$OMP end do
          !$OMP end parallel
        end select
      elseif ( cluster_method .eq. "LMD2" ) then
        select case (nfrg)
        case (1)
          !$OMP parallel default(none)                               & !
          !$OMP private(i)                                           & !
          !$OMP shared(grad,lambda,ixnatm)
          !$OMP do schedule (static)
          do i = 1,ixnatm
            grad(:,i,1) = grad(:,i,1)*lambda
          enddo
          !$OMP end do
          !$OMP end parallel
        case (3)
          rtmp = sqrt(lambda)
          !$OMP parallel default(none)                               & !
          !$OMP private(i)                                           & !
          !$OMP shared(grad,rtmp,lambda,ixnatm)
          !$OMP do schedule (static)
          do i = 1,ixnatm
            grad(:,i,3) = grad(:,i,1)*lambda +                         &
                          grad(:,i,2)*rtmp + grad(:,i,3)
          enddo
          !$OMP end do
          !$OMP end parallel
        end select

      elseif ( psetmp .gt. 1.0d-10 ) then
        rtmp = temps / psetmp
        select case (nfrg)
        case (1)
          grad(:,:,1) = grad(:,:,1) * rtmp
        case (2)
          grad(:,:,2) = grad(:,:,1)*rtmp + grad(:,:,2)
        end select
      endif

      !!! Force check for programing
      if ( first_touch ) then
!        open(unit=111,file="1stp_force.bin",          &
!             form="unformatted",status="replace")
!        write(111)epot
!        write(111)grad(1:3,1:ixnatm,nfrg)
!        close(111) ; stop
        open(unit=111,file="1stp_force.bin",      &
             form="unformatted",status="old")
        allocate(tg(3,ixnatm),Temp_epot(maxene),Tepot(maxene))
        read(111)Temp_epot(1:maxene)
        read(111)tg(1:3,1:ixnatm)
        Tepot(:) = Temp_epot(:)
        Temp_epot(:) = abs(Temp_epot(:) - epot(:))
        tg(1:3,1:ixnatm) = abs( tg(1:3,1:ixnatm) -                   &
                                grad(1:3,1:ixnatm,nfrg))
        rtmp = maxval(tg(1:3,1:ixnatm))
        print*,""
        print*,"***************************************"
        print*,"*        force error calculation      *"
        print*,"***************************************"
        print*,"  Maximum Force error    = ",rtmp
        jjj = maxloc(Temp_epot(2:)) ; j = jjj(1) + 1
        print*,"  Maximum Energy error   = ",Temp_epot(j)
        print*,"                            (",cyenam(j),")"
        print*,"  Threshold value        = ",1.0d-7
        print*,""
        if ( rtmp .gt. 1.0d-7 ) then
          print*,"!! Force error was LARGER than the threshold !!"
          do i = 1,ixnatm
            print*,i,sqrt(dot_product(tg(1:3,i),tg(1:3,i)))
          enddo
        else
          print*,"* Force error was less than the threshold *"
        endif
        write(iprint,*)
        j = 0
        write(iprint,'(10x,3(a15,x))'),"Energy","preEnergy","diff."
        do i = 1,maxene
          if ( iyeflg(i) .eq. 1 ) then
            write(iprint,'("  ",a6," :",3(e15.7,x),$)')cyenam(i),      &
              epot(i),Tepot(i),epot(i)-Tepot(i)
            if ( Temp_epot(i) .gt. 1.0d-7) then
              write(iprint,'(a2)')" *"
            else
              write(iprint,*)
            endif
          endif
        enddo
        if ( maxval(Temp_epot(2:)) .gt. 1.0d-7 ) then
          print*,"!! Energy error was LARGER than the threshold !!"
        else
          print*,"* Energy error was less than the threshold *"
        endif
        close(111)
        first_touch = .false. ; deallocate(tg,Temp_epot,Tepot) ; stop
      endif

      rtmp = -0.5d0 * deltt * convcv
      do i = 1,iynvar
        j = iylvar(i)
        vel(1:3,i) = vel(1:3,i) + rtmp * grad(1:3,j,nfrg) * ifxmass(j)
      enddo

      ! <<<  STOP CENTER OF MASS  >>>
      select case ( idfstp )
      case (1:3)
        call remotv(ixcend(nstpcn),fxmass(1:ixnatm),cord,ier,idfstp)
        if ( ier .ne. 0 ) then
          write(iprint,*)'ERROR> MD '
          write(iprint,*)'   ERROR OCCURES IN REMOTV '
          if ( ier .eq. -1 ) then
            write(iprint,*)'   MASS WEIGHT OF FIRST CHAIN IS ZERO '
          elseif ( ier .eq. -2 ) then
            write(iprint,*)'   MOMENTA OF FIRST CHAIN IS SINGULAR '
          endif
          return
        endif
      end select
 
      ! <<<  SCALING VELOCITY FOR CONSTANT TEMPERATURE  >>>
      if ( idfvsl .eq. 1 ) then
        call scalet(temps,tempc,taut,deltt,iprint,iynvar,lamda)
      else
        lamda = 0.d0
      endif

      ! <<<  CALCULATION OF NO-CONSTRAINED COORDINATE AT TIME DT  >>>
      do i = 1,iynvar
        j = iylvar(i)
        cord(1:3,j) = cord(1:3,j) + deltt * vel(1:3,i)
      enddo
 
      ! <<<  CALCULATION OF CONSTRAINED COORDINATE AT TIME DT  >>>
      ! <<<  CALCULATION OF CONSTRAINED VELOCITY AT TIME (1/2)DT  >>>
      if ( iyfshk .ne. 0 ) then
        select case(ixfbou)
        case (1)
          call shake_PB(cord,cordm1,iprint,ier)
        case default
          call shake_CB(cord,cordm1,iprint,ier)
        end select
        if ( ier .ne. 0 ) return
      endif

!***********************************************

      return
      end subroutine staver 
