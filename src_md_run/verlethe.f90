
      subroutine verlethe(iprint,ier,convcv,deltt,tmpsta,numfre,       &
                          velback,epot,epotcc,idfstp,idfvsl,temps,     &
                          idloop,idhtlp,lmd_vback)

!********************************************************************
!
!     ONE-STEP NUMERICAL INTEGRATION OF MOLECULAR DYNAMICS
!     LEAP-FROG VERLET
!
!     AVAILABLE MOLECULAR DYNAMICS
!       CANONICATL MOLECULAR DYNAMICS WITH AND WITHOUT SHAKE
!
!     NOTE
!       REFERENCE
!         (TOTAL ALGORITHM)
!           J.CHEM.PHYS. , VOL.81 , NO.8(1984) , 3684-3690
!         (SHAKE ALGORITHM)
!           J.COMP.PHYS. , VOL.23(1977) , 327-341
!         (VERLET ALGORITHM)
!           PHYS.REV. , VOL.159(1967) , 98
!         (T-V DISTRIBUTION ALGORITHM)
!           J.CHEM.PHYS. , VOL.81 , NO.8(1984) , 3684-3690
!
!********************************************************************
 
      use COMBAS ; use COMERG ; use COMPSC ; use COMMIS
      use COMMOL ; use PHYCNS ; use COMCMMC ; use COMCMM
      use CALC_TIME
      !$ use omp_lib

      implicit none

      ! Logical unit number for output log
        integer(4),intent(in):: iprint
      ! Condition code (0: NO ERROR, negative: ERROR)
        integer(4),intent(out):: ier
      ! Conversion coefficient from force to velocity
        real(8),intent(in):: convcv
      ! Time step (fsec)
        real(8),intent(in):: deltt
        real(8),intent(in):: tmpsta
      ! Number of freedom
        integer(4),intent(in):: numfre
      ! Current velocity (Angstrom/fsec) 
      ! (in: at T-(1/2)dT, out: at T+(1/2)dT)
        real(8),intent(inout):: velback(3,ixnatm)
        real(8),intent(inout):: lmd_vback
      ! Potential energy at time T (Kcal/mol)
        real(8),intent(inout):: epot(maxene),epotcc(maxene,nfrg)
      ! Flag of stop center of mass
        integer(4),intent(in):: idfstp
      ! Flag of velocity scaling for constant temperature simulation
        integer(4),intent(in):: idfvsl
      ! Temperature for setting (K)
        real(8),intent(in):: temps
      ! Sequential loop number
        integer(4),intent(in):: idloop,idhtlp
        
      ! Coordinate at time T for SHAKE
        real(8):: cordm1(3,ixnatm)
      integer(4):: iatm  ! Atom number
      integer(4):: ivar  ! Free atom number
      real(8):: tempc    ! Temperature at T-(1/2)dt (K)
      real(8):: ek,chi,tempshu,tempf(nfrmx),chif(nfrmx),velt(3,ixnatm),&
                cordback(3,ixnatm),forcec(3,ixnatm),lambda_vt
      integer(4):: i,j,mm,numite,icvara(maxatm),ishk,jshk

      real(8):: dprobl,dew,Tdeltt,lmd_back,lambda_f
      real(8):: rtmp,rtmp2,i_idhtlp,i_tempc,i_tempf(nfrmx)

      integer(4):: jjj(1)
      real(8),parameter:: erth = 1.0d-9
      real(8),allocatable:: tg(:,:),Temp_epot(:),Tepot(:)

!*********************************************************************

      ! <<< Initialization >>>
      ier = 0 ; i_idhtlp = 1.d0 / dble(idhtlp) ; ishk = 0 ; jshk = 1
      Tdeltt = deltt ; lambda_f = 0.d0

      do
        ! <<<  BACK-UP CURRENT COORDIANTE FOR SHAKE  >>>
        !$OMP parallel default (none)                                & !
        !$OMP private(i)                                             & !
        !$OMP shared(ixnatm,iynvar,cordm1,cord,velback,vel)
        !$OMP do schedule (static)
        do i = 1,ixnatm
          cordm1(:,i) = cord(:,i)
        enddo
        !$OMP end do
        !$OMP do schedule (static)
        do i = 1,iynvar
          velback(:,i) = vel(:,i)
        enddo
        !$OMP end do
        !$OMP end parallel
        lmd_back = lambda ; lmd_vback = lambda_v

        ! <<<  CALCULATION OF NO-CONSTRAINED VELOCITY AT T+(1/2)DT  >>>
        !      EPOT IS POTENTIAL ENERGY AT TIME T
        call system_clock(dnergy_tim1)
9999    call dnergy(epot,epotcc)
        call system_clock(dnergy_tim2)
        dnergy_tim = dnergy_tim + dnergy_tim2 - dnergy_tim1
        if ( dnergy_tim2 .le. dnergy_tim1 ) dnergy_tim = dnergy_tim + cm

        ! Calc. weight
        select case ( idfvsl )

        !! Canonical MD
        case (3)

          ! Lambda
          if ( cluster_method .eq. "LMD" ) then
            rtmp = lambda*lambda
            select case (nfrg)
            case (1)
              !$OMP parallel default(none)                           & !
              !$OMP private(i)                                       & !
              !$OMP shared(grad,rtmp,ixnatm)
              !$OMP do schedule (static)
              do i = 1,ixnatm
                grad(:,i,1) = grad(:,i,1)*rtmp
              enddo
              !$OMP end do
              !$OMP end parallel
            case (3)
              !$OMP parallel default(none)                           & !
              !$OMP private(i)                                       & !
              !$OMP shared(grad,rtmp,lambda,ixnatm)
              !$OMP do schedule (static)
              do i = 1,ixnatm
                grad(:,i,3) = grad(:,i,1)*rtmp +                       &
                              grad(:,i,2)*lambda + grad(:,i,3)
              enddo
              !$OMP end do
              !$OMP end parallel
            end select

          ! Lambda square
          elseif ( cluster_method .eq. "LMD2" ) then
            select case (nfrg)
            case (1)
              !$OMP parallel default(none)                           & !
              !$OMP private(i)                                       & !
              !$OMP shared(grad,lambda,ixnatm)
              !$OMP do schedule (static)
              do i = 1,ixnatm
                grad(:,i,1) = grad(:,i,1)*lambda
              enddo
              !$OMP end do
              !$OMP end parallel
            case (3)
              rtmp = sqrt(lambda)
              !$OMP parallel default(none)                           & !
              !$OMP private(i)                                       & !
              !$OMP shared(grad,rtmp,lambda,ixnatm)
              !$OMP do schedule (static)
              do i = 1,ixnatm
                grad(:,i,3) = grad(:,i,1)*lambda +                     &
                              grad(:,i,2)*rtmp + grad(:,i,3)
              enddo
              !$OMP end do
              !$OMP end parallel
            end select

          ! Pseudo temperature for canonical MD
          elseif ( psetmp .gt. eps ) then
            dew = temps / psetmp
            select case (nfrg)
            case (1)
              !$OMP parallel default (none)                          & !
              !$OMP private(i)                                       & !
              !$OMP shared(ixnatm,grad,dew)
              !$OMP do schedule (static)
              do i = 1,ixnatm
                grad(:,i,1) = grad(:,i,1) * dew
              enddo
              !$OMP end do
              !$OMP end parallel
            case (2)
              !$OMP parallel default (none)                          & !
              !$OMP private(i)                                       & !
              !$OMP shared(ixnatm,grad,dew)
              !$OMP do schedule (static)
              do i = 1,ixnatm
                grad(:,i,2) = grad(:,i,1)*dew + grad(:,i,2)
              enddo
              !$OMP end do
              !$OMP end parallel
            end select
          endif

        !! McMD or Adaptive lambda dynamics
        case (4)

          ! normal or partial McMD
          if ( cluster_method(1:3) .ne. "LMD" ) then
            if ( epotcc(1,1) .le. celw ) then
              dprobl = alphalw
            elseif ( epotcc(1,1) .ge. ceup ) then
              dprobl = alphaup
            else
              rtmp = (epotcc(1,1)-celw)/(ceup-celw)
              do i = 1,nwindow
                if ( rtmp .lt. high_v(i) ) exit
              enddo
              rtmp = rtmp - low_v(i)
              dprobl = c(0,i)
              do j = 1,ndeg
                dprobl = dprobl + c(j,i)*rtmp**j
              enddo
            endif
            dew = (gascst/joucal)*1.0d-3
            dew = dew * tempce * dprobl
            select case (nfrg)
            case (1)
              !$OMP parallel default (none)                          & !
              !$OMP private(i)                                       & !
              !$OMP shared(ixnatm,grad,dew)
              !$OMP do schedule (static)
              do i = 1,ixnatm
                grad(:,i,1) = grad(:,i,1) * dew
              enddo
              !$OMP end do
              !$OMP end parallel
            case (2)
              !$OMP parallel default (none)                         & !
              !$OMP private(i)                                       & !
              !$OMP shared(ixnatm,grad,dew)
              !$OMP do schedule (static)
              do i = 1,ixnatm
                grad(:,i,2) = grad(:,i,1)*dew + grad(:,i,2)
              enddo
              !$OMP end do
              !$OMP end parallel
            end select

          ! Adaptive lambda dynamics
          else
            !! Calc. ALSD force
            i = nwindow
            do j = 1,nwindow
              if ( lambda .lt. high_v(j) ) then
                i = j ; exit
              endif
            enddo
            rtmp = lambda - low_v(i)
            dprobl = c(0,i)
            do j = 1,ndeg
              dprobl = dprobl + c(j,i)*rtmp**j
            enddo
            if ( lambda .lt. celw ) then
              rtmp2 = celw - lambda
              dprobl = dprobl + lower*rtmp2
            elseif ( lambda .gt. ceup ) then
              rtmp2 = lambda - ceup
              dprobl = dprobl + upper*rtmp2
            endif
            dew = (gascst/joucal)*1.0d-3*tempce * dprobl

            if ( cluster_method .eq. "LMD" ) then
              !! real space
              rtmp = lambda*lambda
              select case (nfrg)
              case (1)
                !$OMP parallel default(none)                         & !
                !$OMP private(i)                                     & !
                !$OMP shared(grad,rtmp,ixnatm)
                !$OMP do schedule (static)
                do i = 1,ixnatm
                  grad(:,i,1) = grad(:,i,1)*rtmp
                enddo
                !$OMP end do
                !$OMP end parallel
              case (3)
                !$OMP parallel default(none)                         & !
                !$OMP private(i)                                     & !
                !$OMP shared(grad,rtmp,lambda,ixnatm)
                !$OMP do schedule (static)
                do i = 1,ixnatm
                  grad(:,i,3) = grad(:,i,1)*rtmp +                     &
                                grad(:,i,2)*lambda + grad(:,i,3)
                enddo
                !$OMP end do
                !$OMP end parallel
              end select
              !! lambda space
              lambda_f = -2.d0*epotcc(1,1)*lambda - epotcc(1,2) - dew

            elseif ( cluster_method .eq. "LMD2" ) then
              !! real space
              select case (nfrg)
              case (1)
                !$OMP parallel default(none)                         & !
                !$OMP private(i)                                     & !
                !$OMP shared(grad,lambda,ixnatm)
                !$OMP do schedule (static)
                do i = 1,ixnatm
                  grad(:,i,1) = grad(:,i,1)*lambda
                enddo
                !$OMP end do
                !$OMP end parallel
              case (3)
                rtmp = sqrt(lambda)
                !$OMP parallel default(none)                         & !
                !$OMP private(i)                                     & !
                !$OMP shared(grad,rtmp,lambda,ixnatm)
                !$OMP do schedule (static)
                do i = 1,ixnatm
                  grad(:,i,3) = grad(:,i,1)*lambda +                   &
                                grad(:,i,2)*rtmp + grad(:,i,3)
                enddo
                !$OMP end do
                !$OMP end parallel
              end select
              !! lambda square space
              lambda_f = -epotcc(1,1) - 0.5d0*epotcc(1,2)/sqrt(lambda) &
                         - dew
            endif
          endif
        end select

        !!! Force check for programing
        select case (fdebug)
          case(2)
            open(unit=111,file="1stp_force.bin",      &
                 form="unformatted",status="old")
            allocate(tg(3,ixnatm),Temp_epot(maxene),Tepot(maxene))
            read(111)Temp_epot(1:maxene)
            read(111)rtmp2
            read(111)tg(1:3,1:ixnatm)
            Tepot(:) = Temp_epot(:)
            Temp_epot(:) = abs(Temp_epot(:) - epot(:))
            rtmp2 = abs(rtmp2-lambda_f)
            tg(1:3,1:ixnatm) = abs( tg(1:3,1:ixnatm) -                 &
                                    grad(1:3,1:ixnatm,nfrg))
            rtmp = maxval(tg(1:3,1:ixnatm))
            print*,""
            print*,"***************************************"
            print*,"*        force error calculation      *"
            print*,"***************************************"
            print*,"  Maximum Force error    = ",rtmp
            print*,"  Force Error of lambda  = ",rtmp2
            jjj = maxloc(Temp_epot(2:)) ; j = jjj(1) + 1
            print*,"  Maximum Energy error   = ",Temp_epot(j)
            print*,"                            (",cyenam(j),")"
            print*,"  Threshold value        = ",erth
            print*,""
            if ( rtmp .gt. erth ) then
              print*,"Atom list with force larger than the Threshold"
              do i = 1,ixnatm
                rtmp = sqrt(dot_product(tg(1:3,i),tg(1:3,i)))
                if ( rtmp .ge. erth ) print*,i,rtmp
              enddo
              print*,"!! Force error was LARGER than the threshold !!"
            else
              print*,"* Force error was less than the threshold *"
            endif
            write(iprint,*)
            j = 0
            write(iprint,'(10x,3(a15,x))'),"Energy","preEnergy","diff."
            do i = 1,maxene
              if ( iyeflg(i) .eq. 1 ) then
                write(iprint,'("  ",a6," :",3(e15.7,x),$)')cyenam(i),  &
                  epot(i),Tepot(i),epot(i)-Tepot(i)
                if ( Temp_epot(i) .gt. erth) then
                  write(iprint,'(a2)')" *"
                else
                  write(iprint,*)
                endif
              endif
            enddo
            if ( maxval(Temp_epot(2:)) .gt. erth ) then
              print*,"!! Energy error was LARGER than the threshold !!"
            else
              print*,"* Energy error was less than the threshold *"
            endif
            close(111) ; deallocate(tg,Temp_epot) ; stop
          case (3)
            open(unit=111,file="1stp_force.bin",          &
                 form="unformatted",status="replace")
            write(111)epot
            write(111)lambda_f
            write(111)grad(1:3,1:ixnatm,nfrg)
            close(111) ; stop
        end select

        !!! velosity at 1/2 dt (velt) without shake force
        rtmp = 0.5d0*Tdeltt*convcv
        !$OMP parallel default (none)                                & !
        !$OMP private(ivar,iatm)                                     & !
        !$OMP shared(iynvar,iylvar,velt,vel,rtmp,grad,nfrg,ifxmass)
        !$OMP do schedule (static)
        do ivar = 1,iynvar
          iatm = iylvar(ivar)
          velt(1:3,ivar) = vel(1:3,ivar) -                             &
                           rtmp*grad(1:3,iatm,nfrg)*ifxmass(iatm)
        enddo
        !$OMP end do
        !$OMP end parallel
        if ( i_lambda_m .ne. 0.d0 )                                    &
          lambda_vt = lambda_v + rtmp*lambda_f*i_lambda_m

        !!! temperature control without shake force
        if ( idfvsl .ne. 2 ) then

          select case (itmeth)
          case (1)
            if ( iyfshk .ne. 0 ) then
              call caltmp(iynvar,iylvar,velt,fxmass(1:ixnatm),numfre,  &
                          tempc,ek,lambda_vt,lambda_m)
            else
              tempc = temps
            endif
            i_tempc = 1.d0 / tempc

            if ( idhtlp.le.0 .or. idloop.gt.idhtlp ) then
              chi = sqrt(temps*i_tempc)
            else
              chi = sqrt((tmpsta+((temps-tmpsta)*i_idhtlp)*idloop)     &
                    *i_tempc)
            endif

            rtmp = chi * Tdeltt * convcv
            rtmp2 = 2.d0 * chi - 1.d0
            !$OMP parallel default(none)                             & !
            !$OMP private(ivar,iatm)                                 & !
            !$OMP shared(iynvar,iylvar,vel,rtmp2,rtmp,grad,nfrg,ifxmass)
            !$OMP do schedule (static)
            do ivar = 1,iynvar
              iatm = iylvar(ivar)
              vel(1:3,ivar) = rtmp2 * vel(1:3,ivar) -                  &
                              rtmp * grad(1:3,iatm,nfrg) * ifxmass(iatm)
            enddo
            !$OMP end do
            !$OMP end parallel
            if ( i_lambda_m .ne. 0.d0 ) then
              if ( l_temp_control ) then
                lambda_v = rtmp2*lambda_v + rtmp*lambda_f*i_lambda_m
              else
                lambda_v = lambda_v + Tdeltt*convcv*lambda_f*i_lambda_m
              endif
            endif
          case (2)
            if ( iyfshk .ne. 0 ) then
              call caltmp2(iynvar,iylvar,velt,fxmass(1:ixnatm),nfrmx,  &
                           tempf,nfrmol,lambda_vt,lambda_m)
            else
              tempf(:) = temps
            endif
            i_tempf(1:nfrmx) = 1.d0 / tempf(1:nfrmx)

            if ( idhtlp.le.0 .or. idloop.gt.idhtlp ) then
              chif(:)= sqrt(temps*i_tempf(:))
            else
              chif(:) = sqrt((tmpsta + ( (temps-tmpsta) *         &
                                i_idhtlp) * idloop ) * i_tempf(:))
            endif

            !$OMP parallel default(none)                             & !
            !$OMP private(mm,rtmp,rtmp2,i,ivar)                      & !
            !$OMP shared(nfrmx,chif,Tdeltt,convcv,natfr,iatfr,vel,   & !
            !$OMP        grad,nfrg,ifxmass)
            do mm = 1,nfrmx
              rtmp = chif(mm) * Tdeltt * convcv
              rtmp2 = 2.d0 * chif(mm) - 1.d0
              !$OMP do schedule(static)
              do i = 1,natfr(mm)
                ivar = iatfr(i,mm)
                vel(1:3,ivar) = rtmp2 * vel(1:3,ivar) -                &
                              rtmp*grad(1:3,ivar,nfrg)*ifxmass(ivar)
              enddo
              !$OMP end do nowait
            enddo
            !$OMP end parallel
            if ( l_temp_control ) then
              rtmp = chif(ifg_lambda) * Tdeltt * convcv
              rtmp2 = 2.d0 * chif(ifg_lambda) - 1.d0
              lambda_v = rtmp2*lambda_v + rtmp*lambda_f*i_lambda_m
            else
              lambda_v = lambda_v + Tdeltt*convcv*lambda_f*i_lambda_m
            endif
          end select 

        endif

        ! <<<  STOP CENTER OF MASS  >>>
        select case (idfstp)
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

        ! <<<  CALCULATION OF NO-CONSTRAINED COORDINATE AT T+DT  >>>
        !$OMP parallel default (none)                                & !
        !$OMP private(ivar,iatm)                                     & !
        !$OMP shared(iynvar,iylvar,cord,Tdeltt,vel)
        !$OMP do schedule (static)
        do ivar = 1,iynvar
          iatm = iylvar(ivar)
          cord(1:3,iatm) = cord(1:3,iatm) + Tdeltt*vel(1:3,ivar)
        enddo
        !$OMP end do
        !$OMP end parallel
        lambda = lambda + Tdeltt*lambda_v

!       <<<  CALCULATION OF CONSTRAINED COORDINATE AT T+DT  >>>
!       <<<  CALCULATION OF CONSTRAINED VELOCITY AT T+(1/2)DT  >>>
        if ( iyfshk .ne. 0 ) then
          !$OMP parallel default (none)                              & !
          !$OMP private(i)                                           & !
          !$OMP shared(ixnatm,forcec)
          !$OMP do schedule (static)
          do i = 1,ixnatm
            forcec(:,i) = 0.d0
          enddo
          !$OMP end do
          !$OMP end parallel
          lambda_vt = ( lmd_vback + lambda_v ) * 0.5d0
          do numite = 1,icslop
            !$OMP parallel default (none)                            & !
            !$OMP private(i)                                         & !
            !$OMP shared(ixnatm,cordback,cord)
            !$OMP do schedule (static)
            do i = 1,ixnatm
              cordback(1:3,i) = cord(1:3,i)
            enddo
            !$OMP end do
            !$OMP end parallel

            call system_clock(shake_tim1)
            select case(ixfbou)
            case (1)
              call shake_PB(cord,cordm1,iprint,ier)
            case default
              call shake_CB(cord,cordm1,iprint,ier)
            end select
            call system_clock(shake_tim2)
            shake_tim = shake_tim + shake_tim2 - shake_tim1
            if ( shake_tim2 .le. shake_tim1 ) shake_tim = shake_tim + cm
            if ( ier .ne. 0 ) then
              !$OMP parallel default (none)                          & !
              !$OMP private(i)                                       & !
              !$OMP shared(ixnatm,cord,cordm1,vel,velback)
              !$OMP do schedule (static)
              do i = 1,ixnatm
                cord(:,i) = cordm1(:,i)
                vel(:,i) = velback(:,i)
              enddo
              !$OMP end do
              !$OMP end parallel
              lambda = lmd_back ; lambda_v = lmd_vback
              if ( Tdeltt .lt. 0.01d0*deltt ) then
                write(iprint,*)"!! SHAKE ERROR  !!"
                write(iprint,*)"!! PROGRAM STOP !!"
                ier = -9999 ; return
              endif
              write(iprint,*)" SHAKE provision (",idloop," step)"
              write(iprint,*)"        deltT = ",Tdeltt," (fs)"
              ier = 0 ; jshk = jshk * 2 ; ishk = ishk * 2
              Tdeltt = Tdeltt * 0.5d0
              goto 9999
            endif

            rtmp = 1.d0 / Tdeltt
            !$OMP parallel default (none)                            & !
            !$OMP private(ivar,iatm)                                 & !
            !$OMP shared(iynvar,iylvar,vel,cord,cordm1,rtmp,velt,velback)
            !$OMP do schedule (static)
            do ivar = 1,iynvar
              iatm = iylvar(ivar)
              vel(1:3,ivar) = (cord(1:3,iatm)-cordm1(1:3,iatm)) * rtmp
              velt(1:3,ivar) = (velback(1:3,ivar)+vel(1:3,ivar))*0.5d0
            enddo
            !$OMP end do
            !$OMP end parallel

            if ( idfvsl .ne. 2 ) then
              select case (itmeth)
              case (1)
                call caltmp(iynvar,iylvar,velt,fxmass(1:ixnatm),numfre,&
                            tempc,ek,lambda_vt,lambda_m)
              case (2)
                call caltmp2(iynvar,iylvar,velt,fxmass(1:ixnatm),nfrmx,&
                             tempf,nfrmol,lambda_vt,lambda_m)
              end select
            endif

            if ( idhtlp.le.0 .or. idloop.gt.idhtlp ) then
              tempshu = temps
            else
              tempshu = tmpsta+((temps-tmpsta)*i_idhtlp)*idloop
            endif

            select case (itmeth)
            case (1)
              if ( abs(tempc-tempshu) .lt. fcstol*tempshu ) exit
            case (2)
              if ( max(abs(minval(tempf)-tempshu),                     &
                       abs(maxval(tempf)-tempshu)) .lt.                &
                   fcstol*tempshu ) exit
            end select
            if ( numite .eq. icslop ) exit

            rtmp = 1.d0 / (Tdeltt*Tdeltt) ; rtmp2 = 0.5d0*Tdeltt
            !$OMP parallel default (none)                            & !
            !$OMP private(ivar,iatm)                                 & !
            !$OMP shared(iynvar,iylvar,forcec,rtmp,cord,cordback,    & !
            !$OMP        velt,rtmp2,convcv,grad,nfrg,ifxmass,velback)
            !$OMP do schedule (static)
            do ivar = 1,iynvar
              iatm = iylvar(ivar)
              forcec(1:3,iatm) = forcec(1:3,iatm) + rtmp *             &
                                 (cord(1:3,iatm)-cordback(1:3,iatm))
              velt(1:3,ivar) = rtmp2*(-convcv*grad(1:3,iatm,nfrg)*     &
                               ifxmass(iatm) + forcec(1:3,iatm)) +     &
                               velback(1:3,ivar)
            enddo
            !$OMP end do
            !$OMP end parallel

            if ( idfvsl .ne. 2 ) then

              select case (itmeth)
              case (1)
                call caltmp(iynvar,iylvar,velt,fxmass(1:ixnatm),numfre,&
                            tempc,ek,lambda_vt,lambda_m)
                i_tempc = 1.d0 / tempc

                if ( idhtlp.le.0 .or. idloop.gt.idhtlp ) then
                  chi = sqrt(temps*i_tempc)
                else
                  chi = sqrt((tmpsta+((temps-tmpsta)*i_idhtlp)*idloop) &
                               *i_tempc)
                endif

                rtmp = chi * Tdeltt ; rtmp2 = 2.d0 * chi - 1.d0
                !$OMP parallel default(none)                         & !
                !$OMP private(ivar,iatm)                             & !
                !$OMP shared(iynvar,iylvar,vel,convcv,grad,nfrg,     & !
                !$OMP        ifxmass,forcec,velback,rtmp,rtmp2)
                !$OMP do
                do ivar = 1,iynvar
                  iatm = iylvar(ivar)
                  vel(1:3,ivar) = (-convcv*grad(1:3,iatm,nfrg)*        &
                    ifxmass(iatm) + forcec(1:3,iatm))*rtmp +           &
                    velback(1:3,ivar)*rtmp2
                enddo
                !$OMP end do
                !$OMP end parallel
!                if ( i_lambda_m .ne. 0.d0 )         &
!                  lambda_v = rtmp*convcv*lambda_f * i_lambda_m +       &
!                             rtmp2*lmd_vback

              case (2)
                call caltmp2(iynvar,iylvar,velt,fxmass(1:ixnatm),nfrmx,&
                             tempf,nfrmol,lambda_vt,lambda_m)
                i_tempf(1:nfrmx) = 1.d0 / tempf(1:nfrmx)

                if ( idhtlp.le.0 .or. idloop.gt.idhtlp ) then
                  chif(:) = sqrt(temps*i_tempf(:))
                else
                  chif(:) =sqrt((tmpsta+((temps-tmpsta)*i_idhtlp) &
                                 *idloop) * i_tempf(:))
                endif
                
                !$OMP parallel default(none)                         & !
                !$OMP private(mm,rtmp,rtmp2,i,ivar)                  & !
                !$OMP shared(nfrmx,chif,Tdeltt,natfr,iatfr,vel,nfrg, & !
                !$OMP        velback,convcv,grad,ifxmass,forcec)
                do mm = 1,nfrmx
                  rtmp = chif(mm) * Tdeltt
                  rtmp2 = 2.d0 * chif(mm) - 1.d0
                  !$OMP do schedule (static)
                  do i = 1,natfr(mm)
                    ivar = iatfr(i,mm)
                    vel(1:3,ivar) = rtmp2*velback(1:3,ivar) + rtmp *   &
                                    ( -convcv * grad(1:3,ivar,nfrg) *  &
                                      ifxmass(ivar) + forcec(1:3,ivar))
                  enddo
                  !$OMP end do nowait
                enddo
                !$OMP end parallel
!                if ( ifg_lambda .ne. 0 ) then
!                  rtmp = chif(ifg_lambda) * Tdeltt
!                  rtmp2 = 2.d0 * chif(ifg_lambda) - 1.d0
!                  lambda_v = rtmp2*lmd_vback+rtmp*convcv*lambda_f *    &
!                             i_lambda_m
!                endif
              end select

            endif

            ! <<<  STOP CENTER OF MASS  >>>
            select case (idfstp)
            case(1:3)
              call remotv(ixcend(nstpcn),fxmass(1:ixnatm),cordm1,ier,  &
                          idfstp)
              if ( ier .ne. 0 ) then
                write(iprint,*)'ERROR> MD '
                write(iprint,*)'   ERROR OCCURES IN REMOTV '
                if ( ier .eq. -1 ) then
                  write(iprint,*)'   MASS WEIGHT OF FIRST CHAIN IS ZERO'
                elseif ( ier .eq. -2 ) then
                  write(iprint,*)'   MOMENTA OF FIRST CHAIN IS SINGULAR'
                endif
                return
              endif
            end select
 
            ! <<<  CALCULATION OF NO-CONSTRAINED COORDINATE AT T+DT  >>>
            !$OMP parallel default (none)                            & !
            !$OMP private(ivar,iatm)                                 & !
            !$OMP shared(iynvar,iylvar,cord,cordm1,Tdeltt,vel)
            !$OMP do schedule (static)
            do ivar = 1,iynvar
              iatm = iylvar(ivar)
              cord(1:3,iatm) = cordm1(1:3,iatm) + Tdeltt*vel(1:3,ivar)
            enddo
            !$OMP end do
            !$OMP end parallel
            lambda = lmd_back + Tdeltt*lambda_v

          enddo
 
        endif
 
        ! <<<  ELLIPSOIDAL BOUNDARY >>>
        if ( ixfbou .eq. 2 ) then
          call ellpbc(maxatm,ixnchn,ixcend,cordm1,fxmass,velback,      &
                      iynvar,iylvar,icvara,Tdeltt,fxcbou,fxellp)
        endif

        ishk = ishk + 1
        if ( ishk .eq. jshk ) then
          exit
        else
          if ( mod(ishk,2) .eq. 0 ) then
            jshk = jshk / 2 ; ishk = ishk / 2 ; Tdeltt = Tdeltt * 2.d0
          endif
        endif
      enddo

!*********************************

      return
      end subroutine verlethe
