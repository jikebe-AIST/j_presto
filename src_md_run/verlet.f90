
      subroutine verlet(iprint,ier,convcv,deltt,numfre,velback,epot,   &
                        epotcc,idfstp,idfvsl,temps,taut,lamda)

!********************************************************************
!
!     ONE-STEP NUMERICAL INTEGRATION OF MOLECULAR DYNAMICS
!     LEAP-FROG VERLET
!     AVAILABLE MOLECULAR DYNAMICS
!       MICROCANONICAL MOLECULAR DYNAMICS
!       CONSTANT TEMPERATURE AND CONSTANT VOLUME (SUBROUTINE SCALET)
!       CONSTRAINED MOLECULAR DYNAMICS (SUBROUTINE SHAKE)
!          WITH THE ELASTIC COLLISION AT THE ELLIPSOIDAL SURFACE
!
!********************************************************************
 
      use COMBAS ; use COMERG ; use COMCMMC ; use COMCMM ; use CALC_TIME
      !$ use omp_lib

      implicit none

      ! idfstp : flag of stop center of mass
      ! idfvsl : flag of velocity scaling for constant Temp. simulation
        integer(4),intent(in):: iprint,numfre,idfstp,idfvsl
      integer(4),intent(out):: ier
      ! convcv : conversion coefficient from force to velocity (4.184*1.d-4)
      ! taut : relaxation time for constant Temp. simulation (fsec)
        real(8),intent(in):: convcv,deltt,temps,taut
      ! lamda : scaling factor of verocity ( idfvsl .eq. 1)
        real(8),intent(out):: epot(maxene),lamda,velback(3,ixnatm)
      
      ! coordinate at time T for SHAKE
        real(8):: cordm1(3,ixnatm)
      real(8):: tempc,ek,rtmp,epotcc(maxene,nfrg)
      integer(4):: i,j,icvara(maxatm)

!*****************************************

      ier = 0
!     <<<  CALCULATION OF TEMERATURE AT T-(1/2)DT  >>>
      if ( idfvsl .eq. 1 ) then
        call caltmp(iynvar,iylvar,vel,fxmass(1:ixnatm),numfre,tempc,ek,&
                    lambda_v,lambda_m)
      endif
      !$OMP parallel default (none)                                  & !
      !$OMP private(i)                                               & !
      !$OMP shared(cordm1,ixnatm,cord,velback,vel,iynvar)
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

!     <<<  CALCULATION OF NO-CONSTRAINED VELOCITY AT T+(1/2)DT  >>>
!          EPOT IS POTENTIAL ENERGY AT TIME T
      call system_clock(dnergy_tim1)
      call dnergy(epot,epotcc)
      call system_clock(dnergy_tim2)
      dnergy_tim = dnergy_tim + dnergy_tim2 - dnergy_tim1
      if ( dnergy_tim2 .le. dnergy_tim1 ) dnergy_tim = dnergy_tim + cm
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
          !$OMP parallel default (none)                              & !
          !$OMP private(i)                                           & !
          !$OMP shared(ixnatm,grad,rtmp)
          !$OMP do schedule (static)
          do i = 1,ixnatm
            grad(:,i,1) = grad(:,i,1) * rtmp
          enddo
          !$OMP end do
          !$OMP end parallel
        case (2)
          !$OMP parallel default (none)                             & !
          !$OMP private(i)                                           & !
          !$OMP shared(ixnatm,grad,rtmp)
          !$OMP do schedule (static)
          do i = 1,ixnatm
            grad(:,i,2) = grad(:,i,1)*rtmp + grad(:,i,2)
          enddo
          !$OMP end do
          !$OMP end parallel
        end select
      endif

      rtmp = -deltt*convcv
      !$OMP parallel default (none)                                  & !
      !$OMP private(i,j)                                             & !
      !$OMP shared(iynvar,iylvar,vel,rtmp,grad,nfrg,fxmass)
      !$OMP do schedule (static)
      do i = 1,iynvar
        j = iylvar(i)
        vel(1:3,i) = vel(1:3,i) + rtmp*grad(1:3,j,nfrg)/fxmass(j)
      enddo
      !$OMP end do
      !$OMP end parallel
 
!     <<<  STOP CENTER OF MASS  >>>
      if ( idfstp.eq.1 .or. idfstp.eq.2 .or. idfstp.eq.3 ) then
        call remotv(ixcend(nstpcn),fxmass(1:ixnatm),cordm1,ier,idfstp)
        if ( ier .ne. 0 ) then
          write(iprint,*)"ERROR> MD "
          write(iprint,*)"   ERROR OCCURES IN REMOTV "
          if ( ier .eq. -1 ) then
            write(iprint,*)"   MASS WEIGHT OF FIRST CHAIN IS ZERO "
          elseif ( ier .eq. -2 ) then
            write(iprint,*) "   MOMENTA OF FIRST CHAIN IS SINGULAR "
          endif
          return
        endif
      endif

!     <<<  SCALING VELOCITY FOR CONSTANT TEMPERATURE  >>>
      if ( idfvsl .eq. 1 ) then
        call scalet(temps,tempc,taut,deltt,iprint,iynvar,lamda)
      else
        lamda = 0.0D0
      endif
 
!     <<<  CALCULATION OF NO-CONSTRAINED COORDINATE AT T+DT  >>>
      !$OMP parallel default (none)                                  & !
      !$OMP private(i,j)                                             & !
      !$OMP shared(iynvar,iylvar,cord,deltt,vel)
      !$OMP do schedule (static)
      do i = 1,iynvar
        j = iylvar(i)
        cord(1:3,j) = cord(1:3,j) + deltt*vel(1:3,i)
      enddo
      !$OMP end do
      !$OMP end parallel
 
!     <<<  CALCULATION OF CONSTRAINED COORDINATE AT T+DT  >>>
!     <<<  CALCULATION OF CONSTRAINED VELOCITY AT T+(1/2)DT  >>>
      if ( iyfshk .ne. 0 ) then
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
        if ( ier .ne. 0 ) return
        rtmp = 1.d0 / deltt
        !$OMP parallel default (none)                                & !
        !$OMP private(i,j)                                           & !
        !$OMP shared(iynvar,iylvar,vel,cord,cordm1,rtmp)
        !$OMP do schedule (static)
        do i = 1,iynvar
          j = iylvar(i)
          vel(1:3,i) = (cord(1:3,j)-cordm1(1:3,j)) * rtmp
        enddo
        !$OMP end do
        !$OMP end parallel
      endif

!     <<<  ELLIPSOIDAL BOUNDARY >>>
      if ( ixfbou .eq. 2 ) then
        call ellpbc(maxatm,ixnchn,ixcend,cordm1,fxmass,velback,vel,    &
                    iynvar,iylvar,icvara,deltt,fxcbou,fxellp)
      endif

!***********************************************

      return
      end subroutine verlet
