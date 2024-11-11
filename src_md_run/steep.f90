
      subroutine steep(iprint,ier,inloop,inlopl,inupdl,inmons,inflog,  &
                       infbst,fncovg,fnstpi,fnupsl,fndnsl,inuana,      &
                       fnepot)

!*******************************************************************
!
!     STEEPEST DESCENT MINIMIZATION
!
!*******************************************************************
 
      use COMBAS ; use COMERG ; use COMCMMC
      use COMCMM,only: nfrg

      implicit none

      ! inloop: current loop number
      ! inlopl: loop limit
      ! inupdl: update interval
      ! inmons: monitoring inverval
      ! inflog: flag for format of log
      ! infbst: flag for best-fit
      ! inuana: logical unit number for output analysis data
        integer(4),intent(in):: iprint,inloop,inlopl,inupdl,inmons,    &
                                inflog,infbst,inuana
      integer(4),intent(out):: ier
      ! fncovg: conversion criterion (if rmsf <= fncovg, minimization is ended)
      ! fnstpi: initial step length
      ! fnupsl: scale up rate of stemp length
      ! fndnsl: scale down rate of step length
        real(8),intent(in):: fncovg,fnstpi,fnupsl,fndnsl
        real(8),intent(inout):: fnepot(maxene)

      integer(4):: iloop,itmp
      integer(8):: sec
      ! srdir: search direction
        real(8):: stplen,rmsf,rtmp,srdir(3,ixnatm)
 
!**************************************************

      ier = 0 ; iloop = 0 ; stplen = fnstpi*sqrt(dble(iynvar))
      write(iprint,*)"INFORMATION> STEEP (V3.0)"
      write(iprint,*)"             STEEPEST DESCENT MINIMIZATION "
 
!     <<<  MINIMIZATION LOOP  >>>
      do
        ! 1) make search direction
        call calrsf(iynvar,iylvar,grad(1:3,1:ixnatm,nfrg),rmsf,rtmp,   &
                    itmp)
        if ( rmsf .le. fncovg ) then
          write(iprint,*)" "
          write(iprint,*)"INFORMATION> STEEP "
          write(iprint,*)"  MINIMIZATION IS CONVERGED " ; return
        endif
        call mksrdg(iynvar,iylvar,srdir)

        ! 2) one-step line search
        call linsrd(srdir,inloop,fnupsl,fndnsl,stplen,iloop,fnepot,    &
                    iprint,ier)
        if ( ier .ne. 0 ) return

        ! 3) monitoring
        if ( mod(inloop,inmons) .eq. 0 ) then
          call monmin(iprint,inloop,stplen,inflog,infbst,inuana,fnepot)
        endif

        ! 4) update interaction table
        if ( mod(inloop,inupdl).eq.0 ) then
          call celprp(iprint,inflog,ier)
          if ( ier .ne. 0 ) return
        endif 

        ! 5) check cpu time
        call system_clock(sec)
        if ( dble(sec-fxcpus)/dble(cr) .ge. fxcpul ) then
          write(iprint,*)" "
          write(iprint,*)"WARNING> STEEP "
          write(iprint,*)"  CPU TIME IS EXCEEDED THEN STOP MINIMIZATION"
          ier = 1 ; return
        endif
        if ( iloop .eq. inlopl ) then
          write(iprint,*)" "
          write(iprint,*)"INFORMATION> STEEP"
          write(iprint,*)"  MINIMIZATION IS NOT CONVERGED WITHIN LOOP "
          return
        endif
      enddo

!******************************

      return
      end subroutine steep
 
 
!===================================================================


      subroutine linsrd(srdir,inloop,fnupsl,fndnsl,stplen,iloop,       &
                        enenow,iprint,ier)

!*******************************************************************
!
!     ONE STEP LINE SEARCH
!     SEARCH ONLY ONE MINIMIZED STRUCTURE ALONG ONE SEARCH DIRECTION
!
!*******************************************************************
 
      use COMBAS ; use COMERG ; use PHYCNS ; use COMCMM ; use COMCMMC
      !$ use omp_lib

      implicit none    

      integer(4):: inloop,iloop,iprint,ier
      real(8):: srdir(3,ixnatm),fnupsl,fndnsl,stplen,enenow(maxene)

      integer(4):: i,j,ivar
      real(8):: epold,enencc(maxene,nfrg),rtmp
 
!******************************************

      ier = 0 ; epold = enenow(1)
      if ( iyfshk .ne. 0 ) then
        !$OMP parallel default (none)                                & !
        !$OMP private(i)                                             & !
        !$OMP shared(grad,ixnatm,nfrg,cord)
        !$OMP do schedule (static)
        do i = 1,ixnatm
          grad(:,i,nfrg) = cord(:,i)
        enddo
        !$OMP end do
        !$OMP end parallel
      endif

!     <<<  ONE STEP LINE SEARCH  >>>
      !$OMP parallel default (none)                                  & !
      !$OMP private(ivar,i)                                          & !
      !$OMP shared(iynvar,iylvar,cord,stplen,srdir)
      !$OMP do schedule (static)
      do ivar = 1,iynvar
        i = iylvar(ivar)
        cord(1:3,i) = cord(1:3,i) + stplen*srdir(1:3,ivar)
      enddo
      !$OMP end do
      !$OMP end parallel
      if ( iyfshk .ne. 0 ) then
        select case(ixfbou)
        case (1)
          call shake_PB(cord,grad(1:3,1:ixnatm,nfrg),iprint,ier)
        case default
          call shake_CB(cord,grad(1:3,1:ixnatm,nfrg),iprint,ier)
        end select
        if ( ier .gt. 0 ) then
         stplen = stplen * fndnsl
         write(iprint,*)" "
         write(iprint,*)"INFORMATION> STEEP "
         write(iprint,*)"  SHAKE IS NOT CONVERGED "
         write(iprint,*)"  THEN (STEP-LENGTH) = (STEP-LENGTH) * ",fndnsl
         write(iprint,*)" "
         if ( stplen / sqrt(dble(iynvar*3)) .lt. eps ) then
           write(iprint,*)" "
           write(iprint,*)"ERROR> STEEP "
           write(iprint,*)"    STEP LENGTH IS TOO SMALL "
           write(iprint,*)"    AVERAGE STEP LENGTH (ANGSTROMS) : ",    &
                          stplen / sqrt(dble(iynvar*3))
           write(iprint,*)"    LIMIT STEP LENGTH (ANGSTROMS)   : ",eps
           write(iprint,*)" " ; ier = -2
         endif
         return
        elseif ( ier .lt. 0 ) then
          ier = -1 ; return
        endif
      endif

      iloop = iloop + 1 ; inloop = inloop + 1
      call dnergy(enenow,enencc)
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
      endif

      if ( enenow(1) .lt. epold ) then
        stplen = stplen * fnupsl
      else
        stplen = stplen * fndnsl
         if ( stplen / sqrt(dble(iynvar*3)) .lt. eps ) then
           write(iprint,*)" "
           write(iprint,*)"ERROR> STEEP "
           write(iprint,*)"    STEP LENGTH IS TOO SMALL "
           write(iprint,*)"    AVERAGE STEP LENGTH (ANGSTROMS) : ",    &
                          stplen / sqrt(dble(iynvar*3))
           write(iprint,*)"    LIMIT STEP LENGTH (ANGSTROMS)   : ",eps
           write(iprint,*)" " ; ier = -2 ; return
         endif
      endif

!*******************************

      return
      end subroutine linsrd
