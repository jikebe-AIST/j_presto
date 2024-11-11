
      subroutine conj(iprint,ier,inloop,inlopl,inupdl,inmons,inflog,   &
                      infbst,fncovg,fnstpi,fncovl,inllpl,inuana,fnepot)

!*******************************************************************
!
!     CONJUGATE GRADIENT MINIMIZATION
!
!*******************************************************************

      use COMBAS ; use COMERG 
      !$ use omp_lib

      implicit none

      ! Logical unit number for log
        integer(4),intent(in):: iprint
      ! Condition code (0: No error, negative: Error, 
      !                 positive: CPU time is exceeded)
        integer(4),intent(out):: ier
      ! Current loop number
        integer(4),intent(inout):: inloop
      ! Loop limit
        integer(4),intent(in):: inlopl
      ! Update interval
        integer(4),intent(in):: inupdl
      ! Monitoring interval
        integer(4),intent(in):: inmons
      ! Flag of log format
        integer(4),intent(in):: inflog
      ! Flag for best-fit (infbst=1, calc. RMSD of first chain)
        integer(4),intent(in):: infbst
      ! Conversion criterion (if RMSF <= fncovg, minimization is ended)
        real(8),intent(in):: fncovg
      ! Initial step length
        real(8),intent(in):: fnstpi
      ! Conversion criterion of line search
        real(8),intent(in):: fncovl
      ! Loop limit of each line search
        integer(4),intent(in):: inllpl
      ! Logical unit number for output analysis data
        integer(4),intent(in):: inuana
      ! Current potential energy
        real(8),intent(inout):: fnepot(maxene)

      ! Search direction
        real(8):: srdir(3,ixnatm)
      ! Difference vector between current and initial gradient at restart
        real(8):: difgrd(3,ixnatm)
      ! Previous search direction & restart search direction
        real(8):: srdiro(3,ixnatm),srdirr(3,ixnatm)
      ! Step length for line search (angstroms)
        real(8):: stplen
      ! srdirr * difgrd
        real(8):: dtdifg
      ! Current directional derivative (srdir * fngrad)
        real(8):: dirgrd
      ! Flag of search direction (1: Steepest descent direction
      !                           2: Conjugate gradient direction (Normal) 
      !                           3: Conjugate gradient direction (Restart))
        integer(4):: ifsrdr
      ! Search direction number
        integer(4):: insrdr
      ! Restart search direction number
        integer(4):: inrsdr

      real(8):: grado(3,ixnatm)
      logical(4):: onconv,onlopl
      integer(4):: iloop,i
 
!********************************************************************
 
!     <<<  INITILIZATION  >>>
      ier = 0 ; stplen = fnstpi*sqrt(dble(iynvar))
      ifsrdr = 1 ; insrdr = 0 ; inrsdr = 0 ; iloop = 0
      !$OMP parallel default (none)                                  & !
      !$OMP private(i)                                               & !
      !$OMP shared(grado,ixnatm)
      !$OMP do schedule (static)
      do i = 1,ixnatm
        grado(:,i) = 0.d0
      enddo
      !$OMP end do
      !$OMP end parallel
      WRITE(IPRINT,*) 'INFORMATION> CONJ (V1.0)'
      WRITE(IPRINT,*) '       CONJUGATE GRADIENT MINIMIZATION '
 
      if ( iyfshk .ne. 0 ) then
        write(iprint,*)" "
        write(iprint,*)"ERROR> CONJ "
        write(iprint,*)"    SHAKE IS NOT AVAILABLE IN CONJUGATE "//    &
                       "GRADIENT MINIMIZATION "
        write(iprint,*)" "
        ier = -1 ; return
      endif
 
!     <<<  MINIMIZATION LOOP  >>>
      do
        ! Make search direction
        call mksrdc(iprint,ier,ifsrdr,insrdr,inrsdr,difgrd,dtdifg,     &
                    dirgrd,srdir,srdiro,srdirr,grado)
        if ( ier .ne. 0 ) return

        ! Line search
        call linsrc(iprint,ier,inloop,inlopl,inupdl,inmons,inflog,     &
                    infbst,fncovg,fncovl,inllpl,inuana,fnepot,srdir,   &
                    iloop,stplen,dirgrd,onlopl,onconv)
        if ( onlopl ) return
        if ( onconv ) return
        if ( ier .ne. 0 ) return
      enddo
 
!********************************

      return
      end subroutine conj


!===========================================================================
 

      subroutine mksrdc(iprint,ier,ifsrdr,insrdr,inrsdr,difgrd,dtdifg, &
                        dirgrd,srdir,srdiro,srdirr,grado)

      use COMBAS ; use COMERG
      use COMCMMC,only: grad
      use COMCMM,only: nfrg
      !$ use omp_lib

      implicit none

      integer(4):: iprint,ier,ifsrdr,insrdr,inrsdr
      real(8):: dtdifg,dirgrd,grado(3,ixnatm),difgrd(3,ixnatm),        &
                srdir(3,ixnatm),srdiro(3,ixnatm),srdirr(3,ixnatm)

      logical(4):: offdes
      real(8):: rss1,rss2
      integer(4):: ivar,iatm
 
!******************************************************

!     <<<  MAKE SEARCH DIRECTION  >>>
      do
        !! CHECK ORTHOGONALITY BETWEEN CURRENT GRADINET AND PREVIOUS GRADINET
        if ( ifsrdr .eq. 2 ) then
          rss1 = 0.d0 ; rss2 = 0.d0
          do ivar = 1,iynvar
            iatm = iylvar(ivar)
            rss1 = rss1 + grado(1,iatm)*grad(1,iatm,nfrg)              &
                        + grado(2,iatm)*grad(2,iatm,nfrg)              &
                        + grado(3,iatm)*grad(3,iatm,nfrg)
            rss2 = rss2 + grad(1,iatm,nfrg)*grad(1,iatm,nfrg)          &
                        + grad(2,iatm,nfrg)*grad(2,iatm,nfrg)          &
                        + grad(3,iatm,nfrg)*grad(3,iatm,nfrg)
          enddo
          if ( abs(rss1)/rss2 .ge. 0.2d0 ) ifsrdr = 3
        endif

        select case (ifsrdr)
          !! STEEPEST DESCENT DIRECTION
          case (1)
            call mksrdg(iynvar,iylvar,srdir)
            !$OMP parallel default (none)                            & !
            !$OMP private(ivar)                                      & !
            !$OMP shared(srdiro,iynvar,srdir)
            !$OMP do schedule (static)
            do ivar = 1,iynvar
              srdiro(:,ivar) = srdir(:,ivar)
            enddo
            !$OMP end do
            !$OMP end parallel

          !!CONJUGATE GRADINET DIRECTION (NORMAL)
          case (2)
            call mksdc1(iprint,ier,iynvar,iylvar,grado,difgrd,dtdifg,  &
                        srdiro,srdirr,srdir)
            if ( ier .ne. 0 ) return

          !! CONJUGATE GRADINET DIRECTION (RESTART)
          case (3)
            call mksdc2(iprint,ier,iynvar,iylvar,grado,difgrd,dtdifg,  &
                        srdiro,srdirr,srdir)
            if ( ier .ne. 0 ) return
            inrsdr = insrdr
        end select

!       <<<  CHECK SEARCH DIRECTION  >>>
!          CALCULATE INITIAL DIRECTIONAL DERIVATIVE
        call chksrd(iynvar,iylvar,srdir,dirgrd,offdes)
        if ( .not. offdes ) exit
        ifsrdr = 1
      enddo

!     <<<  PREPAR DATA FOR NEXT SEARCH DIRECTION  >>>
      insrdr = insrdr + 1
      if ( ifsrdr .eq. 1 ) then
        ifsrdr = 3
      elseif ( ifsrdr .eq. 3 ) then
        ifsrdr = 2
      endif
      if ( insrdr-inrsdr .eq. iynvar*3 ) ifsrdr = 3
      !$OMP parallel default (none)                                  & !
      !$OMP private(ivar,iatm)                                       & !
      !$OMP shared(iynvar,iylvar,grado,grad,nfrg)
      !$OMP do schedule (static)
      do ivar = 1,iynvar
        iatm = iylvar(ivar)
        grado(1:3,iatm) = grad(1:3,iatm,nfrg)
      enddo
      !$OMP end do
      !$OMP end parallel

!************************************************

      return
      end subroutine mksrdc


!========================================================================        


      subroutine mksdc1(iprint,ier,numvar,listv,grado,difgrd,dtdifg,   &
                        srdiro,srdirr,srdir)

!*******************************************************************
!
!      MAKE SEARCH DIRECTION OF CONJUGATE GRADINET
!        USE SORENSON-WOLFE COEFFICIENT
!        USE BEALE'S RESTART PROCEDURE
!
!*******************************************************************
 
      use COMBAS,only: ixnatm
      use COMCMMC,only: grad ; use COMCMM,only: nfrg
      !$ use omp_lib

      implicit none

      integer(4):: iprint,ier,numvar,listv(ixnatm)
      real(8):: dtdifg,grado(3,ixnatm),difgrd(3,ixnatm),               &
                srdiro(3,ixnatm),srdirr(3,ixnatm),srdir(3,ixnatm)
      real(8):: rss,beta,ganma,difx,dify,difz
      integer(4):: ivar,iatm
 
!*********************************************************
 
!     <<<  CALCULATE SORENSON-WOLFE COEFFICIENT  >>>
      rss = 0.d0 ; beta = 0.d0 ; ganma = 0.d0
      !$OMP parallel default (none)                                  & !
      !$OMP private(ivar,iatm,difx,dify,difz)                        & !
      !$OMP shared(srdiro,srdir,numvar,listv,grad,grado,nfrg)        & !
      !$OMP reduction(+ : beta,rss,ganma,difgrd)
      !$OMP do schedule (static)
      do ivar = 1,numvar
        srdiro(:,ivar) = srdir(:,ivar)
      enddo
      !$OMP end do
      !$OMP do schedule (static)
      do ivar = 1,numvar
        iatm = listv(ivar)
        difx = grad(1,iatm,nfrg) - grado(1,iatm)
        dify = grad(2,iatm,nfrg) - grado(2,iatm)
        difz = grad(3,iatm,nfrg) - grado(3,iatm)
        beta = beta + grad(1,iatm,nfrg)*difx + grad(2,iatm,nfrg)*dify  &
                    + grad(3,iatm,nfrg)*difz
        rss = rss + srdiro(1,ivar)*difx + srdiro(2,ivar)*dify +        &
                    srdiro(3,ivar)*difz
        ganma = ganma + grad(1,iatm,nfrg)*difgrd(1,iatm) +             &
                        grad(2,iatm,nfrg)*difgrd(2,iatm) +             &
                        grad(3,iatm,nfrg)*difgrd(3,iatm)
      enddo
      !$OMP end do
      !$OMP end parallel
        
      if ( rss .eq. 0.d0 ) then
        write(iprint,*)" "
        write(iprint,*)"ERROR> CONJ "
        write(iprint,*)"   STEP LENGTH MAY BE TOO SMALL "
        write(iprint,*)" " ; ier = -1
        return
      endif
      beta = beta / rss
      ganma = ganma / dtdifg

!     <<<  CALCULATE NEW SEARCH DIRECTION  >>>
      rss = 0.d0
      !$OMP parallel default (none)                                  & !
      !$OMP private(ivar,iatm)                                       & !
      !$OMP shared(numvar,listv,srdir,beta,srdiro,ganma,srdirr,grad, & !
      !$OMP        nfrg)                                             & !
      !$OMP reduction (+ : rss)
      !$OMP do schedule (static)
      do ivar = 1,numvar
        iatm = listv(ivar)
        srdir(1,ivar) = beta*srdiro(1,ivar) + ganma*srdirr(1,ivar) -   &
                        grad(1,iatm,nfrg)
        srdir(2,ivar) = beta*srdiro(2,ivar) + ganma*srdirr(2,ivar) -   &
                        grad(2,iatm,nfrg)
        srdir(3,ivar) = beta*srdiro(3,ivar) + ganma*srdirr(3,ivar) -   &
                        grad(3,iatm,nfrg)
        rss = rss + srdir(1,ivar)*srdir(1,ivar)                        &
                  + srdir(2,ivar)*srdir(2,ivar)                        &
                  + srdir(3,ivar)*srdir(3,ivar)
      enddo
      !$OMP end do
      !$OMP end parallel
      rss = 1.d0 / sqrt(rss)
      !$OMP parallel default (none)                                  & !
      !$OMP private(ivar)                                            & !
      !$OMP shared(srdir,rss,numvar)
      !$OMP do schedule (static)
      do ivar = 1,numvar
        srdir(:,ivar) = rss*srdir(:,ivar)
      enddo
      !$OMP end do
      !$OMP end parallel
 
!******************************************

      return
      end subroutine mksdc1
 
 
!========================================================================


      subroutine mksdc2(iprint,ier,numvar,listv,grado,difgrd,dtdifg,   &
                        srdiro,srdirr,srdir)

!*******************************************************************
!
!     MAKE SEARCH DIRECTION OF CONJUGATE GRADINET
!       USE SORENSON-WOLFE COEFFICIENT
!
!*******************************************************************
 
      use COMBAS,only: ixnatm
      use COMCMMC,only: grad ; use COMCMM,only: nfrg
      !$ use omp_lib

      implicit none

      integer(4):: iprint,ier,numvar,listv(ixnatm)
      real(8):: dtdifg,grado(3,ixnatm),difgrd(3,ixnatm),               &
                srdiro(3,ixnatm),srdirr(3,ixnatm),srdir(3,ixnatm)
      real(8):: rss,beta
      integer(4):: ivar,iatm 
 
!*******************************************
 
!     <<<  CALCULATE SORENSON-WOLFE COEFFICIENT  >>>
      rss = 0.d0 ; beta = 0.d0 ; dtdifg = 0.d0
      !$OMP parallel default (none)                                  & !
      !$OMP private(ivar,iatm)                                       & !
      !$OMP shared(srdiro,srdir,srdirr,numvar,listv,difgrd,grad,nfrg,& !
      !$OMP        grado)                                            & !
      !$OMP reduction (+ : beta,rss)
      !$OMP do schedule (static)
      do ivar = 1,numvar
        srdiro(:,ivar) = srdir(:,ivar)
        srdirr(:,ivar) = srdir(:,ivar)
        iatm = listv(ivar)
        difgrd(1:3,iatm) = grad(1:3,iatm,nfrg) - grado(1:3,iatm)
        beta = beta + grad(1,iatm,nfrg)*difgrd(1,iatm) +               &
                      grad(2,iatm,nfrg)*difgrd(2,iatm) +               &
                      grad(3,iatm,nfrg)*difgrd(3,iatm)
        rss = rss + srdiro(1,ivar)*difgrd(1,iatm) +                    &
                    srdiro(2,ivar)*difgrd(2,iatm) +                    &
                    srdiro(3,ivar)*difgrd(3,iatm)
      enddo
      !$OMP end do
      !$OMP end parallel
      dtdifg = rss

      if ( rss .eq. 0.d0 ) then
        write(iprint,*)" "
        write(iprint,*)"ERROR> CONJ "
        write(iprint,*)"   STEP LENGTH MAY BE TOO SMALL "
        write(iprint,*)" " ; ier = -1
        return
      endif
      beta = beta / rss
      
!     <<<  CALCULATE NEW SEARCH DIRECTION  >>>
      rss = 0.d0
      !$OMP parallel default (none)                                  & !
      !$OMP private(ivar,iatm)                                       & !
      !$OMP shared(numvar,listv,srdir,beta,srdiro,grad,nfrg)         & !
      !$OMP reduction (+ : rss)
      !$OMP do schedule (static)
      do ivar = 1,numvar
        iatm = listv(ivar)
        srdir(1:3,ivar) = beta*srdiro(1:3,ivar) - grad(1:3,iatm,nfrg)
        rss = rss + srdir(1,ivar)*srdir(1,ivar) +                      &
                    srdir(2,ivar)*srdir(2,ivar) +                      &
                    srdir(3,ivar)*srdir(3,ivar)
      enddo
      !$OMP end do
      !$OMP end parallel
      rss = 1.d0 / sqrt(rss)
      !$OMP parallel default (none)                                  & !
      !$OMP private(ivar)                                            & !
      !$OMP shared(srdir,numvar,rss)
      !$OMP do schedule (static)
      do ivar = 1,numvar
        srdir(:,ivar) = rss*srdir(:,ivar)
      enddo
      !$OMP end do
      !$OMP end parallel

!****************************************

      return
      end subroutine mksdc2


!================================================================================== 
 

      subroutine chksrd(iynvar,iylvar,srdir,dirgrd,offdes)

      use COMBAS,only: ixnatm
      use COMCMMC,only: grad  ; use COMCMM,only: nfrg
      !$ use omp_lib

      implicit none

      integer(4):: iynvar,iylvar(iynvar)
      real(8):: srdir(3,ixnatm),dirgrd
      logical(4):: offdes
      integer(4):: ivar,iatm
 
!*******************************************************

      dirgrd = 0.d0
      !$OMP parallel default (none)                                  & !
      !$OMP private(iatm,ivar)                                       & !
      !$OMP shared(iynvar,iylvar,grad,nfrg,srdir)                    & !
      !$OMP reduction (+ : dirgrd)
      !$OMP do schedule (static)
      do ivar = 1,iynvar
        iatm = iylvar(ivar)
        dirgrd = dirgrd + grad(1,iatm,nfrg)*srdir(1,ivar) +            &
                          grad(2,iatm,nfrg)*srdir(2,ivar) +            &
                          grad(3,iatm,nfrg)*srdir(3,ivar)
      enddo
      !$OMP end do
      !$OMP end parallel

      if ( dirgrd .ge. 0.d0 ) then
        offdes = .true.
      else
        offdes = .false.
      endif

!**********************************

      return
      end subroutine chksrd


!================================================================================


      subroutine linsrc(iprint,ier,inloop,inlopl,inupdl,inmons,inflog, &
                        infbst,fncovg,fncovl,inllpl,inuana,fnepot,     &
                        srdir,iloop,stplen,dirgrs,onlopl,onconv)

!*******************************************************************
!
!     LINE SEARCH
!
!*******************************************************************
 
      use COMBAS ; use COMERG ; use COMCMMC
      use COMCMM
      !$ use omp_lib

      implicit none

      integer(4):: iprint,ier,i,j,inloop,inlopl,inupdl,inmons,inflog,  &
                   infbst,inllpl,inuana,iloop
      real(8):: fncovg,fncovl,rtmp,fnepot(maxene),trycor(3,ixnatm),    &
                fnepcc(maxene,nfrg),srdir(3,ixnatm),stplen,dirgrs
      logical(4):: onlopl,onconv,onfind

      real(8):: dirgrd,dirgrn,dirgrp,stplnn,stplnp,rss,Ttime
      integer(8):: sec
      integer(4):: ivar,iatm,jloop 
 
!******************************************

      onlopl = .false. ; onconv = .false. ; onfind = .false.
      jloop = 0 ; dirgrn = dirgrs ; stplnn = 0.d0
      !$OMP parallel default (none)                                  & !
      !$OMP private(i)                                               & !
      !$OMP shared(trycor,ixnatm,cord)
      !$OMP do schedule (static)
      do i = 1,ixnatm
        trycor(:,i) = cord(:,i)
      enddo
      !$OMP end do
      !$OMP end parallel

      do
        ! 1) UPDATE INTERACTION TABLE
        if ( iloop.ne.0 .and. mod(inloop,inupdl).eq.0 ) then
          call celprp(iprint,inflog,ier)
          if ( ier .ne. 0 ) return
        endif

        ! 2) MAKE TRIAL STRUCTURE
        !$OMP parallel default (none)                                & !
        !$OMP private(iatm,ivar)                                     & !
        !$OMP shared(iynvar,iylvar,cord,trycor,stplen,srdir)
        !$OMP do schedule (static)
        do ivar = 1,iynvar
          iatm = iylvar(ivar)
          cord(1:3,iatm) = trycor(1:3,iatm) + stplen*srdir(1:3,ivar)
        enddo
        !$OMP end do
        !$OMP end parallel

        ! 3) CALCULATE ENERGY AND GRADIENT OF TRIAL STRUCTURE
        inloop = inloop + 1 ; iloop = iloop + 1 ; jloop = jloop + 1
        call dnergy(fnepot,fnepcc)
        if ( cluster_method .eq. "LMD" ) then
          rtmp = lambda*lambda
          select case (nfrg)
          case (1)
            !$OMP parallel default(none)                             & !
            !$OMP private(i)                                         & !
            !$OMP shared(grad,rtmp,ixnatm)
            !$OMP do schedule (static)
            do i = 1,ixnatm
              grad(:,i,1) = grad(:,i,1)*rtmp
            enddo
            !$OMP end do
            !$OMP end parallel
          case (3)
            !$OMP parallel default(none)                             & !
            !$OMP private(i)                                         & !
            !$OMP shared(grad,rtmp,lambda,ixnatm)
            !$OMP do schedule (static)
            do i = 1,ixnatm
              grad(:,i,3) = grad(:,i,1)*rtmp +                         &
                            grad(:,i,2)*lambda + grad(:,i,3)
            enddo
            !$OMP end do
            !$OMP end parallel
          end select
        elseif ( cluster_method .eq. "LMD2" ) then
          select case (nfrg)
          case (1)
            !$OMP parallel default(none)                             & !
            !$OMP private(i)                                         & !
            !$OMP shared(grad,lambda,ixnatm)
            !$OMP do schedule (static)
            do i = 1,ixnatm
              grad(:,i,1) = grad(:,i,1)*lambda
            enddo
            !$OMP end do
            !$OMP end parallel
          case (3)
            rtmp = sqrt(lambda)
            !$OMP parallel default(none)                             & !
            !$OMP private(i)                                         & !
            !$OMP shared(grad,rtmp,lambda,ixnatm)
            !$OMP do schedule (static)
            do i = 1,ixnatm
              grad(:,i,3) = grad(:,i,1)*lambda +                       &
                            grad(:,i,2)*rtmp + grad(:,i,3)
            enddo
            !$OMP end do
            !$OMP end parallel
          end select
        endif
        
        rss = 0.d0 ; dirgrd = 0.d0
        !$OMP parallel default (none)                                & !
        !$OMP private(ivar,iatm)                                     & !
        !$OMP shared(iynvar,iylvar,grad,srdir,nfrg)                  & !
        !$OMP reduction (+ : dirgrd,rss)
        !$OMP do schedule (static)
        do ivar = 1,iynvar
          iatm = iylvar(ivar)
          dirgrd = dirgrd + grad(1,iatm,nfrg)*srdir(1,ivar) +          &
                            grad(2,iatm,nfrg)*srdir(2,ivar) +          &
                            grad(3,iatm,nfrg)*srdir(3,ivar)
          rss = rss + grad(1,iatm,nfrg)*grad(1,iatm,nfrg) +            &
                      grad(2,iatm,nfrg)*grad(2,iatm,nfrg) +            &
                      grad(3,iatm,nfrg)*grad(3,iatm,nfrg)
        enddo
        !$OMP end do
        !$OMP end parallel
        rss = sqrt(rss/dble(iynvar))

        ! 4) MONITORING
        if ( mod(inloop,inmons) .eq. 0 )                               &
          call monmin(iprint,inloop,stplen,inflog,infbst,inuana,fnepot)

        ! 5) CHECK LOOP LIMIT AND CPU LIMIT
        call system_clock(sec)
        if ( rss .le. fncovg ) then
          write(iprint,*)" "
          write(iprint,*)"INFORMATION> CONJ "
          write(iprint,*)"     MINIMIZATION IS CONVERGED "
          write(iprint,*)" "
          onconv = .true. ; return
        endif
        if ( iloop .ge. inlopl ) then
          write(iprint,*)" "
          write(iprint,*)"INFORMATION> CONJ "
          write(iprint,*)"     MINIMIZATION IS NOT CONVERGED WITHIN"// &
                         " LOOP LIMIT "
          write(iprint,*)" "
          onlopl = .true. ; return
        endif
        if ( jloop .eq. inllpl ) return
        if ( sec .gt. fxcpus ) then
          Ttime = dble(sec-fxcpus)/dble(cr)
        else
          Ttime = dble(cm-fxcpus+sec)/dble(cr)
        endif
        if ( Ttime .gt. fxcpul ) then
          write(iprint,*)" "
          write(iprint,*)"WARNING> CONJ "
          write(iprint,*)"     CPU TIME IS EXCEEDED "
          write(iprint,*)" "
          ier = 1 ; return
        endif
        if ( abs(dirgrd/dirgrs) .le. fncovl ) return

        ! 6) SET NEXT STEP LENGTH
        if ( .not. onfind ) then
          if ( dirgrd .lt. 0.d0 ) then
            stplnn = stplen ; dirgrn = dirgrd
            stplen = stplnn*2.d0
          else
            onfind = .true.
            stplnp = stplen ; dirgrp = dirgrd
            stplen = stplnn - (stplnp-stplnn)/(dirgrp-dirgrn)*dirgrn
          endif
        else
          if ( dirgrd .lt. 0.d0 ) then
            stplnn = stplen ; dirgrn = dirgrd
          else
            stplnp = stplen ; dirgrp = dirgrd
          endif
          stplen = stplnn - (stplnp-stplnn)/(dirgrp-dirgrn)*dirgrn
        endif

      enddo

!***************************************************

      return 
      end subroutine linsrc
