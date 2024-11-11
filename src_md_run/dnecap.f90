
      subroutine dnecap(ecap)

!***********************************************************************
!
!     CALCULATION OF CAP-CONSTRAINED ENERGY AND GRADIENT
!
!***********************************************************************
 
      use COMBAS ; use COMERG ; use COMMIS ; use PHYCNS 
      use COMCMMC,only: cord,grad ; use COMCMM,only: nfrg
      use CALC_TIME
      !$ use omp_lib

      implicit none 

      real(8),intent(out):: ecap    ! CAP constraint energy (KCAL/MOL)

      real(8):: Tcod(3),Tcord(3,nchncap),Tgrad(3,nchncap)!,Tmass(nchncap)
      real(8):: dist,dist2,dist2i,dist3,dist4,tmprat,rr,wrk3,rtmp,Tecap
      real(8):: dx2,dy2,dz2
      integer(4):: i,ichn,iatm,ist,ien,j
      logical(4),save:: fstF = .true.
      logical(4):: ffstF
      integer(4),save,allocatable:: fechn(:,:)
      real(8),save,allocatable:: Tmass(:)

!*********************************************************

      if ( fstF ) then
        fstF = .false. ; ffstF = .true.
        allocate(fechn(2,nchncap),Tmass(nchncap))
        iatm = ixtatb
        OUTER : do i = 1,nchncap
          ichn = ixbfcn(i)
          do
            if ( iatm .gt. ixnatm ) exit OUTER
            if ( ichn .eq. ixachn(iatm) ) then
              if ( ffstF ) then
                fechn(1:2,i) = (/iatm,iatm/) ; ffstF = .false.
              else
                fechn(2,i) = iatm
              endif
            else
              ffstF = .true. ; exit
            endif
            iatm = iatm + 1
          enddo
        enddo OUTER
        !$OMP workshare
        forall(i=1:nchncap)                                            &
          Tmass(i) = 1.d0/sum(fxmass(fechn(1,i):fechn(2,i)))
        !$OMP end workshare
      endif

!     <<<  INITIALIZATION  >>>
      ecap = 0.d0 ; Tecap = 0.d0

!     <<<  CALCULATE MASS-CENTER OF EACH CHAIN  & Tgrad INITIALIZATION >>>
      !$OMP parallel default (none)                                  & !
      !$OMP private(i,ist,ien)                                       & !
      !$OMP shared(nchncap,fechn,Tcord,cord,fxmass,Tmass,fxcbou,Tgrad)
      !$OMP do schedule (static)
      do i = 1,nchncap
        ist = fechn(1,i) ; ien = fechn(2,i)
        Tcord(1,i) = dot_product(cord(1,ist:ien),fxmass(ist:ien))
        Tcord(2,i) = dot_product(cord(2,ist:ien),fxmass(ist:ien))
        Tcord(3,i) = dot_product(cord(3,ist:ien),fxmass(ist:ien))
        Tcord(1:3,i) = Tcord(1:3,i) * Tmass(i) - fxcbou(1:3)
        Tgrad(:,i) = 0.d0
      enddo
      !$OMP end do
      !$OMP end parallel

      if ( capshp.gt.2 .and. iufcap.ne.1 ) then
        print*,"Sorry, AAA" ; stop
      endif

!     <<<  CALCULATE HALF-HARMONIC CAP-ENERGY AND GRADIENT  >>>
      if ( iufcap .eq. 1 ) then

        !! Sphere CAP boundary
        if ( capshp .eq. 1 ) then
          !$OMP parallel default (none)                              & !
          !$OMP private(dist,dist2,rtmp,iatm,Tcod)                   & !
          !$OMP shared(nchncap,Tcord,furcap2,furcap,fukcap,Tgrad,    & !
          !$OMP        natcap,cord,fxcbou,furcap_pro2,furcap_pro,    & !
          !$OMP        fukcap_pro,grad,iamcap,nfrg)                  & !
          !$OMP reduction (+ : ecap)
          !$OMP do schedule (static)
          do i = 1,nchncap
            dist = Tcord(1,i)*Tcord(1,i) + Tcord(2,i)*Tcord(2,i) +     &
                   Tcord(3,i)*Tcord(3,i)
            if ( dist .gt. furcap2 ) then
              dist = sqrt(dist)
              dist2 = dist - furcap
              rtmp = fukcap * dist2
              ecap = ecap + rtmp * dist2
              rtmp = rtmp / dist
              Tgrad(1:3,i) = Tcord(1:3,i) * rtmp
            endif
          enddo
          !$OMP end do nowait
          !$OMP do schedule (static)
          do i = 1,natcap
            iatm = iamcap(i)
            Tcod(1:3) = cord(1:3,iatm) - fxcbou(1:3)
            dist = Tcod(1)*Tcod(1) + Tcod(2)*Tcod(2) + Tcod(3)*Tcod(3)
            if ( dist .gt. furcap_pro2 ) then
              dist = sqrt(dist)
              dist2 = dist - furcap_pro
              rtmp = fukcap_pro * dist2
              ecap = ecap + rtmp * dist2
              rtmp = rtmp / dist
              grad(1:3,iatm,nfrg) = grad(1:3,iatm,nfrg) + rtmp*Tcod(1:3)
            endif
          enddo
          !$OMP end do
          !$OMP end parallel
          ecap = ecap * 0.5d0

        !! Box CAP boundary
        elseif ( capshp .eq. 2 ) then
          !$OMP parallel default (none)                              & !
          !$OMP private(i,dist,iatm,Tcod)                            & !
          !$OMP shared(nchncap,Tcord,celwal,fukcap,natcap,iamcap,    & !
          !$OMP        cord,celwal_pro,grad,Tgrad,fukcap_pro,nfrg)   & !
          !$OMP reduction (+ : ecap,Tecap)
          !$OMP do schedule (static)
          do i = 1,nchncap
            ! X
            if ( Tcord(1,i) .lt. celwal(1) ) then
              dist = Tcord(1,i) - celwal(1)
              ecap = ecap + dist * dist
              Tgrad(1,i) = dist * fukcap
            elseif ( Tcord(1,i) .gt. celwal(2)) then
              dist = Tcord(1,i) - celwal(2)
              ecap = ecap + dist * dist
              Tgrad(1,i) = dist * fukcap
            endif
            ! Y
            if ( Tcord(2,i) .lt. celwal(3) ) then
              dist = Tcord(2,i) - celwal(3)
              ecap = ecap + dist * dist
              Tgrad(2,i) = dist * fukcap
            elseif ( Tcord(2,i) .gt. celwal(4)) then
              dist = Tcord(2,i) - celwal(4)
              ecap = ecap + dist * dist
              Tgrad(2,i) = dist * fukcap
            endif
            ! Z
            if ( Tcord(3,i) .lt. celwal(5) ) then
              dist = Tcord(3,i) - celwal(5)
              ecap = ecap + dist * dist
              Tgrad(3,i) = dist * fukcap
            elseif ( Tcord(3,i) .gt. celwal(6)) then
              dist = Tcord(3,i) - celwal(6)
              ecap = ecap + dist * dist
              Tgrad(3,i) = dist * fukcap
            endif
          enddo
          !$OMP end do nowait
          !$OMP do schedule (static)
          do i = 1,natcap
            iatm = iamcap(i)
            Tcod(1:3) = cord(1:3,iatm)
            ! X
            if ( Tcod(1) .lt. celwal_pro(1) ) then
              dist = Tcod(1) - celwal_pro(1)
              Tecap = Tecap + dist * dist
              grad(1,iatm,nfrg) = grad(1,iatm,nfrg) + dist * fukcap_pro
            elseif ( Tcod(1) .gt. celwal_pro(2) ) then
              dist = Tcod(1) - celwal_pro(2)
              Tecap = Tecap + dist * dist
              grad(1,iatm,nfrg) = grad(1,iatm,nfrg) + dist * fukcap_pro
            endif
            ! Y
            if ( Tcod(2) .lt. celwal_pro(3) ) then
              dist = Tcod(2) - celwal_pro(3)
              Tecap = Tecap + dist * dist
              grad(2,iatm,nfrg) = grad(2,iatm,nfrg) + dist * fukcap_pro
            elseif ( Tcod(2) .gt. celwal_pro(4) ) then
              dist = Tcod(2) - celwal_pro(4)
              Tecap = Tecap + dist * dist
              grad(2,iatm,nfrg) = grad(2,iatm,nfrg) + dist * fukcap_pro
            endif
            ! Z
            if ( Tcod(3) .lt. celwal_pro(5) ) then
              dist = Tcod(3) - celwal_pro(5)
              Tecap = Tecap + dist * dist
              grad(3,iatm,nfrg) = grad(3,iatm,nfrg) + dist * fukcap_pro
            elseif ( Tcod(3) .gt. celwal_pro(6) ) then
              dist = Tcod(3) - celwal_pro(6)
              Tecap = Tecap + dist * dist
              grad(3,iatm,nfrg) = grad(3,iatm,nfrg) + dist * fukcap_pro
            endif
          enddo
          !$OMP end do
          !$OMP end parallel
          ecap = ecap * fukcap
          Tecap = Tecap * fukcap_pro
          ecap = (ecap + Tecap) * 0.5d0

        ! ELLipsoidal CAP boundary
        elseif ( capshp .eq. 3 ) then

          !$OMP parallel default (none)                              & !
          !$OMP private(i,Tcord,dx2,dy2,dz2,dist2,rtmp,dist,iatm,    & !
          !$OMP         Tcod)                                        & !
          !$OMP shared(nchncap,fxellp,Tgrad,natcap,iamcap,cord,      & !
          !$OMP        fxcbou,fxellp_pro,grad,fukcap,fukcap_pro,nfrg)& !
          !$OMP reduction (+ : ecap,Tecap)
          !$OMP do schedule (static)
          do i = 1,nchncap
            dx2 = Tcord(1,i)*Tcord(1,i)
            dy2 = Tcord(2,i)*Tcord(2,i)
            dz2 = Tcord(3,i)*Tcord(3,i)
            dist2 = dx2*fxellp(1) + dy2*fxellp(2) + dz2*fxellp(3)
            if ( dist2 .gt. 1.d0 ) then
              rtmp = 1.d0 - sqrt(1.d0/dist2)
              dist = sqrt(dx2+dy2+dz2)*rtmp
              ecap = ecap + dist * dist
              Tgrad(1:3,i) = Tcord(1:3,i) * rtmp * fukcap
            endif
          enddo
          !$OMP end do nowait
          !$OMP do schedule (static)
          do i = 1,natcap
            iatm = iamcap(i)
            Tcod(1:3) = cord(1:3,iatm) - fxcbou(1:3)
            dx2 = Tcod(1)*Tcod(1)
            dy2 = Tcod(2)*Tcod(2)
            dz2 = Tcod(3)*Tcod(3)
            dist2 = dx2*fxellp_pro(1) + dy2*fxellp_pro(2) +            &
                    dz2*fxellp_pro(3)
            if ( dist2 .gt. 1.d0 ) then
              rtmp = 1.d0 - sqrt(1.d0/dist2)
              dist = sqrt(dx2+dy2+dz2)*rtmp
              Tecap = Tecap + dist * dist
              grad(1:3,iatm,nfrg) = grad(1:3,iatm,nfrg) +              &
                                    Tcod(1:3)*rtmp*fukcap_pro
            endif
          enddo
          !$OMP end do
          !$OMP end parallel
          ecap = ecap * fukcap
          Tecap = Tecap * fukcap_pro
          ecap = (ecap + Tecap) * 0.5d0

        endif
 
!     <<<  CALCULATE HALF-BIQUADRATIC CAP-ENERGY AND GRADIENT  >>>
      elseif ( iufcap .eq. 2 ) then
        dist2 = furcap * furcap ; dist2i = 1.d0 / dist2
        do i = 1,nchncap
          dist = Tcord(1,i)*Tcord(1,i) + Tcord(2,i)*Tcord(2,i) +       &
                 Tcord(3,i)*Tcord(3,i)
          dist3 = (dist - dist2) * dist2i
          dist3 = max(dist3,0.d0)
          ecap = ecap + dist3 * dist3 * dist2
          rtmp = fukcap * dist3
          Tgrad(1:3,i) = Tcord(1:3,i) * rtmp
        enddo
        ecap = ecap * fukcap

        dist2 = furcap_pro2 ; dist2i = 1.d0 / dist2
        do i = 1,natcap
          iatm = iamcap(i)
          Tcod(1:3) = cord(1:3,iatm) - fxcbou(1:3)
          dist = Tcod(1)*Tcod(1) + Tcod(2)*Tcod(2) + Tcod(3)*Tcod(3)
          if ( dist .gt. furcap_pro2 ) then
            dist3 = (dist - dist2) * dist2i
            Tecap = Tecap + dist3 * dist3 * dist2
            rtmp = fukcap_pro * dist3
            grad(1:3,iatm,nfrg) = grad(1:3,iatm,nfrg) + rtmp * Tcod(1:3)
          endif
        enddo
        Tecap = Tecap * fukcap_pro
        ecap = (ecap + Tecap) * 0.25d0

!     <<<  CALCULATION HALF-HARMONIC CAP-CONSTRAINT ENEGY  >>>
!     <<<  & FORCE TO STAY ABOVE XY PLANE  >>>
!              RCRATE = FXCBOU(3)/RCIR (slope)
      elseif ( iufcap .eq. 3 ) then
        do i = 1,nchncap
          wrk3 = Tcord(3,i)
          dist = sqrt(Tcord(1,i)*Tcord(1,i) + Tcord(2,i)*Tcord(2,i) +  &
                       wrk3*wrk3)
          if ( dist.gt.furcap .or. Tcord(3,i).lt.0.d0 ) then !outside of shell
            dist2 = Tcord(1,i)*Tcord(1,i) + Tcord(2,i)*Tcord(2,i)
            if ( dist2.lt.rcir2 .and. Tcord(3,i).lt.0.d0 ) then ! within the cylinder
              ecap = ecap + Tcord(3,i)*Tcord(3,i)
              Tgrad(3,i) = Tcord(3,i) * fukcap
            else
              dist2 = sqrt(dist2)
              rr = -wrk3 / dist2
              if ( rr .lt. rcrate ) then !outside of sphere
                dist2 = dist - furcap
                ecap = ecap + dist2 * dist2
                rtmp = fukcap * dist2 / dist
                Tgrad(1:2,i) = Tcord(1:2,i) * rtmp
                Tgrad(3,i) = wrk3 * rtmp
              else
                rtmp = (dist2 - rcir) / dist2
                Tcord(1:2,i) = Tcord(1:2,i) * rtmp
                dist3 = Tcord(1,i)*Tcord(1,i) + Tcord(2,i)*Tcord(2,i)  &
                      + Tcord(3,i)*Tcord(3,i)
                ecap = ecap + dist3
                Tgrad(1:3,i) = Tcord(1:3,i) * fukcap
              endif
            endif
          endif
        enddo
        ecap = ecap * fukcap

        do i = 1,natcap
          iatm = iamcap(i)
          Tcod(1:2) = cord(1:2,iatm) - fxcbou(1:2)
          wrk3 = cord(3,iatm) - fxcbou(3)
          dist = sqrt(Tcod(1)*Tcod(1) + Tcod(2)*Tcod(2) + wrk3*wrk3)
          if ( dist.gt.furcap_pro .or. cord(3,iatm).lt.plane_pro ) then !outside of shell
            dist2 = Tcod(1)*Tcod(1) + Tcod(2)*Tcod(2)
            Tcod(3) = cord(3,iatm) - plane_pro
            if ( dist2.lt.rcir_pro2 .and.                              &
                 cord(3,iatm).lt.plane_pro) then !within the cylinder
              Tecap = Tecap + Tcod(3) * Tcod(3)
              grad(3,iatm,nfrg) = grad(3,iatm,nfrg) + fukcap_pro*Tcod(3)
            else
              dist2 = sqrt(dist2)
              rr = -wrk3 / dist2
              if ( rr .lt. rcrate_pro ) then !outside of sphere
                dist2 = dist - furcap_pro
                Tecap = Tecap + dist2 * dist2
                rtmp = fukcap_pro * dist2 / dist
                grad(1:2,iatm,nfrg) = grad(1:2,iatm,nfrg)+Tcod(1:2)*rtmp
                grad(3,iatm,nfrg) = grad(3,iatm,nfrg) + wrk3 * rtmp
              else
                rtmp = (dist2 - rcir_pro) / dist2
                Tcod(1:2) = Tcod(1:2) * rtmp
                dist3 = dist2*dist2 + Tcod(3)*Tcod(3)
                Tecap = Tecap + dist3
                grad(1:3,iatm,nfrg) = grad(1:3,iatm,nfrg)              &
                                      + fukcap_pro * Tcod(1:3)
              endif
            endif
          endif
        enddo
        Tecap = Tecap * fukcap_pro
        ecap = (ecap + Tecap) * 0.5d0

!     <<<  CALCULATION HALF-BIQUADRATIC CAP-CONSTRAINT ENEGY  >>>
!     <<<  & FORCE TO STAY ABOVE XY PLANE  >>>
      elseif ( iufcap .eq. 4 ) then

        dist2 = furcap * furcap ; dist2i = 1.d0 / dist2
        do i = 1,nchncap
          wrk3 = Tcord(3,i)
          dist = Tcord(1,i)*Tcord(1,i)+Tcord(2,i)*Tcord(2,i) + wrk3*wrk3
          if ( dist.gt.dist2 .or. Tcord(3,i).lt.0.d0 ) then !outside of shell
            dist4 = dist - wrk3*wrk3
            if ( dist4.lt.rcir2 .and. Tcord(3,i).lt.0.d0 ) then !within the cylinder
              dist3 = Tcord(3,i)*Tcord(3,i) * dist2i
              ecap = ecap + dist3 * dist3 * dist2
              Tgrad(3,i) = Tcord(3,i) * fukcap * dist3
            else
              dist4 = sqrt(dist4)
              rr = -wrk3 / dist4
              if ( rr .lt. rcrate ) then
                dist3 = (dist - dist2) * dist2i
                ecap = ecap + dist3 * dist3 * dist2
                rtmp = fukcap * dist3
                Tgrad(1:2,i) = Tcord(1:2,i) * rtmp
                Tgrad(3,i) = wrk3 * rtmp
              else
                rtmp = (dist4 - rcir) / dist4
                Tcord(1:2,i) = Tcord(1:2,i) * rtmp
                dist = dist4*dist4 + Tcord(3,i)*Tcord(3,i)
                dist3 = (dist - dist2) * dist2i
                ecap = ecap + dist3 * dist3 * dist2
                rtmp = fukcap * dist3
                Tgrad(1:3,i) = Tcord(1:3,i) * rtmp
              endif
            endif
          endif
        enddo
        ecap = ecap * fukcap

        dist2 = furcap_pro2 ; dist2i = 1.d0 / dist2
        do i = 1,natcap
          iatm = iamcap(i)
          Tcod(1:2) = cord(1:2,iatm) - fxcbou(1:2)
          wrk3 = cord(3,iatm) - fxcbou(3)
          dist = Tcod(1)*Tcod(1) + Tcod(2)*Tcod(2) + wrk3*wrk3
          if ( dist.gt.dist2 .or. cord(3,iatm).lt.plane_pro ) then !outside of shell
            dist4 = Tcod(1)*Tcod(1) + Tcod(2)*Tcod(2)
            Tcod(3) = cord(3,iatm) - plane_pro
            if ( dist4 .lt. rcir_pro2 .and.                            &
                 cord(3,iatm).lt.plane_pro ) then !within the cylinder
              dist3 = Tcod(3)*Tcod(3) * dist2i
              Tecap = Tecap + Tcod(3)*Tcod(3) * 2.d0
              grad(3,iatm,nfrg) = grad(3,iatm,nfrg) + Tcod(3)*fukcap_pro
            else
              dist4 = sqrt(dist4)
              rr = -wrk3 / dist4
              if ( rr .lt. rcrate_pro ) then !outside of sphere
                dist3 = (dist - dist2) * dist2i
                Tecap = Tecap + dist3 * dist3 * dist2
                rtmp = fukcap_pro * dist3
                grad(1:2,iatm,nfrg) = grad(1:2,iatm,nfrg)+Tcod(1:2)*rtmp
                grad(3,iatm,nfrg) = grad(3,iatm,nfrg) + wrk3 * rtmp
              else
                rtmp = (dist4 - rcir_pro) / dist4
                Tcod(1:2) = Tcod(1:2) * rtmp
                dist = dist4*dist4 + Tcod(2)*Tcod(3)
                dist3 = (dist - dist2) * dist2i
                Tecap = Tecap + dist3 * dist3 * dist2
                rtmp = fukcap_pro * dist3
                grad(1:3,iatm,nfrg) = grad(1:3,iatm,nfrg)*Tcod(1:3)*rtmp
              endif
            endif
          endif
        enddo
        Tecap = Tecap * fukcap_pro
        ecap = (ecap + Tecap) * 0.25d0

      endif

!     <<<  STORE GRADINET  >>>
      !$OMP parallel default(none)                                   & !
      !$OMP private(i,j,tmprat)                                      & !
      !$OMP shared(nchncap,fechn,fxmass,Tmass,grad,Tgrad,nfrg)
      !$OMP do
      do i = 1,nchncap
      do j = fechn(1,i),fechn(2,i)
        tmprat = fxmass(j) * Tmass(i)
        grad(1:3,j,nfrg) = grad(1:3,j,nfrg) + Tgrad(1:3,i)*tmprat
      enddo
      enddo
      !$OMP end do
      !$OMP end parallel

!************************************************

      return
      end subroutine dnecap
