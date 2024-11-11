
      subroutine ellpbc(maxatm,ixnchn,ixcend,cordm1,fxmass,velback,    &
                        iynvar,iylvar,icvara,deltt,fxcbou,fxellp)

!********************************************************************
!
!     THE ELASTIC COLLISION AT THE ELLIPSOIDAL SURFACE
!     (X-XC)**2/A**2 + (Y-YC)**2/B**2 + (Z-ZC)**2/C**2 = 1
!       [X-Z]C = FXCBOU(1:3), [A-C]  = FXELLP(1:3)
!
!********************************************************************
 
      use COMBAS,only: ixnatm ; use COMCMMC,only: cord,vel

      implicit none

      integer(4):: maxatm,ixnchn,ixcend(ixnchn),iynvar,iylvar(iynvar), &
                   icvara(maxatm)
      real(8):: cordm1(3,ixnatm),fxmass(maxatm),velback(3,ixnatm),     &
                deltt,fxcbou(3),fxellp(3),ixellp(3)

      integer(4):: i,j,ichni,ica,iatm,ivar
      real(8):: cofmas(3),ichnmas,phisrf,oldmas(3),rtmp(3),rtmp2(3),   &
        direct(3),coeffa,coeffb,coeffc,delttp,ellppt(3),ellnrm(3),     &
        fctnrm,velcnt(3),ellort(3),fctort,elltng(3),velcnw(3),velnrm,  &
        veltng,velold(3,ixnatm)

!******************************************
!     <<<  ELLIPSOID BOUNDARY >>>
!     backup the calculated velocity
      velold(1:3,1:iynvar) = vel(1:3,1:iynvar)
      ixellp(1:3) = 1.d0 / fxellp(1:3)

!     find the chains, whose centers of mass are outside.
      do i = 1, ixnchn
        if ( i .eq. 1 ) then
          ichni = 1
        else
          ichni = ixcend(i-1) + 1
        endif
        cofmas(1:3) = 0.d0 ; ichnmas = 1.d0/sum(fxmass(ichni:ixcend(i)))
        do j = ichni,ixcend(i)
          cofmas(1:3) = cofmas(1:3) + cord(1:3,j)*fxmass(j)
        enddo
        cofmas(1:3) = cofmas(1:3) * ichnmas
        phisrf = ((cofmas(1)-fxcbou(1)) * ixellp(1))**2 +              &
                 ((cofmas(2)-fxcbou(2)) * ixellp(2))**2 +              &
                 ((cofmas(3)-fxcbou(3)) * ixellp(3))**2 - 1.d0
        if ( phisrf .gt. 0.d0 ) then
!         find the surface point(ellppt) and delttp
          oldmas(1:3) = 0.d0
          do j = ichni,ixcend(i)
            oldmas(1:3) = oldmas(1:3) + cordm1(1:3,j)*fxmass(j)
          enddo
          oldmas(1:3) = oldmas(1:3) * ichnmas
          direct(1:3) = cofmas(1:3) - oldmas(1:3)
          coeffa = (direct(1)*ixellp(1))**2 + (direct(2)*ixellp(2))**2 &
                 + (direct(3)*ixellp(3))**2
          rtmp(1:3) = oldmas(1:3)-fxcbou(1:3)
          rtmp2(1:3) = rtmp(1:3)*ixellp(1:3)**2
          coeffb = dot_product(rtmp2,direct)
          coeffc = dot_product(rtmp,rtmp2) - 1.d0
          delttp = (-coeffb + sqrt(coeffb*coeffb-coeffa*coeffc))/coeffa
          ellppt(1:3) = oldmas(1:3) + delttp*direct(1:3)

!         set normal vector at ellppt
          ellnrm(1:3) = (ellppt(1:3)-fxcbou(1:3))*ixellp(1:3)**2
          fctnrm = 1.d0 / sqrt(dot_product(ellnrm,ellnrm))
          ellnrm(1:3) = ellnrm(1:3)*fctnrm

!         find free atoms to be re-calculated and set 
!         velocity of center of mass: velcnt
          ica = 0 ; velcnt(1:3) = 0.d0
          do j = 1,iynvar
            iatm = iylvar(j)
            if ( iatm.ge.ichni .and. iatm.le.ixcend(i) ) then
              ica = ica + 1 ; icvara(ica) = j
              velcnt(1:3) = velcnt(1:3) + fxmass(iatm)*vel(1:3,j)
            endif
          enddo
          velcnt(1:3) = velcnt(1:3) * ichnmas

!         reset velocity of the center of mass
          ellort(1) = velcnt(2)*ellnrm(3) - velcnt(3)*ellnrm(2)
          ellort(2) = velcnt(3)*ellnrm(1) - velcnt(1)*ellnrm(3)
          ellort(3) = velcnt(1)*ellnrm(2) - velcnt(2)*ellnrm(1)
          fctort = sqrt(dot_product(ellort,ellort))
          if ( fctort .ne. 0.d0 ) then
            fctort = 1.d0 / fctort ; ellort(1:3) = ellort(1:3)*fctort
            elltng(1) = ellnrm(2)*ellort(3) - ellnrm(3)*ellort(2)
            elltng(2) = ellnrm(3)*ellort(1) - ellnrm(1)*ellort(3)
            elltng(3) = ellnrm(1)*ellort(2) - ellnrm(2)*ellort(1)
            veltng = dot_product(velcnt,elltng)
            velnrm = dot_product(velcnt,ellnrm)
            velcnw(1:3) = -velnrm*ellnrm(1:3) + veltng*elltng(1:3)
          else
            velcnw(1:3) = -velcnt(1:3)
          endif
          do j = 1,ica
            ivar = icvara(j) ; iatm = iylvar(ivar)
            vel(1:3,ivar) = vel(1:3,ivar) - velcnt(1:3) + velcnw(1:3)
            velback(1:3,ivar) = velback(1:3,ivar) + velcnt(1:3) -      &
                                velcnw(1:3)
            cord(1:3,iatm) = delttp*deltt*velold(1:3,ivar) +           &
              (1.d0-delttp)*deltt*vel(1:3,ivar) + cordm1(1:3,iatm)
          enddo
        endif
      enddo

!*****************************************

      return
      end subroutine ellpbc
