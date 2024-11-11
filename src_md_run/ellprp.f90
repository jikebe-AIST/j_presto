
      subroutine ellprp(maxatm,iprint,ixnchn,ixcend,fxmass,fxcbou,     &
                        fxellp,ier)

!********************************************************************
!
!     THE ELASTIC COLLISION AT THE ELLIPSOIDAL SURFACE
!     (X-XC)**2/A**2 + (Y-YC)**2/B**2 + (Z-ZC)**2/C**2 = 1
!     Preparation routine at the initial state.
!     Anychains, the center of mass of which is outside the
!     ellipsoid or sphere, this routine notifies them and stop.
!      ([X-Z]C = FXCBOU(1:3)  [A-C] = FXELLP(1:3))
!
!********************************************************************

      use COMCMMC,only: cord

      implicit none

      integer(4):: maxatm,iprint,ixnchn,ixcend(ixnchn),ier
      real(8):: fxmass(maxatm),fxcbou(3),fxellp(3)

      integer(4):: ichn,ichni,iatm
      real(8):: cofmas(3),chnmas,phisrf

!***************************************

      ier = 0
!     <<<  ELLIPSOID BOUNDARY >>>
!     find the chains, whose centers of mass are outside.
      do ichn = 1,ixnchn
        if ( ichn .eq. 1 ) then
          ichni = 1
        else
          ichni = ixcend(ichn-1) + 1
        endif
        cofmas(1:3) = 0.d0
        do iatm = ichni,ixcend(ichn)
          cofmas(1:3) = cofmas(1:3) + cord(1:3,iatm)*fxmass(iatm)
        enddo
        chnmas = 1.d0 / sum(fxmass(ichni:ixcend(ichn)))
        cofmas(1:3) = cofmas(1:3) * chnmas
        phisrf = ((cofmas(1)-fxcbou(1)) / fxellp(1))**2 +              &
                 ((cofmas(2)-fxcbou(2)) / fxellp(2))**2 +              &
                 ((cofmas(3)-fxcbou(3)) / fxellp(3))**2 - 1.d0
        if ( phisrf .gt. 0.d0 ) then
          if ( ier .eq. 0 ) then
            write(iprint,*)'ERROR> BELOW ATOMS ARE OUTSIDE',           &
                           ' OF THE ELLIPSOID.'
          endif
          ier = ier + 1
          write(iprint,*)'  Chain outside of ellipsoid: ',ichn
          write(iprint,*)'  Outside atoms:'
          write(iprint,'(10(1x,i6))')(iatm,iatm=ichni,ixcend(ichn))
        endif
      enddo

      if ( ier .ne. 0 ) then
        write(iprint,*) '  Total ',ier,' chains are outside.'
        write(iprint,*) ' ' ; ier = -6
      endif

!********************************

      return
      end subroutine ellprp
