
      subroutine dneres_PB(enevcc) ! #SLC2 #PB
      subroutine dneres_CB(enevcc) ! #SLC2 #CB

!***********************************************************************
!
!     CALCULATION OF RESTRAINED ENERGY AND GRADIENT
!       POSITION RESTRAINT ENERGY
!       DISTANCE RESTRAINT ENERGY
!       TORSIONAL RESTRAINT ENERGY (NOT-AVILABLE NOW)
!       CAP-RESTRAINT ENERGY
!
!***********************************************************************
 
      use COMBAS ; use COMERG
      use COMCMM, only: nfrg

      implicit none

      real(8),intent(inout):: enevcc(maxene,nfrg)

!************************************************************

      ! Position
      if ( iyeflg(11) .eq. 1 ) call dposc_PB(enevcc(11,nfrg))  ! #SLC2 #PB
      if ( iyeflg(11) .eq. 1 ) call dposc_CB(enevcc(11,nfrg))  ! #SLC2 #CB
      ! Distance
      if ( iyeflg(12) .eq. 1 ) call ddistc_PB(enevcc(12,nfrg)) ! #SLC2 #PB
      if ( iyeflg(12) .eq. 1 ) call ddistc_CB(enevcc(12,nfrg)) ! #SLC2 #CB
      ! Torsion
      if ( iyeflg(13) .eq. 1 ) call dcdihe_PB(enevcc(13,nfrg)) ! #SLC2 #PB
      if ( iyeflg(13) .eq. 1 ) call dcdihe_CB(enevcc(13,nfrg)) ! #SLC2 #CB
      ! CAP
      if ( iyeflg(14) .eq. 1 ) call dnecap(enevcc(14,nfrg))    ! #SLC2 #CB
      ! GRAV
      if ( iyeflg(15) .eq. 1 ) call dnerep_PB(enevcc(15,nfrg))! #SLC2 #PB
      if ( iyeflg(15) .eq. 1 ) call dnerep_CB(enevcc(15,nfrg))! #SLC2 #CB

!**********************************

      return
      end subroutine dneres_PB ! #SLC2 #PB
      end subroutine dneres_CB ! #SLC2 #CB


!========================================================================

      subroutine dposc_PB(epsc) ! #SLC2 #PB
      subroutine dposc_CB(epsc) ! #SLC2 #CB

!********************************************************************
!
!     SUBROUTINE TO GET POSITIONAL CONSTRAINT ENERGY AND GRADIENT FOR
!     THE FOLLOWING POTENTIALS
!
!********************************************************************

      use COMBAS ; use COMERG ; use COMMIS ; use COMCMMC
      use COMCMM,only: nfrg

      implicit none

      real(8),intent(inout):: epsc

      real(8):: dx,dy,dz,r2,rtmp,work(3)
      integer(4):: i,iatm

!****************************************************
      ! Initialization

!     --- CALCULATE THE POSITIONAL CONSTRAINTS ENERGY AND GRADIENTS ---
!     FUGCNS(I)   : MASS*DELTA OF ATOM IUGPNT(I)
!     BOLTZK : BOLTZMAN CONSTANT ( KCAL/MOL/K )
!     FUCTMP : THE ABSOLUTE TEMPERATURE   ( K )
!     MASS   : MASS OF ATOM I
!     FUWPSC : SCALE FACTOR
!     DELTA  : THE ERROR ESTIMATES OF ATOM I
!              ( /DALTON/A**2 ; MASS WEIGHT OR /A**2 ; NO MASS )
      epsc = 0.d0
      do i = 1,iutatm
        iatm = iugpnt(i)
        dx = cord(1,iatm) - fucord(1,iatm)
        dy = cord(2,iatm) - fucord(2,iatm)
        dz = cord(3,iatm) - fucord(3,iatm)
        dx = dx - fxcell(1) * nint(dx*invcel(1)) ! #SLC2 #PB
        dy = dy - fxcell(2) * nint(dy*invcel(2)) ! #SLC2 #PB
        dz = dz - fxcell(3) * nint(dz*invcel(3)) ! #SLC2 #PB
        r2 = dx*dx + dy*dy + dz*dz
        rtmp = fugcns(i)
        epsc = epsc + rtmp*r2

        work(1)=rtmp*dx ; work(2)=rtmp*dy ; work(3)=rtmp*dz
        grad(1:3,iatm,nfrg) = grad(1:3,iatm,nfrg) + work(1:3)
      enddo
      epsc = epsc * 0.5d0

!***************************************

      return
      end subroutine dposc_PB ! #SLC2 #PB
      end subroutine dposc_CB ! #SLC2 #CB


!============================================================================


      subroutine ddistc_PB(enoe) ! #SLC2 #PB
      subroutine ddistc_CB(enoe) ! #SLC2 #CB

!************************************************************************
!
!     SUBROUTINE TO GET DISTANCE CONSTRAINT ENERGY AND GRADIENT FOR
!       THE FOLLOWING POTENTIALS
!
!      ENOE =  KU*( (<R>**2 - RU**2)/RU**2 )**2  IF <R>.GT.RU
!           =  KL*( (<R>**2 - RL**2)/RL**2 )**2  IF <R>.LT.RL
!           =  0.0              IF (RL.GE.<R>).AND.(<R>.LE.RU)
!
!      <R> : MEAN DISTANCE = [< RIJ**(-6) >]**(-1/6)
!
!      KU  = CONST/(DUPR*DUPR)
!      KL  = CONST/(DLOW*DLOW)
!
!            DUPR   : THE POSITIVE ERROR ESTIMATES
!            DLOW   : THE NEGATIVE ERROR ESTIMATES
!
!----------------------------------------------------------------------C

      use COMBAS ; use COMERG ; use COMMIS ; use PHYCNS ; use COMCMMC
      use COMCMM,only: nfrg

      implicit none

      ! Total distance restraint energy
        real(8),intent(out):: enoe

      integer(4):: i,k,l,itmp,itmp2
      real(8):: dsccns,d(3),totr6,coef,rtmp
      real(8):: const(iutdsc),r0(iutdsc),r6ave(iutdsc),dr2(iutdsc),td(3)

!******************************************************

      enoe = 0.d0
      ! CALCULATE [< R(-6) >]**(-1/6) MEAN DISTANCE
      do i = 1,iutdsc
        totr6 = 0.d0
        do k = 1,iuipar(i)
          itmp = iuipnt(i,k)
          do l = 1,iujpar(i)
            itmp2 = iujpnt(i,l)
            d(1:3) = cord(1:3,itmp) - cord(1:3,itmp2)
            d(1:3) = d(1:3) -                             & ! #SLC2 #PB
                     fxcell(1:3) * nint(d(1:3)*invcel(1:3)) ! #SLC2 #PB
            rtmp = d(1)*d(1) + d(2)*d(2) + d(3)*d(3)
            rtmp = rtmp * rtmp * rtmp
            totr6 = totr6 + 1.d0 / rtmp
          enddo
        enddo
        r6ave(i)= totr6
      enddo

      ! MAIN LOOP TO CALCULATE DISTANCE RESTRAINT ENERGY
      r6ave(1:iutdsc) = r6ave(1:iutdsc) /                              &
                        dble(iuipar(1:iutdsc)*iujpar(1:iutdsc))
      r6ave(1:iutdsc) = r6ave(1:iutdsc)**(-1.d0/6.d0)
      r0(1:iutdsc) = 1.d0 ; const(1:iutdsc) = 0.d0
      dsccns = boltzk*fuctmp*fuwdsc*0.5d0
      do i = 1,iutdsc
        if ( r6ave(i) .gt. furupr(i) ) then
          r0(i) = furupr(i)
          const(i) = dsccns * fudupr(i)
!          const(i) = dsccns / (fudupr(i)*fudupr(i))
        elseif ( r6ave(i) .lt. furlow(i) ) then
          r0(i) = furlow(i)
          const(i) = dsccns * fudlow(i)
!          const(i) = dsccns / (fudlow(i)*fudlow(i))
        endif
        dr2(i) = r6ave(i)*r6ave(i) - r0(i)*r0(i)
        enoe = enoe + const(i)*dr2(i)*dr2(i)
      enddo

      do i = 1,iutdsc
        if ( r6ave(i).ge.furlow(i) .and. r6ave(i).le.furupr(i) ) cycle
        rtmp = r6ave(i) * r6ave(i)
        rtmp = rtmp * rtmp ; rtmp = rtmp * rtmp
        coef = 4.d0*const(i)*dr2(i)*rtmp / (dble(iuipar(i)*iujpar(i)))
        do k = 1,iuipar(i)
          td(1:3) = 0.d0 ; itmp = iuipnt(i,k)
          do l = 1,iujpar(i)
            itmp2 = iujpnt(i,l)
            d(1:3) = cord(1:3,itmp) - cord(1:3,itmp2)
            d(1:3) = d(1:3) -                             & ! #SLC2 #PB
                     fxcell(1:3) * nint(d(1:3)*invcel(1:3)) ! #SLC2 #PB
            rtmp = d(1)*d(1) + d(2)*d(2) + d(3)*d(3)
            rtmp = rtmp * rtmp ; rtmp = rtmp * rtmp
            rtmp = 1.d0 / rtmp
            td(1:3) = td(1:3) + d(1:3) * rtmp
          enddo
          grad(1:3,itmp,nfrg) = grad(1:3,itmp,nfrg) + coef*td(1:3)
        enddo
        do k = 1,iujpar(i)
          td(1:3) = 0.d0 ; itmp = iujpnt(i,k)
          do l = 1,iuipar(i)
            itmp2 = iuipnt(i,l)
            d(1:3) = cord(1:3,itmp) - cord(1:3,itmp2)
            d(1:3) = d(1:3) -                             & ! #SLC2 #PB
                     fxcell(1:3) * nint(d(1:3)*invcel(1:3)) ! #SLC2 #PB
            rtmp = d(1)*d(1) + d(2)*d(2) + d(3)*d(3)
            rtmp = rtmp * rtmp ; rtmp = rtmp * rtmp
            rtmp = 1.d0 / rtmp
            td(1:3) = td(1:3) + d(1:3) * rtmp
          enddo
          grad(1:3,itmp,nfrg) = grad(1:3,itmp,nfrg) + coef*td(1:3)
        enddo
      enddo

!********************************

      return
      end subroutine ddistc_PB ! #SLC2 #PB
      end subroutine ddistc_CB ! #SLC2 #CB


!========================================================================


      subroutine dcdihe_PB(edhc) ! #SLC2 #PB
      subroutine dcdihe_CB(edhc) ! #SLC2 #CB

!**********************************************************************
!
!      THIS SUBROUTINE IS FOR CALCULATION TORSIONAL CONSTRAINT ENERGY
!      AND GRADIENTS.
!
!   ECDIEH =  KU*( PHI - PHIUP )**2  IF PHI.GT.PHIUP
!          =  0.0                    IF(PHI.GE.PHILW).AND.(PHI.LE.PHIUP)
!          =  KL*( PHI - PHILW )**2  IF PHI.LT.PHILW
!
!      KU  =  BOLTZK*FUCTMP*FUWDHC/2*(FUDELU(I)*FUDELU(I))
!      KL  =  BOLTZK*FUCTMP*FUWDHC/2*(FUDELL(I)*FUDELL(I))
!
!            BOLTZK : BOLTZMAN CONSTANT ( KCAL/MOL )
!            FUCTMP : THE ABSOLUTE TEMPERATURE( K )
!            FUWDHC : THE SCALEING FACTOR
!            FUDELU : THE POSITIVE ERROR ESTIMATES ( RAD**2 )
!            FUDELL : THE NEGATIVE ERROR ESTIMATES ( RAD**2 )
!
!**********************************************************************

      use COMBAS ; use COMERG ; use COMMIS ; use PHYCNS ; use COMCMMC
      use COMCMM,only: nfrg

      implicit none

      real(8),intent(inout):: edhc

      integer(4):: i
      real(8):: work(9),wrk0,wrk1,wrk2,wrk3,dx21,dy21,dz21,dx32,dy32,  &
                dz32,dx43,dy43,dz43,px12,py12,pz12,px23,py23,pz23,     &
                px123,py123,pz123,r12,r23,r12r23,s1223,s32123,cosp,phi,&
                const,cons0,cons1,phi0,phi1,delphi,delph0,delph1
 
!*****************************************

      edhc = 0.d0
      do i = 1,iutdhc
        ! CALCULATE VECTOR
        dx21 = cord(1,iudhcp(i,2)) - cord(1,iudhcp(i,1))
        dy21 = cord(2,iudhcp(i,2)) - cord(2,iudhcp(i,1))
        dz21 = cord(3,iudhcp(i,2)) - cord(3,iudhcp(i,1))
        dx32 = cord(1,iudhcp(i,3)) - cord(1,iudhcp(i,2))
        dy32 = cord(2,iudhcp(i,3)) - cord(2,iudhcp(i,2))
        dz32 = cord(3,iudhcp(i,3)) - cord(3,iudhcp(i,2))
        dx43 = cord(1,iudhcp(i,4)) - cord(1,iudhcp(i,3))
        dy43 = cord(2,iudhcp(i,4)) - cord(2,iudhcp(i,3))
        dz43 = cord(3,iudhcp(i,4)) - cord(3,iudhcp(i,3))
        dx21 = dx21 - fxcell(1) * nint(dx21*invcel(1)) ! #SLC2 #PB
        dy21 = dy21 - fxcell(2) * nint(dy21*invcel(2)) ! #SLC2 #PB
        dz21 = dz21 - fxcell(3) * nint(dz21*invcel(3)) ! #SLC2 #PB
        dx32 = dx32 - fxcell(1) * nint(dx32*invcel(1)) ! #SLC2 #PB
        dy32 = dy32 - fxcell(2) * nint(dy32*invcel(2)) ! #SLC2 #PB
        dz32 = dz32 - fxcell(3) * nint(dz32*invcel(3)) ! #SLC2 #PB
        dx43 = dx43 - fxcell(1) * nint(dx43*invcel(1)) ! #SLC2 #PB
        dy43 = dy43 - fxcell(2) * nint(dy43*invcel(2)) ! #SLC2 #PB
        dz43 = dz43 - fxcell(3) * nint(dz43*invcel(3)) ! #SLC2 #PB
        ! VECTOR PRODUCTS OF BONDS 1 AND 2, AND 3
        px12 = dy21*dz32 - dy32*dz21
        py12 = dz21*dx32 - dz32*dx21
        pz12 = dx21*dy32 - dx32*dy21
        px23 = dy32*dz43 - dy43*dz32
        py23 = dz32*dx43 - dz43*dx32
        pz23 = dx32*dy43 - dx43*dy32
        r12 = px12*px12 + py12*py12 + pz12*pz12
        r23 = px23*px23 + py23*py23 + pz23*pz23
        r12r23 = sqrt( r12 * r23 )

        ! IN THE CASE OF BOND ANGLE = 0 OR PI , TORSIONAL ANGLES CAN NOT
        !    BE CALCULATED. HENCE, TORSIONAL ENERGY IS NOT CALCULATED.
        if ( r12r23 .le. 0.d0 ) cycle
          ! CALCULATE COSINE(PHI)
          r12r23 = 1.d0 / r12r23
          ! SCALAR PRODUCT BETWEEN PROD12 AND PROD23
          s1223 = px12*px23 + py12*py23 + pz12*pz23
          ! COSINE PHI
          cosp = min(s1223*r12r23,1.d0)
          cosp = max(cosp,-1.d0)
          ! DEFINE SIGN OF PHI
          !  VECTOR PRODUCT OF PROD12 AND PROD23
          px123 = py12*pz23 - py23*pz12
          py123 = pz12*px23 - pz23*px12
          pz123 = px12*py23 - px23*py12
          ! SCALAR PRODUCT OF VECT32 AND PROD123
          s32123 = dx32*px123 + dy32*py123 + dz32*pz123
          ! DIHEDRAL ANGLE PHI ( 0 < PHI < 2*PI )
          phi = acos(cosp)
          if ( s32123 .lt. 0.d0 ) phi = 2.d0*pi - phi
          ! SWITCHING THE FORM OF A CONSTRAINT POTENTIAL
          ! Problem happens when the angle crosses 0 or 360 deg.
          ! Calculate from the closer wall.
          if ( phi .gt. fucdup(i) ) then
            phi0 = fucdup(i) ; phi1 = fucdlw(i) + 2.d0*pi
            cons0 = boltzk*fuctmp*fuwdhc/(2.d0*fudelu(i)*fudelu(i))
            cons1 = boltzk*fuctmp*fuwdhc/(2.d0*fudell(i)*fudell(i))
          elseif ( phi .lt. fucdlw(i) ) then
            phi0 = fucdup(i) - 2.d0*pi ; phi1 = fucdlw(i)
            cons0 = boltzk*fuctmp*fuwdhc/(2.d0*fudelu(i)*fudelu(i))
            cons1 = boltzk*fuctmp*fuwdhc/(2.d0*fudell(i)*fudell(i))
          else
            phi0 = 0.d0 ; phi1 = 0.d0 ; cons0 = 0.d0 ; cons1 = 0.d0
          endif

          ! CALCULATE TORSION CONSTRAINT ENERGY
          !   TORSION(iudhcp(I,1),iudhcp(I,2),iudhcp(I,3),iudhcp(I,4)))
          delph0 = phi - phi0 ; delph1 = phi1 - phi
          if ( delph0 .lt. delph1 ) then
            delphi = delph0 ; const = cons0
          else
            delphi = -delph1 ; const = cons1
          endif
          edhc = edhc + const*delphi*delphi
 
          ! CALCULATE CONSTANT OF GRADIENT
          !   IF PHI IS ZERO OR PI, THEN SET GRADIENT TO ZERO.
          if ( abs(phi).lt.eps .or. abs(phi-pi).lt.eps ) cycle

          r12 = s1223 / r12 ; r23 = s1223 / r23
          wrk0 = -2.d0*const*delphi/sin(phi)*r12r23
          wrk1 = dx21+dx32 ; wrk2 = dy21+dy32 ; wrk3 = dz21+dz32
       
          ! CALCULATE GRADIENT
          work(1) = dz32*py23-dy32*pz23 - r12*(dz32*py12-dy32*pz12)
          work(2) = dx32*pz23-dz32*px23 - r12*(dx32*pz12-dz32*px12)
          work(3) = dy32*px23-dx32*py23 - r12*(dy32*px12-dx32*py12)
          work(4) = -wrk3*py23 + dz43*py12 + wrk2*pz23 - dy43*pz12 -   &
            r12*(-wrk3*py12 + wrk2*pz12) - r23*(dz43*py23 - dy43*pz23)
          work(5) =  wrk3*px23 - dz43*px12 - wrk1*pz23 + dx43*pz12 -   &
            r12*( wrk3*px12 - wrk1*pz12) - r23*(dx43*pz23 - dz43*px23)
          work(6) = -wrk2*px23 + dy43*px12 + wrk1*py23 - dx43*py12 -   &
            r12*(-wrk2*px12 + wrk1*py12) - r23*(dy43*px23 - dx43*py23)
          work(7) = dz32*py12-dy32*pz12 - r23*(dz32*py23-dy32*pz23)
          work(8) = dx32*pz12-dz32*px12 - r23*(dx32*pz23-dz32*px23)
          work(9) = dy32*px12-dx32*py12 - r23*(dy32*px23-dx32*py23)
          work(1:9) = work(1:9)*wrk0
          grad(1:3,iudhcp(i,1),nfrg) = grad(1:3,iudhcp(i,1),nfrg)      &
                                         + work(1:3)
          grad(1:3,iudhcp(i,2),nfrg) = grad(1:3,iudhcp(i,2),nfrg)      &
                                         + work(4:6)
          grad(1:3,iudhcp(i,3),nfrg) = grad(1:3,iudhcp(i,3),nfrg)      &
                                     - (work(1:3)+work(4:6)+work(7:9))
          grad(1:3,iudhcp(i,4),nfrg) = grad(1:3,iudhcp(i,4),nfrg)      &
                                         + work(7:9)
       enddo

!******************************************

       return
       end subroutine dcdihe_PB ! #SLC2 #PB
       end subroutine dcdihe_CB ! #SLC2 #CB


!===============================================================================


      subroutine dnerep_PB(erep) ! #SLC2 #PB
      subroutine dnerep_CB(erep) ! #SLC2 #CB

      use COMBAS ; use COMCMMC ; use COMCMM ; use COMERG

      implicit none

      real(8),intent(out):: erep

      integer(4):: i,j,iat1,iat2
      real(8):: x,y,z,dx,dy,dz,fx,fy,fz,gx,gy,gz,dd,ddi,r,cc,ee

!***********************************************

      erep = 0.d0
      do j = 1,Nrep(2)
        iat2 = replst(j,2) ; r = fxmass(iat2)*Krep
        x = cord(1,iat2) ; y = cord(2,iat2) ; z = cord(3,iat2)
        do i = 1,ntbrep(j)
          iat1 = itbrep(i,j)
          dx = cord(1,iat1) - x
          dy = cord(2,iat1) - y
          dz = cord(3,iat1) - z
          dx = dx - fxcell(1)*nint(dx*invcel(1))
          dy = dy - fxcell(2)*nint(dy*invcel(2))
          dz = dz - fxcell(3)*nint(dz*invcel(3))
          dd = sqrt(dx*dx+dy*dy+dz*dz)
          ddi = 1.d0 / dd
          dd = max(fycutl-dd,0.d0)
          ee = r*fxmass(iat1)*dd
          erep = erep + ee*dd
          cc = -2.d0*ee*ddi
          fx = dx*cc ; fy = dy*cc ; fz = dz*cc
          gx = gx + fx ; gy = gy + fy ; gz = gz + fz
          grad(1,iat1,nfrg) = grad(1,iat1,nfrg) + fx
          grad(2,iat1,nfrg) = grad(2,iat1,nfrg) + fy
          grad(3,iat1,nfrg) = grad(3,iat1,nfrg) + fz
        enddo
        grad(1,iat2,nfrg) = grad(1,iat2,nfrg) - gx
        grad(2,iat2,nfrg) = grad(2,iat2,nfrg) - gy
        grad(3,iat2,nfrg) = grad(3,iat2,nfrg) - gz
      enddo

!*******************************************

      return
      end subroutine dnerep_PB ! #SLC2 #PB
      end subroutine dnerep_CB ! #SLC2 #CB
