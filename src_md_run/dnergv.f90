
      subroutine dnergv_PB(enevcc) ! #SLC2 #PB
      subroutine dnergv_CB(enevcc) ! #SLC2 #CB

!***********************************************************************
!
!     CALCULATION OF ENERGY AND GRADIENT
!       BOND ENERGY   ANGLE ENERGY   TORSIONAL ENERGY
!       IMPROPER ENERGY
!       ELECTROSTATIC ENERGY (1-4,1-5)
!       VAN DER WAALS ENERGY (1-4,1-5)
!       HYDROGEN BONDED ENERGY (1-5)
!
!       AS ELECTROSTATIC ENERGY , CONSTANT DIELECTRIC CONSTANT
!       AND DISTANCE DEPENDENT DIELECTRIC CONSTANT ARE AVAILABLE
!       PERIODIC BOUNDARY CONDITION IS AVAILABLE
!
!***********************************************************************
 
      use COMBAS ; use COMERG ; use PHYCNS
      use COMCMM, only: nfrg,nlev

      implicit none

      ! Energies between each cluster (kcal/mol)
        real(8),intent(inout):: enevcc(maxene,nfrg)

!************************************************************

!     <<<  BOND ENERGY  >>>
      if ( nbnd .gt. 0 ) call dnebnd_PB(enevcc(2,1:nfrg)) ! #SLC2 #PB
      if ( nbnd .gt. 0 ) call dnebnd_CB(enevcc(2,1:nfrg)) ! #SLC2 #CB
 
!     <<<  ANGLE ENERGY  >>>
      if ( nang .gt. 0 ) call dneang_PB(enevcc(3,1:nfrg)) ! #SLC2 #PB
      if ( nang .gt. 0 ) call dneang_CB(enevcc(3,1:nfrg)) ! #SLC2 #CB
 
!     <<<  1-4 INTERACTION ENERGY  >>>
!         TORSINAL ENERGY
!         IMPROPER TORSIONAL ENERGY
!         1-4 VAN DER WAALS ENERGY
!         1-4 ELECTROSTATIC ENERGY
      if ( iyntor.gt.0 .or. iynimp.gt.0 )                              &
        call dnei14_PB(enevcc(4,1:nfrg),enevcc(5,1:nfrg), & ! #SLC2 #PB
          enevcc(6,1:nfrg),enevcc(7,1:nfrg))                ! #SLC2 #PB
        call dnei14_CB(enevcc(4,1:nfrg),enevcc(5,1:nfrg), & ! #SLC2 #CB
          enevcc(6,1:nfrg),enevcc(7,1:nfrg))                ! #SLC2 #CB

!     <<<  1-5 INTERACTION ENERGY  >>>
!           VAN DER WAALS ENERGY
!           ELECTROSTATIC ENERGY
!           HYDROGEN BONDED ENERGY
      call cmmdo(enevcc(9,1:nfrg),enevcc(8,1:nfrg),enevcc(10,1:nfrg))

!********************************

      return
      end subroutine dnergv_PB ! #SLC2 #PB
      end subroutine dnergv_CB ! #SLC2 #CB


!=======================================================================


      subroutine dnebnd_PB(ebnd) ! #SLC2 #PB
      subroutine dnebnd_CB(ebnd) ! #SLC2 #CB

!***********************************************************************
!
!     CALCULATION OF BOND ENERGY AND GRADIENT
!
!***********************************************************************
 
      use COMBAS ; use COMERG ; use COMCMMC
      use COMCMM, only: nfrg,bndcls
      !$ use omp_lib

      implicit none

      ! Bond energy (kcal/mol)
        real(8),intent(inout):: ebnd(nfrg)

      real(8):: d(3),r12,dif,coef
      real(8):: work(3)
      integer(4):: ibnd,iatm(2),n,i
 
!******************************************

!     <<<  CALCULATION BOND ENERGY AND GRADIENT  >>>
      !$OMP parallel default(none)                                   & !
      !$OMP private(i,ibnd,n,iatm,d,r12,dif,coef,work)               & !
      !$OMP shared(np_bnd,para_bnd,bndcls,iypbnd,cord,fxcell,invcel, & !
      !$OMP        fyqbnd,fyfbnd,grad)                               & !
      !$OMP reduction (+ : ebnd)
      do i = 1,np_bnd
        !$OMP do
        do ibnd = para_bnd(1,i),para_bnd(2,i)
          n = bndcls(ibnd)
          ! Energy
          iatm(1:2) = iypbnd(1:2,ibnd)
          d(1:3) = cord(1:3,iatm(1)) - cord(1:3,iatm(2))
          d(1:3) = d(1:3) -                            & ! #SLC2 #PB
                   fxcell(1:3)*nint(d(1:3)*invcel(1:3))  ! #SLC2 #PB
          r12 = sqrt(d(1)*d(1) + d(2)*d(2) + d(3)*d(3))
          dif = r12 - fyqbnd(ibnd)
          coef = fyfbnd(ibnd) * dif
          ebnd(n) = ebnd(n) + (coef * dif)
          ! Gradient
          coef = 2.d0*coef / r12
          work(1:3) = coef*d(1:3)
          grad(1:3,iatm(1),n) = grad(1:3,iatm(1),n) + work(1:3)
          grad(1:3,iatm(2),n) = grad(1:3,iatm(2),n) - work(1:3)
        enddo
        !$OMP end do
      enddo
      !$OMP end parallel
 
!********************************

      return
      end subroutine dnebnd_PB ! #SLC2 #PB
      end subroutine dnebnd_CB ! #SLC2 #CB


!=======================================================================


      subroutine dneang_PB(eang) ! #SLC2 #PB
      subroutine dneang_CB(eang) ! #SLC2 #CB

!***********************************************************************
!
!     CALCULATION OF ANGLE ENERGY AND GRADIENT
!
!***********************************************************************
 
      use COMBAS ; use COMERG ; use PHYCNS ; use COMCMMC
      use COMCMM, only: nfrg,angcls
      !$ use omp_lib

      implicit none

      ! Angle energy (kcal/mol)
        real(8),intent(inout):: eang(nfrg)

      real(8):: d12(3),d32(3),r12,r32,r12r32,cost,dif,coef,theta,sint
      real(8):: work(6),rtmp
      integer(4):: i,iang,iatm(3),n
      integer(4),allocatable:: icng(:)

!***********************************************

!     <<<  CALCULATION OF ANGLE ENERGY AND GRADIENT  >>>
      !$OMP parallel default(none)                                   & !
      !$OMP private(i,iang,n,iatm,d12,d32,r12,r32,r12r32,cost,theta, & !
      !$OMP         dif,coef,sint,rtmp,work)                         & !
      !$OMP shared(np_ang,para_ang,angcls,iypang,cord,fxcell,invcel, & !
      !$OMP        fyqang,fyfang,grad)                               & !
      !$OMP reduction(+ : eang)
      do i = 1,np_ang
        !$OMP do
        do iang = para_ang(1,i),para_ang(2,i)
          n = angcls(iang)
          iatm(1:3) = iypang(1:3,iang)
          d12(1:3) = cord(1:3,iatm(1)) - cord(1:3,iatm(2))
          d32(1:3) = cord(1:3,iatm(3)) - cord(1:3,iatm(2))
          d12(1:3) = d12(1:3) -                            & ! #SLC2 #PB
                     fxcell(1:3)*anint(d12(1:3)*invcel(1:3)) ! #SLC2 #PB
          d32(1:3) = d32(1:3) -                            & ! #SLC2 #PB
                     fxcell(1:3)*anint(d32(1:3)*invcel(1:3)) ! #SLC2 #PB
          r12 = d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3)
          r32 = d32(1)*d32(1) + d32(2)*d32(2) + d32(3)*d32(3)
          r12r32 = sqrt(1.d0/(r12*r32))
          cost = d12(1)*d32(1) + d12(2)*d32(2) + d12(3)*d32(3)
          r12 = 1.d0 / r12 ; r32 = 1.d0 / r32
          cost = min(cost*r12r32,1.d0) ; cost = max(cost,-1.d0)
          theta = acos(cost)
          ! Energy
          dif = theta - fyqang(iang)
          coef = fyfang(iang) * dif
          eang(n) = eang(n) + coef*dif
          ! Gradient 
          !  SINT .EQ. 0.0D0  THEN TO CALCULATE GRADIENT MUST BE BY
          !  DIFFERENT FORM
          !  NOW , WE REGARD SMALL SINT AS CONSTANT VALUE EPS
          sint = max(sin(theta),eps)
!          if ( sint .lt. eps ) sint = eps
          coef = -2.d0*coef / sint
          rtmp = cost*r12
          work(1:3) = d32(1:3)*r12r32 - d12(1:3)*rtmp
          rtmp = cost*r32
          work(4:6) = d12(1:3)*r12r32 - d32(1:3)*rtmp
          work(1:6) = coef * work(1:6)
          grad(1:3,iatm(1),n) = grad(1:3,iatm(1),n) + work(1:3)
          grad(1:3,iatm(2),n) = grad(1:3,iatm(2),n)                    &
                                     - (work(1:3) + work(4:6))
          grad(1:3,iatm(3),n) = grad(1:3,iatm(3),n) + work(4:6)
        enddo
        !$OMP end do
      enddo
      !$OMP end parallel
 
!************************************
 
      return
      end subroutine dneang_PB ! #SLC2 #PB
      end subroutine dneang_CB ! #SLC2 #CB


!=======================================================================


      subroutine dnei14_PB(etorcn,eimpcn,evancn,eelecn) ! #SLC2 #PB
      subroutine dnei14_CB(etorcn,eimpcn,evancn,eelecn) ! #SLC2 #CB

!***********************************************************************
!
!     CALCULATION OF 1-4 INTERACTION ENERGY AND GRADIENT
!       TORSIONAL ENERGY
!       IMPROPER TORSIONAL ENERGY
!       1-4 VAN DER WAALS ENERGY
!       1-4 ELECTROSTATIC ENERGY
!
!***********************************************************************
 
      use COMBAS ; use COMERG ; use PHYCNS
      use COMCMM, only: nfrg,torcls,impcls

      implicit none

      ! Torsion, improper tors., 1-4 vdW, electrostatic energy (kcal/mol)
        real(8),intent(inout):: etorcn(nfrg),eimpcn(nfrg),evancn(nfrg),&
                                eelecn(nfrg)
 
!***********************************

!     <<<  TORSIONAL ENERGY AND GRADINET  >>>
      if ( ntor .gt. 0 ) then
        call dnetor_PB(ntor,np_tor,para_tor,               & ! #SLC2 #PB
          fyftor(1:ntor),fytrot(1:ntor),fytphs(1:ntor),    & ! #SLC2 #PB
          iyptor(1:4,1:ntor),iytdiv(1:ntor),etorcn,torcls)   ! #SLC2 #PB
        call dnetor_CB(ntor,np_tor,para_tor,               & ! #SLC2 #CB
          fyftor(1:ntor),fytrot(1:ntor),fytphs(1:ntor),    & ! #SLC2 #CB
          iyptor(1:4,1:ntor),iytdiv(1:ntor),etorcn,torcls)   ! #SLC2 #CB
      endif
 
!     <<<  IMPROPER TORSIONAL ENERGY AND GRADINET  >>>
      if ( nimp .gt. 0 ) then
        call dnetor_PB(nimp,np_imp,para_imp,fyfimp(1:nimp),& ! #SLC2 #PB
          fyirot(1:nimp),fyiphs(1:nimp),iypimp(1:4,1:nimp),& ! #SLC2 #PB
          iyidiv(1:nimp),eimpcn,impcls)                      ! #SLC2 #PB
        call dnetor_CB(nimp,np_imp,para_imp,fyfimp(1:nimp),& ! #SLC2 #CB
          fyirot(1:nimp),fyiphs(1:nimp),iypimp(1:4,1:nimp),& ! #SLC2 #CB
          iyidiv(1:nimp),eimpcn,impcls)                      ! #SLC2 #CB
      endif

!     <<<  1-4 VAN DER WAALS AND ELECTROSTATIC ENERGY  >>>
      ! Calc. vdW & electrostatic energy (constant dielectric)
      call dne14_PB(evancn,eelecn) ! #SLC2 #PB
      call dne14_CB(evancn,eelecn) ! #SLC2 #CB

!*************************************

      return
      end subroutine dnei14_PB ! #SLC2 #PB
      end subroutine dnei14_CB ! #SLC2 #CB
 

!=======================================================================


      subroutine dnetor_PB(ntor,np_tor,para_tor,fftor,frot,& ! #SLC2 #PB
                           fphs,iptor,idiv,etor,cls)         ! #SLC2 #PB
      subroutine dnetor_CB(ntor,np_tor,para_tor,fftor,frot,& ! #SLC2 #CB
                           fphs,iptor,idiv,etor,cls)         ! #SLC2 #CB

!***********************************************************************
!
!     CALCULATION OF TORSIONAL ENERGY AND GRADIENT
!
!***********************************************************************

      use COMBAS, only: ixnatm,fxcell,invcel
      use PHYCNS, only: pi,eps
      use COMCMM, only: nfrg
      use COMCMMC
      !$ use omp_lib

      implicit none

      ! Number of torsions
        integer(4),intent(in):: ntor
      ! Number of blocks for torsion calc. parallelization
        integer(4),intent(in):: np_tor
      ! Start and end number of Blocks
        integer(4),intent(in):: para_tor(2,108)
      ! Force constant of torsional energy (kcal/mol)
        real(8),intent(in):: fftor(ntor)
      ! Rotational symmetry number & Phase (radian)
        real(8),intent(in):: frot(ntor),fphs(ntor)
      ! Atom number of torsional pair & Division number
        integer(4),intent(in):: iptor(4,ntor),idiv(ntor)
      ! Torsion energy (kcal/mol)
        real(8),intent(inout):: etor(nfrg)
      ! Torsion cluster
        integer(4),intent(in):: cls(ntor)

      real(8):: d21(3),d32(3),d43(3),p12(3),p23(3),p123(3),r12,r23,    &
                r12r23,s1223,s32123,cosp,sinp,phi,phimod,div,work(9),  &
                wrk0,wrk1,wrk2,wrk3,rtmp
      integer(4):: itor,nrot,iatm(4),n,i

!************************************************

!     <<<  CALCULATION OF TORSIONAL ENERGY AND GRADIENT  >>>
      !$OMP parallel default(none)                                   & !
      !$OMP private(i,itor,iatm,n,d21,d32,d43,p12,p23,r12,r23,r12r23,& !
      !$OMP         s1223,cosp,p123,s32123,phi,phimod,div,sinp,nrot, & !
      !$OMP         wrk0,wrk1,wrk2,wrk3,work)                        & !
      !$OMP shared(np_tor,para_tor,fftor,iptor,cls,cord,fxcell,      & !
      !$OMP        invcel,frot,fphs,idiv,grad)                       & !
      !$OMP reduction(+ : etor)
      do i = 1,np_tor
        !$OMP do
        do itor = para_tor(1,i),para_tor(2,i)
          iatm(1:4) = iptor(1:4,itor)
          ! For cluster
          n = cls(itor)
          d21(1:3) = cord(1:3,iatm(2)) - cord(1:3,iatm(1))
          d32(1:3) = cord(1:3,iatm(3)) - cord(1:3,iatm(2))
          d43(1:3) = cord(1:3,iatm(4)) - cord(1:3,iatm(3))
          d21(1:3) = d21(1:3) -                            & ! #SLC2 #PB
                     fxcell(1:3)*nint(d21(1:3)*invcel(1:3))  ! #SLC2 #PB
          d32(1:3) = d32(1:3) -                            & ! #SLC2 #PB
                     fxcell(1:3)*nint(d32(1:3)*invcel(1:3))  ! #SLC2 #PB
          d43(1:3) = d43(1:3) -                            & ! #SLC2 #PB
                     fxcell(1:3)*nint(d43(1:3)*invcel(1:3))  ! #SLC2 #PB
          p12(1) = d21(2)*d32(3) - d32(2)*d21(3)
          p12(2) = d21(3)*d32(1) - d32(3)*d21(1)
          p12(3) = d21(1)*d32(2) - d32(1)*d21(2)
          p23(1) = d32(2)*d43(3) - d43(2)*d32(3)
          p23(2) = d32(3)*d43(1) - d43(3)*d32(1)
          p23(3) = d32(1)*d43(2) - d43(1)*d32(2)
          r12 = p12(1)*p12(1) + p12(2)*p12(2) + p12(3)*p12(3)
          r23 = p23(1)*p23(1) + p23(2)*p23(2) + p23(3)*p23(3)
          r12r23 = sqrt(r12*r23)

          ! If bond angle = 0 or pi, donot calculate torsion energy
          if ( r12r23 .gt. 0.d0 ) then
            r12r23 = 1.d0 / r12r23
            s1223 = p12(1)*p23(1) + p12(2)*p23(2) + p12(3)*p23(3)
            cosp = s1223 * r12r23
            if ( cosp .gt. 1.d0 ) then
              cosp = 1.d0
            elseif ( cosp .lt. -1.d0 ) then
              cosp = -1.d0
            endif
            p123(1) = p12(2)*p23(3) - p23(2)*p12(3)
            p123(2) = p12(3)*p23(1) - p23(3)*p12(1)
            p123(3) = p12(1)*p23(2) - p23(1)*p12(2)
            s32123 = d32(1)*p123(1) + d32(2)*p123(2) + d32(3)*p123(3)
            phi = acos(cosp)
            if ( s32123 .lt. 0.d0 ) phi = 2.d0*pi - phi

            ! Torsion energy
            phimod = frot(itor)*phi + fphs(itor)
            div = 1.d0 / dble(idiv(itor))
            etor(n) = etor(n) + fftor(itor)*(1.d0+cos(phimod))*div

            ! Gradient
            sinp = sin(phi) ; nrot = nint(frot(itor))
            select case (nrot)
            case (1)
              wrk0 = 1.d0
            case (2)
              wrk0 = 2.d0 * cosp
            case (3)
              wrk0 = -4.d0*sinp*sinp + 3.d0
            case (4)
              wrk0 = (-8.d0*sinp*sinp + 4.d0) * cosp
            case (5)
              wrk0 = (16.d0*sinp*sinp - 20.d0)*sinp*sinp + 5.d0
            case (6)
              wrk0 = cosp*(6.d0 + 32.d0*cosp*cosp*(cosp*cosp-1.d0))
            end select
            wrk0 = cos(fphs(itor))*frot(itor)*fftor(itor)*div*wrk0*    &
                   r12r23

            wrk1 = d21(1)+d32(1)
            wrk2 = d21(2)+d32(2)
            wrk3 = d21(3)+d32(3)
            r12 = s1223 / r12 ; r23 = s1223 / r23

            work(1) = d32(3)*p23(2)-d32(2)*p23(3) -                    &
                      r12*(d32(3)*p12(2)-d32(2)*p12(3))
            work(2) = d32(1)*p23(3)-d32(3)*p23(1) -                    &
                      r12*(d32(1)*p12(3)-d32(3)*p12(1))
            work(3) = d32(2)*p23(1)-d32(1)*p23(2) -                    &
                      r12*(d32(2)*p12(1)-d32(1)*p12(2))
            work(4) = -wrk3*p23(2)+d43(3)*p12(2) +                     &
                       wrk2*p23(3)-d43(2)*p12(3) -                     &
                       r12*(-wrk3*p12(2)+wrk2*p12(3)) -                &
                       r23*(d43(3)*p23(2)-d43(2)*p23(3))
            work(5) = wrk3*p23(1)-d43(3)*p12(1) -                      &
                      wrk1*p23(3)+d43(1)*p12(3) -                      &
                      r12*( wrk3*p12(1)-wrk1*p12(3)) -                 &
                      r23*(d43(1)*p23(3)-d43(3)*p23(1))
            work(6) = -wrk2*p23(1)+d43(2)*p12(1) +                     &
                       wrk1*p23(2)-d43(1)*p12(2) -                     &
                       r12*(-wrk2*p12(1)+wrk1*p12(2)) -                &
                       r23*(d43(2)*p23(1)-d43(1)*p23(2))
            work(7) = d32(3)*p12(2)-d32(2)*p12(3) -                    &
                      r23 * (d32(3)*p23(2)-d32(2)*p23(3))
            work(8) = d32(1)*p12(3)-d32(3)*p12(1) -                    &
                      r23*(d32(1)*p23(3)-d32(3)*p23(1))
            work(9) = d32(2)*p12(1)-d32(1)*p12(2) -                    &
                      r23*(d32(2)*p23(1)-d32(1)*p23(2))
            work(1:9) = work(1:9) * wrk0
            grad(1:3,iatm(1),n) = grad(1:3,iatm(1),n) + work(1:3)
            grad(1:3,iatm(2),n) = grad(1:3,iatm(2),n) + work(4:6)
            grad(1:3,iatm(3),n) = grad(1:3,iatm(3),n) -                &
                              (work(1:3)+work(4:6)+work(7:9))
            grad(1:3,iatm(4),n) = grad(1:3,iatm(4),n) + work(7:9)
          endif
        enddo
        !$OMP end do
      enddo
      !$OMP end parallel

!*****************************************

      return
      end subroutine dnetor_PB ! #SLC2 #PB
      end subroutine dnetor_CB ! #SLC2 #CB


!=======================================================================


      subroutine dne14_PB(evan,eele) ! #SLC2 #PB
      subroutine dne14_CB(evan,eele) ! #SLC2 #CB

!***********************************************************************
!
!     CALCULATION OF 1-4 VAN DER WAALS AND ELECTROSTATIC
!     ENERGY (CONSTANT DIELECTRIC) AND GRADIENT
!
!***********************************************************************
 
      use COMBAS ; use COMERG ; use PHYCNS
      use COMCMM, only: nfrg
      use COMCMMC
      !$ use omp_lib

      implicit none

      ! 1-4 vdW & electrostatic energy (kcal/mol)
        real(8),intent(inout):: evan(nfrg),eele(nfrg)

      real(8):: d14(3),r14,r141,r146,r1412,co6,co12,cc,coef,work(3)
      integer(4):: itor,iatm1,iatm2,i,j,n,ii
 
!************************************************

!     <<<  CALCULATION OF ENERGY AND GRADIENT  >>>
      !$OMP parallel default(none)                                   & !
      !$OMP private(ii,itor,iatm1,iatm2,n,d14,r14,r141,r146,r1412,i, & !
      !$OMP         j,co6,co12,coef,work,cc)                         & !
      !$OMP shared(np_14,para_14,iyp14,i14cls,cord,fxcell,invcel,    & !
      !$OMP        ixatyp,fynbpp,fytvws14,fxchrg,ccoef,fytess14,grad)& !
      !$OMP reduction (+ : evan,eele)
      do ii = 1,np_14
        !$OMP do
        do itor = para_14(1,ii),para_14(2,ii)
          iatm1 = iyp14(1,itor) ; iatm2 = iyp14(2,itor)
          ! For cluster
          n = i14cls(itor)

          d14(1:3) = cord(1:3,iatm1) - cord(1:3,iatm2)
          d14(1:3) = d14(1:3) -                           & ! #SLC2 #PB
                     fxcell(1:3)*nint(d14(1:3)*invcel(1:3)) ! #SLC2 #PB
          r14 = d14(1)*d14(1) + d14(2)*d14(2) + d14(3)*d14(3)
          r14 = 1.d0 / r14
          r141 = sqrt(r14)

          ! 1-4 vdW energy
          r146 = r14 * r14 * r14
          r1412 = r146 * r146
          i = ixatyp(iatm1) ; j = ixatyp(iatm2)
          co6 = fynbpp(1,i,j)*r146
          co12 = fynbpp(2,i,j)*r1412
          evan(n) = evan(n) + fytvws14(itor)*(-co6+co12)
          ! 1-4 vdW gradient
          coef = r14*fytvws14(itor)*(-2.d0*co12+co6)*6.d0
          work(1:3) = coef * d14(1:3)

          ! 1-4 electrostatic energy (constant dielectric)
          !   fytess14 = fxchrg(iatm1)*fxchrg(iatm2)*ccoef*damping_factor
          coef = fytess14(itor) * r141
          eele(n) = eele(n) + coef
          ! 1-4 electrostatic gradient
          coef = coef * (-r14)
          work(1:3) = work(1:3) + coef * d14(1:3)

          grad(1:3,iatm1,n) = grad(1:3,iatm1,n) + work(1:3)
          grad(1:3,iatm2,n) = grad(1:3,iatm2,n) - work(1:3)
        enddo
        !$OMP end do
      enddo
      !$OMP end parallel

!*************************************

      return
      end subroutine dne14_PB ! #SLC2 #PB
      end subroutine dne14_CB ! #SLC2 #CB
