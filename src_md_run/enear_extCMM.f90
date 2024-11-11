
      subroutine enear_extCMM(eelecn)

!*****************************************************************************
!
!     Calculate the near-field exact interactions.
!     The deepest level = 3.
!
!     The vdW interactions are the sum  of the following two atpm pairs:
!       (1) the atom pairs defined by CMM and
!       (2) the atom pairs with shorter pair distances.
!
!      The 1-2, 1-2, & 1-4 interactions are eliminated in
!      the result.  For this elimination, these interactions
!      are calculated.
!
!*****************************************************************************

      use COMBAS ; use COMERG ; use COMCMM ; use COMCMMC
!!      use CALC_TIME
      !$ use omp_lib

      implicit none

!********

      real(8),intent(inout):: eelecn(nfrg)

      integer(4):: i,j,k,n,ii,kk,it1,it2
      real(8):: xi,yi,zi,dx,dy,dz,dd2i,dd3i,dd6i,dd2
      real(8):: v1,v2,cc,ee,a(3)

      real(8):: tENE(nfrg),tENE2,tFOR(1:3)
      integer(4):: num1(1:nfrg)
      integer(4):: ix,iy,iz,jx,jy,jz,iat
      real(8):: q(13),dx2,dy2,dz2,scal,ddi,dd5i,dxy,dyz,dzx,dd5i3,     &
                dd7i5,aa1,aa2,bbx,bby,bbz,cc2x,cc2y,cc2z,cc4,sz_FOR2,  &
                scal_FOR,dd7i,r_ext,rtmpB,rtmpC,rtmpC__2,rtmpC__3

!**************************************
!     nearest cell calc. for external CMM

!!      call system_clock(tim1)
      tENE(1:nfrg) = 0.d0
      ! Select a cell around the boundary
      !$OMP parallel default(none)                                   & !
      !$OMP private(i,ix,iy,iz,num1,j,jx,jy,jz,q,n,k,iat,dx,dy,dz,   & !
      !$OMP         dx2,dy2,dz2,dd2,dd2i,scal,ddi,dd3i,dd5i,dxy,dyz, & !
      !$OMP         dzx,dd5i3,dd7i5,aa1,aa2,tENE2,bbx,bby,bbz,cc2x,  & !
      !$OMP         cc2y,cc2z,cc4,tFOR,sz_FOR2,scal_FOR,dd7i,ee,     & !
      !$OMP         r_ext,rtmpB,rtmpC,rtmpC__2,rtmpC__3)             & !
      !$OMP shared(Nmycel_nEC,mycel_nEC,npcl1,nfrg,Nyrcel_nEC,       & !
      !$OMP        yrcel_nEC,eCMMpot1,ipcl,cord,c1,RR_ext__2,RR_ext,& !
      !$OMP        iRR_ext,iRR_ext__3,iRR_ext__5,i3RR_ext__5,        & !
      !$OMP        RR_ext_2,i5RR_ext__7,tENE,chgmod2,grad)
      !$OMP do reduction(+ : tENE)                                   & !
      !$OMP    schedule (guided)
      do i = 1,Nmycel_nEC
        ix = mycel_nEC(1,i) ; iy = mycel_nEC(2,i) ; iz = mycel_nEC(3,i)
        num1(1:nfrg) = npcl1(1:nfrg,ix,iy,iz)
        ! Select a nearest cell that its exflag1 is .true.
        do j = 1,Nyrcel_nEC(i)
          jx=yrcel_nEC(1,j,i); jy=yrcel_nEC(2,j,i); jz=yrcel_nEC(3,j,i)

          ! Get parameters for external CMM
          q(1:13) = eCMMpot1(1:13,jx,jy,jz)

          ! Calc. external CMM between an atom in the target cell &
          !       the surrounding (or same) cells
          do n = 1,nfrg
            do k = 1,num1(n)
              iat = ipcl(k,n,ix,iy,iz)
              dx = cord(1,iat) - c1(1,jx)
              dy = cord(2,iat) - c1(2,jy)
              dz = cord(3,iat) - c1(3,jz)
              dx2 = dx*dx ; dy2 = dy*dy ; dz2 = dz*dz
              dd2 = dx2 + dy2 + dz2

              ! E & F are scaled, if dd2 < RR_ext__2
              if ( dd2 .lt. RR_ext__2 ) then
                dd2i = 1.d0 / dd2
                scal = sqrt(dd2i*RR_ext__2)
                dx = scal*dx ; dy = scal*dy ; dz = scal*dz
                dx2 = dx*dx ; dy2 = dy*dy ; dz2 = dz*dz

                ddi = iRR_ext
                dd3i = iRR_ext__3
                dd5i = iRR_ext__5
                dxy = dx*dy ; dyz = dy*dz ; dzx = dz*dx
                dd5i3 = i3RR_ext__5 ; dd7i5 = i5RR_ext__7

                !** Energy when dd2 = RR_ext__2
                aa1 = q(2)*dx + q(3)*dy + q(4)*dz
                aa2 = q(5)*dx2 + q(6)*dy2 + q(7)*dz2 +                 &
                      q(11)*dxy + q(12)*dyz + q(13)*dzx
                tENE2 = q(1)*ddi + aa1*dd3i + aa2*dd5i
                !** Force when dd2 = RR_ext__2
                bbx = (q(1)*dx - q(2)) * dd3i
                bby = (q(1)*dy - q(3)) * dd3i
                bbz = (q(1)*dz - q(4)) * dd3i
                cc2x = (q(8)*dx + q(11)*dy + q(13)*dz)*dd5i
                cc2y = (q(11)*dx + q(9)*dy + q(12)*dz)*dd5i
                cc2z = (q(13)*dx + q(12)*dy + q(10)*dz)*dd5i
                cc4 = aa1*dd5i3 + aa2*dd7i5
                tFOR(1) = -bbx + cc2x - cc4*dx
                tFOR(2) = -bby + cc2y - cc4*dy
                tFOR(3) = -bbz + cc2z - cc4*dz
                sz_FOR2 = tFOR(1)*tFOR(1) + tFOR(2)*tFOR(2) +          &
                          tFOR(3)*tFOR(3)
                sz_FOR2 = sqrt(sz_FOR2)

                !** Scale Force
                r_ext = sqrt(dd2)
                rtmpB = sz_FOR2*r_ext
                scal_FOR = r_ext*iRR_ext__3 *                          &
                           ( 6.d0*tENE2*(r_ext-RR_ext) +               &
                             rtmpB*(3.d0*r_ext-RR_ext_2))
                scal_FOR = scal_FOR / sz_FOR2
                tFOR(1:3) = tFOR(1:3) * scal_FOR * chgmod2(iat)
                !** Scale Energy
                rtmpC = r_ext * iRR_ext
                rtmpC__2 = rtmpC * rtmpC ; rtmpC__3 = rtmpC__2 * rtmpC
                tENE2 = (2.d0*tENE2 + rtmpB)*rtmpC__3 -                &
                        (3.d0*tENE2 + rtmpB)*rtmpC__2 + 2.d0 * tENE2

                tENE(n) = tENE(n) + tENE2 * chgmod2(iat)
                grad(1:3,iat,n) = grad(1:3,iat,n) + tFOR(1:3)

              ! E & F are calculated normally, if RR_ext__2 <= dd2
              else
                dd2i = 1.d0 / dd2
                ddi = sqrt(dd2i)
                dd3i = ddi*dd2i
                dd5i = dd3i*dd2i
                dd7i = dd5i*dd2i
                dxy = dx*dy ; dyz = dy*dz ; dzx = dz*dx
                dd5i3 = 3*dd5i ; dd7i5 = 5*dd7i

                !** Energy
                aa1 = q(2)*dx + q(3)*dy + q(4)*dz
                aa2 = q(5)*dx2 + q(6)*dy2 + q(7)*dz2 +                 &
                      q(11)*dxy + q(12)*dyz + q(13)*dzx
                ee = q(1)*ddi + aa1*dd3i + aa2*dd5i
                tENE(n) = tENE(n) + chgmod2(iat)*ee
                !** Force
                bbx = (q(1)*dx - q(2)) * dd3i
                bby = (q(1)*dy - q(3)) * dd3i
                bbz = (q(1)*dz - q(4)) * dd3i
                cc2x = (q(8)*dx + q(11)*dy + q(13)*dz)*dd5i
                cc2y = (q(11)*dx + q(9)*dy + q(12)*dz)*dd5i
                cc2z = (q(13)*dx + q(12)*dy + q(10)*dz)*dd5i
                cc4 = aa1*dd5i3 + aa2*dd7i5
                tFOR(1) = -bbx + cc2x - cc4*dx
                tFOR(2) = -bby + cc2y - cc4*dy
                tFOR(3) = -bbz + cc2z - cc4*dz
                tFOR(1:3) = tFOR(1:3) * chgmod2(iat)
                grad(1:3,iat,n) = grad(1:3,iat,n) + tFOR(1:3)
              endif
            enddo
          enddo
        enddo
      enddo
      !$OMP end do
      !$OMP end parallel

      tENE(1:nfrg) = tENE(1:nfrg) * 0.5d0
      eelecn(1:nfrg) = eelecn(1:nfrg) + tENE(1:nfrg)
!!      call system_clock(tim2)
!!      tmptime = tmptime + tim2 - tim1

!************************************************************

      return
      end subroutine enear_extCMM
