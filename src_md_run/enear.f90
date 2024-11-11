
      subroutine enear_CMM_redu(eelecn,evancn,ehydcn)    ! #SLCT #CMM #REDU
      subroutine enear_CMM_high(eelecn,evancn,ehydcn)    ! #SLCT #CMM #HIGH
      subroutine enear_CMM_doub(eelecn,evancn,ehydcn)    ! #SLCT #CMM #DOUB
      subroutine enear_ACT_redu_PB(eelecn,evancn,ehydcn) ! #SLCT #ACT #REDU #PB
      subroutine enear_ACT_high_PB(eelecn,evancn,ehydcn) ! #SLCT #ACT #HIGH #PB
      subroutine enear_ACT_doub_PB(eelecn,evancn,ehydcn) ! #SLCT #ACT #DOUB #PB
      subroutine enear_ACT_redu_CB(eelecn,evancn,ehydcn) ! #SLCT #ACT #REDU #CB
      subroutine enear_ACT_high_CB(eelecn,evancn,ehydcn) ! #SLCT #ACT #HIGH #CB
      subroutine enear_ACT_doub_CB(eelecn,evancn,ehydcn) ! #SLCT #ACT #DOUB #CB
      subroutine enear_ZD_redu_PB(eelecn,evancn,ehydcn)  ! #SLCT #ZD  #REDU #PB
      subroutine enear_ZD_high_PB(eelecn,evancn,ehydcn)  ! #SLCT #ZD  #HIGH #PB
      subroutine enear_ZD_doub_PB(eelecn,evancn,ehydcn)  ! #SLCT #ZD  #DOUB #PB
      subroutine enear_ZD_redu_CB(eelecn,evancn,ehydcn)  ! #SLCT #ZD  #REDU #CB
      subroutine enear_ZD_high_CB(eelecn,evancn,ehydcn)  ! #SLCT #ZD  #HIGH #CB
      subroutine enear_ZD_doub_CB(eelecn,evancn,ehydcn)  ! #SLCT #ZD  #DOUB #CB

!*****************************************************************************
!
!     Calculate the near-field exact interactions.
!
!     The vdW interactions are the sum of the following two atpm pairs:
!       (1) the atom pairs defined by CMM and
!       (2) the atom pairs with shorter pair distances.
!
!*****************************************************************************

      use COMBAS ; use COMERG ; use COMCMM ; use COMCMMC
      use CALC_TIME
      !$ use omp_lib

      implicit none

!********

      real(8),intent(inout):: evancn(nfrg),eelecn(nfrg),ehydcn(nfrg)

      integer(4):: i,j,ii,kk,it1,it2,itmp
      integer(4):: k,kj     ! #SLC2 #HIGH
      real(8):: xi,yi,zi,dx,dy,dz,dd2i,dd6i,v1,v2,cc,ee,fx,fy,fz,chgi, &
        gwx,gwy,gwz,ee_sum,v_sum
      real(8):: dd2,ddi     ! #SLC2 #ZD
      real(8):: co          ! #SLC2 #ZD #ACT

      ! For vectrization
      integer(4):: j_atoms(i_vec)
      real(8):: xv(i_vec),yv(i_vec),zv(i_vec),chgv(i_vec),vdw6v(i_vec),&
        vdw12v(i_vec),dxv(i_vec),dyv(i_vec),dzv(i_vec),dd2iv(i_vec),   &
        dd6iv(i_vec),eev(i_vec),ccv(i_vec),fxv(i_vec),fyv(i_vec),      &
        fzv(i_vec),v1v(i_vec),v2v(i_vec),gwxv(i_vec),gwyv(i_vec),      &
        gwzv(i_vec),eev_sum(i_vec),vv_sum(i_vec)
      real(8):: dd2v(i_vec),ddiv(i_vec),Eexc(nfrg)  ! #SLC2 #ZD
      real(8):: cov(i_vec)                          ! #SLC2 #ZD #ACT
      real(8),allocatable:: Tgrad(:,:,:)

!**************************************
      ! The initialization is done before calling this subroutine.

!      call system_clock(tim1)
      eelecn(1:nfrg) = 0.d0 ; evancn(1:nfrg) = 0.d0
      allocate(Tgrad(3,ixnatm,nfrg))
      !$OMP parallel default (none)                                  & !
      !$OMP private(i,j)                                             & !
      !$OMP shared(nfrg,ixnatm,Tgrad,Tcord,cord,inum)
      !$OMP do schedule (static) collapse(2)
      do j = 1,nfrg
      do i = 1,ixnatm
        Tgrad(:,i,j) = 0.d0
      enddo
      enddo
      !$OMP end do
      !$OMP do schedule (static)
      do i = 1,ixnatm
        Tcord(1:3,inum(i)) = cord(1:3,i)
      enddo
      !$OMP end do
      !$OMP end parallel

!******************
!  Near field calculation.

      !$OMP parallel default(none)                                   & !
      !$OMP private(ii,it1,xi,yi,zi,i,itmp,j,kk,dx,dy,dz,dd2i,dd6i,  & !
      !$OMP         ee,it2,v1,v2,cc,fx,fy,fz,ee_sum,v_sum,chgi,gwx,  & !
      !$OMP         gwy,gwz,xv,yv,zv,chgv,gwxv,gwyv,gwzv,vdw6v,      & !
      !$OMP         vdw12v,dxv,dyv,dzv,dd2iv,dd6iv,eev,v1v,v2v,ccv,  & !
      !$OMP         fxv,fyv,fzv,eev_sum,vv_sum,j_atoms               & !
      !$OMP         )                                & ! ! #SLCT #CMM #REDU #DOUB
      !$OMP         ,k,kj)                           & ! ! #SLCT #CMM #HIGH
      !$OMP         ,dd2,ddi,dd2v,ddiv,cov,co)       & ! ! #SLCT #ZD #REDU #DOUB #PB #CB
      !$OMP         ,k,kj,dd2,ddi,dd2v,ddiv,cov,co)  & ! ! #SLCT #ZD #HIGH #PB #CB
      !$OMP         ,cov,co)                         & ! ! #SLCT #ACT #REDU #DOUB #PB #CB
      !$OMP         ,k,kj,cov,co)                    & ! ! #SLCT #ACT #HIGH #PB #CB
      !$OMP shared(ixnatm,Tcord,Tixatyp,nfrg,ntb,itb,fynbpp,ntbEL,   & !
      !$OMP        itbEL,i_vec                                       & !
      !$OMP        ,ntbvdW,itbvdW                    & ! ! #SLC2 #CMM
      !$OMP        ,fxcell,invcel                    & ! ! #SLC2 #PB
      !$OMP        ,irlim2                           & ! ! #SLC2 #ZD #ACT
      !$OMP        )                                 & ! ! #SLCT #CMM #ACT #REDU #PB #CB
      !$OMP        ,high_para,np_cell,para_se,Tgrad) & ! ! #SLCT #CMM #ACT #HIGH #PB #CB
      !$OMP        ,Tgrad)                           & ! ! #SLCT #CMM #ACT #DOUB #PB #CB
      !$OMP        ,high_para,np_cell,para_se,Tgrad  & ! ! #SLCT #ZD  #HIGH #PB #CB
      !$OMP        ,Tgrad                            & ! ! #SLCT #ZD  #DOUB #PB #CB
      !$OMP        ,bcoeff,fcoeff,zcore)             & ! ! #SLC2 #ZD
      !$OMP reduction(+ : eelecn,evancn,Tgrad)           ! #SLC2 #REDU
      !$OMP reduction(+ : eelecn,evancn)                 ! #SLC2 #HIGH #DOUB
      do k = 1,high_para                                 ! #SLC2 #HIGH
      do i = nfrg,1,-1
        ee_sum = 0.d0 ; v_sum = 0.d0 ; eev_sum = 0.d0 ; vv_sum = 0.d0
        !$OMP do schedule (static)                   ! #SLC2 #REDU #DOUB
        do ii = 1,ixnatm                             ! #SLC2 #REDU #DOUB
        !$OMP do schedule(dynamic)                   ! #SLC2 #HIGH
        do kj = 1,np_cell(k)                         ! #SLC2 #HIGH
        do ii = para_se(1,kj,k),para_se(2,kj,k)      ! #SLC2 #HIGH
          it1 = Tixatyp(ii)
          gwx = 0.d0 ; gwy = 0.d0 ; gwz = 0.d0
          gwxv = 0.d0 ; gwyv = 0.d0 ; gwzv = 0.d0
          xi = Tcord(1,ii) ; yi = Tcord(2,ii) ; zi = Tcord(3,ii)
          chgi = Tcord(4,ii)

          itmp = ntbEL(ii,i) - mod(ntbEL(ii,i),i_vec)
          ! electrostatic, hand-vectorized
          do j = 1,itmp,i_vec
            j_atoms = itbEL(j:j+i_vec-1,ii,i)
            xv=Tcord(1,j_atoms); yv=Tcord(2,j_atoms);zv=Tcord(3,j_atoms)
            chgv = Tcord(4,j_atoms)
            dxv = xi - xv ; dyv = yi - yv ; dzv = zi - zv
            dxv = dxv - fxcell(1) * nint(dxv*invcel(1))      ! #SLC2 #PB
            dyv = dyv - fxcell(2) * nint(dyv*invcel(2))      ! #SLC2 #PB
            dzv = dzv - fxcell(3) * nint(dzv*invcel(3))      ! #SLC2 #PB
            dd2iv = 1.d0 / (dxv*dxv + dyv*dyv + dzv*dzv)     ! #SLC2 #CMM #ACT
            dd2v = dxv*dxv + dyv*dyv + dzv*dzv               ! #SLC2 #ZD
            dd2iv = 1.d0 / dd2v                              ! #SLC2 #ZD
            cov = sign(0.5d0,dd2iv-irlim2) + 0.5d0           ! #SLC2 #ZD #ACT
            ddiv = sqrt(dd2iv)                               ! #SLC2 #ZD

            eev = chgi * chgv * sqrt(dd2iv)                  ! #SLC2 #CMM
            eev = chgi * chgv * cov                          ! #SLC2 #ZD
            eev = chgi * chgv * sqrt(dd2iv) * cov            ! #SLC2 #ACT
            eev_sum = eev_sum + eev                          ! #SLC2 #CMM #ACT
            eev_sum = eev_sum + eev*(ddiv+bcoeff*dd2v-zcore) ! #SLC2 #ZD

            ccv = eev*dd2iv                                  ! #SLC2 #CMM #ACT
            ccv = eev*(dd2iv*ddiv-fcoeff)                    ! #SLC2 #ZD
            fxv = dxv*ccv ; fyv = dyv*ccv ; fzv = dzv*ccv

            gwxv = gwxv+fxv ; gwyv = gwyv+fyv ; gwzv = gwzv+fzv 
            Tgrad(1,j_atoms,i) = Tgrad(1,j_atoms,i) + fxv    ! #SLC2 #REDU #HIGH
            Tgrad(2,j_atoms,i) = Tgrad(2,j_atoms,i) + fyv    ! #SLC2 #REDU #HIGH
            Tgrad(3,j_atoms,i) = Tgrad(3,j_atoms,i) + fzv    ! #SLC2 #REDU #HIGH
          enddo

          ! electrostatic, remain
          do j = itmp+1,ntbEL(ii,i)
            kk = itbEL(j,ii,i)
            dx=xi-Tcord(1,kk) ; dy=yi-Tcord(2,kk) ; dz=zi-Tcord(3,kk)
            dx = dx - fxcell(1) * nint(dx*invcel(1))     ! #SLC2 #PB
            dy = dy - fxcell(2) * nint(dy*invcel(2))     ! #SLC2 #PB
            dz = dz - fxcell(3) * nint(dz*invcel(3))     ! #SLC2 #PB
            dd2i = 1.d0 / (dx*dx + dy*dy + dz*dz)        ! #SLC2 #CMM #ACT
            dd2 = dx*dx + dy*dy + dz*dz                  ! #SLC2 #ZD
            dd2i = 1.d0 / dd2                            ! #SLC2 #ZD
            co = sign(0.5d0,dd2i-irlim2) + 0.5d0         ! #SLC2 #ZD #ACT
            ddi = sqrt(dd2i)                             ! #SLC2 #ZD

            ee = chgi * Tcord(4,kk) * sqrt(dd2i)         ! #SLC2 #CMM
            ee = chgi * Tcord(4,kk) * co                 ! #SLC2 #ZD
            ee = chgi * Tcord(4,kk) * sqrt(dd2i) * co    ! #SLC2 #ACT
            ee_sum = ee_sum + ee                         ! #SLC2 #CMM #ACT
            ee_sum = ee_sum + ee*(ddi+bcoeff*dd2-zcore)  ! #SLC2 #ZD

            cc = ee*dd2i                                 ! #SLC2 #CMM #ACT
            cc = ee*(dd2i*ddi-fcoeff)                    ! #SLC2 #ZD
            fx = dx*cc ; fy = dy*cc ; fz = dz*cc

            gwx = gwx + fx ; gwy = gwy + fy ; gwz = gwz + fz
            Tgrad(1,kk,i) = Tgrad(1,kk,i) + fx           ! #SLC2 #REDU #HIGH
            Tgrad(2,kk,i) = Tgrad(2,kk,i) + fy           ! #SLC2 #REDU #HIGH
            Tgrad(3,kk,i) = Tgrad(3,kk,i) + fz           ! #SLC2 #REDU #HIGH
          enddo

          itmp = ntb(ii,i) - mod(ntb(ii,i),i_vec)
          ! vdW & electrostatic, hand-vectrized
          do j = 1,itmp,i_vec
            j_atoms = itb(j:j+i_vec-1,ii,i)
            xv=Tcord(1,j_atoms); yv=Tcord(2,j_atoms);zv=Tcord(3,j_atoms)
            chgv = Tcord(4,j_atoms)
            dxv = xi - xv ; dyv = yi - yv ; dzv = zi - zv
            dxv = dxv - fxcell(1) * nint(dxv*invcel(1))      ! #SLC2 #PB
            dyv = dyv - fxcell(2) * nint(dyv*invcel(2))      ! #SLC2 #PB
            dzv = dzv - fxcell(3) * nint(dzv*invcel(3))      ! #SLC2 #PB
            vdw6v = fynbpp(1,Tixatyp(j_atoms),it1)
            vdw12v = fynbpp(2,Tixatyp(j_atoms),it1)
            dd2iv = 1.d0 / (dxv*dxv + dyv*dyv + dzv*dzv)     ! #SLC2 #CMM #ACT
            dd2v = dxv*dxv + dyv*dyv + dzv*dzv               ! #SLC2 #ZD
            dd2iv = 1.d0 / dd2v                              ! #SLC2 #ZD
            cov = sign(0.5d0,dd2iv-irlim2) + 0.5d0           ! #SLC2 #ZD #ACT
            ddiv = sqrt(dd2iv)                               ! #SLC2 #ZD
            dd6iv = dd2iv * dd2iv * dd2iv                    ! #SLC2 #CMM
            dd6iv = dd2iv * dd2iv * dd2iv * cov              ! #SLC2 #ZD #ACT

            eev = chgi * chgv * sqrt(dd2iv)                  ! #SLC2 #CMM
            eev = chgi * chgv * cov                          ! #SLC2 #ZD
            eev = chgi * chgv * sqrt(dd2iv) * cov            ! #SLC2 #ACT
            eev_sum = eev_sum + eev                          ! #SLC2 #CMM #ACT
            eev_sum = eev_sum + eev*(ddiv+bcoeff*dd2v-zcore) ! #SLC2 #ZD

            v1v = vdw6v * dd6iv
            v2v = vdw12v * dd6iv * dd6iv
            vv_sum = vv_sum - v1v + v2v

            ccv = (eev -6.d0*v1v + 12.d0*v2v) * dd2iv        ! #SLC2 #CMM #ACT
            ccv = eev*(dd2iv*ddiv-fcoeff) +                & ! #SLC2 #ZD
                  (-6.d0*v1v+12.d0*v2v)*dd2iv                ! #SLC2 #ZD
            fxv = dxv*ccv ; fyv = dyv*ccv ; fzv = dzv*ccv

            gwxv = gwxv+fxv ; gwyv = gwyv+fyv ; gwzv = gwzv+fzv 
            Tgrad(1,j_atoms,i) = Tgrad(1,j_atoms,i) + fxv    ! #SLC2 #REDU #HIGH
            Tgrad(2,j_atoms,i) = Tgrad(2,j_atoms,i) + fyv    ! #SLC2 #REDU #HIGH
            Tgrad(3,j_atoms,i) = Tgrad(3,j_atoms,i) + fzv    ! #SLC2 #REDU #HIGH
          enddo

          ! vdW & electrostatic, remain
          do j = itmp+1,ntb(ii,i)
            kk = itb(j,ii,i)
            dx=xi-Tcord(1,kk) ; dy=yi-Tcord(2,kk) ; dz=zi-Tcord(3,kk)
            dx = dx - fxcell(1) * nint(dx*invcel(1))    ! #SLC2 #PB
            dy = dy - fxcell(2) * nint(dy*invcel(2))    ! #SLC2 #PB
            dz = dz - fxcell(3) * nint(dz*invcel(3))    ! #SLC2 #PB
            dd2i = 1.d0 / (dx*dx + dy*dy + dz*dz)       ! #SLC2 #CMM #ACT
            dd2 = dx*dx + dy*dy + dz*dz                 ! #SLC2 #ZD
            dd2i = 1.d0 / dd2                           ! #SLC2 #ZD
            co = sign(0.5d0,dd2i-irlim2) + 0.5d0        ! #SLC2 #ZD #ACT
            ddi = sqrt(dd2i)                            ! #SLC2 #ZD
            dd6i = dd2i * dd2i * dd2i                   ! #SLC2 #CMM
            dd6i = dd2i * dd2i * dd2i * co              ! #SLC2 #ZD #ACT

            ee = chgi * Tcord(4,kk) * sqrt(dd2i)        ! #SLC2 #CMM
            ee = chgi * Tcord(4,kk) * co                ! #SLC2 #ZD
            ee = chgi * Tcord(4,kk) * sqrt(dd2i) * co   ! #SLC2 #ACT
            ee_sum = ee_sum + ee                        ! #SLC2 #CMM #ACT
            ee_sum = ee_sum + ee*(ddi+bcoeff*dd2-zcore) ! #SLC2 #ZD

            it2 = Tixatyp(kk)
            v1 = fynbpp(1,it2,it1) * dd6i
            v2 = fynbpp(2,it2,it1) * dd6i * dd6i
            v_sum = v_sum - v1 + v2 

            cc = (ee -6.d0*v1 + 12.d0*v2) * dd2i        ! #SLC2 #CMM #ACT
            cc = ee*(dd2i*ddi-fcoeff) +               & ! #SLC2 #ZD
                 (-6.d0*v1+12.d0*v2)*dd2i               ! #SLC2 #ZD
            fx = dx*cc ; fy = dy*cc ; fz = dz*cc

            gwx = gwx + fx ; gwy = gwy + fy ; gwz = gwz + fz
            Tgrad(1,kk,i) = Tgrad(1,kk,i) + fx          ! #SLC2 #REDU #HIGH
            Tgrad(2,kk,i) = Tgrad(2,kk,i) + fy          ! #SLC2 #REDU #HIGH
            Tgrad(3,kk,i) = Tgrad(3,kk,i) + fz          ! #SLC2 #REDU #HIGH
          enddo

          !  The added vdW interactions for pairs for that  ! #SLC2 #CMM
          !  the pair-distances are shorter than a limit.   ! #SLC2 #CMM
          itmp = ntbvdW(ii,i) - mod(ntbvdW(ii,i),i_vec)     ! #SLC2 #CMM
          ! vdW, hand-vectrized                             ! #SLC2 #CMM
          do j = 1,itmp,i_vec                               ! #SLC2 #CMM
            j_atoms = itbvdW(j:j+i_vec-1,ii,i)              ! #SLC2 #CMM
            xv = Tcord(1,j_atoms) ; yv = Tcord(2,j_atoms)   ! #SLC2 #CMM
            zv = Tcord(3,j_atoms)                           ! #SLC2 #CMM
            vdw6v = fynbpp(1,Tixatyp(j_atoms),it1)          ! #SLC2 #CMM
            vdw12v = fynbpp(2,Tixatyp(j_atoms),it1)         ! #SLC2 #CMM
            dxv = xi - xv ; dyv = yi - yv ; dzv = zi - zv   ! #SLC2 #CMM
            dd2iv = 1.d0 / (dxv*dxv + dyv*dyv + dzv*dzv)    ! #SLC2 #CMM
            dd6iv = dd2iv * dd2iv * dd2iv                   ! #SLC2 #CMM
                                                            ! #SLC2 #CMM
            v1v = vdw6v * dd6iv                             ! #SLC2 #CMM
            v2v = vdw12v * dd6iv * dd6iv                    ! #SLC2 #CMM
            vv_sum = vv_sum - v1v + v2v                     ! #SLC2 #CMM
                                                            ! #SLC2 #CMM
            ccv = (-6.d0 * v1v + 12.d0 * v2v) * dd2iv       ! #SLC2 #CMM
            fxv = dxv*ccv ; fyv = dyv*ccv ; fzv = dzv*ccv   ! #SLC2 #CMM
                                                            ! #SLC2 #CMM
            gwxv = gwxv+fxv ; gwyv = gwyv+fyv               ! #SLC2 #CMM
            gwzv = gwzv+fzv                                 ! #SLC2 #CMM
            Tgrad(1,j_atoms,i) = Tgrad(1,j_atoms,i) + fxv   ! #SLCT #CMM #REDU #HIGH
            Tgrad(2,j_atoms,i) = Tgrad(2,j_atoms,i) + fyv   ! #SLCT #CMM #REDU #HIGH
            Tgrad(3,j_atoms,i) = Tgrad(3,j_atoms,i) + fzv   ! #SLCT #CMM #REDU #HIGH
          enddo                                             ! #SLC2 #CMM
                                                            ! #SLC2 #CMM
          ! vdW, remain                                     ! #SLC2 #CMM
          do j = itmp+1,ntbvdW(ii,i)                        ! #SLC2 #CMM
            kk = itbvdW(j,ii,i)                             ! #SLC2 #CMM
            dx=xi-Tcord(1,kk) ; dy=yi-Tcord(2,kk)           ! #SLC2 #CMM
            dz=zi-Tcord(3,kk)                               ! #SLC2 #CMM
            dd2i = 1.d0 / (dx*dx + dy*dy + dz*dz)           ! #SLC2 #CMM
            dd6i = dd2i * dd2i * dd2i                       ! #SLC2 #CMM
                                                            ! #SLC2 #CMM
            it2 = Tixatyp(kk)                               ! #SLC2 #CMM
            v1 = fynbpp(1,it2,it1) * dd6i                   ! #SLC2 #CMM
            v2 = fynbpp(2,it2,it1) * dd6i * dd6i            ! #SLC2 #CMM
            v_sum = v_sum - v1 + v2                         ! #SLC2 #CMM
                                                            ! #SLC2 #CMM
            cc = (-6.d0*v1 + 12.d0*v2)*dd2i                 ! #SLC2 #CMM
            fx = dx*cc ; fy = dy*cc ; fz = dz*cc            ! #SLC2 #CMM
                                                            ! #SLC2 #CMM
            gwx = gwx + fx ; gwy = gwy + fy                 ! #SLC2 #CMM
            gwz = gwz + fz                                  ! #SLC2 #CMM
            Tgrad(1,kk,i) = Tgrad(1,kk,i) + fx              ! #SLCT #CMM #REDU #HIGH
            Tgrad(2,kk,i) = Tgrad(2,kk,i) + fy              ! #SLCT #CMM #REDU #HIGH
            Tgrad(3,kk,i) = Tgrad(3,kk,i) + fz              ! #SLCT #CMM #REDU #HIGH
          enddo                                             ! #SLC2 #CMM
          
          gwx=gwx+sum(gwxv) ; gwy=gwy+sum(gwyv) ; gwz=gwz+sum(gwzv)
          Tgrad(1,ii,i) = Tgrad(1,ii,i) - gwx
          Tgrad(2,ii,i) = Tgrad(2,ii,i) - gwy
          Tgrad(3,ii,i) = Tgrad(3,ii,i) - gwz
        enddo
        enddo  ! #SLC2 #HIGH
        !$OMP end do nowait

        ee_sum = ee_sum + sum(eev_sum)
        v_sum = v_sum + sum(vv_sum)
        eelecn(i) = eelecn(i) + ee_sum
        evancn(i) = evancn(i) + v_sum
      enddo          ! #SLC2 #HIGH
      !$OMP barrier  ! #SLC2 #HIGH
      enddo
      !$OMP end parallel

!************************************************************

      ! hydrogen bond energy
      if ( iyeflg(10) .eq. 1 ) then
        !$OMP parallel default(none)                                 & !
        !$OMP private(ii,it1,xi,yi,zi,i,itmp,j,kk,dx,dy,dz,dd2i,dd6i,& !
        !$OMP         ee,it2,v1,v2,cc,fx,fy,fz,ee_sum,v_sum,chgi,gwx,& !
        !$OMP         gwy,gwz,xv,yv,zv,chgv,gwxv,gwyv,gwzv,vdw6v,    & !
        !$OMP         vdw12v,dxv,dyv,dzv,dd2iv,dd6iv,eev,v1v,v2v,ccv,& !
        !$OMP         fxv,fyv,fzv,eev_sum,vv_sum,j_atoms             & !
        !$OMP         )                                & ! ! #SLCT #CMM #REDU #DOUB
        !$OMP         ,k,kj)                           & ! ! #SLCT #CMM #HIGH
        !$OMP         ,dd2,ddi,dd2v,ddiv,cov,co)       & ! ! #SLCT #ZD #REDU #DOUB #PB #CB
        !$OMP         ,k,kj,dd2,ddi,dd2v,ddiv,cov,co)  & ! ! #SLCT #ZD #HIGH #PB #CB
        !$OMP         ,cov,co)                         & ! ! #SLCT #ACT #REDU #DOUB #PB #CB
        !$OMP         ,k,kj,cov,co)                    & ! ! #SLCT #ACT #HIGH #PB #CB
        !$OMP shared(ixnatm,Tcord,Tixatyp,nfrg,ntbhyd,itbhyd,fynbpp, & !
        !$OMP        ntbhyd2,itbhyd2,i_vec                           & !
        !$OMP        ,irlim2                           & ! ! #SLC2 #ZD #ACT
        !$OMP        ,bcoeff,fcoeff,zcore              & ! ! #SLC2 #ZD
        !$OMP        ,fxcell,invcel                    & ! ! #SLC2 #PB
        !$OMP        )                                 & ! ! #SLC2 #REDU
        !$OMP        ,high_para,np_cell,para_se,Tgrad) & ! ! #SLC2 #HIGH
        !$OMP        ,Tgrad)                           & ! ! #SLC2 #DOUB
        !$OMP reduction(+ : eelecn,ehydcn,Tgrad)       ! #SLC2 #REDU
        !$OMP reduction(+ : eelecn,ehydcn)             ! #SLC2 #HIGH #DOUB
        do k = 1,high_para                             ! #SLC2 #HIGH
        do i = nfrg,1,-1
          ee_sum = 0.d0 ; v_sum = 0.d0 ; eev_sum = 0.d0 ; vv_sum = 0.d0
          !$OMP do schedule (static)                   ! #SLC2 #REDU #DOUB
          do ii = 1,ixnatm                             ! #SLC2 #REDU #DOUB
          !$OMP do schedule(dynamic)                   ! #SLC2 #HIGH
          do kj = 1,np_cell(k)                         ! #SLC2 #HIGH
          do ii = para_se(1,kj,k),para_se(2,kj,k)      ! #SLC2 #HIGH
            it1 = Tixatyp(ii)
            gwx = 0.d0 ; gwy = 0.d0 ; gwz = 0.d0
            gwxv = 0.d0 ; gwyv = 0.d0 ; gwzv = 0.d0
            xi = Tcord(1,ii) ; yi = Tcord(2,ii) ; zi = Tcord(3,ii)
            chgi = Tcord(4,ii)

            itmp = ntbhyd(ii,i) - mod(ntbhyd(ii,i),i_vec)

            ! hydrogen & electrostatic, hand-vectrized
            do j = 1,itmp,i_vec
              j_atoms = itbhyd(j:j+i_vec-1,ii,i)
              xv=Tcord(1,j_atoms) ; yv=Tcord(2,j_atoms)
              zv=Tcord(3,j_atoms) ; chgv = Tcord(4,j_atoms)
              dxv = xi - xv ; dyv = yi - yv ; dzv = zi - zv
              dxv = dxv - fxcell(1) * nint(dxv*invcel(1))      ! #SLC2 #PB
              dyv = dyv - fxcell(2) * nint(dyv*invcel(2))      ! #SLC2 #PB
              dzv = dzv - fxcell(3) * nint(dzv*invcel(3))      ! #SLC2 #PB
              vdw6v = fynbpp(3,Tixatyp(j_atoms),it1)
              vdw12v = fynbpp(4,Tixatyp(j_atoms),it1)
              dd2iv = 1.d0 / (dxv*dxv + dyv*dyv + dzv*dzv)     ! #SLC2 #CMM #ACT
              dd2v = dxv*dxv + dyv*dyv + dzv*dzv               ! #SLC2 #ZD
              dd2iv = 1.d0 / dd2v                              ! #SLC2 #ZD
              cov = sign(0.5d0,dd2iv-irlim2) + 0.5d0           ! #SLC2 #ZD #ACT
              ddiv = sqrt(dd2iv)                               ! #SLC2 #ZD
              dd6iv = dd2iv*dd2iv
              dd6iv = dd6iv*dd6iv*dd2iv                        ! #SLC2 #CMM
              dd6iv = dd6iv*dd6iv*dd2iv*cov                    ! #SLC2 #ZD #ACT

              eev = chgi * chgv * sqrt(dd2iv)                  ! #SLC2 #CMM
              eev = chgi * chgv * cov                          ! #SLC2 #ZD
              eev = chgi * chgv * sqrt(dd2iv) * cov            ! #SLC2 #ACT
              eev_sum = eev_sum + eev                          ! #SLC2 #CMM #ACT
              eev_sum = eev_sum + eev*(ddiv+bcoeff*dd2v-zcore) ! #SLC2 #ZD

              v1v = vdw6v * dd6iv
              v2v = vdw12v * dd6iv * dd2iv
              vv_sum = vv_sum - v1v + v2v

              ccv = (eev -10.d0*v1v + 12.d0*v2v) * dd2iv       ! #SLC2 #CMM #ACT
              ccv = eev*(dd2iv*ddiv-fcoeff) +                & ! #SLC2 #ZD
                    (-10.d0*v1v + 12.d0*v2v)*dd2iv             ! #SLC2 #ZD
              fxv = dxv*ccv ; fyv = dyv*ccv ; fzv = dzv*ccv

              gwxv = gwxv+fxv ; gwyv = gwyv+fyv ; gwzv = gwzv+fzv 
              Tgrad(1,j_atoms,i) = Tgrad(1,j_atoms,i) + fxv    ! #SLC2 #REDU #HIGH
              Tgrad(2,j_atoms,i) = Tgrad(2,j_atoms,i) + fyv    ! #SLC2 #REDU #HIGH
              Tgrad(3,j_atoms,i) = Tgrad(3,j_atoms,i) + fzv    ! #SLC2 #REDU #HIGH
            enddo

            ! hydrogen & electrostatic, remain
            do j = itmp+1,ntbhyd(ii,i)
              kk = itbhyd(j,ii,i) ; it2 = Tixatyp(kk)
              dx=xi-Tcord(1,kk) ; dy=yi-Tcord(2,kk) ; dz=zi-Tcord(3,kk)
              dx = dx - fxcell(1) * nint(dx*invcel(1))    ! #SLC2 #PB
              dy = dy - fxcell(2) * nint(dy*invcel(2))    ! #SLC2 #PB
              dz = dz - fxcell(3) * nint(dz*invcel(3))    ! #SLC2 #PB
              dd2i = 1.d0 / (dx*dx + dy*dy + dz*dz)       ! #SLC2 #CMM #ACT
              dd2 = dx*dx + dy*dy + dz*dz                 ! #SLC2 #ZD
              dd2i = 1.d0 / dd2                           ! #SLC2 #ZD
              co = sign(0.5d0,dd2i-irlim2) + 0.5d0        ! #SLC2 #ZD #ACT
              ddi = sqrt(dd2i)                            ! #SLC2 #ZD
              dd6i = dd2i * dd2i
              dd6i = dd6i * dd6i * dd2i                   ! #SLC2 #CMM
              dd6i = dd6i * dd6i * dd2i * co              ! #SLC2 #ZD #ACT

              ee = chgi * Tcord(4,kk) * sqrt(dd2i)        ! #SLC2 #CMM
              ee = chgi * Tcord(4,kk) *co                 ! #SLC2 #ZD
              ee = chgi * Tcord(4,kk) * sqrt(dd2i) *co    ! #SLC2 #ACT
              ee_sum = ee_sum + ee                        ! #SLC2 #CMM #ACT
              ee_sum = ee_sum + ee*(ddi+bcoeff*dd2-zcore) ! #SLC2 #ZD


              v1 = fynbpp(3,it2,it1) * dd6i
              v2 = fynbpp(4,it2,it1) * dd6i * dd2i
              v_sum = v_sum - v1 + v2 

              cc = (ee -10.d0*v1 + 12.d0*v2) * dd2i       ! #SLC2 #CMM #ACT
              cc = ee*(dd2i*ddi-fcoeff) +               & ! #SLC2 #ZD
                    (-10.d0*v1 + 12.d0*v2)*dd2i           ! #SLC2 #ZD
              fx = dx*cc ; fy = dy*cc ; fz = dz*cc

              gwx = gwx + fx ; gwy = gwy + fy ; gwz = gwz + fz
              Tgrad(1,kk,i) = Tgrad(1,kk,i) + fx          ! #SLC2 #REDU #HIGH
              Tgrad(2,kk,i) = Tgrad(2,kk,i) + fy          ! #SLC2 #REDU #HIGH
              Tgrad(3,kk,i) = Tgrad(3,kk,i) + fz          ! #SLC2 #REDU #HIGH
            enddo

            !  The added hydrogen interactions for pairs for that  ! #SLC2 #CMM
            !  the pair-distances are shorter than a limit.        ! #SLC2 #CMM
            itmp = ntbhyd2(ii,i) - mod(ntbhyd2(ii,i),i_vec)        ! #SLC2 #CMM
            ! hydrogen, hand-vectrized                             ! #SLC2 #CMM
            do j = 1,itmp,i_vec                                    ! #SLC2 #CMM
              j_atoms = itbhyd2(j:j+i_vec-1,ii,i)                  ! #SLC2 #CMM
              xv=Tcord(1,j_atoms) ; yv=Tcord(2,j_atoms)            ! #SLC2 #CMM
              zv=Tcord(3,j_atoms)                                  ! #SLC2 #CMM
              vdw6v = fynbpp(3,Tixatyp(j_atoms),it1)               ! #SLC2 #CMM
              vdw12v = fynbpp(4,Tixatyp(j_atoms),it1)              ! #SLC2 #CMM
              dxv = xi - xv ; dyv = yi - yv ; dzv = zi - zv        ! #SLC2 #CMM
              dd2iv = 1.d0 / (dxv*dxv + dyv*dyv + dzv*dzv)         ! #SLC2 #CMM
              dd6iv = dd2iv*dd2iv                                  ! #SLC2 #CMM
              dd6iv = dd6iv*dd6iv*dd2iv                            ! #SLC2 #CMM
                                                                   ! #SLC2 #CMM
              v1v = vdw6v * dd6iv                                  ! #SLC2 #CMM
              v2v = vdw12v * dd6iv * dd6iv                         ! #SLC2 #CMM
              vv_sum = vv_sum - v1v + v2v                          ! #SLC2 #CMM
                                                                   ! #SLC2 #CMM
              ccv = (-10.d0 * v1v + 12.d0 * v2v) * dd2iv           ! #SLC2 #CMM
              fxv = dxv*ccv ; fyv = dyv*ccv ; fzv = dzv*ccv        ! #SLC2 #CMM
                                                                   ! #SLC2 #CMM
              gwxv = gwxv+fxv ; gwyv = gwyv+fyv ; gwzv = gwzv+fzv  ! #SLC2 #CMM
              Tgrad(1,j_atoms,i) = Tgrad(1,j_atoms,i) + fxv        ! #SLCT #CMM #REDU #HIGH
              Tgrad(2,j_atoms,i) = Tgrad(2,j_atoms,i) + fyv        ! #SLCT #CMM #REDU #HIGH
              Tgrad(3,j_atoms,i) = Tgrad(3,j_atoms,i) + fzv        ! #SLCT #CMM #REDU #HIGH
            enddo                                                  ! #SLC2 #CMM
                                                                   ! #SLC2 #CMM
            ! hydrogen, remain                                     ! #SLC2 #CMM
            do j = itmp+1,ntbhyd2(ii,i)                            ! #SLC2 #CMM
              kk = itbhyd2(j,ii,i)                                 ! #SLC2 #CMM
              dx = xi-Tcord(1,kk) ; dy = yi-Tcord(2,kk)            ! #SLC2 #CMM
              dz = zi-Tcord(3,kk)                                  ! #SLC2 #CMM
              dd2i = 1.d0 / (dx*dx + dy*dy + dz*dz)                ! #SLC2 #CMM
              dd6i = dd2i*dd2i                                     ! #SLC2 #CMM
              dd6i = dd6i*dd6i*dd2i                                ! #SLC2 #CMM
                                                                   ! #SLC2 #CMM
              it2 = Tixatyp(kk)                                    ! #SLC2 #CMM
              v1 = fynbpp(3,it2,it1) * dd6i                        ! #SLC2 #CMM
              v2 = fynbpp(4,it2,it1) * dd6i * dd6i                 ! #SLC2 #CMM
              v_sum = v_sum - v1 + v2                              ! #SLC2 #CMM
                                                                   ! #SLC2 #CMM
              cc = (-10.d0*v1 + 12.d0*v2)*dd2i                     ! #SLC2 #CMM
              fx = dx*cc ; fy = dy*cc ; fz = dz*cc                 ! #SLC2 #CMM
                                                                   ! #SLC2 #CMM
              gwx = gwx + fx ; gwy = gwy + fy ; gwz = gwz + fz     ! #SLC2 #CMM
              Tgrad(1,kk,i) = Tgrad(1,kk,i) + fx                   ! #SLCT #CMM #REDU #HIGH
              Tgrad(2,kk,i) = Tgrad(2,kk,i) + fy                   ! #SLCT #CMM #REDU #HIGH
              Tgrad(3,kk,i) = Tgrad(3,kk,i) + fz                   ! #SLCT #CMM #REDU #HIGH
            enddo                                                  ! #SLC2 #CMM

            gwx=gwx+sum(gwxv) ; gwy=gwy+sum(gwyv) ; gwz=gwz+sum(gwzv)
            Tgrad(1,ii,i) = Tgrad(1,ii,i) - gwx
            Tgrad(2,ii,i) = Tgrad(2,ii,i) - gwy
            Tgrad(3,ii,i) = Tgrad(3,ii,i) - gwz
          enddo
          enddo     ! #SLC2 #HIGH
          !$OMP end do nowait

          ee_sum = ee_sum + sum(eev_sum)
          v_sum = v_sum + sum(vv_sum)
          eelecn(i) = eelecn(i) + ee_sum
          ehydcn(i) = ehydcn(i) + v_sum
        enddo         ! #SLC2 #HIGH
        !$OMP barrier ! #SLC2 #HIGH
        enddo
        !$OMP end parallel
      endif

!*********************************************************** ! #SLC2 #ZD
                                                             ! #SLC2 #ZD
      ! excess energy calculation                            ! #SLC2 #ZD
      Eexc(:) = 0.d0                                         ! #SLC2 #ZD
      !$OMP parallel default(none)                       & ! ! #SLC2 #ZD
      !$OMP private(i,j,ii,it1,it2,dx,dy,dz,dd2,cc)      & ! ! #SLC2 #ZD
      !$OMP shared(np_excess,para_excess,excls,iexcess,  & ! ! #SLC2 #ZD
      !$OMP        fxcell,invcel,                        & ! ! #SLCT #ZD #REDU #HIGH #DOUB #PB
      !$OMP        cord,chgmod,bcoeff,fcoeff,grad)       & ! ! #SLC2 #ZD
      !$OMP reduction (+ : Eexc)                             ! #SLC2 #ZD
      do i = 1,np_excess                                     ! #SLC2 #ZD
        !$OMP do                                             ! #SLC2 #ZD
        do j = para_excess(1,i),para_excess(2,i)             ! #SLC2 #ZD
          ii = excls(j)                                      ! #SLC2 #ZD
          it1 = iexcess(1,j) ; it2 = iexcess(2,j)            ! #SLC2 #ZD
          dx = cord(1,it1) - cord(1,it2)                     ! #SLC2 #ZD
          dy = cord(2,it1) - cord(2,it2)                     ! #SLC2 #ZD
          dz = cord(3,it1) - cord(3,it2)                     ! #SLC2 #ZD
          dx = dx - fxcell(1) * nint(dx*invcel(1))           ! #SLCT #ZD #REDU #HIGH #DOUB #PB
          dy = dy - fxcell(2) * nint(dy*invcel(2))           ! #SLCT #ZD #REDU #HIGH #DOUB #PB
          dz = dz - fxcell(3) * nint(dz*invcel(3))           ! #SLCT #ZD #REDU #HIGH #DOUB #PB
          dd2 = dx*dx + dy*dy + dz*dz                        ! #SLC2 #ZD
          cc = chgmod(it1) * chgmod(it2)                     ! #SLC2 #ZD
          Eexc(ii) = Eexc(ii) + cc * bcoeff * dd2            ! #SLC2 #ZD
          cc = cc * fcoeff                                   ! #SLC2 #ZD
          dx = dx * cc ; dy = dy * cc ; dz = dz * cc         ! #SLC2 #ZD
          grad(1,it1,ii) = grad(1,it1,ii) + dx               ! #SLC2 #ZD
          grad(2,it1,ii) = grad(2,it1,ii) + dy               ! #SLC2 #ZD
          grad(3,it1,ii) = grad(3,it1,ii) + dz               ! #SLC2 #ZD
          grad(1,it2,ii) = grad(1,it2,ii) - dx               ! #SLC2 #ZD
          grad(2,it2,ii) = grad(2,it2,ii) - dy               ! #SLC2 #ZD
          grad(3,it2,ii) = grad(3,it2,ii) - dz               ! #SLC2 #ZD
        enddo                                                ! #SLC2 #ZD
        !$OMP end do                                         ! #SLC2 #ZD
      enddo                                                  ! #SLC2 #ZD
      !$OMP end parallel                                     ! #SLC2 #ZD
                                                             ! #SLC2 #ZD
!*********************************************************** ! #SLC2 #ZD
                                                             ! #SLC2 #ZD
      !$OMP parallel default(none)                                   & !
      !$OMP private(i,j)                                             & !
      !$OMP shared(grad,ixnatm,nfrg,Tgrad,inum)
      !$OMP do schedule (static) collapse(2)
      do j = 1,nfrg
      do i = 1,ixnatm
        grad(:,i,j) = grad(:,i,j) + Tgrad(:,inum(i),j)
      enddo
      enddo
      !$OMP end do
      !$OMP end parallel
      eelecn(:) = eelecn(:) * 0.5d0              ! #SLC2 #DOUB
      evancn(:) = evancn(:) * 0.5d0              ! #SLC2 #DOUB
      ehydcn(:) = ehydcn(:) * 0.5d0              ! #SLC2 #DOUB
      eelecn(:) = eelecn(:) + dself(:) + Eexc(:) ! #SLC2 #ZD

!      call system_clock(tim2)
!      tmptime = tmptime + tim2 - tim1

!************************************************************

      return
      end subroutine enear_CMM_redu ! #SLCT #CMM #REDU
      end subroutine enear_CMM_high ! #SLCT #CMM #HIGH
      end subroutine enear_CMM_doub ! #SLCT #CMM #DOUB
      end subroutine enear_ACT_redu_PB ! #SLCT #ACT #REDU #PB
      end subroutine enear_ACT_high_PB ! #SLCT #ACT #HIGH #PB
      end subroutine enear_ACT_doub_PB ! #SLCT #ACT #DOUB #PB
      end subroutine enear_ACT_redu_CB ! #SLCT #ACT #REDU #CB
      end subroutine enear_ACT_high_CB ! #SLCT #ACT #HIGH #CB
      end subroutine enear_ACT_doub_CB ! #SLCT #ACT #DOUB #CB
      end subroutine enear_ZD_redu_PB  ! #SLCT #ZD  #REDU #PB
      end subroutine enear_ZD_high_PB  ! #SLCT #ZD  #HIGH #PB
      end subroutine enear_ZD_doub_PB  ! #SLCT #ZD  #DOUB #PB
      end subroutine enear_ZD_redu_CB  ! #SLCT #ZD  #REDU #CB
      end subroutine enear_ZD_high_CB  ! #SLCT #ZD  #HIGH #CB
      end subroutine enear_ZD_doub_CB  ! #SLCT #ZD  #DOUB #CB
