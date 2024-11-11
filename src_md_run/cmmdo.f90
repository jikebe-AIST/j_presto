
      subroutine cmmdo(eelecn,evancn,ehydcn)

!**************************************************************************
!
!     calc. CMM energy & forces
!
!*************************************************************************

      use COMBAS ; use COMERG ; use COMCMM ; use COMCMMC
!      use CALC_TIME
      !$ use omp_lib

      implicit none

      real(8),intent(out):: eelecn(nfrg),evancn(nfrg),ehydcn(nfrg)

      integer(4):: i,k,ii,ix,iy,iz,kk,n,itmp
      integer(4):: num,Tipcl(ipmax)

!**************************************

      ! SIMD check
      if ( SIMD_chk ) call SIMD_check(1)

      ! Exact vdW & ele. Energy calculation intra-near-cells
      select case (iy15m)
      case (1)
        if ( high_para .eq. 0 ) then
          call enear_CMM_redu(eelecn,evancn,ehydcn)
        elseif ( high_para .gt. 0 ) then
          call enear_CMM_high(eelecn,evancn,ehydcn)
        else
          call enear_CMM_doub(eelecn,evancn,ehydcn)
        endif
      case (2)
        select case (ixfbou)
        case (1)
          if ( high_para .eq. 0 ) then
            call enear_ACT_redu_PB(eelecn,evancn,ehydcn)
          elseif ( high_para .gt. 0 ) then
            call enear_ACT_high_PB(eelecn,evancn,ehydcn)
          else
            call enear_ACT_doub_PB(eelecn,evancn,ehydcn)
          endif
        case default
          if ( high_para .eq. 0 ) then
            call enear_ACT_redu_CB(eelecn,evancn,ehydcn)
          elseif ( high_para .gt. 0 ) then
            call enear_ACT_high_CB(eelecn,evancn,ehydcn)
          else
            call enear_ACT_doub_CB(eelecn,evancn,ehydcn)
          endif
        end select
      case (4)
        select case (ixfbou)
        case (1)
          if ( high_para .eq. 0 ) then
            call enear_ZD_redu_PB(eelecn,evancn,ehydcn)
          elseif ( high_para .gt. 0 ) then
            call enear_ZD_high_PB(eelecn,evancn,ehydcn)
          else
            call enear_ZD_doub_PB(eelecn,evancn,ehydcn)
          endif
        case default
          if ( high_para .eq. 0 ) then
            call enear_ZD_redu_CB(eelecn,evancn,ehydcn)
          elseif ( high_para .gt. 0 ) then
            call enear_ZD_high_CB(eelecn,evancn,ehydcn)
          else
            call enear_ZD_doub_CB(eelecn,evancn,ehydcn)
          endif
        end select
      end select

      ! SIMD check
      if ( SIMD_chk ) call SIMD_check(2)

      ! For CMM
      select case (iy15m)
      case (1)
        if ( extCMM_flag ) call enear_extCMM(eelecn)

        ! CMM calculation for inter-cells ele. energy & forces
        if ( itwin .eq. 1 ) then
          !! Calc. among smallest cells
          !! select a target cell (ix,iy,iz)
          !! Energy terms for CMM. and force from CMM cells (gwk)
!          call system_clock(tim1)
          Ecmm(1:nlev,1:nfrg) = 0.d0
          !$OMP parallel default (none)                              & !
          !$OMP private(k,kk)                                        & !
          !$OMP shared(nfrg,ixnatm,gwk)
          !$OMP do schedule (static) collapse(2)
          do k = 1,nfrg
          do kk = 1,ixnatm
            gwk(:,kk,k) = 0.d0
          enddo
          enddo
          !$OMP end do
          !$OMP end parallel

          !$OMP parallel default(none)                               & !
          !$OMP private(k,ii,iz,iy,ix,num,Tipcl,kk,n,itmp)           & !
          !$OMP shared(Nmycel,mycel,nfrag,n_matrix,npcl1,CMMpot1,c1, & !
          !$OMP        ipcl,cord,gwk,chgmod2,nlev,CMMpot2,c2,CMMpot3,& !
          !$OMP        c3,CMMpot4,c4,CMMpot5,c5,i_vec)               & !
          !$OMP reduction (+ : Ecmm)
          do k = 1,nfrag
            !$OMP do schedule (guided)
            do ii = 1,Nmycel(k)
              ix = mycel(1,ii,k) ;iy = mycel(2,ii,k) ;iz = mycel(3,ii,k)
              num = npcl1(k,ix,iy,iz) ; itmp = num - mod(num,i_vec)
              Tipcl(1:num) = ipcl(1:num,k,ix,iy,iz)

              do kk = 1,nfrag
                n = n_matrix(kk,k)
                !!! Loop over the level 1 cells around the target cell
                call compute_particle_multipole(num,itmp,Tipcl,1,      &
                     CMMpot1,c1,ii,kk,k,n,Ecmm)

                !!! Loop over the level 2 cells around the target cell
                if ( nlev .lt. 2 ) cycle
                call compute_particle_multipole(num,itmp,Tipcl,2,      &
                     CMMpot2,c2,ii,kk,k,n,Ecmm)

                !!! Loop over the level 3 cells around the target cell
                if ( nlev .lt. 3 ) cycle
                call compute_particle_multipole(num,itmp,Tipcl,3,      &
                     CMMpot3,c3,ii,kk,k,n,Ecmm)

                !!! Loop over the level 4 cells around the target cell
                if ( nlev .lt. 4 ) cycle
                call compute_particle_multipole(num,itmp,Tipcl,4,      &
                     CMMpot4,c4,ii,kk,k,n,Ecmm)

                !!! Loop over the level 5 cells around the target cell
                if ( nlev .lt. 5 ) cycle
                call compute_particle_multipole(num,itmp,Tipcl,5,      &
                     CMMpot5,c5,ii,kk,k,n,Ecmm)
              enddo

            enddo
            !$OMP end do nowait
          enddo
          !$OMP end parallel

          Ecmm(1:nlev,1:nfrg) = Ecmm(1:nlev,1:nfrg) * 0.5d0
          !$OMP workshare
          forall ( i=1:ixnatm )                                        &
            gwk(1:3,i,1:nfrg) = gwk(1:3,i,1:nfrg) * chgmod2(i)
          !$OMP end workshare

!          call system_clock(tim2)
!          tmptime = tmptime + tim2 - tim1
        endif

        forall ( i=1:nfrg ) eelecn(i) =eelecn(i)+sum(Ecmm(1:nlev,i))
        !$OMP parallel default (none)                                & !
        !$OMP private(i,k)                                           & !
        !$OMP shared(nfrg,ixnatm,grad,gwk)
        !$OMP do schedule (static) collapse(2)
        do k = 1,nfrg
        do i = 1,ixnatm
          grad(:,i,k) = grad(:,i,k) + gwk(:,i,k)
        enddo
        enddo
        !$OMP end do
        !$OMP end parallel
      end select

!************************************************************

      return

!************************************************************


      contains

      subroutine compute_particle_multipole(num,itmp,Tipcl,cell_level, &
                 CMMpot,cell_centers,ii,kk,k,n,tEcmm)

      implicit none

      integer(4),intent(in):: num,itmp,Tipcl(:),cell_level,ii,kk,k,n
      real(8),intent(in):: CMMpot(:,:,:,:,:),cell_centers(:,:)
      real(8),intent(inout):: tEcmm(:,:)

      real(8):: dx,dy,dz,dd2i,ddi,dd3i,dd5i,dxy,dyz,dzx,aa1,aa2,dx2,   &
                dy2,dz2,ee,bbx,bby,bbz,cc2x,cc2y,cc2z,cc4,q(13)
      integer(4):: i,j,iat,jx,jy,jz

      ! hand-vectorization
      real(8):: xv(i_vec),yv(i_vec),zv(i_vec),dxv(i_vec),dyv(i_vec),   &
                dzv(i_vec),dd2iv(i_vec),ddiv(i_vec),dd3iv(i_vec),      &
                dd5iv(i_vec),dxyv(i_vec),dyzv(i_vec),dzxv(i_vec),      &
                aa1v(i_vec),aa2v(i_vec),dx2v(i_vec),dy2v(i_vec),       &
                dz2v(i_vec),eev(i_vec),bbxv(i_vec),bbyv(i_vec),        &
                bbzv(i_vec),cc2xv(i_vec),cc2yv(i_vec),cc2zv(i_vec),    &
                cc4v(i_vec),fx(i_vec),fy(i_vec),fz(i_vec)
      integer(4):: j_atoms(i_vec)

!*************************************

      do i = 1,Nyrcel(ii,cell_level,kk,k)
        jx = yrcel(1,i,ii,cell_level,kk,k)
        jy = yrcel(2,i,ii,cell_level,kk,k)
        jz = yrcel(3,i,ii,cell_level,kk,k)
        !* Get parameters for CMM
        q(1:13) = CMMpot(1:13,kk,jx,jy,jz)

        !!* Calc. CMM between an atom in the target cell & 
        !!*       the surounding cells
        !!* select an atom in the target cell

        !!! hand-vectorization
        do j = 1,itmp,i_vec
          j_atoms(1:i_vec) = Tipcl(j:j+i_vec-1)
          xv = cord(1,j_atoms) ; dxv = xv - cell_centers(1,jx)
          yv = cord(2,j_atoms) ; dyv = yv - cell_centers(2,jy)
          zv = cord(3,j_atoms) ; dzv = zv - cell_centers(3,jz)
          dx2v = dxv*dxv ; dy2v = dyv*dyv ; dz2v = dzv*dzv
          dd2iv = 1.d0 / (dx2v + dy2v + dz2v)
          ddiv = sqrt(dd2iv)
          dd3iv = ddiv*dd2iv
          dd5iv = dd3iv*dd2iv
          dxyv = dxv*dyv ; dyzv = dyv*dzv ; dzxv = dzv*dxv

          !** Energy
          aa1v = q(2)*dxv + q(3)*dyv + q(4)*dzv
          aa2v = q(5)*dx2v + q(6)*dy2v + q(7)*dz2v +                   &
                 q(11)*dxyv + q(12)*dyzv + q(13)*dzxv
          eev = q(1)*ddiv + aa1v*dd3iv + aa2v*dd5iv
          tEcmm(cell_level,n) = tEcmm(cell_level,n) +                  &
                                sum(chgmod2(j_atoms)*eev)

          !** Force
          bbxv = (q(1)*dxv - q(2)) * dd3iv
          bbyv = (q(1)*dyv - q(3)) * dd3iv
          bbzv = (q(1)*dzv - q(4)) * dd3iv
          cc2xv = q(8)*dxv + q(11)*dyv + q(13)*dzv
          cc2yv = q(11)*dxv + q(9)*dyv + q(12)*dzv
          cc2zv = q(13)*dxv + q(12)*dyv + q(10)*dzv
          cc4v = aa1v*3.d0 + aa2v*dd2iv*5.d0
          fx = -bbxv + (cc2xv-cc4v*dxv)*dd5iv
          fy = -bbyv + (cc2yv-cc4v*dyv)*dd5iv
          fz = -bbzv + (cc2zv-cc4v*dzv)*dd5iv
          gwk(1,j_atoms,n) = gwk(1,j_atoms,n) + fx
          gwk(2,j_atoms,n) = gwk(2,j_atoms,n) + fy
          gwk(3,j_atoms,n) = gwk(3,j_atoms,n) + fz
        enddo
        !!! remain of the vectrization
        do j = itmp+1,num
          iat = Tipcl(j)
          dx = cord(1,iat) - cell_centers(1,jx)
          dy = cord(2,iat) - cell_centers(2,jy)
          dz = cord(3,iat) - cell_centers(3,jz)
          dx2 = dx*dx ; dy2 = dy*dy ; dz2 = dz*dz
          dd2i = 1.d0 / (dx2 + dy2 + dz2)
          ddi = sqrt(dd2i)
          dd3i = ddi*dd2i
          dd5i = dd3i*dd2i
          dxy = dx*dy ; dyz = dy*dz ; dzx = dz*dx

          !** Energy
          aa1 = q(2)*dx + q(3)*dy + q(4)*dz
          aa2 = q(5)*dx2 + q(6)*dy2 + q(7)*dz2 +                       &
                q(11)*dxy + q(12)*dyz + q(13)*dzx
          ee = q(1)*ddi + aa1*dd3i + aa2*dd5i
          tEcmm(cell_level,n) = tEcmm(cell_level,n) + chgmod2(iat)*ee

          !** Force
          bbx = (q(1)*dx - q(2)) * dd3i
          bby = (q(1)*dy - q(3)) * dd3i
          bbz = (q(1)*dz - q(4)) * dd3i
          cc2x = q(8)*dx + q(11)*dy + q(13)*dz
          cc2y = q(11)*dx + q(9)*dy + q(12)*dz
          cc2z = q(13)*dx + q(12)*dy + q(10)*dz
          cc4 = aa1*3.d0 + aa2*dd2i*5.d0
          gwk(1,iat,n) = gwk(1,iat,n) - bbx + (cc2x-cc4*dx)*dd5i
          gwk(2,iat,n) = gwk(2,iat,n) - bby + (cc2y-cc4*dy)*dd5i
          gwk(3,iat,n) = gwk(3,iat,n) - bbz + (cc2z-cc4*dz)*dd5i
        enddo
      enddo

!***************************

      return
      end subroutine compute_particle_multipole

      end subroutine cmmdo
