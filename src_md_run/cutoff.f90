
      subroutine cutoff_PB(iprint,ilflag,ier) ! #SLC2 #PB
      subroutine cutoff_CB(iprint,ilflag,ier) ! #SLC2 #CB

!*******************************************************************
!
!     MAKING NEAR TABLE LIST FOR CUT-OFF CALC.
!
!*******************************************************************

      use COMBAS ; use COMERG ; use COMCMM ; use COMCMMC
      !$ use omp_lib

      implicit none

      integer(4),intent(in):: iprint,ilflag
      integer(4),intent(inout):: ier
      integer(4):: ii,iat1,i,j,n,kk,iat2,nn,nchk,itype,jat1,jat2
      real(8):: rtmp

!************************************************

      ier = 0
      call cutatm_PB(iprint,ilflag,ier) ! #SLC2 #PB
      call cutatm_CB(iprint,ilflag,ier) ! #SLC2 #CB
          
!*************************************
!     Check the tables.

      n15mxEL_max = max(n15mxEL_max,maxval(ntbEL(:,:)))
      if ( n15mxEL_max .gt. n15mxEL ) then
        print*,' ??????????????????????????????????????'
        print*,' ? n15mxEL_max is too small. (cutoff) ?'
        print*,' ??????????????????????????????????????'
        print*,' n15mxEL & the max.# detected  =',n15mxEL,' ',         &
               n15mxEL_max
        ier = 51
      endif
      n15mx_max = max(n15mx_max,maxval(ntb(:,:)))
      if ( n15mx_max .gt. n15mx ) then
        print*,' ??????????????????????????????????'
        print*,' ? n15mxEL is too small. (cutoff) ?'
        print*,' ??????????????????????????????????'
        print*,' n15mx & the max.# detected  =',n15mx,' ',n15mx_max
        ier = 51
      endif

!*************************************

      nchk = 0
      ! Remove 1-4 interaction & Rearrangement array ntb & itb
      !$OMP parallel default(none)                                   & !
      !$OMP private(iat1,n,iat2,kk,j,nn,jat1,jat2)                   & !
      !$OMP shared(nchk,ntb,itb,ntab14,itab14,ixnatm,n_matrix,icls,  & !
      !$OMP        itbEL,ntbEL,zero_vdW,inum)
      !$OMP do reduction(+ : nchk)                                   & !
      !$OMP    schedule (auto)
      do iat1 = 1,ixnatm
        jat1 = inum(iat1)
        OUTER : do n = 1,ntab14(iat1)
          iat2 = itab14(n,iat1)
          jat2 = inum(iat2)
          kk = n_matrix(icls(iat1),icls(iat2))
          if ( zero_vdW(iat1) .or. zero_vdW(iat2) ) then
            do j = 1,ntbEL(jat1,kk)
              if ( jat2 .eq. itbEL(j,jat1,kk) ) then
                nn = ntbEL(jat1,kk)
                itbEL(j,jat1,kk) = itbEL(nn,jat1,kk)
                ntbEL(jat1,kk) = nn - 1
                nchk = nchk + 1
                cycle OUTER
              endif
            enddo
          else
            do j = 1,ntb(jat1,kk)
              if ( jat2 .eq. itb(j,jat1,kk) ) then
                nn = ntb(jat1,kk)
                itb(j,jat1,kk) = itb(nn,jat1,kk)
                ntb(jat1,kk) = nn - 1
                nchk = nchk + 1
                cycle OUTER
              endif
            enddo
          endif
        enddo OUTER
      enddo
      !$OMP end do
      !$OMP end parallel

      j = nchk ; if ( high_para .lt. 0 ) j = nchk / 2
      if ( j .ne. nntab14 ) then
        print*,j,nntab14
        print*,"ERROR> CUTOFF"
        print*,"cell size is too small"
        ier = 99 ; return
      endif

      if ( iyeflg(10) .eq. 1 ) then
        ! Separate hydrogen bond energy
        ntbhyd = 0
        !$OMP parallel default(none)                                 & !
        !$OMP private(j,iat1,itype,ii,nn,i,iat2,n)                   & !
        !$OMP shared(nfrg,ixnatm,ntbEL,Tixatyp,itbEL,iyppid,itbhyd,  & !
        !$OMP        ntbhyd,ntb,itb)
        !$OMP do schedule (auto) collapse(2)
        do j = 1,nfrg
        do iat1 = 1,ixnatm
          if ( ntbEL(iat1,j) .ne. 0 ) then
            itype = Tixatyp(iat1)
            ii = 0 ; nn = 0
            do i = 1,ntbEL(iat1,j)
              iat2 = itbEL(i,iat1,j)
              if ( iyppid(itype,Tixatyp(iat2)) .eq. 2 ) then
                nn = nn + 1
                itbhyd(nn,iat1,j) = itbEL(i,iat1,j)
              else
                ii = ii + 1
                itbEL(ii,iat1,j) = itbEL(i,iat1,j)
              endif
            enddo
            ntbhyd(iat1,j) = nn
            ntbEL(iat1,j) = ntbEL(iat1,j) - nn
          endif
          if ( ntb(iat1,j) .ne. 0 ) then
            itype = Tixatyp(iat1)
            ii = 0 ; nn = 0 ; n = ntbhyd(iat1,j)
            do i = 1,ntb(iat1,j)
              iat2 = itb(i,iat1,j)
              if ( iyppid(itype,Tixatyp(iat2)) .eq. 2 ) then
                nn = nn + 1
                itbhyd(n+nn,iat1,j) = itb(i,iat1,j)
              else
                ii = ii + 1
                itb(ii,iat1,j) = itb(i,iat1,j)
              endif
            enddo
            ntbhyd(iat1,j) = nn + n
            ntb(iat1,j) = ntb(iat1,j) - nn
          endif
        enddo
        enddo
        !$OMP end do
        !$OMP end parallel

        ! Check the tables.
        n15mx_max = max(n15mx_max,maxval(ntbhyd(:,:)))
        if ( n15mx_max .gt. n15mx ) then
          print*,' ???????????????????????????????????????????'
          print*,' ? n15mx for ntbhyd is too small. (cutoff) ?'
          print*,' ???????????????????????????????????????????'
          print*,' n15mx & the max.# detected  =',n15mx,' ',n15mx_max
          ier = 51 ; return
        endif
      endif

!****************************************

      return
      end subroutine cutoff_PB ! #SLC2 #PB
      end subroutine cutoff_CB ! #SLC2 #CB
 

!=========================================================================


      subroutine cutatm_PB(iprint,ilflag,ier) ! #SLC2 #PB
      subroutine cutatm_CB(iprint,ilflag,ier) ! #SLC2 #CB

!*******************************************************************
!
!     * ATOM BASE CUT-OFF *
!       IF DISTANCE BETWEEN IATM AND JATM IS LESS EQUAL
!       CUT-OFF DISTANCE  , THEN THIS ATOM PAIR IS INTERACTING
!       ATOM PAIR.
!
!*******************************************************************

      use COMBAS ; use COMERG ; use COMCMM ; use COMCMMC
      !$ use omp_lib

      implicit none

      integer(4),intent(in):: iprint,ilflag
      integer(4),intent(inout):: ier
      integer(4):: ix,iy,iz,nel,ii,iat1,nce,ice,jx,jy,jz,i,j,n,border, &
                   kk,iat2,nn,jj
      integer(4):: num(nfrag),num2(nfrag)
      real(8):: dcod(3),dd,TTcod(3)

!************************************************

      if ( ilflag .eq. 2 ) then
        write(iprint,*)" "
        write(iprint,*)"INFORMATION> CUTOFF (V5.0) "
        write(iprint,*)"             ATOM BASE CUT-OFF IS DONE "
        write(iprint,*)" "
      endif

      !$OMP parallel default (none)                                  & !
      !$OMP private(i,j)                                             & !
      !$OMP shared(nfrg,ixnatm,ntb,ntbEL)
      !$OMP do schedule (static) collapse(2)
      do j = 1,nfrg
      do i = 1,ixnatm
        ntb(i,j) = 0 ; ntbEL(i,j) = 0
      enddo
      enddo
      !$OMP end do
      !$OMP end parallel
      if ( high_para .ge. 0 ) then
        !$OMP parallel default(none)                                   & !
        !$OMP private(iz,iy,ix,num,i,nel,ii,iat1,nce,ice,jx,jy,jz,num2,& !
        !$OMP         j,n,border,kk,iat2,dcod,dd,nn,TTcod)             & !
        !$OMP shared(nv,npcl1,nfrag,cell_nEL,ipcl,nnr,inr,n_matrix,    & !
        !$OMP        fxcell,invcel,                                    & ! ! #SLC2 #PB
        !$OMP        Tcord,rlim2,ntbEL,itbEL,ntb,itb,inum)
        !$OMP do schedule (dynamic) collapse(3)
        do iz = 1,nv(1)
        do iy = 1,nv(1)
        do ix = 1,nv(1)

          num(1:nfrag) = npcl1(1:nfrag,ix,iy,iz)
          if ( all(num .eq. 0 ) ) cycle

          do i = 1,nfrag
            if ( num(i) .eq. 0 ) cycle
            nel = cell_nEL(i,ix,iy,iz)

            ! (1a) Within a cell
            do ii = 1,num(i)
              iat1 = inum(ipcl(ii,i,ix,iy,iz))

              j = i ; n = n_matrix(j,i)
              if ( ii .le. nel ) then
                nn = ntbEL(iat1,n)
                do jj = ii+1,num(i)
                  iat2 = inum(ipcl(jj,j,ix,iy,iz))
                  nn = nn + 1
                  itbEL(nn,iat1,n) = iat2
                enddo
                ntbEL(iat1,n) = nn
              else
                nn = ntb(iat1,n)
                ! i-particle is EL+vdW, and jj>ii+1 is always EL+vdW
                do jj = ii+1,num(i)
                  iat2 = inum(ipcl(jj,j,ix,iy,iz))
                  nn = nn + 1
                  itb(nn,iat1,n) = iat2
                enddo
                ntb(iat1,n) = nn
              endif

              do j = i+1,nfrag
                n = n_matrix(j,i)
                if ( ii .le. nel ) then
                  border = num(j) ! all atoms are EL
                else
                  border = cell_nEL(j,ix,iy,iz)
                endif

                nn = ntbEL(iat1,n)
                do jj = 1,border
                  iat2 = inum(ipcl(jj,j,ix,iy,iz))
                  nn = nn + 1
                  itbEL(nn,iat1,n) = iat2
                enddo
                ntbEL(iat1,n) = nn

                nn = ntb(iat1,n)
                do jj = border+1,num(j)
                  iat2 = inum(ipcl(jj,j,ix,iy,iz))
                  nn = nn + 1
                  itb(nn,iat1,n) = iat2
                enddo
                ntb(iat1,n) = nn
              enddo

              TTcod(1:3) = Tcord(1:3,iat1)
              nce = nnr(ix,iy,iz)
              do ice = 1,nce
                jx = inr(1,ice,ix,iy,iz)
                jy = inr(2,ice,ix,iy,iz)
                jz = inr(3,ice,ix,iy,iz)

                num2(1:nfrag) = npcl1(1:nfrag,jx,jy,jz)
                if ( all(num2 .eq. 0) ) cycle

                do j = 1,nfrag
                  n = n_matrix(j,i)
                  if ( ii .le. nel ) then
                    border = num2(j)
                  else
                    border = cell_nEL(j,jx,jy,jz)
                  endif

                  nn = ntbEL(iat1,n)
                  do kk = 1,border
                    iat2 = inum(ipcl(kk,j,jx,jy,jz))
                    dcod(1:3) = TTcod(1:3) - Tcord(1:3,iat2)
                    dcod(1:3) = dcod(1:3) -                  & ! #SLC2 #PB
                      fxcell(1:3)*nint(dcod(1:3)*invcel(1:3))  ! #SLC2 #PB
                    dd = dcod(1)*dcod(1) + dcod(2)*dcod(2) +           &
                         dcod(3)*dcod(3)
                    if ( dd .le. rlim2 ) then
                      nn = nn + 1
                      itbEL(nn,iat1,n) = iat2
                    endif
                  enddo
                  ntbEL(iat1,n) = nn
                  nn = ntb(iat1,n)
                  do kk = border+1,num2(j)
                    iat2 = inum(ipcl(kk,j,jx,jy,jz))
                    dcod(1:3) = TTcod(1:3) - Tcord(1:3,iat2)
                    dcod(1:3) = dcod(1:3) -                  & ! #SLC2 #PB
                      fxcell(1:3)*nint(dcod(1:3)*invcel(1:3))  ! #SLC2 #PB
                    dd = dcod(1)*dcod(1) + dcod(2)*dcod(2) +           &
                         dcod(3)*dcod(3)
                    if ( dd .le. rlim2 ) then
                      nn = nn + 1
                      itb(nn,iat1,n) = iat2
                    endif
                  enddo
                  ntb(iat1,n) = nn
                enddo
              enddo
            enddo
          enddo
        enddo
        enddo
        enddo
        !$OMP end do
        !$OMP end parallel

      else
        !$OMP parallel default(none)                                   & !
        !$OMP private(iz,iy,ix,num,i,nel,ii,iat1,nce,ice,jx,jy,jz,num2,& !
        !$OMP         j,n,border,kk,iat2,dcod,dd,nn,TTcod)             & !
        !$OMP shared(nv,npcl1,nfrag,cell_nEL,ipcl,nnr,inr,n_matrix,    & !
        !$OMP        fxcell,invcel,                                    & ! ! #SLC2 #PB
        !$OMP        Tcord,rlim2,ntbEL,itbEL,ntb,itb,inum)
        !$OMP do schedule (dynamic) collapse(3)
        do iz = 1,nv(1)
        do iy = 1,nv(1)
        do ix = 1,nv(1)

          num(1:nfrag) = npcl1(1:nfrag,ix,iy,iz)
          if ( all(num .eq. 0 ) ) cycle

          do i = 1,nfrag
            if ( num(i) .eq. 0 ) cycle
            nel = cell_nEL(i,ix,iy,iz)

            ! (1a) Within a cell
            do ii = 1,num(i)
              iat1 = inum(ipcl(ii,i,ix,iy,iz))

              j = i ; n = n_matrix(j,i)
              if ( ii .le. nel ) then
                nn = ntbEL(iat1,n)
                do jj = 1,ii-1
                  iat2 = inum(ipcl(jj,j,ix,iy,iz))
                  nn = nn + 1
                  itbEL(nn,iat1,n) = iat2
                enddo
                do jj = ii+1,num(i)
                  iat2 = inum(ipcl(jj,j,ix,iy,iz))
                  nn = nn + 1
                  itbEL(nn,iat1,n) = iat2
                enddo
                ntbEL(iat1,n) = nn
              else
                nn = ntbEL(iat1,n)
                do jj = 1,nel
                  iat2 = inum(ipcl(jj,j,ix,iy,iz))
                  nn = nn + 1
                  itbEL(nn,iat1,n) = iat2
                enddo
                ntbEL(iat1,n) = nn
                nn = ntb(iat1,n)
                ! i-particle is EL+vdW, and jj>ii+1 is always EL+vdW
                do jj = nel+1,ii-1
                  iat2 = inum(ipcl(jj,j,ix,iy,iz))
                  nn = nn + 1
                  itb(nn,iat1,n) = iat2
                enddo
                do jj = ii+1,num(i)
                  iat2 = inum(ipcl(jj,j,ix,iy,iz))
                  nn = nn + 1
                  itb(nn,iat1,n) = iat2
                enddo
                ntb(iat1,n) = nn
              endif

              do j = 1,i-1
                n = n_matrix(j,i)
                if ( ii .le. nel ) then
                  border = num(j) ! all atoms are EL
                else
                  border = cell_nEL(j,ix,iy,iz)
                endif

                nn = ntbEL(iat1,n)
                do jj = 1,border
                  iat2 = inum(ipcl(jj,j,ix,iy,iz))
                  nn = nn + 1
                  itbEL(nn,iat1,n) = iat2
                enddo
                ntbEL(iat1,n) = nn

                nn = ntb(iat1,n)
                do jj = border+1,num(j)
                  iat2 = inum(ipcl(jj,j,ix,iy,iz))
                  nn = nn + 1
                  itb(nn,iat1,n) = iat2
                enddo
                ntb(iat1,n) = nn
              enddo
              do j = i+1,nfrag
                n = n_matrix(j,i)
                if ( ii .le. nel ) then
                  border = num(j) ! all atoms are EL
                else
                  border = cell_nEL(j,ix,iy,iz)
                endif

                nn = ntbEL(iat1,n)
                do jj = 1,border
                  iat2 = inum(ipcl(jj,j,ix,iy,iz))
                  nn = nn + 1
                  itbEL(nn,iat1,n) = iat2
                enddo
                ntbEL(iat1,n) = nn

                nn = ntb(iat1,n)
                do jj = border+1,num(j)
                  iat2 = inum(ipcl(jj,j,ix,iy,iz))
                  nn = nn + 1
                  itb(nn,iat1,n) = iat2
                enddo
                ntb(iat1,n) = nn
              enddo

              TTcod(1:3) = Tcord(1:3,iat1)
              nce = nnr(ix,iy,iz)
              do ice = 1,nce
                jx = inr(1,ice,ix,iy,iz)
                jy = inr(2,ice,ix,iy,iz)
                jz = inr(3,ice,ix,iy,iz)

                num2(1:nfrag) = npcl1(1:nfrag,jx,jy,jz)
                if ( all(num2 .eq. 0) ) cycle

                do j = 1,nfrag
                  n = n_matrix(j,i)
                  if ( ii .le. nel ) then
                    border = num2(j)
                  else
                    border = cell_nEL(j,jx,jy,jz)
                  endif

                  nn = ntbEL(iat1,n)
                  do kk = 1,border
                    iat2 = inum(ipcl(kk,j,jx,jy,jz))
                    dcod(1:3) = TTcod(1:3) - Tcord(1:3,iat2)
                    dcod(1:3) = dcod(1:3) -                  & ! #SLC2 #PB
                      fxcell(1:3)*nint(dcod(1:3)*invcel(1:3))  ! #SLC2 #PB
                    dd = dcod(1)*dcod(1) + dcod(2)*dcod(2) +           &
                         dcod(3)*dcod(3)
                    if ( dd .le. rlim2 ) then
                      nn = nn + 1
                      itbEL(nn,iat1,n) = iat2
                    endif
                  enddo
                  ntbEL(iat1,n) = nn
                  nn = ntb(iat1,n)
                  do kk = border+1,num2(j)
                    iat2 = inum(ipcl(kk,j,jx,jy,jz))
                    dcod(1:3) = TTcod(1:3) - Tcord(1:3,iat2)
                    dcod(1:3) = dcod(1:3) -                  & ! #SLC2 #PB
                      fxcell(1:3)*nint(dcod(1:3)*invcel(1:3))  ! #SLC2 #PB
                    dd = dcod(1)*dcod(1) + dcod(2)*dcod(2) +           &
                         dcod(3)*dcod(3)
                    if ( dd .le. rlim2 ) then
                      nn = nn + 1
                      itb(nn,iat1,n) = iat2
                    endif
                  enddo
                  ntb(iat1,n) = nn
                enddo
              enddo
            enddo
          enddo
        enddo
        enddo
        enddo
        !$OMP end do
        !$OMP end parallel

      endif

!****************************************

      return
      end subroutine cutatm_PB ! #SLC2 #PB
      end subroutine cutatm_CB ! #SLC2 #CB
