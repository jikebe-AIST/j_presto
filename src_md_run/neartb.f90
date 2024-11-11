
      subroutine neartb(ier)

!****************************************************************************
!
!     Making tables the for near-field exact calc.
!
!     NOTE:
!        The tables are the half of the compleate tables
!        to avoid double counting.
!        This means that each pair is stored in a table only once.
!
!     CONDITION:
!        This subroutine should be called after setting N of atoms in 
!        the system (i.e., ixnatm) and after calculating cell-quantities.
!
!****************************************************************************

      use COMPAR ; use COMCMM ; use COMERG
      use COMBAS ; use COMCMMC,only: cord
      use CALC_TIME
      !$ use omp_lib

      implicit none

      integer(4),intent(inout):: ier

      integer(4):: i,j,n,ii,jj,kk,ix,iy,iz,jx,jy,jz,kx,ky,kz,iat1,iat2
      integer(4):: jat1,jat2,nn,nce,ice,mm,nchk,itype
      real(8):: dcod(3)
      real(8):: dd

      integer(4):: nv1,num(nfrag),num2(nfrag)
      integer(4):: nel,nelvdw,border

!************************************************************
!     Initializing the tables.

!      call system_clock(tim1)
      nv1 = nv(1)
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

!******************************************
!  For the first table.

      if ( high_para .ge. 0 ) then
        !$OMP parallel default(none)                                   & !
        !$OMP private(num,iz,iy,ix,i,ii,n,iat1,j,nn,jj,iat2,nce,ice,jx,& !
        !$OMP         jy,jz,num2,kk,nel,border)                        & !
        !$OMP shared(nfrag,nv1,npcl1,ipcl,n_matrix,ntb,itb,nnr,inr,    & !
        !$OMP        ntbEL,itbEL,cell_nEL,inum)
        !$OMP do schedule (auto) collapse(3)
        do iz = 1,nv1
        do iy = 1,nv1
        do ix = 1,nv1
          num(1:nfrag) = npcl1(1:nfrag,ix,iy,iz)

          if ( all(num .eq. 0) ) cycle

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

              ! (2) Between the nearest cells.
              ! This part is programed to avoid the double count.
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
                    border = num2(j) ! interactions to all atoms are EL
                  else
                    border = cell_nEL(j,jx,jy,jz)
                  endif

                  nn = ntbEL(iat1,n)
                  do kk = 1,border
                    iat2 = inum(ipcl(kk,j,jx,jy,jz))
                    nn = nn + 1
                    itbEL(nn,iat1,n) = iat2
                  enddo
                  ntbEL(iat1,n) = nn
                  nn = ntb(iat1,n)
                  do kk = border+1,num2(j)
                    iat2 = inum(ipcl(kk,j,jx,jy,jz))
                    nn = nn + 1
                    itb(nn,iat1,n) = iat2
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
        !$OMP private(num,iz,iy,ix,i,ii,n,iat1,j,nn,jj,iat2,nce,ice,jx,& !
        !$OMP         jy,jz,num2,kk,nel,border)                        & !
        !$OMP shared(nfrag,nv1,npcl1,ipcl,n_matrix,ntb,itb,nnr,inr,    & !
        !$OMP        ntbEL,itbEL,cell_nEL,inum)
        !$OMP do schedule (auto) collapse(3)
        do iz = 1,nv1
        do iy = 1,nv1
        do ix = 1,nv1
          num(1:nfrag) = npcl1(1:nfrag,ix,iy,iz)

          if ( all(num .eq. 0) ) cycle

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

              ! (2) Between the nearest cells.
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
                    border = num2(j) ! interactions to all atoms are EL
                  else
                    border = cell_nEL(j,jx,jy,jz)
                  endif

                  nn = ntbEL(iat1,n)
                  do kk = 1,border
                    iat2 = inum(ipcl(kk,j,jx,jy,jz))
                    nn = nn + 1
                    itbEL(nn,iat1,n) = iat2
                  enddo
                  ntbEL(iat1,n) = nn
                  nn = ntb(iat1,n)
                  do kk = border+1,num2(j)
                    iat2 = inum(ipcl(kk,j,jx,jy,jz))
                    nn = nn + 1
                    itb(nn,iat1,n) = iat2
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

!******************************************************
!     Make an another table for additional vdW interaction Calc.
!          (for atoms with the dist. less than rlim, and in the 
!           second nearest cells to the central cell. )

      !! In the case that CUT-OFF dist. is shorter than min. box size
      if ( idelay .le. 1 ) goto 800
      !$OMP parallel default (none)                                  & !
      !$OMP private(ix,iy)                                           & !
      !$OMP shared(nfrg,ixnatm,ntbvdW)
      !$OMP do schedule (static) collapse(2)
      do iy = 1,nfrg
      do ix = 1,ixnatm
        ntbvdW(ix,iy) = 0
      enddo
      enddo
      !$OMP end do
      !$OMP end parallel
      !$OMP parallel default(none)                                   & !
      !$OMP private(iz,iy,ix,num,nce,ice,jx,jy,jz,num2,i,ii,iat1,j,n,& !
      !$OMP         jj,iat2,dcod,dd,nn)                              & !
      !$OMP shared(nv1,npcl1,nfrag,nnrv,inrv,ipcl,zero_vdW,n_matrix, & !
      !$OMP        Tcord,rlim2,ntbvdW,itbvdW,cell_nEL,inum)
      !$OMP do schedule (dynamic) collapse(3)
      do iz = 1,nv1
      do iy = 1,nv1
      do ix = 1,nv1
        num(1:nfrag) = npcl1(1:nfrag,ix,iy,iz)
        if ( all(num .eq. 0) ) cycle

        nce = nnrv(ix,iy,iz)
        do ice = 1,nce
          jx = inrv(1,ice,ix,iy,iz)
          jy = inrv(2,ice,ix,iy,iz)
          jz = inrv(3,ice,ix,iy,iz)

          num2(1:nfrag) = npcl1(1:nfrag,jx,jy,jz)
          if ( all(num2 .eq. 0) ) cycle

          do i = 1,nfrag
            do ii = cell_nEL(i,ix,iy,iz)+1,num(i)
              iat1 = inum(ipcl(ii,i,ix,iy,iz))

              do j = 1,nfrag
                n = n_matrix(j,i)
                do jj = cell_nEL(j,jx,jy,jz)+1,num2(j)
                  iat2 = inum(ipcl(jj,j,jx,jy,jz))
                  dcod(1:3) = Tcord(1:3,iat1) - Tcord(1:3,iat2)
                  dd = dcod(1)*dcod(1) + dcod(2)*dcod(2) +           &
                       dcod(3)*dcod(3)
                  if ( dd .lt. rlim2 ) then
                    nn = ntbvdW(iat1,n)
                    nn = nn + 1
                    itbvdW(nn,iat1,n) = iat2
                    ntbvdW(iat1,n) = nn
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      enddo
      enddo
      !$OMP enddo
      !$OMP end parallel

800   continue
!****************
!     Check the tables.
!      print*,"aho ",maxval(ntbEL(1:ixnatm,:)),maxval(ntb(1:ixnatm,:)), &
!        maxval(ntbvdW(1:ixnatm,:))
      n15mxEL_max = max(n15mxEL_max,maxval(ntbEL(:,:)))
      if ( n15mxEL_max .gt. n15mxEL ) then
        print*,' ??????????????????????????????????????'
        print*,' ? n15mxEL_max is too small. (neartb) ?'
        print*,' ??????????????????????????????????????'
        print*,' n15mxEL & the max.# detected  =',n15mxEL,' ',         &
               n15mxEL_max
        ier = 51
      endif
      n15mx_max = max(n15mx_max,maxval(ntb(:,:)))
      if ( n15mx_max .gt. n15mx ) then
        print*,' ??????????????????????????????????'
        print*,' ? n15mxEL is too small. (neartb) ?'
        print*,' ??????????????????????????????????'
        print*,' n15mx & the max.# detected  =',n15mx,' ',n15mx_max
        ier = 51
      endif
      if ( idelay .gt. 1 ) then
        nvdw_max = max(nvdw_max,maxval(ntbvdW(:,:)))
        if ( nvdw_max .gt. nvdw ) then
          print*,' ???????????????????????????????'
          print*,' ? nvdw is too small. (neartb) ?'
          print*,' ???????????????????????????????'
          print*,' nvdw & the max.# detected =',nvdw,'  ',nvdw_max
          ier=51
        endif
      endif
      if ( ier .eq. 51 ) return

!     Remove 1-4 interaction & Rearrangement array ntb & itb
      nchk = 0
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

      ! Remove 1-4 interaction & Rearrangement array ntbvdW & itbvdW
      j = nchk ; if ( high_para .lt. 0 ) j = nchk / 2
      if ( idelay.gt.1 .and. j.ne.nntab14 ) then
        !$OMP parallel default(none)                                 & !
        !$OMP private(iat1,n,iat2,kk,j,nn,jat1,jat2)                 & !
        !$OMP shared(ntbvdW,itbvdW,ntab14,itab14,ixnatm,n_matrix,    & !
        !$OMP        icls,inum,nchk)
        !$OMP do reduction(+ : nchk)                                 & !
        !$OMP    schedule (auto)
        do iat1 = 1,ixnatm
          jat1 = inum(iat1)
          OUTER2 : do n = 1,ntab14(iat1)
            iat2 = itab14(n,iat1)
            jat2 = inum(iat2)
            kk = n_matrix(icls(iat1),icls(iat2))
            do j = 1,ntbvdW(jat1,kk)
              if ( jat2 .eq. itbvdW(j,jat1,kk) ) then
                nn = ntbvdW(jat1,kk) 
                itbvdW(j,jat1,kk) = itbvdW(nn,jat1,kk)
                ntbvdW(jat1,kk) = nn - 1
                nchk = nchk + 1
                cycle OUTER2
              endif
            enddo
          enddo OUTER2
        enddo
        !$OMP end do
        !$OMP end parallel
      endif

      j = nchk ; if ( high_para .lt. 0 ) j = nchk / 2
      if ( j .ne. nntab14 ) then
        print*,"ERROR> NEARTB"
        print*,"CMM cell size is too small"
        ier = 99 ; return
      endif

      ! Separate hydrogen bond energy
      if ( iyeflg(10) .eq. 1 ) then
        ntbhyd = 0 ; ntbhyd2 = 0
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
              if ( iyppid(Tixatyp(iat2),itype) .eq. 2 ) then
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
              if ( iyppid(Tixatyp(iat2),itype) .eq. 2 ) then
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
          print*,' ?????????????????????????????????????????????'
          print*,' ? n15mxEL for ntbhyd is too small. (neartb) ?'
          print*,' ??????????????????????????????????????????????'
          print*,' n15mx & the max.# detected  =',n15mx,' ',n15mx_max
          ier = 51 ; return
        endif

        if ( idelay .gt. 1 ) then
          !$OMP parallel default(none)                               & !
          !$OMP private(j,iat1,itype,ii,nn,i,iat2,n)                 & !
          !$OMP shared(nfrg,ixnatm,ntbvdW,Tixatyp,itbvdW,iyppid,     & !
          !$OMP        itbhyd2,ntbhyd2)
          !$OMP do schedule (auto) collapse(2)
          do j = 1,nfrg
          do iat1 = 1,ixnatm
            if ( ntbvdW(iat1,j) .ne. 0 ) then
              itype = Tixatyp(iat1)
              ii = 0 ; nn = 0 ; n = ntbhyd2(iat1,j)
              do i = 1,ntbvdW(iat1,j)
                iat2 = itbvdW(i,iat1,j)
                if ( iyppid(Tixatyp(iat2),itype) .eq. 2 ) then
                  nn = nn + 1
                  itbhyd2(n+nn,iat1,j) = itbvdW(i,iat1,j)
                else
                  ii = ii + 1
                  itbvdW(ii,iat1,j) = itbvdW(i,iat1,j)
                endif
              enddo
              ntbhyd2(iat1,j) = nn + n
              ntbvdW(iat1,j) = ntbvdW(iat1,j) - nn
            endif
          enddo
          enddo
          !$OMP end do
          !$OMP end parallel

          ! Check the tables.
          if ( idelay .gt. 1 ) then
            nvdw_max = max(nvdw_max,maxval(ntbhyd2(:,:)))
            if ( nvdw_max .gt. nvdw ) then
              print*,' ???????????????????????????????????????????'
              print*,' ? nvdw for ntbhyd2 is too small. (neartb) ?'
              print*,' ???????????????????????????????????????????'
              print*,' nvdw & the max.# detected =',nvdw,'  ',nvdw_max
              ier=51 ; return
            endif
          endif
        endif
      endif

!********************
!     Make CMM table

      ! My cell
      Nmycel(1:nfrag) = 0
      do iz = 1,nv1
      do iy = 1,nv1
      do ix = 1,nv1
        num(1:nfrag) = npcl1(1:nfrag,ix,iy,iz)
        do i = 1,nfrag
          if ( num(i) .ne. 0 ) then
            n = Nmycel(i) + 1
            mycel(1:3,n,i) = (/ix,iy,iz/)
            Nmycel(i) = n
          endif
        enddo
      enddo
      enddo
      enddo

      !$OMP parallel default(none)                                   & !
      !$OMP private(j,ii,ix,iy,iz,nce,i,jx,jy,jz,jj,nn)              & !
      !$OMP shared(nfrag,Nmycel,mycel,ncm1,icm1,npcl1,Nyrcel,yrcel,  & !
      !$OMP        lvup,nlev,ncm2,icm2,npcl2,ncm3,icm3,npcl3,ncm4,   & !
      !$OMP        icm4,npcl4,ncm5,icm5,npcl5)
      do j = 1,nfrag
        ! Level 1
        !$OMP do
        do ii = 1,Nmycel(j)
          ix = mycel(1,ii,j) ; iy = mycel(2,ii,j) ; iz = mycel(3,ii,j)
          nce = ncm1(ix,iy,iz) ; Nyrcel(ii,1,1:nfrag,j) = 0
          do i = 1,nce
            jx = icm1(1,i,ix,iy,iz)
            jy = icm1(2,i,ix,iy,iz)
            jz = icm1(3,i,ix,iy,iz)
            do jj = 1,nfrag
              if ( npcl1(jj,jx,jy,jz) .ne. 0 ) then
                nn = Nyrcel(ii,1,jj,j) + 1
                yrcel(1:3,nn,ii,1,jj,j) = (/jx,jy,jz/)
                Nyrcel(ii,1,jj,j) = nn
              endif
            enddo
          enddo
        enddo
        !$OMP end do nowait

        ! Level 2
        !$OMP do
        do ii = 1,Nmycel(j)
          ix = lvup(mycel(1,ii,j))
          iy = lvup(mycel(2,ii,j))
          iz = lvup(mycel(3,ii,j))
          nce = ncm2(ix,iy,iz) ; Nyrcel(ii,2,1:nfrag,j) = 0
          do i = 1,nce
            jx = icm2(1,i,ix,iy,iz)
            jy = icm2(2,i,ix,iy,iz)
            jz = icm2(3,i,ix,iy,iz)
            do jj = 1,nfrag
              if ( npcl2(jj,jx,jy,jz) .ne. 0 ) then
                nn = Nyrcel(ii,2,jj,j) + 1
                yrcel(1:3,nn,ii,2,jj,j) = (/jx,jy,jz/)
                Nyrcel(ii,2,jj,j) = nn
              endif
            enddo
          enddo
        enddo
        !$OMP end do nowait

        if ( nlev .lt. 3 ) cycle
        ! Level 3
        !$OMP do
        do ii = 1,Nmycel(j)
          ix = lvup(lvup(mycel(1,ii,j)))
          iy = lvup(lvup(mycel(2,ii,j)))
          iz = lvup(lvup(mycel(3,ii,j)))
          nce = ncm3(ix,iy,iz) ; Nyrcel(ii,3,1:nfrag,j) = 0
          do i = 1,nce
            jx = icm3(1,i,ix,iy,iz)
            jy = icm3(2,i,ix,iy,iz)
            jz = icm3(3,i,ix,iy,iz)
            do jj = 1,nfrag
              if ( npcl3(jj,jx,jy,jz) .ne. 0 ) then
                nn = Nyrcel(ii,3,jj,j) + 1
                yrcel(1:3,nn,ii,3,jj,j) = (/jx,jy,jz/)
                Nyrcel(ii,3,jj,j) = nn
              endif
            enddo
          enddo
        enddo
        !$OMP end do nowait

        if ( nlev .lt. 4 ) cycle
        ! Level 4
        !$OMP do
        do ii = 1,Nmycel(j)
          ix = lvup(lvup(lvup(mycel(1,ii,j))))
          iy = lvup(lvup(lvup(mycel(2,ii,j))))
          iz = lvup(lvup(lvup(mycel(3,ii,j))))
          nce = ncm4(ix,iy,iz) ; Nyrcel(ii,4,1:nfrag,j) = 0
          do i = 1,nce
            jx = icm4(1,i,ix,iy,iz)
            jy = icm4(2,i,ix,iy,iz)
            jz = icm4(3,i,ix,iy,iz)
            do jj = 1,nfrag
              if ( npcl4(jj,jx,jy,jz) .ne. 0 ) then
                nn = Nyrcel(ii,4,jj,j) + 1
                yrcel(1:3,nn,ii,4,jj,j) = (/jx,jy,jz/)
                Nyrcel(ii,4,jj,j) = nn
              endif
            enddo
          enddo
        enddo
        !$OMP end do nowait

        if ( nlev .lt. 5 ) cycle
        ! Level 5
        !$OMP do
        do ii = 1,Nmycel(j)
          ix = lvup(lvup(lvup(lvup(mycel(1,ii,j)))))
          iy = lvup(lvup(lvup(lvup(mycel(2,ii,j)))))
          iz = lvup(lvup(lvup(lvup(mycel(3,ii,j)))))
          nce = ncm5(ix,iy,iz) ; Nyrcel(ii,5,1:nfrag,j) = 0
          do i = 1,nce
            jx = icm5(1,i,ix,iy,iz)
            jy = icm5(2,i,ix,iy,iz)
            jz = icm5(3,i,ix,iy,iz)
            do jj = 1,nfrag
              if ( npcl5(jj,jx,jy,jz) .ne. 0 ) then
                nn = Nyrcel(ii,5,jj,j) + 1
                yrcel(1:3,nn,ii,5,jj,j) = (/jx,jy,jz/)
                Nyrcel(ii,5,jj,j) = nn
              endif
            enddo
          enddo
        enddo
        !$OMP end do nowait
      enddo
      !$OMP end parallel

!********************
!     Make table for nearest cell calc. of external CMM

      if ( extCMM_flag ) then
        ! My cell for nearest cell calc. of External Cmm
        Nmycel_nEC = 0
        do iz = 1,nv1
        jz = max(iz-1,1) ; kz = min(iz+1,nv1)
        do iy = 1,nv1
        jy = max(iy-1,1) ; ky = min(iy+1,nv1)
        do ix = 1,nv1
          num(1) = sum(npcl1(1:nfrag,ix,iy,iz))
          if ( num(1) .ne. 0 ) then
          jx = max(ix-1,1) ; kx = min(ix+1,nv1)
          if ( any(exflag1(jx:kx,jy:ky,jz:kz)) ) then
            Nmycel_nEC = Nmycel_nEC + 1
            mycel_nEC(1:3,Nmycel_nEC) = (/ix,iy,iz/)
          endif
          endif
        enddo
        enddo
        enddo

        ! Make your cell for nearest cell calc. of External Cmm
        Nyrcel_nEC(1:Nmycel_nEC) = 0
        !$OMP parallel default(none)                                 & !
        !$OMP private(ii,ix,iy,iz,jz,jy,jx,nn)                       & !
        !$OMP shared(Nmycel_nEC,mycel_nEC,exflag1,Nyrcel_nEC,yrcel_nEC)
        !$OMP do
        do ii = 1,Nmycel_nEC
          ix= mycel_nEC(1,ii) ; iy= mycel_nEC(2,ii) ;iz= mycel_nEC(3,ii)
          do jz = iz-1,iz+1
          do jy = iy-1,iy+1
          do jx = ix-1,ix+1
            if ( exflag1(jx,jy,jz) ) then
              nn = Nyrcel_nEC(ii) + 1
              yrcel_nEC(1:3,nn,ii) = (/jx,jy,jz/)
              Nyrcel_nEC(ii) = nn
            endif
          enddo
          enddo
          enddo
        enddo
        !$OMP end do
        !$OMP end parallel
      endif
!      call system_clock(tim2)
!      tmptime = tmptime + tim2 - tim1

!************************************************************

      return
      end subroutine neartb
