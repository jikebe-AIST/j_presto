
      subroutine incel(iprint)

!*************************************************

      use COMPAR ; use COMCMM ; use COMERG ; use COMBAS ; use COMCMMC
      !$ use omp_lib

      implicit none

      ! unit number of standard output
        integer(4),intent(in):: iprint

      integer(4):: i,j,k,m,ix,iy,iz,jx1,jx2,jy1,jy2,jz1,jz2,itmp,jx,jy,&
        jz,ipair(2),Tpara(2,300),mx1,mx2,my1,my2,mz1,mz2,mxdelay,      &
        mydelay,mzdelay,ii,jj,kk,iix,iiy,iiz
      logical(4):: flgx,flgy,flgz
      logical(4),allocatable:: flg(:),flg2(:)
      integer(4),allocatable:: Tiyp(:,:),Tiytnbf(:),Tcls(:),T14cls(:)
      real(8),allocatable:: Tfyq(:),Tfyf(:),Tfyrot(:),Tfyphs(:),       &
        Tiydiv(:),Tvws(:),Tess(:),dwork(:,:)
      real(8):: rtmp
      real(4):: rix,riy,riz
      logical(4),allocatable:: tcel(:,:,:,:,:,:)

!************************************************************
!     allocation of arrays

      ! Change cell-size for ACT & ZD
      if ( iy15m.eq.2 .or. iy15m.eq.4 ) then
        rlim = fycutl + 1.5d0
        irlim2 = 1.d0 / (fycutl*fycutl)
        cminsiz = rlim * 0.5d0
!        cminsiz = rlim / sqrt(3.d0)
!        cminsiz = int(cminsiz*10.d0) * 0.1d0
      else
        rlim = fycutl
      endif
      rlim2 = rlim * rlim

      ! Change nlev, if the input nlev is too small
      rtmp = maxval(cord(1,:)) - minval(cord(1,:))
      rtmp = max(maxval(cord(2,:))-minval(cord(2,:)),rtmp)
      rtmp = max(maxval(cord(3,:))-minval(cord(3,:)),rtmp)
      if ( dble(2**(nlev+1))*cminsiz .lt. rtmp ) then
        nlev = ceiling(log(rtmp/cminsiz)/log(2.d0)-1.d0)
        if ( nlev .gt. 5 ) then
          write(6,*)"!! CAUTION !! (incel.f90)"
          write(6,*)"This program can apply only nlev <= 5"
          write(6,*)"Program STOP" ; stop
        endif
      endif

      allocate(nv(nlev),siz(nlev))
      forall(i=0:nlev-1) nv(nlev-i) = 4 * 2**(i)
      forall(i=1:nlev) siz(i) = cminsiz*2**(i-1)
      idelay = ceiling(rlim/siz(1))

      if ( iy15m .eq. 1 ) then
        itwin = 1
        allocate(lvdn(2,nv(2)),lvup(nv(1)))
        allocate(gwk(3,ixnatm,nfrg),Ecmm(nlev,nfrg))
        !$OMP parallel default (none)                                & !
        !$OMP private(i,j)                                           & !
        !$OMP shared(gwk,ixnatm,nfrg)
        !$OMP do schedule (static) collapse(2)
        do j = 1,nfrg
        do i = 1,ixnatm
          gwk(:,i,j) = 0.d0
        enddo
        enddo
        !$OMP end do
        !$OMP end parallel
        Ecmm(1:nlev,1:nfrg) = 0.d0
        i = nv(1) * nv(1) * nv(1)
        allocate(Nmycel(nfrag),mycel(3,i,nfrag))
        allocate(Nyrcel(i,nlev,nfrag,nfrag))
        allocate(yrcel(3,jcmx,i,nlev,nfrag,nfrag))
        allocate(CMMpot1(13,nfrag,nv(1),nv(1),nv(1)))
        allocate(CMMpot2(13,nfrag,nv(2),nv(2),nv(2)))
        if ( nlev .ge. 3 ) then
          allocate(CMMpot3(13,nfrag,nv(3),nv(3),nv(3)))
        if ( nlev .ge. 4 ) then
          allocate(CMMpot4(13,nfrag,nv(4),nv(4),nv(4)))
        if ( nlev .ge. 5 ) then
          allocate(CMMpot5(13,nfrag,nv(5),nv(5),nv(5)))
        endif
        endif
        endif
        allocate(npcl2(nfrag,nv(2),nv(2),nv(2)))
        if ( nlev .ge. 3 ) then
          allocate(npcl3(nfrag,nv(3),nv(3),nv(3)))
        if ( nlev .ge. 4 ) then
          allocate(npcl4(nfrag,nv(4),nv(4),nv(4)))
        if ( nlev .ge. 5 ) then
          allocate(npcl5(nfrag,nv(5),nv(5),nv(5)))
        endif
        endif
        endif
      endif

      ! openMP parallel ways
      if ( high_para .eq. -1 ) then
        if ( nlev .le. 3 ) then
          if ( nomp .le. 13 ) then
            high_para = 0
          else
            high_para = 2
          endif
        else
          high_para = 1
        endif
      endif
      write(iprint,*)"     PARALLELAZATION LEVEL of OpenMP : "
      select case (high_para)
      case (0)
        write(iprint,*)"       REDUCTION"
      case (1)
        write(iprint,*)"       HIGH (not using reduction)"
        id2p1 = 2*idelay+1 ; id2p1_2 = id2p1*id2p1
        high_para = id2p1*id2p1*(idelay+1)
        allocate(np_cell(high_para))
        i = nv(1)*nv(1)*nv(1) / high_para + 1
        allocate(para_se(2,i,high_para))
      case (2)
        write(iprint,*)"       DOUBLE"
        high_para = -1
      end select
      if ( i_vec .eq. 0 ) then
        SIMD_chk = .true. ; i_vec = 1
      else
       write(iprint,'(a,i0)')"      SIMD PARALLELAZATION = ",i_vec
        SIMD_chk = .false.
      endif
      write(iprint,*)""

      if ( iy15m .eq. 1 ) then
        rtmp = (4*siz(1))**3
      else
        rtmp = 4.d0 / 3.d0 * 3.14159265358979d0 * (fycutl)**3
      endif
!      if ( high_para .ge. 0 ) i = ceiling(i*0.5)
      i = 1
      if ( high_para .eq. -1 ) i = 2
      allocate(npcl1(nfrag,nv(1),nv(1),nv(1)))
      ipmax = ceiling((siz(1)**3)*COEipmax)*i
      allocate(ipcl(ipmax,nfrag,nv(1),nv(1),nv(1)))
      n15mx = ceiling(COEn15mx * rtmp)*i
      allocate(ntb(ixnatm,nfrg),itb(n15mx,ixnatm,nfrg))
      n15mxEL = ceiling(COEn15mxEL * rtmp)*i
      allocate(ntbEL(ixnatm,nfrg),itbEL(n15mxEL,ixnatm,nfrg))
      nvdw = 0
      if ( iy15m .eq. 1 ) then
        nvdw = ceiling(COEnvdw * rtmp)*i
        allocate(ntbvdW(ixnatm,nfrg),itbvdW(nvdw,ixnatm,nfrg))
        !$OMP parallel default (none)                                & !
        !$OMP private(i,j)                                           & !
        !$OMP shared(ntbvdW,ixnatm,nfrg)
        !$OMP do schedule (static) collapse(2)
        do j = 1,nfrg
        do i = 1,ixnatm
          ntbvdW(i,j) = 0
        enddo
        enddo
        !$OMP end do
        !$OMP end parallel
      endif
      allocate(cell_nEL(nfrag,nv(1),nv(1),nv(1)))
      allocate(inum(ixnatm),Tixatyp(ixnatm),Tcord(4,ixnatm))
      if ( maxval(iyppid) .ge. 2 ) then
        iyeflg(10) = 1
        allocate(ntbhyd(ixnatm,nfrg),itbhyd(n15mx,ixnatm,nfrg))
        allocate(ntbhyd2(ixnatm,nfrg),itbhyd2(nvdw,ixnatm,nfrg))
      endif

!************************************************************
!  General quantities for cells.

      write(iprint,*)"     SIZE OF CELL BOX"
      do i = 1,nlev
        write(iprint,'(a,i0,a,f8.3)')"       LEVEL-",i,": ",siz(i)
      enddo
      write(iprint,'(a,i0)')"      idelay = ",idelay
      write(iprint,*)""

      ! For CMM
      if ( iy15m .eq. 1 ) then
        !  Level lists.
        ! lvdn
        do i = 1,nv(2)
          j = 2 * i
          lvdn(1,i) = j - 1 ; lvdn(2,i) = j
        enddo
        ! lvup
        do i = 1,nv(1)
          j = int(i/2) ; k = i - j * 2
          select case (k)
            case(0)
              lvup(i) = j
            case(1)
              lvup(i) = j + 1
          end select
        enddo

        !  Nearest far cells for each level.
        do i = 1,nlev
          call mk_cm(i)
        enddo

        ! For cells for vdW
        !  Near cells for vdw
        if ( idelay .gt. 1 ) then
          if ( high_para .ge. 0 ) then
!            i = (8*idelay*idelay + 36*idelay + 54) * idelay / 2
            i = (4*idelay*idelay + 6*idelay + 3) * idelay - 13
            allocate(nnrv(nv(1),nv(1),nv(1)),                        &
                     inrv(3,i,nv(1),nv(1),nv(1)))

            do iz = 1,nv(1)
              jz1 = iz ; jz2 = min(nv(1),iz+idelay)
            do iy = 1,nv(1)
              jy1 = max(1,iy-idelay) ; jy2 = min(nv(1),iy+idelay)
            do ix = 1,nv(1)
              jx1 = max(1,ix-idelay) ; jx2 = min(nv(1),ix+idelay)

              itmp = 0
              do i = jz1,jz2
                select case (iz-i)
                  case (-1:1)
                    flgz = .true.
                  case default
                    flgz = .false.
                end select
              do j = jy1,jy2
                select case (iy-j)
                  case (-1:1)
                    flgy = .true.
                  case default
                    flgy = .false.
                end select
              do k = jx1,jx2
                select case (ix-k)
                  case (-1:1)
                    flgx = .true.
                  case default
                    flgx = .false.
                end select

                if ( flgx .and. flgy .and. flgz )  cycle
                if ( (i.eq.iz .and. j.lt.iy) .or.                    &
                   (i.eq.iz .and. j.eq.iy .and. k.le.ix) ) cycle
                itmp = itmp + 1
                inrv(1:3,itmp,ix,iy,iz) = (/k,j,i/)
              enddo
              enddo
              enddo
              nnrv(ix,iy,iz) = itmp
            enddo
            enddo
            enddo
          else
            i = (8*idelay*idelay + 12*idelay + 6) * idelay - 26
            allocate(nnrv(nv(1),nv(1),nv(1)),                        &
                     inrv(3,i,nv(1),nv(1),nv(1)))

            do iz = 1,nv(1)
              jz1 = max(1,iz-idelay) ; jz2 = min(nv(1),iz+idelay)
            do iy = 1,nv(1)
              jy1 = max(1,iy-idelay) ; jy2 = min(nv(1),iy+idelay)
            do ix = 1,nv(1)
              jx1 = max(1,ix-idelay) ; jx2 = min(nv(1),ix+idelay)

              itmp = 0
              do i = jz1,jz2
                select case (iz-i)
                  case (-1:1)
                    flgz = .true.
                  case default
                    flgz = .false.
                end select
              do j = jy1,jy2
                select case (iy-j)
                  case (-1:1)
                    flgy = .true.
                  case default
                    flgy = .false.
                end select
              do k = jx1,jx2
                select case (ix-k)
                  case (-1:1)
                    flgx = .true.
                  case default
                    flgx = .false.
                end select

                if ( flgx .and. flgy .and. flgz )  cycle
                itmp = itmp + 1
                inrv(1:3,itmp,ix,iy,iz) = (/k,j,i/)
              enddo
              enddo
              enddo
              nnrv(ix,iy,iz) = itmp
            enddo
            enddo
            enddo
          endif
        endif
      endif

      !  Coordinates of cells at different levels.
      do i = 1,nlev
        call mk_c(i)
      enddo
      allocate(cg1(3,nv(1))) ; cg1 = c1
      allocate(cg2(3,nv(2))) ; cg2 = c2
      if ( nlev .ge. 3 ) then
        allocate(cg3(3,nv(3))) ; cg3 = c3
      if ( nlev .ge. 4 ) then
        allocate(cg4(3,nv(4))) ; cg4 = c4
      if ( nlev .ge. 5 ) then
        allocate(cg5(3,nv(5))) ; cg5 = c5
      endif
      endif
      endif

      ! Near cells for each level.
      select case (iy15m)
      case (1)
        i = icmx
        if ( high_para .lt. 0 ) i = 2*i
      case default
        select case (ixfbou)
        case (1)
          i = nint(nv(1)*0.5d0)
          mx1 = ceiling(i - 0.5d0*fxcell(1)/siz(1))
          mx2 = ceiling(i + 0.5d0*fxcell(1)/siz(1))
          my1 = ceiling(i - 0.5d0*fxcell(2)/siz(1))
          my2 = ceiling(i + 0.5d0*fxcell(2)/siz(1))
          mz1 = ceiling(i - 0.5d0*fxcell(3)/siz(1))
          mz2 = ceiling(i + 0.5d0*fxcell(3)/siz(1))
          if ( ceiling(0.5d0*fxcell(1)/siz(1))-0.5d0*fxcell(1)/siz(1)  &
               .gt. 0.5d0 ) then
            mxdelay = 1
          else
            mxdelay = 2
          endif
          if ( ceiling(0.5d0*fxcell(2)/siz(1))-0.5d0*fxcell(2)/siz(1)  &
               .gt. 0.5d0 ) then
            mydelay = 1
          else
            mydelay = 2
          endif
          if ( ceiling(0.5d0*fxcell(3)/siz(1))-0.5d0*fxcell(3)/siz(1)  &
               .gt. 0.5d0 ) then
            mzdelay = 1
          else
            mzdelay = 2
          endif
          i = (2*idelay+1+mxdelay)*(2*idelay+1+mydelay)*               &
              (2*idelay+1+mzdelay) - 1
        case default    
          i = (2*idelay+1)**3 - 1
          if ( high_para .ge. 0 ) i = i / 2
        end select
      end select
      allocate(nnr(nv(1),nv(1),nv(1)),inr(3,i,nv(1),nv(1),nv(1)))

      if ( iy15m .eq. 1 ) then
        m = 1
      else
        m = idelay
      endif

      if ( high_para .ge. 0 ) then
        select case (ixfbou)
        case (1)
          allocate(tcel(mx1:mx2,my1:my2,mz1:mz2,mx1:mx2,my1:my2,       &
                        mz1:mz2))
          tcel(:,:,:,:,:,:) = .false.
          iiz = mz2-mz1+1 ; iiy = my2-my1+1 ; iix = mx2-mx1+1
          riz = 1.0/iiz ; riy = 1.0/iiy ; rix = 1.0/iix
          do iz = mz1,mz2
            jz1 = iz ; jz2 = iz + m
            if ( jz2 .ge. mz2 ) jz2 = jz2 + mzdelay
          do iy = my1,my2
            jy1 = iy - m ; jy2 = iy + m
            if ( jy1 .le. my1 ) jy1 = jy1 - mydelay
            if ( jy2 .ge. my2 ) jy2 = jy2 + mydelay
          do ix = mx1,mx2
            jx1 = ix - m ; jx2 = ix + m
            if ( jx1 .le. mx1 ) jx1 = jx1 - mxdelay
            if ( jx2 .ge. mx2 ) jx2 = jx2 + mxdelay
            do i = jz1,jz2
              ii = i - iiz * floor(real(i-mz1)*riz)
            do j = jy1,jy2
              jj = j - iiy * floor(real(j-my1)*riy)
            do k = jx1,jx2
              kk = k - iix * floor(real(k-mx1)*rix)
              tcel(kk,jj,ii,ix,iy,iz) = .true.
            enddo
            enddo
            enddo
            tcel(ix,iy,iz,ix,iy,iz) = .false.
          enddo
          enddo
          enddo
          nnr(:,:,:) = 0
          do iz = mz1,mz2
          do iy = my1,my2
          do ix = mx1,mx2
            itmp = 0
            do i = mz1,mz2
            do j = my1,my2
            do k = mx1,mx2
              if ( tcel(k,j,i,ix,iy,iz) ) then
                itmp = itmp + 1
                inr(1:3,itmp,ix,iy,iz) = (/k,j,i/)
                tcel(ix,iy,iz,k,j,i) = .false.
              endif
            enddo
            enddo
            enddo
            nnr(ix,iy,iz) = itmp
          enddo
          enddo
          enddo
          deallocate(tcel)

        case default
          do iz = 1,nv(1)
            jz1 = iz ; jz2 = min(nv(1),iz+m)
          do iy = 1,nv(1)
            jy1 = max(1,iy-m) ; jy2 = min(nv(1),iy+m)
          do ix = 1,nv(1)
            jx1 = max(1,ix-m) ; jx2 = min(nv(1),ix+m)

            itmp = 0
            do i = jz1,jz2
            do j = jy1,jy2
            do k = jx1,jx2
              if ( (i.eq.iz .and. j.lt.iy) .or.                        &
                   (i.eq.iz .and. j.eq.iy .and. k.le.ix) ) cycle
              itmp = itmp + 1
              inr(1:3,itmp,ix,iy,iz) = (/k,j,i/)
            enddo
            enddo
            enddo
            nnr(ix,iy,iz) = itmp
          enddo
          enddo
          enddo
        end select

      else

        select case (ixfbou)
        case (1)
          allocate(tcel(mx1:mx2,my1:my2,mz1:mz2,mx1:mx2,my1:my2,       &
                        mz1:mz2))
          tcel(:,:,:,:,:,:) = .false.
          iiz = mz2-mz1+1 ; iiy = my2-my1+1 ; iix = mx2-mx1+1
          riz = 1.0/iiz ; riy = 1.0/iiy ; rix = 1.0/iix
          do iz = mz1,mz2
            jz1 = iz - m ; jz2 = iz + m
            if ( jz1 .le. mz1 ) jz1 = jz1 - mzdelay
            if ( jz2 .ge. mz2 ) jz2 = jz2 + mzdelay
          do iy = my1,my2
            jy1 = iy - m ; jy2 = iy + m
            if ( jy1 .le. my1 ) jy1 = jy1 - mydelay
            if ( jy2 .ge. my2 ) jy2 = jy2 + mydelay
          do ix = mx1,mx2
            jx1 = ix - m ; jx2 = ix + m
            if ( jx1 .le. mx1 ) jx1 = jx1 - mxdelay
            if ( jx2 .ge. mx2 ) jx2 = jx2 + mxdelay
            do i = jz1,jz2
              ii = i - iiz * floor(real(i-mz1)*riz)
            do j = jy1,jy2
              jj = j - iiy * floor(real(j-my1)*riy)
            do k = jx1,jx2
              kk = k - iix * floor(real(k-mx1)*rix)
              tcel(kk,jj,ii,ix,iy,iz) = .true.
            enddo
            enddo
            enddo
            tcel(ix,iy,iz,ix,iy,iz) = .false.
          enddo
          enddo
          enddo
          nnr(:,:,:) = 0
          do iz = mz1,mz2
          do iy = my1,my2
          do ix = mx1,mx2
            itmp = 0
            do i = mz1,mz2
            do j = my1,my2
            do k = mx1,mx2
              if ( tcel(k,j,i,ix,iy,iz) ) then
                itmp = itmp + 1
                inr(1:3,itmp,ix,iy,iz) = (/k,j,i/)
              endif
            enddo
            enddo
            enddo
            nnr(ix,iy,iz) = itmp
          enddo
          enddo
          enddo
          deallocate(tcel)

        case default
          do iz = 1,nv(1)
            jz1 = max(1,iz-m) ; jz2 = min(nv(1),iz+m)
          do iy = 1,nv(1)
            jy1 = max(1,iy-m) ; jy2 = min(nv(1),iy+m)
          do ix = 1,nv(1)
            jx1 = max(1,ix-m) ; jx2 = min(nv(1),ix+m)

            itmp = 0
            do i = jz1,jz2
            do j = jy1,jy2
            do k = jx1,jx2
              if ( i.eq.iz .and. j.eq.iy .and. k.eq.ix ) cycle
              itmp = itmp + 1
              inr(1:3,itmp,ix,iy,iz) = (/k,j,i/)
            enddo
            enddo
            enddo
            nnr(ix,iy,iz) = itmp
          enddo
          enddo
          enddo
        end select
      endif

      ! For repustion
      if ( iyeflg(15) .eq. 1 ) then
        allocate(nrepcl(nv(1),nv(1),nv(1),2),                          &
                 irepcl(ipmax,nv(1),nv(1),nv(1),2),ntbrep(Nrep(2)),    &
                 itbrep(n15mx,Nrep(2)),nrepnr(nv(1),nv(1),nv(1)))
        select case (ixfbou)
        case (1)
          i = (2*idelay+1+mxdelay)*(2*idelay+1+mydelay)*               &
              (2*idelay+1+mzdelay) - 1
          allocate(irepnr(3,i,nv(1),nv(1),nv(1)))
          allocate(tcel(mx1:mx2,my1:my2,mz1:mz2,mx1:mx2,my1:my2,       &
                        mz1:mz2))
          tcel(:,:,:,:,:,:) = .false.
          iiz = mz2-mz1+1 ; iiy = my2-my1+1 ; iix = mx2-mx1+1
          riz = 1.0/iiz ; riy = 1.0/iiy ; rix = 1.0/iix
          do iz = mz1,mz2
            jz1 = iz - m ; jz2 = iz + m
            if ( jz1 .le. mz1 ) jz1 = jz1 - mzdelay
            if ( jz2 .ge. mz2 ) jz2 = jz2 + mzdelay
          do iy = my1,my2
            jy1 = iy - m ; jy2 = iy + m
            if ( jy1 .le. my1 ) jy1 = jy1 - mydelay
            if ( jy2 .ge. my2 ) jy2 = jy2 + mydelay
          do ix = mx1,mx2
            jx1 = ix - m ; jx2 = ix + m
            if ( jx1 .le. mx1 ) jx1 = jx1 - mxdelay
            if ( jx2 .ge. mx2 ) jx2 = jx2 + mxdelay
            do i = jz1,jz2
              ii = i - iiz * floor(real(i-mz1)*riz)
            do j = jy1,jy2
              jj = j - iiy * floor(real(j-my1)*riy)
            do k = jx1,jx2
              kk = k - iix * floor(real(k-mx1)*rix)
              tcel(kk,jj,ii,ix,iy,iz) = .true.
            enddo
            enddo
            enddo
            tcel(ix,iy,iz,ix,iy,iz) = .false.
          enddo
          enddo
          enddo
          nrepnr(:,:,:) = 0
          do iz = mz1,mz2
          do iy = my1,my2
          do ix = mx1,mx2
            itmp = 0
            do i = mz1,mz2
            do j = my1,my2
            do k = mx1,mx2
              if ( tcel(k,j,i,ix,iy,iz) ) then
                itmp = itmp + 1
                irepnr(1:3,itmp,ix,iy,iz) = (/k,j,i/)
              endif
            enddo
            enddo
            enddo
            nrepnr(ix,iy,iz) = itmp
          enddo
          enddo
          enddo
          deallocate(tcel)

        case default
          i = (2*idelay+1)**3 - 1
          allocate(irepnr(3,i,nv(1),nv(1),nv(1)))
          do iz = 1,nv(1)
            jz1 = max(1,iz-m) ; jz2 = min(nv(1),iz+m)
          do iy = 1,nv(1)
            jy1 = max(1,iy-m) ; jy2 = min(nv(1),iy+m)
          do ix = 1,nv(1)
            jx1 = max(1,ix-m) ; jx2 = min(nv(1),ix+m)

            itmp = 0
            do i = jz1,jz2
            do j = jy1,jy2
            do k = jx1,jx2
              if ( i.eq.iz .and. j.eq.iy .and. k.eq.ix ) cycle
              itmp = itmp + 1
              irepnr(1:3,itmp,ix,iy,iz) = (/k,j,i/)
            enddo
            enddo
            enddo
            nrepnr(ix,iy,iz) = itmp
          enddo
          enddo
          enddo
        end select

      endif

!******************************

      ! For ZD
      if ( iy15m .eq. 4 ) then
        !! Calc. excess energy
        allocate(dwork(nfrg,ixnatm)) ; dwork(:,:) = 0.d0
        do i = 1,ixnatm
          j = n_matrix(icls(i),icls(i))
          dwork(j,i) = dwork(j,i) + fxchrg(i)
        enddo
        do i = 1,iynbnd
          ipair(1:2) = iypbnd(1:2,i)
          j = n_matrix(icls(ipair(1)),icls(ipair(2)))
          dwork(j,ipair(1)) = dwork(j,ipair(1)) + fxchrg(ipair(2))
          dwork(j,ipair(2)) = dwork(j,ipair(2)) + fxchrg(ipair(1))
        enddo
        do i = 1,iynang
          ipair(1:2) = (/iypang(1,i),iypang(3,i)/)
          j = n_matrix(icls(ipair(1)),icls(ipair(2)))
          dwork(j,ipair(1)) = dwork(j,ipair(1)) + fxchrg(ipair(2))
          dwork(j,ipair(2)) = dwork(j,ipair(2)) + fxchrg(ipair(1))
        enddo  
        do i = 1,iyntor
          if ( iytnbf(i) .ne. 1 ) cycle
          ipair(1:2) = (/iyptor(1,i),iyptor(4,i)/)
!!!          j = n_matrix(icls(iyptor(2,i)),icls(iyptor(3,i)))
          j = n_matrix(icls(ipair(1)),icls(ipair(2)))
          dwork(j,ipair(1)) = dwork(j,ipair(1)) + fxchrg(ipair(2))
          dwork(j,ipair(2)) = dwork(j,ipair(2)) + fxchrg(ipair(1))
        enddo
        rtmp = 1.d0 / fycutl
        zcore = 1.5d0 * rtmp
        fcoeff = rtmp * rtmp * rtmp
        bcoeff = 0.5d0 * fcoeff
        allocate(dself(nfrg))
        forall(i=1:nfrg)                                           &
          dself(i) = sum(fxchrg(1:ixnatm)*dwork(i,1:ixnatm))
        dself(:) = -zcore*dself(:)*0.5d0*ccoef
        deallocate(dwork)

        !! Make an array for excess energy calculation
        nexcess = iynbnd + iynang + iyntor
        allocate(Tiyp(2,nexcess),Tcls(nexcess))
        Tiyp(1:2,1:iynbnd) = iypbnd(1:2,1:iynbnd)
        do i = 1,iynbnd
          Tcls(i) = n_matrix(icls(Tiyp(1,i)),icls(Tiyp(2,i)))
        enddo
        do i = 1,iynang
         Tiyp(1:2,i+iynbnd) = (/iypang(1,i),iypang(3,i)/)
         Tcls(i+iynbnd) =                                              &
           n_matrix(icls(Tiyp(1,i+iynbnd)),icls(Tiyp(2,i+iynbnd)))
        enddo
        j = iynbnd+iynang ; k = 0
        do i = 1,iyntor
          if ( iytnbf(i) .eq. 1 ) then
            k = k + 1
            Tiyp(1:2,j+k) = (/iyptor(1,i),iyptor(4,i)/)
!!!            Tcls(j+k) = n_matrix(icls(iyptor(2,i)),icls(iyptor(3,i)))
            Tcls(j+k) = n_matrix(icls(iyptor(1,i)),icls(iyptor(4,i)))
          endif
        enddo
        nexcess = j + k ; allocate(iexcess(2,nexcess))
        iexcess(1:2,1:nexcess) = Tiyp(1:2,1:nexcess)
        deallocate(Tiyp)
      endif

!*******************************

      if ( nomp .ne. 1 ) then
        ! For bond
        nbnd = iynbnd - iyndbn ; np_bnd = 0 ; k = 0
        allocate(Tiyp(2,nbnd),Tfyq(nbnd),Tfyf(nbnd),bndcls(nbnd),      &
                 flg(ixnatm),flg2(nbnd))
        flg2(:) = .false.
        do
          flg(:) = .true. ; np_bnd = np_bnd + 1
          Tpara(1,np_bnd) = k + 1
          do i = 1,nbnd
            if ( flg2(i) ) cycle
            if ( all(flg(iypbnd(1:2,i))) ) then
              k = k + 1 ; flg2(i) = .true. ; flg(iypbnd(1:2,i)) = .false.
              Tiyp(:,k) = iypbnd(:,i)
              Tfyq(k) = fyqbnd(i) ; Tfyf(k) = fyfbnd(i)
              if ( scale_bond .eq. 1 ) then
                bndcls(k) =n_matrix(icls(iypbnd(1,i)),icls(iypbnd(2,i)))
              else
                bndcls(k) = maxval(n_matrix(:,:))
              endif
            endif
          enddo
          Tpara(2,np_bnd) = k
          if ( k .eq. nbnd ) exit
        enddo
        allocate(para_bnd(2,np_bnd))
        para_bnd(1:2,1:np_bnd) = Tpara(1:2,1:np_bnd)
        iypbnd(1:2,1:nbnd) = Tiyp(1:2,1:nbnd)
        fyqbnd(1:nbnd) = Tfyq(1:nbnd) ; fyfbnd(1:nbnd) = Tfyf(1:nbnd)
        deallocate(Tiyp,flg,flg2,Tfyq,Tfyf)

        ! For ang
        nang = iynang - iyndag ; np_ang = 0 ; k = 0
        allocate(Tiyp(3,nang),Tfyq(nang),Tfyf(nang),angcls(nang),      &
                 flg(ixnatm),flg2(nang))
        flg2(:) = .false.
        do
          flg(:) = .true. ; np_ang = np_ang + 1
          Tpara(1,np_ang) = k + 1
          do i = 1,nang
            if ( flg2(i) ) cycle
            if ( all(flg(iypang(1:3,i))) ) then
              k = k + 1 ; flg2(i) = .true. ; flg(iypang(1:3,i)) = .false.
              Tiyp(:,k) = iypang(:,i)
              Tfyq(k) = fyqang(i) ; Tfyf(k) = fyfang(i)
              if ( scale_angle .eq. 1 ) then
                angcls(k) =n_matrix(icls(iypang(1,i)),icls(iypang(3,i)))
              else
                angcls(k) = maxval(n_matrix(:,:))
              endif
            endif
          enddo
          Tpara(2,np_ang) = k
          if ( k .eq. nang ) exit
        enddo
        allocate(para_ang(2,np_ang))
        para_ang(1:2,1:np_ang) = Tpara(1:2,1:np_ang)
        iypang(1:3,1:nang) = Tiyp(1:3,1:nang)
        fyqang(1:nang) = Tfyq(1:nang) ; fyfang(1:nang) = Tfyf(1:nang)
        deallocate(Tiyp,flg,flg2,Tfyq,Tfyf)

        ! For 1-4 interactions
        n14int = count( iytnbf(1:iyntor).eq.1 ) ; np_14 = 0 ; k = 0
        allocate(Tiyp(2,n14int),Tvws(n14int),Tess(n14int),             &
                 iyp14(2,n14int),fytvws14(n14int),fytess14(n14int),    &
                 i14cls(n14int),flg(ixnatm),flg2(n14int),T14cls(n14int))
        do i = 1,iyntor
          if ( iytnbf(i) .eq. 1 ) then
            k = k + 1
            Tiyp(1,k) = iyptor(1,i) ; Tiyp(2,k) = iyptor(4,i)
            Tvws(k) = fytvws(i) ; Tess(k) = fytess(i)
!!!            T14cls(k) = n_matrix(icls(iyptor(2,i)),icls(iyptor(3,i)))
            T14cls(k) = n_matrix(icls(iyptor(1,i)),icls(iyptor(4,i)))
          endif
        enddo
        k = 0 ; flg2(:) = .false.
        do
          flg(:) = .true. ; np_14 = np_14 + 1
          Tpara(1,np_14) = k + 1
          do i = 1,n14int
            if ( flg2(i) ) cycle
            if ( all(flg(Tiyp(1:2,i))) ) then
              k = k + 1 ; flg2(i) = .true. ; flg(Tiyp(1:2,i)) = .false.
              iyp14(:,k) = Tiyp(:,i) ; fytvws14(k) = Tvws(i)
              fytess14(k) = Tess(i)
              i14cls(k) = T14cls(i)
            endif
          enddo
          Tpara(2,np_14) = k
          if ( k .eq. n14int ) exit
        enddo
        allocate(para_14(2,np_14))
        para_14(1:2,1:np_14) = Tpara(1:2,1:np_14)
!        forall(i=1:n14int)                                             &
!          i14cls(i) = n_matrix(icls(iyp14(1,i)),icls(iyp14(2,i)))
        deallocate(Tiyp,Tvws,Tess,flg,flg2,T14cls)

        ! For tor
        ntor = iyntor - iyndtr ; k = 0
        allocate(Tiyp(4,ntor),Tfyf(ntor),Tfyrot(ntor),Tfyphs(ntor),    &
                 Tiydiv(ntor),flg(ixnatm),flg2(ntor),Tiytnbf(ntor),    &
                 Tvws(ntor),Tess(ntor))
        !! Remove torsion pairs which don't use
        do i = 1,ntor
          if ( abs(fyftor(i)) .ge. epsilon(rtmp) ) then
            k = k + 1
            Tiyp(:,k) = iyptor(:,i)
            Tfyf(k) = fyftor(i) ; Tfyrot(k) = fytrot(i)
            Tfyphs(k) = fytphs(i) ; Tiydiv(k) = iytdiv(i)
            Tiytnbf(k) = iytnbf(i)
            Tvws(k) = fytvws(i) ; Tess(k) = fytess(i)
          endif
        enddo
        ntor = k
        iyptor(1:4,1:ntor) = Tiyp(1:4,1:ntor)
        fyftor(1:ntor) = Tfyf(1:ntor) ; fytrot(1:ntor) = Tfyrot(1:ntor)
        fytphs(1:ntor) = Tfyphs(1:ntor) ; iytdiv(1:ntor) =Tiydiv(1:ntor)
        iytnbf(1:ntor) = Tiytnbf(1:ntor)
        fytvws(1:ntor) = Tvws(1:ntor) ; fytess(1:ntor) = Tess(1:ntor)

        np_tor = 0 ; k = 0 ; flg2(:) = .false.
        do
          flg(:) = .true. ; np_tor = np_tor + 1
          Tpara(1,np_tor) = k + 1
          do i = 1,ntor
            if ( flg2(i) ) cycle
            if ( all(flg(iyptor(1:4,i))) ) then
              k = k + 1 ; flg2(i) = .true. ; flg(iyptor(1:4,i)) = .false.
              Tiyp(:,k) = iyptor(:,i)
              Tfyf(k) = fyftor(i) ; Tfyrot(k) = fytrot(i)
              Tfyphs(k) = fytphs(i) ; Tiydiv(k) = iytdiv(i)
              Tiytnbf(k) = iytnbf(i)
              Tvws(k) = fytvws(i) ; Tess(k) = fytess(i)
            endif
          enddo
          Tpara(2,np_tor) = k
          if ( k .eq. ntor ) exit
        enddo
        allocate(para_tor(2,np_tor))
        para_tor(1:2,1:np_tor) = Tpara(1:2,1:np_tor)
        iyptor(1:4,1:ntor) = Tiyp(1:4,1:ntor)
        fyftor(1:ntor) = Tfyf(1:ntor) ; fytrot(1:ntor) = Tfyrot(1:ntor)
        fytphs(1:ntor) = Tfyphs(1:ntor) ; iytdiv(1:ntor) =Tiydiv(1:ntor)
        iytnbf(1:ntor) = Tiytnbf(1:ntor)
        fytvws(1:ntor) = Tvws(1:ntor) ; fytess(1:ntor) = Tess(1:ntor)
        allocate(torcls(ntor))
        if ( scale_tor .eq. 0 ) then
          torcls(1:ntor) = maxval(n_matrix(:,:))
        else
          forall(i=1:ntor)                                             &
!          torcls(i) = n_matrix(icls(iyptor(1,i)),icls(iyptor(4,i)))
            torcls(i) = n_matrix(icls(iyptor(2,i)),icls(iyptor(3,i)))
          !! peptide bonds
          if ( scale_tor .ge. 2 ) then
            do i = 1,ntor
              j = iyptor(2,i) ; k = iyptor(3,i)
              if ( (cxatmn(j)(1:4).eq."C   " .and.                     &
                    cxatmn(k)(1:4).eq."N   ") .or.                     &
                   (cxatmn(j)(1:4).eq."N   " .and.                     &
                    cxatmn(k)(1:4).eq."C   ") ) then
                torcls(i) = maxval(n_matrix(:,:))
              endif
            enddo
            if ( scale_tor .eq. 3 ) then
              do i = 1,ntor
                j = iyptor(2,i) ; k = iyptor(3,i)
                if ( cxresn(k)(1:3).eq."PRO" .and.                     &
                     cxatmn(j)(1:4).eq."C   " .and.                    &
                     cxatmn(k)(1:4).eq."N   " ) then
                  torcls(i) =                                          &
                    n_matrix(icls(iyptor(2,i)),icls(iyptor(3,i)))
                endif
              enddo
            endif
          endif
        endif
        deallocate(Tiyp,Tfyf,Tfyrot,Tfyphs,Tiydiv,flg,flg2,Tiytnbf,    &
                   Tvws,Tess)

        ! For improper torsion
        nimp = iynimp - iyndip ; k = 0
        allocate(Tiyp(4,nimp),Tfyf(nimp),Tfyrot(nimp),Tfyphs(nimp),    &
                 Tiydiv(nimp),flg(ixnatm),flg2(nimp))
        !! Remove improper torsion pairs which don't use
        do i = 1,nimp
          if ( abs(fyfimp(i)) .ge. epsilon(rtmp) ) then
            k = k + 1
            Tiyp(:,k) = iypimp(:,i)
            Tfyf(k) = fyfimp(i) ; Tfyrot(k) = fyirot(i)
            Tfyphs(k) = fyiphs(i) ; Tiydiv(k) = iyidiv(i)
          endif
        enddo
        nimp = k
        iypimp(1:4,1:nimp) = Tiyp(1:4,1:nimp)
        fyfimp(1:nimp) = Tfyf(1:nimp) ; fyirot(1:nimp) = Tfyrot(1:nimp)
        fyiphs(1:nimp) = Tfyphs(1:nimp) ; iyidiv(1:nimp) =Tiydiv(1:nimp)

        np_imp = 0 ; k = 0 ; flg2(:) = .false.
        do
          flg(:) = .true. ; np_imp = np_imp + 1
          Tpara(1,np_imp) = k + 1
          do i = 1,nimp
            if ( flg2(i) ) cycle
            if ( all(flg(iypimp(1:4,i))) ) then
              k = k + 1 ; flg2(i) = .true. ; flg(iypimp(1:4,i)) = .false.
              Tiyp(:,k) = iypimp(:,i)
              Tfyf(k) = fyfimp(i) ; Tfyrot(k) = fyirot(i)
              Tfyphs(k) = fyiphs(i) ; Tiydiv(k) = iyidiv(i)
            endif
          enddo
          Tpara(2,np_imp) = k
          if ( k .eq. nimp ) exit
        enddo
        allocate(para_imp(2,np_imp))
        para_imp(1:2,1:np_imp) = Tpara(1:2,1:np_imp)
        iypimp(1:4,1:nimp) = Tiyp(1:4,1:nimp)
        fyfimp(1:nimp) = Tfyf(1:nimp) ; fyirot(1:nimp) = Tfyrot(1:nimp)
        fyiphs(1:nimp) = Tfyphs(1:nimp) ; iyidiv(1:nimp) =Tiydiv(1:nimp)
        allocate(impcls(nimp))
        if ( scale_imp .eq. 0 ) then
          impcls(1:nimp) = maxval(n_matrix(:,:))
        else
          forall(i=1:nimp) impcls(i) =                                 &
           n_matrix(minval(icls(iypimp(:,i))),maxval(icls(iypimp(:,i))))
        endif
        deallocate(Tiyp,Tfyf,Tfyrot,Tfyphs,Tiydiv,flg,flg2)

        ! For excess energy of ZD
        if ( iy15m .eq. 4 ) then
          np_excess = 0 ; k = 0
          allocate(Tiyp(3,nexcess),flg(ixnatm),flg2(nexcess))
          flg2(:) = .false.
          do
            flg(:) = .true. ; np_excess = np_excess + 1
            Tpara(1,np_excess) = k + 1
            do i = 1,nexcess
              if ( flg2(i) ) cycle
              if ( all(flg(iexcess(1:2,i))) ) then
                k = k + 1 ; flg2(i) = .true.
                flg(iexcess(1:2,i)) = .false.
                Tiyp(1:2,k) = iexcess(1:2,i)
                Tiyp(3,k) = Tcls(i)
              endif
            enddo
            Tpara(2,np_excess) = k
            if ( k .eq. nexcess ) exit
          enddo
          allocate(para_excess(2,np_excess))
          para_excess(1:2,1:np_excess) = Tpara(1:2,1:np_excess)
          iexcess(1:2,1:nexcess) = Tiyp(1:2,1:nexcess)
          allocate(excls(nexcess))
          excls(1:nexcess) = Tiyp(3,1:nexcess)
          deallocate(Tiyp,flg,flg2,Tcls)
        endif

      ! single
      else

        ! For bond
        nbnd = iynbnd - iyndbn ; np_bnd = 1 ; allocate(para_bnd(2,1))
        para_bnd(1:2,1) = (/1,nbnd/) ; allocate(bndcls(nbnd))
        if ( scale_bond .eq. 0 ) then
          bndcls(1:nbnd) = maxval(n_matrix(:,:))
        else
          forall(i=1:nbnd)                                             &
            bndcls(i) = n_matrix(icls(iypbnd(1,i)),icls(iypbnd(2,i)))
        endif

        ! For angle
        nang = iynang - iyndag ; np_ang = 1 ; allocate(para_ang(2,1))
        para_ang(1:2,1) = (/1,nang/) ; allocate(angcls(nang))
        if ( scale_angle .eq. 0 ) then
          angcls(1:nang) = maxval(n_matrix(:,:))
        else
          forall(i=1:nang)                                             &
            angcls(i) = n_matrix(icls(iypang(1,i)),icls(iypang(3,i)))
        endif

        ! For 1-4 interactions
        n14int = count( iytnbf(1:iyntor).eq.1 ) ; np_14 = 1 ; k = 0
        allocate(para_14(2,1))
        para_14(1:2,1) = (/1,n14int/)
        allocate(iyp14(2,n14int),fytvws14(n14int),fytess14(n14int),    &
                 i14cls(n14int))
        do i = 1,iyntor
          if ( iytnbf(i) .eq. 1 ) then
            k = k + 1
            iyp14(1,k) = iyptor(1,i) ; iyp14(2,k) = iyptor(4,i)
            fytvws14(k) = fytvws(i) ; fytess14(k) = fytess(i)
!!!            i14cls(k) = n_matrix(icls(iyptor(2,i)),icls(iyptor(3,i)))
            i14cls(k) = n_matrix(icls(iyptor(1,i)),icls(iyptor(4,i)))
          endif
        enddo
!        forall(i=1:n14int)                                             &
!          i14cls(i) = n_matrix(icls(iyp14(1,i)),icls(iyp14(2,i)))

        ! For torsion
        ntor = iyntor - iyndtr ; k = 0
        allocate(Tiyp(4,ntor),Tfyf(ntor),Tfyrot(ntor),Tfyphs(ntor),    &
                 Tiydiv(ntor),Tiytnbf(ntor),Tvws(ntor),Tess(ntor))
        !! Remove torsion pairs which don't use
        do i = 1,ntor
          if ( abs(fyftor(i)) .ge. epsilon(rtmp) ) then
            k = k + 1
            Tiyp(:,k) = iyptor(:,i)
            Tfyf(k) = fyftor(i) ; Tfyrot(k) = fytrot(i)
            Tfyphs(k) = fytphs(i) ; Tiydiv(k) = iytdiv(i)
            Tiytnbf(k) = iytnbf(i)
            Tvws(k) = fytvws(i) ; Tess(k) = fytess(i)
          endif
        enddo
        ntor = k
        iyptor(1:4,1:ntor) = Tiyp(1:4,1:ntor)
        fyftor(1:ntor) = Tfyf(1:ntor) ; fytrot(1:ntor) = Tfyrot(1:ntor)
        fytphs(1:ntor) = Tfyphs(1:ntor) ; iytdiv(1:ntor) =Tiydiv(1:ntor)
        iytnbf(1:ntor) = Tiytnbf(1:ntor)
        fytvws(1:ntor) = Tvws(1:ntor) ; fytess(1:ntor) = Tess(1:ntor)

        np_tor = 1 ; allocate(para_tor(2,1))
        para_tor(1:2,1) = (/1,ntor/) ; allocate(torcls(ntor))
        if ( scale_tor .eq. 0 ) then
          torcls(1:ntor) = maxval(n_matrix(:,:))
        else
          forall(i=1:ntor)                                             &
!          torcls(i) = n_matrix(icls(iyptor(1,i)),icls(iyptor(4,i)))
            torcls(i) = n_matrix(icls(iyptor(2,i)),icls(iyptor(3,i)))
          !! peptide bonds
          if ( scale_tor .eq. 1 ) then
            do i = 1,ntor
              j = iyptor(2,i) ; k = iyptor(3,i)
              if ( (cxatmn(j)(1:4).eq."C   " .and.                     &
                    cxatmn(k)(1:4).eq."N   ") .or.                     &
                   (cxatmn(j)(1:4).eq."N   " .and.                     &
                    cxatmn(k)(1:4).eq."C   ") ) then
                torcls(i) = maxval(n_matrix(:,:))
              endif
            enddo
          elseif ( scale_tor .eq. 2 ) then
            do i = 1,ntor
              j = iyptor(2,i) ; k = iyptor(3,i)
              if ( cxresn(k)(1:3) .ne. "PRO" ) cycle
              if ( (cxatmn(j)(1:4).eq."C   " .and.                     &
                    cxatmn(k)(1:4).eq."N   ") .or.                     &
                   (cxatmn(j)(1:4).eq."N   " .and.                     &
                    cxatmn(k)(1:4).eq."C   ") ) then
                torcls(i) = maxval(n_matrix(:,:))
              endif
            enddo
          endif
        endif

        deallocate(Tiyp,Tfyf,Tfyrot,Tfyphs,Tiydiv,Tiytnbf,Tvws,Tess)
        
        ! For improper torsion
        nimp = iynimp - iyndip ; k = 0
        allocate(Tiyp(4,nimp),Tfyf(nimp),Tfyrot(nimp),Tfyphs(nimp),    &
                 Tiydiv(nimp))
        !! Remove improper torsion pairs which don't use
        do i = 1,nimp
          if ( abs(fyfimp(i)) .ge. epsilon(rtmp) ) then
            k = k + 1
            Tiyp(:,k) = iypimp(:,i)
            Tfyf(k) = fyfimp(i) ; Tfyrot(k) = fyirot(i)
            Tfyphs(k) = fyiphs(i) ; Tiydiv(k) = iyidiv(i)
          endif
        enddo
        nimp = k
        iypimp(1:4,1:nimp) = Tiyp(1:4,1:nimp)
        fyfimp(1:nimp) = Tfyf(1:nimp) ; fyirot(1:nimp) = Tfyrot(1:nimp)
        fyiphs(1:nimp) = Tfyphs(1:nimp) ; iyidiv(1:nimp) =Tiydiv(1:nimp)

        np_imp = 1 ; allocate(para_imp(2,1))
        para_imp(1:2,1) = (/1,nimp/) ; allocate(impcls(nimp))
        if ( scale_imp .eq. 0 ) then
          impcls(1:nimp) = maxval(n_matrix(:,:))
        else
          forall(i=1:nimp)                                             &
            impcls(i) = n_matrix(icls(iypimp(1,i)),icls(iypimp(4,i)))
        endif
        deallocate(Tiyp,Tfyf,Tfyrot,Tfyphs,Tiydiv)

        ! For excess energy of ZD
        if ( iy15m .eq. 4 ) then
          np_excess = 1 ; allocate(para_excess(2,1))
          para_excess(1:2,1) = (/1,nexcess/) ; allocate(excls(nexcess))
          excls(1:nexcess) = Tcls(1:nexcess)
          deallocate(Tcls)
        endif

      endif

      do i = 1,n14int
        fytess14(i) = fytess14(i) * ccoef *                            &
                      fxchrg(iyp14(1,i))*fxchrg(iyp14(2,i))
      enddo

!************************************************************

      return
      end subroutine incel


!======================================================================


      subroutine mk_cm(n)

!**************************************************
!
!     Nearest far cells for each level.
!
!**************************************************

      use COMCMM

      implicit none

      integer(4),intent(in):: n

      integer(4),allocatable:: ncm(:,:,:),icm(:,:,:,:,:)
      integer(4):: ix,iy,iz,jx,jy,jz,jx1,jx2,jy1,jy2,jz1,jz2,itmp,kx,  &
                   ky,kz,ix1,ix2,iy1,iy2,iz1,iz2,m,m1
      logical(4):: flgx,flgy,flgz

!*******************************************

      m = nv(n)
      allocate(ncm(m,m,m),icm(3,jcmx,m,m,m))

      if ( n .ne. nlev ) then

        m1 = nv(n+1)
        do iz = 1,m
          jz = lvup(iz) ; iz1 = max(1,iz-1) ; iz2 = min(m,iz+1)
          jz1 = max(1,jz-1) * 2 - 1 ; jz2 = min(m1,jz+1) * 2
        do iy = 1,m
          jy = lvup(iy) ; iy1 = max(1,iy-1) ; iy2 = min(m,iy+1)
          jy1 = max(1,jy-1) * 2 - 1 ; jy2 = min(m1,jy+1) * 2
        do ix = 1,m
          jx = lvup(ix) ; ix1 = max(1,ix-1) ; ix2 = min(m,ix+1)
          jx1 = max(1,jx-1) * 2 - 1 ; jx2 = min(m1,jx+1) * 2

          itmp = 0
          do kz = jz1,jz2
            if ( kz.ge.iz1 .and. kz.le.iz2 ) then
              flgz = .true.
            else
              flgz = .false.
            endif
          do ky = jy1,jy2
            if ( ky.ge.iy1 .and. ky.le.iy2 ) then
              flgy = .true.
            else
              flgy = .false.
            endif
          do kx = jx1,jx2
            if ( kx.ge.ix1 .and. kx.le.ix2 ) then
              flgx = .true.
            else
              flgx = .false.
            endif

            if ( flgx .and. flgy .and. flgz ) cycle
            itmp = itmp + 1
            icm(1:3,itmp,ix,iy,iz) = (/kx,ky,kz/)
          enddo
          enddo
          enddo

          ncm(ix,iy,iz) = itmp
        enddo
        enddo
        enddo

      else

        do iz = 1,m
          iz1 = max(1,iz-1) ; iz2 = min(m,iz+1)
        do iy = 1,m
          iy1 = max(1,iy-1) ; iy2 = min(m,iy+1)
        do ix = 1,m
          ix1 = max(1,ix-1) ; ix2 = min(m,ix+1)

          itmp = 0
          do jz = 1,m
            if ( jz.ge.iz1 .and. jz.le.iz2 ) then
              flgz = .true.
            else
              flgz = .false.
            endif
          do jy = 1,m
            if ( jy.ge.iy1 .and. jy.le.iy2 ) then
              flgy = .true.
            else
              flgy = .false.
            endif
          do jx = 1,m
            if ( jx.ge.ix1 .and. jx.le.ix2 ) then
              flgx = .true.
            else
              flgx = .false.
            endif

            if ( flgx .and. flgy .and. flgz ) cycle
            itmp = itmp + 1
            icm(1:3,itmp,ix,iy,iz) = (/jx,jy,jz/)
          enddo
          enddo
          enddo

          ncm(ix,iy,iz) = itmp
        enddo
        enddo
        enddo

      endif

      select case (n)
        case(1)
          allocate(ncm1(m,m,m),icm1(3,jcmx,m,m,m))
          ncm1 = ncm ; icm1 = icm
        case(2)
          allocate(ncm2(m,m,m),icm2(3,jcmx,m,m,m))
          ncm2 = ncm ; icm2 = icm
        case(3)
          allocate(ncm3(m,m,m),icm3(3,jcmx,m,m,m))
          ncm3 = ncm ; icm3 = icm
        case(4)
          allocate(ncm4(m,m,m),icm4(3,jcmx,m,m,m))
          ncm4 = ncm ; icm4 = icm
        case(5)
          allocate(ncm5(m,m,m),icm5(3,jcmx,m,m,m))
          ncm5 = ncm ; icm5 = icm
      end select

!***************************************

      return
      end subroutine mk_cm


!====================================================================


      subroutine mk_c(n)

!************************************************
!
!     Coordinates of cells at different levels.
!
!************************************************

      use COMCMM

      implicit none

      integer(4),intent(in):: n
      real(8),allocatable:: c(:,:)

      integer(4):: i
      real(8):: r

!****************************************

      allocate(c(3,nv(n)))
      r = -0.5d0*(nv(n)+1)
      forall(i=1:nv(n)) c(1:3,i) = siz(n) * (dble(i)+r)
      select case (n)
        case(1)
          allocate(c1(3,nv(1))) ; c1 = c
        case(2)
          allocate(c2(3,nv(2))) ; c2 = c
        case(3)
          allocate(c3(3,nv(3))) ; c3 = c
        case(4)
          allocate(c4(3,nv(4))) ; c4 = c
        case(5)
          allocate(c5(3,nv(5))) ; c5 = c
      end select

!****************************************

      return
      end subroutine mk_c
