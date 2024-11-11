
      subroutine celpar_PB(ilflag,iprint,ier) ! #SLC2 #PB
      subroutine celpar_CB(ilflag,iprint,ier) ! #SLC2 #CB

!****************************************************************
!
!     Put atoms in the cell-field.
!     The deepest cell-level = 5.
!     The unit-base division is done using geometrical-center
!       of the unit.
!
!****************************************************************

      use COMPAR ; use COMCMM
      use COMBAS ; use COMERG
      use COMCMMC,only: cord
      !$ use omp_lib

      implicit none

      integer(4),intent(in):: ilflag,iprint
      integer(4),intent(inout):: ier

      integer(4):: i,ii,j,jj,nfg,mx,my,mz,num,k1,k2,isum,ix,iy,iz,     &
                   mnx,mny,mnz,mxx,mxy,mxz,ia(1)
      integer(4):: kmax(nfrag),nzero(nfrag),nnonz(nfrag),avn(nfrag)
      real(8):: cord0(3),xx,yy,zz,gx,gy,gz,rtmp,rr
      character(130):: tmp

      logical(4),save:: CMMfstFLG = .true.

!******************************************************
!     Shift the cell coordinates by the field-center.
                                                                   ! #SLC2 #CB
      ! Calc. the system center                                    ! #SLC2 #CB
      if ( .not. extCMM_flag ) then                                ! #SLC2 #CB
        rtmp = 1.d0 / dble(ixnatm)                                 ! #SLC2 #CB
        forall ( i=1:3 ) cord0(i) = sum(cord(i,1:ixnatm)) * rtmp   ! #SLC2 #CB
                                                                   ! #SLC2 #CB
        ! Shift coordinates                                        ! #SLC2 #CB
        forall ( i=1:nv(1) ) c1(:,i) = cg1(:,i) + cord0(:)         ! #SLC2 #CB
        if ( nlev .ge. 2 ) then                                    ! #SLC2 #CB
          forall ( i=1:nv(2) ) c2(:,i) = cg2(:,i) + cord0(:)       ! #SLC2 #CB
          if ( nlev .ge. 3 ) then                                  ! #SLC2 #CB
            forall ( i=1:nv(3) ) c3(:,i) = cg3(:,i) + cord0(:)     ! #SLC2 #CB
            if ( nlev .ge. 4 ) then                                ! #SLC2 #CB
              forall ( i=1:nv(4) ) c4(:,i) = cg4(:,i) + cord0(:)   ! #SLC2 #CB
              if ( nlev .ge. 5 )                                 & ! #SLC2 #CB
                forall ( i=1:nv(5) ) c5(:,i) = cg5(:,i)+cord0(:)   ! #SLC2 #CB
            endif                                                  ! #SLC2 #CB
          endif                                                    ! #SLC2 #CB
        endif                                                      ! #SLC2 #CB
                                                                   ! #SLC2 #CB
      elseif ( CMMfstFLG ) then                                    ! #SLC2 #CB
        ! Shift coordinates                                        ! #SLC2 #CB
        cord0(1:3) = (/celx,cely,celz/)                            ! #SLC2 #CB
        forall ( i=1:nv(1) ) c1(:,i) = c1(:,i) + cord0(:)          ! #SLC2 #CB
        if ( nlev .ge. 2 ) then                                    ! #SLC2 #CB
          forall ( i=1:nv(2) ) c2(:,i) = c2(:,i) + cord0(:)        ! #SLC2 #CB
          if ( nlev .ge. 3 ) then                                  ! #SLC2 #CB
            forall ( i=1:nv(3) ) c3(:,i) = c3(:,i) + cord0(:)      ! #SLC2 #CB
            if ( nlev .ge. 4 ) then                                ! #SLC2 #CB
              forall ( i=1:nv(4) ) c4(:,i) = c4(:,i) + cord0(:)    ! #SLC2 #CB
              if ( nlev .ge. 5 )                                 & ! #SLC2 #CB
                forall ( i=1:nv(5) ) c5(:,i) = c5(:,i) + cord0(:)  ! #SLC2 #CB
            endif                                                  ! #SLC2 #CB
          endif                                                    ! #SLC2 #CB
        endif                                                      ! #SLC2 #CB
        CMMfstFLG = .false.                                        ! #SLC2 #CB
      endif                                                        ! #SLC2 #CB
                                                                   ! #SLC2 #CB
      if ( CMMfstFLG ) then                                        ! #SLC2 #PB
        ! Shift coordinates                                        ! #SLC2 #PB
        forall ( i=1:nv(1) ) c1(:,i) = c1(:,i) + fxcbou(:)         ! #SLC2 #PB
        if ( nlev .ge. 2 ) then                                    ! #SLC2 #PB
          forall ( i=1:nv(2) ) c2(:,i) = c2(:,i) + fxcbou(:)       ! #SLC2 #PB
          if ( nlev .ge. 3 ) then                                  ! #SLC2 #PB
            forall ( i=1:nv(3) ) c3(:,i) = c3(:,i) + fxcbou(:)     ! #SLC2 #PB
            if ( nlev .ge. 4 ) then                                ! #SLC2 #PB
              forall ( i=1:nv(4) ) c4(:,i) = c4(:,i) + fxcbou(:)   ! #SLC2 #PB
              if ( nlev .ge. 5 )                                 & ! #SLC2 #PB
                forall ( i=1:nv(5) ) c5(:,i) = c5(:,i) + fxcbou(:) ! #SLC2 #PB
            endif                                                  ! #SLC2 #PB
          endif                                                    ! #SLC2 #PB
        endif                                                      ! #SLC2 #PB
        CMMfstFLG = .false.                                        ! #SLC2 #PB
      endif                                                        ! #SLC2 #PB

!******************************************************
!     Initialization.

!      call system_clock(tim1)
      zz = siz(1)*0.5d0
      xx = c1(1,1) - zz
      yy = c1(2,1) - zz
      zz = c1(3,1) - zz
      !$OMP parallel default (none)                                  & !
      !$OMP private(i,j)                                             & !
      !$OMP shared(npcl1,nv)
      !$OMP do schedule (static) collapse(2)
      do j = 1,nv(1)
      do i = 1,nv(1)
        npcl1(:,:,i,j) = 0
      enddo
      enddo
      !$OMP end do
      !$OMP end parallel
      rr = 1.d0 / siz(1)

!******************
!     Putting atoms in cells.

      ! CMM or atom-based cut-off
      ! Atom-base
      do i = 1,ncmmatm
      do j = icmmatm(1,i),icmmatm(2,i)
        gx = cord(1,j)-fxcell(1)*floor((cord(1,j)-celwal(1))*invcel(1)) ! #SLC2 #PB
        gy = cord(2,j)-fxcell(2)*floor((cord(2,j)-celwal(3))*invcel(2)) ! #SLC2 #PB
        gz = cord(3,j)-fxcell(3)*floor((cord(3,j)-celwal(5))*invcel(3)) ! #SLC2 #PB
        gx = gx-xx ; gy = gy-yy ; gz = gz-zz                            ! #SLC2 #PB
        gx = cord(1,j)-xx ; gy = cord(2,j)-yy ; gz = cord(3,j)-zz       ! #SLC2 #CB
        mx = int(gx*rr)+1 ; my = int(gy*rr)+1 ; mz = int(gz*rr)+1
        nfg = icls(j)
        num = npcl1(nfg,mx,my,mz) + 1
        npcl1(nfg,mx,my,mz) = num
        ipcl(num,nfg,mx,my,mz) = j
      enddo
      enddo

      ! Residue-base
      do i = 1,ncmmres
      do j = icmmres(1,i),icmmres(2,i)

        k1 = ixrstr(j) ; k2 = ixrend(j)

        rtmp = 1.d0 / dble(k2-k1+1)
        gx = sum(cord(1,k1:k2)) * rtmp
        gy = sum(cord(2,k1:k2)) * rtmp
        gz = sum(cord(3,k1:k2)) * rtmp
        gx = gx-fxcell(1)*floor((gx-celwal(1))*invcel(1)) ! #SLC2 #PB
        gy = gy-fxcell(2)*floor((gy-celwal(3))*invcel(2)) ! #SLC2 #PB
        gz = gz-fxcell(3)*floor((gz-celwal(5))*invcel(3)) ! #SLC2 #PB

        gx = gx - xx ; gy = gy - yy ; gz = gz - zz
        mx = int(gx*rr)+1 ; my = int(gy*rr)+1 ; mz = int(gz*rr)+1
        do jj = k1,k2
          nfg = icls(jj)
          num = npcl1(nfg,mx,my,mz) + 1
          npcl1(nfg,mx,my,mz) = num
          ipcl(num,nfg,mx,my,mz) = jj
        enddo

      enddo
      enddo
!      call system_clock(tim2)
!      tmptime = tmptime + tim2 - tim1

!******************
!  Check.

      j = nv(1) ; avn(:) = 0 ; kmax(:) = 0
      !$OMP parallel default (none)                                  & !
      !$OMP private(i,ix,iy,iz)                                      & !
      !$OMP shared(j,nfrag,npcl1)                                    & !
      !$OMP reduction (+ : avn)                                      & !
      !$OMP reduction (max : kmax)
      !$OMP do schedule (static) collapse (2)
      do iz = 1,j
      do iy = 1,j
      do ix = 1,j
      do i = 1,nfrag
        avn(i) = avn(i) + npcl1(i,ix,iy,iz)
        kmax(i) = max(kmax(i),npcl1(i,ix,iy,iz))
      enddo
      enddo
      enddo
      enddo
      !$OMP end do
      !$OMP end parallel
      isum = sum(avn(1:nfrag))

      ipmax_max = max(ipmax_max,maxval(kmax(:)))
      if ( ipmax_max .gt. ipmax ) then
        print*," ???????????????????????????????????????????"
        print*," ? Too large number of atoms was detected  ?"
        print*," ? in a cell.  Make ipmax bigger. (celpar) ?"
        print*," ???????????????????????????????????????????"
        print*,"   N (detected) =",ipmax_max
        print*,"   ipmax        =",ipmax
        print*,""
        ier = 28 ; return
      endif
        
      if ( isum.ne.ixnatm ) then
        print*," ?????????????????????????????????????????"
        print*," ?    atom out of the field (celpar)     ?"
        print*," ?????????????????????????????????????????"
        print*,"    N (sum) =",isum
        print*,"    N (tot) =",ixnatm
        print*,""
        call outpdb(89,6,"ERROR.pdb",maxatm,ixnatm,cxatmn,cxresn,    &
                    absres,fxmass,fxchrg,10,ixtitl,cxtitl,ier)
        ier = 28 ; return
      endif
!****
      if ( ilflag .eq. 2 ) then
        nzero(1:nfrag) = 0
        !$OMP parallel default (none)                                & !
        !$OMP private(iz,iy,ix,i)                                    & !
        !$OMP shared(j,nfrag,npcl1)                                  & !
        !$OMP reduction (+ : nzero)
        !$OMP do schedule (static) collapse(2)
        do iz = 1,j
        do iy = 1,j
        do ix = 1,j
          do i = 1,nfrag
            if ( npcl1(i,ix,iy,iz) .eq. 0 ) nzero(i) = nzero(i) + 1
          enddo
        enddo
        enddo
        enddo
        !$OMP end do
        !$OMP end parallel
!        forall ( i=1:nfrag ) avn(i) = sum(npcl1(i,1:j,1:j,1:j))
        i = nv(1)*nv(1)*nv(1)
        do ii = 1,nfrag
          nnonz(ii) = i - nzero(ii)
          write(tmp,*)ii
          write(iprint,*)""
          write(iprint,*)"< CLUSTER-"//trim(adjustl(tmp))//" >"
          write(iprint,*)"   No. OF     EMPTY CELLS : ",nzero(ii)
          write(iprint,*)"   No. OF NON-EMPTY CELLS : ",nnonz(ii)
          write(iprint,*)"   MAXIMUM ATOM No. IN A CELL : ",kmax(ii)
          if ( nnonz(ii) .gt. 0 )                                      &
           write(iprint,*)"   AVERAGE ATOM No. OVER NON-EMPTY CELLS "  &
               //": ",dble(avn(ii))/dble(nnonz(ii))
        enddo
      endif

!************************************************************

      return
      end subroutine celpar_PB ! #SLC2 #PB
      end subroutine celpar_CB ! #SLC2 #CB
