
      subroutine inrefbinL

!***********************************************
!
!     INput REF. bin file List for RMSD calc.
!
!***********************************************

      use COMIFN ; use COMVAL

      implicit none

      ! input pdb FILe NAMe
        character(130):: FILnam

      integer(4):: i,j,k,k1,k2,iST,iEN,trseed,ichance,icn
      real(4):: rr,ep
      real(8),allocatable:: tary(:,:,:),tary2(:,:,:),tary3(:)
      integer(4),allocatable:: tiary(:)
      logical(4):: ex
      logical(4):: flg = .false.
      logical(4):: flg2 = .false.

      ! presto output
      integer(4):: istp,iyn15v,iyn15h
      real(4):: sitime,sec,et,ek,temp,lambda,rmsf,rmsd2,ttepotcc(3)
      real(4):: tcod(3,nATMr)

!*******************************

      write(6,'(4x,a)')"+ binary files in "//                          &
        trim(input_reference_binary_list)//" input START"
      write(6,'(a4,$)')"    "

      if ( input_binary_list.eq.input_reference_binary_list ) flg= .true.
      if ( flg .and. (len(trim(input_prob_file)).ne.0 .or.             &
           len(trim(input_weight_file)).ne.0 ) ) flg2 = .true.
      if ( flg2 ) trseed = rseed
      ! Count conformations
      k = 0
      open(unit=1,file=trim(input_reference_binary_list),status="old")
      if ( len(trim(input_weight_file)) .ne. 0 )                       &
        open(unit=3,file=trim(input_weight_file),status="old")
      do
        read(1,*,end=800)FILnam
        if ( len(trim(FILnam)) .eq. 0 ) cycle
        read(1,*,end=800)iST,iEN
        if ( iST .le. 0 ) iST = 1
        if ( iEN .lt. iST ) cycle
        k = k + iEN - iST + 1
      enddo
800   allocate(tary(3,nRMSDf,k),tary2(3,nRMSD,k),tary3(k),tiary(k))

      !! Input binary files
      rewind(1) ; i = 0 ; ichance = 0 ; icn = 0
      do
        read(1,'(a)',end=801)FILnam
        if ( len(trim(FILnam)) .eq. 0 ) cycle
        inquire(file=trim(FILnam),exist=ex)
        if ( .not. ex ) cycle

        read(1,*)iST,iEN
        if ( iST .le. 0 ) iST = 1
        if ( iEN .lt. iST ) cycle

        open(unit=2,file=trim(FILnam),form="unformatted",status="old")
        ! Skip reading
        do j = 1,iST-1
          read(2,end=802)
          read(2,end=802)
        enddo
        ! Reading
        do j = iST,iEN
          read(2,end=802)istp,sitime,sec,et,ek,temp,ep,lambda,iyn15v,  &
                         iyn15h,rmsd2,ttepotcc(1:3)
          read(2,end=802)tcod(1:3,1:nATMr)
          icn = icn + 1

          ! Weighting factor
          if ( flg2 ) then
            if ( len(trim(input_prob_file)) .ne. 0 ) then
              wfac = 0.d0
              if ( GE_method .eq. "MULT" ) then
                k = int((ep-minEc(1))/EbinSZ(1)+1)
                if ( k.ge.1 .and. k.le.nBIN(1) ) wfac = wei(k)
              elseif ( GE_method .eq. "ALSD" ) then
                k = int((ttepotcc(1)-minEc(1))/EbinSZ(1)+1)
                if ( k.ge.1 .and. k.le.nBIN(1) ) then
                  k2 = int((ttepotcc(2)-minEc(2))/EbinSZ(2)+1)
                  if ( k2.ge.1 .and. k2.le.nBIN(2) ) wfac = wei2(k,k2)
                endif
              endif
            elseif ( len(trim(input_weight_file)) .ne. 0 ) then
              write(3,*)k,wfac
            endif
            ! Skip calc. for conformation with very low probability
            ! Select conformation by wfac
            if ( weight_method .eq. "SELECT" ) then
              trseed = 843314861*trseed+453816693
              rr = real(trseed)/2**31 * 0.5d0 + 0.5d0
              if ( wfac.lt.rr ) cycle
            elseif ( weight_method .eq. "WEIGHT" ) then
              if ( wfac .lt. weight_threshold ) cycle
            endif
          else
            wfac = 1.d0
          endif
          ichance = ichance + 1
          if ( flg .and. mod(ichance,interval) .ne. 0 ) cycle
          i = i + 1 ; tiary(i) = icn ; tary3(i) = wfac
          if ( mod(i,100) .eq. 0 ) write(6,'(a1,$)')"*"

          rr = minval(tcod(1,1:nATMr)) ; minX = min(minX,rr)
          rr = maxval(tcod(1,1:nATMr)) ; maxX = max(maxX,rr)
          rr = minval(tcod(2,1:nATMr)) ; minY = min(minY,rr)
          rr = maxval(tcod(2,1:nATMr)) ; maxY = max(maxY,rr)
          rr = minval(tcod(3,1:nATMr)) ; minZ = min(minZ,rr)
          rr = maxval(tcod(3,1:nATMr)) ; maxZ = max(maxZ,rr)
          k1 = 1 ; k2 = 1
          do k = 1,nATMr
            if ( iATMrmsdRF(k1) .eq. k ) then
              tary(1:3,k1,i) = tcod(1:3,k) ; k1 = k1 + 1
            endif
            if ( iATMrmsdR(k2) .eq. k ) then
              tary2(1:3,k2,i) = tcod(1:3,k) ; k2 = k2 + 1
            endif
            if ( k1.gt.nRMSDf .and. k2.gt.nRMSD ) exit
          enddo
        enddo
802     close(2)
      enddo
801   close(1)
      if ( len(trim(input_weight_file)) .ne. 0 ) close(3)

      nREF = i
      write(6,'(4x,a2,i0,a)')"+ ",nREF," ref. PDB file(s) are input"
      write(6,*)
      allocate(RcodRMSDf(3,nRMSDf,nREF),RcodRMSD(3,nRMSD,nREF),        &
               aveRcodF(3,nRMSDf),aveRcod(3,nRMSD),conf_n(nREF),       &
               conf_w(nREF))
      conf_n(1:nREF) = tiary(1:nREF)
      conf_w(1:nREF) = tary3(1:nREF)
      RcodRMSDf(1:3,1:nRMSDf,1:nREF) = tary(1:3,1:nRMSDf,1:nREF)
      RcodRMSD(1:3,1:nRMSD,1:nREF) = tary2(1:3,1:nRMSD,1:nREF)
      aveRcodF = 0.d0 ; aveRcod = 0.d0
      do i = 1,nREF
        forall (j=1:3,k=1:nRMSDf)                                      &
          aveRcodF(j,k) = aveRcodF(j,k) + RcodRMSDf(j,k,i) * conf_w(i)
      enddo
      aveRcodF(:,:) = aveRcodF(:,:) / sum(conf_w)
      do i = 1,nREF
        forall (j=1:3,k=1:nRMSD)                                       &
          aveRcod(j,k) = aveRcod(j,k) + RcodRMSD(j,k,i) * conf_w(i)
      enddo
      aveRcod(:,:) = aveRcod(:,:) / sum(conf_w)

!*******************************

      return
      end subroutine inrefbinL
