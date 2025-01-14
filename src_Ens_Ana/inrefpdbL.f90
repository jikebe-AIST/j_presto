
      subroutine inrefpdbL

!***********************************************
!
!     INput REF. PDB file List for RMSD calc.
!
!***********************************************

      use COMIFN ; use COMVAL

      implicit none

      ! input pdb FILe NAMe
        character(130):: FILnam

      integer(4):: i,j,k,k1,k2,trseed,ichance,icn
      real(4):: rr,ep,ttepotcc(3),tcod(3)
      real(8),allocatable:: tary(:,:,:),tary2(:,:,:),tary3(:)
      integer(4),allocatable:: tiary(:)
      logical(4):: ex
      logical(4):: flg = .false.
      logical(4):: flg2 = .false.
      character(200):: tmp

!*******************************

      write(6,'(4x,a)')"+ Reference PDB files are input "

      if ( input_PDB_list.eq.input_reference_PDB_list ) flg = .true.
      if ( flg .and. (len(trim(input_prob_file)).ne.0 .or.             &
           len(trim(input_weight_file)) .ne. 0 ) ) flg2 = .true.
      if ( flg2 ) then
        trseed = rseed
        open(unit=3,file=trim(input_energy_file),status="old")
      endif
      if ( len(trim(input_weight_file)) .ne. 0 )                       &
        open(unit=4,file=trim(input_weight_file),status="old")

      ! Count number of pdb files
      call system("wc -l "//trim(input_reference_PDB_list)//           &
                  " > "//trim(d_proj2)//".aho.temp")
      open(unit=1,file=trim(d_proj2)//".aho.temp",status="old")
      read(1,*)i
      close(unit=1,status="delete")
      allocate(tary(3,nRMSDf,i),tary2(3,nRMSD,i),tary3(i),tiary(i))

      !! Input PDB files
      open(unit=1,file=trim(input_reference_PDB_list),status="old")
      i = 0 ; ichance = 0 ; icn = 0
      do
        read(1,'(a)',end=800)FILnam
        if ( FILnam .eq. " " ) cycle
        inquire(file=FILnam,exist=ex)
        if ( .not. ex ) cycle
        icn = icn + 1

        ! Weighting factor
        if ( flg2 ) then
          if ( len(trim(input_prob_file)) .ne. 0 ) then
            wfac = 0.d0
            if ( GE_method .eq. "MULT" ) then
              read(3,*)ep
              k = int((ep-minEc(1))/EbinSZ(1)+1)
              if ( k.ge.1 .and. k.le.nBIN(1) ) wfac = wei(k)
            elseif ( GE_method .eq. "ALSD" ) then
              read(3,*)ttepotcc(1:2)
              k = int((ttepotcc(1)-minEc(1))/EbinSZ(1)+1)
              if ( k.ge.1 .and. k.le.nBIN(1) ) then
                k2 = int((ttepotcc(2)-minEc(2))/EbinSZ(2)+1)
                if ( k2.ge.1 .and. k2.le.nBIN(2) ) wfac = wei2(k,k2)
              endif
            endif
          elseif ( len(trim(input_weight_file)) .ne. 0 ) then
            read(4,*)k,wfac
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

        !! Read each PDB files
        write(6,'(6x,a)')trim(FILnam)
        open(unit=2,file=trim(FILnam),status="old") ; k1 = 1 ; k2 = 1
        j = 0
        do
          read(2,'(a)',end=801)tmp
          if ( tmp(1:6).ne."ATOM  " .and. tmp(1:6).ne."HETATM" ) cycle
          j = j + 1
          read(tmp,'(30x,3f8.3)')tcod(1:3)
          minX = min(minX,tcod(1)) ; maxX = max(maxX,tcod(1))
          minY = min(minY,tcod(2)) ; maxY = max(maxY,tcod(2))
          minZ = min(minZ,tcod(3)) ; maxZ = max(maxZ,tcod(3))
          if ( k1 .le. nRMSDf ) then
            if ( iATMrmsdRF(k1) .eq. j ) then
              tary(1:3,k1,i) = tcod(1:3) ; k1 = k1 + 1
            endif
          endif
          if ( k2 .le. nRMSD ) then
            if ( iATMrmsdR(k2) .eq. j ) then
              tary2(1:3,k2,i) = tcod(1:3) ; k2 = k2 + 1
            endif
          endif
        enddo
801     close(unit=2)
      enddo
800   close(1)
      if ( flg2 ) close(3)
      if ( len(trim(input_weight_file)) .ne. 0 ) close(4)

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
        forall (j=1:3,k=1:nRMSD)                                      &
          aveRcod(j,k) = aveRcod(j,k) + RcodRMSD(j,k,i) * conf_w(i)
      enddo
      aveRcod(:,:) = aveRcod(:,:) / sum(conf_w)

!*******************************

      return
      end subroutine inrefpdbL
