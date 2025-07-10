
      subroutine store_element

!*******************************************************
!
!     store element data to variable of input for ToPoLogy
!
!*******************************************************

      use COMINP ; use COMIFN ; use COMVAL

      implicit none

      real(8):: rtmp
      logical(1):: ex

!***************************

      ! input PROJect NaMe (element 1)
      PROJNM = ELEMNT(1)
      write(6,'(2x,a)')"* Input Project name : "
      write(6,'(8x,a)')trim(PROJNM)

      ! Method (element 2) 
      METHOD = ELEMNT(2) 
      if ( METHOD .eq. "CANO" ) then
        write(6,'(2x,a)')"* Method : "   
        write(6,'(8x,a)')"Canonical MD"   
      elseif ( METHOD .eq. "MULT" ) then
        write(6,'(2x,a)')"* Method : "   
        write(6,'(8x,a)')"MultiCanonical MD"   
      elseif ( METHOD(1:4) .eq. "CLMD" ) then
        write(6,'(2x,a)')"* Method : "   
        write(6,'(8x,a)')"Constant lambda molecular dynamics"
      elseif ( METHOD(1:4) .eq. "ALSD" ) then
        write(6,'(2x,a)')"* Method : "   
        write(6,'(8x,a)')"Adaptive lambda square dynamics"   
      else
        write(6,*)"METHOD you input : ",METHOD
        call error(10201)
      endif
      if ( METHOD(5:5) .eq. "2" ) write(6,'(8x,a)')"(lambda square)"

      ! DATa LiST name (element 3)
      DATLST = ELEMNT(3)
      if ( DATLST .eq. " " ) call error(10202)
      inquire(file=trim(DATLST),exist=ex)
      if ( .not. ex ) call error(10203)
      write(6,'(2x,a)')"* Input energy data list name : "
      write(6,'(8x,a)')trim(DATLST)
      
      ! Min & Max reaction cord. range in the current MD (element 4 & 5)
      if ( ELEMNT(4) .eq. " " ) call error(10204)
      if ( ELEMNT(5) .eq. " " ) call error(10205)
      read(ELEMNT(4),*)minE
      read(ELEMNT(5),*)maxE
      if ( METHOD.eq."CANO" .or. METHOD.eq."MULT" ) then
        write(6,'(2x,a)')"* Energy range for sampling : "
        write(6,'(8x,f12.3,a,f12.3,a)')minE," - ",maxE," (kcal/mol)"
      elseif ( METHOD.eq."CLMD" .or. METHOD.eq."ALSD" ) then
        write(6,'(2x,a)')"* Lambda range for sampling : "
        write(6,'(8x,f12.3,a,f12.3,a)')minE," - ",maxE
      elseif ( METHOD.eq."CLMD2" .or. METHOD.eq."ALSD2" ) then
        write(6,'(2x,a)')"* Lambda square range for sampling : "
        write(6,'(8x,f12.3,a,f12.3,a)')minE," - ",maxE
      endif
      if ( minE .ge. maxE ) call error(10206)

      ! current Temperature (element 6)
      if ( ELEMNT(6) .eq. " " ) call error(10207)
      read(ELEMNT(6),*)T
      write(6,'(2x,a)')"* Simulation Temperature : "
      write(6,'(8x,f10.3,a)')T," (K)"
      if ( T .lt. 0 ) call error(10208)

      ! array dimension for fitting (element 7-8)
      read(ELEMNT(7),*)NfitDIM
      read(ELEMNT(8),*)Nwindow
      write(6,'(2x,a)')"* Array dimension for fitting = "
      write(6,'(8x,i0)')NfitDIM
      write(6,'(2x,a)')"* Number of windows = "
      write(6,'(8x,i0)')Nwindow
      if ( (NBIN-1)/Nwindow .lt. 1 ) call error(10209)

      ! Min & Max reaction cord. range in the NeXt MD (element 9 & 10)
      if ( ELEMNT(9) .eq. " " ) then
        nxminE = minE
      else
        read(ELEMNT(9),*)nxminE
      endif
      if ( ELEMNT(10) .eq. " " ) then
        nxmaxE = maxE
      else
        read(ELEMNT(10),*)nxmaxE
      endif
      if ( METHOD.eq."CANO" .or. METHOD.eq."MULT" ) then
        write(6,'(2x,a)')"* Energy range for next simulation : "
        write(6,'(8x,f12.3,a,f12.3,a)')nxminE," - ",nxmaxE," (kcal/mol)"
      elseif ( METHOD.eq."CLMD" .or. METHOD.eq."ALSD" ) then
        write(6,'(2x,a)')"* Lambda range for next simulation : "
        write(6,'(8x,f12.3,a,f12.3,a)')nxminE," - ",nxmaxE
      elseif ( METHOD.eq."CLMD2" .or. METHOD.eq."ALSD2" ) then
        write(6,'(2x,a)')"* Lambda square range for next simulation : "
        write(6,'(8x,f12.3,a,f12.3,a)')nxminE," - ",nxmaxE
      endif
      if ( nxminE .ge. nxmaxE ) call error(10210)
      rtmp = (maxE - minE) / dble(nBIN)
      nBOR = ( (max(nxmaxE,maxE) - min(nxminE,minE)) +                 &
             max(abs(maxE-nxmaxE),abs(minE-nxminE)) ) / rtmp
      nAB = nBOR*2 + nBIN

      ! PREvious fitting data file name (element 11)
      if ( METHOD.eq."MULT" .or. METHOD(1:4).eq."ALSD" ) then
        if ( ELEMNT(11) .eq. " " ) call error(10211)
        PREFIT = ELEMNT(11)
        inquire(file=trim(PREFIT),exist=ex)
        if ( .not. ex ) call error(10212)
        write(6,'(2x,a)')"* Previous fitting data : "
        write(6,'(8x,a)')trim(PREFIT)
      endif

      ! TeMPerture for SaMPling (element 12)
      if ( METHOD .eq. "MULT" ) then
        smpT = -1.d0
        if ( ELEMNT(12) .ne. " " ) then
          read(ELEMNT(12),*)smpT
          write(6,'(2x,a)')"* Temperature for canonical sampling : "
          write(6,'(8x,f10.3,a)')smpT," (K)"
          if ( smpT .lt. 0.d0 ) call error(10213)
        endif
      endif

      ! Next simulation temperature (element 13)
      if ( ELEMNT(13) .eq. " " ) then
        nxT = T
      else
        read(ELEMNT(13),*)nxT
        write(6,'(2x,a)')"* Next simu. Temper. (if you change it) : "
        write(6,'(8x,f10.3,a)')nxT," (K)"
        if ( nxT .lt. 0.d0 ) call error(10214)
      endif

      if ( METHOD(1:4) .eq. "ALSD" ) then
        ! lambda for reweighting (element 14 & 15)
        lmd_rewei = 1.d0
        if ( ELEMNT(14) .ne. " " ) then
          read(ELEMNT(14),*)lmd_rewei
          lmd_rewei2 = lmd_rewei*lmd_rewei
        elseif ( ELEMNT(15) .ne. " " ) then
          read(ELEMNT(15),*)lmd_rewei2
          lmd_rewei = sqrt(lmd_rewei2)
        endif

        ! min & max values for reweighting of ALSD (element 16-23)
        if ( ELEMNT(16).ne." " .and. ELEMNT(17).ne." " .and.           &
             ELEMNT(18).ne." " .and. ELEMNT(19).ne." " .and.           &
             ELEMNT(20).ne." " .and. ELEMNT(21).ne." " ) then
          reweight_flg = .true.
          read(ELEMNT(16),*)minlmd ; read(ELEMNT(17),*)maxlmd
          read(ELEMNT(18),*)minEAA ; read(ELEMNT(19),*)maxEAA
          read(ELEMNT(20),*)minEAB ; read(ELEMNT(21),*)maxEAB
          write(6,'(2x,a,f8.3,a)')"* An ensemble at lambda = ",        &
                                  lmd_rewei," is reweighted"
          write(6,'(2x,a)')"* Range of some values for reweighting ALSD"
          write(6,'(8x,2(a,f15.3))')"LAMBDA : ",minlmd," - ",maxlmd
          write(6,'(8x,2(a,f15.3))')"   EAA : ",minEAA," - ",maxEAA
          write(6,'(8x,2(a,f15.3))')"   EAB : ",minEAB," - ",maxEAB
          if ( minlmd.gt.maxlmd .or. minEAA.gt.maxEAA .or.             &
             minEAB.gt.maxEAB ) call error(10215)

          nlmd = 20 ! 100 = # of bins for lambda
          !nlmd = 40 ! 100 = # of bins for lambda
          slmd = (maxlmd-minlmd) / nlmd
          !sEAA = 4.d0
          sEAA = 1.d0
          rtmp = (maxEAA-minEAA)/sEAA
          nEAA = idint(rtmp)
          if ( abs(rtmp-nEAA) .gt. 0.0001d0 ) nEAA = nEAA + 1
          !sEAB = 4.5d0
          sEAB = 1.d0
          rtmp = (maxEAB-minEAB)/sEAB
          nEAB = idint(rtmp)
          if ( abs(rtmp-nEAB) .gt. 0.0001d0 ) nEAB = nEAB + 1
          nlmd=nlmd+1 ; nEAA=nEAA+1 ; nEAB=nEAB+1
          write(6,'(2x,a)')"* Bin size & number for reweighting ALSD"
          write(6,'(8x,a,f15.3,a3,i6)')"LAMBDA : ",slmd," x ",nlmd
          write(6,'(8x,a,f15.3,a3,i6)')"   EAA : ",sEAA," x ",nEAA
          write(6,'(8x,a,f15.3,a3,i6)')"   EAB : ",sEAB," x ",nEAB
        endif
      endif

      ! Acceleration by flatness & the rate (element 24 & 25)
      if ( ELEMNT(22)(1:1).eq."Y" .or. ELEMNT(22)(1:1).eq."y" ) then
        accele = .true.
      else
        accele = .false.
      endif
      if ( accele .and. (METHOD.eq."MULT" .or. METHOD(1:4).eq."ALSD") ) then
        write(6,'(2x,a)')"* Acceleration by flatness is available"
        read(ELEMNT(23),*)accrat
        write(6,'(8x,a,f8.3)')"scaling factor = ",accrat
      endif

      ! Acceleration for the whole region & the rate (element 26 & 27)
      if ( ELEMNT(24)(1:1).eq."Y" .or. ELEMNT(24)(1:1).eq."y" ) then
        accel2 = .true.
      else
        accel2 = .false.
      endif
      if ( accel2 .and. (METHOD.eq."MULT" .or. METHOD(1:4).eq."ALSD") ) then
        write(6,'(2x,a)')"* Acceleration for the whole region is"//    &
                         " available"
        read(ELEMNT(25),*)accrt2
        write(6,'(8x,a,f8.3)')"scaling factor = ",accrt2
      endif

      ! Acceleration for LAB (element 28)
      if ( ELEMNT(26)(1:1).eq."Y" .or. ELEMNT(26)(1:1).eq."y" ) then
        acclab = .true.
        read(ELEMNT(27),*)accrtl
        if ( accrtl .le. 0.d0 ) accrtl = 1.d0
        read(ELEMNT(28),*)accthl
        if ( accthl .le. 0.d0 ) accthl = 1.d0
      else
        acclab = .false.
      endif
      if ( acclab .and. METHOD(1:4).eq."ALSD" ) then
        write(6,'(2x,a)')"* Acceleration for LAB is available"
        write(6,'(8x,a,f8.3)')"scaling factor = ",accrtl
        write(6,'(8x,a,f8.3,a)')"threshold = ",accthl," sigma"
      endif

      ! Flag of neglect out of the range
      if ( ELEMNT(29)(1:1).eq."Y" .or. ELEMNT(29)(1:1).eq."y" .or.     &
           ELEMNT(29)(1:4).eq."BOTH" .or.                              &
           ELEMNT(29)(1:4).eq."both") then
        write(6,'(2x,a)')"* fitting of dlnP for out of the range is"// &
                         " ignored"
        FneglectH = .true. ; FneglectL = .true.
      elseif ( ELEMNT(29)(1:3).eq."LOW" .or.                           &
               ELEMNT(29)(1:3).eq."low" ) then
        write(6,'(2x,a)')"* fitting of dlnP for out of the LOWER "//   &
                         "range is ignored"
        FneglectH = .false. ; FneglectL = .true.
      elseif ( ELEMNT(29)(1:4).eq."HIGH" .or.                          &
               ELEMNT(29)(1:4).eq."high" ) then
        write(6,'(2x,a)')"* fitting of dlnP for out of the HIGHER "//  &
                         "range is ignored"
        FneglectH = .true. ; FneglectL = .false.
      else
        FneglectH = .false. ; FneglectL = .false.
      endif

      ! force size for out of the range for ALSD
      if ( METHOD(1:4).eq."ALSD" .or. METHOD(1:4).eq."CLMD" ) then
        forcel = 0.d0 ; forceh = 0.d0
        if ( ELEMNT(31) .ne. " " ) read(ELEMNT(31),*)forcel
        if ( ELEMNT(32) .ne. " " ) read(ELEMNT(32),*)forceh
        read(ELEMNT(30),*)forces
        if ( forcel .eq. 0.d0 ) forcel = forces
        if ( forceh .eq. 0.d0 ) forceh = forces
        if ( forcel .eq. forceh ) then
          write(6,'(2x,a,f12.5)')"* force size for out of the range  " &
                         //"for ALSD = ",forcel
        else
          write(6,'(2x,a,f12.5)')"* force size of under the range  " &
                         //"for ALSD = ",forcel
          write(6,'(2x,a,f12.5)')"* force size of over the range  " &
                         //"for ALSD = ",forceh
        endif
      endif

      ! ignore nodata bin flag
      if ( METHOD(1:4).eq."MULT" .or. METHOD(1:4).eq."ALSD") then
        if ( ELEMNT(33)(1:1).eq."Y" .or. ELEMNT(33)(1:1).eq."y" ) then
          igndat = .true.
          write(6,'(2x,a)')"* energy or lambda bins with no data are"//&
                           " ignored in the dlnN fitting"
        else
          igndat = .false.
          write(6,'(2x,a)')"* energy or lambda bins with no data are"//&
                           " used in the dlnN fitting"
        endif 
      endif 

      ! adjust energy for ALSD simulation
      if ( reweight_flg ) then
        if ( ELEMNT(34)(1:1).eq."Y" .or. ELEMNT(34)(1:1).eq."y" ) then
          adjene = .true.
          write(6,'(2x,a)')"* Energy ranges for ALSD reweighting is "//&
            "automatically set."
        endif
      endif

!***************************

      return
      end subroutine store_element
