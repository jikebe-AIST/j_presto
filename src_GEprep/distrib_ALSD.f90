
      subroutine distrib_ALSD

!*********************************************
!
!     DISTRIBution_ADAPTIVE_LAMBDA_DYNAMICS
!
!*********************************************

      use COMIFN ; use COMVAL

      implicit none

      integer(4):: i,j,k,m,ios,ier,itmp,itmp2,iset,Nwindow_old,        &
        ndim_old,iminEAA,imaxEAA,iminEAB,imaxEAB,iEAA,iEAB,iST,iEN
      integer(8):: l,n
      real(8):: EbinSZ,minEc,EE(nAB),rtmp,rtmp2,rtmp3,rtmp4,lnN,       &
        minLMD2,maxLMD2,minEAA2,maxEAA2,minEAB2,maxEAB2,lowd_old,      &
        upd_old,TminE,TmaxE,lower_old,upper_old,low_v_window(Nwindow), &
        high_v_window(Nwindow),lowd,lower,upd,upper,dx,dl
      integer(4),allocatable:: TTPdist(:,:)
      integer(8),allocatable:: Pc(:,:,:)
      real(8),allocatable:: aveEAA(:),aveEAB(:),ELMD(:),EEAA(:),       &
        EEAB(:),flt(:),lmdrng(:,:),eP(:,:),F(:),rF(:),w(:),P(:),lnP(:),&
        COEold(:,:),low_old(:),high_old(:),ttw(:,:,:),weit(:,:),       &
        lnPc(:),tempF(:),COE(:,:),tEE(:),wei(:)
      character(9999):: line
      logical(1):: fwflg

!*****************************
!     Initialization

      EbinSZ = (maxE - minE) / dble(nBIN)
      minEc = minE - EbinSZ*dble(nBOR)
      forall ( i=1:nAB ) EE(i) = minEc + EbinSZ*(dble(i-1)+0.5d0)
      itmp = nBOR+1 ; itmp2 = nBOR+nBIN ; fwflg = accele

      ! Check number of input files
      i = 0
      open(unit=1,file=trim(DATLST),status="old")
      do
        read(1,'(a)',iostat=ios)line
        if ( ios .ne. 0 ) exit
        i = i + 1
      enddo
      close(1) ; i = i / 2
      allocate(flt(i),lmdrng(2,i),eP(nAB,i))

      ! Heap memory
      allocate(aveEAA(nAB),aveEAB(nAB),TTPdist(3,nAB),F(nAB),rF(nAB),  &
               w(nAB))
9999  aveEAA(:) = 0.d0 ; aveEAB(:) = 0.d0 ; l = 0 ; n = 0
      eP(:,:) = 0.d0 ; iset = 0 ; TTPdist(:,:) = 0 ; F(:) = 0.d0
      rF(:) = 0.d0 ; w(:) = 0.d0
      if ( reweight_flg ) then
        allocate(Pc(nlmd,nEAA,nEAB),ELMD(nlmd),EEAA(nEAA),EEAB(nEAB))
        Pc(1:nlmd,1:nEAA,1:nEAB) = 0
        forall ( i=1:nlmd ) ELMD(i) = minlmd + slmd*(dble(i-1)+0.5d0)
        forall ( i=1:nEAA ) EEAA(i) = minEAA + sEAA*(dble(i-1)+0.5d0)
        forall ( i=1:nEAB ) EEAB(i) = minEAB + sEAB*(dble(i-1)+0.5d0)
        minLMD2 = huge(0.d0) ; maxLMD2 = -huge(0.d0)
        minEAA2 = huge(0.d0) ; maxEAA2 = -huge(0.d0)
        minEAB2 = huge(0.d0) ; maxEAB2 = -huge(0.d0)
      endif

!*****************************
!     Input Lambda data

      open(unit=1,file=trim(DATLST),status="old")
      open(unit=4,file=trim(PROJNM)//".bar",status="replace")

      if ( reweight_flg ) then
        call input_lambda_data_weight()
      else
        call input_lambda_data()
      endif

      close(1) ; write(6,*)
      write(4,*)minE,0
      write(4,*)minE,iset+1
      write(4,*)
      write(4,*)maxE,0
      write(4,*)maxE,iset+1
      close(4)

      if ( adjene ) then
        adjene = .false.
        if ( minEAA2.lt.minEAA .or. maxEAA2.gt.maxEAA .or.             &
             minEAB2.lt.minEAB .or. maxEAB2.gt.maxEAB ) then
          minEAA = float(floor(minEAA2/5) * 5)
          maxEAA = ceiling(maxEAA2/5.0) * 5.0
          minEAB = float(floor(minEAB2/5) * 5)
          maxEAB = ceiling(maxEAB2/5.0) * 5.0
          rtmp = (maxEAA-minEAA)/sEAA ; nEAA = idint(rtmp)
          if ( abs(rtmp-nEAA) .gt. 0.0001d0 ) nEAA = nEAA + 1
          rtmp = (maxEAB-minEAB)/sEAB ; nEAB = idint(rtmp)
          if ( abs(rtmp-nEAB) .gt. 0.0001d0 ) nEAB = nEAB + 1
          nEAA = nEAA + 1 ; nEAB = nEAB + 1
          write(6,*)
          write(6,'(2x,a)')"* Bin size recalculation"
          write(6,'(8x,2(a,f15.3))')"LAMBDA : ",minlmd," - ",maxlmd
          write(6,'(8x,2(a,f15.3))')"   EAA : ",minEAA," - ",maxEAA
          write(6,'(8x,2(a,f15.3))')"   EAB : ",minEAB," - ",maxEAB
          write(6,'(2x,a)')"* Bin size & number for reweighting ALSD"
          write(6,'(8x,a,f15.3,a3,i6)')"LAMBDA : ",slmd," x ",nlmd
          write(6,'(8x,a,f15.3,a3,i6)')"   EAA : ",sEAA," x ",nEAA
          write(6,'(8x,a,f15.3,a3,i6)')"   EAB : ",sEAB," x ",nEAB
          deallocate(Pc,ELMD,EEAA,EEAB)
          goto 9999
        endif
      endif
      allocate(P(nAB)) ; P(:) = sum(eP, dim=2)

!*******************************
      ! Make distribution data

      rtmp = log(maxval(P)) ; allocate(lnP(nAB))
      lnP = merge(1.d0, log(P)-rtmp, P.eq.0.d0)
      write(6,'(2x,a,i0,a,i0,a,f8.3,a)')"+ ",n," / ",l," (",           &
           dble(n)/dble(l)*100.d0," (%)) data in the range are counted"
      write(6,*)
      rtmp = sum(exp(lnP(itmp:itmp2)),mask = lnP(itmp:itmp2).ne.1.d0 )
      rtmp = rtmp / dble(nBIN) * 100.d0
      write(6,'(2x,a,f8.3,a)')"+ Flatness              = ",rtmp," (%)"
      ! flatness for individual distribution
      rtmp = sum(flt(1:iset)) / dble(iset)
      write(6,'(2x,a,f8.3,a)')"+ Ave. Flat. for indiv. = ",rtmp," (%)"
      i = minloc(flt(1:iset), dim=1)
      write(6,'(6x,a,f8.3,a,i0,a)')"min. = ",flt(i)," (%) (set ",i,")"
      i = maxloc(flt(1:iset), dim=1)
      write(6,'(6x,a,f8.3,a,i0,a)')"max. = ",flt(i)," (%) (set ",i,")"

      ! lambda range for individual distribution
      rtmp = sum(lmdrng(1,1:iset)) / dble(iset)
      write(6,'(2x,"+ Ave. λ-range (full) for indiv. = ",f8.3)')rtmp
      i = minloc(lmdrng(1,1:iset), dim=1)
      write(6,'(2x,"    min. = ",f8.3," (set ",i0,")")') lmdrng(1,i), i
      i = maxloc(lmdrng(1,1:iset), dim=1)
      write(6,'(2x,"    max. = ",f8.3," (set ",i0,")")') lmdrng(1,i), i
      write(6,*)
      ! lambda range for individual distribution (subrange: lnP ≥ -1.5)
      rtmp = sum(lmdrng(2,1:iset)) / dble(iset)
      write(6,'(2x,"+ Ave. λ-range (sub: lnP ≥ -1.5) for indiv. = ",   &
            f8.3)') rtmp
      i = minloc(lmdrng(2,1:iset), dim=1)
      write(6,'(8x,"min. = ",f8.3," (set ",i0,")")')lmdrng(2,i),i
      i = maxloc(lmdrng(2,1:iset), dim=1)
      write(6,'(8x,"max. = ",f8.3," (set ",i0,")")')lmdrng(2,i),i
      write(6,*)

      ! Output distribution
      open(unit=3,file=trim(PROJNM)//".dist",status="replace")
      do i = 1,nAB
        if ( lnP(i) .ne. 1.d0 ) write(3,*)EE(i),lnP(i)
      enddo
      write(3,*)
      write(3,*)minE,0.d0
      write(3,*)minE,-3.d0
      write(3,*)
      write(3,*)maxE,0.d0
      write(3,*)maxE,-3.d0
      close(3)

      ! Individual distribution
      do i = 1,iset
        rtmp = maxval(eP(:, i))
        if ( rtmp .ne. 0.d0 ) then
          rtmp = log(rtmp)
          where ( eP(:,i) .eq. 0.d0 )
            eP(:,i) = 1.d0
          elsewhere
            eP(:,i) = log(eP(:,i)) - rtmp
          end where
        endif
      end do

      iST = 1 ; iEN = nAB
      do i = 1,nAB
        if ( any(eP(i,:) .ne. 1.d0) ) then
          if ( iST .eq. 1 ) iST = i
          iEN = i
        endif
      enddo
      open(unit=3,file=trim(PROJNM)//".edist",status="replace")
      do i = iST,iEN
        write(3,'(f10.5)',advance='no')EE(i)
        do j = 1,iset
          if ( eP(i,j) .ne. 1.d0 ) then
            write(3,'(f10.5)',advance='no')eP(i,j)
          else
            write(3,'(f10.5)',advance='no')0.d0/0.d0
          endif
        enddo
        write(3,'(a)')''
      enddo
      close(3)

      if ( reweight_flg ) then
        write(6,'(2x,a)')"* Obtained min & max lambda or the square"
        write(6,'(8x,f10.3,a,f10.3)')minLMD2," - ",maxLMD2
        write(6,'(2x,a)')                                              &
          "* Selected min & max lambda or the square for CANO dist."
        write(6,'(8x,f10.3,a,f10.3)')minlmd," - ",maxlmd
        write(6,'(2x,a)')                                              &
            "* min & max EAA & EAB in the selected lrange"
        write(6,'(8x,a,f15.3,a,f15.3)')"EAA : ",minEAA2," - ",maxEAA2
        write(6,'(8x,a,f15.3,a,f15.3)')"EAB : ",minEAB2," - ",maxEAB2
        write(6,*)

        write(6,'(2x,a)')"* N of conf. information for canonical dist."
        write(6,'(8x,a,i0)')  "+                Max : ",maxval(Pc)
        i = count(Pc.ne.0)
        if ( i .gt. 0 ) then
          write(6,'(8x,a,f10.3)')"+                Ave : ",            &
            dble(sum(Pc))/dble(i)
        endif
        write(6,'(8x,a,i0)')  "+ N of non-zero bins : ",i

        ! Input previous dlnP
        open(unit=1,file=trim(PREFIT),status="old")
        read(1,*)Nwindow_old,ndim_old
        allocate(COEold(0:ndim_old,Nwindow_old),low_old(Nwindow_old),  &
                 high_old(Nwindow_old))
        do i = 1,Nwindow_old
          do j = 0,Ndim_old
            read(1,*)COEold(j,i)
          enddo
          read(1,*)low_old(i),high_old(i)
        enddo
        read(1,*)lowd_old,upd_old
        read(1,*)TminE,TmaxE
        read(1,*)
        read(1,*)lower_old,upper_old
        close(1)

        allocate(ttw(nEAA,nEAB,nlmd)) ; ttw(:,:,:) = 0.d0
        m = 1 ; rtmp4 = 1.d0 / (R*T) ; lnN = 0.d0
        do i = 1, nlmd
          rtmp3 = ELMD(i)
          do
            if ( m.ne.Nwindow_old .and. rtmp3.ge.high_old(m) ) then
              m = m + 1 ; cycle
            endif
            exit
          enddo
          rtmp = rtmp3 - low_old(m)
          rtmp2 = COEold(0,m)
          do j = 1, ndim_old
            rtmp2 = rtmp2 + COEold(j,m) * rtmp**j
          enddo
          if ( rtmp3 .lt. TminE ) then
            rtmp2 = rtmp2 + lower_old * (TminE - rtmp3)
          elseif ( rtmp3 .gt. TmaxE ) then
            rtmp2 = rtmp2 + upper_old * (rtmp3 - TmaxE)
          endif
      
          lnN = lnN + rtmp2 * slmd
          do k = 1, nEAB
          do j = 1, nEAA
            select case (METHOD)
            case ("ALSD")
              rtmp = ((rtmp3**2 - lmd_rewei2) * EEAA(j) +  &
                      (rtmp3 - lmd_rewei)   * EEAB(k)) * rtmp4
            case ("ALSD2")
              rtmp = ((rtmp3 - lmd_rewei2) * EEAA(j) +     &
                      (sqrt(rtmp3) - lmd_rewei) * EEAB(k)) * rtmp4
            end select
            ttw(j,k,i) = -(rtmp + lnN)
          enddo
          enddo
        enddo

        rtmp = 0.d0 ; allocate(weit(nEAA,nEAB)) ; weit(:,:) = 0.d0
        if ( any(ttw .ne. 0.d0) ) rtmp = maxval(ttw, mask = ttw.ne.0.d0)
        do i = 1,nlmd
          rtmp2 = dble(sum(Pc(i,1:nEAA,1:nEAB)))
          forall ( j=1:nEAA, k=1:nEAB, ttw(j,k,i).ne.0.d0 )            &
            weit(j,k) = weit(j,k) + exp(ttw(j,k,i)-rtmp)*rtmp2
        enddo
        deallocate(ttw)

        ! lambda distribution
        allocate(lnPc(nlmd)) ; lnPc(:) = 0.d0
        do k = 1,nEAB
        do j = 1,nEAA
          if ( weit(j,k) .eq. 0.d0 ) cycle
          rtmp2 = 1.d0 / weit(j,k)
          do i = 1,nlmd
            rtmp = dble(Pc(i,j,k))
            if ( rtmp .ne. 0.d0 ) lnPc(i) = lnPc(i) + rtmp*rtmp2
          enddo
        enddo
        enddo
        rtmp = log(maxval(lnPc))
        lnPc = merge(1.d0, log(lnPc)-rtmp, lnPc.eq.0.d0)

        open(unit=3,file=trim(PROJNM)//".ldist",status="replace")
        do i = 1,nlmd
          if ( lnPc(i) .ne. 1.d0 ) write(3,*)ELMD(i),lnPc(i)
        enddo
        write(3,*)
        write(3,*)minlmd,0.d0
        write(3,*)minlmd,-3.d0
        write(3,*)
        write(3,*)maxlmd,0.d0
        write(3,*)maxlmd,-3.d0
        close(3)

        ! Output probability file
        weit = merge(1.d0 / weit, weit, weit.ne.0.d0)
        rtmp = 1.d0 / maxval(weit(:,:))
        weit(:,:) = weit(:,:) * rtmp
        rtmp = exp(-7.d0)   ! threshold
        !! Region check
        iminEAA = 1 ; iminEAB = 1 ; imaxEAA = nEAA ; imaxEAB = nEAB
        do i = 1,nEAA
          if ( maxval(weit(i,1:nEAB)) .ge. rtmp ) then
            iminEAA = i ; exit
          endif
        enddo
        do i = nEAA,1,-1
          if ( maxval(weit(i,1:nEAB)) .ge. rtmp ) then
            imaxEAA = i ; exit
          endif
        enddo
        iEAA = imaxEAA - iminEAA + 1
        do i = 1,nEAB
          if ( maxval(weit(1:nEAA,i)) .ge. rtmp ) then
            iminEAB = i ; exit
          endif
        enddo
        do i = nEAB,1,-1
          if ( maxval(weit(1:nEAA,i)) .ge. rtmp ) then
            imaxEAB = i ; exit
          endif
        enddo
        iEAB = imaxEAB - iminEAB + 1
        open(unit=3,file=trim(PROJNM)//".prob",status="replace")
        write(3,*)"ALSD"
        write(3,*)iEAA,iEAB                      ! Number of bins
        write(3,*)sEAA,sEAB                      ! Bin size
        write(3,*)minEAA+sEAA*dble(iminEAA-1),                         &
                  minEAB+sEAB*dble(iminEAB-1)   ! Min value
        do j = iminEAB,imaxEAB
          write(3,*)
          do i = iminEAA,imaxEAA
            if ( weit(i,j) .ge. rtmp ) then
              write(3,*)weit(i,j)
            else
              write(3,*)0.d0
            endif
          enddo
        enddo
        close(3)

        open(unit=3,file=trim(PROJNM)//".cprob",status="replace")
        do j = iminEAB,imaxEAB
          write(3,*)
          rtmp2 = minEAB+dble(j-0.5d0)*sEAB
          do i = iminEAA,imaxEAA
            rtmp3 = minEAA+dble(i-0.5d0)*sEAA
            write(3,'(3(e15.8,x),i8)')                                 &
              rtmp3,rtmp2,weit(i,j)*dble(sum(Pc(:,i,j))),sum(Pc(:,i,j))
          enddo
        enddo
        close(3)

        write(6,'(2x,a)')"* min & max EAA & EAB for reweighting"
        write(6,'(8x,a,f10.3,a,f10.3)')"EAA : ",EEAA(max(iminEAA-1,1)),&
                                       " - ",EEAA(min(imaxEAA+1,nEAA))
        write(6,'(8x,a,f10.3,a,f10.3)')"EAB : ",EEAB(max(iminEAB-1,1)),&
                                       " - ",EEAB(min(imaxEAB+1,nEAB))
        write(6,*)

      endif
  
!*****************************
      ! Calc. Differential Density of States
      rtmp = 1.d0 / (R * T)
      do i = 1,nAB
        if ( P(i) .eq. 0 ) cycle
        rF(i) = -rF(i) / P(i) * rtmp
        F(i) = -F(i) / w(i) * rtmp
      enddo

      ! Output dlnP Function (*.nf)
      allocate(tEE(nAB),tempF(nAB),COE(0:NfitDIM,Nwindow),wei(nAB))
      j = 0
      if ( FneglectL ) then
        iST = nBOR+1
      else
        iST = 1
      endif
      if ( FneglectH ) then
        iEN = nBOR+nBIN
      else
        iEN = nAB-1
      endif
      do i = iST,iEN
        if ( P(i) .eq. 0.d0 ) cycle
        j = j + 1 ; tEE(j) = EE(i) ; tempF(j) = F(i)
        wei(j) = sqrt(P(i))
      enddo

      write(6,'(2x,a)')"+ Fitting for dlnP"
      call lstsquare_final(Nwindow,j,NfitDIM,nxminE,nxmaxE,tEE(1:j),   &
        tempF(1:j),wei(1:j),COE,low_v_window,high_v_window,rtmp,rtmp2, &
        .true.)

      open(unit=1,file=trim(PROJNM)//".nf",status="replace")
      write(1,*)Nwindow,NfitDIM
      do i = 1,Nwindow
        do j = 0,NfitDIM
          write(1,*)COE(j,i)
        enddo
        write(1,*)low_v_window(i),high_v_window(i)
      enddo

      ! Calc. forces for the out of the range
      !! lower side
      lowd = COE(0,1)
      lower = -forcel
      !! upper side
      upd = COE(0,Nwindow) ; rtmp = nxmaxE - low_v_window(Nwindow)
      do i = 1,NfitDIM
        upd = upd + COE(i,Nwindow)*rtmp**i
      enddo
      upper = forceh
      write(1,*)lowd,upd
      write(1,*)nxminE,nxmaxE
      write(1,*)nxT
      write(1,*)lower,upper
      close(1)
      write(6,*)

!*****************************
!     Output dlnP data for plot

      iST = 1 ; iEN = nAB
      open(unit=1,file=trim(PROJNM)//".dlnP",status="replace")
      m = 1
      do i = iST,iEN
        if ( P(i) .eq. 0 ) cycle
        rtmp2 = EE(i)
        do
          if ( m.ne.Nwindow .and. rtmp2.ge.high_v_window(m) ) then
            m = m + 1 ; cycle
          endif
          exit
        enddo
        dx = rtmp2 - low_v_window(m) ; rtmp = COE(0,m)
        do j = 1,NfitDIM
          rtmp = rtmp + COE(j,m)*dx**j
        enddo
        if ( rtmp2 .lt. nxminE ) then
          dl = nxminE - rtmp2
          rtmp3 = rtmp + lower*dl
        elseif ( rtmp2 .gt. nxmaxE ) then
          dl = rtmp2 - nxmaxE
          rtmp3 = rtmp + upper*dl
        else
          rtmp3 = rtmp
        endif
        if ( fwflg ) then
          write(1,'(5(e15.8,x))')rtmp2,F(i),rtmp,rtmp3,rF(i)
        else
          write(1,'(4(e15.8,x))')rtmp2,F(i),rtmp,rtmp3
        endif
      enddo
      write(1,*)
      write(1,*)nxminE,lowd,lowd,lowd
      write(1,*)nxminE,upd,upd,upd
      write(1,*)
      write(1,*)nxmaxE,lowd,lowd,lowd
      write(1,*)nxmaxE,upd,upd,upd
      close(1)

!*****************************
!     Output TTPdist

      open(unit=1,file=trim(PROJNM)//".TTPdist",status="replace")
      do i = 1,nAB
        rtmp2 = EE(i) + 0.5d0*EbinSZ
        if ( TTPdist(1,i).ne.0 ) write(1,*)rtmp2,TTPdist(1:3,i)
      enddo
      close(1)

!*****************************

      return
      contains


!=======================================================================


      subroutine input_lambda_data_weight()

      use COMIFN ; use COMVAL

      implicit none

      integer(4):: ii,jj,k1,k2,k3,iST,iEN,tP(nAB),ttP(nAB)
      real(8):: fw,E,EAA,EAB,ttF(nAB),tF(nAB),tw(nAB)
      logical(1):: onend,ex
      character(9999):: filnam

!********************************************

      tP(:) = 0 ; tF(:) = 0.d0 ; tw(:) = 0.d0
      do
        call infil(1,filnam,iST,iEN,fw,onend,ier)
        if ( ier .ne. 0 ) exit
        inquire(file=trim(filnam),exist=ex)

        if ( ex ) then
          write(6,'(2x,a)')"+ Now input "//trim(filnam)
          if ( fw .ne. 1.d0 ) fwflg = .true.
        
          open(unit=2,file=trim(filnam),status="old")
          do ii = 1,iST-1
            read(2,*,end=800)
          enddo
          ttF(:) = 0.d0 ; ttP(:) = 0
          do ii = iST,iEN
            read(2,*,end=800)E,EAA,EAB
            jj = int((E - minEc) / EbinSZ + 1) ; l = l + 1
            if ( jj.ge.1 .and. jj.le.nAB ) then
              ttP(jj) = ttP(jj) + 1
              aveEAA(jj) = aveEAA(jj) + EAA
              aveEAB(jj) = aveEAB(jj) + EAB
              ttF(jj) = ttF(jj) + 2.d0 * EAA * E + EAB
              if ( jj.ge.itmp .and. jj.le.itmp2 ) n = n + 1
            endif

            k1 = int((E - minlmd) / slmd + 1)
            minLMD2 = min(minLMD2,E) ; maxLMD2 = max(maxLMD2,E)
            if ( k1.ge.1 .and. k1.le.nlmd ) then
              k2 = int((EAA - minEAA) / sEAA + 1)
              k3 = int((EAB - minEAB) / sEAB + 1)
              minEAA2 = min(minEAA2,EAA) ; maxEAA2 = max(maxEAA2,EAA)
              minEAB2 = min(minEAB2,EAB) ; maxEAB2 = max(maxEAB2,EAB)
              if ( k2.ge.1 .and. k2.le.nEAA .and.                      &
                   k3.ge.1 .and. k3.le.nEAB )                          &
                   Pc(k1,k2,k3) = Pc(k1,k2,k3) + 1
            endif
          enddo
800       close(2)

          if ( fw .eq. 1.d0 ) then
            write(6,'(8x,i0,a3,i0)')iST," - ",ii-1
          else
            write(6,'(8x,i0,a3,i0,a5,f8.3,a1)')                        &
              iST," - ",ii-1," ( x ",fw,")"
          endif
          rF(:) = rF(:) + ttF(:)
          tF(:) = tF(:) + ttF(:) * fw
          tw(:) = tw(:) + dble(ttP(:)) * fw
          tP(:) = tP(:) + ttP(:)
        endif
        if ( onend ) call onend_data(tP,tF,tw)
      enddo

      return
      end subroutine input_lambda_data_weight
        

!========================================================================


      subroutine input_lambda_data()

      use COMIFN ; use COMVAL

      implicit none

      integer(4):: ii,jj,iST,iEN,tP(nAB),ttP(nAB)
      real(8):: fw,E,EAA,EAB,ttF(nAB),tF(nAB),tw(nAB)
      logical(1):: onend,ex
      character(9999):: filnam

!********************************************

      tP(:) = 0 ; tF(:) = 0.d0 ; tw(:) = 0.d0
      do
        call infil(1,filnam,iST,iEN,fw,onend,ier)
        if ( ier .ne. 0 ) exit
        inquire(file=trim(filnam),exist=ex)

        if ( ex ) then
          write(6,'(2x,a)')"+ Now input "//trim(filnam)
          if ( fw .ne. 1.d0 ) fwflg = .true.
        
          open(unit=2,file=trim(filnam),status="old")
          do ii = 1,iST-1
            read(2,*,end=800)
          enddo
          ttF(:) = 0.d0 ; ttP(:) = 0
          do ii = iST,iEN
            read(2,*,end=800)E,EAA,EAB ; l = l + 1
            jj = int((E - minEc) / EbinSZ + 1)
            if ( jj.ge.1 .and. jj.le.nAB ) then
              ttP(jj) = ttP(jj) + 1
              aveEAA(jj) = aveEAA(jj) + EAA
              aveEAB(jj) = aveEAB(jj) + EAB
              ttF(jj) = ttF(jj) + 2.d0 * EAA * E + EAB
              if ( jj.ge.itmp .and. jj.le.itmp2 ) n = n + 1
            endif
          enddo
800       close(2)

          if ( fw .eq. 1.d0 ) then
            write(6,'(8x,i0,a3,i0)')iST," - ",ii-1
          else
            write(6,'(8x,i0,a3,i0,a5,f8.3,a1)')                        &
              iST," - ",ii-1," ( x ",fw,")"
          endif
          rF(:) = rF(:) + ttF(:)
          tF(:) = tF(:) + ttF(:) * fw
          tw(:) = tw(:) + dble(ttP(:)) * fw
          tP(:) = tP(:) + ttP(:)
        endif
        if ( onend ) call onend_data(tP,tF,tw)
      enddo

      return
      end subroutine input_lambda_data
        

!========================================================================


      subroutine onend_data(tP,tF,tw)

      use COMIFN ; use COMVAL

      implicit none

      integer(4),intent(inout)::tP(nAB)
      real(8),intent(inout):: tF(nAB),tw(nAB)

      integer(4):: ii,is,ie,is2,ie2
      real(8):: tlnP(nAB)

!****************************************

      iset = iset + 1 ; eP(:,iset) = dble(tP(:))
      rtmp = maxval(eP(:,iset))
      if ( rtmp .ne. 0.d0 ) then
        ! flatness calc.
        rtmp = log(rtmp)
        tlnP = merge(1.d0, log(eP(:,iset))-rtmp, eP(:,iset).eq.0.d0)
        rtmp = sum(exp(tlnP(itmp:itmp2)), mask = tlnP(itmp:itmp2).ne.1.d0)
        rtmp = rtmp / dble(nBIN) * 100.d0
        write(6,'(2x,a,f8.3,a)')"+ Flatness = ",rtmp," (%)"
        flt(iset) = merge(0.01d0, rtmp, rtmp .eq. 0.d0)

        if ( accele ) then
          rtmp = 1.d0/flt(iset)**accrat
          F(:) = F(:) + tF(:) * rtmp ; w(:) = w(:) + tw(:) * rtmp
        else
          F(:) = F(:) + tF(:) ; w(:) = w(:) + tw(:)
        endif

        ! For TTPdist
        do ii = 1,nAB
          if ( tlnP(ii) .le. 0.d0 ) then
            TTPdist(1,ii) = TTPdist(1,ii) + 1
            if ( tlnP(ii) .ge. -1.5d0 ) then
              TTPdist(2,ii) = TTPdist(2,ii) + 1
              if ( tlnP(ii) .ge. -0.5d0 )                              &
                TTPdist(3,ii) = TTPdist(3,ii) + 1
            endif
          endif
        enddo

        do i = 1,nAB
          if ( tlnP(i) .ne. 1.d0 ) then
            is = i ; exit
          endif
        enddo
        do i = is,nAB
          if ( tlnP(i) .ge. -1.5d0 ) then
            is2 = i ; exit
          endif
        enddo
        do i = nAB,1,-1
          if ( tlnP(i) .ne. 1.d0 ) then
            ie = i ; exit
          endif
        enddo
        do i = ie,1,-1
          if ( tlnP(i) .ge. -1.5d0 ) then
            ie2 = i ; exit
          endif
        enddo
        write(4,'(e12.5,x,i0,x,e12.5)')EE(is),iset,EE(is2)
        write(4,'(e12.5,x,i0,x,e12.5)')EE(ie),iset,EE(ie2)
        write(4,*)
        lmdrng(1,iset) = EE(ie) - EE(is) + EbinSZ
        lmdrng(2,iset) = EE(ie2) - EE(is2) + EbinSZ
        write(6,'(2x,a,f8.3,a,f8.3,a)')"+ λ-range = ",                 &
          lmdrng(1,iset)," (",lmdrng(2,iset),")"
        write(6,*)
      endif

      tP(:) = 0 ; tF(:) = 0.d0 ; tw(:) = 0.d0

      return
      end subroutine onend_data


!========================================================================

      end subroutine distrib_ALSD
