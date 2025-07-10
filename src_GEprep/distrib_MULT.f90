
      subroutine distrib_MULT

!*********************************************
!
!     DISTRIBution_MULTICANONICAL
!
!*********************************************

      use COMIFN ; use COMVAL

      implicit none

      ! Energy BIN SiZe
        real(8):: EbinSZ,EsclbinSZ,ENEbinSZ
      ! distribution & lnP
        integer(4):: ttP(nAB),TTPdist(3,nAB)
        real(8):: trP(nAB),tP(nAB),rP(nAB),P(nAB),rlnP(nAB),lnP(nAB),  &
                  tlnP(nAB)
      ! (Temporary) Defferential Density of state
        real(8):: rdlnN(nAB-1),dlnN(nAB),TdlnN(nAB)
      ! MINimum Energy for Count
        real(8):: minEc,minEsclc,minENEc
      ! (Temporary) Energy of each bin (EE is for dlnN)
        real(8):: EE(nAB),Escl(nAB),Escl2(nAB),tEE(nAB),EENE(nBIN),    &
                  EENEscl2(nBIN)
      ! FILE NaMe
        character(9999):: filnam
        integer(4):: iST,iEN
      ! COEfficient
        real(8),allocatable:: COEold(:,:),COE(:,:)
        real(8):: lowd,upd
      ! Flatness calc.
        real(8),allocatable:: flt(:),Erng(:,:)
      ! WEIghting factor
        real(8):: wei(nAB)
      ! Canonical distribution (lnPc)
        real(8):: lnPc(nAB)
      ! MIN & MAX Temperature for sampling range
        real(8):: minT,maxT
      ! Temporary Temperature
        real(8):: TT(nBOR+1:nBOR+1+nBIN)

      integer(4):: i,j,k,l,n,system,itmp,itmp2,ifile,iset,ier,ifirst,  &
                   ilast,ifirst2,ilast2,Nwindow_old,ia(1)
      real(8):: fw,E,rr,ss,TminE,TmaxE,rtmp,rtmp2,rtmp3,rtmp4,dx,      &
                minlnP,re,re2
      real(8),allocatable:: low_old(:),high_old(:)
      real(8):: low_v_window(Nwindow),high_v_window(Nwindow)
      character(130):: tmp,dirname
      logical(1):: onend,fwflg

!*****************************
!     Initialization

      fwflg = .false.
      if ( accele ) fwflg = .true.
      EbinSZ = (maxE - minE) / dble(nBIN)
      EsclbinSZ = 1.d0 / dble(nBIN)
      minEc = minE - EbinSZ*dble(nBOR)
      minEsclc = -EsclbinSZ*dble(nBOR)
      rP(1:nAB) = 0.d0 ; P(1:nAB) = 0.d0 ; l = 0 ; n = 0 ; ifile = 0
      forall ( i=1:nAB ) EE(i) = minEc + EbinSZ*(dble(i-1)+0.5d0)

!      write(dirname,'(a)')trim(PROJNM)//"_DIST"
!      call systemqq("mkdir "//trim(dirname)//" >& /dev/null")

      ! Check number of input files
      call systemqq("wc -l "//trim(DATLST)//" > .aho.temp")
      open(unit=1,file=".aho.temp",status="old")
      read(1,*)i
      close(unit=1,status="delete")
      i = i / 2 ; allocate(flt(i),Erng(2,i))

!*****************************
!     Input Energy data

      itmp = nBOR+1 ; itmp2 = nBOR+nBIN
      trP(1:nAB) = 0.d0 ; tP(1:nAB) = 0.d0 ; iset = 0 ; TTPdist(:,:) = 0
      open(unit=1,file=trim(DATLST),status="old")
      open(unit=4,file=trim(PROJNM)//".bar",status="replace")
      do
        call infil(1,filnam,iST,iEN,fw,onend,ier)
        if ( ier .ne. 0 ) goto 801
        i = system("ls "//trim(filnam)//" >& /dev/null")

        if ( i .eq. 0 ) then
          write(6,'(2x,a)')"+ Now input "//trim(filnam)
          if ( fw .ne. 1.d0 ) fwflg = .true.

          ! Input energy file
          ifile = ifile + 1 ; ttP(1:nAB) = 0
          open(unit=2,file=trim(filnam),status="old")
          do j = 1,iST-1
            read(2,*,end=800)
          enddo
          do j = iST,iEN
            read(2,*,end=800)E ; l = l + 1
            k = int((E - minEc) / EbinSZ) + 1
            if ( k.ge.1 .and. k.le.nAB ) then
              ttP(k) = ttP(k) + 1
              if ( k.ge.itmp .and. k.le.itmp2 ) n = n + 1
            endif
          enddo
800       close(2)

          trP(1:nAB) = trP(1:nAB) + dble(ttP(1:nAB))
          if ( fw .eq. 1.d0 ) then
            write(6,'(8x,i0,a3,i0)')iST," - ",j-1
            tP(1:nAB) = tP(1:nAB) + dble(ttP(1:nAB))
          else
            write(6,'(8x,i0,a3,i0,a5,f8.3,a1)')iST," - ",j-1," ( x ",&
              fw,")"
            tP(1:nAB) = tP(1:nAB) + dble(ttP(1:nAB))*fw
          endif
        endif

        if ( onend ) then
          iset = iset + 1
          rP(1:nAB) = rP(1:nAB) + trP(1:nAB)
          if ( maxval(trP) .ne. 0.d0 ) then
            ! Make individual distribution data
            rtmp = maxval(trP) ; rtmp = log(rtmp)
            where ( trP .eq. 0.d0 )
              tlnP = 1.d0
            elsewhere
              tlnP = log(dble(trP)) - rtmp
            endwhere
            ! flatness for individual distribution
            rtmp = 0.d0
            do k = itmp,itmp2
              if ( tlnP(k) .ne. 1.d0 ) rtmp = rtmp + exp(tlnP(k))
            enddo
            rtmp = rtmp / dble(nBIN) * 100.d0
            write(6,'(2x,a,f8.3,a)')"+ Flatness = ",rtmp," (%)"
            if ( rtmp .eq. 0.d0 ) then
              flt(iset) = 0.01d0
            else
              flt(iset) = rtmp
            endif
            ! For TTPdist
            do k = 1,nAB
              if ( tlnP(k) .le. 0.d0 ) then
                TTPdist(1,k) = TTPdist(1,k) + 1
                if ( tlnP(k) .ge. -1.5d0 ) then
                  TTPdist(2,k) = TTPdist(2,k) + 1
                  if ( tlnP(k) .ge. -0.5d0 )                           &
                    TTPdist(3,k) = TTPdist(3,k) + 1
                endif
              endif
            enddo

            if ( accele ) tP(1:nAB) = tP(1:nAB)/flt(iset)**accrat
            P(1:nAB) = P(1:nAB) + tP(1:nAB)

!            write(tmp,'(a,a1,i0,a5)')trim(dirname)//"/"//trim(PROJNM), &
!                                   "_",iset,".dist"
!            open(unit=3,file=trim(tmp),status="replace")
!            do i = 1,nAB
!              if ( tlnP(i) .ne. 1.d0 ) write(3,*)EE(i),tlnP(i)
!            enddo
!            write(3,*)
!            write(3,*)minE,0.d0
!            write(3,*)minE,-3.d0
!            write(3,*)
!            write(3,*)maxE,0.d0
!            write(3,*)maxE,-3.d0
!            close(3)

            ! Make bar graph
            do i = 1,nAB 
              if ( tlnP(i) .ne. 1.d0 ) then
                ifirst = i ; exit 
              endif
            enddo
            do i = ifirst,nAB
              if ( tlnP(i) .ge. -1.5d0 ) then
                ifirst2 = i ; exit
              endif
            enddo
            do i = nAB,1,-1
              if ( tlnP(i) .ne. 1.d0 ) then
                ilast = i ; exit
              endif
            enddo
            do i = ilast,1,-1
              if ( tlnP(i) .ge. -1.5d0 ) then 
                ilast2 = i ; exit
              endif
            enddo
            write(4,'(e12.5,x,i0,x,e12.5)')EE(ifirst),iset,EE(ifirst2)
            write(4,'(e12.5,x,i0,x,e12.5)')EE(ilast),iset,EE(ilast2)
            write(4,*)
            Erng(1,iset) = (EE(ilast)-EE(ifirst)) + EbinSZ
            Erng(2,iset) = (EE(ilast2)-EE(ifirst2)) + EbinSZ
            write(6,'(2x,a,f12.3,a,f12.3,a)')"+ Erange = ",            &
              Erng(1,iset)," (",Erng(2,iset),")"
            write(6,*)
          endif

          trP(1:nAB) = 0.d0 ; tP(1:nAB) = 0.d0
        endif

      enddo
801   close(1) ; write(6,*)
      write(4,*)minE,0
      write(4,*)minE,iset+1
      write(4,*)
      write(4,*)maxE,0
      write(4,*)maxE,iset+1
      close(4)

!****************************
!     Make distribution data

      !! real P
      rtmp = maxval(rP) ; rtmp = log(rtmp)
      where ( rP .eq. 0.d0 )
        rlnP = 1.d0
      elsewhere
        rlnP = log(rP) - rtmp
      endwhere
      !! P using the weighting factor, fw
      rtmp = maxval(P) ; rtmp = log(rtmp)
      where ( P .eq. 0.d0 )
        lnP = 1.d0
      elsewhere
        lnP = log(P) - rtmp
      endwhere
      write(6,'(2x,a,i0,a,i0,a,f8.3,a)')"+ ",n," / ",l," (",           &
           dble(n)/dble(l)*100.d0," (%)) data in the range are counted"
      write(6,*)

      ! flatness
      rtmp = 0.d0 ; rtmp2 = 0.d0
      do i = itmp,itmp2
        if ( rlnP(i) .ne. 1.d0 ) rtmp = rtmp + exp(rlnP(i))
        if ( lnP(i) .ne. 1.d0 ) rtmp2 = rtmp2 + exp(lnP(i))
      enddo
      rtmp = rtmp / dble(nBIN) * 100.d0
      rtmp2 = rtmp2 / dble(nBIN) * 100.d0
      write(6,'(2x,a,f8.3,a)')"+ Flatness              = ",rtmp," (%)"
      if ( accele ) write(6,'(2x,a,f8.3,a)')                           &
        "+ Flatness (accele)     = ",rtmp2," (%)"
      ! flatness for individual distribution
      rtmp = sum(flt(1:iset)) / dble(iset)
      write(6,'(2x,a,f8.3,a)')"+ Ave. Flat. for indiv. = ",rtmp," (%)"
      ia = minloc(flt(1:iset)) ; i = ia(1)
      write(6,'(2x,a,f8.3,a,i0,a)')"    min. = ",flt(i)," (%) (set ",  &
                                   i,")"
      ia = maxloc(flt(1:iset)) ; i = ia(1)
      write(6,'(2x,a,f8.3,a,i0,a)')"    max. = ",flt(i)," (%) (set ",  &
                                   i,")"
      write(6,*)
      ! E range for individual distribution
      rtmp = sum(Erng(1,1:iset)) / dble(iset)
      write(6,'(2x,a,f12.3)')"+ Ave. Erange for indiv. = ",rtmp
      ia = minloc(Erng(1,1:iset)) ; i = ia(1)
      write(6,'(2x,a,f12.3,a,i0,a)')"    min. = ",Erng(1,i)," (set ",i,&
                                   ")"
      ia = maxloc(Erng(1,1:iset)) ; i = ia(1)
      write(6,'(2x,a,f12.3,a,i0,a)')"    max. = ",Erng(1,i)," (set ",i,&
                                   ")"
      write(6,*)
      ! E range for individual distribution (sub)
      rtmp = sum(Erng(2,1:iset)) / dble(iset)
      write(6,'(2x,a,f12.3)')"+   Ave. Erange for indiv. = ",rtmp
      ia = minloc(Erng(2,1:iset)) ; i = ia(1)
      write(6,'(2x,a,f12.3,a,i0,a)')"      min. = ",Erng(2,i)," (set ",&
                                   i,")"
      ia = maxloc(Erng(2,1:iset)) ; i = ia(1)
      write(6,'(2x,a,f12.3,a,i0,a)')"      max. = ",Erng(2,i)," (set ",&
                                   i,")"
      write(6,*)

!*****************************
!     Output distribution

      open(unit=3,file=trim(PROJNM)//".dist",status="replace")
      if ( fwflg ) then
        do i = 1,nAB
          if ( rlnP(i) .ne. 1.d0 ) write(3,*)EE(i),rlnP(i),lnP(i)
        enddo
        write(3,*)
        write(3,*)minE,0.d0,0.d0
        write(3,*)minE,-3.d0,-3.d0
        write(3,*)
        write(3,*)maxE,0.d0,0.d0
        write(3,*)maxE,-3.d0,-3.d0
      else
        do i = 1,nAB
          if ( lnP(i) .ne. 1.d0 ) write(3,*)EE(i),lnP(i)
        enddo
        write(3,*)
        write(3,*)minE,0.d0
        write(3,*)minE,-3.d0
        write(3,*)
        write(3,*)maxE,0.d0
        write(3,*)maxE,-3.d0
      endif
      close(3)

!      ! Make individual distribution movies
!      open(unit=2,file=".gnu.plt",status="replace")
!      write(2,'(a)')'set title "lambda"'
!      write(2,'(a)')'set data style line'
!      write(2,'(a)')'set terminal png'
!      do i = 1,iset
!        write(tmp,'(a,i0)')trim(dirname)//'/'//trim(PROJNM)//"_",i
!        j = system("ls "//trim(tmp)//".dist >& /dev/null")
!        if ( j .ne. 0 ) cycle
!        write(2,'(a)')'set output "'//trim(tmp)//'.png"'
!        write(2,'(a)')'plot "'//trim(tmp)//'.dist"'
!      enddo
!      write(2,'(a)')'set terminal X11'
!      call systemqq("gnuplot .gnu.plt")
!      close(unit=2,status="delete")

!*****************************
!     Calc. Differential Density of States

      ! Input previous dlnN & Calc. current dlnN
      open(unit=2,file=trim(PREFIT),status="old")
      read(2,*)Nwindow_old,n
      allocate(COEold(0:n,Nwindow_old),low_old(Nwindow_old),         &
                 high_old(Nwindow_old))
      do i = 1,Nwindow_old
        do j = 0,n
          read(2,*)COEold(j,i)
        enddo
        read(2,*)low_old(i),high_old(i)
      enddo
      read(2,*)rr,ss
      read(2,*)TminE,TmaxE
      close(2)

      forall ( i=1:nAB )                                               &
        Escl(i) = (EE(i)+0.5d0*EbinSZ-TminE)/(TmaxE-TminE)

      rtmp = 1.d0 / EbinSZ ; j = 1
      do i = 1,nAB-1
        rtmp3 = Escl(i)
        if ( rtmp3 .lt. 0.d0 ) then
          rtmp2 = rr
        elseif ( rtmp3 .gt. 1.d0 ) then
          rtmp2 = ss
        else
          do
            if ( j.ne.Nwindow_old .and. rtmp3.ge.high_old(j) ) then
              j = j + 1 ; cycle
            endif
            exit
          enddo
          dx = rtmp3 - low_old(j)
          rtmp2 = COEold(0,j)
          do k = 1,n
            rtmp2 = rtmp2 + COEold(k,j)*dx**k
          enddo
        endif
        if ( lnP(i).ne.1.d0 .and. lnP(i+1).ne.1.d0 ) then
          rdlnN(i) = (rlnP(i+1)-rlnP(i))*rtmp + rtmp2
          if ( accel2 .and. rtmp3.ge.0.d0 .and. rtmp3.le.1.d0 ) then
            dlnN(i) = accrt2*(lnP(i+1)-lnP(i))*rtmp + rtmp2
          else
            dlnN(i) = (lnP(i+1)-lnP(i))*rtmp + rtmp2
          endif
        else
          if ( rtmp3.ge.0.d0 .and. rtmp3.le.1.d0 ) then
            rdlnN(i) = rtmp2 ; dlnN(i) = rtmp2
          else
            rdlnN(i) = 0.d0 ; dlnN(i) = 0.d0
          endif
        endif
      enddo
      close(1)

!************************
!     Output dlnN Function (*.nf)

!      ! delete minority data
!      minlnP = -3.5d0
!      do i = 1,nAB
!        if ( lnP(i) .lt. minlnP ) dlnN(i) = 0.d0
!      enddo

      Escl2(1:nAB) = (EE(1:nAB)+0.5d0*EbinSZ-nxminE)/(nxmaxE-nxminE)

!      !! Search min value
!      rtmp = dble(maxval(P))*2.d0
!      do i = 1,nAB
!        if ( P(i)+P(i+1) .ne. 0 ) rtmp = min(rtmp,dble(P(i)+P(i+1)))
!      enddo
!      rtmp = sqrt(rtmp)

      j = 0 ; tEE(1:nAB) = 0.d0 ; TdlnN(1:nAB-1) = 0.d0
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
      if ( igndat ) then
        do i = iST,iEN
          if ( dlnN(i) .eq. 0.d0 ) cycle
          j = j + 1 ; tEE(j) = Escl2(i) ; TdlnN(j) = dlnN(i)
          wei(j) = sqrt(dble(P(i)+P(i+1)))
        enddo
      else
        do i = iST,iEN
          if ( dlnN(i).eq.0.d0 .and.                                   &
               (Escl2(i).lt.TminE .or. Escl2(i).gt.TmaxE) ) cycle
          j = j + 1 ; tEE(j) = Escl2(i) ; TdlnN(j) = dlnN(i)
          wei(j) = max(sqrt(dble(P(i)+P(i+1))),1.d0) ! pseudo count
        enddo
      endif

      open(unit=1,file=trim(PROJNM)//".nf",status="replace")
      write(6,'(2x,a)')"+ Fitting for dlnN"
      allocate(COE(0:NfitDIM,Nwindow))
      call lstsquare_final(Nwindow,j,NfitDIM,0.d0,1.d0,tEE(1:j),       &
        TdlnN(1:j),wei(1:j),COE,low_v_window,high_v_window,re,re2,     &
        .true.)
      write(1,*)Nwindow,NfitDIM
      do i = 1,Nwindow
        do j = 0,NfitDIM
          write(1,*)COE(j,i)
        enddo
        write(1,*)low_v_window(i),high_v_window(i)
      enddo
      lowd = COE(0,1)
      upd = COE(0,Nwindow)
      dx = 1.d0 - low_v_window(Nwindow)
      do i = 1,NfitDIM
        upd = upd + COE(i,Nwindow)*dx**i
      enddo
      write(1,*)lowd,upd
      write(1,*)nxminE,nxmaxE
      write(1,*)nxT
      close(1)
      write(6,*)

!*****************************
!     Output dlnN data for plot

      iST = 1 ; iEN = nAB-1
      open(unit=1,file=trim(PROJNM)//".dlnN",status="replace")
      rtmp3 = 0.5d0*EbinSZ

      if ( fwflg ) then
        j = 1
        do i = iST,iEN
          if ( dlnN(i) .eq. 0.d0 ) cycle
          rtmp2 = Escl2(i)
          do
            if ( j.ne.Nwindow .and. rtmp2.ge.high_v_window(j) ) then
              j = j + 1 ; cycle
            endif
            exit
          enddo
          dx = rtmp2 - low_v_window(j)
          rtmp = COE(0,j)
          do k = 1,NfitDIM
            rtmp = rtmp + COE(k,j)*dx**k
          enddo
          write(1,'(4(e15.8,x))')EE(i)+rtmp3,dlnN(i),rtmp,rdlnN(i)
        enddo
        write(1,*)
        write(1,*)nxminE,lowd,lowd,lowd
        write(1,*)nxminE,upd,upd,upd
        write(1,*)
        write(1,*)nxmaxE,lowd,lowd,lowd
        write(1,*)nxmaxE,upd,upd,upd

      else

        j = 1
        do i = iST,iEN
          if ( dlnN(i) .eq. 0.d0 ) cycle
          rtmp2 = Escl2(i)
          do
            if ( j.ne.Nwindow .and. rtmp2.ge.high_v_window(j) ) then
              j = j + 1 ; cycle
            endif
            exit
          enddo
          dx = rtmp2 - low_v_window(j)
          rtmp = COE(0,j)
          do k = 1,NfitDIM
            rtmp = rtmp + COE(k,j)*dx**k
          enddo
          write(1,'(3(e15.8,x))')EE(i)+rtmp3,dlnN(i),rtmp
        enddo
        write(1,*)
        write(1,*)nxminE,lowd,lowd
        write(1,*)nxminE,upd,upd
        write(1,*)
        write(1,*)nxmaxE,lowd,lowd
        write(1,*)nxmaxE,upd,upd
      endif
      close(1)
 
!*****************************

      if ( smpT.ne.-1.d0 ) then

        ! Calc. canonical distribution (lnPc)
        rtmp = 1.d0 / ( R*smpT ) ; ss = -999999.d0 ; rtmp3 = 0.d0
        j = 1
        do i = 1,nAB
          if ( rlnP(i) .ne. 1.d0 ) then
            rtmp2 = (EE(i)-nxminE)/(nxmaxE-nxminE)
            do
              if ( j.ne.Nwindow .and. rtmp2.ge.high_v_window(j) ) then
                j = j + 1 ; cycle
              endif
              exit
            enddo
            dx = rtmp2 - low_v_window(j)
            ! lnN
            rr = COE(0,j)
            do k = 1,NfitDIM
              rr = rr + COE(k,j)*dx**k
            enddo
            rtmp3 = rtmp3 + rr*EbinSZ
            lnPc(i) = rtmp3 - EE(i) * rtmp
            ss = max(ss,lnPc(i))
          endif
        enddo
        lnPc(1:nAB) = lnPc(1:nAB) - ss

        ! Output Probability file & canonical distribution (lnPc)
        do i = 1,nAB
          if ( rlnP(i).ne.1.d0 .and. lnPc(i).ge.-7.d0 ) then
            itmp = i ; exit
          endif
        enddo
        do i = nAB,1,-1
          if ( rlnP(i).ne.1.d0 .and. lnPc(i).ge.-7.d0 ) then
            itmp2 = i ; exit
          endif
        enddo
        ENEbinSZ = (EE(itmp2)-EE(itmp)+EbinSZ) / dble(nBIN)
        minENEc = EE(itmp) - ENEbinSZ * 0.5d0
        forall ( i=1:nBIN ) EENE(i) = minENEc + ENEbinSZ*(dble(i)-0.5d0)
        forall ( i=1:nBIN )                                            &
          EENEscl2(i) = (EENE(i)-nxminE) / (nxmaxE-nxminE)
        lnPc(1:nAB) = 0.d0 ; rtmp = 1.d0 / ( R*smpT ) ; rtmp3 = 0.d0
        j = 1
        do i = 1,nBIN
          rtmp2 = EENEscl2(i)
          do
            if ( j.ne.Nwindow .and. rtmp2.ge.high_v_window(j) ) then
              j = j + 1 ; cycle
            endif
            exit
          enddo
          dx = rtmp2 - low_v_window(j)
          ! lnN
          rr = COE(0,j)
          do k = 1,NfitDIM
            rr = rr + COE(k,j)*dx**k
          enddo
          rtmp3 = rtmp3 + rr*ENEbinSZ
          lnPc(i) = rtmp3 - EENE(i)*rtmp
        enddo
        ss = maxval(lnPc(1:nBIN))
        lnPc(1:nBIN) = lnPc(1:nBIN) - ss
          
        n = smpT ; write(tmp,'(i0)')n
        tmp = trim(PROJNM)//"_"//trim(tmp)

        !! Output canonical distribution (lnPc)
        open(unit=1,file=trim(tmp)//".cdist",status="replace")
        do i = 1,nBIN
          write(1,*)EENE(i),lnPc(i)
        enddo
        close(1)
        !! Output Probability file
        open(unit=1,file=trim(tmp)//".prob",status="replace")
        write(1,*)"MULT"
        write(1,*)nBIN            ! Number of bins
        write(1,*)ENEbinSZ        ! Bin size
        write(1,*)minENEc         ! Min value
        write(1,*)
        do i = 1,nBIN
          write(1,*)exp(lnPc(i))
        enddo
        close(1)
!        !! Output count file
!        open(unit=1,file=trim(tmp)//".counthead",status="replace")
!        write(1,'(a5)')METHOD
!        write(1,*)nBIN,ENEbinSZ,minENEc
!        close(1)
!        open(unit=1,file=trim(tmp)//".count",status="replace")
!        do i = 1,nBIN
!          write(1,'(x,i0,$)')rP(i)
!        enddo
!        write(1,*)
!        close(1)

!****************

        ! Output current temperature range
        !  ( if dlnN is correct, you can calc the temperature)
        open(unit=1,file=trim(PROJNM)//".et",status="replace")
        j = 1
        ! Temperature for minE
        rtmp = (minE-nxminE) / (nxmaxE-nxminE)
        do
          if ( j.ne.Nwindow .and. rtmp.ge.high_v_window(j) ) then
            j = j + 1 ; cycle
          endif
          exit
        enddo
        dx = rtmp - low_v_window(j)
        rr = COE(0,j)
        do k = 1,NfitDIM
          rr = rr + COE(k,j)*dx**k
        enddo
        minT = 1.d0 / (rr*R)
        write(1,*)minE,minT

        ! Temperature at each energy bin
        do i = nBOR+1,nBOR+nBIN
          rtmp = (EE(i)-nxminE) / (nxmaxE-nxminE)
          do
            if ( j.ne.Nwindow .and. rtmp.ge.high_v_window(j) ) then
              j = j + 1 ; cycle
            endif
            exit
          enddo
          dx = rtmp - low_v_window(j)
          rr = COE(0,j)
          do k = 1,NfitDIM
            rr = rr + COE(k,j)*dx**k
          enddo
          TT(i) = 1.d0 / (rr*R)
          write(1,*)EE(i),TT(i)
        enddo

        ! Temperature for maxE
        rtmp = (maxE-nxminE) / (nxmaxE-nxminE)
        do
          if ( j.ne.Nwindow .and. rtmp.ge.high_v_window(j) ) then
            j = j + 1 ; cycle
          endif
          exit
        enddo
        dx = rtmp - low_v_window(j)
        rr = COE(0,j)
        do k = 1,NfitDIM
          rr = rr + COE(k,j)*dx**k
        enddo
        maxT = 1.d0 / (rr*R)
        write(1,*)maxE,maxT
        close(1)

        write(6,'(2x,a)')"* Current temperature range : "
        write(6,'(8x,f8.3,a,f8.3,a)')minT," - ",maxT," (K)"

        ! SPecific Heat
        open(unit=2,file=trim(PROJNM)//".sph",status="replace")
        do i = nBOR+2,nBOR+nBIN
          rtmp = TT(i) - TT(i-1)
          rr = EbinSZ / rtmp
          ss = rtmp * 0.5d0
          write(2,*)ss,rr
        enddo
        close(2)

      endif

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
      end subroutine distrib_MULT 
