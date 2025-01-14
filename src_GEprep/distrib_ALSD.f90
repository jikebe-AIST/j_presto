
      subroutine distrib_ALSD

!*********************************************
!
!     DISTRIBution_ADAPTIVE_LAMBDA_DYNAMICS
!
!*********************************************

      use COMIFN ; use COMVAL

      implicit none

      ! Energy BIN SiZe
        real(8):: EbinSZ
      ! distribution & lnP
        integer(4):: ttP(nAB),TTPdist(3,nAB)
        real(8):: trP(nAB),tP(nAB),rP(nAB),P(nAB),rlnP(nAB),lnP(nAB),  &
                  tlnP(nAB)
      ! (Temporary) Defferential Density of state
        real(8):: rdlnN(nAB-1),dlnN(nAB-1),TdlnN(nAB-1)
      ! MINimum Energy for Count
        real(8):: minEc
      ! (Temporary) Energy of each bin (EE is for dlnN)
        real(8):: EE(nAB),tEE(nAB),ELMD(nlmd),EEAA(nEAA),EEAB(nEAB)
      ! FILE NaMe
        character(130):: filnam
        integer(4):: iST,iEN
      ! COEfficient
        integer(4):: Nwindow_old,ndim_old
        real(8),allocatable:: COE(:,:),COEold(:,:),low_old(:),         &
                              high_old(:)
        real(8):: lowd,upd,lowd_old,upd_old,TminE,TmaxE,kijyun
      ! Flatness calc.
        real(8),allocatable:: flt(:),lmdrng(:,:)
      ! WEIghting factor
        real(8),allocatable:: wei(:),weit(:,:),tw(:,:,:)
        logical(4),allocatable:: twF(:,:,:)
      ! Canonical distribution (lnPc)
        integer(4),allocatable:: Pc(:,:,:)
        real(8):: lnPc(nlmd)
      ! MIN & MAX Temperature for sampling range
        real(8):: minT,maxT
      ! Temporary Temperature
        real(8):: TT(nBOR+1:nBOR+1+nBIN)
      ! average EAA, EAB
        real(8),allocatable:: aveEAA(:),aveEAB(:),TaveEAA(:),          &
          TaveEAB(:),TcEAAEAB(:),acclabrat(:)
        integer(4),allocatable:: cEAAEAB(:)

      integer(4):: i,j,k,k1,k2,k3,k4,n2,itmp,itmp2,ifile,iset,ier, &
                   ifirst,ilast,ifirst2,ilast2,ii,ia(1)
      integer(8):: n,l
      real(8):: fw,E,EAA,EAB,ENE,rtmp,rtmp2,rtmp3,rtmp4,minlnP,dx,dl,  &
                re,re2
      character(130):: tmp,dirname
      logical(4):: onend,fwflg,ex
      real(8):: minLMD2,maxLMD2,minEAA2,maxEAA2,minEAB2,maxEAB2,avecnt
      real(8):: low_v_window(Nwindow),high_v_window(Nwindow),lnN,      &
                lower,upper,lower_old,upper_old
      integer(4):: iavecnt,maxcnt,iminEAA,imaxEAA,iminEAB,imaxEAB,iEAA,&
                   iEAB
      real(8):: lw(1),hw(1),COEEAAEAB(0:2,1)

!*****************************
!     Initialization

      ! Heap memory
9999  allocate(wei(nAB),weit(nEAA,nEAB),tw(nEAA,nEAB,nlmd))
      allocate(twF(nEAA,nEAB,nlmd),Pc(nlmd,nEAA,nEAB))
      allocate(aveEAA(nAB),aveEAB(nAB),cEAAEAB(nAB))
      allocate(TaveEAA(nAB),TaveEAB(nAB),TcEAAEAB(nAB),acclabrat(nAB))

      ifile = 0 ; fwflg = .false.
      if ( accele ) fwflg = .true.
      EbinSZ = (maxE - minE) / dble(nBIN)
      minEc = minE - EbinSZ*dble(nBOR)
      rP(1:nAB) = 0.d0 ; P(1:nAB) = 0.d0
      l = 0 ; n = 0
      forall ( i=1:nAB ) EE(i) = minEc + EbinSZ*(dble(i-1)+0.5d0)
      aveEAA(:) = 0.d0 ; aveEAB(:) = 0.d0 ; cEAAEAB(:) = 0

      if ( reweight_flg ) then
        Pc(1:nlmd,1:nEAA,1:nEAB) = 0 ; lnPc(1:nlmd) = 0.d0
        forall ( i=1:nlmd ) ELMD(i) = minlmd + slmd*(dble(i-1)+0.5d0)
        forall ( i=1:nEAA ) EEAA(i) = minEAA + sEAA*(dble(i-1)+0.5d0)
        forall ( i=1:nEAB ) EEAB(i) = minEAB + sEAB*(dble(i-1)+0.5d0)
        minLMD2 = 9999999.d0 ; maxLMD2 = -99999999.d0
        minEAA2 = 9999999.d0 ; maxEAA2 = -99999999.d0
        minEAB2 = 9999999.d0 ; maxEAB2 = -99999999.d0
        weit(1:nEAA,1:nEAB) = 0.d0
        tw(1:nEAA,1:nEAB,1:nlmd) = 0.d0
        twF(1:nEAA,1:nEAB,1:nlmd) = .false.
      endif

      ! Check number of input files
      call system("wc -l "//trim(DATLST)//" > .aho.temp")
      open(unit=1,file=".aho.temp",status="old")
      read(1,*)i
      close(unit=1,status="delete")
      i = i / 2 ; allocate(flt(i),lmdrng(2,i))

      ! Input previous dlnP
      open(unit=1,file=trim(PREFIT),status="old")
      read(1,*)Nwindow_old,ndim_old
      allocate(COEold(0:ndim_old,Nwindow_old),low_old(Nwindow_old),    &
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

!*****************************
!     Input Lambda data

      open(unit=1,file=trim(DATLST),status="old")
      open(unit=4,file=trim(PROJNM)//".bar",status="replace")
      itmp = nBOR+1 ; itmp2 = nBOR+nBIN
      trP(1:nAB) = 0.d0 ; tP(1:nAB) = 0.d0 ; iset = 0 ; TTPdist(:,:) = 0

      if ( reweight_flg ) then
        do
          call infil(1,filnam,iST,iEN,fw,onend,ier)
          if ( ier .ne. 0 ) goto 801
          inquire(file=trim(filnam),exist=ex)

          if ( ex ) then
            write(6,'(2x,a)')"+ Now input "//trim(filnam)
            if ( fw .ne. 1.d0 ) fwflg = .true.

            ! Input energy file
            ifile = ifile + 1 ; ttP(1:nAB) = 0
            open(unit=2,file=trim(filnam),status="old")
            do j = 1,iST-1
              read(2,*,end=800)
            enddo
            do j = iST,iEN
              read(2,*,end=800)E,EAA,EAB ; l = l + 1
              k = int((E - minEc) / EbinSZ + 1)
              if ( k.ge.1 .and. k.le.nAB ) then
                ttP(k) = ttP(k) + 1
                aveEAA(k) = aveEAA(k) + EAA
                aveEAB(k) = aveEAB(k) + EAB
                cEAAEAB(k) = cEAAEAB(k) + 1
                if ( k.ge.itmp .and. k.le.itmp2 ) n = n + 1
              endif

              k1 = int((E - minlmd) / slmd + 1)
              minLMD2 = min(minLMD2,E) ; maxLMD2 = max(maxLMD2,E)
              if ( k1.ge.1 .and. k1.le.nlmd ) then
                k2 = int((EAA - minEAA) / sEAA + 1)
                k3 = int((EAB - minEAB) / sEAB + 1)
                minEAA2 = min(minEAA2,EAA) ; maxEAA2 = max(maxEAA2,EAA)
                minEAB2 = min(minEAB2,EAB) ; maxEAB2 = max(maxEAB2,EAB)
                if ( k2.ge.1 .and. k2.le.nEAA .and.                    &
                     k3.ge.1 .and. k3.le.nEAB )                        &
                   Pc(k1,k2,k3) = Pc(k1,k2,k3) + 1
              endif
            enddo
800         close(2)

            trP(1:nAB) = trP(1:nAB) + dble(ttP(1:nAB))
            if ( fw .eq. 1.d0 ) then
              write(6,'(8x,i0,a3,i0)')iST," - ",j-1
              tP(1:nAB) = tP(1:nAB) + dble(ttP(1:nAB))
            else
              write(6,'(8x,i0,a3,i0,a5,f0.0,a1)')iST," - ",j-1," ( x ",&
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
                tlnP = log(trP) - rtmp
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
                    if ( tlnP(k) .ge. -0.5d0 )                         &
                      TTPdist(3,k) = TTPdist(3,k) + 1
                  endif
                endif
              enddo

              if ( accele ) tP(1:nAB) = tP(1:nAB)/flt(iset)**accrat
              P(1:nAB) = P(1:nAB) + tP(1:nAB)

!              write(tmp,'(a,a1,i0,a5)')trim(dirname)//"/"//            &
!                trim(PROJNM),"_",iset,".dist"
!              open(unit=3,file=trim(tmp),status="replace")
!              do i = 1,nAB
!                if ( tlnP(i) .ne. 1.d0 ) write(3,*)EE(i),tlnP(i)
!              enddo
!              write(3,*)
!              write(3,*)minE,0.d0
!              write(3,*)minE,-3.d0
!              write(3,*)
!              write(3,*)maxE,0.d0
!              write(3,*)maxE,-3.d0
!              close(3)

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
              lmdrng(1,iset) = (EE(ilast)-EE(ifirst)) + EbinSZ
              lmdrng(2,iset) = (EE(ilast2)-EE(ifirst2)) + EbinSZ
              write(6,'(2x,a,f8.3,a,f8.3,a)')"+ lrange = ",            &
                lmdrng(1,iset)," (",lmdrng(2,iset),")"
              write(6,*)
            endif

            trP(1:nAB) = 0.d0 ; tP(1:nAB) = 0.d0
          endif
        enddo

      else
        do
          call infil(1,filnam,iST,iEN,fw,onend,ier)
          if ( ier .ne. 0 ) goto 801
          inquire(file=trim(filnam),exist=ex)

          if ( ex ) then
            write(6,'(2x,a)')"+ Now input "//trim(filnam)
            if ( fw .ne. 1.d0 ) fwflg = .true.

            ! Input energy file
            ifile = ifile + 1 ; ttP(1:nAB) = 0
            open(unit=2,file=trim(filnam),status="old")
            do j = 1,iST-1
              read(2,*,end=803)
            enddo
            do j = iST,iEN
              read(2,*,end=803)E,EAA,EAB ; l = l + 1
              k = int((E - minEc) / EbinSZ + 1)

              if ( k.ge.1 .and. k.le.nAB ) then
                ttP(k) = ttP(k) + 1
                aveEAA(k) = aveEAA(k) + EAA
                aveEAB(k) = aveEAB(k) + EAB
                cEAAEAB(k) = cEAAEAB(k) + 1
                if ( k.ge.itmp .and. k.le.itmp2 ) n = n + 1
              endif
            enddo
803         close(2)

            trP(1:nAB) = trP(1:nAB) + dble(ttP(1:nAB))
            if ( fw .eq. 1.d0 ) then
              write(6,'(8x,i0,a3,i0)')iST," - ",j-1
              tP(1:nAB) = tP(1:nAB) + dble(ttP(1:nAB))
            else
              write(6,'(8x,i0,a3,i0,a5,f0.0,a1)')iST," - ",j-1,        &
                " ( x ",fw,")"
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
                tlnP = log(trP) - rtmp
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
                    if ( tlnP(k) .ge. -0.5d0 )                         &
                      TTPdist(3,k) = TTPdist(3,k) + 1
                  endif
                endif
              enddo

              if ( accele ) tP(1:nAB) = tP(1:nAB)/flt(iset)**accrat
              P(1:nAB) = P(1:nAB) + tP(1:nAB)

!              write(tmp,'(a,a1,i0,a5)')trim(dirname)//"/"//            &
!                trim(PROJNM),"_",iset,".dist"
!              open(unit=3,file=trim(tmp),status="replace")
!              do i = 1,nAB
!                if ( tlnP(i) .ne. 1.d0 ) write(3,*)EE(i),tlnP(i)
!              enddo
!              write(3,*)
!              write(3,*)minE,0.d0
!              write(3,*)minE,-3.d0
!              write(3,*)
!              write(3,*)maxE,0.d0
!              write(3,*)maxE,-3.d0
!              close(3)

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
              lmdrng(1,iset) = (EE(ilast)-EE(ifirst)) + EbinSZ
              lmdrng(2,iset) = (EE(ilast2)-EE(ifirst2)) + EbinSZ
              write(6,'(2x,a,f8.3,a,f8.3,a)')"+ lrange = ",            &
                lmdrng(1,iset)," (",lmdrng(2,iset),")"
              write(6,*)
            endif

            trP(1:nAB) = 0.d0 ; tP(1:nAB) = 0.d0
          endif
        enddo
      endif
801   close(1) ; write(6,*)
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
          deallocate(wei,weit,tw,twF,Pc,aveEAA,aveEAB,cEAAEAB,TaveEAA, &
            TaveEAB,TcEAAEAB,acclabrat,flt,lmdrng,COEold,low_old,      &
            high_old)
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
          goto 9999
        endif
      endif

!***************************
!     Output average EAB for each lambda

!      write(6,'(2x,a)')"+ Fitting for lambda-EAB"
      j = 0 ; TaveEAA(:) = 0.d0 ; TaveEAB(:) = 0.d0 ; tEE(1:nAB) = 0.d0
      do i = 1,nAB
        if ( cEAAEAB(i) .eq. 0 ) cycle
        j = j + 1
        rtmp = 1.d0/cEAAEAB(i)
        tEE(j) = EE(i) + 0.5d0*EbinSZ
        aveEAA(i) = aveEAA(i)*rtmp
        TaveEAA(j) = aveEAA(i)
        aveEAB(i) = aveEAB(i)*rtmp
        TaveEAB(j) = aveEAB(i)
        TcEAAEAB(j) = cEAAEAB(i)
      enddo
      call lstsquare_final(1,j,2,nxminE,nxmaxE,tEE(1:j),TaveEAB(1:j),  &
        TcEAAEAB(1:j),COEEAAEAB,lw,hw,re,re2,.false.)

      open(unit=1,file=trim(PROJNM)//".LAB",status="replace")
      rtmp2 = re2/re*accthl ; rtmp3 = maxval(TcEAAEAB(1:j))
      j = 0 ; acclabrat(:) = 0.d0
      if ( acclab ) then
        do i = 1,nAB
          if ( cEAAEAB(i) .eq. 0 ) cycle
          j = j + 1
          dx = tEE(j) - lw(1)
          rtmp = COEEAAEAB(0,1) + COEEAAEAB(1,1)*dx +                  &
                 COEEAAEAB(2,1)*dx**2
          acclabrat(i) = sqrt((TaveEAB(j)-rtmp)**2) / re
          write(1,'(4(f,x,$))')tEE(j),TaveEAA(j),TaveEAB(j),rtmp
          if ( TcEAAEAB(j).ge.rtmp3*0.01d0 .and.                       &
               acclabrat(i).ge.1.d0+rtmp2*sqrt(rtmp3/cEAAEAB(j)) ) then
            P(i) = P(i) * acclabrat(i) * accrtl
            write(1,'(f,$)')TaveEAB(j)
          endif
          write(1,*)
        enddo
      else
        do i = 1,nAB
          if ( cEAAEAB(i) .eq. 0 ) cycle
          j = j + 1
          dx = tEE(j) - lw(1)
          rtmp = COEEAAEAB(0,1) + COEEAAEAB(1,1)*dx +                  &
                 COEEAAEAB(2,1)*dx**2
          write(1,'(4(f,x))')tEE(j),TaveEAA(j),TaveEAB(j),rtmp
        enddo
      endif
      close(1)

!****************************
!     Make distribution data

      !! real P
      rtmp = maxval(rP) ; rtmp = log(rtmp)
      where ( rP .eq. 0.d0 )
        rlnP = 1.d0
      elsewhere
        rlnP = log(rP) - rtmp
      endwhere
      !! P using the weigting factor, fw
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
      ! lambda range for individual distribution
      rtmp = sum(lmdrng(1,1:iset)) / dble(iset)
      write(6,'(2x,a,f8.3)')"+ Ave. lrange for indiv. = ",rtmp
      ia = minloc(lmdrng(1,1:iset)) ; i = ia(1)
      write(6,'(2x,a,f8.3,a,i0,a)')"    min. = ",lmdrng(1,i)," (set ", &
                                   i,")"
      ia = maxloc(lmdrng(1,1:iset)) ; i = ia(1)
      write(6,'(2x,a,f8.3,a,i0,a)')"    max. = ",lmdrng(1,i)," (set ", &
                                   i,")"
      write(6,*)
      ! lambda range for individual distribution (sub)
      rtmp = sum(lmdrng(2,1:iset)) / dble(iset)
      write(6,'(2x,a,f8.3)')"+   Ave. lrange for indiv. = ",rtmp
      ia = minloc(lmdrng(2,1:iset)) ; i = ia(1)
      write(6,'(2x,a,f8.3,a,i0,a)')"      min. = ",lmdrng(2,i),        &
                                   " (set ",i,")"
      ia = maxloc(lmdrng(2,1:iset)) ; i = ia(1)
      write(6,'(2x,a,f8.3,a,i0,a)')"      max. = ",lmdrng(2,i),        &
                                   " (set ",i,")"
      write(6,*)

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

        maxcnt = 0 ; avecnt = sum(Pc) ; iavecnt = 0 ; rtmp4 = 1.d0/(R*T)
        ii = 1 ; lnN = 0.d0

        if ( METHOD .eq. "ALSD" ) then

          do i = 1,nlmd
            rtmp3 = ELMD(i)
            do
              if ( ii.ne.Nwindow_old .and. rtmp3.ge.high_old(ii) ) then
                ii = ii + 1 ; cycle
              endif
              exit
            enddo
            dx = rtmp3 - low_old(ii)
            rtmp2 = COEold(0,ii)
            do l = 1,ndim_old
              rtmp2 = rtmp2 + COEold(l,ii)*dx**l
            enddo
            if ( rtmp3 .lt. TminE ) then
              dl = TminE - rtmp3
              rtmp2 = rtmp2 + lower_old*dl
            elseif ( rtmp3 .gt. TmaxE ) then
              dl = rtmp3 - TmaxE
              rtmp2 = rtmp2 + upper_old*dl
            endif

            lnN = lnN + rtmp2*slmd
            do k = 1,nEAB
            do j = 1,nEAA
              itmp = Pc(i,j,k)
              if ( itmp .ne. 0 ) then
                iavecnt = iavecnt + 1 ; maxcnt = max(maxcnt,itmp)
              endif
              rtmp = ((rtmp3*rtmp3-lmd_rewei2)*EEAA(j)                 &
                    + (rtmp3-lmd_rewei)*EEAB(k))*rtmp4
              tw(j,k,i) = -(rtmp+lnN) ; twF(j,k,i) = .true.
            enddo
            enddo
          enddo

        elseif ( METHOD .eq. "ALSD2" ) then

          do i = 1,nlmd
            rtmp3 = ELMD(i)
            do
              if ( ii.ne.Nwindow_old .and. rtmp3.ge.high_old(ii) ) then
                ii = ii + 1 ; cycle
              endif
              exit
            enddo
            dx = rtmp3 - low_old(ii)
            rtmp2 = COEold(0,ii)
            do l = 1,ndim_old
              rtmp2 = rtmp2 + COEold(l,ii)*dx**l
            enddo
            if ( rtmp3 .lt. TminE ) then
              dl = TminE - rtmp3
              rtmp2 = rtmp2 + lower_old*dl
            elseif ( rtmp3 .gt. TmaxE ) then
              dl = rtmp3 - TmaxE
              rtmp2 = rtmp2 + upper_old*dl
            endif
            lnN = lnN + rtmp2*slmd
            do k = 1,nEAB
            do j = 1,nEAA
              itmp = Pc(i,j,k)
              if ( itmp .ne. 0 ) then
                iavecnt = iavecnt + 1 ; maxcnt = max(maxcnt,itmp)
              endif
              rtmp = ((rtmp3-lmd_rewei2)*EEAA(j)                       &
                    + (sqrt(rtmp3)-lmd_rewei)*EEAB(k))*rtmp4
              tw(j,k,i) = -(rtmp+lnN) ; twF(j,k,i) = .true.
            enddo
            enddo
          enddo

        endif

        rtmp = maxval(tw(1:nEAA,1:nEAB,1:nlmd),mask=twF)
        do i = 1,nlmd
          rtmp2 = dble(sum(Pc(i,1:nEAA,1:nEAB)))
          forall ( j=1:nEAA, k=1:nEAB, twF(j,k,i) )                    &
            weit(j,k) = weit(j,k) + exp(tw(j,k,i)-rtmp)*rtmp2
        enddo

        ! lambda dist.
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
        kijyun = maxval(lnPc) ; kijyun = log(kijyun)
        where ( lnPc .eq. 0.d0 )
          lnPc = 1.d0
        elsewhere
          lnPc = log(lnPc) - kijyun
        endwhere

        write(6,'(2x,a)')"* N of conf. information for canonical dist."
        write(6,'(8x,a,i0)')  "+                Max : ",maxcnt
        write(6,'(8x,a,f10.3)')"+                Ave : ",              &
                               avecnt/dble(iavecnt)
        write(6,'(8x,a,i0)')  "+ N of non-zero bins : ",iavecnt
      endif

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

      if ( reweight_flg ) then
        ! Lambda
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
        ! Output Probability file
        forall ( i=1:nEAA, j=1:nEAB, weit(i,j).ne.0.d0 )               &
          weit(i,j) = 1.d0 / weit(i,j)
        rtmp = maxval(weit(1:nEAA,1:nEAB))
        weit(1:nEAA,1:nEAB) = weit(1:nEAA,1:nEAB) / rtmp
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

!        open(unit=3,file=trim(PROJNM)//".counthead",status="replace")
!        write(3,'(a5)')METHOD
!        write(3,*)nlmd,nEAA,nEAB
!        write(3,*)slmd,sEAA,sEAB
!        write(3,*)minLMD,minEAA,minEAB
!        close(3)
!        open(unit=3,file=trim(PROJNM)//".count",status="replace")
!        do k = 1,nEAB
!        do j = 1,nEAA
!          do i = 1,nlmd
!            write(3,'(x,i0,$)')sum(Pc(i,j,k,:))
!          enddo
!          write(3,*)
!        enddo
!        enddo
!        write(3,*)
!        close(3)

      endif
 
!*****************************
!     Calc. Differential Density of States

      ! Calc. current dlnN
      rtmp = 1.d0 / EbinSZ ; ii = 1
      do i = 1,nAB-1
        rtmp3 = EE(i) + 0.5d0*EbinSZ
        do
          if ( ii.ne.Nwindow_old .and. rtmp3.ge.high_old(ii) ) then
            ii = ii + 1 ; cycle
          endif
          exit
        enddo
        dx = rtmp3 - low_old(ii)
        rtmp2 = COEold(0,ii)
        do j = 1,ndim_old
          rtmp2 = rtmp2 + COEold(j,ii)*dx**j
        enddo
        if ( rtmp3 .lt. TminE ) then
          dl = TminE - rtmp3
          rtmp2 = rtmp2 + lower_old*dl
        elseif ( rtmp3 .gt. TmaxE ) then
          dl = rtmp3 - TmaxE
          rtmp2 = rtmp2 + upper_old*dl
        endif

        if ( lnP(i).ne.1.d0 .and. lnP(i+1).ne.1.d0 ) then
          rdlnN(i) = (rlnP(i+1)-rlnP(i))*rtmp + rtmp2
          if ( accel2 .and. rtmp3.ge.TminE .and. rtmp3.le.TmaxE ) then
            dlnN(i) = accrt2*(lnP(i+1)-lnP(i))*rtmp + rtmp2
          else
            dlnN(i) = (lnP(i+1)-lnP(i))*rtmp + rtmp2
          endif
        else
          if ( rtmp3.ge.TminE .and. rtmp3.le.TmaxE ) then
            rdlnN(i) = rtmp2 ; dlnN(i) = rtmp2
          else
            rdlnN(i) = 0.d0 ; dlnN(i) = 0.d0
          endif
        endif
      enddo

!*********************************
!     Output dlnP Function (*.nf)

      j = 0 ; tEE(1:nAB) = 0.d0 ; TdlnN(1:nAB-1) = 0.d0
      allocate(COE(0:NfitDIM,Nwindow))
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
          j = j + 1 ; tEE(j) = EE(i) ; TdlnN(j) = dlnN(i)
          wei(j) = sqrt(P(i)+P(i+1))
        enddo
      else      
        do i = iST,iEN
          if ( dlnN(i).eq.0.d0 .and.                                   &
               (EE(i).lt.TminE .or. EE(i).gt.TmaxE) ) cycle
          j = j + 1 ; tEE(j) = EE(i) ; TdlnN(j) = dlnN(i)
          wei(j) = max(sqrt(P(i)+P(i+1)),1.d0) ! pseudo count
        enddo
      endif
      tEE(1:j) = tEE(1:j) + 0.5d0*EbinSZ

      write(6,'(2x,a)')"+ Fitting for dlnP"
      call lstsquare_final(Nwindow,j,NfitDIM,nxminE,nxmaxE,tEE(1:j),   &
        TdlnN(1:j),wei(1:j),COE,low_v_window,high_v_window,re,re2,     &
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
      upd = COE(0,Nwindow) ; dx = nxmaxE - low_v_window(Nwindow)
      do i = 1,NfitDIM
        upd = upd + COE(i,Nwindow)*dx**i
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

      iST = 1 ; iEN = nAB-1
      open(unit=1,file=trim(PROJNM)//".dlnP",status="replace")
      ii = 1
      if ( fwflg ) then
        do i = iST,iEN
          if ( dlnN(i) .eq. 0.d0 ) cycle
          rtmp2 = EE(i) + 0.5d0*EbinSZ
          do
            if ( ii.ne.Nwindow .and. rtmp2.ge.high_v_window(ii) ) then
              ii = ii + 1 ; cycle
            endif
            exit
          enddo
          dx = rtmp2 - low_v_window(ii) ; rtmp = COE(0,ii)
          do j = 1,NfitDIM
            rtmp = rtmp + COE(j,ii)*dx**j
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
          write(1,'(5(e15.8,x))')rtmp2,dlnN(i),rtmp,rtmp3,rdlnN(i)
        enddo
        write(1,*)
        write(1,*)nxminE,lowd,lowd,lowd
        write(1,*)nxminE,upd,upd,upd
        write(1,*)
        write(1,*)nxmaxE,lowd,lowd,lowd
        write(1,*)nxmaxE,upd,upd,upd
      else
        do i = iST,iEN
          if ( dlnN(i) .eq. 0.d0 ) cycle
          rtmp2 = EE(i) + 0.5d0*EbinSZ
          do
            if ( ii.ne.Nwindow .and. rtmp2.ge.high_v_window(ii) ) then
              ii = ii + 1 ; cycle
            endif
            exit
          enddo
          dx = rtmp2 - low_v_window(ii) ; rtmp = COE(0,ii)
          do j = 1,NfitDIM
            rtmp = rtmp + COE(j,ii)*dx**j
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
          write(1,'(4(e15.8,x))')rtmp2,dlnN(i),rtmp,rtmp3
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
!     Output TTPdist

      open(unit=1,file=trim(PROJNM)//".TTPdist",status="replace")
      do i = 1,nAB
        rtmp2 = EE(i) + 0.5d0*EbinSZ
        if ( TTPdist(1,i).ne.0 ) write(1,*)rtmp2,TTPdist(1:3,i)
      enddo
      close(1)

!*****************************

      return
      end subroutine distrib_ALSD
