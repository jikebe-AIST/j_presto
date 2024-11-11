
      subroutine distrib_CANO

!*********************************************
!
!     DISTRIBution_CANONICAL
!
!*********************************************

      use COMIFN ; use COMVAL

      implicit none

      ! Energy BIN SiZe
        real(8):: EbinSZ
      ! distribution & lnP
        integer(4):: P(nAB)
        real(8):: lnP(nAB),lnPc(nAB)
      ! (Temporary) Defferential Density of state
        real(8):: dlnN(nAB),TdlnN(nAB)
      ! MINimum Energy for Count
        real(8):: minEc
      ! (Temporary) Energy of each bin (EE is for dlnN)
        real(8):: EE(nAB),Escl2(nAB),tEE(nAB)
      ! FILE NaMe
        character(130):: filnam
        integer(4):: iST,iEN
      ! COEfficient
        real(8),allocatable:: COE(:,:)
        real(8):: lowd,upd
      ! WEIghting factor
        real(8):: wei(nAB),Twei(nAB)
      ! MIN & MAX Temperature for sampling range
        real(8):: minT,maxT

      integer(4):: i,j,k,l,n,system,itmp,itmp2,ifirst,ilast,ifirst2,   &
                   ilast2,ic
      real(8):: E,rtmp,rtmp2,rtmp3,rtmp4,tt,minE2,maxE2,dx,ss,re,re2,  &
                low_v_window(Nwindow),high_v_window(Nwindow)

!*****************************
!     Initialization

      EbinSZ = (maxE - minE) / dble(nBIN)
      minEc = minE - EbinSZ*dble(nBOR)
      l = 0 ; n = 0
      forall ( i=1:nAB ) EE(i) = minEc + EbinSZ*(dble(i-1)+0.5d0)
      lnP(1:nAB) = 0.d0 ; dlnN(1:nAB) = 0.d0 ; wei(1:nAB) = 0.d0
      minE2 = 9999999.d0 ; maxE2 = -9999999.d0

!*****************************
!     Input Energy data

      itmp = nBOR+1 ; itmp2 = nBOR+nBIN ; ic = 0
      open(unit=1,file=trim(DATLST),status="old")
      open(unit=4,file=trim(PROJNM)//".bar",status="replace")
      do
        read(1,'(a)',end=801)filnam
        i = index(filnam,";")
        if ( i .ne. 0 ) filnam = filnam(1:i-1)
        if ( filnam .eq. " " ) cycle
        read(1,*)iST,iEN,tt
        if ( iST .le. 1 ) iST = 1
        i = system("ls "//trim(filnam)//" >& /dev/null")

        if ( i .eq. 0 ) then
          ic = ic + 1
          write(6,'(2x,a)')"+ Now input "//trim(filnam)

          ! Input energy file
          P(1:nAB) = 0
          open(unit=2,file=trim(filnam),status="old")
          do j = 1,iST-1
            read(2,*,end=800)
          enddo
          do j = iST,iEN
            read(2,*,end=800)E ; l = l + 1
            k = int((E - minEc) / EbinSZ) + 1
            minE2 = min(minE2,E) ; maxE2 = max(maxE2,E)
            if ( k.ge.1 .and. k.le.nAB ) then
              P(k) = P(k) + 1
              if ( k.ge.itmp .and. k.le.itmp2 ) n = n + 1
            endif
          enddo
800       close(2)

          ! Make each dist. data
          rtmp = maxval(P) ; rtmp = log(rtmp)
          where ( P .eq. 0 )
            lnP = 1.d0
          elsewhere
            lnP = log(dble(P)) - rtmp
          endwhere

          ! Calc. each differential density of states
          rtmp = 1.d0 / EbinSZ ; rtmp2 = 1.d0/(R*tt)
          do k = 1,nAB-1
            if ( lnP(k).ne.1.d0 .and. lnP(k+1).ne.1.d0 ) then
              rtmp3 = sum(P(k:k+1))
              dlnN(k) = dlnN(k) + rtmp3*((lnP(k+1)-lnP(k))*rtmp + rtmp2)
              wei(k) = wei(k) + rtmp3
            endif
          enddo

          ! Make bar graph
          do i = 1,nAB
            if ( lnP(i) .ne. 1.d0 ) then
              ifirst = i ; exit
            endif
          enddo
          do i = ifirst,nAB
            if ( lnP(i) .ge. -1.5d0 ) then
              ifirst2 = i ; exit
            endif
          enddo
          do i = nAB,1,-1
            if ( lnP(i) .ne. 1.d0 ) then
              ilast = i ; exit
            endif
          enddo
          do i = ilast,1,-1
            if ( lnP(i) .ge. -1.5d0 ) then
              ilast2 = i ; exit
            endif
          enddo
          write(4,'(e12.5,x,i0,x,e12.5)')EE(ifirst),ic,EE(ifirst2)
          write(4,'(e12.5,x,i0,x,e12.5)')EE(ilast),ic,EE(ilast2)
          write(4,*)

          write(6,'(8x,i0,a3,i0)')iST," - ",j-1
        endif

      enddo
801   close(1) ; close(4) ; write(6,*)
      write(6,'(2x,a,i0,a,i0,a,f8.3,a)')"+ ",n," / ",l," (",           &
           dble(n)/dble(l)*100.d0," (%)) data in the range are counted"
      write(6,*)
      write(6,'(2x,a)')"* Obtained min & max E"
      write(6,'(8x,f15.3,a,f15.3)')minE2," - ",maxE2
      write(6,*)

!*****************************
!     Output dlnN Function (*.nf)

      ! Normalization
      forall ( i=1:nAB-1, wei(i).ne.0.d0 ) dlnN(i) = dlnN(i) / wei(i)
      Escl2(1:nAB) = (EE(1:nAB)+0.5d0*EbinSZ-nxminE)/(nxmaxE-nxminE)

      j = 0 ; tEE(1:nAB) = 0.d0 ; TdlnN(1:nAB) = 0.d0 ; Twei(1:nAB) = 0.d0
      do i = 1,nAB
        if ( wei(i) .ne. 0.d0 ) then
          j = j + 1 ; tEE(j) = Escl2(i) ; TdlnN(j) = dlnN(i)
          Twei(j) = sqrt(wei(i))
!          wei(j) = sqrt(dble(P(i)))
        endif
      enddo

      write(6,'(2x,a)')"+ Fitting for dlnN"
      allocate(COE(0:NfitDIM,Nwindow))
      call lstsquare_final(Nwindow,j,NfitDIM,0.d0,1.d0,tEE(1:j),       &
        TdlnN(1:j),Twei(1:j),COE,low_v_window,high_v_window,re,re2,    &
        .true.)

      open(unit=1,file=trim(PROJNM)//".nf",status="replace")
      write(1,*)Nwindow,NfitDIM
      do i = 1,Nwindow
        do j = 0,NfitDIM
          write(1,*)COE(j,i)
        enddo
        write(1,*)low_v_window(i),high_v_window(i)
      enddo
      lowd = COE(0,1) ; upd = COE(0,Nwindow)
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

      open(unit=1,file=trim(PROJNM)//".dlnN",status="replace")
      rtmp3 = 0.5d0*EbinSZ ; j = 1
      do i = 1,nAB-1
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
      close(1)
 
!*****************************
!     Output distribution

      open(unit=3,file=trim(PROJNM)//".cdist",status="replace")
      rtmp2 = 1.d0 / (R*T) ; rtmp4 = 0.d0 ; j = 1 ; ss = -9999999.d0
      do i = 1,nAB
        if ( EE(i) .lt. nxminE ) cycle
        if ( EE(i) .gt. nxmaxE ) exit 
        rtmp3 = (EE(i)-nxminE)/(nxmaxE-nxminE)
        do
          if ( j.ne.Nwindow .and. rtmp3.ge.high_v_window(j) ) then
            j = j + 1 ; cycle
          endif
          exit
        enddo
        dx = rtmp3 - low_v_window(j)
        ! lnN
        rtmp = COE(0,j)
        do k = 1,NfitDIM
          rtmp = rtmp + COE(k,j)*dx**k
        enddo
        rtmp4 = rtmp4 + rtmp*EbinSZ
        lnPc(i) = rtmp4 - EE(i) * rtmp2
        ss = max(ss,lnPc(i))
      enddo
      lnPc(1:nAB) = lnPc(1:nAB) - ss

!      do i = 1,nAB
!        rtmp3 = (EE(i)-nxminE)/(nxmaxE-nxminE)
!        if ( rtmp3 .le. 0.d0 ) then
!          rtmp = rtmp3 * lowd
!        elseif ( EE(i) .ge. 1.d0 ) then
!          rtmp = COE(1)
!          do j = 2,Ndim+1
!            rtmp = rtmp + COE(j) / dble(j)
!          enddo
!          rtmp = rtmp + (rtmp3-1.d0) * upd
!        else
!          rtmp = COE(1)*rtmp3
!          do j = 2,Ndim+1
!            rtmp = rtmp + COE(j)*rtmp3**(j) / dble(j)
!          enddo
!        endif
!        lnPc(i) = rtmp - EE(i)*rtmp2
!      enddo

      do i = 1,nAB
        if ( EE(i) .lt. nxminE ) cycle
        if ( EE(i) .gt. nxmaxE ) exit
        write(3,*)EE(i),lnPc(i)
      enddo
      write(3,*)
      write(3,*)nxminE,0.d0
      write(3,*)nxminE,-3.d0
      write(3,*)
      write(3,*)nxmaxE,0.d0
      write(3,*)nxmaxE,-3.d0
      close(3)

!************************

      return
      end subroutine distrib_CANO
