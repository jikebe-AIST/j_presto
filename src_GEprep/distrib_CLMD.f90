
      subroutine distrib_CLMD

!*********************************************
!
!     DISTRIBution_CANONICAL_LAMBDA_DYNAMICS
!
!*********************************************

      use COMIFN ; use COMVAL

      implicit none

      integer(4)::i,j,k,l,n,ii,iST,iEN
      real(8):: E,EAA,EAB,rtmp,rtmp2,dx,dl,lower,upper,lowd,upd
      integer(4),allocatable:: Tnconf(:),nconf(:)
      real(8),allocatable::tF(:),F(:),Tlambda(:),lambda(:),COE(:,:),   &
        wei(:)
      logical(1),allocatable:: flg(:)

      real(8):: low_v_window(Nwindow),high_v_window(Nwindow)
      character(9999):: line
      logical(1):: ex

!*****************************
!     Input Lambda data & Calc mean force

      ! Check number of input files
      i = 0
      open(unit=1,file=trim(DATLST),status="old")
      do
        read(1,'(a)',iostat=j)line
        if ( j .ne. 0 ) exit
        i = i + 1
      enddo
      close(1) ; i = i / 2
      allocate(tF(i),Tlambda(i),Tnconf(i),flg(i))

      tF(1:i) = 0.d0 ; Tnconf(1:i) = 0 ; n = 0
      open(unit=1,file=trim(DATLST),status="old")
      do
        read(1,'(a)',end=801)line
        i = index(line,";")
        if ( i .ne. 0 ) line = trim(line(1:i-1))
        if ( len(line) .eq. 0 ) cycle
        read(1,*)iST,iEN
        if ( iST .le. 0 ) iST = 1
        inquire(file=trim(line),exist=ex)
        if ( .not. ex ) cycle

        write(6,'(2x,a)')"+ Now input "//trim(line)

        ! Input energy file
        n = n + 1
        open(unit=2,file=trim(line),status="old")
        do i = 1,iST-1
          read(2,*,end=800)
        enddo
        do i = iST,iEN
          read(2,*,end=800)E,EAA,EAB ; l = l + 1
          select case(METHOD)
          case ("CLMD")
            tF(n) = tF(n) - 2.d0*EAA*E - EAB
          case ("CLMD2")
            tF(n) = tF(n) - EAA - 0.5d0*EAB/sqrt(E)
          end select
          Tnconf(n) = Tnconf(n) + 1
        enddo
800     close(2)
        Tlambda(n) = E
        write(6,'(8x,i0,a3,i0)')iST," - ",i-1
      enddo
801   close(1) ; write(6,*)

!****************************

      ! Gather same lambda
      flg(1:n) = .true.
      do i = 1,n-1
        if ( .not. flg(i) ) cycle
        do j = i+1,n
          if ( .not. flg(j) ) cycle
          if ( Tlambda(i) .eq. Tlambda(j) ) then
            tF(i) = tF(i) + tF(j)
            Tnconf(i) = Tnconf(i) + Tnconf(j)
            flg(j) = .false.
          endif
        enddo
      enddo

      ! Calc. mean force
      k = count(flg(1:n)) ; j = 0
      allocate(F(k),lambda(k),nconf(k))
      do i = 1,n
        if ( flg(i) ) then
          j = j + 1
          F(j) = tF(i) / dble(Tnconf(i))
          lambda(j) = Tlambda(i)
          nconf(j) = Tnconf(i)
        endif
      enddo
      ! sort
      tF(1:k) = F(1:k)
      Tlambda(1:k) = lambda(1:k)
      Tnconf(1:k) = nconf(1:k)
      flg(1:k) = .true.
      do i = 1,k
        j = minloc(Tlambda(1:k),mask=flg(1:k),dim=1)
        lambda(i) = Tlambda(j)
        F(i) = tF(j)
        nconf(i) = Tnconf(j)
        flg(j) = .false.
      enddo
      
      rtmp = 1.d0 / (R*T)
      F(1:k) = F(1:k) * rtmp
      write(6,'(2x,a,i0,a,i0,a,f8.3,a)')"+ ",l," data are counted"
      write(6,*)

!*****************************
!     Output dlnP Function (*.nf)

      ! Make weighting factor
      allocate(wei(k),COE(0:NfitDIM,Nwindow))
      wei(1:k) = sqrt(dble(nconf(1:k)))
      write(6,'(2x,a)')"+ Fitting for dlnP"
      call lstsquare_final(Nwindow,k,NfitDIM,nxminE,nxmaxE,lambda(1:k),&
        F(1:k),wei(1:k),COE,low_v_window,high_v_window,rtmp,rtmp2,     &
        .true.)

      open(unit=1,file=trim(PROJNM)//".nf",status="replace")
      write(1,*)Nwindow,NfitDIM
      do i = 1,Nwindow
        do j = 0,NfitDIM
          write(1,*)COE(j,i)
        enddo
        write(1,*)low_v_window(i),high_v_window(i)
      enddo

      ! Calc. parabolic function for the out of the range
      !! lower side
      lowd = COE(0,1)
      lower = -forces
      !! upper side 
      upd = COE(0,Nwindow) ; dx = nxmaxE - low_v_window(Nwindow)
      do i = 1,NfitDIM
        upd = upd + COE(i,Nwindow)*dx**i
      enddo
      upper = forces
      write(1,*)lowd,upd
      write(1,*)nxminE,nxmaxE
      write(1,*)nxT
      write(1,*)lower,upper
      close(1)
      write(6,*)

!*****************************
!     Output dlnP data for plot

      open(unit=1,file=trim(PROJNM)//".dlnP",status="replace")
      j = 1
      do i = 1,k
        do 
          if ( j.ne.Nwindow .and. lambda(i).ge.high_v_window(j) ) then
            j = j + 1 ; cycle
          endif
          exit
        enddo
        dx = lambda(i) - low_v_window(j)
        rtmp = COE(0,j)
        do ii = 1,NfitDIM
          rtmp = rtmp + COE(ii,j)*dx**ii
        enddo
        if ( lambda(i) .lt. nxminE ) then
          dl = nxminE - lambda(i)
          rtmp2 = rtmp + lower*dl
        elseif ( lambda(i) .gt. nxmaxE ) then
          dl = lambda(i) - nxmaxE
          rtmp2 = rtmp + upper*dl
        else
          rtmp2 = rtmp
        endif
        write(1,'(4(e15.8,x))')lambda(i),F(i),rtmp,rtmp2
      enddo

      write(1,*)
      write(1,*)nxminE,lowd,lowd
      write(1,*)nxminE,upd,upd
      write(1,*)
      write(1,*)nxmaxE,lowd,lowd
      write(1,*)nxmaxE,upd,upd
      close(1)
 
!*****************************

      return
      end subroutine distrib_CLMD
