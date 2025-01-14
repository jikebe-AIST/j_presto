
      subroutine Qvalue_calc(nSTR)

!************************************************
!
!     Perform Qvalue calculation
!
!************************************************

      use COMIFN ; use COMVAL
      !$ use omp_lib

      implicit none

      ! STRucture number
        integer(4),intent(in):: nSTR

      integer(4):: i,j,k,l,icn,ires1,ires2,iatm1,iatm2,nc,ist,Ncluster
      real(4):: dcod(3),dist,dsum,tQval,traQ,rtmp
      logical(1):: Qcon_flg(nNCP),ec(nNCP),wc(nNCP),preFtraQcut(nNCP)

!****************************

      Qcon_flg(:) = .false.
      OUTOR : do i = 1,nNCP
        ires1 = iNCP(1,i) ; ires2 = iNCP(2,i)
        do k = iresQ(ires1),iresQ(ires1+1)-1
          iatm1 = iATMQ(k)
          do l = iresQ(ires2),iresQ(ires2+1)-1
            iatm2 = iATMQ(l)
            dcod(1:3) = cod(1:3,iatm1) - cod(1:3,iatm2)
            dcod(1:3) = dcod(1:3)                                      &
                        - bsize(1:3)*nint(dcod(1:3)*ibsize(1:3))
            dcod(1:3) = dcod(1:3)*dcod(1:3)
            rtmp = radQ(k) + radQ(l)
            rtmp = sum(dcod) - rtmp * rtmp
            if ( rtmp .le. 0.0 ) then
              Qcon_flg(i) = .true.
              cycle OUTOR
            endif
          enddo
        enddo
      enddo OUTOR
      icn = count(Qcon_flg(:))
      tQval = real(icn)/real(nNCP)
      Qvalue = Qvalue + tQval*wfac
      Qvalue2 = Qvalue2 + tQval*tQval*wfac
      where ( Qcon_flg ) Qcon = Qcon + wfac

      if ( .not. weakness_NCP_flag ) then
        write(uqv,'(i0,2(x,f))')nSTR,tQval,wfac
      else
        preFtraQcut(:) = FtraQcut(:)
        where(.not.Qcon_flg) FtraQcut = .true.
        icn = count(FtraQcut(:))
        traQ = max(1.0 - real(icn)/real(nNCP),0.0)
        write(uqv,'(i0,3(x,f))')nSTR,tQval,wfac,traQ

        if ( icn .ne. count(preFtraQcut(:)) ) then
          do i = 1,nNCP
            if ( preFtraQcut(i) ) cycle
            if ( .not.FtraQcut(i) ) cycle
            Dmax(i) = max(Dmax(i),traQ)
            where ( .not. FtraQcut ) 
              wc = .true.
            else where
              wc = .false.
            end where
            ec(:) = .false. ; ec(i) = .true. ; nc = 1
            do
              do j = 1,nNCP
                if ( .not. ec(j) ) cycle
                if ( wc(j) ) cycle
                wc(j) = .true. ; ires1 = iNCP(1,j) ; ires2 = iNCP(2,j)
                do k = 1,nNCP
                  if ( wc(k) ) cycle
                  if ( iNCP(1,k).eq.ires1 .or. iNCP(1,k).eq.ires2 .or. &
                       iNCP(2,k).eq.ires1 .or. iNCP(2,k).eq.ires2 )    &
                    ec(k) = .true.
                enddo
              enddo
              l = count(ec(:))
              if ( l .eq. nc ) then
!                where ( ec ) Qpoint(:,i) = Qpoint(:,i) + traQ
                where ( ec ) Qpoint(:,i) = Qpoint(:,i) + 1.d0
                exit
              else
                nc = l
              endif
            enddo
          enddo
        endif

        durability = durability + traQ
        if ( traQ .gt. 0.0 ) traQ_zero = traQ_zero + 1

!        where ( .not. FtraQcut ) 
!          wc = .true.
!        else where
!          wc = .false.
!        end where
!        ist = 1 ; Ncluster = 0
!        !! all
!        do
!          ec(:) = .false.
!          !! cluster
!          CLUSTER : do i = ist,nNCP
!            if ( wc(i) ) cycle
!            if ( FtraQcut(i) ) then
!              ec(i) = .true. ; ist = i + 1 ; nc = 1
!              do
!                do j = i,nNCP
!                  if ( wc(j) ) cycle
!                  if ( ec(j) ) then
!                    ires1 = iNCP(1,j) ; ires2 = iNCP(2,j)
!                    do k = i+1,nNCP
!                      if ( wc(k) ) cycle
!                      if ( iNCP(1,k).eq.ires1 .or. iNCP(1,k).eq.ires2 .or. &
!                           iNCP(2,k).eq.ires1 .or. iNCP(2,k).eq.ires2 ) &
!                        ec(k) = .true.
!                    enddo
!                  endif
!                enddo
!                l = count(ec(:))
!                if ( l .eq. nc ) then
!                  nc = l
!                  do j = i,nNCP
!                    if ( .not. ec(j) ) cycle
!                    if ( .not.preFtraQcut(j) .and. FtraQcut(j) ) then
!                      where (ec) Qpoint(:,j) = Qpoint(:,j) + traQ
!                    endif
!                  enddo
!                  where ( ec ) wc = .true.
!!                  Ncluster = Ncluster + 1
!!                  print*,"Cluster ",Ncluster,nc
!!                  do k = i,nNCP
!!                    if ( ec(k) ) then
!!                      l = iNCP(k)
!!                      print*," * ",con_res_Q(1:2,l)
!!                      Color(k) = dble(nc)
!!                    endif
!!                  enddo
!                  exit CLUSTER
!                else
!                  nc = l
!                endif
!              enddo
!            endif
!          enddo CLUSTER
!          if ( all(wc) ) exit
!        enddo
      endif

!****************************

      return
      end subroutine Qvalue_calc


!==================================================================


      subroutine Qprep

!**********************************************************
!
!     Preparation for Qvalue calc.
!
!**********************************************************

      use COMIFN ; use COMVAL

      implicit none

      integer(4):: icn,i,j,k,l,ires1,ires2,iatm1,iatm2,ierr,ii,jj,    &
        Ncon_res_Q
      integer(4),allocatable:: con_res_Q(:,:),ncon(:)
      real(4):: tcod(3,nQ),dcod(3),rtmp
      character(4):: tATMnm
      character(130):: FILnam,tmp,tmp2
      character(1):: ichn,jchn
      logical(4):: ex

!************************************

      !! Making residue list
      NresQ = 1 ; j = RESnum(iATMQ(1))
      do i = 2,nQ
        if ( RESnum(iATMQ(i)) .ne. j ) then
          NresQ = NresQ + 1 ; j = RESnum(iATMQ(i))
        endif
      enddo
      allocate(iresQ(NresQ+1))
      ii = 1 ; iresQ(1) = 1 ; j = RESnum(iATMQ(1))
      do i = 2,nQ
        if ( RESnum(iATMQ(i)) .ne. j ) then
          ii = ii + 1 ; iresQ(ii) = i ; j = RESnum(iATMQ(i))
        endif
      enddo
      iresQ(NresQ+1) = nQ + 1

      !! Making residue contact pair list
      Ncon_res_Q = 0
      do i = 1,NresQ-1
        ichn = CHN(iATMQ(iresQ(i)))
        do j = i+1,NresQ
          jchn = CHN(iATMQ(iresQ(j)))
          if ( ichn.eq.jchn .and.                                      &
               RESnum(iATMQ(iresQ(j)))-RESnum(iATMQ(iresQ(i))).le.     &
               cut_resnum_Qvalue ) cycle
          Ncon_res_Q = Ncon_res_Q + 1
        enddo
      enddo
      allocate(con_res_Q(2,Ncon_res_Q))
      ii = 0 ; jj = 0
      do i = 1,NresQ-1
        ichn = CHN(iATMQ(iresQ(i)))
        do j = i+1,NresQ
          jchn = CHN(iATMQ(iresQ(j)))
          if ( ichn.eq.jchn .and.                                      &
               RESnum(iATMQ(iresQ(j)))-RESnum(iATMQ(iresQ(i))).le.     &
               cut_resnum_Qvalue ) cycle
          ii = ii + 1
          con_res_Q(1:2,ii) = (/i,j/)
        enddo
      enddo

      ! Make vdW arraies
      allocate(radQ(nQ),radNCP(nQ))
      if ( len(trim(input_topology)) .ne. 0 ) then
        call mk_vdWrad_list(nATM,nQ,iATMQ,RESnm,ATMnm,radQ)
        radNCP(:) = radQ(:)
        radQ(:) = radQ(:) + toleQ*0.5d0
        radNCP(:) = radNCP(:) + toleNCP*0.5d0
      else
        radQ(:) = toleQ * 0.5d0 ; radNCP(:) = toleNCP * 0.5d0
      endif

      ! Search native contact pairs
      open(unit=1,file=trim(input_reference_PDB_list),status="old")
      allocate(ncon(ncon_res_Q)) ; ncon(1:Ncon_res_Q) = 0 ; icn = 0
      do
        read(1,'(a)',end=800)FILnam
        if ( FILnam .eq. " " ) cycle
        inquire(file=FILnam,exist=ex)
        if ( .not. ex ) cycle
        icn = icn + 1

        !! read each PDB files
        open(unit=2,file=trim(FILnam),status="old") ; k = 1 ; j = 0
        do
          read(2,'(a)',end=801)tmp
          if ( tmp(1:6).ne."ATOM  " .and. tmp(1:6).ne."HETATM" ) cycle
          j = j + 1
          if ( j .eq. iATMQ(k) ) then
            read(tmp,'(30x,3f8.3)')tcod(1:3,k)
            k = k + 1
            if ( k .gt. nQ ) exit
          endif
        enddo
801     close(2)

        !! search native contact pairs
        OUTER : do i = 1,Ncon_res_Q
          ires1 = con_res_Q(1,i) ; ires2 = con_res_Q(2,i)
          do iatm1 = iresQ(ires1),iresQ(ires1+1)-1
          do iatm2 = iresQ(ires2),iresQ(ires2+1)-1
            dcod(1:3) = tcod(1:3,iatm1) - tcod(1:3,iatm2)
            dcod(1:3) = dcod(1:3)                                      &
                        - bsize(1:3)*nint(dcod(1:3)*ibsize(1:3))
            dcod(1:3) = dcod(1:3)*dcod(1:3)
            rtmp = radNCP(iatm1) + radNCP(iatm2)
            rtmp = rtmp * rtmp
            if ( sum(dcod) .le. rtmp) then
              ncon(i) = ncon(i) + 1 ; cycle OUTER
            endif
          enddo
          enddo
        enddo OUTER
      enddo
800   close(1)

      !! Determination of the native contact pairs
      open(unit=uqv,file=trim(d_proj)//".Qval",status="replace")
      nNCP = 0
      do i = 1,Ncon_res_Q
        if ( dble(ncon(i))/dble(icn) .ge. NCP_rate ) nNCP = nNCP+1
      enddo
      write(uqv,'(4x,a,i0)')"# The number of native contact pairs is ",&
                            nNCP
      write(uqv,*)
      allocate(iNCP(2,nNCP))
      j = 0
      do i = 1,Ncon_res_Q
        if ( dble(ncon(i))/dble(icn) .ge. NCP_rate ) then
          j = j + 1 ; iNCP(1:2,j) = con_res_Q(1:2,i)
          write(uqv,'(a1,4x,i0,5x,i0,a3,i0)')"#",j,                    &
            RESnum(iATMQ(iresQ(iNCP(1,j))))," - ",                   &
            RESnum(iATMQ(iresQ(iNCP(2,j))))
        endif
      enddo
      write(uqv,*)

      write(6,*)
      write(6,'(4x,a,i0)')"+ Native contact pairs for Qvalue : "&
        ,nNCP
      if ( loglvl.eq."l" .and. nNCP.gt.6 ) then
        write(6,'(8x,a)')"   Res.1  -  Res.2"
        do k = 1,nNCP
          i = iNCP(1,k) ; j = iNCP(2,k)
          write(6,'(i8,a3,i0,a3,i0)')k," : ",RESnum(iATMQ(iresQ(i))),  &
            " - ",RESnum(iATMQ(iresQ(j)))
        enddo
      elseif ( loglvl .eq. "s" ) then
        write(6,'(8x,a)')"   Res.1  -  Res.2"
        do k = 1,3
          i = iNCP(1,k) ; j = iNCP(2,k)
          write(6,'(i8,a3,i0,a3,i0)')k," : ",RESnum(iATMQ(iresQ(i))),  &
            " - ",RESnum(iATMQ(iresQ(j)))
        enddo
        write(6,*)"              ..."
        do k = nNCP-2,nNCP
          i = iNCP(1,k) ; j = iNCP(2,k)
          write(6,'(i8,a3,i0,a3,i0)')k," : ",RESnum(iATMQ(iresQ(i))),  &
            " - ",RESnum(iATMQ(iresQ(j)))
        enddo
      endif

      allocate(Qcon(nNCP)) ; Qcon(:) = 0.d0
      if ( weakness_NCP_flag ) then
        allocate(Qpoint(nNCP,nNCP)) ; Qpoint(:,:) = 0.d0
        allocate(FtraQcut(nNCP)) ; FtraQcut(:) = .false.
        allocate(Dmax(nNCP)) ; Dmax(:) = 0.d0
!        allocate(Color(nNCP)) ; Color(:) = 0.d0
        open(unit=udr,file=trim(d_proj)//".drb",status="replace")
      endif

!************************************

      return
      end subroutine Qprep


!==================================================================


      subroutine Qvalue_final

      use COMVAL
      use COMIFN, only : weakness_NCP_flag
      implicit none

      integer(4):: i,j,ires1,ires2,k,l,iatm1,iatm2,a(1),rank(nNCP)
      real(8):: wNCP(nNCP),rtmp,rtmp2
      logical(1):: flg(nNCP)

      write(uqv,*)
      Qvalue = Qvalue / sumW2 ; Qvalue2 = Qvalue2 / sumW2
      Qvalue2 = sqrt(Qvalue2 - Qvalue*Qvalue)
      write(uqv,*)"# Average Qvalue"
      write(uqv,*)"# ",Qvalue," (",Qvalue2,")"
      close(uqv)

!      rtmp = maxval(Color)
!      Color(:) = Color(:) / rtmp
      Qcon(:) = Qcon(:) / sumW2
 
      if ( weakness_NCP_flag ) then
        wNCP(:) = 0.d0
        do i = 1,nNCP
          wNCP(i) = sum(Qpoint(i,:))
        enddo

        flg(:) = .true. ; i = 0
        do
          i = i + 1
          if ( i .gt. nNCP ) exit
!          a = maxloc(Dmax(:),flg(:))
          a = maxloc(wNCP(:),flg(:))
          flg(a(1)) = .false.
          rank(a(1)) = i
        enddo
        open(unit=1,file=trim(d_proj)//".wNCP",status="replace")
        write(1,*)"# pairs   wNCP rank  Dmax       Qcon"
        do i = 1,nNCP
          ires1 = iNCP(1,i) ; ires2 = iNCP(2,i)
          write(1,'(i0,a1,i0,x,f12.5,x,i5,2(x,f12.5))')                &
            ires1,"-",ires2,wNCP(i),rank(i),Dmax(i),Qcon(i)
        enddo
        close(1)

        rtmp = maxval(wNCP(:)) ; rtmp2 = minval(wNCP(:))
        if ( rtmp-rtmp2 .ne. 0.d0 )                                    &
          wNCP(:) = (wNCP(:)-rtmp2) / (rtmp-rtmp2)

        open(unit=1,file=trim(d_proj)//"_wNCP.bild",status="replace")
        OUTER : do i = 1,nNCP
          ires1 = iNCP(1,i) ; ires2 = iNCP(2,i)
          do k = iresQ(ires1),iresQ(ires1+1)-1
            iatm1 = iATMQ(k)
            if ( ATMnm(iatm1) .eq. " CA " ) then
              do l = iresQ(ires2),iresQ(ires2+1)-1
                iatm2 = iATMQ(l)
                if ( ATMnm(iatm1) .eq. " CA " ) then
!                  rtmp = 1.d0 - Color(i)
!                  rtmp = 1.d0 - wNCP(i)
                  rtmp = 1.d0 - dble(max((31-rank(i)),0))/30.d0
                  write(1,'(a,3f8.3)')".color ",1.d0,rtmp,rtmp
                  write(1,'(a,6f8.3,x,f8.3)')".cylinder ",Rcod(:,iatm1), &
                    Rcod(:,iatm2),0.1d0
                  cycle OUTER
                endif
              enddo
            endif
          enddo
        enddo OUTER
        close(1)

        close(udr)
      endif

!      open(unit=1,file=trim(d_proj)//"_Qscr.pdb",status="replace")
!      do i = istATM,ienATM
!        rtmp = wNCP(RESnum(i))
!        if ( RESnum(i) .lt. 10000 ) then
!          write(1,'(a4,i7,x,a4,x,a4,a1,i4,4x,3f8.3,2f6.3)')      &
!            "ATOM",i,ATMnm(i),RESnm(i),CHN(i),RESnum(i),    &
!            Rcod(1:3,i),0.0,rtmp
!        elseif ( RESnum(i) .lt. 100000 ) then
!          write(1,'(a4,i7,x,a4,x,a4,a1,i5,3x,3f8.3,2f6.3)')      &
!            "ATOM",i,ATMnm(i),RESnm(i),CHN(i),RESnum(i),    &
!            Rcod(1:3,i),0.0,rtmp
!        else
!          write(1,'(a4,i7,x,a4,x,a4,a1,i6,2x,3f8.3,2f6.3)')      &
!            "ATOM",i,ATMnm(i),RESnm(i),CHN(i),RESnum(i),    &
!            Rcod(1:3,i),0.0,rtmp
!        endif
!      enddo
!      close(1)

      return
      end subroutine Qvalue_final
