
      subroutine inbin

!************************************************
!
!     INput BINary files for coordinates
!
!************************************************

      use COMIFN ; use COMVAL

      implicit none

      ! binary FILe NAMe
        character(130):: FILnam,FILnam_old
      ! STart conformation number for binary input
        integer(4):: iST,iiST
      ! ENd conformation number for binary input
        integer(4):: iEN
      ! ICouNt for input conformations
        integer(4):: icn = 0
        integer(4):: ichance = 0
        integer(4):: ichance2 = 0
        integer(4):: TTPncnf = 0
      ! Potential Energy of each structure
        real(4):: ep

      integer(4):: i,j,k,k2,ii
      character(130):: tmp
      real(8):: rr,ranf,each_sumW2
      logical(4):: TTPflg
      logical(1):: ex,Fstream
      logical(1):: remake_cod = .false.
      integer(4):: irmk
      real(4),allocatable:: rmkcod(:,:)
      integer(4),allocatable:: rmkidx(:)  

      ! presto output
      integer(4):: istp,iyn15v,iyn15h,iicn,old_icn,old_ichance2
      real(4):: sitime,sec,et,ek,temp,lambda,rmsd2,ttepotcc(3)

!*********************************

      write(6,'(2x,a)')"+ BINARY files input START"

      if ( output_binary_flag )                                        &
        open(unit=4,file=trim(d_proj)//".cod",status="replace",        &
             form="unformatted")
      open(unit=9,file=trim(d_proj)//".weight",status="replace")
      open(unit=8,file=trim(d_proj)//".TTPncnf",status="replace")
      if ( len(trim(input_weight_file)) .ne. 0 )                       &
        open(unit=10,file=trim(input_weight_file),status="old")
      if ( n_dist_c .ne. 0 )                                           &
        open(unit=udc,file=trim(d_proj)//".distan",status="replace")
      if ( n_ang_c .ne. 0 )                                            &
        open(unit=uac,file=trim(d_proj)//".ang",status="replace")
      if ( n_dih_c .ne. 0 )                                            &
        open(unit=utc,file=trim(d_proj)//".dih",status="replace")
      if ( RMSF_flag ) call RMSF_init
      irmk = count(outF(istATM:ienATM))
      if ( irmk .ne. ienATM-istATM+1 ) then
        remake_cod = .true. ; allocate(rmkcod(3,irmk),rmkidx(irmk))
        k = 0
        do j = istATM,ienATM
          if ( outF(j) ) then
            k = k + 1 ; rmkidx(k) = j
          endif
        enddo
      endif

!*****************************
!     Input BINARY files

      sumW = 0.d0 ; sumW2 = 0.d0 ; old_icn = 0 ; old_ichance2 = 0
      FILnam_old = " "
      open(unit=1,file=trim(input_binary_list),status="old")
      do
        call infil(1,FILnam,iST,iEN,rr,TTPflg,i)
        if ( i .ne. 0 ) goto 800

        ! File existence check
        inquire(file=trim(FILnam),exist=ex)
        if ( .not.ex .or. iEN.lt.iST ) then
          if ( TTPflg ) then
            write(8,*)TTPncnf ; TTPncnf = 0
          endif
          cycle
        endif

        ! Input BINARY file
        if ( trim(FILnam).ne.trim(FILnam_old) .or. iST.lt.iiST ) then
          iiST = 1 ; FILnam_old = FILnam
          inquire(unit=2,opened=ex)
          if ( ex ) close(2)
          j = len_trim(FILnam) ; Fstream = .false.
          if ( adjustl(FILnam(j-4:j)) .eq. ".cods" ) Fstream = .true.
          if ( .not. Fstream ) then
            open(unit=2,file=trim(FILnam),form="unformatted",          &
                 status="old")
          else
            open(unit=2,file=trim(FILnam),form="unformatted",          &
                 access="stream",status="old")
          endif
        endif
        iicn = 0 ; each_sumW2 = 0.d0
        ! Skip reading
        if ( .not. Fstream ) then
          do j = iiST,iST-1
            read(2,end=801)
            read(2,end=801)
          enddo
        else
          do j = iiST,iST-1
            read(2,end=801)istp,sitime,sec,et,ek,temp,ep,lambda,iyn15v,&
                           iyn15h,rmsd2,ttepotcc(1:3)
            read(2,end=801)cod(1:3,1:nATM)
          enddo
        endif
        ! Reading
        do j = iST,iEN
          read(2,end=801)istp,sitime,sec,et,ek,temp,ep,lambda,iyn15v,  &
                         iyn15h,rmsd2,ttepotcc(1:3)
          read(2,end=801)cod(1:3,1:nATM)
          icn = icn + 1 ; iicn = iicn + 1
!          if ( mod(icn,1000) .eq. 0 ) then
!            write(tmp,*)icn
!            write(6,'(4x,a)')"* Now input "//trim(tmp)//" -th conf."
!          endif

          ! Weighting factor
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
            read(10,*)ii,wfac
          else
            wfac = 1.d0
          endif
          ! Skip calc. for conformation with very low probability
          ! Select conformation by wfac
          if ( weight_method .eq. "SELECT" ) then
            rseed = 843314861*rseed+453816693
            rr = real(rseed)/2**31 * 0.5d0 + 0.5d0
            if ( wfac.lt.rr ) cycle
          elseif ( weight_method .eq. "WEIGHT" ) then
            if ( wfac .lt. weight_threshold ) cycle
          endif
          ichance = ichance + 1 ; sumW = sumW + wfac
          if ( mod(ichance,interval) .ne. 0 ) cycle
          if ( weight_method .eq. "SELECT" ) wfac = 1.d0
          ichance2 = ichance2 + 1 ; each_sumW2 = each_sumW2 + wfac
          write(9,'(i0,x,f,x,f8.3,4(x,e12.5))')icn,wfac,lambda,ep,     &
                                              ttepotcc(1:3)
          
          ! NOE intensity calc.
          if ( NMR_flag ) call noe_calc1

          ! Qvalue calculation
          if ( len(trim(atom_spec_Qvalue)).ne.0 ) call Qvalue_calc(icn)

          ! RMSD calculation
          if ( len(trim(RMSD_method)).ne.0 ) call RMSDprep(icn,ichance2)

          ! Making PCA coordinates or Qvalue calculation
          if ( len(trim(atom_spec_PCA)) .ne. 0 ) call PCA_mkcod

          !! distribution analysis
          if ( d_celsiz .gt. 0.d0 ) call distribution

          ! Rg calc.
          if ( len(trim(atom_spec_Rg)) .ne. 0 ) call Rgcalc(icn)

          ! Contact analysis
          if ( len(trim(atom_spec_contact)) .ne. 0 ) call contact(icn)

          if ( n_dist_c .ne. 0 ) call distance_calc
          if ( n_ang_c .ne. 0 ) call angle_calc
          if ( n_dih_c .ne. 0 ) call dihedral_calc

          ! RMSF calc.
          if ( RMSF_flag ) call RMSF_calc

          !! pocket search analysis
          if ( p_celsiz .gt. 0.d0 ) call pocket_search

          ! Output PDB
          if ( output_PDB_flag .or. DSSP_flag ) then
            call outpdb(icn,tmp)
            if ( DSSP_flag ) call DSSP(icn)
            if ( .not. output_PDB_flag ) call system("rm "//trim(tmp))
          endif
          ! Output BINARY
          if ( output_binary_flag ) then
            write(4)istp,sitime,sec,et,ek,temp,ep,lambda,iyn15v,       &
                      iyn15h,rmsd2,ttepotcc(1:3)
            if ( .not. remake_cod ) then
              write(4)cod(1:3,istATM:ienATM)
            else
              rmkcod(1:3,1:irmk) = cod(1:3,rmkidx(1:irmk))
              write(4)rmkcod(1:3,1:irmk)
            endif
          endif
        enddo
801     sumW2 = sumW2 + each_sumW2
        TTPncnf = TTPncnf + ichance2 - old_ichance2
        if ( TTPflg ) then
          write(8,*)TTPncnf ; TTPncnf = 0
          if ( allocated(FtraQcut) ) FtraQcut(:) = .false.
          if ( weakness_NCP_flag ) then
            write(udr,*)durability,traQ_zero
            durability = 0.d0 ; traQ_zero = 0
          endif
        endif
        iiST = iST+iicn
        write(6,'(6x,a)')"* Now input "//trim(FILnam)
        write(6,'(18x,i0,a3,i0)')iST," - ",iiST-1
        write(6,"(12x,a4,i0,a13,i0,a9,f8.3,a4)")"I : ",iicn,           &
          " conf.s, U : ",ichance2-old_ichance2," conf.s (",           &
          each_sumW2/dble(iicn)*100.d0," % )"
        write(6,"(16x,i0,a3,i0,a10)")old_icn+1," - ",icn,"-th conf.s"
        old_icn = icn ; old_ichance2 = ichance2
      enddo
800   close(1)
      inquire(unit=2,opened=ex)
      if ( ex ) close(2)
      if ( output_binary_flag ) close(4)
      if ( DSSP_flag ) call DSSP(0)
      if ( len(trim(RMSD_method)) .ne. 0 ) call RMSD_final
      if ( len(trim(atom_spec_Rg)) .ne. 0 ) call Rg_final
      if ( len(trim(atom_spec_PCA)) .ne. 0 ) call PCA_final
      if ( len(trim(atom_spec_Qvalue)) .ne. 0 ) call Qvalue_final
      if ( len(trim(atom_spec_contact)).ne.0 ) call contact_final(sumW2)
      if ( RMSF_flag ) call RMSF_final

!***
      ! NOE intensity calc.
      if ( NMR_flag ) call noe_calc2

      ! distribution analysis
      if ( d_celsiz .gt. 0.d0 ) call distrib_final
      ! pocket search analysis
      if ( p_celsiz .gt. 0.d0 ) call pocket_search_final

!***
      ! final treatment
      write(6,'(4x,a2,i0,a)')"+ ",icn," conf.s are input"
      if ( len(trim(input_prob_file)).ne.0 .or.                        &
           len(trim(input_weight_file)).ne.0 )                         &
        write(6,'(8x,a,f)')"Sum of weight = ",sumW
      close(9) ; close(8)
      if ( len(trim(input_weight_file)).ne.0 ) close(10)
      write(6,'(4x,a,i0,a)')"+ ",ichance2,                             &
              " conf. were used for analysis or output"
      if ( n_dist_c .ne. 0 ) close(udc)
      if ( n_ang_c .ne. 0 ) close(uac)
      if ( n_dih_c .ne. 0 ) close(utc)

!**********************************

      return
      end subroutine inbin
