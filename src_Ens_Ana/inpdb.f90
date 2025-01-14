
      subroutine inpdb

!****************************************
!
!     INput PDB for coordinates
!
!****************************************

      use COMIFN ; use COMVAL

      implicit none

      ! pdb FILe NAMe
        character(130):: FILnam
      ! Potential Energy of each structure
        real(4):: ep

      integer(4):: icn,i,j,k,k2,ii
      character(130):: tmp
      real(8):: rr,ttepotcc(3)
      logical(1):: remake_cod = .false.
      integer(4):: irmk
      real(4),allocatable:: rmkcod(:,:)
      integer(4),allocatable:: rmkidx(:)

      integer(4):: ichance = 0
      integer(4):: ichance2 = 0

!********************************

      write(6,'(4x,a)')"+ PDB files (in input_PDB_list) input START"

      if ( output_binary_flag )                                        &
        open(unit=4,file=trim(d_proj)//".cod",status="replace",        &
             form="unformatted")
      if ( len(trim(input_prob_file)) .ne. 0 )                         &
        open(unit=9,file=trim(d_proj)//".weight",status="replace")
      if ( len(trim(input_weight_file)) .ne. 0 )                       &
        open(unit=10,file=trim(input_weight_file),status="old")
      if ( n_dist_c .ne. 0 )                                           &
        open(unit=udc,file=trim(d_proj)//".distan",status="replace")
      if ( n_ang_c .ne. " " )                                          &
        open(unit=uac,file=trim(d_proj)//".ang",status="replace")
      if ( n_dih_c .ne. " " )                                          &
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
!     Input PDB files

      open(unit=1,file=trim(input_PDB_list),status="old") ; icn = 0
      ! For weight
      sumW = 0.d0 ; sumW2 = 0.d0
      if ( len(trim(input_prob_file)) .ne. 0 )                         &
          open(unit=3,file=trim(input_energy_file),status="old")
       
      do
        read(1,'(a)',end=800)tmp

        ! Provision for blank columns
        k = index(tmp,";")
        if ( k .ne. 0 ) tmp = tmp(1:k)
        if ( tmp .eq. " " ) cycle

        icn = icn + 1
        FILnam = trim(tmp)
        if ( mod(icn,100) .eq. 0 ) then
          write(6,'(4x,a)')"* Now, input "//trim(FILnam)
        endif

        ! Weighting factor
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
        ichance2 = ichance2 + 1 ; sumW2 = sumW2 + wfac
        if ( len(trim(input_prob_file)) .ne. 0 ) then
          if ( GE_method .eq. "MULT" ) then
            write(9,'(i0,x,f,e12.5)')icn,wfac,ep
          elseif ( GE_method .eq. "ALSD" ) then
            write(9,'(i0,x,f,x,f8.3,3(x,e12.5))')icn,wfac,0.0,         &
                                                ttepotcc(1:3)
          endif
        endif

        ! Read each PDB file
        open(unit=2,file=trim(FILnam),status="old") ; j = 0
        do
          read(2,'(a)',end=801)tmp
          if ( tmp(1:6).ne."ATOM  " .and. tmp(1:6).ne."HETATM" ) cycle
          j = j + 1 ; read(tmp,'(30x,3f8.3)')cod(1:3,j)
        enddo
801     close(2)

        ! NOE reproduction ratio calc.
        if ( NMR_flag ) call noe_calc1

        ! Qvalue calculation
        if ( len(trim(atom_spec_Qvalue)).ne.0 ) call Qvalue_calc(icn)

        ! RMSD calculation
        if ( len(trim(RMSD_method)) .ne. 0 ) call RMSDprep(icn,ichance2)

        ! Making PCA coordinates
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
          write(4)0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0,0.0,ttepotcc(1:3)
          if ( .not. remake_cod ) then
            write(4)cod(1:3,istATM:ienATM)
          else
            rmkcod(1:3,1:irmk) = cod(1:3,rmkidx(1:irmk))
            write(4)rmkcod(1:3,1:irmk)
          endif
        endif
      enddo
800   close(1)
      if ( output_binary_flag ) close(4)
      if ( len(trim(input_prob_file)) .ne. 0 ) close(3)
      if ( DSSP_flag ) call DSSP(0)
      if ( len(trim(RMSD_method)) .ne. 0 ) call RMSD_final
      if ( len(trim(atom_spec_Rg)) .ne. 0 ) call Rg_final
      if ( len(trim(atom_spec_PCA)) .ne. 0 ) call PCA_final
      if ( len(trim(atom_spec_Qvalue)) .ne. 0 ) call Qvalue_final
      if ( len(trim(atom_spec_contact)).ne.0 ) call contact_final(sumW2)
      if ( RMSF_flag ) call RMSF_final

      ! NOE reproduction ratio calc.
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
      close(9)
      if ( len(trim(input_weight_file)).ne.0 ) close(10)
      write(6,'(4x,a,i0,a)')"+ ",ichance2,                             &
              " conf. were used for analysis or output"
      if ( n_dist_c .ne. 0 ) close(udc)
      if ( n_ang_c .ne. 0 ) close(uac)
      if ( n_dih_c .ne. 0 ) close(utc)

!********************************

      return
      end subroutine inpdb
