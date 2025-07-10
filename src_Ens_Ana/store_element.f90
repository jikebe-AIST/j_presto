
      subroutine store_element

!**************************************************
!
!     store element data to variable of input
!
!**************************************************

      use COMINP ; use COMIFN ; use COMVAL

      implicit none

      integer(4):: iERR = 0
      integer(4):: i = 0

      integer(4):: j,k,ii,jj
      character(999):: tmp,tmp1,tmp2,Atmp(4)
      real(8):: rtmp
      logical(4):: ex

      integer(4),allocatable:: tFresT(:),tLresT(:),tFresR(:),tLresR(:),&
        tlist(:)
      character(130),allocatable:: talist(:)
      character(1):: YN
      character(4):: tATMnm
      logical(1),allocatable:: mark(:),markR(:)
      logical(1):: f

!***************************

      ! overwrite_flag
      if ( ELEMNT(1)(1:1).eq."Y" .or. ELEMNT(1)(1:1).eq."y" )          &
        overwrite_flag = .true.

      ! project_name
      project_name = ELEMNT(2)
      write(6,'(2x,a)')"* project_name : "
      write(6,'(6x,a)')trim(project_name)
      d_proj = "Ens_Ana_"//trim(project_name)
      inquire(directory=trim(d_proj),exist=ex)
      if ( ex ) then
        if ( .not. overwrite_flag ) then
          write(6,'(2x,a)')"!! CAUTION !!"
          write(6,'(6x,a)')"Output directory "//trim(d_proj)
          write(6,'(6x,a)')"  has already existed"
          write(6,'(6x,a)')"Do you OVERWRITE it? ( Y or N )"
          open(unit=7,file="/dev/tty",status="replace")
          read(7,*)YN
          close(7)
          write(6,*)
          if ( YN.eq."Y" .or. YN.eq."y" ) then
            write(6,'(2x,a)')"The directory is OVERWRITTEN"
          else
            stop
          endif
        endif
      else
        call system("mkdir "//trim(d_proj))
      endif
      d_proj2 = trim(d_proj)//"/"
      d_proj = trim(d_proj)//"/"//trim(project_name)

      ! log_level
      loglvl = ELEMNT(3)
      if ( loglvl.eq."l" .or. loglvl.eq."L" ) then
        loglvl = "l"
        write(6,'(2x,a)')"* Output log level : long"
      elseif ( loglvl.eq."s" .or. loglvl.eq."S" ) then
        write(6,'(2x,a)')"* Output log level : short"
      else
        write(6,'(2x,a)')"* Output log level : No"
      endif

      ! template_PDB
      write(6,*)
      template_PDB = ELEMNT(4)
      if ( len(trim(template_PDB)) .ne. 0 ) then
        write(6,'(2x,a)')"* template_PDB : "
        write(6,'(6x,a)')trim(template_PDB)
        inquire(file=trim(template_PDB),exist=ex)
        if ( .not. ex ) call error(10201)
        call intmp(template_PDB)
      else
        call error(10202)
      endif
      allocate(mark(nATM))

      ! input_topology
      input_topology = ELEMNT(5)
      if ( len(trim(input_topology)) .ne. 0 ) then
        write(6,*)
        write(6,'(2x,a)')"* input_topology : "
        write(6,'(6x,a)')trim(input_topology)
        inquire(file=trim(input_topology),exist=ex)
        if ( .not. ex ) call error(10201)
        call inptpl(5,6,input_topology,iERR)
      endif

      ! input_PDB_list
      input_PDB_list = ELEMNT(6)
      if ( len(trim(input_PDB_list)) .ne. 0 ) then
        write(6,'(2x,a)')"* input_PDB_list : "
        write(6,'(6x,a)')trim(input_PDB_list)
        inquire(file=trim(input_PDB_list),exist=ex)
        if ( .not. ex ) call error(10201)
      else
        if ( len(trim(ELEMNT(7)//ELEMNT(8)//ELEMNT(9))) .eq. 0 ) then
          input_PDB_list = trim(d_proj2)//".input_PDB_list"
          call system("echo "//trim(template_PDB)//" > "//          &
                      trim(input_PDB_list))
        endif
      endif

      input_binary_list = ELEMNT(7)
      input_otherfmt_list = ELEMNT(8)
      input_restart_list = ELEMNT(9)
      ! input_binary_list
      if ( len(trim(input_binary_list)) .ne. 0 ) then
        write(6,'(2x,a)')"* input_binary_list : "
        write(6,'(6x,a)')trim(input_binary_list)
        inquire(file=trim(input_binary_list),exist=ex)
        if ( .not. ex ) call error(10201)
      ! input_otherfmt_list
!      elseif ( len(trim(input_otherfmt_list)) .ne. 0 ) then
!        write(6,'(2x,a)')"* input_otherfmt_list : "
!        write(6,'(6x,a)')trim(input_otherfmt_list)
!        inquire(file=trim(input_otherfmt_list),exist=ex)
!        if ( .not. ex ) call error(10201)
      ! input_restart_list
      elseif ( len(trim(input_restart_list)) .ne. 0 ) then
        write(6,'(2x,a)')"* input_restart_list : "
        write(6,'(6x,a)')trim(input_restart_list)
        inquire(file=trim(input_restart_list),exist=ex)
        if ( .not. ex ) call error(10201)
      endif

      ! output_PDB_flag
      if ( len(trim(input_PDB_list//input_binary_list//                &
           input_otherfmt_list//input_restart_list)) .ne. 0 ) then
        if ( ELEMNT(10)(1:1).eq."Y" .or. ELEMNT(10)(1:1).eq."y" ) then
          output_PDB_flag = .true.
          d_PDB = trim(d_proj)//"_PDB"
          write(6,'(2x,a)')"* The PDB files are output in ("//         &
                           trim(d_PDB)//" ) directory"
          inquire(directory=trim(d_PDB),exist=ex)
          !! over-write
          if ( ex ) then
            write(6,'(4x,a)')"! CAUTION !"
            write(6,'(4x,a)')"PDB output directory is OVERWRITTEN"
            call systemqq("rm "//trim(d_PDB)//"/*.pdb > /dev/null 2>&1")
          !! Make the directory for PDB files output
          else
            call systemqq("mkdir "//trim(d_PDB))
          endif
          d_PDB = trim(d_PDB)//"/"
        endif
      endif

      ! output_binary_flag
      if ( ELEMNT(11)(1:1).eq."Y" .or. ELEMNT(11)(1:1).eq."y" ) then
        output_binary_flag = .true.
        write(6,'(2x,a)')"* Output a binary file"
        inquire(file=trim(d_proj)//".cod",exist=ex)
        !! over-write
        if ( ex ) then
          write(6,'(4x,a)')"! CAUTION !"
          write(6,'(4x,a)')trim(project_name)//                        &
                           ".cod file is OVERWRITTEN"
          write(6,'(4x,a)')"OK? (if you OK, enter Y)"
          read(7,*)tmp
          if ( tmp(1:1).ne."Y" .and. tmp(1:1).ne."y" ) stop
        endif
      endif

      ! reference_PDB
      if ( len(trim(ELEMNT(12))) .ne. 0 ) then
        reference_PDB = ELEMNT(12)
        write(6,'(2x,a)')"* reference_PDB : "
        write(6,'(6x,a)')trim(reference_PDB)
        if ( len(trim(reference_PDB)) .ne. 0 ) then
          inquire(file=trim(reference_PDB),exist=ex)
          if ( .not. ex ) call error(10201)
          call inref(reference_PDB)
          allocate(markR(nATMr))
        endif
      endif

      ! input reference cord. list
      input_reference_PDB_list = ELEMNT(13)
      input_reference_binary_list = ELEMNT(14)
      !! ref_conf_1_flag
      if ( len(trim(input_reference_PDB_list//                         &
           input_reference_binary_list)) .eq. 0 ) ref_conf_1_flag = .true.
      !! input_reference_PDB_list
      if ( len(trim(input_reference_PDB_list)) .ne. 0 ) then
        if ( len(trim(input_reference_binary_list)) .ne. 0 )           &
          call error(10217)
        write(6,'(2x,a)')"* input_reference_PDB_list : "
        write(6,'(6x,a)')trim(input_reference_PDB_list)
        inquire(file=trim(input_reference_PDB_list),exist=ex)
        if ( .not. ex ) call error(10201)
      !! input_reference_binary_list
      elseif ( len(trim(input_reference_binary_list)) .ne. 0 ) then
        write(6,'(2x,a)')"* input_reference_binary_list : "
        write(6,'(6x,a)')trim(input_reference_binary_list)
        inquire(file=trim(input_reference_binary_list),exist=ex)
        if ( .not. ex ) call error(10201)
      else
        input_reference_PDB_list = trim(d_proj2)//                   &
                                   ".input_reference_PDB_list"
        call system("echo "//trim(reference_PDB)//" > "//            &
                    trim(input_reference_PDB_list))
      endif

      ! PDB_consistency_check
      if ( ELEMNT(15)(1:1).eq."Y" .or. ELEMNT(15)(1:1).eq."y" ) then
        if ( len(trim(input_PDB_list)) .ne. 0 ) then
          write(6,'(2x,a)')"* The consistency of input PDB files are " &
                           //"checked."
          call pdb_chk(template_PDB,input_PDB_list)
        endif
        if ( len(trim(reference_PDB)).ne.0 .and.                       &
          index(input_reference_PDB_list,".input_reference_PDB_list")  &
          .eq.0 ) then
          write(6,'(2x,a)')"* The consistency of reference PDB files " &
                         //"are checked."
          call pdb_chk(reference_PDB,input_reference_PDB_list)
        endif
      else
        write(6,'(2x,a)')"* PDB files coidentity check is SKIPPED"
      endif

      ! atom_spec_output
      atom_spec_output = ELEMNT(16)
      allocate(tlist(nATM),outF(nATM)) ; outF(:) = .false.
      call atom_specifier(len(trim(atom_spec_output)),                 &
        trim(atom_spec_output),nATM,ATMnum,ATMnm,RESnum,RESnm,nCHN,CHN,&
        Rcod,j,tlist)
      write(6,'(2x,a)')"* Atoms output in either PDB or binary files"
      call output_specifier_log(j,loglvl,nATM,ATMnum,ATMnm,RESnum,     &
        RESnm,nCHN,CHN,tlist)
      outF(tlist(1:j)) = .true.
      istATM = tlist(1) ; ienATM = tlist(j)
      fstRES = minval(RESnum(tlist(1:j)))
      fnlRES = maxval(RESnum(tlist(1:j)))
      deallocate(tlist)

      ! analysis & output data interval
      read(ELEMNT(17),*,iostat=iERR)interval
      if ( interval .le. 0 ) interval = 1
      if ( interval .gt. 1 ) then
        write(6,'(2x,a)')"* interval for data analysis : "
        write(6,'(6x,i0)')interval
      endif

      ! input_prob_file
      input_prob_file = ELEMNT(18)
      if ( len(trim(input_prob_file)) .ne. 0 ) then
        write(6,'(2x,a)')"* The conformational ensemble is "//         &
                         "re-weighted"
        write(6,'(2x,a)')"* Probabiligy file : "
        write(6,'(6x,a)')trim(input_prob_file)
        inquire(file=trim(input_prob_file),exist=ex)
        if ( .not. ex ) call error(10201)
      endif

      ! input_weight_file
      input_weight_file = ""
      if ( len(trim(input_prob_file)).eq.0 .and.                       &
           len(trim(ELEMNT(19))).ne.0 ) then
        input_weight_file = ELEMNT(19)
        write(6,'(2x,a)')"* The conformational ensemble is "//         &
                         "re-weighted"
        write(6,'(2x,a)')"* Weight file : "
        write(6,'(6x,a)')trim(input_weight_file)
        inquire(file=trim(input_weight_file),exist=ex)
        if ( .not. ex ) call error(10201)
      endif

      if ( len(trim(input_prob_file)).ne.0 .or.                        &
           len(trim(input_weight_file)).ne.0 ) then
        !! weight way
        if ( ELEMNT(20)(1:6).eq."WEIGHT" .or.                          &
             ELEMNT(20)(1:6).eq."weight" ) then
          weight_method = "WEIGHT"
          !!! weight threshold
          read(ELEMNT(21),*)rtmp
          if ( rtmp.ge.0.d0 .and. rtmp.lt.1.d0 ) weight_threshold = rtmp
          write(6,'(2x,a)')"* weight_threshold for WEIGHT : "
          write(6,'(6x,f)')weight_threshold
        else
          weight_method = "SELECT"
        endif
        write(6,'(2x,a)')"* weight method : "
        write(6,'(6x,a)')weight_method

        !! random_seed
        read(ELEMNT(22),*)rseed
        write(6,'(2x,a)')"* random_seed = "
        write(6,'(6x,i0)')rseed

        !! input_energy_file
        input_energy_file = ELEMNT(23)
        if ( len(trim(input_energy_file)).ne.0 .and.                   &
             len(trim(input_PDB_list)).ne.0 ) then
          write(6,'(2x,a)')"* input_energy_file : "
          write(6,'(6x,a)')trim(input_energy_file)
          inquire(file=trim(input_energy_file),exist=ex)
          if ( .not. ex ) call error(10201)
        endif
        if ( len(trim(input_PDB_list)).ne.0 .and.                      &
             len(trim(input_energy_file)).eq.0 ) call error()
      endif

      !! input parameter for weighting
      if ( len(trim(input_prob_file)) .ne. 0 ) call inweight

      ! periodic boundary condition
      if ( len(trim(ELEMNT(24))) .ne. 0 ) then
        do j = 1,len(trim(ELEMNT(24)))
          YN = ELEMNT(24)(j:j)
          if ( YN.eq."[" .or. YN.eq.":" .or. YN.eq."," .or.            &
               YN.eq."]" )                                             &
            ELEMNT(24) = ELEMNT(24)(:j-1)//" "//ELEMNT(24)(j+1:)
        enddo
        read(ELEMNT(24),*,iostat=iERR)bound(1:2,1:3)
        if ( iERR .ne. 0 ) call error(10207)
        bsize(:) = bound(2,:) - bound(1,:)
        ibsize(:) = 1.0 / bsize(:)
        write(6,'(2x,a)')"* Periodic boundary condition"
        write(6,'(4x,a,f9.3,a3,f9.3)')"x : ",bound(1,1)," - ",bound(2,1)
        write(6,'(4x,a,f9.3,a3,f9.3)')"y : ",bound(1,2)," - ",bound(2,2)
        write(6,'(4x,a,f9.3,a3,f9.3)')"z : ",bound(1,3)," - ",bound(2,3)
        if ( bound(1,1).ge.bound(2,1) .or. bound(1,2).ge.bound(2,2)    &
             .or. bound(1,3).ge.bound(2,3) ) call error(10205)
      else
        bsize(:) = 0.d0 ; ibsize(:) = 0.d0
      endif

!*****

      ! DSSP_flag
      if ( ELEMNT(25)(1:1).eq."Y" .or. ELEMNT(25)(1:1).eq."y" ) then
        DSSP_flag = .true.
      elseif ( ELEMNT(25)(1:1).eq."O" .or. ELEMNT(25)(1:1).eq."o" ) then
        DSSP_flag = .true. ; DSSP_output_flag = .true.
      endif
      if ( DSSP_flag ) then
        call execute_command_line("command -v mkdssp > /dev/null",   &
             exitstat=jj)
        if ( jj .ne. 0 ) then
          write(6,*)"* The DSSP executable file (mkdssp) is not "//  &
            "available, so secondary structure analysis using DSSP"//  &
            " cannot be performed."
          DSSP_flag = .false. ; DSSP_output_flag = .false.
        endif
      endif
      if ( DSSP_flag ) then
        write(6,*)
        write(6,'(2x,a)')"* DSSP execution is performed"
        if ( .not. output_PDB_flag ) then
          d_PDB = trim(d_proj2)//".PDB"
          call system("mkdir "//trim(d_PDB))
          d_PDB = trim(d_PDB)//"/"
        endif

        !! Make directory for DSSP output files
        d_DSSP = trim(d_proj)//"_DSSP"
        inquire(directory=trim(d_DSSP),exist=ex)
        !! The directory has already existed
        if ( ex ) then
          write(6,'(2x,a)')"* ! CAUTION !"
          write(6,'(4x,a)')"DSSP output directory ( "//trim(d_DSSP)//  &
                           ") is OVER-WRITEED"
          call system("rm -f "//trim(d_DSSP)//"/*.dssp > /dev/null")
        !! Make directory for DSSP output files
        else
          write(6,'(2x,a)')"* DSSP files are output in "//             &
                           trim(d_DSSP)//" directory"
          call system("mkdir "//trim(d_DSSP))
        endif
        d_DSSP = trim(d_DSSP)//"/"

        !! Make directory for TOR output files
        d_TOR = trim(d_proj)//"_TOR"
        inquire(directory=trim(d_TOR),exist=ex)
        !! The directory has already existed
        if ( ex ) then
          write(6,'(2x,a)')"* ! CAUTION !"
          write(6,'(4x,a)')"TOR output directory ( "//trim(d_TOR)//    &
                           ") is OVER-WRITEED"
          call system("rm -f "//trim(d_TOR)//"/*.tor > /dev/null")
        !! Make directory for TOR output files
        else
          write(6,'(2x,a)')"* TOR files are output in "//              &
                           trim(d_TOR)//" directory"
          call system("mkdir "//trim(d_TOR))
        endif
        d_TOR = trim(d_TOR)//"/"
      endif

!*****

      ! NMR_flag
      if ( ELEMNT(26)(1:1).eq."Y" .or. ELEMNT(26)(1:1).eq."y" ) then
        NMR_flag = .true.
        write(6,*)
        write(6,'(2x,a)')"* NMR analysis is performed"

        ! input_NOE_file_type
        input_NOE_file_type = ELEMNT(27)
        if ( len(trim(input_NOE_file_type)) .ne. 0 ) then
          write(6,'(2x,a)')"* input_NOE_file_type : "
          write(6,'(6x,a)')trim(input_NOE_file_type)
          if ( input_NOE_file_type .ne. "XPLOR   " .and.               &
               input_NOE_file_type .ne. "DISCOVER" ) call error(10209)

          ! input_NOE_file_name
          input_NOE_file_name = ELEMNT(28)
          if ( len(trim(input_NOE_file_name)) .eq. 0 ) call error(10210)
          write(6,'(2x,a)')"* input_NOE_file_name : "
          write(6,'(6x,a)')trim(input_NOE_file_name)
          inquire(file=trim(input_NOE_file_name),exist=ex)
          if ( .not. ex ) call error(10201)
        endif

        ! output_NOE_file_type
        output_NOE_file_type = ELEMNT(29)
        write(6,'(2x,a)')"* output_NOE_file_type : "
        write(6,'(6x,a)')trim(output_NOE_file_type)
        if ( trim(output_NOE_file_type) .ne. "XPLOR" .and.             &
             trim(output_NOE_file_type) .ne. "DISCOVER" )              &
          call error(10209)

        ! pseudo_Hatom_type_for_output
        read(ELEMNT(30),*)pseudo_Hatom_type_for_output
        if ( pseudo_Hatom_type_for_output .eq. 0 ) then
          write(6,'(2x,a)')"* Pseudo H atom is NOT USED"
        elseif ( pseudo_Hatom_type_for_output .eq. 1 ) then
          write(6,'(2x,a)')"* Pseudo H atom type :"//                  &
               " methyl & aromatic"
        elseif ( pseudo_Hatom_type_for_output .eq. 2 ) then
          write(6,'(2x,a)')"* Pseudo H atom type :"//                  &
               " methyl, aromatic, & prochiral"
        else
          write(6,*)"pseudo_Hatom_type_for_output : ",                 &
                     pseudo_Hatom_type_for_output
          call error(10211)
        endif

        ! input_JCC_file_type
        input_JCC_file_type = ELEMNT(31)
        if ( len(trim(input_JCC_file_type)) .ne. 0 ) then
          write(6,'(2x,a)')"* input_JCC_file_type : "
          write(6,'(6x,a)')trim(input_JCC_file_type)
          if ( trim(input_JCC_file_type) .ne. "XPLOR" .and.       &
               trim(input_JCC_file_type) .ne. "DISCOVER" )        &
            call error(10209)

          ! input_JCC_file_name
          input_JCC_file_name = ELEMNT(32)
          if ( len(trim(input_JCC_file_name)) .eq. 0 ) call error(10212)
          write(6,'(2x,a)')"* input_JCC_file_name : "
          write(6,'(6x,a)')trim(input_JCC_file_name)
          inquire(file=trim(input_JCC_file_name),exist=ex)
          if ( .not. ex ) call error(10201)
        endif

        call mkNOElst ! MaKe NOE LiST for geminal hydrogen atoms
        call mkJCClst ! MaKe JCC LiST
      endif

!*****

      ! Rg calc.
      atom_spec_Rg = ELEMNT(33)
      if ( len(trim(atom_spec_Rg)) .ne. 0 ) then
        write(6,'(2x,a)')"* Radius of gyration calc. is performed"
        write(6,'(4x,a,a)')"atom_spec_Rg : ",trim(atom_spec_Rg)
        allocate(tlist(nATM))
        call atom_specifier(len(trim(atom_spec_Rg)),trim(atom_spec_Rg),&
          nATM,ATMnum,ATMnm,RESnum,RESnm,nCHN,CHN,Rcod,nRg,tlist)
        allocate(iATMrg(nRg)) ; iATMrg(:) = tlist(1:nRg)
        deallocate(tlist)
        call output_specifier_log(size(iATMrg),loglvl,nATM,ATMnum,     &
          ATMnm,RESnum,RESnm,nCHN,CHN,iATMrg)
        call Rg_init()
      endif

!************  Need Reference cord.

      ! RMSD
      RMSD_method = ELEMNT(34) ; RMSD_ref = ELEMNT(35)
      RMSD_fit = ELEMNT(36) ; RMSD_ref_fit = ELEMNT(37)
      if ( len(trim(RMSD_method)) .ne. 0 ) then
        write(6,'(2x,a)')"* RMSD calc. is performed"
        if ( len(trim(RMSD_ref)) .eq. 0 ) RMSD_ref = RMSD_method
        if ( len(trim(RMSD_fit)) .eq. 0 ) RMSD_fit = RMSD_method
        if ( len(trim(RMSD_ref_fit)) .eq. 0 ) RMSD_ref_fit = RMSD_ref

        write(6,'(4x,a)')"atom_spec_RMSD : "//trim(RMSD_method)
        allocate(tlist(nATM))
        call atom_specifier(len(trim(RMSD_method)),trim(RMSD_method),  &
          nATM,ATMnum,ATMnm,RESnum,RESnm,nCHN,CHN,Rcod,nRMSD,tlist)
        allocate(iATMrmsd(nRMSD)) ; iATMrmsd(:) = tlist(1:nRMSD)
        deallocate(tlist)

        if ( RMSD_method .ne. RMSD_ref )                               &
          write(6,'(4x,a)')"atom_spec_RMSD_ref : "//trim(RMSD_ref)
        allocate(tlist(nATMr))
        call atom_specifier(len(trim(RMSD_ref)),trim(RMSD_ref),nATMr,  &
          ATMnumR,ATMnmR,RESnumR,RESnmR,nCHNr,CHNr,refcod,j,tlist)
        if ( j .ne. nRMSD ) call error(10213)
        allocate(iATMrmsdR(nRMSD)) ; iATMrmsdR(:) = tlist(1:nRMSD)
        deallocate(tlist)
        if ( RMSD_method .ne. RMSD_fit )                               &
          write(6,'(4x,a)')"atom_spec_RMSD_fit : "//trim(RMSD_fit)
        allocate(tlist(nATM))
        call atom_specifier(len(trim(RMSD_fit)),trim(RMSD_fit),nATM,   &
          ATMnum,ATMnm,RESnum,RESnm,nCHN,CHN,Rcod,nRMSDf,tlist)
          allocate(iATMrmsdF(nRMSDf)) ; iATMrmsdF(:) = tlist(1:nRMSDf)
          deallocate(tlist)
        if ( RMSD_ref .ne. RMSD_ref_fit )                              &
          write(6,'(4x,a)')"atom_spec_RMSD_ref_fit : "//               &
            trim(RMSD_ref_fit)
        allocate(tlist(nATMr))
        call atom_specifier(len(trim(RMSD_ref_fit)),trim(RMSD_ref_fit),&
          nATMr,ATMnumR,ATMnmR,RESnumR,RESnmR,nCHNr,CHNr,refcod,j,tlist)
        if ( j .ne. nRMSDf ) call error(10213)
        allocate(iATMrmsdRF(nRMSDf)) ; iATMrmsdRF(:) = tlist(1:nRMSDf)
        deallocate(tlist)
        call output_RMSD_log()

        !! RMSD_mode
        RMSD_mode = ELEMNT(38)
        if ( input_reference_PDB_list.ne." " .or.                      &
             input_reference_binary_list.ne." " ) then
          if ( RMSD_mode.eq."AVERAGE" .or.  RMSD_mode.eq."average" ) then
            RMSD_mode = "average" ; ref_conf_1_flag = .true.
            write(6,'(4x,a)')"RMSD_mode : average"
          else
            RMSD_mode = "multi"
            write(6,'(4x,a)')"RMSD_mode : multi"
          endif
        endif
        !! check rotaion flag
        rotation_flag = .true.

        ! interval for timeseries RMSD
        read(ELEMNT(39),*)interval_for_timeseries_RMSD
        if ( interval_for_timeseries_RMSD .gt. 0 ) then
          write(6,'(2x,a)')"* interval for timeseries RMSD : "
          write(6,'(6x,i0)')interval_for_timeseries_RMSD
        else
          interval_for_timeseries_RMSD = 0
        endif

        ! Assign RcodRMSD
        if ( len(trim(input_reference_binary_list)) .ne. 0 ) then
          call inrefbinL
        elseif ( len(trim(input_reference_PDB_list)) .ne. 0 ) then
          call inrefpdbL
        endif
        call RMSD_init()
      endif

!*****

      ! PCA
      atom_spec_PCA = ELEMNT(40)
      if ( len(trim(atom_spec_PCA)) .ne. 0 ) then
        write(6,'(2x,a)')"* Making PCA coordinates is performed"
        write(6,'(4x,a)')"atom_spec_PCA : ",trim(atom_spec_PCA)
        allocate(tlist(nATM))
        call atom_specifier(len(trim(atom_spec_PCA)),                  &
          trim(atom_spec_PCA),nATM,ATMnum,ATMnm,RESnum,RESnm,nCHN,CHN, &
          Rcod,nPCA,tlist)
        allocate(iATMPCA(nPCA)) ; iATMPCA(:) = tlist(1:nPCA)
        deallocate(tlist)
        call output_specifier_log(nPCA,loglvl,nATM,ATMnum,             &
          ATMnm,RESnum,RESnm,nCHN,CHN,iATMPCA)

        PCA_method = ELEMNT(41)
        write(6,'(4x,a)')"PCA_method :"
        if ( PCA_method(1:3).eq."ave" .or.                           &
             PCA_method(1:3).eq."AVE" ) then
          PCA_method = "ave"
          write(6,'(6x,a)')"average distance among residues"
        elseif ( PCA_method(1:3).eq."max" .or.                       &
                 PCA_method(1:3).eq."MAX" ) then
          PCA_method = "max"
          write(6,'(6x,a)')"maximum distance among residues"
        elseif ( PCA_method(1:3).eq."min" .or.                       &
                 PCA_method(1:3).eq."MIN" ) then
          PCA_method = "min"
          write(6,'(6x,a)')"minimum distance among residues"
        elseif ( PCA_method(1:7).eq."contact" .or.                   &
                 PCA_method(1:7).eq."CONTACT" ) then
          PCA_method = "contact"
          write(6,'(6x,a)')"contact among residues"//                &
                           " (0: no contact, 1:contact)"
        elseif ( PCA_method(1:9).eq."intercord" .or.                 &
                 PCA_method(1:9).eq."INTERCORD" ) then
          PCA_method = "intercord"
          write(6,'(6x,a)')"Inter coordinats of the selected atoms"
        else
          PCA_method = "cord"
          write(6,'(6x,a)')"3D coordinates of the selected atoms"
        endif

        if ( PCA_method.ne."cord" .or. PCA_method.ne."intercord" ) then
          !! neighbor_residue_PCA
          read(ELEMNT(42),*)cut_resnum_PCA
          if ( cut_resnum_PCA .lt. 0 ) cut_resnum_PCA = 0
          write(6,'(2x,a,i0)')"* Neighbor residue number = ",     &
                              cut_resnum_PCA

          !! tolerance_PCA
          if ( PCA_method .eq. "contact" ) then
            allocate(radPCA(nPCA))
            read(ELEMNT(43),*)rtmp
            if ( len(trim(input_topology)) .ne. 0 ) then
              if ( rtmp .gt. 0.d0 ) then
                tolePCA = rtmp
              else
                tolePCA = 2.8
              endif
              write(6,'(2x,a,f8.3,a)')                                 &
                "* Tolerance distance for contact judge = ",tolePCA,   &
                " (A)"
              f = .true.
              call mk_vdWrad_list(nATM,nPCA,iATMPCA,RESnm,ATMnm,radPCA,&
                f)
            else
              if ( rtmp .gt. 0.d0 ) then
                tolePCA = rtmp
              else
                tolePCA = 6.0
              endif
              write(6,'(2x,a,f8.3,a)')"* Contact distance = ",tolePCA, &
                " (A)"
              radPCA(:) = tolePCA*0.5d0 ; tolePCA = 0.d0
            endif
          endif
        endif
        call PCAprep() ; call PCA_init()
      endif

!*****

      !! Qvalue
      atom_spec_Qvalue = ELEMNT(44)
      if ( len(trim(atom_spec_Qvalue)) .ne. 0 ) then
        if ( len(trim(reference_PDB)) .eq. 0 ) call error(10217)
        if ( len(trim(input_reference_PDB_list)) .eq. 0 )              &
          call error(10203)
        write(6,'(2x,a)')"* Qvalue analysis is performed"
        write(6,'(4x,a)')"Atom specifier for Qvalue : "//              &
          trim(atom_spec_Qvalue)
        allocate(tlist(nATM))
        if ( nATM.ne.nATMr .or. .not. all(ATMnum.eq.ATMnumR) .or.      &
             .not. all(ATMnm.eq.ATMnmR) .or.                           &
             .not. all(RESnum.eq.RESnumR) .or.                         &
             .not. all(RESnm.eq.RESnmR) .or.                           &
             .not. all(nCHN.eq.nCHNr) .or.                             &
             .not. all(CHN.eq.CHNr)) call error(10218)
        call atom_specifier(len(trim(atom_spec_Qvalue)),               &
          trim(atom_spec_Qvalue),nATM,ATMnum,ATMnm,RESnum,RESnm,nCHN,  &
          CHN,refcod,nQ,tlist)
        allocate(iATMQ(nQ)) ; iATMQ(:) = tlist(1:nQ)
        deallocate(tlist)
        call output_specifier_log(nQ,loglvl,nATM,ATMnum,               &
          ATMnm,RESnum,RESnm,nCHN,CHN,iATMQ)

        !! neighbor_residue_Qvalue
        read(ELEMNT(45),*)cut_resnum_Qvalue
        if ( cut_resnum_Qvalue .lt. 0 ) cut_resnum_Qvalue = 0
        write(6,'(2x,a,i0)')"* Neighbor residue number = ",            &
                            cut_resnum_Qvalue

        !! tolerance_Q
        read(ELEMNT(46),*)rtmp
        if ( len(trim(input_topology)) .ne. 0 ) then
          if ( rtmp .gt. 0.d0 ) then
            toleQ = rtmp
          else
            toleQ = 2.8
          endif
          write(6,'(2x,a,f8.3,a)')                                     &
            "* Tolerant distance for Qvalue contact = ",toleQ," (A)"
          write(6,'(4x,a)')"Radius info. in the input tpl file is used"
        else
          if ( rtmp .gt. 0.d0 ) then
            toleQ = rtmp
          else
            toleQ = 6.0
          endif
          write(6,'(2x,a,f8.3,a)')"* Contact distance for Qvalue ",    &
            toleQ," (A)"
        endif

        !! tolerance_NCP
        read(ELEMNT(47),*)rtmp
        if ( len(trim(input_topology)) .ne. 0 ) then
          if ( rtmp .gt. 0.d0 ) then
            toleNCP = rtmp
          else
            toleNCP = 2.8
          endif
          write(6,'(2x,a,f8.3,a)')                                     &
            "* Tolerant distance for NCP contact = ",toleNCP," (A)"
          write(6,'(4x,a)')"Radius info. in the input tpl file is used"
        else
          if ( rtmp .gt. 0.d0 ) then
            toleNCP = rtmp
          else
            toleNCP = 6.0
          endif
          write(6,'(2x,a,f8.3,a)')"* Contact distance for NCP ",      &
            toleNCP," (A)"
        endif

        !! NCP_rate
        read(ELEMNT(48),*)NCP_rate
        if ( NCP_rate.le.0.0 .or. NCP_rate.lt.1.0 ) NCP_rate = 0.7
        write(6,'(2x,a,f8.3)')                                         &
          "* Threshold ratio for NCP determination = ",NCP_rate

        !! trajectory based Qvalue flag
        if ( ELEMNT(49)(1:1).eq."Y" .or. ELEMNT(49)(1:1).eq."y" ) then
          weakness_NCP_flag = .true.
          write(6,'(2x,a)')"* Weakness native contact pairs are "//    &
                           "searched"
        endif
        call Qprep()
      endif

!*****

      ! cell boundary condition
      if ( len(trim(ELEMNT(50))) .ne. 0 ) then
        do j = 1,len(trim(ELEMNT(50)))
          YN = ELEMNT(50)(j:j)
          if ( YN.eq."[" .or. YN.eq.":" .or. YN.eq."," .or.            &
               YN.eq."]" )                                             &
            ELEMNT(50) = ELEMNT(50)(:j-1)//" "//ELEMNT(50)(j+1:)
        enddo
        read(ELEMNT(50),*,iostat=iERR)cbound(1:2,1:3)
        if ( iERR .ne. 0 ) call error(10207)
        csize(:) = cbound(2,:) - cbound(1,:)
        icsize(:) = 1.0 / csize(:)
        write(6,'(2x,a)')"* Cell boundary for distrib or pocket search"&
          //" analyses"
        write(6,'(4x,a,f9.3,a3,f9.3)')"x : ",cbound(1,1)," - ",        &
          cbound(2,1)
        write(6,'(4x,a,f9.3,a3,f9.3)')"y : ",cbound(1,2)," - ",        &
          cbound(2,2)
        write(6,'(4x,a,f9.3,a3,f9.3)')"z : ",cbound(1,3)," - ",        &
          cbound(2,3)
        if ( cbound(1,1).ge.cbound(2,1) .or. cbound(1,2).ge.cbound(2,2)&
             .or. cbound(1,3).ge.cbound(2,3) ) call error(10205)
      else
        csize(:) = 0.d0 ; icsize(:) = 0.d0
      endif

      ! distribution analysis
      atom_spec_distrib = ELEMNT(51)
      if ( len(trim(atom_spec_distrib)) .ne. 0 ) then
        write(6,'(2x,a)')"* Distribution analysis is performed"
        if ( csize(1) .eq. 0.d0 ) call error(10219)
        write(6,'(4x,a)')"atom_spec_distrib : "//trim(atom_spec_distrib)

        !! distribution cell size
        read(ELEMNT(52),*)rtmp
        if ( rtmp .gt. 0.d0 ) d_celsiz = rtmp
        write(6,'(4x,a,f8.3,a)')"Cell size = ",d_celsiz,"(A)"
        call distrib_init()
      endif

      ! contact analysis
      atom_spec_contact = ELEMNT(53)
      if ( len(trim(atom_spec_contact)) .ne. 0 ) then
        j = index(atom_spec_contact,"_")
        if ( j .eq. 0 ) call error(10215)
        tmp1 = atom_spec_contact(1:j-1) ; tmp2 = atom_spec_contact(j+1:)
        write(6,'(2x,a)')"* Contact analysis is performed"
        write(6,'(4x,a)')"Probe group  : "//trim(tmp1)
        allocate(tlist(nATM))
        call atom_specifier(len(trim(tmp1)),trim(tmp1),nATM,ATMnum,    &
          ATMnm,RESnum,RESnm,nCHN,CHN,Rcod,nprob,tlist)
        allocate(iATMprob(nprob)) ; iATMprob(:) = tlist(1:nprob)
        call output_specifier_log(nprob,loglvl,nATM,ATMnum,     &
          ATMnm,RESnum,RESnm,nCHN,CHN,iATMprob)
        deallocate(tlist) ; allocate(radprob(nprob)) ; f = .false.
        if ( len(trim(input_topology)) .ne. 0 ) then
          call mk_vdWrad_list(nATM,nprob,iATMprob,RESnm,ATMnm,radprob,f)
        else
          call mk_vdWrad_list2(nATM,nprob,iATMprob,RESnm,ATMnm,radprob,&
                               f)
        endif
        write(6,'(4x,a)')"Target group : "//trim(tmp2)
        allocate(tlist(nATM))
        call atom_specifier(len(trim(tmp2)),trim(tmp2),nATM,ATMnum,    &
          ATMnm,RESnum,RESnm,nCHN,CHN,Rcod,ntgt,tlist)
        allocate(iATMtgt(ntgt)) ; iATMtgt(:) = tlist(1:ntgt)
        call output_specifier_log(ntgt,loglvl,nATM,ATMnum,     &
          ATMnm,RESnum,RESnm,nCHN,CHN,iATMtgt)
        deallocate(tlist) ; allocate(radtgt(ntgt)) ; f = .true.
        if ( len(trim(input_topology)) .ne. 0 ) then
          call mk_vdWrad_list(nATM,ntgt,iATMtgt,RESnm,ATMnm,radtgt,f)
        else
          call mk_vdWrad_list2(nATM,ntgt,iATMtgt,RESnm,ATMnm,radtgt,f)
        endif

        read(ELEMNT(54),*)rtmp
        if ( rtmp .gt. 0.d0 ) watrad = rtmp
        write(6,'(4x,a,f8.3,a)')"* Tolerance distance for contact = ", &
          watrad," (A)"
        watrad = watrad * 0.5d0
        if ( ELEMNT(55)(1:1).eq."Y" .or. ELEMNT(55)(1:1).eq."y" ) then
          contact_correlation_flag = .true.
          write(6,'(4x,a)')"Contact_correlation is calculated"
        endif
        call contact_init()
      endif

      ! SOAP calc.
      read(ELEMNT(56),*)nSOAP
      if ( nSOAP.gt.0 .and. nSOAP.lt.fnlRES-fstRES+1 ) then
        if ( mod(nSOAP,2) .eq. 0 ) then
          nSOAP = nSOAP + 1
          write(6,'(4x,a)')"+ window_for_SOAP gets + 1"
          write(6,'(4x,a)')"  because it is an even number"
        endif
        write(6,'(2x,a)')"* SOAP calculation is performed"
        write(6,'(4x,a,i0)')"Window for SOAP = ",nSOAP
      else
        nSOAP = 0
      endif

      ! measurements
      tmp = trim(ELEMNT(57))
      if ( len(trim(tmp)) .ne. 0 ) then
        write(6,'(2x,a)')"* measurements : "//trim(tmp)
        call measure(tmp)
      endif

      ! RMSF calc.
      if ( ELEMNT(58)(1:1).eq."Y" .or. ELEMNT(58)(1:1).eq."y" ) then
        RMSF_flag = .true.
        write(6,'(2x,a)')"* RMSF of each atoms from the average "//    &
          "structure is output in ..."
        write(6,'(8x,a)')trim(d_proj)//"_RMSF.pdb"
      endif

      ! pocket search analysis
      atom_spec_pocket = ELEMNT(59)
      if ( len(trim(atom_spec_pocket)) .ne. 0 ) then
        write(6,'(2x,a)')"* Pocket search analysis is performed"
        if ( csize(1) .eq. 0.d0 ) call error(10219)
        write(6,'(4x,a)')"atom_spec_pocket : "//trim(atom_spec_pocket)
        !! pocket search cell size
        read(ELEMNT(60),*)rtmp
        if ( rtmp .gt. 0.d0 ) p_celsiz = rtmp
        write(6,'(4x,a,f8.3,a)')"Cell size = ",p_celsiz,"(A)"
        !! positive method flag
        if ( ELEMNT(61)(1:1).eq."Y" .or. ELEMNT(61)(1:1).eq."y" ) then
          write(6,'(2x,a)')"* Positive method analysis is performed"
          positive_method_flag = .true.
          read(ELEMNT(62),*)watrad2
          write(6,'(4x,a,f8.3)')"Water radius for pocket search = ",   &
            watrad2
          watrad2 = watrad2 * 2.d0
        endif
        read(ELEMNT(63),*)atom_number_for_pocket_entrance
        if ( atom_number_for_pocket_entrance.lt.1 .or.                 &
             atom_number_for_pocket_entrance.gt.nATM )                 &
             atom_number_for_pocket_entrance = 0
        if ( atom_number_for_pocket_entrance .ne. 0 )                  &
          write(6,'(4x,a,i0)')"atom number for pocket entrance = ",    &
            atom_number_for_pocket_entrance
        call pocket_search_init()
      endif

!************************************

      return
      end subroutine store_element

