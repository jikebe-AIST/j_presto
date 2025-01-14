
      subroutine error(iERR)

!******************************************
!
!     Output error condition
!
!******************************************

      implicit none

      integer(4),intent(in):: iERR

!*********************************

      write(6,*)" "
      write(6,*)" !! ERROR !!"
      write(6,*)" ERROR number = ",iERR
      write(6,*)" "

!*****************************

      ! Main Phase 
      if ( iERR .lt. 10000 ) then
        write(6,*)"phase : whole"

        ! *** main_dssp
        write(6,*)"phase : main"

      ! Input Phase
      elseif ( iERR .lt. 20000 ) then
        write(6,*)"phase : input"

        ! *** asnelm
        if ( iERR .lt. 10200 ) then
          write(6,*)"place : asnelm.f90"

          if ( iERR .eq. 10101 ) then
            write(6,*)"reason : Input file style is strange"
            call system("cat help.txt")
          endif

        ! *** store_element
        elseif ( iERR .lt. 10300 ) then
          write(6,*)"place : store_element.f90"

          if ( iERR .eq. 10201 ) then
            write(6,*)"reason : The above file do not exist"
          elseif ( iERR .eq. 10202 ) then
            write(6,*)"reason : Please input template_PDB"
          elseif ( iERR .eq. 10203 ) then
            write(6,*)"reason : When Qvalue_flag = Y, you should set"//&
                      " input_reference_PDB_list option"
          elseif ( iERR .eq. 10204 ) then
            write(6,*)"reason : You can set only one of 3 input types"
            write(6,*)"         pdb, binary, or otherfmt"
          elseif ( iERR .eq. 10205 ) then
            write(6,*)"reason : The above value(s) you set"//          &
                      " are inadequate"
          elseif ( iERR .eq. 10207 ) then
            write(6,*)"reason : Input boundary condition is strange"
          elseif ( iERR .eq. 10208 ) then
            write(6,*)"reason : please set the above value(s)"
          elseif ( iERR .eq. 10209 ) then
            write(6,*)"reason : You can set it only to ..."
            write(6,*)"         1 : XPLOR"
            write(6,*)"         2 : DISCOVER"
          elseif ( iERR .eq. 10210 ) then
            write(6,*)"reason : Please input input_NOE_file_name"
          elseif ( iERR .eq. 10211 ) then
            write(6,*)"reason : You can set pseudo_Hatom_type_for_"//  &
                      "output only to ..."
            write(6,*)"         1 : methyl & aromatic(low) (default)"
            write(6,*)"         2 : methyl, ethylen & aromatic(high)"
          elseif ( iERR .eq. 10212 ) then
            write(6,*)"reason : Please input input_JCC_file_name"
          elseif ( iERR .eq. 10213 ) then
            write(6,*)"reason : Discrepancies in correspondence between"
            write(6,*)"         selected atoms in the template PDB and "
            write(6,*)"         those in the reference PDB"
          elseif ( iERR .eq. 10214 ) then
            write(6,*)"reason : When you input input_PDB_file and "
            write(6,*)"         input_prob_file, you also must input "
            write(6,*)"         input_energy_file"
          elseif ( iERR .eq. 10215 ) then
            write(6,*)"reason : The atom_spec must be gave as "
            write(6,*)"         [probe group]_[target group]"
          elseif ( iERR .eq. 10216 ) then
            write(6,*)"reason : The above pair is overwrapped with "//&
                      "the other pairs defined previously"
          elseif ( iERR .eq. 10217 ) then
            write(6,*)"reason : When Qvalue_flag = Y, you should set"//&
                      " reference_PDB option"
          elseif ( iERR .eq. 10218 ) then
            write(6,*)"reason : When calculating Qvalue, reference_PDB"
            write(6,*)"         must specify the same system as"
            write(6,*)"         template_PDB"
          elseif ( iERR .eq. 10219 ) then
            write(6,*)"reason : You must set cell_boundary option"
          elseif ( iERR .eq. 10220 ) then
            write(6,*)"reason : The above atom list is inadequate "
          elseif ( iERR .eq. 10220 ) then
            write(6,*)"reason : The above atom list is inadequate "
          endif

        ! *** intmp
        elseif ( iERR .lt. 10400 ) then
          write(6,*)"place : intmp.f90"

          if ( iERR .eq. 10301 ) then
            write(6,*)"reason : In input template PDB, different "//   &
             "residues with the same residue number are listed."
          endif

        ! *** pdbchk
        elseif ( iERR .lt. 10500 ) then
          write(6,*)"place : pdb_chk.f90"

          if ( iERR .eq. 10401 ) then
            write(6,*)"reason : The above file is NOT consistent "// &
                      "with the template (or reference) PDB"
          endif

        ! *** mk_pdb
        elseif ( iERR .lt. 10600 ) then
          write(6,*)"place : mk_pdb.f90"

          if ( iERR .eq. 10501 ) then
            write(6,*)"reason : The above BINARY file do not exist"
          elseif ( iERR .eq. 10502 ) then
            write(6,*)"reason : The above conformation number in  "
            write(6,*)"         BINARY file list is strange"
          endif
          
        ! *** mk_vdWrad_list
        elseif ( iERR .lt. 10700 ) then
          write(6,*)"place : mk_vdWrad_list"

          if ( iERR .eq. 10601 ) then
            write(6,*)"reason : Could not find information with the"//&
             " residue and atom names"
          endif

        ! *** pocket_search_init
        elseif ( iERR .lt. 10800 ) then
          write(6,*)"place : pocket_search_init"

          if ( iERR .eq. 10701 ) then
            write(6,*)"reason : The atom spec must be gave as"
            write(6,*)"         [receptor]_[ligand]"
          endif

        endif

      ! NMR prepare phase
      elseif ( iERR .lt. 30000 ) then
        write(6,*)"phase : NOE file making"

        if ( iERR .lt. 20200 ) then
          write(6,*)"place : mkNOElst.f90"

        elseif ( iERR .lt. 20300 ) then
          write(6,*)"place : Xplor.f90"

        elseif ( iERR .lt. 20400 ) then
          write(6,*)"place : DISCOVER.f90"

        elseif ( iERR .lt. 20500 ) then
          write(6,*)"place : assign_noe.f90"

          if ( iERR .eq. 20401 ) then
            write(6,*)"reason : The above exp. NOE pair don't exist in"
            write(6,*)"         the reference PDB"
          endif

        endif

      ! Analysis phase
      elseif ( iERR .lt. 40000 ) then
        write(6,*)"phase : Analysis"

        ! *** inpdb
        if ( iERR .lt. 30200 ) then
          write(6,*)"place : inpdb.f90"

        elseif ( iERR .lt. 30300 ) then
          write(6,*)"place : inbin.f90"

        endif

      ! Execution phase
      elseif ( iERR .lt. 50000 ) then
        write(6,*)"phase : DSSP execution"

        ! *** dssp_exe
        if ( iERR .lt. 50200 ) then
          write(6,*)"place : dssp_exe.f90"
        endif

      ! Analysis phase
      elseif ( iERR .lt. 60000 ) then
        write(6,*)"phase : DSSP analysis"

        ! *** normal_SSC
        if ( iERR .lt. 50200 ) then
          write(6,*)"place : normal_SSC.f90"

          if ( iERR .eq. 50101 ) then
            write(6,*)"reason : The above DSSP file is strange"
            write(6,*)"         Please check"
          elseif ( iERR .eq. 50102 ) then
            write(6,*)"reason : The above DSSP file do not exist"
          endif

        ! *** weight_SSC
        elseif ( iERR .lt. 50300 ) then
          write(6,*)"place : reweight_analy.f90"

          if ( iERR .eq. 50201 ) then
            write(6,*)"reason : The above DSSP file is strange"
            write(6,*)"         Please check"
          elseif ( iERR .eq. 50202 ) then
            write(6,*)"reason : The above DSSP file do not exist"
          endif

        endif

      ! RMSD calc.
      elseif (iERR .lt. 70000 ) then
        write(6,*)"phase : RMSD calc"

        ! *** preRMSD
        if ( iERR .lt. 60200 ) then
          write(6,*)"place : preRMSD.f90"

        endif

      endif

!*********************

      write(6,*)""

      stop
      end subroutine error
