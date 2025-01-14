
      module COMIFN

!**********************************************
!
!     COMmon data space for Input File Names
!
!**********************************************

      implicit none

      logical(4)::     overwrite_flag = .false.
      character(130):: project_name
      character(130):: template_PDB
      character(130):: input_topology
      character(130):: input_PDB_list
      character(130):: input_binary_list
      character(130):: input_otherfmt_list
      character(130):: input_restart_list
      logical(4)::     output_PDB_flag = .false.
      logical(4)::     output_binary_flag = .false.
      character(130):: reference_PDB = " "
      character(130):: input_reference_PDB_list
      character(130):: input_reference_binary_list
      character(999):: atom_spec_output
      integer(4)::     interval
      character(130):: input_prob_file
      character(130):: input_weight_file
      character(6)::   weight_method
      real(8)::        weight_threshold = 0.d0
      integer(4)::     rseed
      character(130):: input_energy_file
      real(4)::        bsize(3)
      real(4)::        ibsize(3)
      real(4)::        bound(2,3)
      logical(4)::     DSSP_flag = .false.
      logical(4)::     DSSP_output_flag = .false.
      logical(4)::     NMR_flag = .false.
      character(10)::  input_NOE_file_type
      character(130):: input_NOE_file_name
      character(10)::  output_NOE_file_type
      integer(4)::     pseudo_Hatom_type_for_output
      character(8)::   input_JCC_file_type
      character(130):: input_JCC_file_name
      character(999):: atom_spec_Rg
      character(999):: RMSD_method
      character(999):: RMSD_ref
      character(999):: RMSD_fit
      character(999):: RMSD_ref_fit
      character(8)::   RMSD_mode
      character(999):: atom_spec_PCA
      character(10)::  PCA_method
      integer(4)::     cut_resnum_PCA
      real(4)::        tolePCA
      character(999):: atom_spec_Qvalue
      integer(4)::     cut_resnum_Qvalue
      real(4)::        toleQ
      real(4)::        toleNCP
      real(4)::        NCP_rate
      logical(4)::     weakness_NCP_flag = .false.
      character(999):: atom_spec_contact
      integer(4)::     torecon
      logical(4)::     contact_correlation_flag = .false.
      real(4)::        csize(3)
      real(4)::        icsize(3)
      real(4)::        cbound(2,3)
      character(999):: atom_spec_distrib
      real(8)::        d_celsiz
      integer(4)::     nSOAP
      logical(4)::     RMSF_flag = .false.
      character(999):: atom_spec_pocket
      real(8)::        p_celsiz
      logical(4)::     positive_method_flag = .false.
      real(8)::        watrad2
      integer(4)::     atom_number_for_pocket_entrance
      character(1)::   loglvl

!**********************************

      end module COMIFN

