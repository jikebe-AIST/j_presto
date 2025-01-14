
      module COMINP

!*********************************
!
!     COMmon data space for INPut
!
!*********************************

      implicit none

      ! Max number of elements
        integer(4),parameter:: nELM = 63
      ! Control parameter for calc.
        character(999),save:: ELEMNT(nELM) = (/             &
        "N",                                                & ! overwrite_flag
        "a",                                                & ! project_name
        "n",                                                & ! log_level
        " ",                                                & ! template_PDB
        " ",                                                & ! input_topology
        " ",                                                & ! input_PDB_list
        " ",                                                & ! input_binary_list
        " ",                                                & ! input_otherfmt_list
        " ",                                                & ! input_restart_list
        "N",                                                & ! output_PDB_flag
        "N",                                                & ! output_binary_flag
        " ",                                                & ! reference_PDB
        " ",                                                & ! input_reference_PDB_list
        " ",                                                & ! input_reference_binary_list
        "Y",                                                & ! PDB_consistency_check
        "*",                                                & ! atom_spec_output
        "1",                                                & ! interval
        " ",                                                & ! input_prob_file
        " ",                                                & ! input_weight_file
        "SELECT",                                           & ! weight_method
        "0.00001",                                          & ! weight_threshold
        "999999",                                           & ! random_seed
        " ",                                                & ! input_energy_file
        " ",                                                & ! periodic boundary
        "N",                                                & ! DSSP_flag
        "N",                                                & ! NMR_flag
        "XPLOR",                                            & ! input_NOE_file_type
        " ",                                                & ! input_NOE_file_name
        " ",                                                & ! output_NOE_file_type
        "2",                                                & ! pseudo_Hatom_type_for_output
        " ",                                                & ! input_JCC_file_type
        " ",                                                & ! input_JCC_file_name
        " ",                                                & ! atom_spec_Rg
        " ",                                                & ! atom_spec_RMSD
        " ",                                                & ! atom_spec_RMSD_ref
        " ",                                                & ! atom_spec_RMSD_fit
        " ",                                                & ! atom_spec_RMSD_reffit
        "multi",                                            & ! RMSD_mode
        "0",                                                & ! interval_for_timeseries_RMSD
        " ",                                                & ! atom_spec_PCA
        "min",                                              & ! PCA_method
        "1",                                                & ! neighbor_residue_PCA
        "0.d0",                                             & ! tolerance_PCA
        " ",                                                & ! atom_spec_Qvalue
        "1",                                                & ! neighbor_residue_Qvalue
        "0.d0",                                             & ! tolerance_Qvalue
        "0.d0",                                             & ! tolerance_NCP
        "0.7",                                              & ! NCP_rate
        "N",                                                & ! weakness_NCP_flag
        " ",                                                & ! cell_boundary
        " ",                                                & ! atom_spec_distribution
        "3.d0",                                             & ! distribution_cell_size
        " ",                                                & ! atom_spec_contact
        "2.8d0",                                            & ! tolerance_contact
        "N",                                                & ! contact_correlation_flag
        "0",                                                & ! window_for_SOAP
        " ",                                                & ! measurements
        "N",                                                & ! RMSF_flag
        " ",                                                & ! atom_spec_pocket
        "3.d0",                                             & ! pocket_cell_size
        "N",                                                & ! positive_method_flag
        "1.4d0",                                            & ! water_radius_for_pocket
        "0"                                                 & ! atom_number_for_poket_entrance
        /)

      ! Control data flag NAMe of ELeMents
        character(80),save:: namELM(nELM) = (/              &
        "overwrite_flag",                                   &
        "project_name",                                     &
        "log_level",                                        &
        "template_PDB",                                     &
        "input_topology",                                   &
        "input_PDB_list",                                   &
        "input_binary_list",                                &
        "input_otherfmt_list",                              &
        "input_restart_list",                               &
        "output_PDB_flag",                                  &
        "output_binary_flag",                               &
        "reference_PDB",                                    &
        "input_reference_PDB_list",                         &
        "input_reference_binary_list",                      &
        "PDB_consistency_check",                            &
        "atom_spec_output",                                 &
        "interval",                                         &
        "input_prob_file",                                  &
        "input_weight_file",                                &
        "weight_method",                                    &
        "weight_threshold",                                 &
        "random_seed",                                      &
        "input_energy_file",                                &
        "periodic_boundary",                                &
        "DSSP_flag",                                        &
        "NMR_flag",                                         &
        "input_NOE_file_type",                              &
        "input_NOE_file_name",                              &
        "output_NOE_file_type",                             &
        "pseudo_Hatom_type_for_output",                     &
        "input_JCC_file_type",                              &
        "input_JCC_file_name",                              &
        "atom_spec_Rg",                                     &
        "atom_spec_RMSD",                                   &
        "atom_spec_RMSD_ref",                               &
        "atom_spec_RMSD_fit",                               &
        "atom_spec_RMSD_reffit",                            &
        "RMSD_mode",                                        &
        "interval_for_timeseries_RMSD",                     &
        "atom_spec_PCA",                                    &
        "PCA_method",                                       &
        "neighbor_residue_PCA",                             &
        "tolerance_PCA",                                    &
        "atom_spec_Qvalue",                                 &
        "neighbor_residue_Qvalue",                          &
        "tolerance_Qvalue",                                 &
        "tolerance_NCP",                                    &
        "NCP_rate",                                         &
        "weakness_NCP_flag",                                &
        "cell_boundary",                                    &
        "atom_spec_distribution",                           &
        "distribution_cell_size",                           &
        "atom_spec_contact",                                &
        "tolerance_contact",                                &
        "contact_correlation_flag",                         &
        "window_for_SOAP",                                  &
        "measurements",                                     &
        "RMSF_flag",                                        &
        "atom_spec_pocket",                                 &
        "pocket_cell_size",                                 &
        "positive_method_flag",                             &
        "water_radius_for_pocket",                          &
        "atom_number_for_pocket_entrance"                   &
        /)

      end module COMINP

