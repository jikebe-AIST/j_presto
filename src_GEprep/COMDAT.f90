!
!     module COMINP
!     module COMIFN
!     module COMVAL
!
!=======================================================================


      module COMINP

!*********************************
!
!     COMmon data space for INPut
!
!*********************************

      implicit none

      ! Max number of elements
        integer(4),parameter:: nELM = 33

      ! Control parameter for calc.
        character(130),save:: ELEMNT(nELM) = (/                   &
        ! Project name
                    "test",                                       &
        ! Simulation type
                    "",                                           &
        ! ENErgy data LiST name
                    "",                                           &
        ! Min & Max reaction cord. range in the current MD
                    "","",                                        &
        ! Temperature & array dimension for fitting
                    "","4","1",                                   &
        ! Min & Max reaction cord. range in the NeXt MD
                    "","",                                        &
        ! previous dlnN data file name
                    "",                                           &
        ! TeMPerture for SaMPling
                    "","",                                        &
        ! lambda and the square for reweighting
                    "1.d0","1.d0",                                &
        ! min & max values for reweighting
                    "","","","","","",                            &
        ! acceleration & the rate
                    "N","1.d0","N","1.d0","N","1.d0","1.5d0",     &
        ! flag of neglect out of the range
                    "N",                                          &
        ! force size for out of the range for ALSD
                    "1000.d0","","",                              &
        ! ignore nodata bin flag
                    "N"/)

      ! Control data flag NAMe of ELeMents
        character(6),save:: namELM(nELM) = (/                     &
        ! Project name
                     "PROJNM",                                    &
        ! Method
                     "SIMTYP",                                    &
        ! DATa LiST name
                     "CODLST",                                    &
        ! Min & Max reaction cord. range in the current MD
                     "MINVAL","MAXVAL",                           &
        ! Temperature & array dimension for fitting
                     "SIMTEM","FITDIM","NWINDO",                  &
        ! Min & Max reaction cord. range in the NeXt MD
                     "NXMINV","NXMAXV",                           &
        ! PREvious FITting data file name
                     "PREFIT",                                    &
        ! TeMPerture for SaMPling & NeXt simulation Temperature
                     "TMPSMP","NXTEMP",                           &
        ! lambda for reweighting
                     "LAMBDA","LMDSQU",                           &
        ! min & max values for reweighting
                     "MINLMD","MAXLMD","MINEAA","MAXEAA",         &
                     "MINEAB","MAXEAB",                           &
        ! acceleration & the rate
                     "ACCELE","ACCRAT","ACCEL2","ACCRT2","ACCLAB",&
                     "ACCRTL","ACCTHL",                           &
        ! flag of neglect out of the range
                     "NEGLCT",                                    &
        ! force size for out of the range for ALSD
                     "FORCES","FORCEL","FORCEH",                  &
        ! ignore nodata bin flag
                     "IGNDAT"/)

      end module COMINP


!======================================================================


      module COMIFN

!**********************************************
!
!     COMmon data space for Input File Names
!
!**********************************************

      implicit none

      ! PROJect NaMe for data output
        character(130):: PROJNM
      ! DATa LiST name
        character(130):: DATLST
      ! Min & Max Energy range in the current MD
        real(8):: minE,maxE
      ! Method
        character(5):: METHOD
      ! CuRrent Temperature
        real(8):: T
      ! array information for fitting for dlnN
        integer(4):: NfitDIM,Nwindow
      ! Min & Max Energy range in the NeXt MD
        real(8):: nxminE,nxmaxE
      ! PREvious FITting data file name
        character(130):: PREFIT
      ! TeMPerture for SaMPling & Next simulation temperature
        real(8):: smpT,nxT
      ! min & max values for reweighting
        real(8):: minlmd,maxlmd,minEAA,maxEAA,minEAB,maxEAB
      ! Number of bins for reweighting
        integer(4):: nlmd,nEAA,nEAB
      ! Bin size for reweighting
        real(8):: slmd,sEAA,sEAB
      ! acceleration & the rate
        logical(4):: accele,accel2,acclab
        real(8):: accrat,accrt2,accrtl,accthl
      ! lambda for reweighting
        real(8):: lmd_rewei,lmd_rewei2
      ! flag of neglect out of the range
        logical(4):: FneglectH,FneglectL
      ! force size for out of the range for ALSD
        real(8):: forces,forcel,forceh
      ! ignore nodata bins for the dlnN fitting
        logical(4):: igndat

!**********************************

      end module COMIFN


!====================================================================


      module COMVAL

!**************************************************
!
!     COMmon data space for some VALues
!
!**************************************************

      implicit none

      ! Number of BINs
        integer(4),parameter:: nBIN = 100
      ! Number of Bins in Out of the Range
        integer(4):: nBOR
      ! Number of All of the Bins
        integer(4):: nAB
      ! Gas constant
        real(8),parameter:: R = 0.00198721558
      ! reweight flag
        logical(4):: reweight_flg = .false.

!***********************************************

      end module COMVAL
