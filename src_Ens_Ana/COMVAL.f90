
      module COMVAL

!**************************************************
!
!     COMmon data space for some VALues
!
!**************************************************

      implicit none

      ! Directory name
        character(200):: d_proj = ""
        character(200):: d_proj2 = ""
        character(200):: d_PDB = ""
        character(200):: d_DSSP = ""
        character(200):: d_TOR = ""

! File Unit Number
      ! Unit number for Output PDB files
        integer(4):: uop = 7
      ! Unit number for PCA Coordinates
        integer(4):: upc = 12
      ! Unit number for Rg
        integer(4):: urg = 13
      ! Unit number for Qvalue
        integer(4):: uqv = 14
      ! Unit number for durability
        integer(4):: udr = 15
      ! Unit number for ourput RMsd
        integer(4):: urm = 16
      ! Unit number for output time series RMSD
        integer(4):: utr = 17
      ! Unit number for output burial info.
        integer(4):: ubu = 18
      ! Unit number for surface analysis
        integer(4):: usf = 19
      ! Unit number for sscs, SS & tor file
        integer(4):: usc = 20
        integer(4):: uss = 21
        integer(4):: utmp = 22
      ! Unit number for distance & angle output files
        integer(4):: udc = 23
        integer(4):: uac = 24
        integer(4):: utc = 25
      ! Unit number for pocket volume
        integer(4):: upv = 26
      ! Unit number for output burial info. for each atom
        integer(4):: ubu_e = 1000

      ! Number of AToM in template_PDB
        integer(4):: nATM
      ! AToM NUMber array
        integer(4), allocatable:: ATMnum(:)
      ! RESidue name of i-th RESidue
        character(4),allocatable:: aRES(:)
      ! AToM NaMe array
        character(4),allocatable:: ATMnm(:)
      ! RESidue NaMe array
        character(4),allocatable:: RESnm(:)
      ! CHaiN name array
        character(1),allocatable:: CHN(:)
      ! CHaiN Number array
        integer(4),allocatable:: nCHN(:)
      ! RESidue NUMber array
        integer(4),allocatable:: RESnum(:)
      ! FirST and FiNaL number of RESidue for analyses
        integer(4):: fstRES,fnlRES
      ! COorDinates
        real(4),allocatable:: cod(:,:)
      ! COorDinates in template PDB
        real(4),allocatable:: Rcod(:,:)
      ! Occupancy (mass)
        real(8),allocatable:: occupy(:)
      ! B-factor (charge)
        real(8),allocatable:: bfact(:)
      ! Weighting FACtor
        real(8):: wfac
      ! SUMmation of Weighting factor
        real(8):: sumW,sumW2
      ! gas constant
        real(8),parameter:: Rgas = 0.00198721558
      ! min. & max coordinates
        real(4):: maxX = -999999.d0
        real(4):: maxY = -999999.d0
        real(4):: maxZ = -999999.d0
        real(4):: minX =  999999.d0
        real(4):: minY =  999999.d0
        real(4):: minZ =  999999.d0
      ! start & end atom number for output
        integer(4):: istATM,ienATM
      ! Output atom flag
        logical(1),allocatable:: outF(:)

! For weighting
      ! weight method (MULT or ALSD)
        character(4):: GE_method
      ! weight for MULT or ALSD
        real(8),allocatable:: wei(:),wei2(:,:)
      ! number of bins
        integer(4),allocatable:: nBIN(:)
      ! energy bin size
        real(8),allocatable:: EbinSZ(:)
      ! minimum energy
        real(8),allocatable:: minEc(:)

! For PCA analysis
      ! Number of atoms for the atom list
        integer(4):: nPCA
      ! array to keep atom number for the atom list
        integer(4),allocatable:: iATMPCA(:)
      ! Number of residues
        integer(4):: NresPCA
      ! First & Final atom number of residue
        integer(4),allocatable:: iresPCA(:)
      ! Number of residue contacts
        integer(4):: Ncon_res_PCA
      ! Order residue number of contact pairs
        integer(4),allocatable:: con_res_PCA(:,:)
      ! vdW radius
        real(8),allocatable:: radPCA(:)

! For Qvalue
      ! Number of atoms for the atom list
        integer(4):: nQ
      ! array to keep atom number for the atom list
        integer(4),allocatable:: iATMQ(:)
      ! Number of residues
        integer(4):: NresQ
      ! First & Final atom number of residue for Qvalue
        integer(4),allocatable:: iresQ(:)
      ! Number of native contact pairs
        integer(4):: nNCP
      ! Native contact pairs
        integer(4),allocatable:: iNCP(:,:)
      ! Qvalue
        real(4):: Qvalue = 0.0
        real(4):: Qvalue2 = 0.0
      ! vdW radii for Qvalue contact
        real(8),allocatable:: radQ(:),radNCP(:)
        real(8),allocatable:: Qcon(:),Dmax(:),Color(:)
        logical(1),allocatable:: FtraQcut(:)
        real(8),allocatable:: Qpoint(:,:)
        real(4):: durability = 0.0
        integer(4):: traQ_zero = 0

! For distribution analysis
      ! Number of atoms for the atom list
        integer(4):: ndist
      ! array to keep atom number for the atom list
        integer(4),allocatable:: iATMdist(:)
      ! vdW radii for distribution analysis
        real(8),allocatable:: raddist(:)
      ! distribution cell array
        real(8),allocatable:: dcell(:,:,:)
        real(8),allocatable:: dcell_ele(:,:,:)
      ! min & max cell array number
        integer(4):: icminX,icmaxX,icminY,icmaxY,icminZ,icmaxZ
        integer(4):: NdcelX,NdcelY,NdcelZ

! For NOE
      ! Number of Hydrogen AToMs
        integer(4):: NhATM
      ! atom number ID of the i-th HYDrogen atom
        integer(4),allocatable:: hydID(:)
      ! NOE INTensity
        real(8),allocatable:: NOEint(:,:)
      ! Experimental NOE distances & UPper, LOWer Boundary
        real(4),allocatable:: eNOE(:), UPb(:),LOWb(:)
      ! Number of Experimental NOE pairs in the range
        integer(4):: NeNOE
      ! RESidue number of an experimental NOE pairs
        integer(4),allocatable:: iRES(:,:)
      ! AToM name of experimental NOE pairs
        character(4),allocatable:: aATM(:,:)

! For JCC
      ! RESidue number of an Jcc pair
        integer(4),allocatable:: iRESJ(:,:)
      ! AToM name of an Jcc pair
        character(4),allocatable:: aATMJ(:,:)
      ! AVErage dihedral Angle of an jcc pair
        real(4),allocatable:: AveA(:)
      ! DeLTa dihedral Angle from the average jcc one
        real(4),allocatable:: DltA(:)
      ! Number of JCC pairs
        integer(4):: nJCC

! For RMSD calc.
        ! Number of AToM in reference_PDB
        integer(4):: nATMr
      ! RESidue name of i-th RESidue
        character(4),allocatable:: aRESr(:)
      ! AToM NUMber array
        integer(4), allocatable:: ATMnumR(:)
      ! AToM NaMe array
        character(4),allocatable:: ATMnmR(:)
      ! RESidue NaMe array
        character(4),allocatable:: RESnmR(:)
      ! CHaiN name array
        character(1),allocatable:: CHNr(:)
      ! CHaiN number array
        integer(4),allocatable:: nCHNr(:)
      ! RESidue NUMber array
        integer(4),allocatable:: RESnumR(:)
      ! REFerence-COorDinates ( not Rcod )
        real(4),allocatable:: refcod(:,:)

      ! Number of atoms used in RMSD calc.
        integer(4):: nRMSD,nRMSDf
      ! array to keep atom number for RMSD calc.
        integer(4),allocatable:: iATMrmsd(:),iATMrmsdR(:)
        integer(4),allocatable:: iATMrmsdF(:),iATMrmsdRF(:)
      ! Number of ref. conformations
        integer(4):: nREF
      ! Array for the original conformation number & weight
        integer(4),allocatable:: conf_n(:)
        real(8),allocatable:: conf_w(:)
      ! COorDinates of ref. atoms for RMSD calc.
        real(4),allocatable:: RcodRMSDf(:,:,:),RcodRMSD(:,:,:)
      ! AVErage structure COorDinates of the ref. ensemble
        real(4),allocatable:: aveRcodF(:,:),aveRcod(:,:)
      ! average RMSD
        real(8),allocatable:: sumRMSD(:),sumRMSD2(:),sWrmsd(:)
      ! For time series RMSD coordinates
        real(4),allocatable:: codTSRMSDf(:,:,:),codTSRMSD(:,:,:)

      ! Flag whether # of the ref. conf. is 1 or not
        logical(4):: ref_conf_1_flag = .false.
      ! rotation flag for distribution analysis or output conf.
        logical(4):: rotation_flag = .false.
      ! interval for timeseries RMSD
        integer(4):: interval_for_timeseries_RMSD

! For Rg calc.
      ! Number of atoms used in Rg calc.
        integer(4):: nRg
      ! array to keep atom number for Rg calc.
        integer(4),allocatable:: iATMrg(:)
      ! Rg
        real(8):: Rg = 0.d0
        real(8):: Rg2 = 0.d0

! For contact analysis
      ! Number of atoms used in contact analysis
        integer(4):: nprob,ntgt
      ! array to keep atom number for contact analysis
        integer(4),allocatable:: iATMprob(:),iATMtgt(:)
      ! contact ratio
        real(8),allocatable:: rcon(:)
      ! water radius for contact
        real(8):: watrad
      ! Cell size for contact analysis
        real(8):: cellsz_con
      ! Number of particle in a cell
        integer(4):: np_cell
      ! contact correlation
        real(8),allocatable:: con_cor(:,:)
      ! average ASA
        real(8),allocatable:: ave_surface(:)
      ! vdW radii
        real(8),allocatable:: radprob(:),radtgt(:)

! For burial
      ! ave_burial
        real(4),allocatable:: ave_burial(:)

! For distance & angle calc
      ! distance calc. list
        integer(4),allocatable:: dist_c(:,:)
        integer(4):: n_dist_c = 0
      ! angle calc. list
        integer(4),allocatable:: ang_c(:,:)
        integer(4):: n_ang_c = 0
        real(8),parameter:: npi = 3.141592653589793d0
        real(8),parameter:: rpi = 180.d0 / npi
      ! dihedral angle calc. list
        integer(4),allocatable:: dih_c(:,:)
        integer(4):: n_dih_c = 0

! For RMSF calc.
     ! average coordinate vecters & the square
       real(8),allocatable:: acod(:,:)
       real(8),allocatable:: a2cod(:)

! For pocket search analysis
      ! Number of atoms for the atom list
        integer(4):: nrec
        integer(4):: nlig = 0
      ! array to keep atom number for the atom list
        integer(4),allocatable:: iATMrec(:)
        integer(4),allocatable:: iATMlig(:)
      ! vdW radii for pocket search analysis
        real(8),allocatable:: radrec(:)
        real(8),allocatable:: radlig(:)
      ! pocket search cell array
        real(8),allocatable:: pcell(:,:,:)
      ! min & max cell array number
        integer(4):: pcminX,pcmaxX,pcminY,pcmaxY,pcminZ,pcmaxZ
        integer(4):: NpcelX,NpcelY,NpcelZ
      ! pocket atom search
        real(8),allocatable:: pkt_touch(:)

!***********************************************

      end module COMVAL
