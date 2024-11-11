
      module COMCMM

!*********************************************************
!
!     In this file, the quantities for CMM are defined.
!
!*********************************************************

      use COMPAR
!      use COMPAR, only: maxatm, max14n
      ! General parameters.
      ! Used generally, so do not touch.

      ! parameter for SIMD
      integer(4):: i_vec
      logical(4):: SIMD_chk

      ! General quantities for cells. (INPUT)
      integer(4):: nlev = 0             ! The deepest level to be sellected.
      real(8):: cminsiz                 ! minimum cell size
      real(8),allocatable:: siz(:)      ! The cell-size (1D) of each level
      integer(4),allocatable:: nv(:)    ! N of cells (1D) of each level in the full field

      ! Address-Address relationship. (INPUT)
      integer(4),allocatable:: lvdn(:,:) ! (2,nv(1))
      integer(4),allocatable:: lvup(:) ! (nvsecondmax)

      ! MaXimum number of Contact cells (nearest cells), don't touch
      integer(4),parameter:: icmx = 13

      ! Near cells for each level. (INPUT)
      integer(4),allocatable:: nnr(:,:,:)     ! (nv(1),nv(1),nv(1))
      integer(4),allocatable:: inr(:,:,:,:,:) !  (3(x-z),icmx,nv(1),nv(1),nv(1))
      integer(4),allocatable:: nnrv(:,:,:)    ! (nv(1),nv(1),nv(1))
      integer(4),allocatable:: inrv(:,:,:,:,:)!  (3(x-z),mx,nv(1),nv(1),nv(1))

      ! MaXimum number of second nearest cells, don't touch
      integer(4),parameter:: jcmx = 189

      ! Second nearest cells for each level. (INPUT)
      integer(4),allocatable:: ncm1(:,:,:),ncm2(:,:,:),ncm3(:,:,:),    &
                               ncm4(:,:,:),ncm5(:,:,:) ! (nv(?),nv(?),nv(?))
      integer(4),allocatable:: icm1(:,:,:,:,:),icm2(:,:,:,:,:),        &
        icm3(:,:,:,:,:),icm4(:,:,:,:,:),icm5(:,:,:,:,:)
                                     ! (3(x-z),jcmx,nv(?),nv(?),nv(?))

      ! Coordinates of cell-centers at each level. (INPUT)
      real(8),allocatable:: c1(:,:),c2(:,:),c3(:,:),c4(:,:),c5(:,:), & !(3(x-z),nv?)
        cg1(:,:),cg2(:,:),cg3(:,:),cg4(:,:),cg5(:,:) !(3(x-z),nv?)

      ! N of particles in cells in each level (npcl*)
      ! The atom number included in the cells (ipcl*)
      ! This quantity is made only for the deepest level
      integer(4),allocatable:: npcl5(:,:,:,:) ! (nfrag,nv(5),nv(5),nv(5))
      integer(4),allocatable:: npcl4(:,:,:,:) ! (nfrag,nv(4),nv(4),nv(4))
      integer(4),allocatable:: npcl3(:,:,:,:) ! (nfrag,nv(3),nv(3),nv(3))
      integer(4),allocatable:: npcl2(:,:,:,:) ! (nfrag,nv(2),nv(2),nv(2))
      integer(4),allocatable:: npcl1(:,:,:,:) ! (nfrag,nv(1),nv(1),nv(1))
      integer(4),allocatable:: ipcl(:,:,:,:,:) ! (ipmax,nfrag,nv(1),nv(1),nv(1))

      ! CMM table
      integer(4),allocatable:: Nmycel(:)        ! (nfrag)
      integer(1),allocatable:: mycel(:,:,:)     ! (3(xyz),(nv(1))**3,nfrag)
      integer(4),allocatable:: Nyrcel(:,:,:,:)    ! ((nv(1))**3,nlev,nfrag,nfrag)
      integer(1),allocatable:: yrcel(:,:,:,:,:,:) ! (3(xyz),jcmx,(nv(1))**3,nlev,nfrag,nfrag)

      ! Multipoles.
      !! (first dimension, 13 means ...
      !!   1    : mono-pole 
      !!   2-4  : dipole (x,y,z)
      !!   5-7  : quadrapole (xx,yy,zz) )
      !!   8-13 : quadrapole * 2 (xx*2,yy*2,zz*2,xy*2,yz*2,zx*2)
      real(8),allocatable:: CMMpot1(:,:,:,:,:),CMMpot2(:,:,:,:,:),     &
         CMMpot3(:,:,:,:,:),CMMpot4(:,:,:,:,:),CMMpot5(:,:,:,:,:)
                                    !(13,nfrag,nv(?),nv(?),nv(?))

      ! Input parameter from inpara
      integer(4):: nfrg,nfrag
      integer(4),allocatable:: bndcls(:)          !(iynbnd-iyndbn)
      integer(4),allocatable:: angcls(:)          !(iynang-iyndag)
      integer(4),allocatable:: torcls(:)          !(iyntor-iyndtr)
      integer(4),allocatable:: impcls(:)          !(iynimp-iyndip)
      integer(4),allocatable:: excls(:)           !(nexcess)
      integer(4),allocatable:: icls(:)            !(ixnatm)
      integer(4),allocatable:: natfr(:),iatfr(:,:)!(nfrmx),(Natm,nfrmx)
      integer(4):: ncmmatm,ncmmres
      integer(4),allocatable:: icmmatm(:,:),icmmres(:,:) !(2,ncmmatm(or res))
      integer(4),allocatable:: n_matrix(:,:)   !(nfrag,nfrag)

      ! For 1-2, 1-3, 1-4 interaction table
      integer(4):: nntab14     ! Number of 1-4 interactions
      integer(4),allocatable:: ntab14(:) ! (ixnatm)
      integer(4),allocatable:: itab14(:,:) ! (.le.max14n*2,ixnatm)

      ! Modified charge.
      real(8),allocatable:: chgmod(:),chgmod2(:) ! (ixnatm)
      ! Flag for CMM part for twin-range CMM
      integer(4):: itwin
      ! Forces for CMM part used for twin-range CMM.
      real(8),allocatable:: gwk(:,:,:) ! (3,ixnatm,nfrg)
      ! Energy for CMM part used for twin-range CMM.
      real(8),allocatable:: Ecmm(:,:)    ! (nlev,nfrg)

      ! Tables (first tables) for near-field exact calc.
      ! The atom pairs are taken by the procedue of CMM.
      ! ( vdW & electrostatic )
      integer(4),allocatable:: ntb(:,:) ! (ixnatm,nfrg)
      integer(4),allocatable:: itb(:,:,:) ! (n15mx,ixnatm,nfrg)
      ! ( only electrostatic )
      integer(4),allocatable:: ntbEL(:,:) ! (ixnatm,nfrg)
      integer(4),allocatable:: itbEL(:,:,:) ! (n15mxEL,ixnatm,nfrg)

      ! Tables (second tables) for near-field exact calc.
      ! The atom pairs are taken by the pair distances.
      integer(4),allocatable:: ntbvdW(:,:) ! (ixnatm,nfrg)
      integer(4),allocatable:: itbvdW(:,:,:) !(nvdw,ixnatm,nfrg)

      ! Tables (first tables) for near-field exact calc.
      ! The atom pairs are taken by the pair distances.
      integer(4),allocatable:: ntbhyd(:,:) ! (ixnatm,nfrg)
      integer(4),allocatable:: itbhyd(:,:,:) !(n15mx,ixnatm,nfrg)
      integer(4),allocatable:: ntbhyd2(:,:) ! (ixnatm,nfrg)
      integer(4),allocatable:: itbhyd2(:,:,:) !(nvdW,ixnatm,nfrg)

      integer(4),allocatable:: cell_nEL(:,:,:,:) !(nfrag,nv(1),nv(1),nv(1))
      real(8):: rlim,rlim2,irlim2
      integer(4):: idelay,id2p1,id2p1_2
      integer(4),allocatable:: absres(:) !(ixnatm)

      integer(4),allocatable:: inum(:)  ! ixnatm
      integer(4),allocatable:: Tixatyp(:)
      real(8),allocatable:: Tchgmod(:),Tcord(:,:)

!************************************
      ! For external CMM

      ! external CMM flag
      logical(1):: extCMM_flag = .false.

      ! CMM cell center coordinates
      real(4):: celx,cely,celz

      ! external cell calculation flag
      logical(1),allocatable:: exflag1(:,:,:),exflag2(:,:,:),          &
                exflag3(:,:,:),exflag4(:,:,:) ! (nv(?),nv(?),nv(?))

      ! Multipoles.
      !! (first dimension, 16 means ...
      !!   1    : mono-pole 
      !!   2-4  : dipole (x,y,z)
      !!   5-7  : quadrapole (xx,yy,zz) )
      !!   8-13 : quadrapole * 2 (xx*2,yy*2,zz*2,xy*2,yz*2,zx*2)
      real(8),allocatable:: eCMMpot1(:,:,:,:),eCMMpot2(:,:,:,:),       &
          eCMMpot3(:,:,:,:),eCMMpot4(:,:,:,:) !(13,nv(?),nv(?),nv(?))

      ! Some values for nearest cell calc. of external CMM
      real(8):: RR_ext__2,iRR_ext,iRR_ext__3,iRR_ext__5,i3RR_ext__5,   &
                i5RR_ext__7,RR_ext,RR_ext_2
      ! Table for nearest cell calc. of external CMM
      integer(2):: Nmycel_nEC
      integer(1),allocatable:: mycel_nEC(:,:)     ! (3(xyz),(nv(1))**3)
      integer(2),allocatable:: Nyrcel_nEC(:)    ! ((nv(1))**3)
      integer(1),allocatable:: yrcel_nEC(:,:,:) ! (3(xyz),2*icmx,(nv(1))**3)

!************************************

      ! Pseudo temperature (& start pseudo temp.)
      real(8):: psetmp,stapst
      ! Heat loop number for pseudo temp.
      integer(4):: htpslo

      ! For lambda dynamics
      real(8):: lambda,lambda_v,i_lambda_m
      real(8):: lambda_m = 0.d0
      character(4):: cluster_method

      ! No scale energy term flags
      integer(4):: scale_bond = 0
      integer(4):: scale_angle = 0
      integer(4):: scale_tor = 3
      integer(4):: scale_imp = 0

      ! lambda temperature control
      logical(4):: l_temp_control = .false.

!************************************

      ! For Zero Dipole
      real(8):: zcore,fcoeff,bcoeff
      real(8),allocatable:: dself(:) ! nfrg

      ! For repulsion
      integer(4),allocatable:: nrepcl(:,:,:,:),irepcl(:,:,:,:,:)
      integer(4),allocatable:: ntbrep(:),itbrep(:,:)
      integer(4),allocatable:: nrepnr(:,:,:),irepnr(:,:,:,:,:)

!************************************
      end module COMCMM
