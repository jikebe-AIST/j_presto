
      module COMBAS

      use COMPAR

      ! COMMON AREA FOR TOPOLOGY INFORMATION
      integer(4):: ncTPL      ! Number of columns in the TOPOLOGY file

      character(80):: cxtitl(10)               ! TPL. title
      integer(4):: ixtitl,ixfbou,ixcbou,ixtatb !
      integer(4):: nchncap
      logical(1),allocatable:: Tixbfcn(:)
      integer(4),allocatable:: ixbfcn(:)       ! Chain # of i-th for CAP
      real(8):: fxcell(3),invcel(3),fxcbou(3),celwal(6),celwal_pro(6)
      real(8):: fxellp(3),fxellp_pro(3)
      real(8):: ccoef

      integer(4):: ixmolc
      character(40),allocatable:: Tcxmolc(:),cxmolc(:)
      integer(4),allocatable:: Tixsqml(:),ixsqml(:),ixatmm(:),         &
            ixbndm(:),ixangm(:),ixtorm(:),iximpm(:),ixtpcn(:)

      integer(4):: ixnchn,ixnres
      integer(4),allocatable:: ixcend(:)           ! ixnchn
      integer(4),allocatable:: ixrstr(:),ixrend(:) ! ixnres

      integer(4):: ixnatm
      integer(4),allocatable:: ixamol(:),ixachn(:),ixares(:),ixatyp(:),&
        ix14if(:,:),ix14lt(:,:),ixincd(:,:)
      integer(4),allocatable:: ixamol2(:)
      real(8),allocatable:: fxchrg(:),fxmass(:),fxvdwr(:),fxincd(:,:)
      real(8),allocatable:: ifxmass(:)
      character(8),allocatable:: cxatmn(:),cxresn(:)
      character(4),allocatable:: cxatmt(:)

      real(8):: fxcpul
      integer(8):: fxcpus,cr,cm

      logical(4),allocatable:: zero_vdW(:)

      integer(4):: nomp = 1
      integer(4):: high_para
      integer(4),allocatable:: para_se(:,:,:) ! (2,np_cellmax,18)
      integer(4),allocatable:: np_cell(:) ! (18)

      integer(4):: nbnd,np_bnd,nang,np_ang,ntor,np_tor,nimp,np_imp,    &
                   n14int,np_14,nexcess,np_excess
      integer(4),allocatable:: iexcess(:,:) ! (2,nexcess)
      integer(4),allocatable:: para_bnd(:,:),para_ang(:,:),            &
        para_tor(:,:),para_imp(:,:),para_14(:,:),para_excess(:,:)

      integer(4):: nchain
      integer(4),allocatable:: iichain(:,:)

      real(8):: Krep
      integer(4),allocatable:: Nrep(:),replst(:,:)

      integer(4):: nstpcn
      integer(4):: fdebug

      end module COMBAS
