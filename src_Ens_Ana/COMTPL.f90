
      module COMTPL

!**********************************************************************
!
!     For inptpl
!
!**********************************************************************

      implicit none

      ! Circular constant
      real(8),parameter:: pi = 3.141592653589793d0
      ! Radian-degree conversion unit ( pi / 180.0 )
      real(8),parameter:: rad = pi / 180.d0
      ! Avogadro number (/MOL)
      real(8),parameter:: avogad = 6.0221367d+23
      ! Boltzman constant (J/K)
      real(8),parameter:: boltz = 1.380658d-23
      ! Joule - calorie conversion unit (J/CAL)
      real(8),parameter:: joucal = 4.184d0
      real(8),parameter:: boltzk = boltz / joucal * avogad * 0.001d0

      integer(4),parameter:: maxfnc = 5
      integer(4),parameter:: maxnbp = 4
      integer(4),parameter:: max14n = 30

      integer(4):: ncTPL  ! Number of columns in the TOPOLOGY file
      integer(4):: maxatm,maxbnd,maxang,maxtor,maximp,maxtyp
      integer(4):: ixmolc,ixtitl,ixnatm,iynbnd,iynang,iyntor,iynimp,   &
        iynbpf,ixnchn,iyntyp,iynpar(maxfnc),ixnres
      integer(4),allocatable:: Tixsqml(:),ixsqml(:),ixatmm(:),         &
        ixbndm(:),ixangm(:),ixtorm(:),iximpm(:),ixtpcn(:),ixcend(:),   &
        iyppid(:,:),ixachn(:),ixamol(:),ixincd(:,:),ixatyp(:),         &
        iyptor(:,:),iytdiv(:),iytnbf(:),iypimp(:,:),iyidiv(:),         &
        iyinbf(:),ixrstr(:),ixrend(:),ixares(:),ix14if(:,:),           &
        ix14lt(:,:),iypbnd(:,:),iypang(:,:)
      real(8),allocatable:: fyvwme(:),fynbpp(:,:,:),fyvwrd(:),         &
        fy14sv(:),fy14se(:),fxincd(:,:),fyftor(:),fytrot(:),fytphs(:), &
        fytvws(:),fytess(:),fyfimp(:),fyirot(:),fyiphs(:),fyivws(:),   &
        fyiess(:),fxchrg(:),fxmass(:),fxvdwr(:),fyfbnd(:),fyqbnd(:),   &
        fyfang(:),fyqang(:)
      character(40):: cynbpf(maxfnc)
      character(80):: cxtitl
      character(4),allocatable:: cxatmt(:)
      character(8),allocatable:: cxatmn(:),cxresn(:)
      character(40),allocatable:: Tcxmolc(:),cxmolc(:)

!**********************************************************************

      end module COMTPL
