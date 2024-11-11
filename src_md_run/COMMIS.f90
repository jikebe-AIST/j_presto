
      module COMMIS

!************************************************
!
!     COMMON AREA FOR MISCELLANEOUS VARIABLES
!
!************************************************

      use COMPAR

      ! Temperature for constraints
      real(8):: fuctmp

      ! Distance constraints
      integer(4):: iutdsc
      integer(4):: iuipar(maxdst),iujpar(maxdst)
      integer(4):: iuipnt(maxdst,6),iujpnt(maxdst,6)
      real(8):: fuwdsc
      real(8):: fudlow(maxdst),fudupr(maxdst)
      real(8):: furlow(maxdst),furupr(maxdst)

      ! Position constraints
      integer(4):: iutatm
      integer(4),allocatable:: Tiugpnt(:),iugpnt(:) ! iutatm
      real(8):: fuwpsc
      real(8),allocatable:: Tfugcns(:),fugcns(:)    ! iutatm
      real(8),allocatable:: fucord(:,:)

      ! Torsion constraints
      integer(4):: iutdhc
      integer(4):: iudhcp(maxtrc,4),iustor(maxtrc)
      real(8):: fuwdhc
      real(8):: fucdlw(maxtrc),fucdup(maxtrc)
      real(8):: fudell(maxtrc),fudelu(maxtrc)

      ! CAP-constraints
      integer(4):: iufcap,natcap,capshp
      integer(4),allocatable:: iamcap(:)
      real(8):: furcap,furcap2,furcap_pro,furcap_pro2,fukcap,          &
                fukcap_pro,rcir,rcir2,rcrate,plane_pro,rcir_pro2,      &
                rcir_pro,rcrate_pro,CAPbuff

      ! SHAKE
      integer(4):: iutshk,iugshk,iugsk2,iugsk3,iuslop,icslop
      integer(4),allocatable:: iuhshk(:),iuashk(:,:)
      real(8):: fustol,fcstol
      real(8),allocatable:: fudshk(:,:)

      ! Umbrella restraints
      integer(4):: iumbds,iumbdh
      integer(4):: iuipau(maxdst),iujpau(maxdst),iuipnu(maxdst,6)
      integer(4):: iujpnu(maxdst,6),iumdhl(maxtrc,4)
      real(8):: fuwumb
      real(8):: fumbds(maxdst,4),fumbdh(maxtrc,4)

      end module COMMIS
