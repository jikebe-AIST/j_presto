
      module COMERG

!***********************************************
!
!     COMMON AREA FOR ENERGY VARIABLES
!
!***********************************************

      use COMPAR

      ! BOND PARAMETER
      integer(4):: iynbnd,iyndbn
      integer(4),allocatable:: iypbnd(:,:)
      real(8),allocatable:: fyfbnd(:),fyqbnd(:)

      ! ANGLE PARAMETER
      integer(4):: iynang,iyndag
      integer(4),allocatable:: iypang(:,:)
      real(8),allocatable:: fyfang(:),fyqang(:)

      ! TORSION PARAMETER
      integer(4):: iyntor,iyndtr
      integer(4),allocatable:: iyptor(:,:),iytdiv(:),iytnbf(:)
      real(8),allocatable:: fyftor(:),fytrot(:),fytphs(:)
      real(8),allocatable:: fytvws(:),fytess(:)

      ! IMPROPER PARAMETER
      integer(4):: iynimp,iyndip
      integer(4),allocatable:: iypimp(:,:),iyidiv(:),iyinbf(:)
      real(8),allocatable:: fyfimp(:),fyirot(:),fyiphs(:),fyivws(:),   &
                            fyiess(:)

      ! FUNCTION PARAMETER
      integer(4):: iynbpf
      integer(4):: iynpar(maxfnc)
      character(40):: cynbpf(maxfnc)

      ! NONBOND PARAMETER
      integer(4):: iyntyp
      integer(4),allocatable:: iyppid(:,:)
      real(8),allocatable:: fynbpp(:,:,:),fyvwrd(:),fyvwme(:),         &
                            fy14sv(:),fy14se(:)
      real(8):: fydiel
      integer(4),allocatable:: iyp14(:,:) ! (2,n14int)
      real(8),allocatable:: fytvws14(:),fytess14(:) ! (n14int)
      integer(4),allocatable:: i14cls(:) ! (n14int)

      ! INTERACTION TABLE
      integer(4):: iyni12,iyni13,iyni14,iyn15v,iyn15h
      real(8):: fycutl
      integer(4),allocatable:: iyninv(:),iyninh(:),iynine(:)
      integer(4),allocatable:: iyfini(:)   ! ixnres

      ! FREE ATOM INFORMATION
      integer(4):: iynvar,iynfpt,iynfp2,iynfp3,iynfp4,iynfp5,outatm
      integer(4),allocatable:: iytvar(:)
      integer(4),allocatable:: iylvar(:)  ! iynvar

      ! SHAKE INFORMATION
      integer(4):: iyfshk

      ! ENERGY
      integer(4):: iy15m ! 1-5 interaction calc. method
                         ! (1:CMM, 2:ATOM cut-off,3:RESC,4: RESA)
      integer(4):: iyeflg(maxene)
      character(6):: cyenam(maxene)

      end module COMERG
