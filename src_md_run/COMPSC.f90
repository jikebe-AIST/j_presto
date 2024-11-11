
      module COMPSC

!***********************************************************
!
!     COMMON AREA FOR POTENTIAL SCALING
!
!***********************************************************

      use COMPAR

      ! Multicanonical parameters
      integer(4):: nwindow,ndeg,nacntr,ndcntr,ngcntr,ntcntr
      integer(4),allocatable:: ianatm(:)  ! ixnatm
      integer(4),allocatable:: idnatm(:,:)  ! 2,iynbnd
      integer(4),allocatable:: ignatm(:,:)  ! 3,iynang
      integer(4),allocatable:: itnatm(:,:)  ! 4,iyntor

      real(8):: alphalw,alphaup,celw,ceup,tempce,lower,upper
      real(8),allocatable:: c(:,:)             ! (0:ndeg,1:nwindow)
      real(8),allocatable:: low_v(:),high_v(:) ! (nwindow)

      end module COMPSC
