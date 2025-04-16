
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

      end module COMINP


!======================================================================


      module COMIFN

!**********************************************
!
!     COMmon data space for Input File Names
!
!**********************************************

      implicit none

        integer(4):: nfile
        character(130):: PROJNM
        character(999),allocatable:: filelist(:)

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

!***********************************************

      end module COMVAL
