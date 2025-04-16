
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

      ! project & input file list name
        character(130):: PROJNM,AXISFL
        integer(4):: naxis,nfile
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
