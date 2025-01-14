
      subroutine mkJCClst

!****************************************************
!
!     MaKe JCC LiST & Input experimental JCC data
!
!****************************************************

      use COMIFN ; use COMVAL

      implicit none

!****************************
!     Input experimental JCC data, if you need

      if ( input_JCC_file_name .ne. " " ) then

        if ( input_JCC_file_type .eq. "XPLOR" ) then
          call Xplor_jcc ! (Xplor.f90)
        elseif ( input_JCC_file_type .eq. "DISCOVER" ) then
          call DISCOVER_jcc ! (DISCOVER.f90)
        endif
!        call mk_jcc

      endif

!****************************

      return
      end subroutine mkJCClst
