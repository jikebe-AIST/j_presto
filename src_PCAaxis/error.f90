
      subroutine error(iERR)

!******************************************
!
!     Output error condition
!
!******************************************

      implicit none

      integer(4),intent(in):: iERR

!*********************************

      write(6,*)" "
      write(6,*)" !! ERROR !!"
      write(6,*)" ERROR number = ",iERR
      write(6,*)" "

!*****************************

      ! MAIN Phase 
      if ( iERR .lt. 10000 ) then
        write(6,*)"phase : main"

      ! INPUT Phase
      elseif ( iERR .lt. 20000 ) then
        write(6,*)"phase : input"

        ! *** asnelm
        if ( iERR .lt. 10200 ) then
          write(6,*)"place : asnelm.f90"

          if ( iERR .eq. 10101 ) then
            write(6,*)"reason : Input file style is strange"
            call system("cat help.txt")
          endif

        ! *** store_element
        elseif ( iERR .lt. 10300 ) then
          write(6,*)"place : store_element.f90"

          if ( iERR .eq. 10201 ) then
            write(6,*)"reason : Not exist the above file"
          endif

        endif

      endif

!*********************

      write(6,*)""

      stop
      end subroutine error
