
      program main

!*************************************

      use COMIFN,only: METHOD

      implicit none 

      integer(4):: time0,time1,time2,system,count_rate
      character(130):: work
 
!**************************

     write(6,*)'**************************************************'
     write(6,*)'*'
     write(6,*)'*    GEprep (Version 1.0.2)'
     write(6,*)'*'
     write(6,*)'*                         Author : Jinzen Ikebe'
     write(6,*)'*             First Release Date : 2024-11-11'
     write(6,*)'*  Release Date for Current ver. : 2025-05-23'
     write(6,*)'*'
     write(6,*)'**************************************************'
     write(6,*)''

!**************************
!     INPUT phase

      call system_clock(time0,count_rate) 
      call system_clock(time1)
      write(6,'(a)')""
      write(6,'(a)')""
      write(6,'(a)')"  INPUT phase ---------------------"
      write(6,'(a)')""

      ! Input parameters
      call inp  ! ( inpinp.f90 )
      ! Assign parameters
      call store_element ! ( store_element.f90 )

      call system_clock(time2)
      write(6,*)""
      call output_time(6,time1,time2,count_rate,"INPUT")

!***************************
!     distribution phase

      call system_clock(time1)
      write(6,'(a)')""
      write(6,'(a)')""
      write(6,'(a)')"  DISTRIBUTION phase ----------------"
      write(6,'(a)')""

      if ( METHOD .eq. "CANO" ) then
        call distrib_CANO
      elseif ( METHOD .eq. "MULT" ) then
        call distrib_MULT
      elseif ( METHOD(1:4) .eq. "CLMD" ) then
        call distrib_CLMD
      elseif ( METHOD(1:4) .eq. "ALSD" ) then
        call distrib_ALSD
      endif

      call system_clock(time2)
      write(6,*)""
      call output_time(6,time1,time2,count_rate,"DISTRIBUTION")

!***************************

      call system_clock(time2)
      write(6,*)""
      call output_time(6,time0,time2,count_rate,"TOTAL")

      write(6,'(a)')""
      write(6,'(a)')" +++  Program GEprep Normally Ended  +++"
      write(6,'(a)')""

!*****************

      stop
      end program main
