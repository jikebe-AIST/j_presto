
      program main

!**********************************************
!
!     MAIN program for making PCA axis
!
!**********************************************

      implicit none 

      integer(4):: time0,time1,time2,system,iERR,count_rate
 
!**************************

     write(6,*)'**************************************************'
     write(6,*)'*'
     write(6,*)'*    PCAaxis (Version 1.0.1)'
     write(6,*)'*'
     write(6,*)'*                         Author : Jinzen Ikebe'
     write(6,*)'*             First Release Date : 2025-01-28'
     write(6,*)'*  Release Date for Current ver. : 2025-01-28'
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

      call store_element ! ( store_element.f90 )

      call system_clock(time2)
      write(6,*)""
      call output_time(6,time1,time2,count_rate,"INPUT")

!***************************
!     Make PCA axis phase

      call system_clock(time1)
      write(6,'(a)')""
      write(6,'(a)')""
      write(6,'(a)')"  Make PCA axis phase ----------------"
      write(6,'(a)')""

      call mkaxis

      call system_clock(time2)
      write(6,*)""
      call output_time(6,time1,time2,count_rate,"MAKE PCA AXIS")

!***************************

      call system_clock(time2)
      write(6,*)""
      call output_time(6,time0,time2,count_rate,"TOTAL")

      write(6,'(a)')""
      write(6,'(a)')" +++  Program PCAaxis Normally Ended  +++"
      write(6,'(a)')""

!*****************

      stop
      end program main
