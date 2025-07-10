
      program main

!**********************************************
!
!     MAIN program for PCA projection
!
!**********************************************

      implicit none 

      integer(4):: time0,time1,time2,system,iERR,count_rate
 
!**************************

     write(6,*)'**************************************************'
     write(6,*)'*'
     write(6,*)'*    PCAproj (Version 1.0.1)'
     write(6,*)'*'
     write(6,*)'*        Author : Jinzen Ikebe'
     write(6,*)'* First Release : 2025-01-29'
     write(6,*)'*  Current ver. : 2025-01-29'
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

      ! Assign parameters
      call store_element ! ( store_element.f90 )

      call system_clock(time2)
      write(6,*)""
      call output_time(6,time1,time2,count_rate,"INPUT")

!***************************
!     PCA projection phase

      call system_clock(time1)
      write(6,'(a)')""
      write(6,'(a)')""
      write(6,'(a)')"  PCA projection phase ----------------"
      write(6,'(a)')""

      call project

      call system_clock(time2)
      write(6,*)""
      call output_time(6,time1,time2,count_rate,"PCA PROJECTION")

!***************************

      call system_clock(time2)
      write(6,*)""
      call output_time(6,time0,time2,count_rate,"TOTAL")

      write(6,'(a)')""
      write(6,'(a)')" +++  Program PCAproj Normally Ended  +++"
      write(6,'(a)')""

!*****************

      stop
      end program main
