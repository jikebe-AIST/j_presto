
      program ensemble_analysis

!*************************************
!
!     MAIN program for ensemble_analysis
!
!*************************************

      use COMIFN ; use COMVAL
      !$ use omp_lib

      implicit none 

      integer(4)::  time0,time1,time2
      integer(4)::  count_rate     ! count number per one second

      integer(4):: i,j,ii,jj
      logical(4):: ex
      integer(4):: nomp = 1
      character(80):: title
      character(23):: time
      character(40):: user
      integer(4):: lenusr

!**************************

      call system_clock(time0,count_rate) 
      write(6,*)'**************************************************'
      write(6,*)'*'
      write(6,*)'*    Ens_Ana (Version 1.0.0)'
      write(6,*)'*'
      write(6,*)'*                         Author : Jinzen Ikebe'
      write(6,*)'*             First Release Date : 2025-01-14'
      write(6,*)'*  Release Date for Current ver. : 2025-01-14'
      write(6,*)'*'
      write(6,*)'**************************************************'
      write(6,*)''
      !$ nomp = omp_get_max_threads()
      !$ write(6,*)
      !$ write(6,'(a,i3,a)')"     OpenMP thread number = ",nomp
      write(6,'(a)')""
      call infjob(time,user,lenusr)
      title = time//user(1:20)
      write(6,*)title
      write(6,'(a)')""

!**************************
!     Input phase

      call system_clock(time1)
      write(6,'(a)')""
      write(6,'(a)')""
      write(6,'(a)')"  INPUT phase ---------------------"
      write(6,'(a)')""
      call inp  ! ( inpinp.f90 )
      call store_element

      call system_clock(time2)
      write(6,*)""
      call output_time(6,time1,time2,count_rate,"Input")

!***************************
!     Analysis phase
      
      call system_clock(time1)
      write(6,'(a)')""
      write(6,'(a)')""
      write(6,'(a)')"  Analysis phase ------------"
      write(6,'(a)')""

      ! SOAP calc.
      if ( nSOAP .ne. 0 ) call SOAP

      ! Analysis
      if ( len(trim(input_PDB_list)) .ne. 0 ) call inpdb
      if ( len(trim(input_binary_list)) .ne. 0 ) call inbin
!      if ( len(trim(input_otherfmt_list)) .ne. 0 ) call inotherfmt
      if ( len(trim(input_restart_list)) .ne. 0 ) call inrest

      call system_clock(time2)
      write(6,*)""
      call output_time(6,time1,time2,count_rate,"ANALYSIS")
      call system("rm -rf "//trim(d_proj2)//".PDB > /dev/null")
      call system("rm -f "//trim(d_proj2)//".input* > /dev/null")

!***************************

1111  continue
      call system_clock(time2)
      write(6,*)""
      call output_time(6,time0,time2,count_rate,"TOTAL")
      write(6,*)""
      call infjob(time,user,lenusr)
      title = time//user(1:20)
      write(6,*)title

      write(6,'(a)')""
      write(6,'(a)')" +++  Program Ens_Ana Normally Ended  +++"
      write(6,'(a)')""

!*****************

      stop
      end program ensemble_analysis
