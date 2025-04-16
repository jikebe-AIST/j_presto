
      subroutine store_element

!*******************************************************
!
!     store element data to variable of input for ToPoLogy
!
!*******************************************************

      use COMINP ; use COMIFN ; use COMVAL

      implicit none

      integer(4):: iERR = 0
      integer(4):: i,nargs,system
      character(999),allocatable:: argv(:),tlist(:)
      logical(1):: is_arg,ex

!***************************

      nargs = iargc() ; allocate(argv(nargs))
      do i = 1,nargs
        call getarg(i,argv(i))
      enddo
      if ( trim(argv(1)) .ne. "j_presto" ) stop

      ! Help option handling
      do i = 1,nargs
        if ( trim(argv(i)) .eq. "-h" .or. trim(argv(i)) .eq. "--help" ) then
          write(6,'(a)')"Usage: j_presto PCAaxis [options]"
          write(6,'(a)')"Options:"
          write(6,'(a)')"  -h, --help       "//                        &
                        "Show this help message and exit"
          write(6,'(a)')"  -o, --output     "//                        &
                        "Specify output axis file name "//             &
                        '(defalt is "test")'
          write(6,'(a)')"  -a, --axis     PC axis file name "
          write(6,'(a)')"  -i, --input      "//                        &
                        "Specify input PCAcod files (required)"
          write(6,'(a)')"  -n, --num_components      "//               &
                        "N Specifies the number of PC axes for "//     &
                        "projection, ranging from the first (largest"//&
                        ") to the N-th principal component"
          stop
        endif
      enddo

      ! input PROJect NaMe
      PROJNM = "test"
      do i = 1,nargs-1
        if ( trim(argv(i)).eq."-o" .or.                                &
             trim(argv(i)).eq."--outout" ) then
          PROJNM = trim(argv(i+1)) ; exit
        endif
      enddo
      write(6,'(2x,a)')"* Input Project name : "
      write(6,'(8x,a)')trim(PROJNM)

      ! Axis file name
      AXISFL = ""
      do i = 1,nargs-1
        if ( trim(argv(i)).eq."-a" .or.                                &
             trim(argv(i)).eq."--axis" ) then
          AXISFL = trim(argv(i+1)) ; exit
        endif
      enddo
      write(6,'(2x,a)')"* Axis file : "
      write(6,'(8x,a)')trim(AXISFL)
      if ( len(trim(AXISFL)) .eq. 0 ) call error(10201)
      inquire(file=trim(AXISFL),exist=ex)
      if ( .not. ex ) call error(10202)

      ! INPut file names
      allocate(tlist(nargs)) ; tlist(:) = "" ; is_arg = .false.
      nfile = 0
      do i = 1,nargs
        if ( argv(i)(1:1) .eq. "-" ) then
          if ( trim(argv(i)).eq."-i" .or.                              &
               trim(argv(i)).eq."--input" ) then
            is_arg = .true. ; cycle
          else
            is_arg = .false.
          endif
        endif
        if ( is_arg ) then
          nfile = nfile + 1 ; tlist(nfile) = trim(argv(i))
        endif
      enddo
      allocate(filelist(nfile)) ; filelist(:) = tlist(1:nfile)
      write(6,'(2x,a)')"* Input file(s) : "
      do i = 1,nfile
        write(6,'(8x,a)')trim(filelist(i))
        inquire(file=trim(filelist(i)), exist=ex)
        if ( .not. ex ) call error(10201)
      enddo

      ! Number of axes for projection
      naxis = 1
      do i = 1,nargs-1
        if ( trim(argv(i)).eq."-n" .or.                                &
             trim(argv(i)).eq."--num_components" ) then
          read(argv(i+1),*)naxis ; exit
        endif
      enddo
      write(6,'(2x,a)')"* Number of PC axes for projection : "
      write(6,'(8x,i0)')naxis

!***************************

      return
      end subroutine store_element
