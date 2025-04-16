
      program main

      use COMBAS ; use COMCMM, only: nfrg ; use COMCMMC ; use CALC_TIME
      !$ use omp_lib

      implicit none

      integer(4):: efcol,ier,iread,iprint,icol,iersub,iend,i
      integer(8):: sec
      character(8):: phase_name
      character(80):: line,space
      logical(1):: onend,phase_check
      character(80):: title
      character(23):: time
      character(40):: user
      integer(4):: lenusr

!********************************************

!     <<<  INITIALIZATION  >>>

      call system_clock(fxcpus,cr,cm)

      ier = 0 ; space = ' ' ; iread = 5  ; iprint = 6 ; onend = .false.
      iersub = 0 ; phase_check = .false.
      call iniprs
     write(iprint,*)'**************************************************'
     write(iprint,*)'*'
     write(iprint,*)'*    md_run (Version 1.0.1)'
     write(iprint,*)'*'
     write(iprint,*)'*                         Author : Jinzen Ikebe'
     write(iprint,*)'*             First Release Date : 2024-11-11'
     write(iprint,*)'*  Release Date for Current ver. : 2025-01-14'
     write(iprint,*)'*'
     write(iprint,*)'**************************************************'
     write(iprint,*)''
      !$ nomp = omp_get_max_threads()
      !$ write(iprint,'(a,i3,a)')                                    & !
      !$   " *         OpenMP thread number = ",nomp,"          *"
      call infjob(time,user,lenusr)
      title = time//user(1:20)
      write(iprint,*)title

!     <<<  MAIN LOOP  >>>
      do
        read(iread,'(a80)',iostat=iend)line

        if ( iend .eq. 0 ) then
          icol = efcol(line,space,';')
          if ( icol.le.0 .or. index(line(1:icol),'EXE>').eq.0 ) cycle

          !! Phase_name check
          if ( index(line(1:icol),"END") .ne. 0 ) then
            phase_name = "END" ; phase_check = .false.
          elseif ( index(line(1:icol),"INPUT") .ne. 0 ) then
            phase_name = "INPUT" ; phase_check = .true.
          elseif ( index(line(1:icol),"MIN") .ne. 0 ) then
            phase_name = "MINI" ; phase_check = .true.
          elseif ( index(line(1:icol),"MD") .ne. 0 ) then
            phase_name = "MD" ; phase_check = .true.
          elseif ( index(line(1:icol),"ANA") .ne. 0 ) then
            phase_name = "ANALYS" ; phase_check = .true.
          elseif ( index(line(1:icol),"OUTPUT") .ne. 0 ) then
            phase_name = "OUTPUT" ; phase_check = .true.
          else
            phase_check = .false.
          endif
        else
          phase_name = "END" ; phase_check = .false.
        endif

        !! Go to each phase
        if ( phase_check ) then
          write(iprint,*)' '
          write(iprint,*)'++++++++++++++++++++++++++++++++++++++++'
          write(iprint,*)'+                                      +'
          write(iprint,*)'+ INFORMATION> J_PRESTO                +'
          write(iprint,'(x,a12,a8,a20)')                               &
                   '+           ',phase_name,' starts.           +'
          write(iprint,*)'+                                      +'
          write(iprint,*)'++++++++++++++++++++++++++++++++++++++++'
          write(iprint,*)' '

          if ( phase_name .eq. "INPUT" ) then
            call input(iread,iprint,ier,onend)
            allocate(grad(3,ixnatm,nfrg))
          elseif ( phase_name .eq. "MINI" ) then
            call mini(iread,iprint,ier,onend)
          elseif ( phase_name .eq. "MD" ) then
            allocate(vel(3,ixnatm))
            call md(iread,iprint,ier,onend)
          elseif ( phase_name .eq. "OUTPUT" ) then
            call output(iread,iprint,iersub,onend)
          endif
        endif

        !! EXIT loop
        if ( phase_name.eq."END" .or. ier.ne.0 .or. onend .or.         &
             iersub.ne.0 ) exit
      enddo

!===================================================================

      call system_clock(sec)
      write(iprint,*)"* Memory info."
      write(iprint,*)"  n15mxEL & the Max. = ",n15mxEL,n15mxEL_max
      write(iprint,*)"  n15mx   & the Max. = ",n15mx,n15mx_max
      write(iprint,*)"  nvdw    & the Max. = ",nvdw,nvdw_max
      write(iprint,*)"  ipmax   & the Max. = ",ipmax,ipmax_max
      write(iprint,*)' '
      write(iprint,*)'++++++++++++++++++++++++++++++++++++++++'
      write(iprint,*)'+                                      +'
      write(iprint,*)'+ INFORMATION> J_PRESTO                +'
      write(iprint,*)'+       Job has finished.              +'
      write(iprint,*)'+                                      +'
      write(iprint,*)'++++++++++++++++++++++++++++++++++++++++'
      write(iprint,*)' '
      write(iprint,*)' '
      call outcpu(iprint,'TOTAL ',fxcpus,sec,cr,cm)
      write(iprint,*)' '
      call infjob(time,user,lenusr)
      title = time//user(1:20)
      write(iprint,*)title
!      write(6,*)"bakatime= ",dble(tmptime)/dble(cr),"(S)"

      stop
      end program main


!===================================================================


      subroutine iniprs

!********************************************
!
!     This subroutine is for initialization
!
!********************************************

      use COMBAS ; use COMERG ; use COMMIS ; use COMCMM

      implicit none

!*************************

      ! 1) COMBAS
      ixfbou = 0 ; ixcbou = 0 ; fxcpul = 99999999.d0
      celwal(1:6) = 0.d0
      fxcell(1:3) = 40.d0 ; fxcbou(1:3) = 0.d0 ; fxellp(1:3) = 30.d0

      ! 2) COMERG
      iyndbn = 0 ; iyndtr = 0    ; iyndip = 0
      iyfshk = 0 ; fydiel = 1.d0 ; fycutl = 12.d0
      iyeflg(1:9) = 1 ; iyeflg(10:maxene) = 0
      cyenam(1:15) = (/'TOTAL ','BOND  ','ANGLE ','TORS. ','IMPRO.',   &
                       'VDW14 ','ELE14 ','VDW15 ','ELE15 ','HYD.  ',   &
                       'PSR.  ','DSR.  ','DHR.  ','CAP   ','REPUL.'/)
      iyn15v = 0 ; iyn15h = 0

      ! 5) COMMIS
      fuctmp = 300.d0 ; fuwpsc = 1.d0 ; fuwdsc = 1.d0 ; fuwdhc = 1.d0
      furcap = 20.d0  ; CAPbuff = 4.d0 ; furcap_pro = -1.d0
      fukcap = 150.d0 ; fukcap_pro = 150.d0
      iufcap = 1      ; iuslop = 1000  ; icslop = 3
      fustol = 1.d-6  ; fcstol = 1.d-3 ; capshp = 1

      ! 6) COMCMM
      iy15m = 4  ; nfrg = 1 ; itwin = 0
      lambda = 0.d0 ; lambda_v = 0.d0

      return
      end subroutine iniprs
