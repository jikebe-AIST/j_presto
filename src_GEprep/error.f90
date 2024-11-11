
      subroutine error(iERR)

!******************************************
!
!     Output error condition
!
!******************************************

      use COMIFN,only: Nwindow ; use COMVAL,only: nBIN

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
          endif

        ! *** store_element
        elseif ( iERR .lt. 10300 ) then
          write(6,*)"place : store_element.f90"

          if ( iERR .eq. 10201 ) then
            write(6,*)"reason : Please set correct METHOD"
            write(6,*)"         You can input only ..."
            write(6,*)'           "CANO"  : Canonical MD'
            write(6,*)'           "MULT"  : Multicanonical MD'
            write(6,*)                                                 &
              '           "CLMD"  : Constant lambda molecular dynamics'
            write(6,*)                                                 &
              '           "CLMD2" : Constant lambda molecular dynamics'
            write(6,*)                                                 &
              '           "         (lambda square)                   '
            write(6,*)                                                 &
              '           "ALSD"  : Adaptive lambda square dynamics'
            write(6,*)                                                 &
              '           "ALSD2" : Adaptive lambda square dynamics'
            write(6,*)                                                 &
              '           "         (lambda square)                   '
          elseif ( iERR .eq. 10202 ) then
            write(6,*)"reason : Please input 'DATLST'"
          elseif ( iERR .eq. 10203 ) then
            write(6,*)"reason : 'DATLST' file don't exist"
          elseif ( iERR .eq. 10204 ) then
            write(6,*)"reason : Please input 'MINVAL'"
          elseif ( iERR .eq. 10205 ) then
            write(6,*)"reason : Please input 'MAXVAL'"
          elseif ( iERR .eq. 10206 ) then
            write(6,*)"reason : MINVAL > MAXVAL"
            write(6,*)"         Please set MINVAL smaller than MAXVAL"
          elseif ( iERR .eq. 10207 ) then
            write(6,*)"reason : Please set CRTEMP"
          elseif ( iERR .eq. 10208 ) then
            write(6,*)"reason : Please set correct CRTEMP"
          elseif ( iERR .eq. 10209 ) then
            write(6,*)"reason : NWINDO you input is too small"
            write(6,*)"         Number of bins = ",NBIN
            write(6,*)"                 NWINDO = ",Nwindow
          elseif ( iERR .eq. 10210 ) then
            write(6,*)"reason : NXMINV > NXMAXV"
          elseif ( iERR .eq. 10211 ) then
            write(6,*)"reason : Please input 'PREFIT'"
          elseif ( iERR .eq. 10212 ) then
            write(6,*)"reason : 'PREFIT' file don't exist"
          elseif ( iERR .eq. 10213 ) then
            write(6,*)"reason : Please set correct TMPSMP"
          elseif ( iERR .eq. 10214 ) then
            write(6,*)"reason : Please set correct NXTEMP"
          elseif ( iERR .eq. 10215 ) then
            write(6,*)"reason : Please set MIN* < MAX*"
          endif

        ! *** distrib
        elseif ( iERR .lt. 10400 ) then
          write(6,*)"place : distrib.f90"

          if ( iERR .eq. 10301 ) then
            write(6,*)"reason : The above file don't exist"
          endif
        endif

      endif

!*********************

      write(6,*)""

      stop
      end subroutine error
