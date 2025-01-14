
      subroutine mkNOElst

!****************************************************
!
!     MaKe NOE LiST of round robin hydrogen atoms
!          & Input experimental NOE data
!
!****************************************************

      use COMIFN ; use COMVAL

      implicit none

      ! Temporary atom number ID of the i-th HYDrogen atom
        integer(4):: ThydID(nATM)

      integer(4):: i

!****************************
!     Search hydrogen atoms in the objective residue range

      NhATM = 0
      do i = 1,nATM
        ! Remove out of residue range for the analysis
        if ( RESnum(i) .lt. fstRES ) cycle
        if ( RESnum(i) .gt. fnlRES ) exit

        if ( ATMnm(i)(1:1).eq."H" .or. ATMnm(i)(2:2).eq."H" ) then
          NhATM = NhATM + 1
          ThydID(NhATM) = i
!          write(6,'(4x,i7,a3,i0,a9,a4)')                               &
!     &          NhATM," : ",RESnum(i),"-th res. ",ATMnm(i)
        endif
      enddo

      allocate(hydID(NhATM),NOEint(NhATM,NhATM))
      hydID(1:NhATM) = ThydID(1:NhATM)
      NOEint(:,:) = 0.d0
      write(6,*)
      write(6,'(2x,a,i0)')"+ N of hydrogen atoms in the range = ",NhATM

!****************************
!     Input experimental NOE data, if you need

      if ( input_NOE_file_name .ne. " " ) then

        if ( input_NOE_file_type .eq. "XPLOR" ) then
          call Xplor_noe ! (Xplor.f90)
        elseif ( input_NOE_file_type .eq. "DISCOVER" ) then
          call DISCOVER_noe ! (DISCOVER.f90)
        endif

!        call assign_noe
      endif

!****************************

      return
      end subroutine mkNOElst
