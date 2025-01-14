
      subroutine inweight

!*******************************************
!
!     INput parameters for WEIGHTing
!
!*******************************************

      use COMIFN ; use COMVAL

      implicit none

      integer(4):: i,j

!****************************************

      ! Input input_prob_file
      open(unit=1,file=trim(input_prob_file),status="old")
      read(1,*)GE_method
      if ( GE_method .eq. "MULT" ) then
        i = 1
      elseif ( GE_method .eq. "ALSD" ) then
        i = 2
      else
        print*,"input_prob_file is strange" ; stop
      endif
      allocate(nBIN(i),EbinSZ(i),minEc(i))
      read(1,*)nBIN(1:i)
      read(1,*)EbinSZ(1:i)
      read(1,*)minEc(1:i)
      if ( GE_method .eq. "MULT" ) then
        allocate(wei(nBIN(1)))
        read(1,*)
        do i = 1,nBIN(1)
          read(1,*)wei(i)
        enddo
      elseif ( GE_method .eq. "ALSD" ) then
        allocate(wei2(nBIN(1),nBIN(2)))
        do j = 1,nBIN(2)
          read(1,*)
          do i = 1,nBIN(1)
            read(1,*)wei2(i,j)
          enddo
        enddo
      endif
      close(1)
      write(6,'(2x,a)')"* GE method : "//GE_method

!**********************************

      return
      end subroutine inweight
