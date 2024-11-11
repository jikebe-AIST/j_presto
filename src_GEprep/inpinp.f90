!----------------------------------
!
!     subroutine inp
!     subroutine rdline
!     subroutine asnELM
!
!----------------------------------


      subroutine inp

!*********************************
!
!     INPut parameter
!
!*********************************

      use COMINP

      implicit none

      ! Number of DATa In 1 Line
        integer(4):: NDATi1L
      ! LINE in input data
        character(130):: line
      ! INPut DATa separated by spaces
        character(130):: inpdat(100)

!*********************************

      inpdat(1:100) = " "

!     Input parameter
      do 
        ! Input data file
        read(5,'(a)',end=800)line
        if ( line .eq. " " ) cycle
        call rdline(line,";",NDATi1L,inpdat)  ! ( inpinp.f90 )

        ! Assign data to each element
        call asnELM(nELM,namELM,NDATi1L,inpdat,ELEMNT) ! ( inpinp.f90 )
      enddo

800   continue

!********************

      return
      end subroutine inp


!====================================================================


      subroutine rdline(line,CMTout,NDATi1L,inpdat)

!*******************************************************************
!
!     Separate a line to some character blocks deliminated by spaces
!       ( ReaD LINE )
!
!*******************************************************************

      implicit none

      ! input LINE 
        character(*),intent(inout):: line
      ! delimiter for CoMmenT out
        character(1),intent(in):: CMTout
      ! Number of DATa in 1 Line
        integer(4),intent(out):: NDATi1L
      ! INPut DATa separated by spaces
        character(*),intent(inout):: inpdat(100)

      integer(4):: i,ist,ien

!******************************************

      inpdat(:) = " " ; ist = 0 ; ien = 0 ; NDATi1L = 0

      ! Delete line after the delimiter
      i = index(line,CMTout)
      if ( i .ne. 0 ) line = line(1:index(line,CMTout)-1) 

      ! Find space
      do i = 1,len(line)
        if ( index(line(i:i)," ").eq.0 .and. ist.eq.0 ) ist = i
        if ( index(line(i:i)," ").ne.0 .and. ien.eq.0 .and. ist.ne.0 ) &
             ien = i - 1

        ! determination of inpdat
        if ( ist.ne.0 .and. ien.ne.0 ) then
          NDATi1L = NDATi1L + 1
          inpdat(NDATi1L) = line(ist:ien)
          ist = 0 ; ien = 0
        endif

        ! Case the final character in line is not space
        if ( ist.ne.0 .and. i.eq.len(line) .and. ien.eq.0 ) then
          NDATi1L = NDATi1L + 1
          inpdat(NDATi1L) = line(ist:i)
        endif

      enddo

!*******************************

      return
      end subroutine rdline


!=============================================================


      subroutine asnELM(n,namELM,NDATi1L,inpdat,ELEMNT)

!*********************************
!
!     ASsigN ELeMents
!
!*********************************

      implicit none

      ! Number of elements
        integer(4),intent(in):: n
      ! control data flag NAMe of ELeMents
        character(*),intent(in):: namELM(n)
      ! Number of DATa In 1 Line
        integer(4),intent(in):: NDATi1L
      ! INPut DATa separated by spaces
        character(*),intent(in):: inpdat(100)
      ! ELEMeNTs
        character(*),intent(inout):: ELEMNT(n)

      integer(4):: i,j,iERR

!****************************************

      ! Search "=" mark
      do i = 2,NDATi1L-1
        if ( inpdat(i) .eq. "=" ) then
          iERR = 1

          ! Assign elements
          do j = 1,n
            if ( inpdat(i-1) .eq. namELM(j) ) then
              ELEMNT(j) = inpdat(i+1)
              iERR = 0 ; exit
            endif
          enddo

          if ( iERR .eq. 1 ) then
            write(6,*)inpdat(i-1:i+1)
            call error(10101)
          endif

        endif
       enddo

!****************************************

       return
       end subroutine asnELM
