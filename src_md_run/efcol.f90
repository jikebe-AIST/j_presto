
      function efcol(line,space,comlbl)
 
!*******************************************************************
!
!     THIS FUNCTION RETURNS EFFECTIVE COLUMN NUMBER
!       AND CONVERT CHARACTER-CODE
!
!*******************************************************************
 
      implicit none

      ! Effective column number ( line(1:efcol) is meaningful )
        integer(4):: efcol
      ! Input line
        character(80),intent(in):: line
      ! All space line
        character(80),intent(in):: space
      ! Label for comment
        character(1),intent(in):: comlbl

      logical(1):: onnum = .false.
      integer(4):: i,itocol
 
!*******************************************************
 
!     * ATTENTION *
!       IN FACOM OS4/F4-MSP , LINE(73:80) IS LINE-NUMBER
!       IN NUMBERING DATA SET . (LAST 8-BYTES IN EACH LINE ARE
!       LINE-NUMBER)

      do i = 73,80
        if ( index('0123456789',line(i:i)) .eq. 0 ) then
          onnum = .false. ; exit
        endif
      enddo

      if ( onnum ) then
        itocol = 72
      else
        itocol = 80
      endif

      if ( line(1:itocol) .eq. space(1:itocol) ) then
        efcol = -1 ; return
      endif

      efcol = index(line(1:itocol),comlbl) -  1
      if ( efcol .eq. -1 ) efcol = itocol

!**********************************

      return
      end function efcol


!*******************************************************************


      function efcol2(iu,comlbl,ier)

!*******************************************************************
!
!     Return characters in the next effective column
!       except for space
!
!*******************************************************************

      implicit none

      character(80):: efcol2
      ! Unit number for read
        integer(4),intent(in):: iu
      ! Label for comment
        character(1),intent(in):: comlbl
      ! error
        integer(4),intent(out):: ier

      character(80):: space = " "
      character(80):: line
      integer(4):: i

!*******************************************************

      ier = 0
      do
        read(iu,'(a80)',end=800)line
        i = index(line,comlbl)
        if ( i .ne. 0 ) line = line(1:i-1)
        if ( line .eq. space ) cycle
        efcol2 = line
        return
800     ier = 1
        return
      enddo

!**********************************

      return
      end function efcol2
