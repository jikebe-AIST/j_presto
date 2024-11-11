
      subroutine inpcrb(iread,iprint,filena,maxatm,numatm,cord,ier)
 
!*******************************************************************
!
!     THIS SUBROUTINE IS FOR READ BINARY COORDINATE DATA
!
!*******************************************************************
 
      implicit none

      integer(4),intent(in):: iread,iprint,maxatm
      character(80),intent(in):: filena
      real(8),intent(out):: cord(3,maxatm)
      integer(4),intent(out):: numatm,ier
 
      character(23):: time
      character(40):: userna
      integer(4):: iersub
 
!***********************************************************

      ier = 0
      call flopen(iread,filena,20,'NULL',0,iersub)

      if ( iersub .lt. 0 ) then
        write(iprint,*)' '
        write(iprint,*)'ERROR> INPCRB '
        write(iprint,*)'       OPEN ERROR FOR BINARY COORDINATE '
        write(iprint,*)' '
        ier = iersub ; return
      endif

      read(iread)time,userna
      read(iread)numatm
      write(iprint,*)' '
      write(iprint,*)'INFORMATION> INPCRB '
      write(iprint,*)'             INPUT BINARY COORDINATE FILE '
      write(iprint,*)' '
      write(iprint,*)'      TIME            : ',time
      write(iprint,*)'      USERNAME        : ',userna
      write(iprint,*)'      NUMBER OF ATOMS : ',numatm
      write(iprint,*)' '

      if ( numatm .gt. maxatm ) then
        write(iprint,*)'ERROR> INPCRB '
        write(iprint,*)'         NUMBER OF ATOMS IN BINARY COORDINATE'
        write(iprint,*)'         DATA FILE IS OUT OF THE LIMIT '
        write(iprint,*)' '
        write(iprint,*)'  NUMBER OF ATOMS IN BINARY COORD. : ',numatm
        write(iprint,*)'  THE LIMIT IS                     : ',maxatm
        write(iprint,*)' '
        ier = -1750 ; return
      endif
      read(iread)cord(1,1:numatm)
      read(iread)cord(2,1:numatm)
      read(iread)cord(3,1:numatm)
      call flclos(iread,10,ier)

      if ( ier .ne. 0 ) then
        write(iprint,*)'ERROR> INPCRB '
        write(iprint,*)'       FILE CLOSE ERROR '
        write(iprint,*)' '
      endif

!*****************************

      return
      end subroutine inpcrb
