
      subroutine outcrb(iwrite,iprint,filena,ixnatm,ier)
 
!*******************************************************************
!
!     THIS SUBROUTINE IS FOR WRITE BINARY COORDINATE DATA
!
!*******************************************************************
 
      use COMCMMC

      implicit none

      integer(4),intent(in):: iwrite,iprint
      character(80),intent(in):: filena
      integer(4),intent(in):: ixnatm
      integer(4),intent(out):: ier
 
      character(23):: time
      character(40):: userna
      integer(4):: len

!*****************************************

      ier = 0
      write(iprint,*)' '
      write(iprint,*)'INFORMATION> OUTCRB '
      write(iprint,*)'   WRITE BINARY COORDIANTE FILE '
      write(iprint,*)' '
      call flopen(iwrite,filena,22,'NULL',0,ier)
      if ( ier .ne. 0 ) then
        write(iprint,*)'ERROR> OUTCRB '
        write(iprint,*)'       OPEN ERROR OR CLOSE ERROR '
        write(iprint,*)' '
        ier = -1 ; return
      endif

      call infjob(time,userna,len)
 
      write(iwrite)time,userna
      write(iwrite)ixnatm
      write(iwrite)cord(1:3,1:ixnatm)

      call flclos(iwrite,10,ier)
      if ( ier .ne. 0 ) then
        write(iprint,*)'ERROR> OUTCRB '
        write(iprint,*)'       OPEN ERROR OR CLOSE ERROR '
        write(iprint,*)' '
        ier = -1 ; return
      endif

!*********************************

      return
      end subroutine outcrb
