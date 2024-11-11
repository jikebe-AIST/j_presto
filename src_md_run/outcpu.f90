
      subroutine outcpu(iprint,subna,sec1,sec2,cr,cm)

      implicit none

      integer(4),intent(in):: iprint
      integer(8),intent(in):: sec1,sec2,cr,cm
      character(6):: subna
      real(8):: time

!*****************************************

      if ( sec2 .gt. sec1 ) then
        time = dble(sec2-sec1)/dble(cr)
      else
        time = dble(cm-sec1+sec2)/dble(cr)
      endif
        
      write(iprint,'(a18,a6,a5,f12.4,a5)')                             &
        "     CPU TIME FOR ",subna,"   : ",time," (S) "

      return
      end subroutine outcpu 
