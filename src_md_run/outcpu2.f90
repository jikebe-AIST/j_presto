
      subroutine outcpu2(iprint,subna,sec,cr,ftime)

      implicit none

      integer(4),intent(in):: iprint
      integer(8),intent(in):: sec,cr
      character(6):: subna
      real(8):: time,ftime

!*****************************************

      time = dble(sec)/dble(cr)
        
      write(iprint,'(a25,a6,a6,f12.4,a2,f7.4,a5)')                     &
        "                time for ",subna," :    ",time," (",          &
        time/ftime*100.d0," (%))"

      return
      end subroutine outcpu2
