
      subroutine output_time(iprint,time1,time2,count_rate,memo)

!*************************************
!
!     output computational time
!
!*************************************

      implicit none

      integer(4),intent(in):: iprint     ! output unit number
      integer(4),intent(in):: time1      ! initial time
      integer(4),intent(in):: time2      ! final time
      integer(4),intent(in):: count_rate ! count number for 1 sec.
      character(*),intent(in):: memo     ! memo

!*****************************
!     output executed time
      write(iprint,'(a20,a,a3,f,a3)')"  Executed time for ",trim(memo),&
                   " : ",real(dble(time2-time1)/dble(count_rate)),"(s)"
      write(iprint,'(a40)')"  --------------------------------------"

!*****************************

      return
      end subroutine output_time
