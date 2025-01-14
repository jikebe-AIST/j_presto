
      subroutine angle_calc

!**************************************************
!
!     angle calculation
!
!**************************************************

      use COMVAL ; use COMIFN

      implicit none

      integer(4):: i,i1,i2,i3,i4
      real(8):: vec(3,2),d1,d2,ang

!*****************************************

      write(uac,'(f,$)')wfac
      do i = 1,n_ang_c
        i1 = ang_c(1,i) ; i2 = ang_c(2,i)
        i3 = ang_c(3,i) ; i4 = ang_c(4,i)
        vec(1:3,1) = cod(1:3,i2) - cod(1:3,i1)
        vec(1:3,2) = cod(1:3,i4) - cod(1:3,i3)
        vec(1:3,1) = vec(1:3,1) - bsize(1:3) *                         &
                                  anint(vec(1:3,1)*ibsize(1:3))
        vec(1:3,2) = vec(1:3,2) - bsize(1:3) *                         &
                                  anint(vec(1:3,2)*ibsize(1:3))
        d1 = vec(1,1)*vec(1,1) + vec(2,1)*vec(2,1) + vec(3,1)*vec(3,1)
        d1 = sqrt(d1)
        d2 = vec(1,2)*vec(1,2) + vec(2,2)*vec(2,2) + vec(3,2)*vec(3,2)
        d2 = sqrt(d2)
        ang = vec(1,1)*vec(1,2) + vec(2,1)*vec(2,2) + vec(3,1)*vec(3,2)
        ang = acos(ang / (d1*d2)) * rpi
        write(uac,'(x,f8.3,$)')ang
      enddo
      write(uac,*)

!*****************************************

      return
      end subroutine angle_calc
