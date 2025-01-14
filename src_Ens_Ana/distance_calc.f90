
      subroutine distance_calc

!**************************************************
!
!     distance pairs calculation
!
!**************************************************

      use COMVAL ; use COMIFN

      implicit none

      integer(4):: i,i1,i2
      real(8):: x,y,z,d

!*****************************************

      write(udc,'(f,$)')wfac
      do i = 1,n_dist_c
        i1 = dist_c(1,i) ; i2 = dist_c(2,i)
        x = cod(1,i1) - cod(1,i2)
        y = cod(2,i1) - cod(2,i2)
        z = cod(3,i1) - cod(3,i2)
        x = x - bsize(1) * nint(x*ibsize(1))
        y = y - bsize(2) * nint(y*ibsize(2))
        z = z - bsize(3) * nint(z*ibsize(3))
        d = sqrt(x*x + y*y + z*z)
        write(udc,'(x,f8.3,$)')d
      enddo
      write(udc,*)

!*****************************************

      return
      end subroutine distance_calc
