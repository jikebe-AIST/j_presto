
      subroutine dihedral_calc

!**************************************************
!
!     dihedral angle calculation
!
!**************************************************

      use COMVAL ; use COMIFN

      implicit none

      integer(4):: i,i1,i2,i3,i4
      real(8):: d21(3),d32(3),d43(3),p12(3),p23(3),r12,r23,r12r23,     &
                s1223,cosp,p123(3),s32123,phi

!*******************************************************************

      write(utc,'(f,$)')wfac
      do i = 1,n_dih_c
        i1 = dih_c(1,i) ; i2 = dih_c(2,i)
        i3 = dih_c(3,i) ; i4 = dih_c(4,i)
        d21(1:3) = cod(1:3,i2) - cod(1:3,i1)
        d32(1:3) = cod(1:3,i3) - cod(1:3,i2)
        d43(1:3) = cod(1:3,i4) - cod(1:3,i3)
        p12(1) = d21(2)*d32(3) - d32(2)*d21(3)
        p12(2) = d21(3)*d32(1) - d32(3)*d21(1)
        p12(3) = d21(1)*d32(2) - d32(1)*d21(2)
        p23(1) = d32(2)*d43(3) - d43(2)*d32(3)
        p23(2) = d32(3)*d43(1) - d43(3)*d32(1)
        p23(3) = d32(1)*d43(2) - d43(1)*d32(2)
        r12 = p12(1)*p12(1) + p12(2)*p12(2) + p12(3)*p12(3)
        r23 = p23(1)*p23(1) + p23(2)*p23(2) + p23(3)*p23(3)
        r12r23 = 1.d0 / sqrt(r12*r23)
        s1223 = p12(1)*p23(1) + p12(2)*p23(2) + p12(3)*p23(3)
        cosp = s1223 * r12r23
        if ( cosp .gt. 1.d0 ) then
          cosp = 1.d0
        elseif ( cosp .lt. -1.d0 ) then
          cosp = -1.d0
        endif
        phi = acos(cosp)
        p123(1) = p12(2)*p23(3) - p23(2)*p12(3)
        p123(2) = p12(3)*p23(1) - p23(3)*p12(1)
        p123(3) = p12(1)*p23(2) - p23(1)*p12(2)
        s32123 = d32(1)*p123(1) + d32(2)*p123(2) + d32(3)*p123(3)
        if ( s32123 .lt. 0.d0 ) phi = 2.d0*npi - phi
        phi = phi * rpi
        if ( phi .gt. 180.d0 ) phi = phi - 360.d0
        write(utc,'(x,f8.3,$)')phi
      enddo
      write(utc,*)

!*****************************************

      return
      end subroutine dihedral_calc
