
      subroutine SIMD_check(i)

      use COMCMM

      implicit none

      integer(4),intent(in):: i

      integer(4):: tt
      integer(8),save:: t1,t2,cr,cm,mint
      logical(4),save:: fstF = .true.

!*************************************

      select case (i)
      case (1)
        call system_clock(t1)
      case (2)
        call system_clock(t2,cr,cm)
        tt = t2 - t1
        if ( tt .le. 0 ) tt = cm + tt
        if ( fstF ) then
          mint = tt ; fstF = .false.
        else
          if ( tt .lt. mint ) then
            mint = tt ; i_vec = i_vec * 2
          else
            if ( i_vec .ne. 1 ) i_vec = i_vec / 2
            SIMD_chk = .false.
            write(6,*)
            write(6,'(a,i0)')                                         &
              "  SIMD PARALLELAZATION was set to ",i_vec
            write(6,*)
          endif
        endif
      end select

!*************************************

      return
      end subroutine SIMD_check
