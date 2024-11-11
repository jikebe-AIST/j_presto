
      subroutine minver(a,l,m,ierr)

      implicit none

      integer(4),intent(in):: l,m
      real(8),intent(inout):: a(l,m)
      integer(4),intent(out):: ierr

      integer(4):: i,j,k,lr,iw
      integer(4):: iwork(m)
      real(8):: wmax,w,pivot,api
!      real(8):: eps = 1.d-5
      real(8):: eps = 1.d-12
!      real(8):: eps = 1.d-14

!*******************************************************

!      if ( m .lt. 2 ) then
!        ierr = 9991
!        print*,"ERROR (minver) ierr = ",ierr
!        return
!      else
      ierr = 0
!      endif

      forall(i=1:m) iwork(i) = i

      do k = 1,m
        wmax = 0.d0
        do i = k,m
          w = abs(a(i,k))
          if ( w .gt. wmax ) then
            wmax = w ; lr = i
          endif
        enddo
        pivot = a(lr,k) ; api = abs(pivot)
        if ( lr .ne. k ) then
          iw = iwork(k) ; iwork(k) = iwork(lr) ; iwork(lr) = iw
          do j = 1,m
            w = a(k,j) ; a(k,j) = a(lr,j) ; a(lr,j) = w
          enddo
        endif

        ! For regular pivot
        if ( api .gt. eps ) then
          a(k,1:m) = a(k,1:m) / pivot
          do i = 1,m
            if ( i .ne. k ) then
              w = a(i,k)
              if ( w .ne. 0.d0 ) then
                do j = 1,m
                  if ( j .ne. k ) a(i,j) = a(i,j) - w * a(k,j)
                enddo
                a(i,k) = -w / pivot
              endif
            endif
          enddo
          a(k,k) = 1.d0 / pivot

        ! For small pivot
        else
          a(k,1:m) = 0.d0 ; a(1:m,k) = 0.d0
        endif

      enddo

      do i = 1,m
        do while ( iwork(i) .ne. i )
          k = iwork(i) ; iw = iwork(k)
          iwork(k) = iwork(i) ; iwork(i) = iw
          do j = 1,m
            w = a(j,i) ; a(j,i) = a(j,k) ; a(j,k) = w
          enddo
        enddo
      enddo

!******************************

      return
      end subroutine minver
