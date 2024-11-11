
      subroutine simulu(m,n,A,B,X,epsl,iflag,ip,ier)
 
!***********************************************************************
!
!     Solve linear simultaneous equations A*X = B
!       by using LU decomposition
!
!***********************************************************************

      implicit none

      ! Max number of dimension
        integer(4),intent(in):: m
      ! Number of dimension
        integer(4),intent(in):: n
      ! Left hand coefficient matrix, this is LU decomposed
        real(8),intent(inout):: A(m,m)
      ! Right hand coefficient vector
        real(8),intent(in):: B(m)
      ! Solution
        real(8),intent(out):: X(m)
      ! Machine epsilon if diagonal element is less than epsl,
      !   then this matrix is singular
        real(8),intent(in):: epsl
      ! Flag for calc. of LU-decomposition 
      ! (0: NOT LU-decomposition, 1: Apply LU-decomposition)
        integer(4),intent(in):: iflag
      ! Index vector for partial pivotiong (0: (I), 1: (O))
        integer(4),intent(inout):: ip(m)
      ! Condition code (0: NO ERROR, -1: matrix A is singular)
        integer(4),intent(inout):: ier
 
      integer(4):: i,j
      real(8):: t
 
!************************************************

      ier = 0
      ! LU-decomposition
      if ( iflag .eq. 1 ) then
        call ludeco(m,n,A,epsl,ip,ier)
        if ( ier .ne. 0 ) return
      endif
 
      ! Forward substitution
      do i = 1,n
        t = B(ip(i))
        do j = 1,i-1
          t = t - A(ip(i),j)*X(j)
        enddo
        X(i) = T
      enddo

      ! Backward substitution
      do i = n,1,-1
        t = X(i)
        do j = i+1,N
          t = t - A(ip(i),j)*X(j)
        enddo
        X(i) = t * A(IP(i),i)
      enddo

!**************************

      return
      end subroutine simulu

!================================================================ 


      subroutine ludeco(m,n,A,epsl,ip,ier)

      implicit none

      integer(4),intent(in):: m,n
      real(8),intent(inout):: A(m,m)
      real(8),intent(in):: epsl
      integer(4),intent(inout):: ip(m)
      integer(4),intent(inout):: ier

      integer(4):: lv,i,j,k,l
      real(8):: al
 
!*******************************************

      ier = 0
      ! Prepar index vector for partial pivoting
      do k = 1,n
        ip(k) = k
      enddo

      ! LU decompositon
      do k = 1,n
        !! Pivoting
        l = k
        al = abs(A(ip(l),k))
        do i = k+1,n
          if ( abs(A(ip(i),k)) .gt. al ) then
            l = i
            al = abs(A(ip(l),k))
          endif
        enddo
        if ( l .ne. k ) then
          lv = ip(k)
          ip(k) = ip(l)
          ip(l) = lv
        endif

        !!! In case of singular matrix
        if ( abs(A(ip(k),k)) .lt. epsl ) then
          ier = -1 ; return
        endif

        !!! Gauss elimination
        A(ip(k),k) = 1.d0 / A(ip(k),k)
        do i = k+1,n
          A(ip(i),k) = A(ip(i),k) * A(ip(k),k)
          do j = k+1,n
            A(ip(i),j) = A(ip(i),j) - A(ip(i),k)*A(ip(k),j)
          enddo
        enddo

      enddo

!**************************

      return
      end subroutine ludeco 
