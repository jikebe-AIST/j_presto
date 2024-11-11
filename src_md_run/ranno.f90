
      subroutine ranno(n,iseq,mean,sd,x,ier)

!*******************************************************************
!
!     MAKE NORMALIZED RANDOM NUMBER
!     THIS NORMALIZED RANDOM NUMBER IS MADE BY BOX-MULLAR METHOD
!
!*******************************************************************
 
      implicit none

      ! Number of random number
        integer(4),intent(in):: n
      ! Sequence for random number (iseq must be positive number)
        integer(4),intent(inout):: iseq
      ! Mean & SD
        real(4),intent(in):: mean,sd
      ! Random number for 0 - 1
        real(4),intent(out):: x(n)
      ! Condition code
        integer(4),intent(out):: ier

      real(4),parameter:: phi2 = 2.0 * 3.141593
      real(4):: addx(2),dx1,dx2,rtmp,rtmp2
      integer(4):: iersub,i

!**************************************************

      if ( iseq.lt.0 .or. n.lt.1 ) then
        ier = -1 ; return
      else
        ier = 0
      endif

      call ranun(n,iseq,x,iersub)
      if ( mod(n,2) .ne. 0 ) call ranun(2,iseq,addx,iersub)
      if ( iersub .ne. 0 ) then
        ier = -2 ; return
      endif

      if ( mod(n,2) .eq. 0 ) then
        do i = 1,n-1,2
          dx1 = x(i) ; dx2 = x(i+1)
          rtmp = sd * sqrt(-2.0*log(dx1))
          rtmp2 = phi2 * dx2
          x(i) = rtmp * cos(rtmp2) + mean
          x(i+1) = rtmp * sin(rtmp2) + mean
        enddo
      else
        do i = 1,n-2,2
          dx1 = x(i) ; dx2 = x(i+1)
          rtmp = sd * sqrt(-2.0*log(dx1))
          rtmp2 = phi2 * dx2
          x(i) = rtmp * cos(rtmp2) + mean
          x(i+1) = rtmp * sin(rtmp2) + mean
        enddo
        dx1 = x(n) ; dx2 = addx(1)
        x(n) = sd * sqrt(-2.0*log(dx1)) * cos(phi2*dx2) + mean
      endif

!***************************************

      return
      end subroutine ranno
