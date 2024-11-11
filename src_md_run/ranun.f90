
      subroutine ranun (n,iseq,x,ier)
 
!*******************************************************************
!
!     MAKE UNIFORM RANDOM NUMBER
!
!       THIS UNIFORM RANDOM NUMBER IS MADE BY LEHMER'S METHOD
!            ISEQ (INPUT) IS SEED OF RANDOM NUMBER
!            ISEQ(I+1) = MOD( (IA*ISEQ(I))+IB , IC )
!            X(I)      = FLOAT( ISEQ(I+1) ) / FLOAT( IC )
!       FOR PREVENTING OVER-FLOW
!            ISEQ(I) = ID*IU(I) + IL(I)
!            IL(I+1) = MOD( (IA*IL(I)+IB) , ID )
!            IU(I+1) = MOD( (IA*IU(I)+ILS) , IE )
!                ILS = (IA*IL(I)+IB) / ID
!                IE  = IC / ID
!
!*******************************************************************
 
      implicit none

      ! Number of random number
        integer(4),intent(in):: n
      ! Sequence for random number
        integer(4),intent(inout):: iseq
      ! Random number for 0-1
        real(4),intent(out):: x(n)
      integer(4),intent(out):: ier

      ! IA-E must be positive <- they are seed of random number
        integer(4):: ia,ib,ic,id,ie,iu,il,ild,ils,i
      real(8):: rc
 
!***********************************************************

      if ( iseq.lt.0 .or. n.lt.1 ) then
        ier = -1 ; return
      else
        ier = 0
      endif

      ia = 32771 ; ib = 123456781
      ic = 2**30 ; id = 2**15 ; ie = 2**15

      iu = iseq / id
      il = iseq - iu*id
      rc = 1.0 / float(ic)

      do i = 1,n
        ild = il*ia + ib
        ils = ild / id
        il = ild - ils*id
        iu = mod(ia*iu+ils,ie)
        iseq = id*iu + il
        x(i) = float(iseq) * rc
      enddo

!*******************************

      return
      end subroutine ranun
