
      subroutine eigjac(maxdim,ndim,matrix,numlop,covval,wrk,eigvec,ier)

!***********************************************************************
!
!     EIGEN VALUE PROBLEM OF REAL SYMMETRIC MATRIX BY JACOBI METHOD
!       FOR 8 OR LESS DIMENSION MATRIX , JACOBI METHOD IS EFFICINET
!       BUT FOR LARGE DIMENSION MATRIX , JACOBI METHOD IS NOT GOOD
!       HOUSEHOLDER METHOD IS RECOMMENDED FOR LARGE DIMENSION MATRIX
!
!***********************************************************************
 
      implicit none

      integer(4),intent(in):: maxdim,ndim,numlop
      integer(4),intent(out):: ier
      ! convergence criteria
        real(8),intent(in):: covval
      ! eigen vector
        real(8),intent(out):: eigvec(maxdim,maxdim)
      real(8),intent(inout):: matrix(maxdim,maxdim),wrk(maxdim,2)

      integer(4):: i,j,imx,jmx,iloop,ierr
      real(8):: rmxele,rotang,tanval,cosval,sinval
 
!*********************************************
 
!     <<<  INITIAL SETTING  >>>
      if ( ndim .lt. 2 ) then
        ier = -1 ; return
      else
        ier = 0 ; ierr = 1
      endif
      eigvec(1:ndim,1:ndim) = 0.d0
      do i = 1,ndim
        eigvec(i,i) = 1.d0
      enddo

!     <<<  CALCULATION LOOP  >>>
      do iloop = 1,numlop
!       1) SEARCH BIGGEST ELEMENT (CHECK CONVERGENCE)
        rmxele = -1.d0
        do i = 1,ndim-1
        do j = i+1,ndim
          if ( abs(matrix(i,j)) .ge. rmxele ) then
            rmxele = abs( matrix(i,j) ) ; imx = i ; jmx = j
          endif
        enddo
        enddo
        if ( rmxele .le. covval ) then
          ierr = 0 ; exit
        endif

!       2) CALCULATE ROTATION ANGLE
        if ( matrix(imx,imx) .eq. matrix(jmx,jmx) ) then
          cosval = 1.d0 / sqrt(2.d0) ; sinval = cosval
        else
          tanval = 2.d0*matrix(imx,jmx) /                              &
                   (matrix(imx,imx)-matrix(jmx,jmx))
          rotang = atan(tanval) * 0.5d0
          cosval = cos(rotang) ; sinval = sin(rotang) 
        endif

!       3) SYMMETRIC OPERATION
        wrk(1:ndim,1) = matrix(imx,1:ndim)*cosval +                    &
                        matrix(jmx,1:ndim)*sinval
        wrk(1:ndim,2) = -matrix(imx,1:ndim)*sinval +                   &
                        matrix(jmx,1:ndim)*cosval
        wrk(imx,1) = matrix(imx,imx)*cosval*cosval +                   &
                     matrix(jmx,jmx)*sinval*sinval +                   &
                     2.d0*matrix(imx,jmx)*cosval*sinval 
        wrk(jmx,1) = 0.d0 ; wrk(imx,2) = 0.d0
        wrk(jmx,2) = matrix(imx,imx)*sinval*sinval +                   &
                     matrix(jmx,jmx)*cosval*cosval -                   &
                     2.d0*matrix(imx,jmx)*cosval*sinval 
        matrix(imx,1:ndim) = wrk(1:ndim,1)
        matrix(jmx,1:ndim) = wrk(1:ndim,2)
        matrix(1:ndim,imx) = matrix(imx,1:ndim)
        matrix(1:ndim,jmx) = matrix(jmx,1:ndim)

!       4) CALCULATE EIGEN VECTOR
        wrk(1:ndim,1) = eigvec(1:ndim,imx)*cosval +                    &
                          eigvec(1:ndim,jmx)*sinval
        wrk(1:ndim,2) = -eigvec(1:ndim,imx)*sinval +                   &
                          eigvec(1:ndim,jmx)*cosval 
        wrk(1:ndim,1) = eigvec(1:ndim,imx)*cosval +                    &
                        eigvec(1:ndim,jmx)*sinval
        wrk(1:ndim,2) = -eigvec(1:ndim,imx)*sinval +                   &
                        eigvec(1:ndim,jmx)*cosval
        eigvec(1:ndim,imx) = wrk(1:ndim,1)
        eigvec(1:ndim,jmx) = wrk(1:ndim,2)
      enddo

      if ( ierr .ne. 0 ) then
        ier = -2 ; return
      endif

!     <<<  SORTING  >>>
      call sorjac(maxdim,ndim,matrix,eigvec,wrk)

!**********************************

      return
      end subroutine eigjac


!================================================================


      subroutine sorjac(maxdim,ndim,matrix,eigvec,wrk)

      implicit none

      integer(4):: maxdim,ndim
      real(8):: matrix(maxdim,maxdim),eigvec(maxdim,maxdim),           &
                wrk(maxdim,2)

      integer(4):: i,j,isel
      real(8):: rmax

!******************************************

      do i = 1,ndim-1
        rmax = matrix(i,i)
        wrk(1:ndim,1) = eigvec(1:ndim,i)
        do j = i,ndim
          if ( matrix(j,j) .ge. rmax ) then
            isel = j ; rmax = matrix(j,j)
          endif
        enddo
        matrix(isel,isel) = matrix(i,i)
        matrix(i,i) = rmax
        eigvec(1:ndim,i) = eigvec(1:ndim,isel)
        eigvec(1:ndim,isel) = wrk(1:ndim,1)
      enddo

      return
      end subroutine sorjac
