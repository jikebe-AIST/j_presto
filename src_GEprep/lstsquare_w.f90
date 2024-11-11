
      subroutine lstsquare_w(mxNdim,Ndat,X,Y,wei,svAnswer,fNdim)

!***********************************************
!
!     LeaST SQUARE method with Weighing factor
!
!***********************************************

      implicit none

      ! Number of DIMension
        integer(4),intent(in):: mxNdim
      ! Number of DATa
        integer(4),intent(in):: Ndat
      ! Data for fitting
        real(8),intent(in):: X(Ndat),Y(Ndat)
      ! Weighting factor
        real(8),intent(in):: wei(Ndat)
      ! Answer
        real(8),intent(out):: svAnswer(mxNdim+1)
      ! fitted dimension
        integer(4),intent(out):: fNdim

      real(8):: Xsum(mxNdim*2),XYsum(mxNdim),Ysum,Answer(mxNdim+1)
      real(8),allocatable:: A(:,:),B(:),C(:)

      integer(4):: Ndim,i,j,ierr
      real(8):: det,prea0,a0,ex,Resierr,pastResierr,swei,rsum

!*******************************
      ! Initialization
      swei = 1.d0 / sum(wei) ; pastResierr = 999999999.d0
      svAnswer(1:mxNdim+1) = 0.d0

      ! Start fitting cycle
      do Ndim = 2,mxNdim
        Xsum(:) = 0.d0 ; XYsum(:) = 0.d0 ; Ysum = 0.d0
        ! Summation
        do i = 1,Ndat
          do j = 1,Ndim*2
            Xsum(j) = Xsum(j) + X(i)**j * wei(i)
          enddo
          do j = 1,Ndim
            XYsum(j) = XYsum(j) + Y(i) * X(i)**j * wei(i)
          enddo
          Ysum = Ysum + Y(i)*wei(i)
        enddo
        ! Mean calculation
        Xsum(:) = Xsum(:) * swei
        XYsum(:) = XYsum(:) * swei
        Ysum = Ysum * swei
        if ( allocated(A) ) deallocate(A,B,C)
        allocate(A(Ndim,Ndim),B(Ndim),C(Ndim))
        A(1:Ndim,1:Ndim) = 0.d0 ; B(1:Ndim) = 0.d0

        ! Making matrixes A & B
        do i = 1,Ndim
          do j = 1,Ndim
            A(i,j) = Xsum(i+j) - Xsum(i)*Xsum(j)
          enddo
          B(i) = XYsum(i) - Ysum*Xsum(i)
        enddo

        ! Making inverse matrix A
        call Minver(A,Ndim,Ndim,det,ierr)

        ! Making a coefficient matrix
        C = matmul(A,B)
          
        ! Coefficient
        a0 = 0.d0
        do i = 1,Ndim
          prea0 = -C(i) * Xsum(i)
          a0 = a0 + prea0
        enddo
        a0 = a0 + Ysum
        Answer(1) = a0 ; Answer(2:Ndim+1) = C(1:Ndim)

        ! Residual error calculation
        Resierr = 0.d0
        do i = 1,Ndat
          ex = 0.d0
          do j = 1,Ndim+1
            ex = ex + (Answer(j)*X(i)**(j-1))
          enddo
          ex = Y(i) - ex
          ex = ex*ex*wei(i)
          Resierr = Resierr + ex
        enddo
        Resierr = sqrt(Resierr*swei)
        if ( Resierr .lt. pastResierr ) then
          fNdim = Ndim
          pastResierr = Resierr
          svAnswer(1:Ndim+1) = Answer(1:Ndim+1)
          write(6,'(4x,i0,a,f)')Ndim,"-D Fitting Error = ",pastResierr
        endif
      enddo

      if ( fNdim .ne. mxNdim ) svAnswer(fNdim+2:mxNdim+1) = 0.d0
      do i = fNdim,1,-1
        if ( svAnswer(i+1) .ne. 0.d0 ) then
          j = i ; exit
        endif
      enddo
      fNdim = j

      write(6,*)
      write(6,'(4x,a,i0)')"Final fitting dimension = ",fNdim

!***********************************

      return
      end subroutine lstsquare_w
