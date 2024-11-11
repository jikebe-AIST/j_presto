
      subroutine lstsquare(Ndim,Ndat,X,Y,Answer)

!***********************************************
!
!     LeaST SQUARE method
!
!***********************************************

      implicit none

      ! Number of DIMension
        integer(4),intent(in):: Ndim
      ! Number of DATa
        integer(4),intent(in):: Ndat
      ! Data for fitting
        real(8),intent(in):: X(Ndat),Y(Ndat)
      ! Answer
        real(8),intent(out):: Answer(Ndim+1)
 
      real(8):: Xsum(Ndim*2),XYsum(Ndim),Ysum
      real(8):: A(Ndim,Ndim),B(Ndim),C(Ndim)

      integer(4):: i,j,ierr
      real(8):: det,prea0,a0,ex,Resierr,exsum

!*******************************
      ! Initialization
      Xsum(:) = 0.d0 ; XYsum(:) = 0.d0 ; Ysum = 0.d0

      ! Summation
      do i = 1,Ndat
        do j = 1,Ndim*2
          Xsum(j) = Xsum(j) + X(i)**j
        enddo
        do j = 1,Ndim
          XYsum(j) = XYsum(j) + Y(i) * X(i)**j
        enddo
        Ysum = Ysum + Y(i)
      enddo

!     Mean calculation
      Xsum(:) = Xsum(:) / dble(Ndat)
      XYsum(:) = XYsum(:) / dble(Ndat)
      Ysum = Ysum / dble(Ndat)

!     Making matrixes A & B
      do i = 1,Ndim
        do j = 1,Ndim
          A(i,j) = Xsum(i+j) - Xsum(i)*Xsum(j)
        enddo
        B(i) = XYsum(i) - Ysum*Xsum(i)
      enddo

!     Making inverse matrix A
      call Minver(A,Ndim,Ndim,det,ierr)

!     Making a coefficient matrix
      C = matmul(A,B)
          
!     Coefficient
      a0 = 0.d0
      do i = 1,Ndim
        prea0 = -C(i) * Xsum(i)
        a0 = a0 + prea0
      enddo
      a0 = a0 + Ysum
      Answer(1) = a0

      do i = 1,Ndim
        Answer(i+1) = C(i)
      enddo

!     Residual error calculation

      Resierr = 0.d0 ; exsum = 0.d0
      do i = 1,Ndat
        ex = Y(i)
        exsum = exsum + ex*ex
        do j = 1,Ndim+1
          ex = ex - Answer(j)*X(i)**(j-1)
        enddo
        ex = ex*ex
        Resierr = Resierr + ex
      enddo

      exsum = sqrt(Resierr/exsum)*100.d0
      write(6,'(4x,a,f12.3,a)')"Error rate = ",exsum," (%)"

      Resierr = sqrt(Resierr)/dble(Ndat)
      write(6,'(4x,a,e15.6)')"Fitting Error = ",Resierr

!************************************

!      write(Nfil,*)Ndim,Resierr
!      do i = 1,Ndim+1
!        write(Nfil,*)Answer(i)
!      enddo

!      write(Nfil,*)
!      do i = 1,Ndat
!        ex = 0.d0
!        do j = 1,Ndim+1
!          ex = ex + Answer(j)*X(i)**(j-1)
!        enddo
!        write(Nfil,*)X(i),Y(i),ex
!      enddo

!***********************************

      return
      end subroutine lstsquare
