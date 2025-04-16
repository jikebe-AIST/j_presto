
      subroutine mkaxis

!*********************************************************
!
!     make PCA axis
!
!*********************************************************

      use COMIFN

      implicit none

      integer(4):: ifil,i,j,system,ndim,ier,iopt
      integer(4),allocatable:: iflg(:),iwk(:)
      real(4),allocatable:: cod(:)
      real(8),allocatable:: av1(:),av2(:,:),eval(:),evec(:,:),work(:,:)
      real(8):: wei,sumwei,esum,eps

!********************************
!     Input file

      sumwei = 0.d0
      do ifil = 1,nfile
        write(6,'(2x,a)')"+ Now input "//trim(filelist(ifil))
        open(unit=2,file=trim(filelist(ifil)),status="old",            &
             form="unformatted")
        read(2)ndim
        ! Initialization
        if ( .not. allocated(cod) ) then
          allocate(cod(ndim),av1(ndim),av2(ndim,ndim))
          av1(1:ndim) = 0.d0 ; av2(1:ndim,1:ndim) = 0.d0
        endif
        ! Input data
        do
          read(2,end=801)wei,cod(1:ndim)
          sumwei = sumwei + wei
          av1(1:ndim) = av1(1:ndim) + cod(1:ndim)*wei
          do j = 1,ndim
          do i = j,ndim
            av2(i,j) = av2(i,j) + cod(i)*cod(j)*wei
          enddo
          enddo
        enddo
801     continue ; close(2)
      enddo

      ! Calc. variance-covariance matrix
      av1(1:ndim) = av1(1:ndim) / sumwei
      av2(1:ndim,1:ndim) = av2(1:ndim,1:ndim) / sumwei
      do j = 1,ndim
      do i = j,ndim
        av2(i,j) = av2(i,j) - av1(i)*av1(j)
      enddo
      enddo
      do i = 1,ndim-1
      do j = i+1,ndim
        av2(i,j) = av2(j,i)
      enddo
      enddo

      ! Calc. eigenvalues & vectors
      allocate(eval(ndim),evec(ndim,ndim),iflg(ndim),work(ndim,6),     &
               iwk(ndim))
      eps = -1.0 ; iopt = 1
      call def2m(av2,ndim,ndim,ndim,ndim,eps,iopt,eval,evec,iflg,work, &
                 iwk,ier)

      ! Monitering
      esum = 100.d0/sum(eval)
      write(6,*)
      write(6,'(2x,a)')"+ engenvalue & eigenvalue/SUM (%)"
      write(6,*)
      do i = 1,ndim
        write(6,'(8x,i8,x,2(f12.5,x),a)')i,eval(i),eval(i)*esum," (%)"
      enddo

      ! Output data
      open(unit=1,file=trim(PROJNM)//".axis",status="replace",         &
           form="unformatted")
      write(1)ndim
      write(1)av1(1:ndim),eval(1:ndim),evec(1:ndim,1:ndim)

!********************************

      return
      end subroutine mkaxis
