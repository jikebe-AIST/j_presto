
      subroutine project

!*********************************************************
!
!     PCA projection
!
!*********************************************************

      use COMIFN
      !$ use omp_lib

      implicit none

      integer(4):: ifil,i,j,system,ndim
      real(4),allocatable:: cod(:),proj(:)
      real(8),allocatable:: av1(:),eval(:),evec(:,:)
      real(8):: wei

!********************************
!     Input axis

      open(unit=1,file=trim(AXISFL),status="old",form="unformatted")
      read(1)ndim
      allocate(av1(ndim),eval(ndim),evec(ndim,ndim),cod(ndim),         &
               proj(naxis))
      read(1)av1(1:ndim),eval(1:ndim),evec(1:ndim,1:ndim)
      close(1)

!********************************
!     Input file

      open(unit=3,file=trim(PROJNM)//".proj",status="replace")
      do ifil = 1,nfile
        write(6,'(2x,a)')"+ Now input "//trim(filelist(ifil))
        open(unit=2,file=trim(filelist(ifil)),status="old",            &
          form="unformatted")
        read(2)i
        if ( i .ne. ndim ) call error()
        ! Input data
        do
          read(2,end=801)wei,cod(1:ndim)
          cod(1:ndim) = cod(1:ndim) - av1(1:ndim)
          proj(1:naxis) = 0.0
          !$OMP parallel default (none)                            & !
          !$OMP private(i)                                         & !
          !$OMP shared(naxis,ndim,proj,cod,evec)
          !$OMP do
          do i = 1,naxis
            proj(i) = dot_product(cod(1:ndim),evec(1:ndim,i))
          enddo
          !$OMP end do
          !$OMP end parallel
          write(3,'(e12.5,$)')wei
          do i = 1,naxis
            write(3,'(x,e12.5,$)')proj(i)
          enddo
          write(3,*)
        enddo
801     close(2)
      enddo
      close(3)

!********************************

      return
      end subroutine project
