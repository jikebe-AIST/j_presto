
      program transform

!*******************************************************

      implicit none

      integer(4),parameter:: nBIN = 100

      integer(4):: Ndim,i,j
      real(8),allocatable:: coe(:),EE(:),dlnP(:)
      real(8):: lowd,upd,minE,maxE,temp,EbinSZ,rr
      character(130):: inpfil,outfil

!***************************

      read(5,'(a)')inpfil
      read(5,'(a)')outfil
      open(unit=1,file=trim(inpfil),status="old")
      read(1,*)Ndim ; allocate(coe(0:Ndim))
      do i = 0,Ndim
        read(1,*)coe(i)
      enddo
      read(1,*)lowd,upd
      read(1,*)minE,maxE
      read(1,*)temp
      close(1)

      EbinSZ = dble(maxE-minE) / dble(nBIN-1)
      allocate(EE(nBIN),dlnP(nBIN))
      forall ( i=1:nBIN ) EE(i) = minE + EbinSZ*dble(i-1)
      
      do i = 1,nBIN
        rr = coe(0)
        do j = 1,Ndim
          rr = rr + coe(j)*EE(i)**j
        enddo
        dlnP(i) = rr
      enddo

      EbinSZ = 1.d0 / dble(nBIN-1)
      forall ( i=1:nBIN ) EE(i) = EbinSZ*dble(i-1)
      call lstsquare(Ndim,nBIN,EE,dlnP,coe(0:Ndim))

      open(unit=1,file=trim(outfil),status="replace")
      write(1,*)1,Ndim
      do i = 0,Ndim
        write(1,*)coe(i)
      enddo
      write(1,*)0.d0,1.d0
      write(1,*)coe(0),sum(coe(0:Ndim))
      write(1,*)minE,maxE
      write(1,*)temp
      close(1)

!***************************

      stop
      end program transform
