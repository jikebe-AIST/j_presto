
      subroutine calc_pot(Tslv_cod,Tslt_cod,slt_chg,ion_order,ion_chg, &
                          edist,ion)

      implicit none

      integer(4), intent(in):: ion_order(:)
      real(8), intent(in):: Tslv_cod(:),Tslt_cod(:),slt_chg(:),        &
                            ion_chg(:),edist
      integer(4), intent(inout):: ion(:)

      integer(4):: i,j,n,nslv,nslt,iloc(1),iatm,jatm
      real(8):: dist,d
      real(8),allocatable:: slv_cod(:,:),slt_cod(:,:), pot(:)
      logical(1), allocatable:: flag(:)

!**********************************************************************

      ! Make arrays
      nslv = size(Tslv_cod)/3 ; nslt = size(Tslt_cod)/3
      allocate(slv_cod(3,nslv), slt_cod(3,nslt))
      n = 0
      do j = 1, nslv
      do i = 1, 3
        n = n + 1 ; slv_cod(i,j) = Tslv_cod(n)
      enddo
      enddo
      n = 0
      do j = 1, nslt
      do i = 1, 3
        n = n + 1 ; slt_cod(i,j) = Tslt_cod(n)
      enddo
      enddo

      ! Calculate electric fields at each solvent molecule
      allocate(pot(nslv),flag(nslv))
      pot(1:nslv) = 0.d0 ; flag(1:nslv) = .true.
      !$OMP parallel default(none)                                   & !
      !$OMP private(i,j)                                             & !
      !$OMP shared(nslv,nslt,pot,slt_chg,slv_cod,slt_cod)
      !$OMP do
      do i = 1,nslv
      do j = 1,nslt
        pot(i) = pot(i) + slt_chg(j) / dist(slv_cod(:,i),slt_cod(:,j))
      enddo
      enddo
      !$OMP end do
      !$OMP end parallel

      do i = 1,size(ion_order)
        ! Put an ion
        if ( ion_chg(i) .lt. 0.d0 ) then
          iloc = maxloc(pot,mask=flag)
        else
          iloc = minloc(pot,mask=flag)
        endif
        iatm = iloc(1) ; ion(iatm) = ion_order(i) ; flag(iatm) = .false.
        ! Exclude neighboring solvent molecules from exchange with ions
        do jatm = 1,nslv
          if ( flag(jatm) ) then
            d = dist(slv_cod(:,iatm),slv_cod(:,jatm))
            if ( d .le. edist ) flag(jatm) = .false.
            !! Update electric fields
            pot(jatm) = pot(jatm) + ion_chg(i) / d
          endif
        enddo
      enddo

!**********************************************************************

      return
      end subroutine calc_pot


!**********************************************************************


      function dist(cod1, cod2)

      implicit none

      real(8),intent(in):: cod1(3), cod2(3)
      real(8):: dist

      real(8):: dx, dy, dz

!**********************************************************************

      dx = cod1(1) - cod2(1) ; dx = dx * dx
      dy = cod1(2) - cod2(2) ; dy = dy * dy
      dz = cod1(3) - cod2(3) ; dz = dz * dz
      dist = sqrt(dx+dy+dz)

!**********************************************************************

      return
      end function dist
