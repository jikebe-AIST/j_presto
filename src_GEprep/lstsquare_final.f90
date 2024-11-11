
      subroutine lstsquare_final(Nwindow,Ndat,Ndim,minv,maxv,X,Y,wei,  &
                   COE,low_v_window,high_v_window,Resierr,Resierr2,oflg)

!******************************************************
!
!     least square method with weighting factor
!
!******************************************************

      implicit none

      ! Number of windows, data, & dimensions
        integer(4),intent(in):: Nwindow,Ndat,Ndim
      ! Next min & max energy or lambda
        real(8),intent(in):: minv,maxv
      ! Data for fitting
        real(8),intent(in):: X(Ndat),Y(Ndat),wei(Ndat)
      ! Coefficient
        real(8),intent(out):: COE(0:Ndim,Nwindow)
      ! Boundary
        real(8),intent(out):: low_v_window(Nwindow),                   &
                              high_v_window(Nwindow)
      ! Resierr
        real(8),intent(out):: Resierr,Resierr2
      ! Output flag
        logical(4),intent(in):: oflg

      integer(4):: itole,Nmxbin(Nwindow),i,j,k,ii,ierr,maxWIN,imx,imn
      real(8):: Rrange,rtmp,rtmp2,Xsum(0:Ndim*2),XYsum(0:Ndim),swei,   &
                A(Ndim,Ndim),B(Ndim),C(Ndim),det,dx,ex,ex2,            &
                Acon(0:Ndim+2,0:Ndim+2),Bcon(0:Ndim+2),Ccon(0:Ndim+2), &
                Ndata_window(Nwindow)

!*******************************

      ! Initialization
      Rrange = (maxv - minv) / dble(Nwindow)
      rtmp = Rrange + minv ; j = 1 ; rtmp2 = Rrange*0.00001d0
      low_v_window(1) = minv
      do i = 1,Ndat
        if ( rtmp-X(i) .lt. rtmp2 ) then
          high_v_window(j) = X(i)
          Nmxbin(j) = i
          j = j + 1 ; rtmp = Rrange*dble(j) + minv
          if ( j .eq. Nwindow+1 ) exit
          low_v_window(j) = X(i)
        endif
      enddo
      high_v_window(Nwindow) = maxv
      Nmxbin(Nwindow) = Ndat

      ! Tolerance number of bins for fitting
      itole = Ndat / (2*Nwindow)

!*******************************

      ! Search the window included the maximum number of data
      Ndata_window = 0.d0 ; k = 1
      do i = 1,Nwindow
        do j = k,Nmxbin(i)
          Ndata_window(i) = Ndata_window(i) + wei(j)*wei(j)
        enddo
        k = Nmxbin(i) + 1
      enddo
      maxWIN = 1
      do i = 2,Nwindow
        if ( Ndata_window(i) .gt. Ndata_window(maxWIN) ) maxWIN = i
      enddo

      ! normal least square for the maxWIN'th window
      if ( maxWIN .ne. 1 ) then
        ii = Nmxbin(maxWIN-1)
      else
        ii = 1
      endif
      imn = max(1,ii-itole) ; imx = min(Nmxbin(maxWIN)+itole,Ndat)
      swei = 1.d0 / sum(wei(imn:imx))
      Xsum(0:Ndim*2) = 0.d0 ; XYsum(0:Ndim) = 0.d0
      ! Summation
      do i = imn,imx
        dx = X(i)-low_v_window(maxWIN)
        do j = 1,Ndim*2
          Xsum(j) = Xsum(j) + dx**j * wei(i)
        enddo
        do j = 0,Ndim
          XYsum(j) = XYsum(j) + Y(i) * dx**j * wei(i)
        enddo
      enddo
      ! Mean calculation
      Xsum(1:Ndim*2) = Xsum(1:Ndim*2) * swei
      XYsum(0:Ndim) = XYsum(0:Ndim) * swei

      ! Making matrixes A & B
      forall ( i=1:Ndim, j=1:Ndim ) A(i,j) = Xsum(i+j) - Xsum(i)*Xsum(j)
      forall ( i=1:Ndim ) B(i) = XYsum(i) - XYsum(0)*Xsum(i)

      ! Making inverse matrix A
      call minver(A,Ndim,Ndim,det,ierr)

      ! Making a coefficient matrix
      C = matmul(A,B)

      ! Coefficient
      COE(0,maxWIN) = -dot_product(C,Xsum(1:Ndim)) + XYsum(0)
      COE(1:Ndim,maxWIN) = C(1:Ndim)

!***********************************

      ! For the follow windows
      do i = maxWIN+1,Nwindow
        ii = i - 1
        dx = high_v_window(ii) - low_v_window(ii)
        ! Calc. coefficients constrainted
        COE(0,i) = COE(0,ii)
        do j = 1,Ndim
          COE(0,i) = COE(0,i) + COE(j,ii)*dx**j
        enddo
        COE(1,i) = COE(1,ii)
        do j = 2,Ndim
          COE(1,i) = COE(1,i) + dble(j)*COE(j,ii)*dx**(j-1)
        enddo

        imn = max(1,Nmxbin(i-1)-itole) ; imx = min(Nmxbin(i)+itole,Ndat)
        swei = 1.d0 / sum(wei(imn:imx))
        Xsum(0:Ndim*2) = 0.d0 ; XYsum(0:Ndim) = 0.d0
        ! Summation
        do j = imn,imx
          dx = X(j)-low_v_window(i)
          do k = 1,Ndim*2
            Xsum(k) = Xsum(k) + dx**k * wei(j)
          enddo
          do k = 0,Ndim
            XYsum(k) = XYsum(k) + Y(j) * dx**k * wei(j)
          enddo
        enddo
        ! Mean calculation
        Xsum(1:Ndim*2) = Xsum(1:Ndim*2) * swei
        XYsum(0:Ndim) = XYsum(0:Ndim) * swei

        ! Making matrixes A & B
        forall ( j=2:Ndim, k=2:Ndim ) A(j-1,k-1) = Xsum(j+k)
        forall ( j=2:Ndim ) B(j-1) = XYsum(j) - COE(0,i)*Xsum(j)       &
                                              - COE(1,i)*Xsum(j+1)

        ! Making inverse matrix A
        call minver(A(1:Ndim-1,1:Ndim-1),Ndim-1,Ndim-1,det,ierr)

        ! Making a coefficient matrix
        C(1:Ndim-1) = matmul(A(1:Ndim-1,1:Ndim-1),B(1:Ndim-1))

        ! Coefficient
        COE(2:Ndim,i) = C(1:Ndim-1)
      enddo

      ! The previous windows
      do i = maxWIN-1,1,-1
        if ( i .ne. 1 ) then
          ii = Nmxbin(i-1)
        else
          ii = 1
        endif
        imn = max(1,ii-itole) ; imx = min(Nmxbin(i)+itole,Ndat)
        swei = 1.d0 / sum(wei(imn:imx))
        Xsum(0:Ndim*2) = 0.d0 ; XYsum(0:Ndim) = 0.d0
        ! Summation
        do j = imn,imx
          dx = X(j)-low_v_window(i)
          do k = 1,Ndim*2
            Xsum(k) = Xsum(k) + dx**k * wei(j)
          enddo
          do k = 0,Ndim
            XYsum(k) = XYsum(k) + Y(j) * dx**k * wei(j)
          enddo
        enddo
        ! Mean calculation
        Xsum(0) = 1.d0
        Xsum(1:Ndim*2) = Xsum(1:Ndim*2) * swei
        XYsum(0:Ndim) = XYsum(0:Ndim) * swei

        ! Making matrixes A & B
        ii = i + 1
        dx = high_v_window(i) - low_v_window(i)
        forall ( j=0:Ndim, k=0:Ndim )                                  &
          Acon(j,k) = Xsum(j+k)
        forall ( j=0:Ndim ) Acon(j,Ndim+1) = -(dx**j)*swei*0.5d0
        Acon(0,Ndim+2) = 0.d0
        forall ( j=1:Ndim )                                            &
          Acon(j,Ndim+2) = -(dble(j)*dx**(j-1))*swei*0.5d0
        forall ( j=0:Ndim ) Acon(Ndim+1,j) = dx**j
        forall ( j=0:Ndim ) Acon(Ndim+2,j) = dble(j)*dx**(j-1)
        Acon(Ndim+1:Ndim+2,Ndim+1:Ndim+2) = 0.d0

        forall ( j=0:Ndim ) Bcon(j) = XYsum(j)
        Bcon(Ndim+1:Ndim+2) = COE(0:1,ii)
        
        ! Making inverse matrix A
        call minver(Acon,Ndim+3,Ndim+3,det,ierr)

        ! Making a coefficient matrix
        Ccon = matmul(Acon,Bcon)

        ! Coefficient
        COE(0:Ndim,i) = Ccon(0:Ndim)
      enddo

!***************************************************

      ! Residual error calc.
      Resierr = 0.d0 ; Resierr2 = 0.d0
      !! The first window
      do i = 1,Nmxbin(1)
        ex = 0.d0 ; dx = X(i)-low_v_window(1)
        do j = 0,Ndim
          ex = ex + COE(j,1)*dx**j
        enddo
        ex = Y(i) - ex
        ex2 = ex*ex
        ex = ex2*wei(i)
        ex2 = ex*ex2
        Resierr = Resierr + ex
        Resierr2 = Resierr2 + ex2
      enddo
      !! The last window ( Nwindow > 1 )
      if ( Nwindow .gt. 1 ) then
        do i = Nmxbin(Nwindow-1)+1,Ndat
          ex = 0.d0 ; dx = X(i)-low_v_window(Nwindow)
          do j = 0,Ndim
            ex = ex + COE(j,Nwindow)*dx**j
          enddo
          ex = Y(i) - ex
          ex2 = ex*ex
          ex = ex2*wei(i)
          ex2 = ex*ex2
          Resierr = Resierr + ex
          Resierr2 = Resierr2 + ex2
        enddo

        !!! ( Nwindow > 2 )
        if ( Nwindow .gt.2 ) then
          do i = 2,Nwindow-1
          do j = Nmxbin(i-1)+1,Nmxbin(i)
            ex = 0.d0 ; dx = X(j)-low_v_window(i)
            do k = 0,Ndim
              ex = ex + COE(k,i)*dx**k
            enddo
            ex = Y(j) - ex
            ex2 = ex*ex
            ex = ex2*wei(j)
            ex2 = ex*ex2
            Resierr = Resierr + ex
            Resierr2 = Resierr2 + ex2
          enddo
          enddo
        endif
      endif

      swei = 1.d0 / sum(wei)
      Resierr = Resierr*swei
      Resierr2 = sqrt(Resierr2*swei)
      Resierr2 = sqrt(Resierr2-Resierr)
      Resierr = sqrt(Resierr)
      if ( oflg ) write(6,'(4x,a,f)')"Fitting Error = ",Resierr

!********************************************

      return
      end subroutine lstsquare_final
