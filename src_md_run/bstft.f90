
      subroutine bstft(numatm,refcrd,comcrd,weight,iflag,rotmat,cenref,&
                       cencom,rmsd,ier)

!***********************************************************************
!
!     LEAST SQUARE FITTING OF TWO CONFORMATIONS
!        THIS SUBROUTINE USES EIGJAC SUBROUTINE (JACOBI'S DIAGONALIZATIO
!        REFERENCE (ACTA CRYST. (1989).A45,208-210)
!
!***********************************************************************
 
      implicit none

      ! iflag : flag for fitting of comcrd
        integer(4),intent(in):: numatm,iflag
      integer(4),intent(out):: ier
      real(8),intent(in):: refcrd(3,numatm),weight(numatm)
      real(8),intent(out):: rotmat(3,3),cenref(3),cencom(3),rmsd
      real(8),intent(inout):: comcrd(3,numatm)

      integer(4):: i
      real(8):: totmas,dr(3),dc(3),dm(3),dp(3),symmat(4,4),wrk(4,2),   &
                eigvec(4,4),q(4)
 
!****************************************

      if ( numatm .gt. 0 ) then
        ier = 0
      else
        ier = -1 ; return
      endif

!     1) CALCULATION OF MASS CENTER
      cenref(1:3) = 0.d0 ; cencom(1:3) = 0.d0
      totmas = sum(weight(1:numatm))
      do i = 1,numatm
        cenref(1:3) = cenref(1:3) + weight(i)*refcrd(1:3,i)
      enddo
      do i = 1,numatm
        cencom(1:3) = cencom(1:3) + weight(i)*comcrd(1:3,i)
      enddo
      if ( totmas .gt. 0.d0 ) then
        totmas = 1.d0 / totmas
        cenref(1:3) = cenref(1:3) * totmas
        cencom(1:3) = cencom(1:3) * totmas
      else
        ier = -2 ; return
      endif

!     2) CALCULATION OF SYMMETRIC MATRIX
      symmat(1:4,1:4) = 0.d0
      do i = 1,numatm
        dr(1:3) = refcrd(1:3,i) - cenref(1:3)
        dc(1:3) = comcrd(1:3,i) - cencom(1:3)
        dm(1:3) = weight(i) * (dr(1:3) - dc(1:3))
        dp(1:3) = weight(i) * (dr(1:3) + dc(1:3))
        symmat(1,1) = symmat(1,1) + dot_product(dm,dm)
        symmat(1,2) = symmat(1,2) + dp(2)*dm(3) - dm(2)*dp(3)
        symmat(1,3) = symmat(1,3) + dm(1)*dp(3) - dp(1)*dm(3)
        symmat(1,4) = symmat(1,4) + dp(1)*dm(2) - dm(1)*dp(2)
        symmat(2,2) = symmat(2,2) + dp(2)**2 + dp(3)**2 + dm(1)**2
        symmat(2,3) = symmat(2,3) + dm(1)*dm(2) - dp(1)*dp(2)
        symmat(2,4) = symmat(2,4) + dm(1)*dm(3) - dp(1)*dp(3)
        symmat(3,3) = symmat(3,3) + dp(1)**2 + dp(3)**2 + dm(2)**2
        symmat(3,4) = symmat(3,4) + dm(2)*dm(3) - dp(2)*dp(3)
        symmat(4,4) = symmat(4,4) + dp(1)**2 + dp(2)**2 + dm(3)**2
      enddo
      symmat(2,1) = symmat(1,2) ; symmat(3,1) = symmat(1,3)
      symmat(3,2) = symmat(2,3) ; symmat(4,1) = symmat(1,4)
      symmat(4,2) = symmat(2,4) ; symmat(4,3) = symmat(3,4)
 
!     3) DIAGONALIZATION
!           SYMMAT(4,4) IS SMALLEST EIGEN VALUE WHICH GIVES RITATION
!           WHICH MINIMIZE THE SUM OF THE DISTANCES
!           SYMMAT(1,1) IS BIGGEST EIGEN VALUE WHICH GIVES RITATION
!           WHICH MAXMISE THE SUM OF THE DISTANCES
      call eigjac(4,4,symmat,1000,1.0d-10,wrk,eigvec,ier)
      if ( ier .ne. 0 ) then
        ier = -3 ; return
      endif
 
!     4) MAKE ROTATION MATRIX
      q(1:4) = eigvec(1:4,4)
      rotmat(1,1) = q(1)**2 + q(2)**2 - q(3)**2 - q(4)**2
      rotmat(1,2) = 2.d0 * (q(2)*q(3) + q(1)*q(4))
      rotmat(1,3) = 2.d0 * (q(2)*q(4) - q(1)*q(3))
      rotmat(2,1) = 2.d0 * (q(2)*q(3) - q(1)*q(4))
      rotmat(2,2) = q(1)**2 - q(2)**2 + q(3)**2 - q(4)**2
      rotmat(2,3) = 2.d0 * (q(3)*q(4) + q(1)*q(2))
      rotmat(3,1) = 2.d0 * (q(2)*q(4) + q(1)*q(3))
      rotmat(3,2) = 2.d0 * (q(3)*q(4) - q(1)*q(2))
      rotmat(3,3) = q(1)**2 - q(2)**2 - q(3)**2 + q(4)**2

!     5) CALCULATION OF RMSD AND FITTING
      rmsd = 0.d0
      if ( iflag .eq. 1 ) then
        do i = 1,numatm
          dr(1:3) = refcrd(1:3,i) - cenref(1:3)
          dc(1:3) = comcrd(1:3,i) - cencom(1:3)
          dp(1) = dot_product(rotmat(1,1:3),dc)
          dp(2) = dot_product(rotmat(2,1:3),dc)
          dp(3) = dot_product(rotmat(3,1:3),dc)
          comcrd(1:3,i) = dp(1:3) + cenref(1:3)
          dp(1:3) = dr(1:3) - dp(1:3)
          rmsd = rmsd + dot_product(dp,dp)
        enddo
      else
        do i = 1,numatm
          dr(1:3) = refcrd(1:3,i) - cenref(1:3)
          dc(1:3) = comcrd(1:3,i) - cencom(1:3)
          dp(1) = dr(1) - dot_product(rotmat(1,1:3),dc)
          dp(2) = dr(2) - dot_product(rotmat(2,1:3),dc)
          dp(3) = dr(3) - dot_product(rotmat(3,1:3),dc)
          rmsd = rmsd + dot_product(dp,dp)
        enddo
      endif
      rmsd = sqrt(rmsd/dble(numatm))

!*********************************

      return
      end subroutine bstft
