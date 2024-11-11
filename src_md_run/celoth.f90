
      subroutine celoth

!************************************************************
!
!  Get N of atoms in cells at each level.
!
!************************************************************

      use COMPAR ; use COMCMM
      use COMBAS,only: ixnatm

      implicit none

      integer(4):: ix,iy,iz,jx1,jx2,jy1,jy2,jz1,jz2,i

!**************************************

      ! Get N of atoms in level-2 cells from those in level 1.
      !$OMP parallel default (none)                                  & !
      !$OMP private(ix,iy,iz,jx1,jx2,jy1,jy2,jz1,jz2)                & !
      !$OMP shared(nv,lvdn,nfrag,npcl2,npcl1)
      !$OMP do schedule (static)
      do iz = 1,nv(2)
      jz1 = lvdn(1,iz) ; jz2 = lvdn(2,iz)
      do iy = 1,nv(2)
      jy1 = lvdn(1,iy) ; jy2 = lvdn(2,iy)
      do ix = 1,nv(2)
      jx1 = lvdn(1,ix) ; jx2 = lvdn(2,ix)
        forall ( i=1:nfrag )                                           &
          npcl2(i,ix,iy,iz) = sum(npcl1(i,jx1:jx2,jy1:jy2,jz1:jz2))
      enddo
      enddo
      enddo
      !$OMP end do
      !$OMP end parallel

      if ( nlev .lt. 3 ) return
      !  Get N of atoms in level-3 cells from those in level 2.
      do iz = 1,nv(3)
      jz1 = lvdn(1,iz) ; jz2 = lvdn(2,iz)
      do iy = 1,nv(3)
      jy1 = lvdn(1,iy) ; jy2 = lvdn(2,iy)
      do ix = 1,nv(3)
      jx1 = lvdn(1,ix) ; jx2 = lvdn(2,ix)
        forall ( i=1:nfrag )                                           &
          npcl3(i,ix,iy,iz) = sum(npcl2(i,jx1:jx2,jy1:jy2,jz1:jz2))
      enddo
      enddo
      enddo

      if ( nlev .lt. 4 ) return
      !  Get N of atoms in level-4 cells from those in level 3.
      do iz = 1,nv(4)
      jz1 = lvdn(1,iz) ; jz2 = lvdn(2,iz)
      do iy = 1,nv(4)
      jy1 = lvdn(1,iy) ; jy2 = lvdn(2,iy)
      do ix = 1,nv(4)
      jx1 = lvdn(1,ix) ; jx2 = lvdn(2,ix)
        forall ( i=1:nfrag )                                           &
          npcl4(i,ix,iy,iz) = sum(npcl3(i,jx1:jx2,jy1:jy2,jz1:jz2))
      enddo
      enddo
      enddo

      if ( nlev .lt. 5 ) return
      ! Get N of atoms in level-5 cells from those in level 4.
      do iz = 1,nv(5)
      jz1 = lvdn(1,iz) ; jz2 = lvdn(2,iz)
      do iy = 1,nv(5)
      jy1 = lvdn(1,iy) ; jy2 = lvdn(2,iy)
      do ix = 1,nv(5)
      jx1 = lvdn(1,ix) ; jx2 = lvdn(2,ix)
        forall ( i=1:nfrag )                                           &
          npcl5(i,ix,iy,iz) = sum(npcl4(i,jx1:jx2,jy1:jy2,jz1:jz2))
      enddo
      enddo
      enddo

!************************************************************

      return
      end subroutine celoth
