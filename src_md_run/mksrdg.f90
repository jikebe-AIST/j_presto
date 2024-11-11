
      subroutine mksrdg(numvar,listv,srdir)
      
!*******************************************************************
!
!     MAKE SEARCH DIRECTION
!     SEARCH DIRECTION VECTOR IS - (NORMALIZED GRADIENT VECTOR)
!
!*******************************************************************

      use COMBAS,only: ixnatm
      use COMCMMC,only: grad
      use COMCMM,only: nfrg

      implicit none

      integer(4),intent(in):: numvar,listv(ixnatm)
      real(8),intent(inout):: srdir(3,ixnatm)

      integer(4):: ivar,iatm
      real(8):: conv 
 
!*****************************************

      conv = 0.d0
      do ivar = 1,numvar
        iatm = listv(ivar)
        conv = conv + grad(1,iatm,nfrg)*grad(1,iatm,nfrg) +            &
                      grad(2,iatm,nfrg)*grad(2,iatm,nfrg) +            &
                      grad(3,iatm,nfrg)*grad(3,iatm,nfrg)
      enddo
      conv = -1.d0 / sqrt(conv)

      do ivar = 1,numvar
        iatm = listv(ivar)
        srdir(1:3,ivar) = conv*grad(1:3,iatm,nfrg)
      enddo

!********************************

      return
      end subroutine mksrdg
