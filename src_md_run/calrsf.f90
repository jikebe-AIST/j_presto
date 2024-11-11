
      subroutine calrsf(iynvar,iylvar,grad,rmsf,maxf,locf)

!*******************************************************************
!
!     THIS SUBROUTINE IS FOR CALCULATION OF ROOT MEAN
!     SQUARE SUMMANTION OF FORCE
!
!*******************************************************************

      use COMBAS, only: ixnatm

      implicit none

      ! Number of free atoms & List of free atom number
        integer(4),intent(in):: iynvar,iylvar(iynvar)
      ! Gradient
        real(8),intent(in):: grad(3,ixnatm)
      ! R.M.S.F
        real(8),intent(out):: rmsf
      ! Max force 
        real(8),intent(out):: maxf
      ! Atom number, which has max force
        integer(4),intent(out):: locf

      integer(4):: ivar,iatm
      real(8):: rtmp
 
!*********************************************

      rmsf = 0.d0 ; maxf = -9999.d0
      do ivar = 1,iynvar
        iatm = iylvar(ivar)
        rtmp = grad(1,iatm)*grad(1,iatm) +                             &
               grad(2,iatm)*grad(2,iatm) + grad(3,iatm)*grad(3,iatm)
        rmsf = rmsf + rtmp
        if ( rtmp .gt. maxf ) then
          maxf = rtmp ; locf = iatm
        endif
      enddo
      maxf = sqrt(maxf)
      rmsf = sqrt(rmsf / dble(iynvar))

!*********************************

      return
      end subroutine calrsf
