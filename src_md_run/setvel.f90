
      subroutine setvel(istart,gascst,iseed,tempst,iynvar,iylvar,      &
                        maxatm,random,fxmass,ier)
 
!******************************************************
!
!     SET INITIAL VELOCITY FOR START OF MD RUN
!
!******************************************************

      use COMCMMC,only: vel
      !$ use omp_lib

      implicit none

      ! 1 = cold start, 2 = hot start
        integer(4),intent(in):: istart
      ! Gas constant
        real(8),intent(in):: gascst
      ! Random seed
        integer(4),intent(inout):: iseed
      ! Hot start temperature
        real(8),intent(in):: tempst
      integer(4),intent(in):: iynvar,iylvar(iynvar),maxatm
      real(4),intent(out):: random(3,maxatm)
      real(8),intent(in):: fxmass(maxatm)
      integer(4),intent(out):: ier

      integer(4):: i,j
      real(8):: rtmp

!*************************************************

      ier = 0
      if ( istart .eq. 1 ) then
        !$OMP parallel default (none)                                & !
        !$OMP private(i)                                             & !
        !$OMP shared(vel,iynvar)
        !$OMP do schedule (static)
        do i = 1,iynvar
          vel(:,i) = 0.d0
        enddo
        !$OMP end do
        !$OMP end parallel
      elseif ( istart .eq. 2 ) then
        call ranno(3*iynvar,iseed,0.0,1.0,random,ier)
        if ( ier .ne. 0 ) return
        rtmp = sqrt(tempst*gascst*1.0d-7)
        !$OMP parallel default (none)                                & !
        !$OMP private(i,j)                                           & !
        !$OMP shared(iynvar,iylvar,vel,random,rtmp,fxmass)
        !$OMP do schedule (static)
        do i = 1,iynvar
          j = iylvar(i)
          vel(1:3,i) = random(1:3,i) * rtmp / sqrt(fxmass(j))
        enddo
        !$OMP end do
        !$OMP end parallel
      endif

!************************************

      return
      end subroutine setvel
