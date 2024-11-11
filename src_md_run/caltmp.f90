
      subroutine caltmp(iynvar,iylvar,vel,fxmass,numfre,temp,ek,l_vel, &
                        l_mas)

!*********************************************************
!
!     Calculation of current temperature
!
!********************************************************* 

      use COMBAS, only: ixnatm
      use PHYCNS, only: gascst,joucal
      !$ use omp_lib

      implicit none

      integer(4),intent(in):: iynvar,iylvar(iynvar),numfre
      real(8),intent(in):: vel(3,ixnatm),fxmass(ixnatm),l_vel,l_mas
      real(8),intent(out):: temp,ek

      integer(4):: i,iatm
      
!***********************************
 
!     <<<  CALCULATION OF KINETIC ENERGY  >>>
      ek = l_mas*l_vel*l_vel
      !$OMP parallel default(none)                                   & !
      !$OMP private(i,iatm)                                          & !
      !$OMP shared(iynvar,iylvar,fxmass,vel)                         & !
      !$OMP reduction(+ : ek)
      !$OMP do schedule (static)
      do i = 1,iynvar
        iatm = iylvar(i)
        ek = ek + fxmass(iatm) * ( vel(1,iatm)*vel(1,iatm) +           &
             vel(2,iatm)*vel(2,iatm) + vel(3,iatm)*vel(3,iatm) )
      enddo
      !$OMP end do
      !$OMP end parallel
      ek = ek * 1.0d+7
 
!     <<<  CALCULATION OF TEMERATURE  >>>
      temp = ek / (gascst * dble(numfre))
      ek = ek / joucal * 0.0005

!************************************

      return
      end subroutine caltmp


!===================================================================


      subroutine caltmp2(iynvar,iylvar,vel,fxmass,nfrmx,tempf,nfrmol,  &
                         l_vel,l_mas)

!**********************************************************      
!
!     CALCULATION OF CURRENT TEMERATURE
!     FOR EACH MOLECULE (MOLECULE 1 & MOLECULE2).
!
!**********************************************************
 
      use COMBAS, only: ixnatm,ixamol2
      use PHYCNS, only: gascst,joucal
      use COMMOL, only: ifg_lambda
      use COMCMM, only: l_temp_control
      use CALC_TIME
      !$ use omp_lib

      implicit none

      integer(4),intent(in):: iynvar,iylvar(iynvar),nfrmx
      real(8),intent(in):: vel(3,ixnatm),fxmass(ixnatm),nfrmol(nfrmx), &
                           l_vel,l_mas
      real(8),intent(out):: tempf(nfrmx)

      integer(4):: i,iatm,imol

!**********************************************

!     <<<  CALCULATION OF KINETIC ENERGY  >>>
      tempf(1:nfrmx) = 0.d0

      !$OMP parallel default(none)                                   & !
      !$OMP private(i,iatm,imol)                                     & !
      !$OMP shared(iynvar,iylvar,ixamol2,fxmass,vel)                 & !
      !$OMP reduction(+ : tempf)
      !$OMP do schedule (static)
      do i = 1,iynvar
        iatm = iylvar(i) ; imol = ixamol2(iatm)
        tempf(imol) = tempf(imol) + fxmass(iatm) * (                   &
                      vel(1,iatm)*vel(1,iatm) +                        &
                      vel(2,iatm)*vel(2,iatm) +                        &
                      vel(3,iatm)*vel(3,iatm) )
      enddo
      !$OMP end do
      !$OMP end parallel

      if ( l_temp_control )                                            &
        tempf(ifg_lambda) = tempf(ifg_lambda) + l_mas*l_vel*l_vel

!     <<<  CALCULATION OF TEMERATURE  >>>
      i = ixamol2(iynvar)
      tempf(1:i) = tempf(1:i) * 1.0d+7 / (nfrmol(1:i)*gascst)

!****************************************

      return
      end subroutine caltmp2
