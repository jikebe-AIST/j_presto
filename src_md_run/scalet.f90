
      subroutine scalet(temps,tempc,taut,deltt,iprint,iynvar,lamda)

!********************************************************************
!
!     SCALING VELOCITY FOR CONSTANT TEMERATURE
!
!********************************************************************
 
      use COMCMMC,only: vel

      integer(4),intent(in):: iprint,iynvar
      real(8),intent(in):: temps,tempc,taut,deltt
      real(8),intent(out):: lamda
 
!********************************************************
 
!     <<<  CHECK CURRENT TEMPERATURE  >>>
      if ( tempc .le. 0.0D0 ) then
        write(iprint,*)" "
        write(iprint,*)"INFORMATION> SCALET "
        write(iprint,*)"  CURRENT TEMPERATURE IS LESS EQUAL "
        write(iprint,*)"  ZERO "
        write(iprint,*)"  NOW WE DO NOT SCALE VELOCITY FOR "
        write(iprint,*)"  CONSTANT TEMPERATURE "
        write(iprint,*)" " ; return
      endif
 
!     <<<  CALCULATION SCALING FACTOR OF VELOCITY  >>>
      lamda  = 1.0d0 + ( ( deltt/taut ) * ( ( temps/tempc ) - 1.0d0 ) )
      lamda  = sqrt( lamda )
      if ( lamda .ge. 2.0d0 ) then
        write(iprint,*)" "
        write(iprint,*)"INFORMATION> SCALET "
        write(iprint,*)"  SCALING FACTOR OF VELOCITY IS GREATER"
        write(iprint,*)"  EQUAL 2.0 "
        write(iprint,*)"    CURRENT FACTOR IS ",lamda
        write(iprint,*)"    THEN SET FACTOR 2.0 "
        write(iprint,*)" "
        lamda = 2.0D0
      endif
 
!     <<<  SCALING VELOCITY  >>>
      vel(1:3,1:iynvar) = vel(1:3,1:iynvar) * lamda

!****************************************

      return
      end subroutine scalet
