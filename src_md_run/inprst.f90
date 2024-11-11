
      subroutine inprst(iprint,filena,ixnatm,iynvar,et,ek,ep,    &
                        mdloop,curtim,idfvsl,ier)

!********************************************************************
!
!     INPUT RESTART-FILE
!
!********************************************************************
 
      use COMCMM,only: cluster_method,lambda,lambda_v
      use COMCMMC,only: cord,vel

      integer(4),intent(in):: iprint,ixnatm,iynvar,idfvsl
      integer(4),intent(out):: mdloop,ier
      character(*),intent(in):: filena
      ! total,kinetic,potential energy at time T (KCAL/MOL) &
      ! current time (fsec)
      real(8),intent(out):: et,ek,ep,curtim

      character(80):: title
      integer(4):: numatm,numvar
      real(8):: rtmp

!*********************************************
 
      ier = 0
!     <<<  OPEN RESTART FILE  >>>
      call flopen(1,filena,20,'NULL',0,ier)
      if ( ier .ne. 0 ) then
        write(iprint,*)"ERROR> MD "
        write(iprint,*)"   RESTART FILE OPEN ERROR"
        ier = -1 ; return
      endif
 
!     <<<  read RESTART FILE  >>>
      read(1,err=9010,end=9020)title
      read(1,err=9010,end=9020)numatm,numvar
      if ( numatm .ne. ixnatm ) then
        write(iprint,*)"ERROR> MD "
        write(iprint,*)"   NUMBER OF ATOMS "
        write(iprint,*)"      SYSTEM       : ",ixnatm
        write(iprint,*)"      RESTART FILE : ",numatm
        ier = -4 ; return
      endif
      if ( numvar .ne. iynvar ) then
        write(iprint,*)"ERROR> MD "
        write(iprint,*)"   NUMBER OF FREE ATOMS "
        write(iprint,*)"      SYSTEM       : ",iynvar
        write(iprint,*)"      RESTART FILE : ",numvar
        ier = -5 ; return
      endif
      read(1,err=9010,end=9020)mdloop,curtim,et,ek,ep
      read(1,err=9010,end=9020)cord(1:3,1:numatm)
      read(1,err=9010,end=9020)vel(1:3,1:numvar)
      read(1,err=9010)rtmp,lambda_v
      if ( lambda .eq. 0.d0 ) then
        lambda = rtmp
      else
        write(iprint,*)
        write(iprint,*)"    The lambda or the square value for ALSD"
        write(iprint,*)"      is changed"
        write(iprint,'(2(a,f8.3))')"      from ",rtmp," to ",lambda
        write(iprint,*)
      endif
      if ( idfvsl .ne. 4 ) lambda_v = 0.d0
 
!     <<<  OUTPUT INFORMATION  >>>
      write(iprint,*)" "
      write(iprint,*)"INFORMATION> MD "
      write(iprint,*)"   INPUT RESTART FILE "
      write(iprint,'(A80)')title
      write(iprint,*)"    NUMBER OF ATOMS       : ",numatm
      write(iprint,*)"    NUMBER OF FREE ATOMS  : ",numvar
      write(iprint,*)"    LAST MD LOOP          : ",mdloop
      write(iprint,*)"    LAST MD TIME (FSEC)   : ",curtim
      write(iprint,*)"    NEXT STEP ENERGY        "
      write(iprint,*)"         TOTAL ENERGY     : ",et
      write(iprint,*)"         KINETIC ENERGY   : ",ek
      write(iprint,*)"         POTENTIAL ENERGY : ",ep
      if ( cluster_method .eq. "LMD" ) then
        write(iprint,*)"    LAMBDA FOR ALSD       : ",lambda
        write(iprint,*)"    LAMBDA VELOCITY       : ",lambda_v
      elseif ( cluster_method .eq. "LMD2" ) then
        write(iprint,*)"    LMDSQU FOR ALSD       : ",lambda
        write(iprint,*)"    LMDSQU VELOCITY       : ",lambda_v
      endif
      write(iprint,*)" "
 
!     <<<  CLOSE RESTART FILE  >>>
      call flclos(1,10,ier)
      if ( ier .ne. 0 ) then
        write(iprint,*)"ERROR> MD "
        write(iprint,*)"   RESTART FILE CLOSE ERROR"
        ier = -1 ; return
      endif
 
!*********************************

      return
!     <<<  errOR MANAGEMENT  >>>
9010  write(iprint,*)"ERROR> MD "
      write(iprint,*)"   DATA read ERROR IN readING RESTART FILE"
      ier = -2 ; return
9020  write(iprint,*)"ERROR> MD "
      write(iprint,*)"   END OF FILE IS DETECTED DURING"
      write(iprint,*)"   READING RESTART FILE "
      ier = -3 ; return

      end subroutine inprst
