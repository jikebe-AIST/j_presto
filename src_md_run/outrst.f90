
      subroutine outrst(iprint,filena,ixnatm,iynvar,cord,vel,et,&
                        ek,ep,mdloop,curtim,lambda,lambda_v,ier)

!*****************************************************************
!
!     Output restart-file for verlet
!
!*****************************************************************

      implicit none

      ! Logical unit number for input & output
        integer(4),intent(in):: iprint
      ! Restart file name
        character(80),intent(in):: filena
      ! (Max) number of atoms
        integer(4),intent(in):: ixnatm
      ! Number of free atoms
        integer(4),intent(in):: iynvar
      ! Coordinates at time T & verocity at time T-(1/2)dT
        real(8),intent(in):: cord(3,ixnatm),vel(3,ixnatm)
      ! Total, Kinetic, & Potential energy at time T (KCAL/mol)
        real(8),intent(in):: et,ek,ep
      ! Current time (fsec)
        real(8),intent(in):: curtim
      ! Current loop number
        integer(4),intent(in):: mdloop
      ! Lambda & lambda_v for adaptive lambda dynamics
        real(8),intent(in):: lambda,lambda_v
      ! Condition code ( 0: NO-ERROR, -1: File open error, 
      !                 -2: Data write error)
        integer(4),intent(inout):: ier

      character(80):: title
      character(23):: time
      character(40):: user
      integer(4):: lenusr
 
!******************************************

      ier = 0
      call flopen(1,filena,22,'NULL',0,ier)
      if ( ier .ne. 0 ) then
        ier = -1
        write(iprint,*)'ERROR> MD '
        write(iprint,*)'   RESTART FILE OPEN ERROR '
        return
      endif

      ! Make job-information
      call infjob(time,user,lenusr)
      title = '     '//'VERLET'//time//user(1:20)
 
      ! Write restart file
      write(1,err=800)title
      write(1,err=800)ixnatm,iynvar
      write(1,err=800)mdloop,curtim,et,ek,ep
      write(1,err=800)cord(1:3,1:ixnatm)
      write(1,err=800)vel(1:3,1:iynvar)
      write(1,err=800)lambda,lambda_v

      call flclos(1,10,ier)
      if ( ier .ne. 0 ) then 
        ier = -1
        write(iprint,*)'ERROR> MD '
        write(iprint,*)'   RESTART FILE CLOSE ERROR '
        return
      endif

      return

800   ier = -2
      write(iprint,*)"ERROR> MD "
      write(iprint,*)"   DATA WRITE ERROR IN OUTPUT RESTART FILE"
      
!****************************************

      return
      end subroutine outrst
