
      subroutine inmcmd(iprint,iread,ier)

!********************************************************
!
!     Read parameters for multi canonical MD from a file
!
!********************************************************

      use COMBAS ; use COMPSC
      use COMCMM,only: cluster_method

      implicit none

      integer(4),intent(in):: iprint,iread
      integer(4),intent(out):: ier

      integer(4):: i,j

!**********************************

      read(iread,*,err=999)nwindow,ndeg
      allocate(c(0:ndeg,nwindow),low_v(nwindow),high_v(nwindow))
      do i = 1,nwindow
        do j = 0,ndeg
          read(iread,*,err=999)c(j,i)
        enddo
        read(iread,*,err=999)low_v(i),high_v(i)
      enddo
      read(iread,*,err=999)alphalw,alphaup
      read(iread,*,err=999)celw,ceup
      read(iread,*,err=999)tempce
      if ( cluster_method(1:3) .ne. "LMD" ) then
        write(iprint,*)"INFORMATION> INPUT PARAMETERS ",               &
                       "FOR MULTI-CANONICAL MD"
      else
        write(iprint,*)"INFORMATION> INPUT PARAMETERS ",               &
                       "FOR ADAPTIVE LAMBDA SQUARE DYNAMICS"
        read(iread,*,err=999)lower,upper
      endif
!!      write(iprint,*)ndeg
!!      write(iprint,*)c(1:ndeg+1)
!!      write(iprint,*)alphalw,alphaup
!!      write(iprint,*)celw,ceup
!!      write(iprint,*)tempce
      write(iprint,*)
      ier = 0 ; return

999   write(iprint,*)"ERROR> READ FILE ERROR IN INMCMD."
      ier = -1

!*********************************************

      return
      end subroutine inmcmd
