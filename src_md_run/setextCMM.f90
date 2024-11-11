
      subroutine setextCMM(iprint,cieCMn)

!***************************************************
!
!     Input parameter for external CMM
!
!***************************************************

      use COMCMM

      implicit none

      integer(4),intent(in):: iprint
      character(*),intent(in):: cieCMn

      integer(4):: i
      real(4):: r

!*****************************

      write(iprint,*)""
      write(iprint,*)"INFORMATION> INPEXTCMM"

      ! flag for cell shift is changed to .true.
      extCMM_flag = .true.

      open(unit=1,file=cieCMn,status="old",form="unformatted")
      
      read(1)r,i,celx,cely,celz
      write(6,'(6x,a26,i8)')   "External CMM level      : ",i
      write(6,'(6x,a26,f8.3)') "Minimul level cell size : ",r
      write(6,'(6x,a26,3f8.3)')"Cell center coordinate  : ",           &
                               celx,cely,celz
      if ( i .ne. 4 ) then
        write(6,*)"Sorry, this PGM. can move only cell lvl. 4"
        stop
      endif

      if ( r .ne. siz(1) ) then
        write(6,*)" ERROR : The minimum cell size in NAMECE and ",     &
                  "NAMEEC file are different"
        write(6,*)"Cell size in NAMECE : ",siz(1)
      endif
      if ( i .ne. nlev ) then
        write(6,*)" ERROR : nlev in NAMECE and in NAMEEC file are ",   &
                  "different"
        write(6,*)"nlev in NAMECE : ",nlev
      endif

!*************************************

      allocate(exflag1(nv(1),nv(1),nv(1)),exflag2(nv(2),nv(2),nv(2)))
      allocate(exflag3(nv(3),nv(3),nv(3)),exflag4(nv(4),nv(4),nv(4)))
      allocate(eCMMpot1(13,nv(1),nv(1),nv(1)))
      allocate(eCMMpot2(13,nv(2),nv(2),nv(2)))
      allocate(eCMMpot3(13,nv(3),nv(3),nv(3)))
      allocate(eCMMpot4(13,nv(4),nv(4),nv(4)))

!*************************************

      ! external cell calculation flag
      read(1)exflag1(1:nv(1),1:nv(1),1:nv(1))
      read(1)exflag2(1:nv(2),1:nv(2),1:nv(2))
      read(1)exflag3(1:nv(3),1:nv(3),1:nv(3))
      read(1)exflag4(1:nv(4),1:nv(4),1:nv(4))

      ! parameters
      do i = 1,10
        read(1)eCMMpot1(i,1:nv(1),1:nv(1),1:nv(1))
      enddo
      eCMMpot1(11:13,1:nv(1),1:nv(1),1:nv(1)) =                        &
              eCMMpot1(8:10,1:nv(1),1:nv(1),1:nv(1)) * 2.d0
      eCMMpot1(8:10,1:nv(1),1:nv(1),1:nv(1)) =                         &
              eCMMpot1(5:7,1:nv(1),1:nv(1),1:nv(1)) * 2.d0

      do i = 1,10
        read(1)eCMMpot2(i,1:nv(2),1:nv(2),1:nv(2))
      enddo
      eCMMpot2(11:13,1:nv(2),1:nv(2),1:nv(2)) =                        &
              eCMMpot2(8:10,1:nv(2),1:nv(2),1:nv(2)) * 2.d0
      eCMMpot2(8:10,1:nv(2),1:nv(2),1:nv(2)) =                         &
              eCMMpot2(5:7,1:nv(2),1:nv(2),1:nv(2)) * 2.d0

      do i = 1,10
        read(1)eCMMpot3(i,1:nv(3),1:nv(3),1:nv(3))
      enddo
      eCMMpot3(11:13,1:nv(3),1:nv(3),1:nv(3)) =                        &
              eCMMpot3(8:10,1:nv(3),1:nv(3),1:nv(3)) * 2.d0
      eCMMpot3(8:10,1:nv(3),1:nv(3),1:nv(3)) =                         &
              eCMMpot3(5:7,1:nv(3),1:nv(3),1:nv(3)) * 2.d0

      do i = 1,10
        read(1)eCMMpot4(i,1:nv(4),1:nv(4),1:nv(4))
      enddo
      eCMMpot4(11:13,1:nv(4),1:nv(4),1:nv(4)) =                        &
              eCMMpot4(8:10,1:nv(4),1:nv(4),1:nv(4)) * 2.d0
      eCMMpot4(8:10,1:nv(4),1:nv(4),1:nv(4)) =                         &
              eCMMpot4(5:7,1:nv(4),1:nv(4),1:nv(4)) * 2.d0

      close(1)

!*****************************
!     For nearest cell calc.

      ! radius for nearest cell calc. of external CMM
      RR_ext = 1.5d0 * siz(1)
      RR_ext_2 = 3.d0 * siz(1)
      RR_ext__2 = RR_ext * RR_ext
      iRR_ext = 1.d0 / RR_ext
      iRR_ext__3 = iRR_ext * iRR_ext * iRR_ext
      iRR_ext__5 = iRR_ext__3 * iRR_ext * iRR_ext
      i3RR_ext__5 = 3.d0 * iRR_ext__5
      i5RR_ext__7 = 5.d0 * iRR_ext__5 * iRR_ext * iRR_ext

      ! nearest cell calc. table
      i = nv(1)*nv(1)*nv(1)
      allocate(mycel_nEC(3,i))
      allocate(Nyrcel_nEC(i),yrcel_nEC(3,2*icmx+1,i))

!*****************************

      return
      end subroutine setextCMM
