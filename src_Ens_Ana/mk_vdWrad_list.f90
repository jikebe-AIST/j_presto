
      subroutine mk_vdWrad_list(nATM,nlist,ilist,RESnm,ATMnm,vdWrad)

      use COMTPL

      implicit none

      integer(4),intent(in):: nATM,nlist
      integer(4),intent(in):: ilist(nlist)
      character(4),intent(in):: RESnm(nATM)
      character(4),intent(in):: ATMnm(nATM)
      real(8),intent(out):: vdWrad(nlist)

      integer(4):: i,ic,ii,j
      character(4):: res1,atm1,res2,atm2
      logical(1):: check

!***********************************************************************

      ic = 1
      do i = 1,nlist
        check = .false. ; ii = ilist(i)
        res1 = adjustl(RESnm(ii)) ; atm1 = adjustl(ATMnm(ii))
        do j = ic,maxatm
          res2 = adjustl(cxresn(j)) ; atm2 = adjustl(cxatmn(j))
          if ( res1.eq.res2 .and. atm1.eq.atm2 ) then
             vdWrad(i) = fxvdwr(j) ; ic = j + 1 ; check = .true. ; exit
          endif
        enddo
        if ( .not. check ) then
          write(6,*)RESnm(ii),ATMnm(ii) ; call error(10601)
        endif
      enddo
      write(6,'(4x,a)')"Radius info. in the input tpl file is used"

!      do i = 1,nlist
!        ii = ilist(i)
!        print*,i,ii
!        print*,RESnm(ii),ATMnm(ii),vdWrad(i)
!      enddo

!***********************************************************************

      return
      end subroutine mk_vdWrad_list


!========================================================================


      subroutine mk_vdWrad_list2(nATM,nlist,ilist,RESnm,ATMnm,rad,f)

      implicit none

      integer(4),intent(in):: nATM,nlist
      integer(4),intent(in):: ilist(nlist)
      character(4),intent(in):: RESnm(nATM)
      character(4),intent(in):: ATMnm(nATM)
      real(8),intent(out):: rad(nlist)
      logical(1),intent(in):: f

      ! Atom name for vdW radii
        integer(4),parameter:: nANvdW = 11
      ! Atom name array for vdW radii
        character(4),parameter:: ATMnmvdW(nANvdW) = (/                 &
          "CL  ","BR  ","FE  ","C   ","N   ","O   ","S   ","H   ",     &
          "P   ","F   ","I   "/)
      ! Number of atom name characters for vdW radii
        integer(4),parameter:: iANvdW(nANvdW) = (/                     &
          2,2,2,1,1,1,1,1,1,1,1/)
      ! vdW radii array (RvdW = R/2**(1/6))
      real(8),parameter:: vdWrad(nANvdW) = (/                          &
        1.735d0, & ! CL
        1.978d0, & ! BR
        1.26d0,  & ! FE
        1.7d0  , & ! C
        1.625d0, & ! N
        1.48d0 , & ! O
        1.782d0, & ! S
        1.d0   , & ! H
        1.871d0, & ! P
        1.560d0, & ! F
        2.094d0  & ! I
        /)

      integer(4):: i,j,k,l,ierr
      character(4):: tmp

!***********************************************************************

      do i = 1,nlist
        j = ilist(i) ; ierr = 0
        tmp = trim(adjustl(ATMnm(j)))
        do k = 1,nANvdW
          l = iANvdW(k)
          if ( tmp(1:l) .eq. ATMnmvdW(k)(1:l) ) then
            rad(i) = vdWrad(k) ; ierr = 0 ; exit
          endif
          if ( ierr .ne. 0 ) then
            write(6,'(4x,a)')"! ERROR !"
            write(6,'(4x,a)')"In contact_init, the vdW radius "//&
              "parameter for "//tmp//" is NOT defined "
            stop
          endif
        enddo
      enddo
      if ( f ) then
        write(6,'(4x,a)')"The below atom radii are used"
        do k = 1,nANvdW
          write(6,'(8x,a,x,f8.3)')ATMnmvdW(k),vdWrad(k)
        enddo
      endif

!***********************************************************************

      return
      end subroutine mk_vdWrad_list2
