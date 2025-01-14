
      subroutine intmp(FILnam)

!********************************************************
!
!     INput residue number and name for template PDB
!
!********************************************************

      use COMVAL

      implicit none

      ! Input FILe NAMe
      character(*),intent(in):: FILnam

      integer(4):: i,j
      character(1):: c
      character(6):: chkatm
      character(4):: tmp
      character(130):: line
      logical(1),allocatable:: chk(:)

!**********************************
!     assign dimension
      i = 0
      open(unit=1,file=trim(FILnam),status="old")
      do
        read(1,*,end=799)chkatm
        if ( chkatm.eq."ATOM  " .or. chkatm.eq."HETATM" ) i = i + 1
      enddo
799   close(1)
      allocate(ATMnum(i),ATMnm(i),RESnm(i),RESnum(i),CHN(i),cod(3,i),  &
               Rcod(3,i),occupy(i),bfact(i),nCHN(i))
      nATM = i ; occupy = 0.d0 ; bfact = 0.d0

!     Input template PDB file
      open(unit=1,file=trim(FILnam),status="old") ; i = 0
      do
        read(1,'(a)',end=800)line
        if ( line(1:6).ne."ATOM  " .and. line(1:6).ne."HETATM") cycle
        i = i + 1
        read(line,'(6x,i6,a4,x,a4,a1,i8,3f8.3,2f6.3)')ATMnum(i),       &
          ATMnm(i),RESnm(i),CHN(i),RESnum(i),Rcod(1:3,i),occupy(i),    &
          bfact(i)
      enddo
800   close(1)
      fstRES = minval(RESnum) ; fnlRES = maxval(RESnum)

      c = CHN(1) ; j = 1 ; nCHN(1) = j
      do i = 2,nATM
        if ( CHN(i) .ne. c ) then
          j = j + 1 ; c = CHN(i)
        endif
        nCHN(i) = j
      enddo

!***********************************
!     output for check

      write(6,*)
      write(6,'(4x,a)')"+ N of atoms in template_PDB = "
      write(6,'(6x,i0)')nATM
      write(6,'(4x,a)')"+ Residues in template_PDB : "
      write(6,'(6x,i0,a,i0)')fstRES," - ",fnlRES
      allocate(aRES(fstRES:fnlRES),chk(fstRES:fnlRES))
      chk(:) = .false. ; aRES(:) = "    "
      do i = 1,nATM
        tmp = RESnm(i) ; j = RESnum(i)
        if ( chk(j) .and. tmp.ne.RESnm(i) ) call error(10301)
        chk(j) = .true. ; aRES(j) = tmp
      enddo

      write(6,'(4x,a)')"+ Sequence :"
      j = 0
      do i = fstRES,fnlRES
        j = j + 1
        if ( mod(j,10) .eq. 1 ) then
          write(6,'(6x,a4,$)')aRES(i)
        else
          write(6,'(x,a4,$)')aRES(i)
        endif
        if ( mod(j,10) .eq. 0 ) write(6,*)
      enddo
      if ( mod(j,10) .ne. 0 ) write(6,*)

!***********************************

      return
      end subroutine intmp
