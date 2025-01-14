
      subroutine inref(FILnam)

!********************************************************
!
!     INput residue number and name for reference PDB
!
!********************************************************

      use COMVAL

      implicit none

      ! Input FILe NAMe
      character(*),intent(in):: FILnam

      integer(4):: i,j,ist,ien
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
      allocate(ATMnumR(i),ATMnmR(i),RESnmR(i),RESnumR(i),CHNr(i),      &
      nCHNr(i),refcod(3,i))
      nATMr = i

!     Input reference PDB file
      open(unit=1,file=trim(FILnam),status="old") ; i = 0
      do
        read(1,'(a)',end=800)line
        if ( line(1:6).ne."ATOM  " .and. line(1:6).ne."HETATM" ) cycle
        i = i + 1
        read(line,'(6x,i6,a4,x,a4,a1,i8,3f8.3)')ATMnumR(i),ATMnmR(i),  &
          RESnmR(i),CHNr(i),RESnumR(i),refcod(1:3,i)
      enddo
800   close(1)
      ist = minval(RESnumR) ; ien = maxval(RESnumR)

      c = CHNr(1) ; j = 1 ; nCHNr(1) = j
      do i = 2,nATMr
        if ( CHNr(i) .ne. c ) then
          j = j + 1 ; c = CHNr(i)
        endif
        nCHNr(i) = j
      enddo

!***********************************
!     output for check

      write(6,*)
      write(6,'(4x,a)')"+ N of atoms in reference_PDB = "
      write(6,'(6x,i0)')nATMr
      write(6,'(4x,a)')"+ Residues in reference_PDB : "
      write(6,'(6x,i0,a,i0)')ist," - ",ien
      allocate(aRESr(ist:ien),chk(ist:ien))
      chk(:) = .false. ; aRESr(:) = "    "
      do i = 1,nATMr
        tmp = RESnmR(i) ; j = RESnumR(i)
        if ( chk(j) .and. tmp.ne.RESnmR(i) ) call error(10301)
        chk(j) = .true. ; aRESr(j) = tmp
      enddo

      write(6,'(4x,a)')"+ Sequence :"
      j = 0
      do i = ist,ien
        j = j + 1
        if ( mod(j,10) .eq. 1 ) then
          write(6,'(6x,a4,$)')aRESr(i)
        else
          write(6,'(x,a4,$)')aRESr(i)
        endif
        if ( mod(j,10) .eq. 0 ) write(6,*)
      enddo
      if ( mod(j,10) .ne. 0 ) write(6,*)

!***********************************

      return
      end subroutine inref
