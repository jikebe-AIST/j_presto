
      subroutine pdb_chk(REFPDB,PDBLST)

!****************************************
!
!     PDB consistency CHecK
!
!****************************************

      use COMVAL,only: d_proj2,nATMr

      implicit none

      ! Reference PDB file name
        character(*),intent(in):: REFPDB
      ! PDB file name LiST
        character(*),intent(in):: PDBLST

      integer(4):: ic
      character(30):: tmp,list(nATMr)
      character(999):: filnam

!******************************

      write(6,*)""
      write(6,'(2x,a)')"+ PDB consistency check start ("//trim(PDBLST) &
                       //")"

      open(unit=1,file=trim(REFPDB),status="old") ; ic = 0
      do
        read(1,'(a)',end=800)tmp
        if ( tmp(1:6).ne."ATOM  " .and. tmp(1:6).ne."HETATM" ) cycle
        ic = ic + 1 ; list(ic) = tmp
      enddo
800   close(1)

      open(unit=1,file=trim(PDBLST),status="old")
      do
        read(1,'(a)',end=801)filnam
        write(6,'(4x,a)')"Now checking "//trim(filnam)
        open(unit=2,file=trim(filnam),status="old") ; ic = 0
        do
          read(2,'(a)',end=802)tmp
          if ( tmp(1:6).ne."ATOM  " .and. tmp(1:6).ne."HETATM" ) cycle
          ic = ic + 1
          if ( ic .gt. nATMr ) call error(10401)
          if ( tmp .ne. list(ic) ) call error(10401)
        enddo
802     close(2)
        if ( ic .ne. nATMr ) call error(10401)
      enddo
801   close(1)

      write(6,'(2x,a)')"+ PDB coidensity check normally done"
      write(6,*)""

!*************************************************

      return
      end subroutine pdb_chk
