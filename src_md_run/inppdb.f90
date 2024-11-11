
      subroutine inppdb(iread,iprint,filena,maxatm,numatm,atty,        &
     &                  resty,resnum,cord,occupy,bfact,ier)
 
!***********************************************************************
!
!     Read PDB-formatted coordinate data, atom-type, residue-type
!       residue-number coordinate, occupancy, & B-factor
!
!***********************************************************************

      implicit none

      ! Logical unit number for Input & Output
        integer(4),intent(in):: iread,iprint
      ! File name of coordinate
        character(*),intent(in):: filena
      ! MAX number of atoms
        integer(4),intent(in):: maxatm
      ! Number of atoms
        integer(4),intent(out):: numatm
      ! Atom & Residue type of each atom
        character(*),intent(out):: atty(maxatm),resty(maxatm)
      ! Residue number of each atom
        integer(4),intent(out):: resnum(maxatm)
      ! Coordinate, occupancy, & B-factor
        real(8),intent(out)::cord(3,maxatm),occupy(maxatm),bfact(maxatm)
      ! Condition code (0:NO ERROR, -1:FILE OPEN CLOSE ERROR
      !                 -2:DATA READ ERROR, -3: NUMATM .gt. MAXATM)
        integer(4),intent(out):: ier
        
      character(80):: line
      logical(1):: CHNflg = .true.
 
!********************************************

      ! 1) initial setting
      ier = 0 ; numatm = 0
      call flopen(iread,filena,10,'NULL',0,ier)
      if ( ier .ne. 0 ) then
        ier = -1
        write(iprint,*) ' '
        write(iprint,*) 'ERROR> INPPDB '
        write(iprint,*) '    FILE OPEN ERROR'
        write(iprint,*) ' '
        return
      endif
 
      ! 2) read PDB-formatted coordinate file
      do 
        read(iread,'(a80)',end=800,err=910)line
        if ( line(1:4).eq."ATOM" .or. line(1:4).eq."HETA" ) then
          numatm = numatm + 1
          if ( numatm .gt. maxatm ) then
            ier = -3
            write(iprint,*) ' '
            write(iprint,*) 'ERROR> INPPDB '
            write(iprint,*) '    NUMBER OF ATOMS IS GREATER THAN',maxatm
            write(iprint,*) '    NUMBER OF ATOMS IS ',numatm
            write(iprint,*) ' '
            return
          endif
          ! Residue type
          resty(numatm) = line(18:21)
          ! Atom type
          if ( line(13:13) .eq. " " ) then
            atty(numatm) = line(14:16)//" "
          else
            atty(numatm) = line(13:16)
          endif
          ! Residue number
          if ( numatm .eq. 1 ) then
            if ( line(22:22) .eq. " " ) CHNflg = .false.
          endif
          if ( CHNflg ) then
            read(line(23:30),*)resnum(numatm)
          else
            read(line(22:30),*)resnum(numatm)
          endif
          read(line,'(30x,3f8.3,2f6.3)')cord(1:3,numatm),              &
                     occupy(numatm),bfact(numatm)
        endif
      enddo
            
800   if ( numatm .eq. 0 ) then
        ier = -3
        write(iprint,*) ' '
        write(iprint,*) 'ERROR> INPPDB '
        write(iprint,*) '    NUMBER OF ATOMS IS ZERO'
        write(iprint,*) ' '
        return
      endif

      call flclos(iread,10,ier)
      if ( ier .ne. 0 ) then
        ier = -1
        write(iprint,*) ' '
        write(iprint,*) 'ERROR> INPPDB '
        write(iprint,*) '    FILE CLOSE ERROR '
        write(iprint,*) ' '
        return
      endif

      return

!*************

910   ier = -2
      write(iprint,*) ' '
      write(iprint,*) 'ERROR> INPPDB '
      write(iprint,*) '    DATA READ ERROR '
      write(iprint,*) ' '

      return
      end subroutine inppdb
