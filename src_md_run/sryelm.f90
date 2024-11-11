
      subroutine sryelm(iprint,numele,elemnt,line,icol,mxelel,nelel,   &
                        posele,iwork,ier)
 
!*******************************************************************
!
!     THIS SUBROUTINE IS FOR SEARCH ELEMENTS
!
!*******************************************************************
 
      implicit none

      ! Logical unit number of output log
        integer(4),intent(in):: iprint
      ! Number of elements
        integer(4),intent(in):: numele
      ! Elements
        character(6),intent(in):: elemnt(numele)
      ! Input line
        character(80),intent(in):: line
      ! Effective column number
        integer(4),intent(in):: icol
      ! Maximum number of elements in one line
        integer(4),intent(in):: mxelel
      ! Number of elements in one line
        integer(4),intent(out):: nelel
      ! Element number in current line
        integer(4),intent(out):: posele(mxelel)
      ! Work array
        integer(4),intent(inout):: iwork(mxelel)
      ! Condition code (0: NO ERROR, -1: NELEL > MXELEL, -2: UNKNOWN KEYWORD)
        integer(4),intent(out):: ier

      integer(4):: numseg,iele,jcol
 
!**********************************************

      ier = 0 ; nelel = 0 
      ! 1) COUNT NUMBER OF SEGMENTS
      call couseg(line,icol,numseg)

      ! 2) SEARCH POINT OF EACH ELEMEMT IN CURRENT LINE
      do iele = 1,numele
        jcol = index(line(1:icol),elemnt(iele))
        if ( jcol .eq. 1 ) then
          nelel = nelel + 1
          if ( nelel .gt. mxelel ) then
            ier = -1 ; exit
          endif
          iwork(nelel) = jcol ; posele(nelel) = iele
        elseif ( jcol .gt. 1 ) then
          if ( line(jcol-1:jcol-1) .eq. ' ' ) then
            nelel = nelel + 1
            if ( nelel .gt. mxelel ) then
              ier = -1 ; exit
            endif
            iwork(nelel) = jcol ; posele(nelel) = iele
          endif
        endif
      enddo

      if ( ier .eq. -1 ) then
        write(iprint,*)'ERROR> SRYELM '
        write(iprint,*)' NUMBER OF KEYWORDS IN ONE LINE IS TOO LARGE'
        write(iprint,*)'        LIMIT IS ',mxelel
        return
      endif
 
      ! 3) SORTING
      if ( numseg .ne. nelel ) then
        write(iprint,*)'ERROR> SRYELM '
        write(iprint,*)' UNKNOWN KEYWORD IS DETECTED '
        write(iprint,*)line
        ier = -2 ; return
      endif

      if ( nelel .gt. 1 ) call sortel(nelel,posele,iwork)

!***************************************

      return
      end subroutine sryelm
 

!=====================================================================


      subroutine couseg(line,icol,numseg)

!*******************************************************************
!
!     THIS SUBROUTINE IS FOR COUNTING SEGMENTS
!
!*******************************************************************
 
      implicit none

      ! Line
        character(*),intent(in):: line
      ! Effective column number
        integer(4),intent(in):: icol
      ! Number of segments
        integer(4),intent(out):: numseg

      logical(1):: ondat = .false.
      integer(4):: i
 
!*********************************************

      numseg = 0
      do i = 1,icol
        if ( ondat ) then
          if ( line(i:i) .eq. ' ' ) then
            numseg = numseg + 1
            ondat = .false.
          endif
        else
          if ( line(i:i) .ne. ' ' ) ondat = .true.
        endif
      enddo

      if ( line(icol:icol) .ne. ' ' ) numseg = numseg + 1
      numseg = numseg / 2

!*************************************

      return
      end subroutine couseg 
 

!===================================================================


      subroutine sortel(nelel,posele,iwork)

!*******************************************************************
!
!     THIS SUBROUTINE IS FOR SORT ELEMENTS BUBBLE SORT
!
!*******************************************************************

      implicit none

      ! Number of elements in one line
        integer(4),intent(in):: nelel
      ! Element number in current line
        integer(4),intent(inout):: posele(nelel)
      ! Work array
        integer(4),intent(inout):: iwork(nelel)

      integer(4):: i,j(1),jloc,itmp
 
!*******************************************

      if ( nelel .eq. 1 ) return

      do i = 1,nelel-1
        j = minloc(iwork(i:nelel)) ; jloc = j(1) + i - 1
        itmp = iwork(i)
        iwork(i) = iwork(jloc)
        iwork(jloc) = itmp
        itmp = posele(i)
        posele(i) = posele(jloc)
        posele(jloc) = itmp
      enddo

!***************************************

      return
      end subroutine sortel
