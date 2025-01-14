
      subroutine rdfree(line,icol,maxval,numval,inpmod,stcol,encol,    &
                        rval,ival,cval,ier)

!******************************************************************
!
!     READ FREE-FORMAT LINE
!
!     9 NOTE
!    9.1 IF NUMRED (NUMBER OF VALUES IN LINE) IS DIFFERENT FROM
!        NUMVAL , MIN( NUMVAL,NUMRED ) IS USED FOR READING DATA
!         EXAMPLE)
!           IF YOU WANT TO READ TOP 5 DATA IN ONE LINE (THIS LINE
!           CONTAINS MORE 5 DATA) , NUMVAL IS 5
!    9.2 THIS IS NOT AVAILABLE TO CHECK DATA TYPE
!
!******************************************************************
 
      implicit none

      ! Input line
        character(*),intent(in):: line
      ! Effective column number
        integer(4),intent(in):: icol
      ! (Max) number of values
        integer(4),intent(in):: maxval,numval
      ! Input mode (R: real, I: integer, C: character)
        character(1),intent(in):: inpmod(maxval)
      ! Start & end column number
        integer(4),intent(inout):: stcol(maxval),encol(maxval)
      ! values
        real(8),intent(out):: rval(maxval)
        integer(4),intent(out):: ival(maxval)
        character(*),intent(out):: cval(maxval)
      ! Condition code (0: NO ERROR, -1: INPUT PARAM ERROR, -2: READ ERROR)
        integer(4),intent(out):: ier

      integer(4):: kval,nrval,nival,ncval,numred,numchk

!*************************************************************
!     1) CHECK INPUT PARAMETER
      if ( maxval.le.0 .or. numval.le.0 .or. maxval.lt.numval .or.     &
           icol.le.0 ) then
        ier = -1 ; return
      endif
      do kval = 1,numval
        if ( index('RIC',inpmod(kval)) .eq. 0 ) then
          ier = -1 ; return
        endif
      enddo

      ier = 0 ; rval(1:maxval) = 0.d0 ; ival(1:maxval) = 0
      cval(1:maxval) = ' '
 
!     2) ASSIGN COLUMNS
      call ascolm(line,icol,maxval,numred,stcol,encol)
      if ( numred .eq. 0 ) return
      numchk = min(numval,numred)
 
!     3) READ VALUES
      nrval = 0 ; nival = 0 ; ncval = 0
      do kval = 1,numchk
!       3-1) READ REAL*8 VALUES
        if ( inpmod(kval) .eq. "R" ) then
          nrval = nrval + 1
          read(line(stcol(kval):encol(kval)),*,iostat=ier)rval(nrval)
          if ( ier .gt. 0 ) then
            ier = -2 ; return
          endif
!       3-2) READ INTEGER*4 VALUES
        elseif ( inpmod(kval) .eq. "I" ) then
          nival = nival + 1
          read(line(stcol(kval):encol(kval)),*,iostat=ier)ival(nival)
          if ( ier .gt. 0 ) then
            ier = -2 ; return
          endif
!       3-3) READ CHARACTER VALUES
        elseif ( inpmod(kval) .eq. "C" ) then
          ncval = ncval + 1
          cval(ncval)(1:(encol(kval)-stcol(kval)+1)) =                 &
                                     line(stcol(kval):encol(kval))
        endif
      enddo

!****************************************************

      return
      end subroutine rdfree
 

!======================================================================


      subroutine ascolm(line,icol,maxval,numval,stcol,encol) 

      implicit none

      character(*):: line
      integer(4):: icol,maxval,numval,stcol(maxval),encol(maxval)
      integer(4):: kcol
      logical(1):: onval
 
!**********************************************************

      numval = 0 ; stcol(1:maxval) = 0 ; encol(1:maxval) = 0
      onval = .false.

!     * ASSIGN EFFECTIVE COLUMNS *
      do kcol = 1,icol
        if ( line(kcol:kcol).ne." " .and. line(kcol:kcol).ne."," ) then
          if ( .not. onval ) then
            onval = .true. ; numval = numval + 1
            stcol(numval) = kcol
          endif
        else
          if ( onval ) then
            onval = .false. ; encol(numval) = kcol - 1
          endif
        endif
      enddo

      if ( numval.ne.0 .and. encol(numval).eq.0 ) encol(numval) = icol

!****************************************

      return
      end subroutine ascolm 
