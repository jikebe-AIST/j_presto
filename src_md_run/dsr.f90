
      subroutine dsr(infdst,iout,filena,ier)

!*********************************************************************
!
!       THE SUBROUTINE 'DSR' READS DISTANCE CONSTRAINTS DATA,
!       THEN SEARCH ATOM POINTER(ABSOLUTE ATOM NUMBER) FOR DISTANCE
!       CONSTRAINTS ENERGY CALCULATION.
!
!*********************************************************************

      use COMBAS ; use COMMIS 

      implicit none

      ! Logical unit number
        integer(4),intent(in):: infdst,iout
      ! File name
        character(80),intent(in):: filena
      ! Condition code
        integer(4),intent(out):: ier

      integer(4),parameter:: maxval = 20
      integer(4):: inoe,icnt,ichain(2),irsnum(2),icol,efcol,itmp,i,j
      integer(4):: wrk1(maxval),wrk2(maxval),ival(maxval)
      character(1):: inpmod(maxval)
      character(4):: atmnam(2),resnam(2),comand,subcom,output(maxdst)
      character(80):: line,cval(maxval),blank
      real(8):: rval(maxval)

!****************************************************
!     Initialize
      blank = " " ; iuipnt(1:maxdst,1:6) = 0 ; iujpnt(1:maxdst,1:6) = 0
      iuipar(1:maxdst) = 0 ; iujpar(1:maxdst) = 0

      call flopen(infdst,filena,12,'NULL',0,ier)
      if ( ier .lt. 0 ) then
        write(iout,'(a)')"ERROR> DSR"
        write(iout,'(7x,a22,i6)')"OPEN FAILURE IN FLOPEN",infdst
        return
      endif

      ! THE LINE AFTER ';' IS A COMMENT LINE.
      ! IF 'STOP' WAS READ, THEN TERMINATE DISTANCE CONSTRAINTS INPUT.
      inoe = 0
      OUTER : do
        read(infdst,'(a80)',end=800)line
        icol = efcol(line,blank,';')
        if ( icol .le. 0 ) cycle
        inpmod(1:2) = 'C'
        call rdfree(line,icol,maxval,2,inpmod,wrk1,wrk2,rval,ival,cval,&
                    ier)
        if ( ier .lt. 0 ) then
          write(iout,'(a13)')"ERROR> DSR"
          write(iout,'(7x,a20)')"READ ERROR IN RDFREE"
          return
        endif
        comand = cval(1)(1:4) ; subcom = cval(2)(1:4)
        if ( comand .ne. "DSR>" ) cycle
        if ( subcom .eq. "STOP" ) then
          exit
        ! Read distance constraints in free format
        elseif ( subcom .eq. "LIST" ) then
          do
            read(infdst,'(a80)',end=800)line
            icol = efcol(line,blank,';')
            if ( icol .le. 0 ) cycle
            inpmod(1) = 'C'
            call rdfree(line,icol,maxval,1,inpmod,wrk1,wrk2,rval,ival, &
                        cval,ier)
            if ( ier .lt. 0 ) then
              write(iout,'(a13)')"ERROR> DSR"
              write(iout,'(7x,a20)')"READ ERROR IN RDFREE"
              return
            endif
            comand = cval(1)(1:4)
            if ( comand .eq. "END " ) exit OUTER

            inpmod(1:2) = 'I' ; inpmod(3:4) = 'C' ; inpmod(5:6) = 'I'
            inpmod(7:8) = 'C' ; inpmod(9:12) = 'R' ; inpmod(13) = 'C'
            call rdfree(line,icol,maxval,13,inpmod,wrk1,wrk2,rval,ival,&
                        cval,ier)
            if ( ier .lt. 0 ) then
              write(iout,'(a13)')"ERROR> DSR"
              write(iout,'(7x,a20)')"READ ERROR IN RDFREE"
              return
            endif

            inoe = inoe + 1
            if ( inoe .gt. maxdst ) then
              write(iout,'(a)')"ERROR> DSR"
              write(iout,'(7x,a)')"MAXIMUM NUMBER OF INPUT DISTANCES IS"
              write(iout,'(a20,i8)')" EXCEEDED THE LIMIT ",maxdst
              ier = -100 ; return
            endif

            ichain(1) = ival(1) ; ichain(2) = ival(3)
            irsnum(1) = ival(2) ; irsnum(2) = ival(4)
            resnam(1) = cval(1)(1:4) ; resnam(2) = cval(3)(1:4)
            atmnam(1) = cval(2)(1:4) ; atmnam(2) = cval(4)(1:4)
            fudlow(inoe) = rval(1) ; fudupr(inoe) = rval(2)
            furlow(inoe) = rval(3) ; furupr(inoe) = rval(4)
            output(inoe) = cval(5)(1:1)
            if ( rval(1).le.0.d0 .or. rval(2).le.0.d0 ) then
              write(iout,'(a)')"ERROR> DSR"
              write(iout,'(7x,a)')"INPUT DATA ERROR - INPUT"//         &
                                  "POSITIVE VALUES"
              ier = -100 ; return
            endif

            if ( output(inoe)(1:1) .eq. "Y" ) then
              if ( inoe .eq. 1 ) then
                write(iout,*)"INFORMATION> DSR"
                write(iout,'(7x,a)')"READING FOLLOWING "//             &
                                    "DISTANCE RESTRAINTS"
              endif
              write(iout,'(3i5,x,a4,x,a4,x,2i5,x,a4,x,a4,2x,4f6.2)')   &
                  inoe,ichain(1),irsnum(1),resnam(1),atmnam(1),        &
                  ichain(2),irsnum(2),resnam(2),atmnam(2),fudlow(inoe),&
                  fudupr(inoe),furlow(inoe),furupr(inoe)
            endif

            ! Search atom pointer for distance calculations
            do j = 1,2
              icnt = 0 ; itmp = index(atmnam(j),"*")
              if ( itmp .eq. 0 ) then
                itmp = 4
              else
                itmp = itmp - 1
              endif
              do i = 1,ixnatm
                if ( ichain(j).eq.ixachn(i) .and.                      &
                     irsnum(j).eq.ixares(i) .and.                      &
                     resnam(j).eq.cxresn(i) .and.                      &
                     atmnam(j)(1:itmp).eq.cxatmn(i)(1:itmp) ) then
                  icnt = icnt + 1
                  if ( j .eq. 1 ) then
                    iuipnt(inoe,icnt) = i ; iuipar(inoe) = icnt
                  else
                    iujpnt(inoe,icnt) = i ; iujpar(inoe) = icnt
                  endif
                endif
                if ( ixachn(i) .gt. ichain(j) ) exit
              enddo
            enddo

            ! Error finding atom
            if ( iuipar(inoe) .le. 0 ) then
              write(iout,9050)irsnum(1),resnam(1),atmnam(1),inoe
              ier = -100 ; return
            elseif ( iujpar(inoe) .le. 0 ) then
              write(iout,9050)irsnum(2),resnam(2),atmnam(2),inoe
              ier = -100 ; return
            endif
9050        format(/'ERROR> DSR',/,7x,'PROBLEMS FINDING RESIDUE(',  &
                   i4,x,a4' ) OF ATOM( ',a4,' ) : POINTER ',i4)

            ! Dump selected atom information
            if ( output(inoe)(1:1) .eq. "Y" ) then
              write(iout,'(4x,a,6i5,x,a,6i5,a)')                       &
                   "ATOM NUMBER( I-PAIR ;",iuipnt(inoe,1:6),           &
                   "/ J-PAIR : ",iujpnt(inoe,1:6),")"
            endif
          enddo

        endif

      enddo OUTER

800   iutdsc = inoe
      write(iout,*)"INFORMATION> DSR"
      write(iout,'(7x,a,i6)')"TOTAL NUMBER OF DISTANCE RESTRAINTS:",   &
                             iutdsc

      where ( furlow(1:inoe) .lt. 0.d0 ) furlow = 0.d0

      do i = 1,inoe
        if ( furlow(i) .gt. furupr(i) ) then
          write(iout,'(a)')"ERROR> DSR"
          write(iout,'(a,i6,2f8.3)')"       ERRONEOUS DISTANCE "//     &
                "BOUNDS ( LOWER > UPPER ) :",i,furlow(i),furupr(i)
          ier = -100 ; return
        endif
      enddo

      call flclos(infdst,10,ier)
      if ( ier .lt. 0 ) then
        write(iout,*)"INFORMATION> DSR"
        write(iout,'(7x,a,i6)')"TOTAL NUMBER OF DISTANCE RESTRAINTS:", &
                               infdst
        return
      endif

!**********************************************

      return
      end subroutine dsr
