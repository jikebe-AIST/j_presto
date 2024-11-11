
      subroutine setmnt(iimntu,iprint,cimntn,ier)

!**************************************************
!
!  read a file to specify trajectories in MD
!
!**************************************************

      use COMBAS ; use COMPSC ; use COMERG,only: iynbnd,iynang,iyntor

      implicit none

      integer(4):: iimntu,iprint,ier
      character(*):: cimntn

      integer(4):: i,j,k,l,efcol,icol,ich,irs,jresmx,ierr,itemp(4),    &
        ichd(2),irsd(2),ichg(3),irsg(3),icht(4),irst(4),ntotal
      character(80):: line,space
      character(4):: comand,subcom,cat,catd(2),catg(3),catt(4)
      character(1):: yon

      integer(4),parameter:: maxval = 20
      character(1):: inpmod(maxval)
      integer(4):: wrk1(maxval),wrk2(maxval),ival(maxval)
      real(8):: rval(maxval)
      character(80):: cval(maxval)

!**************************************

      ier = 0 ; nacntr = 0 ; ndcntr = 0 ; ngcntr = 0 ; ntcntr = 0
      space = " "
      allocate(ianatm(ixnatm),idnatm(2,iynbnd),ignatm(3,iynang),       &
               itnatm(4,iyntor))
      write(iprint,*)" "
      write(iprint,*)"INFORMATION> SETMNT"
      call flopen(iimntu,cimntn,10,'NULL',0,ier)
      if ( ier .ne. 0 ) then
        write(iprint,*)"ERROR> OPEN FILE ERROR IN SETMNT"
        write(iprint,'(a80)')cimntn ; return
      endif

      OUTER : do
        read(iimntu,'(a80)',err=999,end=800)line
        icol = efcol(line,space,";")
        if ( icol .le. 0 ) cycle
        inpmod(1:2) = "C"
        call rdfree(line,icol,maxval,2,inpmod,wrk1,wrk2,rval,ival,     &
                    cval,ier)
        if ( ier .ne. 0 ) then
          write(iprint,*)"ERROR> READ DATA ERROR IN SETMNT" ; return
        endif
        comand = cval(1)(1:4) ; subcom = cval(2)(1:4)
        if ( comand .ne. "MONI" ) cycle
        if ( subcom .eq. "STOP" ) exit
        if ( subcom .eq. "COOR" ) then
          write(iprint,*)""
          write(iprint,*)" MONITOR ATOM COORDINATES"
          do
            read(iimntu,'(a80)',err=999,end=800)line
            icol = efcol(line,space,";")
            if ( icol .le. 0 ) cycle
            inpmod(1) = "C"
            call rdfree(line,icol,maxval,1,inpmod,wrk1,wrk2,rval,ival, &
                        cval,ier)
            if ( ier .ne. 0 ) then
              write(iprint,*)"ERROR> READ DATA ERROR IN SETMNT" ; return
            endif
            comand = cval(1)(1:4)
            if ( comand .eq. "END " ) cycle OUTER
            inpmod(1:2) = "I" ; inpmod(3:4) = "C"
            call rdfree(line,icol,maxval,4,inpmod,wrk1,wrk2,rval,ival, &
                        cval,ier)
            if ( ier .ne. 0 ) then
              write(iprint,*)"ERROR> READ DATA ERROR IN SETMNT" ; return
            endif
            ich= ival(1) ; irs= ival(2) ; cat= cval(1) ; yon= cval(2)
            if ( ich.gt.0 .and. ich.le.ixnchn ) then
              jresmx = ixares(ixcend(ich))
              if ( irs.gt.0 .and. irs.le.jresmx ) then
                ierr = 1
                do i = 1,ixnatm
                  if ( ich.eq.ixachn(i) .and. irs.eq.ixares(i) .and.   &
                       cat(1:4).eq.cxatmn(i)(1:4) ) then
                    nacntr = nacntr + 1 ; ianatm(nacntr) = i
                    if ( yon(1:1).eq."Y" .or. yon(1:1).eq."y" ) then
                      write(iprint,'(3(i6,x),a4,x,i4,x,a4)')nacntr,i,  &
                        ixachn(i),cxresn(i),ixares(i),cxatmn(i)
                    endif
                    ierr = 0 ; exit
                  endif
                enddo
                if ( ierr .ne. 0 ) then
                  write(iprint,*)"ERROR> SETMNT"
                  write(iprint,*)"      NO TARGET ATOM FOUND."
                  ier = -102 ; return
                endif
              endif
            endif
          enddo

        elseif ( subcom .eq. "DIST" ) then
          write(iprint,*)" "
          write(iprint,*)" MONITOR DISTANCES BETWEEN TWO ATOMS"
          do
            read(iimntu,'(a80)',err=999,end=800)line
            icol = efcol(line,space,";")
            if ( icol .le. 0 ) cycle
            inpmod(1) = "C"
            call rdfree(line,icol,maxval,1,inpmod,wrk1,wrk2,rval,ival, &
                        cval,ier)
            if ( ier .ne. 0 ) then
              write(iprint,*)"ERROR> READ DATA ERROR IN SETMNT" ; return
            endif
            comand = cval(1)(1:4)
            if ( comand .eq. "END " ) cycle OUTER
            inpmod(1:7) = (/"I","I","C","I","I","C","C"/)
            call rdfree(line,icol,maxval,7,inpmod,wrk1,wrk2,rval,ival, &
                        cval,ier)
            if ( ier .ne. 0 ) then
              write(iprint,*)"ERROR> READ DATA ERROR IN SETMNT" ; return
            endif
            ichd(1:2) = (/ival(1),ival(3)/)
            irsd(1:2) = (/ival(2),ival(4)/)
            catd(1:2) = cval(1:2) ; yon = cval(3) ; ierr = 0
            OUTER2 : do i = 1,2
              if ( ichd(i).gt.0 .and. ichd(i).le.ixnchn ) then
                jresmx = ixares(ixcend(ichd(i)))
                if ( irsd(i).gt.0 .and. irsd(i).le.jresmx ) then
                  do j = 1,ixnatm
                    if ( ichd(i).eq.ixachn(j) .and.                    &
                         irsd(i).eq.ixares(j) .and.                    &
                         catd(i)(1:4).eq.cxatmn(j)(1:4) ) then
                      itemp(i) = j ; ierr = ierr + 1 ; cycle OUTER2
                    endif
                  enddo
                endif
              endif
            enddo OUTER2
            ndcntr = ndcntr + 1
            if ( ierr .ne. 2 ) then
              write(iprint,*)"ERROR> SETMNT"
              write(iprint,*)"      NO TARGET ATOM FOUND."
              ier = -112 ; return
            endif
            idnatm(1:2,ndcntr) = itemp(1:2)
            if ( yon(1:1).eq."Y" .or. yon(1:1).eq."y" ) then
              i = idnatm(1,ndcntr) ; j = idnatm(2,ndcntr)
              write(iprint,'(2(2(x,i6),x,a4,x,i4,x,a4))')i,ixachn(i),  &
                cxresn(i),ixares(i),cxatmn(i),j,ixachn(j),cxresn(j),   &
                ixares(j),cxatmn(j)
            endif
          enddo

        elseif ( subcom .eq. "ANGL" ) then
          write(iprint,*)" "
          write(iprint,*)" MONITOR ANGLES FORMED BY THREE ATOMS"
          do
            read(iimntu,'(a80)',err=999,end=800)line
            icol = efcol(line,space,";")
            if ( icol .le. 0 ) cycle
            inpmod(1) = "C"
            call rdfree(line,icol,maxval,1,inpmod,wrk1,wrk2,rval,ival, &
                        cval,ier)
            if ( ier .ne. 0 ) then
              write(iprint,*)"ERROR> READ DATA ERROR IN SETMNT" ; return
            endif
            comand = cval(1)(1:4)
            if ( comand .eq. "END " ) cycle OUTER
            inpmod(1:10) = (/"I","I","C","I","I","C","I","I","C","C"/)
            call rdfree(line,icol,maxval,10,inpmod,wrk1,wrk2,rval,ival,&
                        cval,ier)
            if ( ier .ne. 0 ) then
              write(iprint,*)"ERROR> READ DATA ERROR IN SETMNT" ; return
            endif
            ichg(1:3) = (/ival(1),ival(3),ival(5)/)
            irsg(1:3) = (/ival(2),ival(4),ival(6)/)
            catg(1:3) = cval(1:3) ; yon = cval(4) ; ierr = 0
            OUTER3 : do i = 1,3
              if ( ichg(i).gt.0 .and. ichg(i).le.ixnchn ) then
                jresmx = ixares(ixcend(ichg(i)))
                if ( irsg(i).gt.0 .and. irsg(i).le.jresmx ) then
                  do j = 1,ixnatm
                    if ( ichg(i).eq.ixachn(j) .and.                    &
                         irsg(i).eq.ixares(j) .and.                    &
                         catg(i)(1:4).eq.cxatmn(j)(1:4) ) then
                      itemp(i) = j ; ierr = ierr + 1 ; cycle OUTER3
                    endif
                  enddo
                endif
              endif
            enddo OUTER3
            ngcntr = ngcntr + 1
            if ( ierr .ne. 3 ) then
              write(iprint,*)"ERROR> SETMNT"
              write(iprint,*)"      NO TARGET ATOM FOUND."
              ier = -122 ; return
            endif
            ignatm(1:3,ngcntr) = itemp(1:3)
            if ( yon(1:1).eq."Y" .or. yon(1:1).eq."y" ) then
              i = ignatm(1,ngcntr) ; j = ignatm(2,ngcntr)
              k = ignatm(3,ngcntr)
              write(iprint,'(i6,2(2(x,i6),x,a4,x,i4,x,a4))')ngcntr,    &
                i,ixachn(i),cxresn(i),ixares(i),cxatmn(i),             &
                j,ixachn(j),cxresn(j),ixares(j),cxatmn(j)
              write(iprint,'(6x,2(x,i6),x,a4,x,i4,x,a4)')              &
                k,ixachn(k),cxresn(k),ixares(k),cxatmn(k)
            endif
          enddo

        elseif ( subcom .eq. "TORS" ) then
          write(iprint,*)" "
          write(iprint,*)" MONITOR TORSIONS FORMED BY FOUR ATOMS"
          do
            read(iimntu,'(a80)',err=999,end=800)line
            icol = efcol(line,space,";")
            if ( icol .le. 0 ) cycle
            inpmod(1) = "C"
            call rdfree(line,icol,maxval,1,inpmod,wrk1,wrk2,rval,ival, &
                        cval,ier)
            if ( ier .ne. 0 ) then
              write(iprint,*)"ERROR> READ DATA ERROR IN SETMNT" ; return
            endif
            comand = cval(1)(1:4)
            if ( comand .eq. "END " ) cycle OUTER
            inpmod(1:13) = (/"I","I","C","I","I","C","I","I","C","I",  &
                             "I","C","C"/)
            call rdfree(line,icol,maxval,13,inpmod,wrk1,wrk2,rval,ival,&
                        cval,ier)
            if ( ier .ne. 0 ) then
              write(iprint,*)"ERROR> READ DATA ERROR IN SETMNT" ; return
            endif
            icht(1:4) = (/ival(1),ival(3),ival(5),ival(7)/)
            irst(1:4) = (/ival(2),ival(4),ival(6),ival(8)/)
            catt(1:4) = cval(1:4) ; yon = cval(5) ; ierr = 0
            OUTER4 : do i = 1,4
              if ( icht(i).gt.0 .and. icht(i).le.ixnchn ) then
                jresmx = ixares(ixcend(icht(i)))
                if ( irst(i).gt.0 .and. irst(i).le.jresmx ) then
                  do j = 1,ixnatm
                    if ( icht(i).eq.ixachn(j) .and.                    &
                         irst(i).eq.ixares(j) .and.                    &
                         catt(i)(1:4).eq.cxatmn(j)(1:4) ) then
                      itemp(i) = j ; ierr = ierr + 1 ; cycle OUTER4
                    endif
                  enddo
                endif
              endif
            enddo OUTER4
            ntcntr = ntcntr + 1
            if ( ierr .ne. 4 ) then
              write(iprint,*)"ERROR> SETMNT"
              write(iprint,*)"      NO TARGET ATOM FOUND."
              ier = -132 ; return
            endif
            itnatm(1:4,ntcntr) = itemp(1:4)
            if ( yon(1:1).eq."Y" .or. yon(1:1).eq."y" ) then
              i = itnatm(1,ntcntr) ; j = itnatm(2,ntcntr)
              k = itnatm(3,ntcntr) ; l = itnatm(4,ntcntr)
              write(iprint,'(i6,2(2(x,i6),x,a4,x,i4,x,a4))')ntcntr,    &
                i,ixachn(i),cxresn(i),ixares(i),cxatmn(i),             &
                j,ixachn(j),cxresn(j),ixares(j),cxatmn(j)
              write(iprint,'(6x,2(2(x,i6),x,a4,x,i4,x,a4))')           &
                k,ixachn(k),cxresn(k),ixares(k),cxatmn(k),             &
                l,ixachn(l),cxresn(l),ixares(l),cxatmn(l)
            endif
          enddo

        endif
      enddo OUTER

800   close(iimntu)
      if ( nacntr .gt. 0 ) then
        write(iprint,*)""
        write(iprint,'(i5,a31)')nacntr," COORDINATES WILL BE MONITORED."
      endif
      if ( ndcntr .gt. 0 ) then
        write(iprint,*)""
        write(iprint,'(i5,a31)')ndcntr," DISTANCES   WILL BE MONITORED."
      endif
      if ( ngcntr .gt. 0 ) then
        write(iprint,*)""
        write(iprint,'(i5,a31)')ngcntr," ANGLES      WILL BE MONITORED."
      endif
      if ( ntcntr .gt. 0 ) then
        write(iprint,*)""
        write(iprint,'(i5,a31)')ntcntr," TORSIONS    WILL BE MONITORED."
      endif
      write(iprint,*)""
      ntotal = nacntr + ndcntr + ngcntr + ntcntr
      if ( ntotal .eq. 0 ) then
        write(iprint,'(a)')"  NO VALUES WILL BE MONITORED IN FILE."
      else
        write(iprint,'("  TOTAL ",i5,a)')ntotal,                       &
          " VALUES WILL BE MONITORED IN FILE." 
      endif
      write(iprint,*)""

!************************************

      return

999   write(iprint,*)"ERROR> SETMNT"
      write(iprint,*)"       READ FILE ERROR"
      write(iprint,'(a80)')cimntn ; ier = -1 ; return

      end subroutine setmnt
