
      subroutine inpout(iread,iprint,ier,onend,iotplf,iocrdf,cotpln,   &
                        cocrdn)

!*******************************************************************
!
!     THIS SUBROUTINE IS FOR READ CONTROL DATA FOR IOUTPUT
!
!*******************************************************************

      implicit none

      ! Logical unit number for input & output
        integer(4),intent(in):: iread,iprint
      ! Condition code
        integer(4),intent(out):: ier
      ! Flag for end of control file (end = .true.)
        logical(1),intent(out):: onend
      ! Format type of topology file
      !   1: nowrite, 2: formatted, 3: binary
        integer(4),intent(out):: iotplf
      ! Format type of coordinate file
      !   1: nowrite, 2: PDB, 3: binary
        integer(4),intent(out):: iocrdf
      ! File name of topology & initial coordinate file
        character(80),intent(out):: cotpln,cocrdn

      ! Number of elements in this routine
        integer(4),parameter:: numele = 4
      ! Max number of elements in one line
        integer(4),parameter:: mxelel = 10
      ! Max contents in one element
        integer(4),parameter:: maxcon = 2

      character(6):: elemnt(numele) = (/                               &
               'OTOPOL','OTPLFM','OCOORD','OCRDFM'/)
      character(1):: dataty(numele) = (/                               &
               'D','C','D','C'/)
      character(4):: datcon(maxcon,numele)
      data datcon(1:2,2) /'FORM','BINA'/
      data datcon(1:2,4) /'PDB ','BINA'/
      integer(4):: numspe(numele) = (/0,2,0,2/)

      character(80):: line,space,chaele(numele)
      integer(4):: posele(mxelel),iwork(mxelel),intele(numele)
      logical(1):: onlist

      integer(4),parameter:: maxval = 20
      character(1):: inpmod(maxval),ctmp
      integer(4):: iwork1(maxval),iwork2(maxval),ival(maxval),itmp
      real(8):: rval(maxval)
      character(80):: cval(maxval)

      integer(4):: efcol,numerr,icol,nelel,iele,icon,numval,numint,    &
                   numcha,ierr

!***************************************************************

      ! 0) INITIAL SETTING
      ier = 0 ; numerr = 0 ; onend = .false. ; onlist = .false.
      chaele(1) = " " ; intele(2) = 1
      chaele(3) = " " ; intele(4) = 1
      space = ' '

      ! 1) READ CONTROL DATA FOR INPUT
100   do
        read(iread,'(a80)',iostat=ierr,end=200)line
        icol = efcol(line,space,';')
        if ( icol .le. 0 ) cycle

        if ( index(line(1:icol),'QUIT') .ne. 0 ) then
          exit
        elseif ( index(line(1:icol),'EXE>') .ne. 0 ) then
            backspace(iread) ; exit
        elseif ( index(line(1:icol),'LIST') .ne. 0 ) then
          onlist = .true. ; exit
        endif

        ! 1-1) SERACH ELEMENTS
        call sryelm(iprint,numele,elemnt,line,icol,mxelel,nelel,posele,&
                    iwork,ier)
        if ( ier .ne. 0 ) then
          numerr = numerr + 1
          if ( ier .eq. -2 ) write(iprint,'(10(x,a6))')elemnt(1:numele)
          if ( numerr .ge. 4 ) then
            return
          else
            ier = 0 ; cycle
          endif
        endif
        if ( nelel .le. 0 ) cycle

        ! 1-2) READ DATA OF EACH ELEMENT
        do iele = 1,nelel
          inpmod(2*iele-1) = 'C'
          if ( dataty(posele(iele)) .eq. 'D' ) then
            inpmod(2*iele) = 'C'
          else
            inpmod(2*iele) = dataty(posele(iele))
          endif
        enddo
        numval = nelel * 2
        call rdfree(line,icol,maxval,numval,inpmod,iwork1,iwork2,rval, &
                    ival,cval,ier)
        if ( ier .ne. 0 ) then
          write(iprint,*)'ERROR> INPOUT '
          write(iprint,*)'       DATA TYPE ERROR '
          write(iprint,'(a80)')line
          write(iprint,*)' '
          return
        endif

        ! 1-3) STORE DATA OF EACH ELEMENT
        numint = 0 ; numcha = 0
        do iele = 1,nelel
          ctmp = dataty(posele(iele))
          if ( ctmp .eq. 'I' ) then
            numint = numint + 1 ; numcha = numcha + 1
            intele(posele(iele)) = ival(numint)
          elseif ( ctmp.eq.'C' .or. ctmp.eq.'D' ) then
            numcha = numcha + 2
            if ( ctmp .eq. 'D' ) then
              chaele(posele(iele)) = cval(numcha)
            else
              itmp = posele(iele)
              do icon = 1,numspe(itmp)
                if ( cval(numcha)(1:4) .eq. datcon(icon,itmp) )        &
                  intele(itmp) = icon
              enddo
            endif
          endif
        enddo

      enddo

200   if ( ierr .ne. 0 ) onend = .true.

      ! 2) STORE ELEMENT DATA TO VARIABLE OF INPUT
      cotpln = chaele(1) ; iotplf = intele(2)
      cocrdn = chaele(3) ; iocrdf = intele(4)

      ! 3) OUTPUT CONTROL DATA OF OUTPUT
      if ( cotpln .ne. " " ) then
        write(iprint,*)'      1) TOPOLOGY FILE '
        write(iprint,*)'           FORMAT        : ',datcon(iotplf,2)
        if ( cotpln .ne. space ) then
          write(iprint,*)'           NAME          : '//trim(cotpln)
        endif
      endif
      if ( cocrdn .ne. " " ) then
        write(iprint,*)'      2) COORDINATE FILE '
        write(iprint,*)'           FORMAT        : ',datcon(iocrdf,4)
        if ( cocrdn .ne. space ) then
          write(iprint,*)'           NAME          : '//trim(cocrdn)
        endif
      endif
      write(iprint,*)' '
      if ( onlist ) then
        onlist = .false. ; goto 100
      endif

!**********************************************

      return
      end subroutine inpout
