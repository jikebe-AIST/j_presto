
      subroutine inpinp(iread,iprint,ier,onend,iitplf,iicrdf,iifsor,   &
                        citpln,cicrdn,cishkn,civarn,ciboun,cirefn,     &
                        cipscn,cidscn,cidhcn,cimntn,cieCMn,cirep)

!*******************************************************************
!
!     THIS SUBROUTINE IS FOR READ CONTROL DATA FOR INPUT
!
!*******************************************************************

      use COMCMM ; use COMERG,only: iyeflg,iyfshk

      implicit none

      ! Logical unit number for input & output
        integer(4),intent(in):: iread,iprint
      ! Condition code (0: NO ERROR, negative: ERROR)
        integer(4),intent(out):: ier
      ! .true. detect end of control file
        logical(1),intent(out):: onend
      ! Format type of topology file ( 1: FORMATTED, 2: BINARY)
        integer(4),intent(out):: iitplf
      ! Format type of initial cord. file ( 1: PDB, 2: BINARY)
        integer(4),intent(out):: iicrdf
      ! Flag for set origin (1: NO-transformation, 2: transformation to MASS-center)
        integer(4),intent(out):: iifsor
      ! File name of topology, Initial cord., SHAKE, set variables,
      ! boundary, ref.cord., position, distance, dihedral constraints,
      ! monitoring items, cell data, cell parameters, & external CMM
        character(80),intent(out):: citpln,cicrdn,cishkn,civarn,ciboun,&
          cirefn,cipscn,cidscn,cidhcn,cimntn,cieCMn,cirep

      ! Number of elements in this routine
        integer(4),parameter:: numele = 24
      ! Max number of elements in one line
        integer(4),parameter:: mxelel = 10
      ! Max contents in one element
        integer(4),parameter:: maxcon = 5

      character(6):: elemnt(numele) = (/                               &
        'TOPOLO','TPLFMT','COORDI','CRDFMT','SETSHK','SETVAR','SETBOU',&
        'REFCOO','POSITI','DISTAN','DIHEDR','SETORI','OUTMON','EXTCMM',&
        'REPULS','CELLVL','CELSIZ','CELRES','CLUSTA','CLUSTB','GEPSME',&
        'SCALET','LWEIGH','LTCNTL'/)

      character(1):: dataty(numele) = (/                               &
        'D','C','D','C','D','D','D',  'D','D','D','D','C','D','D',     &
        'D','I','R','D','D','D','C',  'D','R','C'/)

      character(4):: datcon(maxcon,numele)
      data datcon(1:2,2) /'FORM','BINA'/
      data datcon(1:2,4) /'PDB ','BINA'/
      data datcon(1:2,12) /'NO  ','YES '/
      data datcon(1:5,21) /'AND ','OR  ','XOR ','LMD ','LMD2'/
      data datcon(1:2,24) /'NO  ','YES '/

      integer(4):: numspe(numele) = (/                                 &
        0,2,0,2,0,0,0,  0,0,0,0,2,0,0,  0,0,0,0,0,0,5, 0,0,2/)

      character(80):: line,space
      integer(4):: posele(mxelel),iwork(mxelel)
      integer(4):: intele(numele)
      real(8):: reaele(numele)
      character(80):: chaele(numele)
      logical(4):: onlist

      integer(4),parameter:: maxval = mxelel*2
      character(1):: inpmod(maxval)
      integer(4):: iwork1(maxval),iwork2(maxval)
      real(8):: rval(maxval)
      integer(4):: ival(maxval)
      character(80):: cval(maxval)

      integer(4):: efcol,numerr,icol,nelel,numval,numint,numcha,icon,i
      integer(4):: numrea,ierr

!****************************************************************

      ! 0) INITIAL SETTING
      space = " "
      ier = 0 ; numerr = 0 ; onend = .false. ; onlist = .false.
      chaele(1) = " " ; intele(2) = 1
      chaele(3) = " " ; intele(4) = 1
      chaele(5:11) = " " ; intele(12) = 1 ; chaele(13:15) = " "
      intele(16) = 3 ; reaele(17) = 6.d0 ; chaele(18:20) = " "
      intele(21) = 1 ; chaele(22) = "0030" ; reaele(23) = 0.d0
      intele(24) = 1

      ! 1) READ CONTROL DATA FOR INPUT
100   do
        read(iread,'(a80)',iostat=ierr,end=200)line
        icol = efcol(line,space,";")
        if ( icol .le. 0 ) cycle

        if ( index(line(1:icol),"QUIT") .ne. 0 ) then
          exit
        elseif ( index(line(1:icol),"EXE>") .ne. 0 ) then
          backspace(iread) ; exit
        endif

        if ( index(line(1:icol),"LIST") .ne. 0 ) then
          onlist = .true. ; exit
        endif

        ! 1-1 ) SEARCH ELEMENTS
        call sryelm(iprint,numele,elemnt,line,icol,mxelel,nelel,       &
                    posele,iwork,ier)
        if ( ier .ne. 0 ) then
          numerr = numerr + 1
          if ( ier .eq. -2 ) then
            write(iprint,'(10(a1,a6))')" ",elemnt(1:numele)
          endif
          if ( numerr .ge. 4 ) then
            return
          else
            ier = 0 ; cycle
          endif
        endif

        if ( nelel .le. 0 ) cycle

        ! 1-2 ) READ DATA OF EACH ELEMENT
        do i = 1,nelel
          inpmod(2*i-1) = "C"
          if ( dataty(posele(i)) .eq. "D" ) then
            inpmod(2*i) = "C"
          else
            inpmod(2*i) = dataty(posele(i))
          endif
        enddo
        numval = nelel*2
        call rdfree(line,icol,maxval,numval,inpmod,iwork1,iwork2,rval, &
                    ival,cval,ier)
        if ( ier .ne. 0 ) then
          write(iprint,*)'ERROR> INPINP '
          write(iprint,*)'       DATA TYPE ERROR '
          write(iprint,*)' '
          write(iprint,'(a80)')line
          write(iprint,*)' '
          return
        endif

        ! 1-3) STORE DATA OF EACH ELEMENT
        numint = 0 ; numrea = 0 ; numcha = 0
        do i = 1,nelel
          if ( dataty(posele(i)) .eq. "I" ) then
            numint = numint + 1 ; numcha = numcha + 1
            intele(posele(i)) = ival(numint)
          elseif ( dataty(posele(i)) .eq. "R" ) then
            numrea = numrea + 1 ; numcha = numcha + 1
            reaele(posele(i)) = rval(numrea)
          elseif ( dataty(posele(i)) .eq. "C" ) then
            numcha = numcha + 2
            do icon = 1,numspe(posele(i))
              if ( cval(numcha)(1:4) .eq. datcon(icon,posele(i)))      &
                intele(posele(i)) = icon
            enddo
          elseif ( dataty(posele(i)) .eq. "D" ) then
            numcha = numcha + 2
            chaele(posele(i)) = cval(numcha)
          endif
        enddo

      enddo

200   if ( ierr .ne. 0 ) onend = .true.

       citpln = chaele(1) ; iitplf = intele(2) ; cicrdn = chaele(3)
       iicrdf = intele(4) ; cishkn = chaele(5) ; civarn = chaele(6)
       ciboun = chaele(7) ; cirefn = chaele(8) ; cipscn = chaele(9)
       cidscn = chaele(10) ; cidhcn = chaele(11) ; iifsor = intele(12)
       cimntn = chaele(13) ; cieCMn = chaele(14) ; cirep = chaele(15)
       nlev = intele(16) ; cminsiz = reaele(17) ; celres = chaele(18)
       clusta = chaele(19) ; clustb = chaele(20)
       cluster_method = datcon(intele(21),21)
       scaleterm = chaele(22) ; lambda_m = reaele(23)
       if ( intele(24) .eq. 1 ) then
         l_temp_control = .false.
       else
         l_temp_control = .true.
       endif

      ! 3) OUTPUT CONTROL DATA OF INPUT
      if ( citpln .ne. " " ) then
        write(iprint,*)' '
        write(iprint,*)'         TOPOLOGY FILE '
        write(iprint,*)'           FORMAT        : ',datcon(iitplf,2)
        write(iprint,*)'           NAME          : ',trim(citpln)
      endif
      if ( cicrdn .ne. " " ) then
        write(iprint,*)' '
        write(iprint,*)'         INITIAL COORDINATE FILE '
        write(iprint,*)'           FORMAT        : ',datcon(iicrdf,4)
        write(iprint,*)'           NAME          : ',trim(cicrdn)
      endif
      if ( cishkn .ne. " " ) then
        write(iprint,*)' '
        write(iprint,*)'         CONTROL DATA FOR SETTING SHAKE '
        write(iprint,*)'           NAME          : ',trim(cishkn)
        iyfshk = 1
      endif
      if ( civarn .ne. " " ) then
        write(iprint,*)' '
        write(iprint,*)'         CONTROL DATA FOR SETTING VARIABLES '
        write(iprint,*)'           NAME          : ',trim(civarn)
      endif
      if ( ciboun .ne. " " ) then
        write(iprint,*)' '
        write(iprint,*)'         CONTROL DATA FOR SETTING BOUNDARY '
        write(iprint,*)'           NAME          : ',trim(ciboun)
        iyeflg(14) = 1
      endif
      if ( cirefn .ne. " " ) then
        write(iprint,*)' '
        write(iprint,*)'         REFERENCE COORDINATE FILE '
        write(iprint,*)'           NAME          : ',trim(cirefn)
      endif
      if ( cipscn .ne. " " ) then
        write(iprint,*)' '
        write(iprint,*)'         CONTROL DATA FOR POSITION CONS. '
        write(iprint,*)'           NAME          : ',trim(cipscn)
        iyeflg(11) = 1
      endif
      if ( cidscn .ne. " " ) then
        write(iprint,*)' '
        write(iprint,*)'         CONTROL DATA FOR DISTANCE CONS. '
        write(iprint,*)'           NAME          : ',trim(cidscn)
        iyeflg(12) = 1
      endif
      if ( cidhcn .ne. " " ) then
        write(iprint,*)' '
        write(iprint,*)'         CONTROL DATA FOR DIHEDRAL CONS. '
        write(iprint,*)'           NAME          : ',trim(cidhcn)
        iyeflg(13) = 1
      endif
      if ( iifsor .eq. 2 ) then
        write(iprint,*)' '
        write(iprint,*)'         SET ORIGIN      : ',datcon(iifsor,12)
      endif
      if ( cimntn .ne. " " ) then
        write(iprint,*)' '
        write(iprint,*)'         MONITORING TRAJECTORY SELECTION '
        write(iprint,*)'           NAME          : ',trim(cimntn)
      endif
      if ( cieCMn .ne. " " ) then
        write(iprint,*)' '
        write(iprint,*)"         EXTERNAL CMM CALCULATION"
        write(iprint,*)'           NAME          : ',trim(cieCMn)
      endif
      if ( cirep .ne. " " ) then
        write(iprint,*)' '
        write(iprint,*)"         REPULSION CALCULATION"
        write(iprint,*)'           NAME          : ',trim(cirep)
        iyeflg(15) = 1
      endif

      write(iprint,*) ' '
      if ( onlist ) then
        onlist = .FALSE.
        goto 100
      endif

!***********************************************

      return
      end subroutine inpinp
