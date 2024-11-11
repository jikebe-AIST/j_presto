
      subroutine setshk(iread,iprint,filena,ier)
 
!*******************************************************************
!
!     INPUT DATA FILE OF SHAKE AND PREPAR DATA OF SHAKE
!       1) READ CONTROL DATA FOR SHAKE     (INPSHK)
!       2) RE-ARRANGE SHAKE PARAMETER LIST (MODSHK)
!       3) RE-ARRANGE BOND PARAMETER       (REABND)
!       4) RE-ARRANGE ANGLE PARAMETER      (REAANG)
!
!*******************************************************************
 
      use COMBAS ; use COMERG ; use COMMIS

      implicit none

      integer(4),intent(in):: iread,iprint
      character(80),intent(in):: filena
      ! ( 0 : NO ERROR, -1: file open error,
      !  -2 : continuous lines must be less equal 10
      !  -3 : unknown molecule name
      !  -4 : number of atom groups of SHAKE cons. is exceeded
      !  -5 : data type error
      !  -6 : number of atom groups of SHAKE cons. is zero
        integer(4),intent(inout):: ier
 
      integer(4):: istshk,isgshk
 
!******************************************
 
      ier = 0
      if ( iread .gt. 0 ) then
        call inpshk(iread,filena,iprint,ier)
        if ( ier .ne. 0 ) return

        istshk = iutshk ; isgshk = iugshk
        call modshk(iprint,ier)
        if ( ier .ne. 0 ) return

        call reabnd
        call reaang
 
        write(iprint,*)'INFORMATION> INPUT '
        write(iprint,*)'     SHAKE DATA SUMMARY '
        write(iprint,*)'     1) INPUT SHAKE DATA '
        write(iprint,*)'        NUMBER OF SHAKE CONSTRAINTS : ',istshk
        write(iprint,*)'        NUMBER OF SHAKE GROUPS      : ',isgshk
        write(iprint,*)'     2) APPLIED SHAKE DATA '
        write(iprint,*)'        NUMBER OF SHAKE CONSTRAINTS : ',iutshk
        write(iprint,*)'        NUMBER OF SHAKE GROUPS      : ',iugshk
        write(iprint,*)'                 TWO ATOMS          : ',iugsk2
        write(iprint,*)'                 THREE ATOMS        : ',      &
                                            (iugsk3-iugsk2)
        write(iprint,*)'                 FOUR ATOMS         : ',      &
                                            (iugshk-iugsk3)
        write(iprint,*)'     3) BOND PARAMETER '
        write(iprint,*)'        NUMBER OF BONDS             : ',iynbnd
        write(iprint,*)'        NUMBER OF DUMMY BONDS       : ',iyndbn
        write(iprint,*)'     4) ANGLE PARAMETER '
        write(iprint,*)'        NUMBER OF ANGLES            : ',iynang
        write(iprint,*)'        NUMBER OF DUMMY ANGLES      : ',iyndag
        write(iprint,*)' '
      endif
 
      return
      end subroutine setshk


!======================================================================


      subroutine inpshk(iread,filena,iprint,ier)

!*******************************************************************
!
!     INPUT SHAKE CONTROL DATA
!
!*******************************************************************

      use COMBAS ; use COMERG ; use COMMIS 

      implicit none

      integer(4),intent(in):: iread,iprint
      character(80):: filena
      integer(4),intent(inout):: ier

      ! Temporary buffer for store SHAKE control data
        character(800):: buffer
      ! Flag for read molecule name (.true. -> read molecule name)
        logical(1):: onmol
      ! Current molecule number
        integer(4):: molnum
      ! Atom number conversion relative atom number -> absolute one
        integer(4):: iplusa

      integer(4):: efcol,combin,ienbuf,iersub,numshk,imol,icol,icon,igrp
      integer(4):: molgrp(ixmolc)
      character(80):: line,space
 
!*****************************************
 
      ! Initialization
      ier = 0 ; space  = ' ' ; onmol  = .false.
      ienbuf = 0 ; molnum = 0 ; iugshk = 0
      molgrp(1:ixmolc) = 0
 
      call flopen(iread,filena,12,'NULL',0,iersub)
      if ( iersub .ne. 0 ) then
        write(iprint,*)'ERROR> INPSHK '
        write(iprint,*)'       FILE OPEN ERROR  '
        write(iprint,*)' '
        ier = -1 ; return
      endif
 
      ! Reading SHAKE control data
      buffer = ' ' ; ienbuf = 0
      do 
        !! 1) read control data to buffer
        read(iread,'(a80)',end=800)line
        icol = efcol(line,space,';')
        if ( icol .le. 0 ) cycle
        icon = index(line(1:icol),'->') - 1

        if ( icon .ge. 1 ) then
          buffer((ienbuf+1):(ienbuf+icon)) = line(1:icon)
          ienbuf = ienbuf + icon
          cycle
        else
          if ( icon .eq. 0 ) cycle
          buffer((ienbuf+1):(ienbuf+icol)) = line(1:icol)
          ienbuf = ienbuf + icol
        endif
        if ( ienbuf .gt. 800 ) then
          write(iprint,*)'ERROR> INPSHK '
          write(iprint,*)'  CONTINUOUS LINES MUST BE LESS EQUAL 10 '
          write(iprint,*)'  LINES '
          ier = -2 ; return
        endif

        !! 2) read control data from buffer
        if ( index(buffer(1:ienbuf),'SHAKE>') .ne. 0 ) then
          onmol = .true.
          if ( molnum .gt. 0 ) then
            call mkmshk(ixmolc,molgrp,ixsqml,molnum,ixnchn,ixcend,     &
                        maxshk,mxashk,maxequ,iugshk,iuhshk,iuashk,     &
                        fudshk,iprint,ier)
          endif
        else
          if ( onmol ) then
            call srmoln(buffer,ienbuf,ixnchn,ixmolc,cxmolc,ixcend,     &
                        ixsqml,iprint,molnum,iplusa,ier)
            if ( ier .ne. 0 ) return
            onmol = .false.
          else
            call redshk(buffer,ienbuf,iprint,maxshk,mxashk,maxequ,     &
                        iplusa,iugshk,iuhshk,iuashk,fudshk,ier)
            if ( ier .ne. 0 ) return
          endif
        endif

        ienbuf = 0 ; buffer = ' ' ; cycle
      enddo

      !! 3) output information
800   if ( molnum .gt. 0 ) then
        call mkmshk(ixmolc,molgrp,ixsqml,molnum,ixnchn,ixcend,maxshk,  &
                    mxashk,maxequ,iugshk,iuhshk,iuashk,fudshk,iprint,  &
                    ier)
        if ( ier .ne. 0 ) return
      endif
      iutshk = 0
      do igrp = 1,iugshk
        numshk = combin(iuhshk(igrp)+1,2)
        iutshk = iutshk + numshk
      enddo
 
      write(iprint,*)' '
      write(iprint,*)'INFORMATION> INPSHK '
      write(iprint,*)'     NUMBER OF TOTAL SHAKE CONS. : ',iutshk
      write(iprint,*)'     NUMBER OF TOTAL ATOM GROUPS : ',iugshk
      write(iprint,*)'     SHAKE ATOM GROUPS IN EACH MOLECULE '

      do imol = 1,ixmolc
        write(IPRINT,*)'       MOLECULE  ',imol,' : ',molgrp(imol)
      enddo
      write(iprint,*)" "
        
      call flclos(iread,10,ier)
      if ( ier .ne. 0 ) then
        write(iprint,*) 'ERROR> INPSHK '
        write(iprint,*) '       FILE CLOSE ERROR '
        write(iprint,*) ' '
        ier = -1 ; return
      endif

!*****************

      return
      end subroutine inpshk


!=================================================================== 


      subroutine mkmshk(ixmolc,molgrp,ixsqml,molnum,ixnchn,ixcend,     &
                        maxshk,mxashk,maxequ,iugshk,iuhshk,iuashk,     &
                        fudshk,iprint,ier)

!*******************************************************************
!
!     MAKE SHAKE DATA OF ANOTHER CHAIN
!
!*******************************************************************
 
      implicit none

      integer(4),intent(in):: ixmolc
      integer(4),intent(inout):: molgrp(ixmolc)
      integer(4),intent(in):: ixsqml(ixmolc)
      integer(4),intent(in):: molnum,ixnchn,maxshk,mxashk,maxequ
      integer(4),intent(inout):: iugshk
      integer(4),intent(in):: ixcend(ixnchn)
      integer(4),intent(inout):: iuhshk(maxshk),iuashk(mxashk,maxshk)
      real(8),intent(inout):: fudshk(maxequ,maxshk)
      integer(4),intent(in):: iprint
      integer(4),intent(inout):: ier

      integer(4):: combin  ! (function COMBIN)
      integer(4):: lastcn,iminus,numchn,numgrp,iengrp,istgrp,numcon
      integer(4):: imol,ichn,igrp,iplus,itmp1,itmp2
 
!**************************************************

      ! Calculation of number of SHAKE cons. of each molecule
      ier = 0 ; lastcn = 0 ; molgrp(molnum) = iugshk
      if ( molnum .eq. 1 ) then
        iminus = 0
      else
        do imol = 1,molnum-1
          molgrp(molnum) = molgrp(molnum) - molgrp(imol)*ixsqml(imol)
          lastcn = lastcn + ixsqml(imol)
        enddo
        iminus = ixcend(lastcn)
      endif

      ! Make SHAKE constraints of other chains of this molecule 
      numchn = ixsqml(molnum) ; numgrp = molgrp(molnum)
      if ( numchn.ge.2 .and. numgrp.ge.1 ) then
        istgrp = iugshk - numgrp + 1
        iengrp = iugshk
        do ichn = 1,numchn-1
          lastcn = lastcn + 1
          iplus = ixcend(lastcn)
          do igrp = istgrp,iengrp
            iugshk = iugshk + 1
            if ( iugshk .gt. maxshk ) then
              write(iprint,*)'ERROR> INPSHK '
              write(iprint,*)'   NUMBER OF ATOM GROUPS FOR SHAKE'
              write(iprint,*)'   IS EXCEEDED BY LIMIT NUMBER '
              write(iprint,*)' '
              write(iprint,*)'   LIMIT NUMBER IS ',maxshk
              write(iprint,*)' '
              ier = -4 ; return
            endif
            iuhshk(iugshk) = iuhshk(igrp)
            itmp1 = iplus - iminus ; itmp2 = iuhshk(igrp) + 1
            iuashk(1:itmp2,iugshk) = iuashk(1:itmp2,igrp) + itmp1
            numcon = combin(itmp2,2)
            fudshk(1:numcon,iugshk) = fudshk(1:numcon,igrp)
          enddo
        enddo
      endif

!*************************

      return
      end subroutine mkmshk


!=========================================================================

 
      subroutine srmoln(buffer,ienbuf,ixnchn,ixmolc,cxmolc,ixcend,     &
                        ixsqml,iprint,molnum,iplusa,ier)

      implicit none

      character(*),intent(in):: buffer
      integer(4),intent(in):: ienbuf,ixnchn,ixmolc
      character(40),intent(in):: cxmolc(ixmolc)
      integer(4),intent(in):: ixcend(ixnchn),ixsqml(ixmolc)
      integer(4),intent(in):: iprint
      integer(4),intent(inout):: molnum,iplusa,ier

      integer(4):: tpcol,lacol  ! (function)
      integer(4):: itpcol,ilacol,lastcn,imol
      logical(1):: eflg

!********************************************

      eflg = .true.
      do imol = 1,ixmolc
        itpcol = tpcol(cxmolc(imol),40)
        ilacol = lacol(cxmolc(imol),40)
        if ( buffer(1:ienbuf) .eq. cxmolc(imol)(itpcol:ilacol) ) then
          molnum = imol ; eflg = .false. ; exit
        endif
      enddo

      if ( eflg ) then
        molnum = -1 ; ier = -3
        write(iprint,*)'ERROR> INPSHK '
        write(iprint,*)'   UNKNOWN MOLECULE NAME '
        write(iprint,*)buffer(1:ienbuf)
        write(iprint,*)' '
        return
      endif

      if ( molnum .eq. 1 ) then
        iplusa = 0
      else
        lastcn = sum(ixsqml(1:molnum-1))
        iplusa = ixcend(lastcn)
      endif

!*******************

      return
      end subroutine srmoln


!====================================================================== 
 

      subroutine redshk(buffer,ienbuf,iprint,maxshk,mxashk,maxequ,     &
                        iplusa,iugshk,iuhshk,iuashk,fudshk,ier)

      implicit none

      character(*),intent(in):: buffer
      integer(4),intent(in):: ienbuf,iprint,iplusa,maxshk,mxashk,maxequ
      integer(4),intent(inout):: iugshk
      integer(4),intent(inout):: iuhshk(maxshk),iuashk(mxashk,maxshk)
      real(8),intent(inout):: fudshk(maxequ,maxshk)
      integer(4),intent(inout):: ier

      integer(4):: combin ! (function)
      integer(4):: numatm,numcon

      ! Variables for RDFREE
      integer(4):: iwrk1(100),iwrk2(100),ival(100)
      character(1):: inpmod(100)
      real(8):: rval(100)
      character(80):: cval(100)

!*****************************************

      ier = 0
      iugshk = iugshk + 1
      if ( iugshk .gt. maxshk ) then
        write(iprint,*)'ERROR> INPSHK '
        write(iprint,*)'   NUMBER OF ATOM GROUPS FOR SHAKE'
        write(iprint,*)'   IS EXCEEDED BY LIMIT NUMBER '
        write(iprint,*)' '
        write(iprint,*)'   LIMIT NUMBER IS ',maxshk
        write(iprint,*)' '
        ier = -4 ; return
      endif

      inpmod(1) = 'I'
      call rdfree(buffer,ienbuf,100,1,inpmod,iwrk1,iwrk2,rval,ival,    &
                  cval,ier)
      if ( ier .lt. 0 ) then
        write(iprint,*)'ERROR> INPSHK '
        write(iprint,*)'  DATA TYPE ERROR '
        write(iprint,*)buffer(1:ienbuf)
        write(iprint,*)' '
        ier = -5 ; return
      endif
      numatm = ival(1)

      iuhshk(iugshk) = numatm - 1
      numcon = combin(numatm,2)

      inpmod(1:numatm+1) = 'I'
      inpmod(numatm+2:numatm+numcon+1) = 'R'

      call rdfree(buffer,ienbuf,100,numatm+1+numcon,inpmod,iwrk1,iwrk2,&
                  rval,ival,cval,ier)
      if ( ier .lt. 0 ) then
        write(iprint,*)'ERROR> INPSHK '
        write(iprint,*)'  DATA TYPE ERROR '
        write(iprint,*)buffer(1:ienbuf)
        write(iprint,*)' '
        ier = -5 ; return
      endif

      iuashk(1:numatm,iugshk) = ival(2:numatm+1) + iplusa
      fudshk(1:numcon,iugshk) = rval(1:numcon)*rval(1:numcon)

!**************************

      return
      end subroutine redshk 
 

!=============================================================================


      subroutine modshk(iprint,ier)

!****************************************************************
!
!     MODIFY SHAKE CONSTRAINTS DATA
!     REMOVE CONSTRAINTS IF THERE IS FIXED ATOM
!     CHANGE ORDER OF CONSTRAINTS
!
!*******************************************************************
 
      use COMBAS ; use COMERG ; use COMMIS 

      implicit none

      integer(4),intent(in):: iprint
      integer(4),intent(inout):: ier 

      ! Temporal copy of iuhshk & iuashk
        integer(4):: juhshk(maxshk),juashk(mxashk,maxshk)
      ! Temporal copy of fudshk
        real(8):: gudshk(maxequ,maxshk)

      integer(4):: combin ! (function)
      integer(4):: numshk,numatm,numcon,itmp,numdis,igrp,numsk2
      integer(4):: numsk3,numsk4,jugsk2,jugsk3,jugsk4
 
!*******************************************
 
      ier = 0
      ! Remove SHAKE constraints if atom is fixed
      if (ixnatm .ne. iynvar ) then
        numshk = 0
        do igrp = 1,iugshk
          numatm = iuhshk(igrp) + 1
          numcon = sum(iytvar(iuashk(1:numatm,igrp)))
          if ( numcon .eq. numatm ) then
            numshk = numshk + 1
            iuhshk(numshk) = iuhshk(igrp)
            itmp = iuhshk(igrp) + 1
            iuashk(1:itmp,numshk) = iuashk(1:itmp,igrp)
            numdis = combin(itmp,2)
            fudshk(1:numdis,numshk) = fudshk(1:numdis,igrp)
          endif
        enddo

        if ( numshk .eq. 0 ) then
          write(iprint,*)'ERROR> MODSHK '
          write(iprint,*)'  NUMBER OF ATOM GROUPS OF SHAKE IS ZERO '
          write(iprint,*)' '
          ier = -6 ; return
        else
          iugshk = numshk
          iutshk = 0
          do igrp = 1,iugshk
            iutshk = iutshk + combin(iuhshk(igrp)+1,2)
          enddo
        endif
      endif

      ! Change order of constraints
      !      1       -  IUGSK2     ;   TWO PARTICLES
      !  (IUGSK2+1)  -  IUGSK3     ;   THREE PARTICLES
      !  (IUGSK3+1)  -  IUGSHK     ;   FOUR PARTICLES
 
      ! 2-1) Back-up SHAKE constraints
      juhshk(1:iugshk) = iuhshk(1:iugshk)
      do igrp = 1,iugshk
        itmp = juhshk(igrp) + 1
        juashk(1:itmp,igrp) = iuashk(1:itmp,igrp)
        numdis = combin(itmp,2)
        gudshk(1:numdis,igrp) = fudshk(1:numdis,igrp)
      enddo

      ! 2-2) Change order of SHAKE constraints
      numsk2 = 0 ; numsk3 = 0 ; numsk4 = 0
      do igrp = 1,iugshk
        if ( juhshk(igrp) .eq. 1 ) then
          numsk2 = numsk2 + 1
        elseif ( juhshk(igrp) .eq. 2 ) then
          numsk3 = numsk3 + 1
        elseif ( juhshk(igrp) .eq. 3 ) then
          numsk4 = numsk4 + 1
        endif
      enddo
 
      iugsk2 = numsk2 ; iugsk3 = numsk2 + numsk3
      jugsk2 = 0 ; jugsk3 = iugsk2 ; jugsk4 = iugsk3

      do igrp = 1,iugshk
        if ( juhshk(igrp) .eq. 1 ) then
          jugsk2 = jugsk2 + 1
          iuhshk(jugsk2) = juhshk(igrp)
          iuashk(1:2,jugsk2) = juashk(1:2,igrp)
          fudshk(1,jugsk2) = gudshk(1,igrp)
        elseif ( juhshk(igrp) .eq. 2 ) then
          jugsk3 = jugsk3 + 1
          iuhshk(jugsk3) = juhshk(igrp)
          iuashk(1:3,jugsk3) = juashk(1:3,igrp)
          fudshk(1:3,jugsk3) = gudshk(1:3,igrp)
        elseif ( juhshk(igrp) .eq. 3 ) then
          jugsk4 = jugsk4 + 1
          iuhshk(jugsk4) = juhshk(igrp)
          iuashk(1:4,jugsk4) = juashk(1:4,igrp)
          fudshk(1:6,jugsk4) = gudshk(1:6,igrp)
        endif
      enddo

!************************

      return
      end subroutine modshk


!===================================================================


      subroutine reabnd

!*******************************************************************
!
!     RE-ARRANGE BOND PARAMETER FOR SHAKE
!
!     4 TEMPORARY VARIABLES
!       PARLIS     I*4  (MAXEQU,2)   : ATOM PAIR LIST IN EACH SHAKE
!                                      CONSTRAINT
!                                    PLEASE CHECK SHAKE ROUTINE
!
!*******************************************************************
 
      use COMBAS ; use COMERG ; use COMMIS 

      implicit none

      integer(4):: combin
      real(8):: bkfbnd(iynbnd),bkqbnd(iynbnd)
      integer(4):: intwk(2,iynbnd)
      integer(4):: onbond(iynbnd)
      integer(4):: parlis(6,2) ! (maxequ,2)
      data parlis/ 1,2,3,4,3,2,                                        &
                   2,3,1,1,4,4 /
 
      integer(4):: numshk,igrp,ishk,ibnd,iatm1,iatm2,jatm1,jatm2
      integer(4):: jbnd1,jbnd2,nst
 
!************************************************

      ! 1) Back-up bond parameter
      bkfbnd(1:iynbnd) = fyfbnd(1:iynbnd)
      bkqbnd(1:iynbnd) = fyqbnd(1:iynbnd)
      intwk(1:2,1:iynbnd) = iypbnd(1:2,1:iynbnd)
      onbond(1:iynbnd) = 1
 
      ! 2) Search dummy bond under SHAKE constraints
      !! 2 atom SHAKE
      nst = 1
      do igrp = 1,iugsk2
        numshk = combin(iuhshk(igrp)+1,2)
        do ishk = 1,numshk
          iatm1 = iuashk(parlis(ishk,1),igrp)
          iatm2 = iuashk(parlis(ishk,2),igrp)
          do ibnd = nst,iynbnd
            jatm1 = iypbnd(1,ibnd)
            jatm2 = iypbnd(2,ibnd)
            if ( ( iatm1.eq.jatm1 .and. iatm2.eq.jatm2 ) .or.          &
                 ( iatm1.eq.jatm2 .and. iatm2.eq.jatm1 ) ) then
              onbond(ibnd) = 0
              nst = ibnd + 1
              exit
            elseif ( ( iatm1.lt.jatm1 .and. iatm2.lt.jatm2 ) .or.      &
                     ( iatm1.lt.jatm2 .and. iatm2.lt.jatm1 ) ) then
              exit
            endif
          enddo
        enddo
      enddo
      !! 3 atom SHAKE
      nst = 1
      do igrp = iugsk2+1,iugsk3
        numshk = combin(iuhshk(igrp)+1,2)
        do ishk = 1,numshk
          iatm1 = iuashk(parlis(ishk,1),igrp)
          iatm2 = iuashk(parlis(ishk,2),igrp)
          do ibnd = nst,iynbnd
            jatm1 = iypbnd(1,ibnd)
            jatm2 = iypbnd(2,ibnd)
            if ( ( iatm1.eq.jatm1 .and. iatm2.eq.jatm2 ) .or.          &
                 ( iatm1.eq.jatm2 .and. iatm2.eq.jatm1 ) ) then
              onbond(ibnd) = 0
              nst = max(ibnd-1,1) !!(For water molecules)
              exit
            elseif ( ( iatm1.lt.jatm1 .and. iatm2.lt.jatm2 ) .or.      &
                     ( iatm1.lt.jatm2 .and. iatm2.lt.jatm1 ) ) then
              exit
            endif
          enddo
        enddo
      enddo
      !! 4 atom SHAKE
      nst = 1
      do igrp = iugsk3+1,iugshk
        numshk = combin(iuhshk(igrp)+1,2)
        do ishk = 1,numshk
          iatm1 = iuashk(parlis(ishk,1),igrp)
          iatm2 = iuashk(parlis(ishk,2),igrp)
          do ibnd = nst,iynbnd
            jatm1 = iypbnd(1,ibnd)
            jatm2 = iypbnd(2,ibnd)
            if ( ( iatm1.eq.jatm1 .and. iatm2.eq.jatm2 ) .or.          &
                 ( iatm1.eq.jatm2 .and. iatm2.eq.jatm1 ) ) then
              onbond(ibnd) = 0
              nst = ibnd + 1
              exit
            elseif ( ( iatm1.lt.jatm1 .and. iatm2.lt.jatm2 ) .or.      &
                     ( iatm1.lt.jatm2 .and. iatm2.lt.jatm1 ) ) then
              exit
            endif
          enddo
        enddo
      enddo

      iyndbn = iynbnd - sum(onbond(1:iynbnd))
  
      ! 3) Re-arrange bond parameter
      jbnd1 = 0 ; jbnd2 = iynbnd - iyndbn
      do ibnd = 1,iynbnd
        if ( onbond(ibnd) .eq. 1 ) then
          jbnd1 = jbnd1 + 1
          fyfbnd(jbnd1) = bkfbnd(ibnd)
          fyqbnd(jbnd1) = bkqbnd(ibnd)
          iypbnd(1:2,jbnd1) = intwk(1:2,ibnd)
        else
          jbnd2 = jbnd2 + 1
          fyfbnd(jbnd2) = bkfbnd(ibnd)
          fyqbnd(jbnd2) = bkqbnd(ibnd)
          iypbnd(1:2,jbnd2) = intwk(1:2,ibnd)
        endif
      enddo

!***********************************

      return
      end subroutine reabnd
 

!=========================================================================

      subroutine reaang

!*******************************************************************
!
!     RE-ARRANGE ANGLE PARAMETER FOR SHAKE
!
!       PARLIS     I*4  (MAXEQU,2)   : ATOM PAIR LIST IN EACH SHAKE
!                                      CONSTRAINT
!                                    PLEASE CHECK SHAKE ROUTINE
!
!*******************************************************************
 
      use COMBAS ; use COMERG ; use COMMIS 

      implicit none

      integer(4):: combin
      real(8):: bkfang(iynang),bkqang(iynang)
      integer(4):: intwrk(3,iynang),onangl(iynang)
      integer(4):: parlis(6,2) ! (maxequ,2)
      data parlis/ 1,2,3,4,3,2,                                        &
                   2,3,1,1,4,4 /

      integer(4):: numshk,iang,igrp,ishk,iatm1,iatm2,jatm1,jatm2
      integer(4):: jang1,jang2,nst

!***********************************

      ! 1) Back-up angle parameter
      bkfang(1:iynang) = fyfang(1:iynang)
      bkqang(1:iynang) = fyqang(1:iynang)
      intwrk(1:3,1:iynang) = iypang(1:3,1:iynang)
      onangl(1:iynang) = 1

      ! 2) Search dummy angle under SHAKE constraints
      !! 3 atom SHAKE
      nst = 1
      do igrp = iugsk2+1,iugsk3
        numshk = combin(iuhshk(igrp)+1,2)
        do ishk = 1,numshk
          iatm1 = iuashk(parlis(ishk,1),igrp)
          iatm2 = iuashk(parlis(ishk,2),igrp)
          do iang = nst,iynang
            jatm1 = iypang(1,iang)
            jatm2 = iypang(3,iang)
            if ( ( iatm1.eq.jatm1 .and. iatm2.eq.jatm2 ) .or.          &
                 ( iatm1.eq.jatm2 .and. iatm2.eq.jatm1 ) ) then
              onangl(iang) = 0
              nst = max(iang-1,1)  !! (For water molecules)
              exit
            elseif ( ( iatm1.lt.jatm1 .and. iatm2.lt.jatm2 ) .or.      &
                     ( iatm1.lt.jatm2 .and. iatm2.lt.jatm1 ) ) then
              exit
            endif
          enddo
        enddo
      enddo
      !! 4 atom SHAKE
      nst = 1
      do igrp = iugsk3+1,iugshk
        numshk = combin(iuhshk(igrp)+1,2)
        do ishk = 1,numshk
          iatm1 = iuashk(parlis(ishk,1),igrp)
          iatm2 = iuashk(parlis(ishk,2),igrp)
          do iang = nst,iynang
            jatm1 = iypang(1,iang)
            jatm2 = iypang(3,iang)
            if ( ( iatm1.eq.jatm1 .and. iatm2.eq.jatm2 ) .or.          &
                 ( iatm1.eq.jatm2 .and. iatm2.eq.jatm1 ) ) then
              onangl(iang) = 0
              nst = max(iang-2,1) !! (For ACE)
              exit
            elseif ( ( iatm1.lt.jatm1 .and. iatm2.lt.jatm2 ) .or.      &
                     ( iatm1.lt.jatm2 .and. iatm2.lt.jatm1 ) ) then
              exit
            endif
          enddo
        enddo
      enddo

      iyndag = iynang - sum(onangl(1:iynang))

      ! 3) Re-arrange angle parameter 
      jang1 = 0 ; jang2 = iynang - iyndag
      do iang = 1,iynang
        if ( onangl(iang) .eq. 1 ) then
          jang1 = jang1 + 1
          fyfang(jang1) = bkfang(iang) 
          fyqang(jang1) = bkqang(iang) 
          iypang(1:3,jang1) = intwrk(1:3,iang)
        else 
          jang2 = jang2 + 1
          fyfang(jang2) = bkfang(iang) 
          fyqang(jang2) = bkqang(iang) 
          iypang(1:3,jang2) = intwrk(1:3,iang)
        endif
      enddo

!**************************************

      return
      end subroutine reaang
