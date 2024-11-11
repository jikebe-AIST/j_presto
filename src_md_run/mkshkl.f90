
      program mkshkl

!******************************************************************
!
!     THIS PROGRAM IS FOR MAKING SHAKE CONTROL DATA
!
!******************************************************************

      use COMBAS ; use COMERG ; use COMMIS ; use PHYCNS 

      implicit none

      integer(4):: iread,iprint,icard,iwrite,iopt,nummol,imol,ier,     &
                   numgrp
      integer(4),allocatable:: shknum(:)
      integer(4),allocatable:: bnugrp(:),topbnd(:),botbnd(:),topang(:),&
                               botang(:),topchn(:)
      character(80):: file1,file2
 
!**********************************************************

      iread = 5 ; iprint = 6 ; icard = 10 ; iwrite = 20
 
!     1) TERMINAL INPUT
      write(iprint,*)' '
      write(iprint,*)'QUESTION> MKSHKL '
      write(iprint,*)'    ALL BONDS SHAKE OR BONDS AND ANGLES '
      write(iprint,*)'    INCLUDING HYDROGEN ATOMS SHAKE '
      write(iprint,*)' '
      write(iprint,*)' 1) BONDS AND ANGLES INCLUDING HYDROGEN '
      write(iprint,*)'    ATOMS SHAKE (HYDROGEN SHAKE)'
      write(iprint,*)' 2) ALL BONDS SHAKE (FIRST CHAINS)'
      write(iprint,*)'    HYDROGEN SHAKE (OTHER CHAINS) '
      write(iprint,*)' '
      write(iprint,*)'    PLEASE SELECT 1 OR 2 '
      write(iprint,*)' '
      read(iread,*)iopt

      if ( iopt .eq. 1 ) then
        nummol = 0
      elseif ( iopt .eq. 2 ) then
        write(iprint,*)' '
        write(iprint,*)'QUESTION> MKSHKL '
        write(iprint,*)' HOW MANY MOLECULES ARE TREATED '
        write(iprint,*)' AS ALL-BONDS SHAKE'
        write(iprint,*)' '
        write(iprint,*)'  IF YOU SPECIFY 2 THEN FIRST AND SECOND '
        write(iprint,*)'  MOLECULES ARE TREATED AS ALL-BONDS SHAKE '
        write(iprint,*)'  THIRD MOLECULE AND ... ARE TREATED AS '
        write(iprint,*)'  HYDROGEN SHAKE '
        write(iprint,*)' '
        read(iread,*)nummol
      else
        write(iprint,*)' '
        write(iprint,*)'ERROR> MKSHKL '
        write(iprint,*)'     OPTION IS 1 OR 2 '
        write(iprint,*)'     YOUR SPECIFIED VALUE IS ',iopt
        write(iprint,*)' '
        stop
      endif

      write(iprint,*)'  TOPOLOGY FILE (V2.0) <INPUT>  ? '
      read(iread,'(a80)')file1
      write(iprint,*)'  SHAKE CONTROL FILE   <OUTPUT> ? '
      read(iread,'(a80)')file2
 
!     2) read TOPOLOGY FILE
      call inptpl(icard,iprint,file1,ier)
      allocate(bnugrp(ixmolc),topbnd(ixmolc),botbnd(ixmolc),           &
               topang(ixmolc),botang(ixmolc),topchn(ixmolc))
      bnugrp(1:ixmolc) = 0
      if ( ier .ne. 0 ) stop
 
!     3) SEARCH TOP AND BOTTOM BOND NUMBER IN FIRST CHAIN OF EACH
!        MOLECULE
!        SEARCH TOP AND BOTTOM ANGLE NUMBER IN FIRST CHAIN OF EACH
!        MOLECULE
!        SEARCH FIRST CHAIN NUMBER OF EACH MOLECULE
      write(iprint,*)' '
      write(iprint,*)'INFORMATION> MKSHKL '
      write(iprint,*)'  1) SEARCH TOP AND BOTTOM BOND NUMBER IN FIRST'
      write(iprint,*)'  CHAIN OF EACH MOLECULE '
      write(iprint,*)'  2) SEARCH TOP AND BOTTOM ANGLE NUMBER IN FIRST'
      write(iprint,*)'  CHAIN OF EACH MOLECULE '
      write(iprint,*)'  3) SEARCH FIRST CHAIN NUMBER OF EACH MOLECULE '
      write(iprint,*)' '
      call srbnag(maxatm,maxbnd,maxang,ixmolc,iynbnd,iynang,ixamol,    &
                  ixachn,iypbnd,iypang,ixsqml,topbnd,botbnd,topang,    &
                  botang,topchn)
 
!     4) MAKE SHAKE CONTROL DATA (ALL BONDS SHAKE)
!        1-ST MOLECULE TO  nummol-TH MOLECULE
      iugshk = 0
      if ( iopt .eq. 2 ) then
        call mkalls(ixmolc,maxatm,maxbnd,maxshk,mxashk,iynbnd,nummol,  &
                    iugshk,ixamol,ixachn,iypbnd,iuhshk,iuashk,topchn,  &
                    bnugrp,iprint,ier)
        if ( ier .ne. 0 ) stop
      endif

      if ( iopt.eq.1 .or. (iopt.eq.2 .and. ixmolc.gt.nummol) ) then
!       5) CALCULATE NUMBER OF ATOM GROUPS FOR SHAKE IN EACH MOLECULE
        write(iprint,*)' '
        write(iprint,*)'INFORMATION> MKSHKL '
        write(iprint,*)'    CALCULATE NUMBER OF ATOM GROUPS FOR SHAKE '
        write(iprint,*)'    IN EACH MOLECULE '
        write(iprint,*)' '
        allocate(shknum(maxatm))
        call calgrp(ixmolc,maxatm,maxbnd,maxshk,mxashk,ixnatm,iynbnd,  &
                    iugshk,ixamol,ixachn,cxatmn,iypbnd,iuhshk,iuashk,  &
                    nummol,shknum,bnugrp,topchn,iprint,ier)
      endif
 
      write(iprint,*)' '
      write(iprint,*)'INFORMATION> MKSHKL '
      write(iprint,*)'     NUMBER OF ATOM GROUPS FOR SHAKE : ',iugshk
      if ( iugshk .eq. 0 ) then
        write(iprint,*)' '
        write(iprint,*)'ERROR> MKSHKL '
        write(iprint,*)'  ATOM GROUPS OF SHAKE CONSTARINTS ARE NOT'
        write(iprint,*)'  FOUND '
        write(iprint,*)' '
        stop
      endif
      do imol = 1,ixmolc
        if ( imol .eq. 1 ) then
          write(iprint,*)'     MOLECULE ',imol,' - ',bnugrp(imol)
        else
!         If imol th molecule does not have any bond (such as metal ion), 
!         it is not counted in calgrp. So set it properly.
          if ( bnugrp(imol) .eq. 0 ) bnugrp(imol) = bnugrp(imol-1)
          numgrp = bnugrp(imol) - bnugrp(imol-1)
          write(iprint,*)'     MOLECULE ',imol,' - ',numgrp
        endif
      enddo
      write(iprint,*)' '
 
!     6) CALCULATE BOND LENGTH FOR SHAKE
      write(iprint,*)' '
      write(iprint,*)'INFORMATION> MKSHKL '
      write(iprint,*)'    CALCULATE BOND LENGTH FOR SHAKE '
      write(iprint,*)' '
      call calbln(ixmolc,maxatm,maxbnd,maxang,maxshk,mxashk,maxequ,rad,&
                  iugshk,ixamol,cxatmn,cxresn,iypbnd,iypang,fyqbnd,    &
                  fyqang,iuhshk,iuashk,fudshk,topbnd,botbnd,topang,    &
                  botang,iprint,iread)
 
!     7) OUTPUT SHAKE LIST
      write(iprint,*)' '
      write(iprint,*)'INFORMATION> MKSHKL '
      write(iprint,*)'    OUTPUT SHAKE LIST '
      write(iprint,*)' '
      call outshk(ixnchn,maxatm,maxshk,mxashk,maxequ,ixmolc,cxmolc,    &
                  ixcend,ixares,cxatmn,cxresn,iuhshk,iuashk,fudshk,    &
                  bnugrp,topchn,iwrite,iprint,ier,file2)

!****************************************

      stop
      end program mkshkl


!======================================================================
 

      subroutine mkalls(ixmolc,maxatm,maxbnd,maxshk,mxashk,iynbnd,     &
                        nummol,iugshk,ixamol,ixachn,iypbnd,iuhshk,     &
                        iuashk,topchn,bnugrp,iprint,ier)

      implicit none

      integer(4):: ixmolc,maxatm,maxbnd,maxshk,mxashk,iynbnd,nummol,   &
                   iugshk,ixamol(maxatm),ixachn(maxatm),               &
                   iypbnd(2,maxbnd),iuhshk(maxshk),                    &
                   iuashk(mxashk,maxshk),topchn(ixmolc),bnugrp(ixmolc),&
                   iprint,ier
      integer(4):: itmp,ibnd,imol,ichn 

!***************************************************

      do ibnd = 1,iynbnd
        itmp = iypbnd(1,ibnd)
        imol = ixamol(itmp) ; ichn = ixachn(itmp)
        if ( imol.le.nummol .and. ichn.eq.topchn(imol) ) then
          iugshk = iugshk + 1
          if ( iugshk .gt. maxshk ) then
            write(iprint,*)' '
            write(iprint,*)'ERROR> MKSHKL '
            write(iprint,*)'   NUMBER OF ATOM GROUPS ARE EXCEEDED '
            write(iprint,*)'   THE LIMIT'
            write(iprint,*)' '
            write(iprint,*)'   THE LIMIT NUMBER IS ',maxshk
            write(iprint,*)' '
            ier = -1 ; return
          endif
          iuhshk(iugshk) = 1
          iuashk(1:2,iugshk) = iypbnd(1:2,ibnd)
          bnugrp(imol) = iugshk
        endif
      enddo

!********************************

      return
      end subroutine mkalls


!================================================================
 

      subroutine srbnag(maxatm,maxbnd,maxang,ixmolc,iynbnd,iynang,     &
                        ixamol,ixachn,iypbnd,iypang,ixsqml,topbnd,     &
                        botbnd,topang,botang,topchn)
 
      implicit none

      integer(4):: maxatm,maxbnd,maxang,ixmolc,iynbnd,iynang
      integer(4):: ixamol(maxatm),ixachn(maxatm),iypbnd(2,maxbnd),     &
                   iypang(3,maxang),ixsqml(ixmolc),topbnd(ixmolc),     &
                   botbnd(ixmolc),topang(ixmolc),botang(ixmolc),       &
                   topchn(ixmolc)
      integer(4):: imol,ibnd,iatm,ichn,iang,itmp

!************************************************
 
!     <<<  INITILIZATION  >>>
      topbnd(1:ixmolc) = iynbnd + 1 ; topang(1:ixmolc) = iynang + 1
      botbnd(1:ixmolc) = -1 ; botang(1:ixmolc) = -1
      topchn(1) = 1
      do imol = 2,ixmolc
        topchn(imol) = topchn(imol-1) + ixsqml(imol-1)
      enddo
 
!     <<<  SEARCH TOP BOND NUMBER AND BOTTOM BOND NUMBER >>>
      do ibnd = 1,iynbnd
        iatm = iypbnd(1,ibnd)
        imol = ixamol(iatm) ; ichn = ixachn(iatm)
        if ( ichn .eq. topchn(imol) ) then
          topbnd(imol) = min(topbnd(imol),ibnd)
          botbnd(imol) = max(botbnd(imol),ibnd)
        endif
      enddo

!     <<<  SEARCH TOP ANGLE NUMBER AND BOTTOM ANGLE NUMBER >>>
      do iang = 1,iynang
        iatm = iypang(1,iang)
        imol = ixamol(iatm) ; ichn = ixachn(iatm)
        if ( ichn .eq. topchn(imol) ) then
          topang(imol) = min(topang(imol),iang)
          botang(imol) = max(botang(imol),iang)
        endif
      enddo
 
!     <<<  MOLECULE HAS NO BOND OR NO ANGLE  >>>
      itmp = iynbnd + 1
      where ( topbnd .eq. itmp ) topbnd = 0
      itmp = iynang + 1
      where ( topang .eq. itmp ) topang = 0
      where ( botbnd .eq. -1 ) botbnd = 0
      where ( botang .eq. -1 ) botang = 0

!***************************************************

      return
      end subroutine srbnag


!========================================================================
 

      subroutine calgrp(ixmolc,maxatm,maxbnd,maxshk,mxashk,ixnatm,     &
                        iynbnd,iugshk,ixamol,ixachn,cxatmn,iypbnd,     &
                        iuhshk,iuashk,nummol,shknum,bnugrp,topchn,     &
                        iprint,ier)

      implicit none

      integer(4):: ixmolc,maxatm,maxbnd,maxshk,mxashk,ixnatm,iynbnd,   &
                   iugshk,nummol,iprint,ier
      integer(4):: ixamol(maxatm),ixachn(maxatm),iypbnd(2,maxbnd),     &
                   iuhshk(maxshk),iuashk(mxashk,maxshk),shknum(maxatm),&
                   bnugrp(ixmolc),topchn(ixmolc)
      character(8):: cxatmn(maxatm)
      integer(4):: ibnd,iatm(2),imol,ichn,num

!********************************************************

      shknum(1:ixnatm) = 0
!     <<<  SEARCH ATOM GROUPS INCLUDING HYDROGEN ATOMS  >>>
      do ibnd = 1,iynbnd
        iatm(1:2) = iypbnd(1:2,ibnd)
        imol = ixamol(iatm(1)) ; ichn = ixachn(iatm(1))
        if ( imol .le. nummol ) cycle

        if ( ichn .eq. topchn(imol) ) then
          if ( cxatmn(iatm(1))(1:1).eq."H" .and.                       &
               cxatmn(iatm(2))(1:1).ne."H" ) then
            num = shknum(iatm(2))
            if ( num .ne. 0 ) then
              iuhshk(num) = iuhshk(num) + 1
              iuashk(iuhshk(num)+1,num) = iatm(1)
            else
              iugshk = iugshk + 1
              if ( iugshk .gt. maxshk ) then
                ier = -1 ; exit
              endif
              shknum(iatm(2)) = iugshk
              iuashk(1,iugshk) = iatm(2) ; iuashk(2,iugshk) = iatm(1)
              iuhshk(iugshk) = 1
            endif
          elseif ( cxatmn(iatm(2))(1:1).eq."H" .and.                   &
                   cxatmn(iatm(1))(1:1).ne."H" ) then
            num = shknum(iatm(1))
            if ( num .ne. 0 ) then
              iuhshk(num) = iuhshk(num) + 1
              iuashk(iuhshk(num)+1,num) = iatm(2)
            else
              iugshk = iugshk + 1
              if ( iugshk .gt. maxshk ) then
                ier = -1 ; exit
              endif
              shknum(iatm(1)) = iugshk
              iuashk(1:2,iugshk) = iatm(1:2)
              iuhshk(iugshk) = 1
            endif
          endif
          bnugrp(imol) = iugshk
        endif
      enddo

      if ( ier .eq. -1 ) then
        write(iprint,*)' '
        write(iprint,*)'ERROR> MKSHKL '
        write(iprint,*)'  NUMBER OF ATOM GROUPS OF SHAKE CONSTRAINTS'
        write(iprint,*)'  IS EXCEEDED THE LIMIT '
        write(iprint,*)'   '
        write(iprint,*)'     LIMIT NUMBER IS ',maxshk
        write(iprint,*)' '
        return
      endif

!**********************************************

      return
      end subroutine calgrp


!====================================================================


      subroutine calbln(ixmolc,maxatm,maxbnd,maxang,maxshk,mxashk,     &
                        maxequ,rad,iugshk,ixamol,cxatmn,cxresn,iypbnd, &
                        iypang,fyqbnd,fyqang,iuhshk,iuashk,fudshk,     &
                        topbnd,botbnd,topang,botang,iprint,iread)

      implicit none

      integer(4):: ixmolc,maxatm,maxbnd,maxang,maxshk,mxashk,maxequ,   &
                   iugshk,iprint,iread
      integer(4):: ixamol(maxatm),iypbnd(2,maxbnd),iypang(3,maxang),   &
                   iuhshk(maxshk),iuashk(mxashk,maxshk),topbnd(ixmolc),&
                   botbnd(ixmolc),topang(ixmolc),botang(ixmolc)
      real(8):: rad,fyqbnd(maxbnd),fyqang(maxang),fudshk(maxequ,maxshk)
      character(8):: cxatmn(maxatm),cxresn(maxatm)
      integer(4):: igrp,iatmc,imol,istbnd,ienbnd,numcon,iatmh,iatmh1,  &
                   iatmh2,iatmh3,ier,istang,ienang

!*********************************************

      do igrp = 1,iugshk
        iatmc = iuashk(1,igrp) ; imol = ixamol(iatmc)
        istbnd = topbnd(imol) ; ienbnd = botbnd(imol)
        numcon = iuhshk(igrp)
!       1) TWO PARTICLES (ONE BOND)
        if ( numcon .eq. 1 ) then
          iatmh = iuashk(2,igrp)
          call srbnd(istbnd,ienbnd,iatmc,iatmh,maxbnd,iypbnd,fyqbnd,   &
                     fudshk(1,igrp),ier)
          if ( ier .ne. 0 ) then
            call misbnd(cxatmn(iatmc),cxatmn(iatmh),cxresn(iatmc),     &
                        cxresn(iatmh),fudshk(1,igrp),iprint,iread)
          endif

!       2) THREE PARTICLES (THREE BONDS)
        elseif ( numcon .eq. 2 ) then
          iatmh1 = iuashk(2,igrp) ; iatmh2 = iuashk(3,igrp)

!         2-1) IUASHK(IGRP,1) - IUASHK(IGRP,2)
          call srbnd(istbnd,ienbnd,iatmc,iatmh1,maxbnd,iypbnd,fyqbnd,  &
                     fudshk(1,igrp),ier)
          if ( ier .ne. 0 ) then
            call misbnd(cxatmn(iatmc),cxatmn(iatmh1),cxresn(iatmc),    &
                        cxresn(iatmh1),fudshk(1,igrp),iprint,iread)
          endif
 
!         2-2) IUASHK(IGRP,1) - IUASHK(IGRP,3)
          call srbnd(istbnd,ienbnd,iatmc,iatmh2,maxbnd,iypbnd,fyqbnd,  &
                     fudshk(3,igrp),ier)
          if ( ier .ne. 0 ) then
            call misbnd(cxatmn(iatmc),cxatmn(iatmh2),cxresn(iatmc),    &
                        cxresn(iatmh2),fudshk(3,igrp),iprint,iread)
          endif
 
!         2-3) IUASHK(IGRP,2) - IUASHK(IGRP,3)
!              THIS BOND SOMTIMES MISSING
          call srbnd(istbnd,ienbnd,iatmh1,iatmh2,maxbnd,iypbnd,fyqbnd, &
                     fudshk(2,igrp),ier)
          if ( ier .ne. 0 ) then
            istang = topang(imol) ; ienang = botang(imol)
            call srang(istang,ienang,iatmc,iatmh1,iatmh2,istbnd,ienbnd,&
                       cxatmn(iatmc),cxatmn(iatmh1),cxatmn(iatmh2),    &
                       cxresn(iatmc),cxresn(iatmh1),cxresn(iatmh2),    &
                       maxbnd,iypbnd,fyqbnd,maxang,iypang,fyqang,rad,  &
                       fudshk(2,igrp),iprint,iread)
          endif
 
!       3) FOUR PARTICLES (SIX BONDS)
        elseif ( numcon .eq. 3 ) then
          iatmh1 = iuashk(2,igrp) ; iatmh2 = iuashk(3,igrp)
          iatmh3 = iuashk(4,igrp)

!         3-1) IUASHK(IGRP,1) - IUASHK(IGRP,2)
          call srbnd(istbnd,ienbnd,iatmc,iatmh1,maxbnd,iypbnd,fyqbnd,  &
                     fudshk(1,igrp),ier)
          if ( ier .ne. 0 ) then
            call misbnd(cxatmn(iatmc),cxatmn(iatmh1),cxresn(iatmc),    &
                        cxresn(iatmh1),fudshk(1,igrp),iprint,iread)
          endif

!         3-2) IUASHK(IGRP,1) - IUASHK(IGRP,3)
          call srbnd(istbnd,ienbnd,iatmc,iatmh2,maxbnd,iypbnd,fyqbnd,  &
                     fudshk(3,igrp),ier)
          if ( ier .ne. 0 ) then
            call misbnd(cxatmn(iatmc),cxatmn(iatmh2),cxresn(iatmc),    &
                        cxresn(iatmh2),fudshk(3,igrp),iprint,iread)
          endif
 
!         3-3) IUASHK(IGRP,1) - IUASHK(IGRP,4)
          call srbnd(istbnd,ienbnd,iatmc,iatmh3,maxbnd,iypbnd,fyqbnd,  &
                     fudshk(4,igrp),ier)
          if ( ier .ne. 0 ) then
            call misbnd(cxatmn(iatmc),cxatmn(iatmh3),cxresn(iatmc),    &
                        cxresn(iatmh3),fudshk(4,igrp),iprint,iread)
          endif
 
!         3-4) IUASHK(IGRP,2) - IUASHK(IGRP,3)
!              THIS BOND SOMTIMES MISSING
          call srbnd(istbnd,ienbnd,iatmh1,iatmh2,maxbnd,iypbnd,fyqbnd, &
                     fudshk(2,igrp),ier)
          if ( ier .ne. 0 ) then
            istang = topang(imol) ; ienang = botang(imol)
            call srang(istang,ienang,iatmc,iatmh1,iatmh2,istbnd,ienbnd,&
                       cxatmn(iatmc),cxatmn(iatmh1),cxatmn(iatmh2),    &
                       cxresn(iatmc),cxresn(iatmh1),cxresn(iatmh2),    &
                       maxbnd,iypbnd,fyqbnd,maxang,iypang,fyqang,rad,  &
                       fudshk(2,igrp),iprint,iread)
          endif

!         3-5) IUASHK(IGRP,3) - IUASHK(IGRP,4)
!              THIS BOND SOMTIMES MISSING
          call srbnd(istbnd,ienbnd,iatmh2,iatmh3,maxbnd,iypbnd,fyqbnd, &
                     fudshk(5,igrp),ier)
          if ( ier .ne. 0 ) then
            istang = topang(imol) ; ienang = botang(imol)
            call srang(istang,ienang,iatmc,iatmh2,iatmh3,istbnd,ienbnd,&
                       cxatmn(iatmc),cxatmn(iatmh2),cxatmn(iatmh3),    &
                       cxresn(iatmc),cxresn(iatmh2),cxresn(iatmh3),    &
                       maxbnd,iypbnd,fyqbnd,maxang,iypang,fyqang,rad,  &
                       fudshk(5,igrp),iprint,iread)
          endif
 
!         3-6) IUASHK(IGRP,4) - IUASHK(IGRP,2)
!              THIS BOND SOMTIMES MISSING
          call srbnd(istbnd,ienbnd,iatmh1,iatmh3,maxbnd,iypbnd,fyqbnd, &
                     fudshk(6,igrp),ier)
          if ( ier .ne. 0 ) then
            istang = topang(imol) ; ienang = botang(imol)
            call srang(istang,ienang,iatmc,iatmh1,iatmh3,istbnd,ienbnd,&
                       cxatmn(iatmc),cxatmn(iatmh1),cxatmn(iatmh3),    &
                       cxresn(iatmc),cxresn(iatmh1),cxresn(iatmh3),    &
                       maxbnd,iypbnd,fyqbnd,maxang,iypang,fyqang,rad,  &
                       fudshk(6,igrp),iprint,iread)
          endif

        endif

      enddo

!*********************************************

      return
      end subroutine calbln


!============================================================================


      subroutine srbnd(istbnd,ienbnd,iatm1,iatm2,maxbnd,iypbnd,fyqbnd, &
                       bonlen,ier)

      implicit none

      integer(4):: istbnd,ienbnd,iatm1,iatm2,maxbnd,ier
      integer(4):: iypbnd(2,maxbnd)
      real(8):: fyqbnd(maxbnd),bonlen
      integer(4):: ibnd,jatm1,jatm2

!*******************************************************

      do ibnd = istbnd,ienbnd
        jatm1 = iypbnd(1,ibnd) ; jatm2 = iypbnd(2,ibnd)
        if ( (iatm1.eq.jatm1 .and. iatm2.eq.jatm2) .or.                &
             (iatm1.eq.jatm2 .and. iatm2.eq.jatm1) ) then
          bonlen = fyqbnd(ibnd) ; ier = 0 ; return
        endif
      enddo
      ier = - 1

!**********************************

      return
      end subroutine srbnd


!====================================================================
 

      subroutine misbnd(atty1,atty2,resty1,resty2,bonlen,iprint,iread) 

      implicit none

      character(8):: atty1,atty2,resty1,resty2
      real(8):: bonlen
      integer(4):: iprint,iread
 
!**************************************************
 
      write(iprint,*)'WARNING> MKSHKL'
      write(iprint,*)'  BOND LENGTH IS MISSING '
      write(iprint,*)'     BETWEEN ',atty1,' IN ',resty1,' AND ',atty2,&
                     ' IN ',resty2
      write(iprint,*)'  BOND LENGTH ? '
      read(iread,*)bonlen

!***********************************

      return 
      end subroutine misbnd


!========================================================================
 

      subroutine srang(istang,ienang,iatmc,iatm1,iatm2,istbnd,ienbnd,  &
                       attyc,atty1,atty2,restyc,resty1,resty2,maxbnd,  &
                       iypbnd,fyqbnd,maxang,iypang,fyqang,rad,bonlen,  &
                       iprint,iread)

      implicit none

      integer(4):: istang,ienang,iatmc,iatm1,iatm2,istbnd,ienbnd,      &
                   maxbnd,maxang,iprint,iread
      integer(4):: iypbnd(2,maxbnd),iypang(3,maxang)
      real(8):: fyqbnd(maxbnd),fyqang(maxang),rad,bonlen,angle,dist1,  &
                dist2
      character(8):: attyc,atty1,atty2,restyc,resty1,resty2
      integer(4):: iang,jatmc,jatm1,jatm2,ierr,ier

!**********************************************************

      ierr = -1
      do iang = istang,ienang
        jatmc = iypang(2,iang) ; jatm1 = iypang(1,iang)
        jatm2 = iypang(3,iang)
        if ( iatmc.eq.jatmc .and.                                      &
             ( (iatm1.eq.jatm1 .and. iatm2.eq.jatm2) .or.              &
               (iatm1.eq.jatm2 .and. iatm2.eq.jatm1) ) ) then
          angle = fyqang(iang) ; ierr = 0 ; exit
        endif
      enddo

      if ( ierr .eq. -1 ) then
        write(iprint,*)'WARNING> MKSHKL'
        write(iprint,*)'  ANGLE IS MISSING '
        write(iprint,*)'     BETWEEN ',atty1,' IN ',resty1,' - ',attyc,&
                       ' IN ',restyc,' - ',atty2,' IN ',resty2
        write(iprint,*)'  ANGLE (DEGREE) ? '
        read(iread,*)angle
        angle  = angle * rad
      endif
 
      call srbnd(istbnd,ienbnd,iatmc,iatm1,maxbnd,iypbnd,fyqbnd,dist1, &
                 ier)
      if ( ier .ne. 0 ) call misbnd(attyc,atty1,restyc,resty1,dist1,   &
                                    iprint,iread)
      call srbnd(istbnd,ienbnd,iatmc,iatm2,maxbnd,iypbnd,fyqbnd,dist2, &
                 ier)
      if ( ier .ne. 0 ) call misbnd(attyc,atty2,restyc,resty2,dist2,   &
                                    iprint,iread)

      bonlen = dist1*dist1 + dist2*dist2 - 2.d0*dist1*dist2*cos(angle)
      bonlen = sqrt(bonlen)

!*****************************************************

      return
      end subroutine srang


!======================================================================
 

      subroutine outshk(ixnchn,maxatm,maxshk,mxashk,maxequ,ixmolc,     &
                        cxmolc,ixcend,ixares,cxatmn,cxresn,iuhshk,     &
                        iuashk,fudshk,bnugrp,topchn,iwrite,iprint,ier, &
                        filena)

      implicit none

      integer(4):: ixnchn,maxatm,maxshk,mxashk,maxequ,iwrite,   &
                   ixmolc,iprint,ier
      integer(4):: ixcend(ixnchn),ixares(maxatm),iuhshk(maxshk),       &
                   iuashk(mxashk,maxshk),bnugrp(ixmolc),topchn(ixmolc)
      real(8):: fudshk(maxequ,maxshk)
      character(40):: cxmolc(ixmolc)
      character(8):: cxatmn(maxatm),cxresn(maxatm)
      character(80):: filena
      integer(4):: imol,igrp,numgrp,numatm,numdis,istgrp,iminus,iengrp,&
                   lstchn

!*************************************************************

      call flopen(iwrite,filena,12,'ZERO',0,ier)
      if ( ier .ne. 0 ) then
        write(iprint,*)' '
        write(iprint,*)'ERROR> MKSHKL '
        write(iprint,*)'     FILE OPEN ERROR '
        write(iprint,*)' '
        ier = -1 ; return
      endif

      do imol = 1,ixmolc
        if ( imol .eq. 1 ) then
          istgrp = 1 ; iminus = 0
        else
          istgrp = bnugrp(imol-1) + 1 ; lstchn = topchn(imol) - 1
          iminus = ixcend(lstchn)
        endif
        iengrp = bnugrp(imol) ; numgrp = iengrp - istgrp + 1

        if ( numgrp .gt. 0 ) then
          write(iwrite,8000)cxmolc(imol),numgrp
          do igrp = istgrp,iengrp
            numatm = iuhshk(igrp) + 1 ; numdis = (numatm*(numatm-1))/2
            if ( numatm .eq. 2 ) then
              write(iwrite,8100)numatm,iuashk(1:numatm,igrp)-iminus,   &
                  cxatmn(iuashk(1:numatm,igrp)),cxresn(iuashk(1,igrp)),&
                  ixares(iuashk(1,igrp)),fudshk(1:numdis,igrp)
            elseif ( numatm .eq. 3 ) then
              write(iwrite,8200)numatm,iuashk(1:numatm,igrp)-iminus,   &
                  cxatmn(iuashk(1:numatm,igrp)),cxresn(iuashk(1,igrp)),&
                  ixares(iuashk(1,igrp)),fudshk(1:numdis,igrp)
            elseif ( numatm .eq. 4 ) then
              write(iwrite,8300)numatm,iuashk(1:numatm,igrp)-iminus,   &
                  cxatmn(iuashk(1:numatm,igrp)),cxresn(iuashk(1,igrp)),&
                  ixares(iuashk(1,igrp)),fudshk(1:numdis,igrp)
            endif
          enddo
        endif
      enddo
 
8000  format(' '/                                                      &
             'SHAKE> SHAKE '/a40/';  NUMBER OF ATOM GROUPS = ',i10/    &
             ';'/                                                      &
             ';    ORDER OF BOND '/                                    &
             ';      1) NUMBER OF ATOMS = 2 '/                         &
             ';         BOND NUM.   ATOM NUM. '/                       &
             ';           1         1 - 2 '/                           &
             ';      2) NUMBER OF ATOMS = 3 '/                         &
             ';         BOND NUM.   ATOM NUM. '/                       &
             ';           1         1 - 2 '/                           &
             ';           2         2 - 3 '/                           &
             ';           3         1 - 3 '/                           &
             ';      3) NUMBER OF ATOMS = 4 '/                         &
             ';         BOND NUM.   ATOM NUM. '/                       &
             ';           1         1 - 2 '/                           &
             ';           2         2 - 3 '/                           &
             ';           3         1 - 3 '/                           &
             ';           4         1 - 4 '/                           &
             ';           5         3 - 4 '/                           &
             ';           6         2 - 4 '/ )
8100  format(3i5,10x,'  ->  ; ',2a8,16x,'  ',a8,i5/f10.5)
8200  format(4i5, 5x,'  ->  ; ',3a8, 8x,'  ',a8,i5/3f10.5)
8300  format(5i5,    '  ->  ; ',4a8,    '  ',a8,i5/6f10.5)

!**********************************************

      return
      end subroutine outshk
