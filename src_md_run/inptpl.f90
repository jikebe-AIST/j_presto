
      subroutine inptpl(iread,iprint,filena,ier)

!*******************************************************************
!
!     INPUT FORMATTED TOPOLOGY FILE
!
!*******************************************************************

      use COMPAR
      use COMBAS ; use COMERG ; use PHYCNS 

      implicit none

      ! Logical unit number for input & output
        integer(4),intent(in):: iread,iprint
      ! File name of topology file
        character(80),intent(in):: filena
      ! Condition code (0: NO ERROR, NEGATIVE: ERROR)
        integer(4),intent(out):: ier

      integer(4),parameter:: numkey = 10
      character(8000):: buffer
      character(80):: line
      character(4):: keywrd(numkey) = (/                               &
                       'TITL','MOLE','ATOM','BOND','ANGL',             &
                       'TORS','IMPR','FUNC','NONB','LINK'/)
      ! Detect or not keyword
      ! If onkywd is .true., read keyword contents (rdkywd)
        logical(1):: onkywd
      ! Request for reading molecule name
      ! If keyword is atom bond angl tors impr
      ! reqmol is .true., read molecule name (chkmol)
        logical(1)::reqmol
      ! Number of times of reading keywords
      ! timeky(3) is times of atom
        integer(4):: timeky(numkey)
      integer(4):: efcol
      character(80):: space = ' '
 
      integer(4):: i,itmp,ienbuf,icol,iselky,molnum
      logical(4):: Flastlink = .false.

!****************************************************
 
!     <<<  INITIAL SETTING  >>>
      write(iprint,*)'INFORMATION> INPTPL (V3.0) '
      write(iprint,*)'       READ FORMATTED TOPOLOGY FILE '
      ier = 0 ; timeky(1:numkey) = 0 ; iselky = 0
      call  initpl

!     Count columns in the topology file
      call system("wc -l "//trim(filena)//" > aho.temp")
      open(unit=1,file="aho.temp",status="old")
      read(1,*)ncTPL
      close(unit=1,status="delete")
      allocate(Tcxmolc(ncTPL),Tixsqml(ncTPL))

!     <<<  OPEN TOPOLOGY FILE  >>>
!         FORMATTED SEQUENTIAL READONLY        : 10
!         FORMATTED SEQUENTIAL READ/WRITE KEEP : 12 FOR CRAY
      call flopen(iread,filena,12,'NULL',0,ier)
      if ( ier .ne. 0 ) then     
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      TOPOLOGY FILE OPEN ERROR'
        ier = -1 ; return
      endif
 
!     <<<  READ TOPOLOGY FILE   >>>
      ! The first line to get max values for arrays
      read(iread,*)line,maxatm,maxbnd,maxang,maxtor,maximp,maxtyp
      call allocate_arrays
      do
        ienbuf = 0 ; buffer = ' '
 
!       1) READ LINE AND STORE BUFFER
        do
          read(iread,'(a80)',end=3000,err=9100)line
          icol = efcol(line,space,';')
          if ( icol .le. 0 ) cycle
          itmp = index(line(1:icol),'->')
          if ( itmp .ne. 0 ) then
            icol = itmp - 1
            call mkbuff(line,icol,8000,ienbuf,buffer,ier)
            if ( ier .ne. 0 ) goto 9200
            cycle
          else
            call mkbuff(line,icol,8000,ienbuf,buffer,ier)
            if ( ier .ne. 0 ) goto 9200
            exit
          endif
        enddo

!       2) HANDLING BUFFER
!        A) MAKE EACH CHAIN INFORMATION IF MOLECULE HAS
!           TWO OR MORE CHAINS          (MKCHNI)
!        B) SEARCH KEYWORD              (SRKYWD)
!        C) READ EACH KEYWORD CONTENTS  (RDKYWD)
        if ( index(buffer(1:ienbuf),'TPL>') .ne. 0 ) then
          if ( Flastlink ) then
              write(iprint,*)'ERROR> INPTPL '
              write(iprint,*)'      IN TOPOLOGY FILE TPL> LINK MUST BE '
              write(iprint,*)'      APPEARED AT THE END OF THE TOPOLOGY'
              write(iprint,*)'      FILE'
              ier = -4 ; return
          endif
          if ( iselky .eq. 2 ) then
            i = ixmolc
            allocate(ixsqml(i),ixatmm(i),ixbndm(i),ixangm(i),ixtorm(i),&
                     iximpm(i),ixtpcn(i),cxmolc(i))
            ixatmm(1:i) = 0 ; ixbndm(1:i) = 0 ; ixangm(1:i) = 0
            ixtorm(1:i) = 0 ; iximpm(1:i) = 0 ; ixtpcn(1:i) = 0
            cxmolc(1:i) = Tcxmolc(1:i) ; ixsqml(1:i) = Tixsqml(1:i)
            itmp = sum(ixsqml(1:i)) ; allocate(ixcend(itmp))
            deallocate(Tcxmolc,Tixsqml)
          endif
          if ( timeky(3).ne.0 .or. timeky(4).ne.0 .or.                 &
               timeky(5).ne.0 .or. timeky(6).ne.0 .or.                 &
               timeky(7).ne.0 ) then
            call mkchni(iselky,molnum,iprint,ier)
            if ( ier .ne. 0 ) return
          endif

          onkywd = .false.
          do i = 1,numkey
            if ( index(buffer(1:ienbuf),keywrd(i)) .ne. 0 ) then
              iselky = i ; onkywd = .true.
              timeky(iselky) = timeky(iselky) + 1
              if ( iselky .eq. 10 ) Flastlink = .true.
            endif
          enddo
          if ( .not. onkywd ) cycle
          if ( iselky.ge.3 .and. iselky.le.7 ) then
            reqmol = .true.
            if ( ixmolc .eq. 0 ) then
              write(iprint,*)'ERROR> INPTPL '
              write(iprint,*)'      IN TOPOLOGY FILE TPL> MOL MUST BE '
              write(iprint,*)'      APPEARED BEFORE TPL> ATOM TPL> BOND'
              write(iprint,*)'      TPL> TORS TPL> IMPR '
              ier = -4 ; return
            endif
          endif
        else
          if ( .not. onkywd ) cycle
          call rdkywd(buffer,ienbuf,iselky,reqmol,molnum,numkey,timeky,&
                      iprint,ier)
          if ( ier .ne. 0 ) return
        endif

      enddo

!     <<<  MAKE NONBONDED PARAMETER  >>>
3000  call mknopm

!     <<<  SEARCH LAST ATOM NUMBER OF EACH RESIDUE   >>>
      call srres(ier)
      if ( ier .ne. 0 ) return

!     <<<  OUTPUT INFORMATION OF INPTPL  >>>
      call flclos(iread,10,ier)
      if ( ier .ne. 0 ) then
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      TOPOLOGY FILE CLOSE ERROR '
        ier = -1 ; return
      endif
      call infinp(numkey,timeky,iprint,ier)

!**************************************

      return

9100  write(iprint,*)'ERROR> INPTPL '
      write(iprint,*)'      READ ERROR '
      write(iprint,*)LINE
      ier = -2 ; return
9200  write(iprint,*)'ERROR> INPTPL '
      write(iprint,*)'      ONE BUFFER MUST BE LESS EQUAL 8000'
      write(iprint,*)'      CHARACTERS '
      write(iprint,*)'      TOO MANY CONTINUOUS LINES '
      ier = -3 ; return
 
      end subroutine inptpl


!===================================================================


      subroutine initpl

!********************************************************
!
!     INITIALIZE COMBAS COMERG PARAMETERS
!
!********************************************************
 
      use COMBAS ; use COMERG 

      implicit none

!****************************************
 
      ixtitl = 0 ; ixmolc = 0 ; ixnatm = 0 ; iynbnd = 0 ; iynang = 0
      iyntor = 0 ; iynimp = 0 ; iynbpf = 0 ; ixnchn = 0 ; iyntyp = -100
      iynpar(1:maxfnc) = 0 ; cynbpf(1:maxfnc) = ' '
      iyppid(1:maxtyp,1:maxtyp) = 0
      fynbpp(1:maxnbp,1:maxtyp,1:maxtyp) = 0.d0
      fyvwrd(1:maxtyp) = 0.d0 ; fyvwme(1:maxtyp) = 0.d0
      fy14sv(1:maxtyp) = 0.d0 ; fy14se(1:maxtyp) = 0.d0

!     CORRECTION START FOR READ INNER COORDINATES
      ixincd(1:maxatm,1:5) = 0
      fxincd(1:maxatm,1:3) = 0.d0

!***************************

      return
      end subroutine initpl


!===================================================================


      subroutine mkbuff(line,numcol,limbuf,ienbuf,buffer,ier)

      implicit none

      character(*),intent(in):: line
      integer(4),intent(in):: numcol,limbuf
      integer(4),intent(inout):: ienbuf
      character(*),intent(inout):: buffer
      integer(4),intent(out):: ier
 
!*********************************************

      ienbuf = ienbuf + numcol
      if ( ienbuf .gt. limbuf ) then
        ier = -1
      else
        buffer((ienbuf-numcol+1):ienbuf) = line(1:numcol)
        ier = 0
      endif

!*****************************
 
      return
      end subroutine mkbuff


!===================================================================


      subroutine mkchni(iselky,curmol,iprint,ier)

!*******************************************************************
!
!     IF TWO OR MORE CHAINS IN ONE MOLECULE , THEN
!     MAKE ATOM BOND ANGL TORS IMPR INFORMATIONS FOR THESE CHAINS
!
!     CALNMM : MAKE IXATMM OR IXBNDM OR IXANGM OR IXTORM OR IXIMPM
!     IXATMM(IMOL) ; NUMBER OF ATOMS OF IMOL-TH MOLECULE
!     MK***C : MAKE CHAIN INFORMATION
!
!*******************************************************************
 
      use COMBAS ; use COMERG 

      implicit none

      integer(4),intent(in):: iselky,curmol,iprint
      integer(4),intent(out):: ier
 
!*****************************************************

      ier = 0
      select case (iselky)
        case (3)
          call calnmm(ixmolc,curmol,ixsqml,ixnatm,ixatmm)
          ixtpcn(curmol) = ixnchn
          ixcend(ixnchn) = ixnatm
          if ( ixsqml(curmol) .gt. 1 ) call mkatmc(curmol,iprint,ier)
        case (4)
          call calnmm(ixmolc,curmol,ixsqml,iynbnd,ixbndm)
          if ( ixsqml(curmol) .gt. 1 ) call mkbndc(curmol,iprint,ier)
        case (5)
          call calnmm(ixmolc,curmol,ixsqml,iynang,ixangm)
          if ( ixsqml(curmol) .gt. 1 ) call mkangc(curmol,iprint,ier)
        case (6)
          call calnmm(ixmolc,curmol,ixsqml,iyntor,ixtorm)
          if ( ixsqml(curmol) .gt. 1 ) call mktorc(curmol,iprint,ier)
        case (7)
          call calnmm(ixmolc,curmol,ixsqml,iynimp,iximpm)
          if ( ixsqml(curmol) .gt. 1 ) call mkimpc(curmol,iprint,ier)
      end select

!****************************

      return
      end subroutine mkchni


!===================================================================


      subroutine rdkywd(buffer,ienbuf,iselky,reqmol,molnum,numkey,     &
                        timeky,iprint,ier)

!***************************************************************
!
!     READ KEYWORD CONTENTS
!      TPL> TITL TPL> MOLE TPL> FUNC TPL> NONB ARE COMMON 
!      FOR EACH MOLECULE
!      TPL> ATOM TPL> BOND TPL> ANGLE TPL> TORS TPL> IMPR
!           ARE PREPARED FOR EACH MOLECULE
!
!**************************************************************
 
      use COMBAS ; use COMERG 

      implicit none

      character(*),intent(in):: buffer
      integer(4),intent(in):: ienbuf,iselky,numkey,iprint
      integer(4),intent(out):: ier
      integer(4),intent(inout):: timeky(numkey),molnum
      logical(1),intent(inout):: reqmol

      integer(4):: imol
 
!***************************************************************

      ier = 0
      ! <<<  READ EACH KEYWORD CONTENTS  >>>
      select case (iselky)
        case (1)
          call rdtitl(buffer,ienbuf,ixtitl,cxtitl,iprint,ier)
        case (2)
          if ( timeky(2) .gt. 2 ) then
            write(iprint,*)"ERROR> INPTPL"
            write(iprint,*)"   IN TOPOLOGY FILE, TPL> MOL MUST BE "
            write(iprint,*)"   APPEARED ONLY ONCE "
            ier = -1 ; return
          endif
          call rdmole(buffer,ienbuf,iprint,ier)
        case (3:7)
          if ( reqmol ) then
            ! <<<  CHECK MOLECULE NAME  >>>
            reqmol = .false.
            call chkmol(buffer,ienbuf,cxmolc,ixmolc,timeky(iselky),    &
                        timeky(3),molnum,iprint,ier)
            if ( ier .ne. 0 ) return
            if ( iselky .eq. 3 ) then
              if ( molnum .ne. timeky(3) ) then
                write(iprint,*)'ERROR> INPTPL '
                write(iprint,*)'    ORDER OF TPL> ATOM MUST BE SAME '
                write(iprint,*)'    WITH ORDER OF MOLECULES IN TPL> MOL'
                do imol = 1 , ixmolc
                  write(iprint,'(10x,a40)')cxmolc(imol)
                enddo
                ier = -1 ; return
              endif
              ixnchn = ixnchn + 1
            endif
          else
            select case (iselky)
              case (3)
                call rdatom(buffer,ienbuf,iprint,ier)
                if ( ier .ne. 0 ) return
                ixachn(ixnatm) = ixnchn
                ixamol(ixnatm) = timeky(3)
              case (4)
                call rdbond(buffer,ienbuf,molnum,iprint,ier)
              case (5)
                call rdangl(buffer,ienbuf,molnum,iprint,ier)
              case (6)
                call rdtors(buffer,ienbuf,molnum,iprint,ier)
              case (7)
                call rdimpr(buffer,ienbuf,molnum,iprint,ier)
            end select
          endif
        case (8)
          call rdfunc(buffer,ienbuf,iprint,ier)
        case (9)
          call rdnonb(buffer,ienbuf,iprint,ier)
        case (10)
          call rdlink(buffer,ienbuf,iprint,ier)
      end select

!*************************************

      return
      end subroutine rdkywd
 
 
!===================================================================


      subroutine mknopm

      use COMBAS ; use COMERG 

      implicit none

      integer(4):: i,j,ityp1,ityp4,itypf
      real(8):: rtmp,rtmp2,rtmp3

!********************************************************
 
!     <<<  CALCUALTE AMBER-LIKE POTENTIAL PARAMETER  >>>
      if ( index(cynbpf(1),'AMBER') .ne. 0 ) then
        do i = 1,iyntyp
        do j = i,iyntyp
          rtmp = sqrt(fyvwme(i)*fyvwme(j))
          rtmp2 = fyvwrd(i)+fyvwrd(j)
          rtmp2 = rtmp2 * rtmp2
          rtmp2 = rtmp2 * rtmp2 * rtmp2
          fynbpp(1,i,j) = rtmp * rtmp2 * 2.d0
          fynbpp(2,i,j) = rtmp * rtmp2 * rtmp2
          fynbpp(1,j,i) = fynbpp(1,i,j)
          fynbpp(2,j,i) = fynbpp(2,i,j)
          if ( iyppid(i,j) .ne. 2 ) then
            iyppid(i,j) = 1 ; iyppid(j,i) = 1
          endif
        enddo
        enddo
      endif
 
!     <<<  CALCUALTE OPLS-LIKE POTENTIAL PARAMETER  >>>
      if ( index(cynbpf(1),'OPLS') .ne. 0 ) then
        do i = 1,iyntyp
        do j = i,iyntyp
          rtmp = sqrt(fyvwme(i)*fyvwme(j)) * 4.d0
          rtmp2 = fyvwrd(i)*fyvwrd(i)*fyvwrd(i)
          rtmp3 = fyvwrd(j)*fyvwrd(j)*fyvwrd(j)
          fynbpp(1,i,j) = rtmp * rtmp2 * rtmp3
          fynbpp(2,i,j) = rtmp * rtmp2*rtmp2 * rtmp3*rtmp3
          fynbpp(1,j,i) = fynbpp(1,i,j)
          fynbpp(2,j,i) = fynbpp(2,i,j)
          iyppid(i,j) = 1 ; iyppid(j,i) = 1
        enddo
        enddo
      endif

!     <<<  CALCUALTE ECEPP-LIKE POTENTIAL PARAMETER  >>>
      if ( index(cynbpf(1),'ECEPP') .ne. 0 ) then
!       <<<  SET SCALE FOR 1-4 INTERACTION  >>>
        do i = 1,iyntor
          ityp1 = ixatyp(iyptor(1,i))
          ityp4 = ixatyp(iyptor(4,i))
          itypf = iyppid(ityp1,ityp4)
          fytvws(i) = fy14sv(itypf) * dble(iytnbf(i))
          fytess(i) = fy14se(itypf) * dble(iytnbf(i))
        enddo
 
!     <<<  SET SCALE FOR 1-4 INTERACTION  FOR AMBER AND OPLS  >>>
      else
        do i = 1,iyntor
          ityp1 = ixatyp(iyptor(1,i))
          ityp4 = ixatyp(iyptor(4,i))
          fytvws(i) = min(fy14sv(ityp1),fy14sv(ityp4))
          fytess(i) = min(fy14se(ityp1),fy14se(ityp4))
          fytvws(i) = fytvws(i) * dble(iytnbf(i))
          fytess(i) = fytess(i) * dble(iytnbf(i))
        enddo
      endif

!***********************************

      return
      end subroutine  mknopm


!===================================================================


      subroutine srres(ier)

      use COMPAR ; use COMBAS

      integer(4),intent(out):: ier

      integer(4):: ichno,ireso,i
      integer(4):: Tixrstr(ncTPL+ixnchn),Tixrend(ncTPL+ixnchn)

!***********************************************

      ier = 0 ; ixnres = 0 ; Tixrstr(1) = 1
      ichno = ixachn(1) ; ireso = ixares(1)
      do i = 1,ixnatm
        if ( ireso.ne.ixares(i) .or. ichno.ne.ixachn(i) ) then
          ireso = ixares(i) ; ichno = ixachn(i)
          ixnres = ixnres + 1
          Tixrend(ixnres) = i - 1 ; Tixrstr(ixnres+1) = i
        endif
      enddo
      ixnres = ixnres + 1
      Tixrend(ixnres) = ixnatm

      allocate(ixrstr(ixnres),ixrend(ixnres))
      ixrstr(1:ixnres) = Tixrstr(1:ixnres)
      ixrend(1:ixnres) = Tixrend(1:ixnres) 

!***********************************

      return
      end subroutine srres


!===================================================================


      subroutine infinp(numkey,timeky,iprint,ier)

      use COMBAS ; use COMERG 

      implicit none

      integer(4),intent(in):: numkey,timeky(numkey),iprint
      integer(4),intent(out):: ier

      integer(4):: i
 
!**********************************************
 
      ier = 0
      write(iprint,*)' '
      write(iprint,*)'     1) TOTAL NUMBER'
      write(iprint,*)'       NUMBER OF MOLECULES         : ',ixmolc
      write(iprint,*)'       NUMBER OF CHAINS            : ',ixnchn
      write(iprint,*)'       NUMBER OF RESIDUES          : ',ixnres
      write(iprint,*)'       NUMBER OF ATOMS             : ',ixnatm
      write(iprint,*)'       NUMBER OF BONDS             : ',iynbnd
      write(iprint,*)'       NUMBER OF ANGLES            : ',iynang
      write(iprint,*)'       NUMBER OF TORSIONS          : ',iyntor
      write(iprint,*)'       NUMBER OF IMPRO.            : ',iynimp
      write(iprint,*)'       MAXIMUM NUMBER OF FUNC.     : ',iynbpf
      write(iprint,*)'       MAXIMUM NUMBER OF ATOM-TYPE : ',iyntyp
      write(iprint,*)' '
      write(iprint,*)'     2) EACH DATA'
      write(iprint,*)'        A) NUMBER OF CHAINS IN EACH MOLECULE '
      if ( ixmolc .ne. 1 ) then
        do i = 1,ixmolc-1
          write(iprint,*)'         MOLECULE ',i,'   : ',            &
                          ixtpcn(i+1) - ixtpcn(i)
        enddo
      endif
      write(iprint,*)'         MOLECULE ',ixmolc,'   : ',              &
                     ixnchn - ixtpcn(ixmolc) + 1
      write(iprint,*)'        B) NUMBER OF ATOMS IN EACH MOLECULE '
      do i = 1,ixmolc
        write(iprint,*)'         MOLECULE ',i,'   : ',ixatmm(i)
      enddo
      write(iprint,*)'        C) NUMBER OF BONDS IN EACH MOLECULE '
      do i = 1,ixmolc
        write(iprint,*)'         MOLECULE ',i,'   : ',ixbndm(i)
      enddo
      write(iprint,*)'        D) NUMBER OF ANGLES IN EACH MOLECULE '
      do i = 1,ixmolc
        write(iprint,*)'         MOLECULE ',i,'   : ',ixangm(i)
      enddo
      write(iprint,*)'        E) NUMBER OF TORSIONS IN EACH MOLECULE '
      do i = 1,ixmolc
        write(iprint,*)'         MOLECULE ',i,'   : ',ixtorm(i)
      enddo
      write(iprint,*)'        F) NUMBER OF IMPROPER IN EACH MOLECULE '
      do i = 1,ixmolc
        write(iprint,*)'         MOLECULE ',i,'   : ',iximpm(i)
      enddo
 
      if ( timeky(8) .eq. 0 ) then
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      THERE IS NO FUNCTION DATA '
        ier = -1 ; return
      endif
      if ( timeky(9) .eq. 0 ) then
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      THERE IS NO NONBONDED DATA '
        ier = -2 ; return
      endif
 
!***************************************

      return
      end subroutine infinp


!===================================================================


      subroutine calnmm(ixmolc,curmol,ixsqml,num,nummm)

      implicit none

      integer(4),intent(in):: ixmolc,curmol,ixsqml(ixmolc),num
      integer(4),intent(inout):: nummm(ixmolc)
      integer(4):: isum,i
 
!*****************************************************

      if ( curmol .eq. 1 ) then
        nummm(1) = num
      else
        isum = 0
        do i = 1,curmol-1
          isum = isum + nummm(i)*ixsqml(i)
        enddo
        nummm(curmol) = num - isum
      endif

!*******************************

      return
      end subroutine calnmm 


!===================================================================


      subroutine mkatmc(curmol,iprint,ier)

      use COMBAS ; use COMERG 

      implicit none

      integer(4),intent(in):: curmol,iprint
      integer(4),intent(out):: ier 

      integer(4):: istatm,ienatm,i,j,k,nlop,iplus
 
!*******************************************************

      ier = 0
      if ( ixnchn .eq. 1 ) then
        istatm = 1 ; ienatm = ixnatm
      else
        istatm = ixcend(ixnchn-1) + 1
        ienatm = ixnatm
      endif
 
      do i = 1,ixsqml(curmol)-1
        ixnchn = ixnchn + 1
        iplus = ixatmm(curmol) * i

        do j = istatm,ienatm
          ixnatm = ixnatm + 1
          if ( ixnatm .gt. maxatm ) then
            write(iprint,*)'ERROR> INPTPL '
            write(iprint,*)'      NUMBER OF ATOMS IS EXCEEDED THE LIMIT'
            write(iprint,*)'      THE LIMIT IS ',maxatm
            ier = -2 ; return
          endif
          ixamol(ixnatm) = ixamol(j)
          ixachn(ixnatm) = ixnchn
          ixares(ixnatm) = ixares(j)
          ixatyp(ixnatm) = ixatyp(j)
          ix14if(ixnatm,1:3) = ix14if(j,1:3)
          nlop = sum(ix14if(j,1:3))
          if ( nlop .gt. 0 ) then
            do k = 1,nlop
              ix14lt(ixnatm,k) = ix14lt(j,k) + iplus
            enddo
          endif

          ! Correction start for read inner coordinates
          do k = 1,4
            if ( ixincd(j,k) .eq. 0 ) then
              ixincd(ixnatm,k) = 0
            else
              ixincd(ixnatm,k) = ixincd(j,k) + iplus
            endif
          enddo
          ixincd(ixnatm,5) = ixincd(j,5)
          fxincd(ixnatm,1:3) = fxincd(j,1:3)

          ! Correction end for read inner coordinates
          fxchrg(ixnatm) = fxchrg(j) ; fxmass(ixnatm) = fxmass(j)
          fxvdwr(ixnatm) = fxvdwr(j) ; cxatmn(ixnatm) = cxatmn(j)
          cxatmt(ixnatm) = cxatmt(j) ; cxresn(ixnatm) = cxresn(j)
        enddo
 
        ixcend(ixnchn) = ixnatm
      enddo

!*************************************************

      return
      end subroutine mkatmc
 
 
!===================================================================


      subroutine mkbndc(curmol,iprint,ier)

      use COMBAS ; use COMERG 

      implicit none

      integer(4),intent(in)::  curmol,iprint
      integer(4),intent(out):: ier

      integer(4):: i,j,ibndo,ibnds,iplus
 
!****************************************

      ier = 0 ; ibndo = iynbnd - ixbndm(curmol)
      do i = 1,ixsqml(curmol)-1
      do j = 1,ixbndm(curmol)
        ibnds = ibndo + j
        iynbnd = iynbnd + 1
        iplus = ixatmm(curmol) * i
        if ( iynbnd .gt. maxbnd ) then
          write(iprint,*)'ERROR> INPTPL '
          write(iprint,*)'      NUMBER OF BONDS IS EXCEEDED THE LIMIT '
          write(iprint,*)'      THE LIMIT IS ',maxbnd
          ier = -1 ; return
        endif
        iypbnd(1:2,iynbnd) = iypbnd(1:2,ibnds) + iplus
        fyfbnd(iynbnd) = fyfbnd(ibnds) ; fyqbnd(iynbnd) = fyqbnd(ibnds)
      enddo
      enddo

!***********************************

      return
      end subroutine mkbndc
 
 
!===================================================================


      subroutine mkangc(curmol,iprint,ier)

      use COMBAS ; use COMERG 

      implicit none

      integer(4),intent(in):: curmol,iprint
      integer(4),intent(out):: ier

      integer(4):: i,j,iango,iangs,iplus
 
!****************************************

      ier = 0 ; iango = iynang - ixangm(curmol)
      do i = 1,ixsqml(curmol)-1
      do j = 1,ixangm(curmol)
        iynang = iynang + 1
        iangs = iango + j
        iplus = ixatmm(curmol) * i
        if ( iynang .gt. maxang ) then
          write(iprint,*)'ERROR> INPTPL '
          write(iprint,*)'      NUMBER OF ANGLES IS EXCEEDED THE LIMIT'
          write(iprint,*)'      THE LIMIT IS ',maxang
          ier = -1 ; return
        endif
        iypang(1:3,iynang) = iypang(1:3,iangs) + iplus
        fyfang(iynang) = fyfang(iangs) ; fyqang(iynang) = fyqang(iangs)
      enddo
      enddo

!************************************

      return
      end subroutine mkangc
 
 
!===================================================================


      subroutine mktorc(curmol,iprint,ier)

      use COMBAS ; use COMERG 

      implicit none

      integer(4),intent(in):: curmol,iprint
      integer(4),intent(out):: ier

      integer(4):: i,j,itoro,itors,iplus
 
!*******************************************

      ier = 0 ; itoro = iyntor - ixtorm(curmol)
      do i = 1,ixsqml(curmol)-1
      do j = 1,ixtorm(curmol)
        iyntor = iyntor + 1
        itors = itoro + j
        iplus = ixatmm(curmol) * i
        if ( iyntor .gt. maxtor ) then
          write(iprint,*)'ERROR> INPTPL '
          write(iprint,*)'     NUMBER OF TORSIONS IS EXCEEDED THE LIMIT'
          write(iprint,*)'     THE LIMIT IS ',maxtor
          ier = -1 ; return
        endif
        iyptor(1:4,iyntor) = iyptor(1:4,itors) + iplus
        iytdiv(iyntor) = iytdiv(itors) ; fyftor(iyntor) = fyftor(itors)
        fytrot(iyntor) = fytrot(itors) ; fytphs(iyntor) = fytphs(itors)
        iytnbf(iyntor) = iytnbf(itors)
      enddo
      enddo

!*****************************

      return
      end subroutine mktorc
 
 
!===================================================================


      subroutine mkimpc(curmol,iprint,ier)

      use COMBAS ; use COMERG 

      implicit none

      integer(4),intent(in):: curmol,iprint
      integer(4),intent(out):: ier

      integer(4):: i,j,iimpo,iimps,iplus
 
!********************************************

      ier = 0 ; iimpo = iynimp - iximpm(curmol)
      do i = 1,ixsqml(curmol)-1
      do j = 1,iximpm(curmol)
        iynimp = iynimp + 1
        iimps = iimpo + j
        iplus = ixatmm(curmol) * i
        if ( iynimp .gt. maximp ) then
          write(iprint,*)'ERROR> INPTPL '
          write(iprint,*)'      NUMBER OF IMPROPER-TORSIONS IS '
          write(iprint,*)'      EXCEEDED THE LIMIT '
          write(iprint,*)'      THE LIMIT IS ',maximp
          IER = -1 ; return
        endif
        iypimp(1:4,iynimp) = iypimp(1:4,iimps) + iplus
        iyidiv(iynimp) = iyidiv(iimps) ; fyfimp(iynimp) = fyfimp(iimps)
        fyirot(iynimp) = fyirot(iimps) ; fyiphs(iynimp) = fyiphs(iimps)
        iyinbf(iynimp) = iyinbf(iimps)
      enddo
      enddo

!**************************************

      return
      end subroutine mkimpc
 
 
!===================================================================


      subroutine chkmol(buffer,ienbuf,cxmolc,ixmolc,curtim,curmol,     &
                        molnum,iprint,ier)

      implicit none

      character(*),intent(in):: buffer
      integer(4),intent(in):: ienbuf,ixmolc,curtim,iprint
      character(*),intent(in):: cxmolc(ixmolc)
      integer(4),intent(inout):: curmol
      integer(4),intent(out):: molnum,ier

      integer(4):: tpcol,lacol,jtcol,jlcol,i,ktcol,klcol,ichk
 
!*********************************************

      ier = 0 ; molnum = 0 ; ichk = 0
      jtcol = tpcol(buffer,ienbuf) ; jlcol = lacol(buffer,ienbuf)
 
!     * SERACH MOLECULE-NAME *
      do i = 1,ixmolc
        ktcol = tpcol(cxmolc(i),40) ; klcol = lacol(cxmolc(i),40)
        if ( cxmolc(i)(ktcol:klcol) .eq. buffer(jtcol:jlcol) ) then
          molnum = i ; ichk = 1
          if ( molnum .lt. curtim ) then
            write(iprint,*)'ERROR> INPTPL '
            write(iprint,*)'      DATA OF THIS MOLECULE IS ALREADY READ'
            write(iprint,*)'      MOLECULE IS ',buffer(jtcol:jlcol)
            write(iprint,*)' '
            ier = -1 ; return
          endif
          if ( molnum .gt. curmol ) then
            write(iprint,*)'ERROR> INPTPL '
            write(iprint,*)'      NO-ATOM DATA IN THIS MOLECULE '
            write(iprint,*)'      MOLECULE IS ',buffer(jtcol:jlcol)
            write(iprint,*)'      CHECK ORDER OF TPL> ATOM DATA '
            write(iprint,*)'      TPL> ATOM DATA MUST BE THE TOP OF '
            write(iprint,*)'      EACH MOLECULE '
            write(iprint,*)' '
            ier = -2 ; return
          endif
          exit
        endif
      enddo

      if ( ichk .eq. 0 ) then
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      THIS MOLECULE ',buffer(jtcol:jlcol),' IS '
        write(iprint,*)'      UNKNOWN '
        write(iprint,*)' '
        ier = -3 ; return
      endif

!*****************************************

      return
      end subroutine chkmol


!=====================================================


      function tpcol(buffer,ienbuf)

      implicit none

      character(*),intent(in):: buffer
      integer(4),intent(in):: ienbuf
      integer(4):: tpcol,i

!*************************************

      do i = 1,ienbuf
        if ( buffer(i:i) .ne. ' ' ) then
          tpcol = i ; return
        endif
      enddo
      tpcol = ienbuf

!**************************

      return
      end function tpcol


!=========================================================


      function lacol(buffer,ienbuf) 

      implicit none

      character(*),intent(in):: buffer
      integer(4),intent(in):: ienbuf
      integer(4):: lacol,i

!***********************************

      do i = 1,ienbuf
        if ( buffer(i:i) .eq. ' ' ) then
          lacol = i - 1 ; return
        endif
      enddo
      lacol = ienbuf

!******************************

      return
      end function lacol


!===================================================================


      subroutine rdtitl(buffer,ienbuf,ixtitl,cxtitl,iprint,ier)

      implicit none

      character(*),intent(in):: buffer
      character(*),intent(out):: cxtitl(10)
      integer(4),intent(in):: ienbuf,iprint
      integer(4),intent(inout):: ixtitl,ier
 
!*******************************************

      ier = 0
      if ( ienbuf .le. 80 ) then
        ixtitl = ixtitl + 1
        if ( ixtitl .gt. 10 ) then
          write(iprint,*)'ERROR> INPTPL '
          write(iprint,*)'      NUMBER OF TITLE-LINES MUST BE LESS '
          write(iprint,*)'      EQUAL 10 '
          ier = -1 ; return
        endif
        cxtitl(ixtitl) = ' '
        cxtitl(ixtitl)(1:ienbuf) = buffer(1:ienbuf)
      else
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      CHARCATER LENGTH OF ONE TITLE LINE '
        write(iprint,*)'      IS LESS EQUAL 80 '
        write(iprint,*)'      CURRENT LENGTH IS ',ienbuf
        write(iprint,*)'      DO NOT USE -> IN TITLE LINE '
        ier = -2 ; return
      endif

!***********************************

      return
      end subroutine rdtitl
 
 
!===================================================================


      subroutine rdmole(buffer,ienbuf,iprint,ier)

      use COMBAS 

      implicit none

      character(*),intent(in):: buffer
      integer(4),intent(in):: ienbuf,iprint
      integer(4),intent(out):: ier

      character(1):: inpmod(2)
      character(80):: cval(2)
      integer(4):: ival(2),iwork1(2),iwork2(2),ierrd
      real(8):: rval(2)

!******************************************

      ier = 0 ; ixmolc = ixmolc + 1
      inpmod(1) = 'C' ; inpmod(2) = 'I'
      call rdfree(buffer,ienbuf,2,2,inpmod,iwork1,iwork2,rval,ival,    &
                  cval,ierrd)
      if ( ierrd .lt. 0 ) then
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      READ ERROR IN TPL> MOLE '
        ier = -2 ; return
      endif
      Tcxmolc(ixmolc) = cval(1)(1:40)
      if ( ival(1) .lt. 0 ) then
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      NUMBER OF CHAINS IS STRANGE '
        write(iprint,*)'      MOLECULE : ',cxmolc(ixmolc)
        write(iprint,*)'      NUMBER OF CHAINS IS ',ival(1)
        ier = -3 ; return
      endif
      Tixsqml(ixmolc) = ival(1)
 
!***************************************

      return
      end subroutine rdmole
 
 
!===================================================================


      subroutine rdatom(buffer,ienbuf,iprint,ier)

      use COMBAS 

      implicit none

      character(*),intent(in):: buffer
      integer(4),intent(in):: ienbuf,iprint
      integer(4),intent(out):: ier

      integer(4),parameter:: maxval = 60
      character(1):: inpmod(maxval)
      integer(4):: iwork1(maxval),iwork2(maxval),ival(maxval)
      real(8):: rval(maxval)
      character(80):: cval(maxval)
      integer(4):: ierrd,numval,numred,numopt,j
 
!*****************************************************

      ier = 0 ; ixnatm = ixnatm + 1
      if ( ixnatm .gt. maxatm ) then
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      NUMBER OF ATOMS IS EXCEEDED THE LIMIT '
        write(iprint,*)'      LIMIT NUMBER IS ',maxatm
        ier = -1 ; return
      endif
      inpmod(1:2) = 'C' ; inpmod(3) = 'I'   ; inpmod(4) = 'C'
      inpmod(5) = 'I'   ; inpmod(6:8) = 'R' ; inpmod(9:11) = 'I'
      call rdfree(buffer,ienbuf,maxval,11,inpmod,iwork1,iwork2,rval,   &
                  ival,cval,ierrd)
      if ( ierrd .lt. 0 ) then
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      READ ERROR IN TPL> ATOM '
        ier = -2 ; return
      endif

      cxatmn(ixnatm) = cval(1)(1:8) ; cxatmt(ixnatm) = cval(2)(1:4)
      ixatyp(ixnatm) = ival(1)      ; cxresn(ixnatm) = cval(3)(1:8)
      ixares(ixnatm) = ival(2)      ; fxmass(ixnatm) = rval(1)
      fxvdwr(ixnatm) = rval(2)      ; fxchrg(ixnatm) = rval(3)
      ix14if(ixnatm,1:3) = ival(3:5)

      numval = sum(ix14if(ixnatm,1:3))
      if ( numval .ne. 0 ) inpmod(12:11+numval) = 'I'
      numred = numval + 11
      call rdfree(buffer,ienbuf,maxval,numred,inpmod,iwork1,iwork2,    &
                  rval,ival,cval,ierrd)
      if ( ierrd .lt. 0 ) then
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      READ ERROR IN TPL> ATOM '
        ier = -2 ; return
      endif
      if ( numval .ne. 0 )                                             &
        ix14lt(ixnatm,1:numval) = ival(6:numval+5) + ixnatm
      
      !  CORRECTION START FOR READ INNER COORDINATES
      numopt = numval + 11
      inpmod(numopt+1:numopt+4) = 'I'
      inpmod(numopt+5:numopt+7) = 'R'
      numred = numopt + 7

      call rdfree(buffer,ienbuf,maxval,numred,inpmod,iwork1,iwork2,    &
                  rval,ival,cval,ierrd)
      if ( ierrd .lt. 0 ) then
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      READ ERROR IN TPL> ATOM '
        ier = -2 ; return
      endif
      numopt = numval + 5
      do j = 1,4
        if ( ival(numopt+j) .eq. 0 ) then
          ixincd(ixnatm,j) = 0
        else
          ixincd(ixnatm,j) = ival(numopt+j) + ixnatm
        endif
      enddo

      fxincd(ixnatm,1:3) = rval(4:6)
      !  CORRECTION END FOR READ INNER COORDINATES

!*********************************************

      return
      end subroutine rdatom 
 
 
!===================================================================


      subroutine rdbond(buffer,ienbuf,molnum,iprint,ier)

      use COMBAS ; use COMERG 

      implicit none

      character(*),intent(in):: buffer
      integer(4),intent(in):: ienbuf,molnum,iprint
      integer(4),intent(out):: ier

      integer(4):: lastat,iwrk1(10),iwrk2(10),ival(10)
      character(1):: inpmod(10)
      real(8):: rval(10)
      character(80):: cval(10)
 
!**************************************************

      ier = 0 ; iynbnd = iynbnd + 1
      if ( iynbnd .gt. maxbnd ) then
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      NUMBER OF BONDS IS EXCEEDED THE LIMIT '
        write(iprint,*)'      LIMIT NUMBER IS ',maxbnd
        ier = -1 ; return
      endif
      if ( ixtpcn(molnum) .eq. 1 ) then
        lastat = 0
      else
        lastat = ixcend(ixtpcn(molnum)-1)
      endif

      inpmod(1:2) = 'I' ; inpmod(3:4) = 'R'
      call rdfree(buffer,ienbuf,10,4,inpmod,iwrk1,iwrk2,rval,ival,cval,&
                  ier)
      if ( ier .lt. 0 ) then
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      READ ERROR IN TPL> BOND '
        ier = -2 ; return
      endif

      iypbnd(1:2,iynbnd) = ival(1:2)
      fyfbnd(iynbnd) = rval(1) ; fyqbnd(iynbnd) = rval(2)
      iypbnd(1:2,iynbnd) = iypbnd(1:2,iynbnd) + lastat

!****************************************

      return
      end subroutine rdbond 
 
 
!===================================================================


      subroutine rdangl(buffer,ienbuf,molnum,iprint,ier)

      use COMBAS ; use COMERG ; use PHYCNS 

      implicit none

      character(*),intent(in):: buffer
      integer(4),intent(in):: ienbuf,iprint
      integer(4),intent(out):: ier
      integer(4),intent(inout):: molnum

      integer(4):: lastat,iwrk1(10),iwrk2(10),ival(10)
      character(1):: inpmod(10)
      real(8):: rval(10)
      character(80):: cval(10)
 
!*********************************************

      ier = 0 ; iynang = iynang + 1
      if ( iynang .gt. maxang ) then
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      NUMBER OF ANGLES IS EXCEEDED THE LIMIT '
        write(iprint,*)'      THE LIMIT IS ',maxang
        ier = -1 ; return
      endif
      if ( ixtpcn(molnum) .eq. 1 ) then
        lastat = 0
      else
        lastat = ixcend(ixtpcn(molnum)-1)
      endif

      inpmod(1:3) = 'I' ; inpmod(4:5) = 'R'
      call rdfree(buffer,ienbuf,10,5,inpmod,iwrk1,iwrk2,rval,ival,cval,&
                  ier)
      if ( ier .lt. 0 ) then
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      READ ERROR IN TPL> ANGL  '
        ier = -2 ; return
      endif

      iypang(1:3,iynang) = ival(1:3)
      fyfang(iynang) = rval(1) ; fyqang(iynang) = rval(2)
      iypang(1:3,iynang) = iypang(1:3,iynang) + lastat
      fyqang(iynang) = fyqang(iynang) * rad

!****************************************

      return
      end subroutine rdangl

 
!===================================================================


      subroutine rdtors(buffer,ienbuf,molnum,iprint,ier)
 
      use COMBAS ; use COMERG ; use PHYCNS 

      implicit none

      character(*),intent(in):: buffer
      integer(4),intent(in):: ienbuf,iprint
      integer(4),intent(out):: ier
      integer(4),intent(inout):: molnum

      integer(4):: lastat,iwrk1(10),iwrk2(10),ival(10)
      character(1):: inpmod(10)
      real(8):: rval(10)
      character(80):: cval(10)
 
!***********************************************

      ier = 0 ; iyntor = iyntor + 1
      if ( iyntor .gt. maxtor ) then
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      NUMBER OF TORSIONS IS EXCEEDED THE LIMIT '
        write(iprint,*)'      LIMIT NUMBER IS ',maxtor
        ier = -1 ; return
      endif
      if ( ixtpcn(molnum) .eq. 1 ) then
        lastat = 0
      else
        lastat = ixcend(ixtpcn(molnum)-1)
      endif

      inpmod(1:4) = 'I' ; inpmod(5) = 'R' ; inpmod(6) = 'I'
      inpmod(7:8) = 'R' ; inpmod(9) = 'I'
      call rdfree(buffer,ienbuf,10,9,inpmod,iwrk1,iwrk2,rval,ival,cval,&
                  ier)
      if ( ier .lt. 0 ) then
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      READ ERROR IN TPL> TORS '
        ier = -2 ; return
      endif

      iyptor(1:4,iyntor) = ival(1:4)
      fyftor(iyntor) = rval(1) ; iytdiv(iyntor) = ival(5)
      fytrot(iyntor) = rval(2) ; fytphs(iyntor) = rval(3)
      iytnbf(iyntor) = ival(6)
      iyptor(1:4,iyntor) = iyptor(1:4,iyntor) + lastat
      fytphs(iyntor) = fytphs(iyntor) * rad

!*******************************************

      return
      end subroutine rdtors

 
!===================================================================


      subroutine rdimpr(buffer,ienbuf,molnum,iprint,ier)

      use COMBAS ; use COMERG ; use PHYCNS 

      implicit none

      character(*),intent(in):: buffer
      integer(4),intent(in):: ienbuf,iprint
      integer(4),intent(out):: ier
      integer(4),intent(inout):: molnum

      integer(4):: lastat,iwrk1(10),iwrk2(10),ival(10)
      character(1):: inpmod(10)
      real(8):: rval(10)
      character(80):: cval(10)
 
!********************************************

      ier = 0 ; iynimp = iynimp + 1
      if ( iynimp .gt. maximp ) then
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      NUMBER OF IMPROPER IS EXCEEDED THE LIMIT '
        write(iprint,*)'      LIMIT NUMBER IS ',maximp
        ier = -1 ; return
      endif
      if ( ixtpcn(molnum) .eq. 1 ) then
        lastat = 0
      else
        lastat = ixcend(ixtpcn(molnum)-1)
      endif

      inpmod(1:4) = 'I' ; inpmod(5) = 'R' ; inpmod(6) = 'I'
      inpmod(7:8) = 'R' ; inpmod(9) = 'I'
      call rdfree(buffer,ienbuf,10,9,inpmod,iwrk1,iwrk2,rval,ival,cval,&
                  ier)
      if ( ier .lt. 0 ) then
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      READ ERROR IN TPL> IMPRO'
        ier = -2 ; return
      endif

      iypimp(1:4,iynimp) = ival(1:4)
      fyfimp(iynimp) = rval(1) ; iyidiv(iynimp) = ival(5)
      fyirot(iynimp) = rval(2) ; fyiphs(iynimp) = rval(3)
      iyinbf(iynimp) = ival(6)
      iypimp(1:4,iynimp) = iypimp(1:4,iynimp) + lastat
      fyiphs(iynimp) = fyiphs(iynimp) * rad

!******************************************

      return
      end subroutine rdimpr


!===================================================================

      
      subroutine rdfunc(buffer,ienbuf,iprint,ier)

      use COMBAS ; use COMERG 

      implicit none

      character(*),intent(in):: buffer
      integer(4),intent(in):: ienbuf,iprint
      integer(4),intent(out):: ier

      integer(4):: ierrd,ival(3),iwork1(3),iwork2(3)
      character(1):: inpmod(3)
      character(80):: cval(3)
      real(8):: rval(3)
 
!*******************************************

      ier = 0 ; inpmod(1:2) = 'I' ; inpmod(3) = 'C'
      call rdfree(buffer,ienbuf,3,3,inpmod,iwork1,iwork2,rval,ival,    &
                  cval,ierrd)
      if ( ierrd .lt. 0 ) then
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      READ ERROR IN TPL > FUNC '
        ier = -1 ; return
      endif
      if ( ival(1).gt.maxfnc .or. ival(1).le.0 ) then
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      FUNCTION-NUMBER IS OUT OF RANGE '
        write(iprint,*)'      CURRENT IS    ',ival(1)
        write(iprint,*)'      MAX NUMBER IS ',maxfnc
        IER = -2 ; return
      endif
      if ( ival(2).gt.maxnbp .or. ival(2).le.0 ) then
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      NUMBER OF PARAMETERS FOR NONBONDED'
        write(iprint,*)'      INTERACTION IS EXCEEDED THE LIMIT '
        write(iprint,*)'      CURRENT IS   ',ival(2)
        write(iprint,*)'      THE LIMIT IS ',maxnbp
        ier = -3 ; return
      endif

      iynpar(ival(1)) = ival(2)
      cynbpf(ival(1)) = cval(1)(1:40)
      if ( ival(1) .ge. iynbpf ) iynbpf = ival(1)
 
      if ( index(cynbpf(ival(1)),'AMBER').eq.0 .and.                   &
           index(cynbpf(ival(1)),'ECEPP').eq.0 .and.                   &
           index(cynbpf(ival(1)),'OPLS').eq.0 ) then
        write(iprint,*)'ERROR> RDFUNC '
        write(iprint,*)'      AMBER-LIKE POTENTIAL'
        write(iprint,*)'      ECEPP-LIKE POTENTIAL,'
        write(iprint,*)'      OR OPLS-LIKE  POTENTIAL ARE AVILABLE '
        write(iprint,*)' '
        write(iprint,*)'      THERE IS AMBER, ECEPP OR OPLS '
        write(iprint,*)'      IN FUNCTION NAME'
        ier = -4 ; return
      endif
 
!*******************************************

      return
      end subroutine rdfunc
 
 
!===================================================================


      subroutine rdnonb(buffer,ienbuf,iprint,ier)

      use COMBAS ; use COMERG 

      implicit none

      character(*),intent(in):: buffer
      integer(4),intent(in):: ienbuf,iprint
      integer(4),intent(out):: ier

      integer(4):: i1,i2,i3,iwrk1(10),iwrk2(10),ival(10)
      character(1):: inpmod(10)
      character(80):: cval(10)
      real(8):: rval(10)
 
!*************************************************

      ier = 0 ; inpmod(1:3) = 'I'
      call rdfree(buffer,ienbuf,10,3,inpmod,iwrk1,iwrk2,rval,ival,cval,&
                  ier)
      if ( ier .lt. 0 ) then
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      READ ERROR IN TPL> NONB '
        ier = -1 ; return
      endif
      i1 = ival(1) ; i2 = ival(2) ; i3 = ival(3)
      if ( i1.gt.maxtyp .or. i1.le.0 .or. i2.gt.maxtyp ) then
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      ATOM-TYPE NUMBER IS OUT OF RANGE '
        write(iprint,*)'      ATOM-TYPE NUMBER IS ',i1,i2
        write(iprint,*)'      MAX NUMBER IS ',maxtyp
        ier = -2 ; return
      endif
      if ( i3.gt.maxfnc .or. i3.le.0 ) then
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      FUNCTION NUMBER IS OUT OF RANGE '
        write(iprint,*)'      CURRENT IS    ',i3
        write(iprint,*)'      MAX NUMBER IS ',maxfnc
        ier = -3 ; return
      endif
      if ( iynpar(i3) .eq. 0 ) then
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      FUNCTION-NUMBER IS UNKNWON '
        ier = -4 ; return
      endif

!     <<<  READ AMBER-LINK POTENTIAL PARAMETER  >>
!            VDW PARAMETER (FUNCTION-1)    RADIUS + ENERGY(MIN)
!            HYDROGEN-BONDED (FUNCTION-2)  C + D
      if ( index(cynbpf(1),'AMBER') .ne. 0 ) then
        call rdnoam(buffer,ienbuf,i3,maxtyp,maxnbp,iyntyp,iyppid,      &
                    fynbpp,fyvwrd,fyvwme,fy14sv,fy14se,iprint,ier)
        if ( ier .ne. 0 ) return
      endif

!     <<<  READ ECEPP-LINK POTENTIAL PARAMETER  >>
!            VDW PARAMETER (FUNCTION-1)    RADIUS**2 & ENERGY(MIN)
!            HYDROGEN-BONDED (FUNCTION-2)  RADISU**2 & ENERGY(MIN)
      if ( index(cynbpf(1),'ECEPP') .ne. 0 ) then
        call rdnoec(buffer,ienbuf,i3,maxtyp,maxnbp,iyntyp,iyppid,      &
                    fynbpp,fy14sv,fy14se,iprint,ier)
        if ( ier .ne. 0 ) return
      endif

!     <<<  READ OPLS-LINK POTENTIAL PARAMETER  >>>
!            VDW PARAMETER (FUNCTION-1)    RADIUS + ENERGY(MIN)
      if ( index(cynbpf(1),'OPLS') .ne. 0 ) then
        call rdnoop(buffer,ienbuf,i3,maxtyp,iyntyp,fyvwrd,fyvwme,      &
                    fy14sv,fy14se,iprint,ier)
        if ( ier .ne. 0 ) return
      endif

!***************************************

      return
      end subroutine rdnonb
 
 
!===================================================================


      subroutine rdnoam(buffer,ienbuf,funtyp,maxtyp,maxnbp,iyntyp,     &
                        iyppid,fynbpp,fyvwrd,fyvwme,fy14sv,fy14se,     &
                        iprint,ier)

      implicit none

      character(*),intent(in):: buffer
      integer(4),intent(in):: ienbuf,funtyp,maxtyp,maxnbp,iprint
      integer(4),intent(inout):: iyntyp,iyppid(maxtyp,maxtyp)
      integer(4),intent(out):: ier
      real(8),intent(inout):: fynbpp(maxnbp,maxtyp,maxtyp),            &
                            fyvwrd(maxtyp),fyvwme(maxtyp),             &
                            fy14sv(maxtyp),fy14se(maxtyp)

      integer(4):: i1,i2,iwrk1(10),iwrk2(10),ival(10)
      character(1):: inpmod(10)
      character(80):: cval(10)
      real(8):: rval(10)

!******************************************************

      ier = 0
      if ( funtyp.ne.1 .and. funtyp.ne.2 ) then
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      UNKNOWN FUNCTION TYPE '
        write(iprint,*)'      AMBER-TYPE FUNCTION TYPE NUMBER IS'
        write(iprint,*)'      1 OR 2 '
        ier = -1 ; return
      endif

      if ( funtyp .eq. 1 ) then
        inpmod(1:3) = 'I' ; inpmod(4:7) = 'R'
        call rdfree(buffer,ienbuf,10,7,inpmod,iwrk1,iwrk2,rval,ival,   &
                    cval,ier)
        if ( ier .lt. 0 ) then
          write(iprint,*)'ERROR> INPTPL '
          write(iprint,*)'      READ ERROR IN TPL> NONB '
          ier = -2 ; return
        endif
        i1 = ival(1)
        fyvwrd(i1) = rval(1) ; fyvwme(i1) = rval(2)
        fy14se(i1) = rval(3) ; fy14sv(i1) = rval(4)
        iyntyp = max(iyntyp,i1)
      elseif ( funtyp .eq. 2 ) then
        inpmod(1:3) = 'I' ; inpmod(4:5) = 'R'
        call rdfree(buffer,ienbuf,10,5,inpmod,iwrk1,iwrk2,rval,ival,   &
                    cval,ier)
        if ( ier .lt. 0 ) then
          write(iprint,*)'ERROR> INPTPL '
          write(iprint,*)'      READ ERROR IN TPL> NONB '
          ier = -2 ; return
        endif
        i1 = ival(1) ; i2 = ival(2)
        iyntyp = max(iyntyp,i1) ; iyntyp = max(iyntyp,i2)
        iyppid(i1,i2) = ival(3)
        fynbpp(4,i1,i2) = rval(1) ; fynbpp(3,i1,i2) = rval(2)
        iyppid(i2,i1) = ival(3)
        fynbpp(4,i2,i1) = rval(1) ; fynbpp(3,i2,i1) = rval(2)
      endif
 
!***************************************

      return
      end subroutine rdnoam
 

!===================================================================


      subroutine rdnoec(buffer,ienbuf,funtyp,maxtyp,maxnbp,iyntyp,     &
                        iyppid,fynbpp,fy14sv,fy14se,iprint,ier)

      implicit none

      character(*),intent(in):: buffer
      integer(4),intent(in):: ienbuf,funtyp,maxtyp,maxnbp,iprint
      integer(4),intent(inout):: iyntyp,iyppid(maxtyp,maxtyp)
      integer(4),intent(out):: ier
      real(8),intent(inout):: fynbpp(maxnbp,maxtyp,maxtyp),            &
                            fy14sv(maxtyp),fy14se(maxtyp)

      integer(4):: i1,i2,i3,iwrk1(10),iwrk2(10),ival(10)
      real(8):: r1,r2,r3,r4,rval(10),rtmp
      character(1):: inpmod(10)
      character(80):: cval(10)
 
!***********************************************

      ier = 0
      if ( funtyp.ne.1 .and. funtyp.ne.2 ) then
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      UNKNOWN FUNCTION TYPE '
        write(iprint,*)'      ECEPP-TYPE FUNCTION TYPE NUMBER IS'
        write(iprint,*)'      1 OR 2 '
        ier = -1 ; return
      endif
 
!     * READ VDW PARAMETER *
      if ( funtyp .eq. 1 ) then
        inpmod(1:3) = 'I' ; inpmod(4:7) = 'R'
        call rdfree(buffer,ienbuf,10,7,inpmod,iwrk1,iwrk2,rval,ival,   &
                    cval,ier)
        if ( ier .lt. 0 ) then
          write(iprint,*)'ERROR> INPTPL '
          write(iprint,*)'      READ ERROR IN TPL> NONB '
          ier = -2 ; return
        endif
        i1 = ival(1) ; i2 = ival(2) ; i3 = ival(3)
        r1 = rval(1) ; r2 = rval(2) ; r3 = rval(3) ; r4 = rval(4)
        iyppid(i1,i2) = i3 ; iyppid(i2,i1) = i3
        rtmp = r1*r1*r1
        fynbpp(1,i1,i2) = r2 * rtmp * 2.d0
        fynbpp(2,i1,i2) = r2 * rtmp * rtmp
        fynbpp(1:2,i2,i1) = fynbpp(1:2,i1,i2)
        fy14sv(i3) = r3 ; fy14se(i3) = r4
        iyntyp = max(i1,iyntyp) ; iyntyp = max(iyntyp,i2)
 
!     * READ HYDROGEN-BONDED PARAMETER *
      elseif ( funtyp .eq. 2 ) then
        inpmod(1:3) = 'I' ; inpmod(4:7) = 'R'
        call rdfree(buffer,ienbuf,10,7,inpmod,iwrk1,iwrk2,rval,ival,   &
                    cval,ier)
        if ( ier .lt. 0 ) then
          write(iprint,*)'ERROR> INPTPL '
          write(iprint,*)'      READ ERROR IN TPL> NONB '
          ier = -2 ; return
        endif
        i1 = ival(1) ; i2 = ival(2) ; i3 = ival(3)
        r1 = rval(1) ; r2 = rval(2) ; r3 = rval(3) ; r4 = rval(4)
        iyppid(i1,i2) = i3 ; iyppid(i2,i1) = i3
        rtmp = r1*r1*r1
        fynbpp(1,i1,i2) = r2 * rtmp * 2.d0
        fynbpp(3,i1,i2) = r2 * rtmp * r1 * r1 * 2.d0
        fynbpp(4,i1,i2) = r2 * rtmp * rtmp
        fynbpp(3:4,i2,i1) = fynbpp(3:4,i1,i2)
        fy14sv(i3) = r3 ; fy14se(i3) = r4
        iyntyp = max(i1,iyntyp) ; iyntyp = max(iyntyp,i2)
      endif
 
!**************************************

      return
      end subroutine rdnoec

 
!===================================================================


      subroutine rdnoop(buffer,ienbuf,funtyp,maxtyp,iyntyp,fyvwrd,     &
                        fyvwme,fy14sv,fy14se,iprint,ier)

      implicit none

      character(*),intent(in):: buffer
      integer(4),intent(in):: ienbuf,funtyp,maxtyp,iprint
      integer(4),intent(out):: ier
      integer(4),intent(inout):: iyntyp
      real(8),intent(inout):: fyvwrd(maxtyp),fyvwme(maxtyp),           &
                            fy14sv(maxtyp),fy14se(maxtyp)

      integer(4):: i1,iwrk1(10),iwrk2(10),ival(10)
      real(8):: rval(10)
      character(1):: inpmod(10)
      character(80):: cval(10)

!**********************************************************
 
      ier = 0
      if ( funtyp .ne. 1 ) then
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      UNKNOWN FUNCTION TYPE '
        write(iprint,*)'      OPLS-TYPE FUNCTION TYPE NUMBER MUST BE 1 '
        ier = -1 ; return
      endif
      inpmod(1:3) = 'I' ; inpmod(4:7) = 'R'
      call rdfree(buffer,ienbuf,10,7,inpmod,iwrk1,iwrk2,rval,ival,cval,&
                  ier)
      if ( ier .lt. 0 ) then
        write(iprint,*)'ERROR> INPTPL '
        write(iprint,*)'      READ ERROR IN TPL> NONB '
        ier = -2 ; return
      endif

      i1 = ival(1)
      fyvwrd(i1) = rval(1) ; fyvwme(i1) = rval(2)
      fy14se(i1) = rval(3) ; fy14sv(i1) = rval(4)
      iyntyp = max(iyntyp,i1)

!*************************************

      return
      end subroutine rdnoop


!===================================================================


      subroutine rdlink(buffer,ienbuf,iprint,ier)

      use COMBAS ; use COMERG  ; use PHYCNS

      implicit none

      character(*),intent(in):: buffer
      integer(4),intent(in):: ienbuf,iprint
      integer(4),intent(out):: ier

      integer(4):: ierrd,iwork1(10),iwork2(10),ival(10),i,isum,        &
        iat1,iat2
      character(1):: inpmod(10)
      character(80):: cval(10)
      real(8):: rval(10)

!**********************************************************

      ier = 0
      if ( buffer(1:1) .eq. "B" ) then
        inpmod(1) = 'C' ; inpmod(2:3) = 'I' ; inpmod(4:5) = 'R'
        call rdfree(buffer,ienbuf,5,5,inpmod,iwork1,iwork2,rval,ival,  &
                cval,ier)
        iynbnd = iynbnd + 1 ; i = iynbnd
        if ( iynbnd .gt. maxbnd ) then
          write(iprint,*)'ERROR> INPTPL '
          write(iprint,*)'      NUMBER OF BONDS IS EXCEEDED THE LIMIT '
          write(iprint,*)'      THE LIMIT IS ',maxbnd
          ier = -1 ; return
        endif
        iypbnd(1:2,i) = ival(1:2)
        fyfbnd(i) = rval(1) ; fyqbnd(i) = rval(2)

        iat1 = minval(ival(1:2)) ; iat2 = maxval(ival(1:2))
        isum = sum(ix14if(iat1,1:3))
        if ( isum+1 .gt. max14n ) then
          write(iprint,*)'ERROR> INPTPL '
          write(iprint,*)'     IX14IF IS EXCEEDED THE LIMIT (MAX14N)'
          write(iprint,*)'     THE LIMIT IS ',max14n
        endif
        do i = isum,ix14if(iat1,1)+1,-1
          ix14lt(iat1,i+1) = ix14lt(iat1,i)
        enddo
        i = ix14if(iat1,1) + 1
        ix14lt(iat1,i) = iat2
        ix14if(iat1,1) = i

      elseif ( buffer(1:1) .eq. "A" ) then
        inpmod(1) = 'C' ; inpmod(2:4) = 'I' ; inpmod(5:6) = 'R'
        call rdfree(buffer,ienbuf,6,6,inpmod,iwork1,iwork2,rval,ival,  &
                cval,ier)
        iynang = iynang + 1 ; i = iynang
        if ( iynang .gt. maxang ) then
          write(iprint,*)'ERROR> INPTPL '
          write(iprint,*)'      NUMBER OF ANGLES IS EXCEEDED THE LIMIT'
          write(iprint,*)'      THE LIMIT IS ',maxang
          ier = -1 ; return
        endif
        iypang(1:3,i) = ival(1:3)
        fyfang(i) = rval(1) ; fyqang(i) = rval(2) * rad

        iat1 = min(ival(1),ival(3)) ; iat2 = max(ival(1),ival(3))
        isum = sum(ix14if(iat1,1:3))
        if ( isum+1 .gt. max14n ) then
          write(iprint,*)'ERROR> INPTPL '
          write(iprint,*)'     IX14IF IS EXCEEDED THE LIMIT (MAX14N)'
          write(iprint,*)'     THE LIMIT IS ',max14n
        endif
        do i = isum,sum(ix14if(iat1,1:2))+1,-1
          ix14lt(iat1,i+1) = ix14lt(iat1,i)
        enddo
        i = ix14if(iat1,2) + 1
        ix14lt(iat1,ix14if(iat1,1)+i) = iat2
        ix14if(iat1,2) = i

      elseif ( buffer(1:1) .eq. "T" ) then
        inpmod(1) = 'C' ; inpmod(2:5) = 'I' ; inpmod(6) = 'R'
        inpmod(7) = 'I' ; inpmod(8:9) = 'R' ; inpmod(10) = 'I'
        call rdfree(buffer,ienbuf,10,10,inpmod,iwork1,iwork2,rval,ival,&
                cval,ier)
        iyntor = iyntor + 1 ; i = iyntor
        if ( iyntor .gt. maxtor ) then
          write(iprint,*)'ERROR> INPTPL '
          write(iprint,*)'     NUMBER OF TORSIONS IS EXCEEDED THE LIMIT'
          write(iprint,*)'     THE LIMIT IS ',maxtor
          ier = -1 ; return
        endif
        iyptor(1:4,i) = ival(1:4)
        iytdiv(i) = ival(5) ; fyftor(i) = rval(1)
        fytrot(i) = rval(2) ; fytphs(i) = rval(3) * rad
        iytnbf(i) = ival(6)

        if ( iytnbf(i) .eq. 1 ) then
          iat1 = min(ival(1),ival(4)) ; iat2 = max(ival(1),ival(4))
          isum = sum(ix14if(iat1,1:3))
          if ( isum+1 .gt. max14n ) then
            write(iprint,*)'ERROR> INPTPL '
            write(iprint,*)'     IX14IF IS EXCEEDED THE LIMIT (MAX14N)'
            write(iprint,*)'     THE LIMIT IS ',max14n
          endif
          i = ix14if(iat1,3) + 1
          ix14lt(iat1,isum+1) = iat2
          ix14if(iat1,3) = i
        endif

      elseif ( buffer(1:1) .eq. "I" ) then
        inpmod(1) = 'C' ; inpmod(2:5) = 'I' ; inpmod(6) = 'R'
        inpmod(7) = 'I' ; inpmod(8:9) = 'R' ; inpmod(10) = 'I'
        call rdfree(buffer,ienbuf,10,10,inpmod,iwork1,iwork2,rval,ival,&
                cval,ier)
        iynimp = iynimp + 1 ; i = iynimp
        if ( iynimp .gt. maximp ) then
          write(iprint,*)'ERROR> INPTPL '
          write(iprint,*)'      NUMBER OF IMPROPER-TORSIONS IS '
          write(iprint,*)'      EXCEEDED THE LIMIT '
          write(iprint,*)'      THE LIMIT IS ',maximp
          ier = -1 ; return
        endif
        iypimp(1:4,i) = ival(1:4)
        iyidiv(i) = ival(5) ; fyfimp(i) = rval(1)
        fyirot(i) = rval(2) ; fyiphs(i) = rval(3) * rad
        iyinbf(i) = ival(6)
      endif

      if ( ierrd .lt. 0 ) then
        write(iprint,*)'ERROR> INPTPL'
        write(iprint,*)'      READ ERROR IN TPL> LINK '
        ier = -2 ; return
      endif

!************************************************************

      return
      end subroutine rdlink
