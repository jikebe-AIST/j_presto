
      subroutine inpbou(iread,iprint,filena,ier) 

!*******************************************************************
!
!     THIS SUBROUTINE IS FOR READ CONTROL DATA FOR
!       SETTING BOUNDARY CONDITION
!
!*******************************************************************
 
      use COMBAS ; use COMMIS 
 
      implicit none

      ! Logical unit number for read and output log
        integer(4),intent(in):: iread,iprint
      ! File name of control data of CAP
        character(80),intent(in):: filena
      ! Condition code (0: NO ERROR)
        integer(4),intent(out):: ier

      integer(4):: efcol
      character(80):: line,space
      logical(4):: onincl,oncent,onbox,onradi,onoutl,onatom,onbuff,    &
                   onshap,onell,onforc,onforp,onfunc
      integer(4):: iersub,icol,ichn,istchn,j
      
!********************************************************************
 
!     0) INITIAL SETTING
      write(iprint,*)'INFORMATION> INPBOU '
      ier = 0 ; natcap = 0 ; space = " "
      onincl = .false. ; oncent = .false. ; onbox = .false.
      onradi = .false. ; onbuff = .false. ; onoutl = .false.
      onatom = .false. ; onshap = .false. ; onell = .false.
      onforc = .false. ; onforp = .false. ; onfunc = .false.
      allocate(Tixbfcn(ixnchn))
      Tixbfcn(1:ixnchn) = .false.
 
!     1) OPEN CONTROL DATA FOR SETTING BOUNDARY CONDITION
      call flopen(iread,filena,10,'NULL',0,iersub)
      if ( iersub .ne. 0 ) then
        write(iprint,*)'ERROR> INPBOU '
        write(iprint,*)'       FILE OPEN ERROR '
        write(iprint,*)' '
        ier = iersub ; return
      endif
 
!     2) READ CONTROL DATA FOR SETTING VARIABLES
      do
        read(iread,'(a80)',end=200)line
        icol = efcol(line,space,";")
        if ( icol .le. 0 ) cycle
        if ( index(line(1:icol),"BOUND>") .ne. 0 ) then
          onincl = .false. ; oncent = .false. ; onbox = .false.
          onradi = .false. ; onbuff = .false. ; onatom = .false.
          onshap = .false. ; onell  = .false. ; onforc = .false.
          onforp = .false. ; onfunc = .false.
          if ( index(line,"INCL") .ne. 0 ) onincl = .true.
          if ( index(line,"CENT") .ne. 0 ) oncent = .true.
          if ( index(line,"BOX ") .ne. 0 ) onbox  = .true.
          if ( index(line,"RADI") .ne. 0 ) onradi = .true.
          if ( index(line,"BUFF") .ne. 0 ) onbuff = .true.
          if ( index(line,"LIST") .ne. 0 ) onatom = .true.
          if ( index(line,"SHAP") .ne. 0 ) onshap = .true.
          if ( index(line,"ELL ") .ne. 0 ) onell  = .true.
          if ( index(line,"FORC") .ne. 0 ) onforc  = .true.
          if ( index(line,"FORP") .ne. 0 ) onforp  = .true.
          if ( index(line,"FUNC") .ne. 0 ) onfunc  = .true.
          cycle
        endif

!       2-1) READ INCLUDE DATA
        if ( onincl ) then
          call rdbinc(ier,iprint,line,icol,onoutl)
!       2-2) READ CENTER DATA
        elseif ( oncent ) then
          call rdbcen(ier,iprint,line,icol)
!       2-3) READ BOX DATA
        elseif ( onbox ) then
          call rdbbox(ier,iprint,line,icol)
!       2-4) READ RADIUS DATA
        elseif ( onradi ) then
          call rdbrad(ier,iprint,line,icol)
!       2-5) READ RADIUS DATA FOR PROTEIN
        elseif ( onbuff ) then
          call rdbbuff(ier,iprint,line,icol)
!       2-6) READ ATOM INCLUDE DATA
        elseif ( onatom ) then
          call rdbatm(ier,iprint,line,icol,0)
!       2-7) READ CAP SHAPE
        elseif ( onshap ) then
          call rdbshp(ier,iprint,line,icol)
!       2-8) READ ELLipsoid radius
        elseif ( onell ) then
          call rdbell(ier,iprint,line,icol)
!       2-9) READ FORCe constant
        elseif ( onforc ) then
          call rdbforc(ier,iprint,line,icol)
!       2-10) READ FORCe constant for Protein
        elseif ( onforp ) then
          call rdbforp(ier,iprint,line,icol)
!       2-11) READ FUNCtion type for boundary
        elseif ( onfunc ) then
          call rdbfunc(ier,iprint,line,icol)
        endif

        if ( ier .ne. 0 ) return
      enddo

!     3) COUNT NUMBER OF VARIABLES
200   call rdbatm(ier,iprint,line,icol,1)
      nchncap = count(Tixbfcn(1:ixnchn))
      if ( nchncap .eq. 0 ) then
        write(iprint,*)'ERROR> INPBOU '
        write(iprint,*)'    NO-CHAIN IS SELECTED FOR BOUNDARY'
        write(iprint,*)'    CONDITION '
        write(iprint,*)' '
        ier = -1832 ; return
      endif

      ! Make ixbfcn
      allocate(ixbfcn(nchncap))
      j = 0
      do ichn = 1,ixnchn
        if ( Tixbfcn(ichn) ) then
          j = j + 1
          ixbfcn(j) = ichn
        endif
      enddo
      deallocate(Tixbfcn)

      istchn = ixbfcn(1)
      if ( istchn .eq. 1 ) then
        ixtatb = 1
      else
        ixtatb = ixcend(istchn-1) + 1
      endif
      write(iprint,*)'             INCLUDE '
      write(iprint,*)'             NUMBER OF CHAINS   : ',nchncap
      write(iprint,*)'             TOP ATOM NUMBER    : ',ixtatb
      if ( onoutl ) then
        do ichn = 1,nchncap
          write(iprint,'(i10,5x,a40)')                                 &
            ixbfcn(ichn),cxmolc(ixamol(ixcend(ixbfcn(ichn))))
        enddo
      endif
 
!     4) CLOSE CONTROL DATA FOR SETTING BOUNDARY CONDITION
      call flclos(iread,10,iersub)
      if ( iersub .ne. 0 ) then
        write(iprint,*)'ERROR> INPBOU '
        write(iprint,*)'       FILE CLOSE ERROR '
        write(iprint,*)' '
        ier = iersub ; return
      endif

!*******************************

      return
      end subroutine inpbou
 
 
!===========================================================================


      subroutine rdbinc(ier,iprint,line,icol,onoutl)

!*******************************************************************
!
!     THIS SUBROUTINE IS FOR READ INCLUDE DATA
!
!*******************************************************************
 
      use COMBAS 
 
      implicit none

      ! Condition code
        integer(4),intent(out):: ier
      ! Logical unit number for output log
        integer(4),intent(in):: iprint
      ! Input line
        character(*),intent(in):: line
      ! Effective column number
        integer(4),intent(in):: icol
      ! Flag for output list
        logical(1),intent(inout):: onoutl

      ! Work for rdfree
        integer(4),parameter:: maxval = 4
        integer(4):: ival(maxval),iwork1(maxval),iwork2(maxval)
        real(8):: rval(maxval)
        character(80):: cval(maxval) 
        character(1):: inpmod(maxval)

      integer(4):: tpcol,lacol,iersub,istchn,ienchn,jstcol,jencol,imol,&
                   islmol,jstchn,jenchn
 
!*******************************************************************

      ier = 0 ; inpmod(1) = "C" ; inpmod(2:3) = "I" ; inpmod(4) = "C"
      call rdfree(line,icol,maxval,maxval,inpmod,iwork1,iwork2,rval,   &
                  ival,cval,iersub)
      if ( iersub .lt. 0 ) then
        write(iprint,*)'ERROR> INPBOU '
        write(iprint,*)'       DATA TYPE ERROR '
        write(iprint,'(a80)')line
        write(iprint,*)' '
        ier = iersub ; return 
      endif
      if ( index(cval(2),'YES') .ne. 0 ) onoutl = .true.
      istchn = ival(1) ; ienchn = ival(2)
      jstcol = tpcol(cval(1),80) ; jencol = lacol(cval(1),80)
      do imol = 1,ixmolc
        if ( index(cxmolc(imol),cval(1)(jstcol:jencol)) .ne. 0 ) then
          islmol = imol ; exit
        endif
        if ( imol .eq. ixmolc ) then
          write(iprint,*)'ERROR> INPBOU '
          write(iprint,*)'      UNKNOWN MOLECULE NAME '
          write(iprint,*)'    ',cval(1)(jstcol:jencol)
          write(iprint,*)' '
          ier = -1830 ; return
        endif
      enddo
      jstchn = sum(ixsqml(1:islmol-1))

      istchn = max(istchn,1) ; ienchn = min(ienchn,ixsqml(islmol))
      jenchn = jstchn + ienchn ; jstchn = jstchn + istchn

      Tixbfcn(jstchn:jenchn) = .true.

!**************************************

      return
      end subroutine rdbinc
 

!====================================================================

      subroutine rdbcen(ier,iprint,line,icol)

!*******************************************************************
!
!     THIS SUBROUTINE IS FOR READ CENTER DATA
!
!*******************************************************************
 
      use COMBAS ; use COMCMMC
 
      implicit none

      ! Condition code
        integer(4),intent(out):: ier
      ! Logical unit number for output log
        integer(4),intent(in):: iprint
      ! Input line
        character(*),intent(in):: line
      ! Effective column number
        integer(4),intent(in):: icol

      integer(4),parameter:: maxval = 4
      integer(4):: ival(maxval),iwork1(maxval),iwork2(maxval)
      real(8):: rval(maxval)
      character(80):: cval(maxval)
      character(1):: inpmod(maxval)

      real(8):: molwet
      integer(4):: iersub,idchin,istatm,ienatm,iatm,iatsel
 
!*******************************************************************

      ier = 0
!     * READ CHAIN DATA *
      if ( index(line(1:icol),'CHAI') .ne. 0 ) then
        inpmod(1) = "C" ; inpmod(2) = "I"
        call rdfree(line,icol,maxval,2,inpmod,iwork1,iwork2,rval,ival, &
                    cval,iersub)
        if ( iersub .lt. 0 ) then
          write(iprint,*)'ERROR> INPBOU '
          write(iprint,*)'       DATA TYPE ERROR '
          write(iprint,'(a80)')line
          write(iprint,*) ' '
          ier = iersub ; return
        endif
        idchin = ival(1)
        if ( idchin .eq. 1 ) then
          istatm = 1
        else
          istatm = ixcend(idchin-1) + 1
        endif
        ienatm = ixcend(idchin)
        molwet = sum(fxmass(istatm:ienatm))
        molwet = 1.d0 / molwet
        fxcbou(1:3) = 0.d0
        do iatm = istatm,ienatm
          fxcbou(1:3) = fxcbou(1:3) + (fxmass(iatm)*cord(1:3,iatm))
        enddo
        fxcbou(1:3) = fxcbou(1:3) * molwet
 
!     * READ ATOM DATA *
      elseif ( index(line(1:icol),"ATOM") .ne. 0 ) then
        inpmod(1) = "C" ; inpmod(2) = "I"
        call rdfree(line,icol,maxval,2,inpmod,iwork1,iwork2,rval,ival, &
                    cval,iersub)
        if ( iersub .lt. 0 ) then
          write(iprint,*)'ERROR> INPBOU '
          write(iprint,*)'       DATA TYPE ERROR '
          write(iprint,*)line
          write(iprint,*)' '
          ier = iersub ; return
        endif
        iatsel = ival(1)
        if ( iatsel .gt. ixnatm ) then
          write(iprint,*)'ERROR> INPBOU '
          write(iprint,*)'   UNKNOWN ATOM NUMBER '
          write(iprint,*)' '
          ier = -1831 ; return
        endif
        fxcbou(1:3) = cord(1:3,iatsel)
 
!     * READ COORDINATE DATA *
      elseif ( index(line(1:icol),"COOR") .ne. 0 ) then
        inpmod(1) = "C" ; inpmod(2:4) = "R"
        call rdfree(line,icol,maxval,4,inpmod,iwork1,iwork2,rval,ival, &
                    cval,iersub)
        if ( iersub .lt. 0 ) then
          write(iprint,*)'ERROR> INPBOU '
          write(iprint,*)'       DATA TYPE ERROR '
          write(iprint,'(a80)')line
          write(iprint,*)' '
          ier = iersub ; return
        endif
        fxcbou(1:3) = rval(1:3)
      endif

!*****************************************

      return
      end subroutine rdbcen


!====================================================================


      subroutine rdbbox(ier,iprint,line,icol)

!*******************************************************************
!
!     THIS SUBROUTINE IS FOR READ BOX DATA
!
!*******************************************************************
 
      use COMBAS 
 
      implicit none

      ! Condition code
        integer(4),intent(out):: ier
      ! Logical unit numebr for output log
        integer(4),intent(in):: iprint
      ! Input line
        character(*),intent(in):: line
      ! Effective column number
        integer(4),intent(in):: icol

      integer(4),parameter:: maxval = 7
      integer(4):: ival(maxval),iwork1(maxval),iwork2(maxval)
      real(8):: rval(maxval)
      character(80):: cval(maxval)
      character(1):: inpmod(maxval)
      integer(4):: iersub
 
!*******************************************************************

      ier = 0 ; inpmod(1) = "C" ; inpmod(2:7) = "R"
      call rdfree(line,icol,maxval,maxval,inpmod,iwork1,iwork2,rval,   &
                  ival,cval,iersub)
      if ( iersub .lt. 0 ) then
        write(iprint,*)'ERROR> INPBOU '
        write(iprint,*)'       DATA TYPE ERROR '
        write(iprint,'(a80)')line
        write(iprint,*) ' '
        ier = iersub ; return
      endif

      ! Cell size
      if ( cval(1)(1:4) .eq. "SIZE" ) then
        fxcell(1:3) = rval(1:3)
      ! Cell boundary coordinates
      elseif ( cval(1)(1:4) .eq. "COOR" ) then
        celwal(1:6) = rval(1:6)
      else
        write(iprint,*)"ERROR> INPBOU"
        write(iprint,*)"IN BOUND> BOX, YOU MUST SELECT EITHER OPTION"
        write(iprint,*)"1 : SIZE (ASSIGN 3(X,Y,Z) CELL SIZE)"
        write(iprint,*)"2 : COORDINATE (ASSIGN 6(XMIN,XMAX,YMIN,YMAX"//&
                       ",ZMIN,ZMAX) COORDINATES)"
        ier = -1 ; return
      endif

!***************************************

      return
      end subroutine rdbbox


!=========================================================================


      subroutine rdbrad(ier,iprint,line,icol)

!*******************************************************************
!
!     THIS SUBROUTINE IS FOR READ RADIUS DATA
!
!*******************************************************************
 
      use COMBAS ; use COMMIS 
 
      implicit none

      ! Condition code
        integer(4),intent(out):: ier
      ! Logical unit number for output log
        integer(4),intent(in):: iprint
      ! Input line
        character(*),intent(in):: line
      ! Effective column number
        integer(4),intent(in):: icol

      integer(4),parameter:: maxval = 1
      integer(4):: ival(maxval),iwork1(maxval),iwork2(maxval)
      real(8):: rval(maxval)
      character(80):: cval(maxval)
      character(1):: inpmod(maxval)
      integer(4):: iersub
 
!*******************************************************************

      ier = 0 ; inpmod(1) = "R"
      call rdfree(line,icol,maxval,maxval,inpmod,iwork1,iwork2,rval,   &
                  ival,cval,iersub)
      if ( iersub .lt. 0 ) then
        write(iprint,*)'ERROR> INPBOU '
        write(iprint,*)'       DATA TYPE ERROR '
        write(iprint,'(a80)')line
        write(iprint,*)' '
        ier = iersub ; return
      endif
      furcap = rval(1)

!******************************

      return
      end subroutine rdbrad


!====================================================================


      subroutine rdbbuff(ier,iprint,line,icol)

!*******************************************************************
!
!     THIS SUBROUTINE IS FOR READ BUFFer length for solute-solvent areas
!
!*******************************************************************

      use COMBAS ; use COMMIS

      implicit none

      ! Condition code
        integer(4),intent(out):: ier
      ! Logical unit number for output log
        integer(4),intent(in):: iprint
      ! Input line
        character(*),intent(in):: line
      ! Effective column number
        integer(4),intent(in):: icol

      integer(4),parameter:: maxval = 1
      integer(4):: ival(maxval),iwork1(maxval),iwork2(maxval)
      real(8):: rval(maxval)
      character(80):: cval(maxval)
      character(1):: inpmod(maxval)
      integer(4):: iersub

!*******************************************************************

      ier = 0 ; inpmod(1) = "R"
      call rdfree(line,icol,maxval,maxval,inpmod,iwork1,iwork2,rval,   &
                  ival,cval,iersub)
      if ( iersub .lt. 0 ) then
        write(iprint,*)'ERROR> INPBOU '
        write(iprint,*)'       DATA TYPE ERROR '
        write(iprint,'(a80)')LINE
        write(iprint,*)' '
        ier = iersub ; return
      endif
      CAPbuff = rval(1)

!**************************************

      return
      end subroutine rdbbuff


!=======================================================================


      subroutine rdbatm(ier,iprint,line,icol,final_flg)

!************************************************************
!
!      ROUTINE TO SELECT ATOM FOR CAP 
!      AND SEARCH ATOM POINTER( ABSOLUTE ATOM NUMBER )
!
!************************************************************

      use COMBAS ; use COMMIS ; use COMCMM

      implicit none

      integer(4),intent(out):: ier
      integer(4),intent(in):: iprint
      character(*),intent(in):: line
      integer(4),intent(in):: icol
      integer(4),intent(in):: final_flg

      integer(4),parameter:: maxval = 7
      integer(4):: ival(maxval),iwork1(maxval),iwork2(maxval)
      real(8):: rval(maxval)
      character(80):: cval(maxval)
      character(1):: inpmod(maxval)
      integer(4):: ichain,jchain,nrsfst,nrsend,i,j,k,iersub
      character(4):: atmnam,resnam,output,cc
      logical(4),allocatable,save:: bflg(:),yflg(:)

!********************************************************

       ier = 0
       if ( final_flg .eq. 0 ) then
         if ( .not. allocated(bflg) ) then
           allocate(bflg(ixnatm)) ; bflg(:) = .false.
           allocate(yflg(ixnatm)) ; yflg(:) = .false.
         endif
         inpmod(1:4) = "I" ; inpmod(5:7) = "C"
         iwork1(1:maxval) = 1 ; iwork2(1:maxval) = 80

         call rdfree(line,icol,maxval,maxval,inpmod,iwork1,iwork2,rval,&
                     ival,cval,iersub)

!        READ ERROR IN RDFREE
         if ( iersub .lt. 0 ) then
           write(iprint,'(x,a13)')"ERROR> INPBOU"
           write(iprint,'(8x,a20)')"READ ERROR IN RDFREE"
           ier = iersub ; return
         endif

         ichain = ival(1) ; jchain = ival(2)
         nrsfst = ival(3) ; nrsend = ival(4)
         atmnam = cval(1) ; resnam = cval(2) ; output = cval(3)

!        SEARCH ATOM POINTER FOR EACH GROUP
         j = index(resnam,"*") ; if ( j.eq.0 ) j = 5
         k = index(atmnam,"*") ; if ( k.eq.0 ) k = 5
         if ( ichain .eq. 0 ) then
           do i = 1,ixnatm
             if ( absres(i).lt.nrsfst .or. absres(i).gt.nrsend ) cycle
             if ( j.ne.1 .and. resnam(1:j-1).ne.cxresn(i)(1:j-1) ) cycle
             cc = cxatmn(i)
             if ( atmnam.eq."HEAV" .or. atmnam.eq."heav" ) then
               if ( cc(1:1).eq."H" .or. cc(1:1).eq."h" ) cycle
             elseif ( atmnam.eq."SIDE" .or. atmnam.eq."side" ) then
               if ( cc.eq."C   " .or. cc.eq."CA  " .or.                &
                    cc.eq."H   " .or. cc.eq."O   " .or.                &
                    cc.eq."N   " .or. cc.eq."OX  " .or.                &
                    cc.eq."O5' " .or. cc.eq."C5' " .or.                &
                    cc.eq."C4' " .or. cc.eq."C3' " .or.                &
                    cc.eq."O3' " .or. cc.eq."P   " .or.                &
                    cc.eq."H5' " .or. cc.eq."H5''" .or.                &
                    cc.eq."H4' " .or. cc.eq."H3' " .or.                &
                    cc.eq."OP1 " .or. cc.eq."OP2 " ) cycle
             elseif ( atmnam.eq."RING" .or. atmnam.eq."ring" ) then
               if ( cc.eq."C   " .or. cc.eq."CA  " .or.                &
                    cc.eq."H   " .or. cc.eq."O   " .or.                &
                    cc.eq."N   " .or. cc.eq."OX  " .or.                &
                    cc.eq."CB  " ) cycle
             elseif ( ( k.ne.1 .and. atmnam(1:k-1).ne.cc(1:k-1)) ) then
               cycle
             endif
             bflg(i) = .true.
             if ( output(1:1) .eq. "Y" ) yflg(i) = .true.
           enddo
         else
           do i = 1,ixnatm
             if ( ixachn(i).lt.ichain .or. ixachn(i).gt.jchain .or.    &
                  ixares(i).lt.nrsfst .or. ixares(i).gt.nrsend ) cycle
             if ( j.ne.1 .and. resnam(1:j-1).ne.cxresn(i)(1:j-1) ) cycle
             cc = cxatmn(i)
             if ( atmnam.eq."HEAV" .or. atmnam.eq."heav" ) then
               if ( cc(1:1).eq."H" .or. cc(1:1).eq."h" ) cycle
             elseif ( atmnam.eq."SIDE" .or. atmnam.eq."side" ) then
               if ( cc.eq."C   " .or. cc.eq."CA  " .or.                &
                    cc.eq."H   " .or. cc.eq."O   " .or.                &
                    cc.eq."N   " .or. cc.eq."OX  " .or.                &
                    cc.eq."O5' " .or. cc.eq."C5' " .or.                &
                    cc.eq."C4' " .or. cc.eq."C3' " .or.                &
                    cc.eq."O3' " .or. cc.eq."P   " .or.                &
                    cc.eq."H5' " .or. cc.eq."H5''" .or.                &
                    cc.eq."H4' " .or. cc.eq."H3' " .or.                &
                    cc.eq."OP1 " .or. cc.eq."OP2 " ) cycle
             elseif ( atmnam.eq."RING" .or. atmnam.eq."ring" ) then
               if ( cc.eq."C   " .or. cc.eq."CA  " .or.                &
                    cc.eq."H   " .or. cc.eq."O   " .or.                &
                    cc.eq."N   " .or. cc.eq."OX  " .or.                &
                    cc.eq."CB  " ) cycle
             elseif ( ( k.ne.1 .and. atmnam(1:k-1).ne.cc(1:k-1)) ) then
               cycle
             endif
             bflg(i) = .true.
             if ( output(1:1) .eq. "Y" ) yflg(i) = .true.
           enddo
         endif
         write(iprint,'(4i6,3x,2a4)')ichain,jchain,nrsfst,nrsend,      &
                                     atmnam,resnam
       else

!        ALLOCATE ATOM GROUP
         if ( count(bflg) .gt. 0 ) then
           write(iprint,*)                                           &
             "        SELECT FOLLOWING ATOMS FOR CONSTRAINT,"
           do i = 1,ixnatm
             if ( bflg(i) ) then
               natcap = natcap + 1 ; iamcap(natcap) = i
               if (yflg(i)) write(iprint,'(2i6,2x,a4,x,a4,2i6)') &
                 ixachn(i),ixares(i),cxresn(i),cxatmn(i),natcap,i
             endif
           enddo
           write(iprint,*)""
           write(iprint,'(8x,a40,i6)')                                 &
             "TOTAL NUMBER OF ATOMS IN ALL GROUP ARE :",natcap
           write(iprint,*)""
         endif
         deallocate(bflg,yflg)
       endif

!************************************************

       return
       end subroutine rdbatm


!========================================================================


      subroutine rdbshp(ier,iprint,line,icol)

!*******************************************************************
!
!     THIS SUBROUTINE IS FOR READ CAP SHAPE INFO.
!
!*******************************************************************
 
      use COMBAS ; use COMMIS 
 
      implicit none

      ! Condition code
        integer(4),intent(out):: ier
      ! Logical unit number for output log
        integer(4),intent(in):: iprint
      ! Input line
        character(*),intent(in):: line
      ! Effective column number
        integer(4),intent(in):: icol

      integer(4),parameter:: maxval = 1
      integer(4):: ival(maxval),iwork1(maxval),iwork2(maxval)
      real(8):: rval(maxval)
      character(80):: cval(maxval)
      character(1):: inpmod(maxval)
      integer(4):: iersub
 
!*******************************************************************

      ier = 0 ; inpmod(1) = "C"
      call rdfree(line,icol,maxval,maxval,inpmod,iwork1,iwork2,rval,   &
                  ival,cval,iersub)
      if ( iersub .lt. 0 ) then
        write(iprint,*)'ERROR> INPBOU '
        write(iprint,*)'       DATA TYPE ERROR '
        write(iprint,'(a80)')line
        write(iprint,*)' '
        ier = iersub ; return
      endif
      if ( cval(1)(1:4) .eq. "SPHE" ) then
        capshp = 1
      elseif ( cval(1)(1:3) .eq. "BOX" ) then
        capshp = 2
      elseif ( cval(1)(1:3) .eq. "ELL" ) then
        capshp = 3
      else
        write(iprint,*)"ERROR> INPBOU"
        write(iprint,*)"PLEASE INPUT CORRECT SHAPE OPTION"
        write(iprint,*)" ('SPHE','BOX ' or 'ELL ')"
        ier = -1 ; return
      endif

!******************************

      return
      end subroutine rdbshp


!=========================================================================


      subroutine rdbell(ier,iprint,line,icol)

!**********************************************************
!
!     THIS SUBROUTINE IS FOR READ ELLIPSOIDAL CAP RADIUS
!
!**********************************************************

      use COMBAS ; use COMMIS

      implicit none

      ! Condition code
        integer(4),intent(out):: ier
      ! Logical unit number for output log
        integer(4),intent(in):: iprint
      ! Input line
        character(*),intent(in):: line
      ! Effective column number
        integer(4),intent(in):: icol

      integer(4),parameter:: maxval = 4
      integer(4):: ival(maxval),iwork1(maxval),iwork2(maxval)
      real(8):: rval(maxval)
      character(80):: cval(maxval)
      character(1):: inpmod(maxval)
      integer(4):: iersub

!***********************************************

      ier = 0 ; inpmod(1) = "C" ; inpmod(2:4) = "R"
      call rdfree(line,icol,maxval,maxval,inpmod,iwork1,iwork2,rval,   &
                  ival,cval,iersub)
      if ( iersub .lt. 0 ) then
        write(iprint,*)'ERROR> INPBOU '
        write(iprint,*)'       DATA TYPE ERROR '
        write(iprint,'(a80)')line
        write(iprint,*) ' '
        ier = iersub ; return
      endif

      ! Cell size
      if ( cval(1)(1:4) .eq. "SIZE" ) then
        fxellp(1:3) = rval(1:3)
      else
        write(iprint,*)"ERROR> INPBOU"
        write(iprint,*)"IN BOUND> ELL, YOU MUST SELECT ONLY A OPTION"
        write(iprint,*)"1 : SIZE (ASSIGN 3(X,Y,Z) RADIUS SIZE)"
        ier = -1 ; return
      endif

!**************************************************

      return
      end subroutine rdbell


!=========================================================================


      subroutine rdbforc(ier,iprint,line,icol)

!*******************************************************************
!
!     THIS SUBROUTINE IS FOR READ CAP FORCE
!
!*******************************************************************
 
      use COMBAS ; use COMMIS 
 
      implicit none

      ! Condition code
        integer(4),intent(out):: ier
      ! Logical unit number for output log
        integer(4),intent(in):: iprint
      ! Input line
        character(*),intent(in):: line
      ! Effective column number
        integer(4),intent(in):: icol

      integer(4),parameter:: maxval = 1
      integer(4):: ival(maxval),iwork1(maxval),iwork2(maxval)
      real(8):: rval(maxval)
      character(80):: cval(maxval)
      character(1):: inpmod(maxval)
      integer(4):: iersub
 
!*******************************************************************

      ier = 0 ; inpmod(1) = "R"
      call rdfree(line,icol,maxval,maxval,inpmod,iwork1,iwork2,rval,   &
                  ival,cval,iersub)
      if ( iersub .lt. 0 ) then
        write(iprint,*)'ERROR> INPBOU '
        write(iprint,*)'       DATA TYPE ERROR '
        write(iprint,'(a80)')line
        write(iprint,*)' '
        ier = iersub ; return
      endif
      fukcap = rval(1)

!******************************

      return
      end subroutine rdbforc


!====================================================================


      subroutine rdbforp(ier,iprint,line,icol)

!*******************************************************************
!
!     THIS SUBROUTINE IS FOR READ CAP FORCE FOR PROTEIN
!
!*******************************************************************
 
      use COMBAS ; use COMMIS 
 
      implicit none

      ! Condition code
        integer(4),intent(out):: ier
      ! Logical unit number for output log
        integer(4),intent(in):: iprint
      ! Input line
        character(*),intent(in):: line
      ! Effective column number
        integer(4),intent(in):: icol

      integer(4),parameter:: maxval = 1
      integer(4):: ival(maxval),iwork1(maxval),iwork2(maxval)
      real(8):: rval(maxval)
      character(80):: cval(maxval)
      character(1):: inpmod(maxval)
      integer(4):: iersub
 
!*******************************************************************

      ier = 0 ; inpmod(1) = "R"
      call rdfree(line,icol,maxval,maxval,inpmod,iwork1,iwork2,rval,   &
                  ival,cval,iersub)
      if ( iersub .lt. 0 ) then
        write(iprint,*)'ERROR> INPBOU '
        write(iprint,*)'       DATA TYPE ERROR '
        write(iprint,'(a80)')line
        write(iprint,*)' '
        ier = iersub ; return
      endif
      fukcap_pro = rval(1)

!******************************

      return
      end subroutine rdbforp


!====================================================================


      subroutine rdbfunc(ier,iprint,line,icol)

!*******************************************************************
!
!     THIS SUBROUTINE IS FOR READ CAP FUNCTION TYPE
!
!*******************************************************************
 
      use COMBAS ; use COMMIS 
 
      implicit none

      ! Condition code
        integer(4),intent(out):: ier
      ! Logical unit number for output log
        integer(4),intent(in):: iprint
      ! Input line
        character(*),intent(in):: line
      ! Effective column number
        integer(4),intent(in):: icol

      integer(4),parameter:: maxval = 1
      integer(4):: ival(maxval),iwork1(maxval),iwork2(maxval)
      real(8):: rval(maxval)
      character(80):: cval(maxval)
      character(1):: inpmod(maxval)
      integer(4):: iersub
 
!*******************************************************************

      ier = 0 ; inpmod(1) = "R"
      call rdfree(line,icol,maxval,maxval,inpmod,iwork1,iwork2,rval,   &
                  ival,cval,iersub)
      if ( iersub .lt. 0 ) then
        write(iprint,*)'ERROR> INPBOU '
        write(iprint,*)'       DATA TYPE ERROR '
        write(iprint,'(a80)')line
        write(iprint,*)' '
        ier = iersub ; return
      endif
      if ( index(cval(1),"HARM") .ne. 0 ) then
        iufcap = 1
      elseif ( index(cval(1),"BIQU") .ne. 0 ) then
        iufcap = 2
      elseif ( index(cval(1),"HAXY") .ne. 0 ) then
        iufcap = 3
      elseif ( index(cval(1),"BIXY") .ne. 0 ) then
        iufcap = 4
      endif

!******************************

      return
      end subroutine rdbfunc
