
      subroutine rdlstatm(iprint,iread,chkname,outflg,ncou,ilst,ier)

!********************************************************
!
!     Read & search atom info. in list
!
!********************************************************

      use COMBAS ; use COMCMM,only: absres

      implicit none

      integer(4),intent(in):: iprint,iread
      integer(4),intent(out):: ier,ncou,ilst(ixnatm)
      character(6),intent(in):: chkname
      logical(4),intent(in):: outflg

      integer(4),parameter:: maxval = 7
      integer(4):: efcol,icol,ival(maxval),iwork1(maxval),iwork2(maxval)
      real(8):: rval(maxval)
      character(80):: line,space,cval(maxval)
      character(1):: inpmod(maxval)

      integer(4):: i,j,k,iersub,ichn,jchn,ires,jres
      character(4):: catm,cres,output,cc

!*********************************************

      ier = 0 ; iersub = 0 ; space = " " ; ncou = 0
      inpmod(1:4) = "I" ; inpmod(5:7) = "C"

      do
        read(iread,'(a80)',end=800)line
        icol = efcol(line,space,";")
        if ( index(line(1:icol),chkname) .ne. 0 ) then
          backspace(iread) ; exit
        endif
        if ( icol .le. 0 ) cycle

        call rdfree(line,icol,maxval,maxval,inpmod,iwork1,iwork2,rval, &
                    ival,cval,iersub)
        if ( iersub .lt. 0 ) then
          write(iprint,'(x,a7,a6)')"ERROR> ",chkname
          write(iprint,'(8x,a20)')"READ ERROR IN RDFREE"
          ier = iersub ; stop
        endif

        ichn = ival(1) ; jchn = ival(2)
        ires = ival(3) ; jres = ival(4)
        catm = cval(1) ; cres = cval(2) ; output = cval(3)

        ! Search atom pointer
        j = index(cres,"*") ; if ( j.eq.0 ) j = 5
        k = index(catm,"*") ; if ( k.eq.0 ) k = 5
        if ( ichn .eq. 0 ) then
          do i = 1,ixnatm
            if ( absres(i).lt.ires .or. absres(i).gt.jres .or.         &
                 ( j.ne.1 .and. cres(1:j-1).ne.cxresn(i)(1:j-1) ) ) cycle
            cc = cxatmn(i)
            if ( catm.eq."HEAV" .or. catm.eq."heav" ) then
              if (cc(1:1).eq."H" .or. cc(1:1).eq."h") cycle
            elseif ( catm.eq."SIDE" .or. catm.eq."side" ) then
              if ( cc.eq."C   " .or. cc.eq."CA  " .or.                 &
                   cc.eq."H   " .or. cc.eq."O   " .or.                 &
                   cc.eq."N   " .or. cc.eq."OX  " .or.                 &
                   cc.eq."O5' " .or. cc.eq."C5' " .or.                 &
                   cc.eq."C4' " .or. cc.eq."C3' " .or.                 &
                   cc.eq."O3' " .or. cc.eq."P   " .or.                 &
                   cc.eq."H5' " .or. cc.eq."H5''" .or.                 &
                   cc.eq."H4' " .or. cc.eq."H3' " .or.                 &
                   cc.eq."OP1 " .or. cc.eq."OP2 " ) cycle
            elseif ( catm.eq."RING" .or. catm.eq."ring" ) then
              if ( cc.eq."C   " .or. cc.eq."CA  " .or.                 &
                   cc.eq."H   " .or. cc.eq."O   " .or.                 &
                   cc.eq."N   " .or. cc.eq."OX  " .or.                 &
                   cc.eq."CB  " ) cycle
            elseif ( ( k.ne.1 .and. catm(1:k-1).ne.cc(1:k-1)) ) then
              cycle
            endif
            ncou = ncou + 1 ; ilst(ncou) = i
            if ( outflg .and. output(1:1).eq."Y" )                     &
              write(iprint,'(2i8,2x,a4,x,a4,2i8)')0,absres(i),         &
                cxresn(i),cc,ncou,ilst(ncou)
          enddo

        else
          do i = 1,ixnatm
            if ( ixachn(i).lt.ichn .or. ixachn(i).gt.jchn .or.         &
                 ixares(i).lt.ires .or. ixares(i).gt.jres .or.         &
                 ( j.ne.1 .and. cres(1:j-1).ne.cxresn(i)(1:j-1) ) ) cycle
            cc = cxatmn(i)
            if ( catm.eq."HEAV" .or. catm.eq."heav" ) then
              if ( cc(1:1).eq."H" .or. cc(1:1).eq."h" ) cycle
            elseif ( catm.eq."SIDE" .or. catm.eq."side" ) then
              if ( cc.eq."C   " .or. cc.eq."CA  " .or.                 &
                   cc.eq."H   " .or. cc.eq."O   " .or.                 &
                   cc.eq."N   " .or. cc.eq."OX  " .or.                 &
                   cc.eq."O5' " .or. cc.eq."C5' " .or.                 &
                   cc.eq."C4' " .or. cc.eq."C3' " .or.                 &
                   cc.eq."O3' " .or. cc.eq."P   " .or.                 &
                   cc.eq."H5' " .or. cc.eq."H5''" .or.                 &
                   cc.eq."H4' " .or. cc.eq."H3' " .or.                 &
                   cc.eq."OP1 " .or. cc.eq."OP2 " ) cycle
            elseif ( catm.eq."RING" .or. catm.eq."ring" ) then
              if ( cc.eq."C   " .or. cc.eq."CA  " .or.                 &
                   cc.eq."H   " .or. cc.eq."O   " .or.                 &
                   cc.eq."N   " .or. cc.eq."OX  " .or.                 &
                   cc.eq."CB  " ) cycle
            elseif ( ( k.ne.1 .and. catm(1:k-1).ne.cc(1:k-1)) ) then
              cycle
            endif
            ncou = ncou + 1 ; ilst(ncou) = i
            if ( outflg .and. output(1:1).eq."Y" )                     &
              write(iprint,'(2i8,2x,a4,x,a4,2i8)')ixachn(i),ixares(i), &
                cxresn(i),cc,ncou,ilst(ncou)
          enddo
        endif
        write(iprint,'(4i8,3x,2a4,i8)')ichn,jchn,ires,jres,catm,cres,  &
                                       ncou
      enddo

!****************************************

800   return
      end subroutine rdlstatm
