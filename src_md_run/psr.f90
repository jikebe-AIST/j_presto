
      subroutine psr(incons,iout,filena,ier)

!***********************************************************
!
!     ROUTINE TO SELECT ATOM FOR HARMONIC CONSTRAINTS
!     AND SEARCH ATOM POINTER( ABSOLUTE ATOM NUMBER )
!
!***********************************************************

       use COMBAS ; use COMMIS

      implicit none

      integer(4):: incons,iout
      integer(4),intent(out):: ier
      character(80),intent(in):: filena

      integer(4),parameter:: maxval = 20
      integer(4):: wrk1(maxval),wrk2(maxval),ival(maxval)
      character(1):: inpmod(maxval)
      real(8):: rval(maxval)
      character(80):: cval(maxval)

      integer(4):: natgrp,igrp,lnatm,lbatm,icol,efcol,ichain,jchain,   &
        nrsfst,nrsend,icntrs
      real(8):: delta,radmin,radmax
      character(80):: blank,line
      character(4):: comand,subcom,wpmass,atmnam,resnam,output,cntatm

!************************************************************

       ier = 0 ; allocate(Tiugpnt(ixnatm),Tfugcns(ixnatm))
       call flopen(incons,filena,10,'NULL',0,ier)
       if ( ier .ne. 0 ) then
         write(iout,'(a)')" ERROR> PSR"
         write(iout,'(8x,a,i3)')"OPEN FAILURE IN FLOPEN - UNIT ",incons
         return
       endif

       natgrp = 0 ; igrp = 0 ; lnatm = 0 ; lbatm = 0 ; blank = " "
       write(iout,'(a)')" INFORMATION> PSR"
       write(iout,'(8x,a)')"SELECT FOLLOWING ATOMS FOR "//             &
                           "POSITIONAL CONSTRAINTS."

!      MAIN LOOP TO READ PSR CARD
       OUTER : do
         read(incons,'(a80)',end=500)line
         icol = efcol(line,blank,";")
         if ( icol .le. 0 ) cycle
         inpmod(1:2) = "C"
         call rdfree(line,icol,maxval,2,inpmod,wrk1,wrk2,rval,ival,    &
                     cval,ier)
         comand = cval(1)(1:4) ; subcom = cval(2)(1:4)
         if ( comand .ne. "PSR>" ) cycle

         if ( subcom .eq. "STOP" ) then
           exit

         elseif ( subcom .eq. "LIST" ) then
           do
             read(incons,'(a80)',end=500)line
             icol = efcol(line,blank,";")
             if ( icol .le. 0 ) cycle
             inpmod(1) = "C"
             call rdfree(line,icol,maxval,2,inpmod,wrk1,wrk2,rval,ival,&
                         cval,ier)
             comand = cval(1)(1:4)
             if ( comand .eq. "END " ) cycle OUTER
             inpmod(1:4) = "I" ; inpmod(5:6) = "C" ; inpmod(7) = "R"
             inpmod(8:9) = "C"
             call rdfree(line,icol,maxval,9,inpmod,wrk1,wrk2,rval,ival,&
                         cval,ier)
             if ( ier .lt. 0 ) then
               write(iout,'(a)')" ERROR> PSR"
               write(iout,'(8x,a)')"READ ERROR IN RDFREE" ; return
             endif
             if ( rval(1) .le. 0.d0 ) then
               write(iout,'(a)')" ERROR> PSR"
               write(iout,'(8x,a)')"INPUT DATA ERROR. INPUT "//        &
                                   "POSITIVE VALUES."
               ier = -1840 ; return
             endif
             igrp = igrp + 1 ; ichain = ival(1) ; jchain = ival(2)
             nrsfst = ival(3) ; nrsend = ival(4) ; atmnam = cval(1)
             resnam = cval(2) ; delta = rval(1) ; wpmass = cval(3)(1:4)
             output = cval(4)
             ! search atom pointer for each psr
             lbatm = natgrp
             call findat(ichain,jchain,nrsfst,nrsend,atmnam,resnam,    &
                         delta,wpmass,output,natgrp,iout,ier)
             lnatm = natgrp
             write(iout,'(4i6,3x,2a4,f10.2,2x,a4,i6)')ichain,jchain,   &
               nrsfst,nrsend,atmnam,resnam,delta,wpmass,lnatm-lbatm
             if ( ier .lt. 0 ) return
           enddo

         elseif ( subcom .eq. "RADI" ) then
           do
             read(incons,'(a80)',end=500)line
             icol = efcol(line,blank,";")
             if ( icol .le. 0 ) cycle
             inpmod(1) = "C"
             call rdfree(line,icol,maxval,2,inpmod,wrk1,wrk2,rval,ival,&
                         cval,ier)
             comand = cval(1)(1:4)
             if ( comand .eq. "END " ) cycle OUTER
             inpmod(1:2) = "I" ; inpmod(3) = "C" ; inpmod(4:5) = "R"
             inpmod(6) = "C" ; inpmod(7) = "R" ; inpmod(8:9) = "C"
             call rdfree(line,icol,maxval,9,inpmod,wrk1,wrk2,rval,ival,&
                         cval,ier)
             if ( ier .lt. 0 ) then
               write(iout,'(a)')" ERROR> PSR"
               write(iout,'(8x,a)')"READ ERROR IN RDFREE" ; return
             endif
             if ( rval(2) .le. 0.d0 ) then
               write(iout,'(a)')" ERROR> PSR"
               write(iout,'(8x,a)')"INPUT DATA ERROR. INPUT "//        &
                                   "POSITIVE VALUES."
               ier = -1840 ; return
             endif
             igrp = igrp + 1 ; ichain = ival(1) ; icntrs = ival(2)
             cntatm = cval(1) ; atmnam = cval(2) ; radmin = rval(1)
             radmax = rval(2) ; delta = rval(3) ; wpmass = cval(3)(1:4)
             output = cval(4)(1:4)
             ! search atom pointer for each psr ( cut from input radius)
             lbatm = natgrp
             call cutrad(ichain,icntrs,cntatm,radmin,radmax,atmnam,    &
                         delta,wpmass,output,natgrp,iout,ier)
             lnatm = natgrp
             write(iout,'(2i6,3x,a4,x,2f8.2,2x,a4,x,f8.2,x,a4,i6)')    &
               ichain,icntrs,cntatm,radmin,radmax,atmnam,delta,wpmass, &
               lnatm-lbatm
             if ( ier .lt. 0 ) return
           enddo
         endif
       enddo OUTER

500    if ( igrp .le. 0 ) then
         write(iout,'(a)')"ERROR> PSR"
         write(iout,'(8x,a,i3)')"NO PSR CARD READ IN UNIT",incons
         ier = -1841 ; return
       endif
       write(iout,*)
       write(iout,'(8x,a,i6)')                                         &
         "TOTAL NUMBER OF ATOMS IN ALL PSR ARE :",natgrp
       iutatm = natgrp ; allocate(iugpnt(iutatm),fugcns(iutatm))
       iugpnt(1:iutatm) = Tiugpnt(1:iutatm)
       fugcns(1:iutatm) = Tfugcns(1:iutatm)
       deallocate(Tiugpnt,Tfugcns)

!***************************************

       return
       end subroutine psr


!======================================================================


      subroutine findat(ichain,jchain,nrsfst,nrsend,atmnam,resnam,     &
                        delta,wpmass,output,natgrp,iout,ier)

!********************************************************************
!
!     SEARCH ATOM POINTER FOR EACH PSR
!     IF MASS IS SPECIFED IN THE INPUT,THEN 1/DELTA IS MULTIPLIED BY
!     THE MASS OF THE ATOM.
!
!*******************************************************************

       use COMBAS ; use COMMIS ; use COMCMM

      implicit none

      integer(4),intent(in):: ichain,jchain,nrsfst,nrsend,iout
      integer(4),intent(out):: ier
      integer(4),intent(inout):: natgrp
      real(8),intent(in):: delta
      character(4),intent(in):: atmnam,resnam,wpmass,output

      integer(4):: i,j,k
      character(4):: cc

!********************************************
!     SEARCH ATOM POINTER FOR EACH PSR
      ier = 0
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
            if ( cc.eq."C   " .or. cc.eq."CA  " .or. cc.eq."H   " .or. &
                 cc.eq."O   " .or. cc.eq."N   " .or. cc.eq."OX  " .or. &
                 cc.eq."O5' " .or. cc.eq."C5' " .or. cc.eq."C4' " .or. &
                 cc.eq."C3' " .or. cc.eq."O3' " .or. cc.eq."P   " .or. &
                 cc.eq."H5' " .or. cc.eq."H5''" .or. cc.eq."H4' " .or. &
                 cc.eq."H3' " .or. cc.eq."OP1 " .or. cc.eq."OP2 " ) cycle
            elseif ( atmnam.eq."RING" .or. atmnam.eq."ring" ) then
              if ( cc.eq."C   " .or. cc.eq."CA  " .or.                 &
                   cc.eq."H   " .or. cc.eq."O   " .or.                 &
                   cc.eq."N   " .or. cc.eq."OX  " .or.                 &
                   cc.eq."CB  " ) cycle
          elseif ( ( k.ne.1 .and. atmnam(1:k-1).ne.cc(1:k-1)) ) then
            cycle
          endif
          natgrp = natgrp + 1 ; Tiugpnt(natgrp) = i
          if ( wpmass .eq. "MASS" ) then
            Tfugcns(natgrp) = fxmass(i) * delta
          else
            Tfugcns(natgrp) = delta
          endif
          if ( output(1:1) .eq. "Y" )                                  &
            write(iout,'(2i6,2x,a4,x,a4,2i6,f11.4)')0,absres(i),       &
              cxresn(i),cc,natgrp,Tiugpnt(natgrp),Tfugcns(natgrp)
        enddo
      else
        do i = 1,ixnatm
          if ( ixachn(i).lt.ichain .or. ixachn(i).gt.jchain .or.       &
               ixares(i).lt.nrsfst .or. ixares(i).gt.nrsend ) cycle
          if ( j.ne.1 .and. resnam(1:j-1).ne.cxresn(i)(1:j-1) ) cycle
          cc = cxatmn(i)
          if ( atmnam.eq."HEAV" .or. atmnam.eq."heav" ) then
            if ( cc(1:1).eq."H" .or. cc(1:1).eq."h" ) cycle
          elseif ( atmnam.eq."SIDE" .or. atmnam.eq."side" ) then
            if ( cc.eq."C   " .or. cc.eq."CA  " .or. cc.eq."H   " .or. &
                 cc.eq."O   " .or. cc.eq."N   " .or. cc.eq."OX  " .or. &
                 cc.eq."O5' " .or. cc.eq."C5' " .or. cc.eq."C4' " .or. &
                 cc.eq."C3' " .or. cc.eq."O3' " .or. cc.eq."P   " .or. &
                 cc.eq."H5' " .or. cc.eq."H5''" .or. cc.eq."H4' " .or. &
                 cc.eq."H3' " .or. cc.eq."OP1 " .or. cc.eq."OP2 " ) cycle
          elseif ( atmnam.eq."RING" .or. atmnam.eq."ring" ) then
            if ( cc.eq."C   " .or. cc.eq."CA  " .or.                   &
                 cc.eq."H   " .or. cc.eq."O   " .or.                   &
                 cc.eq."N   " .or. cc.eq."OX  " .or.                   &
                 cc.eq."CB  " ) cycle
          elseif ( ( k.ne.1 .and. atmnam(1:k-1).ne.cc(1:k-1)) ) then
            cycle
          endif
          natgrp = natgrp + 1 ; Tiugpnt(natgrp) = i
          if ( wpmass .eq. "MASS" ) then
            Tfugcns(natgrp) = fxmass(i) * delta
          else
            Tfugcns(natgrp) = delta
          endif
          if ( output(1:1) .eq. "Y" )                                  &
            write(iout,'(2i6,2x,a4,x,a4,2i6,f11.4)')ixachn(i),         &
              ixares(i),cxresn(i),cc,natgrp,Tiugpnt(natgrp),           &
              Tfugcns(natgrp)
        enddo
      endif

      if ( natgrp .gt. ixnatm ) then
        write(iout,'(a)')"ERROR> FINDAT"   
        write(iout,'(5x,a,i6,a)')"MAXIMUM PSR ATOM( ",ixnatm,        &
          ") EXCEEDED" ; ier = -1842
      elseif ( natgrp .le. 0 ) then
        write(iout,'(a)')"ERROR> FINDAT"   
        write(iout,'(5x,a)')"ERROR FINDING ATOM IN PSR"   
        write(iout,'(4i6,3x,2a4,f7.2,2x,a4)')ichain,jchain,nrsfst,     &
          nrsend,atmnam,resnam,delta,wpmass ; ier = -1843
      endif

!************************************

      return
      end subroutine findat


!======================================================================


      subroutine cutrad(ichain,icntrs,cntatm,radmin,radmax,atmnam,     &
                        delta,wpmass,output,natgrp,iout,ier)

!**********************************************************************
!
!     SEARCH ATOM POINTER OF ATOMS WITHIN THE INPUT RADIUS FROM THE CENTER
!     IF MASS IS SPECIFED IN THE INPUT,THEN 1/DELTA IS MULTIPLIED BY
!     THE MASS OF THE ATOM.
!
!**********************************************************************

      use COMBAS ; use COMMIS ; use COMCMMC

      implicit none

      integer(4),intent(in):: ichain,icntrs,iout
      integer(4),intent(out):: ier
      integer(4),intent(inout):: natgrp
      real(8),intent(in):: radmin,radmax,delta
      character(4):: cntatm,atmnam,wpmass,output

      integer(4):: i,icount
      real(8):: ccod(3),dcod(3),dist,rmx2,rmn2

!****************************************
!     SEARCH THE CENTER ATOM
      icount = 0
      do i = 1,ixnatm
        if ( ixachn(i).ne.ichain .or. ixares(i).ne.icntrs .or.         &
             cxatmn(i)(1:4).ne.cntatm ) cycle
        icount = icount + 1
        ccod(1:3) = cord(1:3,i)
      enddo
      if ( icount .ne. 1 ) then
        write(iout,'(a)')" ERROR> CUTRAD"
        write(iout,'(5x,a)')"THE CENTER ATOM WAS NOT FOUND"
        ier = -1844 ; return
      endif

!     SEARCH ATOMS WITHIN THE INPUT RADIUS FROM THE CENTER
      rmx2 = radmax*radmax ; rmn2 = radmin*radmin
      do i = 1,ixnatm
        dcod(1:3) = cord(1:3,i) - ccod(1:3)
        dist = dot_product(dcod,dcod)
        if ( dist.gt.rmx2 .or. dist.lt.rmn2 ) cycle
        dist = sqrt(dist)
        if ( atmnam(1:1).eq."*" .or.                                   &
             (atmnam(2:2).eq."*" .and. atmnam(1:1).eq.cxatmn(i)(1:1))  &
             .or.                                                      &
             (atmnam(3:3).eq."*" .and. atmnam(1:2).eq.cxatmn(i)(1:2))  &
             .or. atmnam.eq.cxatmn(i) ) then
          natgrp = natgrp + 1 ; Tiugpnt(natgrp) = i
          if ( wpmass .eq. "MASS" ) then
            Tfugcns(natgrp) = fxmass(i) * delta
          else
            Tfugcns(natgrp) = delta
          endif
          if ( output(1:1) .eq. "Y" )                                  &
            write(iout,'(2i6,2x,a4,x,a4,2i6,2f11.4)')ixachn(i),        &
              ixares(i),cxresn(i),cxatmn(i),natgrp,Tiugpnt(natgrp),    &
              Tfugcns(natgrp),dist
        endif
      enddo

      if ( natgrp .gt. ixnatm ) then
        write(iout,'(a)')" ERROR> CUTRAD"
        write(iout,'(5x,a)')"THE CENTER ATOM WAS NOT FOUND"
        ier = -1842
      elseif ( natgrp .le. 0 ) then
        write(iout,'(a)')" ERROR> CUTRAD"
        write(iout,'(5x,a)')"ERROR FINDING ATOM IN PSR"
        write(iout,'(2i6,3x,a4,2f7.2,2x,f7.2,2(x,a4))')ichain,icntrs,  &
          cntatm,radmin,radmax,delta,atmnam,wpmass
        ier = -1843
      endif

!*********************************************

      return
      end subroutine cutrad
