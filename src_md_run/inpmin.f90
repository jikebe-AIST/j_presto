
      subroutine inpmin(iread,ier,iprint,onend,inmthd,inlopl,inupdl,   &
                        inmons,inflog,infbst,fncovg,fnstpi,fnupsl,     &
                        fndnsl,fncovl,inllpl,inuana,cnnana)

!*******************************************************************
!
!     THIS SUBROUTIne IS FOR READ CONTROL DATA FOR MINMIZATION
!
!*******************************************************************

      use COMBAS ; use COMERG ; use COMMIS ; use COMCMM ; use PHYCNS

      implicit none

      integer(4),intent(in):: iread,iprint
      integer(4),intent(out):: ier
      logical(1),intent(inout):: onend
      integer(4),intent(inout):: inmthd,inlopl,inupdl,inmons,inflog,   &
                                 infbst,inllpl,inuana
      real(8),intent(inout):: fncovg,fnstpi,fnupsl,fndnsl,fncovl
      character(*),intent(inout):: cnnana

      integer(4),parameter:: numele = 55
      integer(4),parameter:: mxelel = 10
      integer(4),parameter:: maxcon = 4

      character(6):: elemnt(numele) = (/                               &
        'METHOD','CPUTIM','LOOPLI','UPDATE','MONITO','LOGFOR','BESTFI',&
        'CONVGR','ISTEPL','UPRATE','DOWNRA','LINESE','CONVLI','NAMEAN',&
        'CAL15M','CUTLEN','DIEVAL','TEMPER','WETPSR','WETDSR','WETDHR',&
        'RADCAP','FORCAP','FUNCAP','SHAKEM','COVSHK','LIMSHK','LXCELL',&
        'LYCELL','LZCELL','BOUNDA','SETCEN','CENTRX','CENTRY','CENTRZ',&
        'ELLIPA','ELLIPB','ELLIPC','RADIUS','FOACAP','RADPRO','CAPSHP',&
        'XMINCL','XMAXCL','YMINCL','YMAXCL','ZMINCL','ZMAXCL','CELBUF',&
        'ELLBUF','LAMBDA','LMDSQU','PARALV','N_SIMD','NSTPCN'/)

      character(1):: dataty(numele) = (/                               &
          'C','R','I','I','I','C','C',   'R','R','R','R','I','R','D',  &
          'C','R','R','R','R','R','R',   'R','R','C','C','R','I','R',  &
          'R','R','C','C','R','R','R',   'R','R','R','R','R','R','C',  &
          'R','R','R','R','R','R','R',   'R','R','R','C','I','I'/)

      character(4):: datcon(maxcon,numele)
      data datcon(1:2,1) /'STEE','CONJ'/
      data datcon(1:2,6) /'SHOR','DETA'/
      data datcon(1:2,7) /'NO  ','YES '/
      data datcon(1:4,15) /'CMM ','ATOM','RESI','ZD  '/
      data datcon(1:4,24) /'HARM','BIQU','HAXY','BIXY'/
      data datcon(1:2,25) /'HBON','ALLB'/
      data datcon(1:3,31) /'NO  ','PBC ','ELLI'/
      data datcon(1:2,32) /'NO  ','YES '/
      data datcon(1:3,42) /'SPHE','BOX ','ELL '/
      data datcon(1:3,53) /'REDU','HIGH','DOUB'/

      integer(4):: numspe(numele) = (/                                 &
        2,0,0,0,0,2,2,  0,0,0,0,0,0,0,  4,0,0,0,0,0,0,  0,0,4,2,0,0,0, &
        0,0,3,2,0,0,0,  0,0,0,0,0,0,3,  0,0,0,0,0,0,0,  0,0,0,3,0,0/)

      character(80):: line,space
      integer(4):: posele(mxelel),iwork(mxelel),intele(numele)
      real(8):: reaele(numele)
      character(80):: chaele(numele)

      logical(1):: onlist

      integer(4),parameter:: maxvl = 20
      character(1):: inpmod(maxvl)
      integer(4):: iwork1(maxvl),iwork2(maxvl),ival(maxvl)
      real(8):: rval(maxvl)
      character(80):: cval(maxvl)

      integer(4):: efcol,numerr,icol,nelel,iele,icon,numval,numint,    &
                   numrea,numcha,ierr
      real(8):: celbuf,ellbuf,rtmp
      
      space  = ' '

!******************************************************************
!     <<<  SET DEFAULT VALUES  >>>

      ier = 0 ; numerr = 0 
      onend = .false. ; onlist = .false.

      intele(1) = 1
      reaele(2) = fxcpul
      intele(3) = 0 ; intele(4) = 10 ; intele(5) = 10
      intele(6) = 1 ; intele(7) = 1
      reaele(8) = 0.1d0 ; reaele(9) = 0.01d0 ; reaele(10) = 1.2d0
      reaele(11) = 0.6d0
      intele(12) = 10
      reaele(13) = 0.1d0
      chaele(14) = " "
      intele(15) = iy15m
      reaele(16) = fycutl ; reaele(17) = fydiel ; reaele(18) = fuctmp
      reaele(19) = fuwpsc ; reaele(20) = fuwdsc ; reaele(21) = fuwdhc
      reaele(22) = furcap ; reaele(23) = fukcap
      intele(24) = iufcap ; intele(25) = iyfshk
      reaele(26) = fustol
      intele(27) = iuslop
      reaele(28:30) = fxcell(1:3)
      intele(31) = ixfbou+1 ; intele(32) = ixcbou+1
      reaele(33:35) = fxcbou(1:3) ; reaele(36:38) = fxellp(1:3)
      reaele(39) = fxellp(1) ; reaele(40) = fukcap_pro
      reaele(41) = CAPbuff
      intele(42) = capshp
      reaele(43:48) = celwal(1:6) ; reaele(49) = CAPbuff
      reaele(50) = CAPbuff ; reaele(51:52) = 0.d0
      intele(53) = 2 ; intele(54) = 0 ; intele(55) = 1

!***********************************************************
!     <<<  READ CONTROL DATA FOR EneRGY MINIMIZATION  >>>

100   do
        read(iread,'(a80)',iostat=ierr,end=200)line
        icol = efcol(line,space,";")
        if ( icol .le. 0 ) cycle

        if ( index(line(1:icol),"QUIT") .ne. 0 ) then
          exit
        elseif ( index(line(1:icol),"EXE>") .ne. 0 ) then
          backspace(iread) ; exit
        elseif ( index(line(1:icol),"LIST") .ne. 0 ) then
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
          inpmod(2*iele-1) = "C"
          if ( dataty(posele(iele)) .eq. "D" ) then
            inpmod(2*iele) = "C"
          else
            inpmod(2*iele) = dataty(posele(iele))
          endif
        enddo
        numval = nelel * 2
        call rdfree(line,icol,maxvl,numval,inpmod,iwork1,iwork2,rval, &
                    ival,cval,ier)
        if ( ier .ne. 0 ) then
          write(iprint,*)'ERROR> INPMIN '
          write(iprint,*)'    DATA TYPE ERROR '
          write(iprint,'(a80)')line
          return
        endif

        ! 1-3) STORE DATA OF EACH ELEMENT
        numint = 0 ; numrea = 0 ; numcha = 0
        do iele = 1,nelel
          if ( dataty(posele(iele)) .eq. "I" ) then
            numint = numint + 1 ; numcha = numcha + 1
            intele(posele(iele)) = ival(numint)
          elseif ( dataty(posele(iele)) .eq. "R" ) then
            numrea = numrea + 1 ; numcha = numcha + 1
            reaele(posele(iele)) = rval(numrea)
          elseif ( dataty(posele(iele)) .eq. "C" ) then
            numcha = numcha + 2
            do icon = 1,numspe(posele(iele))
              if ( cval(numcha)(1:4) .eq. datcon(icon,posele(iele)) ) &
                intele(posele(iele)) = icon
            enddo
          elseif ( dataty(posele(iele)) .eq. "D" ) then
            numcha = numcha + 2
            chaele(posele(iele)) = cval(numcha)
          endif
        enddo

      enddo

200   if ( ierr .ne. 0 ) onend = .true.

!*************************************************************************
!     <<<  STORE ELEMENT DATA TO VARIABLE OF EneRGY MINIMIZATION >>>

      inmthd = intele(1)  ; fxcpul = reaele(2)  ; inlopl = intele(3)
      inupdl = intele(4)  ; inmons = intele(5)  ; inflog = intele(6)
      infbst = intele(7)  ; fncovg = reaele(8)  ; fnstpi = reaele(9)
      fnupsl = reaele(10) ; fndnsl = reaele(11) ; inllpl = intele(12)
      fncovl = reaele(13) ; cnnana = chaele(14) ; iy15m = intele(15)
      fycutl = reaele(16) ; fydiel = reaele(17) ; fuctmp = reaele(18)
      fuwpsc = reaele(19) ; fuwdsc = reaele(20) ; fuwdhc = reaele(21)
      furcap = reaele(22) ; fukcap = reaele(23) ; iufcap = intele(24)
      iyfshk = intele(25) ; fustol = reaele(26) ; iuslop = intele(27)
      fxcell(1:3) = reaele(28:30) ; ixfbou = intele(31)-1
      ixcbou = intele(32)-1 ; fxcbou(1:3) = reaele(33:35)
      fxellp(1:3) = reaele(36:38)
      if ( ixfbou .eq. 3 ) then
        fxellp(1:3) = reaele(39) ; ixfbou = 2
      endif
      fukcap_pro = reaele(40)
      if ( reaele(41) .gt. 0.d0 ) CAPbuff = reaele(41)
      furcap_pro = furcap - CAPbuff
      furcap2 = furcap*furcap ; furcap_pro2 = furcap_pro*furcap_pro
      capshp = intele(42) ; celwal(1:6) = reaele(43:48)
      celbuf = reaele(49) ; ellbuf = reaele(50)
      if ( cluster_method .eq. "LMD" ) then
        if ( reaele(51) .ne. 0.d0 ) then
          lambda = reaele(51)
        elseif ( reaele(52) .ne. 0.d0 ) then
          lambda = sqrt(reaele(52))
        endif
      elseif ( cluster_method .eq. "LMD2" ) then
        if ( reaele(51) .ne. 0.d0 ) then
          lambda = reaele(51)*reaele(51)
        elseif ( reaele(52) .ne. 0.d0 ) then
          lambda = reaele(52)
        endif
      endif
      high_para = intele(53)
      i_vec = intele(54)
      nstpcn = intele(55)

      ! Logical Unit
      inuana = 0
      if ( cnnana .ne. space ) inuana = 30

!****************************************************************************
!     <<<  OUTPUT CONTROL DATA OF EneRGY MINIMIZATION  >>>

      write(iprint,*)'      PARAMETERS FOR MINIMIZATION '
      write(iprint,*)' '
      write(iprint,*)'      1) GEneRAL PARAMETERS  '
      write(iprint,*)'         METHOD            : ',datcon(inmthd,1)
      write(iprint,*)'         CPU LIMIT (S)     : ',fxcpul
      write(iprint,*)'         LOOP LIMIT        : ',inlopl
      write(iprint,*)'         UPDATE            : ',inupdl
      write(iprint,*)'         MONITORING LOG    : ',inmons
      write(iprint,*)'         LOG FORMAT        : ',datcon(inflog,6)
      write(iprint,*)'         CALC. OF RMSD     : ',datcon(infbst,7)
      write(iprint,*)'         CONVERGE (RMSF)   : ',fncovg
      write(iprint,*)'         INITIAL STEP(A)   : ',fnstpi
      if ( inmthd .eq. 1 ) then
        write(iprint,*)'         UP RATE           : ',fnupsl
        write(iprint,*)'         DOWN RATE         : ',fndnsl
      elseif ( inmthd .eq. 2 ) then
        write(iprint,*)'         LINE SEARCH LIMIT : ',inllpl
        write(iprint,*)'         LINE CONV.RATE    : ',fncovl
      endif
      if ( nomp .eq. 1 ) then
        high_para = 0
      else
        if ( high_para .eq. 1 ) then
          high_para = 0
        elseif ( high_para .eq. 2 ) then
          high_para = 1
        elseif ( high_para .eq. 3 ) then
          high_para = 2
        endif
      endif

      if ( inuana .gt. 0 ) then
        write(iprint,*)'      2) ANALYSIS FILE '
        write(iprint,*)'         LOGICAL UNIT      : ',inuana
        write(iprint,'(a80)')cnnana
      endif

      write(iprint,*)'      3) PARAMETERS FOR EneRGY CALCULATION '

      if ( iy15m .eq. 1 ) then
        write(iprint,*)                                                &
              '         CALCULATE NON-BONDED INTERACTIONS BY CMM'
        write(iprint,*)                                                &
          "         VDW CUT-OFF DIST. FOR NO-EXACT CALC. CELL :",fycutl
        if ( ixfbou .eq. 1 ) then
          write(iprint,*)" You cannot choose CMM under periodic "//    &
          "boundary condition." ; stop
        endif
      elseif ( iy15m .eq. 3 ) then
        write(iprint,*)" SORRY!!"
        write(iprint,*)" Residue-based CUT-OFF method is removed."
        stop
      else
        write(iprint,*)'         CUT-OFF LOGIC   : ',datcon(iy15m,15)
        write(iprint,*)'         CUT-OFF DIST.   : ',fycutl
      endif

      write(iprint,*)'         DIELECTRIC VAL.   : ',fydiel

      if ( fuctmp .eq. 0.d0 ) then
        fuctmp = 1.d0 / boltzk 
      elseif ( maxval(iyeflg(11:14)) .eq. 1 ) then
        write(iprint,*) '         TEMP. FOR REST. : ',fuctmp
      endif
      if ( iyeflg(11) .eq. 1 ) then
        write(iprint,*)'         WEIGHT OF POSR  : ',fuwpsc
        fugcns(1:iutatm) = fugcns(1:iutatm)*boltzk*fuctmp*fuwpsc
      endif
      if ( iyeflg(12) .eq. 1 )                                         &
        write(iprint,*)'         WEIGHT OF DISTR : ',fuwdsc
      if ( iyeflg(13) .eq. 1 )                                         &
        write(iprint,*)'         WEIGHT OF DHR   : ',fuwdhc
      if ( iyeflg(14) .eq. 1 ) then
        if ( capshp .eq. 1 ) then
          write(iprint,*)'         CAP SHAPE       :  SPHERE'
          write(iprint,*)'         RADIUS OF CAP   : ',furcap
          write(iprint,*)'         FORCE CONS. CAP : ',fukcap
          write(iprint,*)                                              &
            '         RADIUS OF CAP for specified atoms: ',furcap_pro
          write(iprint,*)                                              &
            '         SQUARE OF THE RADIUS             : ',furcap_pro2
          write(iprint,*)                                              &
          '         FORCE CONS. CAP for specified atoms: ',fukcap_pro
          write(iprint,*)'         FUNC. TYPE(CAP) : ',datcon(iufcap,24)
          write(iprint,*)'          CAP CENTER (X) : ',fxcbou(1)
          write(iprint,*)'          CAP CENTER (Y) : ',fxcbou(2)
          write(iprint,*)'          CAP CENTER (Z) : ',fxcbou(3)

        elseif ( capshp .eq. 2 ) then
          ! Determination of celwal, fxcell, fxcbou
          if ( celwal(1).eq.0.d0 .and. celwal(2).eq.0.d0 ) then
            celwal(1) = fxcbou(1) - fxcell(1) * 0.5d0
            celwal(2) = fxcbou(1) + fxcell(1) * 0.5d0
            celwal(3) = fxcbou(2) - fxcell(2) * 0.5d0
            celwal(4) = fxcbou(2) + fxcell(2) * 0.5d0
            celwal(5) = fxcbou(3) - fxcell(3) * 0.5d0
            celwal(6) = fxcbou(3) + fxcell(3) * 0.5d0
          else
            fxcell(1) = celwal(2) - celwal(1)
            fxcell(2) = celwal(4) - celwal(3)
            fxcell(3) = celwal(6) - celwal(5)
            fxcbou(1) = (celwal(1)+celwal(2)) * 0.5d0
            fxcbou(2) = (celwal(3)+celwal(4)) * 0.5d0
            fxcbou(3) = (celwal(5)+celwal(6)) * 0.5d0
          endif
          celwal_pro(1) = celwal(1) + celbuf
          celwal_pro(2) = celwal(2) - celbuf
          celwal_pro(3) = celwal(3) + celbuf
          celwal_pro(4) = celwal(4) - celbuf
          celwal_pro(5) = celwal(5) + celbuf
          celwal_pro(6) = celwal(6) - celbuf
          write(iprint,*)'         CAP SHAPE       :  BOX'
          write(iprint,*)'         CAP CENTER (X)  : ',fxcbou(1)
          write(iprint,*)'         CAP CENTER (Y)  : ',fxcbou(2)
          write(iprint,*)'         CAP CENTER (Z)  : ',fxcbou(3)
          write(iprint,*)'         CELL SIZE (X)   : ',fxcell(1)
          write(iprint,*)'         CELL SIZE (Y)   : ',fxcell(2)
          write(iprint,*)'         CELL SIZE (Z)   : ',fxcell(3)
          write(iprint,*)'         FORCE CONS. CAP : ',fukcap
          write(iprint,*)'         CELL SIZE for specific atoms (X): ',&
                      fxcell(1)-2.d0*celbuf
          write(iprint,*)'         CELL SIZE for specific atoms (Y): ',&
                      fxcell(2)-2.d0*celbuf
          write(iprint,*)'         CELL SIZE for specific atoms (Z): ',&
                      fxcell(3)-2.d0*celbuf
          write(iprint,*)                                              &
            '         FORCE CONS. CAP for specified atoms: ',fukcap_pro
          write(iprint,*)'         FUNC. TYPE(CAP) : ',datcon(iufcap,24)

        elseif ( capshp .eq. 3 ) then
          write(iprint,*)'         CAP SHAPE       :  ELL'
          write(iprint,*)'         CAP CENTER (X)  : ',fxcbou(1)
          write(iprint,*)'         CAP CENTER (Y)  : ',fxcbou(2)
          write(iprint,*)'         CAP CENTER (Z)  : ',fxcbou(3)
          write(iprint,*)'         CAP RADIUS (X)  : ',fxellp(1)
          write(iprint,*)'         CAP RADIUS (Y)  : ',fxellp(2)
          write(iprint,*)'         CAP RADIUS (Z)  : ',fxellp(3)
          write(iprint,*)'         FORCE CONS. CAP : ',fukcap
          ! Determination of fxellp_pro
          rtmp = minval(fxellp(1:3))
          rtmp = (rtmp - ellbuf)/ rtmp
          fxellp_pro(1:3) = fxellp(1:3) * rtmp
          write(iprint,*)'         CAP RADIUS for specific atoms (X): '&
                      ,fxellp_pro(1)
          write(iprint,*)'         CAP RADIUS for specific atoms (Y): '&
                      ,fxellp_pro(2)
          write(iprint,*)'         CAP RADIUS for specific atoms (Z): '&
                      ,fxellp_pro(3)
          write(iprint,*)                                              &
            '         FORCE CONS. CAP for specified atoms: ',fukcap_pro
          write(iprint,*)'         FUNC. TYPE(CAP) : ',datcon(iufcap,24)
          fxellp(1:3) = 1.d0 / (fxellp(1:3)*fxellp(1:3))
          fxellp_pro(1:3) = 1.d0 / (fxellp_pro(1:3)*fxellp_pro(1:3))
        endif
      endif
      if ( iyeflg(15) .eq. 1 )                                         &
        write(iprint,*)'         COEFFICIENT FOR REPULSION : ',Krep

      if ( ixfbou .eq. 1 ) then
        ! Determination of celwal, fxcell, fxcbou
        if ( celwal(1).eq.0.d0 .and. celwal(2).eq.0.d0 ) then
          celwal(1) = fxcbou(1) - fxcell(1) * 0.5d0
          celwal(2) = fxcbou(1) + fxcell(1) * 0.5d0
          celwal(3) = fxcbou(2) - fxcell(2) * 0.5d0
          celwal(4) = fxcbou(2) + fxcell(2) * 0.5d0
          celwal(5) = fxcbou(3) - fxcell(3) * 0.5d0
          celwal(6) = fxcbou(3) + fxcell(3) * 0.5d0
        else
          fxcell(1) = celwal(2) - celwal(1)
          fxcell(2) = celwal(4) - celwal(3)
          fxcell(3) = celwal(6) - celwal(5)
          fxcbou(1) = (celwal(1)+celwal(2)) * 0.5d0
          fxcbou(2) = (celwal(3)+celwal(4)) * 0.5d0
          fxcbou(3) = (celwal(5)+celwal(6)) * 0.5d0
        endif
        write(iprint,*)'         APPLY PERIODIC BOUNDARY'
        write(iprint,*)'          X AXIS           : ',fxcell(1)
        write(iprint,*)'          Y AXIS           : ',fxcell(2)
        write(iprint,*)'          Z AXIS           : ',fxcell(3)
        if ( ixcbou .ne. 1 ) then
          write(iprint,*)'         CENTER OF BOX IS INPUT'
          write(iprint,*)'          X AXIS           : ',fxcbou(1)
          write(iprint,*)'          Y AXIS           : ',fxcbou(2)
          write(iprint,*)'          Z AXIS           : ',fxcbou(3)
        else
          write(iprint,*)                                              &
          '         CENTER OF BOX IS SET TO THE CENTER OF FIRST CHAIN '
        endif
        if ( minval(fxcell(1:3)) .lt. 2.d0*fycutl ) then
          write(iprint,*)"          The cell size is too small "//     &
            "compared to cut-off dist."
          ier = -1 ; return
        endif

      ! elastic boundary option
      elseif( ixfbou .eq. 2 ) then
        write(iprint,*)'         APPLY ELASTIC ELLIPSOIDAL BOUNDARY'
        write(iprint,*)'          CENTER(X)       : ',fxcbou(1)
        write(iprint,*)'          CENTER(Y)       : ',fxcbou(2)
        write(iprint,*)'          CENTER(Z)       : ',fxcbou(3)
        write(iprint,*)'          ELLIPSOID AXIS;A: ',fxellp(1)
        write(iprint,*)'          ELLIPSOID AXIS;B: ',fxellp(2)
        write(iprint,*)'          ELLIPSOID AXIS;C: ',fxellp(3)
      endif

      if ( iyfshk .ne. 0 ) then
        write(iprint,*)'         APPLY SHAKE CONSTRAINTS'
        write(iprint,*)'          SHAKE METHOD  : ',datcon(iyfshk,25)
        write(iprint,*)'          LOOP LIMIT    : ',iuslop
        write(iprint,*)'          CRI. CONV.    : ',fustol

        if ( iutshk.eq.0 .or. iugshk.eq.0 ) then
          write(iprint,*)' '
          write(iprint,*)'ERROR> INPMIN'
          write(iprint,*)'          SHAKE parameters are defined,'
          write(iprint,*)'          but no SHAKE setting file is read.'
          write(iprint,*)' '
          ier = -1 ; return
        endif
      else
        if ( iutshk.ne.0 .and. iugshk.ne.0 ) then
          write(iprint,*)' '
          write(iprint,*)'ERROR> INPMIN'
          write(iprint,*)'          SHAKE setting file is read,'
          write(iprint,*)                                              &
                '          but no SHAKE parameters are defined. '
          write(iprint,*)' '
          ier = -1 ; return
        endif
      endif

      if ( cluster_method .eq. "LMD" ) then
        write(iprint,*)
        write(iprint,*)"              LAMBDA     : ",lambda
        write(iprint,*)
      elseif ( cluster_method .eq. "LMD2" ) then
        write(iprint,*)
        write(iprint,*)"        LAMBDA SQUARE    : ",lambda
        write(iprint,*)
      endif

      write(iprint,*)' '
      if ( onlist ) then
        onlist = .false. ; goto 100
      endif

      if ( ixfbou .eq. 1 ) then
        invcel(1:3) = 1.d0 / fxcell(1:3)
      else
        invcel(1:3) = 0.d0
      endif

!**************************************************

      return
      end subroutine inpmin
