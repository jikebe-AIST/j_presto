
      subroutine inpmd(iread,ier,iprint,onend,idfrst,idlopl,idupdl,    &
                       idmonl,idflog,fddelt,idfstp,idfvsl,fdtemp,      &
                       fdtaut,idfvel,idseed,fdtmps,idfbst,iducrd,      &
                       iduvel,iduerg,idutrj,idfacc,idfacv,idface,      &
                       idfact,cdnrti,cdnrto,cdncrd,cdnvel,cdnerg,      &
                       cdntrj,idmonc,idmonv,idmone,idmnts,idutte,      &
                       cdntte,idfacm,idhtlp,idmylp,idutce,cdntce,      &
                       idftce,idmoce,idupcm)

!*******************************************************************
!
!     THIS SUBROUTINE IS FOR READ CONTROL DATA FOR MOLECULAR DYNAMICS
!
!*******************************************************************

      use COMBAS ; use COMERG ; use COMMIS ; use COMCMM
      use COMPSC ; use COMMOL ; use PHYCNS

      implicit none

      integer(4),intent(in):: iread,iprint
      logical(1),intent(inout):: onend
      integer(4),intent(inout):: idfrst,idlopl,idupdl,idmonc,idmonv,   &
        idmone,idmnts,idfacc,idfacv,idface,idfact,idfacm,idmonl,idflog,&
        idutce,idftce,idmoce,idupcm,idfstp,idfvsl,idfvel,idseed,idfbst,&
        iducrd,iduvel,iduerg,idutrj,idutte,idmylp,idhtlp,ier
      real(8),intent(inout):: fddelt,fdtemp,fdtaut,fdtmps
      character(*),intent(inout):: cdnrti,cdnrto,cdncrd,cdnvel,cdnerg, &
                                   cdntrj,cdntte,cdntce
      character(80):: cdmpf

      integer(4),parameter:: numele = 85
      integer(4),parameter:: mxelel = 10
      integer(4),parameter:: maxcon = 4

      character(6):: elemnt(numele) = (/                               &
        'RESTAR','CPUTIM','LOOPLI','UPDATE','TIMEST','OUTLOG','LOGFOR',&
        'STOPCE','METHOD','SETTEM','RELAXA','INITIA','RANDOM','STARTT',&
        'BESTFI','NAMERI','NAMERO','NAMECO','NAMEVE','NAMEEN','MNTRCO',&
        'CAL15M','CUTLEN','DIEVAL','TEMPER','WETPSR','WETDSR','WETDHR',&
        'RADCAP','FORCAP','FUNCAP','SHAKEM','COVSHK','LIMSHK','LXCELL',&
        'LYCELL','LZCELL','BOUNDA','SETCEN','CENTRX','CENTRY','CENTRZ',&
        'OUTCOO','OUTVEL','OUTENE','ELLIPA','ELLIPB','ELLIPC','OUTTRJ',&
        'NAMEGE','NAMETR','RADIUS','MNTRVE','MNTREN','MNTRTR','NAMETO',&
        'MNTRTO','HEATLO','DUMMYL','NAMECL','MNTRCL','OUTCLE','TEMPCO',&
        'FOACAP','RADPRO','CMMUPD','CAPSHP','XMINCL','XMAXCL','YMINCL',&
        'YMAXCL','ZMINCL','ZMAXCL','CELBUF','PSETMP','STAPST','HTPSLO',&
        'ELLBUF','LAMBDA','LMDSQU','OUTATM','PARALV','N_SIMD','NSTPCN',&
        'FDEBUG'/)

      character(1):: dataty(numele) = (/                               &
          'C','R','I','I','R','I','C',   'C','C','R','R','C','I','R',  &
          'C','D','D','D','D','D','C',   'C','R','R','R','R','R','R',  &
          'R','R','C','C','R','I','R',   'R','R','C','C','R','R','R',  &
          'I','I','I','R','R','R','I',   'D','D','R','C','C','C','D',  &
          'C','I','I','D','C','I','C',   'R','R','I','C','R','R','R',  &
          'R','R','R','R','R','R','I',   'R','R','R','I','C','I','I',  &
          'C'/)

      character(4):: datcon(maxcon,numele)
      data datcon(1:2,1) /'NO  ','YES '/                ! RESTAR
      data datcon(1:2,7) /'SHOR','DETA'/                ! LOGFOR
      data datcon(1:4,8) /'TRAN','ROTA','BOTH','NO  '/  ! STOPCE
      data datcon(1:4,9) /'CONS','MICR','CANO','GE  '/  ! METHOD
      data datcon(1:2,12) /'ZERO','SET '/               ! INITIA
      data datcon(1:2,15) /'NO  ','YES '/               ! BESTFI
      data datcon(1:4,21) /'NO  ','ASCI','SING','DOUB'/ ! MNTRCO
      data datcon(1:4,22) /'CMM ','ATOM','RESI','ZD  '/ ! CAL15M
      data datcon(1:4,31) /'HARM','BIQU','HAXY','BIXY'/ ! FUNCAP
      data datcon(1:2,32) /'HBON','ALLB'/               ! SHAKEM
      data datcon(1:3,38) /'NO  ','PBC ','ELLI'/        ! BOUNDA
      data datcon(1:2,39) /'NO  ','YES '/               ! SETCEN
      data datcon(1:4,53:55) /'NO  ','ASCI','SING','DOUB',             &
                              'NO  ','ASCI','SING','DOUB',             &
                              'NO  ','ASCI','SING','DOUB'/ ! MNTRVE MNTREN MNTRTR
      data datcon(1:4,57) /'NO  ','ASCI','SING','DOUB'/ ! MNTRTO
      data datcon(1:4,61) /'NO  ','ASCI','SING','DOUB'/ ! MNTRCL
      data datcon(1:2,63) /'NO  ','YES '/               ! TEMPCO
      data datcon(1:3,67) /'SPHE','BOX ','ELL '/        ! CAPSHP
      data datcon(1:3,82) /'REDU','HIGH','DOUB'/        ! PARALV
      data datcon(1:3,85) /'NO  ','YES ','MAKE'/        ! FDEBUG

      integer(4):: numspe(numele) = (/                                 &
        2,0,0,0,0,0,2,  4,4,0,0,2,0,0,  2,0,0,0,0,0,4,  4,0,0,0,0,0,0, &
        0,0,4,2,0,0,0,  0,0,3,2,0,0,0,  0,0,0,0,0,0,0,  0,0,0,4,4,4,0, &
        4,0,0,0,4,0,2,  0,0,0,3,0,0,0,  0,0,0,0,0,0,0,  0,0,0,0,3,0,0, &
        3/)

      character(80):: line,space
      integer(4):: posele(mxelel),iwork(mxelel),intele(numele)
      real(8):: reaele(numele)
      character(80):: chaele(numele)

      logical(1):: onlist

      integer(4),parameter:: maxvl = 20
      character(1):: inpmod(maxvl)
      integer(4):: iwork1(maxvl),iwork2(maxvl)
      real(8):: rval(maxvl)
      integer(4):: ival(maxvl)
      character(80):: cval(maxvl)

      integer(4):: efcol,numerr,icol,nelel,iele,numval,numint,numrea,  &
                   numcha,icon,ierr
      real(8):: rtmp,celbuf,ellbuf

      space = ' '

!********************************************************
!     <<<  SET DEFAULT VALUES  >>>

      ier = 0 ; numerr = 0
      onend = .false. ; onlist = .false.

      intele(1) = 0
      reaele(2) = fxcpul
      intele(3) = 0 ; intele(4) = 10
      reaele(5) = 2.d0
      intele(6) = -1 ; intele(7) = 1 ; intele(8) = 4 ; intele(9) = 3
      reaele(10) = 300.d0 ; reaele(11) = 40.d0
      intele(12) = 2 ; intele(13) = 584287
      reaele(14) = 300.d0
      intele(15) = 1
      chaele(16:20) = ' '
      intele(21) = 3 ; intele(22) = iy15m
      reaele(23) = fycutl ; reaele(24) = fydiel ; reaele(25) = fuctmp
      reaele(26) = fuwpsc ; reaele(27) = fuwdsc ; reaele(28) = fuwdhc
      reaele(29) = furcap ; reaele(30) = fukcap
      intele(31) = iufcap ; intele(32) = iyfshk
      reaele(33) = fustol
      intele(34) = iuslop
      reaele(35:37) = fxcell(1:3)
      intele(38) = ixfbou + 1 ; intele(39) = ixcbou + 1
      reaele(40:42) = fxcbou(1:3)
      intele(43:45) = -1
      reaele(46:48) = fxellp(1:3)
      intele(49) = -1
      chaele(50) = ' ' ; chaele(51) = ' '
      reaele(52) = 10.d0
      intele(53:55) = 1
      chaele(56) = ' '
      intele(57) = 2 ; intele(58) = 0 ; intele(59) = 0
      chaele(60) = ' '
      intele(61) = 1 ; intele(62) = 0 ; intele(63) = 2
      reaele(64) = fukcap_pro ; reaele(65) = furcap_pro
      intele(66) = 5
      intele(67) = capshp
      reaele(68:73) = celwal(1:6) ; reaele(74) = CAPbuff
      reaele(75:76) = 0.d0
      intele(77) = 0
      reaele(78) = CAPbuff
      reaele(79:80) = 0.d0
      intele(81) = iynvar
      intele(82) = 2 ; intele(83) = 0 ; intele(84) = 1
      intele(85) = 1

!********************************************************
!     <<<  READ CONTROL DATA FOR MOLECULAR DYNAMICS  >>>

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

        ! 1-1) Search elements
        call sryelm(iprint,numele,elemnt,line,icol,mxelel,nelel,       &
                    posele,iwork,ier)
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

        ! 1-2) Read data of each element
        do iele = 1,nelel
          inpmod(2*iele-1) = "C"
          if ( dataty(posele(iele)) .eq. "D" ) then
            inpmod(2*iele) = "C"
          else
            inpmod(2*iele) = dataty(posele(iele))
          endif
        enddo
        numval = nelel*2
        call rdfree(line,icol,maxvl,numval,inpmod,iwork1,iwork2,rval,&
                    ival,cval,ier)
        if ( ier .ne. 0 ) then
          write(iprint,*)'ERROR> INPMD '
          write(iprint,*)'    DATA TYPE ERROR ',ier
          write(iprint,'(a80)')line
          return
        endif
        
        ! 1-3) Store data of each element
        numint = 0 ; numrea = 0 ; numcha = 0
        do iele = 1,nelel
          if ( dataty(posele(iele)) .eq. "I" ) then
            numint = numint + 1 ; numcha = numcha + 1
            intele(posele(iele)) = ival(numint)
          elseif ( dataty(posele(iele)) .eq. "R" ) then
            numrea = numrea + 1 ; numcha = numcha + 1
            reaele(posele(iele)) = rval(numrea)
          elseif ( dataty(posele(iele)).eq."C" ) then
            numcha = numcha + 2
            do icon = 1,numspe(posele(iele))
              if ( cval(numcha)(1:4) .eq. datcon(icon,posele(iele)) )  &
                intele(posele(iele)) = icon
            enddo
          elseif ( dataty(posele(iele)).eq."D" ) then
            numcha = numcha + 2
            chaele(posele(iele)) = cval(numcha)
          endif
        enddo

      enddo

200   if ( ierr .ne. 0 ) onend = .true.

!***********************************************************************
!     <<<  STORE ELEMENT DATA TO VARIABLE OF MOLECULAR DYNAMICS  >>>

      idfrst = intele(1) ; fxcpul = reaele(2) ; idlopl = intele(3)
      idupdl = intele(4) ; fddelt = reaele(5) ; idmonl = intele(6)
      idflog = intele(7) ; idfstp = intele(8) ; idfvsl = intele(9)
      fdtemp = reaele(10) ; fdtaut = reaele(11) ; idfvel = intele(12)
      idseed = intele(13) ; fdtmps = reaele(14) ; idfbst = intele(15)
      cdnrti = chaele(16) ; cdnrto = chaele(17) ; cdncrd = chaele(18)
      cdnvel = chaele(19) ; cdnerg = chaele(20) ; idfacc = intele(21)-2
      iy15m = intele(22) ; fycutl = reaele(23) ; fydiel = reaele(24)
      fuctmp = reaele(25) ; fuwpsc = reaele(26) ; fuwdsc = reaele(27)
      fuwdhc = reaele(28) ; furcap = reaele(29) ; fukcap = reaele(30)
      iufcap = intele(31) ; iyfshk = intele(32) ; fustol = reaele(33)
      iuslop = intele(34)
      fxcell(1:3) = reaele(35:37)
      ixfbou = intele(38)-1 ; ixcbou = intele(39)-1
      fxcbou(1:3) = reaele(40:42)
      idmonc = intele(43) ; idmonv = intele(44) ; idmone = intele(45)
      fxellp(1:3) = reaele(46:48)
      idmnts = intele(49) ; cdmpf = chaele(50) ; cdntrj = chaele(51)
      if ( ixfbou .eq. 3 ) then
        fxellp(1:3) = reaele(52) ; ixfbou = 2
      endif
      idfacv = intele(53)-2 ; idface = intele(54)-2
      idfact = intele(55)-2 ; cdntte = chaele(56)
      idfacm = intele(57)-2 ; idhtlp = intele(58) ; idmylp = intele(59)
      cdntce = chaele(60) ; idftce = intele(61)-2 ; idmoce = intele(62)
      itmeth = intele(63) ; fukcap_pro = reaele(64)
      if ( reaele(65) .gt. 0.d0 ) CAPbuff = reaele(65)
      furcap_pro = furcap - CAPbuff
      furcap2 = furcap*furcap ; furcap_pro2 = furcap_pro*furcap_pro
      rcir2 = furcap*furcap - fxcbou(3)*fxcbou(3)
      if ( rcir2 .gt. 0 ) then
        rcir = sqrt(rcir2) ; rcrate = fxcbou(3) / rcir
      else
        rcir = 0.d0 ; rcrate = 0.d0
      endif
      plane_pro = furcap - furcap_pro
      rcir_pro2 = furcap_pro*furcap_pro - (fxcbou(3)-plane_pro)**2
      if ( rcir_pro2 .gt. 0 ) then
        rcir_pro = sqrt(rcir_pro2)
        rcrate_pro = (fxcbou(3)-plane_pro)/rcir_pro
      else
        rcir_pro = 0.d0 ; rcrate_pro = 0.d0
      endif
      idupcm = intele(66)
      if ( iy15m.eq.1 .and. mod(idupdl,idupcm).ne.0 ) then
         write(iprint,*) 'ERROR> INPMD'
         write(iprint,*) '       MOD(UPDATE,CMMUPD) should be 0'
         ierr = -1 ; stop
      endif
      capshp = intele(67) ; celwal(1:6) = reaele(68:73)
      celbuf = reaele(74) ; psetmp = reaele(75) ; stapst = reaele(76)
      htpslo = intele(77) ; ellbuf = reaele(78)
      if ( cluster_method .eq. "LMD" ) then
        if ( reaele(79) .ne. 0.d0 ) then
          lambda = reaele(79)
        elseif ( reaele(80) .ne. 0.d0 ) then
          lambda = sqrt(reaele(80))
        endif
      elseif ( cluster_method .eq. "LMD2" ) then
        if ( reaele(79) .ne. 0.d0 ) then
          lambda = reaele(79)*reaele(79)
        elseif ( reaele(80) .ne. 0.d0 ) then
          lambda = reaele(80)
        endif
      endif
      outatm = intele(81)
      high_para = intele(82)
      i_vec = intele(83)
      nstpcn = intele(84)
      fdebug = intele(85)

      ! Logical Unit
      iducrd = 0 ; iduvel = 0 ; iduerg = 0 ; idutrj = 0 ; idutce = 0
      idutte = 0
      if ( cdncrd .ne. space ) iducrd = 42
      if ( cdnvel .ne. space ) iduvel = 43
      if ( cdnerg .ne. space ) iduerg = 44
      if ( cdntrj .ne. space ) idutrj = 50
      if ( cdntce .ne. space ) idutce = 50
      if ( cdntte .ne. space ) idutte = 59

!***********************************************************************
!     <<<  OUTPUT CONTROL DATA OF MOLECULAR DYNAMICS  >>>

      write(iprint,*)'      PARAMETERS FOR MOLECULAR DYNAMICS '
      write(iprint,*)' '
      write(iprint,*)'      1) GENERAL PARAMETERS  '
      write(iprint,*)'         RESTART         : ',datcon(idfrst,1)
      write(iprint,*)'         CPU LIMIT (S)   : ',fxcpul
      write(iprint,*)'         LOOP LIMIT      : ',idlopl
      write(iprint,*)'         TIME STEP(FSEC) : ',fddelt
      write(iprint,*)'         UPDATE          : ',idupdl
      if ( iy15m .eq. 1 ) then
        write(iprint,*)'         CMM UPDATE      : ',idupcm
        if ( mod(idupdl,idupcm) .ne. 0 ) then
          write(iprint,*)"ERROR> INPMD"
          write(iprint,*)"  CMM UPDATE setting is wrong"
          stop
        endif
      endif
      write(iprint,*)'         MONITORING LOG  : ',idmonl
      write(iprint,*)'         LOG FORMAT      : ',datcon(idflog,7)
      write(iprint,*)'         STOP CENTER     : ',datcon(idfstp,8)
      if ( idfstp .ne. 4 )                                             &
        write(iprint,'(a,i0)')'            #CHN STOPED   : 1-',nstpcn
      write(iprint,*)'         METHOD OF MD    : ',datcon(idfvsl,9)

      select case ( idfvsl )
      case (1)
        write(iprint,*)'           MD WITH CONSTANT TEMPERATURE'
        write(iprint,*)'           (BERENDSEN METHOD)'
        write(iprint,*)'           TEMPERATURE (K) : ',fdtemp
        write(iprint,*)'           RELAXATION(FSEC): ',fdtaut
      case (2)
        write(iprint,*)'           MICRO-CANONICAL MD'
        write(iprint,*)'           TEMPERATURE (K) : ',fdtemp
      case (3)
        write(iprint,*)'           CANONICAL MD'
        write(iprint,*)'           ( HOOVER & EVANS METHOD)'
        write(iprint,*)'           TEMPERATURE (K) : ',fdtemp
      case (4)
        write(iprint,*)'           GENERALIZED ENSEMBLE (McMD or ALSD)'
        write(iprint,*)                                                &
              '           (POTENTIAL SCALING BY LINEAR FUNCTION)'
        write(iprint,*)'           TEMPERATURE (K) : ',fdtemp
      case default
        write(iprint,'("ERROR> INVALID METHOD-ID : ",i3)')idfvsl
        ier = -1 ; return
      end select

      if ( idfrst .eq. 1 ) then
        write(iprint,*)'         INITIAL VELOCITY: ',datcon(idfvel,12)
        if ( idfvel .eq. 2 ) then
          write(iprint,*)'         TEMPERATURE (K) : ',fdtmps
          write(iprint,*)'         RANDOM SEED     : ',idseed
        endif
      endif
      if ( high_para.eq.-1 .and. nomp.eq.1 ) then
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

      write(iprint,*)'         CALC. RMSD      : ',datcon(idfbst,15)

      write(iprint,*)'      2) FILE I/O '
      if ( idfrst .eq. 2 ) then
        write(iprint,*)'         RESTART FILE (INPUT)  '
        write(iprint,'(12x,a4,10x,a2,a)')"NAME",": ",trim(cdnrti)
        write(iprint,*)' '
      else
        write(iprint,*)'         NO RESTART FILE INPUT'
        write(iprint,*)' '
      endif

      if ( cdnrto .ne. space ) then
        write(iprint,*)'         RESTART FILE (OUTPUT) '
        write(iprint,'(12x,a4,10x,a2,a)')"NAME",": ",trim(cdnrto)
        write(iprint,*)' '
      else
        write(iprint,*)'         NO RESTART FILE OUTPUT'
        write(iprint,*)' '
      endif

      if ( idfacc.ge.0 .and. iducrd.gt.0 .and. idmonc.gt.0 ) then
        write(iprint,*)'         COORDINATE FILE (OUTPUT) '
        write(iprint,'(12x,a4,10x,a2,a)')"NAME",": ",trim(cdncrd)
        write(iprint,*)'           LOGICAL UNIT  : ',iducrd
        write(iprint,*)'           OUTPUT F.eq.  : ',idmonc
        write(iprint,*)'           FORMAT OF FILE: ',datcon(idfacc+2,21)
        if ( outatm .ne. iynvar ) then
          write(iprint,*)'             OUTPUT ATOM : 1 - ',outatm
        endif
        write(iprint,*)' '
      else
        write(iprint,*)'         NO COORDINATE FILE OUTPUT'
        write(iprint,*)' '
        idfacc = -1
      endif

      if ( idfacv.ge.0 .and. iduvel.gt.0 .and. idmonv.gt.0 ) then
        write(iprint,*)'         VELOCITY FILE (OUTPUT) '
        write(iprint,'(12x,a4,10x,a2,a)')"NAME",": ",trim(cdnvel)
        write(iprint,*)'           LOGICAL UNIT  : ',iduvel
        write(iprint,*)'           OUTPUT F.eq.  : ',idmonv
        write(iprint,*)'           FORMAT OF FILE: ',datcon(idfacv+2,53)
        if ( outatm .ne. iynvar ) then
          write(iprint,*)'             OUTPUT ATOM : 1 - ',outatm
        endif
        write(iprint,*)' '
      else
        write(iprint,*)'         NO VELOCITY FILE OUTPUT'
        write(iprint,*)' '
        idfacv = -1
      endif

      if ( idface.ge.0 .and. iduerg.gt.0 .and. idmone.gt.0 ) then
        write(iprint,*)'         ENERGY FILE (OUTPUT) '
        write(iprint,'(12x,a4,10x,a2,a)')"NAME",": ",trim(cdnerg)
        write(iprint,*)'           LOGICAL UNIT  : ',iduerg
        write(iprint,*)'           OUTPUT F.eq.  : ',idmone
        write(iprint,*)'           FORMAT OF FILE: ',datcon(idface+2,54)
        write(iprint,*)' '
      else
        write(iprint,*)'         NO ENERGY FILE OUTPUT'
        write(iprint,*)' '
        idface = -1
      endif

      if ( idfact.ge.0 .and. idutrj.gt.0 .and. idmnts.gt.0 ) then
        write(iprint,*)'         STRUCTURAL TRAJECTORY (OUTPUT) '
        write(iprint,'(12x,a4,10x,a2,a)')"NAME",": ",trim(cdntrj)
        write(iprint,*)'           OUTPUT F.eq.  : ',idmnts
        write(iprint,*)'           FORMAT OF FILE: ',datcon(idfact+2,55)
        write(iprint,*)' '
      else
        write(iprint,*)'         NO MON. TRAJECTORY FILE OUTPUT'
        write(iprint,*)' '
        idfact = -1
      endif

      if ( idfacm.ge.0 .and. idutte.gt.0 ) then
        write(iprint,*)'         TOTAL ENERGY AT EVERY STEP (OUTPUT) '
        write(iprint,'(12x,a4,10x,a2,a)')"NAME",": ",trim(cdntte)
        write(iprint,*)'           OUTPUT F.eq.  : EVERY STEP'
        write(iprint,*)'           FORMAT OF FILE: ',datcon(idfacm+2,57)
        write(iprint,*) ' '
      else
        write(iprint,*)                                                &
              '         NO OUTPUT OF TOTAL ENERGY AT EVERY STEP'
        write(iprint,*) ' '
        idfacm = -1
      endif

      if ( idftce.ge.0 .and. idutce.gt.0 .and. idmoce.gt.0 ) then
        write(iprint,*)'         ENERGIES AMONG INDIVIDUAL CLUSTERS',  &
                       ' (OUTPUT) '
        write(iprint,*)'           UNIT          : ',idutce
        if ( cdntce .ne. space ) then
          write(iprint,'(12x,a4,10x,a2,a)')"NAME",": ",trim(cdntce)
        endif
        write(iprint,*)'           OUTPUT F.eq.  : ',idmoce
        write(iprint,*)'           FORMAT OF FILE: ',                  &
                                                  datcon(idftce+2,61)
        write(iprint,*)' '
      else
        write(iprint,*)'         NO OUTPUT OF INDIVIDUAL ENERGY TERMS'
        write(iprint,*)' '
        idftce = -1
      endif

      write(iprint,*)'      5) PARAMETERS FOR ENERGY CALCULATION '

      if ( iy15m .eq. 1 ) then
        write(iprint,*)                                                &
              '         CALCULATE NON-BONDED INTERACTIONS BY CMM'
        write(iprint,*)                                                &
              "         VDW CUT-OFF DIST. FOR NO-EXACT CALC. CELL :",  &
              fycutl
        if ( ixfbou .eq. 1 ) then
          write(iprint,*)" You cannot choose CMM under periodic "//    &
          "boundary condition." ; stop
        endif
      elseif ( iy15m .eq. 3 ) then
        write(iprint,*)" SORRY !!"
        write(iprint,*)" Residue-based CUT-OFF method is removed. "
        stop
      else
        write(iprint,*)'         CUT-OFF LOGIC   : ',datcon(iy15m,22)
        write(iprint,*)'         CUT-OFF DIST.   : ',fycutl
      endif

      write(iprint,*)'         DIELECTRIC VAL. : ',fydiel

      if ( fuctmp .eq. 0.d0 ) then
        fuctmp = 1.d0 / boltzk
      elseif ( maxval(iyeflg(11:14)) .eq. 1 ) then
        write(iprint,*)'         TEMP. FOR REST. : ',fuctmp
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
          write(iprint,*)'         RADIUS OF CAP   : ',furcap
          write(iprint,*)'         FORCE CONS. CAP : ',fukcap
          write(iprint,*)                                              &
            '         RADIUS OF CAP for specified atoms: ',furcap_pro
          write(iprint,*)                                              &
            '         SQUARE OF THE RADIUS             : ',furcap_pro2
          write(iprint,*)                                              &
            '         FORCE CONS. CAP for specified atoms: ',fukcap_pro
          write(iprint,*)'         FUNC. TYPE(CAP) : ',datcon(iufcap,31)
          write(iprint,*)'         CAP CENTER (X)  : ',fxcbou(1)
          write(iprint,*)'         CAP CENTER (Y)  : ',fxcbou(2)
          write(iprint,*)'         CAP CENTER (Z)  : ',fxcbou(3)

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
          write(iprint,*)'         FUNC. TYPE(CAP) : ',datcon(iufcap,31)

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
          write(iprint,*)'         FUNC. TYPE(CAP) : ',datcon(iufcap,31)
          fxellp(1:3) = 1.d0 / (fxellp(1:3)*fxellp(1:3))
          fxellp_pro(1:3) = 1.d0 / (fxellp_pro(1:3)*fxellp_pro(1:3))
        endif
      endif
      if ( iyeflg(15) .eq. 1 )                                         &
        write(iprint,*)'         COEFFICIENT FOR REPULSION : ',Krep

      if ( idfvsl .eq. 3 ) then
        if ( cluster_method .eq. "LMD" ) then
          if ( idfrst .ne. 2 ) then
            write(iprint,*)
            write(iprint,*)'              LAMBDA     : ',lambda
            write(iprint,*)
          endif

        elseif ( cluster_method .eq. "LMD2" ) then
          if ( idfrst .ne. 2 ) then
            write(iprint,*)
            write(iprint,*)'        LAMBDA SQUARE    : ',lambda
            write(iprint,*)
          endif

        elseif ( psetmp .ne. 0.d0 ) then
          write(iprint,*)'         PSEUDO TEMP. is used'
          write(iprint,*)'           PSEUDO TEMP. : ',psetmp,' (K)'
          if ( stapst .ne. 0.d0 ) then
            write(iprint,*)'         START PSEUDO TEMP. : ',stapst,    &
                           " (K)"
            write(iprint,*)"         HEAT LOOP for PSEUDO T : ",htpslo
          endif
        endif
      endif

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

      elseif ( ixfbou .eq. 2 ) then
        write(iprint,*)'         APPLY ELLIPSOIDAL BOUNDARY'
        write(iprint,*)'          CENTER(X)       : ',fxcbou(1)
        write(iprint,*)'          CENTER(Y)       : ',fxcbou(2)
        write(iprint,*)'          CENTER(Z)       : ',fxcbou(3)
        write(iprint,*)'          ELLIPSOID AXIS;A: ',fxellp(1)
        write(iprint,*)'          ELLIPSOID AXIS;B: ',fxellp(2)
        write(iprint,*)'          ELLIPSOID AXIS;C: ',fxellp(3)
      endif

      if ( iyfshk .ne. 0 ) then
        write(iprint,*)'         APPLY SHAKE CONSTRAINTS'
        write(iprint,*)'          SHAKE METHOD  : ',datcon(iyfshk,32)
        write(iprint,*)'          LOOP LIMIT    : ',iuslop
        write(iprint,*)'          CRI. CONV.    : ',fustol

        if ( iutshk.eq.0 .or. iugshk.eq.0 ) then
          write(iprint,*)' '
          write(iprint,*)'ERROR> INPMD'
          write(iprint,*)'          SHAKE parameters are defined,' 
          write(iprint,*)'          but no SHAKE setting file is read. '
          write(iprint,*)' '
          ier = -1 ; return
        endif
      else
        if ( iutshk.ne.0 .and. iugshk.ne.0 ) then
          write(iprint,*)' '
          write(iprint,*)'ERROR> INPMD'
          write(iprint,*)'          SHAKE setting file is read,' 
          write(iprint,*)                                              &
                '          but no SHAKE parameters are defined. '
          write(iprint,*)' '
          ier = -1 ; return
        endif
      endif
     

      if ( onlist ) then
        onlist = .FALSE. ; goto 100
      endif

      ! Parameters for constraint MD method
      if ( idfvsl .ge. 3 ) then
        write(iprint,*)'         HEATLOOP STEPS  : ',idhtlp
        if ( idfvsl .ge. 4 ) then
          write(iprint,*)'         DUMMY CANONICAL LOOPS : ',idmylp
        endif
      endif

      write(iprint,*)' '

      ! READ PARAMETERS FOR GE
      if ( idfvsl .eq. 4 ) then
        call flopen(1,cdmpf,10,'NULL',80,ier)
        if ( ier .ne. 0 ) then
          write(iprint,*)'ERROR> INPMD'
          write(iprint,*)'   FILE OPEN ERROR OF PARAMETER FOR GE.'
          write(iprint,*)'   NAME= ',cdmpf
          ier = -1 ; return
        endif
        call inmcmd(iprint,1,ier)
        close(1)
      endif

      if ( ixfbou .eq. 1 ) then
        invcel(1:3) = 1.d0 / fxcell(1:3)
      else
        invcel(1:3) = 0.d0
      endif

      ! Flag for debug
      select case (fdebug)
        case(2)
          write(iprint,*)' '
          write(iprint,*)'INFORMATION> FLAG for DEBUG = YES'
          write(iprint,*)'  FORCE & ENERGY ERRORS ARE CHECKED'
          write(iprint,*)'  FOR THE FIRST STEP'
          write(iprint,*)' '
        case(3)
          write(iprint,*)' '
          write(iprint,*)'INFORMATION> FLAG for DEBUG = MAKE'
          write(iprint,*)'  FORCE & ENERGY FILE FOR DEBUG'
          write(iprint,*)'  (1stp_force.bin)'
          write(iprint,*)'  IS MADE FOR DEBUG'
          write(iprint,*)' '
      end select

!***********************

      return
      end subroutine inpmd
