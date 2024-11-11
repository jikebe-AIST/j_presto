
      subroutine md(iread,iprint,ier,onend)

!*******************************************************************
!
!     Molecular dynamics verlet method (leap-frog)
!
!*******************************************************************

      use COMBAS ; use COMERG ; use COMMIS ; use COMMOL
      use COMCMMC ; use CALC_TIME
      use COMCMM,only: nfrg,nfrag,icls,natfr,iatfr,cluster_method,    &
                       l_temp_control
      !$ use omp_lib

      implicit none

      ! Logical unit number for input & output
        integer(4),intent(in):: iread,iprint
      ! Condition code (0: No-error, negative: Error, 
      !                 positive: warning (CPU TIME IS EXCEEDED))
        integer(4),intent(inout):: ier 
      ! Flag for end of control data (.true. is End of control data is
      !                                      detected)
        logical(1),intent(inout):: onend

!       FDVELO     R*8  (MAXATM,3)   : TEMPORAY ARRAY FOR VELOCITY
!                       IN VERLET METHOD
!                        VEL(T) = ( VEL(T-0.5DT) + VEL(T+0.5DT) ) / 2
!                       BACK UP VEL(T-0.5DT) TO FDVELO
      real(8):: velo(3,ixnatm)
      real(8):: fdepot(maxene)
      real(8):: fdcvcv,fddelt,fdtemp,fdtaut,fdtmps,fdtims
      character(80):: cdnrti,cdnrto,cdncrd,cdnvel,cdnerg,cdntrj,cdntte
      character(80):: cdntce

      integer(4):: idfrst,idlopl,idupdl,idmonc,idmonv,idmone
      integer(4):: idmonl,idflog,idfstp,idfvsl,idfvel,idseed,idfbst
      integer(4):: iducrd,iduvel,iduerg,idutte
      integer(4):: idnfre,ivar,idfacc,idfacv,idface,idfact
      integer(4):: idfacm,idutrj,idmnts,idhtlp,idmylp,idutce,idftce
      integer(4):: idmoce,idupcm,idloop

      integer(4):: i,j,n,iatm,ifrg,imol
!!      integer(4):: nm1en,nm2en,nm2st,nm2nd
      integer(4):: icon,ilflag,iii(1)
      integer(8):: secst,secen,secst1,secen1
      real(8):: time1,time2,epotcc(maxene,nfrg) 
      real(8),allocatable:: ratio(:)
      logical(4):: flg

!*************************************************************

!     <<<  INITILIZATION  >>>
      call system_clock(secst,cr)
      ier = 0 ; fdcvcv = 4.184D-4 ; idloop = 0 ; fdtims = 0.d0
      write(iprint,*)'INFORMATION> MD (V3.0) '
      write(iprint,*)'     MOLECULAR DYNAMICS  (VERLET METHOD) '
      write(iprint,*)' '

      ! Make inverse of mass
      allocate(ifxmass(ixnatm))
      !$OMP parallel default (none)                                  & !
      !$OMP private(i)                                               & !
      !$OMP shared(ifxmass,ixnatm,fxmass)
      !$OMP do schedule (static)
      do i = 1,ixnatm
        ifxmass(i) = 1.d0 / fxmass(i)
      enddo
      !$OMP end do
      !$OMP end parallel

!     <<<  READ CONTROL DATA  >>>
      call inpmd(iread,ier,iprint,onend,idfrst,idlopl,idupdl,   &
                 idmonl,idflog,fddelt,idfstp,idfvsl,fdtemp,     &
                 fdtaut,idfvel,idseed,fdtmps,idfbst,     &
                 iducrd,iduvel,iduerg,idutrj,idfacc,idfacv,idface,     &
                 idfact,cdnrti,cdnrto,cdncrd,cdnvel,cdnerg,     &
                 cdntrj,idmonc,idmonv,idmone,idmnts,idutte,cdntte,     &
                 idfacm,idhtlp,idmylp,idutce,cdntce,idftce,idmoce,     &
                 idupcm)
      call incel(iprint)
      if ( ier .ne. 0 ) goto 9999

!     <<<  OPEN DATA FILES OF MOLECULAR DYNAMICS TRAJECTORY  >>>
      call opndat(iducrd,iduvel,iduerg,idutrj,idutte,idutce,cdncrd,    &
                  cdnvel,cdnerg,cdntrj,cdntte,cdntce,idfacc,idfacv,    &
                  idface,idfact,idfacm,idftce,ier)
      if ( ier .ne. 0 ) goto 9000
 
!     <<<  CALCULATE THE CENTER OF MASS OF A PROTEIN >>>
      if ( ixfbou.ne.1 .and. ixcbou.eq.1 ) then
        call calcm(ixcend(nstpcn),fxmass(1:ixnatm),fxcbou,ier)
        if ( ier .ne. 0 ) then
          write(iprint,*)'ERROR> MD '
          write(iprint,*)'  MOLECULAR WEIGHT OF FIRST CHAIN IS ZERO '
          goto 9999
        endif
        write(iprint,*)'INFORMATION> MD '
        write(iprint,'(a,i0,a)')                                       &
          '   CENTER OF BOUNDARY IS MASS CENTER OF 1-',nstpcn,' CHAIN '
        write(iprint,*)'             X        : ',fxcbou(1)
        write(iprint,*)'             Y        : ',fxcbou(2)
        write(iprint,*)'             Z        : ',fxcbou(3)
        write(iprint,*)' '
      endif
 
!     <<<  INSPECT THE CHAINS INSIDE OR OUTSIDE THE ELLIPSOID  >>>
      if ( ixfbou .eq. 2 ) then
        call ellprp(maxatm,iprint,ixnchn,ixcend,fxmass,fxcbou,fxellp,  &
                    ier)
        if ( ier .ne. 0 ) goto 9000
      endif
 
!     <<<  SET NUMBER OF FREEDOM  >>>
      idnfre = iynvar*3 - iutshk
      if ( idfstp.eq.1 .or. idfstp.eq.2 ) then
        idnfre = idnfre - 3
      elseif ( idfstp .eq. 3 ) then
        idnfre = idnfre - 6
      endif
      if ( l_temp_control .and. cluster_method(1:3).eq."LMD" .and.     &
           idfvsl.eq.4 ) idnfre = idnfre + 1
      write(iprint,*)'INFORMATION> MD '
      write(iprint,*)'    NUMBER OF FREEDOM OF THIS SYSTEM IS ',idnfre
      if ( idnfre .le. 0 ) then
        ier = -1
        write(iprint,*)'ERROR> MD '
        write(iprint,*)'  NUMBER OF FREEDOM IS TOO SMALL '
        goto 9999
      endif

      ! Make ixamol2 for cluster temperature
      allocate(ixamol2(ixnatm))
      ixamol2(1:ixnatm) = (ixamol(1:ixnatm)-1) * nfrag + 1
      do i = 1,ixnatm
        ixamol2(i) = ixamol2(i) + icls(i) - 1
      enddo
      do i = ixamol(ixnatm)*nfrag,1,-1
        flg = .false.
        do j = 1,ixnatm
          if ( ixamol2(j) .eq. i ) then
            flg = .true. ; exit
          endif
        enddo
        if ( .not. flg ) then
          do j = 1,ixnatm
            if ( ixamol2(j) .ge. i ) ixamol2(j) = ixamol2(j) - 1
          enddo
        endif
      enddo

!     <<<  Get some quantities for each sub-system >>>
!         The following quantities:
!            Number of atoms each sub-system ---> natfr.
      nfrmx = maxval(ixamol2(1:ixnatm))
      allocate(natfr(nfrmx),nfrmol(nfrmx))
      natfr(1:nfrmx) = 0

      do i = 1,ixnatm
        imol = ixamol2(i)
        natfr(imol) = natfr(imol) + 1
      enddo
      allocate(iatfr(maxval(natfr(1:nfrmx)),nfrmx))
      natfr(1:nfrmx) = 0
      do i = 1,ixnatm
        imol = ixamol2(i)
        n = natfr(imol) + 1
        iatfr(n,imol) = i
        natfr(imol) = n
      enddo

!     Number of degrees of freedom for each sub-system.
!     Just giving 3 degrees of freedom to each atom.
      nfrmol(1:nfrmx) = natfr(1:nfrmx)*3.d0

!     Reduce the number of degrees of freedom with SHAKE information.
      do i = 1,iugshk
        imol = ixamol2(iuashk(1,i))
        icon = ( iuhshk(i)+1)*iuhshk(i)*0.5
        nfrmol(imol) = nfrmol(imol) - dble(icon)
      enddo

!     Removal of center of mass motion.
!     Note: the removal is done only for the first chain (fragment).
      select case (idfstp)
      case (1:3)
        do i = 1,ixnatm
          if ( ixamol(1) .eq. ixamol(i) ) then
            iatm = i
          else
            exit
          endif
        enddo
        ifrg = maxval(ixamol2(1:iatm))
        if ( ifrg .ge. 2 ) then
          allocate(ratio(nfrmx))
          do i = 1,ifrg
            ratio(i) = nfrmol(i) / sum(nfrmol(1:ifrg))
          enddo
          if ( idfstp.eq.1 .or. idfstp.eq.2 ) then
            do i = 1,ifrg
              nfrmol(i) = nfrmol(i) - 3.d0*ratio(i)
            enddo
          elseif ( idfstp .eq. 3 ) then
            do i = 1,ifrg
              nfrmol(i) = nfrmol(i) - 6.d0*ratio(i)
            enddo
          endif
          deallocate(ratio)
        else
          if ( idfstp.eq.1 .or. idfstp.eq.2 ) then
            nfrmol(1) = nfrmol(1) - 3.d0
          elseif ( idfstp .eq. 3 ) then
            nfrmol(1) = nfrmol(1) - 6.d0
          endif
        endif
      end select

      ! add a degree of freedom for lambda
      if ( l_temp_control .and. cluster_method(1:3).eq."LMD" .and.     &
           idfvsl.eq.4 ) then
        iii = maxloc(nfrmol(1:nfrmx))
        ifg_lambda = iii(1)
        nfrmol(ifg_lambda) = nfrmol(ifg_lambda) + 1.d0
      else
        ifg_lambda = 0
      endif

!     Atom ordering check.
!     itmeth is temperature control parameter.
!        itmeth=1 ---> whole T control.
!        itmeth=2 ---> independent T control.
      if ( itmeth .eq. 2 ) then
        do ivar = 1 , iynvar
          iatm = iylvar(ivar)
          if ( iatm .ne. ivar ) then
            print*,' ??????????????????????????????????????'
            print*,' ? The independent T control does not ?'
            print*,' ? work in this culculation, because  ?'
            print*,' ? N(atom) =/= N(movable atom).       ?'
            print*,' ??????????????????????????????????????'
            stop
          endif
        enddo
      endif

!     <<<  INITIAL PROCESS OF MOLECULAR DYNAMICS  >>>
      call inimd(iprint,ier,idloop,idhtlp,fdcvcv,fdtims,fddelt,idnfre, &
                 velo,fdepot,idfstp,idfvsl,fdtemp,fdtaut,idfvel,idseed,&
                 fdtmps,idfrst,cdnrti)
        if ( ier .ne. 0 ) goto 9999
 
!     <<<  RUN MOLECULAR DYNAMICS  >>>
      if ( idlopl .gt. 0 ) then
        call system_clock(secst1)
        call runmd(iprint,ier,idloop,idlopl,idhtlp,idupdl,idmonc,      &
                   idmonv,idmone,idmonl,idflog,fdcvcv,fdtims,fddelt,   &
                   idnfre,velo,fdepot,fdtmps,idfstp,idfvsl,fdtemp,     &
                   fdtaut,idfbst,iducrd,iduvel,iduerg,idutrj,idfacc,   &
                   idfacv,idface,idfact,idmnts,idutte,idutce,idftce,   &
                   idmoce,idupcm,cdnrto)
        call system_clock(secen1)
        if ( ier .lt. 0 ) goto 9999
      else
        secst1 = 0.d0 ; secen1 = 0.d0
      endif

!     <<<  FINAL PROCESS OF MOLECULAR DYNAMICS  >>>
      call finmd(iprint,ier,idloop,idhtlp,idflog,fdcvcv,fdtims,fddelt, &
                 idnfre,velo,fdepot,idfstp,idfvsl,fdtemp,fdtaut,fdtmps,&
                 cdnrto)

!     <<<  CLOSE DATA FILES OF MOLECULAR DYNAMICS TRAJECTORY  >>>
      call clsdat(iducrd,iduvel,iduerg,idutrj,idutte,ier)
      if ( ier .ne. 0 ) goto 9000

9999  call system_clock(secen)
      if ( secen .gt. secst ) then
        time1 = dble(secen-secst)/dble(cr)
      else
        time1 = dble(cm-secst+secen)/dble(cr)
      endif
      if ( secen1 .gt. secst1 ) then
        time2 = dble(secen1-secst1)/dble(cr)
      else
        time2 = dble(cm-secst1+secen1)/dble(cr)
      endif
      write(iprint,*)""
      write(iprint,*)"INFORMATION> MD "
      write(iprint,'(12x,a22,f15.5)')'TOTAL CPU TIME (S)  : ',time1
      write(iprint,'(12x,a22,f15.5)')'  MD LOOP           : ',time2
      if ( run_md_tim .ne. 0 ) call outcpu2(iprint,'run_md',           &
        run_md_tim-verlet_tim-celprp_tim,cr,time2)
      if ( verlet_tim .ne. 0 ) call outcpu2(iprint,'verlet',           &
        verlet_tim-dnergy_tim-shake_tim,cr,time2)
      if ( dnergy_tim .ne. 0 ) call outcpu2(iprint,'dnergy',           &
        dnergy_tim,cr,time2)
      if ( shake_tim .ne. 0 ) call outcpu2(iprint,'shake ',            &
        shake_tim,cr,time2)
      if ( celprp_tim .ne. 0 ) call outcpu2(iprint,'celprp',           &
        celprp_tim,cr,time2)
      write(iprint,'(12x,a22,f15.5)')'  OTHERS            : ',         &
        time1-time2
      write(iprint,'(12x,a22,f15.5)')'  SPEED (ns/day)    : ',         &
        dble(idlopl)*fddelt*0.0864d0/time2

!****************************************************

      return

9000  write(iprint,*)'ERROR> MD '
      if ( ier .eq. -1 ) then
        write(iprint,*)'  COORDINATE DATA FILE (BINARY) OPEN ERROR'
        write(iprint,*)'  OR CLOSE ERROR '
      elseif ( ier .eq. -2 ) then
        write(iprint,*)'  VELOCITY DATA FILE OPEN ERROR'
        write(iprint,*)'  OR CLOSE ERROR '
      elseif ( ier .eq. -3 ) then
        write(iprint,*)'  ENERGY DATA FILE OPEN ERROR'
        write(iprint,*)'  OR CLOSE ERROR '
      elseif ( ier .eq. -4 ) then
        write(iprint,*)'  TRAJECTORY DATA FILE OPEN ERROR'
        write(iprint,*)'  OR CLOSE ERROR '
      elseif ( ier .eq. -5 ) then
        write(iprint,*)'  TOTAL ENERGY AT EVERY STEP OPEN ERROR'
        write(iprint,*)'  OR CLOSE ERROR '
      elseif ( ier .eq. -6 ) then
        write(iprint,*)'  SOME CHAINS IN THE INITIAL STRUCTURE'
        write(iprint,*)'  ARE OUTSIDE OF THE ELLIPSOIDE OR SPHERE '
      endif
      goto 9999
 
      end subroutine md


!=================================================================================


      subroutine opndat(iducrd,iduvel,iduerg,idutrj,idutte,idutce,     &
                        cdncrd,cdnvel,cdnerg,cdntrj,cdntte,cdntce,     &
                        idfacc,idfacv,idface,idfact,idfacm,idftce,ier)

      implicit none

      integer(4):: iducrd,iduvel,iduerg,idutrj,idutte,idutce
      character(80):: cdncrd,cdnvel,cdnerg,cdntrj,cdntte,cdntce
      integer(4):: idfacc,idfacv,idface,idfact,idfacm,idftce,ier

!*************************************************

      ier = 0
      if ( cdncrd .ne. " " ) then
        if ( idfacc .gt. 0 ) then
          call flopen(iducrd,cdncrd,22,'NULL',0,ier)
        elseif ( idfacc .eq. 0 ) then
          call flopen(iducrd,cdncrd,12,'NULL',0,ier)
        endif
        if ( ier .ne. 0 ) then
          ier = -1 ; return
        endif
      endif
      if ( cdnvel .ne. " " ) then
        if ( idfacv .gt. 0 ) then
          call flopen(iduvel,cdnvel,22,'NULL',0,ier)
        elseif ( idfacv .eq. 0 ) then
          call flopen(iduvel,cdnvel,12,'NULL',0,ier)
        endif
        if ( ier .ne. 0 ) then
          ier = -2 ; return
        endif
      endif
      if ( cdnerg .ne. " " ) then
        if ( idface .gt. 0 ) then
          call flopen(iduerg,cdnerg,22,'NULL',0,ier)
        elseif ( idface .eq. 0 ) then
          call flopen(iduerg,cdnerg,12,'NULL',0,ier)
        endif
        if ( ier .ne. 0 ) then
          ier = -3 ; return
        endif
      endif
      if ( idutrj .gt. 0 ) then
        if ( idfact .gt. 0 ) then
          call flopen(idutrj,cdntrj,22,'NULL',0,ier)
        elseif ( idfact .eq. 0 ) then
          call flopen(idutrj,cdntrj,12,'NULL',0,ier)
        endif
        if ( ier .ne. 0 ) then
          ier = -4 ; return
        endif
      endif
      if ( idfacm .gt. 0 ) then
        call flopen(idutte,cdntte,22,'NULL',0,ier)
      elseif ( idfacm .eq. 0 ) then
        call flopen(idutte,cdntte,12,'NULL',0,ier)
      endif
      if ( ier .ne. 0 ) then
        ier = -5 ; return
      endif
      if ( idftce .gt. 0 ) then
        call flopen(idutce,cdntce,22,'NULL',0,ier)
      elseif ( idftce .eq. 0 ) then
        call flopen(idutce,cdntce,12,'NULL',0,ier)
      endif
      if ( ier .ne. 0 ) then
        ier = -6 ; return
      endif

!***************************************

      return
      end subroutine opndat
 

!====================================================================


      subroutine inimd(iprint,ier,idloop,idhtlp,fdcvcv,fdtims,fddelt,  &
                       idnfre,velo,fdepot,idfstp,idfvsl,fdtemp,fdtaut, &
                       idfvel,idseed,fdtmps,idfrst,cdnrti)

!*******************************************************************
!
!     INITIAL PROCESS OF MOLECULAR DYNAMICS
!       A) INITIAL START  (INISTA)
!          SET VELOCITY  (SETVEL)
!          ONE STEP MD   (STAVER)
!       B) RESTART        (INIRST)
!          INPUT RESTART FILE (INPRST)
!          CHECK ENERGY CONTINUITY (VERLET)
!
!*******************************************************************
 
      use COMBAS ; use COMERG ; use COMCMMC,only: vel

      implicit none

      integer(4),intent(in):: iprint,idnfre
      integer(4),intent(inout):: ier,idloop,idhtlp,idfstp,idfvsl,      &
                                 idfvel,idseed
      real(8),intent(inout):: fdcvcv,fdtims,fddelt,fdtemp,fdtaut,fdtmps
      real(8),intent(inout):: velo(3,ixnatm),fdepot(maxene)
      integer(4),intent(in):: idfrst
      character(*),intent(in):: cdnrti

!***************************************************

      ier = 0
!     <<<  INITIAL START  >>>
      if ( idfrst .eq. 1 ) then
        call inista(iprint,ier,idloop,fdcvcv,fdtims,fddelt,idnfre,     &
                    fdepot,idfstp,idfvsl,fdtemp,fdtaut,idfvel,idseed,  &
                    fdtmps)
        velo(1:3,1:iynvar) = vel(1:3,1:iynvar)
!     <<<  RESTART  >>>
      elseif ( idfrst .eq. 2 ) then
        call inirst(iprint,ier,idloop,idhtlp,fdcvcv,fdtims,fddelt,     &
                    idnfre,velo,fdepot,idfstp,idfvsl,fdtemp,fdtaut,    &
                    fdtmps,cdnrti)
      endif

!**********************************

      return
      end subroutine inimd
 

!==================================================================


      subroutine runmd(iprint,ier,idloop,idlopl,idhtlp,idupdl,idmonc,  &
                       idmonv,idmone,idmonl,idflog,fdcvcv,fdtims,      &
                       fddelt,idnfre,velo,fdepot,fdtmps,idfstp,idfvsl, &
                       fdtemp,fdtaut,idfbst,iducrd,iduvel,iduerg,      &
                       idutrj,idfacc,idfacv,idface,idfact,idmnts,      &
                       idutte,idutce,idftce,idmoce,idupcm,cdnrto)

!*******************************************************************
!
!      RUNNING MOLECULAR DYNAMICS VERLET METHOD (LEAP FROG)
!
!*******************************************************************
 
      use COMBAS ; use COMERG ; use COMMIS ; use COMCMM
      use COMPSC ; use PHYCNS ; use COMCMMC,only: cord,vel
      use CALC_TIME

      implicit none

      integer(4),intent(in):: iprint,idlopl,idupdl,idutte,idutce,      &
                              idupcm,idnfre,iducrd,iduvel,iduerg,idutrj
      integer(4),intent(inout):: ier,idloop,idhtlp,idmonc,idmonv
      integer(4),intent(inout):: idmone,idmonl,idflog,idftce,idmoce
      integer(4),intent(inout):: idfstp,idfvsl,idfbst,idfacc,idfacv
      integer(4),intent(inout):: idface,idfact,idmnts
      real(8),intent(in):: fddelt
      real(8),intent(inout):: fdtims,fdcvcv,fdtemp,fdtaut,fdtmps
      real(8),intent(inout):: velo(3,ixnatm),fdepot(maxene)
      character(*),intent(in):: cdnrto

      real(8):: fdvelb(3,ixnatm),epotcc(maxene,nfrg),fdcrdm(3,ixnatm)

      integer(4):: isrlop,iloop,i,j
      integer(8):: sec
      real(8):: curtim,lamda,lambda_vt
      logical(4):: mnt_flg
      real(8):: Ttime

!************************************************************

      call system_clock(run_md_tim1)
      ier = 0 ; lamda = 0.d0 ; isrlop = idloop
      write(iprint,*)' '
      write(iprint,*)'INFORMATION> MD '
      write(iprint,*)'             MOLECULAR DYNAMICS LOOP START '

      ! STEP INTEGRAL
      do iloop = 1,idlopl
        idloop = idloop + 1
        curtim = fdtims + fddelt*dble(iloop)
        if ( (idmonc.gt.0 .and. mod(idloop,idmonc).eq.0) .or.          &
             (idmonv.gt.0 .and. mod(idloop,idmonv).eq.0) .or.          &
             (idmone.gt.0 .and. mod(idloop,idmone).eq.0) .or.          &
             (idmnts.gt.0 .and. mod(idloop,idmnts).eq.0) .or.          &
             (idmonl.gt.0 .and. mod(idloop,idmonl).eq.0) ) then
          mnt_flg = .true.
          call system_clock(tim1)
          !$OMP parallel default (none)                              & !
          !$OMP private(i)                                           & !
          !$OMP shared(fdcrdm,cord,ixnatm)
          !$OMP do schedule (static)
          do i = 1,ixnatm
            fdcrdm(:,i) = cord(:,i)
          enddo
          !$OMP end do
          !$OMP end parallel
        else
          mnt_flg = .false.
        endif

        ! Temp.const. or microcanonical
        call system_clock(verlet_tim1)
        if ( idfvsl .le. 2 ) then
          call verlet(iprint,ier,fdcvcv,fddelt,idnfre,fdvelb,fdepot,   &
               epotcc,idfstp,idfvsl,fdtemp,fdtaut,lamda)

        ! canonical (3: cano, 4: McMD)
        else
          call verlethe(iprint,ier,fdcvcv,fddelt,fdtmps,idnfre,fdvelb, &
                        fdepot,epotcc,idfstp,idfvsl,fdtemp,idloop,     &
                        idhtlp,lambda_vt)
        endif
        call system_clock(verlet_tim2)
        verlet_tim = verlet_tim + verlet_tim2 - verlet_tim1
        if ( verlet_tim2 .le. verlet_tim1 ) verlet_tim = verlet_tim + cm

        if ( ier .ne. 0 ) then
          call outpdb(89,6,"ERROR.pdb",maxatm,ixnatm,cxatmn,cxresn,    &
                      absres,fxmass,fxchrg,10,ixtitl,cxtitl,ier)
          j = 0
          do i = 1,maxene
            if ( iyeflg(i) .eq. 1 ) then
              write(iprint,'("  ",a6," : ",e15.7,$)')cyenam(i),fdepot(i)
              j = j + 1
            endif
            if ( j .eq. 3 ) then
              write(iprint,*) ; j = 0
            endif
          enddo
          if ( j .ne. 0 ) write(iprint,*)
          fdtims = fdtims + fddelt*dble(idloop-isrlop) ; return
        endif

!*********************

        ! Output the total energy
        if ( idutte .ne. 0 ) then
          if ( cluster_method(1:3) .eq. "LMD" ) then
            write(idutte,'(2(e15.7,x),$)')lambda,epotcc(1,1)
            if ( nfrg .eq. 1 ) then
              write(idutte,*)
            else
              write(idutte,'(2(e15.7,x))')epotcc(1,2),                 &
                                          sum(epotcc(1,1:nfrg))
            endif

          else
            select case(nfrg)
              case (1)
                write(idutte,'(e15.7)')epotcc(1,1)
              case (2)
                write(idutte,'(2(e15.7,x))')epotcc(1,1),               &
                                            sum(epotcc(1,1:2))
              case (3)
                write(idutte,'(3(e15.7,x))')epotcc(1,1),epotcc(1,2),   &
                                            sum(epotcc(1,1:3))
            end select

          endif
        endif

        ! Output the individual energy terms
        if ( idmoce.gt.0 .and. mod(idloop,idmoce).eq.0 .and.           &
             idutce.gt.0 ) then
          if ( idftce .eq. 0 ) then
            do i = 1,maxene
              write(idutce,'(e15.7,x,e15.7)')epotcc(i,1:nfrg)
            enddo
          elseif ( idftce .eq. 1 ) then
            write(idutce)(sngl(epotcc(i,1:nfrg)),i=1,maxene)
          elseif ( idftce .eq. 2 ) then
            write(idutce)(epotcc(i,1:nfrg),i=1,maxene)
          endif
        endif

        ! Output trajectory and monitoring
        if ( mnt_flg ) then
          !$OMP parallel default (none)                              & !
          !$OMP private(i)                                           & !
          !$OMP shared(velo,iynvar,fdvelb,vel)
          !$OMP do schedule (static)
          do i = 1,iynvar
            velo(:,i) = 0.5d0*(fdvelb(:,i)+vel(:,i))
          enddo
          !$OMP end do
          !$OMP end parallel
          lambda_vt = 0.5d0*(lambda_vt+lambda_v)
          call monmd(iprint,ier,idloop,idflog,idmonc,idmonv,idmone,    &
                     idmnts,idmonl,idnfre,fdcrdm,velo,fdepot,lamda,    &
                     curtim,idfbst,iducrd,iduvel,iduerg,idutrj,idfacc, &
                     idfacv,idface,idfact,idfstp,epotcc(1,1:nfrg),     &
                     lambda_vt)
          if ( ier .ne. 0 ) then
            fdtims = fdtims + fddelt*dble(idloop-isrlop) ; return
          endif

        elseif ( (idmonc.gt.0 .and. mod(idloop+1,idmonc).eq.0) .or.    &
                 (idmonv.gt.0 .and. mod(idloop+1,idmonv).eq.0) .or.    &
                 (idmone.gt.0 .and. mod(idloop+1,idmone).eq.0) .or.    &
                 (idmnts.gt.0 .and. mod(idloop+1,idmnts).eq.0) .or.    &
                 (idmonl.gt.0 .and. mod(idloop+1,idmonl).eq.0) ) then
          !$OMP parallel default (none)                              & !
          !$OMP private(i)                                           & !
          !$OMP shared(velo,iynvar,vel)
          !$OMP do schedule (static)
          do i = 1,iynvar
            velo(:,i) = vel(:,i)
          enddo
          !$OMP end do
          !$OMP end parallel
        endif

        ! Update interaction table and check cpu time
        if ( mod(idloop,idupdl).eq.0 .and. idloop.ne.0 ) then
          call system_clock(celprp_tim1)
          call celprp(iprint,idflog,ier)
          call system_clock(celprp_tim2)
          celprp_tim = celprp_tim + celprp_tim2 - celprp_tim1
          if ( celprp_tim2.le.celprp_tim1 ) celprp_tim=celprp_tim+cm
          if ( ier .ne. 0 ) then
            fdtims = fdtims + fddelt*dble(idloop-isrlop) ; return
          endif
          call system_clock(sec)
          if ( sec .gt. fxcpus ) then
            Ttime = dble(sec-fxcpus)/dble(cr)
          else
            Ttime = dble(cm-fxcpus+sec)/dble(cr)
          endif
          if ( Ttime .ge. fxcpul ) then
            write(iprint,*)" "
            write(iprint,*)"WARNING> MD"
            write(iprint,*)"   CPU TIME IS EXCEEDED"
            fdtims = fdtims + fddelt*dble(idloop-isrlop)
            ier = 1 ; return
          endif
        endif
            
        !! CMM far-field update marker
        if ( iy15m.eq.1 .and. mod(idloop,idupcm).eq.0 ) then
          itwin = 1
        else
          itwin = 0
        endif

        if ( cdnrto.ne." " .and.                                     &
            (idmonl.gt.0 .and. mod(idloop,idmonl).eq.0) ) then
          call midmd(iprint,ier,idloop,idhtlp,fdcvcv,curtim,fddelt,    &
                     idnfre,velo,fdepot,idfstp,idfvsl,fdtemp,fdtaut,   &
                     fdtmps,cdnrto,idlopl)
          if ( ier .ne. 0 ) return
        endif

      enddo

      fdtims = fdtims + fddelt*dble(idloop-isrlop)
      call system_clock(run_md_tim2)
      run_md_tim = run_md_tim + run_md_tim2 - run_md_tim1
      if ( run_md_tim2 .le. run_md_tim1 ) run_md_tim = run_md_tim + cm

!***********************************

      return
      end subroutine runmd


!====================================================================


      subroutine midmd(iprint,ier,idloop,idhtlp,fdcvcv,fdtims,fddelt,  &
                       idnfre,velo,fdepot,idfstp,idfvsl,fdtemp,fdtaut, &
                       fdtmps,cdnrto,idlopl)

!*******************************************************************
!
!     MIDIUM PROCESS OF MOLECULAR DYNAMICS 
!       for writing temporary restart file 
!
!*******************************************************************
 
      use COMBAS ; use COMERG ; use COMMIS ; use COMCMM
      use COMPSC ; use COMCMMC ; use CALC_TIME

      implicit none

      integer(4),intent(in):: iprint,idnfre
      integer(4),intent(in):: idlopl
      integer(4),intent(inout):: ier,idloop,idhtlp
      integer(4),intent(inout):: idfstp,idfvsl
      real(8),intent(in):: fddelt
      real(8),intent(inout):: fdcvcv,fdtims,fdtemp,fdtaut,fdtmps
      real(8),intent(inout):: velo(3,ixnatm),fdepot(maxene)
      character(*),intent(in):: cdnrto

      integer(4):: i
      real(8):: et,ek,temp,lamda,lambda_vt,lambda_old,lambda_vold
      real(8):: fdvelb(3,ixnatm),fdv05p(3,ixnatm),epotcc(maxene,nfrg)
      real(8):: fdcrdm(3,ixnatm)
      character(80):: ctmp
      character(80),save:: Pctmp = " "

!*****************************************************

      ier = 0 ; lamda = 0.d0        

      ! if idloop .eq. idlopl then...
      if ( idloop.eq.idlopl .and. Pctmp.ne." " ) then
        call system("rm "//trim(Pctmp)//" >& /dev/null")
        return
      endif

      ! Back up velocity & coordinates
      lambda_old = lambda
      lambda_vold = lambda_v
      !$OMP parallel default (none)                                  & !
      !$OMP private(i)                                               & !
      !$OMP shared(fdv05p,vel,iylvar,iynvar,fdcrdm,cord,ixnatm)
      !$OMP do schedule (static)
      do i = 1,iynvar
        fdv05p(:,i) = vel(:,i)
      enddo
      !$OMP end do
      !$OMP do schedule (static)
      do i = 1,ixnatm
        fdcrdm(:,i) = cord(:,i)
      enddo
      !$OMP end do
      !$OMP end parallel

!     3) ONE STEP MD

      ! Temp.const. or microcanonical
      call system_clock(verlet_tim1)
      if ( idfvsl .le. 2 ) then
        call verlet(iprint,ier,fdcvcv,fddelt,idnfre,fdvelb,fdepot,     &
                    epotcc,idfstp,idfvsl,fdtemp,fdtaut,lamda)

      ! canonical ( 3: cano, 4: McMD)
      else
        call verlethe(iprint,ier,fdcvcv,fddelt,fdtmps,idnfre,fdvelb,   &
             fdepot,epotcc,idfstp,idfvsl,fdtemp,idloop,idhtlp,lambda_vt)
      endif
      call system_clock(verlet_tim2)
      verlet_tim = verlet_tim + verlet_tim2 - verlet_tim1

      if ( ier .ne. 0 ) return
      !$OMP parallel default (none)                                  & !
      !$OMP private(i)                                               & !
      !$OMP shared(iynvar,cord,fdcrdm,velo,fdvelb,vel,ixnatm)
      !$OMP do schedule (static)
      do i = 1,iynvar
        velo(:,i) = 0.5d0*(fdvelb(:,i)+vel(:,i))
      enddo
      !$OMP end do
      !$OMP do schedule (static)
      do i = 1,ixnatm
        cord(:,i) = fdcrdm(:,i)
      enddo
      !$OMP end do
      !$OMP end parallel
      lambda_vt = 0.5d0*(lambda_vt+lambda_v)

      call caltmp(iynvar,iylvar,velo,fxmass(1:ixnatm),idnfre,temp,ek,  &
                  lambda_vt,lambda_m)
      et = ek + fdepot(1) 

      ! Output restart file
      write(ctmp,'(i0)')idloop
      ctmp = trim(cdnrto)//"_"//trim(ctmp)//"step"
      write(iprint,*)" ( Temp. restart file "//trim(ctmp)//" output )"
      if ( ixfbou .eq. 1 ) call PBcord
      call outrst(iprint,ctmp,ixnatm,iynvar,cord,fdv05p,et,ek,  &
                  fdepot(1),idloop,fdtims,lambda,lambda_vold,ier)
      if ( ier .ne. 0 ) return
      ! Delete past temp_restart file
      if ( Pctmp .ne. " " ) then
        call system("rm "//trim(Pctmp)//" >& /dev/null")
      endif
      Pctmp = ctmp

      !$OMP parallel default (none)                                  & !
      !$OMP private(i)                                               & !
      !$OMP shared(vel,fdv05p,iynvar)
      !$OMP do schedule (static)
      do i = 1,iynvar
        vel(:,i) = fdv05p(:,i)
      enddo
      !$OMP end do
      !$OMP end parallel
      lambda = lambda_old
      lambda_v = lambda_vold

!****************************************

      return
      end subroutine midmd


!======================================================================== 


      subroutine finmd(iprint,ier,idloop,idhtlp,idflog,fdcvcv,fdtims,  &
                       fddelt,idnfre,velo,fdepot,idfstp,idfvsl,fdtemp, &
                       fdtaut,fdtmps,cdnrto)

!*******************************************************************
!
!     FINAL PROCESS OF MOLECULAR DYNAMICS
!
!*******************************************************************
 
      use COMBAS ; use COMERG ; use COMMIS ; use COMCMM
      use COMPSC ; use COMCMMC

      implicit none

      integer(4),intent(in):: iprint,idnfre
      integer(4),intent(inout):: ier,idloop,idflog,idhtlp,idfstp,idfvsl
      real(8),intent(in):: fddelt
      real(8),intent(inout):: fdcvcv,fdtims,fdtemp,fdtaut,fdtmps
      real(8),intent(inout):: velo(3,ixnatm),fdepot(maxene)
      character(*),intent(in):: cdnrto

      integer(4):: ivar,iatm,locf,i,j
      real(8):: et,ek,temp,rmsf,lamda,maxf,lambda_vt,lambda_old,       &
                lambda_vold
      real(8):: fdvelb(3,ixnatm),fdv05p(3,ixnatm),epotcc(maxene,nfrg)
      real(8):: fdcrdm(3,ixnatm)

!*****************************************************

      ier = 0 ; lamda = 0.d0        

      ! Back up velocity
      !$OMP parallel default (none)                                  & !
      !$OMP private(i)                                               & !
      !$OMP shared(fdv05p,iynvar,vel)
      !$OMP do schedule (static)
      do i = 1,iynvar
        fdv05p(:,i) = vel(:,i)
      enddo
      !$OMP end do
      !$OMP end parallel
      lambda_old = lambda
      lambda_vold = lambda_v

      ! Updata interaction table
      call celprp(iprint,idflog,ier)
      if ( ier .ne. 0 ) return

!     3) ONE STEP MD

      !$OMP parallel default (none)                                  & !
      !$OMP private(i)                                               & !
      !$OMP shared(fdcrdm,ixnatm,cord)
      !$OMP do schedule (static)
      do i = 1,ixnatm
        fdcrdm(:,i) = cord(:,i)
      enddo
      !$OMP end do
      !$OMP end parallel

      ! Temp.const. or microcanonical
      if ( idfvsl .le. 2 ) then
        call verlet(iprint,ier,fdcvcv,fddelt,idnfre,fdvelb,fdepot,     &
                    epotcc,idfstp,idfvsl,fdtemp,fdtaut,lamda)

      ! canonical (3: cano, 4: McMD)
      else
        call verlethe(iprint,ier,fdcvcv,fddelt,fdtmps,idnfre,fdvelb,   &
             fdepot,epotcc,idfstp,idfvsl,fdtemp,idloop,idhtlp,lambda_vt)
      endif

      if ( ier .ne. 0 ) return
      !$OMP parallel default (none)                                  & !
      !$OMP private(i)                                               & !
      !$OMP shared(iynvar,ixnatm,cord,fdcrdm,velo,fdvelb,vel)
      !$OMP do schedule (static)
      do i = 1,ixnatm
        cord(:,i) = fdcrdm(:,i)
      enddo
      !$OMP end do
      !$OMP do schedule (static)
      do i = 1,iynvar
        velo(:,i) = 0.5d0*(fdvelb(:,i)+vel(:,i))
      enddo
      !$OMP end do
      !$OMP end parallel
      lambda = lambda_old
      lambda_vt = 0.5d0*(lambda_vt+lambda_v)

      call caltmp(iynvar,iylvar,velo,fxmass(1:ixnatm),idnfre,temp,ek,  &
                  lambda_vt,lambda_m)
      et = ek + fdepot(1)

      ! Output restart file
      if ( cdnrto .ne. " " ) then
        if ( ixfbou .eq. 1 ) call PBcord
        call outrst(iprint,cdnrto,ixnatm,iynvar,cord,fdv05p,et, &
                    ek,fdepot(1),idloop,fdtims,lambda,lambda_vold,ier)
        if ( ier .ne. 0 ) return
        ! Delete temporary restart file
        call system("rm "//trim(cdnrto)//"_*step >& /dev/null")
      endif

      ! Output log
      call calrsf(iynvar,iylvar,grad(1:3,1:ixnatm,nfrg),rmsf,maxf,locf)
      write(iprint,*)""
      write(iprint,'(a17)')" INFORMATION> MD "
      write(iprint,'(a27)')"             FINAL PROCESS "
      write(iprint,'(a25,i15,$)')  " MD LOOP NUMBER        : ",idloop+1
      write(iprint,'(a25,f15.5)')" TIME (PSEC)           : ",          &
                                 (fdtims+fddelt)*0.001d0
      write(iprint,'(a25,e15.7,$)')" TOTAL      (KCAL/MOL) : ",et
      write(iprint,'(a25,e15.7)')  " POTENTIAL  (KCAL/MOL) : ",fdepot(1)
      write(iprint,'(a25,f15.7,$)')" TEMPERATURE (K)       : ",temp
      write(iprint,'(a25,e15.7)')  " KINETIC    (KCAL/MOL) : ",ek
      write(iprint,'(a25,e15.7,$)')" R.M.S.F. (KCAL/MOL*A) : ",rmsf
      write(iprint,'(a25,e15.7)')  " SCALING FACTOR        : ",lamda
      write(iprint,'(a25,e15.7,$)')" MAX F.   (KCAL/MOL*A) : ",maxf
      write(iprint,'(a7,i5,x,a8,i5,x,a8)')" CHAIN ",                   &
           ixachn(locf),cxresn(locf),ixares(locf),cxatmn(locf)

      j = 0
      do i = 1,maxene
        if ( iyeflg(i) .eq. 1 ) then
          write(iprint,'("  ",a6," : ",e15.7,$)')cyenam(i),fdepot(i)
          j = j + 1
        endif
        if ( j .eq. 3 ) then
          write(iprint,*) ; j = 0
        endif
      enddo
      if ( j .ne. 0 ) write(iprint,*)

!****************************************

      return
      end subroutine finmd


!======================================================================== 


      subroutine inista(iprint,ier,idloop,fdcvcv,fdtims,fddelt,idnfre, &
                        fdepot,idfstp,idfvsl,fdtemp,fdtaut,idfvel,     &
                        idseed,fdtmps)
 
      use COMBAS ; use COMERG ; use PHYCNS
      use COMCMMC ; use COMCMM,only: nfrg,lambda_v,lambda_m

      implicit none

      integer(4),intent(in):: iprint,idnfre
      integer(4),intent(inout):: ier,idloop,idfstp,idfvsl,idfvel,idseed
      real(8),intent(in):: fddelt
      real(8),intent(inout):: fdcvcv,fdtims,fdtemp,fdtaut,fdtmps
      real(8),intent(inout):: fdepot(maxene)

      integer(4):: i,iene,locf
      real(8):: ek,ek2,rmsf,lamda,temp,temp1,temp2,maxf
      real(8):: random(3,maxatm)

!*******************************************

      ier = 0 ; idloop = 0 ; fdtims = 0.d0 ; lamda = 0.d0
      ek = 0.d0 ; temp = 0.d0
      ! Set initial velocity
      call setvel(idfvel,gascst,idseed,fdtmps,iynvar,iylvar,maxatm,    &
                  random,fxmass,ier)
      if ( ier .ne. 0 ) then
        write(iprint,*)' '
        write(iprint,*)'ERROR> MD '
        write(iprint,*)'    GENERATION OF RANDOM NUMBER ERROR '
        return
      endif

      ! Remove translational velocity and rotational velocity from 
      !   initial velocity
      if ( idfvel .eq. 2 ) then
        call caltmp(iynvar,iylvar,vel,fxmass(1:ixnatm),idnfre,temp1,ek,&
                    lambda_v,lambda_m)
        call remotv(ixnatm,fxmass(1:ixnatm),cord,ier,3)
        if ( ier .ne. 0 ) then
          write(iprint,*)' '
          write(iprint,*)'ERROR> MD '
          write(iprint,*)'   ERROR IN REMOVING OUTER MOVEMENT '
          write(iprint,*)'   MOMENTA MATRIX MAY BE SINGULAR '
          return
        endif
        call caltmp(iynvar,iylvar,vel,fxmass(1:ixnatm),idnfre,temp2,ek,&
                    lambda_v,lambda_m)
        rmsf = sqrt(fdtmps/temp2)
        !$OMP parallel default (none)                                & !
        !$OMP private(i)                                             & !
        !$OMP shared(vel,iynvar,rmsf)
        !$OMP do schedule (static)
        do i = 1,iynvar
          vel(:,i) = rmsf * vel(:,i)
        enddo
        !$OMP end do
        !$OMP end parallel
        call caltmp(iynvar,iylvar,vel,fxmass(1:ixnatm),idnfre,temp,ek, &
                    lambda_v,lambda_m)
        write(iprint,*)""
        write(iprint,*)"INFORMATION> MD"
        write(iprint,*)"    SET INITIAL VELOCITY "
        write(iprint,*)"        TEMPERATURE  (K) "
        write(iprint,*)                                                &
              "        1) INITIAL SETTING BY RANDOM NUMBER : ",temp1
        write(iprint,*)                                                &
              "        2) AFTER REMOVING OUTER MOTION      : ",temp2
        write(iprint,*)                                                &
              "        3) RE-SCALING                       : ",temp
      endif
 
      ! Make interaction table
      if ( iyeflg(8).eq.1  .or. iyeflg(9).eq.1 .or. iyeflg(10).eq.1 ) then
        call celprp(iprint,2,ier)
        if ( ier .ne. 0 ) return
      endif

      ! Initial process of MD (verlet)
      call staver(iprint,ier,fdcvcv,fddelt,idnfre,fdepot,idfstp,idfvsl,&
                  fdtemp,fdtaut,lamda)
      call calrsf(iynvar,iylvar,grad(1:3,1:ixnatm,nfrg),rmsf,maxf,locf)
      write(iprint,*)' '
      write(Iprint,*)'INFORMATION> MD '
      write(iprint,*)'   INITIAL POTENTIAL ENERGY '
      do iene = 1,maxene
        if ( iyeflg(iene) .eq. 1 ) then
          write(iprint,*)'     ' ,cyenam(iene)," : ",fdepot(iene)
        endif
      enddo
      write(iprint,*)'     RMSF.  : ',rmsf
      write(iprint,*)'    MAX F.  : ',maxf
      write(iprint,*)"    CHAIN   : ",ixachn(locf)," ",cxresn(locf),   &
                                      ixares(locf)," ",cxatmn(locf)

      call caltmp(iynvar,iylvar,vel,fxmass(1:ixnatm),idnfre,temp2,ek2, &
                  lambda_v,lambda_m)
      write(iprint,*)""
      write(iprint,'(a)')" INFORMATION> MD "
      write(iprint,'(a)')"             INITIAL PROCESS "
      write(iprint,'(a23,f10.5,a24,e15.7)')"        TIME (PSEC)  : ",  &
           fdtims," TOTAL     (KCAL/MOL) : ",fdepot(1)+ek
      write(iprint,'(33x,a24,e15.7)')" POTENTIAL (KCAL/MOL) : ",       &
                                     fdepot(1)
      write(iprint,'(33x,a24,e15.7)')" KINETIC   (KCAL/MOL) : ",ek
      write(iprint,'(33x,a24,e15.7)')" TEMPERATURE (K)      : ",temp
      write(iprint,'(a23,f10.5,a24,e15.7)')"        TIME (PSEC)  : ",  &
           fddelt*0.0005d0," KINETIC   (KCAL/MOL) : ",ek2
      write(iprint,'(33x,a24,e15.7)')" TEMPERATURE (K)      : ",temp2

!***************************************

      return
      end subroutine inista
 

!====================================================================
 

      subroutine inirst(iprint,ier,idloop,idhtlp,fdcvcv,fdtims,fddelt, &
                        idnfre,velo,fdepot,idfstp,idfvsl,fdtemp,fdtaut,&
                        fdtmps,cdnrti)

      use COMBAS ; use COMERG ; use COMMIS ; use COMCMM
      use COMPSC ; use PHYCNS ; use COMCMMC

      implicit none

      integer(4),intent(in):: iprint,idnfre
      integer(4),intent(inout):: ier,idloop,idhtlp,idfstp
      integer(4),intent(inout):: idfvsl
      real(8),intent(in):: fddelt
      real(8),intent(inout):: fdcvcv,fdtims,fdtmps,fdtemp,fdtaut
      real(8),intent(inout):: velo(3,ixnatm),fdepot(maxene) 
      character(*),intent(in):: cdnrti

      real(8):: fdcrdm(3,ixnatm)
      real(8):: et,ek,ep,et2,ek2,temp,rmsf,lamda,maxf,lambda_vt,       &
                lambda_old,lambda_vold
      real(8):: fdvelb(3,ixnatm),epotcc(maxene,nfrg)

      integer(4):: i,j,locf
 
!********************************************

      lamda = 0.d0
      ! Input restart file
      call inprst(iprint,cdnrti,ixnatm,iynvar,et,ek,ep,idloop,fdtims,  &
                  idfvsl,ier)
      if ( ier .ne. 0 ) return

      ! Check energy continuity
      if ( iyeflg(8).eq.1  .or. iyeflg(9).eq.1 .or. iyeflg(10).eq.1 ) then
        call celprp(iprint,2,ier)
        if ( ier .ne. 0 ) return
      endif
      !$OMP parallel default (none)                                  & !
      !$OMP private(i)                                               & !
      !$OMP shared(fdcrdm,ixnatm,cord)
      !$OMP do schedule (static)
      do i = 1,ixnatm
        fdcrdm(:,i) = cord(:,i)
      enddo
      !$OMP end do
      !$OMP end parallel
      lambda_old = lambda
      lambda_vold = lambda_v

      ! Temp.const. or microcanonical
      if ( idfvsl .le. 2 ) then
        call verlet(iprint,ier,fdcvcv,fddelt,idnfre,fdvelb,fdepot,     &
                    epotcc,idfstp,idfvsl,fdtemp,fdtaut,lamda)

      ! canonical (3: cano, 4: McMD)
      else
        call verlethe(iprint,ier,fdcvcv,fddelt,fdtmps,idnfre,fdvelb,   &
             fdepot,epotcc,idfstp,idfvsl,fdtemp,idloop,idhtlp,lambda_vt)
      endif

      !$OMP parallel default (none)                                  & !
      !$OMP private(i)                                               & !
      !$OMP shared(velo,iynvar,fdvelb,vel)
      !$OMP do schedule (static)
      do i = 1,iynvar
        velo(:,i) = 0.5d0*(fdvelb(:,i)+vel(:,i))
      enddo
      !$OMP end do
      !$OMP end parallel
      lambda_vt = 0.5d0*(lambda_vt+lambda_v)
      call calrsf(iynvar,iylvar,grad(1:3,1:ixnatm,nfrg),rmsf,maxf,locf)
      call caltmp(iynvar,iylvar,velo,fxmass(1:ixnatm),idnfre,temp,ek2, &
                  lambda_vt,lambda_m)
      et2 = ek2 + fdepot(1)

 
      write(IPRINT,*) ' '
      if ( abs(et2-et) .gt. 0.0001d0 ) then
        write(iprint,*)'WARNING> MD '
        write(iprint,*)'         ENERGY CONTINUITY IS NOT ',           &
                        'SATISFIED AT RESTART PROCESS.'
        write(iprint,*)'     1) START ENERGY '
        write(iprint,*)'         TOTAL        : ',et2
        write(iprint,*)'         KINETIC      : ',ek2
        write(iprint,*)'         POTENTIAL    : ',fdepot(1)
        write(iprint,*)'      2) RESTART FILE ENERGY '
        write(iprint,*)'         TOTAL        : ',et
        write(iprint,*)'         KINETIC      : ',ek
        write(iprint,*)'         POTENTIAL    : ',ep
        write(iprint,*)' '
        write(iprint,*)'   ENERGY IS NOT CONTINUOUS, BECAUSE OF ',    &
                        'BY CHANGING SIMULATION PARAMETERS '
        write(iprint,*)'   FOR  EXAMPLE,'
        write(iprint,*)'        A) CHANGE OF TIME STEP '
        write(iprint,*)'        B) CHANGE OF ENERGY PARAMETER OR '
        write(iprint,*)'           DIELECTRIC CONSTANTi'
        write(iprint,*)'        C) CUT-OFF LOGIC CHANGE'
        write(iprint,*)'        D) BOUNDARY CONDITION CHANGE'
        write(iprint,*)'        E) TEMPERATURE CONTROL CHANGE '
        write(iprint,*)'           eg. CONSTANT TEMP. -> CANONICAL'
        write(iprint,*)'        F) INPUT PARAMETER CHANGE FOR '
        write(iprint,*)'           McMD or ALSD'
        write(iprint,*)'        G) LAMBDA CHANGE FOR ALSD etc.'
      else
        write(iprint,*)'INFORMATION> MD '
        write(iprint,*)'             ENERGY CONTINUITY IS SATISFIED ',&
                        'AT RESTART PROCESS.'
        write(iprint,*)' '
      endif

      call calrsf(iynvar,iylvar,grad(1:3,1:ixnatm,nfrg),rmsf,maxf,locf)
      write(iprint,*)""
      write(iprint,'(a)')" INFORMATION> MD "
      write(iprint,'(a)')"             INITIAL PROCESS (RESTART) "
      write(iprint,'(a25,i15,$)')" MD LOOP NUMBER        : ",idloop+1
      write(iprint,'(a25,f15.5)')" TIME (PSEC)           : ",          &
                                 (fdtims+fddelt)*0.001d0
      write(iprint,'(a25,e15.7,$)')" TOTAL      (KCAL/MOL) : ",et
      write(iprint,'(a25,e15.7)')" POTENTIAL  (KCAL/MOL) : ",fdepot(1)
      write(iprint,'(a25,f15.7,$)')" TEMPERATURE (K)       : ",temp
      write(iprint,'(a25,e15.7)')" KINETIC    (KCAL/MOL) : ",ek
      write(iprint,'(a25,e15.7,$)')" R.M.S.F. (KCAL/MOL*A) : ",rmsf
      write(iprint,'(a25,e15.7)')" SCALING FACTOR        : ",lamda
      write(iprint,'(a25,e15.7,$)')" MAX F.   (KCAL/MOL*A) : ",maxf
      write(iprint,'(a25,i5,x,a8,i5,x,a8)')"          CHAIN        : ",&
           ixachn(locf),cxresn(locf),ixares(locf),cxatmn(locf)

      j = 0
      do i = 1,maxene
        if ( iyeflg(i) .eq. 1 ) then
          write(iprint,'("  ",a6," : ",e15.7,$)')cyenam(i),fdepot(i)
          j = j + 1
        endif
        if ( j .eq. 3 ) then
          write(iprint,*) ; j = 0
        endif
      enddo
      if ( j .ne. 0 ) write(iprint,*)

      !$OMP parallel default (none)                                  & !
      !$OMP private(i)                                               & !
      !$OMP shared(iynvar,ixnatm,cord,fdcrdm,vel,fdvelb)
      !$OMP do schedule (static)
      do i = 1,ixnatm
        cord(:,i) = fdcrdm(:,i)
      enddo
      !$OMP end do
      !$OMP do schedule (static)
      do i = 1,iynvar
        vel(:,i) = fdvelb(:,i)
      enddo
      !$OMP end do
      !$OMP end parallel
      lambda = lambda_old
      lambda_v = lambda_vold

!***********************************************

      return
      end subroutine inirst 


!========================================================================


      subroutine clsdat(iducrd,iduvel,iduerg,idutrj,idutte,ier)

      implicit none

      integer(4),intent(in):: iducrd,iduvel,iduerg,idutrj,idutte
      integer(4),intent(inout):: ier 
 
!***********************************************

      ier = 0
      if ( iducrd .gt. 0 ) then
        call flclos(iducrd,10,ier)
        if ( ier .ne. 0 ) then
          ier = -1 ; return
        endif
      endif
      if ( iduvel .gt. 0 ) then
        call flclos(iduvel,10,ier)
        if ( ier .ne. 0 ) then
          ier = -2 ; return
        endif
      endif
      if ( iduerg .gt. 0 ) then
        call flclos(iduerg,10,ier)
        if ( ier .ne. 0 ) then
          ier = -3 ; return
        endif
      endif
      if ( idutrj .gt. 0 ) then
        call flclos(idutrj,10,ier)
        if ( ier .ne. 0 ) then
          ier = -4 ; return
        endif
      endif
      if ( idutte .gt. 0 ) then
        call flclos(idutte,10,ier)
        if ( ier .ne. 0 ) then
          ier = -5 ; return
        endif
      endif

!********************************

      return
      end subroutine clsdat
