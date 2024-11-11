
      subroutine monmd(iprint,ier,idloop,idflog,idmonc,idmonv,idmone,  &
                       idmnts,idmonl,idnfre,fdcrdm,velo,fdepot,lamda,  &
                       curtim,idfbst,iducrd,iduvel,iduerg,idutrj,      &
                       idfacc,idfacv,idface,idfact,idfstp,tepotcc,     &
                       lambda_vt)

!********************************************************************
!
!     WRITE MD TRAJECTORIES INTO FILES
!
!********************************************************************
!
!       3.2 ARGUMENTS
!       IPRINT (I) I*4               : LOGICAL UNIT NUMBER FOR LOG
!       IER    (O) U*4               : CONDITION CODE
!                      0     ; NO ERROR
!                  NEGATIVE  ; ERROR
!       IDLOOP (I) I*4               : CURRENT MD LOOP NUMBER
!       IDFLOG (I) I*4               : FLAG OF FORMAT OF LOG
!       IDMONC (I) I*4               : OUTPUT INTERVAL FOR COORDINATES
!       IDMONV (I) I*4               : OUTPUT INTERVAL FOR VELOCITIES
!       IDMONE (I) I*4               : OUTPUT INTERVAL FOR ENERGIES
!       IDMNTS (I) I*4               : OUTPUT INTERVAL FOR MON.VALUES
!       IDMONL (I) I*4               : OUTPUT INTERVAL FOR LOG
!            WHEN MOD(IDLOOP,IDMONL) = 0, OUTPUT LOG
!       IDNFRE (I) I*4               : NUMBER OF FREEDOM
!       FDCRDM (I) R*8  (MAXATM,3)   : CURRENT COORDINATE AT TIME CURTIM
!       FDVELO (I) R*8  (MAXATM,3)   : CURRENT VELOCITY AT TIME CURTIM
!       FDGRAD (W) R*8  (MAXATM,3)   :
!                                (I) CURRENT GRADIENT AT TIME CURTIM
!                                (O) CURRENT FREE ATOM CORRDINATE
!       FDEPOT (I) R*8  (MAXENE)     : CURRENT POTENTIAL ENERGY AT TIME
!                                      CURTIM
!       LAMDA  (I) R*8               : SCALING FACTOR FOR CONSTANT
!                                      TEMPEARTURE
!       CURTIM (I) R*8               : CURRENT SIMULATION TIME (FSEC)
!       IDFBST (I) I*4               : FLAG FOR CALCULATION RMSD OF
!                                      FIRST CHAIN
!       IDUCRD (I) I*4               : LOGICAL UNIT NUMBER FOR
!                                      OUTPUT COORDINATE DATA
!       IDUVEL (I) I*4               : LOGICAL UNIT NUMBER FOR
!                                      OUTPUT VELOCITY DATA
!       IDUERG (I) I*4               : LOGICAL UNIT NUMBER FOR
!                                      OUTPUT ENERGY DATA
!       IDUTRJ (I) I*4               : LOGICAL UNIT NUMBER FOR
!                                      OUTPUT MON. TRAJECTORY
!       IDFACC (I) I*4               : FORMAT OF OUTPUT COORDINATES
!                                      IDFACC = -1  NO OUTPUT
!                                      IDFACC = 0   ASCII DATA
!                                      IDFACC = 1   REAL*4 DATA
!                                      IDFACC = 2   REAL*8 DATA
!       IDFACV (I) I*4               : FORMAT OF OUTPUT VELOCITIES
!                                      IDFACV = -1  NO OUTPUT
!                                      IDFACV = 0   ASCII DATA
!                                      IDFACV = 1   REAL*4 DATA
!                                      IDFACV = 2   REAL*8 DATA
!       IDFACE (I) I*4               : FORMAT OF OUTPUT ENERGIES
!                                      IDFACE = -1  NO OUTPUT
!                                      IDFACE = 0   ASCII DATA
!                                      IDFACE = 1   REAL*4 DATA
!                                      IDFACE = 2   REAL*8 DATA
!       IDFACT (I) I*4               : FORMAT OF OUTPUT TRAJECTORIES
!                                      IDFACT = -1  NO OUTPUT
!                                      IDFACT = 0   ASCII DATA
!                                      IDFACT = 1   REAL*4 DATA
!                                      IDFACT = 2   REAL*8 DATA
!       IDFSTP (I) I*4               : FLAG TO STOP CENTER OF MASS
!                                      IDFSTP = 0   NO STOP
!                                      IDFSTP = 1   TRANSLATION STOPS
!                                      IDFSTP = 2   ROTATION STOPS
!                                      IDFSTP = 3   BOTH STOP
!
!     4 TEMPORARY VARIABLES
!       RMSF       R*8               : ROOT MEAN SQUARE FORCE
!       RMSD       R*8               : ROOT MEAN SQUARE DEVIATION
!                                      OF FIRST CHAIN
!       SEC        R*8               : CURRENT CPU TIME (SEC) FROM
!                                      JOB IS STARTED
!       TEMPEM     R*8  (ixmolc)     : TEMPEARTURE (K)
!                                      OF EACH MOLECULE
!       EKEM       R*8  (ixmolc)     : KINETIC ENERGY (KCAL/MOL)
!                                      OF EACH MOLECULE
!       NUMFRM     I*4  (ixmolc)     : NUMBER OF FREEDOM OF EACH
!                                      MOLECULE
!       NUMATF     I*4  (ixmolc)     : NUMBER OF FREE ATOMS OF
!                                      EACH MOLECULE
!
!********************************************************************

      use COMBAS ; use COMERG ; use COMMIS ; use PHYCNS
      use COMCMM,only: nfrg,lambda_v,lambda_m
      use COMCMMC,only: grad

      implicit none

      integer(4):: iprint,ier,idloop,idflog,idmonl,idmonc,idmonv,      &
                   idmone,idnfre,idfbst,iducrd,iduvel,iduerg,idfacc,   &
                   idfacv,idface,idfact,idutrj,idmnts,numvar,idfstp,   &
                   iatm,locf,numatf(ixmolc),numfrm(ixmolc)
      integer(8):: sec
      real(8):: fdcrdm(3,ixnatm),velo(3,ixnatm),fdepot(maxene),lamda,  &
                curtim,et,ek,temp,rmsf,rmsd,maxf,r(3,3),cm1(3),cm2(3), &
                tempem(ixmolc),ekem(ixmolc),fdgradback(3,ixnatm),time, &
                tepotcc(nfrg),lambda_vt

!*********************************************************************

      ier = 0      
      call calrsf(iynvar,iylvar,grad(1:3,1:ixnatm,nfrg),rmsf,maxf,locf)
      call caltmp(iynvar,iylvar,velo(1:3,1:ixnatm),fxmass(1:ixnatm),   &
                  idnfre,temp,ek,lambda_vt,lambda_m)
      if ( idfbst .eq. 2 ) then
        call bstft(ixcend(nstpcn),fucord,fdcrdm,fxmass(1:ixnatm),0,r,  &
                   cm1,cm2,rmsd,ier)
        if ( ier .ne. 0 ) then
          rmsd = 0.d0 ; ier = 0
        endif
      else
        rmsd = 0.d0
      endif

      call system_clock(sec)
      et = ek + fdepot(1)
      time = dble(sec-fxcpus)/dble(cr)
      numvar = 0
      do iatm = 1,ixnatm
        if ( iytvar(iatm) .eq. 1 ) then
          numvar = numvar + 1
          fdgradback(1:3,numvar) = fdcrdm(1:3,iatm)
        endif
      enddo
      call outtra(iducrd,iduvel,iduerg,idutrj,iprint,ier,idfacc,idfacv,&
                  idface,idfact,outatm,maxene,fdgradback,velo,idloop,  &
                  curtim,et,ek,temp,fdepot,rmsf,iyn15v,iyn15h,rmsd,    &
                  time,idmonc,idmonv,idmone,idmnts,idflog,tepotcc)
      if ( ier .ne. 0 ) return

      if ( idmonl.gt.0 .and. mod(idloop,idmonl).eq.0 ) then
        call caletp(maxatm,maxshk,mxashk,iynvar,iugshk,ixamol,fxmass,  &
                    iylvar,iuhshk,iuashk,idfstp,velo,gascst,joucal,    &
                    tempem,ekem,numatf,numfrm)
        call outtlg(iprint,ier,maxene,cxmolc,ixsqml,idflog,idloop,     &
                    curtim,et,ek,fdepot,temp,rmsf,rmsd,lamda,time,     &
                    iyeflg,cyenam,tempem,ekem,numatf,numfrm,maxf,locf)
      endif

!****************************************

      return
      end subroutine monmd


!=======================================================================


      subroutine outtra(iducrd,iduvel,iduerg,idutrj,iprint,ier,idfacc, &
                        idfacv,idface,idfact,outatm,maxene,cord,vel,   &
                        idloop,sitime,et,ek,temp,ep,rmsf,iyn15v,iyn15h,&
                        rmsd,sec,idmonc,idmonv,idmone,idmnts,idflog,   &
                        tepotcc)

      use COMBAS,only: ixnatm
      use COMCMM,only: nfrg

      implicit none

      integer(4):: iducrd,iduvel,iduerg,idutrj,iprint,ier,idfacc,      &
                   idfacv,idface,idfact,outatm,maxene,idloop,iyn15v,   &
                   iyn15h,idmonc,idmonv,idmone,idmnts,idflog
      real(8):: cord(3,ixnatm),vel(3,ixnatm),sitime,et,ek,temp,        &
                ep(maxene),rmsf,rmsd,sec,tepotcc(nfrg)
 
!*********************************************************
 
      ier = 0
!     <<<  OUTPUT COORDINATE DATA  >>>
      if ( idfacc.ge.0 .and. iducrd.gt.0 .and. idmonc.gt.0 .and.       &
           mod(idloop,idmonc).eq.0 ) then
        call outtco(iducrd,idfacc,outatm,cord,idloop,sitime,et,ek,temp,&
                    ep(1),rmsf,iyn15v,iyn15h,rmsd,sec,ier,tepotcc)
        if ( ier .ne. 0 ) then
          write(iprint,*)'ERROR> MD '
          write(iprint,*)'   DATA WRITE ERROR IN TRAJECTORY '
          write(iprint,*)'   COORDIANTE DATA '
          return
        endif
      endif
  
!     <<<  OUTPUT VELOCITY DATA  >>>
      if ( idfacv.ge.0 .and. iduvel.gt.0 .and. idmonv.gt.0 .and.       &
           mod(idloop,idmonv).eq.0 ) then
        call outtve(iduvel,idfacv,outatm,vel,idloop,sitime,et,ek,temp, &
                    ep(1),rmsf,iyn15v,iyn15h,rmsd,sec,ier)
        if ( ier .ne. 0 ) then
          write(iprint,*)'ERROR> MD '
          write(iprint,*)'   DATA WRITE ERROR IN TRAJECTORY '
          write(iprint,*)'   VELOCITY DATA '
          return
        endif
      endif

!     <<<  OUTPUT ENERGY DATA  >>>
      if ( idface.ge.0 .and. iduerg.gt.0 .and. idmone.gt.0 .and.       &
           mod(idloop,idmone).eq.0 ) then
        call outten(iduerg,idface,maxene,idloop,sitime,et,ek,temp,ep,  &
                    rmsf,iyn15v,iyn15h,rmsd,sec,ier)
        if ( ier .ne. 0 ) then
          write(iprint,*)'ERROR> MD '
          write(iprint,*)'   DATA WRITE ERROR IN TRAJECTORY '
          write(iprint,*)'   ENERGY DATA '
          return
        endif
      endif
 
!     <<<  OUTPUT MONITORED TRAJECTORY DATA  >>>
      if ( idfact.ge.0 .and. idutrj.gt.0 .and. idmnts.gt.0 .and.       &
           mod(idloop,idmnts).eq.0 ) then
        call mnstmd(iprint,ier,idfact,cord,idutrj,idloop,idflog)
        if ( ier .ne. 0 ) then
          write(iprint,*)'ERROR> MD '
          write(iprint,*)'   DATA WRITE ERROR IN TRAJECTORY '
          write(iprint,*)'   COORDIANTE DATA '
          return
        endif
      endif

!**************************************

      return
      end subroutine outtra


!=========================================================================
 

      subroutine outtco(ioutco,iflag,outatm,cord,idloop,sitime,et,ek,  &
                        temp,ep,rmsf,iyn15v,iyn15h,rmsd,sec,ier,tepotcc)

      use COMBAS,only: ixnatm,ixfbou
      use COMCMM,only: lambda,nfrg

      implicit none

      integer(4):: ioutco,iflag,outatm,idloop,iyn15v,iyn15h,ier
      real(8):: cord(3,ixnatm),sitime,et,ek,temp,ep,rmsf,rmsd,sec,     &
                tepotcc(nfrg),ttepotcc(3)
 
!*******************************************

      ier = 0 ; ttepotcc(1:3) = 0.d0
      ttepotcc(1:nfrg) = tepotcc(1:nfrg)

      if ( ixfbou .eq. 1 ) call PBcord

      select case ( iflag )
        case (1)
          write(ioutco,err=800)                                        &
!               idloop,sngl(sitime),sngl(sec),sngl(et),sngl(ek),        &
!               sngl(temp),sngl(ep),sngl(rmsf),iyn15v,iyn15h,sngl(rmsd)
               idloop,sngl(sitime),sngl(sec),sngl(et),sngl(ek),        &
               sngl(temp),sngl(ep),sngl(lambda),iyn15v,iyn15h,         &
               sngl(rmsd),sngl(ttepotcc(1:3))
          write(ioutco,err=800)sngl(cord(1:3,1:outatm))

        case (2)
          write(ioutco,err=800)                                        &
!               idloop,sitime,sec,et,ek,temp,ep,rmsf,iyn15v,iyn15h,rmsd
               idloop,sitime,sec,et,ek,temp,ep,lambda,iyn15v,iyn15h,   &
               rmsd,ttepotcc(1:3)
          write(ioutco,err=800)cord(1:3,1:outatm)

        case (0)
          write(ioutco,*,err=800)                                      &
!               idloop,sitime,sec,et,ek,temp,ep,rmsf,iyn15v,iyn15h,rmsd
               idloop,sitime,sec,et,ek,temp,ep,lambda,iyn15v,iyn15h,   &
               rmsd,ttepotcc(1:3)
          write(ioutco,*,err=800)cord(1:3,1:outatm)
      end select

!**********************************

      return
800   ier = -1 ; return
      end subroutine outtco


!===========================================================================


      subroutine outtve(ioutve,iflag,outatm,vel,idloop,sitime,et,ek,   &
                        temp,ep,rmsf,iyn15v,iyn15h,rmsd,sec,ier)
 
      use COMBAS,only: ixnatm

      implicit none

      integer(4):: ioutve,iflag,outatm,idloop,iyn15v,iyn15h,ier
      real(8):: vel(3,ixnatm),sitime,et,ek,temp,ep,rmsf,rmsd,sec
 
!************************************************

      ier = 0
      select case ( iflag )
        case (1)
          write(ioutve,err=800)idloop,sngl(sitime),sngl(sec),sngl(et), &
                               sngl(ek),sngl(temp),sngl(ep),sngl(rmsf),&
                               iyn15v,iyn15h,sngl(rmsd)
          write(ioutve,err=800)sngl(vel(1:3,1:outatm))

        case (2)
          write(ioutve,err=800)idloop,sitime,sec,et,ek,temp,ep,rmsf,   &
                               iyn15v,iyn15h,rmsd
          write(ioutve,err=800)vel(1:3,1:outatm)

        case (0)
          write(ioutve,*,err=800)idloop,sitime,sec,et,ek,temp,ep,rmsf, &
                                 iyn15v,iyn15h,rmsd
          write(ioutve,*,err=800)vel(1:3,1:outatm)
      end select

!*********************************

      return
800   ier = -1 ; return
      end subroutine outtve
 

!==========================================================================


      subroutine outten(iduerg,iflag,maxene,idloop,sitime,et,ek,temp,  &
                        ep,rmsf,iyn15v,iyn15h,rmsd,sec,ier)

      implicit none

      integer(4):: iduerg,iflag,maxene,idloop,iyn15v,iyn15h,ier
      real(8):: sitime,et,ek,temp,ep(maxene),rmsf,rmsd,sec
 
!******************************************

      ier = 0
      select case ( iflag )
        case (1)
          write(iduerg,err=800)idloop,sngl(sitime),sngl(sec)
          write(iduerg,err=800)sngl(et),sngl(ek),sngl(temp),           &
                               sngl(ep(1:maxene)),sngl(rmsf),iyn15v,   &
                               iyn15h,sngl(rmsd)

        case (2)
          write(iduerg,err=800)idloop,sitime,sec
          write(iduerg,err=800)et,ek,temp,ep(1:maxene),rmsf,iyn15v,    &
                               iyn15h,rmsd

        case (0)
          write(iduerg,*,err=800)idloop,sitime,sec
          write(iduerg,*,err=800)et,ek,temp,ep(1:maxene),rmsf,iyn15v,  &
                                 iyn15h,rmsd
      end select

!**************************************

      return
800   ier = -1 ; return
      end subroutine outten


!========================================================================


      subroutine caletp(maxatm,maxshk,mxashk,iynvar,iugshk,ixamol,     &
                        fxmass,iylvar,iuhshk,iuashk,imfstp,vel,gascst, &
                        joucal,temp,ek,numatf,numfre) 
 
      use COMBAS,only: ixnatm,ixmolc

      implicit none

      integer(4):: maxatm,maxshk,mxashk,iynvar,iugshk,ixamol(maxatm),  &
                   iylvar(iynvar),iuhshk(maxshk),iuashk(mxashk,maxshk),&
                   imfstp,numatf(ixmolc),numfre(ixmolc),imol,iatm,ivar,&
                   icon,igrp
      real(8):: fxmass(maxatm),vel(3,ixnatm),gascst,joucal,            &
                temp(ixmolc),ek(ixmolc),rinv
 
!**************************************************
 
!     <<<  INITILIZATION  >>>
      rinv = 0.5d0 / ( joucal * 1.0d+3 )
      ek(1:ixmolc) = 0.d0 ; numfre(1:ixmolc) = 0 ; numatf(1:ixmolc) = 0
 
!     <<<  CALCULATION OF KINETIC ENERGY  >>>
      do ivar = 1,iynvar
        iatm = iylvar(ivar)
        imol = ixamol(iatm)
        ek(imol) = ek(imol) + fxmass(iatm) * ( vel(1,ivar)*vel(1,ivar) &
                   + vel(2,ivar)*vel(2,ivar) + vel(3,ivar)*vel(3,ivar) )
        numfre(imol) = numfre(imol) + 3
        numatf(imol) = numatf(imol) + 1
      enddo

      if ( imfstp.eq.1 .or. imfstp.eq.2 ) then
        numfre(1) = numfre(1) - 3
      elseif ( imfstp .eq. 3 ) then
        numfre(1) = numfre(1) - 6
      endif
 
!     <<<  DESCEND NUMBER OF FREEDOM UNDER SHAKE CONSTARINTS  >>>
      do igrp = 1,iugshk
        iatm = iuashk(1,igrp)
        imol = ixamol(iatm)
        icon = ( iuhshk(igrp)+1) * iuhshk(igrp) / 2
        numfre(imol) = numfre(imol) - icon
      enddo

!     <<<  CALCULATION OF TEMERATURE  >>>
!          EK AT TEMP CALCULATION   JOULE/MOL
!          EK AT FINAL              KCAL/MOL
      do imol = 1,ixmolc
        ek(imol) = ek(imol) * 1.0d+7
        if ( numfre(imol) .ne. 0 ) then
          temp(imol) = ek(imol) / ( gascst * dble(numfre(imol)))
          ek(imol) = ek(imol) * rinv
        else
          temp(imol) = 0.d0 ; ek(imol) = 0.d0
        endif
      enddo
 
!***************************************

      return
      end subroutine caletp


!===================================================================


      subroutine outtlg(iprint,ier,maxene,cxmolc,ixsqml,idflog,curlop, &
                        curtim,et,ek,ep,temp,rmsf,rmsd,lamda,sec,      &
                        iyeflg,cyenam,tempem,ekem,numatf,numfrm,maxf,  &
                        locmaxf)

      use COMBAS,only: ixmolc,ixachn,cxresn,ixares,cxatmn

      implicit none

      integer(4):: iprint,ier,maxene,ixsqml(ixmolc),idflog,curlop,     &
                   iyeflg(maxene),locmaxf,numatf(ixmolc),              &
                   numfrm(ixmolc),i,j
      real(8):: curtim,et,ek,ep(maxene),temp,rmsf,rmsd,lamda,sec,maxf,&
                tempem(ixmolc),ekem(ixmolc)
      character(*):: cxmolc(ixmolc),cyenam(maxene)
 
!*************************************************************

      ier = 0
!!      select case ( idflog )
!!        case (1)
!!          write(iprint,800,err=900)curtim*0.001,et,temp,ep(1)
!!800       format(' '/' TIME (PSEC)           : ',f15.5,                &
!!                   ' TOTAL ENERGY          : ',e15.7/                  &
!!                   ' TEMPERATURE           : ',f15.7,                  &
!!                   ' POTENTIAL             : ',e15.7)
!!        case (2)
          write(iprint,810,err=900)curlop,curtim*0.001,et,ep(1),temp,  &
                ek,rmsf,rmsd,sec,lamda,maxf,ixachn(locmaxf),           &
                cxresn(locmaxf),ixares(locmaxf),cxatmn(locmaxf)
810       format(' '/' MD LOOP NUMBER        : ',i15,                  &
                 ' TIME (PSEC)           : ',f15.5/                    &
                 ' TOTAL      (KCAL/MOL) : ',e15.7,                    &
                 ' POTENTIAL  (KCAL/MOL) : ',e15.7/                    &
                 ' TEMPERATURE (K)       : ',f15.7,                    &
                 ' KINETIC    (KCAL/MOL) : ',e15.7/                    &
                 ' R.M.S.F. (KCAL/MOL*A) : ',e15.7,                    &
                 ' RMSD (ANGSTROMS)      : ',e15.7/                    &
                 ' LAP CPU TIME  (SEC)   : ',f15.7,                    &
                 ' SCALING FACTOR        : ',e15.7/                    &
                 ' MAX FORCE(KCAL/MOL*A) ; ',e15.7,                    &
                 ' CHAIN ',i5," ",a8,i5," ",a8)
          j = 0
          do i = 1,maxene
            if ( iyeflg(i) .eq. 1 ) then
              write(iprint,'("  ",a6," : ",e15.7,$)')cyenam(i),ep(i)
              j = j + 1
            endif
            if ( j .eq. 3 ) then
              write(iprint,*) ; j = 0
            endif
          enddo
          if ( j .ne. 0 ) write(iprint,*)
!!      end select

      if ( ixmolc .ge. 2 ) then
        do i = 1,ixmolc
          write(iprint,'(5x,a20,i5,2i10,2e15.7)',err=900)cxmolc(i),    &
                ixsqml(i),numatf(i),numfrm(i),ekem(i),tempem(i)
        enddo
      endif

!*******************************************

      return
900   write(iprint,*)'ERROR> MD '
      write(iprint,*)'   WRITE ERROR DURING OUTPUT LOG '
      ier = -1 ; return
      end subroutine outtlg
