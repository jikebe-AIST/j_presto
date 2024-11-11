
      subroutine monmin(iprint,inloop,stplen,inflog,infbst,inuana,     &
                        fnepot)

!*******************************************************************
!
!     THIS SUBROUTINE IS FOR MONITORING MINIMIZATION
!     OUTPUT ANALYSIS DATA AND OUTPUT LOG
!
!*******************************************************************
 
      use COMBAS ; use COMERG ; use COMMIS ; use COMCMMC
      use COMCMM,only: nfrg 

      implicit none

      ! inloop : current loop number, inflog : flag of format of log
      ! infbst : flag for best-fit
      ! inuana : logical unit # for output analysis data
        integer(4),intent(in):: iprint,inloop,inflog,infbst,inuana
      ! stplen : step length (ANGSTROMS), fnepot : current potential energy
        real(8),intent(in):: stplen,fnepot(maxene)

      real(8):: rmxfor,rmsf,r(3,3),cm1(3),cm2(3),rmsd,rmsc,time
      integer(4):: iatsel,ier,i,j
      integer(8):: sec
 
!**********************************************
 
!     <<<  CALCULATE RMSF AND RMSD  >>>
      call calrsf(iynvar,iylvar,grad(1:3,1:ixnatm,nfrg),rmsf,rmxfor,   &
                  iatsel)
      if ( infbst .eq. 2 ) then
        call bstft(ixcend(nstpcn),fucord,cord,fxmass(1:ixnatm),0,r,cm1,&
                   cm2,rmsd,ier)
        if ( ier .ne. 0 ) then
          rmsd = 0.d0 ; ier = 0
        endif
      else
        rmsd = 0.d0
      endif
 
      call system_clock(sec)
      time  = dble(sec-fxcpus)/dble(cr)
      RMSC = STPLEN / DBLE(IYNVAR)
 
!     <<<  OUTPUT LOG  >>>
      select case (inflog)
        case (1)
          write(iprint,*)""
          write(iprint,'(a,i15,$)')"  LOOP NUMBER      : ",inloop
          write(iprint,'(a,e15.7)')"  POTENTIAL        : ",fnepot(1)
          write(iprint,'(a,e15.7,$)')"  R.M.S.F.         : ",rmsf
          write(iprint,'(a,e15.7)')"  RMSD             : ",rmsd
        case (2)
          write(iprint,*)""
          write(iprint,'(a,i15,$)')"  MINI LOOP NUMBER      : ",inloop
          write(iprint,'(a,e15.7)')"  POTENTIAL (KCAL/MOL)  : ",       &
                                   fnepot(1)
          write(iprint,'(a,e15.7,$)')"  STEP LENGTH (A)       : ",rmsc
          write(iprint,'(a,e15.7)')"  R.M.S.F. (KCAL/MOL*A) : ",rmsf
          write(iprint,'(a,e15.7,$)')"  RMSD (ANGSTROMS)      : ",rmsd
          write(iprint,'(a,e15.7)')"  LAP CPU TIME (SEC)    : ",time
          write(iprint,'(a,e15.7,$)')"  MAXIMUM FORCE         : ",rmxfor
          write(iprint,'(a8,i5,x,a8,i5,x,a8)')"  CHAIN ",              &
            ixachn(iatsel),cxresn(iatsel),ixares(iatsel),cxatmn(iatsel)
          j = 0 ; write(iprint,'(3x,$)')
          do i = 1,maxene
            if ( iyeflg(i) .eq. 1 ) then
              write(iprint,'(x,a6," : ",e15.7,$)')cyenam(i),fnepot(i)
              j = j + 1
            endif 
            if ( j .eq. 3 ) then
              write(iprint,*) ; j = 0
            endif
          enddo
          if ( j .ne. 0 ) write(iprint,*)
      end select

!     <<<  OUTPUT ANALYSIS DATA  >>>
      if ( inuana .gt. 0 ) then
        write(inuana)inloop,0.d0,time
        write(inuana)stplen,rmsc,0.d0,fnepot(1:maxene),rmsf,iyn15v,    &
                     iyn15h,rmsd
      endif

!*******************************

      return
      end subroutine monmin
