 
      subroutine mini(iread,iprint,ier,onend)

!*******************************************************************
!
!     ENERGY MINIMIZATION
!       STEEPESET DESCENT MINIMIZATION  (ORIGINAL)
!       CONJUGATE GRADIENT MINIMIZATION (ORIGINAL)
!
!*******************************************************************

      use COMBAS ; use COMERG 
      !$ use omp_lib

      implicit none

      integer(4),intent(in):: iread,iprint
      integer(4),intent(inout):: ier
      logical(1),intent(inout):: onend

      integer(4):: inmthd,inlopl,inupdl,inmons,inflog,infbst,inllpl,   &
                   inuana,inloop,i
      real(8):: fncovg,fnstpi,fnupsl,fndnsl,fncovl
      real(8):: fnepot(maxene)
      character(80):: cnnana

      integer(4):: secst,secen,secst1,secen1
 
!**************************************************************
 
      ! <<<  INITILIZATION  >>>
 
      call system_clock(secst)
      ier = 0 ; inloop = 0
      write(iprint,*)' '
      write(iprint,*)'INFORMATION> MINI (V2.0)'
      write(iprint,*)'     ENERGY MINIMIZATION '
      write(iprint,*)' '

      ! Make inverse of mass
      allocate(ifxmass(ixnatm))
      !$OMP parallel default (none)                                  & !
      !$OMP private(i)                                               & !
      !$OMP shared(ifxmass,fxmass,ixnatm)
      !$OMP do schedule (static)
      do i = 1,ixnatm
        ifxmass(i) = 1.d0 / fxmass(i)
      enddo
      !$OMP end do
      !$OMP end parallel
 
      ! <<<  READ CONTROL DATA  >>>
      call inpmin(iread,ier,iprint,onend,inmthd,inlopl,inupdl,inmons,  &
                  inflog,infbst,fncovg,fnstpi,fnupsl,fndnsl,fncovl,    &
                  inllpl,inuana,cnnana)
      call incel(iprint)
      if ( ier .ne. 0 ) goto 9999

      ! <<<  CALCULATE THE CENTER OF MASS OF A PROTEIN >>>
      if ( ixfbou.eq.1 .and. ixcbou.eq.1 ) then
        call calcm(ixcend(nstpcn),fxmass(1:ixnatm),fxcbou,ier)
        if ( ier .ne. 0 ) goto 9999
        write(iprint,*)'INFORMATION> MINI '
        write(iprint,*)                                                &
         '   CENTER OF BOUNDARY IS MASS CENTER OF FIRST CHAIN '
        write(iprint,*)'             X        : ',FXCBOU(1)
        write(iprint,*)'             Y        : ',FXCBOU(2)
        write(iprint,*)'             Z        : ',FXCBOU(3)
        write(iprint,*)' '
      endif
 
      ! <<<  PREPAR ANALYSIS DATA  >>>
      if ( inuana .gt. 0 ) then
        call flopen(inuana,cnnana,21,'NULL',0,ier)
        if ( ier .ne. 0 ) goto 9999
      endif
 
      ! <<<  CALCULATE ENERGY OF START STRUCTURE  >>>
      call inimin(iprint,ier,inloop,infbst,inuana,fnepot)
      if ( ier .ne. 0 ) goto 9999

      ! <<<  RUN MINIMIZATION  >>>
      if ( inlopl .gt. 0 ) then
        call system_clock(secst1)
        call runmin(iprint,ier,inmthd,inloop,inlopl,inupdl,inmons,     &
                    inflog,infbst,fncovg,fnstpi,fnupsl,fndnsl,fncovl,  &
                    inllpl,inuana,fnepot)
        call system_clock(secen1)
      else
        secst1 = 0.d0 ; secen1 = 0.d0
      endif

      ! <<<  FINAL PROCESS OF MINIMIZATION  >>>
      call finmin(iprint,ier,inloop,inflog,infbst,inuana,fnepot)
      if ( inuana.gt.0 ) then
        call flclos(inuana,10,ier)
        if ( ier .ne. 0 ) then
          write(iprint,*)' '
          write(iprint,*)'ERROR> MINI '
          write(iprint,*)'    ANALYSIS FILE OPEN ERROR OR CLOSE ERROR '
          ier = -1
        endif
      endif

9999  call system_clock(secen)
      write(iprint,*)" "
      write(iprint,*)"INFORMATION> MINI "
      write(iprint,'(12x,a22,f15.5)')"TOTAL CPU TIME (S)  : ",         &
                                     dble(secen-secst)/dble(cr)
      write(iprint,'(12x,a22,f15.5)')"  MIN LOOP          : ",         &
                                     dble(secen1-secst1)/dble(cr)
      write(iprint,'(12x,a22,f15.5)')"  OTHERS            : ",         &
                              dble(secen-secst-secen1+secst1)/dble(cr)
 
!****************************

      return
      end subroutine mini


!===================================================================


      subroutine inimin(iprint,ier,inloop,infbst,inuana,fnepot)

!*******************************************************************
!
!     CALCULATION ENERGY OF INITIAL STRUCTURE
!
!*******************************************************************
 
      use COMBAS ; use COMERG ; use COMMIS 
      use COMCMM ; use COMCMMC
      !$ use omp_lib

      implicit none

      integer(4),intent(in):: iprint
      integer(4),intent(inout):: ier,inloop,infbst,inuana
      real(8),intent(inout):: fnepot(maxene)

      real(8):: fnepcc(maxene,nfrg)
 
      integer(4):: sec1,sec2,i,j
      real(8):: rtmp
 
!****************************************
 
      ier = 0
      ! <<<  SHAKE INITIAL STRUCTURE  >>>
      if ( iyfshk .ne. 0 ) then
        !$OMP parallel default (none)                                & !
        !$OMP private(i)                                             & !
        !$OMP shared(grad,ixnatm,nfrg,cord)
        !$OMP do schedule (static)
        do i = 1,ixnatm
          grad(:,i,nfrg) = cord(:,i)
        enddo
        !$OMP end do
        !$OMP end parallel
        call system_clock(sec1)
        select case(ixfbou)
        case (1)
          call shake_PB(cord,grad(1:3,1:ixnatm,nfrg),iprint,ier)
        case default
          call shake_CB(cord,grad(1:3,1:ixnatm,nfrg),iprint,ier)
        end select
        call system_clock(sec2)
        call outcpu(iprint,'SHAKE ',sec1,sec2,cr,cm)
        if ( ier .ne. 0 ) return
      endif

     ! <<<  MAKE 1-5 INTERACTION TABLE  >>>
      call system_clock(sec1)
      call celprp(iprint,2,ier)
      call system_clock(sec2)
      call outcpu(iprint,'CELPRP',sec1,sec2,cr,cm)
      if ( ier .ne. 0 ) return

      ! <<  CALCULATION ENERGY AND GRADIENT OF INITIAL STRRUCTURE  >>
      call system_clock(sec1)
      call dnergy(fnepot,fnepcc)
      if ( cluster_method .eq. "LMD" ) then
        rtmp = lambda*lambda
        select case (nfrg)
        case (1)
          !$OMP parallel default(none)                               & !
          !$OMP private(i)                                           & !
          !$OMP shared(grad,rtmp,ixnatm)
          !$OMP do schedule (static)
          do i = 1,ixnatm
            grad(:,i,1) = grad(:,i,1)*rtmp
          enddo
          !$OMP end do
          !$OMP end parallel
        case (3)
          !$OMP parallel default(none)                               & !
          !$OMP private(i)                                           & !
          !$OMP shared(grad,rtmp,lambda,ixnatm)
          !$OMP do schedule (static)
          do i = 1,ixnatm
            grad(:,i,3) = grad(:,i,1)*rtmp +                           &
                          grad(:,i,2)*lambda + grad(:,i,3)
          enddo
          !$OMP end do
          !$OMP end parallel
        end select
      elseif ( cluster_method .eq. "LMD2" ) then
        select case (nfrg)
        case (1)
          !$OMP parallel default(none)                               & !
          !$OMP private(i)                                           & !
          !$OMP shared(grad,lambda,ixnatm)
          !$OMP do schedule (static)
          do i = 1,ixnatm
            grad(:,i,1) = grad(:,i,1)*lambda
          enddo
          !$OMP end do
          !$OMP end parallel
        case (3)
          rtmp = sqrt(lambda)
          !$OMP parallel default(none)                               & !
          !$OMP private(i)                                           & !
          !$OMP shared(grad,rtmp,lambda,ixnatm)
          !$OMP do schedule (static)
          do i = 1,ixnatm
            grad(:,i,3) = grad(:,i,1)*lambda +                         &
                          grad(:,i,2)*rtmp + grad(:,i,3)
          enddo
          !$OMP end do
          !$OMP end parallel
        end select
      endif

      call system_clock(sec2)
      call outcpu(iprint,'DNERGY',sec1,sec2,cr,cm)
      write(iprint,*)' '
      write(iprint,*)'INFORMATION> MINI '
      write(iprint,*)'    ENERGY OF INITIAL STRUCTURE '
      write(iprint,*)' '
      call otdene(iprint,inloop,infbst,inuana,fnepot)

!**************************************

      return
      end subroutine inimin


!=====================================================================

      subroutine runmin(iprint,ier,inmthd,inloop,inlopl,inupdl,inmons, &
                        inflog,infbst,fncovg,fnstpi,fnupsl,fndnsl,     &
                        fncovl,inllpl,inuana,fnepot)

!*******************************************************************
!
!     RUN ENERGY MINIMIZATION
!       1) STEEPEST DESCENT MINIMIZATION
!       2) CONJUGATE GRADINET MINIMIZATION
!
!*******************************************************************
 
      use COMBAS ; use COMERG

      implicit none

      integer(4),intent(in):: iprint,inmthd
      integer(4),intent(inout):: ier,inloop,inlopl,inupdl,      &
                 inmons,inflog,infbst,inllpl,inuana
      real(8):: fncovg,fnstpi,fnupsl,fndnsl,fncovl
      real(8):: fnepot(maxene)

!*************************************************

      ier = 0
      ! <<<  STEEPEST DESCENT MINIMIZATION  >>>
      if ( inmthd .eq. 1 ) then
        call steep(iprint,ier,inloop,inlopl,inupdl,inmons,inflog,      &
                   infbst,fncovg,fnstpi,fnupsl,fndnsl,inuana,fnepot)
     ! <<<  CONJUGATE GRADIENT MINIMIZATION  >>>
      elseif ( inmthd .eq. 2 ) then
        call conj(iprint,ier,inloop,inlopl,inupdl,inmons,inflog,infbst,&
                  fncovg,fnstpi,fncovl,inllpl,inuana,fnepot)
      endif

!*****************************

      return
      end subroutine runmin


!=====================================================================


      subroutine finmin(iprint,ier,inloop,inflog,infbst,inuana,fnepot)

!*******************************************************************
!
!     FINAL PROCESS OF MINIMIZATION
!
!*******************************************************************
 
      use COMBAS ; use COMERG ; use COMCMM
      use COMCMMC, only: grad

      implicit none

      integer(4),intent(in):: iprint
      integer(4),intent(inout):: ier,inloop,inflog,infbst,inuana
      real(8),intent(inout):: fnepot(maxene)
      real(8):: fnepcc(maxene,nfrg),rtmp
      integer(4):: i,j
 
!****************************************************

      ier = 0
      ! <<<  UPDATE INTERACTION TABLE  >>>
      call celprp(iprint,inflog,ier)
      if ( ier .ne. 0 ) return

      ! <<<  RE-CALCULATE ENERGY OF FINAL STRUCTURE  >>>
      call dnergy(fnepot,fnepcc)
      if ( cluster_method .eq. "LMD" ) then
        rtmp = lambda*lambda
        select case (nfrg)
        case (1)
          !$OMP parallel default(none)                               & !
          !$OMP private(i)                                           & !
          !$OMP shared(grad,rtmp,ixnatm)
          !$OMP do schedule (static)
          do i = 1,ixnatm
            grad(:,i,1) = grad(:,i,1)*rtmp
          enddo
          !$OMP end do
          !$OMP end parallel
        case (3)
          !$OMP parallel default(none)                               & !
          !$OMP private(i)                                           & !
          !$OMP shared(grad,rtmp,lambda,ixnatm)
          !$OMP do schedule (static)
          do i = 1,ixnatm
            grad(:,i,3) = grad(:,i,1)*rtmp +                           &
                          grad(:,i,2)*lambda + grad(:,i,3)
          enddo
          !$OMP end do
          !$OMP end parallel
        end select
      elseif ( cluster_method .eq. "LMD2" ) then
        select case (nfrg)
        case (1)
          !$OMP parallel default(none)                               & !
          !$OMP private(i)                                           & !
          !$OMP shared(grad,lambda,ixnatm)
          !$OMP do schedule (static)
          do i = 1,ixnatm
            grad(:,i,1) = grad(:,i,1)*lambda
          enddo
          !$OMP end do
          !$OMP end parallel
        case (3)
          rtmp = sqrt(lambda)
          !$OMP parallel default(none)                               & !
          !$OMP private(i)                                           & !
          !$OMP shared(grad,rtmp,lambda,ixnatm)
          !$OMP do schedule (static)
          do i = 1,ixnatm
            grad(:,i,3) = grad(:,i,1)*lambda +                         &
                          grad(:,i,2)*rtmp + grad(:,i,3)
          enddo
          !$OMP end do
          !$OMP end parallel
        end select
      endif
 
      ! <<<  OUTPUT MONITORING DATA  >>>
      write(iprint,*)' '
      write(iprint,*)'INFORMATION> MINI '
      write(iprint,*)'   ENERGY OF FINAL STRUCTURE '
      write(iprint,*)' '
      call otdene(iprint,inloop,infbst,inuana,fnepot)

!*************************************

      return
      end subroutine finmin


!===================================================================== 


      subroutine otdene(iprint,inloop,infbst,inuana,fnepot)

!*******************************************************************
!
!     OUTPUT DETAIL ENERGY
!
!*******************************************************************
 
      use COMBAS ; use COMERG ; use COMMIS ; use COMCMMC
      use COMCMM,only: nfrg

      implicit none

      integer(4),intent(in):: iprint
      integer(4),intent(inout):: inloop,infbst,inuana
      real(8),intent(inout):: fnepot(maxene)

      integer(4):: numgrd(8)
      real(8):: rmxfor,rmsf,absfor,work
      real(8):: rotmat(3,3),cenref(3),cencom(3),rmsd 

      integer(4):: iatm,iene,ier,iatsel,ivar
      integer(8):: sec
 
!********************************************************
 
      ! <<<  CALCULATE RMSF AND RMSD  >>>
      call calrsf(iynvar,iylvar,grad(1:3,1:ixnatm,nfrg),rmsf,rmxfor,   &
                  iatsel)
      if ( infbst .eq. 2 ) then
        call bstft(ixcend(nstpcn),fucord,cord,fxmass(1:ixnatm),0,      &
                   rotmat,cenref,cencom,rmsd,ier)
        if ( ier .ne. 0 ) then
          rmsd = 0.d0 ; ier = 0
        endif
      else
        rmsd = 0.d0
      endif
 
      ! <<<  CALCULATE DISTRIBUTION OF GRADIENT  >>>
      numgrd(1:8) = 0
      do ivar = 1,iynvar
        iatm = iylvar(ivar)
        absfor = grad(1,iatm,nfrg)*grad(1,iatm,nfrg) +                 &
                 grad(2,iatm,nfrg)*grad(2,iatm,nfrg) +                 &
                 grad(3,iatm,nfrg)*grad(3,iatm,nfrg)
        if ( absfor .ge. 1.d+6 ) then
          numgrd(1) = numgrd(1) + 1
        elseif ( absfor .ge. 1.d+4 ) then
          numgrd(2) = numgrd(2) + 1
        elseif ( absfor .ge. 1.d+2 ) then
          numgrd(3) = numgrd(3) + 1
        elseif ( absfor .ge. 1.d0 ) then
          numgrd(4) = numgrd(4) + 1
        elseif ( absfor .ge. 1.d-2 ) then
          numgrd(5) = numgrd(5) + 1
        elseif ( absfor .ge. 1.d-4 ) then
          numgrd(6) = numgrd(6) + 1
        elseif ( absfor .ge. 1.d-6 ) then
          numgrd(7) = numgrd(7) + 1
        else
          numgrd(8) = numgrd(8) + 1
        endif
      enddo

      ! <<<  OUTPUT LOG  >>>
      write(iprint,*)'     LOOP       : ',inloop
      do iene = 1,maxene
        if ( cyenam(iene).ne."      " .and. iyeflg(iene).eq.1 ) then
          write(iprint,*)'     ',cyenam(iene),'     : ',fnepot(iene)
        endif
      enddo
      write(iprint,*)'     RMSF       : ',rmsf
      write(iprint,*)'     MAX FORCE  : ',rmxfor,'  CHAIN ',           &
                     ixachn(iatsel),' ',cxresn(iatsel),ixares(iatsel), &
                     ' ',cxatmn(iatsel)
      write(iprint,*)'     RMSD       : ',rmsd
      write(iprint,*)' '
      write(iprint,*)'   DISTRIBUTION OF FORCE  (KCAL/(MOL*A)) '
      write(iprint,*)'     1.0D+3 <=              ',numgrd(1)
      write(iprint,*)'     1.0D+2 <=    < 1.0D+3  ',numgrd(2)
      write(iprint,*)'     1.0D+1 <=    < 1.0D+2  ',numgrd(3)
      write(iprint,*)'     1.0D0  <=    < 1.0D+1  ',numgrd(4)
      write(iprint,*)'     1.0D-1 <=    < 1.0D0   ',numgrd(5)
      write(iprint,*)'     1.0D-2 <=    < 1.0D-1  ',numgrd(6)
      write(iprint,*)'     1.0D-3 <=    < 1.0D-2  ',numgrd(7)
      write(iprint,*)'                  < 1.0D-3  ',numgrd(8)
      write(iprint,*)'   '
 
      ! <<<  OUTPUT ANALYSIS DATA  >>>
      if ( inuana .gt. 0 ) then
        call system_clock(sec)
        work = 0.d0
        write(inuana)inloop,work,dble(sec-fxcpus)/dble(cr)
        write(inuana)work,work,work,fnepot(1:maxene),rmsf,iyn15v,      &
                     iyn15h,rmsd
      endif

!*************************************

      return
      end subroutine otdene
