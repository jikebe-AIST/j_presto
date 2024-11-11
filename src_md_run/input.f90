
      subroutine input(iread,iprint,ier,onend) 

!*******************************************************************
!
!     INPUT DATA SETS
!
!*******************************************************************
 
      use COMBAS ; use COMERG ; use COMMIS ; use COMCMM
      use COMCMMC

      implicit none

      ! Logical unit number for input & output
        integer(4),intent(in):: iread,iprint
      ! Condition code (0: NO ERROR, negative: ERROR)
        integer(4),intent(out):: ier
      ! Flag for end of control data (.true.: end of control data is
      !                                       detected)
        logical(1),intent(inout):: onend

      integer(4):: iitplf,iicrdf,iifsor
      character(80):: citpln,cicrdn,cishkn,civarn,ciboun,cirefn,cipscn,&
                      cidscn,cidhcn,cimntn,cicmpn,cieCMn,cirep

      integer(4),allocatable:: iwork(:)
      real(8),allocatable:: rwork(:)
      character(4),allocatable:: cwork(:)

      integer(8):: sec1,sec2,secst
      integer(4):: i,j
      real(8):: sectpl,seccrd,secoth

!      equivalence (iwork(1),cwork(1))
 
      logical(4),allocatable:: tmp0(:)

!********************************************************

      call system_clock(secst)
      ier = 0 ; sectpl = 0.d0 ; seccrd = 0.d0
      write(iprint,*)' '
      write(iprint,*)'INFORMATION> INPUT (V3.0)'
      write(iprint,*)'             READ SOME DATA FILES '
 
      ! 0) READ INPUT CONTROL DATA
      call inpinp(iread,iprint,ier,onend,iitplf,iicrdf,iifsor,citpln,  &
                  cicrdn,cishkn,civarn,ciboun,cirefn,cipscn,cidscn,    &
                  cidhcn,cimntn,cicmpn,cieCMn,cirep)
      if ( ier .ne. 0 ) return

      ! 1) READ TOPOLOGY FILE
      if ( citpln .ne. " " ) then
        call system_clock(sec1)
        call rdtopo(iitplf,iprint,citpln,ier)
        call system_clock(sec2)
        sectpl = dble(sec2-sec1)/dble(cr)
        if ( ier .ne. 0 ) return
        !! check atoms having vdW parameter = 0
        allocate(tmp0(maxtyp)) ; tmp0(1:maxtyp) = .false.
        do i = 1,maxtyp
          if ( all(abs(fynbpp(:,:,i)) .lt. 1e-16) ) tmp0(i) = .true.
        enddo
        allocate(zero_vdW(ixnatm)) ; zero_vdW(1:ixnatm) = .false.
        do i = 1,ixnatm
          j = ixatyp(i)
          if ( tmp0(j) ) zero_vdW(i) = .true.
        enddo
        deallocate(tmp0)
        allocate(iichain(2,0:ixnatm))
        nchain = 0 ; iichain(1,1) = 1 ; iichain(2,0) = 0
        if ( ixmolc .ne. 1 ) then
          do i = 1,ixmolc-1
          do j = 1,ixtpcn(i+1)-ixtpcn(i)
            nchain = nchain + 1
            iichain(2,nchain) = iichain(2,nchain-1) + ixatmm(i)
          enddo
          enddo
        endif
        do j = 1,ixnchn-ixtpcn(ixmolc)+1
          nchain = nchain + 1
          iichain(2,nchain) = iichain(2,nchain-1) + ixatmm(ixmolc)
        enddo
        do j = 2,nchain
          iichain(1,j) = iichain(2,j-1) + 1
        enddo
      endif

      ! 2) READ INITIAL COORDINATE FILE
      allocate(iwork(maxatm),rwork(maxatm),cwork(maxatm))
      if ( cicrdn .ne. " " ) then
        write(iprint,*)' '
        write(iprint,*)'INFORMATION> INPUT '
        write(iprint,*)'     READ INITIAL STRUCTURE '
        call system_clock(sec1)
        call rdcord(iicrdf,iprint,cicrdn,rwork,iwork,cwork,ier)
        call system_clock(sec2)
        seccrd = dble(sec2-sec1)/dble(cr)
        if ( ier .ne. 0 ) return
      endif
 
      ! 3) READ CONTROL DATA FOR SETTING VARIABLES
      allocate(iylvar(ixnatm))
      if ( civarn .ne. " " ) then
        call inpvar(1,iprint,civarn,ier)
        if ( ier .ne. 0 ) return
      else
        iynvar = ixnatm ; iytvar(1:ixnatm) = 1
        forall ( i=1:ixnatm ) iylvar(i) = i
      endif

      ! 4) SET PARAMETER OF SHAKE
      if ( cishkn .ne. " " ) then
        call setshk(1,iprint,cishkn,ier)
        if ( ier .ne. 0 ) return
      else
        iyndbn = 0 ; iyndag = 0
      endif

      call mkabsr

      ! 5) READ CONTROL DATA FOR BOUNDARY CONDITION
      if ( ciboun .ne. " " ) then
        call inpbou(1,iprint,ciboun,ier)
        if ( ier .ne. 0 ) return
      endif

      ! 6) READ POSITION CONSTRAINTS DATA
      !    READ REFERENCE COORDINATE DATA
      call inppsc(iprint,cipscn,cirefn,rwork,iwork,cwork,ier)
      if ( ier .ne. 0 ) return

      ! 7) READ DISTANCE CONSTRAINT DATA FILE
      if ( cidscn .ne. " " ) then
        call dsr(1,iprint,cidscn,ier)
        if ( ier .ne. 0 ) return
      endif

      ! 8) READ DIHEDRAL CONSTRAINT DATA FILE
      if ( cidhcn .ne. " " ) then
        call dhr(1,iprint,cidhcn,ier)
        if ( ier .ne. 0 ) return
      else
        iyndtr = 0 ; iyndip = 0
      endif

      ! 9) READ CELL and CLUSTER PARAMETER FILE
      call inpcls(iprint,1,cicmpn,ier)
      if ( ier .ne. 0 ) return

      ! 10) READ A FILE TO SPECIFY ITEMS TO BE MONITORED
      if ( cimntn .ne. " " ) then
        call setmnt(1,iprint,cimntn,ier)
        if ( ier .ne. 0 ) return
      endif

      ! 11) READ A external CMM input file
      if ( cieCMn .ne. " " ) then
        if ( iifsor .eq. 2 ) then
          write(iprint,*)
          write(iprint,*)" !!! CAUTION !!!"
          write(iprint,*)"If you want to use external CMM method,"//   &
                         " you canNOT set SETORI option"
          write(iprint,*)"STOP" ; stop
        elseif ( cipscn .eq. " " ) then
          write(iprint,*)
          write(iprint,*)" !!! CAUTION !!!"
          write(iprint,*)"If you want to use external CMM method,"//   &
                         " you must use position constraints"
        endif
        call setextCMM(iprint,cieCMn)
        if ( ier .ne. 0 ) return
      endif

      ! 11.5) READ GRAVITY PARAMETER FILE
      if ( cirep .ne. " " ) then
        call inprep(iprint,cirep,ier) 
        if ( ier .ne. 0 ) return
      endif

      ! 12) TRANSFORMATION
      !     MASS-CENTER OF CURRENT STRUCTURE IS ORIGIN OF COORINATE
      if ( iifsor .eq. 2 ) then
        call tracen(ixnatm,iprint,fxmass(1:ixnatm),fxcbou)
      endif

      ! 13) CALCULATE NUMBER OF INTERACTIONS
      !     CALCULATE TOTAL CHARGE
      call calint(iprint)
 
      ! 14) MAKE IYNINI  (FOR WATER SYSTEM)
      call chkini
 
      ! 15) OUTPUT CPU TIME
      call system_clock(sec2)
      secoth = dble(sec2-secst)/dble(cr)
      write(iprint,*)' '
      write(iprint,*)'INFORMATION> INPUT '
      write(iprint,'(a34,f15.5)')                                      &
            "            TOTAL CPU TIME (S)  : ",secoth
      write(iprint,'(a34,f15.5)')                                      &
            "              READ TPL          : ",sectpl
      write(iprint,'(a34,f15.5)')                                      &
            "              READ CORD         : ",seccrd
      write(iprint,'(a34,f15.5)')                                      &
            "              OTHER             : ",secoth-sectpl-seccrd

!***********************************

      return
      end subroutine input


!========================================================================

      subroutine rdtopo(iitplf,iprint,citpln,ier)

!*******************************************************************
!
!     INPUT TOPOLOGY FILE
!
!*******************************************************************
 
      use COMBAS ; use COMERG ; use PHYCNS 

      implicit none

      ! Format type of topology file ( 1: FORMATTED, 2: BINARY)
        integer(4),intent(in):: iitplf
      ! Logical unit number for output
        integer(4),intent(in):: iprint
      ! File name of topology file
        character(80):: citpln
      ! Condition code
        integer(4):: ier

      real(8):: rtmp
 
!************************************************

      ier = 0
      if ( iitplf .eq. 1 ) then
        call inptpl(1,iprint,citpln,ier)
      elseif ( iitplf .eq. 2 ) then
        call inptpb(1,iprint,citpln,ier)
      endif
      if ( ier .ne. 0 ) return

      rtmp = 1.d0 / 180.d0
      call chkpha(iyntor,iyptor(1:4,1:iyntor),fytphs,rtmp,iprint,ier)
      if ( ier .ne. 0 ) return
      call chkpha(iynimp,iypimp(1:4,1:iynimp),fyiphs,rtmp,iprint,ier)

!**********************************

      return
      end subroutine rdtopo


!=========================================================================
      
      subroutine rdcord(iicrdf,iprint,cicrdn,rwork,iwork,cwork,ier)

!*******************************************************************
!
!     INPUT COORDINATE FILE
!
!*******************************************************************
 
      use COMBAS ; use COMERG ; use COMCMMC

      implicit none

      ! Format type of coordinate file (1: PDB, 2: BINARY)
        integer(4),intent(in):: iicrdf
      ! Logical unit number for output log
        integer(4),intent(in):: iprint
      ! File name of topology file
        character(80),intent(in):: cicrdn
      ! work arrays
        real(8),intent(inout):: rwork(maxatm)
        integer(4),intent(inout):: iwork(maxatm)
        character(*),intent(inout):: cwork(maxatm)
      ! Condition code
        integer(4),intent(out):: ier

      integer(4):: numatm
      real(8):: Tcord(3,maxatm)
 
!***************************************

      ier = 0
      if ( iicrdf .eq. 1 ) then
        call inppdb(1,iprint,cicrdn,maxatm,numatm,cwork,cwork,iwork,   &
                    Tcord,rwork,rwork,ier)
      elseif ( iicrdf .eq. 2 ) then
        call inpcrb(1,iprint,cicrdn,maxatm,numatm,Tcord,ier)
      endif

      if ( ixnatm .ne. numatm ) then
        write(iprint,*)'ERROR> INPUT '
        write(iprint,*)'   NUMBER OF ATOMS IS DIFFERENT '
        write(iprint,*)'      IN TPL   : ',ixnatm
        write(iprint,*)'      IN CRD   : ',numatm
        write(iprint,*)' '
        ier = -1 ; return
      endif
      if ( ier .ne. 0 ) return
      allocate(cord(3,ixnatm))
      cord(1:3,1:ixnatm) = Tcord(1:3,1:ixnatm)
 
      ! CHECK LINEARITY OF TORSION
      call chklin(iyntor,iyptor(1:4,1:iyntor),iprint,ier)
      if ( ier .ne. 0 ) return
      call chklin(iynimp,iypimp(1:4,1:iynimp),iprint,ier)

!*****************************

      return
      end subroutine rdcord
 

!===================================================================


      subroutine inppsc(iprint,cipscn,cirefn,rwork,iwork,cwork,ier)

!*******************************************************************
!
!     READ POSITION RESTRAINTS
!     READ REFERENCE COORDINATE
!
!*******************************************************************
 
      use COMBAS ; use COMMIS

      implicit none

      ! Logical unit number for output
        integer(4),intent(in):: iprint
      ! File name of position restraints data
        character(80),intent(in):: cipscn
      ! File name of reference coordinate file
        character(80),intent(in):: cirefn
      ! Work arrays
        real(8),intent(inout):: rwork(maxatm)
        integer(4),intent(inout):: iwork(maxatm)
        character(*),intent(inout):: cwork(maxatm)
      ! Condition code
        integer(4),intent(out):: ier

      integer(4):: numatm
      real(8):: Tcord(3,maxatm)
 
!************************************************

      if ( cipscn.ne." " .and. cirefn.eq." " ) then
        write(iprint,*)'ERROR> INPUT '
        write(iprint,*)'   REFERENCE COORDINATE MUST BE READ'
        write(iprint,*)' '
        ier = -1 ; return
      else
        ier = 0
      endif
 
      ! 1) READ CONTROL DATA OF POSITION RESTRAINTS
      if ( cipscn .ne. " " ) then
        call psr(1,iprint,cipscn,ier)
        if ( ier .ne. 0 ) return
      endif

      ! 2) READ REFERENCE COORDINATE FOR POSITION CONSTARINTS
      if ( cirefn .ne. " " ) then
        write(iprint,*)' '
        write(iprint,*)'INFORMATION> INPUT '
        write(iprint,*)'     READ REFERECE STRUCTURE '
        call inppdb(1,iprint,cirefn,maxatm,numatm,cwork,cwork,iwork,   &
                    Tcord,rwork,rwork,ier)
        if ( ier .ne. 0 ) return
        allocate(fucord(3,ixnatm))
        fucord(1:3,1:ixnatm) = Tcord(1:3,1:ixnatm)
      endif

!******************************** 

      return
      end subroutine inppsc


!===================================================================


      subroutine tracen(ixnatm,iprint,fxmass,fxcbou)

!*******************************************************************
!
!     TRANSFORMATION
!       MASS-CENTER OF CURREBT STRUCTURE BECOMES ORIGIN
!       OF COORDINATE-SYSTEM
!
!*******************************************************************

      use COMMIS, only: fucord
      use COMCMMC, only: cord

      implicit none

      integer(4),intent(in):: ixnatm,iprint
      real(8),intent(in):: fxmass(ixnatm)
      real(8),intent(inout):: fxcbou(3)

      integer(4):: iatm
      real(8):: molwet,cencod(3)
 
!*****************************************

      molwet = 0.d0 ; cencod(1:3) = 0.d0

      molwet = sum(fxmass(1:ixnatm))
      do iatm = 1,ixnatm
        cencod(1) = cencod(1) + fxmass(iatm)*cord(1,iatm)
        cencod(2) = cencod(2) + fxmass(iatm)*cord(2,iatm)
        cencod(3) = cencod(3) + fxmass(iatm)*cord(3,iatm)
      enddo
      cencod(1:3) = cencod(1:3) / molwet

      cord(1,1:ixnatm) = cord(1,1:ixnatm) - cencod(1)
      cord(2,1:ixnatm) = cord(2,1:ixnatm) - cencod(2)
      cord(3,1:ixnatm) = cord(3,1:ixnatm) - cencod(3)
      fucord(1,1:ixnatm) = fucord(1,1:ixnatm) - cencod(1)
      fucord(2,1:ixnatm) = fucord(2,1:ixnatm) - cencod(2)
      fucord(3,1:ixnatm) = fucord(3,1:ixnatm) - cencod(3)
      fxcbou(1:3) = fxcbou(1:3) - cencod(1:3)

      write(iprint,*)' '
      write(iprint,*)'INFORMATION> INPUT '
      write(iprint,*)'      MASS-CENTER OF THIS SYSTEM '
      write(iprint,*)'      IN OLD COORDINATE SYSTEM '
      write(iprint,'(10x,3f15.7)')cencod(1:3)
      write(iprint,*)'      TOTAL MASS  ',molwet
      write(iprint,*)' '
      write(iprint,*)'      TRANSFORMATION IS DONE '
      write(iprint,*)' '

!********************************

      return
      end subroutine tracen


!==================================================================


      subroutine calint(iprint)

!*******************************************************************
!
!     1) SUM-UP 1-2 1-3 1-4 INTERACTIONS
!     2) SUM-UP TOTAL CHARGE
!     3) SUM-UP NUMBER OF FIXED ATOM PAIRS
!
!*******************************************************************
 
      use COMBAS ; use COMERG 

      implicit none

      integer(4),intent(in):: iprint

      integer(4):: i
      real(8):: totchg
 
!*******************************************

      totchg = sum(fxchrg(1:ixnatm))
      iyni12 = sum(ix14if(1:ixnatm,1))
      iyni13 = sum(ix14if(1:ixnatm,2))
      iyni14 = sum(ix14if(1:ixnatm,3))
      
      if ( ixnatm .ne. iynvar ) then
        iynfpt = ((ixnatm-iynvar)*(ixnatm-iynvar-1)) / 2

        iynfp2 = 0
        do i = 1,iynbnd
          if ( iytvar(iypbnd(1,i)).eq.0 .and.                          &
               iytvar(iypbnd(2,i)).eq.0 ) iynfp2 = iynfp2 + 1
        enddo

        iynfp3 = 0
        do i = 1,iynang
          if ( iytvar(iypang(1,i)).eq.0 .and.                          &
               iytvar(iypang(3,i)).eq.0 ) iynfp3 = iynfp3 + 1
        enddo

        iynfp4 = 0
        do i = 1,iyntor
          if ( iytvar(iyptor(1,i)).eq.0 .and.                          &
               iytvar(iyptor(4,i)).eq.0 .and.                          &
               iytnbf(i).eq.1 ) iynfp4 = iynfp4 + 1
        enddo

        iynfp5 = iynfpt - ( iynfp2 + iynfp3 + iynfp4 )

      else
        iynfpt = 0 ; iynfp2 = 0 ; iynfp3 = 0 ; iynfp4 = 0 ; iynfp5 = 0

      endif

      write(iprint,*)' '
      write(iprint,*)'INFORMATION> INPUT '
      write(iprint,'(a36,i10)')                                        &
           '    NUMBER OF 1-2 INTERACTIONS    : ',iyni12
      write(iprint,'(a36,i10)')                                        &
           '    NUMBER OF 1-3 INTERACTIONS    : ',iyni13
      write(iprint,'(a36,i10)')                                        &
           '    NUMBER OF 1-4 INTERACTIONS    : ',iyni14
      write(iprint,'(a36,f10.5)')                                      &
           '    TOTAL CHARGE (ELECTRON)       : ',totchg

      if ( iynvar .ne. ixnatm ) then
        write(iprint,'(a36,i10)')                                      &
           "    NUMBER OF FIXED ATOM PAIRS    : ",iynfpt
        write(iprint,'(a36,i10)')                                      &
           "              1-2                 : ",iynfp2
        write(iprint,'(a36,i10)')                                      &
           "              1-3                 : ",iynfp3
        write(iprint,'(a36,i10)')                                      &
           "              1-4                 : ",iynfp4
        write(iprint,'(a36,i10)')                                      &
           "              1-5                 : ",iynfp5
      endif

!*************************************

      return
      end subroutine calint


!===================================================================


      subroutine chkini

!*******************************************************************
!
!     MAKE IYFINI DATA
!
!       CONDITION OF IYFINI = 0
!         A) ONE CHAIN HAS ONLY ONE RESIDUE
!         B) THERE ARE ONLY FIXED ATOMS OR FREE ATOMS IN
!            ONE RESIDUE (ONE CHAIN)
!         C) THERE IS NO 1-5 INTERACTION IN ONE RESIDE (ONE CHAIN)
!         D) OPLS UNITED ATOM TYPE  (VDW ONLY)
!       UNDER THIS CONDITION , INTERACTION LIST OF EACH ATOM OF ONE
!       RESIDUE IS SAME WHEN RESIDUE-BASE CUTOFF IS DONE
!
!*******************************************************************
 
      use COMBAS ; use COMERG 

      implicit none

      integer(4):: flgmol(ixmolc)
      integer(4):: ires,ichnst,imol,iatmst,iatmen,iresst,iatm,numi14,  &
                   jatm,jint,numtot
 
!*******************************************

      allocate(iyfini(ixnres))

      ! D) OPLS UNITED ATOM TYPE ?
      if ( index(cynbpf(1),"OPLS") .eq. 0 ) then
        iyfini(1:ixnres) = 1 ; return
      endif

      ! <<<  CHECK MOLECULE INFORMATION  >>>
      ichnst = 1 ; flgmol(1:ixmolc) = 0
      DO_1 : do imol = 1,ixmolc
        if ( ichnst .eq. 1 ) then
          iatmst = 1
        else
          iatmst = ixcend(ichnst-1) + 1
        endif
        iatmen = ixcend(ichnst) ; iresst = ixares(iatmst)

        ! A) IS THERE ONLY ONE RESIDUE IN ONE CHAIN ?
        do iatm = iatmst,iatmen
          if ( iresst .ne. ixares(iatm) ) then
            flgmol(imol) = 1
            ichnst = ichnst + ixsqml(imol)
            cycle DO_1
          endif
        enddo

        ! C) IS THERE NO 1-5 INTERACTION IN ONE RESIDUE ?
        do iatm = iatmst,iatmen-1
          numi14 = sum(ix14if(iatm,1:3))
          DO_2 : do jatm = iatm+1,iatmen
            do jint = 1,numi14
              if ( jatm .eq. ix14lt(iatm,jint) ) cycle DO_2
            enddo
            flgmol(imol) = 1
            ichnst = ichnst + ixsqml(imol)
            cycle DO_1
          enddo DO_2
        enddo

        ichnst = ichnst + ixsqml(imol)
      enddo DO_1

      ! <<<  CHECK FIX/FREE ATOM MIXTURE AND MAKE IYFINI  >>>
      do ires = 1,ixnres
        iatmst = ixrstr(ires) ; iatmen = ixrend(ires)
        imol = ixamol(iatmst)
        if ( flgmol(imol) .eq. 1 ) then
          iyfini(ires) = 1
        ! B) ARE THERE ONLY FIXED ATOMS OR FREE ATOMS IN ONE RESIDUE ?
        else
          numtot = sum(iytvar(iatmst:iatmen))
          if ( numtot.eq.0 .or. numtot.eq.iatmen-iatmst-1 ) then
            iyfini(ires) = 0
          else
            iyfini(ires) = 1
          endif
        endif
      enddo

!*****************************************

      return
      end subroutine chkini
