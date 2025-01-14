!----------------------------------------
!
!     subroutine noe_calc1
!     subroutine noe_calc2
!
!----------------------------------------

      subroutine noe_calc1

!*****************************************
!
!     NOE intensity calculation
!                       (iterative part)
!
!*****************************************

      use COMVAL ; use COMIFN

      implicit none

      integer(4):: i,j,atm1,atm2
      real(8):: rtmp,tcod(3),ttcod(3),d(3)

!***********************************

       do i = 1,NhATM-1
         atm1 = hydID(i)
         tcod(1:3) = cod(1:3,atm1)

         do j = i+1,NhATM
           atm2 = hydID(j)
           ttcod(1:3) = cod(1:3,atm2)
           d(1:3) = tcod(1:3) - ttcod(1:3)
           d(1:3) = d(1:3) - bsize(1:3)*nint(d(1:3)*ibsize(1:3))
           rtmp = d(1)*d(1)+d(2)*d(2)+d(3)*d(3)
           rtmp = rtmp**(-3)
           NOEint(i,j) = rtmp*wfac + NOEint(i,j)
         enddo

       enddo

!***********************************

      return
      end subroutine noe_calc1


!==============================================


      subroutine noe_calc2

!*****************************************
!
!     NOE intensity calculation
!                           (final part)
!
!*****************************************

      use COMIFN ; use COMVAL

      implicit none

      ! OK count, satisfied the NOE distance
        integer(4):: icn1,icn2,icn3,icn4,icn5,icn6
        integer(4):: jcn1,jcn2,jcn3,jcn4,jcn5,jcn6
      ! REProduced NOE distances defined in the Experiment (for Pseudo)
        real(8),allocatable:: eNOErep(:),PeNOErep(:)
      ! Calculated NOE distances from the ensemble
        real(8),allocatable:: cNOE(:,:)
      ! Real calculated NOE intensities from the ensemble (not to consider Pseudo)
        real(8),allocatable:: rNOE(:,:)
      ! Experimental NOE Energy
        real(8),allocatable:: EnoeE(:,:)
      ! COorDinates for each NOE data
        real(8),allocatable:: codNOE(:,:,:)
      ! COorDinates for each NOE data from the Simulation ensemble
        real(8),allocatable:: ScodNOE(:,:)
      ! Geminal Hydrogen AToM NaMe list
        character(4),allocatable:: ghATMnm(:)
      ! RESidue number of Geminal Hydrogen 
        integer(4),allocatable:: GHiRES(:)
      ! AToM name of Geminal Hydrogen
        character(4),allocatable:: GHaATM(:)
      ! (Pre-)Tolerance of NOE distance constraints
        real(8),allocatable:: tole(:,:),Ptole(:)
      ! Same geminal hydrogen atom Check
        logical(1),allocatable:: sc(:)
      ! Number of Geminal Hydrogen AToMs
        integer(4):: NghATM
      ! Color
        real(8):: R,G,B
      ! Count for each NOE data
        integer(4),allocatable:: icn(:)
      ! Count for each NOE data from the simulation ensemble
        integer(4),allocatable:: jcn(:,:)
      ! Count of INTRA, SEQuencial, MEDiam, and LONG noe pairs
        integer(4):: intra,seq,med,long
      
      integer(4):: i,j,k,l,iH,jH,iL,jL
      character(4):: iATMnm,jATMnm,iRESnm
      real(8):: d

!***********************************
      ! Calculate the ensemble average NOE intensities
      NOEint(:,:) = NOEint(:,:) / sumW2

!*************
      ! Calc. NOE dist. for Experimental NOE reproduction result
      if ( input_NOE_file_type .ne. " " ) then
        write(6,'(2x,a)')"+ Exp. NOE reproduction (real) : "
        write(6,'(6x,a)')trim(project_name)//".rrnoe"
        open(unit=1,file=trim(d_proj)//".rrnoe",status="replace")

        write(6,'(2x,a)')"+ Exp. NOE reproduction (pseudo) : "
        write(6,'(6x,a)')trim(project_name)//".rpnoe"
        open(unit=2,file=trim(d_proj)//".rpnoe",status="replace")

        write(6,'(2x,a)')"+ "//                                        &
             "Reproduced Exp. NOE result for chimera (real): "
        write(6,'(6x,a)')trim(project_name)//"_rr.bild"
        open(unit=3,file=trim(d_proj)//"_rr.bild",status="replace")

        write(6,'(2x,a)')"+ "//                                        &
             "Reproduced Exp. NOE result for chimera (pseudo): "
        write(6,'(6x,a)')trim(project_name)//"_rp.bild"
        open(unit=4,file=trim(d_proj)//"_rp.bild",status="replace")

        write(6,'(2x,a)')"+ "//                                        &
         "Reproduced Exp. NOE result (long-range) for chimera (real) : "
        write(6,'(6x,a)')trim(project_name)//"_long.bild"
        open(unit=8,file=trim(d_proj)//"_long.bild",status="replace")


        write(1,*)""
        write(1,*)'pair   <noe>         UPb   <noe>-UPb   NOEenergy'
        write(1,*)""

        write(2,*)""
        write(2,*)'pair   <noe>         UPb   <noe>-UPb   NOEenergy'
        write(2,*)""
   
        allocate(eNOErep(NeNOE),PeNOErep(NeNOE),icn(NeNOE)) 
        eNOErep(:) = 0.d0 ; icn(:) = 0
        allocate(codNOE(3,2,NeNOE)) ; codNOE(:,:,:) = 0.d0
      
        ! Calculation of experimental NOE distances reproduction ratio
        do i = 1,NhATM-1
        iH = hydID(i) ; iATMnm = adjustl(ATMnm(iH))
        do j = i+1,NhATM
        jH = hydID(j) ; jATMnm = adjustl(ATMnm(jH))
          do k = 1,NeNOE
            iL = 4 ; jL = 4
            ! Residue number check
            if ( RESnum(iH).eq.iRES(1,k) .and.                         &
                 RESnum(jH).eq.iRES(2,k) ) then
              ! Atom name check
              iL = len_trim(aATM(1,k))
              if ( index(aATM(1,k),"#") .ne. 0 ) iL = iL - 1
              jL = len_trim(aATM(2,k))
              if ( index(aATM(2,k),"#") .ne. 0 ) jL = jL - 1

              if ( ( aATM(1,k)(1:iL).eq.iATMnm(1:iL) .and.             &
                     aATM(2,k)(1:jL).eq.jATMnm(1:jL) ) ) then
                 eNOErep(k) = eNOErep(k) + NOEint(i,j)
                 icn(k) = icn(k) + 1
                 codNOE(:,1,k) = codNOE(:,1,k) + Rcod(:,iH)
                 codNOE(:,2,k) = codNOE(:,2,k) + Rcod(:,jH)
              endif

            endif
            ! Residue number check
            if ( RESnum(iH).eq.iRES(2,k) .and.                     &
                 RESnum(jH).eq.iRES(1,K) ) then
              ! Atom name check
              iL = len_trim(aATM(2,k))
              if ( index(aATM(2,k),"#") .ne. 0 ) iL = iL - 1
              jL = len_trim(aATM(1,k))
              if ( index(aATM(1,k),"#") .ne. 0 ) jL = jL - 1

              if ( ( aATM(1,k)(1:jL).eq.jATMnm(1:jL) .and.             &
                     aATM(2,k)(1:iL).eq.iATMnm(1:iL) ) ) then
                eNOErep(k) = eNOErep(k) + NOEint(i,j)
                icn(k) = icn(k) + 1
                codNOE(:,1,k) = codNOE(:,1,k) + Rcod(:,jH)
                codNOE(:,2,k) = codNOE(:,2,k) + Rcod(:,iH)
              endif
            endif
            
          enddo

        enddo
        enddo
!        do i = 1,NeNOE-1
!          print*,i,eNOErep(i)
!        enddo
        do i = 1,NeNOE
          codNOE(:,:,i) = codNOE(:,:,i) / dble(icn(i))
        enddo
        PeNOErep(:) = ( eNOErep(:) / dble(icn(:)) )**(-1./6.)
        eNOErep(:) = eNOErep(:)**(-1./6.)
      endif
      
!*****************************
!     Calc. NOE dist. obtained from simulation results
      if ( output_NOE_file_type .ne. " " ) then
        
        ! Make geminal (hydrogen) name list
        allocate(ghATMnm(nATM))
        if ( pseudo_Hatom_type_for_output .eq. 0 ) then
          ghATMnm(:) = adjustl(ATMnm(:))
        else
          call geminal(pseudo_Hatom_type_for_output,nATM,RESnm,        &
                       ATMnm,ghATMnm)
        endif

        ! Remove not important atoms for marge of geminal hydrogen atoms
        allocate(sc(nATM),GHiRES(nATM),GHaATM(nATM))
        sc(:) = .true. ; NghATM = 0
        do i = 1,nATM
          ! Out of the residue range
          if ( RESnum(i).lt.fstRES .or. RESnum(i).gt.fnlRES ) then
            sc(i) = .false. ; cycle
          ! Removed specific hydrogen atoms
          elseif ( ghATMnm(i).eq."H1  "                                &
             .or.  ghATMnm(i).eq."H2  "                                &
             .or.  ghATMnm(i).eq."H3  "                                &
             .or. (RESnm(i)(1:3).eq."ARG".and.ghATMnm(i)(1:2).eq."HH") &
             .or. (RESnm(i)(1:3).eq."ARG".and.ghATMnm(i)(1:2).eq."HE") &
             .or. (RESnm(i)(1:3).eq."ASP".and.ghATMnm(i).eq."HD2 ")    &
             .or. (RESnm(i)(1:3).eq."CYS".and.ghATMnm(i).eq."HG  ")    &
             .or. (RESnm(i)(1:3).eq."GLU".and.ghATMnm(i).eq."HE2 ")    &
             .or. (RESnm(i)(1:3).eq."HIS".and.ghATMnm(i).eq."HD1 ")    &
             .or. (RESnm(i)(1:3).eq."HIS".and.ghATMnm(i).eq."HE2 ")    &
             .or. (RESnm(i)(1:3).eq."LYS".and.ghATMnm(i)(1:2).eq."HZ") &
             .or. (RESnm(i).eq."SER ".and.ghATMnm(i).eq."HG  ")        &
             .or. (RESnm(i).eq."THR ".and.ghATMnm(i).eq."HG1 ")        &
             .or. (RESnm(i).eq."TYR ".and.ghATMnm(i).eq."HH  ") ) then
            sc(i) = .false. ; cycle
          ! Selected hydrogen atoms
          elseif ( sc(i) .and. ghATMnm(i)(1:1).eq."H" ) then
            NghATM = NghATM + 1 
            GHiRES(NghATM) = RESnum(i)
            GHaATM(NghATM) = ghATMnm(i)
          ! duplicated atoms
          else
            sc(i) = .false. ; cycle
          endif

          do j = i+1,nATM
            if ( RESnum(j) .gt. GHiRES(NghATM) ) exit
            if ( (RESnum(j).eq.GHiRES(NghATM)) .and.            &
                 (ghATMnm(j).eq.GHaATM(NghATM)) ) sc(j) = .false.
          enddo
        enddo
        write(6,*)
        write(6,'(4x,a,i0)')                                           &
             "* Number of geminal H atoms in the range for NOE= ",NghATM
        write(6,*)

        ! Calculate NOE intensity simulated from the ensemble
        allocate(jcn(NghATM,NghATM)) ; jcn(:,:) = 0
        allocate(cNOE(NghATM,NghATM)) ; cNOE(:,:) = 0.d0
        allocate(rNOE(NghATM,NghATM)) ; rNOE(:,:) = 0.d0

        do i = 1,NhATM-1
          iH = hydID(i) ; iATMnm = ghATMnm(iH) ; iL = 0
          ! Search gh#
          do j = 1,NghATM
            if ( RESnum(iH).eq.GHiRES(j) .and. iATMnm.eq.GHaATM(j)) then
              iL = j ; exit
            endif
          enddo
          if ( iL .eq. 0 ) cycle

          do j = i+1,NhATM
            jH = hydID(j) ; jATMnm = ghATMnm(jH) ; jL = 0
            ! Search gh#
            do k = 1,NghATM
              if ( RESnum(jH).eq.GHiRES(k) .and. jATMnm.eq.GHaATM(k)) then
                jL = k ; exit
              endif
            enddo
            if ( jL .eq. 0 ) cycle

            cNOE(iL,jL) = cNOE(iL,jL) + NOEint(i,j)
            cNOE(jL,iL) = cNOE(jL,iL) + NOEint(i,j)
            jcn(iL,jL) = jcn(iL,jL) + 1
            jcn(jL,iL) = jcn(jL,iL) + 1
          enddo
        enddo
        
        do i = 1,NghATM
          write(6,'(8x,i,a3,i5,x,a4,x,a4)')                            &
                i," : ",GHiRES(i),aRES(GHiRES(i)),GHaATM(i)
        enddo
        write(6,*)

!***    Determination of NOE distance constraint corrections
        allocate(Ptole(NghATM)) ; Ptole(:) = 0.d0
        do i = 1,NghATM
          iATMnm = GHaATM(i) ; iRESnm = aRES(GHiRES(i))
          if ( index(iATMnm,"#") .ne. 0 ) then
            Ptole(i) = 1.d0
            if ( (iRESnm.eq."VAL " .and. iATMnm.eq."HG# ") .or.        &
                 (iRESnm.eq."LEU " .and. iATMnm.eq."HD# ") .or.        &
                 (iRESnm.eq."PHE " .and. iATMnm.eq."HD# ") .or.        &
                 (iRESnm.eq."PHE " .and. iATMnm.eq."HE# ") .or.        &
                 (iRESnm.eq."TYR " .and. iATMnm.eq."HD# ") .or.        &
                 (iRESnm.eq."TYR " .and. iATMnm.eq."HE# ") ) then
              Ptole(i) = 2.4d0
            endif
          endif
        enddo

        allocate(tole(NghATM,NghATM)) ; tole(:,:) = 0.d0
        do i = 1,NghATM-1
        do j = i+1,NghATM
          ! For intra-residue correction
          if ( GHiRES(i) .eq. GHiRES(j) ) then
            ! Provision for intra-residue correction
            iRESnm = aRES(GHiRES(i)) ; iATMnm = "    " ; jATMnm = "    "
            if ( (GHaATM(i).eq."H   ".or.GHaATM(i).eq."HA  ") .and.  &
                   index(GHaATM(j),"#").ne.0 )  then
               iATMnm = GHaATM(i) ; jATMnm = GHaATM(j)
            elseif ( (GHaATM(j).eq."H   ".or.GHaATM(j).eq."HA  ")    &
                    .and. index(GHaATM(i),"#").ne.0 ) then
               iATMnm = GHaATM(j) ; jATMnm = GHaATM(i)
            endif
            ! assign intra-residue correction
            if ( iATMnm .eq. "    " ) then
              tole(i,j) = Ptole(i) + Ptole(j)
            else
              if ( iRESnm.eq."ALA " .and. iATMnm.eq."H   " .and.       &
                   jATMnm.eq."HB# ") then
                tole(i,j) = 0.6d0 ; cycle
              elseif ( iRESnm.eq."VAL " .and. iATMnm.eq."H   " .and.   &
                       jATMnm.eq."HG# " ) then
                tole(i,j) = 1.7d0 ; cycle
              elseif ( iRESnm.eq."VAL " .and. iATMnm.eq."HA  " .and.   &
                       jATMnm.eq."HG1#" .or. jATMnm.eq."HG2#" ) then
                tole(i,j) = 0.6d0 ; cycle
              elseif ( iRESnm.eq."VAL " .and. iATMnm.eq."HA  " .and.   &
                       jATMnm.eq."HG# " ) then
                tole(i,j) = 0.6d0 ; cycle
              elseif ( iRESnm.eq."ILE " .and. iATMnm.eq."HA  " .and.   &
                      (jATMnm.eq."HG1#" .or. jATMnm.eq."HG2#") ) then
                tole(i,j) = 0.6d0 ; cycle
              elseif ( iRESnm.eq."LEU " .and. iATMnm.eq."H   " .and.   &
                       jATMnm.eq."HB# " ) then
                tole(i,j) = 0.6d0 ; cycle
              elseif ( iRESnm.eq."LEU " .and. iATMnm.eq."HA  " .and.   &
                       jATMnm.eq."HD# " ) then
                tole(i,j) = 1.7d0 ; cycle
              elseif ( iRESnm.eq."SER " .and. iATMnm.eq."H   " .and.   &
                       jATMnm.eq."HB# " ) then
                tole(i,j) = 0.6d0 ; cycle
              elseif ( iRESnm.eq."THR " .and. iATMnm.eq."HA  " .and.   &
                       jATMnm.eq."HG2#" ) then
                tole(i,j) = 0.6d0 ; cycle
              elseif ( (iRESnm(1:3).eq."ASP" .or. iRESnm.eq."ASN ")    &
                  .and. iATMnm.eq."HA  " .and. jATMnm.eq."HB# " ) then
                tole(i,j) = 0.6d0
              elseif ( (iRESnm(1:3).eq."GLU" .or. iRESnm.eq."GLN ")    &
                  .and. iATMnm.eq."H   " .and. jATMnm.eq."HB# " ) then
                tole(i,j) = 0.6d0 ; cycle
              elseif ( (iRESnm(1:3).eq."GLU" .or. iRESnm.eq."GLN " )   &
                  .and. iATMnm.eq."HA  " .and. jATMnm.eq."HG# " ) then
                tole(i,j) = 0.6d0 ; cycle
              elseif ( iRESnm(1:3).eq."LYS" .and. iATMnm.eq."H   "     &
                  .and. jATMnm.eq."HB# " ) then
                tole(i,j) = 0.6d0 ; cycle
              elseif ( iRESnm(1:3).eq."LYS" .and. iATMnm.eq."HA  "     &
                  .and. jATMnm.eq."HG# " ) then
                tole(i,j) = 0.6d0 ; cycle
              elseif ( iRESnm(1:3).eq."ARG" .and. iATMnm.eq."H   "     &
                  .and. jATMnm.eq."HB# " ) then
                tole(i,j) = 0.6d0 ; cycle
              elseif ( iRESnm(1:3).eq."ARG" .and. iATMnm.eq."HA  "     &
                  .and. jATMnm.eq."HG# " ) then
                tole(i,j) = 0.6d0 ; cycle
              elseif ( iRESnm.eq."MET " .and. iATMnm.eq."H   " .and.   &
                       jATMnm.eq."HB# " ) then
                tole(i,j) = 0.6d0 ; cycle
              elseif ( iRESnm.eq."MET " .and. iATMnm.eq."HA  " .and.   &
                       jATMnm.eq."HG# " ) then
                tole(i,j) = 0.6d0 ; cycle
              elseif ( iRESnm(1:3).eq."CYS" .and. iATMnm.eq."H   "     &
                  .and. jATMnm.eq."HB# " ) then
                tole(i,j) = 0.6d0 ; cycle
              elseif ( (iRESnm.eq."PHE " .or. iRESnm.eq."TYR ") .and.  &
                       iATMnm.eq."H   " .and. jATMnm.eq."HB# " ) then
                tole(i,j) = 0.6d0 ; cycle
              elseif ( (iRESnm.eq."PHE " .or. iRESnm.eq."TYR ") .and.  &
                       iATMnm.eq."HA  " .and.                          &
                       (jATMnm.eq."HD# " .or. jATMnm.eq."HE# ") ) then
                tole(i,j) = 2.d0 ; cycle
              elseif ( iRESnm(1:3).eq."HIS" .and. iATMnm.eq."H   "     &
                  .and. jATMnm.eq."HB# " ) then
                tole(i,j) = 0.6d0 ; cycle
              elseif ( iRESnm.eq."TRP " .and. iATMnm.eq."H   " .and.   &
                       jATMnm.eq."HB# " ) then
                tole(i,j) = 0.6d0 ; cycle
              else
                tole(i,j) = Ptole(i) + Ptole(j)
              endif
            endif
          ! For long-range correction
          else
            tole(i,j) = Ptole(i) + Ptole(j)
          endif
        enddo
        enddo

        rNOE(:,:) = cNOE(:,:)
        do i = 1,NghATM-1
        do j = i+1,NghATM
          cNOE(i,j) = cNOE(i,j) / dble(jcn(i,j))
        enddo
        enddo
        
        cNOE(:,:) = cNOE(:,:)**(-1./6.)
       
!       Output NOE restraint list obtainted from simulated ensemble
        if ( output_NOE_file_type .eq. "XPLOR" ) then
          call Xplor_noeout(d_proj,NghATM,cNOE,GHiRES,GHaATM,tole,rNOE)
        endif
 
        ! Make coordinates for calc. NOE
        jcn(:,:) = 0
        allocate(ScodNOE(3,NghATM)) ; ScodNOE(:,:) = 0.d0
        do i = 1,NhATM
          iH = hydID(i) ; iATMnm = ghATMnm(iH) ; iL = 0
          ! Search gh#
          do j = 1,NghATM
            if ( RESnum(iH).eq.GHiRES(j) .and. iATMnm.eq.GHaATM(j)) then
              iL = j ; exit
            endif
          enddo
          if ( iL .eq. 0 ) cycle
          ScodNOE(:,iL) = ScodNOE(:,iL) + Rcod(:,iH)
          jcn(iL,1) = jcn(iL,1) + 1
        enddo
        do i = 1,NghATM
          ScodNOE(:,i) = ScodNOE(:,i) / dble(jcn(i,1))
        enddo

      endif

      NOEint(:,:) = NOEint(:,:)**(-1./6.)

!*************************************
      
      ! Output experimental NOE reproduction result
      if ( input_NOE_file_type .ne. " " ) then
        allocate(EnoeE(NeNOE,2))
        icn1 = 0 ; icn2 = 0 ; icn3 = 0 ; icn4 = 0 ; icn5 = 0 ; icn6 = 0
        jcn1 = 0 ; jcn2 = 0 ; jcn3 = 0 ; jcn4 = 0 ; jcn5 = 0 ; jcn6 = 0
        do i = 1,NeNOE
          ! NOE energy calculation (real)
          if ( eNOErep(i) .gt. UPb(i) ) then
            EnoeE(i,1) = (eNOErep(i) - UPb(i)) * (eNOErep(i) - UPb(i))
          else
            EnoeE(i,1) = 0.d0
          endif
          ! NOE energy calculation (pseudo)
          if ( PeNOErep(i) .gt. UPb(i) ) then
            EnoeE(i,2) = (PeNOErep(i) - UPb(i)) * (PeNOErep(i) - UPb(i))
          else
            EnoeE(i,2) = 0.d0
          endif
        
          ! Output (real)
          if ( eNOErep(i) .le. UPb(i) ) then
            write(1,200)i,eNOErep(i),UPb(i),eNOErep(i)-UPb(i),         &
                 EnoeE(i,1),"OK        "
            icn1 = icn1 + 1
          elseif ( eNOErep(i) .le. UPb(i)+0.5d0 ) then
            write(1,200)i,eNOErep(i),UPb(i),eNOErep(i)-UPb(i),         &
                 EnoeE(i,1),"NG(< 0.5A)"
            icn2 = icn2 + 1
          elseif ( eNOErep(i) .le. UPb(i)+1.d0 ) then
            write(1,200)i,eNOErep(i),UPb(i),eNOErep(i)-UPb(i),         &
                 EnoeE(i,1),"NG(< 1.0A)"
            icn3 = icn3 + 1
          elseif ( eNOErep(i) .le. UPb(i)+2.d0 ) then
            write(1,200)i,eNOErep(i),UPb(i),eNOErep(i)-UPb(i),         &
                 EnoeE(i,1),"NG(< 2.0A)"
            icn4 = icn4 + 1
          elseif ( eNOErep(i) .le. UPb(i)+3.d0 ) then
            write(1,200)i,eNOErep(i),UPb(i),eNOErep(i)-UPb(i),         &
                 EnoeE(i,1),"NG(< 3.0A)"
            icn5 = icn5 + 1
          else
            write(1,200)i,eNOErep(i),UPb(i),eNOErep(i)-UPb(i),         &
                 EnoeE(i,1),"NG(> 3.0A)"
            icn6 = icn6 + 1
          endif
        
          ! Output (pseudo)
          if ( PeNOErep(i) .le. UPb(i) ) then
            write(2,200)i,PeNOErep(i),UPb(i),                         &
                 PeNOErep(i)-UPb(i),EnoeE(i,2),"OK        "
            jcn1 = jcn1 + 1
          elseif ( PeNOErep(i) .le. UPb(i)+0.5d0 ) then
            write(2,200)i,PeNOErep(i),UPb(i),                          &
                 PeNOErep(i)-UPb(i),EnoeE(i,2),"NG(< 0.5A)"
            jcn2 = jcn2 + 1
          elseif ( PeNOErep(i) .le. UPb(i)+1.d0 ) then
            write(2,200)i,PeNOErep(i),UPb(i),                          &
                 PeNOErep(i)-UPb(i),EnoeE(i,2),"NG(< 1.0A)"
            jcn3 = jcn3 + 1
          elseif ( PeNOErep(i) .le. UPb(i)+2.d0 ) then
            write(2,200)i,PeNOErep(i),UPb(i),                          &
                 PeNOErep(i)-UPb(i),EnoeE(i,2),"NG(< 2.0A)"
            jcn4 = jcn4 + 1
          elseif ( PeNOErep(i) .le. UPb(i)+3.d0 ) then
            write(2,200)i,PeNOErep(i),UPb(i),                          &
                 PeNOErep(i)-UPb(i),EnoeE(i,2),"NG(< 3.0A)"
            jcn5 = jcn5 + 1
          else
            write(2,200)i,PeNOErep(i),UPb(i),                          &
                 PeNOErep(i)-UPb(i),EnoeE(i,2),"NG(> 3.0A)"
            jcn6 = jcn6 + 1
          endif
200       format(i5,2x,f6.3,2x,f10.3,2x,2(e10.3,2x),a10)

          ! Output for chimera (real)
          if ( EnoeE(i,1) .eq. 0.d0 ) then
            R = 0.d0 ; G = 0.d0 ; B = 1.d0
          elseif ( EnoeE(i,1) .le. 0.25d0 ) then
            R = 0.d0 ; G = 4.d0*EnoeE(i,1) ; B = 1.d0
          elseif ( EnoeE(i,1) .le. 1.d0 ) then
            R = 0.d0 ; G = 1.d0 ; B = (1.d0-EnoeE(i,1))*4.d0
          elseif ( EnoeE(i,1) .le. 4.d0 ) then
            R = (EnoeE(i,1)-1.d0)/3.d0 ; G = 1.d0 ; B = 0.d0
          elseif ( EnoeE(i,1) .le. 9.d0 ) then
            R = 1.d0 ; G = (9.d0-EnoeE(i,1))/5.d0 ; B = 0.d0
          else
            R = 1.d0 ; G = 0.d0 ; B = 0.d0
          endif
!          if ( EnoeE(i,1) .gt. 1.d0 ) then
          write(3,'(a7,x,3(f8.3,x))')".color ",R,G,B
          write(3,'(a10,x,2(3(f8.3,x),x),x,f8.3)')                     &
               ".cylinder ",codNOE(:,1,i),codNOE(:,2,i),fnlRES*0.002
!          endif
          if ( abs(iRES(1,i)-iRES(2,i)) .gt. 3 ) then
            write(8,'(a7,x,3(f8.3,x))')".color ",R,G,B
            write(8,'(a10,x,2(3(f8.3,x),x),x,f8.3)')                   &
                 ".cylinder ",codNOE(:,1,i),codNOE(:,2,i),fnlRES*0.002
          endif

          ! Output for chimera (pseudo)
          if ( EnoeE(i,2) .eq. 0.d0 ) then
            R = 0.d0 ; G = 0.d0 ; B = 1.d0
          elseif ( EnoeE(i,2) .le. 0.25d0 ) then
            R = 0.d0 ; G = 4.d0*EnoeE(i,2) ; B = 1.d0
          elseif ( EnoeE(i,2) .le. 1.d0 ) then
            R = 0.d0 ; G = 1.d0 ; B = (1.d0-EnoeE(i,2))*4.d0
          elseif ( EnoeE(i,2) .le. 4.d0 ) then
            R = (EnoeE(i,2)-1.d0)/3.d0 ; G = 1.d0 ; B = 0.d0
          elseif ( EnoeE(i,2) .le. 9.d0 ) then
            R = 1.d0 ; G = (9.d0-EnoeE(i,2))/5.d0 ; B = 0.d0
          else
            R = 1.d0 ; G = 0.d0 ; B = 0.d0
          endif
          write(4,'(a7,x,3(f8.3,x))')".color ",R,G,B
          write(4,'(a10,x,2(3(f8.3,x),x),x,f8.3)')                     &
               ".cylinder ",codNOE(:,1,i),codNOE(:,2,i),fnlRES*0.002
        enddo

        write(1,*)""
        write(1,*)"**************************************************"
        write(1,201)"OK         = ",real(icn1)/real(NeNOE)*100.d0," (%)"
        write(1,201)"NG(< 0.5A) = ",real(icn2)/real(NeNOE)*100.d0," (%)"
        write(1,201)"NG(< 1.0A) = ",real(icn3)/real(NeNOE)*100.d0," (%)"
        write(1,201)"NG(< 2.0A) = ",real(icn4)/real(NeNOE)*100.d0," (%)"
        write(1,201)"NG(< 3.0A) = ",real(icn5)/real(NeNOE)*100.d0," (%)"
        write(1,201)"NG(> 3.0A) = ",real(icn6)/real(NeNOE)*100.d0," (%)"
        write(1,*)""
        write(1,'(a13,f14.3)')"NOE energy = ",sum(EnoeE(:,1))
        write(1,'(a13,f14.3,a4)')"Max error = ",                       &
             sqrt(maxval(EnoeE(:,1)))," (A)"
        write(1,'(a13,i)')"Pair numb. = ",maxloc(EnoeE(:,1))
        write(1,*)"***********************************************"
        write(1,*)""
        write(1,'(a)')"       Experimental NOE constraints list"
        write(1,*)
        write(1,*)"#expNOE #res1 atomN1  #res2 atomN2  type      UPb"
        write(1,*)

        write(2,*)""
        write(2,*)"**************************************************"
        write(2,201)"OK         = ",real(jcn1)/real(NeNOE)*100.d0," (%)"
        write(2,201)"NG(< 0.5A) = ",real(jcn2)/real(NeNOE)*100.d0," (%)"
        write(2,201)"NG(< 1.0A) = ",real(jcn3)/real(NeNOE)*100.d0," (%)"
        write(2,201)"NG(< 2.0A) = ",real(jcn4)/real(NeNOE)*100.d0," (%)"
        write(2,201)"NG(< 3.0A) = ",real(jcn5)/real(NeNOE)*100.d0," (%)"
        write(2,201)"NG(> 3.0A) = ",real(jcn6)/real(NeNOE)*100.d0," (%)"
        write(2,*)""
        write(2,'(a13,f14.3)')"NOE energy = ",sum(EnoeE(:,2))
        write(2,'(a13,f14.3,a4)')"Max error = ",                       &
             sqrt(maxval(EnoeE(:,2)))," (A)"
        write(2,'(a13,i)')"Pair numb. = ",maxloc(EnoeE(:,2))
        write(2,*)"***********************************************"
        write(2,*)""
        write(2,'(a)')"       Experimental NOE constraints list"
        write(2,*)
        write(2,*)"#expNOE #res1 atomN1  #res2 atomN2  type      UPb"
        write(2,*)

201     format(a13,f14.3,a4)
      
        intra = 0 ; seq = 0 ; med = 0 ; long = 0
        do i = 1,NeNOE
          j = abs(iRES(1,i)-iRES(2,i))
          if ( j .eq. 0 ) then
            write(1,202)i,iRES(1,i),aRES(iRES(1,i)),aATM(1,i),iRES(2,i)&
                 ,aRES(iRES(2,i)),aATM(2,i)," intra",UPb(i)
            write(2,202)i,iRES(1,i),aRES(iRES(1,i)),aATM(1,i),iRES(2,i)&
                 ,aRES(iRES(2,i)),aATM(2,i)," intra",UPb(i)
            intra = intra + 1
          elseif ( j .le. 1 ) then
            write(1,202)i,iRES(1,i),aRES(iRES(1,i)),aATM(1,i),iRES(2,i)&
                 ,aRES(iRES(2,i)),aATM(2,i),"   seq",UPb(i)
            write(2,202)i,iRES(1,i),aRES(iRES(1,i)),aATM(1,i),iRES(2,i)&
                 ,aRES(iRES(2,i)),aATM(2,i),"   seq",UPb(i)
            seq = seq + 1
          elseif ( j .le. 3 ) then
            write(1,202)i,iRES(1,i),aRES(iRES(1,i)),aATM(1,i),iRES(2,i)&
                 ,aRES(iRES(2,i)),aATM(2,i),"   med",UPb(i)
            write(2,202)i,iRES(1,i),aRES(iRES(1,i)),aATM(1,i),iRES(2,i)&
                 ,aRES(iRES(2,i)),aATM(2,i),"   med",UPb(i)
            med = med + 1
          else
            write(1,202)i,iRES(1,i),aRES(iRES(1,i)),aATM(1,i),iRES(2,i)&
                 ,aRES(iRES(2,i)),aATM(2,i),"  long",UPb(i)
            write(2,202)i,iRES(1,i),aRES(iRES(1,i)),aATM(1,i),iRES(2,i)&
                 ,aRES(iRES(2,i)),aATM(2,i),"  long",UPb(i)
            long = long + 1
          endif
        enddo
202     format(i5,": ",i4,x,a4,x,a4," - ",i4,x,a4,x,a4,a6,2x,f6.3)

        write(1,*)""
        write(1,*)"***********************************************"
        write(1,'(4x,a,i)')"  ALL =",NeNOE
        write(1,'(4x,a,i)')"intra (same res.)             = ",intra
        write(1,'(4x,a,i)')"  seq (adjacent res.)         = ",seq
        write(1,'(4x,a,i)')"  med (2 or 3 neighbor res.)  = ",med
        write(1,'(4x,a,i)')" long (More than 3 res.)      = ",long

        write(2,*)""
        write(2,*)"***********************************************"
        write(2,'(4x,a,i)')"  ALL =",NeNOE
        write(2,'(4x,a,i)')"intra (same res.)             = ",intra
        write(2,'(4x,a,i)')"  seq (adjacent res.)         = ",seq
        write(2,'(4x,a,i)')"  med (2 or 3 neighbor res.)  = ",med
        write(2,'(4x,a,i)')" long (More than 3 res.)      = ",long

        do i = 1,4
          close(i)
        enddo
        close(8)

      endif

!*******************************
      ! Output calculated NOE constraints from the ensemble
      if ( output_NOE_file_type .ne. " " ) then
        write(6,'(2x,a)')                                              &
             "+ "//"NOE from simulated ensemble result for chimera : "
        write(6,'(6x,a)')trim(project_name)//"_c.bild"
        open(unit=1,file=trim(d_proj)//"_c.bild",status="replace")

        do i = 1,NghATM-1
        do j = i+1,NghATM

          if ( cNOE(i,j) .lt. 5.d0 ) then

            ! Color rendaring
            if ( cNOE(i,j) .lt. 2.5d0 ) then
              R = 0.d0 ; G = 0.d0 ; B = 1.d0
            elseif ( cNOE(i,j) .lt. 3.d0 ) then
              R = 0.d0 ; G = 1.d0 ; B = 1.d0
            elseif ( cNOE(i,j) .lt. 4.d0 ) then
              R = 0.d0 ; G = 1.d0 ; B = 0.d0
            elseif ( cNOE(i,j) .lt. 5.d0 ) then
              R = 1.d0 ; G = 1.d0 ; B = 0.d0
            else
              cycle
            endif

            write(1,'(a7,x,3(f8.3,x))')".color ",R,G,B
            write(1,'(a10,x,2(3(f8.3,x),x),x,f8.3)')                   &
                 ".cylinder ",ScodNOE(:,i),ScodNOE(:,j),fnlRES*0.002
          endif

        enddo
        enddo
      endif
      
!************************

      return
      end subroutine noe_calc2
