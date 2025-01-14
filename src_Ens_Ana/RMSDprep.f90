
      subroutine RMSDprep(icn,ichance2)

      use COMIFN ; use COMVAL

      implicit none

      integer(4),intent(in):: icn,ichance2

      integer(4):: i
      integer(4),save:: icnRMSD = 0
      real(8):: RMSD,rtmp

!**************************************

      ! all-against-all
      if ( RMSD_mode .eq. "multi" ) then
        !! same ensemble
        if ( (len(trim(input_binary_list)).ne.0 .and.                  &
            (input_binary_list.eq.input_reference_binary_list)) .or.   &
             (len(trim(input_PDB_list)).ne.0 .and.                     &
            (input_PDB_list.eq.input_reference_PDB_list)) ) then
          do i = ichance2+1,nREF
            call RMSDcalc(nATM,nRMSDf,nRMSD,RcodRMSDf(:,:,i),          &
              RcodRMSD(:,:,i),cod,iATMrmsdF,iATMrmsd,RMSD,rotation_flag)
            write(urm,'(2(i0,x),2(f8.3,x))')icn,conf_n(i),RMSD,        &
                                          wfac*conf_w(i)
            sumRMSD(ichance2) = sumRMSD(ichance2) + RMSD*conf_w(i)
            sumRMSD(i) = sumRMSD(i) + RMSD*wfac
            rtmp = RMSD*RMSD
            sumRMSD2(ichance2) = sumRMSD2(ichance2)+rtmp*conf_w(i)
            sumRMSD2(i) = sumRMSD2(i) + rtmp*wfac
            sWrmsd(ichance2) = sWrmsd(ichance2) + conf_w(i)
            sWrmsd(i) = sWrmsd(i) + wfac
          enddo

        !! query v.s. reference ensemble
        else 
          do i = 1,nREF
            call RMSDcalc(nATM,nRMSDf,nRMSD,RcodRMSDf(:,:,i),          &
              RcodRMSD(:,:,i),cod,iATMrmsdF,iATMrmsd,RMSD,rotation_flag)
            if ( nREF .eq. 1 ) then
              write(urm,'(i0,2(x,f))')icn,RMSD,wfac
            else
              write(urm,'(i0,x,i0,f8.3,x,f)')icn,i,RMSD,wfac
            endif
            sumRMSD(i) = sumRMSD(i) + RMSD * wfac
            sumRMSD2(i) = sumRMSD2(i) + RMSD*RMSD * wfac
            sWrmsd(i) = sWrmsd(i) + wfac
          enddo
        endif

      ! average structure
      elseif ( RMSD_mode .eq. "average" ) then
        call RMSDcalc(nATM,nRMSDf,nRMSD,aveRcodF,aveRcod,cod,iATMrmsdF,&
          iATMrmsd,RMSD,rotation_flag)
        write(urm,'(i0,2(x,f))')icn,RMSD,wfac
        sumRMSD(1) = sumRMSD(1) + RMSD * wfac
        sumRMSD2(1) = sumRMSD2(1) + RMSD*RMSD * wfac
        sWrmsd(1) = sWrmsd(1) + wfac
      endif

      ! Time series RMSD
      if ( interval_for_timeseries_RMSD .ne. 0 ) then
        icnRMSD = icnRMSD + 1
        if ( icn .gt. interval_for_timeseries_RMSD ) then
          call RMSDcalc(nATM,nRMSDf,nRMSD,codTSRMSDf(:,:,icnRMSD),     &
            codTSRMSD(:,:,icnRMSD),cod,iATMrmsdF,iATMrmsd,RMSD,.false.)
          write(utr,'(i0,a3,i0,x,f8.3)')                               &
            icn-interval_for_timeseries_RMSD," - ",icn,RMSD
        endif
        codTSRMSDf(1:3,1:nRMSDf,icnRMSD) = cod(1:3,iATMrmsdF(1:nRMSDf))
        codTSRMSD(1:3,1:nRMSD,icnRMSD) = cod(1:3,iATMrmsd(1:nRMSD))
        if ( icnRMSD .ge. interval_for_timeseries_RMSD ) icnRMSD = 0
      endif

!********************************

      return
      end subroutine RMSDprep


!=================================================================


      subroutine output_RMSD_log()

      use COMVAL ; use COMIFN,only: loglvl

      implicit none

      integer(4):: i,t,r
      character(1):: a,b

!************************************************************

      write(6,'(4x,a,i0)')                                             &
        "N of atoms used for fitting in RMSD calc. : ",nRMSDf
      do i = 1,nRMSDf
        t = iATMrmsdF(i) ; r = iATMrmsdRF(i)
        a = adjustl(ATMnm(t)) ; b = adjustl(ATMnmR(r))
        if ( a .ne. b ) then
           write(6,*)"!! WARNING !!"
           write(6,*)"RMSD fitting has been done on different atom "// &
                    "types"
           write(6,'(2x,i8,a3,3(i0,x,a,x),a3,3(i0,x,a,x))')i," : ",    &
            nCHN(t),trim(adjustl(CHN(t))),RESnum(t),                   &
            trim(adjustl(RESnm(t))),ATMnum(t),trim(adjustl(ATMnm(t))), &
            " - ",                                                     &
            nCHNr(r),trim(adjustl(CHNr(r))),RESnumR(r),                &
            trim(adjustl(RESnmR(r))),ATMnumR(r),trim(adjustl(ATMnmR(r)))
        endif
      enddo
      if ( loglvl.eq."l" .or. nRMSDf.le.6 ) then
        do i = 1,nRMSDf
          t = iATMrmsdF(i) ; r = iATMrmsdRF(i)
          write(6,'(2x,i8,a3,3(i0,x,a,x),a3,3(i0,x,a,x))')i," : ",     &
            nCHN(t),trim(adjustl(CHN(t))),RESnum(t),                   &
            trim(adjustl(RESnm(t))),ATMnum(t),trim(adjustl(ATMnm(t))), &
            " - ",                                                     &
            nCHNr(r),trim(adjustl(CHNr(r))),RESnumR(r),                &
            trim(adjustl(RESnmR(r))),ATMnumR(r),trim(adjustl(ATMnmR(r)))
        enddo
      elseif ( loglvl .eq. "s" ) then
        do i = 1,3
          t = iATMrmsdF(i) ; r = iATMrmsdRF(i)
          write(6,'(2x,i8,a3,3(i0,x,a,x),a3,3(i0,x,a,x))')i," : ",     &
            nCHN(t),trim(adjustl(CHN(t))),RESnum(t),                   &
            trim(adjustl(RESnm(t))),ATMnum(t),trim(adjustl(ATMnm(t))), &
            " - ",                                                     &
            nCHNr(r),trim(adjustl(CHNr(r))),RESnumR(r),                &
            trim(adjustl(RESnmR(r))),ATMnumR(r),trim(adjustl(ATMnmR(r)))
        enddo
        write(6,*)"              ..."
        do i = nRMSDf-2,nRMSDf
          t = iATMrmsdF(i) ; r = iATMrmsdRF(i)
          write(6,'(2x,i8,a3,3(i0,x,a,x),a3,3(i0,x,a,x))')i," : ",     &
            nCHN(t),trim(adjustl(CHN(t))),RESnum(t),                   &
            trim(adjustl(RESnm(t))),ATMnum(t),trim(adjustl(ATMnm(t))), &
            " - ",                                                     &
            nCHNr(r),trim(adjustl(CHNr(r))),RESnumR(r),                &
            trim(adjustl(RESnmR(r))),ATMnumR(r),trim(adjustl(ATMnmR(r)))
        enddo
      endif

      write(6,'(4x,a,i0)')"N of atoms used for RMSD calc. : ",nRMSD
      do i = 1,nRMSD
        t = iATMrmsd(i) ; r = iATMrmsdR(i)
        a = adjustl(ATMnm(t)) ; b = adjustl(ATMnmR(r))
        if ( a .ne. b ) then
           write(6,*)"!! WARNING !!"
           write(6,*)"RMSD calc. has been done on different atom types"
           write(6,'(2x,i8,a3,3(i0,x,a,x),a3,3(i0,x,a,x))')i," : ",    &
            nCHN(t),trim(adjustl(CHN(t))),RESnum(t),                   &
            trim(adjustl(RESnm(t))),ATMnum(t),trim(adjustl(ATMnm(t))), &
            " - ",                                                     &
            nCHNr(r),trim(adjustl(CHNr(r))),RESnumR(r),                &
            trim(adjustl(RESnmR(r))),ATMnumR(r),trim(adjustl(ATMnmR(r)))
        endif
      enddo
      if ( loglvl.eq."l" .or. nRMSDf.le.6 ) then
        do i = 1,nRMSD
          t = iATMrmsd(i) ; r = iATMrmsdR(i)
          write(6,'(2x,i8,a3,3(i0,x,a,x),a3,3(i0,x,a,x))')i," : ",     &
            nCHN(t),trim(adjustl(CHN(t))),RESnum(t),                   &
            trim(adjustl(RESnm(t))),ATMnum(t),trim(adjustl(ATMnm(t))), &
            " - ",                                                     &
            nCHNr(r),trim(adjustl(CHNr(r))),RESnumR(r),                &
            trim(adjustl(RESnmR(r))),ATMnumR(r),trim(adjustl(ATMnmR(r)))
        enddo
      elseif ( loglvl .eq. "s" ) then
        do i = 1,3
          t = iATMrmsd(i) ; r = iATMrmsdR(i)
          write(6,'(2x,i8,a3,3(i0,x,a,x),a3,3(i0,x,a,x))')i," : ",     &
            nCHN(t),trim(adjustl(CHN(t))),RESnum(t),                   &
            trim(adjustl(RESnm(t))),ATMnum(t),trim(adjustl(ATMnm(t))), &
            " - ",                                                     &
            nCHNr(r),trim(adjustl(CHNr(r))),RESnumR(r),                &
            trim(adjustl(RESnmR(r))),ATMnumR(r),trim(adjustl(ATMnmR(r)))
        enddo
        write(6,*)"              ..."
        do i = nRMSD-2,nRMSD
          t = iATMrmsd(i) ; r = iATMrmsdR(i)
          write(6,'(2x,i8,a3,3(i0,x,a,x),a3,3(i0,x,a,x))')i," : ",     &
            nCHN(t),trim(adjustl(CHN(t))),RESnum(t),                   &
            trim(adjustl(RESnm(t))),ATMnum(t),trim(adjustl(ATMnm(t))), &
            " - ",                                                     &
            nCHNr(r),trim(adjustl(CHNr(r))),RESnumR(r),                &
            trim(adjustl(RESnmR(r))),ATMnumR(r),trim(adjustl(ATMnmR(r)))
        enddo
      endif

      return
      end subroutine output_RMSD_log


!=================================================================


      subroutine RMSD_init

      use COMIFN ; use COMVAL

      implicit none

!***************************************

      allocate(sumRMSD(nREF),sumRMSD2(nREF),sWrmsd(nREF))
      sumRMSD = 0.d0 ; sumRMSD2 = 0.d0 ; sWrmsd = 0.d0
      open(unit=urm,file=trim(d_proj)//".rmsd",status="replace")
      if ( RMSD_mode .eq. "multi" ) then
        if ( input_PDB_list .eq. input_reference_PDB_list ) then
          write(urm,*)"# all-against-all RMSDs among the ensemble"
        else
          write(urm,*)"# all-against-all RMSDs among query and "//     &
                      "Ref. Ensemble "
        endif
      elseif ( RMSD_mode .eq. "average" ) then
        write(urm,*)"# RMSDs among query conf. and average str. of"//  &
                    " Ref. conf."
      endif
      write(urm,*)

      if ( interval_for_timeseries_RMSD .ne. 0 ) then
        open(unit=utr,file=trim(d_proj)//".tsrmsd",status="replace")
        allocate(codTSRMSDf(3,nRMSDf,interval_for_timeseries_RMSD))
        allocate(codTSRMSD(3,nRMSD,interval_for_timeseries_RMSD))
      endif

!***********************************

      return
      end subroutine RMSD_init


!=======================================================================


      subroutine RMSD_final

      use COMIFN ; use COMVAL

      implicit none

      integer(4):: i
      real(8):: rtmp,rtmp2,sumRMSD3(nREF)

!****************************************

      write(urm,*)
      where ( sWrmsd .ne. 0.d0 ) sumRMSD = sumRMSD / sWrmsd
      where ( sWrmsd .ne. 0.d0 ) sumRMSD2 = sumRMSD2 / sWrmsd
      sumRMSD3(:) = sqrt(sumRMSD2(:) - sumRMSD(:)*sumRMSD(:))

      if ( nREF.eq.1 .or. RMSD_mode.eq."average" ) then
        write(urm,*)"# Average RMSD"
        write(urm,*)"# ",real(sumRMSD(1))," (",real(sumRMSD3(1)),")"
      else
        write(urm,*)"# Average RMSD against each Ref. conf."
        rtmp = 0.d0 ; rtmp2 = 0.d0
        do i = 1,nREF
          write(urm,*)"# ",i,": ",real(sumRMSD(i))," (",               &
                     real(sumRMSD3(i)),")"
          rtmp = rtmp + sumRMSD(i)
          rtmp2 = rtmp2 + sumRMSD2(i)
        enddo
        rtmp = rtmp / dble(nREF)
        rtmp2 = rtmp2 / dble(nREF)
        rtmp2 = sqrt(rtmp2 - rtmp*rtmp)
        write(urm,*)"# All : ",real(rtmp)," (",real(rtmp2),")"
      endif
      close(urm)
      if ( interval_for_timeseries_RMSD .ne. 0 ) close(utr)

!****************************************

      return
      end subroutine RMSD_final

