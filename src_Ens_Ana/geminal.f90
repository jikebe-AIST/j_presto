
      subroutine geminal(iLvL,nATM,RESnm,ATMnm,ghATMnm)

!*****************************************************
!
!     Make GEMINAL hydrogen atom name list
!
!*****************************************************

      implicit none

      ! pseudo atom LeVeL
        integer(4),intent(in):: iLvL ! ( 1: methyl & aromatic )
                                     ! ( 2: methyl, aromatic, & prochiral )
      ! Number of AToM
        integer(4),intent(in):: nATM
      ! RESidue NaMe array
        character(4),intent(in):: RESnm(nATM)
      ! AToM NaMe array
        character(4),intent(in):: ATMnm(nATM)
      ! Geminal Hydrogen AToM NaMe array
        character(4),intent(inout):: ghATMnm(nATM)

      integer(4):: i
      character(4):: tATMnm(nATM)

!*********************************
!     Exchange atom names from specific to geminal one

      tATMnm(:) = adjustl(ATMnm(:))
      ghATMnm(:) = tATMnm(:)

      if ( iLvL .eq. 1 ) then

        do i = 1,nATM
          ! VAL
          if ( RESnm(i) .eq. "VAL " ) then
            if ( tATMnm(i)(1:3) .eq. "HG1" ) ghATMnm(i) = "HG1#"
            if ( tATMnm(i)(1:3) .eq. "HG2" ) ghATMnm(i) = "HG2#"
          ! ILE
          elseif ( RESnm(i) .eq. "ILE " ) then
            if ( tATMnm(i)(1:3) .eq. "HG2" ) ghATMnm(i) = "HG2#"
            if ( tATMnm(i)(1:2) .eq. "HD" ) ghATMnm(i) = "HD# "
          ! LEU
          elseif ( RESnm(i) .eq. "LEU " ) then
            if ( tATMnm(i)(1:3) .eq. "HD1" ) ghATMnm(i) = "HD1#"
            if ( tATMnm(i)(1:3) .eq. "HD2" ) ghATMnm(i) = "HD2#"
          ! THR
          elseif ( RESnm(i) .eq. "THR " ) then
            if ( tATMnm(i)(1:3) .eq. "HG2" ) ghATMnm(i) = "HG2#"
          ! LYS
          elseif ( RESnm(i)(1:3) .eq. "LYS" ) then
            if ( tATMnm(i)(1:2) .eq. "HZ" ) ghATMnm(i) = "HZ# "
          ! MET
          elseif ( RESnm(i) .eq. "MET " ) then
            if ( tATMnm(i)(1:2) .eq. "HE" ) ghATMnm(i) = "HE# "
          ! PHE
          elseif ( RESnm(i) .eq. "PHE " ) then
            if ( tATMnm(i)(1:2) .eq. "HD" ) ghATMnm(i) = "HD# "
            if ( tATMnm(i)(1:2) .eq. "HE" ) ghATMnm(i) = "HE# "
          ! TYR
          elseif ( RESnm(i) .eq. "TYR " ) then
            if ( tATMnm(i)(1:2) .eq. "HD" ) ghATMnm(i) = "HD# "
            if ( tATMnm(i)(1:2) .eq. "HE" ) ghATMnm(i) = "HE# "
          else
!            write(6,'(a,i0,a)')"RESIDUE name of ",i,"-th atom : "//    &
!     &           RESnm(i)
!            write(6,*)"place : geminal.f90"
!            write(6,*)"reason : The above residue name is unknown"
!            write(6,*)"Program STOP!!"
!            stop
          endif
        enddo

      elseif ( iLvL .eq. 2 ) then

        do i = 1,nATM
          ! GLY
          if ( RESnm(i) .eq. "GLY " ) then
            if ( tATMnm(i)(1:2) .eq. "HA" ) ghATMnm(i) = "HA# "
          ! ALA
          elseif ( RESnm(i) .eq. "ALA " ) then
            if ( tATMnm(i)(1:2) .eq. "HB" ) ghATMnm(i) = "HB# "
          ! VAL
          elseif ( RESnm(i) .eq. "VAL " ) then
            if ( tATMnm(i)(1:2) .eq. "HG" ) ghATMnm(i) = "HG# "
          ! ILE
          elseif ( RESnm(i) .eq. "ILE " ) then
            if ( tATMnm(i)(1:3) .eq. "HG1" ) ghATMnm(i) = "HG1#"
            if ( tATMnm(i)(1:2) .eq. "HD" ) ghATMnm(i) = "HD# "
            if ( tATMnm(i)(1:3) .eq. "HG2" ) ghATMnm(i) = "HG2#"
          ! LEU
          elseif ( RESnm(i) .eq. "LEU " ) then
            if ( tATMnm(i)(1:2) .eq. "HB" ) ghATMnm(i) = "HB# "
            if ( tATMnm(i)(1:2) .eq. "HD" ) ghATMnm(i) = "HD# "
          ! PRO
          elseif ( RESnm(i) .eq. "PRO " ) then
            if ( tATMnm(i)(1:2) .eq. "HB" ) ghATMnm(i) = "HB# "
            if ( tATMnm(i)(1:2) .eq. "HG" ) ghATMnm(i) = "HG# "
            if ( tATMnm(i)(1:2) .eq. "HD" ) ghATMnm(i) = "HD# "
          ! SER
          elseif ( RESnm(i) .eq. "SER " ) then
            if ( tATMnm(i)(1:2) .eq. "HB" ) ghATMnm(i) = "HB# "
          ! THR
          elseif ( RESnm(i) .eq. "THR " ) then
            if ( tATMnm(i)(1:3) .eq. "HG2" ) ghATMnm(i) = "HG2#"
          ! ASP
          elseif ( RESnm(i)(1:3) .eq. "ASP" ) then
            if ( tATMnm(i)(1:2) .eq. "HB" ) ghATMnm(i) = "HB# "
          ! ASN
          elseif ( RESnm(i) .eq. "ASN " ) then
            if ( tATMnm(i)(1:2) .eq. "HB" ) ghATMnm(i) = "HB# "
            if ( tATMnm(i)(1:3) .eq. "HD2" ) ghATMnm(i) = "HD2#"
          ! GLU
          elseif ( RESnm(i)(1:3) .eq. "GLU" ) then
            if ( tATMnm(i)(1:2) .eq. "HB" ) ghATMnm(i) = "HB# "
            if ( tATMnm(i)(1:2) .eq. "HG" ) ghATMnm(i) = "HG# "
          ! GLN
          elseif ( RESnm(i) .eq. "GLN " ) then
            if ( tATMnm(i)(1:2) .eq. "HB" ) ghATMnm(i) = "HB# "
            if ( tATMnm(i)(1:2) .eq. "HG" ) ghATMnm(i) = "HG# "
            if ( tATMnm(i)(1:3) .eq. "HE2" ) ghATMnm(i) = "HE2#"
          ! LYS
          elseif ( RESnm(i)(1:3) .eq. "LYS" ) then
            if ( tATMnm(i)(1:2) .eq. "HB" ) ghATMnm(i) = "HB# "
            if ( tATMnm(i)(1:2) .eq. "HG" ) ghATMnm(i) = "HG# "
            if ( tATMnm(i)(1:2) .eq. "HD" ) ghATMnm(i) = "HD# "
            if ( tATMnm(i)(1:2) .eq. "HE" ) ghATMnm(i) = "HE# "
            if ( tATMnm(i)(1:2) .eq. "HZ" ) ghATMnm(i) = "HZ# "
          ! ARG
          elseif ( RESnm(i)(1:3) .eq. "ARG" ) then
            if ( tATMnm(i)(1:2) .eq. "HB" ) ghATMnm(i) = "HB# "
            if ( tATMnm(i)(1:2) .eq. "HG" ) ghATMnm(i) = "HG# "
            if ( tATMnm(i)(1:2) .eq. "HD" ) ghATMnm(i) = "HD# "
            if ( tATMnm(i)(1:2) .eq. "HH" ) ghATMnm(i) = "HH# "
          ! MET
          elseif ( RESnm(i) .eq. "MET " ) then
            if ( tATMnm(i)(1:2) .eq. "HB" ) ghATMnm(i) = "HB# "
            if ( tATMnm(i)(1:2) .eq. "HG" ) ghATMnm(i) = "HG# "
            if ( tATMnm(i)(1:2) .eq. "HE" ) ghATMnm(i) = "HE# "
          ! CYS
          elseif ( RESnm(i)(1:3) .eq. "CYS" ) then
            if ( tATMnm(i)(1:2) .eq. "HB" ) ghATMnm(i) = "HB# "
          ! PHE
          elseif ( RESnm(i) .eq. "PHE " ) then
            if ( tATMnm(i)(1:2) .eq. "HB" ) ghATMnm(i) = "HB# "
            if ( tATMnm(i)(1:2) .eq. "HD" ) ghATMnm(i) = "HD# "
            if ( tATMnm(i)(1:2) .eq. "HE" ) ghATMnm(i) = "HE# " 
          ! TYR
          elseif ( RESnm(i) .eq. "TYR " ) then
            if ( tATMnm(i)(1:2) .eq. "HB" ) ghATMnm(i) = "HB# "
            if ( tATMnm(i)(1:2) .eq. "HD" ) ghATMnm(i) = "HD# "
            if ( tATMnm(i)(1:2) .eq. "HE" ) ghATMnm(i) = "HE# "
          ! HIS
          elseif ( RESnm(i)(1:3) .eq. "HIS" ) then
            if ( tATMnm(i)(1:2) .eq. "HB" ) ghATMnm(i) = "HB# "
          ! TRP
          elseif ( RESnm(i) .eq. "TRP " ) then
            if ( tATMnm(i)(1:2) .eq. "HB" ) ghATMnm(i) = "HB# "
          else
!            write(6,'(a,i0,a)')"RESIDUE name of ",i,"-th atom : "//    &
!     &           RESnm(i)
!            write(6,*)"place : geminal.f90"
!            write(6,*)"reason : The above residue name is unknown"
!            write(6,*)"Program STOP!!"
!            stop
          endif

        enddo

      endif

!***************************************

      return
      end subroutine geminal
