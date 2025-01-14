!--------------------------------------------------
!
!     subroutine DISCOVER_noe
!     subroutine DISCOVER_jcc
!
!--------------------------------------------------


      subroutine DISCOVER_noe

!**************************************************
!
!     Input *.mr file for DISCOVER
!
!**************************************************

      use COMIFN ; use COMVAL

      implicit none

      ! Temporary RESidue number of an NOE pair
        integer(4),allocatable:: TiRES(:,:)
      ! Temporary AToM & Geminal Hydrogen name of an NOE pair
        character(4),allocatable:: TaATM(:,:)
      ! Temporary UPper & LOWer Bound of an noe pair
        real(4),allocatable:: tUPb(:),tLOWb(:)

      integer(4):: i,j,k,itmp
      real(8):: xxx
      character(4):: tt,atmp
      character(13):: chk
      character(130):: tmp

!****************************************
!     Input *.mr file

!     Case conversion of *.mr file
      call system("cat "//trim(input_NOE_file_name)//                  &
                 " | tr '[a-z]' '[A-Z]' | "//                         &
                 " tr ':' ' ' | tr '_' ' ' > "//trim(d_proj2)//"tmp.mr")

!     Determination of array's dimensions
      call system("wc -l "//trim(d_proj2)//"tmp.mr > "//trim(d_proj2)//&
                  "aho.temp")
      open(unit=1,file=trim(d_proj2)//"aho.temp",status="old")
      read(1,*)i
      close(unit=1,status="delete")
      allocate(TiRES(2,i),TaATM(2,i),tUPb(i),tLOWb(i))
     
!     Input experimental NMR data
      write(6,*)""
      write(6,'(2x,a)')"+ Input experimental NOE data"

      open(unit=1,file=trim(d_proj2)//"tmp.mr",status="old")
      i = 0 ; j = 0
      do
        read(1,'(a13)',end=800)chk
        ! NOE data
        if ( chk .eq. "#NOE DISTANCE" ) then

          do
            read(1,'(a13)',end=800)chk
            if ( chk(1:1).eq."!" .or. chk(1:1).eq." " ) exit
            i = i + 1 ; j = j + 1 ; backspace(1)

            read(1,*)k,tt,TiRES(1,i),TaATM(1,i),k,tt,                  &
                     TiRES(2,i),TaATM(2,i),xxx,tUPb(i),tLOWb(i)

            ! Provision for main chain hydrogen atom name
            if ( TaATM(1,i) .eq. "HN  " ) TaATM(1,i) = "H   "
            if ( TaATM(2,i) .eq. "HN  " ) TaATM(2,i) = "H   "
            ! Remove out of residue range for the analysis
            if ( TiRES(1,i).lt.fstRES .or. TiRES(1,i).gt.fnlRES        &
            .or. TiRES(2,i).lt.fstRES .or. TiRES(2,i).gt.fnlRES ) then
              TiRES(:,i) = 0 ; TaATM(:,i) = "" 
              tUPb(i) = 0.d0 ; tLOWb(i) = 0.d0
              i = i - 1
            endif
            ! Remove strange residue name character (N-terminal)
            if ( TaATM(1,i)(4:4) .eq. "N" ) TaATM(1,i)(4:4) = " "
            if ( TaATM(2,i)(4:4) .eq. "N" ) TaATM(2,i)(4:4) = " "
          enddo

        endif
      enddo
800   close(unit=1,status="delete")
      NeNOE = i
      
!***************************************************
!     Translate atom name (TaATM) to PRESTO format
      
      do i = 1,NeNOE
      do j = 1,2
        itmp = TiRES(j,i) ; atmp = TaATM(j,i)
        if (aRES(itmp).eq."GLY ".and.atmp.eq."HA1 ") TaATM(j,i) = "HA3 "
        if (aRES(itmp).eq."PRO ".and.atmp.eq."HB1 ") TaATM(j,i) = "HB3 "
        if (aRES(itmp).eq."PRO ".and.atmp.eq."HG1 ") TaATM(j,i) = "HG3 "
        if (aRES(itmp).eq."PRO ".and.atmp.eq."HD1 ") TaATM(j,i) = "HD3 "
        if (aRES(itmp).eq."LEU ".and.atmp.eq."HB1 ") TaATM(j,i) = "HB3 "
        if (aRES(itmp).eq."ILE ".and.atmp.eq."HG11") TaATM(j,i) = "HG13"
        if (aRES(itmp).eq."MET ".and.atmp.eq."HB1 ") TaATM(j,i) = "HB3 "
        if (aRES(itmp).eq."MET ".and.atmp.eq."HG1 ") TaATM(j,i) = "HG3 "
        if (aRES(itmp).eq."PHE ".and.atmp.eq."HB1 ") TaATM(j,i) = "HB3 "
        if (aRES(itmp).eq."TYR ".and.atmp.eq."HB1 ") TaATM(j,i) = "HB3 "
        if (aRES(itmp).eq."TRP ".and.atmp.eq."HB1 ") TaATM(j,i) = "HB3 "
        if (aRES(itmp).eq."SER ".and.atmp.eq."HB1 ") TaATM(j,i) = "HB3 "
        if (aRES(itmp)(1:3).eq."CYS".and.atmp.eq."HB1 ")               &
            TaATM(j,i) = "HB3 "
        if (aRES(itmp).eq."ASN ".and.atmp.eq."HB1 ") TaATM(j,i) = "HB3 "
        if (aRES(itmp).eq."GLN ".and.atmp.eq."HB1 ") TaATM(j,i) = "HB3 "
        if (aRES(itmp).eq."GLN ".and.atmp.eq."HG1 ") TaATM(j,i) = "HG3 "
        if (aRES(itmp)(1:3).eq."LYS".and.atmp.eq."HB1 ")               &
            TaATM(j,i) = "HB3 "
        if (aRES(itmp)(1:3).eq."LYS".and.atmp.eq."HG1 ")               &
            TaATM(j,i) = "HG3 "
        if (aRES(itmp)(1:3).eq."LYS".and.atmp.eq."HD1 ")               &
            TaATM(j,i) = "HD3 "
        if (aRES(itmp)(1:3).eq."LYS".and.atmp.eq."HE1 ")               &
            TaATM(j,i) = "HE3 "
        if (aRES(itmp)(1:3).eq."HIS".and.atmp.eq."HB1 ")               &
            TaATM(j,i) = "HB3 "
        if (aRES(itmp)(1:3).eq."ARG".and.atmp.eq."HB1 ")               &
            TaATM(j,i) = "HB3 "
        if (aRES(itmp)(1:3).eq."ARG".and.atmp.eq."HG1 ")               &
            TaATM(j,i) = "HG3 "
        if (aRES(itmp)(1:3).eq."ARG".and.atmp.eq."HD1 ")               &
            TaATM(j,i) = "HD3 "
        if (aRES(itmp)(1:3).eq."ASP".and.atmp.eq."HB1 ")               &
            TaATM(j,i) = "HB3 "
        if (aRES(itmp)(1:3).eq."GLU".and.atmp.eq."HB1 ")               &
            TaATM(j,i) = "HB3 "
        if (aRES(itmp)(1:3).eq."GLU".and.atmp.eq."HG1 ")               &
            TaATM(j,i) = "HG3 "

      enddo
      enddo
      
      allocate(iRES(2,NeNOE),aATM(2,NeNOE),UPb(NeNOE),LOWb(NeNOE))
      allocate(eNOE(NeNOE))
      iRES(:,:) = TiRES(:,1:NeNOE) ; aATM(:,:) = TaATM(:,1:NeNOE)
      UPb(:) = tUPb(1:NeNOE) ; LOWb(:) = tLOWb(1:NeNOE) 
      eNOE(:) = tUPb(1:NeNOE)
      do i = 1,NeNOE
        if ( eNOE(i) .gt. 5.d0 ) eNOE(i) = 5.d0
      enddo

      write(6,'(4x,a,i0,a)')"* ",NeNOE," NOE data in "//          &
                       "the objective res. range are used"
      write(6,'(2x,a)')"+ END input experimental NOE data"

!***************************************************

      return
      end subroutine DISCOVER_noe


!==========================================================


      subroutine DISCOVER_jcc

!**************************************************
!
!     Input JCC data from *.mr file for DISCOVER
!
!**************************************************

      use COMIFN ; use COMVAL

      implicit none

      ! RESidue number of an JCC pair
        integer(4),allocatable:: iRESpre(:,:)
      ! AToM name of an JCC pair
        character(4),allocatable:: aATMpre(:,:)
      ! AVErage dihedral Angle of an jcc pair
        real(4),allocatable:: AveApre(:)
      ! DeLTa dihedral Angle from the average jcc one
        real(4),allocatable:: DltApre(:)

      integer(4):: i,j,k
      real(8):: xx1,xx2
      character(4):: tt
      character(13):: chk
      character(130):: tmp

!******************************
!     Case conversion of *.mr file
      call system("cat "//trim(input_JCC_file_name)//                  &
                 " | tr '[a-z]' '[A-Z]' | "//                         &
                 " tr ':' ' ' | tr '_' ' ' > "//trim(d_proj2)//"tmp.mr")

!     Determination of array's dimensions
      call system("wc -l "//trim(d_proj2)//"tmp.mr > "//trim(d_proj2)//&
                  "aho.temp")
      open(unit=1,file=trim(d_proj2)//"aho.temp",status="old")
      read(1,*)i
      close(unit=1,status="delete")
      allocate(iRESpre(4,i),aATMpre(4,i),AveApre(i),DltApre(i))

!     Reading JCC data
      write(6,*)""
      write(6,'(2x,a)')"+ Reading JCC data from *.mr file for DISCOVER"

      open(unit=1,file=trim(d_proj2)//"tmp.mr",status="old")
      i = 0 ; j = 0
      do
        read(1,'(a13)',end=800)chk
        ! JCC data
        if ( chk .eq. "#NMR_dihedral" ) then

          do
            read(1,'(a13)',end=800)chk
            if ( chk(1:1).eq." " .or. chk(1:1).eq."!" ) exit
            i = i + 1 ; j = j + 1 ; backspace(1)

            read(1,*)k,tt,iRESpre(1,i),aATMpre(1,i),k,tt,              &
                     iRESpre(2,i),aATMpre(2,i),k,tt,                   &
                     iRESpre(3,i),aATMpre(3,i),k,tt,                   &
                     iRESpre(4,i),aATMpre(4,i),xx1,xx2
            AveApre(i) = (xx1+xx2) / 2.d0
            DltApre(i) = (xx2-xx1) / 2.d0
            ! Remove out of residue range for the analysis
            if ( iRESpre(1,i).lt.fstRES .or. iRESpre(1,i).gt.fnlRES    &
            .or. iRESpre(2,i).lt.fstRES .or. iRESpre(2,i).gt.fnlRES    &
            .or. iRESpre(3,i).lt.fstRES .or. iRESpre(3,i).gt.fnlRES    &
            .or. iRESpre(4,i).lt.fstRES .or. iRESpre(4,i).gt.fnlRES ) then
              iRESpre(:,i) = 0 ; aATMpre(:,i) = ""
              AveApre(i) = 0.d0 ; DltApre(i) = 0.d0
              i = i - 1
            endif
            ! Remove strange residue name character
            if ( aATMpre(1,i)(4:4) .eq. "N" ) aATMpre(1,i)(4:4) = " "
            if ( aATMpre(2,i)(4:4) .eq. "N" ) aATMpre(2,i)(4:4) = " "
            
            
          enddo

        endif
      enddo
800   close(unit=1,status="delete")

      nJCC = i
      write(tmp,*)j
      write(6,'(4x,a)')"* "//trim(adjustl(tmp))//" JCC data are input"
      write(tmp,*)nJCC
      write(6,'(4x,a)')"* Only "//trim(adjustl(tmp))//" JCC data in "//&
                       "the objective res. range are used"
      write(6,'(2x,a)')"+ END input a *.mr file"

!     Determinate various matrixes and vectors
      allocate(iRESJ(4,nJCC),aATMJ(4,nJCC),AveA(nJCC),DltA(nJCC))
      iRESJ(1:4,1:nJCC) = iRESpre(1:4,1:nJCC)
      aATMJ(1:4,1:nJCC) = aATMpre(1:4,1:nJCC)
      AveA(1:nJCC) = AveApre(1:nJCC)
      DltA(1:nJCC) = DltApre(1:nJCC)
      deallocate(iRESpre,aATMpre,AveApre,DltApre)

!***************************************************

      return
      end subroutine DISCOVER_jcc
