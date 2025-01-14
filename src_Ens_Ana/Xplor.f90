!------------------------------------------
!
!     subroutine Xplor_noe
!     subroutine Xplor_jcc
!     subroutine Xplor_noeout
!
!------------------------------------------


      subroutine Xplor_noe

!**************************************************
!
!     Input NOE data from *.mr file for X-plor
!
!**************************************************

      use COMIFN ; use COMVAL

      implicit none

      ! Number of DATa In 1 Line
        integer(4):: NDATi1L
      ! INPut DATa separated by spaces
        character(130):: inpdat(100)
      ! Temporary RESidue number of an NOE pair
        integer(4),allocatable:: TiRES(:,:)
      ! Temporary AToM name of an NOE pair
        character(4),allocatable:: TaATM(:,:)
      ! Distance for NOE constraints
        real(4),allocatable:: dNOE(:)
      ! Difference dist. between the noe constraints and the LOW & HIGH bounds
        real(4),allocatable:: dLOW(:),dHIGH(:)

      integer(4):: i,j,k,itmp
      character(4):: atmp
      character(130):: tmp
      logical(1):: chkres,chkatm

!**************************************************
!     Input *.mr file

!     Determination of array's dimensions
      call system("wc -l "//trim(input_NOE_file_name)//" > "//         &
                  trim(d_proj2)//".aho.temp")
      open(unit=1,file=trim(d_proj2)//".aho.temp",status="old")
      read(1,*)i
      close(unit=1,status="delete")
      allocate(TiRES(2,i),TaATM(2,i),dNOE(i),dLOW(i),dHIGH(i))

!     Input experimental NOE data
      write(6,*)""
      write(6,'(2x,a)')"+ Input experimental NOE data"
      
      open(unit=1,file=trim(input_NOE_file_name),status="old")
      i = 0 ; j = 0
      do 
        read(1,'(a)',end=800)tmp
        if ( tmp(1:6).eq."assign" .or. tmp(1:6).eq."ASSIGN" ) then
          i = i + 1 ; j = j + 1 ; chkres = .false. ; chkatm = .false.
          call rdline(tmp,"!",NDATi1L,inpdat) ! ( inpinp.f90 )
          do k = 2,NDATi1L-3
            if ( inpdat(k).eq."(resid" .and. chkres.eq..false. ) then
              read(inpdat(k+1),*)TiRES(1,i) ; chkres = .true.
            elseif ( inpdat(k).eq."name" .and. chkatm.eq..false.) then
              itmp = index(inpdat(k+1)(1:4),")")
              if ( itmp .eq. 0 ) then
                TaATM(1,i) = inpdat(k+1)(1:4)
              else
                TaATM(1,i) = inpdat(k+1)(1:itmp-1)
              endif
              TaATM(1,i) = adjustl(TaATM(1,i)) ; chkatm = .true.
            elseif ( inpdat(k).eq."(resid" .and. chkres ) then
              read(inpdat(k+1),*)TiRES(2,i)
            elseif ( inpdat(k).eq."name" .and. chkatm ) then
              itmp = index(inpdat(k+1)(1:4),")")
              if ( itmp .eq. 0 ) then
                TaATM(2,i) = inpdat(k+1)(1:4)
              else
                TaATM(2,i) = inpdat(k+1)(1:itmp-1)
              endif
              TaATM(2,i) = adjustl(TaATM(2,i))
            endif
          enddo
          read(inpdat(NDATi1L-2),*)dNOE(i)
          read(inpdat(NDATi1L-1),*)dLOW(i)
          read(inpdat(NDATi1L),*)dHIGH(i)
          ! Provision for main chain Hydrogen atom name
          if ( TaATM(1,i) .eq. "HN  " ) TaATM(1,i) = "H   "
          if ( TaATM(2,i) .eq. "HN  " ) TaATM(2,i) = "H   "
          ! Remove out of residue range for the analysis
          if ( TiRES(1,i).lt.fstRES .or. TiRES(1,i).gt.fnlRES .or.     &
               TiRES(2,i).lt.fstRES .or. TiRES(2,i).gt.fnlRES ) then
               TiRES(:,i) = 0 ; TaATM(:,i) = ""
               dNOE(i) = 0.d0 ; dLOW(i) = 0.d0 ; dHIGH(i) = 0.d0
               i = i - 1
          endif
        endif
      enddo
800   close(1)
      NeNOE = i
      write(6,'(4x,a,i0,a)')"* ",j," experimental NOE data are input"
      write(6,'(4x,a,i0,a)')"* Only ",NeNOE," NOE data in "//          &
                       "the objective res. range are used"
      write(6,'(2x,a)')"+ END input experimental NOE data"
          
!     Calculate UP & LOW bound value for NOE constraints
      allocate(iRES(2,NeNOE),aATM(2,NeNOE))
      allocate(eNOE(NeNOE),UPb(NeNOE),LOWb(NeNOE))
      iRES(1:2,1:NeNOE) = TiRES(1:2,1:NeNOE)
      aATM(1:2,1:NeNOE) = TaATM(1:2,1:NeNOE)
      eNOE(1:NeNOE) = dNOE(1:NeNOE)
      UPb(1:NeNOE) = dNOE(1:NeNOE) + dHIGH(1:NeNOE)
      LOWb(1:NeNOE) = dNOE(1:NeNOE) - dLOW(1:NeNOE)

!      do i = 1,NeNOE
!        write(6,201)i,iRES(1,i),aATM(1,i),iRES(2,i),aATM(2,i),dNOE(i),&
!     &              dLOW(i),dHIGH(i)
!      enddo
!201   format(i7,x,i4,x,a4,x,i4,x,a4,x,3(f8.3,x))
      deallocate(TiRES,TaATM,dNOE,dLOW,dHIGH)

!************************
!     Translate atom name (aATM) to PRESTO format

      do i = 1,NeNOE
      do j = 1,2
        itmp = iRES(j,i) ; atmp = aATM(j,i)
        if (aRES(itmp).eq."GLY ".and.atmp.eq."HA1 ") aATM(j,i) = "HA3 "
        if (aRES(itmp).eq."PRO ".and.atmp.eq."HB1 ") aATM(j,i) = "HB3 "
        if (aRES(itmp).eq."PRO ".and.atmp.eq."HG1 ") aATM(j,i) = "HG3 "
        if (aRES(itmp).eq."PRO ".and.atmp.eq."HD1 ") aATM(j,i) = "HD3 "
        if (aRES(itmp).eq."LEU ".and.atmp.eq."HB1 ") aATM(j,i) = "HB3 "
        if (aRES(itmp).eq."ILE ".and.atmp.eq."HG11") aATM(j,i) = "HG13"
        if (aRES(itmp).eq."MET ".and.atmp.eq."HB1 ") aATM(j,i) = "HB3 "
        if (aRES(itmp).eq."MET ".and.atmp.eq."HG1 ") aATM(j,i) = "HG3 "
        if (aRES(itmp).eq."PHE ".and.atmp.eq."HB1 ") aATM(j,i) = "HB3 "
        if (aRES(itmp).eq."TYR ".and.atmp.eq."HB1 ") aATM(j,i) = "HB3 "
        if (aRES(itmp).eq."TRP ".and.atmp.eq."HB1 ") aATM(j,i) = "HB3 "
        if (aRES(itmp).eq."SER ".and.atmp.eq."HB1 ") aATM(j,i) = "HB3 "
        if (aRES(itmp)(1:3).eq."CYS".and.atmp.eq."HB1 ")               &
            aATM(j,i) = "HB3 "
        if (aRES(itmp).eq."ASN ".and.atmp.eq."HB1 ") aATM(j,i) = "HB3 "
        if (aRES(itmp).eq."GLN ".and.atmp.eq."HB1 ") aATM(j,i) = "HB3 "
        if (aRES(itmp).eq."GLN ".and.atmp.eq."HG1 ") aATM(j,i) = "HG3 "
        if (aRES(itmp)(1:3).eq."LYS".and.atmp.eq."HB1 ")               &
            aATM(j,i) = "HB3 "
        if (aRES(itmp)(1:3).eq."LYS".and.atmp.eq."HG1 ")               &
            aATM(j,i) = "HG3 "
        if (aRES(itmp)(1:3).eq."LYS".and.atmp.eq."HD1 ")               &
            aATM(j,i) = "HD3 "
        if (aRES(itmp)(1:3).eq."LYS".and.atmp.eq."HE1 ")               &
            aATM(j,i) = "HE3 "
        if (aRES(itmp)(1:3).eq."HIS".and.atmp.eq."HB1 ")               &
            aATM(j,i) = "HB3 "
        if (aRES(itmp)(1:3).eq."ARG".and.atmp.eq."HB1 ")               &
            aATM(j,i) = "HB3 "
        if (aRES(itmp)(1:3).eq."ARG".and.atmp.eq."HG1 ")               &
            aATM(j,i) = "HG3 "
        if (aRES(itmp)(1:3).eq."ARG".and.atmp.eq."HD1 ")               &
            aATM(j,i) = "HD3 "
        if (aRES(itmp)(1:3).eq."ASP".and.atmp.eq."HB1 ")               &
            aATM(j,i) = "HB3 "
        if (aRES(itmp)(1:3).eq."GLU".and.atmp.eq."HB1 ")               &
            aATM(j,i) = "HB3 "
        if (aRES(itmp)(1:3).eq."GLU".and.atmp.eq."HG1 ")               &
            aATM(j,i) = "HG3 "

      enddo
      enddo
!**************************************

      return
      end subroutine Xplor_noe


!===================================================


      subroutine Xplor_jcc

!**************************************************
!
!     Input JCC data from *.mr file for X-plor
!
!**************************************************

      use COMIFN ; use COMVAL

      implicit none

      ! check ASsigN mark
        character(6):: chkASN
      ! RESidue number of an JCC pair
        integer(4),allocatable:: iRESpre(:,:)
      ! AToM name of an JCC pair
        character(4),allocatable:: aATMpre(:,:)
      ! AVErage dihedral Angle of an jcc pair
        real(4),allocatable:: AveApre(:)

      integer(4):: i,j
      character(130):: tmp

!*****************************
!     Determination of array's dimensions
!      call system("wc -l "//trim(IJCCFN)//" > aho.temp")
!      open(unit=1,file="aho.temp",status="old")
!      read(1,*)i
!      close(unit=1,status="delete")
!      allocate(iRESpre(4,i),aATMpre(4,i),AveApre(i),DltApre(i))

!     Reading JCC data
!      write(6,*)""
!      write(6,'(2x,a)')"+ Reading JCC data from *.mr file for X-plor"

!      open(unit=1,file=trim(IJCCFN),status="old") ; i = 0 ; j = 0
!      do
!        read(1,'(a6)',end=800)chkASN
!        if ( chkASN.ne."assign" .or. chkASN.ne."ASSIGN" ) cycle
!        i = i + 1 ; j = j + 1 ; backspace(1)

        !!! KOKOKARA !!!


!        read(1,'(14x,i3,10x,a4,9x,i3,10x,a4,2x,3(f3.1,x))')            &
!     &       iRESpre(1,i),aATMpre(1,i),iRESpre(2,i),aATMpre(2,i),      &
!     &       dNOE(i),dLOW(i),dHIGH(i)
!        ! Remove out of residue range for the analysis
!        if ( iRESpre(1,i).lt.fstRES .or. iRESpre(1,i).gt.fnlRES .or.   &
!     &       iRESpre(2,i).lt.fstRES .or. iRESpre(2,i).gt.fnlRES ) then
!             iRESpre(:,i) = 0 ; aATMpre(:,i) = ""
!             dNOE(i) = 0.d0 ; dLOW(i) = 0.d0 ; dHIGH(i) = 0.d0
!             i = i - 1
!        endif
!      enddo
!800   close(1)
!      nNOE = i
!      write(tmp,*)j
!      write(6,'(4x,a)')"* "//trim(adjustl(tmp))//" NOE data are input"
!      write(tmp,*)nNOE
!      write(6,'(4x,a)')"* Only "//trim(adjustl(tmp))//" NOE data in "//&
!     &                 "the objective res. range are used"
!      write(6,'(2x,a)')"+ END input a *.mr file"

!     Calculate UP & LOW bound value for NOE constraints
!      allocate(iRES(2,nNOE),aATM(2,nNOE),UPbound(nNOE),LOWbound(nNOE))
!      iRES(1:2,1:nNOE) = iRESpre(1:2,1:nNOE)
!      aATM(1:2,1:nNOE) = aATMpre(1:2,1:nNOE)
!      UPbound(1:nNOE) = dNOE(1:nNOE) + dHIGH(1:nNOE)
!      LOWbound(1:nNOE) = dNOE(1:nNOE) - dLOW(1:nNOE)
!      deallocate(iRESpre,aATMpre,dNOE,dLOW,dHIGH)

!*****************************

      return
      end subroutine Xplor_jcc

!=================================================


      subroutine Xplor_noeout(PROJNM,NghATM,NOE,iRES,aATM,tole,rNOE)

!*************************************
!
!     OUTput reproduced NOE constraints obtained from simulated ensemble
!
!*************************************

      use COMVAL,only: d_proj2

      implicit none

      ! PROJect NaMe
        character(*),intent(in):: PROJNM
      ! Number of Geminal Hydrogen AToM
        integer(4),intent(in):: NghATM
      ! NOE distant
        real(8),intent(in):: NOE(NghATM,NghATM)
      ! RESidue number of i-th geminal hydrogen
        integer(4),intent(in):: iRES(NghATM)
      ! AToM name of i-th geminal hydrogen
        character(4),intent(in):: aATM(NghATM)
      ! Tolerance of NOE distant constraints
        real(8),intent(in):: tole(NghATM,NghATM)
      ! Real NOE intensities not to consider pseudo atom effects
        real(8),intent(in):: rNOE(NghATM,NghATM)

      integer(4):: i,j

!**********************

      write(6,'(2x,a)')"+ Reproduced NOE constraints output file : "
      write(6,'(6x,a)')trim(d_proj2)//trim(PROJNM)//"_calc.mr"

      open(unit=1,file=trim(d_proj2)//".aho.tmp",status="replace")
      do i = 1,NghATM-1
      do j = i+1,NghATM
        if ( NOE(i,j) .lt. 2.5d0 ) then
         write(1,200)iRES(i),adjustr(aATM(i)),iRES(j),adjustr(aATM(j))&
              ,2.5d0,0.7d0,tole(i,j),rNOE(i,j)
        elseif ( NOE(i,j) .lt. 3.d0 ) then
         write(1,200)iRES(i),adjustr(aATM(i)),iRES(j),adjustr(aATM(j))&
              ,3.d0,1.2d0,tole(i,j),rNOE(i,j)
        elseif ( NOE(i,j) .lt. 4.d0 ) then
         write(1,200)iRES(i),adjustr(aATM(i)),iRES(j),adjustr(aATM(j))&
              ,4.d0,2.2d0,tole(i,j),rNOE(i,j)
        elseif ( NOE(i,j) .le. 5.d0 ) then
         write(1,200)iRES(i),adjustr(aATM(i)),iRES(j),adjustr(aATM(j))&
              ,5.d0,3.2d0,tole(i,j),rNOE(i,j)
        endif
      enddo
      enddo
200   format("assign (resid",i4," and name ",a4,") (resid",i4,     &
             " and name ",a4,")",3(x,f3.1),f)
      close(1)

!*** Sort 

      call system("cat "//trim(d_proj2)//".aho.tmp | sort -rnk + 15"// &
                  "| cut -c1-70 > "//trim(d_proj2)//trim(PROJNM)//     &
                  "_calc.mr")
      call system("rm -f "//trim(d_proj2)//".aho.tmp")

!*************************************

      return
      end subroutine Xplor_noeout
