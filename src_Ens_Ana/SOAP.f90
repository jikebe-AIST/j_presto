
      subroutine SOAP

!**************************************************
!
!     SOAP calculation (hydrophobicity)
!
!**************************************************

      use COMIFN ; use COMVAL

      implicit none

      integer(4):: nres,i,j,n_amino,itmp
      integer(4),allocatable:: ir(:)
      real(8):: chg,rtmp,ave_hydro
      real(8),allocatable:: hydro(:)
      character(4),allocatable:: rname(:)
      logical(4),allocatable:: AA_flg(:)

      ! template residue name & hydrophobicity array
      character(3):: Trname(20) = (/"ILE","VAL","LEU","PHE","CYS",     &
        "MET","ALA","GLY","THR","TRP","SER","TYR","PRO","HIS","GLU",   &
        "GLN","ASP","ASN","LYS","ARG"/)
      real(4):: Thydro(20) = (/9.0,8.7,8.3,7.3,7.0,6.4,6.3,4.1,3.8,3.6,&
        3.7,3.2,2.9,1.3,1.0,1.0,1.0,1.0,0.6,0.0/)

!********************************

      open(unit=1,file=trim(d_proj)//".SOAP",status="replace")
      ! Make residue name array
      nres = fnlRES-fstRES+1
      allocate(ir(nres),hydro(nres),rname(nres),AA_flg(nres))
      ir = 0 ; hydro = 0.d0 ; AA_flg = .false.
      do i = fstRES,fnlRES
        rname(i-fstRES+1) = aRES(i)
      enddo

      ! Make residue number array
      do i = 1,nres
        do j = 1,20
          if ( rname(i)(1:3) .eq. Trname(j) ) then
            ir(i) = j ; exit
          endif
        enddo
      enddo
       
      ! Calc. net charge at pH 7
      chg = 0.d0
      do i = 1,nres
        if ( ir(i) .ne. 0 ) then
          AA_flg(i) = .true.
          if ( ir(i).eq.15 .or. ir(i).eq.17 ) then
            chg = chg - 1.d0
          elseif ( ir(i).eq.19 .or. ir(i).eq.20 ) then
            chg = chg + 1.d0
          endif
        endif
      enddo
      n_amino = count(AA_flg)
 
      ! hydrophobicity
      itmp = (nSOAP-1)/2
      hydro(1:itmp) = -1.d0 ; hydro(nres-itmp+1:nres) = -1.d0
      do i = 1+itmp,nres-itmp
        if ( all(AA_flg(i-itmp:i+itmp)) ) then
          do j = i-itmp,i+itmp
            hydro(i) = hydro(i) + Thydro(ir(j))
          enddo
        else
          hydro(i) = -1.d0
        endif
      enddo
      rtmp = 1.d0/(maxval(Thydro)*dble(nSOAP))
      hydro(1:nres) = hydro(1:nres) * rtmp

      ! Output
      itmp = fstRES - 1 ; ave_hydro = 0.d0 ; j = 0
      do i = 1,nres
        if ( hydro(i) .lt. 0.d0 ) then
          write(1,'(i10,x,a4,x,a8)')itmp+i,rname(i),"--------"
        else
          write(1,'(i10,x,a4,x,f8.3)')itmp+i,rname(i),hydro(i)
          ave_hydro = ave_hydro + hydro(i) ; j = j + 1
        endif
      enddo
      ave_hydro = ave_hydro / dble(j)

      write(1,*)
      write(1,'(a,i0,a,i0)')"# Range of residues for this analysis : ",&
                            fstRES," - ",fnlRES
      write(1,'(a,i0)')"# Number of AA residues in the range = ",n_amino
      write(1,'(a,i0)')"# Total charge of amino acids = ",nint(chg)
      ! average charge per a residue
      chg = chg / dble(n_amino)
      write(1,'(a,f8.3)')"# Ave. charge per 1 AA = ",chg
      write(1,'(a,f8.3)')"# Ave. hydrophobicity per 1 AA = ",          &
                         ave_hydro

      ! if this region is IDR or not
      if ( nSOAP .eq. 5 ) then
        rtmp = (chg+1.151d0)/2.785d0
        if ( ave_hydro .lt. rtmp ) then
          write(1,'(a)')"# This region is judged as an ID region"
        else
          write(1,'(a)')"# This region is NOT judged as an ID region"
        endif
          write(1,'(a,f8.3,a)')"# (boundary hydrophobicity = ",rtmp,")"
      endif
      close(1)

!*************************************

      return
      end subroutine SOAP
