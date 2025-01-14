
      subroutine DSSP(icn)

!****************************************************
!
!     DSSP execution and analysis
!
!****************************************************

      use COMIFN ; use COMVAL

      implicit none

      integer(4),intent(in):: icn

      ! Number of column in a DSSP file
        integer(4),save:: nCLM
      ! Number of column I want to neglect
        integer(4),save:: nNEG
      ! Number of residue for the analysis
        integer(4),save:: nRES
      ! Secondary-structure expressed by one character
        character(1):: SS
      ! arrays for counting secondary structure at each residue
        real(4),allocatable,save:: H(:),B(:),E(:),G(:),pI(:),T(:),S(:),&
                                   Bla(:)
      ! Scalor for counting secondary structure of each structure
        integer(4),save:: iH,iB,iE,iG,ipI,iT,iS,iBla
      ! Secondary-Structure contents of whole length peptide chain
        real(4):: AlphaH,Helix3,Pihelix,Bbridge,Estrand,Turn,Bend,Nostr
      ! The size of Large BIN for plot
        integer(4):: Lbin
      ! phi & psi
        real(4):: phi,psi
      integer(4):: i
      real(8):: rtmp
      character(12):: ii
      character(200):: tmp,tmp2,bin
      logical(4),save:: fst_flg = .true.

!*********************************

      if ( icn .gt. 0 ) then
        ! Execute DSSP Program
        write(tmp,'(a,i0,a)') trim(d_PDB),icn,".pdb"
        write(tmp2,'(a,i0,a)') trim(d_DSSP),icn,".dssp"
        call system("mkdssp -i "//trim(tmp)//" -o "//trim(tmp2)//      &
                    " 2> /dev/null")

        ! Initial setting
        if ( fst_flg ) then
          fst_flg = .false.
          ! count column number
          call system("wc -l "//trim(tmp2)//" > "//trim(d_proj2)//    &
            ".baka.temp")
          open(unit=utmp,file=trim(d_proj2)//".baka.temp",status="old")
          read(utmp,*)nCLM
          close(unit=utmp,status="delete")
          ! find # mark
          call system("cat "//trim(tmp2)//" | gawk "//                 &
            "'{ if ( $1 == "//'"#"'//" ) print NR }' > "//             &
            trim(d_proj2)//".baka.temp")
          open(unit=utmp,file=trim(d_proj2)//".baka.temp",status="old")
          read(utmp,*)nNEG
          close(unit=utmp,status="delete")
          nRES = nCLM - nNEG
          open(unit=uss,file=trim(d_proj)//".SS",status="replace")
          open(unit=usc,file=trim(d_proj)//".sscs",status="replace")
          write(usc,'(a7,8(1x,a8))')"# str","AlphaH","3-helix",        &
               "Pi-helix","B-bridge","E-strand","Turn","Bend","No-str"
          do i = 1,nRES
            write(ii,'(i0)')i
            call system ("echo '# phi psi res. #conf. weight' > "//   &
                         trim(d_TOR)//trim(ii)//".tor")
          enddo
          ! Initialize
          allocate(H(nRES),B(nRES),E(nRES),G(nRES),pI(nRES))
          allocate(T(nRES),S(nRES),Bla(nRES))
          H(:) = 0.d0 ; B(:) = 0.d0 ; E(:) = 0.d0 ; G(:) = 0.d0
          pI(:) = 0.d0 ; T(:) = 0.d0 ; S(:) = 0.d0 ; Bla(:) = 0.d0
        endif

        ! analysis
        open(unit=utmp,file=trim(tmp2),status="old")
        do i = 1,nNEG
          read(utmp,*,end=800)
        enddo
        iH = 0 ; iB = 0 ; iE = 0 ; iG = 0 ; ipI = 0
        iT = 0 ; iS = 0 ; iBla = 0
        do i = 1,nRES
          read(utmp,'(16x,a1,86x,2f6.1)',end=800)SS,phi,psi
          if ( phi .eq. 360.0 ) phi = 0.0
          if ( psi .eq. 360.0 ) psi = 0.0
          write(ii,'(i0)')i
          write(bin,'(2(f6.1,x),2(i0,x),f)')phi,psi,i,icn,wfac
          call system ("echo "//trim(bin)//" >> "//trim(d_TOR)//       &
                       trim(ii)//".tor")
          ! count secondary structures
          if ( SS .eq. "H" ) then
            iH = iH + 1 ; H(i) = H(i) + wfac
          elseif ( SS .eq. "B" ) then
            iB = iB + 1 ; B(i) = B(i) + wfac
          elseif ( SS .eq. "E" ) then
            iE = iE + 1 ; E(i) = E(i) + wfac
          elseif ( SS .eq. "G" ) then
            iG = iG + 1 ; G(i) = G(i) + wfac
          elseif ( SS .eq. "I" ) then
            ipI = ipI + 1 ; pI(i) = pI(i) + wfac
          elseif ( SS .eq. "T" ) then
            iT = iT + 1 ; T(i) = T(i) + wfac
          elseif ( SS .eq. "S" ) then
            iS = iS + 1 ; S(i) = S(i) + wfac
          elseif ( SS .eq. " " ) then
            iBla = iBla + 1 ; Bla(i) = Bla(i) + wfac
          endif
          write(uss,'(a1,$)')SS
        enddo
800     write(usc,'(i7,8(1x,i8))')icn,iH,iG,ipI,iB,iE,iT,iS,iBla
        close(utmp)
        write(uss,*)
        if ( .not. DSSP_output_flag ) call system("rm "//trim(tmp2))

      ! Final analysis
      else
        close(usc) ; close(uss)
        ! Make secondary-structure contents at each residue
        H(:) = H(:)/sumW2  ;  B(:) = B(:)/sumW2  ;  E(:) = E(:)/sumW2
        G(:) = G(:)/sumW2  ;  pI(:) = pI(:)/sumW2  ;  T(:) = T(:)/sumW2
        S(:) = S(:)/sumW2  ;  Bla(:) = Bla(:)/sumW2
        open(unit=utmp,file=trim(d_proj)//".sscr",status="replace")
        write(utmp,'(a7,8(1x,a8))')"# res","Alpha H","3-helix",        &
             "Pi-helix","B-bridge","E-strand","Turn","Bend","No-str"
        do i = 1,nRES
          write(utmp,'(i7,8(x,e12.3))')fstRES+i-1,H(i),G(i),pI(i),B(i),&
                                       E(i),T(i),S(i),Bla(i)
        enddo
        close(utmp)

!       Make whole secondary-structure contents
        open(unit=utmp,file=trim(d_proj)//".sscw",status="replace")
        write(utmp,'(a7,8(1x,a8))')"#","Alpha H","3-helix","Pi-helix", &
             "B-bridge","E-strand","Turn","Bend","No-str"
        rtmp = 1.d0 / dble(nRES)
        AlphaH = sum(H)*rtmp ; Helix3 = sum(G)*rtmp
        Pihelix = sum(pI)*rtmp ; Bbridge = sum(B)*rtmp
        Estrand = sum(E)*rtmp ; Turn = sum(T)*rtmp
        Bend = sum(S)*rtmp ; Nostr = sum(Bla)*rtmp
        write(utmp,'(7x,8(x,e12.3))')                                  &
              AlphaH,Helix3,Pihelix,Bbridge,Estrand,Turn,Bend,Nostr
        close(utmp)

      ! determination of bin size
      if ( mod(nRES,5).eq.0 ) then
        Lbin = nRES / 5
      else
        Lbin = ( nRES / 5 ) + 1
      endif
      write(bin,*)Lbin

      open(unit=utmp,file=trim(d_proj)//".plt",status="replace")
        write(utmp,*)"# Secondary structure contents at each residue"
        write(utmp,*)""
        write(utmp,*)"set size 0.48"
        write(utmp,*)"set size ratio 0.7"
        write(tmp,*)fstRES
        write(tmp2,*)fnlRES
        write(utmp,*)"set xtics "//trim(adjustl(tmp))//","//           &
                  trim(adjustl(bin))//","//trim(adjustl(tmp2))
        write(utmp,*)"set mxtics "//trim(adjustl(bin))
        write(tmp,*)fstRES-1
        write(tmp2,*)fnlRES+1
        write(utmp,*)"set xrange ["//trim(adjustl(tmp))//":"//         &
                  trim(adjustl(tmp2)),"]"
        write(utmp,*)"set ytics 0,0.2,1"
        write(utmp,*)"set mytics 2"
        write(utmp,*)"set yrange [0:1.05]"
        write(utmp,*)'set format y "%3.1f"'
        write(utmp,*)'set xlabel "residue number" font "Helvetica"'
        write(utmp,*)'set ylabel "frequency" font "Helvetica"'
        write(utmp,*)'set key outside'
        write(utmp,*)"set style data lines"
        write(utmp,*)""
        write(utmp,'(a,$)')'plot "'//trim(d_proj)//".sscr"             &
                  //'" u 1:2 lw 3 title "Alpha H"'                     &
                  //', "" u 1:3 lw 3 title "3-helix"'                  &
                  //', "" u 1:4 lw 3 title "Pi-helix"'                 &
                  //', "" u 1:5 lw 3 title "B-bridge"'                 &
                  //', "" u 1:6 lw 3 title "E-strand"'                 &
                  //', "" u 1:7 lw 3 title "Turn"'                     &
                  //', "" u 1:8 lw 3 title "Bend"'                     &
                  //', "" u 1:9 lw 3 title "No-str"'
        write(utmp,*)""
        write(utmp,*)"set terminal postscript eps enhanced color"
        write(utmp,*)'set output "md.eps"'
        write(utmp,*)"replot"
        write(utmp,*)"set terminal x11"
        write(utmp,*)""
        write(utmp,*)"! convert -density 300 md.eps "//trim(d_proj)//  &
                  "_ER.tif"
        write(utmp,*)""
        write(utmp,*)""
        write(utmp,*)"# Secondary structure contents of whole protein"
        write(utmp,*)""
        write(utmp,*)"reset"
        write(utmp,*)"set size 0.48"
        write(utmp,*)"set style histogram rowstacked"
        write(utmp,*)"set style fill solid"
        write(utmp,*)'set format y "%3.1f"'
        write(utmp,*)'set xlabel ""'
        write(utmp,*)'set xtics 11,2,12'
        write(utmp,*)'set mxtics 0'
        write(utmp,*)'set xrange [-1:2]'
        write(utmp,*)'set ylabel "frequency"'
        write(utmp,*)'set ytics 0,0.2,1'
        write(utmp,*)'set mytics 4'
        write(utmp,*)'set yrange [0:1]'
        write(utmp,*)'set grid noxtics ytics mytics'
        write(utmp,*)''
        write(utmp,'(a,$)')'plot "'//trim(d_proj)//".sscw"       &
                  //'" u 1 title"Alpha H" with histograms'             &
                  //', "" u 2 title "3-helix" with histograms'         &
                  //', "" u 3 title "Pi-helix" with histograms'        &
                  //', "" u 4 title "B-bridge" with histograms'        &
                  //', "" u 5 title "E-strand" with histograms'        &
                  //', "" u 6 title "Turn" with histograms'            &
                  //', "" u 7 title "Bend" with histograms'            &
                  //', "" u 8 title "Nostr" with histograms'
        write(utmp,*)""
        write(utmp,*)"set terminal postscript eps enhanced color"
        write(utmp,*)'set output "md.eps"'
        write(utmp,*)"replot"
        write(utmp,*)"set terminal x11"
        write(utmp,*)""
        write(utmp,*)"! convert -density 300 md.eps "//trim(d_proj)    &
                  //"_WL.tif"
        write(utmp,*)""
        write(utmp,'(4x,a)')"!rm md.eps"
        write(utmp,*)""

        close(utmp)

!        call system("gnuplot < "//trim(d_proj)//".plt")
        write(6,*)""
        write(6,'(2x,a)')"+ The following files are output"
        write(6,*)""
        write(6,'(4x,a)')"* SS at each residue : "
        write(6,'(8x,a)')trim(project_name)//".SS"
        write(6,'(4x,a)')"* SS contents at each residue : "
        write(6,'(8x,a)')trim(project_name)//".sscr"
        write(6,'(4x,a)')"* SS contents of whole protein : "
        write(6,'(8x,a)')trim(project_name)//".sscw"
        write(6,'(4x,a)')"* SS contents of each structure : "
        write(6,'(8x,a)')trim(project_name)//".sscs"
        write(6,'(4x,a)')"* Batch file for gnuplot : "
        write(6,'(8x,a)')trim(project_name)//".plt"
        write(6,'(4x,a)')"* TIF file of SSC at each residue : "
        write(6,'(8x,a)')trim(project_name)//"_ER.tif"
        write(6,'(4x,a)')"* TIF file of SSC of whole protein : "
        write(6,'(8x,a)')trim(project_name)//"_WL.tif"
      endif

!*******************************************

      return
      end subroutine DSSP
