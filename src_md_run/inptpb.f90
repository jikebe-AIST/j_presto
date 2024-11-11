
      subroutine inptpb(iread,iprint,filena,ier)

!*******************************************************************
!
!     INPUT BINARY TOPOLOGY FILE
!
!*******************************************************************
 
      use COMBAS ; use COMERG 

      implicit none

      integer(4):: iread,iprint,ier
      character(80):: filena

      integer(4):: i,j,k
 
!************************************************

      ier = 0 ; iynpar(1:maxfnc) = 0 ; cynbpf(1:maxfnc) = " "
      iyppid(1:maxtyp,1:maxtyp) = 0
      fynbpp(1:maxnbp,1:maxtyp,1:maxtyp) = 0.d0
      write(iprint,*)" "
      write(iprint,*)"INFORMATION> INPTPB (V2.0)"
      write(iprint,*)"       READ UNFORMATTED TOPOLOGY FILE"
      write(iprint,*)" "
 
!     <<<  OPEN TOPOLOGY FILE  >>>
!         UNFORMATTED SEQUENTIAL READONLY
      call flopen(iread,filena,20,'NULL',0,ier)
      if ( ier .ne. 0 ) then
        write(iprint,*)"ERROR> INPTPB "
        write(iprint,*)"       TOPOLOGY FILE OPEN ERROR"
        write(iprint,*)" " ; ier = -1 ; return
      endif

!     <<<  READ TOPOLOGY FILE START  >>>
 
!     * READ TITLE *
      read(iread)ixtitl
      select case ( ixtitl )
        case (1:10)
          do i = 1,ixtitl
            read(iread)cxtitl(i)
          enddo

        case (11:)
          write(iprint,*)"WARNING> INPTPB "
          write(iprint,*)"       NUMBER OF TITLES IS EXCEEDED "
          write(iprint,*)"       THE LIMIT NUMBER "
          write(iprint,*)" "
          write(iprint,*)"       THE LIMIT IS 10 "
          write(iprint,*)" " ; ier = -2
      end select

!     * READ MOLECULES *
      read(iread)ixmolc ; i = ixmolc
      allocate(ixsqml(i),ixatmm(i),ixbndm(i),ixangm(i),ixtorm(i))
      allocate(iximpm(i),ixtpcn(i),cxmolc(i))
      if ( ixmolc .lt. 0 ) then
        write(iprint,*)"ERROR> INPTPB "
        write(iprint,*)"     NUMBER OF MOLECULES IS OUT OF ORDER "
        write(iprint,*)"       CURRENT IS   ",ixmolc
        write(iprint,*)" " ; ier = -3 ; return
      endif
      do i = 1,ixmolc
        read(iread)cxmolc(i),ixsqml(i)
      enddo

!     * READ ATOMS *
      read(iread)ixnatm
      if ( ixnatm.gt.maxatm .or. ixnatm.lt.0 ) then
        write(iprint,*)"ERROR> INPTPB "
        write(iprint,*)"     NUMBER OF ATOMS IS OUT OF ORDER "
        write(iprint,*)"       CURRENT IS   ",ixnatm
        write(iprint,*)"       THE LIMIT IS ",maxatm
        write(iprint,*)" " ; ier = -4 ; return
      endif
      do i = 1,ixnatm
        read(iread)ixamol(i),ixachn(i),ixares(i),ixatyp(i),            &
          ix14if(i,1:3),fxchrg(i),fxmass(i),fxvdwr(i),cxatmn(i),       &
          cxatmt(i),cxresn(i)
        j = sum(ix14if(i,1:3))
        if ( j .gt. 0 ) read(iread)ix14lt(i,1:j)
      enddo

      read(iread)ixnres
      if ( ixnres .lt. 0 ) then
        write(iprint,*)"ERROR> INPTPB "
        write(iprint,*)"     NUMBER OF RESIDUES IS OUT OF ORDER "
        write(iprint,*)"       CURRENT IS   ",ixnres
        write(iprint,*)" " ; ier = -5 ; return
      endif
      do i = 1,ixnres
        read(iread)ixrstr(i),ixrend(i)
      enddo

      read(iread)ixnchn
      if ( ixnchn .lt. 0 ) then      
        write(iprint,*)"ERROR> INPTPB "
        write(iprint,*)"     NUMBER OF CHAINS IS OUT OF ORDER "
        write(iprint,*)"       CURRENT IS   ",ixnchn
        write(iprint,*)" " ; ier = -6 ; return
      endif
      do i = 1,ixnchn
        read(iread)ixcend(i)
      enddo

!     * READ BONDS *
      read(iread)iynbnd
      if ( iynbnd .gt. maxbnd ) then
        write(iprint,*)"ERROR> INPTPB "
        write(iprint,*)"     NUMBER OF BONDS IS EXCEEDED THE LIMIT"
        write(iprint,*)"       CURRENT IS   ",iynbnd
        write(iprint,*)"       THE LIMIT IS ",maxbnd
        write(iprint,*)" " ; ier = -7 ; return
      endif
      do i = 1,iynbnd
        read(iread)iypbnd(1:2,i)
        read(iread)fyfbnd(i),fyqbnd(i)
      enddo

!     * READ ANGLES *
      read(iread)iynang
      if ( iynang .gt. maxang ) then
        write(iprint,*)"ERROR> INPTPB "
        write(iprint,*)"     NUMBER OF ANGLES IS EXCEEDED THE LIMIT "
        write(iprint,*)"       CURRENT IS   ",iynang
        write(iprint,*)"       THE LIMIT IS ",maxang
        write(iprint,*)" " ; ier = -8 ; return
      endif
      do i = 1,iynang
        read(iread)iypang(1:3,i)
        read(iread)fyfang(i),fyqang(i)
      enddo

!     * READ TORSIONS *
      read(iread)iyntor
      if ( iyntor .gt. maxtor ) then
        write(iprint,*)"ERROR> INPTPB "
        write(iprint,*)"     NUMBER OF TORSIONS IS EXCEEDED THE LIMIT "
        write(iprint,*)"       CURRENT IS   ",iyntor
        write(iprint,*)"       THE LIMIT IS ",maxtor
        write(iprint,*)" " ; ier = -9 ; return
      endif
      do i = 1,iyntor
        read(iread)iyptor(1:4,i),iytdiv(i)
        read(iread)fyftor(i),fytrot(i),fytphs(i),fytvws(i),fytess(i),  &
                   iytnbf(i)
      enddo

!     * READ IMPROPER-TORSIONS *
      read(iread)iynimp
      if ( iynimp .gt. maximp ) then
        write(iprint,*)"ERROR> INPTPB "
        write(iprint,*)"     NUMBER OF IMPROPER IS EXCEEDED THE LIMIT"
        write(iprint,*)"       CURRENT IS   ",iynimp
        write(iprint,*)"       THE LIMIT IS ",maximp
        write(iprint,*)" " ; ier = -10 ; return
      endif
      do i = 1,iynimp
        read(iread)iypimp(1:4,i),iyidiv(i)
        read(iread)fyfimp(i),fyirot(i),fyiphs(i),fyivws(i),fyiess(i),  &
                   iyinbf(i)
      enddo

!     * READ FUNCTIONS *
      read(iread)iynbpf
      if ( iynbpf.gt.maxfnc .or. iynbpf.lt.0 ) then
        write(iprint,*)"ERROR> INPTPB "
        write(iprint,*)"     MAXIMUM NUMBER OF FUNCTION NUMBER IS "
        write(iprint,*)"     OUT OF ORDER "
        write(iprint,*)"       CURRENT IS   ",iynbpf
        write(iprint,*)"       THE LIMIT IS ",maxfnc
        write(iprint,*)" " ; ier = -11 ; return
      endif
      do i = 1,iynbpf
        read(iread)iynpar(i),cynbpf(i)
      enddo

!     * READ NONBONDS *
      read(iread)iyntyp
      if ( iyntyp.gt.maxtyp .or. iyntyp.lt.0 ) then
        write(iprint,*)"ERROR> INPTPB "
        write(iprint,*)"     MAXIMUM NUMBER OF ATOM TYPE NUMBER IS "
        write(iprint,*)"     OUT OF OREDER "
        write(iprint,*)"       CURRENT IS   ",iyntyp
        write(iprint,*)"       THE LIMIT IS ",maxtyp
        write(iprint,*)" " ; ier = -12 ; return
      endif
      do i = 1,iyntyp
      do j = 1,iyntyp
        read(iread)iyppid(i,j)
        if ( iyppid(i,j) .gt. 0 ) then
          k = iynpar(iyppid(i,j))
          if ( k.gt.maxnbp .or. k.lt.0 ) then
            write(iprint,*)"ERROR> INPTPB "
            write(iprint,*)"       NUMBER OF PARAMETERS FOR NONBONDED "
            write(iprint,*)"       FUNCTION IS OUT OF ORDER "
            write(iprint,*)"       CURRENT IS   ",k
            write(iprint,*)"       THE LIMIT IS ",maxnbp
            write(iprint,*)" " ; ier = -13 ; return
          endif
          read(iread)fynbpp(1:k,i,j)
        endif
      enddo
      enddo
      do i = 1,iyntyp
        read(iread)fyvwrd(i),fyvwme(i),fy14sv(i),fy14se(i)
      enddo

!     <<<  write INFORMATION  >>>
      write(iprint,*)"       NUMBER OF MOLECULES         : ",ixmolc
      write(iprint,*)"       NUMBER OF CHAINS            : ",ixnchn
      write(iprint,*)"       NUMBER OF RESIDUES          : ",ixnres
      write(iprint,*)"       NUMBER OF ATOMS             : ",ixnatm
      write(iprint,*)"       NUMBER OF BONDS             : ",iynbnd
      write(iprint,*)"       NUMBER OF ANGLES            : ",iynang
      write(iprint,*)"       NUMBER OF TORSIONS          : ",iyntor
      write(iprint,*)"       NUMBER OF IMPRO.            : ",iynimp
      write(iprint,*)"       MAXIMUM NUMBER OF FUNC.     : ",iynbpf
      write(iprint,*)"       MAXIMUM NUMBER OF ATOM-TYPE : ",iyntyp
      write(iprint,*)" "
 
!     <<<  CLOSE TOLPOLOGY FILE  >>>
      call flclos(iread,10,ier)
      if ( ier .ne. 0 ) then
        write(iprint,*)"ERROR> INPTPB "
        write(iprint,*)"       TOPOLOGY FILE CLOSE ERROR"
        write(iprint,*)" " ; ier = -1
      endif

!*****************************************

      return
      end subroutine inptpb
