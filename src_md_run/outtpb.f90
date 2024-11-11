
      subroutine outtpb(iwrite,iprint,filena,ier)

!*******************************************************************
!
!     OUTPUT BINARY TOPOLOGY FILE
!
!*******************************************************************
 
      use COMBAS ; use COMERG 

      implicit none

      integer(4):: iwrite,iprint,ier
      character(80):: filena

      integer(4):: i,j,i14n,numnbp

!******************************************
 
      IER = 0
      write(iprint,*)" "
      write(iprint,*)"INFORMATION> OUTTPB (V2.0)"
      write(iprint,*)"        write BINARY TOPOLOGY FILE "
      write(iprint,*)" "
 
!     <<<  OPEN TOPOLOGY FILE  >>>
      call flopen(iwrite,filena,22,'NULL',0,ier)
      if ( ier .ne. 0 ) then
        write(iprint,*)"ERROR> OUTTPB "
        write(iprint,*)"     TOPOLOGY FILE OPEN ERROR"
        write(iprint,*)" " ; ier = -1 ; return
      endif

!     <<< write TOPOLOGY FILE START  >>>
!     * write TITLE *
      write(iwrite)ixtitl
      if ( ixtitl .gt. 0 ) then
        do i = 1,ixtitl
          write(iwrite)cxtitl(i)
        enddo
      endif
 
!     * write MOLECULES *
      write(iwrite)ixmolc
      if ( ixmolc .gt. 0 ) then
        do i = 1,ixmolc
          write(iwrite)cxmolc(i),ixsqml(i)
        enddo
      endif
 
!     * write ATOMS *
      write(iwrite)ixnatm
      if ( ixnatm .gt. 0 ) then
        do i = 1,ixnatm
          write(iwrite)ixamol(i),ixachn(i),ixares(i),ixatyp(i),        &
            ix14if(i,1:3),fxchrg(i),fxmass(i),fxvdwr(i),cxatmn(i),     &
            cxatmt(i),cxresn(i)
          i14n = sum(ix14if(i,1:3))
          if ( i14n .gt. 0 ) write(iwrite)ix14lt(i,1:i14n)
        enddo
      endif
  
      write(iwrite)ixnres
      if ( ixnres .gt. 0 ) then
        do i = 1,ixnres
          write(iwrite)ixrstr(i),ixrend(i)
        enddo
      endif
      write(iwrite)ixnchn
      if ( ixnchn .gt. 0 ) then
        do i = 1,ixnchn
          write(iwrite)ixcend(i)
        enddo
      endif
 
!     * write BONDS *
      write(iwrite)iynbnd
      if ( iynbnd .gt. 0 ) then
        do i = 1,iynbnd
          write(iwrite)iypbnd(1:2,i)
          write(iwrite)fyfbnd(i),fyqbnd(i)
        enddo
      endif
 
!     * write ANGLES *
      write(iwrite)iynang
      if ( iynang .gt. 0 ) then
        do i = 1,iynang
          write(iwrite)iypang(1:3,i)
          write(iwrite)fyfang(i),fyqang(i)
        enddo
      endif
 
!     * write TORSIONS *
      write(iwrite)iyntor
      if ( iyntor .gt. 0 ) then
        do i = 1,iyntor
          write(iwrite)iyptor(1:4,i),iytdiv(i)
          write(iwrite)fyftor(i),fytrot(i),fytphs(i),fytvws(i),        &
                       fytess(i),iytnbf(i)
        enddo
      endif
 
!     * write IMPROPER-TORSIONS *
      write(iwrite)iynimp
      if ( iynimp .gt. 0 ) then
        do i = 1,iynimp
          write(iwrite)iypimp(1:4,i),iyidiv(i)
          write(iwrite)fyfimp(i),fyirot(i),fyiphs(i),fyivws(i),        &
                       fyiess(i),iyinbf(i)
        enddo
      endif
 
!     * write FUNCTIONS *
      write(iwrite)iynbpf
      if ( iynbpf .gt. 0 ) then
        do i = 1,iynbpf
          write(iwrite)iynpar(i),cynbpf(i)
        enddo
      endif
 
!     * write NONBONDS *
      write(iwrite)iyntyp
      if ( iyntyp .gt. 0 ) then
        do i = 1,iyntyp
        do j = 1,iyntyp
          write(iwrite)iyppid(i,j)
          if ( iyppid(i,j) .gt. 0 ) then
            numnbp = iynpar(iyppid(i,j))
            if ( numnbp .gt. 0 ) then
              write(iwrite)fynbpp(1:numnbp,i,j)
            endif
          endif
        enddo
        enddo

        do i = 1,iyntyp
          write(iwrite)fyvwrd(i),fyvwme(i),fy14sv(i),fy14se(i)
        enddo
      endif

      call flclos(iwrite,10,ier)
      if ( ier .ne. 0 ) then
        write(iprint,*)"ERROR> OUTTPB "
        write(iprint,*)"     TOPOLOGY FILE CLOSE ERROR"
        write(iprint,*)" " ; ier = -1 ; return
      endif

!*****************************

      return
      end subroutine outtpb
