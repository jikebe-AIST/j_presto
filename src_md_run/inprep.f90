
      subroutine inprep(iprint,cirep,ier)

      use COMBAS ; use COMERG

      implicit none

      integer(4),intent(in):: iprint
      character(80),intent(in):: cirep
      integer(4),intent(out):: ier

      integer(4):: nlist,icol,i,efcol,n,Tlst(ixnatm),Treplst(ixnatm,2)
      character(80):: space
      character(130):: tmp
      logical(4):: ex

!****************************************************

      ier = 0 ; space = " "
      inquire(file=trim(cirep),exist=ex)
      if ( .not. ex ) then
        write(iprint,*)"ERROR> REPULSION FILE OPEN ERROR IN INPGRAV"
        write(iprint,'(a80)')cirep ; ier = -1 ; return
      endif

      write(iprint,*)
      write(iprint,*)"INFORMATION> INPREP"
      write(iprint,*)"        SELECT FOLLOWING ATOMS FOR REPUL. CALC."
      open(unit=1,file=trim(cirep),status="old")
      nlist = 0 ; Krep = 0.0 ; allocate(Nrep(2))
      do
        read(1,'(a)',end=800)tmp
        icol = efcol(tmp,space,";")
        tmp = tmp(1:icol)
        i = index(tmp,"INPREP")
        if ( i .eq. 0 ) cycle
        if ( index(tmp,"KREP") .ne. 0 ) then
          do
            read(1,'(a)',end=800)tmp
            icol = efcol(tmp,space,";")
            tmp = tmp(1:icol)
            if ( tmp .eq. " " ) cycle
            read(tmp,*)Krep ; exit
          enddo
        elseif ( index(tmp,"LIST") .ne. 0 ) then
          nlist = nlist + 1
          if ( nlist .ge. 3 ) then
            write(iprint,*)"ERROR> INPREP"
            write(iprint,*)"You can input only two lists in inprep."
            ier = -1 ; return
          endif
          write(iprint,'(a,i0)')"      LIST ",nlist
          call rdlstatm(iprint,1,"INPREP",.true.,n,Tlst,ier)
          Nrep(nlist) = n ; Treplst(1:n,nlist) = Tlst(1:n)
          write(iprint,*)
        endif
      enddo
800   close(1)
      if ( Krep.eq.0.d0 .or. nlist.ne.2 ) iyeflg(15) = 0
      allocate(replst(maxval(Nrep(1:2)),2))
      if ( Nrep(1) .ge. Nrep(2) ) then
        replst(1:Nrep(1),1) = Treplst(1:Nrep(1),1)
        replst(1:Nrep(2),2) = Treplst(1:Nrep(2),2)
      else
        replst(1:Nrep(1),2) = Treplst(1:Nrep(1),1)
        replst(1:Nrep(2),1) = Treplst(1:Nrep(2),2)
        i = Nrep(1) ; Nrep(1) = Nrep(2) ; Nrep(2) = i
      endif
      write(iprint,*)
      write(iprint,*)"INFORMATION> INPREP"
      write(iprint,*)"        COEFFICIENT FOR REPULSION =  ",Krep
      Krep = Krep * 0.0001d0

!*******************************************************

      return
      end subroutine inprep
