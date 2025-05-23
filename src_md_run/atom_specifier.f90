
      subroutine atom_specifier(n,specifiers,natm,anum,aname,rnum,     &
                                rname,cnum,cname,cod,icount,flist)

      implicit none

      integer(4),intent(in):: n
      character(n),intent(in):: specifiers
      integer(4),intent(in):: natm
      integer(4),intent(in):: anum(natm)
      character(*),intent(in):: aname(natm)
      integer(4),intent(in):: rnum(natm)
      character(*),intent(in):: rname(natm)
      integer(4),intent(in):: cnum(natm)
      character(*),intent(in):: cname(natm)
      real(4),intent(in):: cod(3,natm)
      integer(4),intent(out):: icount
      integer(4),intent(out):: flist(natm)

      integer(4):: ii,jj,kk,ll
      integer(4),allocatable:: list(:),sum_list(:)
      character(n),allocatable:: specifier(:)
      logical(1):: fmark(natm)

!***********************************************************************

      ! Separate '+' sign
      kk = 0 ; jj = 1
      do
        ii = index(specifiers(kk+1:),"+")
        if ( ii .ne. 0 ) then
          jj = jj + 1 ; kk = ii
        else
          exit
        endif
      enddo
      allocate(specifier(jj),sum_list(natm*jj))
      kk = 0 ; jj = 0
      do
        ii = index(specifiers(kk+1:),"+")
        if ( ii .ne. 0 ) then
          jj = jj + 1
          specifier(jj) = specifiers(kk+1:kk+ii-1)
          kk = ii
        else
          jj = jj + 1
          specifier(jj) = specifiers(kk+ii+1:)
          exit
        endif
      enddo

      kk = 1
      do ii = 1,jj
        call separate_specifier(len(trim(specifier(ii))),             &
          trim(specifier(ii)),fmark)
        ll = count(fmark) ; allocate(list(ll))
        call make_specifier_list(ll,natm,list,fmark)
        sum_list(kk:kk+ll-1) = list(1:ll) ; kk = kk + ll
        deallocate(list)
      enddo
      icount = kk - 1
      if ( icount .eq. 0 ) then
        write(6,*)"!! ERROR !!"
        write(6,*)"None of atoms were selected under the specified"//" &
                  selection criteria "//trim(specifiers)//"."
        stop
      endif
      flist(1:icount) = sum_list(1:icount)

      return
      contains

!***********************************************************************

      recursive subroutine separate_specifier(nlen,sp,mark)

      integer(4),intent(in):: nlen
      character(nlen),intent(in):: sp
      logical(1),intent(inout):: mark(natm)

      integer(4):: i,j,ic,ispa
      character(nlen),allocatable:: spa(:,:)
      logical(1):: check,nmark(natm),logic

!***********************************************************************

      ic = 1 ; check = .false.
      do i = 1,nlen
        if ( sp(i:i) .eq. "(" ) then
          check = .true.
        elseif ( sp(i:i) .eq. ")" ) then
          check = .false.
        endif
        if ( check ) cycle
        if ( sp(i:i) .eq. "&" ) then
          ic = ic + 1
        elseif ( sp(i:i) .eq. "|" ) then
          ic = ic + 1
        endif
      enddo
      ispa = ic

      allocate(spa(2,ic)) ; spa(1,1) = ""
      ic = 1 ; j = 1 ; check = .false.
      do i = 1,nlen
        if ( sp(i:i) .eq. "(" ) then
          check = .true.
        elseif ( sp(i:i) .eq. ")" ) then
          check = .false.
        endif
        if ( check ) cycle
        if ( sp(i:i) .eq. "&" ) then
          spa(2,ic) = sp(j:i-1) ; ic = ic + 1
          spa(1,ic) = "&" ; j = i + 1
        elseif ( sp(i:i) .eq. "|" ) then
          spa(2,ic) = sp(j:i-1) ; ic = ic + 1
          spa(1,ic) = "|" ; j = i + 1
        endif
      enddo
      spa(2,ic) = sp(j:nlen)

      do i = 1,ispa
        nmark(:) = .false.
        logic = ( spa(2,i)(1:1) .ne. "!" )
        if ( .not. logic ) spa(2,i) = spa(2,i)(2:)
        if ( trim(spa(2,i)(1:1)).ne."(" )  then
          call check_specifier(spa(1,i)(1:1),len(trim(spa(2,i))),      &
                               trim(spa(2,i)),nmark)
        else
          j = len(trim(spa(2,i)))
          call separate_specifier(j-2,spa(2,i)(2:j-1),nmark)
        endif

        if ( .not. logic ) nmark = .not. nmark
        if ( len(trim(spa(1,i))) .eq. 0 ) then
          mark(:) = nmark(:)
        elseif ( trim(spa(1,i)) .eq. "&" ) then
          mark = mark .and. nmark
        elseif ( trim(spa(1,i)) .eq. "|" ) then
          mark = mark .or. nmark
        endif
      enddo

      return
      end subroutine separate_specifier

!***********************************************************************


      subroutine check_specifier(kigo,nlen,sp,mark)

      character(1),intent(in):: kigo
      integer(4),intent(in):: nlen
      character(nlen),intent(in):: sp
      logical(1),intent(inout):: mark(natm)

      integer(4):: i
      character(nlen):: cond
      character(nlen),allocatable:: spa(:)

!***********************************************************************

      cond = sp
      if ( cond(1:1) .eq. "[" ) then
        cond = cond(2:len(trim(cond))-1)
        if ( index(cond,":") .ne. 0 ) then
          call range_selector(len(trim(cond)),cond,mark)
!        else
!          call cod_selector
        endif
      else
        call atom_selector(len(trim(cond)),cond,mark)
      endif

      return
      end subroutine check_specifier


!***********************************************************************


      subroutine atom_selector(nlen,cond,mark)

      integer(4),intent(in):: nlen
      character(nlen),intent(in):: cond
      logical(1),intent(inout):: mark(natm)

      integer(4):: i,j
      character(nlen):: chains,residues,atoms
      logical(1):: rmark(natm),amark(natm),cmark(natm)

!***********************************************************************

      chains = "*" ; residues = "*" ; atoms = "*"
      i = index(cond,":") ; j = index(cond(i+1:),":")
      if ( j .ne. 0 ) then
        chains = cond(:i-1) ; residues = cond(i+1:i+j-1)
        atoms = cond(i+j+1:)
      elseif ( i .ne. 0 ) then
        residues = cond(:i-1) ; atoms = cond(i+1:)
      else
        atoms = cond
      endif

      if ( chains .eq. "*" ) then
        cmark(:) = .true.
      else
        cmark(:) = .false.
        call cond_checker("c",len(trim(chains)),chains,cmark)
      endif
      if ( residues .eq. "*" ) then
        rmark(:) = .true.
      else
        rmark(:) = .false.
        call cond_checker("r",len(trim(residues)),residues,rmark)
      endif
      if ( atoms .eq. "*" ) then
        amark(:) = .true.
      else
        amark(:) = .false.
        call cond_checker("a",len(trim(atoms)),atoms,amark)
      endif

      where(cmark .and. rmark .and. amark) mark = .true.

      return
      end subroutine atom_selector


!***********************************************************************


      subroutine cond_checker(knd,nlen,txt,mark)

      character(1),intent(in):: knd
      integer(4),intent(in):: nlen
      character(nlen),intent(in):: txt
      logical(1),intent(inout):: mark(natm)
      integer(4):: i,j
      character(nlen),allocatable:: array(:)
      logical(1):: check

      call name_separator(nlen,txt,",",array)
      j = nlen
      do i = 1,size(array)
        call digit_check(nlen,array(i),check)
        if ( check ) then
          if ( knd.eq."c" ) call number_checker(j,array(i),cnum,mark)
          if ( knd.eq."r" ) call number_checker(j,array(i),rnum,mark)
          if ( knd.eq."a" ) call number_checker(j,array(i),anum,mark)
        else
          if ( knd.eq."c") call name_checker(j,array(i),cname,mark)
          if ( knd.eq."r") call name_checker(j,array(i),rname,mark)
          if ( knd.eq."a") call name_checker(j,array(i),aname,mark)
        endif
      enddo

      return
      end subroutine cond_checker


!***********************************************************************


      subroutine name_separator(nlen,txt,sep,array)

      integer(4),intent(in):: nlen
      character(nlen),intent(in):: txt
      character(1),intent(in):: sep
      character(nlen),intent(inout),allocatable:: array(:)
      integer(4):: i,j,ic

      j = 0 ; ic = 0
      do
        i = index(txt(j+1:),sep)
        ic = ic + 1
        if ( i .eq. 0 ) exit
        j = j + i
      enddo
      if ( ic .eq. 0 ) return
      allocate(array(ic))

      j = 0 ; ic = 0
      do
        i = index(txt(j+1:),sep)
        ic = ic + 1
        if ( i .eq. 0 ) then
          array(ic) = txt(j+1:) ; exit
        else
          array(ic) = txt(j+1:j+i-1)
          j = j + i
        endif
      enddo

      return
      end subroutine name_separator


!***********************************************************************


      subroutine digit_check(nlen,txt,digit)

      integer(4),intent(in):: nlen
      character(nlen),intent(in):: txt
      logical(1),intent(out):: digit
      integer(4):: i
      character(1):: c

      digit = .true.
      do i = 1,len(trim(txt))
        c = txt(i:i)
        if ( c.ne."<" .and. .not. (c.ge."0" .and. c.le."9") ) then
          digit = .false. ; exit
        endif
      enddo

      return
      end subroutine digit_check


!***********************************************************************


      subroutine number_checker(nlen,txt,num,mark)

      integer(4),intent(in):: nlen
      character(nlen),intent(in):: txt
      integer(4),intent(in):: num(natm)
      logical(1),intent(inout):: mark(natm)
      character(nlen),allocatable:: cond2(:)
      integer(4):: i,j

      call name_separator(nlen,txt,"<",cond2)
      if ( size(cond2) .eq. 1 ) then
        read(cond2(1),*)i ; where(num .eq. i) mark = .true.
      else
        read(cond2(1),*)i ; read(cond2(2),*)j
        where(num.ge.i .and. num.le.j) mark = .true.
      endif

      return
      end subroutine number_checker


!***********************************************************************


      subroutine name_checker(nlen,txt,nam,mark)

      integer(4),intent(in):: nlen
      character(nlen),intent(in):: txt
      character(*),intent(in):: nam(natm)
      logical(1),intent(inout):: mark(natm)
      integer(4):: i,j,k
      character(999):: tmp

      i = index(txt,"*")
      if ( i .eq. 0 ) then
        do j = 1,natm
          tmp = trim(adjustl(nam(j)))
          if ( trim(txt) .eq. trim(adjustl(nam(j))) ) mark(j) = .true.
        enddo
      elseif ( i .eq. 1 ) then
        mark(:) = .true.
      else
        do j = 1,natm
          tmp = trim(adjustl(nam(j)))
          if ( txt(:i-1) .eq. tmp(:i-1) ) mark(j) = .true.
        enddo
      endif

      return
      end subroutine name_checker


!***********************************************************************


      subroutine range_selector(nlen,cond,mark)

      integer(4),intent(in):: nlen
      character(nlen),intent(in):: cond
      logical(1),intent(inout):: mark(natm)

      integer(4):: i,j
      character(nlen):: x,y,z
      logical(1):: xmark(natm),ymark(natm),zmark(natm)

!***********************************************************************

      i = index(cond,":") ; j = index(cond(i+1:),":")
      if ( i.eq.0 .or. j.eq.0 ) then
        write(6,*)"!! ERROR !!"
        write(6,*)'The range selector must be given as "[x:y:z]".'
        stop
      endif
      x = cond(:i-1) ; y = cond(i+1:i+j-1) ; z = cond(i+j+1:)

      if ( len(trim(x)) .eq. 0 ) then
        xmark(:) = .true.
      else
        xmark(:) = .false.
        call range_checker(1,len(trim(x)),trim(x),xmark)
      endif
      if ( len(trim(y)) .eq. 0 ) then
        ymark(:) = .true.
      else
        ymark(:) = .false.
        call range_checker(2,len(trim(y)),trim(y),ymark)
      endif
      if ( len(trim(z)) .eq. 0 ) then
        zmark(:) = .true.
      else
        zmark(:) = .false.
        call range_checker(3,len(trim(z)),trim(z),zmark)
      endif

      mark(:) = .false.
      where(xmark .and. ymark .and. zmark) mark = .true.

      return
      end subroutine range_selector


!***********************************************************************


      subroutine range_checker(ix,nlen,txt,mark)

      integer(4),intent(in):: ix
      integer(4),intent(in):: nlen
      character(nlen),intent(in):: txt
      logical(1),intent(inout):: mark(natm)
      integer(4):: i,j
      real(4):: rmin,rmax
      character(nlen),allocatable:: array(:),array2(:)

      call name_separator(nlen,txt,",",array) ; mark(:) = .false.
      do i = 1,size(array)
        call name_separator(nlen,array(i),"<",array2)
        if ( len(trim(array2(1))) .eq. 0 ) then
          rmin = minval(cod(ix,:))
        else
          read(array2(1),*)rmin
        endif
        if ( len(trim(array2(2))) .eq. 0 ) then
          rmax = maxval(cod(ix,:))
        else
          read(array2(2),*)rmax
        endif
        do j = 1,natm
          if ( cod(ix,j).ge.rmin .and. cod(ix,j).le.rmax ) then
            mark(j) = .true.
          endif
        enddo
        deallocate(array2)
      enddo

      return
      end subroutine range_checker


!***********************************************************************


      end subroutine atom_specifier


!***********************************************************************


      subroutine make_specifier_list(n,natm,list,mark)

      integer(4),intent(in):: n,natm
      integer(4),intent(inout):: list(n)
      logical(1),intent(in):: mark(natm)

      integer(4):: i,ic

!***********************************************************************

      ic = 0
      do i = 1,natm
        if ( mark(i) ) then
          ic = ic + 1 ; list(ic) = i
        endif
      enddo

      return
      end subroutine make_specifier_list


!***********************************************************************


      subroutine output_specifier_log(n,loglvl,natm,anum,aname,rnum,   &
                                      rname,cnum,cname,l)

      integer(4),intent(in):: n
      character(1),intent(in):: loglvl
      integer(4),intent(in):: natm
      integer(4),intent(in):: anum(natm)
      character(*),intent(in):: aname(natm)
      integer(4),intent(in):: rnum(natm)
      character(*),intent(in):: rname(natm)
      integer(4),intent(in):: cnum(natm)
      character(*),intent(in):: cname(natm)
      integer(4),intent(in):: l(n)

      integer(4):: i,j
      character(999):: txt

!***********************************************************************

      write(6,'(a,i0)')"    Total number of atoms selected = ",n
      if ( loglvl.eq."l" .or. n.le.6 ) then
        do i = 1,n
          j = l(i)
          write(txt,'(3(i0,x,a,x))')cnum(j),trim(adjustl(cname(j))),   &
            rnum(j),trim(adjustl(rname(j))),anum(j),                   &
            trim(adjustl(aname(j)))
          write(6,'(2x,i8,a3,a)')i," : ",trim(txt)
        enddo
      elseif ( loglvl .eq. "s" ) then
        do i = 1,3
          j = l(i)
          write(txt,'(3(i0,x,a,x))')cnum(j),trim(adjustl(cname(j))),   &
            rnum(j),trim(adjustl(rname(j))),anum(j),                   &
            trim(adjustl(aname(j)))
          write(6,'(2x,i8,a3,a)')i," : ",trim(txt)
        enddo
        write(6,*)"              ..."
        do i = n-2,n
          j = l(i)
          write(txt,'(3(i0,x,a,x))')cnum(j),trim(adjustl(cname(j))),   &
            rnum(j),trim(adjustl(rname(j))),anum(j),                   &
            trim(adjustl(aname(j)))
          write(6,'(2x,i8,a3,a)')i," : ",trim(txt)
        enddo
      endif

      return
      end subroutine output_specifier_log
