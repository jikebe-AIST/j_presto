

      subroutine measure(txt)

      use COMVAL ; use COMIFN,only: loglvl

      implicit none

      character(999),intent(in):: txt

      integer(4):: i,j,k,ic,list(nATM),iseg(4)
      character(999):: tmp,part,seg,tmp2

!**********************************************************************

      n_dist_c = 0 ; n_ang_c = 0 ; n_dih_c = 0 ; tmp = txt
      do
        i = index(tmp,";")
        if ( i .eq. 0 ) then
          part = tmp ; tmp = ""
        else
          part = tmp(:i-1) ; tmp = tmp(i+1:)
        endif
        ic = 0
        do
          j = index(part,"_")
          if ( j .eq. 0 ) exit
          ic = ic + 1
          part = part(j+1:)
        enddo
        select case(ic)
          case(1)
            n_dist_c = n_dist_c + 1
          case(2)
            n_ang_c = n_ang_c + 1
          case(3)
            n_dih_c = n_dih_c + 1
        end select
        if ( len(trim(tmp)) .eq. 0 ) exit
      enddo
      allocate(dist_c(2,n_dist_c),ang_c(4,n_ang_c),dih_c(4,n_dih_c))
      n_dist_c = 0 ; n_ang_c = 0 ; n_dih_c = 0 ; tmp = txt
      do
        i = index(tmp,";")
        if ( i .eq. 0 ) then
          part = tmp ; tmp = ""
        else
          part = tmp(:i-1) ; tmp = tmp(i+1:)
        endif
        ic = 0
        do
          j = index(part,"_")
          if ( j .eq. 0 ) then
            seg = part ; part = ""
            call atom_specifier(len(trim(seg)),trim(seg),nATM,   &
              ATMnum,ATMnm,RESnum,RESnm,nCHN,CHN,Rcod,k,list)
            if ( k .ne. 1 ) exit
            iseg(ic+1) = list(1)
          else
            ic = ic + 1
            seg = part(:j-1) ; part = part(j+1:)
            call atom_specifier(len(trim(seg)),trim(seg),nATM,   &
              ATMnum,ATMnm,RESnum,RESnm,nCHN,CHN,Rcod,k,list)
            if ( k .ne. 1 ) cycle
            iseg(ic) = list(1)
          endif
          if ( len(trim(part)) .eq. 0 ) exit
        enddo
        select case(ic)
          case(1)
            n_dist_c = n_dist_c + 1 ; dist_c(:,n_dist_c) = iseg(1:2)
          case(2)
            n_ang_c = n_ang_c + 1
            ang_c(:,n_ang_c) = (/iseg(3), iseg(1), iseg(2), iseg(3)/)
          case(3)
            n_dih_c = n_dih_c + 1 ; dih_c(:,n_dih_c) = iseg(1:4)
        end select
        if ( len(trim(tmp)) .eq. 0 ) exit
      enddo

      if ( n_dist_c .gt. 0 ) then
        write(6,'(2x,a)')"* distance analysis is performed"
        write(6,'(4x,a,i0)')"Total number of distance pairs = ",n_dist_c
        if ( loglvl.eq."l" .or. n_dist_c.le.6 ) then
          do i = 1,n_dist_c
            do j = 1,2
              k = dist_c(j,i)
              write(tmp2,'(3(i0,x,a,x))')nCHN(k),trim(CHN(k)),         &
                RESnum(k),trim(RESnm(k)),ATMnum(k),trim(adjustl(ATMnm(k)))
              if ( j .ne. 1 ) then 
                tmp = trim(tmp)//" - ["//trim(tmp2)//"]"
              else
                tmp = "["//trim(tmp2)//"]"
              endif
            enddo
            write(6,'(i8,a3,a)')i," : ",trim(tmp)
          enddo
        elseif ( loglvl .eq. "s" ) then 
          do i = 1,3
            do j = 1,2
              k = dist_c(j,i)
              write(tmp2,'(3(i0,x,a,x))')nCHN(k),trim(CHN(k)),         &
                RESnum(k),trim(RESnm(k)),ATMnum(k),trim(adjustl(ATMnm(k)))
              if ( j .ne. 1 ) then 
                tmp = trim(tmp)//" - ["//trim(tmp2)//"]"
              else
                tmp = "["//trim(tmp2)//"]"
              endif
            enddo
            write(6,'(i8,a3,a)')i," : ",trim(tmp)
          enddo
          write(6,*)"              ..."
          do i = n_dist_c-2,n_dist_c
            do j = 1,2
              k = dist_c(j,i)
              write(tmp2,'(3(i0,x,a,x))')nCHN(k),trim(CHN(k)),         &
                RESnum(k),trim(RESnm(k)),ATMnum(k),trim(adjustl(ATMnm(k)))
              if ( j .ne. 1 ) then 
                tmp = trim(tmp)//" - ["//trim(tmp2)//"]"
              else
                tmp = "["//trim(tmp2)//"]"
              endif
            enddo
            write(6,'(i8,a3,a)')i," : ",trim(tmp)
          enddo
        endif
      endif
      if ( n_ang_c .gt. 0 ) then
        write(6,'(2x,a)')"* angle analysis is performed"
        write(6,'(4x,a,i0)')"Total number of angle pairs = ",n_ang_c
        if ( loglvl.eq."l" .or. n_ang_c.le.6 ) then
          do i = 1,n_ang_c
            do j = 1,3
              k = ang_c(j,i)
              write(tmp2,'(3(i0,x,a,x))')nCHN(k),trim(CHN(k)),         &
                RESnum(k),trim(RESnm(k)),ATMnum(k),trim(adjustl(ATMnm(k)))
              if ( j .ne. 1 ) then 
                tmp = trim(tmp)//" - ["//trim(tmp2)//"]"
              else
                tmp = "["//trim(tmp2)//"]"
              endif
            enddo
            write(6,'(i8,a3,a)')i," : ",trim(tmp)
          enddo
        elseif ( loglvl .eq. "s" ) then 
          do i = 1,3
            do j = 1,3
              k = ang_c(j,i)
              write(tmp2,'(3(i0,x,a,x))')nCHN(k),trim(CHN(k)),         &
                RESnum(k),trim(RESnm(k)),ATMnum(k),trim(adjustl(ATMnm(k)))
              if ( j .ne. 1 ) then 
                tmp = trim(tmp)//" - ["//trim(tmp2)//"]"
              else
                tmp = "["//trim(tmp2)//"]"
              endif
            enddo
            write(6,'(i8,a3,a)')i," : ",trim(tmp)
          enddo
          write(6,*)"              ..."
          do i = n_ang_c-2,n_ang_c
            do j = 1,3
              k = ang_c(j,i)
              write(tmp2,'(3(i0,x,a,x))')nCHN(k),trim(CHN(k)),         &
                RESnum(k),trim(RESnm(k)),ATMnum(k),trim(adjustl(ATMnm(k)))
              if ( j .ne. 1 ) then 
                tmp = trim(tmp)//" - ["//trim(tmp2)//"]"
              else
                tmp = "["//trim(tmp2)//"]"
              endif
            enddo
            write(6,'(i8,a3,a)')i," : ",trim(tmp)
          enddo
        endif
      endif
      if ( n_dih_c .gt. 0 ) then
        write(6,'(2x,a)')"* dihedral angle analysis is performed"
        write(6,'(4x,a,i0)')"Total number of dihedral angle pairs = ", &
          n_dih_c
        if ( loglvl.eq."l" .or. n_dih_c.le.6 ) then
          do i = 1,n_dih_c
            do j = 1,4
              k = dih_c(j,i)
              write(tmp2,'(3(i0,x,a,x))')nCHN(k),trim(CHN(k)),         &
                RESnum(k),trim(RESnm(k)),ATMnum(k),trim(adjustl(ATMnm(k)))
              if ( j .ne. 1 ) then 
                tmp = trim(tmp)//" - ["//trim(tmp2)//"]"
              else
                tmp = "["//trim(tmp2)//"]"
              endif
            enddo
            write(6,'(i8,a3,a)')i," : ",trim(tmp)
          enddo
        elseif ( loglvl .eq. "s" ) then 
          do i = 1,3
            do j = 1,4
              k = dih_c(j,i)
              write(tmp2,'(3(i0,x,a,x))')nCHN(k),trim(CHN(k)),         &
                RESnum(k),trim(RESnm(k)),ATMnum(k),trim(adjustl(ATMnm(k)))
              if ( j .ne. 1 ) then 
                tmp = trim(tmp)//" - ["//trim(tmp2)//"]"
              else
                tmp = "["//trim(tmp2)//"]"
              endif
            enddo
            write(6,'(i8,a3,a)')i," : ",trim(tmp)
          enddo
          write(6,*)"              ..."
          do i = n_dih_c-2,n_dih_c
            do j = 1,4
              k = dih_c(j,i)
              write(tmp2,'(3(i0,x,a,x))')nCHN(k),trim(CHN(k)),         &
                RESnum(k),trim(RESnm(k)),ATMnum(k),trim(adjustl(ATMnm(k)))
              if ( j .ne. 1 ) then 
                tmp = trim(tmp)//" - ["//trim(tmp2)//"]"
              else
                tmp = "["//trim(tmp2)//"]"
              endif
            enddo
            write(6,'(i8,a3,a)')i," : ",trim(tmp)
          enddo
        endif
      endif

!**********************************************************************

      return
      end subroutine measure
