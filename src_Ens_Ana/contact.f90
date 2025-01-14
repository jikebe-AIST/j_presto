
      subroutine contact(icn)

!*******************************************************
!
!     Do contact analysis
!
!*******************************************************

      use COMIFN ; use COMVAL
!      use morph_f
      !$ use omp_lib

      implicit none

      integer(4),intent(in):: icn

      real(4):: tcod(3,nprob)
      real(4):: rtmp
      integer(4):: iminX,imaxX,iminY,imaxY,iminZ,imaxZ
      integer(4):: iiminX,iimaxX,iiminY,iimaxY,iiminZ,iimaxZ
      integer(4),allocatable:: grid(:,:,:,:)
      logical(4):: con_flg(ntgt)
      integer(4):: i,iatm,n,x,y,z,j,ii,icod(3),icod2(3),nn(ntgt)
      real(4):: tcod2(3),dcod(3),d
      real(8),allocatable:: ASAcod(:,:),ASAvdW(:)
      real(8):: tASA(1),ASA(ntgt),rtmp2

!*********************************

      ! make tcod
      do i = 1,nprob
        iatm = iATMprob(i)
        tcod(1:3,i) = cod(1:3,iatm)
      enddo

      ! make grid
      rtmp = 1.0 / cellsz_con
      iminX = nint(minval(tcod(1,1:nprob))*rtmp)
      imaxX = nint(maxval(tcod(1,1:nprob))*rtmp)
      iminY = nint(minval(tcod(2,1:nprob))*rtmp)
      imaxY = nint(maxval(tcod(2,1:nprob))*rtmp)
      iminZ = nint(minval(tcod(3,1:nprob))*rtmp)
      imaxZ = nint(maxval(tcod(3,1:nprob))*rtmp)
800   allocate(grid(0:np_cell,iminX:imaxX,iminY:imaxY,iminZ:imaxZ))
      grid(0,:,:,:) = 0

      do i = 1,nprob
        icod(1:3) = nint(tcod(1:3,i)*rtmp)
        n = grid(0,icod(1),icod(2),icod(3)) + 1
        if ( n .gt. np_cell ) then
          np_cell = np_cell + 1
          deallocate(grid) ; goto 800
        endif
        grid(n,icod(1),icod(2),icod(3)) = i
        grid(0,icod(1),icod(2),icod(3)) = n
      enddo

!***************************************

      ! Contact, burial & ASA analysis
      !$OMP parallel default (none)                                  & !
      !$OMP private(i,tcod2,icod2,iiminX,iimaxX,iiminY,iimaxY,iiminZ,& !
      !$OMP         iimaxZ,j,ASAcod,ASAvdW,ii,rtmp2,x,y,z,iatm,dcod, & !
      !$OMP         d)                                               & !
      !$OMP shared(ntgt,con_flg,nn,cod,iATMtgt,rtmp,iminX,imaxX,     & !
      !$OMP        iminY,imaxY,iminZ,imaxZ,grid,ASA,radtgt,watrad,   & !
      !$OMP        tcod,radprob)
      !$OMP do schedule (guided)
      do i = 1,ntgt
        con_flg(i) = .false. ; nn(i) = 0
        tcod2(1:3) = cod(1:3,iATMtgt(i))
        icod2(1:3) = nint(tcod2(1:3)*rtmp)
        iiminX = max(icod2(1)-1,iminX)
        iimaxX = min(icod2(1)+1,imaxX)
        iiminY = max(icod2(2)-1,iminY)
        iimaxY = min(icod2(2)+1,imaxY)
        iiminZ = max(icod2(3)-1,iminZ)
        iimaxZ = min(icod2(3)+1,imaxZ)
        j = sum(grid(0,iiminX:iimaxX,iiminY:iimaxY,iiminZ:iimaxZ))
        if ( j .eq. 0 ) then
          ASA(i) = 4.d0*npi*(radtgt(i)+watrad)*(radtgt(i)+watrad)
        else
          allocate(ASAcod(3,0:j),ASAvdW(0:j))
          ASAcod(1:3,0) = tcod2(1:3)
          ASAvdW(0) = radprob(i)
          ii = 0 ; rtmp2 = ASAvdW(0)+watrad*2.d0
          do z = iiminZ,iimaxZ
          do y = iiminY,iimaxY
          do x = iiminX,iimaxX
            do j = 1,grid(0,x,y,z)
              iatm = grid(j,x,y,z)
              dcod(1:3) = tcod2(1:3) - tcod(1:3,iatm)
              d = dcod(1)*dcod(1) + dcod(2)*dcod(2) + dcod(3)*dcod(3)
              if ( d .lt. 0.00001 ) cycle
              if ( d .lt. (radprob(iatm)+rtmp2)**2 ) then
                ii = ii + 1 ; con_flg(i) = .true.
                ASAcod(1:3,ii) = tcod(1:3,iatm)
                ASAvdW(ii) = radprob(iatm)
              endif
            enddo
          enddo
          enddo
          enddo
          nn(i) = ii
#ifdef ASA
          call calc_morphs_selective(ii+1,ASAcod(1:3,0:ii),            &
               ASAvdW(0:ii),watrad,1,(/1/),ASA(i))
#endif
          deallocate(ASAcod,ASAvdW)
        endif
      enddo
      !$OMP end do
      !$OMP end parallel

      where(con_flg) rcon(:) = rcon(:) + wfac

      ! correlation
      if ( contact_correlation_flag ) then
        do i = 1,ntgt-1
          if ( con_flg(i) ) then
            do j = i+1,ntgt
              if ( con_flg(j) ) then
                con_cor(j,i) = con_cor(j,i) + wfac
              else
                con_cor(j,i) = con_cor(j,i) - wfac
              endif
            enddo
          else
            do j = i+1,ntgt
              if ( con_flg(j) ) then
                con_cor(j,i) = con_cor(j,i) - wfac
              else
                con_cor(j,i) = con_cor(j,i) + wfac
              endif
            enddo
          endif
        enddo
      endif

!**************************************

#ifdef ASA
      ! Output
      write(usf,'(i0,2(x,f),2(x,i0))')icn,wfac,sum(ASA),sum(nn),       &
                                      count(con_flg)
      if ( ntgt .le. 100 ) then
        do i = 1,ntgt
          write(ubu_e+i,'(i0,x,f,2(x,f8.3),x,i0)')icn,wfac,ASA(i), &
            100.d0*ASA(i)/(4.d0*npi*(vdWcon2(i)+watrad)**2),nn(i)
        enddo
      endif
      ave_burial(:) = ave_burial(:) + dble(nn(:))*wfac
      ave_surface(:) = ave_surface(:) + ASA(:)*wfac
#endif

!**************************************

      return
      end subroutine contact


!=========================================================================


      subroutine contact_init()

      use COMIFN ; use COMVAL

      implicit none

      integer(4):: i,j,k,ierr
      real(8):: rtmp
      character(4):: tATMnm
      character(130):: tmp,tmp2
      logical(4):: ex

      allocate(rcon(ntgt))
      if ( contact_correlation_flag ) then
        allocate(con_cor(ntgt,ntgt)) ; con_cor(:,:) = 0.d0
      endif
      cellsz_con = maxval(radprob)+maxval(radtgt)+watrad*2.d0
      np_cell = int(0.2*(cellsz_con**3))+1

#ifdef ASA
      open(unit=ubu,file=trim(d_proj)//".bur",status="replace")
      open(unit=usf,file=trim(d_proj)//".sur",status="replace")
      allocate(ave_burial(nCON2),ave_surface(nCON2))
      ave_burial(:) = 0.0 ; ave_surface(:) = 0.d0

      if ( ntgt .le. 100 ) then
        inquire(file=trim(d_proj)//"_BUR",exist=ex)
        !! over-write
        if ( ex ) then
          write(6,'(4x,a)')"! CAUTION !"
          write(6,'(4x,a)')"The directory for burial info. ( "//       &
            trim(d_proj)//"_BUR ) is OVERWRITEN"
          call system("rm -rf "//trim(d_proj)//                        &
                      "_BUR/*.bur 2> /dev/null")
        !! Make the directory for sur files
        else
          call system("mkdir "//trim(d_proj)//"_BUR 2> /dev/null")
        endif

        write(tmp2,'(a)')trim(d_proj)//"_BUR/"
        do i = 1,ntgt
          j = iATMtgt(i)
          tATMnm = trim(adjustl(ATMnm(j)))
          write(tmp,'(a,i0,a)')trim(tmp2),RESnum(j),"_"//              &
                trim(adjustl(RESnm(j)))//"_"//tATMnm
          open(unit=i+ubu_e,file=trim(tmp)//".bur",status="replace")
        enddo
      endif
#endif

      return
      end subroutine contact_init


!=========================================================================


      subroutine contact_final(swei)

      use COMIFN ; use COMVAL

      implicit none

      real(8),intent(in):: swei

      integer(4),parameter:: nbin_per = 101
      integer(4):: i,j,k,ist,ic,max_iburi
      real(8):: rtmp,rtmp2
      real(4):: afac,bfac
      character(4):: tATMnm
      character(130):: tmp,tmp2
 
      integer(4),allocatable:: iburi(:)
      real(4),allocatable:: wburi(:),proba(:),ASA(:),ASApercent(:)

!*********************************************

      ! Contact
      open(unit=1,file=trim(d_proj)//"_con.pdb",status="replace")
      j = 1
      rcon(:) = rcon(:) / sumW2
      rtmp2 = maxval(rcon)
      if ( rtmp2 .gt. 0.d0 ) then
        rtmp2 = 1.d0 / rtmp2
      else
        rtmp2 = 1.d0
      endif
      do i = istATM,ienATM
        if ( j .gt. ntgt ) then
          rtmp = 0.d0
        else
          if ( i .eq. iATMtgt(j) ) then
            rtmp = rcon(j) ; j = j + 1
          else
            rtmp = 0.d0
          endif
        endif
        if ( RESnum(i) .lt. 10000 ) then
          write(1,'(a4,i7,x,a4,x,a4,a1,i4,4x,3f8.3,2f6.3)')          &
               "ATOM",i,ATMnm(i),RESnm(i),CHN(i),                &
               RESnum(i),Rcod(1:3,i),rtmp*rtmp2,rtmp
        elseif ( RESnum(i) .lt. 100000 ) then
          write(1,'(a4,i7,x,a4,x,a4,a1,i5,3x,3f8.3,2f6.3)')          &
               "ATOM",i,ATMnm(i),RESnm(i),CHN(i),                &
               RESnum(i),Rcod(1:3,i),rtmp*rtmp2,rtmp
        else
          write(1,'(a4,i7,x,a4,x,a4,a1,i6,2x,3f8.3,2f6.3)')          &
               "ATOM",i,ATMnm(i),RESnm(i),CHN(i),                &
               RESnum(i),Rcod(1:3,i),rtmp*rtmp2,rtmp
        endif
      enddo
      close(1)

      !! correlation
      if ( contact_correlation_flag ) then
        con_cor = con_cor / sumW2
        forall(i=1:ntgt) con_cor(i,i) = 1.d0
        do i = 1,ntgt-1
        do j = i,ntgt
          con_cor(i,j) = con_cor(j,i)
        enddo
        enddo
        open(unit=1,file=trim(d_proj)//".con_cor",status="replace")
        do i = 1,ntgt
          do j = 1,ntgt
            write(1,*)i,j,con_cor(j,i)
          enddo
          write(1,*)
        enddo
        close(1) 
      endif

#ifdef ASA
      ! Burial & ASA
      ave_burial(:) = ave_burial(:) / swei
      ave_surface(:) = ave_surface(:) / swei
      do i = 1,ntgt
        j = iATMtgt(i)
!        write(tmp,'(i0,a)')RESnum(j),"_"//trim(adjustl(RESnm(j)))//"_" &
!          //trim(adjustl(ATMnm(j)))
        write(ubu,'(i0,2(x,a),4(x,f8.3))')RESnum(j),RESnm(j),ATMnm(j), &
          ave_surface(i),                                              &
          100.d0*ave_surface(i)/(4.d0*npi*(vdWcon2(i)+watrad)**2),     &
          ave_burial(i),rcon(i)
      enddo

      write(usf,'(a,e12.5)')"# average ASA                = ",         &
                           sum(ave_surface)
      write(usf,'(a,e12.5)')"# average burial             = ",         &
                           sum(ave_burial)
      write(usf,'(a,e12.5)')"# average number of contacts = ",sum(rcon)

      close(ubu) ; close(usf)
      if ( ntgt .le. 100 ) then
        do i = ubu_e+1,ubu_e+ntgt
          close(i)
        enddo
        !! make burial graph
        write(tmp2,'(a)')trim(d_proj)//"_BUR/"
        do i = 1,ntgt
          j = iATMtgt(i)
          tATMnm = trim(adjustl(ATMnm(j)))
          write(tmp,'(a,i0,a)')trim(tmp2),RESnum(j),"_"//              &
               trim(adjustl(RESnm(j)))//"_"//tATMnm
          !!! count columns
          call systemqq("wc -l "//trim(tmp)//".bur > .ahocount")
          open(unit=i+ubu_e,file=".ahocount",status="old")
          read(i+ubu_e,*)ic
          allocate(iburi(ic),wburi(ic),ASA(ic),ASApercent(ic))
          close(unit=i+ubu_e,status="delete")
          !!! read burial data       
          open(unit=i+ubu_e,file=trim(tmp)//".bur",status="old")
          do j = 1,ic
            read(i+ubu_e,*)k,wburi(j),ASA(j),ASApercent(j),iburi(j)
          enddo
          close(i+ubu_e)
          !!! make graph
          max_iburi = maxval(iburi(:))
          allocate(proba(0:max_iburi)) ; proba(:) = 0.0
          do j = 1,ic
            k = iburi(j)
            proba(k) = proba(k) + wburi(j)
          enddo
          rtmp = sum(proba) ; proba(:) = proba(:) / rtmp
          open(unit=i+ubu_e,file=trim(tmp)//".graph",status="replace")
          do j = 0,max_iburi
            write(i+ubu_e,*)j,proba(j)
          enddo
          close(i+ubu_e)
          deallocate(iburi,proba)

          allocate(proba(0:100)) ; proba(:) = 0.0
          do j = 1,ic
            k = nint(ASApercent(j))
            proba(k) = proba(k) + wburi(j)
          enddo
          proba(:) = proba(:) / rtmp
          open(unit=i+ubu_e,file=trim(tmp)//".sur_graph",              &
               status="replace")
          do j = 0,100
            write(i+ubu_e,*)real(j),proba(j)
          enddo
          close(i+ubu_e)
          deallocate(wburi,proba,ASA,ASApercent)
        enddo

      endif

      open(unit=1,file=trim(d_proj)//"_bur.pdb",status="replace")
      ist = 1
      do i = istATM,ienATM
        afac = 0.0 ; bfac = 0.0
        do j = ist,ntgt
          if ( iATMtgt(j) .eq. i ) then
            afac = ave_surface(j)/(4.d0*npi*(vdWcon2(j)+watrad)**2)
            bfac = ave_burial(j) ; ist = j+1 ; exit
          endif
        enddo
        if ( RESnum(i) .lt. 10000 ) then
          write(1,'(a4,i7,x,a4,x,a4,a1,i4,4x,3f8.3,2f6.3)')"ATOM",i,   &
            ATMnm(i),RESnm(i),CHN(i),RESnum(i),Rcod(1:3,i),afac,bfac
        elseif ( RESnum(i) .lt. 100000 ) then
          write(1,'(a4,i7,x,a4,x,a4,a1,i5,3x,3f8.3,2f6.3)')"ATOM",i,   &
            ATMnm(i),RESnm(i),CHN(i),RESnum(i),Rcod(1:3,i),afac,bfac
        else
          write(1,'(a4,i7,x,a4,x,a4,a1,i6,2x,3f8.3,2f6.3)')"ATOM",i,   &
            ATMnm(i),RESnm(i),CHN(i),RESnum(i),Rcod(1:3,i),afac,bfac
        endif
      enddo
#endif

!***************************************

      return
      end subroutine contact_final
