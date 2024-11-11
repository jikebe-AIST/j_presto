
      subroutine celprp(iprint,ilflag,ier)

!******************************************
!     
!     Prepare cell information
!
!******************************************

      use COMPAR ; use COMBAS ; use COMERG ; use COMCMM ; use COMCMMC
      !$ use omp_lib

      implicit none

      integer(4),intent(in):: iprint
      integer(4),intent(in):: ilflag
      integer(4),intent(inout):: ier

      integer(4):: ix,iy,iz,num(nfrag),i,j,nel,nelvdw,ii,iat1,nv1,     &
                   ipcl_el(ipmax),ipcl_elvdw(ipmax),k,ix1,iy1,iz1,jy,  &
                   jz,jx

!**************************************

      ier = 0
      if ( ilflag .eq. 2 ) then
        write(iprint,*)" "
        write(iprint,*)"INFORMATION> CELPRP"
        write(iprint,*)"             CELLS ARE PREPARED."
        write(iprint,*)" "
      endif

      select case (ixfbou)
      case (1)
        call celpar_PB(ilflag,iprint,ier)
      case default
        call celpar_CB(ilflag,iprint,ier)
      end select
      if ( ier .ne. 0 ) then
        write(iprint,*)"ERROR> CELPRP"
        write(iprint,*)"    CELSET ERROR IN CELPAR HAS OCCURED."
        return
      endif

!***********************

      ! Reorder ipcl so that atoms w/o EL and with EL are contiguous
      ! within each group.
      nv1 = nv(1)
      !$OMP parallel default(none)                                   & !
      !$OMP private(iz,iy,ix,num,i,nel,nelvdw,ii,iat1,ipcl_el,       & !
      !$OMP         ipcl_elvdw)                                      & !
      !$OMP shared(nv1,npcl1,nfrag,ipcl,zero_vdW,cell_nEL)
      !$OMP do collapse(3)
      do iz = 1,nv1
      do iy = 1,nv1
      do ix = 1,nv1
        num(1:nfrag) = npcl1(1:nfrag,ix,iy,iz)

        do i = 1,nfrag
          ! The order of ipcl is not important
          ! Thus we split ipcl to two parts: EL, and EL+vdW
          nel = 0 ; nelvdw = 0
          do ii = 1,num(i)
            iat1 = ipcl(ii,i,ix,iy,iz)
            if ( zero_vdW(iat1) ) then
              nel = nel + 1 ; ipcl_el(nel) = iat1
            else
              nelvdw = nelvdw + 1 ; ipcl_elvdw(nelvdw) = iat1
            endif
          enddo
          ipcl(1:nel,i,ix,iy,iz) = ipcl_el(1:nel)
          ipcl(nel+1:num(i),i,ix,iy,iz) = ipcl_elvdw(1:nelvdw)
          cell_nEL(i,ix,iy,iz) = nel
        enddo
      enddo
      enddo
      enddo
      !$OMP end do
      !$OMP end parallel

      ! Make inum array
      ii = 0
      do iz = 1,nv1,2
        iz1 = iz+1
      do iy = 1,nv1,2
        iy1 = iy+1
      do ix = 1,nv1,2
        ix1 = ix+1
        do jz = iz,iz1
        do jy = iy,iy1
        do jx = ix,ix1
        do i = 1,nfrag
        do j = 1,npcl1(i,jx,jy,jz)
          ii = ii + 1
          inum(ipcl(j,i,jx,jy,jz)) = ii
        enddo
        enddo
        enddo
        enddo
        enddo
      enddo
      enddo
      enddo

      !! For high-level parallelization
      if ( high_para .gt. 0 ) then
        np_cell(:) = 0
        do k = 1,high_para
          iz1 = (k-1) / id2p1_2
          jz = iz1*id2p1_2
          iy1 = (k-1-jz) / id2p1
          jy = jz + iy1*id2p1
          ix1 = k-jy ; iy1 = iy1+1 ; iz1 = iz1+1

          do iz = iz1,nv1,idelay+1
          do iy = iy1,nv1,id2p1
          do ix = ix1,nv1,id2p1
            if ( all(npcl1(:,ix,iy,iz) .eq. 0) ) cycle
            np_cell(k) = np_cell(k) + 1
            CHK1 : do i = 1,nfrag
            do j = 1,npcl1(i,ix,iy,iz)
              para_se(1,np_cell(k),k) = inum(ipcl(j,i,ix,iy,iz))
              exit CHK1
            enddo
            enddo CHK1
            CHK2 : do i = nfrag,1,-1
            do j = npcl1(i,ix,iy,iz),1,-1
              para_se(2,np_cell(k),k) = inum(ipcl(j,i,ix,iy,iz))
              exit CHK2
            enddo
            enddo CHK2
          enddo
          enddo
          enddo
        enddo
      endif

      !$OMP parallel default (none)                                  & !
      !$OMP private(i)                                               & !
      !$OMP shared(Tixatyp,Tcord,ixnatm,inum,ixatyp,cord,chgmod)
      !$OMP do schedule (static)
      do i = 1,ixnatm
        Tixatyp(inum(i)) = ixatyp(i)
        Tcord(1:3,inum(i)) = cord(1:3,i)
        Tcord(4,inum(i)) = chgmod(i)
      enddo
      !$OMP end do 
      !$OMP end parallel

      if ( iy15m .eq. 1 ) then
        call celoth
        call celpol(ilflag,iprint)
        call neartb(ier)
      else
        select case(ixfbou)
        case (1)
          call cutoff_PB(iprint,ilflag,ier)
        case default
          call cutoff_CB(iprint,ilflag,ier)
        end select
      endif

      if ( iyeflg(15) .eq. 1 ) call repprp

      if ( ier .ne. 0 ) then
        write(iprint,*)"ERROR> CELPRP"
        write(iprint,*)"  NEARTB ERROR IN CMMPRP HAS OCCURED."
        return
      endif

!**************************

      return
      end subroutine celprp
