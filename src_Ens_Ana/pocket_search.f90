
      subroutine pocket_search

      use COMIFN ; use COMVAL

      implicit none

      integer(4):: i,j,iary(3),ix,jx,iy,jy,iz,jz,mx,my,mz,ii,cx1,cx2,  &
        cy1,cy2,cz1,cz2,bx,by,bz,iii,jjj
      real(8):: c,rdist,r,r2,rtmp,tcod(3),rary(3),dx,dy,dz,ttcod(3)
      logical(4),allocatable:: cflg(:,:,:),cflg2(:,:,:),hole(:,:,:),   &
        hole2(:,:,:),hole3(:,:,:)
      logical(4):: tflg,ltc(nrec),ptc(nrec)

!**********************************************

      allocate(cflg(pcminX:pcmaxX,pcminY:pcmaxY,pcminZ:pcmaxZ))
      allocate(hole(pcminX:pcmaxX,pcminY:pcmaxY,pcminZ:pcmaxZ))
      ! Put atoms in a grid
      cflg = .false. ; rdist = 1.d0 / p_celsiz
      if ( bsize(1) .ne. 0.d0 ) then
        do i = 1,nrec
          j = iATMrec(i)
          ttcod(1:3) = cod(1:3,j) - bsize(1:3) *                       &
                      floor((cod(1:3,j)-bound(1,1:3))*ibsize(1:3))
          r = radrec(i)
          cx1 = 0 ; cx2 = 0 ; cy1 = 0 ; cy2 = 0 ; cz1 = 0 ; cz2 = 0
          if ( ttcod(1)-bound(1,1) .lt. r ) cx2 = 1
          if ( bound(2,1)-ttcod(1) .lt. r ) cx1 = -1
          if ( ttcod(2)-bound(1,2) .lt. r ) cy2 = 1
          if ( bound(2,2)-ttcod(2) .lt. r ) cy1 = -1
          if ( ttcod(3)-bound(1,3) .lt. r ) cz2 = 1
          if ( bound(2,3)-ttcod(3) .lt. r ) cz1 = -1
          do bz = cz1,cz2
            tcod(3) = ttcod(3) + bz*bsize(3)
          do by = cy1,cy2
            tcod(2) = ttcod(2) + by*bsize(2)
          do bx = cx1,cx2
            tcod(1) = ttcod(1) + bx*bsize(1)
            call pocket_sub
          enddo
          enddo
          enddo
        enddo
      else
        do i = 1,nrec
          j = iATMrec(i)
          tcod(1:3) = cod(1:3,j)
          r = radrec(i)
          call pocket_sub
        enddo
      endif

      ! Search pocket
      hole = .true.
      do iz = pcminZ,pcmaxZ
      do iy = pcminY,pcmaxY
      do ix = pcminX,pcmaxX
        if ( cflg(ix,iy,iz) ) hole(ix,iy,iz) = .false.
      enddo
      enddo
      enddo
      !! X
      do iz = pcminZ,pcmaxZ
      do iy = pcminY,pcmaxY
        do ix = pcminX,pcmaxX
          if ( cflg(ix,iy,iz) ) exit
          hole(ix,iy,iz) = .false.
        enddo
        do ix = pcmaxX,pcminX,-1
          if ( cflg(ix,iy,iz) ) exit
          hole(ix,iy,iz) = .false.
        enddo
      enddo
      enddo
      !! Y
      do iz = pcminZ,pcmaxZ
      do ix = pcminX,pcmaxX
        do iy = pcminY,pcmaxY
          if ( cflg(ix,iy,iz) ) exit
          hole(ix,iy,iz) = .false.
        enddo
        do iy = pcmaxY,pcminY,-1
          if ( cflg(ix,iy,iz) ) exit
          hole(ix,iy,iz) = .false.
        enddo
      enddo
      enddo
      !! Z
      do iy = pcminY,pcmaxY
      do ix = pcminX,pcmaxX
        do iz = pcminZ,pcmaxZ
          if ( cflg(ix,iy,iz) ) exit
          hole(ix,iy,iz) = .false.
        enddo
        do iz = pcmaxZ,pcminZ,-1
          if ( cflg(ix,iy,iz) ) exit
          hole(ix,iy,iz) = .false.
        enddo
      enddo
      enddo

      if ( nlig .gt. 0 ) then
        ! Put ligand atoms in a grid
        allocate(cflg2(pcminX:pcmaxX,pcminY:pcmaxY,pcminZ:pcmaxZ))
        allocate(hole2(pcminX:pcmaxX,pcminY:pcmaxY,pcminZ:pcmaxZ))
        cflg2 = .false.
        do i = 1,nlig
          j = iATMlig(i)
          tcod(1:3) = cod(1:3,j)
          r = radlig(i)
          ix = max(pcminX,nint((tcod(1)-r)*rdist))
          jx = min(pcmaxX,nint((tcod(1)+r)*rdist))
          iy = max(pcminY,nint((tcod(2)-r)*rdist))
          jy = min(pcmaxY,nint((tcod(2)+r)*rdist))
          iz = max(pcminZ,nint((tcod(3)-r)*rdist))
          jz = min(pcmaxZ,nint((tcod(3)+r)*rdist))
          rary(1:3) = tcod(1:3)*rdist
          iary(1:3) = nint(rary(1:3))
          r2 = r*r
          do mz = iz,jz
            ii = mz - iary(3)
            select case (ii)
            case (:-1)
              dz = (dble(mz)-rary(3)+0.5d0)*p_celsiz
            case (0)
              dz = 0.d0
            case (1:)
              dz = (dble(mz)-rary(3)-0.5d0)*p_celsiz
            end select
            dz = dz * dz
          do my = iy,jy
            ii = my - iary(2)
            select case (ii)
            case (:-1)
              dy = (dble(my)-rary(2)+0.5d0)*p_celsiz
            case (0)
              dy = 0.d0
            case (1:)
              dy = (dble(my)-rary(2)-0.5d0)*p_celsiz
            end select
            dy = dy * dy
          do mx = ix,jx
            ii = mx - iary(1)
            select case (ii)
            case (:-1)
              dx = (dble(mx)-rary(1)+0.5d0)*p_celsiz
            case (0)
              dx = 0.d0
            case (1:)
              dx = (dble(mx)-rary(1)-0.5d0)*p_celsiz
            end select
            dx = dx * dx
            if ( dx+dy+dz .le. r2 ) cflg2(mx,my,mz) = .true.
          enddo
          enddo
          enddo
        enddo

        hole2 = cflg2 ; i = 0
        do
          do iz = pcminZ,pcmaxZ
          do iy = pcminY,pcmaxY
          do ix = pcminX,pcmaxX
            if ( hole2(ix,iy,iz) ) then
              do jz = iz-1,pcminZ,-1
                if ( .not. hole(ix,iy,jz) ) exit
                hole2(ix,iy,jz) = .true.
              enddo
              do jz = iz+1,pcmaxZ
                if ( .not. hole(ix,iy,jz) ) exit
                hole2(ix,iy,jz) = .true.
              enddo
              do jy = iy-1,pcminY,-1
                if ( .not. hole(ix,jy,iz) ) exit
                hole2(ix,jy,iz) = .true.
              enddo
              do jy = iy+1,pcmaxY
                if ( .not. hole(ix,jy,iz) ) exit
                hole2(ix,jy,iz) = .true.
              enddo
              do jx = ix-1,pcminX,-1
                if ( .not. hole(jx,iy,iz) ) exit
                hole2(jx,iy,iz) = .true.
              enddo
              do jx = ix+1,pcmaxX
                if ( .not. hole(jx,iy,iz) ) exit
                hole2(jx,iy,iz) = .true.
              enddo
            endif
          enddo
          enddo
          enddo
          j = count(hole2)
          if ( i .eq. j ) then
            exit
          else
            i = j
          endif
        enddo
        where(cflg2) hole2 = .false.

        ! If you want to remove pocket entrance region
        if ( atom_number_for_pocket_entrance .ne. 0 ) then
          allocate(hole3(pcminX:pcmaxX,pcminY:pcmaxY,pcminZ:pcmaxZ))
          hole3 = .false.
          !! Find pocket entrance
          j = atom_number_for_pocket_entrance
          tcod(1:3) = cod(1:3,j)
          ix = nint(tcod(1)*rdist)
          iy = nint(tcod(2)*rdist)
          iz = nint(tcod(3)*rdist)
          !! Find the clothest hole2
          i = 1
          OUTER2 : do
            do jz = max(pcminZ,iz-i),max(pcminZ,iz+i)
            do jy = max(pcminY,iy-i),max(pcminY,iy+i)
            do jx = max(pcminx,ix-i),max(pcminX,ix+i)
              if ( hole2(jx,jy,jz) ) then
                mx = jx ; my = jy ; mz = jz ; exit OUTER2
              endif
            enddo
            enddo
            enddo
            i = i + 1
          enddo OUTER2
          !! Find pocket entrance
          hole3(mx,my,mz) = .true.
          do jz = mz-1,pcminZ,-1
            if ( .not. hole2(mx,my,jz) ) exit
            hole3(mx,my,jz) = .true.
          enddo
          do jz = mz+1,pcmaxZ
            if ( .not. hole2(mx,my,jz) ) exit
            hole3(mx,my,jz) = .true.
          enddo
          do jy = my-1,pcminY,-1
            if ( .not. hole2(mx,jy,mz) ) exit
            hole3(mx,jy,mz) = .true.
          enddo
          do jy = my+1,pcmaxY
            if ( .not. hole2(mx,jy,mz) ) exit
            hole3(mx,jy,mz) = .true.
          enddo
          do jx = mx-1,pcminX,-1
            if ( .not. hole2(jx,my,mz) ) exit
            hole3(jx,my,mz) = .true.
          enddo
          do jx = mx+1,pcmaxX
            if ( .not. hole2(jx,my,mz) ) exit
            hole3(jx,my,mz) = .true.
          enddo
          i = count(hole3)
          do
            do iz = pcminZ,pcmaxZ
            do iy = pcminY,pcmaxY
            do ix = pcminX,pcmaxX
              if ( hole3(ix,iy,iz) ) then
                do jz = iz-1,pcminZ,-1
                  if ( .not. hole2(ix,iy,jz) ) exit
                  hole3(ix,iy,jz) = .true.
                enddo
                do jz = iz+1,pcmaxZ
                  if ( .not. hole2(ix,iy,jz) ) exit
                  hole3(ix,iy,jz) = .true.
                enddo
                do jy = iy-1,pcminY,-1
                  if ( .not. hole2(ix,jy,iz) ) exit
                  hole3(ix,jy,iz) = .true.
                enddo
                do jy = iy+1,pcmaxY
                  if ( .not. hole2(ix,jy,iz) ) exit
                  hole3(ix,jy,iz) = .true.
                enddo
                do jx = ix-1,pcminX,-1
                  if ( .not. hole2(jx,iy,iz) ) exit
                  hole3(jx,iy,iz) = .true.
                enddo
                do jx = ix+1,pcmaxX
                  if ( .not. hole2(jx,iy,iz) ) exit
                  hole3(jx,iy,iz) = .true.
                enddo
              endif
            enddo
            enddo
            enddo
            j = count(hole3)
            if ( i .eq. j ) then
              exit
            else
              i = j
            endif
          enddo
          where(hole3) hole2 = .false.
        endif

        hole = hole2
      endif

      ! Pocket atom search
      if ( positive_method_flag ) then
        rdist = 1.d0 / p_celsiz
        ltc = .false. ; ptc = .false.
        if ( bsize(1) .ne. 0.d0 ) then
          do i = 1,nrec
            j = iATMrec(i)
            ttcod(1:3) = cod(1:3,j) - bsize(1:3) *                     &
                        floor((cod(1:3,j)-bound(1,1:3))*ibsize(1:3))
            r = radrec(i)
            cx1 = 0 ; cx2 = 0 ; cy1 = 0 ; cy2 = 0 ; cz1 = 0 ; cz2 = 0
            if ( ttcod(1)-bound(1,1) .lt. r ) cx2 = 1
            if ( bound(2,1)-ttcod(1) .lt. r ) cx1 = -1
            if ( ttcod(2)-bound(1,2) .lt. r ) cy2 = 1
            if ( bound(2,2)-ttcod(2) .lt. r ) cy1 = -1
            if ( ttcod(3)-bound(1,3) .lt. r ) cz2 = 1
            if ( bound(2,3)-ttcod(3) .lt. r ) cz1 = -1
            do bz = cz1,cz2
              tcod(3) = ttcod(3) + bz*bsize(3)
            do by = cy1,cy2
              tcod(2) = ttcod(2) + by*bsize(2)
            do bx = cx1,cx2
              tcod(1) = ttcod(1) + bx*bsize(1)
              call pocket_sub2
            enddo
            enddo
            enddo
          enddo
        else
          do i = 1,nrec
            j = iATMrec(i)
            tcod(1:3) = cod(1:3,j)
            r = radrec(i)
            call pocket_sub2
          enddo
        endif
        where(ptc) pkt_touch = pkt_touch + wfac
!        where(ltc) pkt_touch = pkt_touch - wfac
      endif

      where(hole) pcell = pcell + wfac
      r = p_celsiz
      write(upv,*)wfac,dble(count(hole))*r*r*r

!***********************************************

      return

!***********************************************

      contains

        subroutine pocket_sub

        ix = max(pcminX,nint((tcod(1)-r)*rdist))
        jx = min(pcmaxX,nint((tcod(1)+r)*rdist))
        iy = max(pcminY,nint((tcod(2)-r)*rdist))
        jy = min(pcmaxY,nint((tcod(2)+r)*rdist))
        iz = max(pcminZ,nint((tcod(3)-r)*rdist))
        jz = min(pcmaxZ,nint((tcod(3)+r)*rdist))
        rary(1:3) = tcod(1:3)*rdist
        iary(1:3) = nint(rary(1:3))
        r2 = r*r
        do mz = iz,jz
          ii = mz - iary(3)
          select case (ii)
          case (:-1)
            dz = (dble(mz)-rary(3)+0.5d0)*p_celsiz
          case (0)
            dz = 0.d0
          case (1:)
            dz = (dble(mz)-rary(3)-0.5d0)*p_celsiz
          end select
          dz = dz * dz
        do my = iy,jy
          ii = my - iary(2)
          select case (ii)
          case (:-1)
            dy = (dble(my)-rary(2)+0.5d0)*p_celsiz
          case (0)
            dy = 0.d0
          case (1:)
            dy = (dble(my)-rary(2)-0.5d0)*p_celsiz
          end select
          dy = dy * dy
        do mx = ix,jx
          ii = mx - iary(1)
          select case (ii)
          case (:-1)
            dx = (dble(mx)-rary(1)+0.5d0)*p_celsiz
          case (0)
            dx = 0.d0
          case (1:)
            dx = (dble(mx)-rary(1)-0.5d0)*p_celsiz
          end select
          dx = dx * dx
          if ( dx+dy+dz .le. r2 ) cflg(mx,my,mz) = .true.
        enddo
        enddo
        enddo

        return
        end subroutine pocket_sub

!*******************

        subroutine pocket_sub2

!        if ( ptc(i) ) return
!        ! ligand check
!        do ix = 1,nlig
!          iy = iATMlig(ix)
!          if ( RESnum(j) .eq. RESnum(iy) ) then
!            ptc(i) = .false. ; return
!          endif
!        enddo
!        if ( ltc(i) ) return
        ! ligand touch check
        do ix = 1,nlig
          iy = iATMlig(ix)
          ttcod(1:3) = cod(1:3,iy)
          dx = ttcod(1) - tcod(1)
          dy = ttcod(2) - tcod(2)
          dz = ttcod(3) - tcod(3)
          rtmp = radlig(ix)+r+watrad2
          if ( dx*dx+dy*dy+dz*dz .le. rtmp*rtmp ) then
            ltc(i) = .true. ; ptc(i) = .false.
!            jz = RESnum(j)
!            do iz = i+1,nrec
!              if ( RESnum(iATMrec(iz)) .ne. jz ) exit
!              ltc(iz) = .true. ; ptc(iz) = .false.
!            enddo
!            do iz = i-1,1,-1
!              if ( RESnum(iATMrec(iz)) .ne. jz ) exit
!              ltc(iz) = .true. ; ptc(iz) = .false.
!            enddo
            return
          endif
        enddo

        ix = max(pcminX,nint((tcod(1)-r)*rdist))
        jx = min(pcmaxX,nint((tcod(1)+r)*rdist))
        iy = max(pcminY,nint((tcod(2)-r)*rdist))
        jy = min(pcmaxY,nint((tcod(2)+r)*rdist))
        iz = max(pcminZ,nint((tcod(3)-r)*rdist))
        jz = min(pcmaxZ,nint((tcod(3)+r)*rdist))
        rary(1:3) = tcod(1:3)*rdist
        iary(1:3) = nint(rary(1:3))
        r2 = r*r
        OUTER : do mz = iz,jz
          ii = mz - iary(3)
          select case (ii)
          case (:-1)
            dz = (dble(mz)-rary(3)+0.5d0)*p_celsiz
          case (0)
            dz = 0.d0
          case (1:)
            dz = (dble(mz)-rary(3)-0.5d0)*p_celsiz
          end select
          dz = dz * dz
        do my = iy,jy
          ii = my - iary(2)
          select case (ii)
          case (:-1)
            dy = (dble(my)-rary(2)+0.5d0)*p_celsiz
          case (0)
            dy = 0.d0
          case (1:)
            dy = (dble(my)-rary(2)-0.5d0)*p_celsiz
          end select
          dy = dy * dy
        do mx = ix,jx
          ii = mx - iary(1)
          select case (ii)
          case (:-1)
            dx = (dble(mx)-rary(1)+0.5d0)*p_celsiz
          case (0)
            dx = 0.d0
          case (1:)
            dx = (dble(mx)-rary(1)-0.5d0)*p_celsiz
          end select
          dx = dx * dx
          if ( dx+dy+dz .le. r2 ) then
            if ( any(hole(max(pcminX,mx-1):min(pcmaxX,mx+1),           &
                          max(pcminY,my-1):min(pcmaxY,my+1),           &
                          max(pcminZ,mz-1):min(pcmaxZ,mz+1))) ) then
              ptc(i) = .true.
!              iii = RESnum(j)
!              do jjj = i+1,nrec
!                if ( RESnum(iATMrec(jjj)) .ne. iii ) exit
!                if ( .not. ltc(jjj) ) ptc(jjj) = .true.
!              enddo
!              do jjj = i-1,1,-1
!                if ( RESnum(iATMrec(jjj)) .ne. iii ) exit
!                if ( .not. ltc(jjj) ) ptc(jjj) = .true.
!              enddo
              exit OUTER
            endif
          endif
        enddo
        enddo
        enddo OUTER
!        pkt_touch(i) = pkt_touch(i) + wfac*c

        return
        end subroutine pocket_sub2

!*******************************************

      end subroutine pocket_search


!====================================================================


      subroutine pocket_search_init

      use COMIFN ; use COMVAL

      implicit none

      integer(4):: i,j,tlist(nATM)
      character(999):: tmp,tmp1,tmp2
      logical(1):: f

!*********************************************

      open(unit=upv,file=trim(d_proj)//".pktvol",status="replace")

      j = index(atom_spec_pocket,"_")
      if ( j .eq. 0 ) call error(10701)
      tmp1 = atom_spec_pocket(1:j-1) ; tmp2 = atom_spec_pocket(j+1:)
      write(6,'(4x,a)')"Receptor group : "//trim(tmp1)
      call atom_specifier(len(trim(tmp1)),trim(tmp1),nATM,ATMnum,ATMnm,&
        RESnum,RESnm,nCHN,CHN,Rcod,nrec,tlist)
      allocate(iATMrec(nrec)) ; iATMrec(:) = tlist(1:nrec)
      call output_specifier_log(nrec,loglvl,nATM,ATMnum,ATMnm,RESnum,  &
        RESnm,nCHN,CHN,iATMrec)
      allocate(radrec(nrec)) ; f = .false.
      if ( len(trim(input_topology)) .ne. 0 ) then
        call mk_vdWrad_list(nATM,nrec,iATMrec,RESnm,ATMnm,radrec,f)
      else
        call mk_vdWrad_list2(nATM,nrec,iATMrec,RESnm,ATMnm,radrec,f)
      endif
      write(6,'(4x,a)')"Ligand group : "//trim(tmp2)
      call atom_specifier(len(trim(tmp2)),trim(tmp2),nATM,ATMnum,ATMnm,&
        RESnum,RESnm,nCHN,CHN,Rcod,nlig,tlist)
      allocate(iATMlig(nlig)) ; iATMlig(:) = tlist(1:nlig)
      call output_specifier_log(nlig,loglvl,nATM,ATMnum,ATMnm,RESnum,  &
        RESnm,nCHN,CHN,iATMlig)
      allocate(radlig(nlig)) ; f = .true.
      if ( len(trim(input_topology)) .ne. 0 ) then
        call mk_vdWrad_list(nATM,nlig,iATMlig,RESnm,ATMnm,radlig,f)
      else
        call mk_vdWrad_list2(nATM,nlig,iATMlig,RESnm,ATMnm,radlig,f)
      endif

      pcminX = int(bound(1,1)/p_celsiz)-1
      pcminY = int(bound(1,2)/p_celsiz)-1
      pcminZ = int(bound(1,2)/p_celsiz)-1
      pcmaxX = int(bound(2,1)/p_celsiz)+1
      pcmaxY = int(bound(2,2)/p_celsiz)+1
      pcmaxZ = int(bound(2,3)/p_celsiz)+1
      NpcelX = pcmaxX - pcminX + 1
      NpcelY = pcmaxY - pcminY + 1
      NpcelZ = pcmaxZ - pcminZ + 1
      allocate(pcell(pcminX:pcmaxX,pcminY:pcmaxY,pcminZ:pcmaxZ))
      pcell = 0.d0

      if ( positive_method_flag ) then
        allocate(pkt_touch(nrec)) ; pkt_touch(:) = 0.d0
      endif

!*********************************************

      return
      end subroutine pocket_search_init


!====================================================================


      subroutine pocket_search_final

      use COMIFN ; use COMVAL

      implicit none

      integer(4):: i,j,k,l,ia(1)
      real(8),allocatable:: Tpkt_touch(:)
      logical(4),allocatable:: flg(:),flg2(:),flg3(:)
      real(8):: thre = 0.5d0
      real(8),allocatable:: Nscore(:)
      integer(4),allocatable:: cresatm(:)

!*********************************************

      pcell = pcell / sumW2
      open(unit=1,file=trim(d_proj)//"_pkt.xplor",status="replace")
      write(1,*)
      write(1,'(i8)')1
      write(1,'(a)')trim(project_name)//"_pkt"
      write(1,'(9i8)')NpcelX,pcminX,pcmaxX,NpcelY,pcminY,pcmaxY,       &
                      NpcelZ,pcminZ,pcmaxZ
      write(1,'(6e12.5)')dble(NpcelX)*p_celsiz,                &
        dble(NpcelY)*p_celsiz,dble(NpcelZ)*p_celsiz,   &
        90.d0,90.d0,90.d0
      write(1,'(a3)')"ZYX"
      do i = pcminZ,pcmaxZ
        write(1,'(i8)')i
       write(1,'(6e12.5)')pcell(pcminX:pcmaxX,pcminY:pcmaxY,i)
      enddo
      close(1) ; close(upv)

      if ( positive_method_flag ) then
        allocate(Tpkt_touch(istATM:ienATM))
        Tpkt_touch(:) = 0.d0
        pkt_touch(:) = pkt_touch(:) / sumW2
        ! ligand check
        do i = 1,nrec
          j = iATMrec(i)
          do k = 1,nlig
            if ( iATMlig(k) .eq. j ) then
              pkt_touch(i) = 0.d0 ; exit
            endif
          enddo
        enddo
        open(unit=1,file=trim(d_proj)//"_posm.pdb",status="replace")
        do i = 1,nrec
          j = iATMrec(i)
          if ( j.lt.istATM .or. j.gt.ienATM ) cycle
          Tpkt_touch(j) = pkt_touch(i)
        enddo
        do i = istATM,ienATM
          if ( RESnum(i) .lt. 10000 ) then
            write(1,'(a4,i7,x,a4,x,a4,a1,i4,4x,3f8.3,2f6.3)')          &
              "ATOM",i,ATMnm(i),RESnm(i),CHN(i),RESnum(i),             &
              cod(1:3,i),0.d0,Tpkt_touch(i)
          elseif ( RESnum(i) .lt. 100000 ) then
            write(1,'(a4,i7,x,a4,x,a4,a1,i5,3x,3f8.3,2f6.3)')          &
              "ATOM",i,ATMnm(i),RESnm(i),CHN(i),RESnum(i),             &
              cod(1:3,i),0.d0,Tpkt_touch(i)
          else
            write(1,'(a4,i7,x,a4,x,a4,a1,i6,2x,3f8.3,2f6.3)')          &
              "ATOM",i,ATMnm(i),RESnm(i),CHN(i),RESnum(i),             &
              cod(1:3,i),0.d0,Tpkt_touch(i)
          endif
        enddo
        close(1)

        open(unit=1,file=trim(d_proj)//"_posm.dat",status="replace")
        allocate(flg(nrec),flg2(nrec)) ; flg(1:nrec) = .true.
        do i = 1,nrec
          j = iATMrec(i)
          if ( pkt_touch(i) .eq. 0.d0 )                                &
!          if ( pkt_touch(i).eq.0.d0 .or. ATMnm(j).eq." N  " .or.       &
!               ATMnm(j).eq." CA " .or. ATMnm(j).eq." C  " .or.         &
!               ATMnm(j).eq." O  " .or. ATMnm(j).eq." H  " )            &
             flg(i) = .false.
        enddo
        do
          if ( count(flg(:)) .eq. 0 ) exit
          ia = maxloc(pkt_touch(:),flg(:))
          i = ia(1) ; flg2(:) = .false.
          j = iATMrec(i)
          do k = 1,nrec
            l = iATMrec(k)
            if ( RESnum(l).eq.RESnum(j) .and. flg(k) ) flg2(k) = .true.
          enddo
          do
            if ( count(flg2(:)) .eq. 0 ) exit
            ia = maxloc(pkt_touch(:),flg2(:))
            k = ia(1) ; flg(k) = .false. ; flg2(k) = .false.
            l = iATMrec(k)
            if ( pkt_touch(k) .ge. thre ) then
              write(1,'(i0,x,a4,x,a4,x,f6.3,x,a1)')RESnum(l),RESnm(l), &
                ATMnm(l),pkt_touch(k),"o"
            elseif ( pkt_touch(k) .le. -thre ) then
              write(1,'(i0,x,a4,x,a4,x,f6.3,x,a1)')RESnum(l),RESnm(l), &
                ATMnm(l),pkt_touch(k),"x"
            else
              write(1,'(i0,x,a4,x,a4,x,f6.3)')RESnum(l),RESnm(l),      &
                ATMnm(l),pkt_touch(k)
            endif
          enddo
          write(1,*)
        enddo
        close(1)

        open(unit=1,file=trim(d_proj)//"_posm.scr",status="replace")
        allocate(flg3(fstRES:fnlRES)) ; flg3 = .false.
        flg(1:nrec) = .true.
        allocate(Nscore(fstRES:fnlRES)) ; Nscore = 0.d0
        allocate(cresatm(fstRES:fnlRES)) ; cresatm = 0
        do i = 1,nrec
          j = iATMrec(i)
          if ( pkt_touch(i) .eq. 0.d0 ) then
!          if ( pkt_touch(i).eq.0.d0 .or. ATMnm(j).eq." N  " .or.       &
!               ATMnm(j).eq." CA " .or. ATMnm(j).eq." C  " .or.         &
!               ATMnm(j).eq." O  " .or. ATMnm(j).eq." H  " ) then
            flg(i) = .false.
          else
            Nscore(RESnum(j)) = Nscore(RESnum(j)) + pkt_touch(i)
            flg3(RESnum(j)) = .true.
          endif
          cresatm(RESnum(j)) = cresatm(RESnum(j)) + 1
        enddo
!        Nscore(:) = Nscore(:) / dble(cresatm(:))
        do
          if ( count(flg(:)) .eq. 0 ) exit
          ia = maxloc(Nscore(fstRES:fnlRES),flg3(fstRES:fnlRES))       &
               + fstRES - 1
          i = ia(1) ; flg2(:) = .false. ; flg3(i) = .false.
          write(1,'(i0,x,a4,x,f12.3)')i,aRES(i),Nscore(i)
          do j = 1,nrec
            l = iATMrec(j)
            if ( i.eq.RESnum(l) .and. flg(j) ) flg2(j) = .true.
          enddo
          do
            if ( count(flg2(:)) .eq. 0 ) exit
            ia = maxloc(pkt_touch(:),flg2(:))
            k = ia(1) ; flg(k) = .false. ; flg2(k) = .false.
            l = iATMrec(k)
            write(1,'(5x,a4,x,f6.3)')ATMnm(l),pkt_touch(k)
          enddo
          write(1,*)
        enddo
        close(1)
      endif

!***********************************************

      return
      end subroutine pocket_search_final
