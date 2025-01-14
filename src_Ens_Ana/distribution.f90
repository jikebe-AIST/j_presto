
      subroutine distribution

      use COMIFN ; use COMVAL ; use COMTPL,only: fxchrg

      implicit none

      integer(4):: i,j,iary(3),ix,jx,iy,jy,iz,jz,mx,my,mz,ii,cx1,cx2,  &
        cy1,cy2,cz1,cz2,bx,by,bz
      real(8):: rdist,r,tcod(3),ttcod(3),rary(3),dx,dy,dz
      logical(4),allocatable:: cflg(:,:,:)

!**********************************************

      allocate(cflg(icminX:icmaxX,icminY:icmaxY,icminZ:icmaxZ))
      cflg = .false. ; rdist = 1.d0 / d_celsiz
      if ( bsize(1) .ne. 0.d0 ) then
        do i = 1,ndist
          j = iATMdist(i)
          ttcod(1:3) = cod(1:3,j) - bsize(1:3) *                       &
                      floor((cod(1:3,j)-bound(1,1:3))*ibsize(1:3))
          r = raddist(i)
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
            call dist_sub
          enddo
          enddo
          enddo
          rary(1:3) = ttcod(1:3)*rdist
          iary(1:3) = nint(rary(1:3))
          if ( iary(1).lt.icminX .or. iary(1).gt.icmaxX .or.           &
               iary(2).lt.icminY .or. iary(2).gt.icmaxY .or.           &
               iary(3).lt.icminZ .or. iary(3).gt.icmaxZ ) cycle
          dcell_ele(iary(1),iary(2),iary(3)) =                         &
            dcell_ele(iary(1),iary(2),iary(3)) + fxchrg(j)*wfac
        enddo
      else
        do i = 1,ndist
          j = iATMdist(i)
          tcod(1:3) = cod(1:3,j)
          r = raddist(i)
          call dist_sub
          if ( iary(1).lt.icminX .or. iary(1).gt.icmaxX .or.           &
               iary(2).lt.icminY .or. iary(2).gt.icmaxY .or.           &
               iary(3).lt.icminZ .or. iary(3).gt.icmaxZ ) cycle
          dcell_ele(iary(1),iary(2),iary(3)) =                         &
            dcell_ele(iary(1),iary(2),iary(3)) + fxchrg(j)*wfac
        enddo
      endif
      where(cflg) dcell = dcell + wfac
      deallocate(cflg)

!***********************************************

      return

!***********************************************

      contains

        subroutine dist_sub

        ix = max(icminX,nint((tcod(1)-r)*rdist))
        jx = min(icmaxX,nint((tcod(1)+r)*rdist))
        iy = max(icminY,nint((tcod(2)-r)*rdist))
        jy = min(icmaxY,nint((tcod(2)+r)*rdist))
        iz = max(icminZ,nint((tcod(3)-r)*rdist))
        jz = min(icmaxZ,nint((tcod(3)+r)*rdist))
        rary(1:3) = tcod(1:3)*rdist
        iary(1:3) = nint(rary(1:3))
        r = r*r
        do mz = iz,jz
          ii = mz - iary(3)
          select case (ii)
          case (:-1)
            dz = (dble(mz)-rary(3)+0.5d0)*d_celsiz
          case (0)
            dz = 0.d0
          case (1:)
            dz = (dble(mz)-rary(3)-0.5d0)*d_celsiz
          end select
        do my = iy,jy
          ii = my - iary(2)
          select case (ii)
          case (:-1)
            dy = (dble(my)-rary(2)+0.5d0)*d_celsiz
          case (0)
            dy = 0.d0
          case (1:)
            dy = (dble(my)-rary(2)-0.5d0)*d_celsiz
          end select
        do mx = ix,jx
          ii = mx - iary(1)
          select case (ii)
          case (:-1)
            dx = (dble(mx)-rary(1)+0.5d0)*d_celsiz
          case (0)
            dx = 0.d0
          case (1:)
            dx = (dble(mx)-rary(1)-0.5d0)*d_celsiz
          end select
          if ( dx*dx+dy*dy+dz*dz .le. r ) cflg(mx,my,mz) = .true.
        enddo
        enddo
        enddo

        return
        end subroutine dist_sub

!***********************************************

      end subroutine distribution


!====================================================================


      subroutine distrib_init

      use COMIFN ; use COMVAL

      implicit none

      integer(4):: i,j,k,l,tlist(nATM),ierr
      character(4):: tATMnm
      logical(1):: f = .true.

!*********************************************

      call atom_specifier(len(trim(atom_spec_distrib)),                &
        trim(atom_spec_distrib),nATM,ATMnum,ATMnm,RESnum,RESnm,nCHN,   &
        CHN,Rcod,ndist,tlist)
      allocate(iATMdist(ndist)) ; iATMdist(:) = tlist(1:ndist)
      call output_specifier_log(ndist,loglvl,nATM,ATMnum,ATMnm,RESnum, &
        RESnm,nCHN,CHN,iATMdist)

      ! Make vdW arraies
      allocate(raddist(ndist))
      if ( len(trim(input_topology)) .ne. 0 ) then
        call mk_vdWrad_list(nATM,ndist,iATMdist,RESnm,ATMnm,raddist,f)
        write(6,'(4x,a)')"Radius info. in the input tpl file is used"
      else
        call mk_vdWrad_list(nATM,ndist,iATMdist,RESnm,ATMnm,raddist,f)
      endif

      icminX = int(cbound(1,1)/d_celsiz)-1
      icminY = int(cbound(1,2)/d_celsiz)-1
      icminZ = int(cbound(1,3)/d_celsiz)-1
      icmaxX = int(cbound(2,1)/d_celsiz)+1
      icmaxY = int(cbound(2,2)/d_celsiz)+1
      icmaxZ = int(cbound(2,3)/d_celsiz)+1
      NdcelX = icmaxX - icminX + 1
      NdcelY = icmaxY - icminY + 1
      NdcelZ = icmaxZ - icminZ + 1
      allocate(dcell(icminX:icmaxX,icminY:icmaxY,icminZ:icmaxZ))
      allocate(dcell_ele(icminX:icmaxX,icminY:icmaxY,icminZ:icmaxZ))
      dcell = 0.d0 ; dcell_ele = 0.d0

!*********************************************

      return
      end subroutine distrib_init


!====================================================================


      subroutine distrib_final

      use COMIFN ; use COMVAL

      implicit none

      integer(4):: i,j,k

!*********************************************

      dcell = dcell / sumW2 ; dcell_ele = dcell_ele / sumW2
      open(unit=1,file=trim(d_proj)//"_dist.xplor",status="replace")
      write(1,*)
      write(1,'(i8)')1
      write(1,'(a)')trim(project_name)//"_dist"
      write(1,'(9i8)')NdcelX,icminX,icmaxX,NdcelY,icminY,icmaxY,       &
                      NdcelZ,icminZ,icmaxZ
      write(1,'(6e12.5)')dble(NdcelX)*d_celsiz,                  &
        dble(NdcelY)*d_celsiz,dble(NdcelZ)*d_celsiz,       &
        90.d0,90.d0,90.d0
      write(1,'(a3)')"ZYX"
      do i = icminZ,icmaxZ
        write(1,'(i8)')i
       write(1,'(6e12.5)')dcell(icminX:icmaxX,icminY:icmaxY,i)
      enddo
      close(1)

      open(unit=1,file=trim(d_proj)//"_dist.dx",status="replace")
      write(1,'(a,3(x,i0))')"object 1 class gridpositions counts",     &
        NdcelX,NdcelY,NdcelZ
      write(1,'(a,3(x,f8.3))')"origin",                                &
        (dble(icminX)-0.5d0)*d_celsiz,                           &
        (dble(icminY)-0.5d0)*d_celsiz,                           &
        (dble(icminZ)-0.5d0)*d_celsiz
      write(1,'(a,3(x,f8.3))')"delta",d_celsiz,0.d0,0.d0
      write(1,'(a,3(x,f8.3))')"delta",0.d0,d_celsiz,0.d0
      write(1,'(a,3(x,f8.3))')"delta",0.d0,0.d0,d_celsiz
      write(1,'(a,3(x,i0))')"object 2 class gridconnections counts"    &
        ,NdcelX,NdcelY,NdcelZ
      write(1,'(a,i0,a)')"object 3 class array type double rank "//    &
        "0 items ",NdcelX*NdcelY*NdcelZ," data follows"
      write(1,'(3(e12.5,x))')(((dcell_ele(i,j,k),                      &
        k=icminZ,icmaxZ),j=icminY,icmaxY),i=icminX,icmaxX)
      write(1,'(a)')'attribute "dep" string "positions"'
      write(1,'(a)')'object "regular positions regular '//             &
       'connections" class field'
      write(1,'(a)')'component "positions" value 1'
      write(1,'(a)')'component "connections" value 2'
      write(1,'(a)')'component "data" value 3'
      close(1)

!***********************************************

      return
      end subroutine distrib_final
