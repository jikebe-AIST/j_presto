
      subroutine celpol(ilflag,iprint)

!***********************************************************
!
!     Calculate CMM parameters
!
!************************************************************

      use COMPAR ; use COMCMM
      use COMBAS,only: ixnatm,fxchrg
      use COMCMMC,only: cord
      !$ use omp_lib

      implicit none

      integer(4),intent(in):: ilflag,iprint

      integer(4):: ncou(nfrag)
      integer(4):: i,ii,jj,ix,iy,iz,num,iat,mx,my,mz
      real(8):: qz,rx,ry,rz,rxx,ryy,rzz,rxy,ryz,rzx,cc
      real(8):: gx,gy,gz,dx,dy,dz,dxx,dyy,dzz,qq,xx,yy,zz
      real(8):: ax,ay,az,xax,yay,zaz

!************************************************************
!  Generating level 1

      ncou(1:nfrag) = 0
      !$OMP parallel default(none)                                   & !
      !$OMP private(iz,iy,ix,i,num,qz,rx,ry,rz,rxx,ryy,rzz,rxy,ryz,  & !
      !$OMP         rzx,gx,gy,gz,ii,iat,cc,dx,dy,dz,dxx,dyy,dzz,jj,  & !
      !$OMP         mx,my,mz,xx,yy,zz,ax,ay,az,xax,yay,zaz,qq)       & !
      !$OMP shared(nv,nfrag,npcl1,c1,ipcl,fxchrg,cord,CMMpot1,lvdn,  & !
      !$OMP        npcl2,CMMpot2,npcl3,CMMpot3,npcl4,CMMpot4,npcl5,  & !
      !$OMP        CMMpot5,nlev,c2,c3,c4,c5,ilflag,iprint)           & !
      !$OMP reduction (+ : ncou)
      !$OMP do collapse(3) schedule (static)
      do iz = 1,nv(1)
      do iy = 1,nv(1)
      do ix = 1,nv(1)
        do i = 1,nfrag
          num = npcl1(i,ix,iy,iz)
          if ( num .ne. 0 ) then
            qz = 0.d0
            rx = 0.d0  ; ry = 0.d0  ; rz = 0.d0
            rxx = 0.d0 ; ryy = 0.d0 ; rzz = 0.d0
            rxy = 0.d0 ; ryz = 0.d0 ; rzx = 0.d0
            gx = c1(1,ix) ; gy = c1(2,iy) ; gz = c1(3,iz)

            do ii = 1,num
              iat = ipcl(ii,i,ix,iy,iz)
              cc = fxchrg(iat)

              ! Mono-pole.
              qz = qz + cc
              ! Di-pole.
              dx = cord(1,iat) - gx
              dy = cord(2,iat) - gy
              dz = cord(3,iat) - gz
              rx = rx + cc*dx ; ry = ry + cc*dy ; rz = rz + cc*dz
              ! Quadru-pole.
              dxx = dx*dx ; dyy = dy*dy ; dzz = dz*dz

              rxx = rxx + cc*(2*dxx - dyy  - dzz)
              ryy = ryy + cc*(2*dyy - dzz  - dxx)
              rzz = rzz + cc*(2*dzz - dxx  - dyy)
              rxy = rxy + cc*dx*dy
              ryz = ryz + cc*dy*dz
              rzx = rzx + cc*dz*dx
            enddo

            CMMpot1(1:13,i,ix,iy,iz) =                                 &
              (/qz,rx,ry,rz,rxx*0.5d0,ryy*0.5d0,rzz*0.5d0,             &
                rxx,ryy,rzz,rxy*3.d0,ryz*3.d0,rzx*3.d0/)
            if ( qz .ne. 0.d0 ) ncou(i) = ncou(i) + 1

          else
            CMMpot1(:,i,ix,iy,iz) = 0.d0
          endif
        enddo
      enddo
      enddo
      enddo
      !$OMP enddo

!****************
!  Level 1 ---> Level 2.

      if ( nlev .lt. 2 ) goto 9999
      !$OMP do collapse(3) schedule (static)
      do iz = 1,nv(2)
      do iy = 1,nv(2)
      do ix = 1,nv(2)
        do i = 1,nfrag
          num = npcl2(i,ix,iy,iz)
          if ( num .ne. 0 ) then
            qz = 0.d0
            rx = 0.d0  ; ry = 0.d0  ; rz = 0.d0
            rxx = 0.d0 ; ryy = 0.d0 ; rzz = 0.d0
            rxy = 0.d0 ; ryz = 0.d0 ; rzx = 0.d0
   
            do jj = 1,8

              select case (jj)
                case (1)
                  mx = lvdn(1,ix) ; my = lvdn(1,iy) ; mz = lvdn(1,iz)
                case (2)
                  mx = lvdn(2,ix) ; my = lvdn(1,iy) ; mz = lvdn(1,iz)
                case (3)
                  mx = lvdn(1,ix) ; my = lvdn(2,iy) ; mz = lvdn(1,iz)
                case (4)
                  mx = lvdn(2,ix) ; my = lvdn(2,iy) ; mz = lvdn(1,iz)
                case (5)
                  mx = lvdn(1,ix) ; my = lvdn(1,iy) ; mz = lvdn(2,iz)
                case (6)
                  mx = lvdn(2,ix) ; my = lvdn(1,iy) ; mz = lvdn(2,iz)
                case (7)
                  mx = lvdn(1,ix) ; my = lvdn(2,iy) ; mz = lvdn(2,iz)
                case (8)
                  mx = lvdn(2,ix) ; my = lvdn(2,iy) ; mz = lvdn(2,iz)
              end select

              num = npcl1(i,mx,my,mz)
              if ( num .eq. 0 ) cycle

              qq = CMMpot1(1,i,mx,my,mz)
              qz = qz + qq

              xx = c1(1,mx) - c2(1,ix)
              yy = c1(2,my) - c2(2,iy)
              zz = c1(3,mz) - c2(3,iz)

              ax = CMMpot1(2,i,mx,my,mz)
              ay = CMMpot1(3,i,mx,my,mz)
              az = CMMpot1(4,i,mx,my,mz)

              rx = rx + ax + qq*xx
              ry = ry + ay + qq*yy
              rz = rz + az + qq*zz

              dxx = xx*xx ; dyy = yy*yy ; dzz = zz*zz
              xax = xx*ax ; yay = yy*ay ; zaz = zz*az
              rxx = rxx + CMMpot1(5,i,mx,my,mz)                        &
                        + 2*xax - yay - zaz                            &
                        + qq*(dxx - 0.5d0*(dyy + dzz))
              ryy = ryy + CMMpot1(6,i,mx,my,mz)                        &
                        + 2*yay - zaz - xax                            &
                        + qq*(dyy - 0.5d0*(dzz + dxx))
              rzz = rzz + CMMpot1(7,i,mx,my,mz)                        &
                        + 2*zaz - xax - yay                            &
                        + qq*(dzz - 0.5d0*(dxx + dyy))
              rxy = rxy + CMMpot1(11,i,mx,my,mz) * 0.5d0               &
                        + (yy*ax + xx*ay + xx*yy*qq)*1.5d0
              ryz = ryz + CMMpot1(12,i,mx,my,mz) * 0.5d0               &
                        + (zz*ay + yy*az + yy*zz*qq)*1.5d0
              rzx = rzx + CMMpot1(13,i,mx,my,mz) * 0.5d0               &
                        + (xx*az + zz*ax + zz*xx*qq)*1.5d0
            enddo

            CMMpot2(1:13,i,ix,iy,iz) =                                 &
              (/qz,rx,ry,rz,rxx,ryy,rzz,rxx*2.d0,ryy*2.d0,rzz*2.d0,    &
                rxy*2.d0,ryz*2.d0,rzx*2.d0/)
          else
            CMMpot2(:,i,ix,iy,iz) = 0.d0
          endif
        enddo
      enddo
      enddo
      enddo
      !$OMP end do

!*******************************
!  Level 2 ---> Level 3.

      if ( nlev .lt. 3 ) goto 9999
      !$OMP do collapse(3) schedule(static)
      do iz = 1,nv(3)
      do iy = 1,nv(3)
      do ix = 1,nv(3)
        do i = 1,nfrag
          num = npcl3(i,ix,iy,iz)
          if ( num .ne. 0 ) then
            qz = 0.d0
            rx = 0.d0  ; ry = 0.d0  ; rz = 0.d0
            rxx = 0.d0 ; ryy = 0.d0 ; rzz = 0.d0
            rxy = 0.d0 ; ryz = 0.d0 ; rzx = 0.d0
   
            do jj = 1,8

              select case (jj)
                case (1)
                  mx = lvdn(1,ix) ; my = lvdn(1,iy) ; mz = lvdn(1,iz)
                case (2)
                  mx = lvdn(2,ix) ; my = lvdn(1,iy) ; mz = lvdn(1,iz)
                case (3)
                  mx = lvdn(1,ix) ; my = lvdn(2,iy) ; mz = lvdn(1,iz)
                case (4)
                  mx = lvdn(2,ix) ; my = lvdn(2,iy) ; mz = lvdn(1,iz)
                case (5)
                  mx = lvdn(1,ix) ; my = lvdn(1,iy) ; mz = lvdn(2,iz)
                case (6)
                  mx = lvdn(2,ix) ; my = lvdn(1,iy) ; mz = lvdn(2,iz)
                case (7)
                  mx = lvdn(1,ix) ; my = lvdn(2,iy) ; mz = lvdn(2,iz)
                case (8)
                  mx = lvdn(2,ix) ; my = lvdn(2,iy) ; mz = lvdn(2,iz)
              end select

              num = npcl2(i,mx,my,mz)
              if ( num .eq. 0 ) cycle

              qq = CMMpot2(1,i,mx,my,mz)
              qz = qz + qq

              xx = c2(1,mx) - c3(1,ix)
              yy = c2(2,my) - c3(2,iy)
              zz = c2(3,mz) - c3(3,iz)

              ax = CMMpot2(2,i,mx,my,mz)
              ay = CMMpot2(3,i,mx,my,mz)
              az = CMMpot2(4,i,mx,my,mz)

              rx = rx + ax + qq*xx
              ry = ry + ay + qq*yy
              rz = rz + az + qq*zz

              dxx = xx*xx ; dyy = yy*yy ; dzz = zz*zz
              xax = xx*ax ; yay = yy*ay ; zaz = zz*az
              rxx = rxx + CMMpot2(5,i,mx,my,mz)                        &
                        + 2*xax - yay - zaz                            &
                        + qq*(dxx - 0.5d0*(dyy + dzz))
              ryy = ryy + CMMpot2(6,i,mx,my,mz)                        &
                        + 2*yay - zaz - xax                            &
                        + qq*(dyy - 0.5d0*(dzz + dxx))
              rzz = rzz + CMMpot2(7,i,mx,my,mz)                        &
                        + 2*zaz - xax - yay                            &
                        + qq*(dzz - 0.5d0*(dxx + dyy))
              rxy = rxy + CMMpot2(11,i,mx,my,mz) * 0.5d0               &
                        + (yy*ax + xx*ay + xx*yy*qq)*1.5d0
              ryz = ryz + CMMpot2(12,i,mx,my,mz) * 0.5d0               &
                        + (zz*ay + yy*az + yy*zz*qq)*1.5d0
              rzx = rzx + CMMpot2(13,i,mx,my,mz) * 0.5d0               &
                        + (xx*az + zz*ax + zz*xx*qq)*1.5d0
            enddo

            CMMpot3(1:13,i,ix,iy,iz) =                                 &
              (/qz,rx,ry,rz,rxx,ryy,rzz,rxx*2.d0,ryy*2.d0,rzz*2.d0,    &
                rxy*2.d0,ryz*2.d0,rzx*2.d0/)
          else
            CMMpot3(:,i,ix,iy,iz) = 0.d0
          endif
        enddo
      enddo
      enddo
      enddo
      !$OMP end do

!****************
!  Level 3 ---> Level 4.

      if ( nlev .lt. 4 ) goto 9999
      !$OMP do collapse(3) schedule(static)
      do iz = 1,nv(4)
      do iy = 1,nv(4)
      do ix = 1,nv(4)
        do i = 1,nfrag
          num = npcl4(i,ix,iy,iz)
          if ( num .ne. 0 ) then
            qz = 0.d0
            rx = 0.d0  ; ry = 0.d0  ; rz = 0.d0
            rxx = 0.d0 ; ryy = 0.d0 ; rzz = 0.d0
            rxy = 0.d0 ; ryz = 0.d0 ; rzx = 0.d0

            do jj = 1,8

              select case (jj)
                case (1)
                  mx = lvdn(1,ix) ; my = lvdn(1,iy) ; mz = lvdn(1,iz)
                case (2)
                  mx = lvdn(2,ix) ; my = lvdn(1,iy) ; mz = lvdn(1,iz)
                case (3)
                  mx = lvdn(1,ix) ; my = lvdn(2,iy) ; mz = lvdn(1,iz)
                case (4)
                  mx = lvdn(2,ix) ; my = lvdn(2,iy) ; mz = lvdn(1,iz)
                case (5)
                  mx = lvdn(1,ix) ; my = lvdn(1,iy) ; mz = lvdn(2,iz)
                case (6)
                  mx = lvdn(2,ix) ; my = lvdn(1,iy) ; mz = lvdn(2,iz)
                case (7)
                  mx = lvdn(1,ix) ; my = lvdn(2,iy) ; mz = lvdn(2,iz)
                case (8)
                  mx = lvdn(2,ix) ; my = lvdn(2,iy) ; mz = lvdn(2,iz)
              end select

              num = npcl3(i,mx,my,mz)
              if ( num .eq. 0 ) cycle

              qq = CMMpot3(1,i,mx,my,mz)
              qz = qz + qq

              xx = c3(1,mx) - c4(1,ix)
              yy = c3(2,my) - c4(2,iy)
              zz = c3(3,mz) - c4(3,iz)

              ax = CMMpot3(2,i,mx,my,mz)
              ay = CMMpot3(3,i,mx,my,mz)
              az = CMMpot3(4,i,mx,my,mz)

              rx = rx + ax + qq*xx
              ry = ry + ay + qq*yy
              rz = rz + az + qq*zz

              dxx = xx*xx ; dyy = yy*yy ; dzz = zz*zz
              xax = xx*ax ; yay = yy*ay ; zaz = zz*az
              rxx = rxx + CMMpot3(5,i,mx,my,mz)                        &
                        + 2*xax - yay - zaz                            &
                        + qq*(dxx - 0.5d0*(dyy + dzz))
              ryy = ryy + CMMpot3(6,i,mx,my,mz)                        &
                        + 2*yay - zaz - xax                            &
                        + qq*(dyy - 0.5d0*(dzz + dxx))
              rzz = rzz + CMMpot3(7,i,mx,my,mz)                        &
                        + 2*zaz - xax - yay                            &
                        + qq*(dzz - 0.5d0*(dxx + dyy))
              rxy = rxy + CMMpot3(11,i,mx,my,mz) * 0.5d0               &
                        + (yy*ax + xx*ay + xx*yy*qq)*1.5d0
              ryz = ryz + CMMpot3(12,i,mx,my,mz) * 0.5d0               &
                        + (zz*ay + yy*az + yy*zz*qq)*1.5d0
              rzx = rzx + CMMpot3(13,i,mx,my,mz) * 0.5d0               &
                        + (xx*az + zz*ax + zz*xx*qq)*1.5d0
            enddo

            CMMpot4(1:13,i,ix,iy,iz) =                                 &
              (/qz,rx,ry,rz,rxx,ryy,rzz,rxx*2.d0,ryy*2.d0,rzz*2.d0,    &
                rxy*2.d0,ryz*2.d0,rzx*2.d0/)
          else
            CMMpot4(:,i,ix,iy,iz) = 0.d0
          endif
        enddo
      enddo
      enddo
      enddo
      !$OMP end do

!****************
!  Level 4 ---> Level 5.

      if ( nlev .lt. 5 ) goto 9999
      !$OMP do collapse(3) schedule(static)
      do iz = 1,nv(5)
      do iy = 1,nv(5)
      do ix = 1,nv(5)
        do i = 1,nfrag
          num = npcl5(i,ix,iy,iz)
          if ( num .ne. 0 ) then
            qz = 0.d0
            rx = 0.d0  ; ry = 0.d0  ; rz = 0.d0
            rxx = 0.d0 ; ryy = 0.d0 ; rzz = 0.d0
            rxy = 0.d0 ; ryz = 0.d0 ; rzx = 0.d0

            do jj = 1,8

              select case (jj)
                case (1)
                  mx = lvdn(1,ix) ; my = lvdn(1,iy) ; mz = lvdn(1,iz)
                case (2)
                  mx = lvdn(2,ix) ; my = lvdn(1,iy) ; mz = lvdn(1,iz)
                case (3)
                  mx = lvdn(1,ix) ; my = lvdn(2,iy) ; mz = lvdn(1,iz)
                case (4)
                  mx = lvdn(2,ix) ; my = lvdn(2,iy) ; mz = lvdn(1,iz)
                case (5)
                  mx = lvdn(1,ix) ; my = lvdn(1,iy) ; mz = lvdn(2,iz)
                case (6)
                  mx = lvdn(2,ix) ; my = lvdn(1,iy) ; mz = lvdn(2,iz)
                case (7)
                  mx = lvdn(1,ix) ; my = lvdn(2,iy) ; mz = lvdn(2,iz)
                case (8)
                  mx = lvdn(2,ix) ; my = lvdn(2,iy) ; mz = lvdn(2,iz)
              end select

              num = npcl4(i,mx,my,mz)
              if ( num .eq. 0 ) cycle

              qq = CMMpot4(1,i,mx,my,mz)
              qz = qz + qq

              xx = c4(1,mx) - c5(1,ix)
              yy = c4(2,my) - c5(2,iy)
              zz = c4(3,mz) - c5(3,iz)

              ax = CMMpot4(2,i,mx,my,mz)
              ay = CMMpot4(3,i,mx,my,mz)
              az = CMMpot4(4,i,mx,my,mz)

              rx = rx + ax + qq*xx
              ry = ry + ay + qq*yy
              rz = rz + az + qq*zz

              dxx = xx*xx ; dyy = yy*yy ; dzz = zz*zz
              xax = xx*ax ; yay = yy*ay ; zaz = zz*az
              rxx = rxx + CMMpot4(5,i,mx,my,mz)                        &
                        + 2*xax - yay - zaz                            &
                        + qq*(dxx - 0.5d0*(dyy + dzz))
              ryy = ryy + CMMpot4(6,i,mx,my,mz)                        &
                        + 2*yay - zaz - xax                            &
                        + qq*(dyy - 0.5d0*(dzz + dxx))
              rzz = rzz + CMMpot4(7,i,mx,my,mz)                        &
                        + 2*zaz - xax - yay                            &
                        + qq*(dzz - 0.5d0*(dxx + dyy))
              rxy = rxy + CMMpot4(11,i,mx,my,mz) * 0.5d0               &
                        + (yy*ax + xx*ay + xx*yy*qq)*1.5d0
              ryz = ryz + CMMpot4(12,i,mx,my,mz) * 0.5d0               &
                        + (zz*ay + yy*az + yy*zz*qq)*1.5d0
              rzx = rzx + CMMpot4(13,i,mx,my,mz) * 0.5d0               &
                        + (xx*az + zz*ax + zz*xx*qq)*1.5d0
            enddo

            CMMpot5(1:13,i,ix,iy,iz) =                                 &
              (/qz,rx,ry,rz,rxx,ryy,rzz,rxx*2.d0,ryy*2.d0,rzz*2.d0,    &
                rxy*2.d0,ryz*2.d0,rzx*2.d0/)
          else
            CMMpot5(:,i,ix,iy,iz) = 0.d0
          endif
        enddo
      enddo
      enddo
      enddo
      !$OMP end do
9999  continue
      !$OMP end parallel
      if ( ilflag .eq. 2 ) then
        write(iprint,*)"  "
        write(iprint,*)" No. OF NON-NEUTRAL CELLS AT LEVEL-1"
        do i = 1,nfrag
          write(iprint,'(a12,i1,a3,i0)')"    CLUSTER-",i," : ",ncou(i)
        enddo
        write(iprint,*)"  "
      endif

!*************************************
!     For external CMM

!!      if ( extCMM_flag ) then
!!        CMMpot5(1:13,1:nv(1),1:nv(1),1:nv(1),nfrag) =                      &
!!     &     CMMpot5(1:13,1:nv(1),1:nv(1),1:nv(1),nfrag) +                   &
!!     &     eCMMpot5(1:13,1:nv(1),1:nv(1),1:nv(1))
!!        CMMpot4(1:13,1:nv(2),1:nv(2),1:nv(2),nfrag) =                   &
!!     &     CMMpot4(1:13,1:nv(2),1:nv(2),1:nv(2),nfrag) +                &
!!     &     eCMMpot4(1:13,1:nv(2),1:nv(2),1:nv(2))
!!        CMMpot3(1:13,1:nv(3),1:nv(3),1:nv(3),nfrag) =                   &
!!     &     CMMpot3(1:13,1:nv(3),1:nv(3),1:nv(3),nfrag) +                &
!!     &     eCMMpot3(1:13,1:nv(3),1:nv(3),1:nv(3))
!!        CMMpot2(1:13,1:nv(4),1:nv(4),1:nv(4),nfrag) =                   &
!!     &     CMMpot2(1:13,1:nv(4),1:nv(4),1:nv(4),nfrag) +                &
!!     &     eCMMpot2(1:13,1:nv(4),1:nv(4),1:nv(4))
!!        CMMpot1(1:13,1:nv(5),1:nv(5),1:nv(5),nfrag) =                   &
!!     &     CMMpot1(1:13,1:nv(5),1:nv(5),1:nv(5),nfrag) +                &
!!     &     eCMMpot1(1:13,1:nv(5),1:nv(5),1:nv(5))
!!      endif
!!      call system_clock(tim2)
!!      tmptime = tmptime + tim2 - tim1

!*************************************

      return
      end subroutine celpol
