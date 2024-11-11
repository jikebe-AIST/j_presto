
      subroutine repprp

      use COMBAS ; use COMCMM ; use COMCMMC

      implicit none

      integer(4):: i,j,n,mx,my,mz,ii,jj,iat1,iat2,ice,nce,ix,iy,iz
      real(8):: xx,yy,zz,rr,gx,gy,gz,d(3)

!*****************************************************

      ! Put atom in the grid
      zz = siz(1)*0.5d0
      xx = c1(1,1) - zz
      yy = c1(2,1) - zz
      zz = c1(3,1) - zz
      nrepcl(1:nv(1),1:nv(1),1:nv(1),1:2) = 0
      rr = 1.d0 / siz(1)

      do j = 1,Nrep(1)
        i = replst(j,1)
        gx = cord(1,i)-fxcell(1)*floor((cord(1,i)-celwal(1))*invcel(1))
        gy = cord(2,i)-fxcell(2)*floor((cord(2,i)-celwal(3))*invcel(2))
        gz = cord(3,i)-fxcell(3)*floor((cord(3,i)-celwal(5))*invcel(3))
        gx = gx-xx ; gy = gy-yy ; gz = gz-zz
        mx = int(gx*rr)+1 ; my = int(gy*rr)+1 ; mz = int(gz*rr)+1
        n = nrepcl(mx,my,mz,1) + 1
        nrepcl(mx,my,mz,1) = n
        irepcl(n,mx,my,mz,1) = i
      enddo
      do j = 1,Nrep(2)
        i = replst(j,2)
        gx = cord(1,i)-fxcell(1)*floor((cord(1,i)-celwal(1))*invcel(1))
        gy = cord(2,i)-fxcell(2)*floor((cord(2,i)-celwal(3))*invcel(2))
        gz = cord(3,i)-fxcell(3)*floor((cord(3,i)-celwal(5))*invcel(3))
        gx = gx-xx ; gy = gy-yy ; gz = gz-zz
        mx = int(gx*rr)+1 ; my = int(gy*rr)+1 ; mz = int(gz*rr)+1
        n = nrepcl(mx,my,mz,2) + 1
        nrepcl(mx,my,mz,2) = n
        irepcl(n,mx,my,mz,2) = j
      enddo

      ntbrep(:) = 0
      do mz = 1,nv(1)
      do my = 1,nv(1)
      do mx = 1,nv(1)
        do i = 1,nrepcl(mx,my,mz,2)
          ii = irepcl(i,mx,my,mz,2)
          iat1 = replst(ii,2)
          n = ntbrep(ii)
          do j = 1,nrepcl(mx,my,mz,1)
            jj = irepcl(j,mx,my,mz,1)
            n = n + 1
            itbrep(n,ii) = jj
          enddo
          nce = nrepnr(mx,my,mz)
          do ice = 1,nce
            ix = irepnr(1,ice,mx,my,mz)
            iy = irepnr(2,ice,mx,my,mz)
            iz = irepnr(3,ice,mx,my,mz)
            do j = 1,nrepcl(ix,iy,iz,1)
              iat2 = irepcl(j,ix,iy,iz,1)
              d(1:3) = cord(1:3,iat2) - cord(1:3,iat1)
              d(1:3) = d(1:3) - fxcell(1:3)*nint(d(1:3)-invcel(1:3))
              rr = d(1)*d(1) + d(2)*d(2) + d(3)*d(3)
              if ( rr .le. rlim2 ) then
                n = n + 1
                itbrep(n,ii) = iat2
              endif
            enddo
          enddo
          ntbrep(ii) = n
        enddo
      enddo
      enddo
      enddo

!**********************************************

      return
      end subroutine repprp
