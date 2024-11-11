
      subroutine PBcord

      use COMBAS ; use COMCMMC

      implicit none

      integer(4):: i,j,ist,ien
      real(8):: tcod(3),dcod(3)

!****************************************************

      do i = 1,nchain
        ist = iichain(1,i) ; ien = iichain(2,i)
        tcod(1) = sum(cord(1,ist:ien))
        tcod(2) = sum(cord(2,ist:ien))
        tcod(3) = sum(cord(3,ist:ien))
        tcod(1:3) = tcod(1:3) / dble(ien-ist+1)
        dcod(1) = - fxcell(1) * floor((tcod(1)-celwal(1))*invcel(1))
        dcod(2) = - fxcell(2) * floor((tcod(2)-celwal(3))*invcel(2))
        dcod(3) = - fxcell(3) * floor((tcod(3)-celwal(5))*invcel(3))
        do j = ist,ien
          cord(1:3,j) = cord(1:3,j) + dcod(1:3)
        enddo
      enddo

!****************************************************

      return
      end subroutine PBcord
