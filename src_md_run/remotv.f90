
      subroutine remotv(numatm,mass,cordp,ier,idfstp)
 
!********************************************************************
!
!     REMOVE OUTER-FREEDOM'S VELOCITY
!       REMOVE TRANSLATIONAL VELOCITY OF MASS-CENTER
!       REMOVE ROTATIONAL VELOCITY OF MASS-CENTER
!
!********************************************************************

      use COMBAS,only: ixnatm
      use COMCMMC,only: vel
 
      implicit none

      ! (MAX) number of atoms
        integer(4),intent(in):: numatm
      ! Mass, coordinate & velocity of atoms
        real(8),intent(in):: mass(ixnatm)
      ! past Coordinate
        real(8),intent(in):: cordp(3,ixnatm)
      ! Condition code (0: NO ERROR, -1: MOLECULAR WEIGHT IS ZERO
      !                              -2: MOMENTA MATRIX IS SINGULAR)
        integer(4),intent(out):: ier
      ! Flag for stop center of mass of first chain
      ! (1: STOP TRANSLATIONAL MOVEMENT, 2: STOP ROTATIONAL MOVEMENT
      !  3: STOP TRANSLATIONAL AND ROTATIONAL MOVEMENT)
        integer(4),intent(in):: idfstp

      real(8):: dcord(3,numatm)
!     TOTMAS     R*8               : MOLECULE WEIGHT
!     CORCEN     R*8               : CENTER OF MASS
!     TVELCN     R*8               : TRANSLATIONAL VELOCITY
!     AVELCN     R*8               : ROTATIONAL VELOCITY
!     IXX,IYY,IZZ,IXY,IYZ,IZX      : MOMENTA
!     INVMAT     R*8               : INVERSE MATRIX OF MOMENTA
      real(8):: totmas,corcen(3),tvelcn(3),avelcn(3),ixx,iyy,izz,ixy,  &
                iyz,izx,ixy2,iyz2,izx2,invmat(3,3),detmat,dx,dy,dz
      integer(4):: iatm

!******************************************
 
      ier = 0
 
!     <<<  CALCULATION OF TRANSLATIONAL VELOCITY OF CENTER OF MASS  >>>
      totmas = sum(mass(1:numatm))
      if ( totmas .le. 0.d0 ) then
        ier = -1 ; return
      else
        totmas = 1.d0 / totmas
      endif

      corcen(1:3) = 0.d0 ; tvelcn(1:3) = 0.d0
      if ( idfstp.eq.1 .or. idfstp.eq.3 ) then
        do iatm = 1,numatm
          corcen(1:3) = corcen(1:3) + mass(iatm)*cordp(1:3,iatm)
          tvelcn(1:3) = tvelcn(1:3) + mass(iatm)*vel(1:3,iatm)
        enddo
        corcen(1:3) = corcen(1:3) * totmas
        tvelcn(1:3) = tvelcn(1:3) * totmas
      else
        do iatm = 1,numatm
          corcen(1:3) = corcen(1:3) + mass(iatm)*cordp(1:3,iatm)
        enddo
        corcen(1:3) = corcen(1:3) * totmas
      endif

      do iatm = 1,numatm
        dcord(1:3,iatm) = cordp(1:3,iatm) - corcen(1:3)
      enddo

      if ( idfstp.eq.2 .or. idfstp.eq.3 ) then
 
!       <<<  CALCULATION OF MOMENTA  >>>
        ixx = 0.d0 ; iyy = 0.d0 ; izz = 0.d0
        ixy = 0.d0 ; iyz = 0.d0 ; izx = 0.d0
        do iatm = 1,numatm
          dx = dcord(1,iatm) ; dy = dcord(2,iatm) ; dz = dcord(3,iatm)
          ixx = ixx + mass(iatm) * (dy*dy + dz*dz)
          iyy = iyy + mass(iatm) * (dz*dz + dx*dx)
          izz = izz + mass(iatm) * (dx*dx + dy*dy)
          ixy = ixy - mass(iatm) * dx*dy
          iyz = iyz - mass(iatm) * dy*dz
          izx = izx - mass(iatm) * dz*dx
        enddo
 
!       <<<  CALCULATION OF INVERSE MATRIX OF MOMENTA  >>>
        ixy2 = ixy*ixy ; iyz2 = iyz*iyz ; izx2 = izx*izx
        detmat = ixx*iyy*izz + 2.d0*ixy*iyz*izx                        &
                 -iyz2*ixx - izx2*iyy -ixy2*izz
        if ( detmat .eq. 0.d0 ) then
          ier = -2 ; return
        else
          detmat = 1.d0 / detmat
          invmat(1,1) = (iyy*izz - iyz2) * detmat
          invmat(1,2) = (iyz*izx - izz*ixy) * detmat
          invmat(1,3) = (ixy*iyz - izx*iyy) * detmat
          invmat(2,1) = invmat(1,2)
          invmat(2,2) = (izz*ixx - izx2) * detmat
          invmat(2,3) = (izx*ixy - ixx*iyz) * detmat
          invmat(3,1) = invmat(1,3)
          invmat(3,2) = invmat(2,3)
          invmat(3,3) = (ixx*iyy - ixy2) * detmat
        endif

!       <<<  CALCULATION OF ANGULAR VELOCITY  >>>
        avelcn(1:3) = 0.d0
        do iatm = 1,numatm
          dx = dcord(1,iatm) ; dy = dcord(2,iatm) ; dz = dcord(3,iatm)
          ixx = dy*vel(3,iatm) - dz*vel(2,iatm)
          iyy = dz*vel(1,iatm) - dx*vel(3,iatm)
          izz = dx*vel(2,iatm) - dy*vel(1,iatm)
          avelcn(1:3) = avelcn(1:3) + mass(iatm) *                     &
            (invmat(1:3,1)*ixx + invmat(1:3,2)*iyy + invmat(1:3,3)*izz)
        enddo

      endif
 
!     <<<  REMOVE TRANSLATIONAL AND ROTATIONAL VELOCITY  >>>
      do iatm = 1,numatm
        dx = dcord(1,iatm) ; dy = dcord(2,iatm) ; dz = dcord(3,iatm)
        vel(1,iatm) = vel(1,iatm)-avelcn(2)*dz+avelcn(3)*dy-tvelcn(1)
        vel(2,iatm) = vel(2,iatm)-avelcn(3)*dx+avelcn(1)*dz-tvelcn(2)
        vel(3,iatm) = vel(3,iatm)-avelcn(1)*dy+avelcn(2)*dx-tvelcn(3)
      enddo

!************************************

      return
      end subroutine remotv       
