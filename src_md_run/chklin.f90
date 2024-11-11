
      subroutine chklin(numtor,atmlis,iprint,ier)

      use COMCMMC

      implicit none

      integer(4),intent(in):: numtor
      integer(4),intent(in):: atmlis(4,numtor)
      integer(4),intent(in):: iprint
      integer(4),intent(out):: ier
 
      real(8):: dx21,dy21,dz21,dx32,dy32,dz32,dx43,dy43,dz43
      real(8):: px12,py12,pz12,px23,py23,pz23,r12,r23,r12r23
      integer(4):: itor,iatm(4)
 
!*************************************

      ier = 0
      do itor = 1,numtor
        iatm(1:4) = atmlis(1:4,itor)
        dx21   = cord(1,iatm(2)) - cord(1,iatm(1))
        dy21   = cord(2,iatm(2)) - cord(2,iatm(1))
        dz21   = cord(3,iatm(2)) - cord(3,iatm(1))
        dx32   = cord(1,iatm(3)) - cord(1,iatm(2))
        dy32   = cord(2,iatm(3)) - cord(2,iatm(2))
        dz32   = cord(3,iatm(3)) - cord(3,iatm(2))
        dx43   = cord(1,iatm(4)) - cord(1,iatm(3))
        dy43   = cord(2,iatm(4)) - cord(2,iatm(3))
        dz43   = cord(3,iatm(4)) - cord(3,iatm(3))
        px12   = dy21*dz32 - dy32*dz21
        py12   = dz21*dx32 - dz32*dx21
        pz12   = dx21*dy32 - dx32*dy21
        px23   = dy32*dz43 - dy43*dz32
        py23   = dz32*dx43 - dz43*dx32
        pz23   = dx32*dy43 - dx43*dy32
        r12    = px12*px12 + py12*py12 + pz12*pz12
        r23    = px23*px23 + py23*py23 + pz23*pz23
        R12R23 = r12 * r23
        if ( r12r23 .eq. 0.0D0 ) then
          write(iprint,'(5x,i5,a15,4i7,a)')itor," - TH TORSION (",     &
            iatm(1:4)," ) IS LINEAR THEN CAN NOT CALCULATE TORSIONAL"//&
            " ANGLE"
          ier = -1
        endif
      enddo

!************************************

      return
      end subroutine chklin
