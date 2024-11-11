
      subroutine caldih(dx1,dy1,dz1,dx2,dy2,dz2,dx3,dy3,dz3,dx4,dy4,   &
                        dz4,pi,phi)

!*******************************************************************
!
!          CALCULATE DIHEDRAL ANGLE
!
!*******************************************************************

      ! (x,y,z)-coordinates of atom(1-4)
      real(8),intent(in):: dx1,dy1,dz1,dx2,dy2,dz2,dx3,dy3,dz3,dx4,dy4,&
                           dz4,pi
      ! dihedral angle (degree)
      real(8),intent(out):: phi
     
      real(8):: dx21,dy21,dz21,dx32,dy32,dz32,dx43,dy43,dz43,px12,py12,&
                pz12,px23,py23,pz23,r12,r23,r12r23,s1223,cosp,px123,   &
                py123,pz123,s32123
 
!******************************************

      dx21 = dx2 - dx1 ; dy21 = dy2 - dy1 ; dz21 = dz2 - dz1
      dx32 = dx3 - dx2 ; dy32 = dy3 - dy2 ; dz32 = dz3 - dz2
      dx43 = dx4 - dx3 ; dy43 = dy4 - dy3 ; dz43 = dz4 - dz3

      px12 = dy21*dz32 - dy32*dz21
      py12 = dz21*dx32 - dz32*dx21
      pz12 = dx21*dy32 - dx32*dy21
      px23 = dy32*dz43 - dy43*dz32
      py23 = dz32*dx43 - dz43*dx32
      pz23 = dx32*dy43 - dx43*dy32
 
      r12  = px12*px12 + py12*py12 + pz12*pz12
      r23  = px23*px23 + py23*py23 + pz23*pz23
      r12r23    = sqrt( r12 * r23 )

!     IF BOND ANGLE = 0 OR PI , TORSIONAL ANGLE CAN'T BE CALCULATED
      if ( r12r23 .gt. 0.0D0 ) then
        r12r23 = 1.0D0 / r12r23
        s1223 = px12*px23 + py12*py23 + pz12*pz23
        cosp = s1223 * r12r23
        if ( cosp .gt. 1.d0 ) cosp = 1.d0
        if ( cosp .lt. -1.d0 ) cosp = -1.d0
 
        px123     = py12*pz23 - py23*pz12
        py123     = pz12*px23 - pz23*px12
        pz123     = px12*py23 - px23*py12
        s32123    = dx32*px123 + dy32*py123 + dz32*pz123
        phi = acos(cosp)
        if ( s32123 .lt. 0.d0 ) phi = -phi
      else
        phi = 0.d0
      endif

      phi = phi * ( 180.d0 / pi )

!*******************************************

      return
      end subroutine caldih
