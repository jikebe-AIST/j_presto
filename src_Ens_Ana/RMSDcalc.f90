
      subroutine RMSDcalc(nATM,nRMSD,nRMSD2,Rcod,Rcod2,cod,iATMrmsd,   &
                          iATMrmsd2,RMSD,flag)

!*************************************************
!
!     Calculate RMSD between two conformations
!
!*************************************************

      implicit none

      ! Number of atoms
        integer(4),intent(in):: nATM
      ! Numer of atoms for RMSD calc.
        integer(4),intent(in):: nRMSD,nRMSD2
      ! (Reference) COorDinates
        real(4),intent(in):: Rcod(3,nRMSD),Rcod2(3,nRMSD2)
        real(4),intent(inout):: cod(3,nATM)
      ! array to keep AToM number for RMSD calc.
        integer(4),intent(in):: iATMrmsd(nRMSD),iATMrmsd2(nRMSD2)
      ! RMSD
        real(8),intent(out):: RMSD
      ! Flag for output PDB
        logical(4),intent(in):: flag

      ! Temporary (Reference) COorDinates only for RMSD calc.
        real(4):: tRcod(3,nRMSD),tcod(3,nRMSD),tcod2(3,nRMSD2),        &
                  tRcod2(3,nRMSD2)
      ! (Reference) coordinate of Center Of Mass
        real(4):: cof(3),Rcof(3),cof2(3)
      ! Rotation Matrix
        real(8):: RM(3,3)

      integer(4):: i
      real(8):: ttcod(3)
     
!*******************
!***This part was changed by IKEBE

      ! Selection atoms used for RMSD calc.
      tcod(1:3,1:nRMSD) = cod(1:3,iATMrmsd(1:nRMSD))
      tcod2(1:3,1:nRMSD2) = cod(1:3,iATMrmsd2(1:nRMSD2))
      tRcod(:,:) = Rcod(:,:) ; tRcod2 = Rcod2

!***

      ! Shift conformations to the origin
      call shift(nRMSD,tRcod,Rcof) !(internal subroutine)
      call shift(nRMSD,tcod,cof) !(internal subroutine)
      ! Rotate the molecule & calc. RMSD
      call rot(nRMSD,tcod,tRcod,RM)

      ! Shift conformation
      forall(i=1:3) tcod2(i,1:nRMSD2) = tcod2(i,1:nRMSD2) - cof(i)
      forall(i=1:3) tRcod2(i,1:nRMSD2) = tRcod2(i,1:nRMSD2) - Rcof(i)
      ! Rotate
      do i = 1,nRMSD2
        ttcod(1:3) = tcod2(1:3,i)
        tcod2(1,i) = RM(1,1)*ttcod(1)+RM(1,2)*ttcod(2)+RM(1,3)*ttcod(3)
        tcod2(2,i) = RM(2,1)*ttcod(1)+RM(2,2)*ttcod(2)+RM(2,3)*ttcod(3)
        tcod2(3,i) = RM(3,1)*ttcod(1)+RM(3,2)*ttcod(2)+RM(3,3)*ttcod(3)
      enddo
      ! calc. RMSD
      RMSD = 0.d0
      do i = 1,nRMSD2
        ttcod(1:3) = tcod2(1:3,i)-tRcod2(1:3,i)
        RMSD = dot_product(ttcod,ttcod) + RMSD
      enddo
      RMSD = sqrt( RMSD / dble(nRMSD2) )

      ! Make superimposed coordinates for Output PDB
      if ( flag ) then
        ! translation to (0,0,0)
        forall(i=1:3) cod(i,1:nATM) = cod(i,1:nATM) - cof(i)

        ! Rotate
        do i = 1,nATM
          ttcod(1:3) = cod(1:3,i)
          cod(1,i) = RM(1,1)*ttcod(1)+RM(1,2)*ttcod(2)+RM(1,3)*ttcod(3)
          cod(2,i) = RM(2,1)*ttcod(1)+RM(2,2)*ttcod(2)+RM(2,3)*ttcod(3)
          cod(3,i) = RM(3,1)*ttcod(1)+RM(3,2)*ttcod(2)+RM(3,3)*ttcod(3)
        enddo

        ! translation to center of mass of ref. conf.
        do i = 1,nATM
          cod(1:3,i) = cod(1:3,i) + Rcof(1:3)
        enddo
      endif

!*******************

      return
      end subroutine RMSDcalc


!=======================================================


      subroutine shift(nATM,cod,cof)

!*************************************
!
!     SHIFT mass center to the origin
!
!*************************************

      implicit none

      ! Number of atoms
        integer(4),intent(in):: nATM
      ! COorDinate
        real(4),intent(inout):: cod(3,nATM)
      ! coodinate of Center Of Mass
        real(4),intent(out):: cof(3)

      integer(4):: i

!*****************************

      ! Calc. of center of mass
      forall(i=1:3) cof(i) = sum(cod(i,1:nATM)) / dble(nATM)
      
      ! Shift conformation
      forall(i=1:3) cod(i,1:nATM) = cod(i,1:nATM) - cof(i)

!*****************************

      return
      end subroutine shift


!======================================================


      subroutine rot(nATM,cod,Rcod,trns)

!*******************************************
!
!     ROTate the molecule
!
!*******************************************

      implicit none

      ! Number of atoms
        integer(4),intent(in):: nATM
      ! (Reference) COorDinates
        real(4),intent(in):: Rcod(3,nATM)
        real(4),intent(inout):: cod(3,nATM)
!      ! RMSD
!        real(8),intent(out):: RMSD
      ! Rotation Matrix
        real(8),intent(out):: trns(3,3)

      integer(4):: i,j,m,n,k1,k2,k3
      real(8):: w(3,3),s(6),ev(3),uv(3,3),work(3,3),tr(3,3),tcod(3)
      real(8):: trs,wsq(3),dett,det

!****************************
      
      ! calc. w matrix
      w(:,:) = 0.d0
      do n = 1,3
      do m = 1,3
        do i = 1,nATM
          w(m,n) = w(m,n) + Rcod(m,i)*cod(n,i)
        enddo
      enddo
      enddo
      w(:,:) = w(:,:) / dble(nATM)

      ! symmetrize s = w'*w
      j = 0
      do m = 1,3
      do n = 1,m
        j = j + 1
        s(j) = w(1,m)*w(1,n) + w(2,m)*w(2,n) + w(3,m)*w(3,n)
      enddo
      enddo

      ! diagonalize s
      call diag(s,ev,uv) !(internal subroutine)
      wsq(:) = sqrt(ev(:))
      
      ! 2**3 possible combination for (w'*w)**(-1/2)
      do k1 = 1,-1,-2
      do k2 = 1,-1,-2
      do k3 = 1,-1,-2

        do n = 1,3
        do m = 1,3
          work(m,n) = uv(m,1) * uv(n,1) / wsq(1) * k1 +              &
                      uv(m,2) * uv(n,2) / wsq(2) * k2 +              &
                      uv(m,3) * uv(n,3) / wsq(3) * k3
        enddo
        enddo

        ! transformation matrix: tr = w*(w'w)**(-1/2)
        do n = 1,3
        do m = 1,3
          tr(m,n) = w(m,1)*work(1,n) + w(m,2)*work(2,n) +            &
                    w(m,3)*work(3,n)
        enddo
        enddo

        ! trace(tr*w')=trace((w'*w)**(1/2))
        trs = wsq(1)*k1 + wsq(2)*k2 + wsq(3)*k3

        ! get matrix: -tr*w = -w*(w'*w)**(1/2)*w'
        j = 0
        do m = 1,3
          do n = 1,m
            j = j + 1
            s(j) = -tr(1,m)*w(1,n) - tr(2,m)*w(2,n) - tr(3,m)*w(3,n)
          enddo
          s(j) = s(j) + trs
        enddo

        ! get curvatures:
        ! diagonalize trace((w'*w)**(1/2))-w*(w'*w)**(1/2)*w'
        call diag(s,ev,work)
        dett = det(tr)

        if ( ev(1).ge.0.d0 .and. ev(2).ge.0.d0 .and. ev(3).ge.0.d0   &
             .and. dett.ge.0.d0 ) then
          trns(1:3,1:3) = tr(1:3,1:3) ; goto 100
        endif

      enddo
      enddo
      enddo

100   continue

!******************

!      ! rotate the molecule
!      do i = 1,nATM
!        tcod(1:3) = cod(1:3,i)
!        cod(1,i) = trns(1,1)*tcod(1)+trns(1,2)*tcod(2)+trns(1,3)*tcod(3)
!        cod(2,i) = trns(2,1)*tcod(1)+trns(2,2)*tcod(2)+trns(2,3)*tcod(3)
!        cod(3,i) = trns(3,1)*tcod(1)+trns(3,2)*tcod(2)+trns(3,3)*tcod(3)
!      enddo

!      ! calc. RMSD
!      RMSD = 0.d0
!      do i = 1,nATM
!        RMSD = (cod(1,i)-Rcod(1,i))*(cod(1,i)-Rcod(1,i)) +          &
!               (cod(2,i)-Rcod(2,i))*(cod(2,i)-Rcod(2,i)) +          &
!               (cod(3,i)-Rcod(3,i))*(cod(3,i)-Rcod(3,i)) + RMSD
!      enddo
!      RMSD = sqrt( RMSD / dble(nATM) )
      
!****************************

      return
      end subroutine rot


!=======================================================


      subroutine diag(s,e,v)

      implicit none

      integer(4),parameter:: worklen = 3*3-1
      real(8):: s(6),e(3),v(3,3),wk(3),work(worklen)
      integer(4):: i,j,k,ind

!************************

      k = 0 ; ind = -1
      do i = 1,3
      do j = 1,i
        k = k + 1
        v(i,j) = s(k)
        v(j,i) = s(k)
      enddo
      enddo

      ! LAPACK call
      call dsyev("Vector","Up",3,v,3,e,work,worklen,ind)
      ! LAPACK returns eigenvalues in an ascending order
      ! Reorder to a descending order
      work(1:3) = v(1:3,3) ; v(1:3,3) = v(1:3,1)
      v(1:3,1) = work(1:3)
      work(1) = e(3) ; e(3) = e(1) ; e(1) = work(1)

      ! or NUMPAC version
!      call hoqrvd(v,3,3,e,wk,1d-7,ind)
      if (ind.ne.0) write(6,*)"####### error ind = ",ind

!*************************

      return
      end subroutine diag


!=======================================================


      function det(a)

      implicit none
      real(8):: a(3,3),det

      det = a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1) +         &
            a(1,3)*a(3,2)*a(2,1) - a(1,2)*a(2,1)*a(3,3) -         &
            a(1,3)*a(2,2)*a(3,1) - a(1,1)*a(2,3)*a(3,2)

      return
      end function det
