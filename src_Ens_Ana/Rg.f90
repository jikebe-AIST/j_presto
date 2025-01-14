
      subroutine Rgcalc(icn)

!****************************************************
!
!     Calc. Rg
!
!****************************************************

      use COMIFN ; use COMVAL

      implicit none

      integer(4),intent(in):: icn
      integer(4):: i,iatm
      real(8)::com(3),tcod(3,nRg),tRg,ttRg(nRg),maxd

!****************************************************

      ! calc. geometrical center of mass
      do i = 1,nRg
        iatm = iATMrg(i)
        tcod(1:3,i) = cod(1:3,iatm)
      enddo
      forall(i=1:3) com(i) = sum(tcod(i,1:nRg)) / dble(nRg)

      ! calc. Rg
      forall(i=1:3) tcod(i,1:nRg) = tcod(i,1:nRg) - com(i)
      tcod(1:3,1:nRg) = tcod(1:3,1:nRg)*tcod(1:3,1:nRg)
      forall(i=1:nRg) ttRg(i) = sum(tcod(1:3,i))
      maxd = sqrt(maxval(ttRg))
      tRg = sqrt(sum(ttRg)/dble(nRg))
      write(urg,'(i0,2(x,f8.3),x,f)')icn,tRg,maxd,wfac
      Rg = Rg + tRg*wfac ; Rg2 = Rg2 + tRg*tRg*wfac

!*****************************

      return
      end subroutine Rgcalc


!=========================================================================


      subroutine Rg_init

      use COMIFN ; use COMVAL
      implicit none

      open(unit=urg,file=trim(d_proj)//".Rg",status="replace")
      write(urg,*)"# radius of gyration ("//trim(atom_spec_Rg)//")"
      return
      end subroutine Rg_init


!=================================================================


      subroutine Rg_final

      use COMVAL
      implicit none

      write(urg,*)
      Rg = Rg / sumW2 ; Rg2 = Rg2 / sumW2 ; Rg2 = sqrt(Rg2 - Rg*Rg)
      write(urg,*)"# Average Rg"
      write(urg,*)"# ",Rg," (",Rg2,")"
      close(urg)
      return
      end subroutine Rg_final
