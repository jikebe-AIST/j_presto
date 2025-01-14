
      subroutine RMSF_init

      use COMIFN ; use COMVAL

      implicit none

!****************************************************

      allocate(acod(1:3,istATM:ienATM),a2cod(istATM:ienATM))
      acod(:,:) = 0.d0 ; a2cod(:) = 0.d0

!****************************************************

      return
      end subroutine RMSF_init


!==============================================================


      subroutine RMSF_calc

      use COMIFN ; use COMVAL

      implicit none

      integer(4):: i

!******************************************************

      forall(i=istATM:ienATM) acod(:,i) = acod(:,i) + cod(:,i)*wfac
      do i = istATM,ienATM
        a2cod(i) = a2cod(i) + wfac*(cod(1,i)*cod(1,i) +                &
          cod(2,i)*cod(2,i) + cod(3,i)*cod(3,i))
      enddo

!******************************************************

      return
      end subroutine RMSF_calc


!===============================================================


      subroutine RMSF_final

      use COMIFN ; use COMVAL

      implicit none

      integer(4):: i
      real(8):: r

!*******************************************************

      r = 1.d0 / sumW2
      acod(:,:) = acod(:,:) * r
      a2cod(:) = a2cod(:) * r

      open(unit=1,file=trim(d_proj)//"_RMSF.pdb",status="replace")
      open(unit=2,file=trim(d_proj)//"_ave.pdb",status="replace")
      do i = istATM,ienATM
        if ( .not. outF(i) ) cycle
        r = a2cod(i) - (acod(1,i)*acod(1,i) + acod(2,i)*acod(2,i) +    &
                        acod(3,i)*acod(3,i))
        r = sqrt(r)
        if ( RESnum(i) .lt. 10000 ) then
          write(1,'(a4,i7,x,a4,x,a4,a1,i4,4x,3f8.3,2f6.3)')      &
            "ATOM",i,ATMnm(i),RESnm(i),CHN(i),RESnum(i),    &
            Rcod(1:3,i),0.0,r
          write(2,'(a4,i7,x,a4,x,a4,a1,i4,4x,3f8.3,2f6.3)')      &
            "ATOM",i,ATMnm(i),RESnm(i),CHN(i),RESnum(i),    &
            acod(1:3,i),0.0,r
        elseif ( RESnum(i) .lt. 100000 ) then
          write(1,'(a4,i7,x,a4,x,a4,a1,i5,3x,3f8.3,2f6.3)')      &
            "ATOM",i,ATMnm(i),RESnm(i),CHN(i),RESnum(i),    &
            Rcod(1:3,i),0.0,r
          write(2,'(a4,i7,x,a4,x,a4,a1,i5,3x,3f8.3,2f6.3)')      &
            "ATOM",i,ATMnm(i),RESnm(i),CHN(i),RESnum(i),    &
            acod(1:3,i),0.0,r
        else
          write(1,'(a4,i7,x,a4,x,a4,a1,i6,2x,3f8.3,2f6.3)')      &
            "ATOM",i,ATMnm(i),RESnm(i),CHN(i),RESnum(i),    &
            Rcod(1:3,i),0.0,r
          write(2,'(a4,i7,x,a4,x,a4,a1,i6,2x,3f8.3,2f6.3)')      &
            "ATOM",i,ATMnm(i),RESnm(i),CHN(i),RESnum(i),    &
            acod(1:3,i),0.0,r
        endif
      enddo
      close(1) ; close(2)

!*******************************************************

      return
      end subroutine RMSF_final
