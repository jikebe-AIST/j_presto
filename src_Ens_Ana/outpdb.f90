
      subroutine outpdb(icn,tmp)

!*********************************************
!
!     output pdb files
!
!*********************************************

      use COMIFN ; use COMVAL

      implicit none

      integer(4),intent(in):: icn
      character(130),intent(out):: tmp

      integer(4):: i
      character(1):: a1

!***************************************

      write(tmp,'(a,i0,a)')trim(d_PDB),icn,".pdb"
      open(unit=uop,file=tmp,status="replace")
      do i = istATM,ienATM
        if ( .not. outF(i) ) cycle
        if ( RESnum(i) .lt. 10000 ) then
          a1 = adjustl(ATMnm(i))
          write(uop,'(a4,i7,x,a4,x,a4,a1,i4,4x,3f8.3,23x,a1,a2)')      &
          "ATOM",i,ATMnm(i),RESnm(i),CHN(i),RESnum(i),cod(1:3,i),a1,"  "
        elseif ( RESnum(i) .lt. 100000 ) then
          write(uop,'(a4,i7,x,a4,x,a4,a1,i4,4x,3f8.3,23x,a1,a2)')      &
          "ATOM",i,ATMnm(i),RESnm(i),CHN(i),RESnum(i),cod(1:3,i),a1,"  "
        else
          write(uop,'(a4,i7,x,a4,x,a4,a1,i4,4x,3f8.3,23x,a1,a2)')      &
          "ATOM",i,ATMnm(i),RESnm(i),CHN(i),RESnum(i),cod(1:3,i),a1,"  "
        endif
      enddo
      close(uop)

!*******************************************

      return
      end subroutine outpdb
