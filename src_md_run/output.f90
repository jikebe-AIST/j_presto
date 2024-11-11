
      subroutine output(iread,iprint,ier,onend)
 
!*******************************************************************
!
!     THIS SUBROUTINE IS FOR OUTPUT SOME DATA FILES
!
!*******************************************************************
 
      use COMBAS ; use COMCMM,only: absres
 
      implicit none

      integer(4),intent(in):: iread,iprint
      integer(4),intent(out):: ier
      logical(1),intent(inout):: onend

      integer(4):: iotplf,iocrdf
      character(80):: cotpln,cocrdn
      integer(8):: sec1,sec2

!*******************************************

      ier = 0 
      write(iprint,*)" "
      write(iprint,*)"INFORMATION> OUTPUT (V2.0)"
      write(iprint,*)"             WRITE SOME DATA FILES "
      write(iprint,*)" "
 
!     ** READ INPUT CONTROL DATA **
      call inpout(iread,iprint,ier,onend,iotplf,iocrdf,cotpln,cocrdn)
      if ( ier .ne. 0 ) return
 
!     ** WRITE TOPOLOGY FILE **
      if ( cotpln .ne. " " ) then
        call system_clock(sec1,cr)
        if ( iotplf .eq. 1 ) then
          call outtpl(1,iprint,cotpln,ier)
          call system_clock(sec2)
          call outcpu(iprint,'OUTTPL',sec1,sec2,cr,cm)
        elseif ( iotplf .eq. 2 ) then
          call outtpl(1,iprint,cotpln,ier)
          call system_clock(sec2)
          call outcpu(iprint,'OUTTPB',sec1,sec2,cr,cm)
        endif
        write(iprint,*)""
        if ( ier .ne. 0 ) return
      endif

!     ** WRITE CURRENT STRUCTURE DATA **
      if ( cocrdn .ne. " " ) then
        call system_clock(sec1,cr)
        if ( ixfbou .eq. 1 ) call PBcord
        if ( iocrdf .eq. 1 ) then
          call outpdb(1,iprint,cocrdn,maxatm,ixnatm,cxatmn,cxresn,     &
                      absres,fxmass,fxchrg,10,ixtitl,cxtitl,ier)
          call system_clock(sec2)
          call outcpu(iprint,'OUTPDB',sec1,sec2,cr,cm)
        elseif ( iocrdf .eq. 2 ) then
          call outcrb(1,iprint,cocrdn,ixnatm,ier)
          call system_clock(sec2)
          call outcpu(iprint,'OUTCRB',sec1,sec2,cr,cm)
        endif
        write(iprint,*)""
        if ( ier .ne. 0 ) return
      endif

!*********************************************

      return
      end subroutine output
