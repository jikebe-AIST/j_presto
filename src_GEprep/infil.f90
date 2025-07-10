
      subroutine infil(iu,filnam,iST,iEN,fw,onend,ier)

!****************************************

      implicit none

      integer(4),intent(in):: iu
      character(9999),intent(out):: filnam
      integer(4),intent(out):: iST,iEN,ier
      real(8),intent(out):: fw
      logical(1),intent(out):: onend

      integer(4):: i,ic,iend
      character(9999):: line

!**********************************

      ic = 0 ; fw = 1.d0 ; onend = .true. ; ier = 0
      do
        read(iu,'(a)',end=800)line
        i = index(line,";")
        if ( i .ne. 0 ) line = line(1:i-1)
        if ( line .eq. " " ) cycle

        ic = ic + 1
        select case (ic)
          case (1)
            filnam = line
          case (2)
            i = index(line,"->")
            if ( i .ne. 0 ) then
              onend = .false. ; line = line(1:i-1)
            endif
            line = adjustl(line)
            iend = index(line," ") - 1
            read(line(1:iend),*)iST
            line = adjustl(line(iend+1:))
            iend = index(line," ") - 1
            read(line(1:iend),*)iEN
            line = adjustl(line(iend+1:))
            if ( line .ne. " " ) read(line,*)fw
            exit
        end select
      enddo

      if ( iST .le. 0 ) iST = 1
      if ( fw .lt. 0.d0 ) fw = 1.d0
      return

800   ier = -1

!************************************

      return
      end subroutine infil
