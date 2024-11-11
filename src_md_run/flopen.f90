
      subroutine flopen(iounit,filena,idacc,cblank,lrec,ier)

!*******************************************************************
!
!     THIS SUBROUTINE IS FOR OPEN FILE FOR UNIX
!
!*******************************************************************

      implicit none

      ! Logical unit number
        integer(4),intent(in):: iounit
      ! File name
        character(80),intent(in):: filena
      ! ID number for access
      !   10: formatted sequential readonly keep
      !   11: formatted sequential writeonly keep
      !   12: formatted sequential read/write keep
      !   13: formatted sequential read/write delete
      !   20: unformatted sequential readonly keep
      !   21: unformatted sequential writeonly keep
      !   22: unformatted sequential read/write keep
      !   23: unformatted sequential read/write delete
      !   30: unformatted direct readonly keep
      !   31: unformatted direct writeonly keep
      !   32: unformatted direct read/write keep
      !   33: unformatted direct read/write delete
        integer(4),intent(in):: idacc
      ! 'NULL' regard blank as NULL, 'ZERO' regard blank as zero
      !   In case of 10, 12, 13, this option is available
        character(4),intent(in):: cblank
      ! Record length (In case of 30 - 33, this option is available)
        integer(4),intent(in):: lrec
      ! Condition code (0: NO ERROR, -100: OPEN ERROR, -101: PARAMETER ERROR)
        integer(4),intent(out):: ier

      character(80):: space = ' '

!*****************************************************

      ! * CHECK PARAMETER *
      if ( iounit.lt.0 .or. iounit.ge.100 ) then
        ier = -101 ; return
      else
        ier = 0
      endif

      if ( cblank.ne.'NULL' .and. cblank.ne.'ZERO' ) then
        print*,"ahoaho " ; stop
      endif

      ! * OPEN FILE IN UNIX *
      select case (idacc)
        case (10)
          if ( filena .eq. space ) then
            open(unit=iounit,form='formatted',access='sequential',     &
                 status='old',blank=cblank,err=999)
          else
            open(unit=iounit,form='formatted',access='sequential',     &
                 file=filena,status='old',blank=cblank,err=999)
          endif
        case (11)
          if ( filena .eq. space ) then
            open(unit=iounit,form='formatted',access='sequential',     &
                 status='new',err=999)
          else
            open(unit=iounit,form='formatted',access='sequential',     &
                 file=filena,status='new',err=999)
          endif
        case (12)
          if ( filena .eq. space ) then
            open(unit=iounit,form='formatted',access='sequential',     &
                 status='unknown',blank=cblank,err=999)
          else
            open(unit=iounit,form='formatted',access='sequential',     &
                 file=filena,status='unknown',blank=cblank,err=999)
          endif
        case (13)
          open(unit=iounit,form='formatted',access='sequential',       &
               status='scratch',blank=cblank,err=999)
        case (20)
          if ( filena .eq. space ) then
            open(unit=iounit,form='unformatted',access='sequential',   &
                 status='old',err=999)
          else
            open(unit=iounit,form='unformatted',access='sequential',   &
                 file=filena,status='old',err=999)
          endif
        case (21)
          if ( filena .eq. space ) then
            open(unit=iounit,form='unformatted',access='sequential',   &
                 status='new',err=999)
          else
            open(unit=iounit,form='unformatted',access='sequential',   &
                 file=filena,status='new',err=999)
          endif
        case (22)
          if ( filena .eq. space ) then
            open(unit=iounit,form='unformatted',access='sequential',   &
                 status='unknown',err=999)
          else
            open(unit=iounit,form='unformatted',access='sequential',   &
                 file=filena,status='unknown',err=999)
          endif
        case (23)
          open(unit=iounit,form='unformatted',access='sequential',     &
               status='scratch',err=999)
        case (30:33)
          if ( lrec .le. 0 ) then
            ier = -101 ; return
          endif
          select case (idacc)
            case (30)
              if ( filena .eq. space ) then
                open(unit=iounit,form='unformatted',access='direct',   &
                     status='old',recl=lrec,err=999)
              else
                open(unit=iounit,form='unformatted',access='direct',   &
                     file=filena,status='old',recl=lrec,err=999)
              endif
            case (31)
              if ( filena .eq. space ) then
                open(unit=iounit,form='unformatted',access='direct',   &
                     status='new',recl=lrec,err=999)
              else
                open(unit=iounit,form='unformatted',access='direct',   &
                     file=filena,status='new',recl=lrec,err=999)
              endif
            case (32)
              if ( filena .eq. space ) then
                open(unit=iounit,form='unformatted',access='direct',   &
                     status='unknown',recl=lrec,err=999)
              else
                open(unit=iounit,form='unformatted',access='direct',   &
                     file=filena,status='unknown',recl=lrec,err=999)
              endif
            case (33)
              open(unit=iounit,form='unformatted',access='direct',   &
                   status='scratch',recl=lrec,err=999)
          end select
        case default
          ier = -101
      end select

!*********************************************

      return

999   ier = -100 ; return
      end subroutine flopen
