
      subroutine infjob(timnow,user,len)    

      implicit none

      ! Current date and time
        character(23),intent(out):: timnow
      ! User name
        character(*),intent(out):: user
      ! User name character length
        integer(4),intent(out):: len

      character(24):: today

!*****************************************

      call  fdate(today)
      timnow = today(1:20)//"'"//today(23:24)
      call  getlog(user)

      len = index(adjustl(user)," ") - 1

!*****************************

      return
      end subroutine infjob
