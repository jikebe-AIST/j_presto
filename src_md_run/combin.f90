
      function combin(n,r)

!*******************************************************************
!
!     THIS FUNCTION IS FOR CALCULATION COMBINATION
!
!*******************************************************************

      implicit none

      integer(4),intent(in):: n,r
      integer(4):: combin,nfact,nrfact,rfact,factor 
 
!*****************************************

      nfact = factor(n)
      nrfact = factor(n-r)
      rfact = factor(r)
      combin = nfact / (nrfact*rfact) 

!*****************************************

      return
      end function combin
