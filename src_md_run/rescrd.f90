
      subroutine rescrd(numatm,numvar,onfree,cord)

!*******************************************************************
!
!     RESTORE COORDINATE
!         IN BINARY COORDINATE DATA , THERE ARE ONLY FREE ATOMS
!         INFORMATION
!         THIS SUBROUTINE MURGES FREE AND FIX ATOM INFORMATION
!
!*******************************************************************

      use COMPAR,only: maxatm ; use COMMIS

      implicit none

      integer(4),intent(in):: numatm,numvar,onfree(maxatm)
      real(8),intent(out):: cord(maxatm,3)

      integer(4):: numpnt,i 
 
!************************************************

      numpnt = numvar + 1 
      do i = numatm,1,-1
        if ( onfree(i) .eq. 1 ) then
          numpnt = numpnt - 1
          cord(i,1:3) = cord(numpnt,1:3)
        else
          cord(i,1:3) = fucord(1:3,i)
        endif
      enddo

!**************************************

      return
      end subroutine rescrd
