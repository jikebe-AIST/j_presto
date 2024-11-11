
      module COMPAR

!************************************************

      implicit none

      integer(4):: maxatm,maxbnd,maxang,maxtor,maximp,maxshk

      integer(4),parameter:: max14n = 30

      ! For interaction table (emprical coefficients)
      real(8),parameter:: COEn15mxEL = 0.15d0
      real(8),parameter:: COEn15mx = 0.15d0
      real(8),parameter:: COEnvdw = 0.025d0
      real(8),parameter:: COEipmax = 0.25d0

      integer(4):: n15mxEL,n15mx,nvdw,ipmax
      integer(4):: n15mxEL_max = 0
      integer(4):: n15mx_max = 0
      integer(4):: nvdw_max = 0
      integer(4):: ipmax_max = 0

      ! MAX N of atom TYPes for vdW param., FuNCtion kinds, 
      !  and Non-Bond Parameters
      integer(4):: maxtyp
!      integer(4),parameter:: maxtyp = 40
      integer(4),parameter:: maxfnc = 5
      integer(4),parameter:: maxnbp = 4

      ! For SHAKE
      integer(4),parameter:: mxashk = 4
      integer(4),parameter:: maxequ = (mxashk*(mxashk-1))/2

      integer(4),parameter:: maxtrc = 200
      integer(4),parameter:: maxdst = 3000
      integer(4),parameter:: maxene = 15

      end module COMPAR
