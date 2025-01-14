
      subroutine allocate_arrays

!***********************************************************************

      use COMTPL

      implicit none

!***********************************************************************

      ! For COMBAS
      allocate(ixamol(maxatm),ixachn(maxatm),ixares(maxatm),           &
        ixatyp(maxatm),ix14if(maxatm,3),ix14lt(maxatm,max14n),         &
        ixincd(maxatm,5),fxchrg(maxatm),fxmass(maxatm),fxvdwr(maxatm), &
        fxincd(maxatm,3),cxatmn(maxatm),cxresn(maxatm),cxatmt(maxatm))

      ! For COMERG
      allocate(iypbnd(2,maxbnd),fyfbnd(maxbnd),fyqbnd(maxbnd))
      allocate(iypang(3,maxang),fyfang(maxang),fyqang(maxang))
      allocate(iyptor(4,maxtor),iytdiv(maxtor),iytnbf(maxtor),         &
        fyftor(maxtor),fytrot(maxtor),fytphs(maxtor),fytvws(maxtor),   &
        fytess(maxtor))
      allocate(iypimp(4,maximp),iyidiv(maximp),iyinbf(maximp),         &
        fyfimp(maximp),fyirot(maximp),fyiphs(maximp),fyivws(maximp),   &
        fyiess(maximp))
      allocate(iyppid(maxtyp,maxtyp),fynbpp(maxnbp,maxtyp,maxtyp),     &
        fyvwrd(maxtyp),fyvwme(maxtyp),fy14sv(maxtyp),fy14se(maxtyp))

!***********************************************************************

      return
      end subroutine allocate_arrays
