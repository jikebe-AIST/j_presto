
      subroutine inpcls(iprint,iicmpu,cicmpn,ier)

!*************************************************************
!
!     read cluster-parameters and set parameters
!
!*************************************************************

      use COMPAR,only: max14n
      use COMBAS
      use COMERG
      use COMCMM ; use PHYCNS,only: eekcal
      !$ use omp_lib

      implicit none

      ! Logical unit number for output log and read
        integer(4),intent(in):: iprint,iicmpu
      ! Input file name
        character(*),intent(in):: cicmpn
      ! Condition code
        integer(4),intent(out):: ier

      integer(4):: efcol,icol,ncou,ilst(ixnatm)
      real(8):: fac,tmp
      character(80):: line,space,line2,efcol2

      integer(4):: i,j,k,n,ia,ib,isum,pint,n14,curn
      ! heap memory
      integer(4),allocatable:: Tntab14(:,:),Titab14(:,:)
      integer(4):: Tnntab14
      integer(4),allocatable:: Ticmmres(:,:),Ticmmatm(:,:)
      logical(4):: fstflg,ex
      logical(4),allocatable:: resflg(:)
      allocate(Tntab14(2,ixnatm),Titab14(max14n,ixnatm))

!******************************************

!     0) Initital setting
      allocate(icls(ixnatm))
      !$OMP parallel default (none)                                  & !
      !$OMP private(i)                                               & !
      !$OMP shared(icls,ixnatm,Titab14)
      !$OMP do schedule (static)
      do i = 1,ixnatm
        icls(i) = 1
        Titab14(:,i) = 0
      enddo
      !$OMP end do
      !$OMP end parallel
      ier = 0 ; space = " "
      nlev = 3 ; cminsiz = 6.d0
      ncmmatm = 1 ; ncmmres = 0
      allocate(icmmatm(2,1)) ; icmmatm(1:2,1) = (/1,ixnatm/)
      nfrg = 1
      cluster_method = "AND"

      inquire(file=trim(cicmpn),exist=ex)
      if ( ex ) then
        ! 1) Open file
        call flopen(iicmpu,cicmpn,10,'NULL',80,ier)
        if ( ier .ne. 0 ) then
          write(iprint,*)"ERROR> FILE OPEN ERROR IN INPCLS"
          write(iprint,'(a80)')cicmpn
          return
        endif

        ! 2) Read control data
        do
          read(iicmpu,'(a80)',end=800)line
          icol = efcol(line,space,";")
          if ( icol .le. 0 ) cycle
          if ( index(line(1:icol),"INPCLS>") .ne. 0 ) then

            ! CMM LEVEL
            if ( index(line(1:icol),"CMMLVL").ne.0 .or.                &
                 index(line(1:icol),"CELLVL").ne.0 ) then
              line = efcol2(iicmpu,";",ier)
              if ( ier .ne. 0 ) then
                write(iprint,*)"ERROR> INPCLS"
                write(iprint,*)"CMMLVL or CELLVL is inappropriate"
                return
              endif
              read(line,*)nlev
              if ( nlev.lt.3 .or. nlev.gt.5 ) then
                write(iprint,*)"ERROR> INPCLS"
                write(iprint,*)"You can only set CMM level = 3 ~ 5"
                stop
              endif

            ! CELL SIZE
            elseif ( index(line(1:icol),"CMMCSZ").ne.0 .or.            &
                     index(line(1:icol),"CELSIZ").ne.0 ) then
              line = efcol2(iicmpu,";",ier)
              if ( ier .ne. 0 ) then
                write(iprint,*)"ERROR> INPCLS"
                write(iprint,*)"CMMCSZ or CELSIZ is inappropriate"
                return
              endif
              read(line,*)cminsiz
              if ( cminsiz .le. 0.d0 ) cminsiz = 6.d0

            ! RESIDUE CALCULATION LIST FOR CMM
            elseif ( index(line(1:icol),"RESLIST") .ne. 0 ) then
              write(iprint,*)
              write(iprint,*)"INFORMATION> INPCLS"
              write(iprint,*)"        SELECT FOLLOWING ATOMS FOR "//   &
                             "RESIDUE CELL CEPARATION FOR CMM"
              call rdlstatm(iprint,iicmpu,"INPCLS",.true.,ncou,ilst,ier)
              if ( ier .ne. 0 ) then
                write(iprint,*)"ERROR> INPCLS"
                write(iprint,*)"RESLIST is inappropriate"
                return
              endif

              ! Make cmm residue list
              if ( ncou .ne. 0 ) then
                n = absres(ixnatm)
                allocate(resflg(n)) ; resflg(1:n) = .false.
                do i = 1,ncou
                  resflg(absres(ilst(i))) = .true.
                enddo

                allocate(Ticmmres(2,n),Ticmmatm(2,n))
                ncmmatm = 0 ; ncmmres = 0
                ia = 1 ; ib = 1
                do i = 1,n-1
                  if ( resflg(i) .neqv. resflg(i+1) ) then
                    if ( resflg(i) ) then
                      ncmmres = ncmmres + 1
                      Ticmmres(1:2,ncmmres) = (/ia,ib/)
                    else
                      ncmmatm = ncmmatm + 1
                      Ticmmatm(1:2,ncmmatm) = (/ixrstr(ia),ixrend(ib)/)
                    endif
                    ia = i+1
                  endif
                  ib = i+1
                enddo
                ! final loop provision
                if ( resflg(n) ) then
                  ncmmres = ncmmres + 1
                  Ticmmres(1:2,ncmmres) = (/ia,ib/)
                else
                  ncmmatm = ncmmatm + 1
                  Ticmmatm(1:2,ncmmatm) = (/ixrstr(ia),ixrend(ib)/)
                endif

                deallocate(icmmatm)
                allocate(icmmatm(2,ncmmatm),icmmres(2,ncmmres))
                icmmatm(1:2,1:ncmmatm) = Ticmmatm(1:2,1:ncmmatm)
                icmmres(1:2,1:ncmmres) = Ticmmres(1:2,1:ncmmres)
                deallocate(Ticmmatm,Ticmmres)
              endif

            ! CLUSTER LIST 2
            elseif ( index(line(1:icol),"CLSTLIST2") .ne. 0 ) then
              write(iprint,*)
              write(iprint,*)"INFORMATION> INPCLS"
              write(iprint,*)"        SELECT FOLLOWING ATOMS FOR "//   &
                             "CLUSTER2"
              call rdlstatm(iprint,iicmpu,"INPCLS",.true.,ncou,ilst,ier)
              if ( ier .ne. 0 ) then
                write(iprint,*)"ERROR> INPCLS"
                write(iprint,*)"CLSTLIST2 is inappropriate"
                return
              endif
              do i = 1,ncou
                icls(ilst(i)) = 5
              enddo

            ! CLUSTER LIST
            elseif ( index(line(1:icol),"CLSTLIST") .ne. 0 ) then
              write(iprint,*)
              write(iprint,*)"INFORMATION> INPCLS"
              write(iprint,*)"        SELECT FOLLOWING ATOMS FOR "//   &
                             "CLUSTER"
              call rdlstatm(iprint,iicmpu,"INPCLS",.true.,ncou,ilst,ier)
              if ( ier .ne. 0 ) then
                write(iprint,*)"ERROR> INPCLS"
                write(iprint,*)"CLSTLIST is inappropriate"
                return
              endif
              do i = 1,ncou
                icls(ilst(i)) = 4
              enddo

            ! CLUSTER METHOD
            elseif ( index(line(1:icol),"CLSTMETHOD") .ne. 0 ) then
              line = efcol2(iicmpu,";",ier)
              if ( ier .ne. 0 ) then
                write(iprint,*)"ERROR> INPCLS"
                write(iprint,*)"CLSTMETHOD is inappropriate"
                return
              endif
              read(line,*)cluster_method
              if ( cluster_method(1:3).ne."LMD" .and.                  &
                   cluster_method.ne."AND" .and.                       &
                   cluster_method.ne."OR"  .and.                       &
                   cluster_method.ne."XOR" ) then
                write(iprint,*)" !! CAUTION !!"
                write(iprint,*)" CLSTMETHOD you input is strange."
                write(iprint,*)" CLSTMETHOD is changed to 'LMD'"
                cluster_method = "LMD"
              endif

            ! SCALED_TERM
            elseif ( index(line(1:icol),"SCALED_TERM") .ne. 0 ) then
              line = efcol2(iicmpu,";",ier)
              if ( ier .ne. 0 ) then
                write(iprint,*)"ERROR> INPCLS"
                write(iprint,*)"SCALED_TERM is inappropriate"
                return
              endif
              read(line(1:1),*)scale_bond
              read(line(2:2),*)scale_angle
              read(line(3:3),*)scale_tor
              read(line(4:4),*)scale_imp

            ! lambda weight
            elseif ( index(line(1:icol),"LAMBDA_WEIGHT") .ne. 0 ) then
              line = efcol2(iicmpu,";",ier)
              if ( ier .ne. 0 ) then
                write(iprint,*)"ERROR> INPCLS"
                write(iprint,*)"LAMBDA_WEIGHT is inappropriate"
                return
              endif
              read(line,*)tmp
              if ( tmp .gt. 0.d0 ) lambda_m = tmp

            ! lambda temp control
            elseif ( index(line(1:icol),"LAMBDA_TEMP_CONTROL")         &
                     .ne. 0) then
              line = efcol2(iicmpu,";",ier)
              if ( ier .ne. 0 ) then
                write(iprint,*)"ERROR> INPCLS"
                write(iprint,*)"LAMBDA_TEMP_CONTROL is inappropriate"
                return
              endif
              if ( line(1:1).eq."Y" .or. line(1:1).eq."y" )            &
                l_temp_control = .true.
            endif
          endif
        enddo
      endif
800   call flclos(iicmpu,10,ier)

!**********************
!     3) Set variables

      nfrag = maxval(icls(1:ixnatm))
      if ( nfrag .eq. 4 ) then
        nfrag = 2
      elseif ( nfrag .eq. 5 ) then
        nfrag = 3
      endif

      if ( nfrag .ne. 1 ) then
        if ( cluster_method(1:3) .eq. "LMD" ) then
          nfrg = 3
        else
          nfrg = 2
        endif
      endif

      if ( nfrag .eq. 2 ) then
        where ( icls .eq. 1 ) icls = 2
        where ( icls .eq. 4 ) icls = 1
      elseif ( nfrag .eq. 3 ) then
        where ( icls .eq. 1 ) icls = 3
        where ( icls .eq. 4 ) icls = 1
        where ( icls .eq. 5 ) icls = 2
      endif

      ! Monitor the contents
      write(iprint,*)
      write(iprint,*)"INFORMATION> INPCLS"
      write(iprint,*)"     CONTENTS OF CLUSTERING PARAMETERS"
      write(iprint,*)
      write(iprint,*)"     NUMBER OF CLUSTERS TO BE SEPARATED : ",nfrag
      if ( nfrag .ne. 1 ) then
        write(iprint,*)"                  CLUSTER SCALED METHOD : ",   &
                     cluster_method
      endif
      if ( cluster_method(1:3).eq."LMD" .and. nfrag.ne.1 .and.         &
           nfrag.ne.2 ) then
        write(iprint,*)
        write(iprint,*)"INFORMATION> INPCLS"
        write(iprint,*)"     IF YOU USE CLSTMETHOD = 'LMD' or 'LMD2', "
        write(iprint,*)"       YOU CAN'T SET CLSTLIST2 "
        stop
      endif
      
      ! Prepare modified partial charge of atoms &
      ! Set Titab14 (table for atoms with 1-4 relation)
      ccoef = eekcal/fydiel ; fac = sqrt(ccoef)
      allocate(chgmod(ixnatm),chgmod2(ixnatm))
      !$OMP parallel default (none)                                  & !
      !$OMP private(i)                                               & !
      !$OMP shared(chgmod,chgmod2,ixnatm,fac,ccoef,fxchrg,Tntab14)
      !$OMP do schedule (static)
      do i = 1,ixnatm
        chgmod(i) = fac * fxchrg(i)
        chgmod2(i) = ccoef * fxchrg(i)
        Tntab14(:,i) = 0
      enddo
      !$OMP end do
      !$OMP end parallel

      Tnntab14 = 0
      do ia = 1,ixnatm
        isum = sum(ix14if(ia,1:3))
        if ( isum .gt. 0 ) then
          fstflg = .true.
          do ib = 1,isum
            if ( ix14lt(ia,ib) .gt. ia ) then
              if ( fstflg ) then
                Tnntab14 = Tnntab14 + 1
                Tntab14(1,Tnntab14) = ia
                fstflg = .false.
              endif
              pint = Tntab14(2,Tnntab14) + 1
              Titab14(pint,Tnntab14) = ix14lt(ia,ib)
              Tntab14(2,Tnntab14) = pint
            endif
          enddo
        endif
      enddo

      ! Determination of itab14 array size
      pint = maxval(Tntab14(2,1:Tnntab14))
      if ( pint .gt. max14n ) then
        print*,"INFORMATION> SETCMM"
        print*,"?????????????????????????????"
        print*,"? Detected: ntab14 > max14n ?"
        print*,"?????????????????????????????"
        print*,"max14n & the Max. = ",max14n,pint
        print*," Change max14n in COMPAR"
        ier = ier + 1 ; return
      endif

      allocate(ntab14(ixnatm),itab14(2*pint,ixnatm))
      ! Number of 1-4 interactions
      nntab14 = sum(Tntab14(2,1:Tnntab14))

      ! Construct symmetrized 1-4 table
      itab14(:,:) = 0 ; ntab14(:) = 0
      do i = 1,Tnntab14
        ia = Tntab14(1,i) ; n14 = Tntab14(2,i)
        itab14(1:n14,ia) = Titab14(1:n14,i)
        ntab14(ia) = n14
      enddo
      do i = 1,Tnntab14
        ia = Tntab14(1,i) ; n14 = Tntab14(2,i)
        do j = 1,n14
          ib = Titab14(j,i)
          ntab14(ib) = ntab14(ib) + 1
          itab14(ntab14(ib),ib) = ia
        enddo
      enddo
      ! Sort 1-4 table for later use
      do ia = 1,ixnatm
        curn = ntab14(ia)
        ! insertion sort, this should be efficient enough, as the pint ~ 15
        do j= 2,curn
          n = itab14(j,ia)
          do k = j-1,0,-1
            if ( k .eq. 0 ) then
              itab14(1,ia) = n ; exit
            elseif ( itab14(k,ia) .lt. n ) then
              itab14(k+1,ia) = n ; exit
            endif
            itab14(k+1,ia) = itab14(k,ia)
          enddo
        enddo
      enddo

      ! Make nfrag - nfrg relation matrix
      allocate(n_matrix(nfrag,nfrag))
      if ( nfrag .eq. 1 ) then
        n_matrix(1,1) = 1
      elseif ( nfrag .eq. 2 ) then
        if ( cluster_method(1:3) .eq. "LMD" ) then
          n_matrix(1,1) = 1 ; n_matrix(2,2) = 3
          n_matrix(1,2) = 2 ; n_matrix(2,1) = 2
        elseif ( cluster_method .eq. "AND" ) then
          n_matrix(1:2,1:2) = 2 ; n_matrix(1,1) = 1
        elseif ( cluster_method .eq. "OR" ) then
          n_matrix(1:2,1:2) = 1 ; n_matrix(2,2) = 2
        elseif ( cluster_method .eq. "XOR" ) then
          n_matrix(1,1) = 2 ; n_matrix(2,2) = 2
          n_matrix(1,2) = 1 ; n_matrix(2,1) = 1
        endif
      elseif ( nfrag .eq. 3 ) then
        if ( cluster_method .eq. "AND" ) then
          n_matrix(1:3,1:3) = 2
          n_matrix(1,1) = 1 ; n_matrix(2,2) = 1
        elseif ( cluster_method .eq. "OR" ) then
          n_matrix(1:3,1:3) = 2
          n_matrix(1,1) = 1 ; n_matrix(1,2) = 1 ; n_matrix(2,1) = 1
        elseif ( cluster_method .eq. "XOR" ) then
          n_matrix(1:3,1:3) = 2
          n_matrix(1,2) = 1 ; n_matrix(2,1) = 1
        endif
      endif

      ! Monitor scaled_terms
      if ( .not.(nfrg.eq.1 .and. cluster_method.eq."AND") ) then
        write(iprint,*)
        write(iprint,*)"  SCALED_ENERGY_TERMS INFORMATION"
        write(iprint,*)"    * vdW"
        write(iprint,*)"    * Electrostatic"
        if ( scale_bond.eq.1 ) write(iprint,*)"    * Bond"
        if ( scale_angle.eq.1 ) write(iprint,*)"    * Angle"
        if ( scale_tor.ge.1 ) then
          write(iprint,*)"    * Torsion angle"
          if ( scale_tor.eq.2 ) then
            write(iprint,*)"       except for peptide bonds"
          elseif ( scale_tor.eq.3 ) then
            write(iprint,*)"       except for peptide bonds, but "//   &
                 "scale the XXX-PRO peptide bonds"
          endif
        endif
        if ( scale_imp.eq.1 ) write(iprint,*)                          &
          "    * Improper torsion angle"
        write(iprint,*)
      endif

      ! Calc. lambda mass
      i_lambda_m = 0.d0
      if ( cluster_method(1:3) .eq. "LMD" ) then
        if ( lambda_m .eq. 0.d0 ) then
          do i = 1,ixnatm
            if ( icls(i) .eq. 1 ) lambda_m = lambda_m + fxmass(i)
          enddo
        endif
        i_lambda_m = 1.d0 / lambda_m
        write(iprint,*)
        write(iprint,*)"  LAMBDA WEIGHT : ",lambda_m
        write(iprint,*)
      endif

      ! lambda temperature control
      if ( l_temp_control ) then
        write(iprint,*)"  LAMBDA TEMPERATURE CONTROL : ON"
      else
        write(iprint,*)"  LAMBDA TEMPERATURE CONTROL : OFF"
      endif

!**********************

      return
      end subroutine inpcls
