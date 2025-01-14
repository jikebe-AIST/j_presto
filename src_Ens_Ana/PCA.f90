
      subroutine PCA_mkcod

!************************************************
!
!     Make inter-COoDinate for PCA analysis
!
!************************************************

      use COMIFN ; use COMVAL

      implicit none

      integer(4):: i,j,k,l,icn,ires1,ires2,iatm1,iatm2,n
      real(4):: dcod(3),dist,dsum,rad1,rad2,thre
      real(4),allocatable:: PCAcod(:)

!****************************

      if ( trim(PCA_method) .eq. "cord" ) then
        n = nPCA * 3
      elseif ( trim(PCA_method) .eq. "intercord" ) then
        n = (nPCA*(nPCA-1))/2
      else
        n = Ncon_res_PCA
      endif
      allocate(PCAcod(n))

      if ( trim(PCA_method) .eq. "contact" ) then
        OUTER : do i = 1,Ncon_res_PCA
          ires1 = con_res_PCA(1,i) ; ires2 = con_res_PCA(2,i)
          do j = iresPCA(ires1),iresPCA(ires1+1)-1
            iatm1 = iATMPCA(j) ; rad1 = radPCA(j)
            do k = iresPCA(ires2),iresPCA(ires2+1)-1
              iatm2 = iATMPCA(k) ; rad2 = radPCA(k)
              dcod(1:3) = cod(1:3,iatm1) - cod(1:3,iatm2)
              dcod(1:3) = dcod(1:3) - bsize(1:3) *                     &
                                      nint(dcod(1:3)*ibsize(1:3))
              dcod(1:3) = dcod(1:3)*dcod(1:3)
              thre = rad1 + rad2 + tolePCA ; thre = thre*thre
              if ( sum(dcod) .le. thre ) then
                PCAcod(i) = 1.0 ; cycle OUTER
              endif
            enddo
          enddo
          PCAcod(i) = 0.0
        enddo OUTER

      elseif ( trim(PCA_method) .eq. "ave" ) then
        do i = 1,Ncon_res_PCA
          ires1 = con_res_PCA(1,i) ; ires2 = con_res_PCA(2,i)
          dsum = 0.0
          do j = iresPCA(ires1),iresPCA(ires1+1)-1
            iatm1 = iATMPCA(j)
            do k = iresPCA(ires2),iresPCA(ires2+1)-1
              iatm2 = iATMPCA(k)
              dcod(1:3) = cod(1:3,iatm1) - cod(1:3,iatm2)
              dcod(1:3) = dcod(1:3) - bsize(1:3) *                     &
                                      nint(dcod(1:3)*ibsize(1:3))
              dcod(1:3) = dcod(1:3)*dcod(1:3)
              dsum = dsum + sqrt(sum(dcod))
            enddo
          enddo
          dsum = dsum / real((iresPCA(ires1+1)-iresPCA(ires1))*        &
                                 (iresPCA(ires2+1)-iresPCA(ires2)))
          PCAcod(i) = dsum
        enddo

      elseif ( trim(PCA_method) .eq. "max" ) then
        do i = 1,Ncon_res_PCA
          ires1 = con_res_PCA(1,i) ; ires2 = con_res_PCA(2,i)
          dist = 0.0
          do j = iresPCA(ires1),iresPCA(ires1+1)-1
            iatm1 = iATMPCA(j)
            do k = iresPCA(ires2),iresPCA(ires2+1)-1
              iatm2 = iATMPCA(k)
              dcod(1:3) = cod(1:3,iatm1) - cod(1:3,iatm2)
              dcod(1:3) = dcod(1:3) - bsize(1:3) *                     &
                                      nint(dcod(1:3)*ibsize(1:3))
              dcod(1:3) = dcod(1:3)*dcod(1:3)
              dist = max(dist,sum(dcod))
            enddo
          enddo
          PCAcod(i) = sqrt(dist)
        enddo

      elseif ( trim(PCA_method) .eq. "min" ) then
        do i = 1,Ncon_res_PCA
          ires1 = con_res_PCA(1,i) ; ires2 = con_res_PCA(2,i)
          dist = 99999999.0
          do j = iresPCA(ires1),iresPCA(ires1+1)-1
            iatm1 = iATMPCA(j)
            do k = iresPCA(ires2),iresPCA(ires2+1)-1
              iatm2 = iATMPCA(k)
              dcod(1:3) = cod(1:3,iatm1) - cod(1:3,iatm2)
              dcod(1:3) = dcod(1:3) - bsize(1:3) *                     &
                                      nint(dcod(1:3)*ibsize(1:3))
              dcod(1:3) = dcod(1:3)*dcod(1:3)
              dist = min(dist,sum(dcod))
            enddo
          enddo
          PCAcod(i) = sqrt(dist)
        enddo

      elseif ( trim(PCA_method) .eq. "cord" ) then
        do i = 1,nPCA
          PCAcod((i-1)*3+1:i*3) = cod(1:3,iATMPCA(i))
        enddo
      
      elseif ( trim(PCA_method) .eq. "intercord" ) then
        k = 0
        do i = 1,nPCA-1
          iatm1 = iATMPCA(i)
          do j = i+1,nPCA
            iatm2 = iATMPCA(j)
            k = k + 1
            dcod(1:3) = cod(1:3,iatm1) - cod(1:3,iatm2)
            dcod(1:3) = dcod(1:3) - bsize(1:3) *                       &
                                    nint(dcod(1:3)*ibsize(1:3))
            dcod(1:3) = dcod(1:3)*dcod(1:3)
            PCAcod(k) = sqrt(sum(dcod(1:3)))
          enddo
        enddo
      endif

      write(upc)wfac,PCAcod(1:n)
      deallocate(PCAcod)

!****************************

      return
      end subroutine PCA_mkcod


!==================================================================


      subroutine PCAprep()

      use COMVAL ; use COMIFN

      integer(4):: i,j,k,ii,jj
      character(1):: ichn,jchn

!**********************************************

      if ( trim(PCA_method).ne."cord" .and.                            &
           trim(PCA_method).ne."intercord" ) then
        ! Making residue list
        NresPCA = 1 ; j = RESnum(iATMPCA(1))
        do i = 2,nPCA
          if ( RESnum(iATMPCA(i)) .ne. j ) then
            NresPCA = NresPCA + 1 ; j = RESnum(iATMPCA(i))
          endif
        enddo
        allocate(iresPCA(NresPCA+1))
        k = 1 ; iresPCA(1) = 1 ; j = RESnum(iATMPCA(1))
        do i = 2,nPCA
          if ( RESnum(iATMPCA(i)) .ne. j ) then
            k = k + 1 ; iresPCA(k) = i ; j = RESnum(iATMPCA(i))
          endif
        enddo
        iresPCA(NresPCA+1) = nPCA + 1

        ! Making residue contact pair list
        Ncon_res_PCA = 0
        do i = 1,NresPCA-1
          ichn = CHN(iATMPCA(iresPCA(i)))
          do j = i+1,NresPCA
            jchn = CHN(iATMPCA(iresPCA(j)))
            if ( ichn.eq.jchn .and.                                    &
              RESnum(iATMPCA(iresPCA(j)))-RESnum(iATMPCA(iresPCA(i)))  &
              .le.cut_resnum_PCA ) cycle
            Ncon_res_PCA = Ncon_res_PCA + 1
          enddo
        enddo
        allocate(con_res_PCA(2,Ncon_res_PCA))
        ii = 0 ; jj = 0
        do i = 1,NresPCA-1
          ichn = CHN(iATMPCA(iresPCA(i)))
          do j = i+1,NresPCA
            jchn = CHN(iATMPCA(iresPCA(j)))
            if ( ichn.eq.jchn .and.                                    &
              RESnum(iATMPCA(iresPCA(j)))-RESnum(iATMPCA(iresPCA(i)))  &
              .le.cut_resnum_PCA ) cycle
            ii = ii + 1
            con_res_PCA(1:2,ii) = (/i,j/)
          enddo
        enddo

        write(6,*)
        write(6,'(4x,a,i0)')"+ Residue pairs for PCA internal cord. : "&
          ,Ncon_res_PCA
        if ( loglvl.eq."l" .and. Ncon_res_PCA.gt.6 ) then
          write(6,'(8x,a)')"   Res.1  -  Res.2"
          do k = 1,Ncon_res_PCA
            i = con_res_PCA(1,k) ; j = con_res_PCA(2,k)
            write(6,'(i8,a3,i0,a3,i0)')k," : ",                     &
              RESnum(iATMPCA(iresPCA(i)))," - ",                       &
              RESnum(iATMPCA(iresPCA(j)))
          enddo
        elseif ( loglvl .eq. "s" ) then
          write(6,'(8x,a)')"   Res.1  -  Res.2"
          do k = 1,3
            i = con_res_PCA(1,k) ; j = con_res_PCA(2,k)
            write(6,'(i8,a3,i0,a3,i0)')k," : ",                     &
              RESnum(iATMPCA(iresPCA(i)))," - ",                       &
              RESnum(iATMPCA(iresPCA(j)))
          enddo
          write(6,*)"              ..."
          do k = Ncon_res_PCA-2,Ncon_res_PCA
            i = con_res_PCA(1,k) ; j = con_res_PCA(2,k)
            write(6,'(i8,a3,i0,a3,i0)')k," : ",                     &
              RESnum(iATMPCA(iresPCA(i)))," - ",                       &
              RESnum(iATMPCA(iresPCA(j)))
          enddo
        endif
      endif

      return
      end subroutine PCAprep


!==================================================================


      subroutine PCA_init()

      use COMIFN ; use COMVAL
      implicit none

      open(unit=upc,file=trim(d_proj)//".PCAcod",status="replace",     &
           form="unformatted")
      if ( trim(PCA_method) .eq. "cord" ) then
        write(upc)nPCA*3
      elseif ( trim(PCA_method) .eq. "intercord" ) then
        write(upc)nPCA*(nPCA-1)/2
      else
        write(upc)Ncon_res_PCA
      endif
      return
      end subroutine PCA_init


!===================================================================


      subroutine PCA_final()

      use COMVAL
      implicit none

      close(upc)
      return
      end subroutine PCA_final
