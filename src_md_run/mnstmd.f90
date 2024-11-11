
      subroutine mnstmd(iprint,ier,idfact,cord,idutrj,idloop,idflog )

      use COMBAS ; use COMERG ; use COMPSC ; use PHYCNS

      implicit none

      integer(4):: iprint,ier,idutrj,idloop,idfact,idflog
      real(8):: cord(3,ixnatm)

      integer(4):: maxbuf,nseq,i
      real(8):: rd2deg,tcod(3,3),r12,r32,r23,px12,py12,pz12,px23,py23, &
                pz23,px123,py123,pz123,s32123,cosp
      real(8),allocatable:: buf(:)

!************************************************

      if ( idutrj .le. 0 ) then
        ier = -1 ; return
      else
        ier = 0 ; rd2deg = 180.d0 / pi ; nseq = 0
        maxbuf = 3*nacntr+ndcntr+ngcntr+ntcntr
        allocate(buf(maxbuf))
      endif

!     calculate coordinates information
      do i = 1,nacntr
        buf(nseq+1:nseq+3) = cord(1:3,ianatm(i))
        nseq = nseq + 3
      enddo

!     distance information
      do i = 1,ndcntr
        nseq = nseq + 1
        tcod(1:3,1) = cord(1:3,idnatm(1,i))-cord(1:3,idnatm(2,i))
        buf(nseq) = sqrt(dot_product(tcod(1:3,1),tcod(1:3,1)))
      enddo

!     angle information
      do i = 1,ngcntr
        tcod(1:3,1) = cord(1:3,ignatm(1,i))-cord(1:3,ignatm(2,i))
        tcod(1:3,2) = cord(1:3,ignatm(3,i))-cord(1:3,ignatm(2,i))
        r12 = dot_product(tcod(1:3,1),tcod(1:3,1))
        r32 = dot_product(tcod(1:3,2),tcod(1:3,2))
        if ( r12*r32 .ne. 0.d0 ) then
          cosp = dot_product(tcod(1:3,1),tcod(1:3,2)) / sqrt(r12*r32)
          if ( cosp .ge. 1.d0 ) then
            cosp = 1.d0
          elseif ( cosp .le. -1.d0 ) then
            cosp = -1.d0
          endif
          nseq = nseq + 1
          buf(nseq) = acos(cosp)*rd2deg
        else
          nseq = nseq + 1
          buf(nseq) = 9999.d0
        endif
      enddo

!     torsion information
      do i = 1,ntcntr
        tcod(1:3,1) = cord(1:3,itnatm(2,i))-cord(1:3,itnatm(1,i))
        tcod(1:3,2) = cord(1:3,itnatm(3,i))-cord(1:3,itnatm(2,i))
        tcod(1:3,3) = cord(1:3,itnatm(4,i))-cord(1:3,itnatm(3,i))
        px12 = tcod(2,1)*tcod(3,2) - tcod(2,2)*tcod(3,1)
        py12 = tcod(3,1)*tcod(1,2) - tcod(3,2)*tcod(1,1)
        pz12 = tcod(1,1)*tcod(2,2) - tcod(1,2)*tcod(2,1)
        px23 = tcod(2,2)*tcod(3,3) - tcod(2,3)*tcod(3,2)
        py23 = tcod(3,2)*tcod(1,3) - tcod(3,3)*tcod(1,2)
        pz23 = tcod(1,2)*tcod(2,3) - tcod(1,3)*tcod(2,2)
        r12 = px12**2 + py12**2 + pz12**2
        r23 = px23**2 + py23**2 + pz23**2
        if ( r12*r23 .ne. 0.d0 ) then
          cosp = (px12*px23 + py12*py23 + pz12*pz23) / sqrt(r12*r23)
          if ( cosp .ge. 1.d0 ) then
            cosp = 1.d0
          elseif ( cosp .le. -1.d0 ) then
            cosp = -1.d0
          endif
          px123 = py12*pz23 - py23*pz12
          py123 = pz12*px23 - pz23*px12
          pz123 = px12*py23 - px23*py12
          s32123 = tcod(1,2)*px123 + tcod(2,2)*py123 + tcod(3,2)*pz123
          nseq = nseq + 1 ; buf(nseq) = acos(cosp)*rd2deg
          if ( s32123 .lt. 0.d0 ) buf(nseq) = -buf(nseq)
        else
          nseq = nseq + 1 ; buf(nseq) = 9999.d0
        endif
      enddo    

!     write information to file in binary
      if ( nseq .gt. 0 ) then
        select case (idfact)
        case (1)
          write(idutrj,err=90)idloop,nacntr,ndcntr,ngcntr,ntcntr,nseq
          write(idutrj,err=90)sngl(buf(1:nseq))
        case (2)
          write(idutrj,err=90)idloop,nacntr,ndcntr,ngcntr,ntcntr,nseq
          write(idutrj,err=90)buf(1:nseq)
        case default
          write(idutrj,*,err=90)idloop,nacntr,ndcntr,ngcntr,ntcntr,nseq
          write(idutrj,*,err=90)buf(1:nseq)
        end select
      endif

      if ( idflog .eq. 2 ) then
        write(iprint,*)' '
        write(iprint,*)'INFORMATION> MNSTMD'
        write(iprint,*)'   ',nseq,'     STRUCTURAL PARAMETERS'
        write(iprint,*)'     AT ',idloop,' TH STEP WERE WRITTEM,'
        select case (idfact)
        case (1)
          write(iprint,*)'     WITH SINGLE PRECISION.'
        case (2)
          write(iprint,*)'     WITH DOUBLE PRECISION.'
        case default
          write(iprint,*)'     BY ASCII DATA.'
        end select
      endif

!*******************************

      return
90    write(iprint,*)' '
      write(iprint,*)'ERROR> MNSTMD'
      write(iprint,*)'       IN WRITING STRUCTURAL PARAMETERS'
      write(iprint,*)'       AT ',idloop,' TH STEP.'
      ier = -1 ; return

      end subroutine mnstmd
