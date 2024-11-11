
      subroutine outpdb(iwrite,iprint,filena,maxatm,numatm,atty,resty, &
                        resnum,occupy,bfact,maxcom,numcom,commnt,ier)
 
!***********************************************************************
!
!     Write PDB-formatted coordinate data Atom type, Residue type,
!       Residue number, Coordinate, Occupancy, & B-factor
!
!***********************************************************************

      use COMCMMC,only: cord

      implicit none

      ! Logical unit number for input and output
        integer(4),intent(in):: iwrite,iprint
      ! File name of coordinate
        character(*),intent(in):: filena
      ! Max number of atoms & number of atoms
        integer(4),intent(in):: maxatm,numatm
      ! Atom & Residue type of each atom
        character(*),intent(in):: atty(maxatm),resty(maxatm)
      ! Residue number of each atom
        integer(4),intent(in):: resnum(maxatm)
      ! Coordinate, Occupancy, & B-factor
        real(8),intent(in):: occupy(maxatm),bfact(maxatm)
      ! (Max) number of comment lines
        integer(4),intent(in):: maxcom
        integer(4),intent(inout):: numcom
      ! Comments
        character(*),intent(inout):: commnt(maxcom)
      ! Condition code ( 0: NO ERROR, -1:Number of atom s is not valid
      !                 -2: File open or close error
      !                 -3: data write error)
        integer(4),intent(out):: ier

      character(23):: curtim
      character(80):: userna
      integer(4):: len,icom,iatm
 
!*****************************************************

      if ( numatm.le.0 .or. numatm.gt.maxatm ) then
        ier = -1
        write(iprint,*)" "
        write(iprint,*)"ERROR> OUTPDB "
        write(iprint,*)"  NUMBER OF ATOMS IS NOT VALID "
        write(iprint,*)"    NUMBER OF ATOMS IS ",numatm
        write(iprint,*)" "
        return
      else
        ier = 0
      endif

      call flopen(iwrite,filena,12,'NULL',0,ier)
      if ( ier .ne. 0 ) then
        ier = -2
        write(iprint,*)" "
        write(iprint,*)"ERROR> OUTPDB "
        write(iprint,*)"  FILE OPEN OR CLOSE ERROR "
        write(iprint,*)" "
        return
      endif

      if ( numcom .lt. maxcom-1 ) then
        call infjob(curtim,userna,len)
        numcom = numcom + 1
        commnt(numcom) = " "
        commnt(numcom)(1:23) = curtim
        numcom = numcom + 1
        commnt(numcom) = " "
        commnt(numcom)(1:len) = userna(1:len)
      endif

      do icom = 1,numcom
        write(iwrite,'(a6,4x,a70)')"REMARK",commnt(icom)(1:70)
      enddo
 
      do iatm = 1,numatm
        if( atty(iatm)(4:4) .eq. ' ' ) then
          if ( resnum(iatm) .lt. 100000 ) then
            write(iwrite,8100,err=9200)                                &
                  iatm,atty(iatm)(1:3),resty(iatm)(1:4),resnum(iatm),  &
                  cord(1:3,iatm),occupy(iatm),bfact(iatm)
          else
            write(iwrite,8101,err=9200)                                &
                  iatm,atty(iatm)(1:3),resty(iatm)(1:4),resnum(iatm),  &
                  cord(1:3,iatm),occupy(iatm),bfact(iatm)
          endif
        else
          if ( resnum(iatm) .lt. 100000 ) then
            write(iwrite,8200,err=9200)                                &
                  iatm,atty(iatm)(1:4),resty(iatm)(1:4),resnum(iatm),  &
                  cord(1:3,iatm),occupy(iatm),bfact(iatm)
          else
            write(iwrite,8201,err=9200)                                &
                  iatm,atty(iatm)(1:4),resty(iatm)(1:4),resnum(iatm),  &
                  cord(1:3,iatm),occupy(iatm),bfact(iatm)
          endif
        endif
      enddo

8100  format('ATOM',i7,2x,a3,1x,a4,i5,4x,3f8.3,2f6.2)
8101  format('ATOM',i7,2x,a3,1x,a4,i6,3x,3f8.3,2f6.2)
8200  format('ATOM',i7,1x,a4,1x,a4,i5,4x,3f8.3,2f6.2)
8201  format('ATOM',i7,1x,a4,1x,a4,i6,3x,3f8.3,2f6.2)

      call flclos(iwrite,10,ier)
      if ( ier .ne. 0 ) then 
        ier = -2
        write(iprint,*)" "
        write(iprint,*)"ERROR> OUTPDB "
        write(iprint,*)"  FILE CLOSE ERROR "
        write(iprint,*)" "
        return
      endif

      return

!*******************************

9200  ier = -3
      write(iprint,*)" "
      write(iprint,*)"ERROR> OUTPDB "
      write(iprint,*)"  DATA WRITE ERROR "
      write(iprint,*)" "

      return
      end subroutine outpdb
