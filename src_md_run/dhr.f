       SUBROUTINE DHR(INCONS,IOUT,FILENA,IER)
C----------------------------------------------------------------------
C
C $SUBROUTINE DHR T.NAKAI 1989.10.23 VERSION1.3
C
C      ROUTINE TO SELECT ATOMS( I-J-K-L ) FOR DIHEDRAL CONSTRAINTS
C      AND SEARCH ATOM POINTERS( ABSOLUTE ATOM NUMBER ).
C
C    PERI  2-ND   T.NAKAI  1989. 5.31. VERSION 1.0
C                                6.12. VERSION 1.1
C                               10.16. VERSION 1.2
C                               10.23. VERSION 1.3
C
C   ( ARGUMENTS )
C
C      INCONS I*4  (I) LOGICAL UNIT NUMBER
C      IOUT   I*4  (I) LOGICAL UNIT NUMBER FOR OUTPUT DIAGNOSTICS
C      FILENA A*80 (I) FILE NAME
C      IER    I*4  (O) ERROR CODE
C
C----------------------------------------------------------------------
C  ERROR CODE
C----------------------------------------------------------------------
C      -1860 : INPUT DATA ERROR ( INPUT POSITIVE VALUE )
C      -1861 : ERROR IN FINDIH
C      -1862 : ERROR IN CHKDIH
C      -1863 : ERROR IN RAGTOR
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C  INCLUDE FILES
C----------------------------------------------------------------------

      use COMBAS ; use COMMIS ; use PHYCNS

C----------------------------------------------------------------------
C  VARIABLES
C----------------------------------------------------------------------
C    LOCAL VARIABLES FOR RDFREE
C
      implicit none

       INTEGER*4  MAXVAL
       PARAMETER(MAXVAL = 20)
       INTEGER*4    WRK1(MAXVAL),WRK2(MAXVAL),IVAL(MAXVAL)
       INTEGER*4    ICOL,EFCOL
       CHARACTER*1  INPMOD(MAXVAL)
       REAL*8       RVAL(MAXVAL)
       CHARACTER*80 CVAL(MAXVAL)
C
       INTEGER*4    INCONS, IOUT, IER
       INTEGER*4    ITDIHC, ICHAIN, IRESNM(4)
       CHARACTER*4  COMAND, SUBCOM, OUTPUT, ATMNAM(4)
       CHARACTER*80 LINE, FILENA, BLANK
       INTEGER*4    IERSUB, K
C----------------------------------------------------------------------
C  SET INITIAL PARAMETERS
C
       IER = 0
       ITDIHC = 0
       BLANK = ' '
C
C  OPEN UNIT 'INCONS' FOR READING DIHEDRAL CONSTRAINTS
C
       CALL FLOPEN(INCONS,FILENA,10,'NULL',0,IERSUB)
C
C  OPEN ERROR IN FLOPEN
C
       IF(IERSUB.LT.0) THEN
         WRITE(6,990) INCONS
         IER = IERSUB
         RETURN
       ENDIF
C
C  READ DIHEDRAL CONSTRAINTS ( FREE FORMAT )
C  MAIN LOOP TO READ DIHEDRAL CONSTRAINTS
C
 1000  CONTINUE
       READ(INCONS,'(A80)',END=500) LINE
       ICOL = EFCOL(LINE,BLANK,';')
       IF(ICOL.LE.0) GOTO 1000
       INPMOD(1) = 'C'
       INPMOD(2) = 'C'
C----------------------------------------------------------------------
       CALL RDFREE(LINE,ICOL,MAXVAL,2,INPMOD,WRK1,WRK2,
     +             RVAL,IVAL,CVAL,IERSUB)
       IF(IERSUB.LT.0) THEN
          WRITE(IOUT,970)
          IER = IERSUB
          RETURN
       ENDIF
       COMAND = CVAL(1)(1:4)
       SUBCOM = CVAL(2)(1:4)
       IF(COMAND.NE.'DHR>') GOTO 1000
C
C      WILL KEEP READING UNITL 'STOP' IS READ
C
       IF(SUBCOM.EQ.'STOP') GOTO 500
C
C      ALL ATOM PAIRS IN THIS CARD WILL BE PUT IN THE SAME GROUP AND
C      WILL KEEP READING UNITL 'END' IS READ
C
       IF(SUBCOM.EQ.'LIST') THEN
  100     CONTINUE
          READ(INCONS,'(A80)',END=500) LINE
          ICOL = EFCOL(LINE,BLANK,';')
          IF(ICOL.LE.0) GOTO 100
          INPMOD(1) = 'C'
          CALL RDFREE(LINE,ICOL,MAXVAL,2,INPMOD,WRK1,WRK2,
     +                RVAL,IVAL,CVAL,IERSUB)
          IF(IERSUB.LT.0) THEN
             WRITE(IOUT,970)
             IER = IERSUB
             RETURN
          ENDIF
          COMAND = CVAL(1)(1:4)
          IF(COMAND.EQ.'END ') GOTO 1000
C
          INPMOD(1) = 'I'
          INPMOD(2) = 'I'
          INPMOD(3) = 'C'
          INPMOD(4) = 'I'
          INPMOD(5) = 'C'
          INPMOD(6) = 'I'
          INPMOD(7) = 'C'
          INPMOD(8) = 'I'
          INPMOD(9) = 'C'
          INPMOD(10)= 'R'
          INPMOD(11)= 'R'
          INPMOD(12)= 'R'
          INPMOD(13)= 'R'
          INPMOD(14)= 'C'
C
          CALL RDFREE(LINE,ICOL,MAXVAL,14,INPMOD,WRK1,WRK2,
     +                RVAL,IVAL,CVAL,IERSUB)
          IF(IERSUB.LT.0) THEN
             WRITE(IOUT,970)
             IER = IERSUB
             RETURN
          ENDIF
C
C  INPUT DATA ERROR IN 'INCONS'
C
          IF(RVAL(1).LE.0.0D0) THEN
             WRITE(IOUT,960)
             IER = -1860
             RETURN
          ENDIF
C
          ITDIHC    = ITDIHC + 1
          ICHAIN    = IVAL(1)
          IRESNM(1) = IVAL(2)
          ATMNAM(1) = CVAL(1)
          IRESNM(2) = IVAL(3)
          ATMNAM(2) = CVAL(2)
          IRESNM(3) = IVAL(4)
          ATMNAM(3) = CVAL(3)
          IRESNM(4) = IVAL(5)
          ATMNAM(4) = CVAL(4)
          FUDELL(ITDIHC) = RVAL(1)
          FUDELU(ITDIHC) = RVAL(2)
          FUCDLW(ITDIHC) = RVAL(3)
          FUCDUP(ITDIHC) = RVAL(4)
          OUTPUT    = CVAL(5)
C
C    CONVERT INPUT DIHEDRAL ANGLE FROM DEGREE TO RADIAN
C
          IF(FUCDUP(ITDIHC).LT.0.0D0) THEN
            FUCDUP(ITDIHC) = ( FUCDUP(ITDIHC) + 360.0D0 )*RAD
          ELSE
            FUCDUP(ITDIHC) = FUCDUP(ITDIHC)*RAD
          ENDIF
 
          IF(FUCDLW(ITDIHC).LT.0.0D0) THEN
            FUCDLW(ITDIHC) = ( FUCDLW(ITDIHC) + 360.0D0 )*RAD
          ELSE
            FUCDLW(ITDIHC) = FUCDLW(ITDIHC)*RAD
          ENDIF
C  DEBUG
C         WRITE(IOUT,890) ((IRESNM(I),ATMNAM(I)),I=1,4),FUDELL(ITDIHC),
C    +                    FUDELU(ITDIHC),FUCDLW(ITDIHC),FUCDUP(ITDIHC)
C 890  FORMAT(1X,4(I6,1X,A4),4F8.3)
C
C  SEARCH ATOM POINTERS
C
          CALL FINDIH(ICHAIN,IRESNM,ATMNAM,OUTPUT,ITDIHC,IOUT,IERSUB)
          IF(IERSUB.LT.0) THEN
             WRITE(IOUT,950)
             IER = -1861
             RETURN
          ENDIF
C
          GOTO 100
C
       ELSE IF(SUBCOM.EQ.'NUMB') THEN
  200     CONTINUE
          READ(INCONS,'(A80)',END=500) LINE
          ICOL = EFCOL(LINE,BLANK,';')
          IF(ICOL.LE.0) GOTO 200
          INPMOD(1) = 'C'
          CALL RDFREE(LINE,ICOL,MAXVAL,2,INPMOD,WRK1,WRK2,
     +                RVAL,IVAL,CVAL,IERSUB)
          IF(IERSUB.LT.0) THEN
             WRITE(IOUT,970)
             IER = IERSUB
             RETURN
          ENDIF
          COMAND = CVAL(1)(1:4)
          IF(COMAND.EQ.'END ') GOTO 1000
C
          INPMOD(1) = 'I'
          INPMOD(2) = 'I'
          INPMOD(3) = 'I'
          INPMOD(4) = 'I'
          INPMOD(5) = 'R'
          INPMOD(6) = 'R'
          INPMOD(7) = 'R'
          INPMOD(8) = 'R'
          INPMOD(9) = 'C'
C
          CALL RDFREE(LINE,ICOL,MAXVAL,9,INPMOD,WRK1,WRK2,
     +                RVAL,IVAL,CVAL,IERSUB)
          IF(IERSUB.LT.0) THEN
             WRITE(IOUT,970)
             IER = IERSUB
             RETURN
          ENDIF
C
C  INPUT DATA ERROR IN 'INCONS'
C
          IF(RVAL(1).LE.0.0D0) THEN
             WRITE(IOUT,960)
             IER = -1860
             RETURN
          ENDIF
C
          ITDIHC = ITDIHC + 1
          IUDHCP(ITDIHC,1) = IVAL(1)
          IUDHCP(ITDIHC,2) = IVAL(2)
          IUDHCP(ITDIHC,3) = IVAL(3)
          IUDHCP(ITDIHC,4) = IVAL(4)
          FUDELL(ITDIHC)   = RVAL(1)
          FUDELU(ITDIHC)   = RVAL(2)
          FUCDLW(ITDIHC)   = RVAL(3)
          FUCDUP(ITDIHC)   = RVAL(4)
          OUTPUT = CVAL(1)
C
C    CONVERT INPUT DIHEDRAL ANGLE FROM DEGREE TO RADIAN
C
          IF(FUCDUP(ITDIHC).LT.0.0D0) THEN
            FUCDUP(ITDIHC) = ( FUCDUP(ITDIHC) + 360.0D0 )*RAD
          ELSE
            FUCDUP(ITDIHC) = FUCDUP(ITDIHC)*RAD
          ENDIF
 
          IF(FUCDLW(ITDIHC).LT.0.0D0) THEN
            FUCDLW(ITDIHC) = ( FUCDLW(ITDIHC) + 360.0D0 )*RAD
          ELSE
            FUCDLW(ITDIHC) = FUCDLW(ITDIHC)*RAD
          ENDIF
C  DEBUG
C         WRITE(IOUT,880) (IUDHCP(ITDIHC,I),I=1,4),FUDELL(ITDIHC),
C    +                    FUDELU(ITDIHC),FUCDLW(ITDIHC),FUCDUP(ITDIHC)
C 880  FORMAT(1X,4I6,4F8.3)
C
C  CHECK WHETHER DIHEDRAL ANGLE IS EXIST OR NOT IN THE MOLECULE
C
          CALL CHKDIH(ITDIHC,IOUT,IERSUB)
C
          IF(IERSUB.LT.0) THEN
             WRITE(IOUT,910)
             IER = -1862
             RETURN
          ENDIF
C
C  DUMP ATOM POINTERS ; SWITCHING IS AVAILABLE BY "OUTPUT"
C
          IF(OUTPUT(1:1).EQ.'Y') THEN
             IF(ITDIHC.LE.1) WRITE(IOUT,980)
             WRITE(IOUT,920) ITDIHC,(IUDHCP(ITDIHC,K),K=1,4),
     +                       FUDELL(ITDIHC),FUDELU(ITDIHC),
     +                       FUCDLW(ITDIHC)/RAD,FUCDUP(ITDIHC)/RAD
          ENDIF
C
          GOTO 200
C
       ENDIF
       GOTO 1000
C
C  END OF MAIN LOOP
C
  500  CONTINUE
C
C DEBUG
C         WRITE(IOUT,940) (I,(IUDHCP(I,J),J=1,4),I=1,ITDIHC)
          WRITE(IOUT,930) ITDIHC
C
          IUTDHC = ITDIHC
 
C
C  REARRANGE TORSION PAIR ATOM ARRAY
C
       CALL RAGTOR(IOUT,IERSUB)
 
       IF(IERSUB.LT.0) THEN
          WRITE(IOUT,900)
          IER = -1863
          RETURN
       ENDIF
 
C----------------------------------------------------------------------
C  FORMAT
C----------------------------------------------------------------------
  990  FORMAT(1X,'ERROR> DHR',/,5X,'OPEN FAILURE IN FLOPEN - UNIT ',
     +        I3)
  980  FORMAT(1X,'INFORMATION> DHR',/,5X,'SELECT FOLLOWING ATOMS FOR'
     +           ' DIHEDRAL CONSTRAINTS.')
  970  FORMAT(1X,'ERROR> DHR',/,5X,'READ ERROR IN RDFREE')
  960  FORMAT(1X,'ERROR> DHR',/,5X,'INPUT DATA ERROR. INPUT '
     +           'POSITIVE VALUES.')
  950  FORMAT(1X,'ERROR> DHR',/,5X,'ERROR FOUND IN FIDIH.')
  940  FORMAT( (1H ,2(5I6)))
  930  FORMAT(1X,'INFORMATION> DHR',/,5X,
     +           'TOTAL NUMBER OF DIHEDRAL CONSTRAINTS ARE :',I6)
  920  FORMAT(5I6,4F9.3)
  910  FORMAT(1X,'ERROR> DHR',/,5X,'ERROR FOUND IN CHKDIH.')
  900  FORMAT(1X,'ERROR> DHR',/,5X,'ERROR FOUND IN RAGTOR.')
 
       CALL  FLCLOS( INCONS , 10 , IER )
 
       RETURN
       END
C----------------------------------------------------------------------
       SUBROUTINE  FINDIH(ICHAIN,IRESNM,ATMNAM,OUTPUT,ITDIHC,IOUT,IER)
C----------------------------------------------------------------------
C
C $SUBROUTINE FINDIH T.NAKAI 1989.10.30 VERSION1.2
C
C  SEARCH ATOM POINTERS( I-J-K-L ) FOR DIHEDRAL CONSTRAINTS
C
C    PERI 2-ND  T.NAKAI  1989. 6.12 VERSION 1.0
C                        1989.10.30 VERSION 1.2
C
C  ( ARGUMENTS )
C
C     ICHAIN   I*4  (I)  THE CHAIN NUMBER FOR SELECTION
C     IRESNM   I*4  (I)  THE RESIDUE NUMBER ( I-J-K-L )
C     ATMNAM   A*4  (I)  THE ATOM NAME ( I-J-K-L )
C     OUTPUT   A*4  (I)  OPTION TO OUTPUT INFORMATION
C     ITDIHC   I*4  (M)  TOTAL NUMBER OF DIHEDRAL CONSTRAINTS
C     IOUT     I*4  (I)  LOGICAL UNIT NUMBER TO OUTPUT INFORMATIONS
C     IER      I*4  (O)  ERROR CODE
C
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C  INCLUDE FILES
C----------------------------------------------------------------------

      use COMBAS ; use COMERG ; use COMMIS ; use PHYCNS

C----------------------------------------------------------------------
C  VARIABLES
C----------------------------------------------------------------------
      implicit none

       INTEGER*4    ITDIHC, ICHAIN, IRESNM(4)
       CHARACTER*4  OUTPUT, ATMNAM(4)
       INTEGER*4    IER, IOUT, IFIND
       INTEGER*4    IERSUB, I, J, K
C----------------------------------------------------------------------
C  SEARCH ATOM POINTERS OF DIHEDRAL ANGLES ( I-J-K-L )
C
C  DEBUG
C      WRITE(IOUT,890) ITDIHC,ICHAIN,((IRESNM(J),ATMNAM(J)),J=1,4)
C 890  FORMAT(1X,2I6,4(I6,1X,A4))
C
       DO 100 J=1,4
          IFIND = 0
          DO 200 I=1,IXNATM
             IF(IXACHN(I).NE.ICHAIN) GOTO 200
             IF( (IXARES(I).EQ.IRESNM(J)) .AND.
     +           (CXATMN(I).EQ.ATMNAM(J)) ) THEN
                 IFIND = IFIND + 1
                 IUDHCP(ITDIHC,J) = I
             ENDIF
  200     CONTINUE
          IF(IFIND.NE.1) WRITE(IOUT,990) ITDIHC,ICHAIN,IRESNM(J),
     +                                   ATMNAM(J)
  100  CONTINUE
C
C  DUMP ATOM POINTERS ; SWITCHING IS AVAILABLE BY "OUTPUT"
C
       IF(OUTPUT(1:1).EQ.'Y') THEN
          IF(ITDIHC.LE.1) WRITE(IOUT,960)
          WRITE(IOUT,980) ITDIHC,(IUDHCP(ITDIHC,K),K=1,4),
     +                    FUDELL(ITDIHC),FUDELU(ITDIHC),
     +                    FUCDLW(ITDIHC)/RAD,FUCDUP(ITDIHC)/RAD
       ENDIF
C
C  CHECK WHETHER DIHEDRAL ANGLE IS EXIST OR NOT IN THE MOLECULE
C
       CALL CHKDIH(ITDIHC,IOUT,IERSUB)
C
       IF(IERSUB.LT.0) THEN
          WRITE(IOUT,970)
          IER = -1862
          RETURN
       ENDIF
C----------------------------------------------------------------------
C  FORMAT
C----------------------------------------------------------------------
  990  FORMAT(1X,'ERROR> FINDIH',/,5X,'ERROR FINDING ATOM OF DIHEDRAL',
     +           ' ANGLE.',3I6,1X,A4)
  980  FORMAT(5I6,4F9.3)
  970  FORMAT(1X,'ERROR> FINDIH',/,5X,'ERROR FOUND IN CHKDIH.')
  960  FORMAT(1X,'INFORMATION> DHR',/,5X,'SELECT FOLLOWING ATOMS FOR'
     +           ' DIHEDRAL CONSTRAINTS.')
       RETURN
       END
C----------------------------------------------------------------------
       SUBROUTINE  CHKDIH(ITDIHC,IOUT,IER)
C----------------------------------------------------------------------
C
C $SUBROUTINE CHKDIH T.NAKAI 1989.06.12 VERSION1.1
c  Modified by K.Teraishi 1993.8.12
c  To cope with the multiply defined torsion in tpl file.
C
C  CHECK WHETHER DIHEDRAL ANGLE IS EXIST OR NOT IN THE MOLECULE
C
C  ( ARGUMENTS )
C
C     ITDIHC   I*4  (M)  TOTAL NUMBER OF DIHEDRAL CONSTRAINTS
C     IOUT     I*4  (I)  LOGICAL UNIT NUMBER TO OUTPUT INFORMATIONS
C     IER      I*4  (O)  ERROR CODE
C
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C  INCLUDE FILES
C----------------------------------------------------------------------

      use COMBAS ; use COMERG ; use COMMIS

C----------------------------------------------------------------------
C  VARIABLES
C----------------------------------------------------------------------

      implicit none

       INTEGER*4    ITDIHC, IER, IOUT, ICHECK
C#930107
       INTEGER*4    I
C----------------------------------------------------------------------
C  CHECK DIHEDRAL ANGLES
C
       ICHECK = 0
       DO 300 I=1,IYNTOR
          IF( (IYPTOR(1,I).EQ.IUDHCP(ITDIHC,1)) .AND.
     +        (IYPTOR(2,I).EQ.IUDHCP(ITDIHC,2)) .AND.
     +        (IYPTOR(3,I).EQ.IUDHCP(ITDIHC,3)) .AND.
     +        (IYPTOR(4,I).EQ.IUDHCP(ITDIHC,4)) ) THEN
c             ICHECK = ICHECK + 1
              ICHECK = 1
          ELSE IF( (IYPTOR(4,I).EQ.IUDHCP(ITDIHC,1)) .AND.
     +             (IYPTOR(3,I).EQ.IUDHCP(ITDIHC,2)) .AND.
     +             (IYPTOR(2,I).EQ.IUDHCP(ITDIHC,3)) .AND.
     +             (IYPTOR(1,I).EQ.IUDHCP(ITDIHC,4)) ) THEN
c             ICHECK = ICHECK + 1
              ICHECK = 1
          ENDIF
  300  CONTINUE
C
C  ERROR/RETURN
C
       IF(ICHECK.NE.1) THEN
          WRITE(IOUT,970) ITDIHC
          IER = -1862
          RETURN
       ENDIF
C----------------------------------------------------------------------
C  FORMAT
C----------------------------------------------------------------------
  970  FORMAT(1X,'ERROR> CHKDIH',/,5X,'NO SUCH DIHEDRAL ANGLE WAS',
     +           ' FOUND IN MOLECULE.',3I6)
       RETURN
       END
C----------------------------------------------------------------------
      SUBROUTINE  RAGTOR(IOUT,IER)
C----------------------------------------------------------------------
C
C $SUBROUTINE RAGTOR T.NAKAI 1989.10.19 VERSION1.1
C
C  REARRANGE TORSION PAIR ATOM ARRAY
C
C   PERI 2-ND  T.NAKAI  1989.10.16  VERSION 1.0
C                       1989.10.19  VERSION 1.1
C
C  ( ARGUMENTS )
C
C     IOUT     I*4  (I)  LOGICAL UNIT NUMBER TO OUTPUT INFORMATIONS
C     IER      I*4  (O)  ERROR CODE
C
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C  INCLUDE FILES
C----------------------------------------------------------------------

      use COMBAS ; use COMERG ; use COMMIS

C----------------------------------------------------------------------
C  LOCAL VARIABLES
C----------------------------------------------------------------------
      implicit none

      INTEGER*4  ICDIHE, IRAGTR, ICOUNT
      INTEGER*4  IATOM, JATOM, KATOM, LATOM, IOUT, IER
      INTEGER*4  ISKIP(MAXTOR),JSKIP(MAXTOR)
C#930107
      INTEGER*4  I, J
C----------------------------------------------------------------------
 
      ICDIHE = 0
      IRAGTR = 0
      ICOUNT = 0
 
      DO 50 I=1,MAXTOR
         JSKIP(I) = 0
   50 CONTINUE
 
      DO 100 I=1,IUTDHC
         IATOM = IUDHCP(I,1)
         JATOM = IUDHCP(I,2)
         KATOM = IUDHCP(I,3)
         LATOM = IUDHCP(I,4)
         DO 200 J=1,IYNTOR
C           IF((IYPTOR(1,J).EQ.IATOM).AND.(IYPTOR(2,J).EQ.JATOM).AND.
C    +         (IYPTOR(3,J).EQ.KATOM).AND.(IYPTOR(4,J).EQ.LATOM)) THEN
            IF((IYPTOR(2,J).EQ.JATOM).AND.(IYPTOR(3,J).EQ.KATOM)) THEN
               ICDIHE = ICDIHE + 1
               ISKIP(ICDIHE) = J
            ENDIF
  200    CONTINUE
  100 CONTINUE
 
      DO 300 I=1,IYNTOR
         DO 350 J=1,ICDIHE
            IF(I.EQ.ISKIP(J)) THEN
              JSKIP(I) = -1
            ENDIF
  350   CONTINUE
  300 CONTINUE
 
      IF(IYNTOR+ICDIHE.GT.MAXTOR) THEN
         WRITE(IOUT,910) IYNTOR+ICDIHE-MAXTOR
         IER = -1910
         RETURN
      ENDIF
 
      DO 150 I=1,IYNTOR
         IF(JSKIP(I).LT.0) THEN
            IRAGTR = IRAGTR + 1
            IYPTOR(1,IYNTOR+IRAGTR) = IYPTOR(1,I)
            IYPTOR(2,IYNTOR+IRAGTR) = IYPTOR(2,I)
            IYPTOR(3,IYNTOR+IRAGTR) = IYPTOR(3,I)
            IYPTOR(4,IYNTOR+IRAGTR) = IYPTOR(4,I)
            IYTDIV(IYNTOR+IRAGTR) = IYTDIV(I)
            IYTNBF(IYNTOR+IRAGTR) = IYTNBF(I)
            FYFTOR(IYNTOR+IRAGTR) = FYFTOR(I)
            FYTROT(IYNTOR+IRAGTR) = FYTROT(I)
            FYTPHS(IYNTOR+IRAGTR) = FYTPHS(I)
            FYTVWS(IYNTOR+IRAGTR) = FYTVWS(I)
            FYTESS(IYNTOR+IRAGTR) = FYTESS(I)
        ENDIF
  150 CONTINUE
 
      IF(ICDIHE.NE.IRAGTR) THEN
         WRITE(IOUT,900)
         IER = - 1900
         RETURN
      ENDIF
 
C---------------------------------------------------------------------C
C DEBUG
C     WRITE(IOUT,*) ' ICDIHE : ITUDHC ', ICDIHE,IUTDHC
C     WRITE(IOUT,*) ' IYNTOR : IRAGTR ', IYNTOR,IRAGTR
C     WRITE(IOUT,800) ( (IUDHCP(I,J),J=1,4),I=1,IUTDHC )
C     WRITE(IOUT,810) ( IYPTOR(1,I),IYPTOR(2,I),IYPTOR(3,I),IYPTOR(4,I),
C    +                  IYTDIV(I),IYTNBF(I),
C    +                  FYFTOR(I),FYTROT(I),FYTPHS(I),FYTVWS(I),
C    +                  FYTESS(I),I=1,IYNTOR )
C 800 FORMAT(1H ,4I6)
C 810 FORMAT(1H ,6I6,5F8.3)
C---------------------------------------------------------------------C
 
C  REARRANGE ARRAY
 
      DO 400 I=1,IYNTOR+IRAGTR
         IF(JSKIP(I).GE.0) THEN
            ICOUNT = ICOUNT + 1
            IYPTOR(1,ICOUNT) = IYPTOR(1,I)
            IYPTOR(2,ICOUNT) = IYPTOR(2,I)
            IYPTOR(3,ICOUNT) = IYPTOR(3,I)
            IYPTOR(4,ICOUNT) = IYPTOR(4,I)
            IYTDIV(ICOUNT) = IYTDIV(I)
            IYTNBF(ICOUNT) = IYTNBF(I)
            FYFTOR(ICOUNT) = FYFTOR(I)
            FYTROT(ICOUNT) = FYTROT(I)
            FYTPHS(ICOUNT) = FYTPHS(I)
            FYTVWS(ICOUNT) = FYTVWS(I)
            FYTESS(ICOUNT) = FYTESS(I)
         ENDIF
  400 CONTINUE
C---------------------------------------------------------------------C
C DEBUG
C     WRITE(IOUT,*) ' ICDIHE : ITUDHC ', ICDIHE,IUTDHC
C     WRITE(IOUT,*) ' IYNTOR : IRAGTR ', IYNTOR,IRAGTR
C     WRITE(IOUT,800) ( (IUDHCP(I,J),J=1,4),I=1,IUTDHC )
C     WRITE(IOUT,810) ( IYPTOR(1,I),IYPTOR(2,I),IYPTOR(3,I),IYPTOR(4,I),
C    +                  IYTDIV(I),IYTNBF(I),
C    +                  FYFTOR(I),FYTROT(I),FYTPHS(I),FYTVWS(I),
C    +                  FYTESS(I),I=1,IYNTOR )
C 800 FORMAT(1H ,4I6)
C 810 FORMAT(1H ,6I6,5F8.3)
C---------------------------------------------------------------------C
 
      WRITE(IOUT,920)
      WRITE(IOUT,930) ( IYPTOR(1,I),IYPTOR(2,I),IYPTOR(3,I),IYPTOR(4,I),
     +                  I=IYNTOR+1,IYNTOR+IRAGTR )
 
      IYNDTR = IRAGTR
 
  900 FORMAT(/' ERROR> RAGTOR',/
     +        '        NUMBER OF DIHEDRAL CONSTRAINTS IS MISMATCHED.'/)
  910 FORMAT(/' ERROR> RAGTOR',/
     +        '        MAXTOR IS EXCEEDED BY ',I6/)
  920 FORMAT(/' INFORMATION> RAGTOR',/
     +        '        FOLLOWING DIHEDRALS SHOUD BE RESTRICTED.'/)
  930 FORMAT(10X,4I6)
      RETURN
      END
