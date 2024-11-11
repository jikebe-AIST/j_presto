 
C*******************************************************************
C
C     SUBROUTINE  INPVAR
C                 THIS SUBROUTINE IS FOR READ CONTROL DATA FOR
C                 SETTING VARIABLES AND MAKE VARIABLE TABLE
C
C                 PERI  2ND. DEP.  K.MORIKAMI 1989. 4. 3 VERSION 1.0
C                                             1989. 5.17 VERSION 1.1
C                                             1989. 5.24 VERSION 2.0
C
C*******************************************************************
 
      SUBROUTINE  INPVAR( IREAD,IPRINT,FILENA,IER )

      use COMBAS ; use COMERG
 
C===================================================================
C
C     ** ARGUNENTS **
C              (I) INPUT ARGUMENT
C              (O) OUTPUT ARGUMENT
C              (M) MODIFY ARGUMENT
C       IREAD   I*4       (I)   LOGICAL UNIT NUMBER FOR READ
C                               CONTROL DATA OF SETTING VARIABLES
C       FILENA  A*80      (I)   FILE NAME OF CONTROL DATA OF CAP
C       IPRINT  I*4       (I)   LOGICAL UNIT NUMBER FOR OUTPUT LOG
C       IER     I*4       (O)   CONDITION CODE
C                                0    ;  NO-ERROR
C
C===================================================================

      implicit none
 
      INTEGER*4     IREAD,IPRINT,IER
      CHARACTER*80  FILENA
 
      INTEGER*4     MAXVAL
      PARAMETER     ( MAXVAL = 9 )
      INTEGER*4     IWORK1(MAXVAL),IWORK2(MAXVAL),IVAL(MAXVAL)
      CHARACTER*1   INPMOD(MAXVAL)
      CHARACTER*80  CVAL(MAXVAL)
      REAL*8        RVAL(MAXVAL)
      INTEGER*4     EFCOL
      CHARACTER*80  LINE,SPACE
      LOGICAL*4     ONLIST,ONRADI,ONOUTL
      INTEGER*4     IERSUB, IATM, ICOL
 
C===================================================================
 
      IER    = 0
      ONLIST = .FALSE.
      ONRADI = .FALSE.
      ONOUTL = .FALSE.
      SPACE  = ' '
      DO 100 IATM = 1 , IXNATM
        IYTVAR(IATM) = 1
 100  CONTINUE
 
C     * OPEN CONTROL DATA FOR SETTING VARIABLES *
 
      CALL  FLOPEN( IREAD,FILENA,10,'NULL',0,IERSUB )
      IF( IERSUB .NE. 0 ) THEN
        WRITE(IPRINT,*) ' '
        WRITE(IPRINT,*) 'ERROR> INPVAR '
        WRITE(IPRINT,*) '       FILE OPEN ERROR '
        WRITE(IPRINT,*) ' '
        IER = IERSUB
        RETURN
      ENDIF
 
C     * READ CONTROL DATA FOR SETTING VARIABLES *
 
 200  CONTINUE
        READ(IREAD,'(A80)',END=300) LINE
        ICOL = EFCOL( LINE,SPACE,';' )
        IF( ICOL .LE. 0 ) GOTO 200
        IF( INDEX(LINE(1:ICOL),'SETVAR>') .NE. 0 ) THEN
          ONLIST = .FALSE.
          ONRADI = .FALSE.
          IF( INDEX(LINE(1:ICOL),'LIST') .NE. 0 ) ONLIST = .TRUE.
          IF( INDEX(LINE(1:ICOL),'RADI') .NE. 0 ) ONRADI = .TRUE.
          GOTO 200
        ENDIF
 
C       * READ DIRECT LIST INPUT *
 
        IF( ONLIST ) THEN
          CALL  RDVLIS( MAXVAL,IER,IPRINT,LINE,ICOL,
     +                  IVAL,IWORK1,IWORK2,RVAL,CVAL,INPMOD,
     +                  ONOUTL )
          IF( IER .NE. 0 ) RETURN
          GOTO 200
        ENDIF
 
C       * READ RADIUS DATA *
 
        IF( ONRADI ) THEN
          CALL  RDVRAD( MAXVAL,IER,IPRINT,LINE,ICOL,
     +                  IVAL,IWORK1,IWORK2,RVAL,CVAL,INPMOD,
     +                  ONOUTL )
          IF( IER .NE. 0 ) RETURN
          GOTO 200
        ENDIF
 
C     * COUNT NUMBER OF VARIABLES *
 
 300  CONTINUE
        IYNVAR = 0
        DO 310 IATM = 1 , IXNATM
          IF( IYTVAR(IATM) .EQ. 1 ) THEN
            IYNVAR         = IYNVAR + 1
            IYLVAR(IYNVAR) = IATM
          ENDIF
 310    CONTINUE
        IF( IYNVAR .LE. 0 ) THEN
          WRITE(IPRINT,*) ' '
          WRITE(IPRINT,*) 'ERROR> INPVAR '
          WRITE(IPRINT,*) '       NUMBER OF VARIABLES IS ZERO '
          WRITE(IPRINT,*) ' '
          IER = -1821
          RETURN
        ELSE
          WRITE(IPRINT,*) ' '
          WRITE(IPRINT,*) 'INFORMATION> INPVAR '
          WRITE(IPRINT,*) '        NUMBER OF FREE ATOMS IS ',IYNVAR
          WRITE(IPRINT,*) ' '
        ENDIF
C
        IF( ONOUTL ) THEN
          WRITE(IPRINT,*) '   LIST OF FREE ATOM '
          WRITE(IPRINT,*) ' '
          DO 320 IATM = 1 , IXNATM
            IF( IYTVAR(IATM) .EQ. 1 ) THEN
              WRITE(IPRINT,'(I7,5X,A8,I5,A8,I5)')
     +        IATM,CXATMN(IATM),IXARES(IATM),CXRESN(IATM),IXACHN(IATM)
            ENDIF
 320      CONTINUE
        ENDIF
 
      CALL  FLCLOS( IREAD , 10 , IER )
 
      RETURN
      END
 
C*******************************************************************
C
C     SUBROUTINE  RDVLIS
C                 THIS SUBROUTINE IS FOR READ DIRECT LIST INPUT
C
C                 PERI  2ND. DEP.  K.MORIKAMI 1989. 5.24 VERSION 1.0
C
C*******************************************************************
 
      SUBROUTINE  RDVLIS( MAXVAL,IER,IPRINT,LINE,ICOL,
     +                    IVAL,IWORK1,IWORK2,RVAL,CVAL,INPMOD,
     +                    ONOUTL )

      use COMBAS ; use COMERG 
 
C===================================================================
C
C     ** ARGUNENTS **
C              (I) INPUT ARGUMENT
C              (O) OUTPUT ARGUMENT
C              (M) MODIFY ARGUMENT
C              (W) WORK ARGUMENT
C       MAXVAL  I*4       (I)   MAX NUMBER OF VARIABLES
C       IER     I*4       (O)   CONDITION CODE
C       IPRINT  I*4       (I)   LOGICAL UNIT NUMBER FOR OUTPUT LOG
C       LINE    A*(*)     (I)   INPUT LINE
C       ICOL    I*4       (I)   EFFECTIVE COLUMN NUMBER
C       IVAL    I*4       (W)   WORK FOR RDFREE
C       IWORK1  I*4       (W)   WORK FOR RDFREE
C       IWORK2  I*4       (W)   WORK FOR RDFREE
C       RVAL    R*8       (W)   WORK FOR RDFREE
C       CVAL    A*80      (W)   WORK FOR RDFREE
C       INPMOD  A*1       (W)   WORK FOR RDFREE
C       ONOUTL  L*4       (O)   FLAG FOR OUTPUT LIST
C
C*******************************************************************

      implicit none
 
      INTEGER*4       MAXVAL,IER,IPRINT,ICOL
      INTEGER*4       IVAL(MAXVAL),IWORK1(MAXVAL),IWORK2(MAXVAL)
      REAL*8          RVAL(MAXVAL)
      CHARACTER*(*)   LINE
      CHARACTER*80    CVAL(MAXVAL)
      CHARACTER*1     INPMOD(MAXVAL)
      LOGICAL*4       ONOUTL
      INTEGER*4       IERSUB, ISTCHN, IENCHN, ISTRES, IENRES, ICHN, 
     +                JCOL, ISTATM, IENATM, IATM
 
C*******************************************************************
 
      IER       = 0
      INPMOD(1) = 'C'
      INPMOD(2) = 'I'
      INPMOD(3) = 'I'
      INPMOD(4) = 'I'
      INPMOD(5) = 'I'
      INPMOD(6) = 'C'
      INPMOD(7) = 'C'
      CALL  RDFREE( LINE,ICOL,MAXVAL,7,INPMOD,IWORK1,IWORK2,
     +              RVAL,IVAL,CVAL,IERSUB )
      IF( IERSUB .LT. 0 ) THEN
        WRITE(IPRINT,*) ' '
        WRITE(IPRINT,*) 'ERROR> INPVAR '
        WRITE(IPRINT,*) '       DATA TYPE ERROR '
        WRITE(IPRINT,*) ' '
        WRITE(IPRINT,'(A80)') LINE
        WRITE(IPRINT,*) ' '
        IER = IERSUB
        RETURN
      ENDIF
      IF( INDEX(CVAL(3),'YES') .NE. 0 ) ONOUTL = .TRUE.
 
C     * SEARCH ATOMS TO BE FIXED/FREE *
 
      ISTCHN = IVAL(1)
      IENCHN = IVAL(2)
      IF( ( ISTCHN .GT. IENCHN ) .OR.
     +    ( ISTCHN .GT. IXNCHN )       ) THEN
        WRITE(IPRINT,*) ' '
        WRITE(IPRINT,*) 'ERROR> INPVAR '
        WRITE(IPRINT,*) '    CHAIN NUMBER IS ILLEGAL '
        WRITE(IPRINT,*) ' '
        IER = -1820
        RETURN
      ENDIF
      ISTCHN = MAX( ISTCHN , 1      )
      IENCHN = MIN( IENCHN , IXNCHN )
      ISTRES = IVAL(3)
      IENRES = IVAL(4)
      JCOL = INDEX(CVAL(2)(1:8),'*') - 1
      DO 100 ICHN = ISTCHN , IENCHN
        IF( ICHN .EQ. 1 ) THEN
          ISTATM = 1
        ELSE
          ISTATM = IXCEND(ICHN-1) + 1
        ENDIF
        IENATM = IXCEND(ICHN)
        DO 110 IATM = ISTATM , IENATM
          IF( ( ( IXARES(IATM) .GE. ISTRES ) .OR.
     +          ( ISTRES .LT. 0            )       ) .AND.
     +        ( ( IXARES(IATM) .LE. IENRES ) .OR.
     +          ( IENRES .LT. 0            )       )       ) THEN
            IF( JCOL .EQ. -1 ) THEN
              IF( CXATMN(IATM) .EQ. CVAL(2)(1:8) ) THEN
                IF( CVAL(1)(1:3) .EQ. 'FIX'  ) IYTVAR(IATM) = 0
                IF( CVAL(1)(1:4) .EQ. 'FREE' ) IYTVAR(IATM) = 1
               ENDIF
            ENDIF
            IF( JCOL .EQ. 0 ) THEN
                IF( CVAL(1)(1:3) .EQ. 'FIX'  ) IYTVAR(IATM) = 0
                IF( CVAL(1)(1:4) .EQ. 'FREE' ) IYTVAR(IATM) = 1
            ENDIF
            IF( JCOL .GE. 1 ) THEN
              IF( INDEX(CXATMN(IATM),CVAL(2)(1:JCOL)) .EQ. 1 ) THEN
                IF( CVAL(1)(1:3) .EQ. 'FIX'  ) IYTVAR(IATM) = 0
                IF( CVAL(1)(1:4) .EQ. 'FREE' ) IYTVAR(IATM) = 1
              ENDIF
            ENDIF
          ENDIF
 110    CONTINUE
 100  CONTINUE
 
      RETURN
      END
 
C*******************************************************************
C
C     SUBROUTINE  RDVRAD
C                 THIS SUBROUTINE IS FOR READ RADIUS DATA
C
C                 PERI  2ND. DEP.  K.MORIKAMI 1989. 5.24 VERSION 1.0
C
C*******************************************************************
 
      SUBROUTINE  RDVRAD( MAXVAL,IER,IPRINT,LINE,ICOL,
     +                    IVAL,IWORK1,IWORK2,RVAL,CVAL,INPMOD,
     +                    ONOUTL )

      use COMBAS ; use COMERG 
 
C===================================================================
C
C     ** ARGUNENTS **
C              (I) INPUT ARGUMENT
C              (O) OUTPUT ARGUMENT
C              (M) MODIFY ARGUMENT
C              (W) WORK ARGUMENT
C       MAXVAL  I*4       (I)   MAX NUMBER OF VARIABLES
C       IER     I*4       (O)   CONDITION CODE
C       IPRINT  I*4       (I)   LOGICAL UNIT NUMBER FOR OUTPUT LOG
C       LINE    A*(*)     (I)   INPUT LINE
C       ICOL    I*4       (I)   EFFECTIVE COLUMN NUMBER
C       IVAL    I*4       (W)   WORK FOR RDFREE
C       IWORK1  I*4       (W)   WORK FOR RDFREE
C       IWORK2  I*4       (W)   WORK FOR RDFREE
C       RVAL    R*8       (W)   WORK FOR RDFREE
C       CVAL    A*80      (W)   WORK FOR RDFREE
C       INPMOD  A*1       (W)   WORK FOR RDFREE
C       ONOUTL  L*4       (O)   FLAG FOR OUTPUT LIST
C
C*******************************************************************

      use COMCMMC

      implicit none
 
      INTEGER*4       MAXVAL,IER,IPRINT,ICOL
      INTEGER*4       IVAL(MAXVAL),IWORK1(MAXVAL),IWORK2(MAXVAL)
      REAL*8          RVAL(MAXVAL)
      CHARACTER*(*)   LINE
      CHARACTER*80    CVAL(MAXVAL)
      CHARACTER*1     INPMOD(MAXVAL)
      LOGICAL*4       ONOUTL
 
      REAL*8          XCENT,YCENT,ZCENT,RADIN,RADOUT,DIST
      CHARACTER*8     ATTYSL
      INTEGER*4       IERSUB, IATM, JCOL
 
C*******************************************************************
 
      IER       = 0
 
C     * READ ATOM DATA *
 
      IF( INDEX(LINE(1:ICOL),'ATOM') .NE. 0 ) THEN
        INPMOD(1) = 'C'
        INPMOD(2) = 'C'
        INPMOD(3) = 'I'
        INPMOD(4) = 'R'
        INPMOD(5) = 'R'
        INPMOD(6) = 'C'
        INPMOD(7) = 'C'
        CALL  RDFREE( LINE,ICOL,MAXVAL,7,INPMOD,IWORK1,IWORK2,
     +                RVAL,IVAL,CVAL,IERSUB )
        IF( IERSUB .LT. 0 ) THEN
          WRITE(IPRINT,*) ' '
          WRITE(IPRINT,*) 'ERROR> INPVAR '
          WRITE(IPRINT,*) '       DATA TYPE ERROR '
          WRITE(IPRINT,*) ' '
          WRITE(IPRINT,'(A80)') LINE
          WRITE(IPRINT,*) ' '
          IER = IERSUB
          RETURN
        ENDIF
        IF( INDEX(CVAL(4),'YES') .NE. 0 ) ONOUTL = .TRUE.
        IATM   = IVAL(1)
        XCENT  = CORD(1,IATM)
        YCENT  = CORD(2,IATM)
        ZCENT  = CORD(3,IATM)
        RADIN  = RVAL(1)
        RADOUT = RVAL(2)
        ATTYSL = CVAL(3)(1:8)
      ENDIF
 
C     * READ COORDINATE DATA *
 
      IF( INDEX(LINE(1:ICOL),'COOR') .NE. 0 ) THEN
        INPMOD(1) = 'C'
        INPMOD(2) = 'C'
        INPMOD(3) = 'R'
        INPMOD(4) = 'R'
        INPMOD(5) = 'R'
        INPMOD(6) = 'R'
        INPMOD(7) = 'R'
        INPMOD(8) = 'C'
        INPMOD(9) = 'C'
        CALL  RDFREE( LINE,ICOL,MAXVAL,9,INPMOD,IWORK1,IWORK2,
     +                RVAL,IVAL,CVAL,IERSUB )
        IF( IERSUB .LT. 0 ) THEN
          WRITE(IPRINT,*) ' '
          WRITE(IPRINT,*) 'ERROR> INPVAR '
          WRITE(IPRINT,*) '       DATA TYPE ERROR '
          WRITE(IPRINT,*) ' '
          WRITE(IPRINT,'(A80)') LINE
          WRITE(IPRINT,*) ' '
          IER = IERSUB
          RETURN
        ENDIF
        IF( INDEX(CVAL(4),'YES') .NE. 0 ) ONOUTL = .TRUE.
        XCENT  = RVAL(1)
        YCENT  = RVAL(2)
        ZCENT  = RVAL(3)
        RADIN  = RVAL(4)
        RADOUT = RVAL(5)
        ATTYSL = CVAL(3)(1:8)
      ENDIF
 
C     * SEARCH ATOMS FREE/FIX *
 
      JCOL   = INDEX(ATTYSL,'*') - 1
      RADIN  = RADIN**2
      RADOUT = RADOUT**2
      DO 100 IATM = 1 , IXNATM
        DIST = ( CORD(1,IATM) - XCENT )**2 +
     +         ( CORD(2,IATM) - YCENT )**2 +
     +         ( CORD(3,IATM) - ZCENT )**2
        IF( ( DIST .GE. RADIN ) .AND. ( DIST .LE. RADOUT ) ) THEN
          IF( JCOL .EQ. -1 ) THEN
            IF( ATTYSL .EQ. CXATMN(IATM) ) THEN
              IF( CVAL(1)(1:3) .EQ. 'FIX'  ) IYTVAR(IATM) = 0
              IF( CVAL(1)(1:4) .EQ. 'FREE' ) IYTVAR(IATM) = 1
            ENDIF
          ENDIF
          IF( JCOL .EQ. 0 ) THEN
              IF( CVAL(1)(1:3) .EQ. 'FIX'  ) IYTVAR(IATM) = 0
              IF( CVAL(1)(1:4) .EQ. 'FREE' ) IYTVAR(IATM) = 1
          ENDIF
          IF( JCOL .GE. 1 ) THEN
            IF( INDEX(CXATMN(IATM),ATTYSL(1:JCOL)) .EQ. 1 ) THEN
              IF( CVAL(1)(1:3) .EQ. 'FIX'  ) IYTVAR(IATM) = 0
              IF( CVAL(1)(1:4) .EQ. 'FREE' ) IYTVAR(IATM) = 1
            ENDIF
          ENDIF
        ENDIF
 100  CONTINUE
 
      RETURN
      END
