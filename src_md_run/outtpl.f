                                                                        
C*******************************************************************    
C                                                                       
C     1 PROGRAM NAME                                                    
C          OUTTPL                                                       
C                                                                       
C     2 FUNCTION                                                        
C          OUTPUT FORMATTED TOPOLOGY FILE                               
C                                                                       
C     3 CALLING SEQUENCE                                                
C       3.1 INTERFACE                                                   
C          CALL  OUTTPL( IWRITE,IPRINT,FILENA,IER )                     
C       3.2 ARGUMENTS                                                   
C       IWRITE (I) I*4               : LOGICAL UNIT NUMBER FOR          
C                                      OUTPUT TOPOLOGY FILE             
C       IPRINT (I) I*4               : LOGICAL UNIT NUMBER FOR          
C                                      OUTPUT LOG                       
C       FILENA (I) A*80              : FILE NAME OF TOPOLOGY FILE       
C       IER    (O) I*4               : CONDITION CODE                   
C                         0   ;   NO ERROR                              
C                        -1   ;   FILE OPEN ERROR                       
C                                                                       
C     5 USED INCLUDE FILE                                               
C          COMBAS  COMERG  PHYCNS                                       
C                                                                       
C     6 USED SUBROUTINES                                                
C       6.1 OUTER       :   FLOPEN  FLCLOS                              
C                                                                       
C     7 AUTHOR                                                          
C       K.MORIKAMI                            1988.11.21 VERSION 1.0    
C                                             1989. 5.31 VERSION 2.0    
C                                             1990. 8. 8 VERSION 3.0    
C                                             1991. 5.28 VERSION 3.1    
C       T.NAKAI                               1991.11.20 VERSION 3.2    
C                                                                       
C*******************************************************************    
                                                                        
      SUBROUTINE  OUTTPL( IWRITE , IPRINT , FILENA , IER )              

      use COMBAS ; use COMERG ; use PHYCNS

      implicit none
                                                                        
      INTEGER*4     IWRITE,IPRINT,IER                                   
      CHARACTER*80  FILENA                                              
                                                                        
      INTEGER*4     ISTATM(ixmolc),IENATM(ixmolc)                       
      CHARACTER*80  LINE                                                
      INTEGER*4     IMOL, ITITL, NUMATM, IATM, I14N, ILOOP, J14NEN,
     +              K14N, J14NST, NUMBND, IBND, NUMANG, IANG, NLOOP,
     +              NUMTOR, ITOR, NUMIMP, IIMP, IFNC, ITYP, JTYP
C for inner coordinat
	  INTEGER*4     INCD(4),j

C*******************************************************************    
                                                                        
      IER = 0                                                           
      WRITE(IPRINT,*) ' '                                               
      WRITE(IPRINT,*) 'INFORMATION> OUTTPL (V3.0)'                      
      WRITE(IPRINT,*) '     WRITE FORMATTED TOPOLOGY FILE '             
      WRITE(IPRINT,*) ' '                                               
                                                                        
C     <<<  OPEN TOPOLOGY FILE  >>>                                      
                                                                        
      CALL  FLOPEN( IWRITE,FILENA,12,'NULL',0,IER )                     
      IF ( IER .NE. 0 ) GOTO 9000                                       
                                                                        
C     <<<  SEARCH TOP ATOM NUMBER OF FIRST CHAIN OF EACH MOLECULE  >>>  
                                                                        
      DO 1000 IMOL = 1 , IXMOLC                                         
        IF( IXTPCN(IMOL) .EQ. 1 ) THEN                                  
          ISTATM(IMOL) = 1                                              
        ELSE                                                            
          ISTATM(IMOL) = IXCEND( IXTPCN(IMOL)-1 ) + 1                   
        ENDIF                                                           
        IENATM(IMOL)   = IXCEND( IXTPCN(IMOL) )                         
 1000 CONTINUE                                                          
                                                                        
C     <<<  WRITE TOPOLOGY FILE  >>>                                     
                                                                        
C     1) WRITE TITLES                                                   
                                                                        
      WRITE(IWRITE,8100)                                                
      IF( IXTITL .EQ. 0 ) THEN                                          
        WRITE(IWRITE,8110)                                              
      ELSE                                                              
        DO 1100 ITITL = 1 , IXTITL                                      
          WRITE(IWRITE,'(A80)') CXTITL(ITITL)                           
 1100   CONTINUE                                                        
      ENDIF                                                             
                                                                        
 8100 FORMAT(' '/'TPL> TITLE')                                          
 8110 FORMAT(    '    NO-TITLE ')                                       
                                                                        
C     2) WRITE MOLECULES                                                
                                                                        
      WRITE(IWRITE,8200) IXMOLC                                         
      DO 1200 IMOL = 1 , IXMOLC                                         
        WRITE(IWRITE,8210) CXMOLC(IMOL),IXSQML(IMOL)                    
 1200 CONTINUE                                                          
                                                                        
 8200 FORMAT(' '/'TPL> MOLECULES '/';NUMBER OF MOLECULES : ',I5)        
 8210 FORMAT( A40,10X,I10 )                                             
                                                                        
C     3) WRITE ATOMS                                                    
                                                                        
      DO 1300 IMOL = 1 , IXMOLC                                         
                                                                        
        WRITE(IWRITE,8300) CXMOLC(IMOL) ,                               
     +                     ( IENATM(IMOL) - ISTATM(IMOL) + 1 )          
        NUMATM = 0                                                      
        DO 1310 IATM = ISTATM(IMOL) , IENATM(IMOL)                      
          NUMATM = NUMATM + 1                                           
          I14N = IX14IF(IATM,1) + IX14IF(IATM,2) + IX14IF(IATM,3)       
          IF( I14N .NE. 0 ) THEN                                        
             WRITE(IWRITE,8310) CXATMN(IATM),CXATMT(IATM),IXATYP(IATM), 
     +             CXRESN(IATM),IXARES(IATM),FXMASS(IATM),FXVDWR(IATM), 
     +             FXCHRG(IATM),                                        
     +             IX14IF(IATM,1),IX14IF(IATM,2),IX14IF(IATM,3),NUMATM  
            IF( MOD( I14N , 10 ) .EQ. 0 ) THEN                          
              NLOOP = I14N / 10                                         
            ELSE                                                        
              NLOOP = ( I14N/10 ) + 1                                   
            ENDIF                                                       
            DO 1311 ILOOP = 1 , NLOOP                                   
              LINE  = ' '                                               
              J14NST = 10*ILOOP - 9                                     
              IF( ILOOP .EQ. NLOOP ) THEN                               
                J14NEN      = I14N                                      
              ELSE                                                      
                J14NEN      = 10*ILOOP                                  
C               LINE(72:73) = '->'                                      
              ENDIF                                                     
              LINE(72:73) = '->' 
              WRITE(LINE(1:70),8320)                                    
     +       ( ( IX14LT(IATM,K14N) - IATM ) , K14N = J14NST , J14NEN )  
              WRITE(IWRITE,'(A80)') LINE                                
 1311       CONTINUE                                                    
          ELSE                                                          
C           WRITE(IWRITE,8330) CXATMN(IATM),CXATMT(IATM),IXATYP(IATM),  
            WRITE(IWRITE,8310) CXATMN(IATM),CXATMT(IATM),IXATYP(IATM),  
     +            CXRESN(IATM),IXARES(IATM),FXMASS(IATM),FXVDWR(IATM),  
     +            FXCHRG(IATM),                                         
     +            IX14IF(IATM,1),IX14IF(IATM,2),IX14IF(IATM,3),NUMATM   
          ENDIF                                                         
C
C  OUTPUT INFORMATION FOR INNER COORDINATES
C                           1992. 3.21 T.NAKAI
          DO 1320 J=1,4
            IF(IXINCD(IATM,J).NE.0) THEN
              INCD(J)=IXINCD(IATM,J) - IATM
            ELSE
              INCD(J)=0
            END IF
 1320     CONTINUE
          WRITE(IWRITE,8340) (INCD(J),J=1,4),(FXINCD(IATM,J),J=1,3)
C
 1310   CONTINUE                                                        
 1300 CONTINUE                                                          
                                                                        
 8300 FORMAT(' '/'TPL> ATOMS '/A40/'; NUMBER OF ATOMS = ',I10)          
 8310 FORMAT(A8,A4,I5,1X,A8,I5,3F10.4,3I3,' -> ;',I5)                   
 8330 FORMAT(A8,A4,I5,1X,A8,I5,3F10.4,3I3,'    ;',I5)                   
 8320 FORMAT(10I7)                                                      
 8340 FORMAT(4I7,3X,3F10.4)                                             
                                                                        
C     4) WRITE BONDS                                                    
                                                                        
      DO 1400 IMOL = 1 , IXMOLC                                         
        NUMBND = 0                                                      
        DO 1410 IBND = 1 , IYNBND                                       
          IF( ( IXAMOL(IYPBND(1,IBND)) .EQ. IMOL ) .AND.                
     +        ( IXACHN(IYPBND(1,IBND)) .EQ. IXTPCN(IMOL) )      ) THEN  
            NUMBND = NUMBND + 1                                         
          ENDIF                                                         
 1410   CONTINUE                                                        
        IF ( NUMBND .NE. 0 ) THEN                                       
          WRITE(IWRITE,8400) CXMOLC(IMOL),NUMBND                        
          NUMBND = 0                                                    
          DO 1420 IBND = 1 , IYNBND                                     
          IF( ( IXAMOL(IYPBND(1,IBND)) .EQ. IMOL ) .AND.                
     +        ( IXACHN(IYPBND(1,IBND)) .EQ. IXTPCN(IMOL) )      ) THEN  
              NUMBND = NUMBND + 1                                       
              WRITE(IWRITE,8410)                                        
     +                ( IYPBND(1,IBND) - ISTATM(IMOL) + 1 ) ,           
     +                ( IYPBND(2,IBND) - ISTATM(IMOL) + 1 ) ,           
     +                FYFBND(IBND),FYQBND(IBND),NUMBND                  
            ENDIF                                                       
 1420     CONTINUE                                                      
        ENDIF                                                           
 1400 CONTINUE                                                          
                                                                        
 8400 FORMAT(' '/'TPL> BONDS '/A40/'; NUMBER OF BONDS = ',I10)          
 8410 FORMAT(2I10,2F15.7,' ; ',I10)                                     
                                                                        
C     5) WRITE ANGLES                                                   
                                                                        
      DO 1500 IMOL = 1 , IXMOLC                                         
        NUMANG = 0                                                      
        DO 1510 IANG = 1 , IYNANG                                       
          IF( ( IXAMOL(IYPANG(1,IANG)) .EQ. IMOL ) .AND.                
     +        ( IXACHN(IYPANG(1,IANG)) .EQ. IXTPCN(IMOL) )      ) THEN  
            NUMANG = NUMANG + 1                                         
          ENDIF                                                         
 1510   CONTINUE                                                        
        IF ( NUMANG .NE. 0 ) THEN                                       
          WRITE(IWRITE,8500) CXMOLC(IMOL),NUMANG                        
          NUMANG = 0                                                    
          DO 1520 IANG = 1 , IYNANG                                     
            IF( ( IXAMOL(IYPANG(1,IANG)) .EQ. IMOL ) .AND.              
     +          ( IXACHN(IYPANG(1,IANG)) .EQ. IXTPCN(IMOL) )  ) THEN    
              NUMANG = NUMANG + 1                                       
              WRITE(IWRITE,8510)                                        
     +                 ( IYPANG(1,IANG) - ISTATM(IMOL) + 1 ) ,          
     +                 ( IYPANG(2,IANG) - ISTATM(IMOL) + 1 ) ,          
     +                 ( IYPANG(3,IANG) - ISTATM(IMOL) + 1 ) ,          
     +                 FYFANG(IANG),(FYQANG(IANG)/RAD),                 
     +                 NUMANG                                           
            ENDIF                                                       
 1520     CONTINUE                                                      
        ENDIF                                                           
 1500 CONTINUE                                                          
                                                                        
 8500 FORMAT(' '/'TPL> ANGLES '/A40/'; NUMBER OF ANGLES = ',I10)        
 8510 FORMAT(3I10,2F15.7,' ; ',I10)                                     
                                                                        
C     6) WRITE TORSIONS                                                 
                                                                        
      DO 1600 IMOL = 1 , IXMOLC                                         
        NUMTOR = 0                                                      
        DO 1610 ITOR = 1 , IYNTOR                                       
          IF( ( IXAMOL(IYPTOR(1,ITOR)) .EQ. IMOL ) .AND.                
     +        ( IXACHN(IYPTOR(1,ITOR)) .EQ. IXTPCN(IMOL) )      ) THEN  
            NUMTOR = NUMTOR + 1                                         
          ENDIF                                                         
 1610   CONTINUE                                                        
        IF ( NUMTOR .NE. 0 ) THEN                                       
          WRITE(IWRITE,8600) CXMOLC(IMOL),NUMTOR                        
          NUMTOR = 0                                                    
          DO 1620 ITOR = 1 , IYNTOR                                     
            IF( ( IXAMOL(IYPTOR(1,ITOR)) .EQ. IMOL ) .AND.              
     +          ( IXACHN(IYPTOR(1,ITOR)) .EQ. IXTPCN(IMOL) ) ) THEN     
              NUMTOR = NUMTOR + 1                                       
              WRITE(IWRITE,8610)                                        
     +                ( IYPTOR(1,ITOR) - ISTATM(IMOL) + 1 ) ,           
     +                ( IYPTOR(2,ITOR) - ISTATM(IMOL) + 1 ) ,           
     +                ( IYPTOR(3,ITOR) - ISTATM(IMOL) + 1 ) ,           
     +                ( IYPTOR(4,ITOR) - ISTATM(IMOL) + 1 ) ,           
     +                FYFTOR(ITOR),IYTDIV(ITOR),                        
     +                INT(FYTROT(ITOR)),( FYTPHS(ITOR)/RAD ),           
     +                IYTNBF(ITOR),NUMTOR                               
            ENDIF                                                       
 1620     CONTINUE                                                      
        ENDIF                                                           
 1600 CONTINUE                                                          
                                                                        
 8600 FORMAT(' '/'TPL> TORSIONS '/A40/'; NUMBER OF TORSIONS = ',I10)    
 8610 FORMAT(4I8,F10.4,2I5,F10.4,I5,' ; ',I6 )                          
                                                                        
C     7) WRITE IMPROPER-TORSIONS                                        
C                                                                       
      DO 1700 IMOL = 1 , IXMOLC                                         
        NUMIMP = 0                                                      
        DO 1710 IIMP = 1 , IYNIMP                                       
          IF( ( IXAMOL(IYPIMP(1,IIMP)) .EQ. IMOL ) .AND.                
     +        ( IXACHN(IYPIMP(1,IIMP)) .EQ. IXTPCN(IMOL) )      ) THEN  
            NUMIMP = NUMIMP + 1                                         
          ENDIF                                                         
 1710   CONTINUE                                                        
        IF ( NUMIMP .NE. 0 ) THEN                                       
          WRITE(IWRITE,8700) CXMOLC(IMOL),NUMIMP                        
          NUMIMP = 0                                                    
          DO 1720 IIMP = 1 , IYNIMP                                     
            IF( ( IXAMOL(IYPIMP(1,IIMP)) .EQ. IMOL ) .AND.              
     +          ( IXACHN(IYPIMP(1,IIMP)) .EQ. IXTPCN(IMOL) )     ) THEN 
              NUMIMP = NUMIMP + 1                                       
              WRITE(IWRITE,8610)                                        
     +              ( IYPIMP(1,IIMP) - ISTATM(IMOL) + 1 ) ,             
     +              ( IYPIMP(2,IIMP) - ISTATM(IMOL) + 1 ) ,             
     +              ( IYPIMP(3,IIMP) - ISTATM(IMOL) + 1 ) ,             
     +              ( IYPIMP(4,IIMP) - ISTATM(IMOL) + 1 ) ,             
     +              FYFIMP(IIMP),IYIDIV(IIMP),                          
     +              INT(FYIROT(IIMP)),( FYIPHS(IIMP)/RAD ),             
     +              IYINBF(IIMP),NUMIMP                                 
            ENDIF                                                       
 1720     CONTINUE                                                      
        ENDIF                                                           
 1700 CONTINUE                                                          
                                                                        
 8700 FORMAT(' '/'TPL> IMPROPER-TORSIONS '/A40/                         
     +             '; NUMBER OF IMPROPER-TORSIONS = ',I10)              
                                                                        
C     8) WRITE FUNCTIONS                                                
                                                                        
      WRITE(IWRITE,8800)                                                
      DO 1800 IFNC = 1 , IYNBPF                                         
        WRITE(IWRITE,8810) IFNC,IYNPAR(IFNC),CYNBPF(IFNC)               
 1800 CONTINUE                                                          
                                                                        
 8800 FORMAT(' '/'TPL> FUNCTIONS ')                                     
 8810 FORMAT(2I10,10X,A40)                                              
                                                                        
C     9) WRITE NONBONDS                                                 
                                                                        
      WRITE(IWRITE,8900) IYNTYP                                         
                                                                        
C     9-1) WRITE AMBER-LIKE NONBONDED PARAMETERS                        
                                                                        
      IF( INDEX(CYNBPF(1),'AMBER') .NE. 0 ) THEN                        
        DO 1900 ITYP = 1 , IYNTYP                                       
          WRITE(IWRITE,8910) ITYP,FYVWRD(ITYP),FYVWME(ITYP),            
     +                       FY14SE(ITYP),FY14SV(ITYP)                  
 1900   CONTINUE                                                        
        DO 1910 ITYP = 1 , IYNTYP                                       
          DO 1920 JTYP = 1 , IYNTYP                                     
            IF( IYPPID(ITYP,JTYP) .EQ. 2 ) THEN                         
              WRITE(IWRITE,8920) ITYP,JTYP,IYPPID(ITYP,JTYP),           
     +                           FYNBPP(4,ITYP,JTYP),                   
     +                           FYNBPP(3,ITYP,JTYP)                    
            ENDIF                                                       
 1920     CONTINUE                                                      
 1910   CONTINUE                                                        
      ENDIF                                                             
                                                                        
C     9-2) WRITE OPLS-LIKE NONBONDED PARAMETERS                         
                                                                        
      IF( INDEX(CYNBPF(1),'OPLS') .NE. 0 ) THEN                         
        DO 1930 ITYP = 1 , IYNTYP                                       
          WRITE(IWRITE,8910) ITYP,FYVWRD(ITYP),FYVWME(ITYP),            
     +                       FY14SE(ITYP),FY14SV(ITYP)                  
 1930   CONTINUE                                                        
      ENDIF                                                             
                                                                        
 8900 FORMAT(' '/'TPL> NONBONDS '/'; MAXINUM TYPE NUMBER = ',I10)       
 8910 FORMAT(I5,'    0','   1',4F15.7)                                  
 8920 FORMAT(3I5,2F15.7)                                                
                                                                        
C     <<<  CLOSE FILE  >>>                                              
                                                                        
      CALL  FLCLOS( IWRITE , 10 , IER )                                 
      IF( IER .NE. 0 ) GOTO 9000                                        
                                                                        
 9999 RETURN                                                            
                                                                        
 9000 CONTINUE                                                          
        WRITE(IPRINT,*) 'ERROR> OUTTPL '                                
        WRITE(IPRINT,*) '    FILE OPEN ERROR OR CLOSE ERROR '           
        WRITE(IPRINT,*) ' '                                             
        IER = -1                                                        
        GOTO 9999                                                       
      END                                                               
