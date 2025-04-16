
      SUBROUTINE DEF2M (A,N,NA,NE,NV,EPS,IOPT,E,V,IFLG,WK,IWK,IER)
C
C------------------------------------------------------------------------------
C
C     NAME                - DEF2M : DOUBLE PRECISION
C
C     USAGE               - CALL DEF2M(A,N,NA,NV,EPS,IOPT,E,V,IFLG,
C                                                        WK,IWK,IER)
C
C     FUNCTION            - BY THE HOUSEHOLDER-BISECTION METHOD AND THE
C                             INVERSE ITERATION METHOD, WE OBTAIN THE  
C                             EIGENVALUES AND EIGENVECTORS OF A REAL SYM-
C                             METRIC MATRIX.
C
C     ARGUMENTS  A(NA,N)  - INPUT. REAL SYMMETRIC MATRIX.
C                N        - INPUT. ORDER N OF THE MATRIX.
C                             (NE<=N<=NA)
C                NA       - INPUT. NUMBER OF THE ROWS OF A, V IN THE
C                             DIMENSION STATEMENT OF A MAIN PROGRAM.
C                NE       - INPUT. NUMBER OF EIGENVALUES TO BE OBTAINED. 
C                             (1<=NV<=NE)
C                NV       - INPUT. NUMBER OF EIGENVECTORS TO BE OBTAINED.
C                EPS      - INPUT. CRITERION FOR CONVERGENCE.
C                             WHEN EPS<0.0 IS GIVEN,EPS IS SET TO    
C                             THE STANDARD VALUE.
C                IOPT     - INPUT. OPTION PARAMETER.
C                             IF IOPT=1, EIGENVALUES AND EIGENVECTORS ARE
C                                        COMPUTED IN DESCENDING ORDER.
C                             IF IOPT=2, EIGENVALUES AND EIGENVECTORS ARE
C                                        COMPUTED IN ASCENDING ORDER.
C                             IF IOPT=3, EIGENVALUES ARE COMPUTED IN
C                                        DESCENDING ORDER.
C                             IF IOPT=4, EIGENVALUES ARE COMPUTED IN
C                                        ASCENDING ORDER.
C                E(NE)    - OUTPUT. THE K-TH EIGENVALUE ENTERS INTO E(K)
C                             (K=1,2,....,NE).
C                V(NA,NV) - OUTPUT. THE K-TH COLUMN BECOMES THE EIGEN-
C                             VECTOR CORRESPONDING TO THE EIGENVALUE
C                             E(K) (K=1,2,...,NV).
C                IFLG(NV) - OUTPUT. 
C                             IF IFLG(K)=1, THE K-TH EIGENVECTOR WAS
C                                           OBTAIND.
C                             IF IFLG(K)=0, THE K-TH EIGENVECTOR DID NOT
C                                           CONVERGE.
C                WK(N,6)  - WORK AREA.
C                IWK(N)   - WORK AREA.
C                IER      - ERROR INDICATOR.
C                             IF IER=   0, NO ERROR WAS DETECTED.
C                             IF IER=1000, THE INVERSE ITERATION DID NOT
C                                          CONVERGE. IFLG WAS SET TO ZERO.
C                             IF IER=2000, PARAMETER ERROR.
C                                          (NV<1,NV>NE,NE>N,N>NA,
C                                           IOPT<1 OR IOPT>4).
C
C     REQD. SUBPROGRAM    - DET2M,DETAM,DETBM,SUERM,SURTM
C
C     STATUS              - S-1511-1 05-02
C
C     HISTORY             - DATE.    1980. 3
C                                    1980.11
C                                    1981.10
C                                    1986. 6
C                                    1991. 3 JST J.Takeuchi
C
C------------------------------------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
C     GENERIC
      DIMENSION A(NA,N),E(NE),V(NA,NV),IFLG(NV),WK(N,6),NAME(2),IWK(N)
      DATA ZERO,ONE,PT5,UEPS/0.0D0,1.0D0,0.5D0,1.0D-12/
C      DATA ZERO,ONE,PT5,UEPS/0.0D0,1.0D0,0.5D0,Z990/  1991. 3/4 JST J.Takeuchi
C                                               ^^^^Kishu Izon Teisuu zyokyo
C      DATA XMIN /Z341/  1991. 3/4 JST J.Takeuchi
C                 ^^^^Kishu Izon Teisuu zyokyo
      DATA NAME(1),NAME(2)/'DEF','2M  '/

D     integer loop,subloop

C
C     CHECK OF THE PARAMETERS
C

D     do loop = 1,N
D        do subloop = 1,NA
D           write(6,*)'A(',subloop,',',loop,')=',A(subloop,loop)
D        end do
D     end do

      IF (NV.LT.1 .OR. NV.GT.NE)    GO TO 9000
      IF (N.LT.NE .OR. N.GT.NA)     GO TO 9000
      IF (IOPT.LT.1 .OR. IOPT.GT.4) GO TO 9000
      IER=0
      NM1=N-1
      NM2=N-2
C      TOL=XMIN/UEPS
C      Seido hensuu no settei wo teisuu ni sita.
      TOL=1.0D-12
C
C     REDUCE TO TRIDIAGONAL FORM BY HOUSEHOLDER'S METHOD
C
      IF (NM2) 10,130,20
C
C     WHEN N=1
C
   10 E(1)=A(1,1)
      IF (EPS.LT.ZERO) EPS=ABS(A(1,1))*UEPS
      IF (IOPT.GE.3) GO TO 9999
      V(1,1)=ONE
      IFLG(1)=1
      GO TO 9999
C
C
C     WHEN N=3,4,...
C
   20 DO 120 K=1,NM2
         K1=K+1
         S=ZERO
         DO 30 I=K1,N
            S=S+A(I,K)*A(I,K)
   30    CONTINUE
         WK(K,1)=ZERO
         IF (S.LE.TOL) GO TO 110
         SR=SQRT(S)
         A1=A(K1,K)
         IF (A1.LT.ZERO) SR=-SR
         WK(K,1)=-SR
         R=ONE/(S+A1*SR)
         A(K1,K)=A1+SR
         DO 40 I=K1,N
            WK(I,1)=ZERO
   40    CONTINUE
         DO 60 J=K1,NM1
            T=A(J,K)
            S=A(J,J)*T
            T=T*R
            JP1=J+1
            DO 50 I=JP1,N
               S=S+A(I,J)*A(I,K)
               WK(I,1)=WK(I,1)+T*A(I,J)
   50       CONTINUE
            WK(J,1)=WK(J,1)+S*R
   60    CONTINUE
         WK(N,1)=WK(N,1)+A(N,N)*A(N,K)*R
         S=ZERO
         DO 70 I=K1,N
            S=S+A(I,K)*WK(I,1)
   70    CONTINUE
         T=-PT5*S*R
         DO 80 I=K1,N
            WK(I,1)=WK(I,1)+T*A(I,K)
   80    CONTINUE
         DO 100 J=K1,N
            WJ1=-WK(J,1)
            AJK=-A(J,K)
            DO 90 I=J,N
               A(I,J)=A(I,J)+WJ1*A(I,K)
               A(I,J)=A(I,J)+AJK*WK(I,1)
   90       CONTINUE
  100    CONTINUE
  110    WK(K,2)=A(K,K)
  120 CONTINUE
  130 WK(NM1,2)=A(NM1,NM1)
      WK(N,2)=A(N,N)
      WK(NM1,1)=A(N,NM1)
C
C     CALL THE ROUTINE TO COMPUTE EIGENVALUES AND EIGENVECTORS FOR
C        TRIDIAGONAL MATRIX
C
      CALL DET2M(WK(1,2),N,WK(1,1),NE,NV,EPS,IOPT,E,V,NA,IFLG,
     *            WK(1,3),IWK,IER)
      IF (IOPT.GE.3) GO TO 9999
C
C     BACKWARD TRANSFORMATION
C
      IF (N.EQ.2) GO TO 9999
      DO 140 I=1,NM2
         WK(I,1)=-WK(I,1)*A(I+1,I)
  140 CONTINUE
      DO 180 JJ=1,NV
         DO 170 K=1,NM2
            J=NM1-K
            IF (WK(J,1).EQ.ZERO) GO TO 170
            S=ZERO
            JP1=J+1
            DO 150 I=JP1,N
               S=S+A(I,J)*V(I,JJ)
  150       CONTINUE
            IF (S.EQ.ZERO) GO TO 170
            T=-S/WK(J,1)
            DO 160 I=JP1,N
               V(I,JJ)=V(I,JJ)+T*A(I,J)
  160       CONTINUE
  170    CONTINUE
  180 CONTINUE
C
C     NORMAL END
C

D     do loop = 1,NE
D        write(6,*)'E(',loop,')=',E(loop)
D        do subloop = 1,NA
D           write(6,*)' V(',subloop,',',loop,')=',V(subloop,loop)
D        end do
D     end do

      GO TO 9999
C
C     PARAMETER ERROR
C
 9000 CALL SUERM(NAME,2000,IER)
 9999 CALL SURTM(NAME,IER)
      RETURN
      END
                 
       


      SUBROUTINE DET2M (D,N,D1,NE,NV,EPS,IOPT,E,V,NV1,IFLG,WK,IWK,IER)
C
C------------------------------------------------------------------------------
C
C     NAME                - DET2M : DOUBLE PRECISION
C
C     USAGE               - CALL DET2M(D,N,D1,NE,NV,EPS,IOPT,E,V,NV1,
C                                                      IFLG,WK,IWK,IER)
C
C     FUNCTION            - BY THE BISECTION METHOD AND INVERSE ITERATION
C                             METHOD, WE OBTAIN THE EIGENVALUES AND
C                             EIGENVECTORS OF A REAL SYMMETRIC TRAIDIAGONAL
C                             MATRIX.
C
C     ARGUMENTS  D(N)     - INPUT. DIAGONAL ELEMENTS OF A REAL SYMMETRIC
C                             TRIDIAGONAL MATRIX.
C                N        - INPUT. ORDER N OF THE MATRIX.
C                             (NE<=N<=NV1)
C                D1(N)    - INPUT. SUB-DIAGONAL ELEMENTS.
C                             THE K-TH SUB-DIAGONAL ELEMENT IS GIVEN TO
C                             D1(K) (K=1,2,...,N-1).
C                NE       - INPUT. NUMBER OF EIGENVALUES TO BE OBTAINED. 
C                             (1<=NV<=NE)
C                NV       - INPUT. NUMBER OF EIGENVECTORS TO BE OBTAINED.
C                EPS      - INPUT. CRITERION FOR CONVERGENCE OF EIGEN-
C                             VALUE.
C                             WHEN EPS<0.0 IS GIVEN,EPS IS SET TO    
C                             THE STANDARD VALUE.
C                IOPT     - INPUT. OPTION PARAMETER.
C                             IF IOPT=1, EIGENVALUES AND EIGENVECTORS ARE
C                                        COMPUTED IN DESCENDING ORDER.
C                             IF IOPT=2, EIGENVALUES AND EIGENVECTORS ARE
C                                        COMPUTED IN ASCENDING ORDER.
C                             IF IOPT=3, EIGENVALUES ARE COMPUTED IN
C                                        DESCENDING ORDER.
C                             IF IOPT=4, EIGENVALUES ARE COMPUTED IN
C                                        ASCENDING ORDER.
C                E(NE)    - OUTPUT. THE K-TH EIGENVALUE ENTERS INTO E(K)
C                             (K=1,2,....,NE).
C                V(NV1,NV)- OUTPUT. THE K-TH COLUMN BECOMES THE EIGEN-
C                             VECTOR CORRESPONDING TO THE EIGENVALUE
C                             E(K) (K=1,2,...,NV).
C                NV1      - INPUT. NUMBER OF THE ROWS OF V IN THE
C                             DIMENSION STATEMENT OF A MAIN PROGRAM.    
C                IFLG(NV) - OUTPUT. 
C                             IF IFLG(K)=1, THE K-TH EIGENVECTOR WAS
C                                           OBTAIND.
C                             IF IFLG(K)=0, THE K-TH EIGENVECTOR DID NOT
C                                           CONVERGE.
C                WK(N,4)  - WORK AREA.
C                IWK(N)   - WORK AREA.
C                IER      - ERROR INDICATOR.
C                             IF IER=   0, NO ERROR WAS DETECTED.
C                             IF IER=1000, THE INVERSE ITERATION DID NOT
C                                          CONVERGE. IFLG WAS SET TO ZERO.
C                             IF IER=2000, PARAMETER ERROR.
C                                          (NV<1,NV>NE,NE>N,N>NV1,
C                                           IOPT<1 OR IOPT>4).
C
C     REQD. SUBPROGRAM    - DETAM,DETBM,SUERM,SURTM
C
C     STATUS              - S-1511-1 05-02
C
C     HISTORY             - DATE.    1980. 3
C                                    1980.11
C                                    1981.10
C                                    1986. 6
C
C------------------------------------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
C     GENERIC
      DIMENSION D(N),D1(N),E(NE),V(NV1,NV),IFLG(NV),WK(N,4),
     *          IWK(N),NAME(2)
      DATA NAME(1),NAME(2)/'DET','2M  '/
C
C     CHECK OF THE PARAMETERS
C
      IF (NV.LT.1 .OR. NV.GT.NE)    GO TO 9000
      IF (N.LT.NE .OR. N.GT.NV1)    GO TO 9000
      IF (IOPT.LT.1 .OR. IOPT.GT.4) GO TO 9000
      IER=0
      IF (IOPT.EQ.1 .OR. IOPT.EQ.3) MOPT=1
      IF (IOPT.EQ.2 .OR. IOPT.EQ.4) MOPT=2
      EPS1=EPS
C
C     CALL THE ROUTINE TO COMPUTE EIGENVALUES
C
      CALL DETAM(D,N,D1,NE,EPS,MOPT,E,WK,IER)
      IF (IOPT.GT.2) GO TO 9999
C
C     CALL THE ROUTINE TO COMPUTE EIGENVECTORS
C
      CALL DETBM(D,N,D1,E,NV,EPS1,V,NV1,IFLG,WK(1,1),WK(1,2),WK(1,3),
     *           WK(1,4),IWK,IER)
C
C     NORMAL END
C
      GO TO 9999
C
C     PARAMETER ERROR
C
 9000 CALL SUERM(NAME,2000,IER)
 9999 CALL SURTM(NAME,IER)
      RETURN
      END
                 
       


      SUBROUTINE DETAM (D,N,D1,NE,EPS,IOPT,E,WK,IER)
C
C------------------------------------------------------------------------------
C
C     NAME                - DETAM : DOUBLE PRECISION
C
C     USAGE               - CALL DETAM(D,N,D1,NE,EPS,IOPT,E,WK,IER)
C
C     FUNCTION            - BY THE BISECTION METHOD WE OBTAIN THE EIGEN-
C                             VALUES OF A REAL SYMMETRIC TRAIDIAGONAL
C                             MATRIX.
C
C     ARGUMENTS  D(N)     - INPUT. DIAGONAL ELEMENTS OF A REAL SYMMETRIC
C                             TRIDIAGONAL MATRIX.
C                N        - INPUT. ORDER N OF THE MATRIX.
C                             (1<=NE<=N)
C                D1(N)    - INPUT. SUB-DIAGONAL ELEMENTS.
C                             THE K-TH SUB-DIAGONAL ELEMENT IS GIVEN TO
C                             D1(K) (K=1,2,...,N-1).
C                NE       - INPUT. NUMBER OF EIGENVALUES TO BE OBTAINED. 
C                EPS      - INPUT. CRITERION FOR CONVERGENCE OF EIGEN-
C                             VALUE.
C                             WHEN EPS<0.0 IS GIVEN,EPS IS SET TO    
C                             THE STANDARD VALUE.
C                IOPT     - INPUT. OPTION PARAMETER.
C                             IF IOPT=1, EIGENVALUES ARE COMPUTED IN
C                                        DESCENDING ORDER.
C                             IF IOPT=2, EIGENVALUES ARE COMPUTED IN
C                                        ASCENDING ORDER.
C                E(NE)    - OUTPUT. THE K-TH EIGENVALUE ENTERS INTO E(K)
C                             (K=1,2,....,NE).
C                WK(N)    - WORK AREA.
C                IER      - ERROR INDICATOR.
C                             IF IER=   0, NO ERROR WAS DETECTED.
C                             IF IER=2000, PARAMETER ERROR.
C                                          (NE<1,NE>N,N<1,
C                                           IOPT<1 OR IOPT>2).
C
C     REQD. SUBPROGRAM    - SUERM,SURTM
C
C     STATUS              - S-1511-1 05-02
C
C     HISTORY             - DATE.    1980. 3
C                                    1980.11
C                                    1983.12
C                                    1986. 6
C
C------------------------------------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
C     GENERIC
      DIMENSION D(N),D1(N),E(NE),WK(N),NAME(2)
      DATA ZERO,TWO,PT5,UEPS/0.0D0,2.0D0,0.5D0,1.0D-12/
C      DATA ZERO,TWO,PT5,UEPS/0.0D0,2.0D0,0.5D0,Z990/1991. 3/4 JST J.Takeuchi
C                                               ^^^^kishu izon no hensuu jokyo
      DATA ONE/1.0D0/
      DATA NAME(1),NAME(2)/'DET','AM  '/
C
C     CHECK OF THE PARAMETERS
C
      IF (NE.LT.1 .OR. NE.GT.N)     GO TO 9000
      IF (IOPT.LT.1 .OR. IOPT.GT.2) GO TO 9000
      IER=0
      IF (N.NE.1) GO TO 10
      E(1)=D(1)
      R=ABS(D(1))
      IF (EPS.LT.ZERO) EPS=R*UEPS
      GO TO 9999
C
   10 NM1=N-1
      XH=MAX(D(1)+ABS(D1(1)),D(N)+ABS(D1(NM1)))
      XL=MIN(D(1)-ABS(D1(1)),D(N)-ABS(D1(NM1)))
      IF (N.EQ.2) GO TO 30
      DO 20 I=2,NM1
         T=ABS(D1(I-1))+ABS(D1(I))
         IF (D(I)+T.GT.XH) XH=D(I)+T
         IF (D(I)-T.LT.XL) XL=D(I)-T
   20 CONTINUE
   30 R=MAX(ABS(XH),ABS(XL))
      T=R*UEPS
      XH=XH+T
      XL=XL-T
      IF (EPS.LT.ZERO) EPS=R*UEPS
      DO 40 I=1,NM1
         WK(I)=D1(I)*D1(I)
   40 CONTINUE
      F=XH
      T=XL
      IF (IOPT.EQ.1) GO TO 45
      F=XL
      T=XH
   45 DO 50 I=1,NE
         E(I)=T
   50 CONTINUE
      UEPSR=ONE/UEPS
      DO 120 K=1,NE
         DD=E(K)
   60    T=PT5*(DD+F)
         TFE=TWO*(ABS(DD)+ABS(F))*UEPS+EPS
         IF (ABS(DD-F).LE.TFE) GO TO 115
         J=0
         Q=D(1)-T
         IF (Q.GE.ZERO) J=J+1
         DO 90 I=2,N
            IF (Q.NE.ZERO) GO TO 70
            QQ=ABS(D1(I-1))*UEPSR
            GO TO 80
C
   70       QQ=WK(I-1)/Q
   80       Q=D(I)-T-QQ
C            IF (Q.GE.ZERO) J=N-J 1991. 3/4 JST J.Takeuchi
C            Gen'zaino tiiki ni aru koyuuti no kazu no ruikei wo okonatteiru
C           konotame, koyuuti no kazu (hensuu J) no kasan to sita.
            IF (Q.GE.ZERO) J=J+1
   90    CONTINUE
         IF (IOPT.EQ.2) J=N-J
         IF (J.GE.K) GO TO 100
         F=T
         GO TO 60
C
  100    DD=T
         M=MIN(J,NE)
         DO 110 I=K,M
            E(I)=T
  110    CONTINUE
         GO TO 60
C
  115    E(K)=T
  120 CONTINUE   
C
C     NORMAL END
C
      GO TO 9999
C
C     PARAMETER ERROR
C
 9000 CALL SUERM(NAME,2000,IER)
 9999 CALL SURTM(NAME,IER)
      RETURN
      END

      SUBROUTINE DETBM (D,N,D1,E,NV,EPS,V,NV1,IFLG,WK1,WK2,WK3,WK4,
     *                  IWK,IER)
C
C------------------------------------------------------------------------------
C
C     NAME                - DETBM : DOUBLE PRECISION
C
C     USAGE               - CALL DETBM(D,N,D1,E,NV,EPS,V,NV1,IFLG,
C                                        WK1,WK2,WK3,WK4,IWK,IER)
C
C     FUNCTION            - BY THE INVERSE ITERATION METHOD, WE OBTAIN THE
C                             EIGENVECTORS OF A REAL SYMMETRIC TRAIDIAGONAL
C                             MATRIX.
C
C     ARGUMENTS  D(N)     - INPUT. DIAGONAL ELEMENTS OF A REAL SYMMETRIC
C                             TRIDIAGONAL MATRIX.
C                N        - INPUT. ORDER N OF THE MATRIX.
C                             (1<=NV<=N<=NV1)
C                D1(N)    - INPUT. SUB-DIAGONAL ELEMENTS.
C                             THE K-TH SUB-DIAGONAL ELEMENT IS GIVEN TO
C                             D1(K) (K=1,2,...,N-1).
C                E(NV)    - INPUT. THE K-TH EIGENVALUE IS GIVEN TO E(K)
C                             (K=1,2,...,NV).
C                NV       - INPUT. NUMBER OF EIGENVECTORS TO BE OBTAINED. 
C                EPS      - INPUT. CRITERION FOR CONVERGENCE.
C                V(NV1,NV)- OUTPUT. THE K-TH COLUMN BECOMES THE EIGEN-
C                             VECTOR CORRESPONDING TO THE EIGENVALUE
C                             E(K) (K=1,2,...,NV)
C                NV1      - INPUT. NUMBER OF THE ROWS OF V IN THE
C                             DIMENSION STATMENT OF A MAIN PROGRAM.
C                IFLG(NV) - OUTPUT.
C                             IF IFLG(K)=1, THE K-TH EIGENVECTOR WAS
C                                           OBTAINED.
C                             IF IFLG(K)=0, THE K-TH EIGENVECTOR DID NOT
C                                           CONVERGE.
C                WK1(N)   - WORK AREA.
C                WK2(N)   - WORK AREA.
C                WK3(N)   - WORK AREA.
C                WK4(N)   - WORK AREA.
C                IWK(N)   - WORK AREA.                
C                IER      - ERROR INDICATOR.
C                             IF IER=   0, NO ERROR WAS DETECTED.
C                             IF IER=1000, THE INVERSE ITERATION DID NOT
C                                          CONVERGE. IFLG WAS SET TO ZERO.
C                             IF IER=2000, PARAMETER ERROR.
C                                          (NV<1,NV>N.N>1 OR N>NV1).
C
C     REQD. SUBPROGRAM    - SUERM,SURTM
C
C     STATUS              - S-1511-1 05-02
C
C     HISTORY             - DATE.    1980. 3
C                                    1980.11
C                                    1981.10
C                                    1983.12
C                                    1985. 1
C                                    1986. 6
C
C------------------------------------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
C     GENERIC
      DIMENSION D(N),D1(N),E(NV),V(NV1,NV),WK1(N),WK2(N),WK3(N),WK4(N),
     *          IWK(N),IFLG(NV),NAME(2)
      DATA ZERO,ONE,DCM3/0.0D0,1.0D0,1.0D-3/
      DATA UEPS/1.0D-12/
C      DATA UEPS/Z990/1991. 3/4 JST J.Takeuchi
C                ^^^^Kishu izon no Teisuu zyokyo
      DATA NAME(1),NAME(2)/'DET','BM  '/
C
C     CHECK OF THE PARAMETERS
C
      IF (NV.LT.1 .OR. NV.GT.N)     GO TO 9000
      IF (N.GT.NV1)                 GO TO 9000
      IER=0
      NM1=N-1
      NM2=N-2
      IF (N.NE.1) GO TO 10
      V(1,1)=ONE
      IFLG(1)=1
      GO TO 9999
C
C     DECIDE EPS1 TO REPLACE ZERO PIVOT IN DECOMPOSITION AND
C                 TO MODIFY CLOSE EIGENVALUES
C            EPS2 TO CHECK FOR GROUPING
C            EPS4 TO SET THE INITIAL VECTOR
C
   10 TNORM=ABS(D(N))
      DO 30 I=1,NM1
         TNORM=TNORM+ABS(D(I))+ABS(D1(I))
   30 CONTINUE
      IF (TNORM.EQ.ZERO) TNORM=ONE
      FN=FLOAT(N)
      EPS1=MAX(TNORM*UEPS,EPS)
      EPS2=TNORM*DCM3
      EPS3=FN*EPS1
      EPS4=SQRT(FN)*EPS1
      EPS5=FN*EPS4
      NGRP=0
      IPNT=1
      DO 330 JJ=1,NV
         ITER=1
         NSCNT=0
         NSROW=0
         EJJ=E(JJ)
         IF (JJ.EQ.1) GO TO 60
C
C     CHECK IF THE EIGENVALUE IS CLOSE TO THE PREVIOUS ONE OR NOT
C
         T=ABS(EJJ-E(JJ-1))
         IF (T.GT.EPS1) GO TO 40
         NGRP=NGRP+1
         EJJ=EJJ+FLOAT(NGRP)*EPS1
         GO TO 60
C
   40    IF (T.LT.EPS2) GO TO 50
         NGRP=0
         GO TO 60
C
   50    NGRP=NGRP+1
C
C     INITIALIZE THE VECTOR
C
   60    DO 70 I=1,N
            V(I,JJ)=EPS4
   70    CONTINUE
C
C     TRIANGULAR DECOMPOSITION WITH INTERCHANGES
C
         DO 80 I=1,N
            WK2(I)=D(I)-EJJ
            WK3(I)=D1(I)
   80    CONTINUE
         DO 100 I=1,NM1
            IF (ABS(D1(I)).LE.ABS(WK2(I))) GO TO 90
            IWK(I)=1
            PIV=D1(I)
            S=ONE/PIV
            WK1(I)=WK2(I)*S
            WK2(I)=S
            S=WK2(I+1)
            WK2(I+1)=WK3(I)-S*WK1(I)
            WK3(I)=S
            WK4(I)=WK3(I+1)
            WK3(I+1)=-WK4(I)*WK1(I)
            GO TO 95
C
   90       IWK(I)=0
            IF (WK2(I).EQ.ZERO) WK2(I)=EPS1
            PIV=WK2(I)
            WK2(I)=ONE/PIV
            WK1(I)=D1(I)*WK2(I)
            WK2(I+1)=WK2(I+1)-WK1(I)*WK3(I)
            WK4(I)=ZERO
   95       IF (ABS(PIV).GT.EPS3) GO TO 100
            IF (NSCNT.EQ.NGRP) NSROW=I
            NSCNT=NSCNT+1
  100    CONTINUE
         IF (WK2(N).EQ.ZERO) WK2(N)=EPS1
         IF (ABS(WK2(N)).GT.EPS3) GO TO 105
         IF (NSCNT.EQ.NGRP) NSROW=N
  105    IF (NSROW.NE.0) V(NSROW,JJ)=V(NSROW,JJ)+EPS4
         WK2(N)=ONE/WK2(N)
C
C     BACKWARD SUBSTITUTION
C
  110    V(N,JJ)=V(N,JJ)*WK2(N)
         V(NM1,JJ)=(V(NM1,JJ)-WK3(NM1)*V(N,JJ))*WK2(NM1)
         IF (N.EQ.2) GO TO 130
         DO 120 K=1,NM2
            I=NM1-K
            V(I,JJ)=(V(I,JJ)-WK3(I)*V(I+1,JJ)-WK4(I)*V(I+2,JJ))*WK2(I)
  120    CONTINUE
  130    VNORM=ZERO
         DO 140 I=1,N
            VNORM=VNORM+ABS(V(I,JJ))
  140    CONTINUE
         IF (ITER.EQ.1) GO TO 230
         IF (NGRP.EQ.0) GO TO 200
C
C     ORTHOGONALIZE WITH RESPECT TO PREVIOUS MEMBERS
C        OF THE GROUP
C
         JS=JJ-NGRP
         JE=JJ-1
         DO 170 J=JS,JE
            SUM=ZERO
            DO 150 I=1,N
               SUM=SUM+V(I,JJ)*V(I,J)
  150       CONTINUE
            IF (SUM.EQ.ZERO) GO TO 170
            S=-SUM
            DO 160 I=1,N
               V(I,JJ)=V(I,JJ)+S*V(I,J)
  160       CONTINUE
  170    CONTINUE
C
C     CHECK FOR THE CONVERGENCE
C
         VNORM=ZERO
         DO 190 I=1,N
            VNORM=VNORM+ABS(V(I,JJ))
  190    CONTINUE
  200    IF (VNORM.GE.ONE) GO TO 300
         IF (ITER.LT.5) GO TO 230
C
C     NO CONVERGENCE
C
         CALL SUERM(NAME,1000,IER)
         IFLG(JJ)=0
         GO TO 305
C
  230    IF (VNORM.NE.ZERO) GO TO 260
C
C     CHANGE THE INITIAL VECTOR
C
         V(IPNT,JJ)=EPS5
         IPNT=IPNT+1
         IF (IPNT.GT.N) IPNT=1
         GO TO 280
C
C
C     NORMALIZE THE INITIAL VECTOR
C
  260    S=EPS5/VNORM
         DO 270 I=1,N
            V(I,JJ)=V(I,JJ)*S
  270    CONTINUE
C
C     FORWARD SUBSITTUTION
C
  280    DO 290 I=1,NM1
            T=V(I+1,JJ)
            IF (IWK(I).EQ.0) GO TO 285
            T=V(I,JJ)
            V(I,JJ)=V(I+1,JJ)
  285       V(I+1,JJ)=T-WK1(I)*V(I,JJ)
  290    CONTINUE
         ITER=ITER+1
         GO TO 110
C
C
C     NORMALIZE AS NORM=1
C
  300    IFLG(JJ)=1
  305    SUM=ZERO
         DO 310 I=1,N
            SUM=SUM+V(I,JJ)*V(I,JJ)
  310    CONTINUE
         IF (SUM.EQ.ZERO) GO TO 330
         S=ONE/SQRT(SUM)
         DO 320 I=1,N
            V(I,JJ)=V(I,JJ)*S
  320    CONTINUE
  330 CONTINUE               
C
C     NORMAL END
C
      GO TO 9999
C
C
C     PARAMETER ERROR
C
 9000 CALL SUERM(NAME,2000,IER)
 9999 CALL SURTM(NAME,IER)
      RETURN
      END
                 
       


      SUBROUTINE SUERM(NAME,ISET,IER)
C
C------------------------------------------------------------------------------
C
C     NAME                - SUERM
C                            ALIAS (SUSTM,SUCRM,SURTM)
C
C     USAGE               - CALL SUERM(NAME,ISET,IER)
C                         - CALL SUSTM(LEVEL,NFILE)
C                         - CALL SUCRM(KLEVEL,KFILE)
C                         - CALL SURTM(NAME,IER)
C
C     FUNCTION            - SUERM : WRITE ERROR CODE AND
C                                     ERROR INFORMATION.
C                           SUSTM : SET REFERENCE NUMBER FOR
C                                     ERROR INFORMATION OUTPUT DATASET.
C                           SUCRM : OUTPUT THE VALUE OF COMMON BLOCK
C                                     SUCOM.
C                           SURTM : WRITE ENDING MESSARGE.
C                             SUSTM,SUCRM AND SURTM ARE ENTRY NAME OF
C                               SUERM.
C
C     REQD. SUBPROGRAM    - NONE.
C
C     STATUS              - S-1511-1 05-03
C
C     HISTORY             - DATE.    1980. 3
C                                    1986. 6
C                                    1988.11
C
C------------------------------------------------------------------------------
C
      EXTERNAL SUCO3
      DIMENSION NAME(2)
      COMMON /SUCOM/ MLEVEL,NOUT
      IER=ISET
      ILEVEL=MLEVEL
      IF (ILEVEL.GE.11) ILEVEL=ILEVEL-10
      IF (IER.LT.1000 .OR. ILEVEL.EQ.4) GO TO 999
      IF (ILEVEL.NE.2) GO TO 10
        IF (IER.LT.2000) GO TO 999
        GO TO 30
C
   10 IF (ILEVEL.NE.3) GO TO 20
        IF (IER.LT.3000) GO TO 999
        GO TO 40
C
   20 IF (IER.GE.2000) GO TO 30
C        WRITE(NOUT,100) NAME,IER
        GO TO 999
C
   30 IF (IER.GE.3000) GO TO 40
C        WRITE(NOUT,200) NAME,IER
        GO TO 999
C
C   40 WRITE(NOUT,300) NAME,IER
   40 CONTINUE
      GO TO 999
C
      ENTRY SUSTM(LEVEL,NFILE)
C
      MLEVEL=1
      IF (LEVEL.GE.1 .AND. LEVEL.LE.4) MLEVEL=LEVEL
      IF (LEVEL.GE.11 .AND. LEVEL.LE.14) MLEVEL=LEVEL
      NOUT=NFILE
      GO TO 999
C
      ENTRY SUCRM(KLEVEL,NFILE)
C
      KLEVEL=MLEVEL
      KFILE=NOUT
      GO TO 999
C
      ENTRY SURTM(NAME,IER)
C
C      IF (MLEVEL.GE.11) WRITE(NOUT,400) NAME,IER
  999 RETURN
C
  100 FORMAT(1H0,'*** WARNING ERROR   FROM MSLII ROUTINE ',2A4,
     *       ' ( IER = ',I6,' ) ***')
  200 FORMAT(1H0,'*** PARAMETER ERROR FROM MSLII ROUTINE ',2A4,
     *       ' ( IER = ',I6,' ) ***')
  300 FORMAT(1H0,'*** TERMINAL ERROR  FROM MSLII ROUTINE ',2A4,
     *       ' ( IER = ',I6,' ) ***')
  400 FORMAT(1H0,'*** RETURN MESSAGE  FROM MSLII ROUTINE ',2A4,
     *       ' ( IER = ',I6,' ) ***')
      END
      BLOCK DATA SUCO3
      COMMON /SUCOM/ MLEVEL,NOUT
      DATA MLEVEL,NOUT /1,6/
      END  
