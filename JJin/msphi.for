        PROGRAM MSPHI
C
C       ======================================================
C       Purpose: This program computes the modified spherical 
C                Bessel functions of the first kind in(x) and 
C                in'(x) using subroutine SPHI
C       Input :  x --- Argument of in(x)
C                n --- Order of in(x) ( 0 � n � 250 )
C       Output:  SI(n) --- in(x)
C                DI(n) --- in'(x)
C       Example: x = 10.0
C                  n          in(x)               in'(x)
C                --------------------------------------------
C                  0     .1101323287D+04     .9911909633D+03
C                  1     .9911909633D+03     .9030850948D+03
C                  2     .8039659985D+03     .7500011637D+03
C                  3     .5892079640D+03     .5682828129D+03
C                  4     .3915204237D+03     .3934477522D+03
C                  5     .2368395827D+03     .2494166741D+03
C       ======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION SI(0:250),DI(0:250)
        WRITE(*,*)'Please enter n and x '
        READ(*,*)N,X
        WRITE(*,30)N,X
        IF (N.LE.10) THEN
           NS=1
        ELSE
           WRITE(*,*) 'Please enter order step Ns'
           READ(*,*) NS
        ENDIF
        CALL SPHI(N,X,NM,SI,DI)
        WRITE(*,*)
        WRITE(*,*)'  n          in(x)               in''(x)'
        WRITE(*,*)'--------------------------------------------'
        DO 10 K=0,NM,NS
10         WRITE(*,20)K,SI(K),DI(K)
20      FORMAT(1X,I3,2D20.10)
30      FORMAT(3X,'Nmax =',I3,',     ','x =',F6.1)
        END


        SUBROUTINE SPHI(N,X,NM,SI,DI)
C
C       ========================================================
C       Purpose: Compute modified spherical Bessel functions
C                of the first kind, in(x) and in'(x)
C       Input :  x --- Argument of in(x)
C                n --- Order of in(x) ( n = 0,1,2,... )
C       Output:  SI(n) --- in(x)
C                DI(n) --- in'(x)
C                NM --- Highest order computed
C       Routines called:
C                MSTA1 and MSTA2 for computing the starting
C                point for backward recurrence
C       ========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION SI(0:N),DI(0:N)
        NM=N
        IF (DABS(X).LT.1.0D-100) THEN
           DO 10 K=0,N
              SI(K)=0.0D0
10            DI(K)=0.0D0
           SI(0)=1.0D0
           DI(1)=0.333333333333333D0
           RETURN
        ENDIF
        SI(0)=DSINH(X)/X
        SI(1)=-(DSINH(X)/X-DCOSH(X))/X
        SI0=SI(0)
        IF (N.GE.2) THEN
           M=MSTA1(X,200)
           IF (M.LT.N) THEN
              NM=M
           ELSE
              M=MSTA2(X,N,15)
           ENDIF
           F0=0.0D0
           F1=1.0D0-100
           DO 15 K=M,0,-1
              F=(2.0D0*K+3.0D0)*F1/X+F0
              IF (K.LE.NM) SI(K)=F
              F0=F1
15            F1=F
           CS=SI0/F
           DO 20 K=0,NM
20            SI(K)=CS*SI(K)
        ENDIF
        DI(0)=SI(1)
        DO 25 K=1,NM
25         DI(K)=SI(K-1)-(K+1.0D0)/X*SI(K)
        RETURN
        END


        INTEGER FUNCTION MSTA1(X,MP)
C
C       ===================================================
C       Purpose: Determine the starting point for backward  
C                recurrence such that the magnitude of    
C                Jn(x) at that point is about 10^(-MP)
C       Input :  x     --- Argument of Jn(x)
C                MP    --- Value of magnitude
C       Output:  MSTA1 --- Starting point   
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        N0=INT(1.1*A0)+1
        F0=ENVJ(N0,A0)-MP
        N1=N0+5
        F1=ENVJ(N1,A0)-MP
        DO 10 IT=1,20             
           NN=N1-(N1-N0)/(1.0D0-F0/F1)                  
           F=ENVJ(NN,A0)-MP
           IF(ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
 10        F1=F
 20     MSTA1=NN
        RETURN
        END


        INTEGER FUNCTION MSTA2(X,N,MP)
C
C       ===================================================
C       Purpose: Determine the starting point for backward
C                recurrence such that all Jn(x) has MP
C                significant digits
C       Input :  x  --- Argument of Jn(x)
C                n  --- Order of Jn(x)
C                MP --- Significant digit
C       Output:  MSTA2 --- Starting point
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        HMP=0.5D0*MP
        EJN=ENVJ(N,A0)
        IF (EJN.LE.HMP) THEN
           OBJ=MP
           N0=INT(1.1*A0)
        ELSE
           OBJ=HMP+EJN
           N0=N
        ENDIF
        F0=ENVJ(N0,A0)-OBJ
        N1=N0+5
        F1=ENVJ(N1,A0)-OBJ
        DO 10 IT=1,20
           NN=N1-(N1-N0)/(1.0D0-F0/F1)
           F=ENVJ(NN,A0)-OBJ
           IF (ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
10         F1=F
20      MSTA2=NN+10
        RETURN
        END

        REAL*8 FUNCTION ENVJ(N,X)
        DOUBLE PRECISION X
        ENVJ=0.5D0*DLOG10(6.28D0*N)-N*DLOG10(1.36D0*X/N)
        RETURN
        END
