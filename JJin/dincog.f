        SUBROUTINE DINCOG(A,X,GIN,GIM,GIP)
C
C       ===================================================
C       Purpose: Compute the incomplete gamma function
C                r(a,x), ג(a,x) and P(a,x)
C       Input :  a   --- Parameter ( a ף 170 )
C                x   --- Argument 
C       Output:  GIN --- r(a,x)
C                GIM --- ג(a,x)
C                GIP --- P(a,x)
C       Routine called: GAMMA for computing ג(x)
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        XAM=-X+A*DLOG(X)
        IF (XAM.GT.700.0.OR.A.GT.170.0) THEN
           WRITE(*,*)'a and/or x too large'
           STOP
        ENDIF
        IF (X.EQ.0.0) THEN
           GIN=0.0
           CALL GAMMA(A,GA)
           GIM=GA
           GIP=0.0
        ELSE IF (X.LE.1.0+A) THEN
           S=1.0D0/A
           R=S
           DO 10 K=1,60
              R=R*X/(A+K)
              S=S+R
              IF (DABS(R/S).LT.1.0D-15) GO TO 15
10         CONTINUE
15         GIN=DEXP(XAM)*S
           CALL GAMMA(A,GA)
           GIP=GIN/GA
           GIM=GA-GIN
        ELSE IF (X.GT.1.0+A) THEN
           T0=0.0D0
           DO 20 K=60,1,-1
              T0=(K-A)/(1.0D0+K/(X+T0))
20         CONTINUE
           GIM=DEXP(XAM)/(X+T0)
           CALL GAMMA(A,GA)
           GIN=GA-GIM
           GIP=1.0D0-GIM/GA
        ENDIF
        END


c        SUBROUTINE GAMMA(X,GA)
C
C       ==================================================
C       Purpose: Compute gamma function ג(x)
C       Input :  x  --- Argument of ג(x)
C                       ( x is not equal to 0,-1,-2,תתת)
C       Output:  GA --- ג(x)
C       ==================================================
C
c        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c        DIMENSION G(26)
c        PI=3.141592653589793D0
c        IF (X.EQ.INT(X)) THEN
c           IF (X.GT.0.0D0) THEN
c              GA=1.0D0
c              M1=X-1
c              DO 10 K=2,M1
c10               GA=GA*K
c           ELSE
c              GA=1.0D+300
c           ENDIF
c        ELSE
c           IF (DABS(X).GT.1.0D0) THEN
c              Z=DABS(X)
c              M=INT(Z)
c              R=1.0D0
c              DO 15 K=1,M
c15               R=R*(Z-K)
c              Z=Z-M
c           ELSE
c              Z=X
c           ENDIF
c           DATA G/1.0D0,0.5772156649015329D0,
c     &          -0.6558780715202538D0, -0.420026350340952D-1,
c     &          0.1665386113822915D0,-.421977345555443D-1,
c     &          -.96219715278770D-2, .72189432466630D-2,
c     &          -.11651675918591D-2, -.2152416741149D-3,
c     &          .1280502823882D-3, -.201348547807D-4,
c     &          -.12504934821D-5, .11330272320D-5,
c     &          -.2056338417D-6, .61160950D-8,
c     &          .50020075D-8, -.11812746D-8,
c     &          .1043427D-9, .77823D-11,
c     &          -.36968D-11, .51D-12,
c     &          -.206D-13, -.54D-14, .14D-14, .1D-15/
c           GR=G(26)
c           DO 20 K=25,1,-1
c20            GR=GR*Z+G(K)
c           GA=1.0D0/(GR*Z)
c           IF (DABS(X).GT.1.0D0) THEN
c              GA=GA*R
c              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
c           ENDIF
c        ENDIF
c        RETURN
c        END
