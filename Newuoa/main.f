C
C     The Chebyquad test problem (Fletcher, 1965) for N = 2,4,6 and 8,
C     with NPT = 2N+1.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(10),W(10000)
      IPRINT=2
      MAXFUN=5000
      RHOEND=1.0D-6
      DO 30 N=2,8,2
      NPT=2*N+1
      DO 10 I=1,N
   10 X(I)=DFLOAT(I)/DFLOAT(N+1)
      RHOBEG=0.2D0*X(1)
      PRINT 20, N,NPT
   20 FORMAT (//4X,'Results with N =',I2,' and NPT =',I3)
      CALL NEWUOA (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W)
   30 CONTINUE
      STOP
      END

      SUBROUTINE CALFUN (N,X,F)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),Y(10,10)
      DO 10 J=1,N
      Y(1,J)=1.0D0
   10 Y(2,J)=2.0D0*X(J)-1.0D0
      DO 20 I=2,N
      DO 20 J=1,N
   20 Y(I+1,J)=2.0D0*Y(2,J)*Y(I,J)-Y(I-1,J)
      F=0.0D0
      NP=N+1
      IW=1
      DO 40 I=1,NP
      SUM=0.0D0
      DO 30 J=1,N
   30 SUM=SUM+Y(I,J)
      SUM=SUM/DFLOAT(N)
      IF (IW .GT. 0) SUM=SUM+1.0D0/DFLOAT(I*I-2*I)
      IW=-IW
   40 F=F+SUM*SUM
      RETURN
      END



