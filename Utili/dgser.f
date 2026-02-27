c
c  Funzione gamma incompleta.
c
c  Nota: per valori di A < 1, si e' usata la formula di riflessione.
c
      SUBROUTINE DGSER(DGAMSER,A,X,DGLN)
      DOUBLE PRECISION A,X,DGLN,PI,Z,ZZ,dgammln,ap,dgamser,sum,del
      PARAMETER (ITMAX=100,EPS=3.D-7)
      external dgammln
      pi=4.d0*datan(1.d0)
      if(a.eq.1.d0) dgln=dlog(1.d0)
      if (a.gt.0.d0.and.a.lt.1.d0) then
        z=1.d0-a
        zz=1+z
        DGLN=dlog(pi*z)-DGAMMLN(zz)-dlog(dsin(pi*z))
      else if(a.gt.1.) then
        dgln=dgammln(a)
      endif
      IF(X.LE.0.d0)THEN
        IF(X.LT.0.d0) stop
        DGAMSER=0.d0
        RETURN
      ENDIF
      AP=A
      SUM=1.d0/A
      DEL=SUM
      DO 11 N=1,ITMAX
        AP=AP+1.d0
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF(DABS(DEL).LT.DABS(SUM)*EPS)GO TO 1
11    CONTINUE
      STOP('A too large, ITMAX too small')
1     DGAMSER=SUM*DEXP(-X+A*DLOG(X)-DGLN)
      RETURN
      END
