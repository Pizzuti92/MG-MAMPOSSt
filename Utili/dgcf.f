      SUBROUTINE DGCF(DGAMMCF,A,X,DGLN)
      double precision dgammcf,a,x,dgln,dgammln
      double precision a0,a1,b0,b1,fac,an,ana,anf,g,gold
      PARAMETER (ITMAX=100,EPS=3.d-7)
      GLN=DGAMMLN(A)
      GOLD=0.d0
      A0=1.d0
      A1=X
      B0=0.d0
      B1=1.d0
      FAC=1.d0
      DO 11 N=1,ITMAX
        AN=dFLOAT(N)
        ANA=AN-A
        A0=(A1+A0*ANA)*FAC
        B0=(B1+B0*ANA)*FAC
        ANF=AN*FAC
        A1=X*A0+ANF*A1
        B1=X*B0+ANF*B1
        IF(A1.NE.0.d0)THEN
          FAC=1.d0/A1
          G=B1*FAC
          IF(ABS((G-GOLD)/G).LT.EPS)GO TO 1
          GOLD=G
        ENDIF
11    CONTINUE
      stop('A too large, ITMAX too small')
1     dGAMMCF=dEXP(-X+A*dLOG(X)-dGLN)*G
      RETURN
      END
