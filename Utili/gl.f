C
c
c  calcola l'integrale usando la tecnica di Gauss-Legendre
c
c  il numero di punti usati viene aumentato fino a che il risultato non
c  raggiunge la stabilita` richiesta.
c
      subroutine gl(f,ax,bx,eps,ris,ifl)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (nmax=500)
      real*8 x(nmax), wx(nmax)
      ifl=0
c
c  inizializza il numero di punti ed ilvalore del risultato
c
      ris0=-1.d+10
      n=10
c
c  ciclo di integrazione
c
 100  continue
      if(n.ge.nmax) then
        print*,' n > nmax !'
        ifl=1
        ris=0.d0
        return
      endif
      call gauleg(ax,bx,x,wx,n)
      sum=0.d0
      do 1 i=1,n
        xi=x(i)
        sum=sum+wx(i)*f(xi)
   1  continue
      err=dabs(sum-ris0)
      if (err.ge.dabs(eps*sum)) then
        n=n+2
        ris0=sum
        goto 100
      else
c        print*,' nr. pt. = ',n
        ris=sum
      endif
      return
      end 
c
c
c
      SUBROUTINE GAULEG(X1,X2,X,W,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      REAL*8 X1,X2,X(N),W(N)
      PARAMETER (EPS=3.D-14)
      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)
      DO 12 I=1,M
        Z=DCOS(3.141592654D0*(I-.25D0)/(N+.5D0))
1       CONTINUE
          P1=1.D0
          P2=0.D0
          DO 11 J=1,N
            P3=P2
            P2=P1
            P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
11        CONTINUE
          PP=N*(Z*P1-P2)/(Z*Z-1.D0)
          Z1=Z
          Z=Z1-P1/PP
        IF(DABS(Z-Z1).GT.EPS)GO TO 1
        X(I)=XM-XL*Z
        X(N+1-I)=XM+XL*Z
        W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
        W(N+1-I)=W(I)
12    CONTINUE
      RETURN
      END
C

