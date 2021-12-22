      SUBROUTINE linmin(p,xi,n,fret)
      INTEGER*4 n,NMAX
      DOUBLE PRECISION fret,p(n),xi(n),TOL
      PARAMETER (NMAX=50,TOL=1.d-4)
CU    USES brent,f1dim,mnbrak
      INTEGER j,ncom
      DOUBLE PRECISION ax,bx,fa,fb,fx,xmin,xx,pcom(NMAX),xicom(NMAX)
     *,brent
      COMMON /f1com/ pcom,xicom,ncom
      EXTERNAL f1dim
      ncom=n
      do 11 j=1,n
        pcom(j)=p(j)
        xicom(j)=xi(j)
11    continue
      ax=0.d0
      xx=1.d0
      call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
      fret=brent(ax,xx,bx,f1dim,TOL,xmin)
      do 12 j=1,n
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software )!0.d0
