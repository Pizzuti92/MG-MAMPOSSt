      FUNCTION f1dim(x)
      INTEGER*4 NMAX
      double precision f1dim,func,x
      PARAMETER (NMAX=50)
CU    USES func
      INTEGER*4 j,ncom
      double precision pcom(NMAX),xicom(NMAX),xt(NMAX)
      COMMON /f1com/ pcom,xicom,ncom
      do 11 j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
11    continue
      f1dim=func(xt)
      return
      END
