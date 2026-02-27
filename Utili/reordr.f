c
c
      SUBROUTINE REORDR(N,ARRAY,INDX)
      implicit real*8 (a-h,l-m,o-z)
      DIMENSION ARRAY(N), INDX(N)
      PARAMETER (NMAX=20000)
      REAL*8    HLP(NMAX)
      DO I=1,N
         HLP(I)=ARRAY(INDX(I))
      END DO         
      DO I=1,N
         ARRAY(I)=HLP(I)
      END DO         
      RETURN
      END 
