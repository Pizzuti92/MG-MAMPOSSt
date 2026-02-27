c
c     Inverse hyperbolic cosine function
c
      
      double precision function dacosh(z)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (pig=3.1415926535897932d0)
      dacosh=dlog(z+dsqrt(z*z-1.d0))
      return
      end
