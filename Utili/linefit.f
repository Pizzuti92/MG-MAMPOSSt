      subroutine linefit(x,y,n,b0,b1)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension x(n),y(n)

C Least square fitting of a straight line to data points

      sumx=0.d0
      sumy=0.d0
      sumsqx=0.d0
      sumxy=0.d0
      do i=1,n
         sumx=sumx+x(i)
         sumy=sumy+y(i)
         sumsqx=sumsqx+x(i)*x(i)
         sumxy=sumxy+x(i)*y(i)
      enddo
      deno=n*sumsqx-sumx*sumx
      b1=(n*sumxy-sumx*sumy)/deno  ! slope
      b0=(sumsqx*sumy-sumx*sumxy)/deno ! intercept
c      write(*,*)'Slope, Intercept= ',b1,b0
      return
      end
