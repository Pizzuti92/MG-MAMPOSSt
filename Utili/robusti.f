c
c
      subroutine robusti(xin,n,ibwt,c,s,ifail)

C     Routine that uses robust techniques to estimate the central location C
C     and the spread S, for a distribution of N sorted values X.
C     Based on the work of Beers, Flynn and Gebhardt, AJ 100, 32, 1990.
C     When necessary, use is made of Numerical Recipes routines.

      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      parameter (pi=3.1415926535897932,nmax=20000)

      dimension xin(n),x(nmax),xred(nmax)

      data      t1,t2 /6.0,9.0/

      ifail=0

C     Reset output variables

      c=0.0
      s=0.0

C     Copy input data into work array

      do i=1,n
         x(i)=xin(i)
         c=c+xin(i)
      end do
      c=c/float(n)

C     If N<3 calculation of spread is impossible

      if (n.lt.3) return

C     Calculate the median M, given the sorted distribution of X

      call sort(n,x)
      n2=n/2
      if(2*n2.eq.n)then
        xm=0.5*(x(n2)+x(n2+1))
      else
        xm=x(n2+1)
      endif

C     Calculate the Median Absolute Deviation

      do i=1,n
         xred(i)=abs(x(i)-xm)
      end do                 
      call sort(n,xred)
      n2=n/2
      if(2*n2.eq.n)then
         xmad=0.5*(xred(n2)+xred(n2+1))
      else
         xmad=xred(n2+1)
      endif
      if (xmad.eq.0.0) then
cc         write(*,*) 'robust: no differentiation in data!'
         ifail=1
         return
      end if

C     Calculate the biweight location estimator

      fact=1.0/(t1*xmad)
      sum1=0.0
      sum2=0.0
      do i=1,n
         u=(x(i)-xm)*fact
         if (abs(u).lt.1.0) then
            uhlp=(1.0-u**2)**2
            sum2=sum2+uhlp
            sum1=sum1+(x(i)-xm)*uhlp
         end if
      end do
      cbi=xm+sum1/sum2
      c=cbi

C     Calculate the biweight scale if requested

      if (ibwt.gt.0.5) then
         fact=1.0/(t2*xmad)
         sum1=0.0
         sum2=0.0
         do i=1,n
            u=(x(i)-xm)*fact
            if (abs(u).lt.1.0) then
               uhlp1=1.0-u**2
               uhlp2=uhlp1**4
               uhlp3=1.0-5.0*u**2
               sum1=sum1+(x(i)-xm)**2*uhlp2
               sum2=sum2+uhlp1*uhlp3
            end if
         end do
         s=sqrt(float(n))*sqrt(sum1)/abs(sum2)
      else

C     ... and use a gapper algorithm otherwise

         sum=0.0
         do i=1,n-1
            g=x(i+1)-x(i)
            w=float(i*(n-i))
            sum=sum+w*g
         end do
         s=sqrt(pi)/(n*(n-1))*sum
      end if

      return
      end
