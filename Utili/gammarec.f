c
c     Function gammarec: uses recursive relation
c                Gamma(x)=Gamma(x+1)/x
c     when x<0 to call dgamma(arg) with arg>0
c

      function gammarec(x)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      nit=-nint(x-0.5)

      if (nit.le.0) then

         gammarec=dgamma(x)
         return

      else

         prod=1.
         do i=1,nit
            prod=prod*(x+i-1.)
         enddo
         gammarec=dgamma(x+nit)/prod
         return

      endif

      end

