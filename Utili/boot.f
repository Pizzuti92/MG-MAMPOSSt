c
c     bootstrap y -> yboo
c
      subroutine boot(y,w,yboo,wboo,ny)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension y(ny),w(ny),yboo(ny),wboo(ny),r(25000)
      idum = 12234
      do j=1,ny
         r(j) = ran2(idum)
         jj=nint((ny-1)*r(j)+1)
         yboo(j)=y(jj)
         wboo(j)=w(jj)
      enddo
      return
      end

