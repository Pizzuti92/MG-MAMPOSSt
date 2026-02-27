	subroutine solvecs(x,y,n,y0,xsol,nsol)
c
c SOLVECS fits cubic splines to a set of 1D data, and returns
c the solution
c
	implicit none
	integer n, nsol
	real*8 x(n), y(n), xsol(nsol), y0
	integer ic, ier, i
	real*8 c(100,3), solve3, c1, c2, c3
c
c Obtain cubic spline coefficients
c
	ic = 100
	call icsccu(x,y,n,c,ic,ier)
c
c Search intervals for possible solutions
c
	nsol = 0
	do i = 1, n-1
c		print *,' i x yi yip1 = ', i, x(i), y(i), y(i+1)
		if ((y(i) .le. y0 .and. y(i+1) .ge. y0) .or.
     &		 (y(i) .ge. y0 .and. y(i+1) .le. y0)) then
			nsol = nsol + 1
			c1 = c(i,1)
			c2 = c(i,2)
			c3 = c(i,3)

			xsol(nsol) = x(i) + solve3(y(i),y0,c1,c2,c3,ic)
		endif
	enddo
c
	return
	end
c
c
	double precision function solve3(y1,y0,c1,c2,c3,ic)
c
	implicit none
	integer ic
	real*8 y1, y0, c1, c2, c3
	real*8 a, b, c, q, r, aa, bb, discrim
c
	a = c2/c3
	b = c1/c3
	c = (y1 - y0)/c3
c
	q = (a*a - (b+b+b))/9.d0
	r = (2.d0*a*a*a-9.d0*a*b+27.d0*c)/54.d0
c
	discrim = r*r-q*q*q
	if (discrim .lt. 0.d0) then
		print *,' SOLVE3: Q R discrim = ', q, r, discrim
		print *,' SOLVE3: a b c = ', a, b, c
		stop
	endif

	if (r .gt. 0.d0) then
		aa = -1.d0*(r + dsqrt(discrim))**(1.d0/3.d0)
	elseif (r .lt. 0.d0) then
		aa = (-1.d0*r + dsqrt(discrim))**(1.d0/3.d0)
	else
		aa = 0.d0
	endif

	if (aa .ne. 0.d0) then
		bb = q/aa
	else
		bb = 0.d0
	endif

	solve3 = aa + bb - a/3.d0

	return
	end



