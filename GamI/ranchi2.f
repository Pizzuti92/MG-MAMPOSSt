	double precision function ranchi2(nu,seed)
c
	implicit none
	integer nu, seed
	real*8 aerr, xmin, xmax, ran2, gammln, rtsafe
	real*8 nuover2, gamfac, q
	common /area/ nuover2, gamfac, q
	external franchi2
c
	if (nu .lt. 2) then
		print *,' RANCHI2: nu must be greater or equal to 2'
		stop
	endif
c
	aerr = 1.d-4
	xmin = 0.d0
	xmax = dacos(0.d0)
c
	q = ran2(seed)
	nuover2 = 0.5d0*dfloat(nu)
	gamfac = dexp(gammln(nuover2))
	ranchi2 = dtan(rtsafe(franchi2,xmin,xmax,aerr))
c
	return
	end
c
c
	subroutine franchi2(x,y,dy)
c
	implicit none
	real*8 x, y, dy
c
	real*8 nuover2, gamfac, q
	common /area/ nuover2, gamfac, q
c
	real*8 gammp, hafchi2, eps
c
	eps = 1.d-4
	if ((dacos(0.d0)-x) .lt. eps) then
		y = 1.d0 - q
		dy = 0.d0
	elseif (x .lt. eps) then
		y = (0.5d0*x)**nuover2/(nuover2*gamfac) - q
		if (nuover2 .gt. 1.d0) then
			dy = 0.5**nuover2*x**(nuover2-1.d0)/gamfac
		elseif (nuover2 .eq. 1.d0) then
			dy = 0.5d0
		endif
	else
		hafchi2 = 0.5d0*dtan(x)
		y = gammp(nuover2,hafchi2) - q
		dy = 0.5d0*dexp(-hafchi2)*hafchi2**(nuover2-1.d0)/
     &	 dcos(x)**2.d0/gamfac
	endif
c
	return
	end

