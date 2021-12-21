	double precision function ellip3(phi,ak,tol)
c
	implicit real*8 (a-h,l,m,o-z)
c
	pi = 3.1415926d0
	a = 1.d0
	b = dsqrt(1.d0-ak*ak)
	c = ak
	phidum = phi
	ipass = 0
c
111	continue
c
	print *,' ELL3: ipass a b c phi =',ipass,a,b,c,phidum
	if (dabs(c) .gt. tol) then
		aold = a
		bold = b
		cold = c
		phiold = phidum
		a = .5d0*(aold+bold)
		b = dsqrt(aold*bold)
		c = .5d0*(aold-bold)
		phidum = phiold + datan(bold/aold*dtan(phiold))
		if (phidum .lt. phiold) phidum = phidum + pi
		tanphi = dtan(phiold)
		arg = bold/aold*tanphi
		print *,' tanphi arg =',tanphi,arg
		ipass = ipass + 1
		ellip3 = phidum/(2.d0**ipass*a)
		print *,' ipass ell =',ipass,ellip3
		go to 111
	endif
c
	ellip3 = phidum/(2.d0**ipass*a)
	return
	end
