	double precision function ellip4(phi,ak,tol)
c
c ELLIP4 computes incomplete elliptical integral of the first kind F(phi
c by the method of the ascending Landen transformation (e.g. Abramowitz
c Stegun 1970, Section 17.5).
c
	implicit real*8 (a-h,l,m,o-z)
c
	qartpi = 0.7853981633d0
	alpha = dasin(ak)
	phidum = phi
	prdsin = 1.d0
c
111	continue
c
	osinal = dsin(alpha)
	phiold = phidum
	alpha = dacos(2.d0/(1.+osinal)-1.d0)
	phidum = .5d0*(phiold+dasin(osinal*dsin(phiold)))
	prdsin = prdsin*dsin(alpha)
	ellip4 = dsqrt(prdsin/ak)*dlog(dtan(qartpi+.5d0*phidum))
c
	arg = (phidum-phiold)/phiold
	if (dabs(arg) .gt. tol) go to 111
c
	return
	end
c
c
c
