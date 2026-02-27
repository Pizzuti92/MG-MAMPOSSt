	double precision function ellip2(iell,phi,ak,n)
c
c ELLIP2 computes incomplete elliptical integrals of the forst kind,
c F(phi,k), where phi is close to pi/2, by
c F = K(k) - simson(f,phi,k,n) where f = 1/sqrt(1-k2sin2x)
c n is the number of knots for the Simpson routine (1 is sufficient,
c except when k is close to 1: 11 is sufficient for k = .999)
c
	implicit real*8 (a-h,o-z)
	common /area/ akk, iellip
	double precision mmdelk, mmdele
	external felip
	pi = 3.1415926d0
	hafpi = .5d0*pi
	akk = ak
	iellip = iell
	if (iell .eq. 0) compl = mmdelk(2,ak,ier)
	if (iell .eq. 1) compl = mmdele(2,ak,ier)
	if (ier .eq. 130) then
		print *,' FUNC: elliptical integration failed because ...'
		print *,' k = ',ak
		stop
	endif
	ellip2 = compl - simson(felip,phi,hafpi,n)
	return
	end
c
	double precision function felip(x)
c
	implicit real*8 (a-h,o-z)
	common /area/ ak, iell
	g = dsqrt(1.d0-ak*ak*dsin(x)*dsin(x))
	if (iell .eq. 0) felip = 1.d0/g
	if (iell .eq. 1) felip = g
	return
	end
c
c
c
