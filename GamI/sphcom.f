	subroutine sphcom(theta, phi, m, n, theta0, phi0)
c
c SPHCOM returns the center of mass in spherical coordinates
c (assuming spherical coordinates in radians).
c
	implicit none
	integer n
	real*8 theta(n), phi(n), m(n), theta0, phi0
	integer i
	real*8 sum0, sum1, sum2, mm, th, ph, sum3, denom, pi, arg2, arg3
	real*8 arccos
c
	pi = dacos(-1.d0)
c
c Find sum in cartesian coordinates of weighted vectors.
c
	sum0 = 0.d0
	sum1 = 0.d0
	sum2 = 0.d0
	sum3 = 0.d0
	do 100 i = 1, n
		th = theta(i)
		ph = phi(i)
		mm = m(i)
		sum0 = sum0 + mm
		sum1 = sum1 + mm*dsin(th)
		sum2 = sum2 + mm*dcos(th)*dsin(ph)
		sum3 = sum3 + mm*dcos(th)*dcos(ph)
100	continue
c
c Find factor denom not necessarily equal to sum0, which forces resultan
c vector to have unit magnitude.
c
	denom = dsqrt(sum1*sum1 + sum2*sum2 + sum3*sum3)
c
c If resultant vector is smack right in the center of the sphere, then c
c and exit.
c
	if (denom .eq. 0.d0) then
		print *,' SPHCOM: cannot find weighted centroid'
		do 120 i = 1, n
			print *,' i theta phi weight =',i,theta(i),phi(i),m(i)
120		continue
		return
	endif
c
c Compute theta0
c
	theta0 = dasin(sum1/denom)
c
c If theta0 is of (hopefully only a little bit from extremum permissible
c then set at correct extrema (poles) and set phi0 to 0 (because phi0 is
c undefined).
c
	if (theta0 .gt. pi/2.d0) then
		theta0 = pi/2
		phi0 = 0.d0
	elseif (theta0 .lt. -1.d0*pi/2.d0) then
		theta0 = -1.d0*pi/2
		phi0 = 0.d0
	elseif (dabs(theta0) .eq. pi/2.d0) then
		phi0 = 0.d0
	else
		arg2 = sum2/denom/dcos(theta0)
		arg3 = sum3/denom/dcos(theta0)
		if (arg2 .ge. 0.d0) then
			phi0 = arccos(arg3)
		else
			phi0 = pi + pi - arccos(arg3)
		endif
	endif
c
	return
	end
c
c
c
