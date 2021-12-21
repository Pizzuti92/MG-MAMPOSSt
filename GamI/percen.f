	double precision function percen(x,n,q)
c
c PERCEN determines the 100*qth percentile for the sorted array x.
c
	integer n
	real*8 x(n), q
	integer iz
	real*8 z, dz
c
	z = 1.d0 + q*dfloat(n-1)
	iz = int(z)
	dz = z - dfloat(iz)
	percen = (1.d0-dz)*x(iz) + dz*x(iz+1)
c
	return
	end
c
c
c
