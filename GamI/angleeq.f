	double precision function angleeq(alpha1,delta1,alpha2,delta2)
c
c ANGLEEQ determines the angle in degrees between two points whose equatorial
c coordinates are specified (in degrees).
c
	implicit none
	real*8 delta1, alpha1, delta2, alpha2
	real*8 degree, eps, arg
c	real*8 arccos
c
	degree = dacos(-1.d0)/180.d0
	eps = 1.d-8
c
	if (delta1 .eq. delta2 .and. alpha1 .eq. alpha2) then
		angleeq = 0.d0
	else
		arg = dcos(delta1*degree)*dcos(delta2*degree)
     &	 *dcos((alpha1-alpha2)*degree)
     &	 + dsin(delta1*degree)*dsin(delta2*degree)
		if (dabs(arg) .le. 1.d0) then
			angleeq = dacos(arg)/degree
		elseif (arg .ge. -1.d0-eps) then
			angleeq = 180.d0
		elseif (arg .le. 1.d0-eps) then
			angleeq = 0.d0
		else
			print *, ' ANGLEEQ: cos(angle) = ', arg
			stop
		endif
	endif
c
	return
	end
c
c
c
