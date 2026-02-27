	subroutine eqtogal(alpha,delta,lrad,brad)
c
c EQTOGAL converts RA and DEC in radians to lII and bII in radians.
c
	implicit none
	real*8 alpha, delta, lrad, brad
	real*8 pi, degree, eps, eps1, epsm1, ang1, ang2, ang3
	real*8 sinang1, sinang3, cosang1, cosang3, sindel, cosdel
	real*8 sinb, cosb, cosl, sinl
c
	pi = dacos(-1.d0)
	degree = pi/180.d0
	eps = 1.d-3
	eps1 = 1.d0+eps
	epsm1 = -1.d0-eps
	ang1 = 62.6d0*degree	! angle NGP,NCP
	ang2 = 282.25d0*degree	! 
	ang3 = 33.d0*degree
	sinang1 = dsin(ang1)
	cosang1 = dcos(ang1)
	sinang3 = dsin(ang3)
	cosang3 = dcos(ang3)
	sindel = dsin(delta)
	cosdel = dcos(delta)
	sinb = sindel*dcos(ang1) - cosdel*dsin(alpha-ang2)*sinang1
	if (sinb .ge. 1.d0 .and. sinb .le. eps1) then
		brad = pi/2.d0
	elseif (sinb .le. -1.d0 .and. sinb .ge. epsm1) then
		brad = -pi/2.d0
	elseif (sinb .gt. -1.d0 .and. sinb .lt. 1.d0) then
		brad = dasin(sinb)
	else
		print *,' EQTOGAL: sinb = ', sinb
		stop
	end if
c
	cosb = dcos(brad)
	if (cosb .lt. 0.d0) then
		print *,' EQTOGAL: cosb = ', cosb
		stop
	elseif (cosb .eq. 0.d0) then
		lrad = 0.d0
	else
		cosl = (cosdel*dcos(alpha-ang2)*cosang3
     &		 - cosdel*dsin(alpha-ang2)*cosang1*sinang3
     &		 - sindel*sinang1*sinang3)
     &		 /cosb
		sinl = (cosdel*dsin(alpha-ang2)*cosang1*cosang3
     &		 + sindel*sinang1*cosang3
     &		 + cosdel*dcos(alpha-ang2)*sinang3)
     &		 /cosb
		if (sinl .ge. 0.d0) then
			lrad = dacos(cosl)
		else
			lrad = pi + pi - dacos(cosl)
		end if
	end if
c
	return
	end

