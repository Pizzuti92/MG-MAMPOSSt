	double precision function angle(theta1,phi1,theta2,phi2)
c
c ANGLE determines the angle in radians between two points whose spherical
c coordinates are specified (in radians). theta = 0 at equator.
c
	implicit none
	real*8 theta1, phi1, theta2, phi2
	real*8 arccos
c
	if (theta1 .eq. theta2 .and. phi1 .eq. phi2) then
		angle = 0.d0
	else
c		angle = arccos(dcos(theta1)*dcos(theta2)*dcos(phi1-phi2)
		angle = dacos(dcos(theta1)*dcos(theta2)*dcos(phi1-phi2)
     &		 + dsin(theta1)*dsin(theta2))
	endif
c
	return
	end
c
c
c
