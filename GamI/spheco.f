	subroutine spheco(rah,ram,ras,ded,dem,des,sgn,phi,theta)
c
c SPHECO determines the spherical coordinates theta and phi (in degrees)
c from RA and DEC in hours mins secs, degrees, mins, secs.
c
	implicit none
	integer rah, ded, ram, dem, des
	real*8 ras
	character*1 sgn
	real*8 theta, phi
c
	theta = dfloat(ded) + dfloat(dem)/60.d0 + dfloat(des)/3600.d0
	if (sgn .eq. '-') theta = -1.d0*theta
c
	phi = 15.d0*(dfloat(rah) + dfloat(ram)/60.d0 + ras/3600.d0)
c
	return
	end
c
c
c
