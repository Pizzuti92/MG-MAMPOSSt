	subroutine eqtosph(rah,ram,ras,sgn,ded,dem,des,alpha,delta)
c
	implicit none
	integer rah, ram, ded, dem
	real*8 ras, des, alpha, delta
	character*1 sgn
	real*8 degree, pi
c
	pi = dacos(-1.d0)
	degree = pi/180.d0
c
	alpha = 15.d0*(dfloat(rah) + dfloat(ram)/60.d0 + ras/3600.d0)
     &  *degree
	delta = (dfloat(ded) + dfloat(dem)/60.d0 + des/3600.d0)*degree
	if (sgn .eq. '-') then
		delta = -1.d0*delta
	endif
c
	return
	end
