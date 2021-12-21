	double precision function arccos(x)
c
c ARCCOS computes the arc-cosine and allows the argument to have absolut
c value slightly larger than unity.
c
	implicit none
	real*8 x
	real*8 EPS, dx
	parameter (EPS=1.d-8)
c
	dx = dabs(x)
	if (dx .le. 1.d0) then
		arccos = dacos(x)
	elseif (dx .le. (1.d0+eps)) then
		arccos = dacos(x/dx)
	else
		print *,' ARCCOS: argument = ',x,' out of range'
		stop
	endif
c
	return
	end
c
c
c
