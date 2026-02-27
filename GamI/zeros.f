	subroutine zeros(fzero,xmin,xmax,delx,itol,ipass,zero,izero,
     &    islope,ier)
c
c ZEROS uses ZBRENT to look for all zeros of a 1D function f(x) between
c xmin and xmax, scanning through with a resolution of delx: with a rela
c tolerance of 10**(-itol) and using a maximum of ipass passes per subin
c If delx < 0, the function is searched in logarithmic increments instea
c of linear ones.
c So two zeros that are closer than delx in linear space (or -delx in lo
c space) could be missed.
c The routine returns the izero zeros in the array 'zero' with for each
c flag indicating whether the slope is positive (islope = 1) or negative
c (islope = -1).
c
c Arguments:
c ---------
c
c f: double precision function declared external in the calling program
c
c xmin (I): lower limit of interval to test
c
c xmax (I): upper limit of interval to test
c
c delx (I): size of internal subroutine-interval to test
c delx < 0 -> -1*size in logarithmic decades
c
c itol (I): relative tolerance for finding each zero
c
c ipass (I): maximum number of passes for ZBRENT to find solution
c
c zero (O): array of zeros returned
c
c izero (O): number of zeros returned
c
c islope (O): array of integers indicating positive or negative slope at
c
c ier (O): maximum ZBRENT error
c
	implicit none
	integer izero
	integer itol, ipass, ier
	integer islope(izero)
	real*8 xmin, xmax, delx, zero(izero)
	integer i, imax, izero0, isgn1, isgn2
	real*8 x1, x2, lx1, lx2, y1, y2, eps, fzero
	external fzero
c
	izero0 = izero
	izero = 0
	ier = 0
c
	if (delx .gt. 0.d0) then
		imax = int((xmax-xmin)/delx) + 1
	elseif (delx .lt. 0.d0) then
		imax = int(dlog10(xmax/xmin)/(-1.d0*delx)) + 1
	else
		print *,' ZEROS: cannot handle delx = 0'
		stop
	endif
c
	do 100 i = 1, imax
		if (delx .gt. 0.d0) then
			x1 = xmin + delx*dfloat(i-1)
			x2 = x1 + delx
		else
			lx1 = dlog10(xmin) - delx*dfloat(i-1)
			lx2 = lx1 - delx
			x1 = 10.d0**lx1
			x2 = 10.d0**lx2
		endif
c
		if (x1 .gt. xmax) goto 100
		y1 = fzero(x1)
		y2 = fzero(x2)
		if (y1 .ne. 0.d0) then
			isgn1 = nint(y1/dabs(y1))
		else
			isgn1 = 0
		endif
		if (y2 .ne. 0.d0) then
			isgn2 = nint(y2/dabs(y2))
		else
			isgn2 = 0
		endif
		if ((isgn1*isgn2) .gt. 0) goto 100
c
		eps = 0.d0
		call zbrent(fzero,eps,itol,x1,x2,ipass,ier)
		if (ier .gt. 0) then
			print *,' ZEROS: x1 x2 ier = ', x1,x2,ier
c			goto 100
		endif
		izero = izero + 1
		if (izero .le. izero0) then
			zero(izero) = x2
			if (y1 .lt. 0.d0) then
				islope(izero) = 1
			else
				islope(izero) = -1
			endif
		else
			print *,' ZEROS: more than ', izero0, ' zeros'
			print *,' exiting now'
			go to 111
		endif
100	continue
c
111	continue
c
	return
	end
