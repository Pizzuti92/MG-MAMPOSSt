	double precision function median (x,imax)
c
c MEDIAN computes the median value of a set of imax <= 1000 real numbers
c
	implicit none
	integer imax
	real*8 x(imax)
	integer i, iarg
	real*8 xdum(1000000)
	if (imax .gt. 1000000) then
		print *,'MEDIAN cannot handle ',imax,' elements'
		stop
	endif
c
c Sort the array (in increasing order).
c
	do 100 i = 1, imax
		xdum(i) = x(i)
100	continue

	call sort(imax,xdum)	! Numerical Recipes QUICKSORT
				! Numerical recipes SELECT is not faster
c
c Determine median (depending on whether imax is even or odd).
c
	if (mod(imax,2) .eq. 0) then
		iarg = imax/2
		median = 0.5d0*(xdum(iarg)+xdum(iarg+1))
	else
		iarg = (imax+1)/2
		median = xdum(iarg)
	endif
c
	return
	end
c
c
c
