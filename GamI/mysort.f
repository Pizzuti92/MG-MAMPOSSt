	subroutine mysort(n,x)
c
c MYSORT sorts in increasing order a real*8 array x of length n
c and, according to the value of n, uses  the optimal NUMERICAL RECIPES
c
	integer n
	real*8 x(n)
	integer NMAX
	parameter (NMAX=1000000)
c
	if (n .le. 1) then
c		do nothing!
		return
	endif
c
	if (n .gt. 1 .and. n .lt. 50) then
		call piksrt(n,x)
	elseif (n .lt. 1000)  then
		call shell(n,x)
	elseif (n .lt. NMAX) then
c		print *,' about to call SORT'
		call sort(n,x)
c		print *,' past SORT'
	else
		print *,' MYSORT: n = ', n, ' ... is too large or <= 0'
		stop
	endif
c
	return
	end
c
c
c
