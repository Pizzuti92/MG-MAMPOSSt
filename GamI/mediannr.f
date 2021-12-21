	double precision function mediannr(nmax,array)
c
c MEDIANNR is fast way to obtain medians as per Numerical recipes 2nd edition
c
	integer nmax
	real*8 array(nmax)
	integer n, j, nhalf, nhalf1
	real*8 tmparray(50000)
c
	do n = 1, nmax
		tmparray(n) = array(n)
	enddo

	if (mod(nmax,2) .eq. 1) then
		nhalf = (nmax+1)/2
		mediannr = select(nhalf,nmax,tmparray)		
	else
		nhalf = nmax/2
		nhalf1 = nhalf + 1
		mediannr = 0.5d0*(select(nhalf,nmax,tmparray) + 
     &		 select(nhalf1,nmax,tmparray))
	endif
c
	return
	end
c
c
c
