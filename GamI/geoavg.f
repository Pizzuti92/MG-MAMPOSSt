	double precision function geoavg(x,imax)
c
c GEOAVG computes the geometric average of a positive definite array.
c This subprogram stops if it encounters any non-positive element in the
c array.
c
	implicit real*8 (a-h,l,m,o-z)
	real*8 x(imax)
	if (imax .le. 0) then
		print *,'Non-positive array length: program will stop ...'
		print *,'array length =',imax
		stop
	endif
	sum = 0.d0
	do 100 i = 1, imax
		if (x(i) .le. 0.) then
			print *,'Non-positive element in the array'	
			print *,'Cannot, proceed with geometric average'
			print *,'Program will stop ...'
			print *,'array element #',i,' = ', x(i)
			i1 = i-1
			print *,'previous elts.: ',(x(ii),ii=1,i1)
			stop
		endif
		sum = sum + dlog(x(i))
100	continue
	geoavg = dexp(sum/dfloat(imax))
c
	return
	end
c
c
c
