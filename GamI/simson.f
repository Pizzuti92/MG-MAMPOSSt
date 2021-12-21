	double precision function simson(f,a,b,n)
c
c SIMSON is a Simpson integral evaluator.
c a and b are the integration limits, and n is the number of knots.
c
	implicit real*8 (a-h,o-z)
	external f
c
	h = .5d0*(b-a)/float(n)
	x = a
	sum = f(a)
	do 100 i = 1, n-1
		x = x + h
		sum = sum + 4.d0*f(x)
		x = x + h
		sum = sum + 2.d0*f(x)
100	continue
	x = x + h
	sum = sum + 4.d0*f(x)
	sum = sum + f(b)
c
	simson = sum*h/3.d0
	return
	end
c
c
c
