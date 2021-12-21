	double precision function dotprd(x,y,n)
c
c DOTPRD evaluates the dot-product of the two n-vectors x and y.
c I.e. the vectors can have arbitrary though equal dimensions.
c
	real*8 x(n), y(n)
	dotprd = 0.d0
	do 100 i = 1, n
		dotprd = dotprd + x(i)*y(i)
100	continue
c
	return
	end
c
c
c
