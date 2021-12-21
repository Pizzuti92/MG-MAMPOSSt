	double precision function erfun2(x)
c
c ERFUN2 computes the fraction of vectors with magnitudes greater than
c x, if the projections of the vectors have a Maxwellian distribution.
c This is:
c             erfun2 = erf(x) - x erf'(x)
c
	implicit real*8(a-h,l,m,o-z)
	data third/.3333333333333d0/,fifth/.2d0/,fortnh/7.1428571429d-2/
	data ono54/.0185185185186d0/
	data tosrpi/1.12837916710d0/,fosrpi/2.25675833419d0/
	data p/0.3275911d0/
	data a1/0.254829592d0/,a2/-0.284496736d0/,a3/1.421413741d0/
	data a4/-1.453152027d0/,a5/1.061405429d0/
c
	x2 = x*x
	t = 1.d0 + p*x
c
	if (x .ge. 6.d0) then
		erfun2 = 1.d0
	elseif (x .gt. 0.4d0) then
		t2 = t*t
		t4 = t2*t2
		eminx2 = dexp(-x2)
		erfun2 = 1.d0-eminx2*(a1/t+a2/t2+a3/(t*t2)+a4/t4
     &		+a5/(t4*t)+tosrpi*x)
	elseif (x .ge. 0.0d0) then
		x4 = x2*x2
		erfun2 = fosrpi*x2*x*(third-fifth*x2+fortnh*x4-ono54*x4*x2)
	else
		print *,'ERFUN2: x = ',x,' < 0'
		stop
	endif
c
	return
	end
c
c
c
