	real function bulelp(x,kc,a,b,itol)
c
c BULELP is Bulirsch's incomplete elliptical integral solver.
c
	implicit real*8 (a-h,o-z)
	real*8 kc, m
	common /BULL/ l, ll
c
	pi = 3.1415926d0
	ca = 10.d0**(-itol)
	cb = 10.d0**(-2*(itol+1))
c
	if (x .eq. 0.d0) then
		bulelp = 0.d0
	elseif (kc .eq. 0.d0) then
		print *,' cannot handle kc = 0, stopping ...'
		stop
	else
		c = x*x
		d = 1.d0+c
		p = dsqrt((1.d0+kc*kc*c)/d)
		d = x/d
		c = d/(2.d0*p)
		z = a-b
		ai = a
		aa = (b+a)/2.d0
		bb = b
		y = dabs(1.d0/x)
		f = 0.d0
		l = 0
		m = 1.d0
		akc = dabs(kc)
		ll = 0
c
111		continue
c
		ll = ll+1
		bb = ai*akc+bb
		e = m*akc
		g = e/p
		d = f*g+d
		f = c
		ai = aa
		p = g+p
		c = (d/p+c)/2.d0
		g = m
		m = akc+m
		aa = (bb/m+aa)/2.d0
		y = -e/y+y
		if (y .eq. 0.d0) y = dsqrt(e)*cb
		if (dabs(g-akc) .gt. ca*g) then
			akc = dsqrt(e)*2.d0
			l = l*2
			if (y .lt. 0.d0) l = l+1
			go to 111
		endif
		if (y .lt. 0.d0) l = l+1
		e = (datan(m/y)+pi*l)*aa/m
		if (x .lt. 0.d0) e = -e
		bulelp = e+c*z
	endif
	return
c
	end
