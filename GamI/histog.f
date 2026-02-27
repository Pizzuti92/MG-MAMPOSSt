	subroutine histog(x,nmax,bin0,binf,jbin,bin,isc,an)
c
c HISTOG evaluates histograms and returns either the differential histog
c (for isc = 0), or the normalized integrated histogram (isc = 1).
c
	implicit real*8 (a-h,l,m,o-z)
	real*8 x(nmax), an(jbin), bin(jbin)
c
	delb = (binf-bin0)/float(jbin)
	d2 = 0.5d0*delb
	bin2 = bin0 
	anint = 0.d0
	do 300 j = 1, jbin
		bin1 = bin2
		bin2 = bin1 + delb
		bin(j) = bin1 + d2
		ann = 0.d0
		do 301 n = 1, nmax
			if (x(n) .ge. bin1 .and. x(n) .lt. bin2) then
				ann = ann + 1.
			endif
301		continue
c
		anint = anint + ann
		anint0 = anint/float(nmax)
c
c Total values in bins for isc = 0, relative integral for isc = 1.
c
		if (isc .eq. 0) an(j) = ann
		if (isc .eq. 1) an(j) = anint0
300	continue
c
	return
	end
c
c
c
