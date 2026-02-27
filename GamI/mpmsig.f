	subroutine mpmsig(x,n,xmin,xminus,xmed,xplus,xmax)
	implicit none
	integer n
	real*8 x(n), xmin, xminus, xmed, xplus, xmax
	integer i, NMAX
	parameter(NMAX=10000)
	real*8 xdum(NMAX), SIG1, SIG2, xfac
	parameter(SIG2 = 0.84134d0)	
	parameter(SIG1 = 1.d0-SIG2)
c
	if (n .eq. 0) then
		xmin = -1.d0
		xminus = -1.d0
		xmed = -1.d0
		xplus = -1.d0
		xmax = -1.d0
		return
	endif
c
	do 100 i = 1, n
		xdum(i) = x(i)
100	continue
	call qsort(xdum,n,.true.)
c
	xmin = xdum(1)
	xmax = xdum(n)
c
	xmed = xfac(xdum,n,.5d0)
	xminus = xfac(xdum,n,SIG1)
	xplus = xfac(xdum,n,SIG2)
c
	return
	end
c
c
	double precision function xfac(xdum,n,fac)
	implicit none
	integer n
	real*8 xdum(n), fac
	integer i
	real*8 j, jrem
c
	j = 1.d0 + fac*dfloat(n-1)
	i = int(j)
	jrem = j - i
c
	if (dabs(jrem-0.5d0) .lt. 0.5d0) then
		xfac = .5d0*(xdum(i)+xdum(i+1))
	elseif (jrem .gt. 0.5d0) then
		xfac = xdum(i+1)
	else
		xfac = xdum(i)
	endif
c
	return
	end
