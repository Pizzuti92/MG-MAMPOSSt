	character*24 function dtoa(x,n)
c
c DTOA converts a real*8 to a character string with n numbers after the
c (format = g15.n)
c
	implicit none
	integer n
	real*8 x
	character*24 dummy
c
	if (n .eq. 0) then
		write(dummy,90) x
	elseif (n .eq. 1) then
		write(dummy,91) x
	elseif (n .eq. 2) then
		write(dummy,92) x
	elseif (n .eq. 3) then
		write(dummy,93) x
	elseif (n .eq. 4) then
		write(dummy,94) x
	elseif (n .eq. 5) then
		write(dummy,95) x
	elseif (n .eq. 6) then
		write(dummy,96) x
	else
		print *,'  n = ', n, ' ... is out of range'
		stop
	endif
c
90	format('''',1pg22.0,'''')
91	format('''',1pg22.1,'''')
92	format('''',1pg22.2,'''')
93	format('''',1pg22.3,'''')
94	format('''',1pg22.4,'''')
95	format('''',1pg22.5,'''')
96	format('''',1pg22.6,'''')
c
	read(dummy,*) dtoa
c
	return
	end
