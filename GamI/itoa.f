	character*12 function itoa(n)
c
c ITOA converts an integer to a character string
c
	implicit none
	integer n
	character*12 dummy
c
	write(dummy,90) n
90	format('''',i10,'''')
	read(dummy,*) itoa
c
	return
	end
