	integer function atoi(string)
c
c ATOI converts a string into an integer*4 number.
c
	implicit none
	character*80 string
	character*80 dummy
c
	write(dummy,90) string
90	format(a80)
	read(dummy,*) atoi
c
	return
	end
c
c
c
