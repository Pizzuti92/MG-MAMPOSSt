	double precision function atod(string)
c
c ATOD converts a string into a real*8 number.
c
	implicit none
	character*80 string
	character*80 dummy
c
	write(dummy,90) string
90	format(a80)
	read(dummy,*) atod
c
	return
	end
c
c
c
