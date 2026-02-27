	integer function binser(array,value,imax)
c
c BINSER performs a binary search of 'value' on the increase-sorted arra
c 'array' of size 'imax'. BINSER returns the index of the first array el
c whose value is larger or equal to 'value'.
c
	implicit none
	integer imax
	real*8 array(imax), value
	integer i1, i2, i
c
	if (array(imax) .lt. value) then
		binser = imax+1
		return
	elseif (array(1) .gt. value) then
		binser = 0
		return
	endif
	i1 = 1
	i2 = imax
c
111	continue
	i = (i1+i2)/2
	if (array(i) .ge. value) then
		if (array(i-1) .lt. value) then
			binser = i
			return
		endif
		i2 = i
		go to 111
	else
		i1 = i
		go to 111
	endif
c
	return
	end
c
c
c

		
