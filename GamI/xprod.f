	subroutine xprod(x, y, cross)
c
c XPROD calculates the cross-product of vectors x and y.
c
	real*8 x(3), y(3), cross(3)
c
	cross(1) = x(2)*y(3) - x(3)*y(2)
	cross(2) = x(3)*y(1) - x(1)*y(3)
	cross(3) = x(1)*y(2) - x(2)*y(1)
c
	return
	end
c
c
c
