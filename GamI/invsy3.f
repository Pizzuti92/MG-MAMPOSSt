	subroutine invsy3(aa, yy, xx)
c
c INVSY3 solves the system of 3 linear equations Ax=y, where A is a
c real symetric 3-3 matrix:
c
c      ~        * a(1)   a(4)   a(5) *           * d1   a1   a2 *
c      A   =   *  a(4)   a(2)   a(6)  *    =    *  a1   d2   a3  *
c               * a(5)   a(6)   a(3) *           * a2   a3   d3 *
c
	implicit none
	real*8 aa(6), yy(3), xx(3)
	real*8 a1, a2, a3, d1, d2, d3, y1, y2, y3
	real*8 del, del1, del2, del3
c
	a1 = aa(4)
	a2 = aa(5)
	a3 = aa(6)
	d1 = aa(1)
	d2 = aa(2)
	d3 = aa(3)
	y1 = yy(1)
	y2 = yy(2)
	y3 = yy(3)
c
	del =d1*d2*d3+a1*a2*a3+a1*a2*a3 -d1*a3*a3 -d2*a2*a2-d3*a1*a1
	del1 = y1*d2*d3+y2*a2*a3+y3*a1*a3-y1*a3*a3-d2*a2*y3-d3*a1*y2
	del2 = d1*y2*d3+y1*a2*a3+a2*a1*y3-d1*y3*a3-y2*a2*a2-d3*y1*a1
	del3 = d1*d2*y3+a1*y1*a3+a2*a1*y2-d1*y2*a3-d2*y1*a2-y3*a1*a1
	xx(1) = del1/del
	xx(2) = del2/del
	xx(3) = del3/del
c
	return
	end
c
c
c
