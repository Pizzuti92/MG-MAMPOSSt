program main

  use root_module, wp => root_module_rk

  implicit none

  real(wp) :: x, f
  integer :: iflag

  call root_scalar('bisection',func,-9.0_wp,31.0_wp,x,f,iflag)

  write(*,*) 'f(',x,') = ', f
  write(*,*) 'iflag    = ', iflag

contains

  function func(x) result(f)

  implicit none

  real(wp),intent(in) :: x
  real(wp) :: f

  f = -200.0_wp * x * exp(-3.0_wp*x)

  end function func

end program main



c       subroutine zero_rc ( a, b, t, arg, status, value )

c!*****************************************************************************80
c!
c!! zero_rc() seeks a root of a function F(X) using reverse communication.
c!
c!  Discussion:
c!
c!    The interval [A,B] must be a change of sign interval for F.
c!    That is, F(A) and F(B) must be of opposite signs.  Then
c!    assuming that F is continuous implies the existence of at least
c!    one value C between A and B for which F(C) = 0.
c!
c!    The location of the zero is determined to within an accuracy
c!    of 4 * EPSILON * abs ( C ) + 2 * T.
c!
c!    The routine is a revised version of the Brent zero finder 
c!    algorithm, using reverse communication.
c!
c!    Thanks to Thomas Secretin for pointing out a transcription error in the
c!    setting of the value of P, 11 February 2013.
c!
c!  Licensing:
c!
c!    This code is distributed under the GNU LGPL license. 
c!
c!  Modified:
c!
c!    13 July 2021
c!
c!  Author:
c!
c!    John Burkardt
c!
c!  Reference:
c!
c!    Richard Brent,
c!    Algorithms for Minimization Without Derivatives,
c!    Dover, 2002,
c!    ISBN: 0-486-41998-3,
c!    LC: QA402.5.B74.
c!
c!  Input:
c!
c!    real ( kind = rk ) A, B, the endpoints of the change of sign interval.
c!
c!    real ( kind = rk ) T, a positive error tolerance.
c!
c!    integer STATUS, used to communicate between the user 
c!    and the routine.  The user only sets STATUS to zero on the first call, 
c!    to indicate that this is a startup call.
c!
c!    real ( kind = rk ) VALUE, the function value at ARG, as requested
c!    by the routine on the previous call.
c!
c!  Output:
c!
c!    real ( kind = rk ) ARG, the currently considered point.  For the next call,
c!    the user is requested to evaluate the function at ARG, and return
c!    the value in VALUE.  On return with STATUS zero, ARG is the routine's
c!    estimate for the function's zero.
c!
c!    integer STATUS, used to communicate between the user 
c!    and the routine.  The routine returns STATUS positive to request 
c!    that the function be evaluated at ARG, or returns STATUS as 0, to 
c!    indicate that the iteration is complete and that ARG is the estimated zero.
c!
c       implicit none

c       integer, parameter :: rk = kind ( 1.0D+00 )

c       real ( kind = rk ) a
c       real ( kind = rk ) arg
c       real ( kind = rk ) b
c       real ( kind = rk ), save :: c
c       real ( kind = rk ), save :: d
c       real ( kind = rk ), save :: e
c       real ( kind = rk ), save :: fa
c       real ( kind = rk ), save :: fb
c       real ( kind = rk ), save :: fc
c       real ( kind = rk ) m
c       real ( kind = rk ) p
c       real ( kind = rk ) q
c       real ( kind = rk ) r
c       real ( kind = rk ) s
c       real ( kind = rk ), save :: sa
c       real ( kind = rk ), save :: sb
c        integer status
c       real ( kind = rk ) t
c       real ( kind = rk ) tol
c       real ( kind = rk ) value
c!
c!  Input STATUS = 0.
c!  Initialize, request F(A).
c!
c       if ( status == 0 ) then

c        sa = a
c        sb = b
c        e = sb - sa
c        d = e

c        status = 1
c        arg = a
c        return
c!
c!  Input STATUS = 1.
c!  Receive F(A), request F(B).
c!
c        else if ( status == 1 ) then

c        fa = value

c        status = 2
c        arg = sb
c        return
c!
c!  Input STATUS = 2
c!  Receive F(B).
c!
c        else if ( status == 2 ) then

c            fb = value

c            if ( 0.0D+00 < fa * fb ) then
c            status = -1
c            return
c            end if

c            c = sa
c            fc = fa

c        else

c        fb = value

c        if ( ( 0.0D+00 < fb .and. 0.0D+00 < fc ) .or. &
c            ( fb <= 0.0D+00 .and. fc <= 0.0D+00 ) ) then
c        c = sa
c        fc = fa
c        e = sb - sa
c        d = e
c        end if

c        end if
c!
c!  Compute the next point at which a function value is requested.
c!
c        if ( abs ( fc ) < abs ( fb ) ) then

c        sa = sb
c        sb = c
c        c = sa
c        fa = fb
c        fb = fc
c        fc = fa

c        end if

c        tol = 2.0D+00 * epsilon ( sb ) * abs ( sb ) + t
c        m = 0.5D+00 * ( c - sb )

c        if ( abs ( m ) <= tol .or. fb == 0.0D+00 ) then
c            status = 0
c            arg = sb
c            return
c        end if

c        if ( abs ( e ) < tol .or. abs ( fa ) <= abs ( fb ) ) then

c            e = m
c            d = e

c        else

c            s = fb / fa

c            if ( sa == c ) then

c                p = 2.0D+00 * m * s
c                q = 1.0D+00 - s

c            else

c                q = fa / fc
c                r = fb / fc
c                p = s * ( 2.0D+00 *m*q*(q-r) -(sb-sa)*(r-1.0D+00 ) )
c              q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )

c            end if

c            if ( 0.0D+00 < p ) then
c                q = - q
c            else
c                p = - p
c            end if

c            s = e
c            e = d

c        if ( 2.0D+00 * p < 3.0D+00 * m * q - abs ( tol * q ) .and. &
c            p < abs ( 0.5D+00 * s * q ) ) then
c            d = p / q
c        else
c            e = m
c            d = e
c        end if

c        end if

c        sa = sb
c        fa = fb

c        if ( tol < abs ( d ) ) then
c            sb = sb + d
c        else if ( 0.0D+00 < m ) then
c            sb = sb + tol
c        else
c            sb = sb - tol
c        end if

c            arg = sb
c        status = status + 1

c        return
c        end
