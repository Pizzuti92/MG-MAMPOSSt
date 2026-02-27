        real*8 function lfdeviate(alpha,xmin,seed)
c
c LFDEVIATE generates random x = L/L* from Schechter LF
c from xmin L* to 10 L*.
c
c compile with -lnr
c
        implicit none
        integer seed
        real*8 alpha, xmin
        real*8 xmax, alpha1, q, ran2
        real ( kind = 8 ) r_uniform_01
c
        if (alpha .eq. -1.d0) then
            print *,' LFDEVIATE: alpha must not equal -1'
            stop
        endif
        alpha1 = 1.d0 + alpha
        xmax = 10.d0

        lfdeviate = -1.d0
        do while (lfdeviate .lt. 0.d0)
            q = r_uniform_01(seed)
            lfdeviate = (q*(xmax**alpha1-xmin**alpha1) 
     &		 + xmin**alpha1)**(1.d0/alpha1)
c
            q = r_uniform_01(seed)
            if (q .gt. exp(-lfdeviate)) then
                lfdeviate = -1.d0
            endif
        enddo
c
        return
        end
c
c
c
      function r_uniform_01 ( seed )

!*****************************************************************************80

!! R8_UNIFORM_01 returns a unit pseudorandom R8.

!  Discussion:

!    This routine implements the recursion

!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )

!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.

!    If the initial seed is 12345, then the first three computations are

!      Input     Output      R8_UNIFORM_01
!      SEED      SEED

!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702

!  Licensing:

!    This code is distributed under the GNU LGPL license.

!  Modified:

!    31 May 2007

!  Author:

!    John Burkardt

!  Reference:

!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.

!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.

!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.

!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, 1969, pages 136-143.

!  Parameters:

!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.
!    On output, SEED has been updated.

!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.

       implicit none

       integer ( kind = 4 ) k
       real ( kind = 8 ) r_uniform_01
       integer ( kind = 4 ) seed

       k = seed / 127773

       seed = 16807 * ( seed - k * 127773 ) - k * 2836

       if ( seed < 0 ) then
        seed = seed + 2147483647
       end if

!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!

       r_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

       return
      end



