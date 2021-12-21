DOUBLE PRECISION FUNCTION D9LGIC (A, X, ALX)
!
!! D9LGIC computes the log complementary incomplete Gamma function ...
!            for large X and for A  <=  X.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7E
!***TYPE      DOUBLE PRECISION (R9LGIC-S, D9LGIC-D)
!***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, LARGE X,
!             LOGARITHM, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Compute the log complementary incomplete gamma function for large X
! and for A  <=  X.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900720  Routine changed from user-callable to subsidiary.  (WRB)
!***END PROLOGUE  D9LGIC
  DOUBLE PRECISION A, X, ALX, EPS, FK, P, R, S, T, XMA, XPA, D1MACH
  SAVE EPS
  DATA EPS / 0.D0 /
!***FIRST EXECUTABLE STATEMENT  D9LGIC
  if (EPS == 0.D0) EPS = 0.5D0*D1MACH(3)
!
  XPA = X + 1.0D0 - A
  XMA = X - 1.D0 - A
!
  R = 0.D0
  P = 1.D0
  S = P
  DO 10 K=1,300
    FK = K
    T = FK*(A-FK)*(1.D0+R)
    R = -T/((XMA+2.D0*FK)*(XPA+2.D0*FK)+T)
    P = R*P
    S = S + P
    if (ABS(P) < EPS*S) go to 20
 10   CONTINUE
  call XERMSG ('SLATEC', 'D9LGIC', &
     'NO CONVERGENCE IN 300 TERMS OF CONTINUED FRACTION', 1, 2)
!
 20   D9LGIC = A*ALX - X + LOG(S/XPA)
!
  return
end
