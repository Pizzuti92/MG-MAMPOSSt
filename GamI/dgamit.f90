  DOUBLE PRECISION FUNCTION DGAMIT (A, X)
!
!! DGAMIT calculates Tricomi's form of the incomplete Gamma function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7E
!***TYPE      DOUBLE PRECISION (GAMIT-S, DGAMIT-D)
!***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB,
!             SPECIAL FUNCTIONS, TRICOMI
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
!   Evaluate Tricomi's incomplete Gamma function defined by
!
!   DGAMIT = X**(-A)/GAMMA(A) * integral from 0 to X of EXP(-T) *
!              T**(A-1.)
!
!   for A  >  0.0 and by analytic continuation for A  <=  0.0.
!   GAMMA(X) is the complete gamma function of X.
!
!   DGAMIT is evaluated for arbitrary real values of A and for non-
!   negative values of X (even though DGAMIT is defined for X  <
!   0.0), except that for X = 0 and A  <=  0.0, DGAMIT is infinite,
!   which is a fatal error.
!
!   The function and both arguments are DOUBLE PRECISION.
!
!   A slight deterioration of 2 or 3 digits accuracy will occur when
!   DGAMIT is very large or very small in absolute value, because log-
!   arithmic variables are used.  Also, if the parameter  A  is very
!   close to a negative integer (but not a negative integer), there is
!   a loss of accuracy, which is reported if the result is less than
!   half machine precision.
!
!***REFERENCES  W. Gautschi, A computational procedure for incomplete
!                 gamma functions, ACM Transactions on Mathematical
!                 Software 5, 4 (December 1979), pp. 466-481.
!               W. Gautschi, Incomplete gamma functions, Algorithm 542,
!                 ACM Transactions on Mathematical Software 5, 4
!                 (December 1979), pp. 482-489.
!***ROUTINES CALLED  D1MACH, D9GMIT, D9LGIC, D9LGIT, DGAMR, DLGAMS,
!                    DLNGAM, XERCLR, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
!***END PROLOGUE  DGAMIT
  DOUBLE PRECISION A, X, AEPS, AINTA, ALGAP1, ALNEPS, ALNG, ALX, &
    BOT, H, SGA, SGNGAM, SQEPS, T, D1MACH, DGAMR, D9GMIT, D9LGIT, &
    DLNGAM, D9LGIC
  LOGICAL FIRST
  SAVE ALNEPS, SQEPS, BOT, FIRST
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DGAMIT
  if (FIRST) THEN
     ALNEPS = -LOG (D1MACH(3))
     SQEPS = SQRT(D1MACH(4))
     BOT = LOG (D1MACH(1))
  end if
  FIRST = .FALSE.
!
  if (X  <  0.D0) call XERMSG ('SLATEC', 'DGAMIT', 'X IS NEGATIVE' &
     , 2, 2)
!
  if (X /= 0.D0) ALX = LOG (X)
  SGA = 1.0D0
  if (A /= 0.D0) SGA = SIGN (1.0D0, A)
  AINTA = AINT (A + 0.5D0*SGA)
  AEPS = A - AINTA
!
  if (X > 0.D0) go to 20
  DGAMIT = 0.0D0
  if (AINTA > 0.D0 .OR. AEPS /= 0.D0) DGAMIT = DGAMR(A+1.0D0)
  return
!
 20   if (X > 1.D0) go to 30
  if (A >= (-0.5D0) .OR. AEPS /= 0.D0) call DLGAMS (A+1.0D0, ALGAP1, &
    SGNGAM)
  DGAMIT = D9GMIT (A, X, ALGAP1, SGNGAM, ALX)
  return
!
 30   if (A < X) go to 40
  T = D9LGIT (A, X, DLNGAM(A+1.0D0))
  if (T < BOT) call XERCLR
  DGAMIT = EXP (T)
  return
!
 40   ALNG = D9LGIC (A, X, ALX)
!
! EVALUATE DGAMIT IN TERMS OF LOG (DGAMIC (A, X))
!
  H = 1.0D0
  if (AEPS == 0.D0 .AND. AINTA <= 0.D0) go to 50
!
  call DLGAMS (A+1.0D0, ALGAP1, SGNGAM)
  T = LOG (ABS(A)) + ALNG - ALGAP1
  if (T > ALNEPS) go to 60
!
  if (T > (-ALNEPS)) H = 1.0D0 - SGA * SGNGAM * EXP(T)
  if (ABS(H) > SQEPS) go to 50
!
  call XERCLR
  call XERMSG ('SLATEC', 'DGAMIT', 'RESULT LT HALF PRECISION', 1, &
     1)
!
 50   T = -A*ALX + LOG(ABS(H))
  if (T < BOT) call XERCLR
  DGAMIT = SIGN (EXP(T), H)
  return
!
 60   T = T - A*ALX
  if (T < BOT) call XERCLR
  DGAMIT = -SGA * SGNGAM * EXP(T)
  return
!
end
