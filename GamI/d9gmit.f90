DOUBLE PRECISION FUNCTION D9GMIT (A, X, ALGAP1, SGNGAM, ALX)
!
!! D9GMIT computes Tricomi's incomplete Gamma function for small arguments.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7E
!***TYPE      DOUBLE PRECISION (R9GMIT-S, D9GMIT-D)
!***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, SMALL X,
!             SPECIAL FUNCTIONS, TRICOMI
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Compute Tricomi's incomplete gamma function for small X.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DLNGAM, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   890911  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900720  Routine changed from user-callable to subsidiary.  (WRB)
!***END PROLOGUE  D9GMIT
  DOUBLE PRECISION A, X, ALGAP1, SGNGAM, ALX, AE, AEPS, ALGS, ALG2, &
    BOT, EPS, FK, S, SGNG2, T, TE, D1MACH, DLNGAM
  LOGICAL FIRST
  SAVE EPS, BOT, FIRST
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  D9GMIT
  if (FIRST) THEN
     EPS = 0.5D0*D1MACH(3)
     BOT = LOG (D1MACH(1))
  end if
  FIRST = .FALSE.
!
  if (X  <=  0.D0) call XERMSG ('SLATEC', 'D9GMIT', &
     'X SHOULD BE GT 0', 1, 2)
!
  MA = A + 0.5D0
  if (A < 0.D0) MA = A - 0.5D0
  AEPS = A - MA
!
  AE = A
  if (A < (-0.5D0)) AE = AEPS
!
  T = 1.D0
  TE = AE
  S = T
  DO 20 K=1,200
    FK = K
    TE = -X*TE/FK
    T = TE/(AE+FK)
    S = S + T
    if (ABS(T) < EPS*ABS(S)) go to 30
 20   CONTINUE
  call XERMSG ('SLATEC', 'D9GMIT', &
     'NO CONVERGENCE IN 200 TERMS OF TAYLOR-S SERIES', 2, 2)
!
 30   if (A >= (-0.5D0)) ALGS = -ALGAP1 + LOG(S)
  if (A >= (-0.5D0)) go to 60
!
  ALGS = -DLNGAM(1.D0+AEPS) + LOG(S)
  S = 1.0D0
  M = -MA - 1
  if (M == 0) go to 50
  T = 1.0D0
  DO 40 K=1,M
    T = X*T/(AEPS-(M+1-K))
    S = S + T
    if (ABS(T) < EPS*ABS(S)) go to 50
 40   CONTINUE
!
 50   D9GMIT = 0.0D0
  ALGS = -MA*LOG(X) + ALGS
  if (S == 0.D0 .OR. AEPS == 0.D0) go to 60
!
  SGNG2 = SGNGAM * SIGN (1.0D0, S)
  ALG2 = -X - ALGAP1 + LOG(ABS(S))
!
  if (ALG2 > BOT) D9GMIT = SGNG2 * EXP(ALG2)
  if (ALGS > BOT) D9GMIT = D9GMIT + EXP(ALGS)
  return
!
 60   D9GMIT = EXP (ALGS)
  return
!
end
