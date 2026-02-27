DOUBLE PRECISION FUNCTION D9LGIT (A, X, ALGAP1)
!
!! D9LGIT computes the logarithm of Tricomi's incomplete Gamma function ...
!  with Perron's continued fraction for large X and A  >=  X.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7E
!***TYPE      DOUBLE PRECISION (R9LGIT-S, D9LGIT-D)
!***KEYWORDS  FNLIB, INCOMPLETE GAMMA FUNCTION, LOGARITHM,
!             PERRON'S CONTINUED FRACTION, SPECIAL FUNCTIONS, TRICOMI
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Compute the log of Tricomi's incomplete gamma function with Perron's
! continued fraction for large X and for A  >=  X.
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
!***END PROLOGUE  D9LGIT
  DOUBLE PRECISION A, X, ALGAP1, AX, A1X, EPS, FK, HSTAR, P, R, S, &
    SQEPS, T, D1MACH
  LOGICAL FIRST
  SAVE EPS, SQEPS, FIRST
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  D9LGIT
  if (FIRST) THEN
     EPS = 0.5D0*D1MACH(3)
     SQEPS = SQRT(D1MACH(4))
  end if
  FIRST = .FALSE.
!
  if (X  <=  0.D0 .OR. A  <  X) call XERMSG ('SLATEC', 'D9LGIT', &
     'X SHOULD BE GT 0.0 AND LE A', 2, 2)
!
  AX = A + X
  A1X = AX + 1.0D0
  R = 0.D0
  P = 1.D0
  S = P
  DO 20 K=1,200
    FK = K
    T = (A+FK)*X*(1.D0+R)
    R = T/((AX+FK)*(A1X+FK)-T)
    P = R*P
    S = S + P
    if (ABS(P) < EPS*S) go to 30
 20   CONTINUE
  call XERMSG ('SLATEC', 'D9LGIT', &
     'NO CONVERGENCE IN 200 TERMS OF CONTINUED FRACTION', 3, 2)
!
 30   HSTAR = 1.0D0 - X*S/A1X
  if (HSTAR  <  SQEPS) call XERMSG ('SLATEC', 'D9LGIT', &
     'RESULT LESS THAN HALF PRECISION', 1, 1)
!
  D9LGIT = -X - ALGAP1 - LOG(HSTAR)
  return
!
end
