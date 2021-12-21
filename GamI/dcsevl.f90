  DOUBLE PRECISION FUNCTION DCSEVL (X, CS, N)
!
!! DCSEVL evaluates a Chebyshev series.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C3A2
!***TYPE      DOUBLE PRECISION (CSEVL-S, DCSEVL-D)
!***KEYWORDS  CHEBYSHEV SERIES, FNLIB, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
!  Evaluate the N-term Chebyshev series CS at X.  Adapted from
!  a method presented in the paper by Broucke referenced below.
!
!       Input Arguments --
!  X    value at which the series is to be evaluated.
!  CS   array of N terms of a Chebyshev series.  In evaluating
!       CS, only half the first coefficient is summed.
!  N    number of terms in array CS.
!
!***REFERENCES  R. Broucke, Ten subroutines for the manipulation of
!                 Chebyshev series, Algorithm 446, Communications of
!                 the A.C.M. 16, (1973) pp. 254-256.
!               L. Fox and I. B. Parker, Chebyshev Polynomials in
!                 Numerical Analysis, Oxford University Press, 1968,
!                 page 56.
!***ROUTINES CALLED  D1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900329  Prologued revised extensively and code rewritten to allow
!           X to be slightly outside interval (-1,+1).  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DCSEVL
  DOUBLE PRECISION B0, B1, B2, CS(*), ONEPL, TWOX, X, D1MACH
  LOGICAL FIRST
  SAVE FIRST, ONEPL
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DCSEVL
  if (FIRST) ONEPL = 1.0D0 + D1MACH(4)
  FIRST = .FALSE.
  if (N  <  1) call XERMSG ('SLATEC', 'DCSEVL', &
     'NUMBER OF TERMS  <=  0', 2, 2)
  if (N  >  1000) call XERMSG ('SLATEC', 'DCSEVL', &
     'NUMBER OF TERMS  >  1000', 3, 2)
  if (ABS(X)  >  ONEPL) call XERMSG ('SLATEC', 'DCSEVL', &
     'X OUTSIDE THE INTERVAL (-1,+1)', 1, 1)
!
  B1 = 0.0D0
  B0 = 0.0D0
  TWOX = 2.0D0*X
  DO 10 I = 1,N
     B2 = B1
     B1 = B0
     NI = N + 1 - I
     B0 = TWOX*B1 - B2 + CS(NI)
   10 CONTINUE
!
  DCSEVL = 0.5D0*(B0-B2)
!
  return
end
