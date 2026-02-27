  DOUBLE PRECISION FUNCTION DGAMR (X)
!
!! DGAMR computes the reciprocal of the Gamma function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7A
!***TYPE      DOUBLE PRECISION (GAMR-S, DGAMR-D, CGAMR-C)
!***KEYWORDS  FNLIB, RECIPROCAL GAMMA FUNCTION, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DGAMR(X) calculates the double precision reciprocal of the
! complete Gamma function for double precision argument X.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DGAMMA, DLGAMS, XERCLR, XGETF, XSETF
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900727  Added EXTERNAL statement.  (WRB)
!***END PROLOGUE  DGAMR
  DOUBLE PRECISION X, ALNGX, SGNGX, DGAMMA
  EXTERNAL DGAMMA
!***FIRST EXECUTABLE STATEMENT  DGAMR
  DGAMR = 0.0D0
  if (X <= 0.0D0 .AND. AINT(X) == X) RETURN
!
  call XGETF (IROLD)
  call XSETF (1)
  if (ABS(X) > 10.0D0) go to 10
  DGAMR = 1.0D0/DGAMMA(X)
  call XERCLR
  call XSETF (IROLD)
  return
!
 10   call DLGAMS (X, ALNGX, SGNGX)
  call XERCLR
  call XSETF (IROLD)
  DGAMR = SGNGX * EXP(-ALNGX)
  return
!
end
