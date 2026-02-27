subroutine DLGAMS (X, DLGAM, SGNGAM)
!
!! DLGAMS computes the logarithm of the absolute value of the Gamma function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7A
!***TYPE      DOUBLE PRECISION (ALGAMS-S, DLGAMS-D)
!***KEYWORDS  ABSOLUTE VALUE OF THE LOGARITHM OF THE GAMMA FUNCTION,
!             FNLIB, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DLGAMS(X,DLGAM,SGNGAM) calculates the double precision natural
! logarithm of the absolute value of the Gamma function for
! double precision argument X and stores the result in double
! precision argument DLGAM.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DLNGAM
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  DLGAMS
  DOUBLE PRECISION X, DLGAM, SGNGAM, DLNGAM
!***FIRST EXECUTABLE STATEMENT  DLGAMS
  DLGAM = DLNGAM(X)
  SGNGAM = 1.0D0
  if (X > 0.D0) RETURN
!
  INT = MOD (-AINT(X), 2.0D0) + 0.1D0
  if (INT == 0) SGNGAM = -1.0D0
!
  return
end
