function INITDS (OS, NOS, ETA)
!
!! INITDS determines the number of terms needed in an orthogonal ...
!            polynomial series so that it meets a specified accuracy.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C3A2
!***TYPE      DOUBLE PRECISION (INITS-S, INITDS-D)
!***KEYWORDS  CHEBYSHEV, FNLIB, INITIALIZE, ORTHOGONAL POLYNOMIAL,
!             ORTHOGONAL SERIES, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
!  Initialize the orthogonal series, represented by the array OS, so
!  that INITDS is the number of terms needed to insure the error is no
!  larger than ETA.  Ordinarily, ETA will be chosen to be one-tenth
!  machine precision.
!
!             Input Arguments --
!   OS     double precision array of NOS coefficients in an orthogonal
!          series.
!   NOS    number of coefficients in OS.
!   ETA    single precision scalar containing requested accuracy of
!          series.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891115  Modified error message.  (WRB)
!   891115  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  INITDS
  DOUBLE PRECISION OS(*)
!***FIRST EXECUTABLE STATEMENT  INITDS
  if (NOS  <  1) call XERMSG ('SLATEC', 'INITDS', &
     'Number of coefficients is less than 1', 2, 1)
!
  ERR = 0.
  DO 10 II = 1,NOS
    I = NOS + 1 - II
    ERR = ERR + ABS(REAL(OS(I)))
    if (ERR > ETA) go to 20
   10 CONTINUE
!
   20 if (I  ==  NOS) call XERMSG ('SLATEC', 'INITDS', &
     'Chebyshev series too short for specified accuracy', 1, 1)
  INITDS = I
!
  return
end
