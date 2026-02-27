program factorizzazione_check
  implicit none
  real :: r, s, rb, b0, bi
  real :: f, Rr, Ss, prodotto
  real, parameter :: tol = 1.0e-6

  ! Parametri
  rb = 2.0
  b0 = 1.5
  bi = 0.7

  ! Valori di prova
  r = 1.2
  s = 0.8

  ! Calcolo della funzione originale
  f = exp(bi * (-exp((r**2)/(rb**2)) + exp((s**2)/(rb**2)))) &
      * (s/r)**(2.0*b0) &
      * ((rb + s)/(r + rb))**(-2.0*b0 + 2.0*bi)

  ! Calcolo della parte in r
  Rr = exp(-bi * exp((r**2)/(rb**2))) * r**(-2.0*b0) * (r + rb)**(-2.0*bi + 2.0*b0)

  ! Calcolo della parte in s
  Ss = exp(bi * exp((s**2)/(rb**2))) * s**(2.0*b0) * (s + rb)**(-2.0*b0 + 2.0*bi)

  ! Verifica
  prodotto = Rr * Ss

  print*, "f(r,s)        = ", f
  print*, "R(r) * S(s)    = ", prodotto
  if (abs(f - prodotto) < tol) then
     print*, "CHECK OK: la fattorizzazione è corretta entro la tolleranza."
  else
     print*, "ERRORE: la fattorizzazione NON torna!"
  end if

end program factorizzazione_check
