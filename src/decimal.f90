MODULE decimal
  IMPLICIT NONE
  !
  ! Define precision
  !
  INTEGER, PARAMETER :: dp=KIND(1.0d0)
  REAL(kind=dp), PARAMETER :: precision = EPSILON(1.0_dp)
  real(kind=dp), parameter :: PI = 4 * atan (1.0_dp)  
END MODULE decimal
