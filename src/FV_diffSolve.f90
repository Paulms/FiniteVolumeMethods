MODULE FV_diffSolve
USE decimal
USE FVTypes
USE scheme_utils

IMPLICIT NONE
CONTAINS

subroutine solve(problem, time_scheme, algorithm, CFL)
  INTEGER                   :: time_scheme
  INTEGER                   :: N, M, percentage
  REAL(kind = dp)           :: dt, limit, Tend, CFL
  REAL(kind = dp)           :: tt
  REAL(kind = dp), ALLOCATABLE  :: tmp(:,:), uold(:,:),rhs(:,:),rhs2(:,:)
  REAL(kind=dp)             :: tiempo1, tiempo2
  type(CLS1DDiuffsionProblem) :: problem
  CLASS(FVDiff1DAlgorithm)     :: algorithm
  CHARACTER(LEN=32)             :: message

  !Get variables
  Tend = problem%Tend
  N = problem%mesh%N
  M = problem%M
  dt = 0.0_dp
  tt = 0.0_dp
  ! Setup progress message
  tiempo1 = 0.0_dp;tiempo2 = 0.0_dp
  percentage = 0
  limit = Tend/5
  !Init variables
  ALLOCATE(uold(N,M), tmp(N,M))
  ALLOCATE(rhs(N,M),rhs2(N,M))

  uold = 0.0_dp; rhs = 0.0_dp; rhs2 = 0.0_dp
  tmp = 0.0_dp
  message = algorithm%GetStartMessage()
  Print *, message
  CALL cpu_TIME(tiempo1)
  tt = 0.0_dp
  problem%uu = problem%u0
   DO WHILE (tt <= Tend)
    uold = problem%uu
    dt = cdt(uold, problem, CFL)
    !Scheme with forward Euler
    IF (time_scheme == FORWARD_EULER) THEN
      CALL algorithm%update(rhs, uold, dt,problem)
      problem%uu = uold + dt*rhs
      !Scheme with SSPRK22
    ELSE IF (time_scheme == SSPRK22) THEN
      !FIRST Step
      CALL algorithm%update(rhs, uold, dt,problem)
      problem%uu = 0.5*(uold + dt*rhs)
      !Second Step
      CALL algorithm%update(rhs, uold+dt*rhs, dt,problem)
      problem%uu = problem%uu + 0.5*(uold + dt*rhs)
    ELSE IF (time_scheme == SSPRRK33) THEN
      CALL algorithm%update(rhs, uold, dt,problem)
      tmp = uold + dt*rhs
      CALL algorithm%update(rhs, tmp, dt,problem)
      tmp = (3*uold + tmp + dt*rhs) / 4
      CALL algorithm%update(rhs, tmp, dt/2.0,problem)
      problem%uu = (uold + 2*tmp + 2*dt*rhs) / 3
    ELSE
      print *, "ERROR time integraton method not available..."
      stop
    END IF
    !Print progress
    IF (tt > limit) THEN
      percentage = percentage + 20
      limit = limit + Tend/5
      print *, percentage, "% completed"
    END IF
    tt = tt + dt
  END DO
  print *, "completed..."
  CALL cpu_TIME(tiempo2)
  PRINT*,'tiempo de CPU = ',tiempo2-tiempo1
  DEALLOCATE(uold, rhs, tmp,rhs2)
end subroutine solve

END MODULE FV_diffSolve