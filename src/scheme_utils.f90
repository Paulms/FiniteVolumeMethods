MODULE scheme_utils
USE decimal
USE FVTypes
IMPLICIT NONE
CONTAINS

function cdtdiff(u, prob, CFL) result(dt)
  type(CLS1DDiffusionProblem) :: prob
  REAL(kind = dp), intent(in)  :: u(:,:), CFL
  REAL(kind = dp)              :: dt, dx
  REAL(kind = dp)              :: maxp, maxpB
  character jobvl, jobvr
  real (kind = dp)              :: WR(size(u,2)), WL(size(u,2)), work(1000)
  real (kind = dp)              :: VL(size(u,2),size(u,2)), VR(size(u,2),size(u,2))
  REAL(kind = dp), ALLOCATABLE :: B(:,:), J(:,:)
  INTEGER :: M, N, info, ldvl, ldvr, lwork, i
  ! =========================
  dx = prob%mesh%dx
  maxp = 0
  maxpB = 0
  N = size(u,1)
  M = size(u,2)
  jobvl = 'N'
  jobvr = 'N'
  LDVL = M; LDVR = M; lwork = -1
  WR = 0.0_dp; WL = 0.0_dp; VL = 0.0_dp; VR = 0.0_dp
  ALLOCATE(B(M,M), J(M,M))
  J = 0.0_dp; B = 0.0_dp
  do i = 1,N
    J = prob%Jf(u(i,:))
    call dgeev (jobvl, jobvr, M, J, M, WR, WL, VL, LDVL, VR, LDVR, WORK, 4*M, INFO)
    maxp = max(maxp, maxval(abs(WR)))
    B = prob%K(u(i,:))
    call dgeev (jobvl, jobvr, M, B, M, WR, WL, VL, LDVL, VR, LDVR, WORK, 4*M, INFO)
    maxpB = max(maxpB, maxval(abs(WR)))
    IF (info /= 0) THEN
      PRINT *, "Error while computing eigenvalues"
      STOP
    END IF
  end do
  dt = CFL/(1/dx*maxp+1/(2*dx**2)*maxpB)
  DEALLOCATE(B, J)
end function cdtdiff

function  fluxp(u, JacF) result(pp)
  REAL(kind = dp), intent(in)  :: u(:)
  REAL(kind = dp) :: pp
  character jobvl, jobvr
  real (kind = dp)              :: WR(size(u,1)), WL(size(u,1)), work(1000)
  real (kind = dp)              :: VL(size(u,1),size(u,1)), VR(size(u,1),size(u,1))
  REAL(kind = dp), ALLOCATABLE :: J(:,:)
  INTEGER :: M, N, info, ldvl, ldvr, lwork

    interface
      function JacF(phi) result(Jf)
          import  :: dp
          real(kind = dp), intent(in) :: phi(:)
          real(kind = dp) :: Jf(SIZE(phi,1),SIZE(phi,1))
      end function JacF
    end interface

  M = size(u,1)
  jobvl = 'N'
  jobvr = 'N'
  LDVL = M; LDVR = M; lwork = -1
  WR = 0.0_dp; WL = 0.0_dp; VL = 0.0_dp; VR = 0.0_dp
  WORK = 0.0_dp
  ALLOCATE(J(M,M))
  J = 0.0_dp
  J = JacF(u)
  call dgeev (jobvl, jobvr, M, J, M, WR, WL, VL, LDVL, VR, LDVR, WORK, 4*M, INFO)
  pp = maxval(abs(WR))
  DEALLOCATE(J)
end function fluxp
END MODULE scheme_utils