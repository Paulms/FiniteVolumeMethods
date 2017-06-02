! Test Problem 3
! Vehicular traffic problem
! Test based on example 1 of:
! Bürger, Mulet, Villada, A difussion Corrected Multiclass LWR Traffic model
! with Anticipation Lengths and Reaction times, Advances in Applied Mathematics
! and mechanics, 2013
MODULE test3
USE decimal
USE plot
USE FVTypes
USE KT_scheme
USE EC_scheme
USE FV_diffSolve
IMPLICIT NONE
PUBLIC test3_run
PRIVATE
REAL(kind = dp), PARAMETER :: phic = exp(-7.0_dp/exp(1.0_dp))
REAL(kind = dp), PARAMETER :: Vmax(4) = [60.0,55.0,50.0,45.0]
REAL(kind = dp), PARAMETER :: Lmin = 0.03_dp
REAL(kind = dp), PARAMETER :: CC = exp(1.0_dp)/7.0_dp
REAL(kind = dp), PARAMETER :: tau = 1e-3_dp
REAL(kind = dp), PARAMETER :: kappa = 1e-3_dp
REAL(kind = dp), PARAMETER :: wc = 50.0_dp

CONTAINS
subroutine test3_run()
  REAL(kind = dp)           :: Tend
  INTEGER                   :: N, M, bdtype
  REAL(kind = dp)           :: dx
  REAL(kind = dp)           :: CFL, L
  REAL(kind = dp), ALLOCATABLE    :: xx(:)
  REAL(kind = dp), ALLOCATABLE    :: uinit(:,:)
  INTEGER                         :: i
  CHARACTER(LEN=32)               :: name           ! File name to save plot data
  REAL(kind = dp), ALLOCATABLE    :: results(:,:)
  CHARACTER(LEN=8), ALLOCATABLE  :: names(:)
  type(Uniform1DMesh) :: mesh
  type(CLS1DDiffusionProblem) :: prob
  type(KT1DAlgorithm)         :: KtAlg
  type(ESJPe1DAlgorithm)       :: ESPJeAlg
  ! Zero variables
  Tend = 0.0_dp; dx = 0.0_dp; CFL = 0.0_dp
  ! Initialize variables
  Tend = 0.2_dp      ! Final Time
  CFL = 0.25_dp
  L = 10.0
  M = 4
  bdtype = PERIODIC
  !Run numerical schemes
  N = 100          ! Number of nodes
  CALL setup_problem(0.0_dp, L, N, M, mesh, uinit, bdtype)
  CALL prob%Initialize(mesh, uinit, M, Tend, Flux, JacF, BB)
  name = 'test_3_100_ini'
  ALLOCATE(results(N, M+1), names(M+1))
  names = ['X       ', 'y1      ','y2      ', 'y3      ', 'y4      ']
  results(:,1) = mesh%x
  results(:,2:5) = uinit(:,1:4)
  CALL save_matrix(results, names, name, 0)
  CALL solve(prob, SSPRK22, KtAlg, CFL)
  results(:,2:5) = prob%uu(:,1:4)
  name = 'test_3_100_kt'
  CALL save_matrix(results, names, name, 0)

  CALL ESPJeAlg%Initialize(Nfluxe, KKNe, ve, 0.0_dp)  
  CALL solve(prob, SSPRK22, ESPJeAlg, CFL)
  results(:,2:5) = prob%uu(:,1:4)
  name = 'test_3_100_es'
  CALL save_matrix(results, names, name, 0)
  
  DEALLOCATE(results, names, uinit)
end subroutine test3_run

SUBROUTINE setup_problem(xinit, xend, N, M, mesh, uinit, bdtype)
  INTEGER, INTENT(IN)             :: N,M
  REAL(kind = dp)                 :: L, xinit, xend
  REAL(kind = dp), ALLOCATABLE    :: uinit(:,:)
  INTEGER                         :: i, j, bdtype
  type(Uniform1DMesh)             :: mesh
    ! Allocate memory
  call mesh%Initialize(N, xinit, xend, bdtype)
  ALLOCATE(uinit(N,M))
  uinit = 0.0_dp;
  ! Initial Conditions
  DO j = 1, N
    if (0.0 < mesh%x(j) .and. mesh%x(j) <= 0.1) then
      uinit(j,:) = 10.0_dp*mesh%x(j)*[0.2,0.3,0.2,0.3]
    elseif (0.1 < mesh%x(j) .and. mesh%x(j) <=0.9) then
      uinit(j,:) = [0.2,0.3,0.2,0.3]
    elseif (0.9<mesh%x(j) .and. mesh%x(j)<=1) then
      uinit(j,:) = -10.0_dp*(mesh%x(j)-1)*[0.2,0.3,0.2,0.3]
    else
      uinit(j,:) = 0.0_dp
    end if
  END DO
END SUBROUTINE   

FUNCTION flux(phi) RESULT(ff)
  REAL(kind = dp), INTENT(IN)  :: phi(:)
  REAL(kind = dp)              :: ff(SIZE(phi,1))
  ff = 0.0_dp
  ff = VV(sum(phi))*phi*Vmax
END FUNCTION flux

FUNCTION JacF(phi) RESULT(F)
  REAL(kind = dp), INTENT(IN)  :: phi(:)
  REAL(kind = dp)              :: F(SIZE(phi,1),SIZE(phi,1))
  INTEGER                      :: M, i, j
  REAL(kind = dp)              :: Vphi, VPphi
  F = 0.0_dp
  M = SIZE(phi,1)
  Vphi = VV(sum(phi))
  VPphi = VP(sum(phi))
  do i =  1,M
    do j = 1,M
      if (i == j) then
        F(i,j)=Vmax(i)*(Vphi + phi(i)*VPphi)
      else
        F(i,j)=Vmax(i)*(phi(i)*VPphi)
      end if
    end do
  end do
END FUNCTION JacF

function BB(phi) result(k)
  real(kind = dp), intent(in) :: phi(:)
  real(kind = dp) :: k(SIZE(phi,1), SIZE(phi,1))
  REAL(kind = dp) :: VPphi, phiVm
  INTEGER         :: M, i
  k = 0.0_dp
  if (sum(phi) > phic) then
    M = size(phi,1)
    VPphi = VP(sum(phi))
    phiVm = sum(phi*Vmax)/M
    do i = 1,M
     k(i,i) = -VPphi*Lmin*phiVm
    end do
  end if
end function BB

function VV(phi)
  real(kind = dp), intent(in) :: phi
  REAL(kind = dp) :: VV
  VV = 1.0_dp
  if (phi > phic) then
    VV = -CC*log(phi)
  end if
end function VV

function VP(phi)
  real(kind = dp), intent(in) :: phi
  REAL(kind = dp) :: VP
  VP = 0.0_dp
  if (phi > phic) then
    VP = -CC*1/phi
  end if
end function VP

! Functions for EC_scheme
!k(u) numerico
FUNCTION KKN(ul, ur) RESULT(kk)
  !Numerical viscosity matrix
  REAL(kind = dp), INTENT(IN)  :: ul(:)
  REAL(kind = dp), INTENT(IN)  :: ur(:)
  REAL(kind = dp)              :: kk(SIZE(ul,1),SIZE(ul,1))
  kk = 0.0_dp
  !kk = 0.5*(BB(ul)+BB(ur))
END FUNCTION KKN

FUNCTION NFlux(ul, ur) RESULT(ff)
  !Entropy conservative flux
  REAL(kind = dp), INTENT(IN)  :: ul(:)
  REAL(kind = dp), INTENT(IN)  :: ur(:)
  REAL(kind = dp)              :: ff(SIZE(ul,1)), V, sul, sur
  REAL(kind = dp)              :: wl(SIZE(ul,1)), wr(SIZE(ul,1))
  REAL(kind = dp)              :: lul, lur
  INTEGER                      :: M,i
  M = size(ul,1)
  ff = 0.0_dp; wl = 0.0_dp; wr = 0.0_dp; lul = 0.0_dp; lur = 0.0_dp
  V = 0.0;sul = 0.0_dp; sur = 0.0_dp
  wl = vl(ul); wr = vl(ur)
  sul = sum(ul); sur = sum(ur)
  lul = log(sul); lur = log(sur)
  if (sul <= phic .and. sur <= phic) then
    V = (sur-sul)
  else
    V = - CC*((sur*lur-sul*lul) - (sur-sul))
  end if
  do i = 1,M
    if (abs(ur(i)-ul(i)) > kappa) then
      ff(i) = V*(Vmax(i)/(wr(i)-wl(i)))
    else
      ff(i) = VV(0.5*(sul+sur))*Vmax(i)*0.5*(ul(i)+ur(i))
    end if
  end do
END FUNCTION NFlux

! functions in entropy variables
FUNCTION vl(phi)
  REAL(kind = dp), INTENT(IN) :: phi(:)
  REAL(kind = dp)             :: vl(SIZE(phi,1))
  INTEGER                     :: i, M
  vl = 0.0_dp
  M = size(phi,1)
  do i=1,M
    if (phi(i) < tau) then
      vl(i) = -wc
    else
      vl(i) = max(log(phi(i)),-wc)
    end if
  end do
END FUNCTION vl
FUNCTION ve(phi)
  REAL(kind = dp), INTENT(IN)  :: phi(:)
  REAL(kind = dp)              :: ve(SIZE(phi,1))
  ve = 0.0_dp
  ve = vl(phi)/Vmax
END FUNCTION ve
!k(u) numerico
FUNCTION KKNe(wl, wr) RESULT(kk)
  !Numerical viscosity matrix
  REAL(kind = dp), INTENT(IN)  :: wl(:)
  REAL(kind = dp), INTENT(IN)  :: wr(:)
  REAL(kind = dp)              :: kk(SIZE(wl,1),SIZE(wl,1))
  INTEGER                      :: M, i
  REAL(kind = dp)              :: meanW(SIZE(wl,1)), VPw, wVm

  kk = 0.0_dp
  meanW = 0.0_dp
  if (sum(exp(Vmax*wl)) > phic .and. sum(exp(Vmax*wr)) > phic) then
    M = size(wl,1)
    meanW = 0.5*((exp(Vmax*wl))+(exp(Vmax*wr)))
    VPw = VP(sum(meanW))
    wVm = sum(meanW*Vmax)/M
    do i = 1,M
     kk(i,i) = -VPw*Lmin*wVm
    end do
  end if
END FUNCTION KKNe

FUNCTION NFluxe(wl, wr) RESULT(ff)
  !Entropy conservative flux
  REAL(kind = dp), INTENT(IN)  :: wl(:)
  REAL(kind = dp), INTENT(IN)  :: wr(:)
  REAL(kind = dp)              :: ff(SIZE(wl,1))
  ff = 0.0_dp
  ff = Nflux(exp(Vmax*wl),exp(Vmax*wr))
END FUNCTION NFluxe

END MODULE test3