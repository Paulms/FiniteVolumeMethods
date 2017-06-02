! Test Problem 2
! Test based on example a of:
! Abeynaike et al, The experimental measurement and modelling of sedimentation
! and creaming for glycerol/biodiesel droplet dispersions, 2012

! 20% glycerol-rich phase

MODULE test2
USE decimal
USE plot
USE FVTypes
USE KT_scheme
USE EC_scheme
IMPLICIT NONE
PUBLIC test2_run
PRIVATE
INTEGER, PARAMETER         :: M = 8
REAL(kind = dp), PARAMETER :: phio(M) = [0.0, 0.006, 0.018, 0.048, 0.08, 0.042, 0.006, 0.0] 
REAL(kind = dp), PARAMETER :: dd(M) = [2.0, 23.0, 34.0, 50.0, 70.0, 100.0, 150.0, 200.0]*1e-6 !m
REAL(kind=dp), PARAMETER   :: phic = 0.0_dp
REAL(kind=dp), PARAMETER   :: rhod = 1090.0_dp !Kg/m^3
REAL(kind=dp), PARAMETER   :: rhoc = 880.0_dp  !Kg/m^3
REAL(kind=dp), PARAMETER   :: gr = 9.81_dp   !m/s^2
REAL(kind=dp), PARAMETER   :: muc = 6.5e-3_dp !Pa/s
REAL(kind = dp), PARAMETER :: ut(M) = (rhod - rhoc)*gr/18*muc*1e4*dd**2 !m/s
REAL(kind =dp), PARAMETER  :: nn = 4.65_dp   ! index of Richardson and Zaki
REAL(kind=dp), PARAMETER   :: Ddo = 1.0e-7_dp  !m/s^2
CONTAINS
subroutine test2_run()
  REAL(kind = dp)           :: Tend
  INTEGER                   :: N
  REAL(kind = dp)           :: dx
  REAL(kind = dp)           :: CFL, L
  REAL(kind = dp), ALLOCATABLE    :: xx(:)
  REAL(kind = dp), ALLOCATABLE    :: uinit(:,:), uu(:,:)
  INTEGER                         :: i
  CHARACTER(LEN=32)               :: name           ! File name to save plot data
  REAL(kind = dp), ALLOCATABLE    :: results(:,:)
  CHARACTER(LEN=8), ALLOCATABLE  :: names(:)
  ! Zero variables
  Tend = 0.0_dp; dx = 0.0_dp; CFL = 0.0_dp
  ! Initialize variables
  Tend = 300.0_dp      ! Final Time
  CFL = 0.3_dp
  L = 20e-3_dp
  !Run numerical schemes
  N = 200          ! Number of nodes
  CALL setup_problem(dx, L, N, M, xx, uu, uinit)
  name = 'test_2_500_3_ini'
  ALLOCATE(results(N, M+1), names(M+1))
  names = ['X       ', 'y1      ','y2      ', 'y3      ', 'y4      ','y5      ','y6      ','y7      ','y8      ']
  results(:,1) = xx
  results(:,2:M+1) = uinit(:,1:M)
  CALL save_matrix(results, names, name, 0)
  !CALL KT(SSPRK22, uu, dx, CFL, Tend, ZERO_FLUX, 1.0_dp, Flux, JacF, BB)
  results(:,2:M+1) = uu(:,1:M)
  name = 'test_2_500_3_kt'
  CALL save_matrix(results, names, name, 0)
  
  !uu = uinit
  !CALL Entropy_NonConservative(SSPRK22, .false., uu, Tend, dx, CFL, NFlux, KKN, JacF, BB, 0.0_dp, PERIODIC)
  !results(:,2:5) = uu(:,1:4)
  !name = 'test_1_800_ec'
  !CALL save_matrix(results, names, name, 0)

  DEALLOCATE(results, uu, names, uinit, xx)
end subroutine test2_run

SUBROUTINE setup_problem(dx, L, N, M, xx, uu, uinit)
  INTEGER, INTENT(IN)             :: N,M
  REAL(kind = dp)                 :: dx, L
  REAL(kind = dp), ALLOCATABLE    :: xx(:)
  REAL(kind = dp), ALLOCATABLE    :: uinit(:,:), uu(:,:)
  INTEGER                         :: i, j
  dx = L/N
    ! Allocate memory
  ALLOCATE(xx(N), uu(N,M))
  xx = 0.0_dp; uu = 0.0_dp
  xx = (/(i*dx+dx/2,i=0, N-1)/)      ! Location of grid points

  ! Initial Conditions
  DO j = 1, N
    uu(j,:) = phio
  END DO
  !Save initial condition
  ALLOCATE(uinit(N,M))
  uinit = 0.0_dp;   
  uinit = uu
END SUBROUTINE   

FUNCTION flux(phi) RESULT(ff)
  REAL(kind = dp), INTENT(IN)  :: phi(:)
  REAL(kind = dp)              :: ff(SIZE(phi,1))
  ff = 0.0_dp
  ff = VV(sum(phi))*phi*ut
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
        F(i,j)=ut(i)*(Vphi + phi(i)*VPphi)
      else
        F(i,j)=ut(i)*(phi(i)*VPphi)
      end if
    end do
  end do
END FUNCTION JacF

function BB(phi) result(k)
  real(kind = dp), intent(in) :: phi(:)
  real(kind = dp) :: k(SIZE(phi,1), SIZE(phi,1))
  REAL(kind = dp) :: sphi
  INTEGER         :: M, i
  k = 0.0_dp
  if (sum(phi) > phic) then
    M = size(phi,1)
    sphi = sum(phi)
    do i = 1,M
     k(i,i) = ddo*(1-sphi)**nn
    end do
  end if
end function BB

function VV(phi)
  real(kind = dp), intent(in) :: phi
  REAL(kind = dp) :: VV
  VV = 1.0_dp
  !if (phi > phic) then
    VV = (1-phi)**nn
  !end if
end function VV

function VP(phi)
  real(kind = dp), intent(in) :: phi
  REAL(kind = dp) :: VP
  VP = 0.0_dp
  !if (phi > phic) then
    VP = -nn*(1-phi)**(nn-1)
  !end if
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
  REAL(kind = dp)              :: ff(SIZE(ul,1))
  INTEGER                      :: M
  M = size(ul,1)
  ff = 0.0_dp
  ff = 0.5*(ur+ul)*0.5*(VV(sum(ur))+VV(sum(ur)))*ut
END FUNCTION NFlux

END MODULE test2