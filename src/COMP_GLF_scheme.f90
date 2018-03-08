! Component Wise global Lax-Friedrichs scheme
! Based on:
!  Raimund Bürger , Rosa Donat , Pep Mulet , Carlos A. Vega, 
!  On the implementation of WENO schemes for a class of polydisperse sedimentation 
!  models, Journal of Computational Physics, v.230 n.6, p.2322-2344, 
!  March, 2011  [doi>10.1016/j.jcp.2010.12.019]
MODULE CGLF_scheme
USE decimal
USE FVTypes
USE scheme_utils
USE eno_weno
USE limiters

IMPLICIT NONE
  ! Numerical function types
  abstract interface
    function AbstractAlphaF(u, prob) result(alpha)
      import  :: dp
      import  :: CLS1DDiffusionProblem 
      type(CLS1DDiffusionProblem) :: prob
      real(kind = dp), intent(in) :: u(:,:)
      real(kind = dp) :: alpha
    end function AbstractAlphaF
  end interface
  !Algorithm Type
  !
  ! Explicit Component Wise Global Lax Friedrichs algorithm + diffusion
  type, public, extends(FVDiff1DAlgorithm) :: COMP_GLF_1DAlgorithm
    procedure(AbstractAlphaF),pointer,nopass :: alphaf
    real(kind = dp) :: alpha
  contains
    procedure :: update_dt => update_dt_CGLF
    procedure :: update => update_CGLF
    procedure :: getStartMessage => start_message_CGLF
    procedure :: initialize => init_CGLF
  end type COMP_GLF_1DAlgorithm

CONTAINS
  !Not Optional
  subroutine init_CGLF(alg, alphaf)
      CLASS(COMP_GLF_1DAlgorithm)    :: alg
      procedure(AbstractAlphaF) :: alphaf
      !========================
      alg%alphaf => alphaf
      alg%alpha = 0.0_dp
  end subroutine init_CGLF
  function start_message_CGLF(alg) result(message)
      CLASS(COMP_GLF_1DAlgorithm)  :: alg
      CHARACTER(LEN=32)             :: message
      message = "Starting Comp-GLF scheme ..."
  end function start_message_CGLF
  ! Update time based on CFL conditions
  function update_dt_CGLF (alg, u, CFL) result(dt)
    CLASS(COMP_GLF_1DAlgorithm)  :: alg
    REAL(kind = dp), intent(in)  :: u(:,:), CFL
    REAL(kind = dp)              :: dt, dx
    alg%alpha = alg%alphaf(u, alg%problem)
    dx = alg%problem%mesh%dx
    dt = CFL*dx/alg%alpha
  end function update_dt_CGLF

  !!!!!!!!!!! Main methods, used to update solution in each time integration step
  subroutine update_CGLF(alg, rhs, uold, dt)
      CLASS(COMP_GLF_1DAlgorithm)  :: alg
      real(kind = dp), intent(in) :: uold(:,:)
      real(kind = dp)             :: rhs(:,:), dt
      REAL(kind = dp)               :: dx, alpha
      INTEGER                       :: N, j,M, boundary, order, k, i
      REAL(kind = dp), ALLOCATABLE  :: uu(:,:), hh(:,:), pp(:,:)
      REAL(kind = dp), ALLOCATABLE  :: fminus(:,:), fplus(:,:), crj(:,:), Du(:,:)

      !==================
      N = alg%problem%mesh%N
      M = alg%problem%M
      dx = alg%problem%mesh%dx
      boundary = alg%problem%mesh%bdtype
      alpha = alg%alpha
      order = 3; k = 2  !order = k + 1

      ALLOCATE(uu(0:N+1,M), crj(order+1, order))
      ALLOCATE(hh(N+1,M), pp(N+1,M))
      uu = 0.0_dp; hh = 0.0_dp; pp = 0.0_dp; crj = 0.0_dp

      ! Compute Weno coefficients
      crj = unif_crj(order)

      ! Add ghost cells
      uu(1:N,:) = uold
      if (boundary == PERIODIC) then
        uu(0,:) = uold(N,:); uu(N+1,:)=uold(1,:)
      else
        uu(0,:) = uold(1,:); uu(N+1,:)=uold(N,:)
      end if

      ! Global Lax Frierichs splitting
      ALLOCATE(fminus(-k:N+k+1,M), fplus(-k:N+k+1,M))
      fminus = 0.0_dp; fplus = 0.0_dp
      !$omp parallel
      !$omp do
      do j = 1,N
        fminus(j,:) = 0.5_dp*(alg%problem%f(uu(j,:))-alpha*uu(j,:))
        fplus(j,:) = 0.5_dp*(alg%problem%f(uu(j,:))+alpha*uu(j,:))
      end do
      !$omp end do
      !$omp end parallel
      do j = 1,(k+1)
        if (boundary == ZERO_FLUX) then
          fminus(1-j,:) = fminus(1,:); fminus(N+j,:) = fminus(N,:)
          fplus(1-j,:) = fplus(1,:); fplus(N+j,:) = fplus(N,:)
        else !PERIODIC
          fminus(1-j,:) = fminus(N-j+1,:); fminus(N+j,:) = fminus(j,:)
          fplus(1-j,:) = fplus(N-j+1,:); fplus(N+j,:) = fplus(j,:)
        end if
      end do

      ! limit slope
      ALLOCATE(Du(0:N+1,M))
      Du = 0.0_dp
      do i = 1,M
        do j = 1,N
          Du(j,i) = minmod((uu(j,i)-uu(j-1,i)),(uu(j+1,i)-uu(j-1,i))/2,(uu(j+1,i)-uu(j,i)))
        end do
      end do

      ! WENO 5 reconstruction
      !$omp parallel
      !$omp do
      do j = 0,N
        do i = 1,M
          ! Diffusion
          hh(j+1,i) = sum(WENO_pm_rec(fminus(j-k+1:j+k+1,i),fplus(j-k:j+k,i),order*2-1, crj))
          pp(j+1,:) = 0.5*MATMUL(alg%problem%K(uu(j+1,:))+alg%problem%K(uu(j,:)),Du(j,:)/dx)
        end do
      end do
      !$omp end do
      !$omp end parallel

      !Compute Numeric Flux + Diffusion term
      if (boundary == ZERO_FLUX) then
        hh(1,:)=0.0_dp; pp(1,:)=0.0_dp
        hh(N+1,:)=0.0_dp; pp(N+1,:)=0.0_dp
      end if
      do j = 1,N
        rhs(j,:) = - 1/dx * (hh(j+1,:)-hh(j,:)-(pp(j+1,:)-pp(j,:)))
      end do
      DEALLOCATE(uu, hh, pp, fminus, fplus, Du)
  end subroutine update_CGLF

END MODULE CGLF_scheme