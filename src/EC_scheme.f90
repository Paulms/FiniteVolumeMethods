! Based on
! Jerez, Pares, Entropy stable schemes for degenerate convection-difusion
! equations. 2017. Society for Industrial and Applied Mathematics. SIAM. Vol. 55.
! No. 1. pp. 240-264
MODULE EC_scheme
USE decimal
USE FVTypes
USE scheme_utils
IMPLICIT NONE
  ! Numerical function types
  abstract interface
    function AbstractNflux(ul, ur) result(f)
        import  :: dp
        real(kind = dp), intent(in) :: ul(:), ur(:)
        real(kind = dp) :: f(SIZE(ul,1))
    end function AbstractNflux
    function AbstractNMatDiff(ul, ur) result(k)
        import :: dp
        real(kind = dp), intent(in) :: ul(:), ur(:)
        real(kind = dp) :: k(SIZE(ul,1), SIZE(ul,1))
    end function AbstractNMatDiff
    function AbstractEntropyV(u) result(w)
      import  :: dp
      real(kind = dp), intent(in) :: u(:)
      real(kind = dp) :: w(SIZE(u,1))
    end function AbstractEntropyV
  end interface
  !Algorithm Type
  ! Pares Scheme using original variables
  type, public, extends(FVDiff1DAlgorithm) :: ESJP1DAlgorithm
    real(kind=dp)           ::  epsilon = 0.0
    procedure(AbstractNflux),pointer,nopass :: Nf
    procedure(AbstractNMatDiff),pointer,nopass :: Nk
  contains
    procedure :: update_dt => update_dt_ES
    procedure :: update => update_ES
    procedure :: getStartMessage => start_message_ES
    procedure :: initialize => init_ES
  end type ESJP1DAlgorithm
  ! Pares Scheme using entropy variables
  type, public, extends(FVDiff1DAlgorithm) :: ESJPe1DAlgorithm
    real(kind=dp)           ::  epsilon = 0.0
    procedure(AbstractNflux),pointer,nopass :: Nf
    procedure(AbstractNMatDiff),pointer,nopass :: Nk
    procedure(AbstractEntropyV),pointer,nopass :: ve
  contains
    procedure :: update_dt => update_dt_EES
    procedure :: update => update_EES
    procedure :: getStartMessage => start_message_EES
    procedure :: initialize => init_EES
  end type ESJPe1DAlgorithm

CONTAINS
    !Not Optional
    subroutine init_ES(alg, Nflux, NDiffMat, epsilon)
      CLASS(ESJP1DAlgorithm)    :: alg
      procedure(AbstractNFLux)  :: NFlux
      procedure(AbstractNMatDiff) :: NDiffMat
      real(kind=dp)           ::epsilon
      !========================
      alg%epsilon = epsilon
      alg%Nf => NFlux
      alg%Nk => NDiffMat
  end subroutine init_ES
  subroutine init_EES(alg, Nflux, NDiffMat, entropyV, epsilon)
      CLASS(ESJPe1DAlgorithm)    :: alg
      procedure(AbstractNFLux)  :: NFlux
      procedure(AbstractNMatDiff) :: NDiffMat
      procedure(AbstractEntropyV) :: entropyV
      real(kind=dp)           ::epsilon
      !========================
      alg%epsilon = epsilon
      alg%Nf => NFlux
      alg%Nk => NDiffMat
      alg%ve => entropyV
  end subroutine init_EES
  function start_message_ES(alg) result(message)
      CLASS(ESJP1DAlgorithm)  :: alg
      CHARACTER(LEN=32)             :: message
      message = "Starting Entropy Stable - Non conservative diffusion ..."
  end function start_message_ES
  function start_message_EES(alg) result(message)
      CLASS(ESJPe1DAlgorithm)  :: alg
      CHARACTER(LEN=32)             :: message
      message = "Starting Entropy Stable (ev) - Non conservative diffusion ..."
  end function start_message_EES
  
  ! Update time based on CFL condition
  function update_dt_ES (alg, u, CFL) result(dt)
    CLASS(ESJP1DAlgorithm)  :: alg
    REAL(kind = dp), intent(in)  :: u(:,:), CFL
    REAL(kind = dp)              :: dt
    dt = cdtdiff(u, alg%problem, CFL)
  end function update_dt_ES
    function update_dt_EES (alg, u, CFL) result(dt)
    CLASS(ESJPe1DAlgorithm)  :: alg
    REAL(kind = dp), intent(in)  :: u(:,:), CFL
    REAL(kind = dp)              :: dt
    dt = cdtdiff(u, alg%problem, CFL)
  end function update_dt_EES

  !!!!!!!!!!! Main methods, used to update solution in each time integration step
    subroutine update_ES(alg, rhs, uold, dt)
      CLASS(ESJP1DAlgorithm)  :: alg
      real(kind = dp), intent(in) :: uold(:,:)
      real(kind = dp)             :: rhs(:,:), dt
      REAL(kind = dp)               :: dx, epsilon
      INTEGER                       :: N, j,M, boundary
      REAL(kind = dp), ALLOCATABLE  :: uu(:,:), hh(:,:), pp(:,:)
      !==================
      N = alg%problem%mesh%N
      M = alg%problem%M
      dx = alg%problem%mesh%dx
      epsilon = alg%epsilon
      boundary = alg%problem%mesh%bdtype

      ALLOCATE(uu(0:N+1,M))
      ALLOCATE(hh(N+1,M), pp(N+1,M))
      uu = 0.0_dp; hh = 0.0_dp; pp = 0.0_dp

      uu(1:N,:) = uold
      if (boundary == PERIODIC) then
        uu(0,:) = uold(N,:); uu(N+1,:)=uold(1,:)
      else
        uu(0,:) = uold(1,:); uu(N+1,:)=uold(N,:)
      end if

      ! Numerical Fluxes
      do j = 1,(N+1)
        hh(j,:) = alg%Nf(uu(j-1,:), uu(j,:))
      end do
      ! Diffusion
      do j = 1,(N+1)
        pp(j,:) = 1/dx*MATMUL(alg%Nk(uu(j-1,:), uu(j,:)),(uu(j,:)-uu(j-1,:)))+&
        epsilon*1.0_dp/dx*(uu(j+1,:)-2*uu(j,:)+uu(j-1,:))
      end do

      !Compute Numeric Flux + Diffusion term
      if (boundary == ZERO_FLUX) then
        hh(1,:)=0.0; pp(1,:)=0.0
        hh(N+1,:)=0.0; pp(N+1,:)=0.0
      end if
      do j = 1,N
        rhs(j,:) = - 1/dx * (hh(j+1,:)-hh(j,:)-(pp(j+1,:)-pp(j,:)))
      end do
      DEALLOCATE(uu, hh, pp)

  end subroutine update_ES

    subroutine update_EES(alg, rhs, uold, dt)
      CLASS(ESJPe1DAlgorithm)  :: alg
      real(kind = dp), intent(in) :: uold(:,:)
      real(kind = dp)             :: rhs(:,:), dt
      REAL(kind = dp)               :: dx, epsilon
      INTEGER                       :: N, j,M, boundary
      REAL(kind = dp), ALLOCATABLE  :: uu(:,:), vv(:,:), hh(:,:), pp(:,:)
      !==================
      N = alg%problem%mesh%N
      M = alg%problem%M
      dx = alg%problem%mesh%dx
      epsilon = alg%epsilon
      boundary = alg%problem%mesh%bdtype

      ALLOCATE(uu(0:N+1,M), vv(0:N+1,M))
      ALLOCATE(hh(N+1,M), pp(N+1,M))
      uu = 0.0_dp; vv = 0.0_dp; hh = 0.0_dp; pp = 0.0_dp

      uu(1:N,:) = uold
      if (boundary == PERIODIC) then
        uu(0,:) = uold(N,:); uu(N+1,:)=uold(1,:)
      else
        uu(0,:) = uold(1,:); uu(N+1,:)=uold(N,:)
      end if

      do j = 0,(N+1)
        vv(j,:) = alg%ve(uu(j,:))
      end do

      ! Numerical Fluxes
      do j = 1,(N+1)
        hh(j,:) = alg%Nf(vv(j-1,:), vv(j,:))
      end do
      ! Diffusion
      do j = 1,(N+1)
        pp(j,:) = 1/dx*MATMUL(alg%Nk(vv(j-1,:), vv(j,:)),vv(j,:)-vv(j-1,:))+&
        epsilon*1.0_dp/dx*(vv(j,:)-vv(j-1,:))
      end do

      !Compute Numeric Flux + Diffusion term
      if (boundary == ZERO_FLUX) then
        hh(1,:)=0.0; pp(1,:)=0.0
        hh(N,:)=0.0; pp(N,:)=0.0
      end if
      do j = 1,N
        rhs(j,:) = - 1/dx * (hh(j+1,:)-hh(j,:)-(pp(j+1,:)-pp(j,:)))
      end do
      DEALLOCATE(uu, vv, hh, pp)

  end subroutine update_EES
END MODULE EC_scheme
