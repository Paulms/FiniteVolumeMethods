! Lax-Friedrichs scheme
MODULE LF_scheme
USE decimal
USE FVTypes
USE scheme_utils
IMPLICIT NONE
  ! Numerical function types
  abstract interface
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
  !
  type, public, extends(FVDiff1DAlgorithm) :: LLFDiff1DAlgorithm
    procedure(AbstractNMatDiff),pointer,nopass :: Nk
    procedure(AbstractEntropyV),pointer,nopass :: ve
  contains
    procedure :: update_dt => update_dt_LLFd
    procedure :: update => update_LLFd
    procedure :: getStartMessage => start_message_LLFd
    procedure :: initialize => init_LLFd
  end type LLFDiff1DAlgorithm

CONTAINS
    !Not Optional
  subroutine init_LLFd(alg, NDiffMat, entropyV)
      CLASS(LLFDiff1DAlgorithm)    :: alg
      procedure(AbstractNMatDiff) :: NDiffMat
      procedure(AbstractEntropyV) :: entropyV
      !========================
      alg%Nk => NDiffMat
      alg%ve => entropyV
  end subroutine init_LLFd
  function start_message_LLFd(alg) result(message)
      CLASS(LLFDiff1DAlgorithm)  :: alg
      CHARACTER(LEN=32)             :: message
      message = "Starting Local Lax Friedrichs + Non conservative diffusion ..."
  end function start_message_LLFd
  ! Update time based on CFL conditions
  function update_dt_LLFd (alg, u, CFL) result(dt)
    CLASS(LLFDiff1DAlgorithm)  :: alg
    REAL(kind = dp), intent(in)  :: u(:,:), CFL
    REAL(kind = dp)              :: dt
    dt = cdtdiff(u, alg%problem, CFL)
  end function update_dt_LLFd

  !!!!!!!!!!! Main methods, used to update solution in each time integration step
  subroutine update_LLFd(alg, rhs, uold, dt)
      CLASS(LLFDiff1DAlgorithm)  :: alg
      real(kind = dp), intent(in) :: uold(:,:)
      real(kind = dp)             :: rhs(:,:), dt, alphar, alphal
      REAL(kind = dp)               :: dx
      INTEGER                       :: N, j,M, boundary, ss
      REAL(kind = dp), ALLOCATABLE  :: uu(:,:), vv(:,:), hh(:,:), pp(:,:), alphak(:)
      !==================
      N = alg%problem%mesh%N
      M = alg%problem%M
      dx = alg%problem%mesh%dx
      boundary = alg%problem%mesh%bdtype

      ALLOCATE(uu(0:N+1,M), vv(0:N+1,M), alphak(N+1))
      ALLOCATE(hh(N+1,M), pp(N+1,M))
      uu = 0.0_dp; vv = 0.0_dp; hh = 0.0_dp; pp = 0.0_dp; alphak = 0.0_dp

      uu(1:N,:) = uold
      if (boundary == PERIODIC) then
        uu(0,:) = uold(N,:); uu(N+1,:)=uold(1,:)
      else
        uu(0,:) = uold(1,:); uu(N+1,:)=uold(N,:)
      end if

      do j = 0,(N+1)
        vv(j,:) = alg%ve(uu(j,:))
      end do

      ! Local Lax Friedrichs for convex flux
      alphal = fluxp(uu(0,:), alg%problem%Jf)
      do j = 1,(N+1)
        alphar = fluxp(uu(j,:), alg%problem%Jf)
        alphak(j) = MAX(alphal, alphar)
        alphal = alphar
      end do
      ! Numerical Fluxes
      do j = 1,(N+1)
        hh(j,:) = 0.5*(alg%problem%f(uu(j-1,:))+alg%problem%f(uu(j,:))-alphak(j)*(uu(j,:)-uu(j-1,:)))
      end do
      ! Diffusion
      do j = 1,(N+1)
        pp(j,:) = 1/dx*MATMUL(alg%Nk(vv(j-1,:), vv(j,:)),vv(j,:)-vv(j-1,:))
      end do

      !Compute Numeric Flux + Diffusion term
      if (boundary == ZERO_FLUX) then
        hh(1,:)=0.0; pp(1,:)=0.0
        hh(N,:)=0.0; pp(N,:)=0.0
      end if
      do j = 1,N
        rhs(j,:) = - 1/dx * (hh(j+1,:)-hh(j,:)-(pp(j+1,:)-pp(j,:)))
      end do
      DEALLOCATE(uu, vv, hh, pp, alphak)

  end subroutine update_LLFd

END MODULE LF_scheme