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
  type, public, extends(FVDiff1DAlgorithm) :: ESJP1DAlgorithm
    real(kind=dp)           ::  epsilon = 0.0
    procedure(AbstractNflux),pointer,nopass :: Nf
    procedure(AbstractNMatDiff),pointer,nopass :: Nk
  contains
    procedure :: update => update_ES
    procedure :: getStartMessage => start_message_ES
    procedure :: initialize => init_ES
  end type ESJP1DAlgorithm

  type, public, extends(FVDiff1DAlgorithm) :: ESJPe1DAlgorithm
    real(kind=dp)           ::  epsilon = 0.0
    procedure(AbstractNflux),pointer,nopass :: Nf
    procedure(AbstractNMatDiff),pointer,nopass :: Nk
    procedure(AbstractEntropyV),pointer,nopass :: ve
  contains
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

  subroutine update_ES(alg, rhs, uold, dt, prob)
      CLASS(ESJP1DAlgorithm)  :: alg
      real(kind = dp), intent(in) :: uold(:,:)
      real(kind = dp) :: rhs(:,:), dt
      type(CLS1DDiuffsionProblem) :: prob
      REAL(kind = dp)               :: dx, epsilon

      INTEGER                       :: N, j,M, boundary, ss
      REAL(kind = dp), ALLOCATABLE  :: uleft(:), uright(:), kleft(:,:), kright(:,:)

      !==================
      N = prob%mesh%N
      M = prob%M
      dx = prob%mesh%dx
      epsilon = alg%epsilon
      boundary = prob%mesh%bdtype

      ALLOCATE(uleft(M), uright(M), kleft(M,M), kright(M,M))
      uleft = 0.0_dp; uright = 0.0_dp; kleft = 0.0_dp;kright = 0.0_dp
      uleft = uold(1,:); uright = uold(N,:)
      if (boundary == PERIODIC) then
        uleft = uold(N,:)
        uright = uold(1,:)
      end if

      j = 1
      kleft = alg%Nk(uleft,uold(j,:))
      kright = alg%Nk(uold(j,:),uold(j+1,:))
      rhs(j,:) = -1.0_dp/dx*(alg%Nf(uold(j,:), uold(j+1,:))-alg%Nf(uleft, uold(j,:))) +&
      1.0_dp/dx**2*(MATMUL(kright,(uold(j+1,:)-uold(j,:))) - MATMUL(kleft,(uold(j,:)-uleft))) +&
      epsilon*1.0_dp/dx**2*(uold(j+1,:)-2*uold(j,:)+uleft)
      DO j = 2,(N-1)
      kleft = alg%Nk(uold(j-1,:),uold(j,:))
      kright = alg%Nk(uold(j,:),uold(j+1,:))
        rhs(j,:) = - 1.0_dp/dx*(alg%Nf(uold(j,:), uold(j+1,:))-alg%Nf(uold(j-1,:), uold(j,:))) +&
        1.0_dp/dx**2*(MATMUL(kright,(uold(j+1,:)-uold(j,:))) - MATMUL(kleft,(uold(j,:)-uold(j-1,:))))+&
        epsilon*1.0_dp/dx**2*(uold(j+1,:)-2*uold(j,:)+uold(j-1,:))
      END DO
      j = N
      kleft = alg%Nk(uold(j-1,:),uold(j,:))
      kright = alg%Nk(uold(j,:),uright)
      rhs(j,:) = -1.0_dp/dx*(alg%Nf(uold(j,:), uright)-alg%Nf(uold(j-1,:), uold(j,:))) +&
      1.0_dp/dx**2*(MATMUL(kright,(uright-uold(j,:))) - MATMUL(kleft,(uold(j,:)-uold(j-1,:))))+&
      epsilon*1.0_dp/dx**2*(uright-2*uold(j,:)+uold(j-1,:))
      DEALLOCATE(uleft, uright, kleft, kright)

  end subroutine update_ES

  subroutine update_EES(alg, rhs, uold, dt, prob)
      CLASS(ESJPe1DAlgorithm)  :: alg
      real(kind = dp), intent(in) :: uold(:,:)
      real(kind = dp) :: rhs(:,:), dt
      type(CLS1DDiuffsionProblem) :: prob
      REAL(kind = dp)               :: dx, epsilon

      INTEGER                       :: N, j,M, boundary, ss
      REAL(kind = dp), ALLOCATABLE  :: vold(:,:),vleft(:), vright(:), kleft(:,:), kright(:,:)

      !==================
      N = prob%mesh%N
      M = prob%M
      dx = prob%mesh%dx
      epsilon = alg%epsilon
      boundary = prob%mesh%bdtype

      ALLOCATE(vold(N,M),vleft(M), vright(M), kleft(M,M), kright(M,M))
      vleft = 0.0_dp; vright = 0.0_dp; kleft = 0.0_dp;kright = 0.0_dp
      vleft = uold(1,:); vright = uold(N,:); vold = 0.0_dp

      do j = 1,N
        vold(j,:) = alg%ve(uold(j,:))
      end do

      if (boundary == PERIODIC) then
        vleft = vold(N,:)
        vright = vold(1,:)
      end if

      j = 1
      kleft = alg%Nk(vleft,vold(j,:))
      kright = alg%Nk(vold(j,:),vold(j+1,:))
      rhs(j,:) = -1.0_dp/dx*(alg%Nf(vold(j,:), vold(j+1,:))-alg%Nf(vleft, vold(j,:))) +&
      1.0_dp/dx**2*(MATMUL(kright,(vold(j+1,:)-vold(j,:))) - MATMUL(kleft,(vold(j,:)-vleft))) +&
      epsilon*1.0_dp/dx**2*(vold(j+1,:)-2*vold(j,:)+vleft)
      DO j = 2,(N-1)
      kleft = alg%Nk(vold(j-1,:),vold(j,:))
      kright = alg%Nk(vold(j,:),vold(j+1,:))
        rhs(j,:) = - 1.0_dp/dx*(alg%Nf(vold(j,:), vold(j+1,:))-alg%Nf(vold(j-1,:), vold(j,:))) +&
        1.0_dp/dx**2*(MATMUL(kright,(vold(j+1,:)-vold(j,:))) - MATMUL(kleft,(vold(j,:)-vold(j-1,:))))+&
        epsilon*1.0_dp/dx**2*(vold(j+1,:)-2*vold(j,:)+vold(j-1,:))
      END DO
      j = N
      kleft = alg%Nk(vold(j-1,:),vold(j,:))
      kright = alg%Nk(vold(j,:),vright)
      rhs(j,:) = -1.0_dp/dx*(alg%Nf(vold(j,:), vright)-alg%Nf(vold(j-1,:), vold(j,:))) +&
      1.0_dp/dx**2*(MATMUL(kright,(vright-vold(j,:))) - MATMUL(kleft,(vold(j,:)-vold(j-1,:))))+&
      epsilon*1.0_dp/dx**2*(vright-2*vold(j,:)+vold(j-1,:))
      DEALLOCATE(vold, vleft, vright, kleft, kright)

  end subroutine update_EES

END MODULE EC_scheme
