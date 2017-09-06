! Second-Order upwind central scheme
! Based on:
! Kurganov A., Noelle S., Petrova G., Semidiscrete Central-Upwind schemes
! for hyperbolic Conservation Laws and Hamilton-Jacobi Equations. SIAM. Sci Comput,
! Vol 23, No 3m pp 707-740. 2001

MODULE CU_scheme
USE decimal
USE FVTypes
USE scheme_utils
USE limiters
IMPLICIT NONE
  !Problem Type
  type, public, extends(FVDiff1DAlgorithm) :: CU1DAlgorithm
    real(kind=dp)           ::  theta = 1.0
  contains
    procedure :: update_dt => update_dt_CU
    procedure :: update => update_CU
    procedure :: getStartMessage => start_message_CU
    procedure :: initialize => init_CU
  end type CU1DAlgorithm

CONTAINS
    !Optional
    subroutine init_CU(alg, theta)
      REAL (kind=dp)    :: theta
      CLASS(CU1DAlgorithm)    :: alg
      !========================
      alg%theta = theta
  end subroutine init_CU
  function start_message_CU(alg) result(message)
      CLASS(CU1DAlgorithm)  :: alg
      CHARACTER(LEN=32)             :: message
      message = "Start 2o Central-Upwind Scheme..."
  end function start_message_CU
  
  ! Update time based on CFL condition
  function update_dt_CU (alg, u, CFL) result(dt)
    CLASS(CU1DAlgorithm)  :: alg
    REAL(kind = dp), intent(in)  :: u(:,:), CFL
    REAL(kind = dp)              :: dt
    dt = cdtdiff(u, alg%problem, CFL)
  end function update_dt_CU
  
  !!!!!!!!!!! Main methods, used to update solution in each time integration step
  subroutine update_CU(alg, rhs, uold, dt)
      CLASS(CU1DAlgorithm)  :: alg
      real(kind = dp), intent(in) :: uold(:,:)
      real(kind = dp) :: rhs(:,:), dt
      REAL(kind = dp)               :: dx, theta, lambdm(2), lambdp(2)
      REAL(kind = dp), ALLOCATABLE  :: uu(:,:), Du(:,:), uminus(:,:), uplus(:,:), aam(:), aap(:)
      REAL(kind = dp), ALLOCATABLE  :: hh(:,:), pp(:,:), Du_ap(:,:)
      INTEGER                       :: N, M, i, j ,boundary
      REAL(kind = dp)               :: lmb
      !==================
      N = alg%problem%mesh%N
      M = alg%problem%M
      dx = alg%problem%mesh%dx
      theta = alg%theta
      boundary = alg%problem%mesh%bdtype

      !Add ghost cells
      ALLOCATE(uu(0:N+1,M))
      uu = 0.0_dp
      uu(1:N,:) = uold
      if (boundary == PERIODIC) then
        uu(0,:) = uold(N,:); uu(N+1,:)=uold(1,:)
      else
        uu(0,:) = uold(1,:); uu(N+1,:)=uold(N,:)
      end if

      !Compute diffusion
      lmb = dt/dx

      ! update vector
      ! 1. Reconstruct approximate derivatives
      ALLOCATE(Du(0:N+1,M), uminus(N+1,M), uplus(N+1,M))
      Du = 0.0_dp; uminus = 0.0_dp; uplus = 0.0_dp
      do i = 1,M
        do j = 1,N
          Du(j,i) = minmod(theta*(uu(j,i)-uu(j-1,i)),(uu(j+1,i)-uu(j-1,i))/2,theta*(uu(j+1,i)-uu(j,i)))
        end do
      end do
      if (boundary == PERIODIC) then
        Du(0,:) = Du(N,:); Du(N+1,:)=Du(1,:)
      else
        Du(0,:) = Du(1,:); Du(N+1,:)=Du(N,:)
      end if
      ! Local speeds of propagation
      uminus = uu(0:N,:)+0.5*Du(0:N,:)
      uplus = uu(1:N+1,:)-0.5*Du(1:N+1,:)
      ALLOCATE(aam(N+1),aap(N+1))
      aam = 0.0_dp; aap = 0.0_dp; lambdm = 0.0_dp; lambdp = 0.0_dp
      do j = 1,(N+1)
        lambdm = minmaxl(uminus(j,:), alg%problem%Jf)
        lambdp = minmaxl(uplus(j,:), alg%problem%Jf)
        aam(j)=min(0.0_dp,min(lambdm(1), lambdp(1)))
        aap(j)=max(0.0_dp,max(lambdm(2), lambdp(2)))
      end do
  
      ! Numerical Fluxes
      ALLOCATE(hh(N+1,M), pp(N+1,M), Du_ap(0:N+1,M))
      hh = 0.0_dp; pp = 0.0_dp; Du_ap = 0.0_dp
      do j = 1,(N+1)
        if (abs(aap(j)-aam(j))<1.0e-8_dp) then
          hh(j,:) = 0.5_dp*(alg%problem%f(uminus(j,:))+alg%problem%f(uplus(j,:)))
        else
          hh(j,:) = (aap(j)*alg%problem%f(uminus(j,:))-aam(j)*alg%problem%f(uplus(j,:)))/(aap(j)-aam(j)) + &
          (aap(j)*aam(j))/(aap(j)-aam(j))*(uplus(j,:) - uminus(j,:))
        end if
      end do
      !Du_ap = (uu(1:N+1,:)-uu(0:N,:))/dx
      Du_ap = Du/dx
      ! Diffusion
      do j = 1,N+1
        pp(j,:) = 0.5*MATMUL(alg%problem%K(uu(j,:))+alg%problem%K(uu(j-1,:)),Du_ap(j,:))
      end do

      !Compute Numeric Flux + Diffusion term
      if (boundary == ZERO_FLUX) then
        hh(1,:)=0.0; pp(1,:)=0.0
        hh(N+1,:)=0.0; pp(N+1,:)=0.0
      end if
      do j = 1,N
        rhs(j,:) = - 1/dx * (hh(j+1,:)-hh(j,:)-(pp(j+1,:)-pp(j,:)))
      end do

      DEALLOCATE(uu, uminus, uplus, aam, aap, Du)
      DEALLOCATE(hh, pp, Du_ap)
  end subroutine update_CU
END MODULE CU_scheme