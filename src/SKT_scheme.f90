! Semidiscrete KT Central-Scheme: Second-Order
! Based on:
! Kurganov, Tadmor. New High Resolution Central Schemes for Non Linear Conservation
! Laws and Convection-Difussion Equations. Journal of Comp Physics 160, pp 241-282. 2000

MODULE SKT_scheme
USE decimal
USE FVTypes
USE scheme_utils
USE limiters
IMPLICIT NONE
  !Problem Type
  type, public, extends(FVDiff1DAlgorithm) :: SKT1DAlgorithm
    real(kind=dp)           ::  theta = 1.0
  contains
    procedure :: update_dt => update_dt_SKT
    procedure :: update => update_SKT
    procedure :: getStartMessage => start_message_SKT
    procedure :: initialize => init_SKT
  end type SKT1DAlgorithm

CONTAINS
    !Optional
    subroutine init_SKT(alg, theta)
      REAL (kind=dp)    :: theta
      CLASS(SKT1DAlgorithm)    :: alg
      !========================
      alg%theta = theta
  end subroutine init_SKT
  function start_message_SKT(alg) result(message)
      CLASS(SKT1DAlgorithm)  :: alg
      CHARACTER(LEN=32)             :: message
      message = "Starting Kurganov-Tadmor Scheme..."
  end function start_message_SKT
  
  ! Update time based on CFL condition
  function update_dt_SKT (alg, u, CFL) result(dt)
    CLASS(SKT1DAlgorithm)  :: alg
    REAL(kind = dp), intent(in)  :: u(:,:), CFL
    REAL(kind = dp)              :: dt
    dt = cdtdiff(u, alg%problem, CFL)
  end function update_dt_SKT
  
  !!!!!!!!!!! Main methods, used to update solution in each time integration step
  subroutine update_SKT(alg, rhs, uold, dt)
      CLASS(SKT1DAlgorithm)  :: alg
      real(kind = dp), intent(in) :: uold(:,:)
      real(kind = dp) :: rhs(:,:), dt
      REAL(kind = dp)               :: dx, theta
      REAL(kind = dp), ALLOCATABLE  :: uu(:,:), Du(:,:), uminus(:,:), uplus(:,:), aa(:)
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
      ALLOCATE(aa(N+1))
      aa = 0.0_dp
      do j = 1,(N+1)
        aa(j)=max(fluxp(uminus(j,:), alg%problem%Jf),fluxp(uplus(j,:), alg%problem%Jf))
      end do
  
      ! Numerical Fluxes
      ALLOCATE(hh(N+1,M), pp(N+1,M), Du_ap(0:N+1,M))
      hh = 0.0_dp; pp = 0.0_dp; Du_ap = 0.0_dp
      do j = 1,(N+1)
        hh(j,:) = 0.5*(alg%problem%f(uplus(j,:))+alg%problem%f(uminus(j,:))) - aa(j)/2*(uplus(j,:) - uminus(j,:))
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
        hh(N,:)=0.0; pp(N,:)=0.0
      end if
      do j = 1,N
        rhs(j,:) = - 1/dx * (hh(j+1,:)-hh(j,:)-(pp(j+1,:)-pp(j,:)))
      end do

      DEALLOCATE(uu, uminus, uplus, aa, Du)
      DEALLOCATE(hh, pp, Du_ap)
  end subroutine update_SKT
END MODULE SKT_scheme