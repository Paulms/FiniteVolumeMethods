! Dissipation Reduced Central upwind Scheme: Second-Order
! Based on:
! Kurganov A., Lin C., On the reduction of Numerical Dissipation in Central-Upwind
! Schemes, Commun. Comput. Phys. Vol 2. No. 1, pp 141-163, Feb 2007.

MODULE DRCU_scheme
USE decimal
USE FVTypes
USE scheme_utils
USE limiters
IMPLICIT NONE
  !Problem Type
  type, public, extends(FVDiff1DAlgorithm) :: DRCU1DAlgorithm
    real(kind=dp)           ::  theta = 1.0
  contains
    procedure :: update_dt => update_dt_DRCU
    procedure :: update => update_DRCU
    procedure :: getStartMessage => start_message_DRCU
    procedure :: initialize => init_DRCU
  end type DRCU1DAlgorithm

CONTAINS
    !Optional
    subroutine init_DRCU(alg, theta)
      REAL (kind=dp)    :: theta
      CLASS(DRCU1DAlgorithm)    :: alg
      !========================
      alg%theta = theta
  end subroutine init_DRCU
  function start_message_DRCU(alg) result(message)
      CLASS(DRCU1DAlgorithm)  :: alg
      CHARACTER(LEN=32)             :: message
      message = "Start 2o DR Central-Upwind Scheme..."
  end function start_message_DRCU
  
  ! Update time based on CFL condition
  function update_dt_DRCU (alg, u, CFL) result(dt)
    CLASS(DRCU1DAlgorithm)  :: alg
    REAL(kind = dp), intent(in)  :: u(:,:), CFL
    REAL(kind = dp)              :: dt
    dt = cdtdiff(u, alg%problem, CFL)
  end function update_dt_DRCU
  
  !!!!!!!!!!! Main methods, used to update solution in each time integration step
  subroutine update_DRCU(alg, rhs, uold, dt)
      CLASS(DRCU1DAlgorithm)  :: alg
      real(kind = dp), intent(in) :: uold(:,:)
      real(kind = dp) :: rhs(:,:), dt
      REAL(kind = dp)               :: dx, theta, lambdm(2), lambdp(2)
      REAL(kind = dp), ALLOCATABLE  :: uu(:,:), Du(:,:), uminus(:,:), uplus(:,:), aam(:), aap(:)
      REAL(kind = dp), ALLOCATABLE  :: hh(:,:), pp(:,:), Du_ap(:,:), flm(:), flp(:), wint(:), qj(:)
      INTEGER                       :: N, M, i, j ,k, boundary
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
      ALLOCATE(flm(M), flp(M), wint(M), qj(M))
      ALLOCATE(hh(N+1,M), pp(N+1,M), Du_ap(0:N+1,M))
      hh = 0.0_dp; pp = 0.0_dp; Du_ap = 0.0_dp
      flm = 0.0_dp; flp = 0.0_dp; wint = 0.0_dp; qj = 0.0_dp
      do j = 1,(N+1)
        if (abs(aap(j)-aam(j))<1.0e-8_dp) then
          hh(j,:) = 0.0_dp
        else
          flm = alg%problem%f(uminus(j,:))
          flp = alg%problem%f(uplus(j,:))
          wint = 1/(aap(j)-aam(j))*(aap(j)*uplus(j,:)-aam(j)*uminus(j,:)-(flp-flm))
          do k = 1,M
            qj(k) = minmod2((uplus(j,k)-wint(k))/(aap(j)-aam(j)),(wint(k)-uminus(j,k))/(aap(j)-aam(j)))
          end do 
          hh(j,:) = (aap(j)*flm-aam(j)*flp)/(aap(j)-aam(j)) + &
          (aap(j)*aam(j))*((uplus(j,:) - uminus(j,:))/(aap(j)-aam(j))-qj)
        end if
      end do
      !Du_ap = (uu(1:N+1,:)-uu(0:N,:))/dx
      Du_ap = Du/dx
      ! Diffusion
      do j = 1,N+1
        pp(j,:) = 0.5*MATMUL(alg%problem%K(uu(j,:))+alg%problem%K(uu(j-1,:)),Du_ap(j-1,:))
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
      DEALLOCATE(hh, pp, Du_ap, flm, flp, wint, qj)
  end subroutine update_DRCU
END MODULE DRCU_scheme
