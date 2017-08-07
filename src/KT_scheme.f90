!!!!!!!!!!!!!!!!!!!!!!!!!11 Kurganov Tadmor Centered Scheme
! Complete Discrete form
! Base on:
! Stefan Berres, Raimund Burger, Kenneth H. Karlsen and Elmer M,
! Strongly Degenerate Parabolic-Hyperbolic Systems Modeling Polydisperse
! Sedimentation with Compression
MODULE KT_scheme
USE decimal
USE FVTypes
USE scheme_utils
USE limiters
IMPLICIT NONE
  !Problem Type
  type, public, extends(FVDiff1DAlgorithm) :: KT1DAlgorithm
    real(kind=dp)           ::  theta = 1.0
  contains
    procedure :: update_dt => update_dt_KT
    procedure :: update => update_KT
    procedure :: getStartMessage => start_message_Kt
    procedure :: initialize => init_KT
  end type KT1DAlgorithm

CONTAINS
    !Optional
    subroutine init_KT(alg, theta)
      REAL (kind=dp)    :: theta
      CLASS(KT1DAlgorithm)    :: alg
      !========================
      alg%theta = theta
  end subroutine init_KT
  function start_message_Kt(alg) result(message)
      CLASS(KT1DAlgorithm)  :: alg
      CHARACTER(LEN=32)             :: message
      message = "Starting Kurganov-Tadmor Scheme..."
  end function start_message_Kt
  
  ! Update time based on CFL condition
  function update_dt_KT (alg, u, CFL) result(dt)
    CLASS(KT1DAlgorithm)  :: alg
    REAL(kind = dp), intent(in)  :: u(:,:), CFL
    REAL(kind = dp)              :: dt
    dt = cdtdiff(u, alg%problem, CFL)
  end function update_dt_KT
  
  !!!!!!!!!!! Main methods, used to update solution in each time integration step
  subroutine update_KT(alg, rhs, uold, dt)
      CLASS(KT1DAlgorithm)  :: alg
      real(kind = dp), intent(in) :: uold(:,:)
      real(kind = dp) :: rhs(:,:), dt
      REAL(kind = dp)               :: dx, theta
      REAL(kind = dp), ALLOCATABLE  :: uu(:,:), Du(:,:), uminus(:,:), uplus(:,:), aa(:)
      REAL(kind = dp), ALLOCATABLE  :: u_l(:,:), u_r(:,:), Df_l(:,:), Df_r(:,:), Psir(:,:), Psi(:,:)
      REAL(kind = dp), ALLOCATABLE  :: phi_l(:,:), phi_r(:,:), Fphir(:,:), Fphil(:,:), DPsi(:,:)
      REAL(kind = dp), ALLOCATABLE  :: hh(:,:), pp(:,:), Du_ap(:,:), Ful(:), Fulm(:), Fulp(:), Fur(:), Furm(:), Furp(:)
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
      ! 1. slopes
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
  
        ! Flux slopes
      ALLOCATE(u_l(N+1,M), u_r(N+1,M))
      u_l = 0.0_dp; u_r = 0.0_dp
      do i = 1,M
        do j = 1,(N+1)
          u_l(j,i) = uu(j-1,i) + (0.5-lmb*aa(j))*Du(j-1,i)
          u_r(j,i) = uu(j,i) - (0.5-lmb*aa(j))*Du(j,i)
        end do
      end do
      ALLOCATE(Df_l(N+1,M), Df_r(N+1,M))
      Df_l = 0.0_dp; Df_r = 0.0_dp
      ALLOCATE(Ful(M), Fulm(M), Fulp(M), Fur(M), Furm(M), Furp(M))
      Ful=0.0_dp; Fulm=0.0_dp; Fulp=0.0_dp; Fur=0.0_dp; Furm=0.0_dp; Furp=0.0_dp
      do j = 2,N
        Ful = alg%problem%f(u_l(j,:)); Fulm = alg%problem%f(u_l(j-1,:)); Fulp = alg%problem%f(u_l(j+1,:))
        Fur = alg%problem%f(u_r(j,:)); Furm = alg%problem%f(u_r(j-1,:)); Furp = alg%problem%f(u_r(j+1,:))
        do i = 1, M
          Df_l(j,i) = minmod(theta*(Ful(i)-Fulm(i)),(Fulp(i)-Fulm(i))/2,theta*(Fulp(i)-Ful(i)))
          Df_r(j,i) = minmod(theta*(Fur(i)-Furm(i)),(Furp(i)-Furm(i))/2,theta*(Furp(i)-Fur(i)))
        end do
      end do

      ! Predictor solution values
      ALLOCATE(phi_l(N+1,M), phi_r(N+1,M))
      phi_l = 0.0_dp; phi_r = 0.0_dp
      phi_l = u_l - lmb/2*Df_l
      phi_r = u_r - lmb/2*Df_r

      ! Aproximate cell averages
      ALLOCATE(Psir(N+1,M), Psi(0:N+1,M), Fphir(N+1,M), Fphil(N+1,M))
      Psir = 0.0_dp; Psi = 0.0_dp; Fphir = 0.0_dp; Fphil = 0.0_dp
      do j = 1,N+1
        Fphir(j,:) = alg%problem%f(phi_r(j,:))
        Fphil(j,:) = alg%problem%f(phi_l(j,:))
        if (abs(aa(j)) > 1.0e-6_dp) then
          Psir(j,:) = 0.5*(uu(j-1,:)+uu(j,:))+(1-lmb*aa(j))/4*(Du(j-1,:)-Du(j,:))-1/(2*aa(j))*(Fphir(j,:)-Fphil(j,:))
        else
          Psir(j,:) = 0.5*(uu(j-1,:)+uu(j,:))
        end if
      end do
      do j = 1,N
        Psi(j,:) = uu(j,:) - lmb/2*(aa(j+1)-aa(j))*Du(j,:)-lmb/(1-lmb*(aa(j+1)+aa(j)))*(Fphil(j+1,:)-Fphir(j,:))
      end do
      if (boundary == PERIODIC) then
        Psi(0,:) = Psi(N,:); Psi(N+1,:)=Psi(1,:)
      else
        Psi(0,:) = Psi(1,:); Psi(N+1,:)=Psi(N,:)
      end if

      ! Discrete derivatives
      ALLOCATE(DPsi(N+1,M))
      DPsi = 0.0_dp
      do j = 2,N
        do i = 1,M
          DPsi(j,i)= 2/dx*minmod(theta*(psir(j,i)-Psi(j-1,i))/(1+lmb*(aa(j)-aa(j-1))), &
          (Psi(j,i)-Psi(j-1,i))/(2+lmb*(2*aa(j)-aa(j-1)-aa(j+1))), &
          theta*(Psi(j,i)-psir(j,i))/(1+lmb*(aa(j)-aa(j+1))))
        end do
      end do

      ! Numerical Fluxes
      ALLOCATE(hh(N+1,M), pp(N+1,M), Du_ap(0:N+1,M))
      hh = 0.0_dp; pp = 0.0_dp; Du_ap = 0.0_dp
      do j = 1,(N+1)
        hh(j,:) = 0.5*(Fphir(j,:)+Fphil(j,:))-0.5*(uu(j,:)-uu(j-1,:))*aa(j)+&
        aa(j)*(1-lmb*aa(j))/4*(Du(j,:)+Du(j-1,:)) + lmb*dx/2*(aa(j))**2*DPsi(j,:)
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

      DEALLOCATE(uu, uminus, uplus, aa, u_l, u_r, Df_l, Df_r, phi_l, phi_r, fphir, fphil, DPsi, Du)
      DEALLOCATE(hh, pp, Du_ap, Ful, Fulm, Fulp, Fur, Furm, Furp)
  end subroutine update_KT
END MODULE KT_scheme