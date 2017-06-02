 !!!!!!!!!!!!!!!!!!!!!!!!!11 Kurganov Tadmor Centered Scheme
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
    procedure :: update => update_KT
    procedure :: getStartMessage => start_message_Kt
    procedure :: initialize => init_KT
  end type KT1DAlgorithm

CONTAINS
    !Optional
    subroutine init_KT(alg, theta)
      CLASS(KT1DAlgorithm)    :: alg
      real(kind=dp)           ::theta
      !========================
      alg%theta = theta
  end subroutine init_KT
  function start_message_Kt(alg) result(message)
      CLASS(KT1DAlgorithm)  :: alg
      CHARACTER(LEN=32)             :: message
      message = "Starting Kurganov-Tadmor Scheme..."
  end function start_message_Kt

  subroutine update_KT(alg, rhs, uold, dt, prob)
      CLASS(KT1DAlgorithm)  :: alg
      real(kind = dp), intent(in) :: uold(:,:)
      real(kind = dp) :: rhs(:,:), dt
      type(CLS1DDiffusionProblem) :: prob
      REAL(kind = dp)               :: dx, theta
      REAL(kind = dp), ALLOCATABLE  :: uo(:,:), Du(:,:), uminus(:,:), uplus(:,:), aa(:)
      REAL(kind = dp), ALLOCATABLE  :: u_l(:,:), u_r(:,:), Df_l(:,:), Df_r(:,:), Psir(:,:), Psi(:,:)
      REAL(kind = dp), ALLOCATABLE  :: phi_l(:,:), phi_r(:,:), Fphir(:,:), Fphil(:,:), DPsi(:,:)
      REAL(kind = dp), ALLOCATABLE  :: hh(:,:), pp(:,:), Du_ap(:,:), Ful(:), Fulm(:), Fulp(:), Fur(:), Furm(:), Furp(:)
      REAL(kind = dp), ALLOCATABLE  :: hhleft(:), hhright(:), ppleft(:), ppright(:)
      INTEGER                       :: N, M, i, j ,boundary
      INTEGER                       :: ss
      REAL(kind = dp)               :: lmb
      !==================
      N = prob%mesh%N
      M = prob%M
      dx = prob%mesh%dx
      theta = alg%theta
      boundary = prob%mesh%bdtype

      ss = 0
      IF (boundary == PERIODIC) THEN
        ! Ghost Cells
        ss = 1
        N = N + 2
        ALLOCATE(uo(N,M))
        uo(2:N-1,:) = uold(:,:)
        uo(1,:) = uold(N-2,:)
        uo(N,:) = uold(1,:)
      ELSE
        ALLOCATE(uo(N,M))
        uo = 0.0_dp
        uo = uold
      END IF

      !Compute diffusion
      lmb = dt/dx

      ! update vector
      ! 1. slopes
      ALLOCATE(Du(N,M), uminus(N-1,M), uplus(N-1,M))
      Du = 0.0_dp; uminus = 0.0_dp; uplus = 0.0_dp
      do i = 1,M
        do j = 2,(N-1)
          Du(j,i) = minmod(theta*(uo(j,i)-uo(j-1,i)),(uo(j+1,i)-uo(j-1,i))/2,theta*(uo(j+1,i)-uo(j,i)))
        end do
      end do
      ! Local speeds of propagation
      uminus = uo(1:N-1,:)+0.5*Du(1:N-1,:)
      uplus = uo(2:N,:)-0.5*Du(2:N,:)
      ALLOCATE(aa(N-1))
      aa = 0.0_dp
      do j = 1,(N-1)
        aa(j)=max(fluxp(uminus(j,:), prob%Jf),fluxp(uplus(j,:), prob%Jf))
      end do
  
        ! Flux slopes
      ALLOCATE(u_l(N-1,M), u_r(N-1,M))
      u_l = 0.0_dp; u_r = 0.0_dp
      do i = 1,M
        do j = 2,N
          u_l(j-1,i) = uo(j-1,i) + (0.5-lmb*aa(j-1))*Du(j-1,i)
          u_r(j-1,i) = uo(j,i) - (0.5-lmb*aa(j-1))*Du(j,i)
        end do
      end do
      ALLOCATE(Df_l(N-1,M), Df_r(N-1,M))
      Df_l = 0.0_dp; Df_r = 0.0_dp
      ALLOCATE(Ful(M), Fulm(M), Fulp(M), Fur(M), Furm(M), Furp(M))
      Ful=0.0_dp; Fulm=0.0_dp; Fulp=0.0_dp; Fur=0.0_dp; Furm=0.0_dp; Furp=0.0_dp
      do j = 2,(N-2)
        Ful = prob%f(u_l(j,:)); Fulm = prob%f(u_l(j-1,:)); Fulp = prob%f(u_l(j+1,:))
        Fur = prob%f(u_r(j,:)); Furm = prob%f(u_r(j-1,:)); Furp = prob%f(u_r(j+1,:))
        do i = 1, M
          Df_l(j,i) = minmod(theta*(Ful(i)-Fulm(i)),(Fulp(i)-Fulm(i))/2,theta*(Fulp(i)-Ful(i)))
          Df_r(j,i) = minmod(theta*(Fur(i)-Furm(i)),(Furp(i)-Furm(i))/2,theta*(Furp(i)-Fur(i)))
        end do
      end do

      ! Predictor solution values
      ALLOCATE(phi_l(N-1,M), phi_r(N-1,M))
      phi_l = 0.0_dp; phi_r = 0.0_dp
      phi_l = u_l - lmb/2*Df_l
      phi_r = u_r - lmb/2*Df_r

      ! Aproximate cell averages
      ALLOCATE(Psir(N-1,M), Psi(N,M), Fphir(N-1,M), Fphil(N-1,M))
      Psir = 0.0_dp; Psi = 0.0_dp; Fphir = 0.0_dp; Fphil = 0.0_dp
      do j = 1,(N-1)
        Fphir(j,:) = prob%f(phi_r(j,:))
        Fphil(j,:) = prob%f(phi_l(j,:))
        if (abs(aa(j)) > 1.0e-6_dp) then
          Psir(j,:) = 0.5*(uo(j,:)+uo(j+1,:))+(1-lmb*aa(j))/4*(Du(j,:)-Du(j+1,:))-1/(2*aa(j))*(Fphir(j,:)-Fphil(j,:))
        else
          Psir(j,:) = 0.5*(uo(j,:)+uo(j+1,:))
        end if
      end do
      do j = 2,(N-1)
        Psi(j,:) = uo(j,:) - lmb/2*(aa(j)-aa(j-1))*Du(j,:)-lmb/(1-lmb*(aa(j)+aa(j-1)))*(Fphil(j,:)-Fphir(j-1,:))
      end do

      ! Discrete derivatives
      ALLOCATE(DPsi(N-1,M))
      DPsi = 0.0_dp
      do j = 2,(N-2)
        do i = 1,M
          DPsi(j,i)= 2/dx*minmod(theta*(psir(j,i)-Psi(j,i))/(1+lmb*(aa(j)-aa(j-1))), &
          (Psi(j+1,i)-Psi(j,i))/(2+lmb*(2*aa(j)-aa(j-1)-aa(j+1))), &
          theta*(Psi(j+1,i)-psir(j,i))/(1+lmb*(aa(j)-aa(j+1))))
        end do
      end do

      ! Numerical Fluxes
      ALLOCATE(hh(N-1,M), pp(N-1,M), Du_ap(N-1,M))
      hh = 0.0_dp; pp = 0.0_dp; Du_ap = 0.0_dp
      do j = 1,(N-1)
        hh(j,:) = 0.5*(Fphir(j,:)+Fphil(j,:))-0.5*(uo(j+1,:)-uo(j,:))*aa(j)+&
        aa(j)*(1-lmb*aa(j))/4*(Du(j+1,:)+Du(j,:)) + lmb*dx/2*(aa(j))**2*DPsi(j,:)
      end do
      Du_ap = (uo(2:N,:)-uo(1:N-1,:))/dx
      !Du_ap = Du/dx
      ! Diffusion
      do j = 1,N-1
        pp(j,:) = 0.5*MATMUL(prob%K(uo(j+1,:))+prob%K(uo(j,:)),Du_ap(j,:))
      end do
      ALLOCATE(hhleft(M),ppleft(M),hhright(M),ppright(M))
      hhleft = 0.0_dp; hhright = 0.0_dp; ppleft = 0.0_dp; ppright = 0.0_dp
      if (boundary == PERIODIC) then
        hhleft = hh(1,:); ppleft = pp(1,:)
        hhright = hh(N-1,:); ppright = pp(N-1,:)
      end if
      j = 1 + ss
      rhs(j-ss,:) = - 1/dx * (hh(j,:) -hhleft - (pp(j,:)-ppleft))
      do j = (2+ss),(N-1-ss)
        rhs(j-ss,:) = - 1/dx * (hh(j,:)-hh(j-1,:)-(pp(j,:)-pp(j-1,:)))
      end do
      j = N-ss
      rhs(j-ss,:) =  -1/dx*(hhright-hh(j-1,:)-(ppright - pp(j-1,:)))

      DEALLOCATE(uo, uminus, uplus, aa, u_l, u_r, Df_l, Df_r, phi_l, phi_r, fphir, fphil, DPsi, Du)
      DEALLOCATE(hh, pp, Du_ap, Ful, Fulm, Fulp, Fur, Furm, Furp, hhleft, hhright, ppleft, ppright)
  end subroutine update_KT
END MODULE KT_scheme