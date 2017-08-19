MODULE eno_weno
USE decimal
IMPLICIT NONE
contains
function unif_crj(k) result(coefs)
  INTEGER k,i,j,l,m,r,q
  REAL(kind = dp), ALLOCATABLE  :: coefs(:,:)
  REAL(kind = dp)               :: psum, numer, numer_prod, denom

  ALLOCATE(coefs(k+1,k))
  coefs = 0.0_dp

  if (k == 1) then
    coefs = reshape([1.0_dp, 1.0_dp], shape(coefs))
  end if

  do i = 1,(k+1)
    r = i-2
    do j = 1,k
      psum = 0.0_dp
      do m = j,k
        numer = 0.0_dp
        do l = 0,k
          if (l/=m) then
            numer_prod = 1.0_dp
            do q = 0,k
              if (q /= m .and. q /= l) then
                numer_prod = numer_prod*(r-q+1)
              end if
            end do
            numer = numer + numer_prod
          end if
        end do
        denom = 1.0_dp
        do l = 0,k
          if (l/=m) then
            denom = denom * (m-l)
          end if
        end do
        psum = psum + numer/denom
      end do
      coefs(i,j) = psum
    end do
  end do
end function unif_crj

!Eno reconstruction for uniform mesh
function ENO_urec(dx, vloc, k, crj) result(vv)
  INTEGER           :: k, i, r, m, N, j
  REAL(kind = dp)   :: vloc(:), crj(:,:)
  REAL(kind = dp)   :: dx
  REAL(kind = dp)   :: vv(2), vl, vr
  REAL(kind = dp), ALLOCATABLE  :: vdiffs(:,:)

  N = size(vloc,1)
  ALLOCATE(vdiffs(k,N))
  vl = 0.0_dp ; vr = 0.0_dp; vv = 0.0_dp
  if (N /= 2*k-1) then
    print *, "dimension of vloc is not consistent with order ",k," ENO"
    stop
  end if
  ! Calculation of divided differences
  vdiffs(1,:) = vloc
  do i = 2,k
    vdiffs(i,1:N-i+1) = (vdiffs(i-1,2:N-i+2)-vdiffs(i-1,1:N-i+1))/dx 
  end do
  ! Calculation of Stencil
  r = 0
  if (k < 2) then 
    vv = [vloc(1),vloc(1)]
  else
    do m = 2,k
      if (abs(vdiffs(m,k-r-1)) < abs(vdiffs(m,k-r))) then
        r = r + 1
      end if
    end do
    do j = 0,(k-1)
      vl = vl + vloc(k-r+j)*crj(r+1,j+1)
      vr = vr + vloc(k-r+j)*crj(r+2,j+1)
    end do
    vv = [vl, vr]
  end if
  deallocate(vdiffs)
end function ENO_urec

! Eno reconstruction for unestructured mesh
function ENO_rec(xloc, vloc, k, crj) result(vv)
  INTEGER           :: k, i, r, m, N, j
  REAL(kind = dp)   :: xloc(:), vloc(:), crj(:,:)
  REAL(kind = dp)   :: dx, vl, vr
  REAL(kind = dp)   :: vv(2)
  REAL(kind = dp), ALLOCATABLE  :: vdiffs(:,:)

  N = size(vloc,1)
  ALLOCATE(vdiffs(k,N))
  vl = 0.0_dp; vr = 0.0_dp
  if (N /= 2*k-1) then
    print *, "dimension of vloc is not consistent with order ",k," ENO"
    stop
  end if
  ! Calculation of divided differences
  vdiffs(1,:) = vloc
  do i = 2,k
    vdiffs(i,1:N-i+1) = (vdiffs(i-1,2:N-i+2)-vdiffs(i-1,1:N-i+1))/(xloc(i:N)-xloc(1:N-i+1))
  end do
  ! Calculation of Stencil
  r = 0
  if (k < 2) then 
    vv = [vloc(1),vloc(1)]
  else
    do m = 2,k
      if (abs(vdiffs(m,k-r-1)) < abs(vdiffs(m,k-r))) then
        r = r + 1
      end if
    end do
    do j = 0,(k-1)
      vl = vl + vloc(k-r+j)*crj(r+1,j+1)
      vr = vr + vloc(k-r+j)*crj(r+2,j+1)
    end do
    vv = [vl, vr]
  end if
  deallocate(vdiffs)
end function ENO_rec

! WEno reconstruction for uniform mesh
! Order available: 1, 3, 5
subroutine get_betak(vloc, k, betak, dr)
  REAL(kind = dp)   :: vloc(:)
  INTEGER           :: k
  REAL(kind = dp)   :: betak(:), dr(:)

  betak = 0.0_dp; dr = 0.0_dp;
  if (k==2) then
      dr = [2/3,1/3]
      betak = [(vloc(3)-vloc(2))**2, (vloc(2)-vloc(1))**2]
  else if (k==3) then
      dr = [3.0_dp/10.0_dp, 3.0_dp/5.0_dp, 1/10.0_dp]
      betak = [13/12.0_dp*(vloc(3)-2*vloc(4)+vloc(5))**2 + 1/4.0_dp*(3*vloc(3)-4*vloc(4)+vloc(5))**2, &
      13/12.0_dp*(vloc(2)-2*vloc(3)+vloc(4))**2 + 1/4.0_dp*(vloc(2)-vloc(4))**2, &
      13/12.0_dp*(vloc(1)-2*vloc(2)+vloc(3))**2 + 1/4.0_dp*(3*vloc(3)-4*vloc(2)+vloc(1))**2]
  else
    print *, "WENO reconstruction of order $order is not implemented yet!"
    stop
  end if
end subroutine get_betak

function WENO_urec(vloc, order, crj) result (vv)
  REAL(kind = dp)   :: vloc(:), crj(:,:)
  INTEGER           :: order, r, i, k, N
  REAL(kind = dp)   :: vv(2), epsilon, vl, vr
  REAL(kind = dp), allocatable   :: alphal(:), alphar(:), omegal(:), omegar(:)
  REAL(kind = dp), allocatable   :: betak(:), dr(:), ulr(:), urr(:)

  epsilon = 1e-12_dp
  vl = 0.0_dp; vr = 0.0_dp
  N = size(vloc,1)
  if (N /= order) then
    print *, "dimension of vloc is not consistent with order ",order," WENO"
  end if
  k = floor((order + 1)/2.0_dp)
  ! Special case k = 1
  if (k == 1) then
    vl = vloc(1); vr = vloc(1)
    vv = [vl, vr]
  else
    ! Apply WENO procedure
    ALLOCATE(alphal(k), alphar(k), omegal(k), omegar(k))
    ALLOCATE(betak(k), dr(k))
    alphal = 0.0_dp; alphar = 0.0_dp; omegal = 0.0_dp; omegar = 0.0_dp
    betak = 0.0_dp; dr = 0.0_dp

    ! Compute k values of xl and xr based on different stencils
    ALLOCATE(ulr(k), urr(k))
    ulr = 0.0_dp; urr = 0.0_dp
    do r = 0,(k-1)
        do i=0,(k-1)
            urr(r+1) = urr(r+1) + crj(r+2,i+1)*vloc(k-r+i);
            ulr(r+1) = ulr(r+1) + crj(r+1,i+1)*vloc(k-r+i);
        end do
    end do

    ! Set up WENO coefficients do dif,erent orders - 2k-1
    CALL get_betak(vloc, k, betak, dr)

    ! Compute alpha parameters
    do r = 1,k
        alphar(r) = dr(r)/(epsilon+betak(r))**2;
        alphal(r) = dr(k+1-r)/(epsilon+betak(r))**2;
    end do

    ! Compute wENO weights parameters
    do r = 1,k
        omegal(r) = alphal(r)/sum(alphal);
        omegar(r) = alphar(r)/sum(alphar);
    end do

    ! Compute cell interface values
    do r = 1,k
        vl = vl + omegal(r)*ulr(r);
        vr = vr + omegar(r)*urr(r);
    end do
    vv = [vl, vr]
    DEALLOCATE(alphal, alphar, omegal, omegar, betak, dr, ulr, urr)
  end if
end function WENO_urec

function WENO_pm_rec(vmloc, vploc, order, crj) result(vv)
  real(kind = dp) :: vmloc(:), vploc(:), epsilon, crj(:,:), vv(2)
  INTEGER         :: order, N, k, r, i
  REAL(kind = dp), allocatable   :: alphal(:), alphar(:), omegal(:), omegar(:)
  REAL(kind = dp), allocatable   :: betark(:), betalk(:), dr(:), ulr(:), urr(:)
  REAL(kind = dp) :: vl, vr

  epsilon = 1e-12_dp
  vl = 0.0_dp; vr = 0.0_dp
  N = size(vmloc,1)
  if (N /= order) then
    print *, "dimension of vloc is not consistent with order ",order," WENO"
  end if
  k = floor((order + 1)/2.0_dp)
  ! Special case k = 1
  if (k == 1) then
    vl = vmloc(1); vr = vploc(1)
    vv = [vl, vr]
  else
    ! Apply WENO procedure
    ALLOCATE(alphal(k), alphar(k), omegal(k), omegar(k))
    ALLOCATE(betark(k), betalk(k), dr(k))
    alphal = 0.0_dp; alphar = 0.0_dp; omegal = 0.0_dp; omegar = 0.0_dp
    betark = 0.0_dp; betalk = 0.0_dp; dr = 0.0_dp

    ! Compute k values of xl and xr based on different stencils
    ALLOCATE(ulr(k), urr(k))
    ulr = 0.0_dp; urr = 0.0_dp

    do r=0,(k-1)
        do i=0,k-1
            urr(r+1) = urr(r+1) + crj(r+2,i+1)*vploc(k-r+i);
            ulr(r+1) = ulr(r+1) + crj(r+1,i+1)*vmloc(k-r+i);
        end do
    end do

    ! Set up WENO coefficients do dif,erent orders - 2k-1
    CALL get_betak(vploc, k, betark, dr)
    CALL get_betak(vmloc, k, betalk, dr)

    ! Compute alpha parameters
    do r=1,k
        alphar(r) = dr(r)/(epsilon+betark(r))**2
        alphal(r) = dr(k+1-r)/(epsilon+betalk(r))**2
    end do

    ! Compute wENO weights parameters
    do r=1,k
        omegal(r) = alphal(r)/sum(alphal)
        omegar(r) = alphar(r)/sum(alphar)
    end do

    ! Compute cell interface values
    do r=1,k
        vl = vl + omegal(r)*ulr(r);
        vr = vr + omegar(r)*urr(r);
    end do
    vv = [vl, vr]
    DEALLOCATE(alphal, alphar, omegal, omegar, betalk, betark, dr, ulr, urr)
  end if
end function WENO_pm_rec


! Mapped WEno reconstruction do uni,om me,h
! Reference:
!   A. Henrick, T. Aslam, J. Powers, Mapped weighted essentially non-oscillatory
!   schemes: Achiving optimal order  donear critical points

! Mapping function
function gk(omega, dr) result(g)
  real(kind = dp) :: omega(:), dr(:)
  real(kind = dp) :: g(size(omega,1))
  INTEGER         :: i
  g = 0.0_dp
  do i = 1,size(omega,1)
    g(i) = omega(i)*(dr(i)+dr(i)**2-3*dr(i)*omega(i)+omega(i)**2)/(dr(i)**2+omega(i)*(1-2*dr(i)))
  end do
end function gk

function MWENO_urec(vloc, order, crj) result(vv)
  real(kind = dp) :: vloc(:), crj(:,:), vv(2)
  INTEGER         :: order, N, k, r, i
  real(kind = dp) :: epsilon, vl, vr
  real(kind = dp), allocatable :: alphal(:), alphar(:), omegal(:), omegar(:), alphaml(:), alphamr(:), omegaml(:), omegamr(:)
  real(kind = dp), allocatable :: betak(:), dr(:), drl(:), ulr(:), urr(:)

  epsilon = 1e-12_dp
  vl = 0.0_dp; vr = 0.0_dp; vv = 0.0_dp
  N = size(vloc,1)
  if (N /= order) then
    print *,"dimension of vloc is not consistent with order $order WENO"
  end if
  k = floor((order + 1)/2.0_dp)
  ! Special case k = 1
  if (k == 1) then
    vl = vloc(1); vr = vloc(1)
    vv = [vl, vr]
  else
    ! Apply WENO procedure
    ALLOCATE(alphal(k), alphar(k), omegal(k), omegar(k), alphaml(k), alphamr(k), omegaml(k), omegamr(k))
    ALLOCATE(betak(k), dr(k), ulr(k), urr(k))
    alphal = 0.0_dp; alphar = 0.0_dp; omegal = 0.0_dp; omegar = 0.0_dp
    alphaml = 0.0_dp; alphamr = 0.0_dp; omegaml = 0.0_dp; omegamr = 0.0_dp
    betak = 0.0_dp; dr = 0.0_dp; drl = 0.0_dp; ulr = 0.0_dp; urr = 0.0_dp

    ! Compute k values of xl and xr based on different stencils
    do r=0,(k-1)
        do i=0,k-1
            urr(r+1) = urr(r+1) + crj(r+2,i+1)*vloc(k-r+i)
            ulr(r+1) = ulr(r+1) + crj(r+1,i+1)*vloc(k-r+i)
        end do
    end do

    ! Set up WENO coefficients for different orders 2k-1
    if (k==2) then
        dr = [2.0_dp/3.0_dp, 1.0_dp/3.0_dp]
        betak = [(vloc(3)-vloc(2))**2, (vloc(2)-vloc(1))**2]
    else if (k==3) then
        dr = [3.0_dp/10.0_dp, 3.0_dp/5.0_dp, 1.0_dp/10.0_dp]
        betak = [13.0_dp/12.0_dp*(vloc(3)-2.0_dp*vloc(4)+vloc(5))**2 + 1/4.0_dp*(3*vloc(3)-4*vloc(4)+vloc(5))**2, &
        13.0_dp/12.0_dp*(vloc(2)-2*vloc(3)+vloc(4))**2 + 1/4.0_dp*(vloc(2)-vloc(4))**2, &
        13.0_dp/12.0_dp*(vloc(1)-2*vloc(2)+vloc(3))**2 + 1/4.0_dp*(3.0_dp*vloc(3)-4*vloc(2)+vloc(1))**2]
    else
      print *, "WENO reconstruction of order ", order," is not implemented yet!"
      stop
    end if

    ! Compute alpha parameters
    do r=1,k
        alphar(r) = dr(r)/(epsilon+betak(r))**2;
        alphal(r) = dr(k+1-r)/(epsilon+betak(r))**2;
    end do

    ! Compute wENO weights parameters
    do r=1,k
        omegal(r) = alphal(r)/sum(alphal);
        omegar(r) = alphar(r)/sum(alphar);
    end do

    ! Compute alpha mapped parameters
    drl = dr(k : 1 : -1)
    alphamr = gk(omegar, dr)
    alphaml = gk(omegal, drl)

    ! Compute mapped wENO weights parameters
    do r=1,k
        omegaml(r) = alphaml(r)/sum(alphaml);
        omegamr(r) = alphamr(r)/sum(alphamr);
    end do

    ! Compute cell interface values
    do r=1,k
        vl = vl + omegaml(r)*ulr(r);
        vr = vr + omegamr(r)*urr(r);
    end do
    vv = [vl, vr]
    DEALLOCATE(alphal, alphar, omegal, omegar, alphaml, alphamr, omegaml, omegamr)
    DEALLOCATE(betak, dr, drl, ulr, urr)
  end if
end function MWENO_urec

function MWENO_pm_rec(vmloc,vploc,order, crj) result(vv)
  real(kind = dp), intent(in) :: vmloc(:), vploc(:), crj(:,:)
  INTEGER         :: order, N, k, i, r
  real(kind = dp) :: epsilon, vl, vr, vv(2)
  real(kind = dp), allocatable :: alphal(:), alphar(:), omegal(:), omegar(:), alphaml(:), alphamr(:), omegaml(:), omegamr(:)
  real(kind = dp), allocatable :: betark(:), betalk(:), dr(:), drl(:), ulr(:), urr(:)
  epsilon = 1e-12
  vl = 0.0_dp; vr = 0.0_dp; vv = 0.0_dp
  N = size(vmloc,1)
  if (N /= order) then
    print *, "dimension of vloc is not consistent with order ",order," WENO"
  end if
  k = floor((order + 1)/2.0_dp)
  ! Special case k = 1
  if (k == 1) then
    vl = vmloc(1); vr = vploc(1)
    vv = [vl, vr]
  else
    ! Apply WENO procedure
    ALLOCATE(alphal(k), alphar(k), omegal(k), omegar(k), alphaml(k), alphamr(k), omegaml(k), omegamr(k))
    ALLOCATE(betark(k), betalk(k), dr(k), drl(k), ulr(k), urr(k))
      alphal = 0.0_dp; alphar = 0.0_dp; omegal = 0.0_dp; omegar = 0.0_dp
    alphaml = 0.0_dp; alphamr = 0.0_dp; omegaml = 0.0_dp; omegamr = 0.0_dp
    betark = 0.0_dp; betalk = 0.0_dp; dr = 0.0_dp; ulr = 0.0_dp; urr = 0.0_dp

    ! Compute k values of xl and xr based on different stencils
    do r=0,(k-1)
        do i=0,k-1
            urr(r+1) = urr(r+1) + crj(r+2,i+1)*vploc(k-r+i)
            ulr(r+1) = ulr(r+1) + crj(r+1,i+1)*vmloc(k-r+i)
        end do
    end do

    ! Set up WENO coefficients do dif,erent orders - 2k-1
    CALL get_betak(vploc, k, betark, dr)
    CALL get_betak(vmloc, k, betalk, dr)

    ! Compute alpha parameters
    do r=1,k
        alphar(r) = dr(r)/(epsilon+betark(r))**2;
        alphal(r) = dr(k+1-r)/(epsilon+betalk(r))**2;
    end do

    ! Compute wENO weights parameters
    do r=1,k
        omegal(r) = alphal(r)/sum(alphal);
        omegar(r) = alphar(r)/sum(alphar);
    end do

    ! Compute alpha mapped parameters
    drl = dr(k : 1 : -1)
    alphamr = gk(omegar,dr)
    alphaml = gk(omegal,drl)

    ! Compute mapped wENO weights parameters
    do r=1,k
        omegaml(r) = alphaml(r)/sum(alphaml);
        omegamr(r) = alphamr(r)/sum(alphamr);
    end do

    ! Compute cell interface values
    do r=1,k
        vl = vl + omegaml(r)*ulr(r);
        vr = vr + omegamr(r)*urr(r);
    end do
    vv = [vl,vr]
    DEALLOCATE(alphal, alphar, omegal, omegar, alphaml, alphamr, omegaml, omegamr)
    DEALLOCATE(betark, betalk, dr, drl, ulr, urr)
  end if
end function MWENO_pm_rec


END MODULE eno_weno