MODULE limiters
USE decimal
USE FVTypes
IMPLICIT NONE
CONTAINS

function minmod(a,b,c) result(mm)
  REAL(kind=dp), intent(in) :: a,b,c
  REAL(kind=dp) :: mm
  if (a > 0.0_dp .and. b > 0.0_dp .and. c > 0.0_dp) then
    mm = minval([a,b,c])
  elseif (a < 0.0_dp .and. b < 0.0_dp .and. c < 0.0_dp) then 
    mm = maxval([a,b,c])
  else
    mm = 0.0_dp
  end if
end function minmod
function minmod2(a,b) result(mm)
  REAL(kind=dp), intent(in) :: a,b
  REAL(kind=dp) :: mm
  if (a > 0.0_dp .and. b > 0.0_dp) then
    mm = min(a,b)
  elseif (a < 0.0_dp .and. b < 0.0_dp) then 
    mm = max(a,b)
  else
    mm = 0.0_dp
  end if
end function minmod2
END MODULE limiters