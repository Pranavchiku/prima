module xinbd_mod

implicit none

private
public :: xinbd


contains


function xinbd(xbase, step, xl, xu, sl, su) result(x)
!--------------------------------------------------------------------------------------------------!
! This function sets X to XBASE + STEP, paying careful attention to the following bounds.
! 1. XBASE is a point between XL and XU (guaranteed);
! 2. STEP is a step between SL and SU (may be with rounding errors);
! 3. SL = XL - XBASE, SU = XU - XBASE;
! 4. X should be between XL and XU.
!--------------------------------------------------------------------------------------------------!
! Common modules
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, EPS, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : trueloc, linspace
use, non_intrinsic :: infnan_mod, only : is_finite

implicit none

! Inputs
real(RP), intent(in) :: xbase(:)
real(RP), intent(in) :: step(:)
real(RP), intent(in) :: xl(:)
real(RP), intent(in) :: xu(:)
real(RP), intent(in) :: sl(:)
real(RP), intent(in) :: su(:)
integer(IK) :: i
logical, allocatable :: mask(:)
integer(IK), allocatable :: tmp(:)
integer(IK), allocatable :: tmp2(:)

! Outputs
real(RP) :: x(size(xbase))

! Local variables
character(len=*), parameter :: srname = 'XINBD'
integer(IK) :: n
real(RP) :: s(size(xbase))

! Sizes
n = int(size(xbase), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(all(is_finite(xbase)), 'SIZE(XBASE) == N, XBASE is finite', srname)
    call assert(size(xl) == n .and. size(xu) == n, 'SIZE(XL) == N == SIZE(XU)', srname)
    call assert(all(xbase >= xl .and. xbase <= xu), 'XL <= XBASE <= XU', srname)
    call assert(size(sl) == n .and. size(su) == n, 'SIZE(SL) == N == SIZE(SU)', srname)
    call assert(all(step + 1.0E2_RP * EPS * max(ONE, abs(step)) >= sl .and. &
        & step - 1.0E2_RP * EPS * max(ONE, abs(step)) <= su), 'SL <= STEP <= SU', srname)
end if

!====================!
! Calculation starts !
!====================!

s = max(sl, min(su, step))
x = max(xl, min(xu, xbase + s))
allocate(mask(size(s)))
! mask = trueloc(s <= sl)
mask = s<=sl
allocate(tmp(count(mask)))
tmp = trueloc(mask)
! x(trueloc(s <= sl)) = xl(trueloc(s <= sl))
! x(tmp) = xl(tmp)
do i = 1, count(mask)
    x(tmp(i)) = xl(tmp(i))
end do

mask = s >= su
allocate(tmp2(count(mask)))
tmp2 = trueloc(mask)
! x(trueloc(s >= su)) = xu(trueloc(s >= su))
! x(tmp) = xu(tmp)
do i = 1, count(mask)
    x(tmp2(i)) = xu(tmp2(i))
end do

! deallocate(tmp)
! allocate(mask(count(s>=su)))
! ! x(trueloc(s >= su)) = xu(trueloc(s >= su))
! mask = trueloc(s >= su)
! x(mask) = xu(mask)
! deallocate(mask)

!====================!
!  Calculation ends  !
!====================!

if (DEBUGGING) then
    call assert(size(x) == n .and. all(x >= xl .and. x <= xu), 'SIZE(X) == N, XL <= X <= XU', srname)
    call assert(all(x <= xl .or. step > sl), 'X == XL if STEP <= SL', srname)
    call assert(all(x >= xu .or. step < su), 'X == XU if STEP >= SU', srname)
end if

end function xinbd


end module xinbd_mod
