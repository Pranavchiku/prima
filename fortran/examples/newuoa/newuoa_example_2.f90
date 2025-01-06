!--------------------------------------------------------------------------------------------------!
! This is an example to illustrate the usage of the solver.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code.
!
! Started: July 2020
!
! Last Modified: Friday, March 15, 2024 PM02:06:18
!--------------------------------------------------------------------------------------------------!


!-------------------------------- THE MODULE THAT IMPLEMENTS CALFUN -------------------------------!
module calfun_mod

implicit none
private
public :: RP, IK, calfun
integer, parameter :: RP = kind(0.0D0)
integer, parameter :: IK = kind(0)
! N.B.: We assume that PRIMA_REAL_PRECISION = 64 (double precision) and PRIMA_INTEGER_KIND = 0
! (default kind). Revise RP and IK if this is not the case.

contains

subroutine calfun(x, f)
! The Chebyquad test problem (Fletcher, 1965)
implicit none

real(RP), intent(in) :: x(:)
real(RP), intent(out) :: f

integer :: i, n
real(RP) :: y(size(x) + 1, size(x) + 1), tmp

n = size(x)

y(1:n, 1) = 1.0_RP
y(1:n, 2) = 2.0_RP * x - 1.0_RP
do i = 2, n
    y(1:n, i + 1) = 2.0_RP * y(1:n, 2) * y(1:n, i) - y(1:n, i - 1)
end do

f = 0.0_RP
do i = 1, n + 1
    tmp = sum(y(1:n, i)) / real(n, RP)
    if (modulo(i, 2) /= 0) then
        tmp = tmp + 1.0_RP / real(i * i - 2 * i, RP)
    end if
    f = f + tmp * tmp
end do
end subroutine calfun

end module calfun_mod

module newuoa_mod
!--------------------------------------------------------------------------------------------------!
! NEWUOA_MOD is a module providing the reference implementation of Powell's NEWUOA algorithm in
!
! M. J. D. Powell, The NEWUOA software for unconstrained optimization without derivatives, In Large-
! Scale Nonlinear Optimization, eds. G. Di Pillo and M. Roma, 255--297, Springer, New York, 2006
!
! NEWUOA approximately solves
!
!   min F(X),
!
! where X is a vector of variables that has N components and F is a real-valued objective function.
! It tackles the problem by a trust region method that forms quadratic models by interpolation.
! There can be some freedom in the interpolation conditions, which is taken up by minimizing the
! Frobenius norm of the change to the second derivative of the quadratic model, beginning with a
! zero matrix.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on the NEWUOA paper and Powell's code, with
! modernization, bug fixes, and improvements.
!
! Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: July 2020
!
! Last Modified: Thursday, February 22, 2024 PM03:30:12
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: newuoa


contains


subroutine newuoa(calfun, x, &
    & f, nf, rhobeg, rhoend, ftarget, maxfun, npt, iprint, eta1, eta2, gamma1, gamma2, &
    & xhist, fhist, maxhist, callback_fcn, info)
! Common modules
use, non_intrinsic :: consts_mod, only : DEBUGGING
use, non_intrinsic :: consts_mod, only : MAXFUN_DIM_DFT
use, non_intrinsic :: consts_mod, only : RHOBEG_DFT, RHOEND_DFT, FTARGET_DFT, IPRINT_DFT
use, non_intrinsic :: consts_mod, only : RP, IK, TWO, HALF, TEN, TENTH, EPS
use, non_intrinsic :: debug_mod, only : assert, warning
use, non_intrinsic :: evaluate_mod, only : moderatex
use, non_intrinsic :: history_mod, only : prehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite, is_posinf
use, non_intrinsic :: memory_mod, only : safealloc
use, non_intrinsic :: pintrf_mod, only : OBJ, CALLBACK
use, non_intrinsic :: preproc_mod, only : preproc
use, non_intrinsic :: string_mod, only : num2str

! Solver-specific modules
use, non_intrinsic :: newuob_mod, only : newuob

implicit none

! Compulsory arguments
procedure(OBJ) :: calfun  ! N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER
real(RP), intent(inout) :: x(:)

! Optional inputs
procedure(CALLBACK), optional :: callback_fcn
integer(IK), intent(in), optional :: iprint
integer(IK), intent(in), optional :: maxfun
integer(IK), intent(in), optional :: maxhist
integer(IK), intent(in), optional :: npt
real(RP), intent(in), optional :: eta1
real(RP), intent(in), optional :: eta2
real(RP), intent(in), optional :: ftarget
real(RP), intent(in), optional :: gamma1
real(RP), intent(in), optional :: gamma2
real(RP), intent(in), optional :: rhobeg
real(RP), intent(in), optional :: rhoend

! Optional outputs
integer(IK), intent(out), optional :: info
integer(IK), intent(out), optional :: nf
real(RP), intent(out), optional :: f
real(RP), intent(out), optional, allocatable :: fhist(:)
real(RP), intent(out), optional, allocatable :: xhist(:, :)

! Local variables
character(len=*), parameter :: solver = 'NEWUOA'
character(len=*), parameter :: srname = 'NEWUOA'
integer(IK) :: info_loc
integer(IK) :: iprint_loc
integer(IK) :: maxfun_loc
integer(IK) :: maxhist_loc
integer(IK) :: n
integer(IK) :: nf_loc
integer(IK) :: nhist
integer(IK) :: npt_loc
real(RP) :: eta1_loc
real(RP) :: eta2_loc
real(RP) :: f_loc
real(RP) :: ftarget_loc
real(RP) :: gamma1_loc
real(RP) :: gamma2_loc
real(RP) :: rhobeg_loc
real(RP) :: rhoend_loc
real(RP), allocatable :: fhist_loc(:)
real(RP), allocatable :: xhist_loc(:, :)

! Sizes
n = int(size(x), kind(n))

! Replace any NaN in X by ZERO and Inf/-Inf in X by REALMAX/-REALMAX.
x = moderatex(x)

! Read the inputs.

! If RHOBEG is present, then RHOBEG_LOC is a copy of RHOBEG; otherwise, RHOBEG_LOC takes the default
! value for RHOBEG, taking the value of RHOEND into account. Note that RHOEND is considered only if
! it is present and it is VALID (i.e., finite and positive). The other inputs are read similarly.
if (present(rhobeg)) then
    rhobeg_loc = rhobeg
elseif (present(rhoend)) then
    ! Fortran does not take short-circuit evaluation of logic expressions. Thus it is WRONG to
    ! combine the evaluation of PRESENT(RHOEND) and the evaluation of IS_FINITE(RHOEND) as
    ! "IF (PRESENT(RHOEND) .AND. IS_FINITE(RHOEND))". The compiler may choose to evaluate the
    ! IS_FINITE(RHOEND) even if PRESENT(RHOEND) is false!
    if (is_finite(rhoend) .and. rhoend > 0) then
        rhobeg_loc = max(TEN * rhoend, RHOBEG_DFT)
    else
        rhobeg_loc = RHOBEG_DFT
    end if
else
    rhobeg_loc = RHOBEG_DFT
end if

if (present(rhoend)) then
    rhoend_loc = rhoend
elseif (rhobeg_loc > 0) then
    rhoend_loc = max(EPS, min((RHOEND_DFT / RHOBEG_DFT) * rhobeg_loc, RHOEND_DFT))
else
    rhoend_loc = RHOEND_DFT
end if

if (present(ftarget)) then
    ftarget_loc = ftarget
else
    ftarget_loc = FTARGET_DFT
end if

if (present(maxfun)) then
    maxfun_loc = maxfun
else
    maxfun_loc = MAXFUN_DIM_DFT * n
end if

if (present(npt)) then
    npt_loc = npt
elseif (maxfun_loc >= n + 3_IK) then  ! Take MAXFUN into account if it is valid.
    npt_loc = min(maxfun_loc - 1_IK, 2_IK * n + 1_IK)
else
    npt_loc = 2_IK * n + 1_IK
end if

if (present(iprint)) then
    iprint_loc = iprint
else
    iprint_loc = IPRINT_DFT
end if

if (present(eta1)) then
    eta1_loc = eta1
elseif (present(eta2)) then
    if (eta2 > 0 .and. eta2 < 1) then
        eta1_loc = max(EPS, eta2 / 7.0_RP)
    end if
else
    eta1_loc = TENTH
end if

if (present(eta2)) then
    eta2_loc = eta2
elseif (eta1_loc > 0 .and. eta1_loc < 1) then
    eta2_loc = (eta1_loc + TWO) / 3.0_RP
else
    eta2_loc = 0.7_RP
end if

if (present(gamma1)) then
    gamma1_loc = gamma1
else
    gamma1_loc = HALF
end if

if (present(gamma2)) then
    gamma2_loc = gamma2
else
    gamma2_loc = TWO
end if

if (present(maxhist)) then
    maxhist_loc = maxhist
else
    maxhist_loc = maxval([maxfun_loc, n + 3_IK, MAXFUN_DIM_DFT * n])
end if

! Preprocess the inputs in case some of them are invalid.
call preproc(solver, n, iprint_loc, maxfun_loc, maxhist_loc, ftarget_loc, rhobeg_loc, rhoend_loc, &
    & npt=npt_loc, eta1=eta1_loc, eta2=eta2_loc, gamma1=gamma1_loc, gamma2=gamma2_loc)

! Further revise MAXHIST_LOC according to MAXHISTMEM, and allocate memory for the history.
! In MATLAB/Python/Julia/R implementation, we should simply set MAXHIST = MAXFUN and initialize
! FHIST = NaN(1, MAXFUN), XHIST = NaN(N, MAXFUN) if they are requested; replace MAXFUN with 0 for
! the history that is not requested.
call prehist(maxhist_loc, n, present(xhist), xhist_loc, present(fhist), fhist_loc)


!-------------------- Call NEWUOB, which performs the real calculations. --------------------------!
if (present(callback_fcn)) then
    call newuob(calfun, iprint_loc, maxfun_loc, npt_loc, eta1_loc, eta2_loc, ftarget_loc, gamma1_loc, &
        & gamma2_loc, rhobeg_loc, rhoend_loc, x, nf_loc, f_loc, fhist_loc, xhist_loc, info_loc, callback_fcn)
else
    call newuob(calfun, iprint_loc, maxfun_loc, npt_loc, eta1_loc, eta2_loc, ftarget_loc, gamma1_loc, &
        & gamma2_loc, rhobeg_loc, rhoend_loc, x, nf_loc, f_loc, fhist_loc, xhist_loc, info_loc)
end if
!--------------------------------------------------------------------------------------------------!


! Write the outputs.

if (present(f)) then
    f = f_loc
end if

if (present(nf)) then
    nf = nf_loc
end if

if (present(info)) then
    info = info_loc
end if

! Copy XHIST_LOC to XHIST if needed.
if (present(xhist)) then
    nhist = min(nf_loc, int(size(xhist_loc, 2), IK))
    !----------------------------------------------------!
    call safealloc(xhist, n, nhist)  ! Removable in F2003.
    !----------------------------------------------------!
    xhist = xhist_loc(:, 1:nhist)
    ! N.B.:
    ! 0. Allocate XHIST as long as it is present, even if the size is 0; otherwise, it will be
    ! illegal to enquire XHIST after exit.
    ! 1. Even though Fortran 2003 supports automatic (re)allocation of allocatable arrays upon
    ! intrinsic assignment, we keep the line of SAFEALLOC, because some very new compilers (Absoft
    ! Fortran 21.0) are still not standard-compliant in this respect.
    ! 2. NF may not be present. Hence we should NOT use NF but NF_LOC.
    ! 3. When SIZE(XHIST_LOC, 2) > NF_LOC, which is the normal case in practice, XHIST_LOC contains
    ! GARBAGE in XHIST_LOC(:, NF_LOC + 1 : END). Therefore, we MUST cap XHIST at NF_LOC so that
    ! XHIST contains only valid history. For this reason, there is no way to avoid allocating
    ! two copies of memory for XHIST unless we declare it to be a POINTER instead of ALLOCATABLE.
end if
! F2003 automatically deallocate local ALLOCATABLE variables at exit, yet we prefer to deallocate
! them immediately when they finish their jobs.
deallocate (xhist_loc)

! Copy FHIST_LOC to FHIST if needed.
if (present(fhist)) then
    nhist = min(nf_loc, int(size(fhist_loc), IK))
    !--------------------------------------------------!
    call safealloc(fhist, nhist)  ! Removable in F2003.
    !--------------------------------------------------!
    fhist = fhist_loc(1:nhist)  ! The same as XHIST, we must cap FHIST at NF_LOC.
end if
deallocate (fhist_loc)

! If MAXFHIST_IN >= NF_LOC > MAXFHIST_LOC, warn that not all history is recorded.
if ((present(xhist) .or. present(fhist)) .and. maxhist_loc < nf_loc) then
    call warning(solver, 'Only the history of the last '//num2str(maxhist_loc)//' function evaluation(s) is recorded')
end if

! Postconditions
if (DEBUGGING) then
    call assert(nf_loc <= maxfun_loc, 'NF <= MAXFUN', srname)
    call assert(size(x) == n .and. .not. any(is_nan(x)), 'SIZE(X) == N, X does not contain NaN', srname)
    nhist = min(nf_loc, maxhist_loc)
    if (present(xhist)) then
        call assert(size(xhist, 1) == n .and. size(xhist, 2) == nhist, 'SIZE(XHIST) == [N, NHIST]', srname)
        call assert(.not. any(is_nan(xhist)), 'XHIST does not contain NaN', srname)
    end if
    if (present(fhist)) then
        call assert(size(fhist) == nhist, 'SIZE(FHIST) == NHIST', srname)
        call assert(.not. any(is_nan(fhist) .or. is_posinf(fhist)), 'FHIST does not contain NaN/+Inf', srname)
        call assert(.not. any(fhist < f_loc), 'F is the smallest in FHIST', srname)
    end if
end if

end subroutine newuoa


end module newuoa_mod

!---------------------------------------- THE MAIN PROGRAM ----------------------------------------!
program newuoa_exmp

! ! The following line makes the solver available.
 use newuoa_mod, only : newuoa

! ! The following line specifies which module provides CALFUN.
 use calfun_mod, only : RP, IK, calfun

 implicit none

 integer, parameter :: n = 6
 integer :: i, nf, info
 real(RP) :: f, x(n)

! ! The following lines illustrates how to call the solver to solve the Chebyquad problem.
 x = [(real(i, RP) / real(n + 1, RP), i=1, n)]  ! Define the starting point.
 call newuoa(calfun, x, f)  ! This call will not print anything.

! ! In addition to the compulsory arguments, the following illustration specifies also RHOBEG and
! ! IPRINT, which are optional. All the unspecified optional arguments (RHOEND, MAXFUN, etc.) will
! ! take their default values coded in the solver.
 x = [(real(i, RP) / real(n + 1, RP), i=1, n)]  ! Define the starting point.
 call newuoa(calfun, x, f, rhobeg=0.2_RP * x(1), iprint=1_IK, nf=nf, info=info)

end program newuoa_exmp
