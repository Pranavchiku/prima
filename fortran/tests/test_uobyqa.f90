module recursive_mod
implicit none
private
public :: recursive_fun1

contains

subroutine chrosen(x, f)
use, non_intrinsic :: consts_mod, only : RP
implicit none
real(RP), intent(in) :: x(:)
real(RP), intent(out) :: f
integer :: n
real(RP), parameter :: alpha = 4.0_RP
n = size(x)
f = sum((x(1:n - 1) - 1.0_RP)**2 + alpha * (x(2:n) - x(1:n - 1)**2)**2); 
end subroutine chrosen

subroutine recursive_fun1(x, f)
use, non_intrinsic :: consts_mod, only : RP
use, non_intrinsic :: uobyqa_mod, only : uobyqa
implicit none
real(RP), intent(in) :: x(:)
real(RP), intent(out) :: f
real(RP) :: x_loc(size(x))
x_loc = x
call uobyqa(chrosen, x_loc, f)
end subroutine recursive_fun1

end module recursive_mod


module test_solver_mod
!--------------------------------------------------------------------------------------------------!
! This module tests UOBYQA on a few simple problems.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: September 2021
!
! Last Modified: Thursday, March 07, 2024 PM02:41:31
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: test_solver


contains


subroutine test_solver(probs, mindim, maxdim, dimstride, nrand, randseed, testdim)

use, non_intrinsic :: consts_mod, only : RP, IK, TWO, TEN, ZERO, REALMAX
use, non_intrinsic :: memory_mod, only : safealloc
use, non_intrinsic :: noise_mod, only : noisy, noisy_calfun, orig_calfun
use, non_intrinsic :: param_mod, only : MINDIM_DFT, MAXDIM_DFT, DIMSTRIDE_DFT, NRAND_DFT, RANDSEED_DFT
use, non_intrinsic :: prob_mod, only : PNLEN, PROB_T, construct, destruct
use, non_intrinsic :: rand_mod, only : setseed, rand, randn
use, non_intrinsic :: string_mod, only : strip, istr
use, non_intrinsic :: uobyqa_mod, only : uobyqa

implicit none

character(len=PNLEN), intent(in), optional :: probs(:)
integer(IK), intent(in), optional :: mindim
integer(IK), intent(in), optional :: maxdim
integer(IK), intent(in), optional :: dimstride
integer(IK), intent(in), optional :: nrand
integer, intent(in), optional :: randseed
character(len=*), intent(in), optional :: testdim

character(len=*), parameter :: bigprob = 'bigprob'
character(len=*), parameter :: solname = 'uobyqa'
character(len=:), allocatable :: testdim_loc
character(len=PNLEN) :: probname
character(len=PNLEN) :: probs_loc(100)
integer :: randseed_loc
integer :: rseed
integer(IK) :: dimstride_loc
integer(IK) :: iprint
integer(IK) :: iprob
integer(IK) :: irand
integer(IK) :: maxdim_loc
integer(IK) :: maxfun
integer(IK) :: maxhist
integer(IK) :: mindim_loc
integer(IK) :: n
integer(IK) :: nprobs
integer(IK) :: npt
integer(IK) :: nrand_loc
integer(IK), parameter :: bign = 80_IK
integer(IK), parameter :: largen = 160_IK
real(RP) :: f
real(RP) :: ftarget
real(RP) :: rhobeg
real(RP) :: rhoend
real(RP), allocatable :: fhist(:)
real(RP), allocatable :: x(:)
real(RP), allocatable :: xhist(:, :)
type(PROB_T) :: prob

if (present(probs)) then
    nprobs = int(size(probs), kind(nprobs))
    probs_loc(1:nprobs) = probs
else
    nprobs = 5
    probs_loc(1:nprobs) = ['chebyquad', 'chrosen  ', 'trigsabs ', 'trigssqs ', 'vardim   ']
end if

if (present(mindim)) then
    mindim_loc = mindim
else
    mindim_loc = MINDIM_DFT
end if

if (present(maxdim)) then
    maxdim_loc = maxdim
else
    maxdim_loc = MAXDIM_DFT
end if

if (present(dimstride)) then
    dimstride_loc = dimstride
else
    dimstride_loc = DIMSTRIDE_DFT
end if

if (present(nrand)) then
    nrand_loc = nrand
else
    nrand_loc = NRAND_DFT
end if

if (present(randseed)) then
    randseed_loc = randseed
else
    randseed_loc = RANDSEED_DFT
end if

if (present(testdim)) then
    testdim_loc = testdim
else
    testdim_loc = 'small'
end if


! Test the big problem
if (testdim_loc == 'big' .or. testdim_loc == 'large') then
    probname = bigprob
    n = merge(bign, largen, testdim_loc == 'big')
    call construct(prob, probname, n)
    do irand = 1, 1  ! The test is expensive
        rseed = int(sum(istr(solname)) + sum(istr(probname)) + n + irand + RP + randseed_loc)
        call setseed(rseed)
        iprint = 2_IK
        npt = (n + 2_IK) * (n + 1_IK) / 2_IK
        if (int(npt) + 2000 > huge(0_IK)) then
            maxfun = huge(0_IK)
        else
            maxfun = npt + int(2000.0_RP * rand(), IK)
        end if
        maxhist = maxfun
        ftarget = -REALMAX
        rhobeg = noisy(prob % Delta0)
        rhoend = max(1.0E-6_RP, rhobeg * 10.0_RP**(6.0_RP * rand() - 6.0_RP))
        call safealloc(x, n) ! Not all compilers support automatic allocation yet, e.g., Absoft.
        x = noisy(prob % x0)
        orig_calfun => prob % calfun

        print '(/A, I0, A, I0, A, I0)', strip(probname)//': N = ', n, ', MAXFUN = ', maxfun, ', Random test ', irand
        call uobyqa(noisy_calfun, x, f, rhobeg=rhobeg, rhoend=rhoend, maxfun=maxfun, &
            & maxhist=maxhist, fhist=fhist, xhist=xhist, ftarget=ftarget, iprint=iprint)

        deallocate (x)
        ! nullify (orig_calfun)
    end do
    ! DESTRUCT deallocates allocated arrays/pointers and nullify the pointers. Must be called.
    call destruct(prob)  ! Destruct the testing problem.
else
    do iprob = 1, nprobs
        probname = probs_loc(iprob)
        do n = mindim_loc, maxdim_loc, dimstride_loc
            call construct(prob, probname, n)  ! Construct the testing problem.

            do irand = 1, nrand_loc
                ! Initialize the random seed using N, IRAND, RP, and RANDSEED_LOC. Do not include IK so
                ! that the results for different IK are the same.
                rseed = int(sum(istr(solname)) + sum(istr(probname)) + n + irand + RP + randseed_loc)
                call setseed(rseed)
                iprint = int(sign(min(3.0_RP, 1.5_RP * abs(randn())), randn()), kind(iprint))
                maxfun = int(2.0E2_RP * rand() * real(n, RP), kind(maxfun))
                if (rand() <= 0.2) then
                    maxfun = 0
                end if
                maxhist = int(TWO * rand() * real(max(10_IK * n, maxfun), RP), kind(maxhist))
                if (rand() <= 0.2) then
                    maxhist = 0
                end if
                if (rand() <= 0.2) then
                    ftarget = -TEN**abs(TWO * randn())
                elseif (rand() <= 0.2) then  ! Note that the value of rand() changes.
                    ftarget = REALMAX
                else
                    ftarget = -REALMAX
                end if

                rhobeg = noisy(prob % Delta0)
                rhoend = max(1.0E-6_RP, rhobeg * 10.0_RP**(6.0_RP * rand() - 5.0_RP))
                if (rand() <= 0.2) then
                    rhoend = rhobeg
                elseif (rand() <= 0.2) then  ! Note that the value of rand() changes.
                    rhobeg = ZERO
                end if
                call safealloc(x, n) ! Not all compilers support automatic allocation yet, e.g., Absoft.
                x = noisy(prob % x0)
                orig_calfun => prob % calfun

                print '(/A, I0, A, I0)', strip(probname)//': N = ', n, ', Random test ', irand
                call uobyqa(noisy_calfun, x, f, rhobeg=rhobeg, rhoend=rhoend, maxfun=maxfun, &
                    & maxhist=maxhist, fhist=fhist, xhist=xhist, ftarget=ftarget, iprint=iprint)

                deallocate (x)
                ! nullify (orig_calfun)
            end do

            ! DESTRUCT deallocates allocated arrays/pointers and nullify the pointers. Must be called.
            call destruct(prob)  ! Destruct the testing problem.
        end do
    end do
end if


! Test recursive call.
! The depth of the recursion is 2. The first recursion is in RECURSIVE_FUN1, and the second is in
! RECURSIVE_FUN2. RECURSIVE_FUN1(Y) is defined by minimizing the CHROSEN function with Y being
! the starting point. RECURSIVE_FUN2(Y) is defined by RECURSIVE_FUN1 in a similar way. Note
! that RECURSIVE_FUN1 is essentially a constant function.
n = 3_IK
print '(/A, I0)', 'Testing recursive call of '//solname//' on a problem with N = ', n
call safealloc(x, n)
x = randn(n)
call uobyqa(recursive_fun2, x, f, iprint=2_IK)
deallocate (x)

contains

subroutine recursive_fun2(x_internal, f_internal)
use, non_intrinsic :: recursive_mod, only : recursive_fun1
implicit none
real(RP), intent(in) :: x_internal(:)
real(RP), intent(out) :: f_internal
real(RP) :: x_loc(size(x_internal))
x_loc = x_internal
call uobyqa(recursive_fun1, x_loc, f_internal)
end subroutine recursive_fun2

end subroutine test_solver


end module test_solver_mod
