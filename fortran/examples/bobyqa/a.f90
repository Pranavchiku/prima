!--------------------------------------------------------------------------------------------------!
! This is an example to illustrate the usage of the solver.
!
! The objective function is trivial. This is intentional, as the focus is how to use the API.
!--------------------------------------------------------------------------------------------------!


!------------------------- THE MODULE THAT IMPLEMENTS CALFUN, CALLBACK_FCN ------------------------!
module calfun_mod2

    implicit none
    private
    public :: IK, RP, calfun, callback_fcn
    integer, parameter :: RP = kind(0.0D0)
    integer, parameter :: IK = kind(0)
    ! N.B.: We assume that PRIMA_REAL_PRECISION = 64 (double precision) and PRIMA_INTEGER_KIND = 0
    ! (default kind). Revise RP and IK if this is not the case.
    
    contains
    
    ! Objective function
    subroutine calfun(x, f)
    implicit none
    
    ! Inputs
    real(RP), intent(in) :: x(:)
    
    ! Outputs
    real(RP), intent(out) :: f
    
    f = (x(1) - 5.0_RP)**2 + (x(2) - 4.0_RP)**2
    
    end subroutine calfun
    
    ! Callback function
    subroutine callback_fcn(x, f, nf, tr, cstrv, nlconstr, terminate)
    implicit none
    real(RP), intent(in) :: x(:)
    real(RP), intent(in) :: f
    integer(IK), intent(in) :: nf
    integer(IK), intent(in) :: tr
    real(RP), intent(in), optional :: cstrv
    real(RP), intent(in), optional :: nlconstr(:)
    logical, intent(out), optional :: terminate
    
    if (.false.) print *, cstrv      ! Suppress compiler warning about unused variable
    if (.false.) print *, nlconstr   ! Suppress compiler warning about unused variable
    
    write (*, '("Best point so far: x = [", F6.4, ", ", F6.4, "], f = ", F6.3, ", nf = ", I0, ", tr = ", I0, "")') &
        & x(1), x(2), f, nf, tr
    
    terminate = .false.
    
    end subroutine callback_fcn
    
    end module calfun_mod2
    
    ! TODO:
    ! 1. Improve RESCUE so that it accepts an [XNEW, FNEW] that is not interpolated yet, or even accepts
    ! [XRESERVE, FRESERVE], which contains points that have been evaluated.
    !
    module bobyqb_mod2
    !--------------------------------------------------------------------------------------------------!
    ! This module performs the major calculations of BOBYQA.
    !
    ! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the BOBYQA paper.
    !
    ! N.B. (Zaikun 20230312): In Powell's code, the strategy concerning RESCUE is a bit complex.
    !
    ! 1. Suppose that a trust-region step D is calculated. Powell's code sets KNEW_TR before evaluating
    ! F at the trial point XOPT+D, assuming that the value of F at this point is not better than the
    ! current FOPT. With this KNEW_TR, the denominator of the update is calculated. If this denominator
    ! is sufficiently large, then evaluate F at XOPT+D, recalculate KNEW_TR if the function value turns
    ! out better than FOPT, and perform the update to include XOPT+D in the interpolation. If the
    ! denominator is not sufficiently large, then RESCUE is called, and another trust-region step is
    ! taken immediately after, discarding the previously calculated trust-region step D.
    !
    ! 2. Suppose that a geometry step D is calculated. Then KNEW_GEO must have been set before. Powell's
    ! code then calculates the denominator of the update. If the denominator is sufficiently large, then
    ! evaluate F at XOPT+D, and perform the update. If the denominator is not sufficiently large, then
    ! RESCUE is called; if RESCUE does not evaluate F at any new point (allowed by Powell's code but not
    ! ours), then take a new geometry step, or else take a trust-region step, discarding the previously
    ! calculated geometry step D in both cases.
    !
    ! 3. If it turns out necessary to call RESCUE again, but no new function value has been evaluated
    ! after the last RESCUE, then Powell's code will terminate.
    !
    ! Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    !
    ! Started: February 2022
    !
    ! Last Modified: Sunday, March 31, 2024 PM07:31:32
    !--------------------------------------------------------------------------------------------------!
    
    implicit none
    private
    public :: bobyqb
    
    
    contains
    
    
    subroutine bobyqb(calfun, iprint, maxfun, npt, eta1, eta2, ftarget, gamma1, gamma2, rhobeg, rhoend, &
        & xl, xu, x, nf, f, fhist, xhist, info, callback_fcn)
    !--------------------------------------------------------------------------------------------------!
    ! This subroutine performs the major calculations of BOBYQA.
    !
    ! IPRINT, MAXFUN, MAXHIST, NPT, ETA1, ETA2, FTARGET, GAMMA1, GAMMA2, RHOBEG, RHOEND, XL, XU, X, NF,
    ! F, FHIST, XHIST, and INFO are identical to the corresponding arguments in subroutine BOBYQA.
    !
    ! XBASE holds a shift of origin that should reduce the contributions from rounding errors to values
    !   of the model and Lagrange functions.
    ! SL and SU hold XL - XBASE and XU - XBASE, respectively.
    ! XOPT is the displacement from XBASE of the best vector of variables so far (i.e., the one provides
    !   the least calculated F so far). XOPT satisfies SL(I) <= XOPT(I) <= SU(I), with appropriate
    !   equalities when XOPT is on a constraint boundary. FOPT = F(XOPT + XBASE). However, we do not
    !   save XOPT and FOPT explicitly, because XOPT = XPT(:, KOPT) and FOPT = FVAL(KOPT), which is
    !   explained below.
    ! [XPT, FVAL, KOPT] describes the interpolation set:
    ! XPT contains the interpolation points relative to XBASE, each COLUMN for a point; FVAL holds the
    !   values of F at the interpolation points; KOPT is the index of XOPT in XPT.
    ! [GOPT, HQ, PQ] describes the quadratic model: GOPT will hold the gradient of the quadratic model
    !   at XBASE + XOPT; HQ will hold the explicit second order derivatives of the quadratic model; PQ
    !   will contain the parameters of the implicit second order derivatives of the quadratic model.
    ! [BMAT, ZMAT] describes the matrix H in the BOBYQA paper (eq. 2.7), which is the inverse of
    !   the coefficient matrix of the KKT system for the least-Frobenius norm interpolation problem:
    ! ZMAT will hold a factorization of the leading NPT*NPT submatrix of H, the factorization being
    !   OMEGA = ZMAT*ZMAT^T, which provides both the correct rank and positive semi-definiteness. BMAT
    !   will hold the last N ROWs of H except for the (NPT+1)th column. Note that the (NPT + 1)th row
    !   and column of H are not saved as they are unnecessary for the calculation.
    ! D is reserved for trial steps from XOPT. It is chosen by subroutine TRSBOX or GEOSTEP. Usually
    !   XBASE + XOPT + D is the vector of variables for the next call of CALFUN.
    !--------------------------------------------------------------------------------------------------!
    
    ! Common modules
    use, non_intrinsic :: checkexit_mod, only : checkexit
    use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, TEN, TENTH, REALMAX, DEBUGGING
    use, non_intrinsic :: debug_mod, only : assert!, wassert, validate
    use, non_intrinsic :: evaluate_mod, only : evaluate
    use, non_intrinsic :: history_mod, only : savehist, rangehist
    use, non_intrinsic :: infnan_mod, only : is_nan, is_finite, is_posinf
    use, non_intrinsic :: infos_mod, only : INFO_DFT, SMALL_TR_RADIUS, MAXTR_REACHED, DAMAGING_ROUNDING,&
        & NAN_INF_MODEL, CALLBACK_TERMINATE
    use, non_intrinsic :: linalg_mod, only : norm
    use, non_intrinsic :: message_mod, only : retmsg, rhomsg, fmsg
    use, non_intrinsic :: pintrf_mod, only : OBJ, CALLBACK
    use, non_intrinsic :: powalg_mod, only : quadinc, calden, calvlag!, errquad
    use, non_intrinsic :: ratio_mod, only : redrat
    use, non_intrinsic :: redrho_mod, only : redrho
    use, non_intrinsic :: shiftbase_mod, only : shiftbase
    use, non_intrinsic :: xinbd_mod, only : xinbd
    
    ! Solver-specific modules
    use, non_intrinsic :: geometry_bobyqa_mod, only : geostep, setdrop_tr
    use, non_intrinsic :: initialize_bobyqa_mod, only : initxf, initq, inith
    use, non_intrinsic :: rescue_mod, only : rescue
    use, non_intrinsic :: trustregion_bobyqa_mod, only : trsbox, trrad
    use, non_intrinsic :: update_bobyqa_mod, only : updatexf, updateq, tryqalt, updateh
    
    implicit none
    
    ! Inputs
    procedure(OBJ) :: calfun  ! N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER
    procedure(CALLBACK), optional :: callback_fcn
    integer(IK), intent(in) :: iprint
    integer(IK), intent(in) :: maxfun
    integer(IK), intent(in) :: npt
    real(RP), intent(in) :: eta1
    real(RP), intent(in) :: eta2
    real(RP), intent(in) :: ftarget
    real(RP), intent(in) :: gamma1
    real(RP), intent(in) :: gamma2
    real(RP), intent(in) :: rhobeg
    real(RP), intent(in) :: rhoend
    real(RP), intent(in) :: xl(:)  ! XL(N)
    real(RP), intent(in) :: xu(:)  ! XU(N)
    
    ! In-outputs
    real(RP), intent(inout) :: x(:)  ! X(N)
    
    ! Outputs
    integer(IK), intent(out) :: info
    integer(IK), intent(out) :: nf
    real(RP), intent(out) :: f
    real(RP), intent(out) :: fhist(:)  ! FHIST(MAXFHIST)
    real(RP), intent(out) :: xhist(:, :)  ! XHIST(N, MAXXHIST)
    
    ! Local variables
    character(len=*), parameter :: solver = 'BOBYQA'
    character(len=*), parameter :: srname = 'BOBYQB'
    integer(IK) :: ij(2, max(0_IK, int(npt - 2 * size(x) - 1, IK)))
    integer(IK) :: itest
    integer(IK) :: k
    integer(IK) :: knew_geo
    integer(IK) :: knew_tr
    integer(IK) :: kopt
    integer(IK) :: maxfhist
    integer(IK) :: maxhist
    integer(IK) :: maxtr
    integer(IK) :: maxxhist
    integer(IK) :: n
    integer(IK) :: subinfo
    integer(IK) :: tr
    logical :: accurate_mod
    logical :: adequate_geo
    logical :: bad_trstep
    logical :: close_itpset
    logical :: improve_geo
    logical :: reduce_rho
    logical :: rescued
    logical :: shortd
    logical :: small_trrad
    logical :: terminate
    logical :: trfail
    logical :: ximproved
    real(RP) :: bmat(size(x), npt + size(x))
    real(RP) :: crvmin
    real(RP) :: d(size(x))
    real(RP) :: delbar
    real(RP) :: delta
    real(RP) :: den(npt)
    real(RP) :: distsq(npt)
    real(RP) :: dnorm
    real(RP) :: dnorm_rec(2)  ! Powell's implementation: DNORM_REC(3)
    real(RP) :: ebound
    real(RP) :: fval(npt)
    real(RP) :: gamma3
    real(RP) :: gopt(size(x))
    real(RP) :: hq(size(x), size(x))
    real(RP) :: moderr
    real(RP) :: moderr_rec(size(dnorm_rec))
    real(RP) :: pq(npt)
    real(RP) :: qred
    real(RP) :: ratio
    real(RP) :: rho
    real(RP) :: sl(size(x))
    real(RP) :: su(size(x))
    real(RP) :: vlag(npt + size(x))
    real(RP) :: xbase(size(x))
    real(RP) :: xdrop(size(x))
    real(RP) :: xosav(size(x))
    real(RP) :: xpt(size(x), npt)
    real(RP) :: zmat(npt, npt - size(x) - 1)
    real(RP), parameter :: trtol = 1.0E-2_RP  ! Convergence tolerance of trust-region subproblem solver
    
    ! Sizes.
    n = int(size(x), kind(n))
    maxxhist = int(size(xhist, 2), kind(maxxhist))
    maxfhist = int(size(fhist), kind(maxfhist))
    maxhist = int(max(maxxhist, maxfhist), kind(maxhist))
    
    ! Preconditions
    if (DEBUGGING) then
        call assert(abs(iprint) <= 3, 'IPRINT is 0, 1, -1, 2, -2, 3, or -3', srname)
        call assert(n >= 1, 'N >= 1', srname)
        call assert(npt >= n + 2, 'NPT >= N+2', srname)
        call assert(maxfun >= npt + 1, 'MAXFUN >= NPT+1', srname)
        call assert(eta1 >= 0 .and. eta1 <= eta2 .and. eta2 < 1, '0 <= ETA1 <= ETA2 < 1', srname)
        call assert(gamma1 > 0 .and. gamma1 < 1 .and. gamma2 > 1, '0 < GAMMA1 < 1 < GAMMA2', srname)
        call assert(rhobeg >= rhoend .and. rhoend > 0, 'RHOBEG >= RHOEND > 0', srname)
        call assert(size(xl) == n .and. size(xu) == n, 'SIZE(XL) == N == SIZE(XU)', srname)
        call assert(all(rhobeg <= (xu - xl) / TWO), 'RHOBEG <= MINVAL(XU-XL)/2', srname)
        call assert(all(is_finite(x)), 'X is finite', srname)
        call assert(all(x >= xl .and. (x <= xl .or. x - xl >= rhobeg)), 'X == XL or X - XL >= RHOBEG', srname)
        call assert(all(x <= xu .and. (x >= xu .or. xu - x >= rhobeg)), 'X == XU or XU - X >= RHOBEG', srname)
        call assert(maxhist >= 0 .and. maxhist <= maxfun, '0 <= MAXHIST <= MAXFUN', srname)
        call assert(size(xhist, 1) == n .and. maxxhist * (maxxhist - maxhist) == 0, &
            & 'SIZE(XHIST, 1) == N, SIZE(XHIST, 2) == 0 or MAXHIST', srname)
        call assert(maxfhist * (maxfhist - maxhist) == 0, 'SIZE(FHIST) == 0 or MAXHIST', srname)
    end if
    
    !====================!
    ! Calculation starts !
    !====================!
    
    ! Initialize XBASE, XPT, SL, SU, FVAL, and KOPT, together with the history, NF, and IJ.
    call initxf(calfun, iprint, maxfun, ftarget, rhobeg, xl, xu, x, ij, kopt, nf, fhist, fval, &
        & sl, su, xbase, xhist, xpt, subinfo)
    
    ! Report the current best value, and check if user asks for early termination.
    terminate = .false.
    if (present(callback_fcn)) then
        !call callback_fcn(xbase + xpt(:, kopt), fval(kopt), nf, 0_IK, terminate=terminate)
        if (terminate) then
            subinfo = CALLBACK_TERMINATE
        end if
    end if
    
    ! Initialize X and F according to KOPT.
    x = xinbd(xbase, xpt(:, kopt), xl, xu, sl, su)  ! In precise arithmetic, X = XBASE + XOPT.
    f = fval(kopt)
    
    ! Finish the initialization if INITXF completed normally and CALLBACK did not request termination;
    ! otherwise, do not proceed, as XPT etc may be uninitialized, leading to errors or exceptions.
    if (subinfo == INFO_DFT) then
        ! Initialize [BMAT, ZMAT], representing inverse of KKT matrix of the interpolation system.
        call inith(ij, xpt, bmat, zmat)
    
        ! Initialize the quadratic represented by [GOPT, HQ, PQ], so that its gradient at XBASE+XOPT is
        ! GOPT; its Hessian is HQ + sum_{K=1}^NPT PQ(K)*XPT(:, K)*XPT(:, K)'.
        call initq(ij, fval, xpt, gopt, hq, pq)
        if (.not. (all(is_finite(gopt)) .and. all(is_finite(hq)) .and. all(is_finite(pq)))) then
            subinfo = NAN_INF_MODEL
        end if
    end if
    
    ! Check whether to return due to abnormal cases that may occur during the initialization.
    if (subinfo /= INFO_DFT) then
        info = subinfo
        ! Arrange FHIST and XHIST so that they are in the chronological order.
        call rangehist(nf, xhist, fhist)
        ! Print a return message according to IPRINT.
        !call retmsg(solver, info, iprint, nf, f, x)
        ! Postconditions
        if (DEBUGGING) then
            call assert(nf <= maxfun, 'NF <= MAXFUN', srname)
            call assert(size(x) == n .and. .not. any(is_nan(x)), 'SIZE(X) == N, X does not contain NaN', srname)
            call assert(all(x >= xl) .and. all(x <= xu), 'XL <= X <= XU', srname)
            call assert(.not. (is_nan(f) .or. is_posinf(f)), 'F is not NaN/+Inf', srname)
            call assert(size(xhist, 1) == n .and. size(xhist, 2) == maxxhist, 'SIZE(XHIST) == [N, MAXXHIST]', srname)
            call assert(.not. any(is_nan(xhist(:, 1:min(nf, maxxhist)))), 'XHIST does not contain NaN', srname)
            ! The last calculated X can be Inf (finite + finite can be Inf numerically).
            do k = 1, min(nf, maxxhist)
                call assert(all(xhist(:, k) >= xl) .and. all(xhist(:, k) <= xu), 'XL <= XHIST <= XU', srname)
            end do
            call assert(size(fhist) == maxfhist, 'SIZE(FHIST) == MAXFHIST', srname)
            call assert(.not. any(is_nan(fhist(1:min(nf, maxfhist))) .or. is_posinf(fhist(1:min(nf, maxfhist)))), &
                & 'FHIST does not contain NaN/+Inf', srname)
            call assert(.not. any(fhist(1:min(nf, maxfhist)) < f), 'F is the smallest in FHIST', srname)
        end if
        return
    end if
    
    ! Set some more initial values.
    ! We must initialize RATIO. Otherwise, when SHORTD = TRUE, compilers may raise a run-time error that
    ! RATIO is undefined. But its value will not be used: when SHORTD = FALSE, its value will be
    ! overwritten; when SHORTD = TRUE, its value is used only in BAD_TRSTEP, which is TRUE regardless of
    ! RATIO. Similar for KNEW_TR.
    ! No need to initialize SHORTD unless MAXTR < 1, but some compilers may complain if we do not do it.
    rho = rhobeg
    delta = rho
    ebound = ZERO
    rescued = .false.
    shortd = .false.
    trfail = .false.
    ratio = -ONE
    dnorm_rec = REALMAX
    moderr_rec = REALMAX
    knew_tr = 0
    knew_geo = 0
    itest = 0
    
    ! If DELTA <= GAMMA3*RHO after an update, we set DELTA to RHO. GAMMA3 must be less than GAMMA2. The
    ! reason is as follows. Imagine a very successful step with DENORM = the un-updated DELTA = RHO.
    ! Then TRRAD will update DELTA to GAMMA2*RHO. If GAMMA3 >= GAMMA2, then DELTA will be reset to RHO,
    ! which is not reasonable as D is very successful. See paragraph two of Sec. 5.2.5 in
    ! T. M. Ragonneau's thesis: "Model-Based Derivative-Free Optimization Methods and Software".
    ! According to test on 20230613, for BOBYQA, this Powellful updating scheme of DELTA works better
    ! than setting directly DELTA = MAX(NEW_DELTA, RHO).
    gamma3 = max(ONE, min(0.75_RP * gamma2, 1.5_RP))
    
    ! MAXTR is the maximal number of trust-region iterations. Here, we set it to HUGE(MAXTR) so that
    ! the algorithm will not terminate due to MAXTR. However, this may not be allowed in other languages
    ! such as MATLAB. In that case, we can set MAXTR to 10*MAXFUN, which is unlikely to reach because
    ! each trust-region iteration takes 1 or 2 function evaluations unless the trust-region step is short
    ! or fails to reduce the trust-region model but the geometry step is not invoked.
    ! N.B.: Do NOT set MAXTR to HUGE(MAXTR), as it may cause overflow and infinite cycling in the DO
    ! loop. See
    ! https://fortran-lang.discourse.group/t/loop-variable-reaching-integer-huge-causes-infinite-loop
    ! https://fortran-lang.discourse.group/t/loops-dont-behave-like-they-should
    maxtr = huge(maxtr) - 1_IK  !!MATLAB: maxtr = 10 * maxfun;
    info = MAXTR_REACHED
    
    ! Begin the iterative procedure.
    ! After solving a trust-region subproblem, we use three boolean variables to control the workflow.
    ! SHORTD: Is the trust-region trial step too short to invoke a function evaluation?
    ! IMPROVE_GEO: Should we improve the geometry?
    ! REDUCE_RHO: Should we reduce rho?
    ! BOBYQA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
    do tr = 1, maxtr
        ! Generate the next trust region step D.
        call trsbox(delta, gopt, hq, pq, sl, su, trtol, xpt(:, kopt), xpt, crvmin, d)
        dnorm = min(delta, norm(d))
        shortd = (dnorm <= HALF * rho)  ! `<=` works better than `<` in case of underflow.
    
        ! Set QRED to the reduction of the quadratic model when the move D is made from XOPT. QRED
        ! should be positive. If it is nonpositive due to rounding errors, we will not take this step.
        qred = -quadinc(d, xpt, gopt, pq, hq)  ! QRED = Q(XOPT) - Q(XOPT + D)
        trfail = (.not. qred > 1.0E-6 * rho**2)  ! QRED is tiny/negative or NaN.
    
        ! When D is short, make a choice between reducing RHO and improving the geometry depending
        ! on whether or not our work with the current RHO seems complete. RHO is reduced if the
        ! errors in the quadratic model at the recent interpolation points compare favourably
        ! with predictions of likely improvements to the model within distance HALF*RHO of XOPT.
        ! Why do we reduce RHO when SHORTD is true and the entries of MODERR_REC and DNORM_REC are all
        ! small? The reason is well explained by the BOBYQA paper in the paragraphs surrounding
        ! (6.8)--(6.11). Roughly speaking, in this case, a trust-region step is unlikely to decrease the
        ! objective function according to some estimations. This suggests that the current trust-region
        ! center may be an approximate local minimizer up to the current "resolution" of the algorithm.
        ! When this occurs, the algorithm takes the view that the work for the current RHO is complete,
        ! and hence it will reduce RHO, which will enhance the resolution of the algorithm in general.
        if (shortd .or. trfail) then
            delta = TENTH * delta
            if (delta <= gamma3 * rho) then
                delta = rho  ! Set DELTA to RHO when it is close to or below.
            end if
            ! Evaluate EBOUND. It will be used as a bound to test if the entries of MODERR_REC are small.
            ebound = errbd(crvmin, d, gopt, hq, moderr_rec, pq, rho, sl, su, xpt(:, kopt), xpt)
        else
            ! Calculate the next value of the objective function.
            x = xinbd(xbase, xpt(:, kopt) + d, xl, xu, sl, su)  ! X = XBASE + XOPT + D without rounding.
            call evaluate(calfun, x, f)
            nf = nf + 1_IK
            rescued = .false.  ! Set RESCUED to FALSE after evaluating F at a new point.
    
            ! Print a message about the function evaluation according to IPRINT.
            !call fmsg(solver, 'Trust region', iprint, nf, delta, f, x)
            ! Save X, F into the history.
            call savehist(nf, x, xhist, f, fhist)
    
            ! Check whether to exit
            subinfo = checkexit(maxfun, nf, f, ftarget, x)
            if (subinfo /= INFO_DFT) then
                info = subinfo
                exit
            end if
    
            ! Update DNORM_REC and MODERR_REC.
            ! DNORM_REC records the DNORM of the recent function evaluations with the current RHO.
            dnorm_rec = [dnorm_rec(2:size(dnorm_rec)), dnorm]
            ! MODERR is the error of the current model in predicting the change in F due to D.
            ! MODERR_REC records the prediction errors of the recent models with the current RHO.
            moderr = f - fval(kopt) + qred
            moderr_rec = [moderr_rec(2:size(moderr_rec)), moderr]
    
            ! Calculate the reduction ratio by REDRAT, which handles Inf/NaN carefully.
            ratio = redrat(fval(kopt) - f, qred, eta1)
    
            ! Update DELTA. After this, DELTA < DNORM may hold.
            delta = trrad(delta, dnorm, eta1, eta2, gamma1, gamma2, ratio)
            if (delta <= gamma3 * rho) then
                delta = rho  ! Set DELTA to RHO when it is close to or below.
            end if
    
            ! Is the newly generated X better than current best point?
            ximproved = (f < fval(kopt))
    
            ! Call RESCUE if rounding errors have damaged the denominator corresponding to D.
            ! RESCUE is invoked sometimes though not often after a trust-region step, and it does
            ! improve the performance, especially when pursing high-precision solutions.
            vlag = calvlag(kopt, bmat, d, xpt, zmat)
            den = calden(kopt, bmat, d, xpt, zmat)
            if (ximproved .and. any(den > maxval(vlag(1:npt)**2))) then
                ! Below are some alternatives conditions for calling RESCUE. They perform fairly well.
                ! !if (.false.) then  ! Do not call RESCUE at all.
                ! !if (ximproved .and. .not. any(den > 0.25_RP * maxval(vlag(1:npt)**2))) then
                ! !if (ximproved .and. .not. any(den > HALF * maxval(vlag(1:npt)**2))) then
                ! !if (.not. any(den > HALF * maxval(vlag(1:npt)**2))) then  ! Powell's code.
                ! !if (.not. any(den > maxval(vlag(1:npt)**2))) then
                if (rescued) then
                    info = DAMAGING_ROUNDING  ! The last RESCUE did not improve the situation.
                    exit
                end if
                call rescue(calfun, solver, iprint, maxfun, delta, ftarget, xl, xu, kopt, nf, fhist, &
                    & fval, gopt, hq, pq, sl, su, xbase, xhist, xpt, bmat, zmat, subinfo)
                ! if (subinfo /= INFO_DFT) then
                !     info = subinfo
                !     exit
                ! end if
                ! rescued = .true.
                ! dnorm_rec = REALMAX
                ! moderr_rec = REALMAX
    
                ! ! RESCUE shifts XBASE to the best point before RESCUE. Update D, MODERR, and XIMPROVED.
                ! ! Do NOT calculate QRED according to this D, as it is not really a trust region step.
                ! ! Note that QRED will be used afterward for defining IMPROVE_GEO and REDUCE_RHO.
                ! d = max(sl, min(su, d)) - xpt(:, kopt)
                ! moderr = f - fval(kopt) - quadinc(d, xpt, gopt, pq, hq)
                ! ximproved = (f < fval(kopt))
            end if
    
            ! Set KNEW_TR to the index of the interpolation point to be replaced with XOPT + D.
            ! KNEW_TR will ensure that the geometry of XPT is "good enough" after the replacement.
            ! knew_tr = setdrop_tr(kopt, ximproved, bmat, d, delta, rho, xpt, zmat)
    
            ! Update [BMAT, ZMAT] (representing H in the BOBYQA paper), [GQ, HQ, PQ] (the quadratic
            ! model), and [FVAL, XPT, KOPT, FOPT, XOPT] so that XPT(:, KNEW_TR) becomes XOPT + D. If
            ! KNEW_TR = 0, the updating subroutines will do essentially nothing, as the algorithm
            ! decides not to include XOPT + D into XPT.
            ! if (knew_tr > 0) then
            !     xdrop = xpt(:, knew_tr)
            !     xosav = xpt(:, kopt)
            !     call updateh(knew_tr, kopt, d, xpt, bmat, zmat)
            !     call updatexf(knew_tr, ximproved, f, max(sl, min(su, xosav + d)), kopt, fval, xpt)
            !     call updateq(knew_tr, ximproved, bmat, d, moderr, xdrop, xosav, xpt, zmat, gopt, hq, pq)
            !     ! Try whether to replace the new quadratic model with the alternative model, namely the
            !     ! least Frobenius norm interpolant.
            !     call tryqalt(bmat, fval - fval(kopt), ratio, sl, su, xpt(:, kopt), xpt, zmat, itest, gopt, hq, pq)
            !     if (.not. (all(is_finite(gopt)) .and. all(is_finite(hq)) .and. all(is_finite(pq)))) then
            !         info = NAN_INF_MODEL
            !         exit
            !     end if
            ! end if
        end if  ! End of IF (SHORTD .OR. TRFAIL). The normal trust-region calculation ends.
    
    
        !----------------------------------------------------------------------------------------------!
        ! Before the next trust-region iteration, we may improve the geometry of XPT or reduce RHO
        ! according to IMPROVE_GEO and REDUCE_RHO, which in turn depend on the following indicators.
        ! N.B.: We must ensure that the algorithm does not set IMPROVE_GEO = TRUE at infinitely many
        ! consecutive iterations without moving XOPT or reducing RHO. Otherwise, the algorithm will get
        ! stuck in repetitive invocations of GEOSTEP. To this end, make sure the following.
        ! 1. The threshold for CLOSE_ITPSET is at least DELBAR, the trust region radius for GEOSTEP.
        ! Normally, DELBAR <= DELTA <= the threshold (In Powell's UOBYQA, DELBAR = RHO < the threshold).
        ! 2. If an iteration sets IMPROVE_GEO = TRUE, it must also reduce DELTA or set DELTA to RHO.
    
        ! ACCURATE_MOD: Are the recent models sufficiently accurate? Used only if SHORTD is TRUE.
        !accurate_mod = all(abs(moderr_rec) <= ebound) .and. all(dnorm_rec <= rho)
        ! CLOSE_ITPSET: Are the interpolation points close to XOPT?
        !distsq = sum((xpt - spread(xpt(:, kopt), dim=2, ncopies=npt))**2, dim=1)
        !!MATLAB: distsq = sum((xpt - xpt(:, kopt)).^2)  % Implicit expansion
        ! close_itpset = all(distsq <= max(delta**2, (TEN * rho)**2))
        ! ! Below are some alternative definitions of CLOSE_ITPSET.
        ! ! N.B.: The threshold for CLOSE_ITPSET is at least DELBAR, the trust region radius for GEOSTEP.
        ! ! !close_itpset = all(distsq <= max((TWO * delta)**2, (TEN * rho)**2))  ! Powell's code.
        ! ! !close_itpset = all(distsq <= 4.0_RP * delta**2)  ! Powell's NEWUOA code.
        ! ! !close_itpset = all(distsq <= max(delta**2, 4.0_RP * rho**2))  ! Powell's LINCOA code.
        ! ! ADEQUATE_GEO: Is the geometry of the interpolation set "adequate"?
        ! ! N.B. (Zaikun 20240314): Even if RESCUE has just been called (RESCUED = TRUE), the geometry may
        ! ! still be inadequate/improvable if XPT contains points far away from XOPT.
        ! adequate_geo = (shortd .and. accurate_mod) .or. close_itpset
        ! ! SMALL_TRRAD: Is the trust-region radius small? This indicator seems not impactive in practice.
        ! small_trrad = (max(delta, dnorm) <= rho)  ! Powell's code. See also (6.7) of the BOBYQA paper.
        ! !small_trrad = (delsav <= rho)  ! Behaves the same as Powell's version. DELSAV = unupdated DELTA.
    
        ! ! IMPROVE_GEO and REDUCE_RHO are defined as follows.
        ! ! N.B.: If SHORTD is TRUE at the very first iteration, then REDUCE_RHO will be set to TRUE.
        ! ! Powell's code does not have TRFAIL in BAD_TRSTEP; it terminates if TRFAIL is TRUE.
    
        ! ! BAD_TRSTEP (for IMPROVE_GEO): Is the last trust-region step bad?
        ! bad_trstep = (shortd .or. trfail .or. ratio <= eta1 .or. knew_tr == 0)
        ! improve_geo = bad_trstep .and. .not. adequate_geo  ! See the text above (6.7) of the BOBYQA paper.
        ! ! BAD_TRSTEP (for REDUCE_RHO): Is the last trust-region step bad?
        ! bad_trstep = (shortd .or. trfail .or. ratio <= 0 .or. knew_tr == 0)
        ! reduce_rho = bad_trstep .and. adequate_geo .and. small_trrad  ! See (6.7) of the BOBYQA paper.
        ! Zaikun 20221111: What if RESCUE has been called? Is it still reasonable to use RATIO?
        ! Zaikun 20221127: If RESCUE has been called, then KNEW_TR may be 0 even if RATIO > 0.
    
        ! Equivalently, REDUCE_RHO can be set as follows. It shows that REDUCE_RHO is TRUE in two cases.
        ! !bad_trstep = (shortd .or. trfail .or. ratio <= 0 .or. knew_tr == 0)
        ! !reduce_rho = (shortd .and. accurate_mod) .or. (bad_trstep .and. close_itpset .and. small_trrad)
    
        ! With REDUCE_RHO properly defined, we can also set IMPROVE_GEO as follows.
        ! !bad_trstep = (shortd .or. trfail .or. ratio <= eta1 .or. knew_tr == 0)
        ! !improve_geo = bad_trstep .and. (.not. reduce_rho) .and. (.not. close_itpset)
    
        ! With IMPROVE_GEO properly defined, we can also set REDUCE_RHO as follows.
        ! !bad_trstep = (shortd .or. trfail .or. ratio <= 0 .or. knew_tr == 0)
        ! !reduce_rho = bad_trstep .and. (.not. improve_geo) .and. small_trrad
    
        ! BOBYQA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
        !call assert(.not. (improve_geo .and. reduce_rho), 'IMPROVE_GEO and REDUCE_RHO are not both TRUE', srname)
        !
        ! If SHORTD or TRFAIL is TRUE, then either IMPROVE_GEO or REDUCE_RHO is TRUE unless CLOSE_ITPSET
        ! is TRUE but SMALL_TRRAD is FALSE.
        !call assert((.not. (shortd .or. trfail)) .or. (improve_geo .or. reduce_rho .or. &
        !    & (close_itpset .and. .not. small_trrad)), 'If SHORTD or TRFAIL is TRUE, then either &
        !    & IMPROVE_GEO or REDUCE_RHO is TRUE unless CLOSE_ITPSET is TRUE but SMALL_TRRAD is FALSE', srname)
        !----------------------------------------------------------------------------------------------!
    
    
        ! Since IMPROVE_GEO and REDUCE_RHO are never TRUE simultaneously, the following two blocks are
        ! exchangeable: IF (IMPROVE_GEO) ... END IF and IF (REDUCE_RHO) ... END IF.
    
        ! Improve the geometry of the interpolation set by removing a point and adding a new one.
        if (improve_geo) then
            ! XPT(:, KNEW_GEO) will become XOPT + D below. KNEW_GEO /= KOPT unless there is a bug.
            knew_geo = int(maxloc(distsq, dim=1), kind(knew_geo))
    
            ! Set DELBAR, which will be used as the trust-region radius for the geometry-improving
            ! scheme GEOSTEP. Note that DELTA has been updated before arriving here.
            delbar = max(min(TENTH * sqrt(maxval(distsq)), delta), rho)  ! Powell's code
            !delbar = rho  ! Powell's UOBYQA code
            !delbar = max(min(TENTH * sqrt(maxval(distsq)), HALF * delta), rho)  ! Powell's NEWUOA code
            !delbar = max(TENTH * delta, rho)  ! Powell's LINCOA code
    
            ! Find D so that the geometry of XPT will be improved when XPT(:, KNEW_GEO) becomes XOPT + D.
            d = geostep(knew_geo, kopt, bmat, delbar, sl, su, xpt, zmat)
    
            ! Call RESCUE if rounding errors have damaged the denominator corresponding to D.
            ! 1. This does make a difference, yet RESCUE seems not invoked often after a geometry step.
            ! 2. In Powell's implementation, it may happen that RESCUE only recalculates [BMAT, ZMAT]
            ! without introducing any new point into XPT. In that case, GEOSTEP will have to be called
            ! after RESCUE, without which the code may encounter an infinite cycling. We have modified
            ! RESCUE so that it introduces at least one new point into XPT and there is no need to call
            ! GEOSTEP afterward. This improves the performance a bit and simplifies the flow of the code.
            ! 3. It is tempting to incorporate XOPT+D into the interpolation even if RESCUE is called.
            ! However, this cannot be done without recalculating KNEW_GEO, as XPT has been changed by
            ! RESCUE, so that it is invalid to replace XPT(:, KNEW_GEO) with XOPT+D anymore. With a new
            ! KNEW_GEO, the step D will become improper as it was chosen according to the old KNEW_GEO.
            ! vlag = calvlag(kopt, bmat, d, xpt, zmat)
            ! den = calden(kopt, bmat, d, xpt, zmat)
            if (den(knew_geo) > HALF * vlag(knew_geo)**2) then
                if (rescued) then
                    info = DAMAGING_ROUNDING  ! The last RESCUE did not improve the situation.
                    exit
                end if
                call rescue(calfun, solver, iprint, maxfun, delta, ftarget, xl, xu, kopt, nf, fhist, &
                    & fval, gopt, hq, pq, sl, su, xbase, xhist, xpt, bmat, zmat, subinfo)
                if (subinfo /= INFO_DFT) then
                    info = subinfo
                    exit
                end if
                rescued = .true.
                dnorm_rec = REALMAX
                moderr_rec = REALMAX
            else
                ! Calculate the next value of the objective function.
                x = xinbd(xbase, xpt(:, kopt) + d, xl, xu, sl, su)  ! X = XBASE + XOPT + D without rounding.
                call evaluate(calfun, x, f)
                nf = nf + 1_IK
                rescued = .false.  ! Set RESCUED to FALSE after evaluating F at a new point.
    
                ! Print a message about the function evaluation according to IPRINT.
                !call fmsg(solver, 'Geometry', iprint, nf, delbar, f, x)
                ! Save X, F into the history.
                call savehist(nf, x, xhist, f, fhist)
    
                ! Check whether to exit
                subinfo = checkexit(maxfun, nf, f, ftarget, x)
                if (subinfo /= INFO_DFT) then
                    info = subinfo
                    exit
                end if
    
                ! Update DNORM_REC and MODERR_REC.
                ! DNORM_REC records the DNORM of the recent function evaluations with the current RHO.
                ! Powell's code does not update DNORM. Therefore, DNORM is the length of the last
                ! trust-region trial step, inconsistent with MODERR_REC. The same problem exists in NEWUOA.
                dnorm = min(delbar, norm(d))
                dnorm_rec = [dnorm_rec(2:size(dnorm_rec)), dnorm]
                ! MODERR is the error of the current model in predicting the change in F due to D.
                ! MODERR_REC records the prediction errors of the recent models with the current RHO.
                moderr = f - fval(kopt) - quadinc(d, xpt, gopt, pq, hq)  ! QRED = Q(XOPT) - Q(XOPT + D)
                moderr_rec = [moderr_rec(2:size(moderr_rec)), moderr]
    
                ! Is the newly generated X better than current best point?
                ximproved = (f < fval(kopt))
    
                ! Update [BMAT, ZMAT] (represents H in the BOBYQA paper), [FVAL, XPT, KOPT, FOPT, XOPT],
                ! and [GQ, HQ, PQ] (the quadratic model), so that XPT(:, KNEW_GEO) becomes XOPT + D.
                xdrop = xpt(:, knew_geo)
                xosav = xpt(:, kopt)
                call updateh(knew_geo, kopt, d, xpt, bmat, zmat)
                call updatexf(knew_geo, ximproved, f, max(sl, min(su, xosav + d)), kopt, fval, xpt)
                call updateq(knew_geo, ximproved, bmat, d, moderr, xdrop, xosav, xpt, zmat, gopt, hq, pq)
                if (.not. (all(is_finite(gopt)) .and. all(is_finite(hq)) .and. all(is_finite(pq)))) then
                    info = NAN_INF_MODEL
                    exit
                end if
            end if
        end if  ! End of IF (IMPROVE_GEO). The procedure of improving geometry ends.
    
        ! The calculations with the current RHO are complete. Enhance the resolution of the algorithm
        ! by reducing RHO; update DELTA at the same time.
        if (reduce_rho) then
            if (rho <= rhoend) then
                info = SMALL_TR_RADIUS
                exit
            end if
            delta = max(HALF * rho, redrho(rho, rhoend))
            rho = redrho(rho, rhoend)
            ! Print a message about the reduction of RHO according to IPRINT.
            !call rhomsg(solver, iprint, nf, delta, fval(kopt), rho, xbase + xpt(:, kopt))
            ! DNORM_REC and MODERR_REC are corresponding to the recent function evaluations with
            ! the current RHO. Update them after reducing RHO.
            dnorm_rec = REALMAX
            moderr_rec = REALMAX
        end if  ! End of IF (REDUCE_RHO). The procedure of reducing RHO ends.
    
        ! Shift XBASE if XOPT may be too far from XBASE.
        ! Powell's original criteria for shifting XBASE is as follows.
        ! 1. After a trust region step that is not short, shift XBASE if SUM(XOPT**2) >= 1.0E3*DNORM**2.
        ! In this case, it seems quite important for the performance to recalculate QRED.
        ! 2. Before a geometry step, shift XBASE if SUM(XOPT**2) >= 1.0E3*DELBAR**2.
        if (sum(xpt(:, kopt)**2) >= 1.0E3_RP * delta**2) then
            ! Other possible criteria: SUM(XOPT**2) >= 1.0E4*DELTA**2, SUM(XOPT**2) >= 1.0E4*RHO**2.
            sl = min(sl - xpt(:, kopt), ZERO)
            su = max(su - xpt(:, kopt), ZERO)
            call shiftbase(kopt, xbase, xpt, zmat, bmat, pq, hq)
            xbase = max(xl, min(xu, xbase))
        end if
    
        ! Report the current best value, and check if user asks for early termination.
        if (present(callback_fcn)) then
            !call callback_fcn(xbase + xpt(:, kopt), fval(kopt), nf, tr, terminate=terminate)
            if (terminate) then
                info = CALLBACK_TERMINATE
                exit
            end if
        end if
    
    end do  ! End of DO TR = 1, MAXTR. The iterative procedure ends.
    
    ! Return from the calculation, after trying the Newton-Raphson step if it has not been tried yet.
    if (info == SMALL_TR_RADIUS .and. shortd .and. dnorm > TENTH * rhoend .and. nf < maxfun) then
        x = xinbd(xbase, xpt(:, kopt) + d, xl, xu, sl, su)  ! In precise arithmetic, X = XBASE + XOPT + D.
        call evaluate(calfun, x, f)
        nf = nf + 1_IK
        ! Print a message about the function evaluation according to IPRINT.
        ! Zaikun 20230512: DELTA has been updated. RHO is only indicative here. TO BE IMPROVED.
        !call fmsg(solver, 'Trust region', iprint, nf, rho, f, x)
        ! Save X, F into the history.
        call savehist(nf, x, xhist, f, fhist)
    end if
    
    ! Choose the [X, F] to return: either the current [X, F] or [XBASE + XOPT, FOPT].
    if (fval(kopt) < f .or. is_nan(f)) then
        x = xinbd(xbase, xpt(:, kopt), xl, xu, sl, su)  ! In precise arithmetic, X = XBASE + XOPT.
        f = fval(kopt)
    end if
    
    ! Arrange FHIST and XHIST so that they are in the chronological order.
    call rangehist(nf, xhist, fhist)
    
    ! Print a return message according to IPRINT.
    call retmsg(solver, info, iprint, nf, f, x)
    
    !====================!
    !  Calculation ends  !
    !====================!
    
    ! Postconditions
    if (DEBUGGING) then
        call assert(nf <= maxfun, 'NF <= MAXFUN', srname)
        call assert(size(x) == n .and. .not. any(is_nan(x)), 'SIZE(X) == N, X does not contain NaN', srname)
        call assert(all(x >= xl) .and. all(x <= xu), 'XL <= X <= XU', srname)
        call assert(.not. (is_nan(f) .or. is_posinf(f)), 'F is not NaN/+Inf', srname)
        call assert(size(xhist, 1) == n .and. size(xhist, 2) == maxxhist, 'SIZE(XHIST) == [N, MAXXHIST]', srname)
        call assert(.not. any(is_nan(xhist(:, 1:min(nf, maxxhist)))), 'XHIST does not contain NaN', srname)
        ! The last calculated X can be Inf (finite + finite can be Inf numerically).
        do k = 1, min(nf, maxxhist)
            call assert(all(xhist(:, k) >= xl) .and. all(xhist(:, k) <= xu), 'XL <= XHIST <= XU', srname)
        end do
        call assert(size(fhist) == maxfhist, 'SIZE(FHIST) == MAXFHIST', srname)
        call assert(.not. any(is_nan(fhist(1:min(nf, maxfhist))) .or. is_posinf(fhist(1:min(nf, maxfhist)))), &
            & 'FHIST does not contain NaN/+Inf', srname)
        call assert(.not. any(fhist(1:min(nf, maxfhist)) < f), 'F is the smallest in FHIST', srname)
    end if
    
    end subroutine bobyqb
    
    
    function errbd(crvmin, d, gopt, hq, moderr_rec, pq, rho, sl, su, xopt, xpt) result(ebound)
    !--------------------------------------------------------------------------------------------------!
    ! This function defines EBOUND, which will be used as a bound to test whether the errors in recent
    ! models are sufficiently small. See the elaboration on pages 30--31 of the BOBYQA paper, in the
    ! paragraphs surrounding (6.8)--(6.11).
    !--------------------------------------------------------------------------------------------------!
    
    ! Common modules
    use, non_intrinsic :: consts_mod, only : RP, IK, HALF, DEBUGGING
    use, non_intrinsic :: debug_mod, only : assert
    use, non_intrinsic :: infnan_mod, only : is_finite
    use, non_intrinsic :: linalg_mod, only : matprod, diag, issymmetric, trueloc
    use, non_intrinsic :: powalg_mod, only : hess_mul
    
    implicit none
    
    ! Inputs
    real(RP), intent(in) :: crvmin
    real(RP), intent(in) :: d(:)
    real(RP), intent(in) :: gopt(:)
    real(RP), intent(in) :: hq(:, :)
    real(RP), intent(in) :: moderr_rec(:)
    real(RP), intent(in) :: pq(:)
    real(RP), intent(in) :: rho
    real(RP), intent(in) :: sl(:)
    real(RP), intent(in) :: su(:)
    real(RP), intent(in) :: xopt(:)
    real(RP), intent(in) :: xpt(:, :)
    
    ! Outputs
    real(RP) :: ebound
    
    ! Local variables
    character(len=*), parameter :: srname = 'ERRBD'
    integer(IK) :: n
    integer(IK) :: npt
    integer(IK) :: i
    integer(IK) :: temp(size(d))
    real(RP) :: bfirst(size(d))
    real(RP) :: bsecond(size(d))
    real(RP) :: gnew(size(d))
    real(RP) :: xnew(size(d))
    
    ! Sizes
    n = int(size(xpt, 1), kind(n))
    npt = int(size(xpt, 2), kind(npt))
    
    ! Preconditions
    if (DEBUGGING) then
        call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2', srname)
        call assert(crvmin >= 0, 'CRVMIN >= 0', srname)
        call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N, D is finite', srname)
        call assert(size(gopt) == n, 'SIZE(GOPT) == N', srname)
        call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is n-by-n and symmetric', srname)
        call assert(size(pq) == npt, 'SIZE(PQ) == NPT', srname)
        call assert(rho > 0, 'RHO > 0', srname)
        call assert(size(sl) == n .and. size(su) == n, 'SIZE(SL) == N == SIZE(SU)', srname)
        call assert(size(xopt) == n .and. all(is_finite(xopt)), 'SIZE(XOPT) == N, XOPT is finite', srname)
        call assert(all(is_finite(xpt)), 'XPT is finite', srname)
        call assert(all(xopt >= sl .and. xopt <= su), 'SL <= XOPT <= SU', srname)
        !call assert(all(xpt >= spread(sl, dim=2, ncopies=npt) .and. &
        !    & xpt <= spread(su, dim=2, ncopies=npt)), 'SL <= XPT <= SU', srname)
    end if
    
    !====================!
    ! Calculation starts !
    !====================!
    
    xnew = xopt + d
    gnew = gopt + hess_mul(d, xpt, pq, hq)
    !bfirst = maxval(abs(moderr_rec))
    ! bfirst(trueloc(xnew <= sl)) = gnew(trueloc(xnew <= sl)) * rho
    ! bfirst(trueloc(xnew >= su)) = -gnew(trueloc(xnew >= su)) * rho
    ! bsecond = HALF * (diag(hq) + matprod(xpt**2, pq)) * rho**2
    end function errbd
    
    
    end module bobyqb_mod2
    
    
    module bobyqa_mod2
    !--------------------------------------------------------------------------------------------------!
    ! BOBYQA_MOD is a module providing the reference implementation of Powell's BOBYQA algorithm in
    !
    ! M. J. D. Powell, The BOBYQA algorithm for bound constrained optimization without derivatives,
    ! Technical Report DAMTP 2009/NA06, Department of Applied Mathematics and Theoretical Physics,
    ! Cambridge University, Cambridge, UK, 2009
    !
    ! BOBYQA approximately solves
    !
    !   min F(X) subject to XL <= X <= XU,
    !
    ! where X is a vector of variables that has N components, and F is a real-valued objective function.
    ! XL and XU are a pair of N-dimensional vectors indicating the lower and upper bounds of X. The
    ! algorithm assumes that XL < XU entrywise. It tackles the problem by applying a trust region method
    ! that forms quadratic models by interpolation. There is usually some freedom in the interpolation
    ! conditions, which is taken up by minimizing the Frobenius norm of the change to the second
    ! derivative of the model, beginning with the ZERO matrix. The values of the variables are
    ! constrained by upper and lower bounds. The arguments of the subroutine are as follows.
    !
    ! Coded by Zaikun ZHANG (www.zhangzk.net) based on the BOBYQA paper and Powell's code, with
    ! modernization, bug fixes, and improvements.
    !
    ! Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    !
    ! Started: February 2022
    !
    ! Last Modified: Thursday, February 22, 2024 PM03:30:31
    !--------------------------------------------------------------------------------------------------!
    
    implicit none
    private
    public :: bobyqa
    
    
    contains
    
    
    subroutine bobyqa(calfun, x, &
        & f, xl, xu, &
        & nf, rhobeg, rhoend, ftarget, maxfun, npt, iprint, &
        & eta1, eta2, gamma1, gamma2, xhist, fhist, maxhist, honour_x0, callback_fcn, info)
    !--------------------------------------------------------------------------------------------------!
    ! Among all the arguments, only CALFUN and X are obligatory. The others are OPTIONAL and you can
    ! neglect them unless you are familiar with the algorithm. Any unspecified optional input will take
    ! the default value detailed below. For instance, we may invoke the solver as follows.
    !
    ! ! First define CALFUN and X, and then do the following.
    ! call bobyqa(calfun, x, f)
    !
    ! or
    !
    ! ! First define CALFUN, X, and XL, and then do the following.
    ! call bobyqa(calfun, x, f, xl = xl, rhobeg = 1.0D0, rhoend = 1.0D-6)
    !
    ! See examples/bobyqa_exmp.f90 for a concrete example.
    !
    ! A detailed introduction to the arguments is as follows.
    ! N.B.: RP and IK are defined in the module CONSTS_MOD. See consts.F90 under the directory named
    ! "common". By default, RP = kind(0.0D0) and IK = kind(0), with REAL(RP) being the double-precision
    ! real, and INTEGER(IK) being the default integer. For ADVANCED USERS, RP and IK can be defined by
    ! setting PRIMA_REAL_PRECISION and PRIMA_INTEGER_KIND in common/ppf.h. Use the default if unsure.
    !
    ! CALFUN
    !   Input, subroutine.
    !   CALFUN(X, F) should evaluate the objective function at the given REAL(RP) vector X and set the
    !   value to the REAL(RP) scalar F. It must be provided by the user, and its definition must conform
    !   to the following interface:
    !   !-------------------------------------------------------------------------!
    !    subroutine calfun(x, f)
    !    real(RP), intent(in) :: x(:)
    !    real(RP), intent(out) :: f
    !    end subroutine calfun
    !   !-------------------------------------------------------------------------!
    !
    ! X
    !   Input and output, REAL(RP) vector.
    !   As an input, X should be an N dimensional vector that contains the starting point, N being the
    !   dimension of the problem. As an output, X will be set to an approximate minimizer.
    !
    ! F
    !   Output, REAL(RP) scalar.
    !   F will be set to the objective function value of X at exit.
    !
    ! XL, XU
    !   Input, REAL(RP) vectors, default: XL = [], XU = [].
    !   XL is the lower bound for X. Its size is either N or 0, the latter signifying that X has no
    !   lower bound. Any entry of XL that is NaN or below -BOUNDMAX will be taken as -BOUNDMAX, which
    !   effectively means there is no lower bound for the corresponding entry of X. The value of
    !   BOUNDMAX is 0.25*HUGE(X), which is about 8.6E37 for single precision and 4.5E307 for double
    !   precision. XU is similar.
    !   N.B.:
    !   1. It is required that XU - XL > 2*EPSILON(X), which is about 2.4E-7 for single precision and
    !   4.5E-16 for double precision. Otherwise, the solver will return after printing a warning.
    !   2. Why don't we set BOUNDMAX to REALMAX? Because we want to avoid overflow when calculating
    !   XU - XL and when defining/updating SU and SL. This is not a problem in MATLAB/Python/Julia/R.
    !
    ! NF
    !   Output, INTEGER(IK) scalar.
    !   NF will be set to the number of calls of CALFUN at exit.
    !
    ! RHOBEG, RHOEND
    !   Inputs, REAL(RP) scalars, default: RHOBEG = 1, RHOEND = 10^-6. RHOBEG and RHOEND must be set to
    !   the initial and final values of a trust-region radius, both being positive and RHOEND <= RHOBEG.
    !   Typically RHOBEG should be about one tenth of the greatest expected change to a variable, and
    !   RHOEND should indicate the accuracy that is required in the final values of the variables.
    !
    ! FTARGET
    !   Input, REAL(RP) scalar, default: -Inf.
    !   FTARGET is the target function value. The algorithm will terminate when a point with a function
    !   value <= FTARGET is found.
    !
    ! MAXFUN
    !   Input, INTEGER(IK) scalar, default: MAXFUN_DIM_DFT*N with MAXFUN_DIM_DFT defined in the module
    !   CONSTS_MOD (see common/consts.F90). MAXFUN is the maximal number of calls of CALFUN.
    !
    ! NPT
    !   Input, INTEGER(IK) scalar, default: 2N + 1.
    !   NPT is the number of interpolation conditions for each trust region model. Its value must be in
    !   the interval [N+2, (N+1)(N+2)/2]. Powell commented that "the value NPT = 2*N+1 being recommended
    !   for a start ... much larger values tend to be inefficient, because the amount of routine work of
    !   each iteration is of magnitude NPT**2, and because the achievement of adequate accuracy in some
    !   matrix calculations becomes more difficult. Some excellent numerical results have been found in
    !   the case NPT=N+6 even with more than 100 variables." And "choices that exceed 2*N+1 are not
    !   recommended" by Powell.
    !
    ! IPRINT
    !   Input, INTEGER(IK) scalar, default: 0.
    !   The value of IPRINT should be set to 0, 1, -1, 2, -2, 3, or -3, which controls how much
    !   information will be printed during the computation:
    !   0: there will be no printing;
    !   1: a message will be printed to the screen at the return, showing the best vector of variables
    !      found and its objective function value;
    !   2: in addition to 1, each new value of RHO is printed to the screen, with the best vector of
    !      variables so far and its objective function value;
    !   3: in addition to 2, each function evaluation with its variables will be printed to the screen;
    !   -1, -2, -3: the same information as 1, 2, 3 will be printed, not to the screen but to a file
    !      named BOBYQA_output.txt; the file will be created if it does not exist; the new output will
    !      be appended to the end of this file if it already exists.
    !   Note that IPRINT = +/-3 can be costly in terms of time and/or space.
    !
    ! ETA1, ETA2, GAMMA1, GAMMA2
    !   Input, REAL(RP) scalars, default: ETA1 = 0.1, ETA2 = 0.7, GAMMA1 = 0.5, and GAMMA2 = 2.
    !   ETA1, ETA2, GAMMA1, and GAMMA2 are parameters in the updating scheme of the trust-region radius
    !   detailed in the subroutine TRRAD in trustregion.f90. Roughly speaking, the trust-region radius
    !   is contracted by a factor of GAMMA1 when the reduction ratio is below ETA1, and enlarged by a
    !   factor of GAMMA2 when the reduction ratio is above ETA2. It is required that 0 < ETA1 <= ETA2
    !   < 1 and 0 < GAMMA1 < 1 < GAMMA2. Normally, ETA1 <= 0.25. It is NOT advised to set ETA1 >= 0.5.
    !
    ! XHIST, FHIST, MAXHIST
    !   XHIST: Output, ALLOCATABLE rank 2 REAL(RP) array;
    !   FHIST: Output, ALLOCATABLE rank 1 REAL(RP) array;
    !   MAXHIST: Input, INTEGER(IK) scalar, default: MAXFUN
    !   XHIST, if present, will output the history of iterates, while FHIST, if present, will output the
    !   history function values. MAXHIST should be a nonnegative integer, and XHIST/FHIST will output
    !   only the history of the last MAXHIST iterations. Therefore, MAXHIST = 0 means XHIST/FHIST will
    !   output nothing, while setting MAXHIST = MAXFUN requests XHIST/FHIST to output all the history.
    !   If XHIST is present, its size at exit will be [N, min(NF, MAXHIST)]; if FHIST is present, its
    !   size at exit will be min(NF, MAXHIST).
    !
    !   IMPORTANT NOTICE:
    !   Setting MAXHIST to a large value can be costly in terms of memory for large problems.
    !   MAXHIST will be reset to a smaller value if the memory needed exceeds MAXHISTMEM defined in
    !   CONSTS_MOD (see consts.F90 under the directory named "common").
    !   Use *HIST with caution! (N.B.: the algorithm is NOT designed for large problems).
    !
    ! HONOUR_X0
    !  Input, LOGICAL scalar, default: it is .false. if RHOBEG is present and 0 < RHOBEG < Inf, and it
    !  is .true. otherwise. HONOUR_X0 indicates whether to respect the user-defined X0 or not.
    !  BOBYQA requires that the distance between X0 and the inactive bounds is at least RHOBEG. X0 or
    !  RHOBEG is revised if this requirement is not met. If HONOUR_X0 == TRUE, revise RHOBEG if needed;
    !  otherwise, revise X0 if needed. See the PREPROC subroutine for more information.
    !
    ! CALLBACK_FCN
    !   Input, function to report progress and optionally request termination.
    !
    ! INFO
    !   Output, INTEGER(IK) scalar.
    !   INFO is the exit flag. It will be set to one of the following values defined in the module
    !   INFOS_MOD (see common/infos.f90):
    !   SMALL_TR_RADIUS: the lower bound for the trust region radius is reached;
    !   FTARGET_ACHIEVED: the target function value is reached;
    !   MAXFUN_REACHED: the objective function has been evaluated MAXFUN times;
    !   MAXTR_REACHED: the trust region iteration has been performed MAXTR times (MAXTR = 2*MAXFUN);
    !   NAN_INF_MODEL: NaN or Inf occurs in the model;
    !   NAN_INF_X: NaN or Inf occurs in X;
    !   DAMAGING_ROUNDING: the rounding error becomes damaging;
    !   NO_SPACE_BETWEEN_BOUNDS: there is not enough space between some lower and upper bounds, namely
    !   one of the difference XU(I)-XL(I) is less than 2*RHOBEG.
    !   !--------------------------------------------------------------------------!
    !   The following case(s) should NEVER occur unless there is a bug.
    !   NAN_INF_F: the objective function returns NaN or +Inf;
    !   TRSUBP_FAILED: a trust region step failed to reduce the model.
    !   !--------------------------------------------------------------------------!
    !--------------------------------------------------------------------------------------------------!
    
    ! Common modules
    use, non_intrinsic :: consts_mod, only : RP, IK, TWO, HALF, TEN, TENTH, EPS, BOUNDMAX, DEBUGGING
    use, non_intrinsic :: consts_mod, only : RHOBEG_DFT, RHOEND_DFT, FTARGET_DFT, MAXFUN_DIM_DFT, IPRINT_DFT
    use, non_intrinsic :: debug_mod, only : assert, warning
    use, non_intrinsic :: evaluate_mod, only : moderatex
    use, non_intrinsic :: history_mod, only : prehist
    use, non_intrinsic :: infnan_mod, only : is_nan, is_finite, is_posinf
    use, non_intrinsic :: infos_mod, only : NO_SPACE_BETWEEN_BOUNDS
    use, non_intrinsic :: linalg_mod, only : trueloc
    use, non_intrinsic :: memory_mod, only : safealloc
    use, non_intrinsic :: pintrf_mod, only : OBJ, CALLBACK
    use, non_intrinsic :: preproc_mod, only : preproc
    use, non_intrinsic :: string_mod, only : num2str
    
    ! Solver-specific modules
    use, non_intrinsic :: bobyqb_mod2, only : bobyqb
    
    implicit none
    
    ! Compulsory arguments
    procedure(OBJ) :: calfun  ! N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER
    real(RP), intent(inout) :: x(:)  ! X(N)
    
    ! Optional inputs
    procedure(CALLBACK), optional :: callback_fcn
    integer(IK), intent(in), optional :: iprint
    integer(IK), intent(in), optional :: maxfun
    integer(IK), intent(in), optional :: maxhist
    integer(IK), intent(in), optional :: npt
    logical, intent(in), optional :: honour_x0
    real(RP), intent(in), optional :: eta1
    real(RP), intent(in), optional :: eta2
    real(RP), intent(in), optional :: ftarget
    real(RP), intent(in), optional :: gamma1
    real(RP), intent(in), optional :: gamma2
    real(RP), intent(in), optional :: rhobeg
    real(RP), intent(in), optional :: rhoend
    real(RP), intent(in), optional :: xl(:)  ! XL(N)
    real(RP), intent(in), optional :: xu(:)  ! XU(N)
    
    ! Optional outputs
    integer(IK), intent(out), optional :: info
    integer(IK), intent(out), optional :: nf
    real(RP), intent(out), optional :: f
    real(RP), intent(out), optional, allocatable :: fhist(:)  ! FHIST(MAXFHIST)
    real(RP), intent(out), optional, allocatable :: xhist(:, :)  ! XHIST(N, MAXXHIST)
    
    ! Local variables
    character(len=*), parameter :: solver = 'BOBYQA'
    character(len=*), parameter :: srname = 'BOBYQA'
    integer(IK) :: info_loc
    integer(IK) :: iprint_loc
    integer(IK) :: k
    integer(IK) :: maxfun_loc
    integer(IK) :: maxhist_loc
    integer(IK) :: n
    integer(IK) :: i
    integer(IK) :: temp(size(x))
    integer(IK) :: nf_loc
    integer(IK) :: nhist
    integer(IK) :: npt_loc
    logical :: has_rhobeg
    logical :: honour_x0_loc
    real(RP) :: eta1_loc
    real(RP) :: eta2_loc
    real(RP) :: f_loc
    real(RP) :: ftarget_loc
    real(RP) :: gamma1_loc
    real(RP) :: gamma2_loc
    real(RP) :: rhobeg_loc
    real(RP) :: rhoend_loc
    real(RP) :: xl_loc(size(x))
    real(RP) :: xu_loc(size(x))
    real(RP), allocatable :: fhist_loc(:)  ! FHIST_LOC(MAXFHIST)
    real(RP), allocatable :: xhist_loc(:, :)  ! XHIST_LOC(N, MAXXHIST)
    
    ! Sizes
    n = int(size(x), kind(n))
    
    ! Preconditions
    if (DEBUGGING) then
        call assert(n >= 1, 'N >= 1', srname)
        if (present(xl)) then
            call assert(size(xl) == n .or. size(xl) == 0, 'SIZE(XL) == N unless XL is empty', srname)
        end if
        if (present(xu)) then
            call assert(size(xu) == n .or. size(xu) == 0, 'SIZE(XU) == N unless XU is empty', srname)
        end if
    end if
    
    ! Read the inputs
    
    xl_loc = -BOUNDMAX
    if (present(xl)) then
        if (size(xl) > 0) then
            xl_loc = xl
        end if
    end if
    
    temp = trueloc(is_nan(xl_loc) .or. xl_loc > BOUNDMAX)
    do i = 1,size(temp)
        xl_loc(temp(i)) = BOUNDMAX
    end do
    
    xu_loc = BOUNDMAX
    if (present(xu)) then
        if (size(xu) > 0) then
            xu_loc = xu
        end if
    end if
    temp = trueloc(is_nan(xu_loc) .or. xu_loc > BOUNDMAX)
    do i = 1,size(temp)
        xu_loc(temp(i)) = BOUNDMAX
    end do
    ! The solver requires that MINVAL(XU-XL) >= 2*RHOBEG, and we return if MINVAL(XU-XL) < 2*EPS.
    ! It would be better to fix the variables at (XU+XL)/2 wherever XU and XL almost equal, as is done
    ! in the MATLAB/Python interface of the solvers. In Fortran, this is doable using internal functions,
    ! but we choose not to implement it in the current version.
    if (any(xu_loc - xl_loc < TWO * EPS)) then
        if (present(info)) then
            info = NO_SPACE_BETWEEN_BOUNDS
        end if
        call warning(solver, 'There is no space between the lower and upper bounds of variable '// &
            & num2str(minval(trueloc(xu_loc - xl_loc < TWO * EPS)))//'. The solver cannot continue')
        return
    end if
    
    x = max(xl_loc, min(xu_loc, moderatex(x)))
    
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
    
    has_rhobeg = present(rhobeg)
    honour_x0_loc = .true.
    if (present(honour_x0)) then
        honour_x0_loc = honour_x0
    else if (has_rhobeg) then
        ! HONOUR_X0 is FALSE if user provides a valid RHOBEG. Is this the best choice?
        honour_x0_loc = (.not. (is_finite(rhobeg) .and. rhobeg > 0))
    end if
    
    
    ! Preprocess the inputs in case some of them are invalid. It does nothing if all inputs are valid.
    call preproc(solver, n, iprint_loc, maxfun_loc, maxhist_loc, ftarget_loc, rhobeg_loc, rhoend_loc, &
        & npt=npt_loc, eta1=eta1_loc, eta2=eta2_loc, gamma1=gamma1_loc, gamma2=gamma2_loc, &
        & has_rhobeg=has_rhobeg, honour_x0=honour_x0_loc, xl=xl_loc, xu=xu_loc, x0=x)
    
    ! Further revise MAXHIST_LOC according to MAXHISTMEM, and allocate memory for the history.
    ! In MATLAB/Python/Julia/R implementation, we should simply set MAXHIST = MAXFUN and initialize
    ! FHIST = NaN(1, MAXFUN), XHIST = NaN(N, MAXFUN)
    ! if they are requested; replace MAXFUN with 0 for the history that is not requested.
    call prehist(maxhist_loc, n, present(xhist), xhist_loc, present(fhist), fhist_loc)
    
    
    !-------------------- Call BOBYQB, which performs the real calculations. --------------------------!
    if (present(callback_fcn)) then
        call bobyqb(calfun, iprint_loc, maxfun_loc, npt_loc, eta1_loc, eta2_loc, ftarget_loc, &
            & gamma1_loc, gamma2_loc, rhobeg_loc, rhoend_loc, xl_loc, xu_loc, x, nf_loc, f_loc, &
            & fhist_loc, xhist_loc, info_loc, callback_fcn)
    else
        call bobyqb(calfun, iprint_loc, maxfun_loc, npt_loc, eta1_loc, eta2_loc, ftarget_loc, &
            & gamma1_loc, gamma2_loc, rhobeg_loc, rhoend_loc, xl_loc, xu_loc, x, nf_loc, f_loc, &
            & fhist_loc, xhist_loc, info_loc)
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
    
    ! If NF_LOC > MAXHIST_LOC, warn that not all history is recorded.
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
    
        if (present(xl)) then
            if (size(xl) == size(x)) then
                call assert(all(x >= xl), 'X >= XL', srname)
                if (present(xhist)) then
                    do k = 1, nhist
                        call assert(all(xhist(:, k) >= xl), 'XHIST >= XL', srname)
                    end do
                end if
            end if
        end if
    
        if (present(xu)) then
            if (size(xu) == size(x)) then
                call assert(all(x <= xu), 'X <= XU', srname)
                if (present(xhist)) then
                    do k = 1, nhist
                        call assert(all(xhist(:, k) <= xu), 'XHIST <= XU', srname)
                    end do
                end if
            end if
        end if
    
        if (present(fhist)) then
            call assert(size(fhist) == nhist, 'SIZE(FHIST) == NHIST', srname)
            call assert(.not. any(is_nan(fhist) .or. is_posinf(fhist)), 'FHIST does not contain NaN/+Inf', srname)
            call assert(.not. any(fhist < f_loc), 'F is the smallest in FHIST', srname)
        end if
    end if
    
    end subroutine bobyqa
    
    
    end module bobyqa_mod2
    
    
    
    
    !---------------------------------------- THE MAIN PROGRAM ----------------------------------------!
    program bobyqa_exmp
    
    !  The following line makes the solver available.
     use bobyqa_mod2, only : bobyqa
     ! The following line specifies which module provides CALFUN and CALLBACK_FCN.
     use calfun_mod2, only : RP, IK, calfun, callback_fcn
     implicit none
     integer, parameter :: n = 2
     integer :: nf, info
     real(RP) :: f, x(n), x0(n), lb(n), ub(n)
     ! Define the starting point.
     x0 = 0.0_RP
     ! Define the lower and upper bounds. We define an upper bound that will be active
     ! in order to demonstrate the usage of bounds.
     lb = -1.0_RP
     ub = 4.5_RP
     ! The following lines illustrates how to call the solver.
     x = x0
     call bobyqa(calfun, x, f, lb, ub)  ! This call will not print anything.
     ! In addition to the compulsory arguments, the following illustration specifies also RHOBEG and
     ! IPRINT, which are optional. All the unspecified optional arguments (RHOEND, MAXFUN, etc.) will
     ! take their default values coded in the solver.
    !  x = x0
    ! call bobyqa(calfun, x, f, lb, ub, rhobeg=1.0_RP, iprint=1_IK, nf=nf, info=info, callback_fcn=callback_fcn)
    
    end program bobyqa_exmp
    