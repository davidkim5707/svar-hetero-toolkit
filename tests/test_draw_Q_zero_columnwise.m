function test_draw_Q_zero_columnwise()
%TEST_DRAW_Q_ZERO_COLUMNWISE  Self-test for draw_Q_zero_columnwise.m (zero path).
%
%   Builds a synthetic cache, forms zero constraints with build_zero_constraint,
%   and checks that the drawn rotation Q (a) is orthonormal, (b) satisfies the
%   zeros EXACTLY (IRF elements ~ 0), (c) satisfies the impact-sign screen, and
%   (d) is reproduced by gaussian_to_Q(X). Covers single- and two-shock zeros.
%
%   Run:  test_draw_Q_zero_columnwise
%   Requires: build_zero_constraint.m, gaussian_to_Q_zero.m  (zero path).
%             find_admissible_Q_columnwise.m is only needed for the no-zero
%             delegation, which this test does not exercise.

    rng(7);
    ny = 7; nReg = 4; horizon = 6;
    tol = 1e-10;

    A0 = randn(ny);
    A0(1:ny+1:end) = abs(diag(A0)).' + 1;
    A0inv = A0 \ eye(ny);
    lam = 0.5 + rand(ny, nReg);

    By = 0.2 * randn(ny);
    Psi = zeros(ny, ny, horizon); Psi(:, :, 1) = eye(ny);
    for h = 2:horizon, Psi(:, :, h) = By * Psi(:, :, h-1); end

    cache = struct('tvA', false, 'nRegimes', nReg, 'nlambda', nReg, ...
                   'A0inv', A0inv, 'Psi', Psi, 'lambda', lam);

    srd = 1;                                   % sign_regime_dependent
    check_full_fn = @(Q) true;                 % isolate the construction

    fprintf('\n== test_draw_Q_zero_columnwise ==\n');

    % ---------- Case A: single zero on shock 1 + impact signs on shock 1 ----
    zrA = struct('variable_idx', 2, 'shock_idx', 1, 'horizon', 0, 'regime_idx', 4);
    zcA = build_zero_constraint(zrA, cache, srd, ny, [], false);   % per-regime zeros, the case these checks evaluate

    impact_signs = zeros(ny);
    impact_signs(1, 1) = +1;   % FFR up to MP
    impact_signs(5, 1) = -1;   % GDP down to MP

    [X, Q, ok, no] = draw_Q_zero_columnwise(ny, check_full_fn, A0inv, lam, ...
                        impact_signs, srd, zcA, 200, 10000);
    assert(ok, 'Case A: no admissible Q found');
    report('A (single zero)', Q, X, A0inv, lam, Psi, impact_signs, srd, ...
           zcA, no, tol);
    assert(zcA.single_column, 'Case A should be single_column');

    % ---------- Case B: two zeros on two DIFFERENT shocks -------------------
    zrB = struct('variable_idx', {2, 7}, 'shock_idx', {1, 3}, ...
                 'horizon', {0, 2}, 'regime_idx', {4, 1});   % one impact, one h=2
    zcB = build_zero_constraint(zrB, cache, srd, ny, [], false);   % per-regime zeros, the case these checks evaluate

    impact_signs_B = zeros(ny);
    impact_signs_B(1, 1) = +1;          % MP
    impact_signs_B(2, 3) = -1;          % SOMA down to QT2

    [XB, QB, okB, noB] = draw_Q_zero_columnwise(ny, check_full_fn, A0inv, lam, ...
                        impact_signs_B, srd, zcB, 200, 10000);
    assert(okB, 'Case B: no admissible Q found');
    report('B (two-shock zeros)', QB, XB, A0inv, lam, Psi, impact_signs_B, srd, ...
           zcB, noB, tol);
    assert(~zcB.single_column, 'Case B should NOT be single_column');

    fprintf('ALL TESTS PASSED\n\n');
end


function report(tag, Q, X, A0inv, lam, Psi, impact_signs, srd, zc, n_outer, tol)
    ny = size(Q, 1);

    % (a) orthonormality
    orth = max(max(abs(Q*Q.' - eye(ny))));

    % (b) zeros: IRF element for each restriction in zc.meta
    zmax = 0;
    for k = 1:numel(zc.meta)
        m = zc.meta(k);
        v = irf_element(Q.', A0inv, lam, Psi, m.variable_idx, m.shock_idx, ...
                        m.horizon, m.regime_idx, srd);
        zmax = max(zmax, abs(v));
    end

    % (c) impact-sign screen (regime 1 is enough to display; screen used all)
    sgn_ok = true;
    for j = 1:ny
        if any(impact_signs(:, j) ~= 0)
            for r = 1:size(lam, 2)
                Dh = diag(sqrt(max(lam(:, r), 0)));
                smat = (A0inv * Dh) * Q.';          % srd == 1 convention
                if srd == 0, smat = (A0inv * Q.') * Dh; end
                s = smat(:, j);
                idx = impact_signs(:, j) ~= 0;
                if any(impact_signs(idx, j) .* s(idx) <= 0), sgn_ok = false; end
            end
        end
    end

    % (d) X reproduces Q via gaussian_to_Q_zero (the map the ESS uses)
    Qx = gaussian_to_Q_zero(X, zc);
    xrep = max(max(abs(Qx - Q)));

    fprintf(['  %-22s outer=%d  |QQ''-I|=%.1e  max|zero IRF|=%.1e  ' ...
             'signs_ok=%d  |g2Q(X)-Q|=%.1e\n'], ...
             tag, n_outer, orth, zmax, sgn_ok, xrep);
    assert(orth < tol,  [tag ': not orthonormal']);
    assert(zmax < tol,  [tag ': zero restriction violated']);
    assert(sgn_ok,      [tag ': impact-sign screen violated']);
    assert(xrep < 1e-8, [tag ': gaussian_to_Q(X) != Q']);
end


function val = irf_element(P, A0inv, lam, Psi, i0, j0, H, r, srd)
% P has columns = drawn vectors (P = Q'); convention smat=(A0inv*Dh)*P (srd=1).
    Dh = diag(sqrt(max(lam(:, r), 0)));
    if srd == 0, smat = (A0inv * P) * Dh; else, smat = (A0inv * Dh) * P; end
    if H == 0, PsiH = eye(size(P, 1)); else, PsiH = Psi(:, :, H + 1); end
    IRF = PsiH * smat;
    val = IRF(i0, j0);
end