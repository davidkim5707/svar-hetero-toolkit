function test_ess_draw_Q_zero_columnwise()
%TEST_ESS_DRAW_Q_ZERO_COLUMNWISE  Self-test for the zero-confined ESS.
%
%   Pipeline: draw_Q_zero_columnwise (init) -> ess_draw_Q_zero_columnwise (many
%   sweeps). Checks, for single- and two-shock zero patterns:
%     (i)   seed consistency:  gaussian_to_Q_zero(X_init, zc) == Q_init
%     (ii)  every swept Q stays orthonormal              (~1e-15)
%     (iii) every swept Q satisfies the zeros EXACTLY    (~1e-15)
%     (iv)  every accepted Q satisfies the impact-sign screen
%     (v)   the chain actually MOVES (explores the manifold)
%
%   Run:  test_ess_draw_Q_zero_columnwise
%   Requires: build_zero_constraint.m, draw_Q_zero_columnwise.m,
%             gaussian_to_Q_zero.m, ess_draw_Q_zero_columnwise.m.
%             (Zero path only; the no-zero delegation is not exercised.)

    rng(11);
    ny = 7; nReg = 4; horizon = 6;
    tol = 1e-10;

    A0 = randn(ny); A0(1:ny+1:end) = abs(diag(A0)).' + 1;
    A0inv = A0 \ eye(ny);
    lam = 0.5 + rand(ny, nReg);
    By = 0.2 * randn(ny);
    Psi = zeros(ny, ny, horizon); Psi(:, :, 1) = eye(ny);
    for h = 2:horizon, Psi(:, :, h) = By * Psi(:, :, h-1); end
    cache = struct('tvA', false, 'nRegimes', nReg, 'nlambda', nReg, ...
                   'A0inv', A0inv, 'Psi', Psi, 'lambda', lam);

    srd        = 1;
    n_sweeps   = 40;
    max_shrink = 30;

    fprintf('\n== test_ess_draw_Q_zero_columnwise ==\n');

    % ---------- Case A: single zero on shock 1 -----------------------------
    zrA = struct('variable_idx', 2, 'shock_idx', 1, 'horizon', 0, 'regime_idx', 4);
    zcA = build_zero_constraint(zrA, cache, srd, ny, [], false);   % per-regime zeros, the case these checks evaluate
    sgnA = zeros(ny); sgnA(1, 1) = +1; sgnA(5, 1) = -1;
    run_case('A (single zero)', cache, zcA, sgnA, srd, n_sweeps, max_shrink, tol);

    % ---------- Case B: two zeros on two different shocks ------------------
    zrB = struct('variable_idx', {2, 7}, 'shock_idx', {1, 3}, ...
                 'horizon', {0, 2}, 'regime_idx', {4, 1});
    zcB = build_zero_constraint(zrB, cache, srd, ny, [], false);   % per-regime zeros, the case these checks evaluate
    sgnB = zeros(ny); sgnB(1, 1) = +1; sgnB(2, 3) = -1;
    run_case('B (two-shock zeros)', cache, zcB, sgnB, srd, n_sweeps, max_shrink, tol);

    fprintf('ALL TESTS PASSED\n\n');
end


function run_case(tag, cache, zc, impact_signs, srd, n_sweeps, max_shrink, tol)
    ny    = size(cache.A0inv, 1);
    A0inv = cache.A0inv; lam = cache.lambda; Psi = cache.Psi;

    check_fn = @(Q) sign_screen(Q, A0inv, lam, impact_signs, srd);

    % --- init ---
    [X, Q, ok, ~] = draw_Q_zero_columnwise(ny, check_fn, A0inv, lam, ...
                        impact_signs, srd, zc, 200, 10000);
    assert(ok, [tag ': init failed']);

    % (i) seed consistency
    seed_err = max(max(abs(gaussian_to_Q_zero(X, zc) - Q)));
    assert(seed_err < tol, [tag ': seed does not reproduce Q']);

    Q_init = Q;
    worst_orth = max(max(abs(Q*Q.' - eye(ny))));
    worst_zero = zero_residual(Q, zc, A0inv, lam, Psi, srd);
    worst_sgn  = 0;

    % --- sweeps ---
    for sweep = 1:n_sweeps
        [X, Q, ~] = ess_draw_Q_zero_columnwise(X, check_fn, zc, max_shrink);

        worst_orth = max(worst_orth, max(max(abs(Q*Q.' - eye(ny)))));
        worst_zero = max(worst_zero, zero_residual(Q, zc, A0inv, lam, Psi, srd));
        if ~sign_screen(Q, A0inv, lam, impact_signs, srd), worst_sgn = 1; end
    end

    moved = norm(Q - Q_init, 'fro');

    fprintf(['  %-22s sweeps=%d  seed_err=%.1e  worst|QQ''-I|=%.1e  ' ...
             'worst|zero|=%.1e  sign_viol=%d  moved=%.3f\n'], ...
             tag, n_sweeps, seed_err, worst_orth, worst_zero, worst_sgn, moved);

    assert(worst_orth < tol, [tag ': orthonormality lost during sweeps']);
    assert(worst_zero < tol, [tag ': zero restriction violated during sweeps']);
    assert(worst_sgn == 0,   [tag ': accepted a sign-violating Q']);
    assert(moved > 1e-2,     [tag ': chain did not move (no exploration)']);
end


function ok = sign_screen(Q, A0inv, lam, impact_signs, srd)
% Strict impact-sign screen across all regimes (no sign flipping), matching the
% per-column screen used by draw_Q_zero_columnwise.
    ny = size(Q, 1); nReg = size(lam, 2); ok = true;
    P = Q.';
    for j = 1:ny
        if ~any(impact_signs(:, j) ~= 0), continue; end
        for r = 1:nReg
            Dh = diag(sqrt(max(lam(:, r), 0)));
            if srd == 0, smat = (A0inv * P) * Dh; else, smat = (A0inv * Dh) * P; end
            s   = smat(:, j);
            idx = impact_signs(:, j) ~= 0;
            if any(impact_signs(idx, j) .* s(idx) <= 0), ok = false; return; end
        end
    end
end


function z = zero_residual(Q, zc, A0inv, lam, Psi, srd)
% Max |IRF element| over all zero restrictions in zc.meta.
    z = 0; P = Q.';
    for k = 1:numel(zc.meta)
        m  = zc.meta(k);
        Dh = diag(sqrt(max(lam(:, m.regime_idx), 0)));
        if srd == 0, smat = (A0inv * P) * Dh; else, smat = (A0inv * Dh) * P; end
        if m.horizon == 0, PsiH = eye(size(P, 1)); else, PsiH = Psi(:, :, m.horizon + 1); end
        IRF = PsiH * smat;
        z = max(z, abs(IRF(m.variable_idx, m.shock_idx)));
    end
end