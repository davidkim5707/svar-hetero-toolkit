function test_build_zero_constraint()
%TEST_BUILD_ZERO_CONSTRAINT  Standalone self-test for build_zero_constraint.m.
%
%   Verifies the core identity  c' * p_{j0} = IRF_H^{(r*)}(i0,j0)  by building a
%   synthetic cache, drawing p_{j0} in null(c'), completing it to an orthogonal
%   P, and confirming the targeted IRF element is machine-zero (while a generic
%   unconstrained column is not). Covers both sign_regime_dependent conventions
%   and several horizons/regimes. Also checks the metadata fields (single_column,
%   col_order, feasible, rank guard).
%
%   Run:  test_build_zero_constraint
%   Requires build_zero_constraint.m on the path.

    rng(0);
    ny = 7; nReg = 4; horizon = 6;
    tol = 1e-10;

    % --- synthetic structural pieces (mirror a build_theta_cache cache) ---
    A0 = randn(ny);
    A0(1:ny+1:end) = abs(diag(A0)).' + 1;          % keep invertible / well-cond
    A0inv = A0 \ eye(ny);
    lam = 0.5 + rand(ny, nReg);                    % lambda(:,r)

    By = 0.2 * randn(ny);                          % toy 1-lag companion (stable-ish)
    Psi = zeros(ny, ny, horizon);
    Psi(:, :, 1) = eye(ny);                        % impact
    for h = 2:horizon, Psi(:, :, h) = By * Psi(:, :, h-1); end

    cache = struct('tvA', false, 'nRegimes', nReg, 'nlambda', nReg, ...
                   'A0inv', A0inv, 'Psi', Psi, 'lambda', lam);

    fprintf('\n== test_build_zero_constraint ==\n');
    fprintf('core identity: p_{j0} in null(c'')  =>  IRF_H(i0,j0) == 0\n');

    cases = [ 2 1 0 4 ;     % i0 j0 H r
              4 3 0 4 ;
              2 1 2 1 ;
              6 1 3 2 ];
    worst = 0;
    for srd = [0 1]
        for ci = 1:size(cases, 1)
            i0 = cases(ci,1); j0 = cases(ci,2); H = cases(ci,3); r = cases(ci,4);

            zr = struct('variable_idx', i0, 'shock_idx', j0, ...
                        'horizon', H, 'regime_idx', r);
            zc = build_zero_constraint(zr, cache, srd, ny, [], false);   % per-regime zeros, the case these checks evaluate
            assert(zc.has_zero && zc.single_column, 'single_column flag wrong');
            assert(isequal(zc.restricted_cols, j0), 'restricted_cols wrong');
            assert(zc.col_order(1) == j0, 'col_order must build restricted col first');

            c = zc.Z{j0}.';                        % ny x 1 (single zero)

            % p_{j0} in null(c'): random, project out c, normalize
            p = randn(ny, 1); p = p - c * (c' * p); p = p / norm(p);

            % orthonormal P with column j0 = (+/-) p
            M = randn(ny); M(:, 1) = p;
            [Qo, ~] = qr(M);                       % Qo(:,1) = +/- p
            P = Qo; P(:, [1 j0]) = P(:, [j0 1]);   % move constrained col to position j0

            val      = irf_element(P, A0inv, lam, Psi, i0, j0, H, r, srd);
            Pc       = qr_rand(ny);                % generic unconstrained rotation
            val_ctrl = irf_element(Pc, A0inv, lam, Psi, i0, j0, H, r, srd);

            worst = max(worst, abs(val));
            flag  = '';
            if abs(val) > tol, flag = '   <-- FAIL'; end
            fprintf(['  srd=%d (i0=%d,j0=%d,H=%d,r=%d): |IRF|=%.2e   ' ...
                     'control|IRF|=%.2e%s\n'], ...
                     srd, i0, j0, H, r, abs(val), abs(val_ctrl), flag);
            assert(abs(val) < tol, 'targeted IRF element is not zero');
        end
    end
    fprintf('  worst |IRF| over all cases: %.2e  (tol %.1e)\n', worst, tol);

    % --- metadata: two zeros on two DIFFERENT shocks ---
    zr2 = struct('variable_idx', {2, 7}, 'shock_idx', {1, 3}, ...
                 'horizon', {0, 0}, 'regime_idx', {4, 4});
    zc2 = build_zero_constraint(zr2, cache, 1, ny);
    assert(~zc2.single_column, 'two-column case must not be single_column');
    assert(isequal(zc2.restricted_cols, [1 3]), 'restricted_cols (two-col) wrong');
    assert(zc2.feasible, 'two single zeros must be feasible');
    fprintf('  two-shock case: single_column=%d feasible=%d restricted=[%s]\n', ...
            zc2.single_column, zc2.feasible, num2str(zc2.restricted_cols));

    % --- empty case ---
    z0 = build_zero_constraint([], cache, 1, ny);
    assert(~z0.has_zero && isequal(z0.col_order, 1:ny), 'empty case wrong');
    fprintf('  empty case: has_zero=%d\n', z0.has_zero);

    % --- duplicate rows: two identical restrictions on one shock collapse to one ---
    zdup = struct('variable_idx', {2, 2}, 'shock_idx', {1, 1}, ...
                  'horizon', {0, 0}, 'regime_idx', {4, 4});
    zd = build_zero_constraint(zdup, cache, 1, ny);
    assert(zd.nz(1) == 1 && size(zd.Z{1}, 1) == 1, 'duplicate restriction should collapse to one row');
    fprintf('  duplicate restriction: collapsed to nz = %d\n', zd.nz(1));

    fprintf('ALL TESTS PASSED\n\n');
end


function val = irf_element(P, A0inv, lam, Psi, i0, j0, H, r, srd)
% IRF_H(i0,j0) under the code's convention:
%   s.r.d.=1: smat = (A0inv*Dh)*P ;  s.r.d.=0: smat = (A0inv*P)*Dh
    Dh = diag(sqrt(max(lam(:, r), 0)));
    if srd == 0
        smat = (A0inv * P) * Dh;
    else
        smat = (A0inv * Dh) * P;
    end
    if H == 0, PsiH = eye(size(P, 1)); else, PsiH = Psi(:, :, H + 1); end
    IRF = PsiH * smat;
    val = IRF(i0, j0);
end


function Q = qr_rand(n)
    [Q, ~] = qr(randn(n));
end