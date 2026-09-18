function zc = build_zero_constraint(ZeroRestrictions, cache, sign_regime_dependent, ny, validate, zero_on_average)
%BUILD_ZERO_CONSTRAINT  Zero IRF restrictions -> per-shock linear constraints
%   on the columns of the drawn rotation P (= Q'), for the Arias, Rubio-Ramirez
%   & Waggoner (2018) null-space draw + volume-element reweighting.
%
%   zc = build_zero_constraint(ZeroRestrictions, cache, sign_regime_dependent, ny)
%   zc = build_zero_constraint(..., validate)
%   zc = build_zero_constraint(..., validate, zero_on_average)
%
% ============================================================================
% ZERO ON THE AVERAGE (REPRESENTATIVE) IRF  -- zero_on_average (default true)
%   The stored "average" IRF (output.sign_irfs_temp, built by computeOutputs in
%   SignRestrictionCheck) uses the REGIME-1 structure with lambda = 1:
%       IRF_avg = Psi_1 * A0inv_1 * Q'        (no diag(sqrt lambda_r) scaling).
%   With sign_regime_dependent = true the per-regime transmission instead
%   carries diag(sqrt lambda_r), so a zero imposed on a single regime's
%   transmission is EXACTLY zero in that regime's IRF but is in general
%   NONZERO in the average IRF (the lambda scaling does not cancel).
%
%   zero_on_average = true (DEFAULT) builds the constraint against the SAME
%   transmission that produces the stored average IRF, i.e. lambda is dropped
%   (Mmat = Psi_H * A0inv) regardless of sign_regime_dependent. The zero then
%   holds EXACTLY in output.sign_irfs_temp, and because lambda is gone the
%   constraint is the SAME for every regime -> one restriction suffices (no
%   need to repeat it per regime), which also keeps the admissible set wide.
%   The SIGN restrictions are unaffected and are still screened in EVERY regime
%   by the caller (checkSignRestrictions / the impact-sign screen in
%   draw_Q_zero_columnwise), so signs stay regime-robust while the zero is an
%   average-IRF (representative) restriction.
%
%   zero_on_average = false reproduces the old per-regime behaviour: the
%   constraint uses regime-r transmission (lambda included when
%   sign_regime_dependent = true), so the zero holds in that regime's IRF but
%   not necessarily in the stored average.
%
% ============================================================================
% REGIME INDEX CONVENTION  (corrected)
%   build_theta_cache reports TWO regime counts:
%     cache.nRegimes : number of regime-specific A0 blocks (1 when tvA is OFF,
%                      since A0 is common; > 1 only under tvA).
%     cache.nlambda  : number of lambda (variance-ratio) columns = the number
%                      of heteroskedasticity regimes, present even when
%                      tvA is OFF.
%   A zero restriction's regime_idx selects which regime's TRANSMISSION carries
%   the zero. The transmission in regime r is  Psi_r * A0inv_r * diag(sqrt
%   lambda_r) ; A0inv_r is regime-specific only under tvA (else the common A0),
%   while lambda_r is always the r-th lambda column. So regime_idx ranges over
%   the LAMBDA regimes (1..cache.nlambda), NOT 1..cache.nRegimes. A0 slicing
%   below uses min(r, nRegimes-effectively) via local_slice*, which return the
%   common A0 when tvA is OFF, so regime_idx > 1 on a non-tvA run is valid and
%   simply selects that lambda column (the A0 piece is common).
%
%   For s.r.d. = 0 the lambda scalar drops out of the constraint entirely, so
%   on a non-tvA run EVERY regime_idx gives the SAME constraint row (A0inv and
%   Psi are common); the choice is then immaterial and any 1..nlambda is fine.
% ============================================================================
% WHAT THIS DOES
%   A zero restriction "shock j0 has no effect on variable i0 at horizon H in
%   regime r*" is a LINEAR restriction on the j0-th column p_{j0} of P (= Q'):
%
%       IRF_H^{(r*)}(i0,j0) = c' * p_{j0} = 0,
%       c = ( Psi_{r*}^{(H)} * L_{r*} )' * e_{i0}         (s.r.d. = 1)
%       c = ( Psi_{r*}^{(H)} * A0^{-1}_{r*} )' * e_{i0}   (s.r.d. = 0)
%
%   Psi^{(0)} = I; H is the ECONOMIC horizon (0 = impact), at cache.Psi(:,:,H+1).
%   Rows c' are stacked per shock j into Z{j}; the column draw builds p_j in the
%   null space of M_j = [ p_{built before j} ; Z{j} ].
%
% ============================================================================
% INPUT SCHEMA  (options.ZeroRestrictions); see header of the original.
% OUTPUT  zc (struct); see header of the original.
% ============================================================================

    if nargin < 5 || isempty(validate), validate = true; end
    if nargin < 6 || isempty(zero_on_average), zero_on_average = true; end

    % ---- Empty / no-zero fast path ----------------------------------------
    if isempty(ZeroRestrictions) || numel(ZeroRestrictions) == 0
        zc = local_empty(ny);
        return
    end

    % ---- Cache geometry ----------------------------------------------------
    tvA  = isfield(cache, 'tvA') && cache.tvA;
    if ~isfield(cache, 'nRegimes') || isempty(cache.nRegimes)
        error('build_zero_constraint:noCache', 'cache.nRegimes missing.');
    end
    nReg = cache.nRegimes;
    nlam = cache.nlambda;
    if ~isfield(cache, 'A0inv') || isempty(cache.A0inv)
        error('build_zero_constraint:noCache', ...
              'cache.A0inv is empty; pass a cache built by build_theta_cache.');
    end

    % Effective regime ceiling for regime_idx: the number of TRANSMISSION
    % regimes available, which is the lambda-column count (>= the A0-block
    % count). Under tvA nReg == nlam; under non-tvA nReg == 1 < nlam.
    nReg_eff = max(nReg, nlam);

    nR = numel(ZeroRestrictions);

    nz   = zeros(1, ny);
    Z    = cell(1, ny);
    meta = repmat(struct('variable_idx', [], 'shock_idx', [], ...
                         'horizon', [], 'regime_idx', []), nR, 1);

    % ---- Resolve / validate each restriction; evaluate c -------------------
    for k = 1:nR
        rk = ZeroRestrictions(k);

        i0 = local_getfield(rk, 'variable_idx', []);
        j0 = local_getfield(rk, 'shock_idx',    []);
        H  = local_getfield(rk, 'horizon',      0);
        if nReg_eff > 1
            r = local_getfield(rk, 'regime_idx', []);
        else
            r = local_getfield(rk, 'regime_idx', 1);
        end

        if isempty(i0) || isempty(j0)
            error('build_zero_constraint:missingField', ...
                  'ZeroRestrictions(%d): variable_idx and shock_idx are required.', k);
        end
        local_checkint(i0, 1, ny,  k, 'variable_idx');
        local_checkint(j0, 1, ny,  k, 'shock_idx');
        local_checkint(H,  0, Inf, k, 'horizon');
        if isempty(r)
            error('build_zero_constraint:missingRegime', ...
                  ['ZeroRestrictions(%d): regime_idx is required when there is ' ...
                   'more than one regime.'], k);
        end
        % regime_idx is validated against the TRANSMISSION-regime count
        % (lambda columns), NOT cache.nRegimes (the A0-block count).
        local_checkint(r, 1, nReg_eff, k, 'regime_idx');
        if r > nlam
            error('build_zero_constraint:regimeLambda', ...
                  ['ZeroRestrictions(%d): regime_idx = %d exceeds the number of ' ...
                   'lambda columns (%d).'], k, r, nlam);
        end

        % --- regime transmission pieces ---
        % zero_on_average = true: match the STORED average IRF, which uses the
        %   regime-1 A0inv and lambda = 1 (computeOutputs: smat = A0_1 \ Q').
        %   So drop lambda and use regime-1 A0inv regardless of sign_regime_
        %   dependent and regardless of the restriction's regime_idx.
        % zero_on_average = false: per-regime transmission (regime r,
        %   lambda included when sign_regime_dependent = true).
        if zero_on_average
            A0inv_r   = local_slice3(cache.A0inv, 1, tvA);       % regime-1 A0inv
            sqrtlam_r = ones(ny, 1);                             % lambda = 1 (average)
            r_eff     = 1;                                       % horizon slice uses regime 1
        else
            A0inv_r   = local_slice3(cache.A0inv, r, tvA);
            sqrtlam_r = sqrt(max(cache.lambda(:, r), 0));        % ny x 1 (r-th lambda col)
            r_eff     = r;
        end

        % --- H-step propagator (economic horizon H -> Psi index H+1) ---
        if H == 0
            PsiH = eye(ny);
        else
            Psi_r = local_slice4(cache.Psi, r_eff, tvA);          % ny x ny x horizon
            if isempty(Psi_r), navail = 0; else, navail = size(Psi_r, 3); end
            if navail < H + 1
                error('build_zero_constraint:horizonRange', ...
                      ['ZeroRestrictions(%d): horizon H = %d needs cached ' ...
                       'propagator slice Psi(:,:,%d), but only %d slice(s) are ' ...
                       'stored. Raise the cache horizon (sr_max_horizon) so it ' ...
                       'covers the zero-restriction horizons.'], k, H, H+1, navail);
            end
            PsiH = Psi_r(:, :, H + 1);
        end

        % --- constraint row c' (see header) ---
        % Under zero_on_average the sqrtlam_r = 1 vector makes both branches
        % identical (lambda drops out); the average IRF is exactly Psi_H*A0inv*q.
        if zero_on_average || sign_regime_dependent == 0
            Mmat = PsiH * A0inv_r;                            % Lambda drops out
        else
            Lr   = A0inv_r .* sqrtlam_r.';                    % A0inv_r * diag(sqrtlam_r)
            Mmat = PsiH * Lr;
        end
        c = Mmat(i0, :).';                                    % ny x 1 ; c' p_{j0} = IRF_H(i0,j0)

        nc = norm(c);
        if ~(nc > 0) || ~all(isfinite(c))
            error('build_zero_constraint:degenerate', ...
                  ['ZeroRestrictions(%d): constraint vector is zero or non-finite ' ...
                   '(degenerate restriction at this draw).'], k);
        end
        c = c / nc;                                           % unit-norm (null space unchanged)

        Z{j0}  = [Z{j0}; c.'];                                % stack as a row
        nz(j0) = nz(j0) + 1;

        meta(k).variable_idx = i0;
        meta(k).shock_idx    = j0;
        meta(k).horizon      = H;
        meta(k).regime_idx   = r;
    end

    % ---- Collapse duplicate constraint rows per shock ---------------------
    % Under zero_on_average the same (variable, horizon) restriction imposed on
    % several regimes produces IDENTICAL rows (lambda dropped, A0inv common),
    % so a restriction repeated across all regimes collapses to ONE row.
    % Remove duplicate (and sign-flipped duplicate) rows so Z{j} is full row
    % rank and the null-space construction in draw_Q_zero_columnwise gets the
    % true constraint count. Rows are already unit-norm.
    for j = find(nz > 0)
        Zj = Z{j};
        keep = true(size(Zj, 1), 1);
        for a = 2:size(Zj, 1)
            for b = 1:a-1
                if keep(b)
                    d = norm(Zj(a, :) - Zj(b, :));
                    s = norm(Zj(a, :) + Zj(b, :));   % sign-flipped duplicate
                    if min(d, s) < 1e-9
                        keep(a) = false; break
                    end
                end
            end
        end
        Z{j}  = Zj(keep, :);
        nz(j) = sum(keep);
    end

    % ---- Optional: full-row-rank check on each Z{j} (ARRW requirement) -----
    if validate
        for j = find(nz > 0)
            rk_j = rank(Z{j}, 1e-9);
            if rk_j < nz(j)
                error('build_zero_constraint:rankDeficient', ...
                      ['Shock %d carries %d zero restrictions whose constraint ' ...
                       'vectors are linearly dependent (rank %d). ARRW require ' ...
                       'Z_j full row rank; drop the redundant restriction(s).'], ...
                      j, nz(j), rk_j);
            end
        end
    end

    % ---- Assemble output ---------------------------------------------------
    zc.has_zero        = true;
    zc.nz              = nz;
    zc.restricted_cols = find(nz > 0);
    zc.Z               = Z;
    zc.single_column   = (isscalar(zc.restricted_cols));
    zc.meta            = meta;

    % build order: most-restricted shock first (MATLAB sort is stable, so ties
    % keep ascending index order); unrestricted shocks (nz = 0) land last.
    [~, zc.col_order]  = sort(nz, 'descend');

    % feasibility: dim null(M_k) = ny - (k-1) - nz(col_order(k)) >= 1
    %           => nz(col_order(k)) <= ny - k.
    nz_ord = nz(zc.col_order);
    zc.feasible = all(nz_ord <= (ny - (1:ny)));
    if validate && ~zc.feasible
        kbad = find(nz_ord > (ny - (1:ny)), 1);
        warning('build_zero_constraint:infeasible', ...
                ['Zero pattern over-restricts: at build step %d the column carries ' ...
                 '%d zeros but only %d free direction(s) remain (need <= %d). ' ...
                 'No rotation can satisfy all zeros.'], ...
                kbad, nz_ord(kbad), ny - kbad, ny - kbad);
    end
end


% =====================================================================
% LOCAL HELPERS
% =====================================================================

function zc = local_empty(ny)
    zc.has_zero        = false;
    zc.nz              = zeros(1, ny);
    zc.restricted_cols = [];
    zc.Z               = cell(1, ny);
    zc.single_column   = false;
    zc.col_order       = 1:ny;
    zc.feasible        = true;
    zc.meta            = repmat(struct('variable_idx', [], 'shock_idx', [], ...
                                       'horizon', [], 'regime_idx', []), 0, 1);
end


function v = local_getfield(s, name, default)
    if isfield(s, name) && ~isempty(s.(name))
        v = s.(name);
    else
        v = default;
    end
end


function local_checkint(x, lo, hi, k, name)
    ok = isscalar(x) && isfinite(x) && (x == floor(x)) && (x >= lo) && (x <= hi);
    if ~ok
        if isinf(hi), rng = sprintf('>= %d', lo); else, rng = sprintf('in %d..%d', lo, hi); end
        error('build_zero_constraint:badField', ...
              'ZeroRestrictions(%d): %s must be an integer %s.', k, name, rng);
    end
end


function M = local_slice3(A, r, tvA)
% A0inv: ny x ny (non-tvA) or ny x ny x nRegimes (tvA)
    if tvA, M = A(:, :, r); else, M = A; end
end


function M = local_slice4(A, r, tvA)
% By / Psi: 3-D (non-tvA) or 4-D (tvA); empty stays empty.
    if isempty(A)
        M = A;
    elseif tvA
        M = A(:, :, :, r);
    else
        M = A;
    end
end