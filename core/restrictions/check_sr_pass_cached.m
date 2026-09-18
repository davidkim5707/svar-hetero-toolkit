function [pass, omega_hat, om_diag] = check_sr_pass_cached(Q, cache, var, ...
    SignRestrictions, NarrativeRestrictions, horizon, ...
    sign_regime_dependent, penalty_offdiagonal_on, ny, nrr_omega_opts)
%CHECK_SR_PASS_CACHED  Fast Q-only sign-restriction check (tvA-aware).
%
%    Optional 10th input nrr_omega_opts and optional
%   outputs omega_hat / om_diag. When the caller requests omega_hat
%   (nargout >= 2) AND supplies nrr_omega_opts AND the draw PASSES all
%   restrictions, the AD-RR narrative weight omega(theta, Q) is estimated by
%   estimate_nrr_omega at the exact objects (fsign, fsign_cell,
%   irTemp_narrative_cell) this function passes to
%   check_narrative_regime_avg -- so the simulated event is byte-identical
%   to the retained-draw check. Hot-loop calls with one output are
%   UNCHANGED: no omega computation, identical RNG stream, identical result.
%
%   When cache.tvA is true, A0inv / By / Psi carry an extra trailing regime
%   dimension and every regime r is checked with its OWN transmission
%   A0inv(:,:,r), By(:,:,:,r), Psi(:,:,:,r). The shared rotation Q must
%   satisfy the restrictions in EVERY regime. When cache.tvA is false the
%   transmission is common (legacy heteroskedastic identification) and only
%   lambda varies; behaviour is then identical to the previous version.
%
%   Regime count for the loop is cache.nlambda (build_theta_cache expands it
%   to the A0-regime count under tvA, with lambda = ones for noLmd runs).
%
%   [FSIGN-CONSISTENCY PATCH 07/2026] The rotation Q is common across
%   regimes, so each restricted shock must satisfy its sign restrictions
%   with ONE orientation. The regime loop now rejects any draw whose
%   per-regime flip fsign_i disagrees with the baseline flip fsign on a
%   shock where both are active. Allowing the flip to vary by regime would
%   accept rotations outside the admissible set (a shock cannot carry one
%   orientation under the average transmission and the opposite one in some
%   regime). The step-(1) impact screen pins the orientation of every
%   REGIME check (a flipped column fails the as-is impact sign there), but
%   the baseline check at the average transmission carries no screen, so a
%   draw whose average IRF satisfies the signs only after a flip was
%   previously accepted with an orientation that disagrees across objects.
%   The guard rejects exactly those draws, which makes the effective
%   admissible set "signs hold as-is under the average transmission and in
%   every regime" -- the paper's set.

    %  omega defaults; only computed on the storage-stage call path.
    omega_hat  = 1;
    om_diag    = struct();
    want_omega = (nargout >= 2) && (nargin >= 10) && ~isempty(nrr_omega_opts);

    tvA      = isfield(cache, 'tvA') && cache.tvA;
    lambda   = cache.lambda;
    nlambda  = cache.nlambda;
    impact_s = cache.impact_signs;
    accshk   = cache.accepted_shock;

    use_fast = isfield(cache, 'SR_parsed') && ~isempty(cache.SR_parsed);

    % =================================================================
    % (1) FAST IMPACT SCREEN  (no IRF computation)
    % =================================================================
    has_imp_gate = isfield(cache, 'impact_signs_byregime') && ~isempty(cache.impact_signs_byregime);
    if ~isempty(impact_s) && any(impact_s(:) ~= 0)
        for m = 1:nlambda
            if has_imp_gate
                impact_m = cache.impact_signs_byregime{min(m, numel(cache.impact_signs_byregime))};
            else
                impact_m = impact_s;
            end
            mask_m = (impact_m ~= 0);
            if ~any(mask_m(:)), continue; end

            A0inv_m = local_slice3(cache.A0inv, m, tvA);
            Dhalf   = diag(sqrt(max(lambda(:, m), 0)));
            if sign_regime_dependent == 0
                smat_m = (A0inv_m * Q') * Dhalf;
            else
                smat_m = (A0inv_m * Dhalf) * Q';
            end
            prod = impact_m .* smat_m;
            if any(prod(mask_m) <= 0)
                pass = false;
                return
            end
        end
    end

    % =================================================================
    % (2) FULL IRF + SIGN-RESTRICTION CHECK
    % =================================================================
    A0inv_1 = local_slice3(cache.A0inv, 1, tvA);
    By_1    = local_slice4(cache.By,  1, tvA);
    Psi_1   = local_slice4(cache.Psi, 1, tvA);

    smat_baseline = A0inv_1 * Q';
    ir0 = local_compute_ir(smat_baseline, Psi_1, By_1, horizon);

    if isempty(SignRestrictions)
        fsign = ones(1, ny);
        irTemp_narrative_cell = cell(1, nlambda);
        fsign_cell            = cell(1, nlambda);
        for i = 1:nlambda
            A0inv_i = local_slice3(cache.A0inv, i, tvA);
            By_i    = local_slice4(cache.By,  i, tvA);
            Psi_i   = local_slice4(cache.Psi, i, tvA);
            smat_regime = local_regime_smat(A0inv_i, Q, lambda(:, i), sign_regime_dependent);
            irTemp_narrative_cell{i} = local_compute_ir(smat_regime, Psi_i, By_i, horizon);
            fsign_cell{i} = fsign;
        end
    else
        % --- Baseline check (regime 1, lambda = 1) ---
        if use_fast
            [ok0, fsign0] = checkrestrictions_per_shock_fast(cache.SR_parsed, ir0, 1e-12);
        else
            [ok0, fsign0, ~] = checkrestrictions_per_shock(SignRestrictions, ir0, 1e-12, true);
        end
        if ~ok0
            pass = false;
            return
        end
        fsign = fsign0;

        if nlambda == 1
            irTemp_narrative_cell    = cell(1, 1);
            fsign_cell               = cell(1, 1);
            irTemp_narrative_cell{1} = ir0 .* reshape(fsign0, 1, 1, []);
            fsign_cell{1}            = fsign0;
        else
            irTemp_narrative_cell = cell(1, nlambda);
            fsign_cell            = cell(1, nlambda);
            for i = 1:nlambda
                A0inv_i = local_slice3(cache.A0inv, i, tvA);
                By_i    = local_slice4(cache.By,  i, tvA);
                Psi_i   = local_slice4(cache.Psi, i, tvA);

                smat_regime = local_regime_smat(A0inv_i, Q, lambda(:, i), sign_regime_dependent);
                iri_raw     = local_compute_ir(smat_regime, Psi_i, By_i, horizon);

                if isfield(cache,'SR_parsed_byregime') && ~isempty(cache.SR_parsed_byregime)
                    SRP_i = cache.SR_parsed_byregime{min(i, numel(cache.SR_parsed_byregime))};
                else
                    SRP_i = cache.SR_parsed;
                end

                if use_fast
                    [ok_i, fsign_i] = checkrestrictions_per_shock_fast(SRP_i, iri_raw, 1e-12);
                else
                    [ok_i, fsign_i, ~] = checkrestrictions_per_shock(SignRestrictions, iri_raw, 1e-12, true);
                end

                if isfield(cache,'SR_parsed_byregime') && ~isempty(cache.SR_parsed_byregime) && ~isempty(fsign)
                    gated = (fsign_i(:).' == 0);
                    fsign_i(gated) = fsign(gated);
                end

                % [FSIGN-CONSISTENCY] one orientation per shock across the
                % baseline (average transmission) and every regime.
                incons = (fsign_i(:).' ~= 0) & (fsign(:).' ~= 0) & ...
                         (fsign_i(:).' ~= fsign(:).');
                if any(incons)
                    pass = false;
                    return
                end

                penalty = local_offdiag_penalty(Q, lambda(:, i), accshk, penalty_offdiagonal_on);

                if ~(ok_i && (penalty < 1))
                    pass = false;
                    return
                end

                fsign_cell{i}            = fsign_i;
                irTemp_narrative_cell{i} = iri_raw .* reshape(fsign_i, 1, 1, []);
            end
        end
    end

    % =================================================================
    % (3) NARRATIVE RESTRICTIONS  (point-anchor, optional)
    % =================================================================
    if ~isempty(NarrativeRestrictions)
        narr_ok = local_narrative_check(NarrativeRestrictions, cache, var, Q, ...
            fsign, fsign_cell, nlambda, irTemp_narrative_cell);
        if ~narr_ok
            pass = false;
            return
        end
    end

    % =================================================================
    % (4) REGIME / WINDOW-AVERAGE NARRATIVE RESTRICTIONS
    % =================================================================
    has_nrr = isfield(cache, 'NarrativeRegimeRestrictions') && ...
        ~isempty(cache.NarrativeRegimeRestrictions);
    if has_nrr
        narr_reg_ok = check_narrative_regime_avg(Q, cache, var, ...
            cache.NarrativeRegimeRestrictions, fsign_cell, fsign, ...
            irTemp_narrative_cell);
        if ~narr_reg_ok
            pass = false;
            return
        end
    end

    % =================================================================
    % (5) AD-RR omega weight at the surviving draw (storage stage
    %     only). Uses the SAME fsign / fsign_cell / irTemp_narrative_cell
    %     objects step (4) just used, so the simulated event is identical
    %     to the retained-draw check.
    % =================================================================
    has_nrr_point = ~isempty(NarrativeRestrictions);
    if want_omega && (has_nrr || has_nrr_point)
        pass_fns = {};
        if has_nrr_point   % point anchors (AD-RR NSR 4/5 style)
            pass_fns{end+1} = @(xi_sim) local_narrative_check(NarrativeRestrictions, ...
                cache, var, Q, fsign, fsign_cell, nlambda, irTemp_narrative_cell, xi_sim);
        end
        if has_nrr         % window-average anchors
            pass_fns{end+1} = @(xi_sim) check_narrative_regime_avg(Q, cache, var, ...
                cache.NarrativeRegimeRestrictions, fsign_cell, fsign, ...
                irTemp_narrative_cell, xi_sim);
        end
        if has_nrr
            NRR_for_dates = cache.NarrativeRegimeRestrictions;
        else
            NRR_for_dates = [];
        end
        use_fast_omega = ~isfield(nrr_omega_opts, 'fast') || nrr_omega_opts.fast;
        if use_fast_omega   % [PERF] vectorized omega, validated against pass_fns
            [omega_hat, om_diag] = estimate_nrr_omega_fast(Q, cache, var, ...
                NarrativeRestrictions, NRR_for_dates, fsign, fsign_cell, ...
                irTemp_narrative_cell, nrr_omega_opts, pass_fns);
        else
            [sim_dates, sim_regimes] = local_collect_narrative_dates( ...
                NarrativeRestrictions, NRR_for_dates, nlambda);
            [omega_hat, om_diag] = estimate_nrr_omega(Q, cache, var, ...
                sim_dates, sim_regimes, pass_fns, nrr_omega_opts);
        end
    end

    pass = true;
end


% =====================================================================
% LOCAL HELPERS  (unchanged)
% =====================================================================

function M = local_slice3(A, r, tvA)
    if tvA
        M = A(:, :, r);
    else
        M = A;
    end
end


function M = local_slice4(A, r, tvA)
    if isempty(A)
        M = A;
    elseif tvA
        M = A(:, :, :, r);
    else
        M = A;
    end
end


function ir = local_compute_ir(smat, Psi_r, By_r, horizon)
    if ~isempty(Psi_r) && size(Psi_r, 3) >= horizon
        [ny_local, n_shk] = size(smat);
        ir = zeros(ny_local, horizon, n_shk);
        for h = 1:horizon
            ir(:, h, :) = Psi_r(:, :, h) * smat;
        end
    else
        ir = impulsdtrf(By_r, smat, horizon);
    end
end


function smat = local_regime_smat(A0inv, Q, lambda_i, sign_regime_dependent)
    Dhalf = diag(sqrt(max(lambda_i, 0)));
    if sign_regime_dependent == 0
        smat = (A0inv * Q') * Dhalf;
    else
        smat = (A0inv * Dhalf) * Q';
    end
end


function penalty = local_offdiag_penalty(Q, lambda_i, accepted_shock, penalty_on)
    if ~penalty_on
        penalty = 0;
        return
    end
    Qlambda = Q * diag(lambda_i) * Q';
    K = length(accepted_shock);
    is_pen = false(K, 1);
    for k = 1:K
        j = accepted_shock(k);
        own_var = abs(Qlambda(j, j));
        cross   = abs(Qlambda(j, :));
        cross(j) = 0;
        is_pen(k) = own_var < 3 * max(cross);
    end
    penalty = sum(is_pen);
end


function ok = local_narrative_check(NR, cache, var, Q, fsign, fsign_cell, nlambda, ir_cell, xi_override)
%   [XI PATCH 07/2026] Labeled shocks are eps_tilde = Q * xi with
%   xi = Lambda^{-1/2} * eps. udraw_true stores eps on the original scale,
%   so each restricted date t is standardized by its OWN regime's
%   sqrt(lambda) before rotation. This restores sum_k HD_t(k) = u_{v,t}
%   (the regime IRFs already carry Dhalf) and evaluates the sign check on
%   the model's labeled shock. xi_override (optional, T x ny STANDARDIZED
%   shocks; non-restricted rows may be NaN) is used by estimate_nrr_omega.
    num_r = length(NR.time_index);

    use_override = (nargin >= 9) && ~isempty(xi_override);
    u_true = var.udraw_true;

    function v_t = local_v_row(t_, r_idx, fs)
        lam_r = sqrt(max(cache.lambda(:, r_idx), 0))';
        if any(lam_r <= 0), v_t = []; return; end
        if use_override
            xi_t = xi_override(t_, :);
        else
            xi_t = u_true(t_, :) ./ lam_r;
        end
        v_t = (xi_t * Q') .* reshape(fs, 1, []);
    end

    sign_pass = false(num_r, 1);
    for i = 1:num_r
        r_idx = local_resolve_regime(NR, i, nlambda);
        if r_idx < 1 || r_idx > nlambda || ~isfinite(r_idx)
            sign_pass(i) = false; continue
        end

        if nlambda > 1
            fs = fsign_cell{r_idx};
        else
            fs = fsign;
        end
        if isempty(fs)
            sign_pass(i) = false; continue
        end

        t_  = NR.time_index(i);
        s_  = NR.shock_index(i);
        sg  = NR.sign(i);
        v_t = local_v_row(t_, r_idx, fs);
        if isempty(v_t)
            sign_pass(i) = false; continue
        end
        val = v_t(s_);

        if sg > 0
            sign_pass(i) = val > 0;
        elseif sg < 0
            sign_pass(i) = val < 0;
        else
            sign_pass(i) = true;
        end
    end

    if ~all(sign_pass)
        ok = false; return
    end

    dom_pass = true(num_r, 1);
    for i = 1:num_r
        if isfield(NR, 'type') && length(NR.type) >= i
            type_i = NR.type{i};
        else
            type_i = 'dominance';
        end

        if strcmpi(type_i, 'sign')
            continue
        end

        r_idx = local_resolve_regime(NR, i, nlambda);
        if r_idx < 1 || r_idx > numel(ir_cell) || ~isfinite(r_idx)
            dom_pass(i) = false; continue
        end

        ir_reg = ir_cell{r_idx};
        if nlambda > 1
            fs_d = fsign_cell{r_idx};
        else
            fs_d = fsign;
        end

        t_ = NR.time_index(i);
        s_ = NR.shock_index(i);
        if isfield(NR, 'variable_index') && length(NR.variable_index) >= i
            v_idx = NR.variable_index(i);
        else
            v_idx = 1;
        end

        v_t = local_v_row(t_, r_idx, fs_d);
        if isempty(v_t)
            dom_pass(i) = false; continue
        end
        HD  = getHDs_fast(ir_reg, v_t, v_idx);
        tgt = abs(HD(s_));

        switch lower(type_i)
            case 'dominance_max'
                others_idx  = setdiff(1:length(HD), s_);
                max_other   = max(abs(HD(others_idx)));
                dom_pass(i) = (tgt > max_other);
            case 'dominance_min'
                others_idx  = setdiff(1:length(HD), s_);
                min_other   = min(abs(HD(others_idx)));
                dom_pass(i) = (tgt < min_other);
            otherwise
                oth         = sum(abs(HD)) - tgt;
                dom_pass(i) = (tgt > oth);
        end
    end

    ok = all(dom_pass);
end


function r_idx = local_resolve_regime(NR, i, nlambda)
    if nlambda == 1
        r_idx = 1; return
    end
    if isfield(NR, 'regime_indices') && ~isempty(NR.regime_indices) && ...
            length(NR.regime_indices) >= i
        r_idx = NR.regime_indices(i);
    elseif isfield(NR, 'regime_index') && ~isempty(NR.regime_index)
        r_idx = NR.regime_index;
    else
        r_idx = 1;
    end
end


function [dates, regimes] = local_collect_narrative_dates(NR_point, NRR_regime, nlambda)
%LOCAL_COLLECT_NARRATIVE_DATES  Union of all narrative-restricted dates with
% their regime indices, for the omega simulator.
    dates = []; regimes = [];
    if ~isempty(NR_point)
        for i = 1:length(NR_point.time_index)
            dates(end+1, 1)   = NR_point.time_index(i);            %#ok<AGROW>
            regimes(end+1, 1) = local_resolve_regime(NR_point, i, nlambda); %#ok<AGROW>
        end
    end
    if ~isempty(NRR_regime)
        for k = 1:numel(NRR_regime)
            w = NRR_regime(k).window(:);
            dates   = [dates;   w];                                 %#ok<AGROW>
            regimes = [regimes; repmat(NRR_regime(k).regime_idx, numel(w), 1)]; %#ok<AGROW>
        end
    end
end