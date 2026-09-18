function [SignCheck, accepted_shock, irf_save, Aplus, Qdraw, x1_rotated, stability, lambda, By, fsign, fsign_cell] = ...
    SignRestrictionCheck(x1, var1, SignRestrictions, NarrativeRestrictions, lags, lcA0, constant, ny, horizon, ...
    sign_regime_dependent, penalty_offdiagonal_on, lcLmd, fix_first, Qdraw, SR_parsed, NarrativeRegimeRestrictions, ...
    enforce_stability, lcA0_tv, SR_parsed_byregime)
%--------------------------------------------------------------------------
% SignRestrictionCheck  (tvA + PARTIAL-tvA AWARE)
%   Verification-pass check for sign + narrative restrictions, consistent
%   with check_sr_pass_cached.m.
%
%   Time-varying A0 (full or partial):
%   New trailing arg lcA0_tv (subset of lcA0 that is regime-specific):
%     []          -> legacy: full tvA / non-tvA inferred from numel(x1)
%     given       -> partial tvA: A0 = [common once] + [tv per regime],
%                    seedx layout [ a_cm ; a_tv x nReg ; lambda ].
%   A0 is reconstructed via the shared helper unpackA0_tvA so this routine
%   agrees with build_theta_cache / bvar_posterior_tvA / the driver exactly.
%
%   Under tvA the shared Q must satisfy the restrictions in EVERY regime;
%   lambda is expanded to ny x nRegimes (ones for a noLmd run). The returned
%   By is 4-D (ny x ny x lags x nRegimes) under tvA. irf_save is the regime-1
%   representative IRF (per-regime IRFs are reconstructed in the driver).
%
%   Homoskedastic-shock handling:
%   Two changes so a homoskedastic shock (all-false lcLmd row -> lambda fixed
%   to 1 across regimes) is handled consistently with the cache / posterior:
%     (1) extractLambda: zero ONLY rows carrying free lambda (not all rows),
%         so homoskedastic rows keep the ones-init and pass the positivity
%         gate as a flat row of 1s. The old wholesale zeroing left them at 0,
%         which the "any(lambda(:) <= 0)" validation below would REJECT.
%     (2) rotateLambdaParameters: write only the FREE lambda cells (column-
%         major, matching the seed packing), instead of always writing
%         ny*(nlambda-1) entries. With homoskedastic shocks (or fix_first)
%         the old loop overran/misaligned the "_wQ" diagnostic lambda block.
%
%   x1_rotated bakes Q into A0 for the diagnostic "_wQ" posterior. Under
%   PARTIAL tvA the rotation A0_rot = Q*A0_r generally breaks the common
%   structure (a common element rotates into different values per regime), so
%   x1_rotated keeps the REDUCED layout by taking the regime-1 rotation for
%   the common block and the per-regime rotation for the tv block. This is a
%   representative diagnostic only (xrot feeds the "_wQ" outputs, never the
%   IRF/FEVD reconstruction).
%
%   STABILITY (enforce_stability, 17th arg, default true): per-regime under
%   tvA (reject if ANY regime companion has a root of modulus >= 1).
%
%   Requires (standalone): unpackA0_tvA.m
%--------------------------------------------------------------------------

    % === Input Validation ===
    if nargin < 13 || isempty(fix_first)
        fix_first = false;
    end
    fix_first = logical(fix_first);

    if nargin < 15
        SR_parsed = [];
    end
    if nargin < 16
        NarrativeRegimeRestrictions = [];
    end
    if nargin < 17 || isempty(enforce_stability)
        enforce_stability = true;
    end
    if nargin < 18
        lcA0_tv = [];                 % [] -> legacy detection
    end
    if nargin < 19 || isempty(SR_parsed_byregime)
        SR_parsed_byregime = [];      % [] -> no regime gate (full set everywhere)
    end

    % === Detect tvA, regime count and the common/tv partition ===
    nA        = sum(lcA0(:));
    nLmd_free = sum(lcLmd(:));

    if isempty(lcA0_tv)
        % legacy: full tvA or non-tvA inferred from length
        nA_blocks = round((numel(x1) - nLmd_free) / nA);
        if nA_blocks < 1, nA_blocks = 1; end
        nReg    = nA_blocks;
        lcA0_tv = lcA0;            % full tvA (nReg == 1 -> non-tvA)
        nA_cm   = 0;  nA_tv = nA;
    else
        % partial tvA: explicit common/tv split
        lcA0_cm = lcA0 & ~lcA0_tv;
        nA_cm   = sum(lcA0_cm(:));
        nA_tv   = sum(lcA0_tv(:));
        if nA_tv > 0
            nReg = round((numel(x1) - nLmd_free - nA_cm) / nA_tv);
        else
            nReg = 1;
        end
        if nReg < 1, nReg = 1; end
    end
    tvA      = nReg > 1;
    A0offset = nA_cm + nA_tv * nReg;     % lambda block start (full/partial/non-tvA)

    % === Extract VAR Coefficients (A0 / By per regime under tvA) ===
    [A0, Aplus, By] = extractVARCoefficients(x1, var1, lcA0, lcA0_tv, ny, lags, tvA, nReg);

    % === Extract Lambda (expanded to nReg under tvA) ===
    [lambda, nlambda] = extractLambda(x1, A0offset, lcLmd, ny, fix_first, tvA, nReg);

    % === Lambda Validation ===
    if any(lambda(:) <= 0) || any(~isfinite(lambda(:)))
        [SignCheck, accepted_shock, irf_save, Aplus, Qdraw, x1_rotated, stability, fsign, fsign_cell] = ...
            deal(0, [], [], [], [], [], 0, [], []);
        By = []; lambda = [];
        return
    end

    % === Stability Check (per regime under tvA; optional rejection) ===
    [isStable, stability] = checkStability(By, ny, lags, tvA, nReg);
    if enforce_stability && ~isStable
        [SignCheck, accepted_shock, irf_save, Aplus, Qdraw, x1_rotated, fsign, fsign_cell] = ...
            deal(0, [], [], [], [], [], [], []);
        return
    end
    if ~enforce_stability
        stability = 1;
    end

    % === Identify Shocks from Restrictions ===
    accepted_shock = identifyShocks(SignRestrictions, NarrativeRestrictions);

    % === Check Sign Restrictions (per-regime transmission under tvA) ===
    [SignCheck, fsign, smat, irTemp_narrative_cell, fsign_cell] = ...
        checkSignRestrictions(SignRestrictions, A0, Qdraw, By, lambda, nlambda, ...
        horizon, sign_regime_dependent, penalty_offdiagonal_on, accepted_shock, ny, SR_parsed, tvA, ...
        SR_parsed_byregime);

    % === Check Narrative Restrictions (point-anchor) ===
    if ~isempty(NarrativeRestrictions) && SignCheck
        SignCheck = checkNarrativeRestrictionsWrapper(NarrativeRestrictions, var1, Qdraw, ...
            fsign, fsign_cell, nlambda, irTemp_narrative_cell, lambda);
    end

    % === Check Regime/Window-Average Narrative Restrictions ===
    if SignCheck && ~isempty(NarrativeRegimeRestrictions)
        cache_for_NRR.lambda  = lambda;
        cache_for_NRR.nlambda = nlambda;
        SignCheck = check_narrative_regime_avg(Qdraw, cache_for_NRR, var1, ...
            NarrativeRegimeRestrictions, fsign_cell, fsign, irTemp_narrative_cell);
    end

    % === Compute Outputs if Restrictions Satisfied ===
    if SignCheck
        [irf_save, x1_rotated] = computeOutputs(x1, A0, Qdraw, By, smat, lambda, ...
            lcA0, lcA0_tv, lcLmd, nlambda, ny, horizon, accepted_shock, fsign, ...
            sign_regime_dependent, tvA, nReg);
    else
        irf_save = [];
        x1_rotated = [];
    end
end


%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function [A0, Aplus, By] = extractVARCoefficients(x1, var1, lcA0, lcA0_tv, ny, lags, tvA, nReg)
    % Aplus = common structural lag block (rfvar3 estimates one B; A0 carries
    % the regime under tvA). A0 is rebuilt via the shared partial-tvA unpacker.
    Btrim = var1.Bdraw(1:lags*ny, :);
    Aplus = NaN(ny, ny, lags);
    for iLag = 1:lags
        rows = (iLag - 1) * ny + (1:ny);
        Aplus(:, :, iLag) = Btrim(rows, :)';
    end

    A0stack = unpackA0_tvA(x1, ny, lcA0, lcA0_tv, nReg);   % ny x ny x nReg

    if tvA
        A0 = A0stack;
        By = zeros(ny, ny, lags, nReg);
        for r = 1:nReg
            for iLag = 1:lags
                By(:, :, iLag, r) = A0stack(:, :, r) \ Aplus(:, :, iLag);
            end
        end
    else
        A0 = A0stack(:, :, 1);
        By = zeros(ny, ny, lags);
        for iLag = 1:lags
            By(:, :, iLag) = A0 \ Aplus(:, :, iLag);
        end
    end
end


function [lambda, nlambda] = extractLambda(x1, offA0, lcLmd, ny, fix_first, tvA, nReg)
    if all(lcLmd(:) == 0)
        lambda  = ones(ny, 1);
        nlambda = 1;
    else
        nlambda = size(lcLmd, 2);
        lambda  = ones(ny, nlambda);

        if fix_first
            lcLmd(:, 1)  = false;
            lambda(:, 1) = 1;
        else
            % HOMOSKEDASTIC-SHOCK FIX: zero ONLY rows carrying free lambda.
            % Homoskedastic shocks (all-false lcLmd row) keep the ones-init so
            % they reconstruct to a flat row of 1s below and pass the
            % positivity gate. The old "lambda(:,:) = 0" left them at 0, which
            % the "any(lambda(:) <= 0)" validation in the caller would reject.
            free_rows = any(lcLmd, 2);
            lambda(free_rows, :) = 0;
        end

        nFree = sum(lcLmd(:));
        if nFree > 0
            lambda(lcLmd) = x1(offA0 + (1:nFree));
        end

        if fix_first
            if nlambda >= 3
                lambda(:, end) = (nlambda - 1) - sum(lambda(:, 2:end-1), 2);
            else
                lambda(:, end) = 1;
            end
        else
            lambda(:, end) = nlambda - sum(lambda(:, 1:end-1), 2);
        end
    end

    % Expand to A0-regime count under tvA (noLmd -> ones across regimes)
    if tvA && nlambda < nReg
        lambda  = ones(ny, nReg);
        nlambda = nReg;
    end
end


function [isStable, stability] = checkStability(By, ny, lags, tvA, nReg)
    if tvA
        isStable = true;
        for r = 1:nReg
            if ~local_stab(By(:, :, :, r), ny, lags)
                isStable = false; break
            end
        end
    else
        isStable = local_stab(By, ny, lags);
    end
    stability = double(isStable);
end


function ok = local_stab(By, ny, lags)
    Phi = reshape(By, ny, ny * lags);
    if lags > 1
        Companion = [Phi; eye(ny*(lags-1)), zeros(ny*(lags-1), ny)];
    else
        Companion = Phi;
    end
    ok = all(abs(eig(Companion)) < 1);
end


function accepted_shock = identifyShocks(SignRestrictions, NarrativeRestrictions)
    shock_ids = [];
    if ~isempty(SignRestrictions)
        shock_token_cells = cellfun(@(s) regexp(s, ',\s*(\d+)\)', 'tokens', 'once'), ...
            SignRestrictions, 'UniformOutput', false);
        shock_ids = cellfun(@(c) sscanf(c{1}, '%d'), shock_token_cells);
    end
    if ~isempty(NarrativeRestrictions)
        shock_ids = [shock_ids, NarrativeRestrictions.shock_index(:)'];
    end
    accepted_shock = unique(shock_ids);
end


function [SignCheck, fsign, smat, irTemp_narrative_cell, fsign_cell] = ...
    checkSignRestrictions(SignRestrictions, A0, Qdraw, By, lambda, nlambda, ...
    horizon, sign_regime_dependent, penalty_offdiagonal_on, accepted_shock, ny, SR_parsed, tvA, ...
    SR_parsed_byregime)

    if nargin < 12, SR_parsed = []; end
    if nargin < 13, tvA = false;    end
    if nargin < 14 || isempty(SR_parsed_byregime), SR_parsed_byregime = []; end
    use_fast = ~isempty(SR_parsed);

    SignCheck = 0;
    fsign = [];

    A0_1 = local_a0slice(A0, 1, tvA);
    By_1 = local_byslice(By, 1, tvA);
    smat = A0_1 \ Qdraw';                 % baseline: regime 1, lambda = 1

    irTemp_narrative_cell = cell(1, nlambda);
    fsign_cell            = cell(1, nlambda);

    % --- No sign restrictions ---
    if isempty(SignRestrictions)
        fsign = ones(1, ny);
        for i = 1:nlambda
            A0_i = local_a0slice(A0, i, tvA);
            By_i = local_byslice(By, i, tvA);
            smat_regime = computeRegimeSmat(A0_i, Qdraw, lambda(:, i), sign_regime_dependent);
            irTemp_narrative_cell{i} = impulsdtrf(By_i, smat_regime, horizon);
            fsign_cell{i} = fsign;
        end
        SignCheck = 1;
        return
    end

    % ============================================================
    % (A) Baseline check (regime 1, lambda = 1)
    % ============================================================
    ir0 = impulsdtrf(By_1, smat, horizon);

    if use_fast
        [ok0, fsign0] = checkrestrictions_per_shock_fast(SR_parsed, ir0, 1e-12);
    else
        [ok0, fsign0, ~] = checkrestrictions_per_shock(SignRestrictions, ir0, 1e-12, true);
    end
    if ~ok0
        return
    end
    fsign = fsign0;

    if nlambda == 1
        irTemp_narrative_cell{1} = ir0 .* reshape(fsign0, 1, 1, []);
        fsign_cell{1} = fsign0;
        SignCheck = 1;
        return
    end

    % ============================================================
    % (B) Multi-regime: each regime uses its OWN transmission + fsign.
    %     REGIME GATE: in regime i, enforce only the restrictions in
    %     SR_parsed_byregime{i} (full set when no gate provided). The
    %     baseline check (A) above always uses the FULL set, so a gated-out
    %     shock keeps a clean, regime-agnostic sign orientation.
    % ============================================================
    for i = 1:nlambda
        A0_i = local_a0slice(A0, i, tvA);
        By_i = local_byslice(By, i, tvA);

        smat_regime = computeRegimeSmat(A0_i, Qdraw, lambda(:, i), sign_regime_dependent);
        iri_raw     = impulsdtrf(By_i, smat_regime, horizon);

        % regime-gated parsed set (fall back to the full set if no gate)
        if isempty(SR_parsed_byregime)
            SRP_i = SR_parsed;
        else
            SRP_i = SR_parsed_byregime{min(i, numel(SR_parsed_byregime))};
        end

        if use_fast
            [ok_i, fsign_i] = checkrestrictions_per_shock_fast(SRP_i, iri_raw, 1e-12);
        else
            % slow path is NOT gated (only used when SR_parsed is empty).
            [ok_i, fsign_i, ~] = checkrestrictions_per_shock(SignRestrictions, iri_raw, 1e-12, true);
        end

        % shocks gated out of this regime keep the baseline orientation
        % (irrelevant for R1/R2 here -- no narrative anchors there -- but
        %  protects the per-regime fsign_cell/HD machinery if one is added).
        if ~isempty(SR_parsed_byregime) && ~isempty(fsign)
            gated = (fsign_i(:).' == 0);
            fsign_i(gated) = fsign(gated);
        end

        penalty = computeOffDiagonalPenalty(Qdraw, lambda(:, i), accepted_shock, penalty_offdiagonal_on);

        if ~(ok_i && (penalty < 1))
            return
        end

        fsign_cell{i}            = fsign_i;
        irTemp_narrative_cell{i} = iri_raw .* reshape(fsign_i, 1, 1, []);
    end

    SignCheck = 1;
end


function M = local_a0slice(A0, r, tvA)
    if tvA, M = A0(:, :, r); else, M = A0; end
end


function M = local_byslice(By, r, tvA)
    if tvA, M = By(:, :, :, r); else, M = By; end
end


function smat = computeRegimeSmat(A0, Qdraw, lambda_i, sign_regime_dependent)
    Dhalf = diag(sqrt(max(lambda_i, 0)));
    if sign_regime_dependent == 0
        smat = (A0 \ Qdraw') * Dhalf;
    else
        smat = (A0 \ Dhalf) * Qdraw';
    end
end


function penalty = computeOffDiagonalPenalty(Qdraw, lambda_i, accepted_shock, penalty_offdiagonal_on)
    if ~penalty_offdiagonal_on
        penalty = 0;
        return
    end

    Qlambda = Qdraw * diag(lambda_i) * Qdraw';
    is_penalty = false(length(accepted_shock), 1);

    for k = 1:length(accepted_shock)
        j = accepted_shock(k);
        own_var = abs(Qlambda(j, j));
        cross_var = abs(Qlambda(j, :));
        cross_var(j) = 0;
        is_penalty(k) = own_var < 3 * max(cross_var);
    end

    penalty = sum(is_penalty);
end


function SignCheck = checkNarrativeRestrictionsWrapper(NarrativeRestrictions, var1, Qdraw, ...
    fsign, fsign_cell, nlambda, irTemp_narrative_cell, lambda)
%   [XI PATCH 07/2026] Labeled shocks are eps_tilde = Q * xi with
%   xi = Lambda^{-1/2} * eps; udraw_true stores eps on the original scale.
%   Each restricted date is standardized by its OWN regime's sqrt(lambda)
%   before rotation, mirroring local_narrative_check in check_sr_pass_cached
%   EXACTLY so the storage-stage recheck agrees with the hot loop.

    u_true = var1.udraw_true;
    num_restrictions = length(NarrativeRestrictions.time_index);
    sign_checks_passed = false(num_restrictions, 1);

    for nar_idx = 1:num_restrictions
        regime_idx = getRegimeIndex(NarrativeRestrictions, nar_idx, nlambda);

        if isempty(regime_idx) || ~isfinite(regime_idx) || ...
                regime_idx < 1 || regime_idx > nlambda
            sign_checks_passed(nar_idx) = false;
            continue
        end

        if nlambda > 1
            fs = fsign_cell{regime_idx};
        else
            fs = fsign;
        end
        if isempty(fs)
            sign_checks_passed(nar_idx) = false;
            continue
        end

        t_idx = NarrativeRestrictions.time_index(nar_idx);
        s_idx = NarrativeRestrictions.shock_index(nar_idx);
        sgn   = NarrativeRestrictions.sign(nar_idx);

        lam_r = sqrt(max(lambda(:, regime_idx), 0))';
        if any(lam_r <= 0)
            sign_checks_passed(nar_idx) = false;
            continue
        end
        v_t = ((u_true(t_idx, :) ./ lam_r) * Qdraw') .* reshape(fs, 1, []);
        val = v_t(s_idx);
        if sgn > 0
            sign_checks_passed(nar_idx) = val > 0;
        elseif sgn < 0
            sign_checks_passed(nar_idx) = val < 0;
        else
            sign_checks_passed(nar_idx) = true;   % sgn == 0: no sign restriction
        end
    end

    NarrativeCheck = all(sign_checks_passed);

    dominance_checks_passed = checkDominance(NarrativeRestrictions, var1, Qdraw, ...
        fsign, fsign_cell, nlambda, irTemp_narrative_cell, lambda);

    SignCheck = NarrativeCheck && all(dominance_checks_passed);
end


function dominance_checks_passed = checkDominance(NarrativeRestrictions, var1, Qdraw, ...
    fsign, fsign_cell, nlambda, irTemp_narrative_cell, lambda)
%   [XI PATCH 07/2026] Same standardization as the wrapper above: the HD is
%   getHDs_fast(regime IRF with Dhalf, standardized-and-rotated shock row),
%   which restores sum_k HD_t(k) = u_{v,t}.
    u_true = var1.udraw_true;
    num_restrictions = length(NarrativeRestrictions.time_index);
    dominance_checks_passed = true(num_restrictions, 1);

    for nar_idx = 1:num_restrictions

        if isfield(NarrativeRestrictions, 'type') && ...
                length(NarrativeRestrictions.type) >= nar_idx
            type_i = NarrativeRestrictions.type{nar_idx};
        else
            type_i = 'dominance';
        end

        if strcmpi(type_i, 'sign')
            continue
        end

        regime_idx = getRegimeIndex(NarrativeRestrictions, nar_idx, nlambda);

        if isempty(regime_idx) || ~isfinite(regime_idx) || ...
                regime_idx < 1 || regime_idx > numel(irTemp_narrative_cell)
            dominance_checks_passed(nar_idx) = false;
            continue
        end

        if nlambda > 1
            irTemp_regime = irTemp_narrative_cell{regime_idx};
            fs_d = fsign_cell{regime_idx};
        else
            irTemp_regime = irTemp_narrative_cell{1};
            fs_d = fsign;
        end

        t_idx = NarrativeRestrictions.time_index(nar_idx);
        s_idx = NarrativeRestrictions.shock_index(nar_idx);
        v_idx = getVariableIndex(NarrativeRestrictions, nar_idx);

        lam_r = sqrt(max(lambda(:, regime_idx), 0))';
        if any(lam_r <= 0)
            dominance_checks_passed(nar_idx) = false;
            continue
        end
        v_t = ((u_true(t_idx, :) ./ lam_r) * Qdraw') .* reshape(fs_d, 1, []);
        HD                 = getHDs_fast(irTemp_regime, v_t, v_idx);
        abs_contrib_target = abs(HD(s_idx));

        switch lower(type_i)
            case 'dominance_max'
                others_idx = setdiff(1:length(HD), s_idx);
                max_other  = max(abs(HD(others_idx)));
                dominance_checks_passed(nar_idx) = (abs_contrib_target > max_other);
            case 'dominance_min'
                others_idx = setdiff(1:length(HD), s_idx);
                min_other  = min(abs(HD(others_idx)));
                dominance_checks_passed(nar_idx) = (abs_contrib_target < min_other);
            otherwise   % 'dominance' (legacy strict, sum-based)
                abs_contrib_others = sum(abs(HD)) - abs_contrib_target;
                dominance_checks_passed(nar_idx) = ...
                    (abs_contrib_target > abs_contrib_others);
        end
    end
end


function regime_idx = getRegimeIndex(NarrativeRestrictions, nar_idx, nlambda)
    if nlambda == 1
        regime_idx = 1;
        return
    end

    if isfield(NarrativeRestrictions, 'regime_indices') && ...
            ~isempty(NarrativeRestrictions.regime_indices) && ...
            length(NarrativeRestrictions.regime_indices) >= nar_idx
        regime_idx = NarrativeRestrictions.regime_indices(nar_idx);
    elseif isfield(NarrativeRestrictions, 'regime_index') && ...
            ~isempty(NarrativeRestrictions.regime_index)
        regime_idx = NarrativeRestrictions.regime_index;
    else
        error(['NarrativeRestrictions requires either ''regime_indices'' ' ...
               '(per-anchor vector) or ''regime_index'' (single scalar) ' ...
               'when nlambda > 1. Neither field was provided.']);
    end
end


function v_idx = getVariableIndex(NarrativeRestrictions, nar_idx)
    if isfield(NarrativeRestrictions, 'variable_index')
        if length(NarrativeRestrictions.variable_index) >= nar_idx
            v_idx = NarrativeRestrictions.variable_index(nar_idx);
        else
            v_idx = NarrativeRestrictions.variable_index(1);
        end
    else
        v_idx = 1;
    end
end


function [irf_save, x1_rotated] = computeOutputs(x1, A0, Qdraw, By, smat, lambda, ...
    lcA0, lcA0_tv, lcLmd, nlambda, ny, horizon, accepted_shock, fsign, sign_regime_dependent, tvA, nReg)

    nA      = sum(lcA0(:));
    n_shock = length(accepted_shock);

    % --- Representative (regime 1) IRF: smat = A0_1 \ Qdraw' (lambda = 1) ---
    irf_save = zeros(ny, horizon, n_shock);
    By_1     = local_byslice(By, 1, tvA);
    irTemp   = impulsdtrf(By_1, smat, horizon);
    for j = 1:n_shock
        s = accepted_shock(j);
        irf_save(:, :, j) = fsign(s) * irTemp(:, :, s);
    end

    % --- Rotate A0 into x1_rotated, preserving the (partial-)tvA layout ---
    %   Under partial tvA the rotation A0_rot = Q*A0_r breaks the common
    %   structure (a common element rotates differently per regime), so the
    %   common block uses the regime-1 rotation as a representative value
    %   (xrot feeds the diagnostic "_wQ" posterior only). For full tvA
    %   (lcA0_tv == lcA0, nA_cm == 0) this reduces to the legacy per-regime
    %   rotation; for non-tvA to the single rotation.
    x1_rotated = x1;
    if tvA
        lcA0_cm = lcA0 & ~lcA0_tv;
        nA_cm   = sum(lcA0_cm(:));
        nA_tv   = sum(lcA0_tv(:));

        % common block: regime-1 rotation (representative)
        if nA_cm > 0
            A01_rot = Qdraw * A0(:, :, 1);
            x1_rotated(1:nA_cm) = A01_rot(lcA0_cm);
        end
        % tv block: per-regime rotation
        for r = 1:nReg
            A0r_rot = Qdraw * A0(:, :, r);
            x1_rotated(nA_cm + (r-1)*nA_tv + (1:nA_tv)) = A0r_rot(lcA0_tv);
        end
        offset_start = nA_cm + nA_tv * nReg;
    else
        A0_rotated = Qdraw * A0;
        x1_rotated(1:nA) = A0_rotated(lcA0);
        offset_start = nA;
    end

    % --- Rotate lambda parameters only when free lambda exists ---
    if sum(lcLmd(:)) > 0
        x1_rotated = rotateLambdaParameters(x1_rotated, Qdraw, lambda, offset_start, lcLmd, ...
            nlambda, ny, sign_regime_dependent);
    end

    x1_rotated = x1_rotated(:);
end


function x1_rotated = rotateLambdaParameters(x1_rotated, Qdraw, lambda, offset_start, lcLmd, ...
    nlambda, ny, sign_regime_dependent)
% HOMOSKEDASTIC-SHOCK FIX: write ONLY the free lambda cells (in column-major
% order, matching the seed packing seedLmd(lcLmd) and constructLambda's
% lmd(lcLmd) = seedx(...)). The previous version always wrote ny*(nlambda-1)
% entries, which overran / misaligned the lambda block of x1_rotated whenever
% free cells were fewer than that (homoskedastic shocks, or fix_first). This
% only affects the diagnostic "_wQ" posterior outputs, never the IRFs.
%
% In the normal case (all rows free in regimes 1..nlambda-1) find(lcLmd)
% reproduces the old cell set and column-major order exactly, so behaviour is
% byte-identical there.

    % rotated lambda diagonal, per regime
    rot_diag = zeros(ny, nlambda);
    for r = 1:nlambda
        if sign_regime_dependent == 0
            rot_diag(:, r) = lambda(:, r);
        else
            Lr = Qdraw * diag(lambda(:, r)) * Qdraw';
            rot_diag(:, r) = diag(Lr);
        end
    end

    free_idx = find(lcLmd);     % column-major linear indices of free lambda cells
    x1_rotated(offset_start + (1:numel(free_idx))) = rot_diag(free_idx);
end