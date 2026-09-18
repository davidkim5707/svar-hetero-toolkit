function [omega_hat, om_diag] = estimate_nrr_omega_fast(Q, cache, var, ...
    NR_point, NRR_reg, fsign, fsign_cell, ir_cell, opts, pass_fns_exact)
%ESTIMATE_NRR_OMEGA_FAST  Vectorized AD-RR omega at fixed (theta, Q).
%
% Same estimand as estimate_nrr_omega (v3): omega = Pr[all narrative checks
% pass | fresh standardized shocks xi on the restricted dates]. The slow
% path calls the exact checkers once PER SIMULATION (2M function calls for
% M = 2,000 and 1,000 stored draws); this path evaluates all M simulations
% with one matmul per restricted date plus (M x ny) array comparisons,
% replicating the checkers' impact-only HD semantics EXACTLY:
%
%     v^(m)  = (xi^(m) * Q') .* fs_r          per date, regime r
%     HD(j)  = ir_cell{r}(v_idx, 1, j) * v(j) (single-row getHDs_fast)
%
% SAFETY GATE: when opts.validate is true (set it on the first stored draw),
% the first K = min(200, M) simulations are replayed through the EXACT
% checker handles (pass_fns_exact, xi_override path) and the per-simulation
% indicators are asserted identical. Any divergence throws, so the fast
% path can never silently change the estimand.
%
% [XI HOIST 07/2026] The xi sample is drawn ONCE per (nW, ny, M, law, nu)
% configuration and cached in a persistent variable, then REUSED across all
% stored draws (common random numbers). Validity: omega(theta,Q) =
% E_xi[indicator] with the SAME xi law at every draw, so each omega_hat
% remains an unbiased MC estimate at its own (theta, Q); reuse only
% correlates the MC errors across draws, which is harmless for the
% self-normalized importance weights (and reduces the variance of RELATIVE
% weights, since the common error component partially cancels in ratios).
% The validation gate replays the SAME cached xi through the exact
% checkers, so it is unaffected. Disable with opts.hoist_xi = false to
% recover per-draw fresh sampling. NOTE: the persistent cache survives
% within a MATLAB session -- call `clear estimate_nrr_omega_fast` at the
% top of the driver so every run starts from a fresh, seed-determined xi.
% chi2rnd is replaced by the distribution-identical 2*randg(nu/2,.)
% (chi2(nu) = 2*Gamma(nu/2,1)); this changes the RNG stream consumption but
% not the estimand.
%
% Mirrored semantics (from local_narrative_check / check_narrative_regime_avg):
%   point:  sign (sgn = 0 auto-pass), type 'sign' skips dominance,
%           'dominance_max' / 'dominance_min' / legacy 'dominance' (strict,
%           sum-based); regime via regime_indices(i) -> regime_index -> 1;
%           variable_index default 1; invalid regime => restriction fails
%           (all simulations fail).
%   window: focal sign on window mean, dominance on mean |HD| with modes
%           'on'/'max', 'min', 'strict', anti-dominance 'min';
%           other_idx / anti_other_idx overrides honored.
%
% PRODUCT FORM (opts.product_form, default true): shocks are iid across
% dates, so restrictions touching DISJOINT date sets define independent
% events and omega = prod over groups of omega_g. Estimating each group's
% marginal and multiplying targets the SAME omega with far smaller variance
% than the joint indicator -- essential when the joint event is rare (the
% fiscal strict-dominance case, where the joint estimator floors at 1/M).
% Groups are connected components of restrictions linked by shared dates.
% The joint per-simulation indicator is still produced for the validation
% gate, which is unchanged.

    if nargin < 9  || isempty(opts), opts = struct(); end
    if nargin < 10, pass_fns_exact = {}; end
    M         = get_field(opts, 'M', 2000);
    shock_law = get_field(opts, 'shock_law', 'student');
    nu        = get_field(opts, 'nu', 5.703233415698);
    validate  = get_field(opts, 'validate', false);
    prodform  = get_field(opts, 'product_form', true);
    momchk    = get_field(opts, 'moment_check', false);

    om_diag = struct('n_pass', 0, 'M', M, 'floored', false, 'validated', false);
    nlambda = cache.nlambda;
    [Tpl, ny] = size(var.udraw_true);

    % ---------- collect restricted dates with regimes ----------
    dates = []; regs = [];
    has_pt = ~isempty(NR_point);
    if has_pt
        for i = 1:length(NR_point.time_index)
            dates(end+1,1) = NR_point.time_index(i);                 %#ok<AGROW>
            regs(end+1,1)  = resolve_regime(NR_point, i, nlambda);   %#ok<AGROW>
        end
    end
    has_win = ~isempty(NRR_reg);
    if has_win
        for e = 1:numel(NRR_reg)
            w = NRR_reg(e).window(:);
            dates = [dates; w];  regs = [regs; repmat(NRR_reg(e).regime_idx, numel(w), 1)]; %#ok<AGROW>
        end
    end
    if isempty(dates), omega_hat = 1; return; end
    [t_uniq, iu] = unique(dates, 'stable');
    r_uniq = regs(iu);
    for k = 1:numel(t_uniq)
        if numel(unique(regs(dates == t_uniq(k)))) > 1
            error('estimate_nrr_omega_fast: date %d maps to two regimes.', t_uniq(k));
        end
    end
    if max(t_uniq) > Tpl || min(t_uniq) < 1
        error('estimate_nrr_omega_fast: narrative date outside udraw_true range.');
    end
    nW = numel(t_uniq);
    row_of = containers.Map(t_uniq, 1:nW);

    % ---------- draw all standardized shocks once, reuse across draws ------
    % [XI HOIST 07/2026] see header. Keyed on (nW, ny, M, law, nu) so any
    % configuration change forces a fresh sample within the session.
    persistent XI_cache XI_key
    hoist = get_field(opts, 'hoist_xi', true);
    key   = [nW, ny, M, double(strcmpi(shock_law, 'student')), nu];
    if hoist && ~isempty(XI_cache) && isequal(XI_key, key)
        XI = XI_cache;
    else
        XI = randn(nW, ny, M);
        if strcmpi(shock_law, 'student')
            if nu <= 2, error('estimate_nrr_omega_fast: nu <= 2.'); end
            % chi2(nu) = 2 * Gamma(nu/2, 1): identical law, cheaper than chi2rnd.
            XI = XI ./ sqrt(2 * randg(nu/2, nW, ny, M) / nu);
        end
        if hoist
            XI_cache = XI;
            XI_key   = key;
        end
    end

    % ---------- rotate + flip: one matmul per date ----------
    V = cell(nW, 1);                       % each (M x ny), flipped
    ok_setup = true;
    for k = 1:nW
        fs = get_fs(r_uniq(k), fsign, fsign_cell, nlambda);
        if isempty(fs), ok_setup = false; break; end
        V{k} = (reshape(XI(k, :, :), ny, M).' * Q.') .* reshape(fs, 1, []);
    end

    % per-restriction indicators (R x M) and their date sets, for grouping
    IND = {};          % each: M x 1 logical
    DSET = {};         % each: vector of rows (dates) the restriction touches
    if ~ok_setup
        IND{end+1} = false(M, 1);  DSET{end+1} = 1:nW;
    end

    % ---------- point restrictions ----------
    if has_pt && ok_setup
        for i = 1:length(NR_point.time_index)
            ind = true(M, 1);
            k_i = row_of(NR_point.time_index(i));
            r_i = resolve_regime(NR_point, i, nlambda);
            if ~isfinite(r_i) || r_i < 1 || r_i > max(nlambda, 1) || r_i ~= r_uniq(k_i)
                IND{end+1} = false(M, 1); DSET{end+1} = k_i; continue %#ok<AGROW>
            end
            Vk = V{k_i};
            s_ = NR_point.shock_index(i);
            sg = NR_point.sign(i);
            if sg > 0
                ind = ind & (Vk(:, s_) > 0);
            elseif sg < 0
                ind = ind & (Vk(:, s_) < 0);
            end
            type_i = 'dominance';
            if isfield(NR_point, 'type') && length(NR_point.type) >= i
                type_i = NR_point.type{i};
            end
            if strcmpi(type_i, 'sign')
                IND{end+1} = ind; DSET{end+1} = k_i; continue %#ok<AGROW>
            end
            if r_i > numel(ir_cell) || isempty(ir_cell{r_i})
                IND{end+1} = false(M, 1); DSET{end+1} = k_i; continue %#ok<AGROW>
            end
            v_idx = 1;
            if isfield(NR_point, 'variable_index')
                vi = NR_point.variable_index;
                if length(vi) >= i, v_idx = vi(i); else, v_idx = vi(1); end
            end
            irrow = reshape(ir_cell{r_i}(v_idx, 1, :), 1, []);
            absHD = abs(Vk .* irrow);                      % M x ny
            tgt   = absHD(:, s_);
            oth   = setdiff(1:ny, s_);
            switch lower(type_i)
                case 'dominance_max'
                    ind = ind & (tgt > max(absHD(:, oth), [], 2));
                case 'dominance_min'
                    ind = ind & (tgt < min(absHD(:, oth), [], 2));
                otherwise                                   % legacy strict
                    ind = ind & (tgt > sum(absHD, 2) - tgt);
            end
            IND{end+1} = ind; DSET{end+1} = k_i; %#ok<AGROW>
        end
    end

    % ---------- window-average restrictions ----------
    if has_win && ok_setup
        for e = 1:numel(NRR_reg)
            ind = true(M, 1);
            res = NRR_reg(e);  win = res.window(:)';  r = res.regime_idx;
            if isempty(win) || ~isfinite(r) || r < 1 || r > max(nlambda, 1)
                IND{end+1} = false(M, 1); DSET{end+1} = 1:nW; continue %#ok<AGROW>
            end
            rows = arrayfun(@(t) row_of(t), win);
            has_f = isfield(res, 'shock_idx') && ~isempty(res.shock_idx);
            has_a = isfield(res, 'anti_shock_idx') && ~isempty(res.anti_shock_idx);
            mode  = 'off';
            if has_f && isfield(res, 'dominance') && ~isempty(res.dominance)
                mode = lower(res.dominance);
            end
            do_f = has_f && ~strcmp(mode, 'off');
            do_a = has_a && isfield(res, 'anti_dominance') && strcmpi(res.anti_dominance, 'on');
            if has_f                                        % focal sign on window mean
                j = res.shock_idx;
                mv = zeros(M, 1);
                for k = rows, mv = mv + V{k}(:, j); end
                ind = ind & (res.sign * (mv / numel(rows)) > 0);
            end
            if do_f || do_a
                if r > numel(ir_cell) || isempty(ir_cell{r})
                    IND{end+1} = false(M, 1); DSET{end+1} = rows; continue %#ok<AGROW>
                end
                v_idx = res.variable_idx;
                irrow = reshape(ir_cell{r}(v_idx, 1, :), 1, []);
                mHD = zeros(M, ny);
                for k = rows, mHD = mHD + abs(V{k} .* irrow); end
                mHD = mHD / numel(rows);
                if do_f
                    j = res.shock_idx;
                    if isfield(res, 'other_idx') && ~isempty(res.other_idx), oth = res.other_idx;
                    else, oth = setdiff(1:ny, j); end
                    switch mode
                        case {'on', 'max'}, ind = ind & (mHD(:, j) > max(mHD(:, oth), [], 2));
                        case 'min',         ind = ind & (mHD(:, j) < min(mHD(:, oth), [], 2));
                        case 'strict',      ind = ind & (mHD(:, j) > sum(mHD(:, oth), 2));
                        otherwise, error('estimate_nrr_omega_fast: dominance mode "%s".', mode);
                    end
                end
                if do_a
                    a = res.anti_shock_idx;
                    if isfield(res, 'anti_other_idx') && ~isempty(res.anti_other_idx), ao = res.anti_other_idx;
                    else, ao = setdiff(1:ny, a); end
                    ind = ind & (mHD(:, a) < min(mHD(:, ao), [], 2));
                end
            end
            IND{end+1} = ind; DSET{end+1} = rows; %#ok<AGROW>
        end
    end

    % ---------- joint indicator (validation + non-product fallback) ----------
    R = numel(IND);
    pass = true(M, 1);
    for ri = 1:R, pass = pass & IND{ri}; end

    % ---------- validation gate: replay through the EXACT checkers ----------
    if validate && ~isempty(pass_fns_exact)
        K = min(get_field(opts, 'validate_K', 200), M);
        xi_sim = NaN(Tpl, ny);
        for m = 1:K
            xi_sim(t_uniq, :) = XI(:, :, m);
            ok = true;
            for c = 1:numel(pass_fns_exact)
                if ~pass_fns_exact{c}(xi_sim), ok = false; break; end
            end
            if ok ~= pass(m)
                error(['estimate_nrr_omega_fast: VALIDATION FAILED at sim %d ' ...
                       '(fast = %d, exact = %d). Falling back to the exact ' ...
                       'path is required; do not trust fast omega.'], m, pass(m), ok);
            end
        end
        om_diag.validated = true;
        om_diag.validate_K = K;
    end

    % ---------- grouped omega ----------
    if prodform && R > 1
        % union-find over shared dates -> independent groups
        gid = 1:R;
        for a = 1:R
            for b = a+1:R
                if ~isempty(intersect(DSET{a}, DSET{b}))
                    gid(gid == gid(b)) = gid(a);
                end
            end
        end
        gs = unique(gid);
        omega_hat = 1;  gvec = zeros(1, numel(gs));  gfloor = false;
        for gi = 1:numel(gs)
            ig = true(M, 1);
            for ri = find(gid == gs(gi)), ig = ig & IND{ri}; end
            pg = mean(ig);
            if pg < 1 / M, pg = 1 / M; gfloor = true; end
            gvec(gi) = pg;  omega_hat = omega_hat * pg;
        end
        om_diag.group_omegas = gvec;
        om_diag.n_groups     = numel(gs);
        om_diag.floored      = gfloor;
        om_diag.n_pass       = sum(pass);
    else
        om_diag.n_pass = sum(pass);
        omega_hat = om_diag.n_pass / M;
        if omega_hat < 1 / M, omega_hat = 1 / M; om_diag.floored = true; end
    end

    if momchk                                              % same shape check as v3
        lamW = sqrt(max(cache.lambda(:, r_uniq), 0))';
        xi_real = var.udraw_true(t_uniq, :) ./ lamW;
        om_diag.xi_std      = std(xi_real, 0, 1);
        om_diag.xi_kurtosis = kurtosis(xi_real, 0, 1);
        om_diag.n_rows      = nW;
    end
end

function r = resolve_regime(NR, i, nlambda)
    if nlambda == 1, r = 1; return; end
    if isfield(NR, 'regime_indices') && ~isempty(NR.regime_indices) && length(NR.regime_indices) >= i
        r = NR.regime_indices(i);
    elseif isfield(NR, 'regime_index') && ~isempty(NR.regime_index)
        r = NR.regime_index;
    else
        r = 1;
    end
end

function fs = get_fs(r, fsign, fsign_cell, nlambda)
    if nlambda > 1
        if r >= 1 && r <= numel(fsign_cell), fs = fsign_cell{r}; else, fs = []; end
    else
        fs = fsign;
    end
end

function v = get_field(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end