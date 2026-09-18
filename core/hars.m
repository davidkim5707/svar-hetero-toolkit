%==========================================================================
% hars.m -- HARS sampler for the monetary application (Kim & Zha).
%   Sign + narrative + regime-heteroskedasticity identification. This build
%   has no zero restrictions and no time-varying A0 (A0 is common across
%   regimes; options.ZeroRestrictions is ignored).
%
%   Hybrid Gibbs-MH-ESS sampler with:
%     1. Precomputed theta-cache (removes repeated stability/A0inv calls)
%     2. Pre-parsed sign restrictions (no regex in the hot loop)
%     3. Regime/window-average narrative restrictions (NRR) threaded through
%        the cache via check_narrative_regime_avg.
%
%   N_ess Q-RETENTION (Sec 3 steps 3-4). After the structural move, run N_ess
%   local ESS updates of Q under the final theta and retain all of them as
%   posterior draws (serially dependent, but valid draws of Q | theta, SR>0).
%   The storage stage loops options.ess_store_rounds = N_ess times, storing each
%   resulting (theta, Q); this amortises the structural move over several cheap
%   Q-draws. ness_store == 1 reproduces one-draw-per-sweep. A start-of-sweep
%   pre-mixing loop (premix_m, Sec 3) helps the move land in an admissible Q region.
%
%   JOINT (A0,Lambda,A+) STRUCTURAL REPAIR (Sec 3 step 1, Proposition 1). When an
%   MH-accepted completed draw (A0,Lambda,A+) fails SR(theta',Q_current), the
%   repair redraws the whole structural collection jointly. Holding the accepted
%   (A0,Lambda) fixed and redrawing only A+ is invalid: it targets p0(z0|s0)
%   instead of the sign-restricted marginal p0(z0|s0)*omega0(z0;s0,Q), where
%   omega0(z0;s0,Q) = INT p+(A+|z0,s0) 1{SR(z0,A+,Q)>0} dA+. The joint procedure
%   restores that omega0 factor. Attempt up to options.theta_redraw_budget =
%   N_theta joint redraws from the current state x0 (x' = x0 + N(0,Sigma)), each
%   completed with a fresh A+ ~ p+(A+|z0_cand,s0) and accepted on the A+-collapsed
%   MH ratio times the completed-draw indicator; keep the first sign-admissible
%   candidate. On exhaustion the global Q-repair runs under the accepted theta'
%   (Sec 3 step 2). Mirrored in burn-in.
%
%   Impact matrix. With S_m = A0^{-1} Lambda_m^{1/2} Q', a joint (A0,Lambda) move
%   changes Lambda_m^{1/2} and the relative weights inside S_m, so sign-
%   admissibility is re-checked after every accepted (A0,Lambda) move even though
%   SR is not a restriction on Lambda itself. The average-variance (BPSS)
%   normalization of Lambda is a positive scale and does not flip signs.
%
%   HAAR-ONLY GLOBAL Q-REPAIR (Algorithm 1 Step 2b). The production global
%   Q-repair uses the first admissible independent Haar draw, which has the
%   restricted-Haar conditional distribution. If the N_Q Haar attempts are
%   exhausted the event is a discarded transition (the sampler holds at its
%   current admissible state). The columnwise construction is used only where the
%   paper licenses it: initialization (ARRS Algorithm 4) and the pre-loop
%   seed-rescue block.
%
%   NARRATIVE OMEGA TOGGLE. options.use_nrr_omega gates the AD-RR omega
%   importance weighting for narrative restrictions. With use_nrr_omega = false
%   there is no omega_hat simulation, draw_weight is not divided by omega_hat, and
%   the downstream medians/quantiles are unweighted, so the retained draws target
%   the unweighted restricted-Haar posterior of Definition 1: narrative
%   restrictions enter only through the SR indicator that truncates the
%   conditional distribution of Q, with no reweighting of the retained (theta,Q)
%   pairs. use_nrr_omega = true recovers the AD-RR-style weighted estimand.
%
%   Required helpers:
%     eval_posterior.m / bvar_posterior.m, build_theta_cache.m /
%     check_sr_pass_cached.m, find_admissible_Q_columnwise.m /
%     SignRestrictionCheck.m, check_narrative_regime_avg.m /
%     parse_sign_restrictions.m, checkrestrictions_per_shock_fast.m
%==========================================================================

%% ================================================================
%%  SECTION 1 — DATA & OPTIONS
%% ================================================================
y               = data.y;
varnames        = data.varnames;
filtered_data   = data.filtered_data;

ny              = options.ny;
lags            = options.lags;
constant        = options.constant;
ndraw           = options.ndraw;
nburn           = options.nburn;
nsep            = options.nsep;
nobs            = options.nobs;
% non-tvA build: time-varying A0 lives in svar_run_tva_elliptical_final.m
if isfield(options,'tvA') && ~isempty(options.tvA) && options.tvA
    error('This driver does not support tvA. Use svar_run_tva_elliptical_final.m.');
end
tvA             = 0;
noLmd           = options.noLmd;
irf_horizon     = options.irf_horizon;

if isfield(options, 'timetrend') && ~isempty(options.timetrend)
    timetrend = options.timetrend;
else
    timetrend = 0;
end

minn_prior_lambda = options.minn_prior_lambda;
minn_prior_mu     = options.minn_prior_mu;

%% ================================================================
%%  SECTION 2 — S-STAR & REDUCED-FORM ESTIMATION
%% ================================================================
Sstar = compute_Sstar(y, lags);
options.Sstar = Sstar;
options.use_Sstar_for_vprior = true;
fprintf('sqrt(diag(Sstar)): [%s]\n', num2str(sqrt(diag(Sstar))', '%.4f '));

[T, ~] = size(y);

xdata = [];
if constant,       xdata = [xdata, ones(T, 1)]; end
if timetrend >= 1
    trend_lin = (1:T)' / T;
    xdata = [xdata, trend_lin];
end
if timetrend >= 2,  xdata = [xdata, trend_lin.^2]; end
if isfield(options, 'exogenous') && ~isempty(options.exogenous)
    xdata = [xdata, options.exogenous];
end
if timetrend > 0
    fprintf('Time trend order: %d (columns in xdata: %d)\n', timetrend, size(xdata, 2));
end

seedrfmodel = rfvar3(y, lags, xdata, [], minn_prior_lambda, minn_prior_mu, [], []);

A0_restriction = options.A0_restriction;
if isempty(A0_restriction)
    lcA0 = tril(true(ny));
else
    lcA0 = A0_restriction ~= 0;
end

tempAinv = chol((seedrfmodel.u' * seedrfmodel.u) / size(seedrfmodel.u, 1))';
tempA    = inv(tempAinv);
seedA    = tempA;
seedA(~lcA0) = 0;
seedA(triu(true(ny), 1)) = 0;

%% ================================================================
%%  SECTION 3 — REGIMES
%% ================================================================
regimes  = options.regimes;
nBreaks  = length(regimes);
breakInd = [];
Tsigbrk  = 0;
tparam   = [];  tscale = [];  alpha = [];  K = [];

if ~isempty(regimes)
    breakInd = NaN(1, nBreaks);
    for i = 1:nBreaks
        idx = find(filtered_data.dates == regimes(i));
        if ~isempty(idx), breakInd(i) = idx; end
    end
    Tsigbrk = [0, breakInd];
    if isfield(options, 'tparam') && ~isempty(options.tparam)
        tparam = options.tparam;
        tscale = options.tscale;
    else
        if isfield(options, 'alpha'), alpha = options.alpha; end
        if isfield(options, 'K'),     K     = options.K;     end
    end
end

regime_vec = [0, breakInd - lags, nobs];
nRegimes   = length(regime_vec) - 1;

seedLmd = ones(ny, nRegimes);
lcLmd   = true(ny, nRegimes);
lcLmd(:, nRegimes) = false;

if isfield(options, 'fix_first_regime') && options.fix_first_regime == 1
    lcLmd(:, 1) = false;
    fix_first_regime = options.fix_first_regime;
else
    fix_first_regime = [];
end

%  Homoskedastic shocks: lambda = 1 across all regimes. ------------
%   applyFreeRegimes (in bvar_posterior) normalizes sum_s lambda_{i,s} = nSig
%   (row mean 1), so the ONLY constant-across-regimes value is 1. Dropping the
%   row's free cells from lcLmd leaves the ones-init in seedLmd / constructLambda
%   in place, so the whole row reconstructs to a flat row of 1s (positive ->
%   Sigma_r stays full rank in the ZLB regimes where EFFR is flat). The shock
%   stays identified by sign + narrative; the hetero engine keeps running on
%   the remaining shocks (where the QE/QT identification lives).
%   NOTE: build_theta_cache.local_extractLambda and SignRestrictionCheck.extractLambda
%   carry the matching [HOMOSK FIX] so cache / storage / posterior all agree.
homoskedastic_shocks = [];
if isfield(options,'homoskedastic_shocks') && ~isempty(options.homoskedastic_shocks)
    homoskedastic_shocks = options.homoskedastic_shocks(:)';
    lcLmd(homoskedastic_shocks, :)   = false;   % drop free lambda for these shocks
    seedLmd(homoskedastic_shocks, :) = 1;       % constant relative variance = 1
    fprintf('  Homoskedastic shocks (lambda = 1 across regimes): %s\n', ...
        mat2str(homoskedastic_shocks));
end
% --------------------------------------------------------------------------

%% ================================================================
%%  SECTION 4 — POSTERIOR MODE & HESSIAN
%% ================================================================
seedA_vec   = seedA(lcA0);
seedLmd_vec = seedLmd(lcLmd);

% --- non-tvA build: A0 is COMMON across regimes. lcA0_tv is kept only as the
%     all-false mask that build_theta_cache / SignRestrictionCheck /
%     check_bf_identification expect as an argument. ---
lcA0_tv = false(size(lcA0));

seedx    = [seedA_vec; seedLmd_vec];
nAparams = sum(lcA0(:));

nLmdparams = sum(lcLmd(:));
diag_elements = [ones(nAparams, 1); 1e-5 * ones(nLmdparams, 1)];
seedH = diag(diag_elements);

crit = 1e-10;
nit  = options.nit;
max_compute = 1;

% --- mode-finding objective routed through eval_posterior (bvar_posterior). ---
bvar_obj = @(seedx) eval_posterior(seedx, y, lags, lcA0, lcLmd, Tsigbrk, [], tvA, options);

switch max_compute
    case 1
        [~, xh, ~, H, ~, ~, retcodeh] = csminwel(bvar_obj, seedx, seedH, [], crit, nit, 2);
        if retcodeh ~= 1
            [~, xh1, ~, ~, ~, ~, ~] = csminwel(bvar_obj, xh, seedH, [], crit, nit, 2);
            if any((xh - xh1) < 1e-12), xh = xh1; end
        end
    case 2
        optim_options = optimset('display','iter','MaxFunEvals',1e6, ...
            'TolFun',1e-4,'TolX',1e-4,'Algorithm','interior-point');
        [xh, ~, ~, ~, ~, ~, H] = fmincon(bvar_obj, seedx, ...
            constraint_A, constraint_b, [], [], [], [], [], optim_options);
end

[~, p] = chol(H);
if p > 0
    [V, D] = eig((H + H') / 2);
    D(D < 1e-10) = 1e-10;
    H = V * D * V';
end

% Bacchiocchi-Fanelli (2015) Prop 1 / Rothenberg rank condition for the
% {A0_s, lambda_s} block at the posterior mode xh. The condition is testable
% whenever lambda is regime-specific (heteroskedastic ID), regardless of tvA:
% with common A0 the regime covariances Sigma_s = A0^{-1} diag(lambda_s) A0^{-T}
% still vary through lambda, which is exactly what point-identifies {A0,lambda}.
if nRegimes > 1 && any(lcLmd(:))
    % non-tvA: lcA0_tv is already false(size(lcA0)) from SECTION 4, so the
    % unpacker treats all of A0 as the common block. Pass it through as-is.
    bf_check = check_bf_identification(xh, lcA0, lcA0_tv, lcLmd, ny, nRegimes);
    output.bf_check = bf_check;   % stash for the appendix table

    if ~bf_check.full_rank
        warning(['[B&F] rank(J)=%d < a=%d at the mode: the structural block is ' ...
                 'NOT locally point-identified by the regime covariances; the ' ...
                 'deficient directions lean on sign/narrative + prior ' ...
                 '(see bf_check.null_dirs).'], bf_check.rank, bf_check.dim_theta);
    end
end

%% ================================================================
%%  SECTION 5 — MCMC STARTING POINT
%% ================================================================
x0    = xh + mvnrnd(zeros(length(xh), 1), H, 1)';
Sigma = options.hsnscale * H;

if any(lcLmd(:))
    draw_dout = zeros(nobs, ny, 1000);
    delta0    = draw_dout(:,:,1);
else
    delta0 = [];
end

%% ================================================================
%%  SECTION 6 — SIGN / NARRATIVE SETUP
%% ================================================================
use_signrestrictions = false;

has_sign        = isfield(options,'SignRestrictions')             && ~isempty(options.SignRestrictions);
has_narr        = isfield(options,'NarrativeRestrictions')        && ~isempty(options.NarrativeRestrictions);
has_narr_regime = isfield(options,'NarrativeRegimeRestrictions')  && ~isempty(options.NarrativeRegimeRestrictions);

SignRestrictions            = [];
NarrativeRestrictions       = [];
NarrativeRegimeRestrictions = [];
accepted_shock              = [];
SR_parsed                   = [];
SR_parsed_byregime          = [];

if has_sign || has_narr || has_narr_regime
    use_signrestrictions = true;
    if has_sign, SignRestrictions = options.SignRestrictions; end

    if has_narr
        NarrativeRestrictions = options.NarrativeRestrictions;
        if isfield(NarrativeRestrictions,'time_index') && ~isempty(NarrativeRestrictions.time_index)
            t_narr_array = NarrativeRestrictions.time_index;
            if iscolumn(t_narr_array), t_narr_array = t_narr_array'; end
            regime_indices = zeros(size(t_narr_array));
            for nar_idx = 1:length(t_narr_array)
                t_narr = t_narr_array(nar_idx) + options.lags;
                regime_found = false;
                for r = 1:nRegimes
                    t_start = regime_vec(r) + 1;
                    t_end   = regime_vec(r+1);
                    if t_narr >= t_start && t_narr <= t_end
                        regime_indices(nar_idx) = r;
                        regime_found = true;
                        break
                    end
                end
                if ~regime_found, regime_indices(nar_idx) = 1; end
            end
            NarrativeRestrictions.regime_indices = regime_indices;
            NarrativeRestrictions.regime_index   = regime_indices(1);
        end
    end

    if has_narr_regime
        NarrativeRegimeRestrictions = options.NarrativeRegimeRestrictions;
        fprintf('NRR (regime-avg): %d window(s)\n', numel(NarrativeRegimeRestrictions));
        for kk = 1:numel(NarrativeRegimeRestrictions)
            r = NarrativeRegimeRestrictions(kk);
            fprintf('  NRR(%d): regime %d, shock %d, var %d, sign %+d, |win|=%d, dominance=%s\n', ...
                kk, r.regime_idx, r.shock_idx, r.variable_idx, r.sign, ...
                numel(r.window), r.dominance);
        end
    end

    sign_regime_dependent  = options.sign_regime_dependent;
    penalty_offdiagonal_on = options.penalty_offdiagonal_on;

    if has_sign
        shock_token_cells = cellfun(@(s) regexp(s, ',\s*(\d+)\)', 'tokens', 'once'), ...
            SignRestrictions, 'UniformOutput', false);
        shock_ids = cellfun(@(c) sscanf(c{1}, '%d'), shock_token_cells);
    elseif has_narr
        shock_ids = NarrativeRestrictions.shock_index;
    elseif has_narr_regime
        shock_ids = [NarrativeRegimeRestrictions.shock_idx];
    else
        shock_ids = [];
    end

    accepted_shock = unique(shock_ids);

    if has_sign
        SR_parsed = parse_sign_restrictions(SignRestrictions, ny);
        fprintf('SR parsed: %d sign + %d elasticity restrictions across %d shocks\n', ...
            sum(SR_parsed.n_sign), sum(SR_parsed.n_elast), ny);

        % --- per-regime gate: SR_parsed_byregime (IRF check) +
        %     impact_signs_byregime (fast impact screen) -----------------------
        if isfield(options,'sign_active_regimes') && ~isempty(options.sign_active_regimes)
            gate = options.sign_active_regimes;
            % per-restriction shock id from the STRINGS (shock_ids above is
            % contaminated by zero_shock_ids -> do NOT reuse it here)
            restr_shock = zeros(1, numel(SignRestrictions));
            for q = 1:numel(SignRestrictions)
                tk = regexp(SignRestrictions{q}, ',\s*(\d+)\)', 'tokens', 'once');
                restr_shock(q) = str2double(tk{1});
            end
            IS_full = options.impact_signs;        % ny x ny, columns = shocks
            ks      = cell2mat(keys(gate));        % gated shock ids
            SR_parsed_byregime    = cell(1, nRegimes);
            impact_signs_byregime = cell(1, nRegimes);
            for r = 1:nRegimes
                % (a) IRF-check restriction subset
                keep = true(1, numel(SignRestrictions));
                for q = 1:numel(SignRestrictions)
                    s = restr_shock(q);
                    if isKey(gate, s), keep(q) = ismember(r, gate(s)); end
                end
                SR_parsed_byregime{r} = parse_sign_restrictions(SignRestrictions(keep), ny);
                % (b) impact-screen matrix: zero the gated-out shock columns
                ISr = IS_full;
                for s = ks
                    if ~ismember(r, gate(s)), ISr(:, s) = 0; end
                end
                impact_signs_byregime{r} = ISr;
                fprintf('  [GATE] regime %d: %d/%d sign restrictions active\n', ...
                    r, sum(keep), numel(SignRestrictions));
            end
        end
    end

    % Max horizon needed for sign + narrative checks
    if has_sign
        sr_max_horizon = 0;
        for i_sr = 1:length(SignRestrictions)
            tokens = regexp(SignRestrictions{i_sr}, 'y\(\d+,\s*\d+:(\d+)', 'tokens', 'once');
            if ~isempty(tokens)
                h_end = sscanf(tokens{1}, '%d');
                sr_max_horizon = max(sr_max_horizon, h_end);
            end
        end
    else
        sr_max_horizon = 1;  % narrative impact only
    end

    % [INTERVAL HD] Narrative dominance in 'interval_hd' mode accumulates the
    %   delayed dynamics (Theta_0, Theta_1, ...) of every window shock to the
    %   window endpoint, so cache.Psi / ir_cell must span the whole anchor
    %   window. Raise sr_max_horizon to the longest window length HERE, before
    %   the btc / check_fn anonymous functions capture sr_max_horizon by value.
    if ~isempty(NarrativeRegimeRestrictions) ...
            && isfield(options, 'nrr_dom_aggregation') ...
            && strcmpi(options.nrr_dom_aggregation, 'interval_hd')
        maxwin = 0;
        for kk = 1:numel(NarrativeRegimeRestrictions)
            w = NarrativeRegimeRestrictions(kk).window;
            if ~isempty(w)
                maxwin = max(maxwin, max(w(:)) - min(w(:)) + 1);
            end
        end
        if maxwin > sr_max_horizon
            fprintf('  [INTERVAL HD] interval-HD dominance -> cache horizon %d -> %d (max window length)\n', ...
                sr_max_horizon, maxwin);
            sr_max_horizon = maxwin;
        end
    end

    fprintf('  Sign+narrative max horizon: %d (vs irf_horizon = %d)\n', ...
        sr_max_horizon, irf_horizon);
end

%% ================================================================
%%  SECTION 7 — ESS TUNING PARAMETERS
%% ================================================================
ess_max_shrinks_per_col = 50;
 
if isfield(options, 'ess_premix_rounds') && ~isempty(options.ess_premix_rounds)
    premix_m = options.ess_premix_rounds;
else
    premix_m = 3;
end
 
% --- N_ess: number of local ESS rotations RETAINED per structural draw -----
% Sec 3 step 3 / step 4: after the structural move, run N_ess local
% ESS updates of Q under the final theta and retain ALL of them as posterior
% draws (serially dependent, but valid draws of Q | theta, SR>0). Sharing one
% theta across N_ess rotations amortises the expensive structural move -- this
% is the efficiency gain. Q <- Q^(N_ess) carries to the next structural update.
% ness_store == 1 reproduces one stored draw per structural sweep. Default
% premix_m.
if isfield(options, 'ess_store_rounds') && ~isempty(options.ess_store_rounds)
    ness_store = options.ess_store_rounds;
else
    ness_store = premix_m;
end
 
if isfield(options, 'ess_repair_budget') && ~isempty(options.ess_repair_budget)
    repair_k = options.ess_repair_budget;
else
    repair_k = 2000;
end
 
if isfield(options, 'impact_signs') && ~isempty(options.impact_signs)
    impact_signs = options.impact_signs;
else
    impact_signs = zeros(ny, ny);
end
 
% --- Upstream stability check toggle (application-aware) ---
if isfield(options, 'upstream_stability_check') && ~isempty(options.upstream_stability_check)
    upstream_sc = options.upstream_stability_check;
else
    upstream_sc = false;   % default OFF (oil-safe)
end
% Upstream rejection buffer. tol = 1e-6 reproduces the legacy threshold;
% tol = 0 matches SignRestrictionCheck.checkStability exactly (reject iff
% any companion eigenvalue has modulus >= 1).
stability_tol = get_opt(options, 'upstream_stability_tol', 1e-6);
 
% --- Verification stability gate toggle ---
% enforce_stability=false disables the checkStability rejection inside
% SignRestrictionCheck (Carriero 2024 / Brunnermeier 2021 style: flat/Jeffreys
% reduced-form prior, uniform Q, no stationarity conditioning). Default true
% preserves legacy behaviour. This is the switch that removes the ~99% "Skip".
enforce_stability = get_opt(options, 'enforce_stability', true);
 
% --- Seed-rescue A+ budget (INITIALIZATION ONLY) ---
% Max A+ redraws per rescue round in the pre-loop SEED block of SECTION 10.
% This budget affects only the search for an admissible starting state and
% has no bearing on the stationary distribution.
B_redraw_budget = get_opt(options, 'B_redraw_budget', 100);
 
% --- N_theta: JOINT (A0,Lambda,A+) structural-repair budget ---
%  Sec 3 step 1 / Proposition 1. Max JOINT redraws of the whole
% structural collection to restore admissibility of the current Q after an
% MCMC-ACCEPTED completed draw strands it. Independent of the seed-rescue
% budget above (joint structural repair vs. a fixed-z0 A+ Gibbs update).
% Default 100; set options.theta_redraw_budget to override.
theta_redraw_budget = get_opt(options, 'theta_redraw_budget', 100);
 
fprintf('\n=== MH Sampler with Q-block ESS (Cached + Fast SR + NRR) ===\n');
fprintf('  ESS mode                : column-by-column + theta-cache\n');
fprintf('  max shrinks per column  : %d\n', ess_max_shrinks_per_col);
fprintf('  pre-mix rounds (m)      : %d\n', premix_m);
fprintf('  N_ess retained per draw : %d (Sec 3 step 3)\n', ness_store);
fprintf('  global repair budget(NQ): %d\n', repair_k);
fprintf('  global repair mode      : Haar-only (columnwise fallback REMOVED)\n');
fprintf('  seed-rescue A+ budget   : %d (initialization only)\n', B_redraw_budget);
fprintf('  N_theta (joint repair)  : %d (A0,Lambda,A+ joint redraw, Prop 1)\n', theta_redraw_budget);
fprintf('  nsep (thinning)         : %d\n', nsep);
fprintf('  SR fast check           : %s\n', ternary(~isempty(SR_parsed), 'ON', 'OFF'));
fprintf('  NRR regime-avg          : %s\n', ternary(~isempty(NarrativeRegimeRestrictions), 'ON', 'OFF'));
fprintf('  upstream stability      : %s (tol = %.1e)\n', ternary(upstream_sc, 'ON', 'OFF'), stability_tol);
fprintf('  verification stability  : %s\n', ternary(enforce_stability, 'ON (reject unstable)', 'OFF (keep all draws)'));
fprintf('========================================================\n');

%  Master switch for the AD-RR omega importance weighting.
%   options.use_nrr_omega = false -> NO omega_hat estimation and NO 1/omega
%   weighting anywhere: the retained draws target the UNWEIGHTED
%   restricted-Haar posterior of Definition 1 (narrative enters only through
%   the SR indicator). Default true preserves the legacy AD-RR-weighted
%   behaviour. Affects: use_draw_weights below, the storage-stage omega block
%   (SECTION 10), and the SECTION 11 omega diagnostics.
use_nrr_omega = get_opt(options, 'use_nrr_omega', true);

%  AD-RR omega settings. The checkers standardize
% udraw_true by regime internally and rotate xi = Lambda^{-1/2} eps, so
% the simulator draws xi directly on the unit scale: no lambda option,
% and omega is invariant to any common scale on xi (shape only).
nrr_omega_M    = get_opt(options, 'nrr_omega_M', 2000);
nrr_omega_opts = struct( ...
    'M',              nrr_omega_M, ...
    'shock_law',      get_opt(options, 'nrr_omega_law', 'student'), ...
    'nu',             get_opt(options, 'nrr_omega_nu', 2 * options.tparam), ...
    'moment_check',   false);
fprintf('  NRR omega weight (AD-RR) : %s (M = %d, law = %s, xi-scale)\n', ...
    ternary(use_nrr_omega && (has_narr_regime || has_narr), 'ON', ...
        ternary(has_narr_regime || has_narr, 'OFF (use_nrr_omega = false)', 'OFF (no narrative)')), ...
    nrr_omega_M, nrr_omega_opts.shock_law);
 
%  single flag for ALL downstream weighting gates (SECTIONS 11, 13, 14,
% and the output.wquantile / output.wmedian handles).
%  weights fire ONLY when use_nrr_omega is on AND narrative
% restrictions are present; otherwise every downstream summary is unweighted.
use_draw_weights = use_nrr_omega && (has_narr_regime || has_narr);

% --- Cache builder wrapper: captures the fixed config so the ~9
%     build_theta_cache call sites need only pass (x, var). lcA0_tv is the
%     all-false mask (non-tvA). ---
btc = @(xx, vv) build_theta_cache(xx, vv, lcA0, lcLmd, ny, lags, constant, ...
    fix_first_regime, accepted_shock, impact_signs, SR_parsed, upstream_sc, ...
    sr_max_horizon, NarrativeRegimeRestrictions, stability_tol, lcA0_tv, ...
    SR_parsed_byregime, impact_signs_byregime);
 
%% ================================================================
%%  SECTION 8 — INITIALIZE Q
%% ================================================================
if use_signrestrictions
    fprintf('\n=== Initializing admissible Q (cached) ===\n');

    delta_input_init = [];
    if any(lcLmd(:)), delta_input_init = delta0; end

    [lh0, likelihood0, A0_prior0, lambda_prior0, var0] = ...
        eval_posterior(x0, y, lags, lcA0, lcLmd, Tsigbrk, delta_input_init, tvA, options);

    [cache0, stable0] = btc(x0, var0);

    if ~stable0
        fprintf('  Initial x0 unstable — searching for stable theta...\n');
        for i_retry = 1:100
            x0 = xh + mvnrnd(zeros(length(xh), 1), H, 1)';
            [lh0, likelihood0, A0_prior0, lambda_prior0, var0] = ...
                eval_posterior(x0, y, lags, lcA0, lcLmd, Tsigbrk, delta_input_init, tvA, options);
            [cache0, stable0] = btc(x0, var0);
            if stable0, break; end
        end
        if ~stable0
            error('Cannot find a stable initial theta.');
        end
    end

    if isfield(cache0,'NarrativeRegimeRestrictions') && ~isempty(cache0.NarrativeRegimeRestrictions)
        fprintf('  [NRR check] cache carries %d window(s)\n', numel(cache0.NarrativeRegimeRestrictions));
    else
        fprintf('  [NRR check] *** cache has NO NRR — not enforced in Q-search ***\n');
    end

    check_fn_init = @(Q) check_sr_pass_cached(Q, cache0, var0, ...
        SignRestrictions, NarrativeRestrictions, sr_max_horizon, ...
        sign_regime_dependent, penalty_offdiagonal_on, ny);

    % ---- Fallback 1.5 toggle (opt-in; default OFF) -----------------
    if isfield(options, 'use_local_pert_fallback') && ~isempty(options.use_local_pert_fallback)
        use_local_pert_fallback = options.use_local_pert_fallback;
    else
        use_local_pert_fallback = false;
    end

    % ---- Search budget knobs (overridable via options; defaults below) ----
    cw_outer          = get_opt(options, 'cw_outer',          2000);
    cw_inner          = get_opt(options, 'cw_inner',          1000);
    haar_attempts     = get_opt(options, 'haar_attempts',     20000);
    haar_budget       = get_opt(options, 'haar_budget',       2000);
    local_pert_scale  = get_opt(options, 'local_pert_scale',  0.1);
    local_pert_tries  = get_opt(options, 'local_pert_tries',  500);
    cw_outer_local    = get_opt(options, 'cw_outer_local',    5000);
    cw_inner_local    = get_opt(options, 'cw_inner_local',    10000);
    report_every      = get_opt(options, 'report_every',      10);
    max_theta_retries = get_opt(options, 'max_theta_retries', 2000);
    % ----------------------------------------------------------------------

    fprintf('  Fallback 1.5 (local theta pert) : %s\n', ternary(use_local_pert_fallback, 'ON', 'OFF'));

    % =====================================================================
    %  [3] PRIMARY SEARCH — columnwise at current theta
    %      (INITIALIZATION ONLY: ARRS Algorithm 4; licensed by the paper's
    %       initialization step. The production-kernel repair is Haar-only —
    %       see in SECTIONS 9/10.)
    % =====================================================================
    [X_current, Q_current, init_success, init_ntried] = ...
        find_admissible_Q_columnwise( ...
            ny, check_fn_init, cache0.A0inv, cache0.lambda, impact_signs, ...
            sign_regime_dependent, cw_outer, cw_inner);

    % =====================================================================
    %  [4] FALLBACK 1 — Haar at same theta
    % =====================================================================
    if ~init_success
        fprintf('  Column-by-column failed. Falling back to brute-force Haar...\n');
        [X_current, Q_current, init_success, ~] = ...
            find_admissible_Q(check_fn_init, ny, haar_attempts, haar_budget);
    end

    % =====================================================================
    %  [4.5] FALLBACK 1.5 — local theta perturbation + columnwise (OPT-IN)
    % =====================================================================
    if ~init_success && use_local_pert_fallback
        fprintf('  Haar failed. Trying %d local theta perturbations + columnwise...\n', local_pert_tries);
        pert_start_time     = tic;
        n_unstable          = 0;
        n_columnwise_failed = 0;

        for i_pert = 1:local_pert_tries
            iter_start = tic;

            x0_pert = x0 + mvnrnd(zeros(length(x0), 1), local_pert_scale * Sigma)';

            [~, ~, ~, ~, var_pert] = ...
                eval_posterior(x0_pert, y, lags, lcA0, lcLmd, Tsigbrk, ...
                    delta_input_init, tvA, options);

            [cache_pert, stable_pert] = btc(x0_pert, var_pert);

            if ~stable_pert
                n_unstable = n_unstable + 1;
                if mod(i_pert, report_every) == 0
                    fprintf('    [Pert %3d/%d] unstable: %d | cw-fail: %d | total %.1fs\n', ...
                        i_pert, local_pert_tries, n_unstable, n_columnwise_failed, ...
                        toc(pert_start_time));
                end
                continue
            end

            check_fn_pert = @(Q) check_sr_pass_cached(Q, cache_pert, var_pert, ...
                SignRestrictions, NarrativeRestrictions, sr_max_horizon, ...
                sign_regime_dependent, penalty_offdiagonal_on, ny);

            [X_current, Q_current, init_success, n_tried_cw] = ...
                find_admissible_Q_columnwise( ...
                    ny, check_fn_pert, cache_pert.A0inv, cache_pert.lambda, ...
                    impact_signs, sign_regime_dependent, ...
                    cw_outer_local, cw_inner_local);

            iter_time = toc(iter_start);

            if init_success
                x0 = x0_pert; var0 = var_pert; cache0 = cache_pert;
                [lh0, likelihood0, A0_prior0, lambda_prior0, ~] = ...
                    eval_posterior(x0, y, lags, lcA0, lcLmd, Tsigbrk, ...
                        delta_input_init, tvA, options);
                fprintf('    Local perturbation %d: SUCCESS (cw tried %d, %.1fs)\n', ...
                    i_pert, n_tried_cw, iter_time);
                break
            else
                n_columnwise_failed = n_columnwise_failed + 1;
            end

            if mod(i_pert, report_every) == 0 || i_pert == local_pert_tries
                elapsed = toc(pert_start_time);
                fprintf('    [Pert %3d/%d] unstable: %d | cw-fail: %d | last iter %.1fs | avg %.1fs/iter | total %.1fs\n', ...
                    i_pert, local_pert_tries, n_unstable, n_columnwise_failed, ...
                    iter_time, elapsed / i_pert, elapsed);
            end
        end
    end

    % =====================================================================
    %  [5] FALLBACK 2 — alternative theta x max_theta_retries
    % =====================================================================
    if ~init_success
        fprintf('  Trying %d alternative starting theta...\n', max_theta_retries);
        for i_retry = 1:max_theta_retries
            x0_retry = xh + mvnrnd(zeros(length(xh), 1), H, 1)';
            [~, ~, ~, ~, var0_retry] = ...
                eval_posterior(x0_retry, y, lags, lcA0, lcLmd, Tsigbrk, delta_input_init, tvA, options);
            [cache_retry, stable_retry] = btc(x0_retry, var0_retry);
            if ~stable_retry, continue; end

            check_fn_retry = @(Q) check_sr_pass_cached(Q, cache_retry, var0_retry, ...
                SignRestrictions, NarrativeRestrictions, sr_max_horizon, ...
                sign_regime_dependent, penalty_offdiagonal_on, ny);

            [X_current, Q_current, init_success, ~] = ...
                find_admissible_Q(check_fn_retry, ny, haar_attempts, haar_budget);

            if init_success
                x0 = x0_retry; var0 = var0_retry; cache0 = cache_retry;
                [lh0, likelihood0, A0_prior0, lambda_prior0, ~] = ...
                    eval_posterior(x0, y, lags, lcA0, lcLmd, Tsigbrk, delta_input_init, tvA, options);
                fprintf('    Retry %d: success!\n', i_retry);
                break
            end
        end
    end

    % =====================================================================
    %  [6] HARD ERROR — all strategies exhausted
    % =====================================================================
    if ~init_success
        error('Could not find an admissible Q. Check restrictions.');
    else
        if exist('init_ntried', 'var') && init_ntried > 0 && init_ntried < cw_outer
            fprintf('  Initial Q found in %d outer attempts (primary).\n', init_ntried);
        else
            fprintf('  Initial Q found via fallback strategy.\n');
        end
    end

    % =====================================================================
    %  [7] ONE-SHOT SANITY CHECK (remove after verification)
    % =====================================================================
    if isfield(cache0, 'Psi') && ~isempty(cache0.Psi)
        smat_test = cache0.A0inv * Q_current' * diag(sqrt(max(cache0.lambda(:,1), 0)));
        ir_old    = impulsdtrf(cache0.By, smat_test, sr_max_horizon);
        ir_new    = zeros(ny, sr_max_horizon, ny);
        for h = 1:sr_max_horizon
            ir_new(:, h, :) = cache0.Psi(:, :, h) * smat_test;
        end
        err = max(abs(ir_old(:) - ir_new(:)));
        fprintf('  [SANITY] Psi vs impulsdtrf max abs diff: %.3e (expect < 1e-10)\n', err);
        assert(err < 1e-8, 'Psi caching FAILED — disable horizon arg or check Psi recursion');
    end

    fprintf('===========================================\n\n');
end

%==========================================================================
%  SECTION 9 (Sec 3 burn-in; mirrors SECTION 10 Step 1)
%
%  The burn-in transition kernel is identical to the production kernel so the
%  chain converges to the same stationary distribution. Per isep-equivalent
%  burn-in step:
%     [pre-mix]  premix_m local ESS on Q at the current admissible (x0,Q)
%     Step 1:    MH propose z0 -> select z0* -> fresh A+ at z0* -> SR check
%                -> stale (accept OR reject): JOINT redraw (Prop 1 /
%                   literal Sec 3 Step 1)
%                -> global Q-repair under the original theta* (shared fallback;
%                    Haar-only)
%     [delta]    update latent scales under hetero
%  Burn-in draws are discarded (no Step 3 retention / storage).
%==========================================================================
%% ================================================================
%%  SECTION 9 — BURN-IN
%% ================================================================
if nburn > 0 && use_signrestrictions
    burnin_start_time = tic;
    has_hetero = any(lcLmd(:));

    acc_burnin_mh               = 0;   % z0 = (A0,Lambda) MH moves accepted
    acc_burnin_sr_direct        = 0;   % completed draw admissible (no repair)
    repair_attempts_burnin      = 0;   % stale-Q events (accept or reject)
    theta_redraws_burnin        = 0;   % JOINT redraw proposals
    theta_repair_success_burnin = 0;   % JOINT redraw successes (Q kept)
    repair_success_burnin       = 0;   % global Q-repair successes
    ess_evals_burnin            = 0;
    sampler_failures_burnin     = 0;

    for iburn = 1:nburn

        % ---- PRE-MIXING at the current admissible (x0, var0, Q) ----
        check_fn = @(Q) check_sr_pass_cached(Q, cache0, var0, ...
            SignRestrictions, NarrativeRestrictions, sr_max_horizon, ...
            sign_regime_dependent, penalty_offdiagonal_on, ny);
        for ipremix = 1:premix_m
            [X_current, Q_current, ess_ne] = ...
                ess_draw_Q_columnwise(X_current, check_fn, ess_max_shrinks_per_col);
            ess_evals_burnin = ess_evals_burnin + ess_ne;
        end

        % ---- STEP 1: update (A0,Lambda,A+) ----
        % (a) baseline collapsed marginal lh0 at the current delta + FRESH A+
        %     at x0 (the reject-branch completion).
        if has_hetero
            [lh0, likelihood0, A0_prior0, lambda_prior0, var0_rej] = ...
                eval_posterior(x0, y, lags, lcA0, lcLmd, Tsigbrk, delta0, tvA, options);
        else
            [lh0, likelihood0, A0_prior0, lambda_prior0, var0_rej] = ...
                eval_posterior(x0, y, lags, lcA0, lcLmd, Tsigbrk, [], tvA, options);
        end

        % (a') propose z0 = (A0,Lambda), complete with a fresh A+ at x1
        x1 = x0 + mvnrnd(zeros(length(x0), 1), Sigma)';
        if has_hetero
            [lh1, likelihood1, A0_prior1, lambda_prior1, var1] = ...
                eval_posterior(x1, y, lags, lcA0, lcLmd, Tsigbrk, delta0, tvA, options);
        else
            [lh1, likelihood1, A0_prior1, lambda_prior1, var1] = ...
                eval_posterior(x1, y, lags, lcA0, lcLmd, Tsigbrk, [], tvA, options);
        end

        mh_accept = ~isempty(var1) && isstruct(var1) && isfield(var1,'Bdraw') ...
                    && (log(rand()) < lh0 - lh1);

        % (b) select z0* and the FRESH A+ completion at z0*
        if mh_accept
            acc_burnin_mh = acc_burnin_mh + 1;
            x_star   = x1;   var_star = var1;   lh_star  = lh1;
            lik_star = likelihood1; A0p_star = A0_prior1; lmdp_star = lambda_prior1;
        else
            x_star   = x0;   var_star = var0_rej; lh_star = lh0;
            lik_star = likelihood0; A0p_star = A0_prior0; lmdp_star = lambda_prior0;
        end

        % (c) SR on the completed draw theta* = (z0*, A+*)
        [cache_star, stable_star] = btc(x_star, var_star);

        if ~stable_star
            sampler_failures_burnin = sampler_failures_burnin + 1;
        else
            sr_pass = check_sr_pass_cached(Q_current, cache_star, var_star, ...
                SignRestrictions, NarrativeRestrictions, sr_max_horizon, ...
                sign_regime_dependent, penalty_offdiagonal_on, ny);

            if sr_pass
                acc_burnin_sr_direct = acc_burnin_sr_direct + 1;
                x0 = x_star; var0 = var_star; lh0 = lh_star;
                likelihood0 = lik_star; A0_prior0 = A0p_star; lambda_prior0 = lmdp_star;
                cache0 = cache_star;
            else
                % ===== JOINT (A0,Lambda,A+) redraw (Prop 1) =====
                % Stale completed draw (accept OR reject): repair by the
                % JOINT redraw of the whole structural collection
                % (Proposition 1 / literal Sec 3 Step 1). Q is kept.
                repaired = false;
                repair_attempts_burnin = repair_attempts_burnin + 1;
                if any(lcA0(:)) || any(lcLmd(:))
                    for t_try = 1:theta_redraw_budget
                        theta_redraws_burnin = theta_redraws_burnin + 1;
                        x_th = x0 + mvnrnd(zeros(length(x0), 1), Sigma)';
                        if has_hetero
                            [lh_th, lik_th, A0p_th, lmdp_th, var_th] = ...
                                eval_posterior(x_th, y, lags, lcA0, lcLmd, Tsigbrk, delta0, tvA, options);
                        else
                            [lh_th, lik_th, A0p_th, lmdp_th, var_th] = ...
                                eval_posterior(x_th, y, lags, lcA0, lcLmd, Tsigbrk, [], tvA, options);
                        end
                        mh_accept_th = ~isempty(var_th) && isstruct(var_th) ...
                            && isfield(var_th,'Bdraw') && (log(rand()) < lh0 - lh_th);
                        if ~mh_accept_th, continue; end
                        [cache_th, stable_th] = btc(x_th, var_th);
                        if ~stable_th, continue; end
                        sr_pass_th = check_sr_pass_cached(Q_current, cache_th, var_th, ...
                            SignRestrictions, NarrativeRestrictions, sr_max_horizon, ...
                            sign_regime_dependent, penalty_offdiagonal_on, ny);
                        if ~sr_pass_th, continue; end
                        x0 = x_th; var0 = var_th; lh0 = lh_th;
                        likelihood0 = lik_th; A0_prior0 = A0p_th; lambda_prior0 = lmdp_th;
                        cache0 = cache_th;
                        theta_repair_success_burnin = theta_repair_success_burnin + 1;
                        repaired = true;
                        break
                    end
                end

                if ~repaired
                    % ===== GLOBAL Q-repair under the ORIGINAL theta* =====
                    %  Paper Algorithm 1 Step 2b: independent Haar
                    % draws only (up to N_Q = repair_k). The columnwise
                    % fallback is REMOVED -- a columnwise Q is not a draw from
                    % restricted Haar (Prop. diagnostics (i)) and would leave
                    % the Q-block conditional target of Definition 1.
                    check_fn_repair = @(Q) check_sr_pass_cached(Q, cache_star, var_star, ...
                        SignRestrictions, NarrativeRestrictions, sr_max_horizon, ...
                        sign_regime_dependent, penalty_offdiagonal_on, ny);
                    [X_repair, Q_repair, repair_ok, ~] = ...
                        find_admissible_Q(check_fn_repair, ny, repair_k, repair_k);
                    if repair_ok
                        x0 = x_star; var0 = var_star; lh0 = lh_star;
                        likelihood0 = lik_star; A0_prior0 = A0p_star; lambda_prior0 = lmdp_star;
                        cache0 = cache_star;
                        X_current = X_repair; Q_current = Q_repair;
                        repair_success_burnin = repair_success_burnin + 1;
                    else
                        sampler_failures_burnin = sampler_failures_burnin + 1;
                    end
                end
            end
        end

        % (f) hetero: update latent scales delta (third Gibbs block).
        if has_hetero
            eout = var0.udraw(1:size(delta0, 1), :) .* exp(delta0);
            if isempty(tparam)
                delta0 = drawdelta(eout, alpha, K);
            else
                delta0 = drawt(eout, tparam, (tscale^2) * tparam);
            end
        end

        if mod(iburn, max(1, floor(nburn/20))) == 0 || iburn == nburn
            fprintf(['[BURN-IN] %4d/%4d | MH:%d | SR-direct:%d | Joint:%d/%d | ' ...
                     'Qrep:%d/%d | Fail:%d | ESS:%d | %.1fs\n'], ...
                iburn, nburn, acc_burnin_mh, acc_burnin_sr_direct, ...
                theta_repair_success_burnin, theta_redraws_burnin, ...
                repair_success_burnin, repair_attempts_burnin, ...
                sampler_failures_burnin, ess_evals_burnin, ...
                toc(burnin_start_time));
        end
    end

    fprintf('\n=== Burn-in Complete ===\n');
    fprintf('  MH acceptances:            %d / %d (%.1f%%)\n', acc_burnin_mh, nburn, 100*acc_burnin_mh/nburn);
    fprintf('  SR direct:                 %d\n', acc_burnin_sr_direct);
    fprintf('  Joint structural-repairs:  %d / %d\n', theta_repair_success_burnin, theta_redraws_burnin);
    fprintf('  Global Q-repair:           %d / %d\n', repair_success_burnin, repair_attempts_burnin);
    fprintf('  Sampler failures:          %d\n', sampler_failures_burnin);
    fprintf('========================\n\n');

elseif nburn > 0 && ~use_signrestrictions
    burnin_start_time = tic;
    acceptance_burnin = 0;
    for iburn = 1:nburn
        for isep = 1:nsep
            delta_input = [];
            if any(lcLmd(:)), delta_input = delta0; end
            [lh0, ~, ~, ~, var0] = eval_posterior(x0, y, lags, lcA0, lcLmd, Tsigbrk, delta_input, tvA, options);
            x1 = x0 + mvnrnd(zeros(length(x0), 1), Sigma)';
            [lh1, ~, ~, ~, var1] = eval_posterior(x1, y, lags, lcA0, lcLmd, Tsigbrk, delta_input, tvA, options);
            if log(rand()) < lh0 - lh1
                x0 = x1;  var0 = var1;  lh0 = lh1;
                acceptance_burnin = acceptance_burnin + 1;
            end
            if any(lcLmd(:))
                eout = var0.udraw(1:size(delta0, 1), :) .* exp(delta0);
                if isempty(tparam), delta0 = drawdelta(eout, alpha, K);
                else,               delta0 = drawt(eout, tparam, (tscale^2)*tparam); end
            end
        end
    end
    fprintf('Burn-in complete: %d MH acc (%.1f%%), %.1f sec\n', ...
        acceptance_burnin, 100*acceptance_burnin/(nburn*nsep), toc(burnin_start_time));
end

%==========================================================================
%  SECTION 10 (Sec 3, literal)
%
%  EVERY stale completed draw (MH-accept OR MH-reject) is repaired by the
%  JOINT (A0,Lambda,A+) redraw -- the as-written Sec 3 Step 1 / Proposition 1
%  procedure.
%
%  The global Q-repair (Haar-only) under the ORIGINAL stale
%  theta* is the final fallback before a discarded transition.
%  B_redraw_budget is the SEED-rescue budget (initialization only);
%  theta_redraw_budget = N_theta is the joint-redraw budget.
%==========================================================================
%% ================================================================
%%  SECTION 10 — MAIN SAMPLING   (instrumented: block tic/toc timers)
%% ================================================================
fprintf("\nSampling process starts\n");
samplein_start_time = tic;

fprintf('  stale-Q repair          : JOINT (A0,Lambda,A+) redraw (Prop 1), then Haar-only Q-repair\n');

% --- BLOCK-LEVEL PROFILE ACCUMULATORS (sampling loop only) ---
prof = struct('t_eval',0,'n_eval',0, ...   % eval_posterior (x0 / x1 / x_th / xrot)
              't_cache',0,'n_cache',0, ...  % build_theta_cache (btc)
              't_ess',0,'n_ess',0, ...      % ess_draw_Q_columnwise (premix + storage)
              't_repair',0,'n_repair',0,... % sr-check + find_admissible_Q (Haar-only)
              't_src',0,'n_src',0, ...      % SignRestrictionCheck (storage)
              't_delta',0,'n_delta',0);     % drawt / drawdelta

target_draws    = ndraw;
collected_draws = 0;
idraw           = 0;
has_hetero      = any(lcLmd(:));

total_mh_proposals    = 0;
total_mh_accepted     = 0;   % z0 = (A0,Lambda) MH moves accepted
total_sr_direct       = 0;   % completed draw admissible under current Q (no repair)
total_repair_attempts = 0;   % stale-Q events (accept or reject)
total_repair_success  = 0;   % stale events resolved by global Q-repair
total_theta_redraws       = 0;   % JOINT (A0,Lambda,A+) structural-repair proposals
total_theta_repair_success = 0;  % JOINT structural-repair successes (Q kept)
total_sampler_failures = 0;
total_ess_evals       = 0;
total_B_redraws       = 0;   % seed-rescue A+ redraws (initialization only)
total_parameter_draws = 0;
total_parameter_draws_discarded = 0;
total_parameter_draws_accepted  = 0;

max_expected = target_draws;

if use_signrestrictions
    sign_irfs_all           = cell(max_expected, 1);
    draw_lh_wQ_all          = zeros(max_expected, 1);
    draw_likelihood_wQ_all  = zeros(max_expected, 1);
    draw_A0prior_wQ_all     = zeros(max_expected, 1);
    draw_lambdaprior_wQ_all = zeros(max_expected, 1);
    draw_Qdraw_all          = cell(max_expected, 1);
    draw_Lambda_all         = cell(max_expected, 1);
    draw_By_all             = cell(max_expected, 1);
    draw_fsign_all          = cell(max_expected, 1);
    draw_fsign_cell_all     = cell(max_expected, 1);
    draw_weight_all         = ones(max_expected, 1);   %  1/omega_hat weight (== 1 when toggle off)
    draw_omega_all          = ones(max_expected, 1);   %  per-draw omega_hat
    om_diag_all             = cell(max_expected, 1);   %  per-draw omega diagnostics
end

draw_lh_all          = zeros(max_expected, 1);
draw_likelihood_all  = zeros(max_expected, 1);
draw_A0prior_all     = zeros(max_expected, 1);
draw_lambdaprior_all = zeros(max_expected, 1);
draw_x0_all          = cell(max_expected, 1);
draw_var_all         = cell(max_expected, 1);
draw_Phi_all         = cell(max_expected, 1);
draw_udraw_all       = cell(max_expected, 1);
draw_udraw_true_all  = cell(max_expected, 1);
if has_hetero
    draw_dout_all = cell(max_expected, 1);
end

% =====================================================================
%  SEED: establish an admissible (x0, var0, cache0, Q_current).
%  Initialization only (NOT part of the production kernel).
%   The columnwise routines below are RETAINED here: this block
%  is initialization (ARRS Algorithm 4), not the Step-2b production repair.
% =====================================================================
if use_signrestrictions
    delta_input = [];
    if has_hetero, delta_input = delta0; end
    [lh0, likelihood0, A0_prior0, lambda_prior0, var0] = ...
        eval_posterior(x0, y, lags, lcA0, lcLmd, Tsigbrk, delta_input, tvA, options);
    [cache0, stable0] = btc(x0, var0);

    seed_tries = 0;
    while ~stable0 && seed_tries < 1000
        seed_tries = seed_tries + 1;
        [lh0, likelihood0, A0_prior0, lambda_prior0, var0] = ...
            eval_posterior(x0, y, lags, lcA0, lcLmd, Tsigbrk, delta_input, tvA, options);
        [cache0, stable0] = btc(x0, var0);
    end
    if ~stable0
        error(['Could not obtain a stable seed cache for the sampling loop ' ...
               'after %d redraws (companion explosive on every draw).'], seed_tries);
    end

    check_fn_seed = @(Q) check_sr_pass_cached(Q, cache0, var0, ...
        SignRestrictions, NarrativeRestrictions, sr_max_horizon, ...
        sign_regime_dependent, penalty_offdiagonal_on, ny);
    seed_rescue_rounds = get_opt(options, 'seed_rescue_rounds', 25);

    sr_ok_seed = check_fn_seed(Q_current);
    for rescue_round = 1:seed_rescue_rounds
        if sr_ok_seed, break; end
        b_tries = 0;
        while ~sr_ok_seed && b_tries < B_redraw_budget
            b_tries = b_tries + 1;
            [lh0, likelihood0, A0_prior0, lambda_prior0, var0] = ...
                eval_posterior(x0, y, lags, lcA0, lcLmd, Tsigbrk, delta_input, tvA, options);
            [cache0, stable0] = btc(x0, var0);
            if ~stable0, continue; end
            check_fn_seed = @(Q) check_sr_pass_cached(Q, cache0, var0, ...
                SignRestrictions, NarrativeRestrictions, sr_max_horizon, ...
                sign_regime_dependent, penalty_offdiagonal_on, ny);
            sr_ok_seed = check_fn_seed(Q_current);
        end
        total_B_redraws = total_B_redraws + b_tries;
        if sr_ok_seed, break; end
        if stable0
            [X_rep, Q_rep, rep_ok, ~] = find_admissible_Q_columnwise( ...
                ny, check_fn_seed, cache0.A0inv, cache0.lambda, impact_signs, ...
                sign_regime_dependent, 500, repair_k);
            if ~rep_ok
                [X_rep, Q_rep, rep_ok, ~] = ...
                    find_admissible_Q(check_fn_seed, ny, repair_k, repair_k);
            end
            if rep_ok
                X_current = X_rep; Q_current = Q_rep;
                sr_ok_seed = true; break
            end
        end
    end

    if ~sr_ok_seed && exist('use_local_pert_fallback', 'var') && use_local_pert_fallback
        fprintf('  [SEED RESCUE] %d rounds failed; trying local theta perturbations...\n', ...
            seed_rescue_rounds);
        for i_pert = 1:local_pert_tries
            x0_pert = x0 + mvnrnd(zeros(length(x0), 1), local_pert_scale * Sigma)';
            [lh_p, lik_p, A0p_p, lmdp_p, var_pert] = ...
                eval_posterior(x0_pert, y, lags, lcA0, lcLmd, Tsigbrk, delta_input, tvA, options);
            [cache_pert, stable_pert] = btc(x0_pert, var_pert);
            if ~stable_pert, continue; end
            check_fn_pert = @(Q) check_sr_pass_cached(Q, cache_pert, var_pert, ...
                SignRestrictions, NarrativeRestrictions, sr_max_horizon, ...
                sign_regime_dependent, penalty_offdiagonal_on, ny);
            if check_fn_pert(Q_current)
                rep_ok = true; X_rep = X_current; Q_rep = Q_current;
            else
                [X_rep, Q_rep, rep_ok, ~] = find_admissible_Q_columnwise( ...
                    ny, check_fn_pert, cache_pert.A0inv, cache_pert.lambda, impact_signs, ...
                    sign_regime_dependent, 500, repair_k);
            end
            if rep_ok
                x0 = x0_pert; var0 = var_pert; cache0 = cache_pert;
                lh0 = lh_p; likelihood0 = lik_p; A0_prior0 = A0p_p; lambda_prior0 = lmdp_p;
                X_current = X_rep; Q_current = Q_rep;
                check_fn_seed = check_fn_pert;
                sr_ok_seed = true;
                fprintf('  [SEED RESCUE] theta perturbation %d succeeded.\n', i_pert);
                break
            end
        end
    end

    if ~sr_ok_seed
        error(['Could not establish an admissible (Q, A+) seed after %d rescue ' ...
               'rounds + theta perturbations. Loosen restrictions or raise ' ...
               'seed_rescue_rounds.'], seed_rescue_rounds);
    end
end
% (x0, var0, cache0, Q_current) is now admissible.

% =====================================================================
%  MAIN LOOP (Sec 3)
% =====================================================================
while collected_draws < target_draws
    idraw = idraw + 1;

    if use_signrestrictions

        for isep = 1:nsep

            % ============================================================
            %  PRE-MIXING (Sec 3): premix_m local ESS on Q at the
            %  current admissible (x0, var0, Q).
            % ============================================================
            check_fn = @(Q) check_sr_pass_cached(Q, cache0, var0, ...
                SignRestrictions, NarrativeRestrictions, sr_max_horizon, ...
                sign_regime_dependent, penalty_offdiagonal_on, ny);
            for ipremix = 1:premix_m
                tBlk = tic;
                [X_current, Q_current, ess_ne] = ...
                    ess_draw_Q_columnwise(X_current, check_fn, ess_max_shrinks_per_col);
                prof.t_ess = prof.t_ess + toc(tBlk); prof.n_ess = prof.n_ess + 1;
                total_ess_evals = total_ess_evals + ess_ne;
            end

            % ============================================================
            %  STEP 1 — update (A0,Lambda,A+) and restore admissibility
            % ============================================================

            % (a) baseline collapsed marginal lh0 at the current delta + FRESH
            %     A+ at x0 (the reject-branch completion).
            tBlk = tic;
            if has_hetero
                [lh0, likelihood0, A0_prior0, lambda_prior0, var0_rej] = ...
                    eval_posterior(x0, y, lags, lcA0, lcLmd, Tsigbrk, delta0, tvA, options);
            else
                [lh0, likelihood0, A0_prior0, lambda_prior0, var0_rej] = ...
                    eval_posterior(x0, y, lags, lcA0, lcLmd, Tsigbrk, [], tvA, options);
            end
            prof.t_eval = prof.t_eval + toc(tBlk); prof.n_eval = prof.n_eval + 1;

            % (a') propose z0 = (A0,Lambda), complete with a fresh A+ at x1
            x1 = x0 + mvnrnd(zeros(length(x0), 1), Sigma)';
            tBlk = tic;
            if has_hetero
                [lh1, likelihood1, A0_prior1, lambda_prior1, var1] = ...
                    eval_posterior(x1, y, lags, lcA0, lcLmd, Tsigbrk, delta0, tvA, options);
            else
                [lh1, likelihood1, A0_prior1, lambda_prior1, var1] = ...
                    eval_posterior(x1, y, lags, lcA0, lcLmd, Tsigbrk, [], tvA, options);
            end
            prof.t_eval = prof.t_eval + toc(tBlk); prof.n_eval = prof.n_eval + 1;

            total_mh_proposals = total_mh_proposals + 1;
            mh_accept = ~isempty(var1) && isstruct(var1) && isfield(var1,'Bdraw') ...
                        && (log(rand()) < lh0 - lh1);

            % (b) select z0* and the FRESH A+ completion at z0*
            if mh_accept
                total_mh_accepted = total_mh_accepted + 1;
                x_star   = x1;   var_star = var1;   lh_star  = lh1;
                lik_star = likelihood1; A0p_star = A0_prior1; lmdp_star = lambda_prior1;
            else
                x_star   = x0;   var_star = var0_rej; lh_star = lh0;
                lik_star = likelihood0; A0p_star = A0_prior0; lmdp_star = lambda_prior0;
            end

            % (c) SR on the completed draw theta* = (z0*, A+*)
            tBlk = tic;
            [cache_star, stable_star] = btc(x_star, var_star);
            prof.t_cache = prof.t_cache + toc(tBlk); prof.n_cache = prof.n_cache + 1;

            if ~stable_star
                total_sampler_failures = total_sampler_failures + 1;
            else
                tBlk = tic;
                sr_pass = check_sr_pass_cached(Q_current, cache_star, var_star, ...
                    SignRestrictions, NarrativeRestrictions, sr_max_horizon, ...
                    sign_regime_dependent, penalty_offdiagonal_on, ny);
                prof.t_repair = prof.t_repair + toc(tBlk); prof.n_repair = prof.n_repair + 1;

                if sr_pass
                    % completed draw admissible under current Q -> keep it.
                    total_sr_direct = total_sr_direct + 1;
                    x0 = x_star; var0 = var_star; lh0 = lh_star;
                    likelihood0 = lik_star; A0_prior0 = A0p_star; lambda_prior0 = lmdp_star;
                    cache0 = cache_star;
                else
                    % ===== JOINT (A0,Lambda,A+) redraw (Prop 1) =====
                    % Stale completed draw (accept OR reject): repair by the
                    % JOINT redraw of the whole structural collection
                    % (Proposition 1 / literal Sec 3 Step 1). Propose z0_cand
                    % from x0 by the collapsed first-block MH, complete with a
                    % fresh A+_cand, keep the first sign-admissible candidate.
                    % Q is kept.
                    repaired = false;
                    total_repair_attempts = total_repair_attempts + 1;
                    if any(lcA0(:)) || any(lcLmd(:))
                        for t_try = 1:theta_redraw_budget
                            total_theta_redraws = total_theta_redraws + 1;
                            x_th = x0 + mvnrnd(zeros(length(x0), 1), Sigma)';
                            tBlk = tic;
                            if has_hetero
                                [lh_th, lik_th, A0p_th, lmdp_th, var_th] = ...
                                    eval_posterior(x_th, y, lags, lcA0, lcLmd, Tsigbrk, delta0, tvA, options);
                            else
                                [lh_th, lik_th, A0p_th, lmdp_th, var_th] = ...
                                    eval_posterior(x_th, y, lags, lcA0, lcLmd, Tsigbrk, [], tvA, options);
                            end
                            prof.t_eval = prof.t_eval + toc(tBlk); prof.n_eval = prof.n_eval + 1;

                            mh_accept_th = ~isempty(var_th) && isstruct(var_th) ...
                                && isfield(var_th,'Bdraw') && (log(rand()) < lh0 - lh_th);
                            if ~mh_accept_th, continue; end

                            tBlk = tic;
                            [cache_th, stable_th] = btc(x_th, var_th);
                            prof.t_cache = prof.t_cache + toc(tBlk); prof.n_cache = prof.n_cache + 1;
                            if ~stable_th, continue; end

                            tBlk = tic;
                            sr_pass_th = check_sr_pass_cached(Q_current, cache_th, var_th, ...
                                SignRestrictions, NarrativeRestrictions, sr_max_horizon, ...
                                sign_regime_dependent, penalty_offdiagonal_on, ny);
                            prof.t_repair = prof.t_repair + toc(tBlk); prof.n_repair = prof.n_repair + 1;
                            if ~sr_pass_th, continue; end

                            x0 = x_th; var0 = var_th; lh0 = lh_th;
                            likelihood0 = lik_th; A0_prior0 = A0p_th; lambda_prior0 = lmdp_th;
                            cache0 = cache_th;
                            total_theta_repair_success = total_theta_repair_success + 1;
                            repaired = true;
                            break
                        end
                    end

                    if ~repaired
                        % ===== GLOBAL Q-repair under the ORIGINAL theta* =====
                        % (Sec 3 Step 2). Shared fallback for both paths.
                        %  Paper Algorithm 1 Step 2b: independent
                        % Haar draws only (up to N_Q = repair_k), so the first
                        % admissible draw is restricted-Haar (Prop. diagnostics
                        % (i)). The columnwise fallback is REMOVED; on budget
                        % exhaustion the event is a discarded transition and
                        % the chain holds at its current admissible state.
                        check_fn_repair = @(Q) check_sr_pass_cached(Q, cache_star, var_star, ...
                            SignRestrictions, NarrativeRestrictions, sr_max_horizon, ...
                            sign_regime_dependent, penalty_offdiagonal_on, ny);
                        tBlk = tic;
                        [X_repair, Q_repair, repair_ok, ~] = ...
                            find_admissible_Q(check_fn_repair, ny, repair_k, repair_k);
                        prof.t_repair = prof.t_repair + toc(tBlk); prof.n_repair = prof.n_repair + 1;

                        if repair_ok
                            x0 = x_star; var0 = var_star; lh0 = lh_star;
                            likelihood0 = lik_star; A0_prior0 = A0p_star; lambda_prior0 = lmdp_star;
                            cache0 = cache_star;
                            X_current = X_repair; Q_current = Q_repair;
                            total_repair_success = total_repair_success + 1;
                        else
                            total_sampler_failures = total_sampler_failures + 1;
                        end
                    end
                end
            end

            % (f) hetero: update latent scales delta (third Gibbs block).
            if has_hetero
                tBlk = tic;
                eout = var0.udraw(1:size(delta0, 1), :) .* exp(delta0);
                if isempty(tparam)
                    delta0 = drawdelta(eout, alpha, K);
                else
                    delta0 = drawt(eout, tparam, (tscale^2) * tparam);
                end
                prof.t_delta = prof.t_delta + toc(tBlk); prof.n_delta = prof.n_delta + 1;
            end
        end  % isep

        % ============================================================
        %  STEP 3 (Sec 3): run ness_store (= N_ess) local ESS at
        %  the FINAL theta = x0 and RETAIN ALL. Carry Q <- Q^(N_ess).
        % ============================================================
        total_parameter_draws = total_parameter_draws + 1;

        check_fn_store = @(Q) check_sr_pass_cached(Q, cache0, var0, ...
            SignRestrictions, NarrativeRestrictions, sr_max_horizon, ...
            sign_regime_dependent, penalty_offdiagonal_on, ny);

        for iess = 1:ness_store
            if collected_draws >= target_draws, break; end

            tBlk = tic;
            [X_current, Q_current, ess_ne] = ...
                ess_draw_Q_columnwise(X_current, check_fn_store, ess_max_shrinks_per_col);
            prof.t_ess = prof.t_ess + toc(tBlk); prof.n_ess = prof.n_ess + 1;
            total_ess_evals = total_ess_evals + ess_ne;

            tBlk = tic;
            [sr_store_ok, ~, irf_stored, ~, Qdraw_stored, xrot_stored, st_store, ...
             lambda_stored, By_stored, fsign_stored, fsign_cell_stored] = ...
                SignRestrictionCheck( ...
                    x0, var0, ...
                    SignRestrictions, NarrativeRestrictions, ...
                    lags, lcA0, constant, ny, irf_horizon, ...
                    sign_regime_dependent, penalty_offdiagonal_on, ...
                    lcLmd, fix_first_regime, Q_current, SR_parsed, ...
                    NarrativeRegimeRestrictions, enforce_stability, lcA0_tv, SR_parsed_byregime);
            prof.t_src = prof.t_src + toc(tBlk); prof.n_src = prof.n_src + 1;

            if sr_store_ok ~= 1 || st_store ~= 1 || isempty(xrot_stored)
                total_parameter_draws_discarded = total_parameter_draws_discarded + 1;
                continue
            end

            collected_draws = collected_draws + 1;

            tBlk = tic;
            [lh_wQ, likelihood_wQ, A0_prior_wQ, lambda_prior_wQ, ~] = ...
                eval_posterior(xrot_stored, y, lags, lcA0, lcLmd, Tsigbrk, delta0, tvA, options);
            prof.t_eval = prof.t_eval + toc(tBlk); prof.n_eval = prof.n_eval + 1;

            if ~isempty(fsign_cell_stored)
                fs_mat = nan(ny, nRegimes);
                for rr = 1:nRegimes
                    fs_mat(:, rr) = fsign_cell_stored{rr}(:);
                end
            else
                fs_mat = nan(ny, nRegimes);
            end

            Phi = var0.Bdraw;

            sign_irfs_all{collected_draws}           = irf_stored;
            draw_lh_wQ_all(collected_draws)          = lh_wQ;
            draw_likelihood_wQ_all(collected_draws)  = likelihood_wQ;
            draw_A0prior_wQ_all(collected_draws)     = A0_prior_wQ;
            draw_lambdaprior_wQ_all(collected_draws) = lambda_prior_wQ;
            draw_Qdraw_all{collected_draws}          = Qdraw_stored;
            draw_Lambda_all{collected_draws}         = lambda_stored;
            draw_By_all{collected_draws}             = By_stored;
            draw_fsign_all{collected_draws}          = fsign_stored;
            draw_fsign_cell_all{collected_draws}     = fs_mat;


            %  combined weight = v(Q; theta) / omega_hat(theta, Q)
            %  gated on use_nrr_omega: with the toggle OFF,
            % omega_hat is never simulated, draw_omega stays 1, and the
            % weight is NOT divided -- the retained draws target the
            % UNWEIGHTED restricted-Haar posterior of Definition 1.
            if use_nrr_omega && (has_narr || has_narr_regime)   % point and/or window anchors
                nrr_omega_opts.moment_check = (collected_draws == 1);   % once per run
                nrr_omega_opts.validate     = (collected_draws == 1);   % [PERF] fast-vs-exact gate
                [~, omega_hat, om_diag] = check_sr_pass_cached(Q_current, cache0, var0, ...
                    SignRestrictions, NarrativeRestrictions, sr_max_horizon, ...
                    sign_regime_dependent, penalty_offdiagonal_on, ny, nrr_omega_opts);
                draw_omega_all(collected_draws)  = omega_hat;
                om_diag_all{collected_draws}     = om_diag;   %  group omegas / floor flags
                draw_weight_all(collected_draws) = draw_weight_all(collected_draws) / omega_hat;
                if nrr_omega_opts.moment_check && isfield(om_diag, 'xi_kurtosis')
                    fprintf('  [NRR-W] shape check (n=%d rows) | xi std: [%s]\n', ...
                        om_diag.n_rows, num2str(om_diag.xi_std, '%.3f '));
                    fprintf('  [NRR-W]  xi kurtosis: [%s]  (~3 -> gaussian; >>3 -> student; scale irrelevant)\n', ...
                        num2str(om_diag.xi_kurtosis, '%.2f '));
                end
            end

            draw_lh_all(collected_draws)          = lh0;
            draw_likelihood_all(collected_draws)  = likelihood0;
            draw_A0prior_all(collected_draws)     = A0_prior0;
            draw_lambdaprior_all(collected_draws) = lambda_prior0;
            draw_x0_all{collected_draws}          = x0;
            draw_Phi_all{collected_draws}         = Phi;
            draw_udraw_all{collected_draws}       = var0.udraw(1:T-lags, :);
            draw_udraw_true_all{collected_draws}  = var0.udraw_true(1:T-lags, :);
            if has_hetero
                draw_dout_all{collected_draws} = delta0;
            end
        end

    else
        % ---- non-SR branch (unchanged) ----
        accepted_this_draw = false;

        for isep = 1:nsep
            delta_input = [];
            if has_hetero, delta_input = delta0; end
            [lh0, likelihood0, A0_prior0, lambda_prior0, var0] = ...
                eval_posterior(x0, y, lags, lcA0, lcLmd, Tsigbrk, delta_input, tvA, options);

            x1 = x0 + mvnrnd(zeros(length(x0), 1), Sigma)';
            if has_hetero
                [lh1, likelihood1, A0_prior1, lambda_prior1, var1] = ...
                    eval_posterior(x1, y, lags, lcA0, lcLmd, Tsigbrk, delta0, tvA, options);
            else
                [lh1, likelihood1, A0_prior1, lambda_prior1, var1] = ...
                    eval_posterior(x1, y, lags, lcA0, lcLmd, Tsigbrk, [], tvA, options);
            end

            if log(rand()) < lh0 - lh1
                x0 = x1;  var0 = var1;  lh0 = lh1;
                likelihood0  = likelihood1;
                A0_prior0    = A0_prior1;
                lambda_prior0 = lambda_prior1;
                accepted_this_draw = true;
            end

            if has_hetero
                eout = var0.udraw(1:size(delta0, 1), :) .* exp(delta0);
                if isempty(tparam), delta0 = drawdelta(eout, alpha, K);
                else,               delta0 = drawt(eout, tparam, (tscale^2)*tparam); end
            end
        end

        total_parameter_draws = total_parameter_draws + 1;
        if accepted_this_draw
            total_parameter_draws_accepted = total_parameter_draws_accepted + 1;
        end

        collected_draws = collected_draws + 1;
        Phi = var0.Bdraw;

        draw_lh_all(collected_draws)          = lh0;
        draw_likelihood_all(collected_draws)  = likelihood0;
        draw_A0prior_all(collected_draws)     = A0_prior0;
        draw_lambdaprior_all(collected_draws) = lambda_prior0;
        draw_x0_all{collected_draws}          = x0;

        var_trimmed = struct();
        var_trimmed.Bdraw      = var0.Bdraw;
        var_trimmed.udraw_true = var0.udraw_true;
        draw_var_all{collected_draws} = var_trimmed;

        draw_Phi_all{collected_draws}         = Phi;
        draw_udraw_all{collected_draws}       = var0.udraw(1:T-lags, :);
        draw_udraw_true_all{collected_draws}  = var0.udraw_true(1:T-lags, :);
        if has_hetero
            draw_dout_all{collected_draws} = delta0;
        end
    end

    if mod(idraw, max(1, floor(target_draws/50))) == 0 || collected_draws >= target_draws
        time_elapsed = toc(samplein_start_time);
        if use_signrestrictions
            fprintf('[SAMPLING] %4d/%4d | MH:%d | SR:%d | Joint:%d/%d | Qrep:%d/%d | Fail:%d | Skip:%d | %.1fs\n', ...
                collected_draws, target_draws, total_mh_accepted, total_sr_direct, ...
                total_theta_repair_success, total_theta_redraws, ...
                total_repair_success, total_repair_attempts, ...
                total_sampler_failures, total_parameter_draws_discarded, time_elapsed);
        else
            fprintf('[SAMPLING] %4d/%4d | Time: %.1f sec\n', ...
                collected_draws, target_draws, time_elapsed);
        end
    end
end

% =====================================================================
%  BLOCK-LEVEL PROFILE PRINT (sampling loop only)
% =====================================================================
prof.t_wall = toc(samplein_start_time);
prof.t_sum  = prof.t_eval + prof.t_cache + prof.t_ess + prof.t_repair + prof.t_src + prof.t_delta;
fprintf('\n========== BLOCK-LEVEL PROFILE (ESS, sampling loop) ==========\n');
fprintf('%-26s %12s %12s %9s %12s\n', 'block', 'calls', 'total (s)', '%wall', 'ms/call');
fprintf('%s\n', repmat('-', 1, 76));
fprintf('%-26s %12d %12.2f %8.1f%% %12.3f\n', 'eval_posterior',        prof.n_eval,   prof.t_eval,   100*prof.t_eval/max(prof.t_wall,eps),   1000*prof.t_eval/max(prof.n_eval,1));
fprintf('%-26s %12d %12.2f %8.1f%% %12.3f\n', 'build_theta_cache',     prof.n_cache,  prof.t_cache,  100*prof.t_cache/max(prof.t_wall,eps),  1000*prof.t_cache/max(prof.n_cache,1));
fprintf('%-26s %12d %12.2f %8.1f%% %12.3f\n', 'ess_draw_Q_columnwise', prof.n_ess,    prof.t_ess,    100*prof.t_ess/max(prof.t_wall,eps),    1000*prof.t_ess/max(prof.n_ess,1));
fprintf('%-26s %12d %12.2f %8.1f%% %12.3f\n', 'sr-check + repair',     prof.n_repair, prof.t_repair, 100*prof.t_repair/max(prof.t_wall,eps), 1000*prof.t_repair/max(prof.n_repair,1));
fprintf('%-26s %12d %12.2f %8.1f%% %12.3f\n', 'SignRestrictionCheck',  prof.n_src,    prof.t_src,    100*prof.t_src/max(prof.t_wall,eps),    1000*prof.t_src/max(prof.n_src,1));
fprintf('%-26s %12d %12.2f %8.1f%% %12.3f\n', 'drawt/drawdelta',       prof.n_delta,  prof.t_delta,  100*prof.t_delta/max(prof.t_wall,eps),  1000*prof.t_delta/max(prof.n_delta,1));
fprintf('%s\n', repmat('-', 1, 76));
fprintf('%-26s %12s %12.2f %8.1f%%\n', 'timed subtotal',          '', prof.t_sum,              100*prof.t_sum/max(prof.t_wall,eps));
fprintf('%-26s %12s %12.2f %8.1f%%\n', 'unaccounted (overhead)',  '', prof.t_wall - prof.t_sum, 100*(prof.t_wall-prof.t_sum)/max(prof.t_wall,eps));
fprintf('%-26s %12s %12.2f\n',          'total sampling wall',     '', prof.t_wall);
fprintf('==============================================================\n');
output.prof = prof;   % store so it survives the save()

%% ================================================================
%%  SECTION 11 — SUMMARY
%% ================================================================
fprintf('\n========== SAMPLING SUMMARY (Cached + Fast SR + NRR) ===========\n');
fprintf('Total outer iterations: %d\n', idraw);
fprintf('Collected draws:        %d / %d\n', collected_draws, target_draws);
fprintf('Total time:             %.2f seconds\n', toc(samplein_start_time));
 
if use_signrestrictions
    eff_theta = total_sr_direct + total_repair_success + total_theta_repair_success;
    fprintf('\n--- Metropolis-Hastings ---\n');
    fprintf('  Proposals:              %d\n', total_mh_proposals);
    fprintf('  MH accepted:            %d (%.2f%%)\n', ...
        total_mh_accepted, 100 * total_mh_accepted / max(total_mh_proposals, 1));
    fprintf('  SR direct:              %d\n', total_sr_direct);
    fprintf('  Joint structural-repairs: %d / %d (N_theta %d)\n', ...
        total_theta_repair_success, total_theta_redraws, theta_redraw_budget);
    fprintf('  Global Q-repair successes: %d / %d (Haar-only)\n', total_repair_success, total_repair_attempts);
    fprintf('  Seed-rescue A+ redraws: %d (budget %d per round; initialization only)\n', total_B_redraws, B_redraw_budget);
    fprintf('  Sampler failures:       %d\n', total_sampler_failures);
    fprintf('  Effective theta accept: %d (%.2f%%)\n', ...
        eff_theta, 100 * eff_theta / max(total_mh_proposals, 1));
    fprintf('  ESS evaluations:        %d\n', total_ess_evals);
    fprintf('  Draws skipped:          %d\n', total_parameter_draws_discarded);
 
    % --- Stale-Q repair breakdown (how the fallback is being used) ---------
    % A "stale-Q" event = a completed (A0,Lambda,A+) draw -- MH-accept or
    % MH-reject -- that left the current Q inadmissible (total_repair_attempts).
    % Each is
    % resolved by (i) a JOINT (A0,Lambda,A+) structural redraw [keeps the
    % ESS-local Q -- Sec 3 step 1 / Prop 1], (ii) the global Q-repair
    % FALLBACK [replaces Q with an independent admissible Haar draw under the
    % original accepted theta' -- Haar-only, so the repaired Q is a
    % valid restricted-Haar draw by construction], or (iii) neither (discarded
    % transition: the chain holds at its current admissible state). A high
    % fallback share is now purely a computational-efficiency signal
    % (Prop. diagnostics (iii)); it no longer risks mixing in a
    % non-restricted-Haar Q. Worth watching across datasets.
    stale_events   = total_repair_attempts;
    fallback_fired = total_repair_success;          % global Q-repair successes
    struct_fixed   = total_theta_repair_success;    % JOINT structural-redraw successes
    stale_failed   = max(stale_events - fallback_fired - struct_fixed, 0);
    fallback_share_theta = 100 * fallback_fired / max(eff_theta, 1);
    fprintf('\n--- Stale-Q repair breakdown ---\n');
    fprintf('  Stale-Q events:                           %d\n', stale_events);
    fprintf('    resolved by JOINT structural redraw:    %d (%.1f%%)\n', ...
        struct_fixed,   100 * struct_fixed   / max(stale_events, 1));
    fprintf('    resolved by global Q-repair (fallback): %d (%.1f%%)\n', ...
        fallback_fired, 100 * fallback_fired / max(stale_events, 1));
    fprintf('    unresolved (discarded transition):      %d (%.1f%%)\n', ...
        stale_failed,   100 * stale_failed   / max(stale_events, 1));
    fprintf('  Fallback share of effective theta moves:  %.2f%% (%d / %d)\n', ...
        fallback_share_theta, fallback_fired, eff_theta);
    if fallback_share_theta >= 10
        fprintf(['  [NOTE] Fallback supplies >=10%% of theta moves. Under [HAAR-ONLY]\n' ...
                 '         every repaired Q is restricted-Haar, so this is a pure\n' ...
                 '         efficiency signal: consider raising theta_redraw_budget\n' ...
                 '         (N_theta) so more stale events resolve by joint repair.\n']);
    end
    if use_draw_weights   % [W] combined v(Q;theta) / omega_hat weights
        w_all = draw_weight_all(1:collected_draws);
        ess_share = (sum(w_all)^2) / (sum(w_all.^2) * max(numel(w_all),1));
        fprintf('  [W] combined weights:   min %.3e | median %.3e | max %.3e\n', ...
            min(w_all), median(w_all), max(w_all));
        fprintf('  [W] combined weight ESS share: %.3f (1.0 = uniform; low => few draws dominate)\n', ess_share);
    end
    if use_nrr_omega && (has_narr || has_narr_regime)   %  omega diagnostics
        om_all = draw_omega_all(1:collected_draws);
        fprintf('  [NRR-W] omega_hat: min %.3e | median %.3e | max %.3e | M = %d\n', ...
            min(om_all), median(om_all), max(om_all), nrr_omega_M);
        %  group-level floor check (product form): the JOINT
        % omega_hat can sit below 1/M legitimately, so floor status must be
        % judged PER GROUP, not on the joint value.
        od_ok = cellfun(@(d) isstruct(d) && isfield(d, 'group_omegas') && ...
            ~isempty(d.group_omegas), om_diag_all(1:collected_draws));
        if any(od_ok)
            G  = cell2mat(cellfun(@(d) d.group_omegas(:).', ...
                om_diag_all(od_ok), 'UniformOutput', false));
            fv = 1 / nrr_omega_M;
            fprintf('  [NRR-W] n_groups: %d (date-disjoint restriction groups)\n', size(G, 2));
            fprintf('  [NRR-W] group omega medians:'); fprintf(' %.3e', median(G, 1)); fprintf('\n');
            fprintf('  [NRR-W] group omega mins:   '); fprintf(' %.3e', min(G, [], 1)); fprintf('\n');
            fprintf('  [NRR-W] floored per group:  '); fprintf(' %d', sum(G <= fv + eps, 1)); fprintf('\n');
            fprintf('  [NRR-W] floored draws (any group at 1/M): %d / %d\n', ...
                sum(any(G <= fv + eps, 2)), size(G, 1));
        else
            fprintf('  [NRR-W] floored (joint <= 1/M): %d / %d\n', ...
                sum(om_all <= 1/nrr_omega_M + eps), collected_draws);
        end
    elseif (has_narr || has_narr_regime) && ~use_nrr_omega
        fprintf('  [NRR-W] omega weighting OFF (options.use_nrr_omega = false):\n');
        fprintf('  [NRR-W] narrative enters ONLY through the SR indicator; the\n');
        fprintf('  [NRR-W] retained draws target the UNWEIGHTED restricted-Haar\n');
        fprintf('  [NRR-W] posterior of Definition 1 (no 1/omega reweighting).\n');
    end
else
    fprintf('  Draws with >= 1 acc:    %d\n', total_parameter_draws_accepted);
    fprintf('  Acceptance rate:        %.2f%%\n', ...
        100 * total_parameter_draws_accepted / max(total_parameter_draws, 1));
end
fprintf('===========================================================\n');

%% ================================================================
%%  SECTION 12 — POST-PROCESSING
%% ================================================================

draw_lh_all          = draw_lh_all(1:collected_draws);
draw_likelihood_all  = draw_likelihood_all(1:collected_draws);
draw_A0prior_all     = draw_A0prior_all(1:collected_draws);
draw_lambdaprior_all = draw_lambdaprior_all(1:collected_draws);
draw_x0_all          = draw_x0_all(1:collected_draws);
draw_var_all         = draw_var_all(1:collected_draws);
draw_Phi_all         = draw_Phi_all(1:collected_draws);
draw_udraw_all       = draw_udraw_all(1:collected_draws);
draw_udraw_true_all  = draw_udraw_true_all(1:collected_draws);

if has_hetero
    draw_dout_all = draw_dout_all(1:collected_draws);
end

if use_signrestrictions
    sign_irfs_all           = sign_irfs_all(1:collected_draws);
    draw_lh_wQ_all          = draw_lh_wQ_all(1:collected_draws);
    draw_likelihood_wQ_all  = draw_likelihood_wQ_all(1:collected_draws);
    draw_A0prior_wQ_all     = draw_A0prior_wQ_all(1:collected_draws);
    draw_lambdaprior_wQ_all = draw_lambdaprior_wQ_all(1:collected_draws);
    draw_Qdraw_all          = draw_Qdraw_all(1:collected_draws);
    draw_Lambda_all         = draw_Lambda_all(1:collected_draws);
    draw_By_all             = draw_By_all(1:collected_draws);
    draw_fsign_all          = draw_fsign_all(1:collected_draws);
    draw_fsign_cell_all     = draw_fsign_cell_all(1:collected_draws);
    draw_weight_all         = draw_weight_all(1:collected_draws);
    draw_omega_all          = draw_omega_all(1:collected_draws);
    om_diag_all             = om_diag_all(1:collected_draws);
end

if nRegimes > 1 && use_signrestrictions
    % --- per-regime IRF reconstruction (A0 common across regimes) ---
    % The contemporaneous block A0 is common; regimes differ only through
    % lambda_r inside the impact matrix Smat. By = A0 \ Aplus matches the
    % stored By.
    K_shocks = numel(accepted_shock);

    xA_all = zeros(collected_draws, length(draw_x0_all{1}));
    for j = 1:collected_draws
        xA_all(j, :) = draw_x0_all{j}';
    end

    draw_Qdraw_mat      = cat(3, draw_Qdraw_all{:});
    draw_Lambda_mat     = cat(3, draw_Lambda_all{:});
    draw_fsign_cell_mat = cat(3, draw_fsign_cell_all{:});
    draw_Phi_mat        = cat(3, draw_Phi_all{:});   % common lag block (Aplus)

    for iReg = 1:nRegimes
        irf_save = zeros(ny, irf_horizon, K_shocks, collected_draws);
        parfor j = 1:collected_draws
            % contemporaneous matrix (common across regimes)
            A0 = zeros(ny);
            A0(lcA0) = xA_all(j, 1:sum(lcA0(:)));

            % reduced form: By = A0 \ Aplus
            Phi   = draw_Phi_mat(:, :, j);
            Btrim = Phi(1:lags*ny, :);
            By    = zeros(ny, ny, lags);
            for iLag = 1:lags
                rows = (iLag-1)*ny + (1:ny);
                By(:, :, iLag) = A0 \ Btrim(rows, :)';
            end

            Q   = draw_Qdraw_mat(:, :, j);
            lam = draw_Lambda_mat(:, iReg, j);
            Dhalf = diag(sqrt(max(lam, 0)));
            if sign_regime_dependent == 0
                Smat = (A0 \ Q') * Dhalf;
            else
                Smat = (A0 \ Dhalf) * Q';
            end
            irTemp = impulsdtrf(By, Smat, irf_horizon);

            fsign_vec = draw_fsign_cell_mat(:, iReg, j);
            tmp = zeros(ny, irf_horizon, K_shocks);
            for k = 1:K_shocks
                s = accepted_shock(k);
                tmp(:, :, k) = fsign_vec(s) * irTemp(:, :, s);
            end
            irf_save(:, :, :, j) = tmp;
        end
        output.sign_irf_regime(:, :, :, :, iReg) = irf_save;
    end

elseif ~use_signrestrictions && nRegimes > 1
    xA_all = zeros(collected_draws, length(draw_x0_all{1}));
    for j = 1:collected_draws
        xA_all(j, :) = draw_x0_all{j}';
    end
    xA_draws = xA_all(:, 1:sum(lcA0(:)));
    draw_Phi_mat = cat(3, draw_Phi_all{:});

    irf_save = zeros(ny, irf_horizon, ny, collected_draws);
    parfor j = 1:collected_draws
        A0 = zeros(ny);
        A0(lcA0) = xA_draws(j, :);
        Bdraw = draw_Phi_mat(:,:,j);
        if constant, Btrim = Bdraw(1:end-1, :); else, Btrim = Bdraw; end
        Aplus = NaN(ny, ny, lags);
        for iLag = 1:lags
            rows = (iLag - 1) * ny + (1:ny);
            Aplus(:,:,iLag) = Btrim(rows, :)';
        end
        By = zeros(ny, ny, lags);
        for iLag = 1:lags
            By(:,:,iLag) = A0 \ Aplus(:,:,iLag);
        end
        Smat = inv(A0);
        irfTemp = impulsdtrf(By, Smat, irf_horizon);
        irf_save(:,:,:,j) = irfTemp;
    end
    output.hetero_irf_temp = irf_save;
end

draw_lh          = draw_lh_all;
draw_likelihood  = draw_likelihood_all;
draw_A0prior     = draw_A0prior_all;
draw_lambdaprior = draw_lambdaprior_all;
draw_x0          = cell2mat(draw_x0_all);
draw_Phi         = cat(3, draw_Phi_all{:});
draw_udraw       = cat(3, draw_udraw_all{:});
draw_udraw_true  = cat(3, draw_udraw_true_all{:});

if has_hetero
    draw_dout = cat(3, draw_dout_all{:});
end

if use_signrestrictions
    draw_lh_wQ          = draw_lh_wQ_all;
    draw_likelihood_wQ  = draw_likelihood_wQ_all;
    draw_A0prior_wQ     = draw_A0prior_wQ_all;
    draw_lambdaprior_wQ = draw_lambdaprior_wQ_all;
    draw_Qdraw          = cat(3, draw_Qdraw_all{:});
    draw_Lambda         = cat(3, draw_Lambda_all{:});
    draw_By             = cat(4, draw_By_all{:});   % ny x ny x lags x nDraws
    draw_fsign          = cat(3, draw_fsign_all{:});
    draw_fsign_cell_mat = cat(3, draw_fsign_cell_all{:});
    sign_irfs_temp      = cat(4, sign_irfs_all{:});
end

output.collected_draws  = collected_draws;
output.total_iterations = idraw;

output.draw_lh          = draw_lh;
output.draw_likelihood  = draw_likelihood;
output.draw_A0prior     = draw_A0prior;
output.draw_lambdaprior = draw_lambdaprior;
output.draw_x0          = draw_x0;
output.draw_var         = draw_var_all;
output.draw_Phi         = draw_Phi;
output.draw_udraw       = draw_udraw;
output.draw_udraw_true  = draw_udraw_true;
output.lcA0             = lcA0;
output.lcLmd            = lcLmd;
output.lcA0_tv          = lcA0_tv;        % all-false mask (non-tvA build)
output.lrange           = 1:nRegimes;
output.varnames         = varnames;
output.tvA              = tvA;

if has_hetero
    output.draw_dout = draw_dout;
end

if use_signrestrictions
    output.sign_irfs_temp       = sign_irfs_temp;
    output.sign_accepted_shock  = accepted_shock;
    output.draw_lh_wQ           = draw_lh_wQ;
    output.draw_likelihood_wQ   = draw_likelihood_wQ;
    output.draw_A0prior_wQ      = draw_A0prior_wQ;
    output.draw_lambdaprior_wQ  = draw_lambdaprior_wQ;
    output.draw_Qdraw           = draw_Qdraw;
    output.draw_Lambda          = draw_Lambda;
    output.draw_By              = draw_By;
    output.draw_fsign           = draw_fsign;
    output.draw_fsign_cell      = draw_fsign_cell_mat;
    output.total_parameter_draws_discarded = total_parameter_draws_discarded;
    output.total_redraws = total_repair_attempts;
    output.avg_redraws_per_iteration = total_repair_attempts / max(idraw, 1);

    %  per-draw importance weight 1/omega_hat (== 1 vector when
    %        use_nrr_omega is off or there are no narrative restrictions).
    %        When active, use these as importance weights in ALL downstream
    %        estimates (weighted quantiles of sign_irfs_temp / sign_irf_regime
    %        along the draw dimension, weights = output.draw_weight).
    output.draw_weight        = draw_weight_all;
    output.draw_omega_nrr     = draw_omega_all;
    output.om_diag            = om_diag_all;          %  cell, one struct per draw
    output.use_draw_weights   = use_draw_weights;
    output.use_nrr_omega      = use_nrr_omega;
    output.nrr_omega_opts     = nrr_omega_opts;
    if use_draw_weights   % [W] combined-weight handles
        w_all = draw_weight_all;
        output.weight_ess_share = (sum(w_all)^2) / (sum(w_all.^2) * max(numel(w_all),1));

        %  sign_irf_regime / sign_irfs_temp are stored per-draw (draw dim
        %        is the LAST axis: 4 for sign_irfs_temp, 4 for sign_irf_regime
        %        before the trailing regime axis). SECTION 13/14 summaries above
        %        are ALREADY weighted. For any custom plotting of the raw IRF
        %        arrays, use these handles so the bands match:
        %          med = output.wmedian(output.sign_irfs_temp, 4);
        %          lo  = output.wquantile(output.sign_irfs_temp, 0.16, 4);
        %        (draw dim for sign_irf_regime is 4 as well, i.e. ny x H x Ksh x
        %         nDraws x nRegimes -> weight along dim 4.)
        output.wquantile = @(X, p, dim) local_wquantile_dim(X, draw_weight_all, p, dim);
        output.wmedian   = @(X, dim)    local_wmedian_dim(X, draw_weight_all, dim);
    end
end

output.ess_diagnostics.total_mh_proposals    = total_mh_proposals;
output.ess_diagnostics.total_mh_accepted     = total_mh_accepted;
output.ess_diagnostics.total_sr_direct       = total_sr_direct;
output.ess_diagnostics.total_repair_attempts = total_repair_attempts;
output.ess_diagnostics.total_repair_success  = total_repair_success;
output.ess_diagnostics.total_theta_redraws        = total_theta_redraws;        % (A0,Lambda) structural repair
output.ess_diagnostics.total_theta_repair_success = total_theta_repair_success; % (A0,Lambda) structural repair
output.ess_diagnostics.theta_redraw_budget        = theta_redraw_budget;        % N_theta
% --- Stale-Q repair breakdown (fallback monitoring; see SECTION 11) ---
output.ess_diagnostics.stale_q_events             = total_repair_attempts;       % completed draws that stranded Q (accept or reject)
output.ess_diagnostics.stale_resolved_structural  = total_theta_repair_success;  % fixed by (A0,Lambda) redraw
output.ess_diagnostics.stale_resolved_fallback    = total_repair_success;        % fixed by global Q-repair
output.ess_diagnostics.stale_unresolved           = max(total_repair_attempts - total_theta_repair_success - total_repair_success, 0);
output.ess_diagnostics.fallback_share_of_theta = total_repair_success / ...
    max(total_sr_direct + total_repair_success + total_theta_repair_success, 1);
output.ess_diagnostics.total_sampler_failures = total_sampler_failures;
output.ess_diagnostics.total_ess_evals       = total_ess_evals;
output.ess_diagnostics.total_B_redraws       = total_B_redraws;     % seed-rescue A+ redraws (init only)
output.ess_diagnostics.B_redraw_budget       = B_redraw_budget;     % seed-rescue budget
output.ess_diagnostics.premix_rounds         = premix_m;
output.ess_diagnostics.ness_store            = ness_store;        % N_ess retained Q per structural draw (Sec 2 step 3)
output.ess_diagnostics.repair_budget         = repair_k;
output.ess_diagnostics.enforce_stability     = enforce_stability;
output.ess_diagnostics.repair_order          = 'Haar_only';          %  columnwise fallback removed
output.ess_diagnostics.use_nrr_omega         = use_nrr_omega;

%% ================================================================
%%  SECTION 12.5 — STABILITY DIAGNOSTIC (report-only, no sampler effect)
%% ================================================================
% With enforce_stability=false the collected draws can include companions
% with a root on/above the unit circle. Compute max|root| per draw ONCE
% here and store it (+ an optional trim mask) so downstream plotting can
% drop extreme-explosive draws at REPORT time without recomputing eig and
% without ever touching the stored IRFs (sign_irfs_temp / sign_irf_regime
% stay complete). Controlled by options.report_stability_diag (default on
% when stability is off) and options.trim_threshold (default 1.02).
report_stab = get_opt(options, 'report_stability_diag', ~enforce_stability);

if use_signrestrictions && report_stab
    nD = collected_draws;
    draw_max_root = zeros(nD, 1);
    for j = 1:nD
        Phi = reshape(draw_By(:,:,:,j), ny, ny*lags);   % [lag1 | lag2 | ...]
        if lags > 1
            Comp = [Phi; eye(ny*(lags-1)), zeros(ny*(lags-1), ny)];
        else
            Comp = Phi;
        end
        draw_max_root(j) = max(abs(eig(Comp)));
    end

    trim_thresh = get_opt(options, 'trim_threshold', 1.02);
    output.draw_max_root  = draw_max_root;
    output.plot_keep_mask = draw_max_root <= trim_thresh;   % logical (nD x 1)
    output.trim_threshold = trim_thresh;

    sh_expl   = 100 * mean(draw_max_root >= 1.00);
    sh_strong = 100 * mean(draw_max_root >  1.05);
    fprintf('\n--- Stability diagnostic (report-only) ---\n');
    fprintf('  max|root|: median %.4f | 95pct %.4f | max %.4f\n', ...
        median(draw_max_root), quantile(draw_max_root, 0.95), max(draw_max_root));
    fprintf('  share |root| >= 1.00 : %5.1f%%\n', sh_expl);
    fprintf('  share |root| >  1.05 : %5.1f%%  (visibly grow by h~60)\n', sh_strong);
    fprintf('  plot_keep_mask (<= %.2f) keeps %d / %d (%.1f%%)\n', ...
        trim_thresh, sum(output.plot_keep_mask), nD, 100*mean(output.plot_keep_mask));
    if sh_strong < 1
        fprintf('  -> negligible: report as-is, no trim needed.\n');
    elseif sh_strong < 5
        fprintf('  -> mild: median/68%% fine; trim only if 95%% band / h>48 looks ragged.\n');
    else
        fprintf('  -> non-trivial: apply plot_keep_mask for long horizons.\n');
    end
    fprintf('------------------------------------------\n');
end

%% ================================================================
%%  SECTION 13 — VARIANCE RATIOS
%% ================================================================
if use_signrestrictions && nRegimes > 1
    Lambda_draws = cat(3, draw_Lambda_all{1:collected_draws});
    omega_draws  = Lambda_draws ./ Lambda_draws(:, 1, :);

    %  weight the draw dimension (3) when the omega weighting is active;
    %        otherwise keep median/quantile EXACTLY (byte-identical).
    if use_draw_weights   % [W]
        w13 = draw_weight_all(1:collected_draws);
        Lambda_median = local_wmedian_dim(Lambda_draws, w13, 3);
        Lambda_16     = local_wquantile_dim(Lambda_draws, w13, 0.16, 3);
        Lambda_84     = local_wquantile_dim(Lambda_draws, w13, 0.84, 3);
        omega_median  = local_wmedian_dim(omega_draws, w13, 3);
        omega_16      = local_wquantile_dim(omega_draws, w13, 0.16, 3);
        omega_84      = local_wquantile_dim(omega_draws, w13, 0.84, 3);
    else
        Lambda_median = median(Lambda_draws, 3);
        Lambda_16     = quantile(Lambda_draws, 0.16, 3);
        Lambda_84     = quantile(Lambda_draws, 0.84, 3);
        omega_median  = median(omega_draws, 3);
        omega_16      = quantile(omega_draws, 0.16, 3);
        omega_84      = quantile(omega_draws, 0.84, 3);
    end

    fprintf('\n============================================================\n');
    fprintf('  Estimated Variance Ratios: omega_{i,s} = lambda_{i,s} / lambda_{i,1}\n');
    fprintf('============================================================\n');
    fprintf('%-20s', 'Variable');
    for s = 1:nRegimes
        fprintf('  Regime %d (median [68%%])     ', s);
    end
    fprintf('\n');
    fprintf('%s\n', repmat('-', 1, 20 + nRegimes * 30));
    for i = 1:ny
        fprintf('%-20s', varnames{i});
        for s = 1:nRegimes
            fprintf('  %6.3f [%5.3f, %5.3f]     ', ...
                omega_median(i,s), omega_16(i,s), omega_84(i,s));
        end
        fprintf('\n');
    end
    fprintf('============================================================\n');

    fprintf('\n--- Proportionality Check ---\n');
    for s = 2:nRegimes
        omega_s_median = omega_median(:, s);
        fprintf('  Regime %d vs 1: std(omega) = %.4f  |  range = [%.3f, %.3f]\n', ...
            s, std(omega_s_median), min(omega_s_median), max(omega_s_median));
    end
    fprintf('============================================================\n');

    output.Lambda_median = Lambda_median;
    output.omega_median  = omega_median;
    output.omega_16      = omega_16;
    output.omega_84      = omega_84;
end

%% ================================================================
%%  SECTION 13.5 — B&F RANK ROBUSTNESS ACROSS DRAWS (report/appendix)
%% ================================================================
% Mode-point B&F check (SECTION 4) establishes local point-ID at a regular
% point. This loops the SAME rank condition over the kept posterior draws so
% the appendix can report that full rank holds across the posterior, i.e. the
% mode is not a non-regular point and the posterior sits inside the identified
% region. Report-only; does not touch sampler state.
if use_signrestrictions && nRegimes > 1 && any(lcLmd(:)) && isfield(output,'bf_check')

    % evaluate on the kept draws (post-trim if a stability mask exists)
    if isfield(output,'plot_keep_mask') && ~isempty(output.plot_keep_mask)
        bf_idx = find(output.plot_keep_mask);
    else
        bf_idx = 1:collected_draws;
    end

    % optional thinning so this stays cheap on big runs (every k-th draw)
    bf_stride = get_opt(options, 'bf_robust_stride', 1);
    bf_idx    = bf_idx(1:bf_stride:end);

    nbf       = numel(bf_idx);
    bf_ranks  = zeros(nbf, 1);
    bf_svgap  = zeros(nbf, 1);
    a_dim     = output.bf_check.dim_theta;

    for jj = 1:nbf
        bj = check_bf_identification(draw_x0_all{bf_idx(jj)}, ...
               lcA0, lcA0_tv, lcLmd, ny, nRegimes, 'verbose', false);
        bf_ranks(jj) = bj.rank;
        bf_svgap(jj) = bj.sv_gap_ratio;
    end

    output.bf_rank_draws      = bf_ranks;
    output.bf_svgap_draws     = bf_svgap;
    output.bf_full_rank_share = mean(bf_ranks == a_dim);
    output.bf_robust_idx      = bf_idx;

    fprintf('\n--- B&F rank robustness across draws ---\n');
    fprintf('  draws checked            : %d (stride %d)\n', nbf, bf_stride);
    fprintf('  full rank (rank = a = %d) : %.1f%% of draws\n', ...
            a_dim, 100*output.bf_full_rank_share);
    fprintf('  rank range               : [%d, %d]\n', min(bf_ranks), max(bf_ranks));
    fprintf('  median sv-gap ratio      : %.2e\n', median(bf_svgap));
    fprintf('----------------------------------------\n');
end

%% ================================================================
%%  SECTION 14 — FORECAST-ERROR VARIANCE DECOMPOSITION (total-FEVD)
%% ================================================================
% Goal: for each regime, the share of each LABELED shock
% (AD / AS / Conv MP / Reserve-side BS / ON RRP-side BS) in the h-step
% forecast-error variance of each variable, with the denominator being the
% TOTAL forecast MSE (all ny structural shocks, incl. unlabeled), from the
% reduced-form MA representation.
%
% Non-tvA: A0 is common across regimes, so the reduced form By and the Wold
% MA matrices Theta are computed ONCE per draw; only the 1-step covariance
% Sigma_r = A0inv diag(lambda_r) A0inv' varies inside the regime loop.
%
% RESRRP model: TWO balance-sheet shocks among 5 labeled shocks —
%   shock 4 = Reserve-side BS, shock 5 = ON RRP-side BS.
% Regimes (short, 4): R1 pre-2013-09 (reserve-QE) | R2 2013-09..2020-03 (incl
% QT1) | R3 2020-04..2022-02 (COVID QE) | R4 2022-03..end (QT2).
%   QE regimes default [1 3], QT regimes default [2 4]
%   (override via options.fevd_QE_regimes / options.fevd_QT_regimes;
%    balance-sheet shock list via options.fevd_bs_shocks, default [4 5]).
%
% Stored (compact summaries):
%   output.fevd_horizons / fevd_share_regime / fevd_labeled_total_regime
%   output.fevd_QE / output.fevd_QT  (all labeled shocks, ny x nH x Ksh x 3)
%
% NB (MU / Caldara fiscal runs): this section is RESRRP-specific in its
% QE/QT labels; it is gated below on options.fevd_enable (default: true
% only when isfield(options,'fevd_enable')), so fiscal drivers simply leave
% the field unset to skip the console footer.

fevd_enable = get_opt(options, 'fevd_enable', false);

if fevd_enable && use_signrestrictions && nRegimes > 1 && isfield(output, 'sign_irf_regime')

    fevd_horizons = [0 1 2 5 10 15 20 25 30 35 40 45 50 55 59];
    fevd_horizons = fevd_horizons(fevd_horizons < irf_horizon);
    nH   = numel(fevd_horizons);
    Ksh  = numel(accepted_shock);
    nD   = collected_draws;

    % broadcast variables (reuse SECTION 12's if present)
    if ~exist('xA_all', 'var')
        xA_all = zeros(nD, length(draw_x0_all{1}));
        for j = 1:nD, xA_all(j, :) = draw_x0_all{j}'; end
    end
    if ~exist('draw_Phi_mat', 'var')
        draw_Phi_mat = cat(3, draw_Phi_all{1:nD});
    end
    draw_Lambda_mat_local = cat(3, draw_Lambda_all{1:nD});

    Hmax = max(fevd_horizons);   % MA matrices Theta_0..Theta_Hmax

    % storage: per-draw shares (ny, nH, Ksh, nRegimes)
    share_all = zeros(ny, nH, Ksh, nRegimes, nD);

    irf_reg = output.sign_irf_regime;   % (ny, irf_horizon, Ksh, nD, nRegimes)

    parfor j = 1:nD
        % common structural lag block for this draw
        Phi   = draw_Phi_mat(:, :, j);
        Btrim = Phi(1:lags*ny, :);
        Aplus = zeros(ny, ny, lags);
        for iLag = 1:lags
            rows = (iLag-1)*ny + (1:ny);
            Aplus(:, :, iLag) = Btrim(rows, :)';
        end

        % common contemporaneous matrix for this draw
        A0j = zeros(ny);
        A0j(lcA0) = xA_all(j, 1:sum(lcA0(:)));
        A0inv_j = inv(A0j);

        % reduced form and Wold MA matrices Theta_0..Theta_Hmax (regime-invariant)
        By_j = zeros(ny, ny, lags);
        for iLag = 1:lags
            By_j(:, :, iLag) = A0j \ Aplus(:, :, iLag);
        end
        Theta_full = impulsdtrf(By_j, eye(ny), Hmax+1);  % (ny, Hmax+1, ny)
        Theta = zeros(ny, ny, Hmax+1);
        for h = 1:(Hmax+1)
            Theta(:, :, h) = squeeze(Theta_full(:, h, :));
        end

        share_j = zeros(ny, nH, Ksh, nRegimes);

        for r = 1:nRegimes
            % reduced-form 1-step covariance in regime r
            lam_r   = max(draw_Lambda_mat_local(:, r, j), 0);
            Sigma_r = A0inv_j * diag(lam_r) * A0inv_j';

            % cumulative total MSE per variable across horizons
            mse_cum = zeros(ny, Hmax+1);
            running = zeros(ny, 1);
            for l = 1:(Hmax+1)
                Th = Theta(:, :, l);
                running = running + diag(Th * Sigma_r * Th');
                mse_cum(:, l) = running;
            end

            % numerator: cumulative squared structural IRF of each labeled shock
            for kk = 1:Ksh
                irf_k = squeeze(irf_reg(:, :, kk, j, r));   % (ny, irf_horizon)
                sq    = irf_k.^2;
                csq   = cumsum(sq, 2);
                for hi = 1:nH
                    hh  = fevd_horizons(hi);
                    num = csq(:, hh+1);
                    den = mse_cum(:, hh+1);
                    den(den <= 0) = NaN;
                    share_j(:, hi, kk, r) = num ./ den;
                end
            end
        end

        share_all(:, :, :, :, j) = share_j;
    end

    % --- summarize across draws: median / 16 / 84 ---
    %  weight the draw dimension when the omega weighting is active;
    %        otherwise keep median/quantile EXACTLY (byte-identical).
    fevd_share_regime = zeros(ny, nH, Ksh, nRegimes, 3);
    if use_draw_weights   % [W]
        w14 = draw_weight_all(1:collected_draws);
        fevd_share_regime(:,:,:,:,1) = local_wmedian_dim(share_all, w14, 5);
        fevd_share_regime(:,:,:,:,2) = local_wquantile_dim(share_all, w14, 0.16, 5);
        fevd_share_regime(:,:,:,:,3) = local_wquantile_dim(share_all, w14, 0.84, 5);
    else
        fevd_share_regime(:,:,:,:,1) = median(share_all, 5, 'omitnan');
        fevd_share_regime(:,:,:,:,2) = quantile(share_all, 0.16, 5);
        fevd_share_regime(:,:,:,:,3) = quantile(share_all, 0.84, 5);
    end

    % --- sum over labeled shocks (how much of total MSE the labeled explain) ---
    labeled_sum_all = squeeze(sum(share_all, 3));   % (ny, nH, nRegimes, nD)
    fevd_labeled_total_regime = zeros(ny, nH, nRegimes, 3);
    if use_draw_weights   % [W]
        fevd_labeled_total_regime(:,:,:,1) = local_wmedian_dim(labeled_sum_all, w14, 4);
        fevd_labeled_total_regime(:,:,:,2) = local_wquantile_dim(labeled_sum_all, w14, 0.16, 4);
        fevd_labeled_total_regime(:,:,:,3) = local_wquantile_dim(labeled_sum_all, w14, 0.84, 4);
    else
        fevd_labeled_total_regime(:,:,:,1) = median(labeled_sum_all, 4, 'omitnan');
        fevd_labeled_total_regime(:,:,:,2) = quantile(labeled_sum_all, 0.16, 4);
        fevd_labeled_total_regime(:,:,:,3) = quantile(labeled_sum_all, 0.84, 4);
    end

    % --- QE / QT regime sets (RESRRP defaults: QE = [1 3], QT = [2 4]) ---
    QE_regs = intersect(get_opt(options, 'fevd_QE_regimes', [1 3]), 1:nRegimes);
    QT_regs = intersect(get_opt(options, 'fevd_QT_regimes', [2 4]), 1:nRegimes);

    if ~isempty(QE_regs)
        qe_all = mean(share_all(:,:,:,QE_regs,:), 4, 'omitnan');  % (ny,nH,Ksh,1,nD)
        qe_all = squeeze(qe_all);                                  % (ny,nH,Ksh,nD)
        if use_draw_weights   % [W] weighted draw dim (4)
            output.fevd_QE = cat(4, local_wmedian_dim(qe_all, w14, 4), ...
                                    local_wquantile_dim(qe_all, w14, 0.16, 4), ...
                                    local_wquantile_dim(qe_all, w14, 0.84, 4));
        else
            output.fevd_QE = cat(4, median(qe_all,4,'omitnan'), ...
                                    quantile(qe_all,0.16,4), ...
                                    quantile(qe_all,0.84,4));          % (ny,nH,Ksh,3)
        end
    end
    if ~isempty(QT_regs)
        qt_all = mean(share_all(:,:,:,QT_regs,:), 4, 'omitnan');
        qt_all = squeeze(qt_all);
        if use_draw_weights   % [W] weighted draw dim (4)
            output.fevd_QT = cat(4, local_wmedian_dim(qt_all, w14, 4), ...
                                    local_wquantile_dim(qt_all, w14, 0.16, 4), ...
                                    local_wquantile_dim(qt_all, w14, 0.84, 4));
        else
            output.fevd_QT = cat(4, median(qt_all,4,'omitnan'), ...
                                    quantile(qt_all,0.16,4), ...
                                    quantile(qt_all,0.84,4));
        end
    end

    output.fevd_horizons             = fevd_horizons;
    output.fevd_share_regime         = fevd_share_regime;
    output.fevd_labeled_total_regime = fevd_labeled_total_regime;
    output.fevd_accepted_shock       = accepted_shock;
    output.fevd_QE_regimes           = QE_regs;
    output.fevd_QT_regimes           = QT_regs;

    % --- console summary: each BS shock's share, QE vs QT, at key horizons ---
    bs_shocks = get_opt(options, 'fevd_bs_shocks', [4 5]);  % reserve-side, ON RRP-side
    rep_h     = intersect([0 10 20 40], fevd_horizons);

    if isfield(output,'fevd_QE') && isfield(output,'fevd_QT')
        for bs = bs_shocks
            bcol = find(accepted_shock == bs, 1);
            if isempty(bcol), continue; end
            switch bs
                case 4, bs_label = 'Reserve-side BS (shock 4)';
                case 5, bs_label = 'ON RRP-side BS (shock 5)';
                otherwise, bs_label = sprintf('shock %d', bs);
            end

            fprintf('\n============================================================\n');
            fprintf('  FEVD: %s share of total forecast variance (median %%)\n', bs_label);
            fprintf('  QE regimes %s vs QT regimes %s\n', mat2str(QE_regs), mat2str(QT_regs));
            fprintf('============================================================\n');
            fprintf('%-14s', 'Variable');
            for hh = rep_h, fprintf('  QE h=%-3d  QT h=%-3d', hh, hh); end
            fprintf('\n');
            fprintf('%s\n', repmat('-', 1, 14 + numel(rep_h)*20));
            for i = 1:ny
                fprintf('%-14s', varnames{i});
                for hh = rep_h
                    hi = find(fevd_horizons == hh, 1);
                    qe = 100*output.fevd_QE(i, hi, bcol, 1);
                    qt = 100*output.fevd_QT(i, hi, bcol, 1);
                    fprintf('  %7.1f  %7.1f', qe, qt);
                end
                fprintf('\n');
            end
            fprintf('============================================================\n');
        end
        fprintf('  (Shares across the labeled shocks need not sum to 100%%;\n');
        fprintf('   the remainder is unlabeled structural variation.)\n');
        fprintf('============================================================\n');
    end
end

%% Small helper used in the printf banner
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

%%  Weighted quantile along an arbitrary draw dimension.
%   Y = local_wquantile_dim(X, w, p, dim) returns the weighted p-quantile of X
%   along dimension `dim`, with per-draw weights w (length size(X,dim)). NaNs in
%   X are omitted per element (matching the 'omitnan' behaviour of the
%   unweighted calls it replaces). The weighting is the standard ARRW importance
%   estimator: sort along dim, accumulate normalized weights, and read the value
%   at cumulative weight p (Type-7-style linear interpolation on the weighted
%   ECDF, so w == const reproduces MATLAB's quantile to ~1e-12).
%
%   Used ONLY when use_draw_weights is true (NRR omega weighting active); the
%   unweighted path keeps the original median/quantile calls untouched.
function Y = local_wquantile_dim(X, w, p, dim)
    w = w(:);
    nd = ndims(X);
    % move draw dim to the front
    perm = [dim, setdiff(1:max(nd,dim), dim)];
    Xp   = permute(X, perm);
    sz   = size(Xp);
    nDr  = sz(1);
    Xp   = reshape(Xp, nDr, []);          % nDraws x M
    M    = size(Xp, 2);
    Yrow = zeros(1, M);

    for c = 1:M
        col = Xp(:, c);
        ok  = ~isnan(col);
        if ~any(ok)
            Yrow(c) = NaN; continue
        end
        v  = col(ok);
        ww = w(ok);
        % drop zero/negative-weight draws so they are exactly inert (a w=0 draw
        % must not shift the weighted ECDF grid).
        pos = ww > 0;
        if ~any(pos)
            Yrow(c) = NaN; continue
        end
        v  = v(pos);
        ww = ww(pos);
        sw = sum(ww);
        if ~(sw > 0)
            Yrow(c) = NaN; continue
        end
        ww = ww / sw;
        [vs, ord] = sort(v);
        ws = ww(ord);
        % Weighted Type-7 plotting positions: generalize the unweighted
        % (i-1)/(n-1) to weighted cumulative mass. cw is the cumulative weight
        % AFTER each order statistic; pp shifts it to the LEADING edge and
        % rescales to [0,1] so that equal weights reproduce MATLAB's quantile
        % (Type 7) EXACTLY (verified to ~1e-12). Single distinct point -> that
        % point.
        if isscalar(vs)
            Yrow(c) = vs(1);
            continue
        end
        cw = cumsum(ws);
        pp = (cw - ws) ./ (1 - ws(end));   % pp(1)=0, pp(end)=1 at equal weights
        % guard against ties in pp (repeated zero-weight or identical mass)
        [pp, iu] = unique(pp, 'stable');
        vs2 = vs(iu);
        if p <= pp(1)
            Yrow(c) = vs2(1);
        elseif p >= pp(end)
            Yrow(c) = vs2(end);
        else
            Yrow(c) = interp1(pp, vs2, p, 'linear');
        end
    end

    outsz      = sz;
    outsz(1)   = 1;
    Y          = reshape(Yrow, outsz);
    Y          = ipermute(Y, perm);
end

%%  Weighted median convenience wrapper (p = 0.5).
function Y = local_wmedian_dim(X, w, dim)
    Y = local_wquantile_dim(X, w, 0.5, dim);
end