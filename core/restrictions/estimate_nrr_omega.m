function [omega_hat, om_diag] = estimate_nrr_omega(Q, cache, var, ...
    sim_dates, sim_regimes, pass_fns, opts)
%ESTIMATE_NRR_OMEGA  (v3) AD-RR (2018) narrative admissibility probability
% at fixed (theta, Q):
%
%     omega(theta, Q) = Pr[ ALL narrative checks pass | fresh model shocks
%                           on the narrative-restricted dates ]
%
% v3 works on the STANDARDIZED shock scale, matching the
% checkers: the labeled shocks are eps_tilde = Q * xi with
% xi = Lambda^{-1/2} * eps, and the checkers standardize udraw_true by
% regime internally. This function therefore simulates xi directly --
% iid across dates and shocks, symmetric, unit scale -- and injects it via
% the checkers' xi_override argument. No lambda enters the simulation
% (the include_lambda option is GONE), and omega is invariant to any common
% scale factor on xi: both the sign event and the |HD| dominance ranking
% are homogeneous of degree zero/one in a common scaling. Only the SHAPE
% of the law matters (student vs gaussian), and only mildly.
%
% INPUTS
%   Q            : rotation of the retained draw
%   cache, var   : theta cache and shock container (sizes only)
%   sim_dates    : nW x 1 dates (rows of udraw_true) to redraw
%   sim_regimes  : nW x 1 regime index per date (kept for API stability and
%                  the moment check; NOT used in the simulation)
%   pass_fns     : cell of handles, each @(xi_sim) -> logical; xi_sim is a
%                  T x ny matrix of standardized shocks, NaN outside the
%                  restricted rows (a tripwire: any checker touching an
%                  unrestricted row fails loudly)
%   opts         : .M (default 2000), .shock_law 'student'|'gaussian'
%                  (default 'student'), .nu (default 5.703233415698),
%                  .moment_check (default false)
%
% OUTPUTS
%   omega_hat : max(n_pass / M, 1/M)
%   om_diag   : .n_pass, .M, .floored, optional moment-check fields

    if nargin < 7 || isempty(opts), opts = struct(); end
    M            = get_field(opts, 'M', 2000);
    shock_law    = get_field(opts, 'shock_law', 'student');
    nu           = get_field(opts, 'nu', 5.703233415698);
    moment_check = get_field(opts, 'moment_check', false);

    om_diag = struct('n_pass', 0, 'M', M, 'floored', false);

    if isempty(pass_fns) || isempty(sim_dates)
        omega_hat = 1; return
    end

    [Tpl, ny] = size(var.udraw_true);

    % ---- dedupe dates, reject regime conflicts ----
    sim_dates   = sim_dates(:);
    sim_regimes = sim_regimes(:);
    [t_uniq, iu] = unique(sim_dates, 'stable');
    r_uniq = sim_regimes(iu);
    for k = 1:numel(t_uniq)
        rk = sim_regimes(sim_dates == t_uniq(k));
        if numel(unique(rk)) > 1
            error('estimate_nrr_omega: date %d maps to two regimes.', t_uniq(k));
        end
    end
    if max(t_uniq) > Tpl || min(t_uniq) < 1
        error('estimate_nrr_omega: narrative date outside udraw_true range.');
    end
    nW = numel(t_uniq);

    use_t = strcmpi(shock_law, 'student');
    if use_t && nu <= 2
        error('estimate_nrr_omega: nu <= 2 has no finite variance.');
    end

    % ---- Monte Carlo over fresh STANDARDIZED shocks ----
    n_pass = 0;
    xi_sim = NaN(Tpl, ny);           % NaN tripwire outside restricted rows
    XI = randn(nW, ny, M);           % [PERF] all randoms in ONE call
    if use_t                         % raw t_nu; common scale is irrelevant
        XI = XI ./ sqrt(chi2rnd(nu, nW, ny, M) / nu);
    end
    for m = 1:M
        xi_sim(t_uniq, :) = XI(:, :, m);

        ok = true;
        for c = 1:numel(pass_fns)
            if ~pass_fns{c}(xi_sim), ok = false; break; end
        end
        n_pass = n_pass + ok;
    end

    om_diag.n_pass = n_pass;
    omega_hat      = n_pass / M;
    if omega_hat < 1 / M
        omega_hat       = 1 / M;
        om_diag.floored = true;
    end

    % ---- one-time shape check on realized standardized shocks ----
    % Standardize the restricted rows by their own regime and inspect the
    % law: kurtosis ~ 3 -> gaussian; >> 3 -> student. Common scale does not
    % affect omega, so only the shape matters here.
    if moment_check
        lamW   = sqrt(max(cache.lambda(:, r_uniq), 0))';   % nW x ny
        xi_real = var.udraw_true(t_uniq, :) ./ lamW;
        om_diag.xi_std      = std(xi_real, 0, 1);
        om_diag.xi_kurtosis = kurtosis(xi_real, 0, 1);
        om_diag.n_rows      = nW;
        om_diag.note = ['Shape check on regime-standardized restricted rows ' ...
            '(needs enough rows to be informative; omega is scale-invariant, ' ...
            'so only kurtosis matters: ~3 -> gaussian, >>3 -> student).'];
    end
end

function v = get_field(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end