function [log_posterior, log_likelihood, A0_prior_logprob, lplmd, var] = ...
    bvar_posterior_thetero(seedx, y, lags, lcA0, lcLmd, Tsigbrk, oweights, options)
%--------------------------------------------------------------------------
% BVAR_POSTERIOR_THETERO
%   Threshold-augmented heteroskedastic BVAR posterior (option i).
%
%   Identical to BVAR_POSTERIOR, plus a Qian-Marin-Veiga (2026) style level
%   shift in the (orthogonalized) structural log-variances:
%
%       log var_i(t) = log lambda_{i, s(t)} + delta_i * d_t,
%
%   where s(t) is the (unchanged) calendar regime, d_t in {0,1} is a
%   pre-determined stress indicator (options.d_t, length nobs), and delta_i
%   is a free shift for each structural shock listed in options.delta_shocks.
%   delta is appended to seedx AFTER the A0 and lambda blocks:
%        seedx = [ A0(lcA0) ; lambda(lcLmd) ; delta(delta_shocks) ].
%   With options.thetero off / no delta_shocks this reduces to bvar_posterior.
%
%   delta enters the LIKELIHOOD only (through sigpar.delta / sigpar.d_t, which
%   rfvar3 folds into lmdseries). It does NOT enter the prior-normalisation
%   block (computePriorLikelihood), which is correct: delta is a sample-period
%   variance feature, not part of the Minnesota dummy normaliser.
%
% Extra options fields (set in the driver when options.thetero == 1):
%   options.d_t             : nobs x 1 binary stress indicator (pre-determined)
%   options.delta_shocks    : row/col vector of structural-shock indices (1..ny)
%                             carrying a free stress shift (e.g. 1:ny, or a subset)
%   options.delta_prior_var : scalar c3 in delta_i ~ N(0, c3)
%
% Outputs: same 5 as bvar_posterior.
%--------------------------------------------------------------------------

    %% === Input Validation ===
    if nargin < 7
        oweights = [];
    end

    %% === Setup ===
    [T, ny] = size(y);
    presample = 0;
    lmdmean = 0;

    %% === Construct A0 Matrix and Prior ===
    [A0, A0_prior_logprob] = constructA0(seedx, ny, lcA0, options);

    %% === Construct Lambda Matrix and Prior ===
    [lmd, lplmd, isValid] = constructLambda(seedx, ny, lcA0, lcLmd, Tsigbrk, options);

    if ~isValid
        log_posterior = 1e5;
        log_likelihood = -1e10;
        var = [];
        return;
    end

    %% === Prepare Lambda for Likelihood ===
    llmd = prepareLogLambda(lmd, lmdmean);

    %% === Threshold stress shift delta (option i) ===
    [delta_full, delta_prior_logprob] = ...
        constructDeltaShift(seedx, ny, lcA0, lcLmd, options);

    %% === Setup Sigma Parameters ===
    sigpar.A0      = A0;
    sigpar.lmd     = llmd;
    sigpar.Tsigbrk = [Tsigbrk, T];
    sigpar.delta   = delta_full;        % ny x 1, zeros on non-loaded shocks
    sigpar.d_t     = options.d_t;       % nobs x 1 binary indicator

    % transmission threshold (tthr): pass the high-inflation indicator + which
    % lags are state-dependent into rfvar3 (it interacts those lags). Orthogonal
    % to delta (variance) and to the hetero ID of B.
    if isfield(options,'tthr') && ~isempty(options.tthr) && options.tthr
        sigpar.d_inf   = options.d_inf;
        sigpar.sd_lags = options.sd_lags;
        % optional shrinkage prior on the increment Gamma (Gamma ~ N(0,c^2) in
        % the A0 frame); rfvar3 adds Theil dummy rows. [] / Inf -> flat (no-op).
        if isfield(options,'gamma_prior_sd') && ~isempty(options.gamma_prior_sd)
            sigpar.gamma_sd = options.gamma_prior_sd;
        end
    end

    %% === Construct Prior Dummy Observations ===
    [ydum, xdum, pbreaks, xdata, nx] = constructPriorDummies(y, lags, ny, options);

    %% === Compute Posterior Likelihood ===
    [log_dnsty, var] = computePosteriorLikelihood(y, ydum, xdum, xdata, ...
        lags, T, ny, nx, pbreaks, sigpar, oweights, options);

    %% === Compute Prior Likelihood (delta deliberately excluded here) ===
    log_prior = computePriorLikelihood(y, ydum, xdum, xdata, lags, ny, nx, ...
        pbreaks, presample, llmd, sigpar, A0);

    %% === Combine into Final Posterior ===
    log_likelihood = log_dnsty - log_prior;
    log_posterior = -sum(log_likelihood, 'all') - A0_prior_logprob ...
                    - lplmd - delta_prior_logprob;
end


%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function [delta_full, delta_prior_logprob] = ...
    constructDeltaShift(seedx, ny, lcA0, lcLmd, options)
%--------------------------------------------------------------------------
% Unpack the stress-shift block from seedx (after A0 and lambda) and compute
% its Normal(0, c3) prior. delta_full is ny x 1 with the free shifts placed
% on the rows in options.delta_shocks and zeros elsewhere.
%--------------------------------------------------------------------------
    nA   = sum(lcA0(:));
    nLmd = sum(lcLmd(:));

    ds = options.delta_shocks(:)';          % structural shocks carrying a shift
    nd = numel(ds);

    delta_vec = seedx(nA + nLmd + (1:nd));   % free shift params (real)

    delta_full = zeros(ny, 1);
    delta_full(ds) = delta_vec;

    % delta_i ~ N(0, c3)  (Qian's Delta_d prior)
    c3 = options.delta_prior_var;
    delta_prior_logprob = -0.5 * sum(delta_vec.^2) / c3 ...
                          - nd * (0.5 * log(2 * pi) + 0.5 * log(c3));
end


function [A0, A0_prior_logprob] = constructA0(seedx, ny, lcA0, options)
%--------------------------------------------------------------------------
% Construct A0 matrix and compute its prior log probability
%--------------------------------------------------------------------------
    A0 = zeros(ny, ny);
    nA = sum(lcA0(:));
    A0(lcA0) = seedx(1:nA);

    switch options.a0_prior_setting
        case 0
            % Normal prior: A0 ~ N(100*I, 200^2*I)
            A0_mean = 100;
            A0_std = 200;
            A0_prior_logprob = -0.5 * sum(((A0 - eye(ny) * A0_mean).^2) / (A0_std^2), 'all') ...
                               - ny^2 * (log(2 * pi)/2 + log(A0_std));
        case 1
            % Improper prior: p(A0) proportional to |det(A0)|^(-ny)
            A0_prior_logprob = -ny * log(abs(det(A0)));
        otherwise
            % Carriero et al. (2024) style prior using Sstar
            A0_prior_logprob = computeCarrieroPrior(A0, ny, options.Sstar);
    end
end


function A0_prior_logprob = computeCarrieroPrior(A0, ny, Sstar)
%--------------------------------------------------------------------------
% Compute Carriero et al. (2024) style prior for A0
%--------------------------------------------------------------------------
    A0_prior_logprob = 0;
    for i = 1:ny
        var_i = Sstar(i, i);
        A0_prior_logprob = A0_prior_logprob ...
            - 0.5 * sum(A0(i, :).^2) / var_i ...
            - ny * (0.5 * log(2 * pi) + 0.5 * log(var_i));
    end
end


function [lmd, lplmd, isValid] = constructLambda(seedx, ny, lcA0, lcLmd, Tsigbrk, options)
%--------------------------------------------------------------------------
% Construct lambda matrix (regime-specific variances) and compute prior
%--------------------------------------------------------------------------
    nA = sum(lcA0(:));
    nLmd = sum(lcLmd(:));
    nSig = length(Tsigbrk);

    fix_first = isfield(options, 'fix_first_regime') && logical(options.fix_first_regime);

    if nLmd == 0
        % No heteroskedasticity
        lmd = ones(ny, nSig);
        lplmd = 0;
        isValid = true;
        return;
    end

    % Initialize lambda matrix
    lmd = ones(ny, nSig);
    lmd(lcLmd) = seedx(nA + (1:nLmd));

    if fix_first
        [lmd, lplmd, isValid] = applyFixedFirstRegime(lmd, ny, nSig);
    else
        [lmd, lplmd, isValid] = applyFreeRegimes(lmd, ny, nSig);
    end
end


function [lmd, lplmd, isValid] = applyFixedFirstRegime(lmd, ny, nSig)
%--------------------------------------------------------------------------
% Apply constraint: first regime fixed to 1, Dirichlet prior on remaining
%--------------------------------------------------------------------------
    lmd(:, 1) = 1;

    if nSig >= 3
        partial = sum(lmd(:, 2:(nSig-1)), 2);
        lmd(:, nSig) = (nSig - 1) - partial;
    else
        lmd(:, nSig) = 1;
    end

    if any(lmd(:, 2:end) <= 0, 'all')
        lplmd = -1e5;
        isValid = false;
        return;
    end

    lpL = sum(log(lmd(:, 2:end)) - log(nSig - 1), 1) - gammaln(2 * (nSig - 1));
    lplmd = sum(lpL) - (nSig - 2) * log(nSig - 1);
    isValid = true;
end


function [lmd, lplmd, isValid] = applyFreeRegimes(lmd, ny, nSig)
%--------------------------------------------------------------------------
% Apply constraint: row sums equal nSig, Dirichlet prior on all regimes
%--------------------------------------------------------------------------
    rowSums = sum(lmd(:, 1:(nSig-1)), 2);
    lmd(:, nSig) = nSig - rowSums;

    if any(lmd(:) <= 0)
        lplmd = -1e5;
        isValid = false;
        return;
    end

    lpL = sum(log(lmd) - log(nSig), 1) - gammaln(2 * nSig);
    lplmd = sum(lpL) - (nSig - 1) * log(nSig);
    isValid = true;
end


function llmd = prepareLogLambda(lmd, lmdmean)
%--------------------------------------------------------------------------
% Convert lambda to log scale and optionally append mean
%--------------------------------------------------------------------------
    llmd = -log(lmd);

    if isvector(llmd)
        lmdbar = llmd;
    elseif lmdmean
        lmdbar = mean(llmd, 2);
    else
        lmdbar = zeros(size(llmd, 1), 1);
    end

    llmd = [llmd, lmdbar];
end


function [ydum, xdum, pbreaks, xdata, nx] = constructPriorDummies(y, lags, ny, options)
%--------------------------------------------------------------------------
% Construct Minnesota prior dummy observations
%--------------------------------------------------------------------------
    [T, ~] = size(y);

    mnprior.tight = options.minn_prior_tau;
    mnprior.decay = options.minn_prior_decay;
    mnprior.unit_root_ = options.unitroot;

    urprior.lambda = options.minn_prior_lambda;
    urprior.mu = options.minn_prior_mu;

    vprior = setupVariancePrior(options, ny, y, lags);

    xdata = [];
    if options.constant
        xdata = [xdata, ones(T, 1)];
    end

    if isfield(options, 'timetrend') && ~isempty(options.timetrend) && options.timetrend >= 1
        trend_lin = (1:T)' / T;
        xdata = [xdata, trend_lin];
        if options.timetrend >= 2
            xdata = [xdata, trend_lin.^2];
        end
    end

    if isfield(options, 'exogenous') && ~isempty(options.exogenous)
        xdata = [xdata, options.exogenous];
    end

    nx = size(xdata, 2);
    if nx > 0
        xbar = mean(xdata(1:lags, :), 1);
    else
        xbar = [];
    end

    ybar = mean(y(1:lags, :), 1);

    if isfield(options, 'unitroot') && ~isempty(options.unitroot)
        nstat = options.unitroot;
    else
        nstat = ones(ny, 1);
    end

    if isempty(mnprior.tight) && isempty(mnprior.decay) && ...
       isempty(urprior.lambda) && isempty(urprior.mu) && ...
       (isempty(vprior.w) || vprior.w == 0)
        ydum = [];
        xdum = [];
        pbreaks = [];
    else
        [ydum, xdum, pbreaks] = varprior(ny, nx, lags, mnprior, vprior, ...
                                          urprior, ybar, xbar, nstat);
    end
end


function vprior = setupVariancePrior(options, ny, y, lags)
%--------------------------------------------------------------------------
% Setup variance prior (Sims-Zha convention).
%--------------------------------------------------------------------------
    if isfield(options, 'vprior_sig') && ~isempty(options.vprior_sig)
        sig = options.vprior_sig(:);
        if numel(sig) ~= ny
            error('setupVariancePrior: options.vprior_sig must have ny=%d entries.', ny);
        end
        vprior.sig = sig;
    elseif isfield(options, 'use_Sstar_for_vprior') && options.use_Sstar_for_vprior ...
            && isfield(options, 'Sstar') && ~isempty(options.Sstar)
        vprior.sig = sqrt(diag(options.Sstar));
    else
        vprior.sig = computeUnivariateARStd(y, lags);
    end
    vprior.w = options.minn_prior_omega;
end


function sigvec = computeUnivariateARStd(y, lags)
%--------------------------------------------------------------------------
% Per-variable AR(p) residual std (Sims-Zha convention).
%--------------------------------------------------------------------------
    [T, ny] = size(y);
    sigvec  = zeros(ny, 1);
    Teff    = T - lags;

    for j = 1:ny
        Yj = y(lags+1:T, j);
        Xj = ones(Teff, 1);
        for l = 1:lags
            Xj = [Xj, y(lags+1-l : T-l, j)]; %#ok<AGROW>
        end
        bj  = Xj \ Yj;
        res = Yj - Xj * bj;
        sigvec(j) = max(std(res), 1e-6);
    end
end


function [log_dnsty, var] = computePosteriorLikelihood(y, ydum, xdum, xdata, ...
    lags, T, ny, nx, pbreaks, sigpar, oweights, options)
%--------------------------------------------------------------------------
% Compute posterior log density using rfvar3
%--------------------------------------------------------------------------
    var = rfvar3([y; ydum], lags, [xdata; xdum], [T; T + pbreaks], [], [], sigpar, oweights);

    Tu = size(var.u, 1);
    lmdllh = 0.5 * sum(var.lmdseries, 'all');
    detTerm = log(abs(det(sigpar.A0)));

    llh = -0.5 * sum(var.u(:).^2) + Tu * (-ny * log(2 * pi) / 2 + detTerm) + lmdllh;

    nX = lags * ny + nx;

    log_dnsty = llh + 0.5 * sum(var.logdetxxi, 'all') + ny * nX * log(2 * pi) / 2;
end


function log_prior = computePriorLikelihood(y, ydum, xdum, xdata, lags, ny, nx, ...
    pbreaks, presample, llmd, sigpar, A0)
%--------------------------------------------------------------------------
% Compute prior log density (presample + dummy observations).
% NOTE: priorSigpar carries NO delta/d_t, so the stress shift correctly does
% not enter the Minnesota prior normaliser.
%--------------------------------------------------------------------------
    Tp = presample + lags;

    priorSigpar.A0 = A0;
    priorSigpar.lmd = [llmd(:, 1), llmd(:, end)];
    priorSigpar.Tsigbrk = [0, Tp];

    if isempty(xdata)
        xprior = xdum;
    else
        xprior = [xdata(1:Tp, :); xdum];
    end

    varp = rfvar3([y(1:Tp, :); ydum], lags, xprior, [Tp; Tp + pbreaks], [], [], priorSigpar, []);

    Tup = size(varp.u, 1);
    lmdllhp = 0.5 * sum(varp.lmdseries, 'all');
    detPriorA0 = log(abs(det(A0)));

    llhp = -0.5 * sum(varp.u(:).^2) - Tup * (ny * log(2 * pi) / 2 - detPriorA0) + lmdllhp;

    nX = lags * ny + nx;

    normalizer = 0.5 * sum(varp.logdetxxi, 'all') + ny * nX * log(2 * pi) / 2;
    log_prior = llhp + normalizer;
end