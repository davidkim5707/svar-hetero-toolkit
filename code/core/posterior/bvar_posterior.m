function [log_posterior, log_likelihood, A0_prior_logprob, lplmd, var] = ...
    bvar_posterior(seedx, y, lags, lcA0, lcLmd, Tsigbrk, oweights, options)
%--------------------------------------------------------------------------
% BVAR_POSTERIOR
%   Computes the posterior density for a Bayesian SVAR with heteroskedasticity
%
% Inputs:
%   seedx    : parameter vector (A0 elements followed by lambda elements)
%   y        : T x ny matrix of endogenous variables
%   lags     : number of VAR lags
%   lcA0     : ny x ny logical mask for free parameters in A0
%   lcLmd    : ny x nSig logical mask for free parameters in lambda
%   Tsigbrk  : vector of regime break points
%   oweights : observation weights (optional)
%   options  : structure with prior settings and hyperparameters
%
%   TTHR: if options.tthr is on, the high-inflation indicator
%       options.d_inf and the state-dependent lags options.sd_lags are passed
%       into rfvar3 via sigpar, which interacts those lags (state-dependent
%       conditional-mean dynamics A_j). This is regressor-side only and leaves
%       the heteroskedasticity identification of B = A0^{-1} untouched.
%
% Outputs:
%   log_posterior    : negative log posterior (for minimization)
%   log_likelihood   : log likelihood value
%   A0_prior_logprob : log prior probability of A0
%   lplmd            : log prior probability of lambda
%   var              : structure with VAR estimation results
%
% Based on R code originally written by Karthik Sastry (SVAR toolkit)
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

    %% === Setup Sigma Parameters ===
    sigpar.A0 = A0;
    sigpar.lmd = llmd;
    sigpar.Tsigbrk = [Tsigbrk, T];

    % transmission threshold (tthr): pass high-inflation indicator + state-
    % dependent lags into rfvar3 (regressor-side; hetero ID of B untouched).
    if isfield(options,'tthr') && ~isempty(options.tthr) && options.tthr
        sigpar.d_inf   = options.d_inf;
        sigpar.sd_lags = options.sd_lags;
        % optional shrinkage prior on the increment Gamma (Gamma ~ N(0,c^2) in
        % the A0 frame); rfvar3 adds Theil dummy rows. [] / Inf -> flat (no-op).
        if isfield(options,'gamma_prior_sd') && ~isempty(options.gamma_prior_sd)
            sigpar.gamma_sd = options.gamma_prior_sd;
        end
        % optional restriction of state-dependence to chosen EQUATIONS (rows),
        % e.g. [GDP GDPDEF]; others are forced state-invariant by rfvar3.
        if isfield(options,'sd_eqs') && ~isempty(options.sd_eqs)
            sigpar.sd_eqs = options.sd_eqs;
        end
    end

    %% === Construct Prior Dummy Observations ===
    [ydum, xdum, pbreaks, xdata, nx] = constructPriorDummies(y, lags, ny, options);

    %% === Compute Posterior Likelihood ===
    [log_dnsty, var] = computePosteriorLikelihood(y, ydum, xdum, xdata, ...
        lags, T, ny, nx, pbreaks, sigpar, oweights, options);

    %% === Compute Prior Likelihood ===
    log_prior = computePriorLikelihood(y, ydum, xdum, xdata, lags, ny, nx, ...
        pbreaks, presample, llmd, sigpar, A0);

    %% === Combine into Final Posterior ===
    log_likelihood = log_dnsty - log_prior;
    log_posterior = -sum(log_likelihood, 'all') - A0_prior_logprob - lplmd;
end


%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

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

    % Enforce sum constraint on remaining regimes
    if nSig >= 3
        partial = sum(lmd(:, 2:(nSig-1)), 2);
        lmd(:, nSig) = (nSig - 1) - partial;
    else
        lmd(:, nSig) = 1;
    end

    % Positivity check
    if any(lmd(:, 2:end) <= 0, 'all')
        %disp('Warning: Non-positive lambda created with fix_first_regime; rejecting draw.');
        lplmd = -1e5;
        isValid = false;
        return;
    end

    % Dirichlet(2) prior on regimes 2:end
    lpL = sum(log(lmd(:, 2:end)) - log(nSig - 1), 1) - gammaln(2 * (nSig - 1));
    lplmd = sum(lpL) - (nSig - 2) * log(nSig - 1);
    isValid = true;
end


function [lmd, lplmd, isValid] = applyFreeRegimes(lmd, ny, nSig)
%--------------------------------------------------------------------------
% Apply constraint: row sums equal nSig, Dirichlet prior on all regimes
%--------------------------------------------------------------------------
    % Enforce arithmetic-average constraint
    rowSums = sum(lmd(:, 1:(nSig-1)), 2);
    lmd(:, nSig) = nSig - rowSums;

    % Positivity check
    if any(lmd(:) <= 0)
        % disp('Warning: Non-positive lambda; rejecting draw.');
        lplmd = -1e5;
        isValid = false;
        return;
    end

    % Dirichlet(2) prior
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
%
%   setupVariancePrior now takes y and lags so that vprior.sig defaults to
%   the univariate AR(lags) residual std per variable (Sims-Zha convention).
%   The earlier default of 0.01 was wildly off-scale for most real datasets
%   and made the Minnesota prior essentially flat.
%--------------------------------------------------------------------------
    [T, ~] = size(y);

    % Setup hyperparameters
    mnprior.tight = options.minn_prior_tau;
    mnprior.decay = options.minn_prior_decay;
    mnprior.unit_root_ = options.unitroot;

    urprior.lambda = options.minn_prior_lambda;
    urprior.mu = options.minn_prior_mu;

    % vprior.sig now scales with each variable's residual std
    vprior = setupVariancePrior(options, ny, y, lags);

    % Setup exogenous data: [constant, trend, trend^2, user exogenous]
    xdata = [];
    if options.constant
        xdata = [xdata, ones(T, 1)];
    end

    % Time trend support:
    %   timetrend = 0 : no trend (default)
    %   timetrend = 1 : linear trend
    %   timetrend = 2 : linear + quadratic trend
    if isfield(options, 'timetrend') && ~isempty(options.timetrend) && options.timetrend >= 1
        trend_lin = (1:T)' / T;
        xdata = [xdata, trend_lin];
        if options.timetrend >= 2
            xdata = [xdata, trend_lin.^2];
        end
    end

    % User-supplied exogenous variables
    if isfield(options, 'exogenous') && ~isempty(options.exogenous)
        xdata = [xdata, options.exogenous];
    end

    nx = size(xdata, 2);
    if nx > 0
        xbar = mean(xdata(1:lags, :), 1);
    else
        xbar = [];
    end

    % Mean of first lags observations
    ybar = mean(y(1:lags, :), 1);

    % Use the user-defined unitroot vector instead of hardcoded ones
    if isfield(options, 'unitroot') && ~isempty(options.unitroot)
        nstat = options.unitroot; 
    else
        nstat = ones(ny, 1); % Default fallback if you forget to set it
    end

    % Generate dummy observations
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
%
%   Default vprior.sig is now computed from univariate AR(p) residual stds
%   per variable. This makes the Minnesota prior scale-invariant: the dummy
%   magnitude tight*sig(j) tracks each variable's natural scale, so prior
%   precision on coefficients is (tight)^2 regardless of variable units.
%
%   Override options (priority order):
%     1. options.vprior_sig                  — user-supplied ny×1 vector
%     2. options.Sstar + options.use_Sstar_for_vprior — sqrt(diag(Sstar))
%     3. (default)                            — univariate AR(lags) residual std
%
%   Earlier default was 0.01*ones(ny,1), which gave essentially flat prior
%   on datasets with residual stds far from 0.01.
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
        % Default: univariate AR(lags) residual std per variable
        vprior.sig = computeUnivariateARStd(y, lags);
    end
    vprior.w = options.minn_prior_omega;
end


function sigvec = computeUnivariateARStd(y, lags)
%--------------------------------------------------------------------------
% COMPUTEUNIVARIATEARSTD  Per-variable AR(p) residual std (Sims-Zha convention).
%
%   For each column j of y, fit a univariate AR(lags) with constant by OLS
%   and return std of residuals. Used to scale Minnesota dummies so that
%   prior precision on coefficients is invariant across variables with
%   different scales.
%
%   Floor at 1e-6 to avoid degenerate sig=0 on near-zero columns.
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

    % nx already includes constant + trend + exogenous columns
    nX = lags * ny + nx;

    log_dnsty = llh + 0.5 * sum(var.logdetxxi, 'all') + ny * nX * log(2 * pi) / 2;
end


function log_prior = computePriorLikelihood(y, ydum, xdum, xdata, lags, ny, nx, ...
    pbreaks, presample, llmd, sigpar, A0)
%--------------------------------------------------------------------------
% Compute prior log density (presample + dummy observations)
%--------------------------------------------------------------------------
    Tp = presample + lags;

    % Prior sigma parameters
    priorSigpar.A0 = A0;
    priorSigpar.lmd = [llmd(:, 1), llmd(:, end)];
    priorSigpar.Tsigbrk = [0, Tp];

    % Prior x data
    if isempty(xdata)
        xprior = xdum;
    else
        xprior = [xdata(1:Tp, :); xdum];
    end

    % Compute prior VAR
    varp = rfvar3([y(1:Tp, :); ydum], lags, xprior, [Tp; Tp + pbreaks], [], [], priorSigpar, []);

    % Prior likelihood components
    Tup = size(varp.u, 1);
    lmdllhp = 0.5 * sum(varp.lmdseries, 'all');
    detPriorA0 = log(abs(det(A0)));

    llhp = -0.5 * sum(varp.u(:).^2) - Tup * (ny * log(2 * pi) / 2 - detPriorA0) + lmdllhp;

    % nx already includes constant + trend + exogenous columns
    nX = lags * ny + nx;

    normalizer = 0.5 * sum(varp.logdetxxi, 'all') + ny * nX * log(2 * pi) / 2;
    log_prior = llhp + normalizer;
end