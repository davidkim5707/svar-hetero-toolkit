function [log_posterior, log_likelihood, A0_prior_logprob, lplmd, var] = ...
    bvar_posterior_tvA(seedx, y, lags, lcA0, lcLmd, Tsigbrk, oweights, options)
%--------------------------------------------------------------------------
% BVAR_POSTERIOR_TVA   *** PARTIAL-tvA UPDATED ***
%   Time-varying-A0 counterpart of bvar_posterior. A0 is regime-specific, but
%   now supports a PARTIAL split: a chosen subset of free A0 elements (lcA0_tv,
%   via options.tvA_rows / options.lcA0_tv -> resolveTvMask) is regime-specific
%   while the rest is COMMON across regimes. seedx layout:
%
%       seedx = [ a_cm (nA_cm) ; a_tv (nA_tv x nSig, regime-major) ; lmd_free ]
%
%   With lcA0_tv == lcA0 (no partial spec) this reduces exactly to the
%   full-tvA layout [ vec(A0_1); ...; vec(A0_nSig); lmd_free ].
%
%   I/O matches bvar_posterior so eval_posterior can call either with the
%   same 5-output signature. Priors, Minnesota dummies, trend/exogenous
%   handling and the A0 prior family are ALL driven by `options` exactly as
%   in bvar_posterior (no hardcoded hyperparameters).
%
%
%   The regime-mean A0 slice, the freqs-based det term and the [first, mean]
%   prior-regime structure are preserved from the original tvA routine, so
%   the underlying rfvar3 3D-A0 / var.freqs behaviour is unchanged.
%
%   Requires (standalone): resolveTvMask.m, unpackA0_tvA.m
%--------------------------------------------------------------------------

    if nargin < 7, oweights = []; end
    if nargin < 8
        error('bvar_posterior_tvA: options is required (8th argument).');
    end

    [T, ny] = size(y);
    presample = 0;
    lmdmean   = 0;
    nSig      = length(Tsigbrk);    % number of regimes (= nRegimes)

    %% === Partial-tvA partition (single source of truth) ===
    lcA0_tv = resolveTvMask(lcA0, ny, options);   % subset of lcA0 (== lcA0 -> full tvA)
    lcA0_cm = lcA0 & ~lcA0_tv;
    nA_cm   = sum(lcA0_cm(:));
    nA_tv   = sum(lcA0_tv(:));
    nA0blk  = nA_cm + nA_tv * nSig;               % total A0-block length in seedx

    %% === Regime-specific A0 matrices and prior ===
    [A, A0_prior_logprob] = constructA0_tvA(seedx, ny, nSig, lcA0_cm, lcA0_tv, options);

    %% === Lambda matrix and prior (A0 block is nA0blk long) ===
    [lmd, lplmd, isValid] = constructLambda_tvA(seedx, ny, nA0blk, nSig, lcLmd, options);
    if ~isValid
        log_posterior  = 1e5;
        log_likelihood = -1e10;
        var = [];
        return;
    end

    %% === Log-lambda (+ mean column) ===
    llmd = prepareLogLambda(lmd, lmdmean);

    %% === Append regime-mean A0 slice (Sastry convention) ===
    A_aug = cat(3, A, mean(A, 3));   % ny x ny x (nSig+1)

    sigpar.A0      = A_aug;
    sigpar.lmd     = llmd;
    sigpar.Tsigbrk = [Tsigbrk, T];

    %% === Minnesota prior dummies (options-driven) ===
    [ydum, xdum, pbreaks, xdata, nx] = constructPriorDummies(y, lags, ny, options);

    %% === Posterior likelihood ===
    [log_dnsty, var] = computePosteriorLikelihood_tvA(y, ydum, xdum, xdata, ...
        lags, T, ny, nx, pbreaks, sigpar, oweights);

    %% === Prior likelihood ===
    log_prior = computePriorLikelihood_tvA(y, ydum, xdum, xdata, lags, ny, nx, ...
        pbreaks, presample, llmd, A_aug);

    %% === Assemble ===
    log_likelihood = log_dnsty - log_prior;
    log_posterior  = -sum(log_likelihood, 'all') - A0_prior_logprob - lplmd;
end


%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function [A, A0_prior_logprob] = constructA0_tvA(seedx, ny, nSig, lcA0_cm, lcA0_tv, options)
%--------------------------------------------------------------------------
% Build regime-stacked A0 (ny x ny x nSig) under PARTIAL tvA and form the
% A0 log prior. Common free elements (lcA0_cm) are shared across regimes;
% tv free elements (lcA0_tv) are regime-specific. The matrices are built by
% the shared unpacker so every consumer reconstructs A0 identically.
%
% The prior is NESTED (two layers):
%   (i)  per-regime level prior  : sum_s a0LogPrior(A_s)
%        NOTE: the common block is identical across regimes, so for the
%        element-separable settings (0 Normal, 2 Sstar) it is counted nSig
%        times -> an nSig-times tighter level prior on common elements. With
%        T >> 1 and the weak default A0 prior this shift is negligible
%        (likelihood-dominated). For exact cross-model MDD, count the common
%        block once instead.
%   (ii) cross-regime shrinkage  : -1/(2 tau^2) * sum_s ||A_s - Abar||^2 over
%        the TV elements only (common deviations are zero by construction).
%
% Controlled by options.tvA_shrink_tau (default Inf):
%   tau = Inf    -> shrinkage term vanishes; the no-pooling endpoint. For
%                   partial tvA this is the recommended default, since pooling
%                   is now STRUCTURAL (common vs free), not soft.
%   tau finite   -> the regime-specific tv elements are additionally pulled
%                   toward their common cross-regime mean.
%
% Scale standardizing the deviations (per A0 row i), in priority order:
%   options.tvA_shrink_scale  (scalar -> all rows; ny-vector -> per row)
%   sqrt(diag(options.Sstar)) (matches the Sstar A0 prior row scaling)
%   ones(ny,1)                (raw: tau then carries A0-element units)
%--------------------------------------------------------------------------
    lcA0 = lcA0_cm | lcA0_tv;
    A    = unpackA0_tvA(seedx, ny, lcA0, lcA0_tv, nSig);   % ny x ny x nSig

    % --- (i) per-regime level prior (full A_s) ---
    A0_prior_logprob = 0;
    for s = 1:nSig
        A0_prior_logprob = A0_prior_logprob + a0LogPrior(A(:, :, s), ny, options);
    end

    % --- (ii) cross-regime shrinkage on tv elements (common dev = 0) ---
    tau = Inf;
    if isfield(options, 'tvA_shrink_tau') && ~isempty(options.tvA_shrink_tau)
        tau = options.tvA_shrink_tau;
    end

    if nSig >= 2 && isfinite(tau) && tau > 0
        rowscale = resolveShrinkScale(options, ny);    % ny x 1, positive
        Sstd     = repmat(rowscale(:), 1, ny);         % element (i,j) scale = rowscale(i)
        Abar     = mean(A, 3);
        for s = 1:nSig
            dev = (A(:, :, s) - Abar) ./ Sstd;
            A0_prior_logprob = A0_prior_logprob - 0.5 * sum(dev(lcA0_tv).^2) / tau^2;
        end
    end
end


function rowscale = resolveShrinkScale(options, ny)
%--------------------------------------------------------------------------
% Per-row scale used to standardize cross-regime A0 deviations, so that
% options.tvA_shrink_tau is a unitless pooling strength. See constructA0_tvA.
%--------------------------------------------------------------------------
    if isfield(options, 'tvA_shrink_scale') && ~isempty(options.tvA_shrink_scale)
        sc = options.tvA_shrink_scale;
        if isscalar(sc)
            rowscale = sc * ones(ny, 1);
        else
            rowscale = sc(:);
            if numel(rowscale) ~= ny
                error('constructA0_tvA: options.tvA_shrink_scale must be scalar or an ny=%d vector.', ny);
            end
        end
    elseif isfield(options, 'Sstar') && ~isempty(options.Sstar)
        rowscale = sqrt(diag(options.Sstar));
    else
        rowscale = ones(ny, 1);
    end

    rowscale = max(rowscale(:), 1e-8);   % guard against zero/degenerate scale
end


function lp = a0LogPrior(A0, ny, options)
%--------------------------------------------------------------------------
% Single-regime A0 log prior. Same families as bvar_posterior.constructA0,
% so set options.a0_prior_setting identically to your non-tvA run.
%--------------------------------------------------------------------------
    switch options.a0_prior_setting
        case 0   % Normal: A0 ~ N(100*I, 200^2*I)
            A0_mean = 100; A0_std = 200;
            lp = -0.5 * sum(((A0 - eye(ny) * A0_mean).^2) / (A0_std^2), 'all') ...
                 - ny^2 * (log(2 * pi)/2 + log(A0_std));
        case 1   % Improper: p(A0) ∝ |det(A0)|^(-ny)
            lp = -ny * log(abs(det(A0)));
        otherwise % Carriero et al. (2024) via Sstar
            lp = 0;
            for i = 1:ny
                var_i = options.Sstar(i, i);
                lp = lp - 0.5 * sum(A0(i, :).^2) / var_i ...
                        - ny * (0.5 * log(2 * pi) + 0.5 * log(var_i));
            end
    end
end


function [lmd, lplmd, isValid] = constructLambda_tvA(seedx, ny, nA0blk, nSig, lcLmd, options)
%--------------------------------------------------------------------------
% Same constraints / Dirichlet prior as bvar_posterior, but the lambda block
% sits AFTER the (partial-tvA) A0 block of length nA0blk = nA_cm + nA_tv*nSig.
% (Legacy full tvA: nA0blk = nA*nSig.)
%--------------------------------------------------------------------------
    nLmd      = sum(lcLmd(:));
    fix_first = isfield(options, 'fix_first_regime') && logical(options.fix_first_regime);

    if nLmd == 0
        lmd = ones(ny, nSig); lplmd = 0; isValid = true; return;
    end

    lmd = ones(ny, nSig);
    lmd(lcLmd) = seedx(nA0blk + (1:nLmd));

    if fix_first
        [lmd, lplmd, isValid] = applyFixedFirstRegime(lmd, ny, nSig);
    else
        [lmd, lplmd, isValid] = applyFreeRegimes(lmd, ny, nSig);
    end
end


function [lmd, lplmd, isValid] = applyFixedFirstRegime(lmd, ny, nSig)
%--------------------------------------------------------------------------
% First regime fixed to 1, Dirichlet prior on remaining (matches non-tvA).
%--------------------------------------------------------------------------
    lmd(:, 1) = 1;

    if nSig >= 3
        partial = sum(lmd(:, 2:(nSig-1)), 2);
        lmd(:, nSig) = (nSig - 1) - partial;
    else
        lmd(:, nSig) = 1;
    end

    if any(lmd(:, 2:end) <= 0, 'all')
        lplmd = -1e5; isValid = false; return;
    end

    lpL = sum(log(lmd(:, 2:end)) - log(nSig - 1), 1) - gammaln(2 * (nSig - 1));
    lplmd = sum(lpL) - (nSig - 2) * log(nSig - 1);
    isValid = true;
end


function [lmd, lplmd, isValid] = applyFreeRegimes(lmd, ny, nSig)
%--------------------------------------------------------------------------
% Row sums equal nSig, Dirichlet(2) prior on all regimes (matches non-tvA).
%--------------------------------------------------------------------------
    rowSums = sum(lmd(:, 1:(nSig-1)), 2);
    lmd(:, nSig) = nSig - rowSums;

    if any(lmd(:) <= 0)
        lplmd = -1e5; isValid = false; return;
    end

    lpL = sum(log(lmd) - log(nSig), 1) - gammaln(2 * nSig);
    lplmd = sum(lpL) - (nSig - 1) * log(nSig);
    isValid = true;
end


function llmd = prepareLogLambda(lmd, lmdmean)
%--------------------------------------------------------------------------
% Convert lambda to log scale and append mean column (Sims-Zha convention).
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
% Minnesota prior dummy observations (identical to bvar_posterior).
%--------------------------------------------------------------------------
    [T, ~] = size(y);

    mnprior.tight      = options.minn_prior_tau;
    mnprior.decay      = options.minn_prior_decay;
    mnprior.unit_root_ = options.unitroot;

    urprior.lambda = options.minn_prior_lambda;
    urprior.mu     = options.minn_prior_mu;

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
        ydum = []; xdum = []; pbreaks = [];
    else
        [ydum, xdum, pbreaks] = varprior(ny, nx, lags, mnprior, vprior, ...
                                          urprior, ybar, xbar, nstat);
    end
end


function vprior = setupVariancePrior(options, ny, y, lags)
%--------------------------------------------------------------------------
% Variance prior, Sims-Zha convention (identical to bvar_posterior).
%   1. options.vprior_sig                       — user-supplied ny x 1
%   2. options.Sstar + use_Sstar_for_vprior     — sqrt(diag(Sstar))
%   3. default                                  — univariate AR(lags) std
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
% Per-variable AR(lags) residual std, floored at 1e-6 (matches bvar_posterior).
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


function [log_dnsty, var] = computePosteriorLikelihood_tvA(y, ydum, xdum, xdata, ...
    lags, T, ny, nx, pbreaks, sigpar, oweights)
%--------------------------------------------------------------------------
% Posterior log density. detTerm averages log|det A0(regime_t)| over obs via
% var.freqs (regime-varying A0).
%--------------------------------------------------------------------------
    var = rfvar3([y; ydum], lags, [xdata; xdum], [T; T + pbreaks], [], [], sigpar, oweights);

    Tu      = size(var.u, 1);
    lmdllh  = 0.5 * sum(var.lmdseries, 'all');
    detTerm = regimeMeanLogDet(sigpar.A0, var.freqs);

    llh = -0.5 * sum(var.u(:).^2) + Tu * (-ny * log(2 * pi) / 2 + detTerm) + lmdllh;
    nX  = lags * ny + nx;
    log_dnsty = llh + 0.5 * sum(var.logdetxxi, 'all') + ny * nX * log(2 * pi) / 2;
end


function log_prior = computePriorLikelihood_tvA(y, ydum, xdum, xdata, lags, ny, nx, ...
    pbreaks, presample, llmd, A_aug)
%--------------------------------------------------------------------------
% Prior log density. Prior regimes = [first, mean] A0 slices (as in the
% original tvA routine); det averaged over varp.freqs.
%--------------------------------------------------------------------------
    Tp = presample + lags;

    priorSigpar.A0      = cat(3, A_aug(:, :, 1), A_aug(:, :, end));  % [regime 1, mean]
    priorSigpar.lmd     = [llmd(:, 1), llmd(:, end)];
    priorSigpar.Tsigbrk = [0, Tp];

    if isempty(xdata)
        xprior = xdum;
    else
        xprior = [xdata(1:Tp, :); xdum];
    end

    varp = rfvar3([y(1:Tp, :); ydum], lags, xprior, [Tp; Tp + pbreaks], [], [], priorSigpar, []);

    Tup        = size(varp.u, 1);

    if Tup == 0
        % Degenerate prior sub-sample: Tp = lags yields no usable observations
        % (e.g. flat prior, presample = 0, no Minnesota dummies). rfvar3 then
        % returns an empty sample -> smpl = [], freqs = [], logdetxxi = NaN.
        % The data-dependent prior terms all vanish, so return only the
        % dimension constant. This shifts log_posterior by at most an additive
        % constant, which is irrelevant for the mode finder and MH accept/reject.
        nX = lags * ny + nx;
        log_prior = ny * nX * log(2 * pi) / 2;
        return;
    end

    lmdllhp    = 0.5 * sum(varp.lmdseries, 'all');
    detPriorA0 = regimeMeanLogDet(priorSigpar.A0, varp.freqs);

    llhp = -0.5 * sum(varp.u(:).^2) - Tup * (ny * log(2 * pi) / 2 - detPriorA0) + lmdllhp;
    nX   = lags * ny + nx;
    normalizer = 0.5 * sum(varp.logdetxxi, 'all') + ny * nX * log(2 * pi) / 2;
    log_prior  = llhp + normalizer;
end


function d = regimeMeanLogDet(Aset, freqs)
%--------------------------------------------------------------------------
% Mean over observations of log|det A0(regime_t)|. Aset is ny x ny x nSlices,
% freqs is the per-observation regime index returned by rfvar3.
%
% When there are no usable observations (freqs empty -- e.g. the prior
% sub-sample has Tp = lags so rfvar3 yields 0 rows), fall back to the
% unweighted slice mean. This keeps d FINITE: the det term is multiplied by
% Tup = 0 in that case, and a NaN here would poison the result via 0*NaN.
%--------------------------------------------------------------------------
    nS   = size(Aset, 3);
    dets = zeros(nS, 1);
    for s = 1:nS
        dets(s) = log(abs(det(Aset(:, :, s))));
    end
    if isempty(freqs)
        d = mean(dets);
    else
        d = mean(dets(freqs));
    end
end