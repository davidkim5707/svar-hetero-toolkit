function [lh, likelihood, A0_prior, lambda_prior, var] = ...
    eval_posterior(x, y, lags, lcA0, lcLmd, Tsigbrk, delta0, tvA, options)
%==========================================================================
% Evaluate the log posterior and associated output structure for a BVAR model.
%
% Selects the posterior routine by heteroskedasticity (lcLmd), whether A0 is
% time-varying (tvA), and whether the threshold stress shift is on
% (options.thetero). The thetero branch mirrors the non-tvA hetero branch but
% calls bvar_posterior_thetero, which adds a Qian-style additive shift
% delta_i*d_t to the (orthogonalized) structural log-variances.
%
% Inputs:
%   x        : current parameter draw (vectorized; tvA stacks A0 per regime;
%              thetero appends delta after the A0 and lambda blocks)
%   y        : T x ny data matrix
%   lags     : VAR lags
%   lcA0     : ny x ny logical mask for free elements of A0
%   lcLmd    : ny x nSig logical mask for regime-specific variances
%   Tsigbrk  : regime breakpoints
%   delta0   : log-volatility draw (used when lcLmd is nonzero)
%   tvA      : true if A0 is time-varying across regimes
%   options  : prior settings / hyperparameters (options.thetero toggles the
%              threshold stress shift; thetero and tvA are mutually exclusive)
%
% Outputs:
%   lh, likelihood, A0_prior, lambda_prior, var
%==========================================================================

thetero = isfield(options,'thetero') && ~isempty(options.thetero) && options.thetero;
tthr    = isfield(options,'tthr')    && ~isempty(options.tthr)    && options.tthr;

% tvA takes dispatch priority below, so enforce the documented mutual
% exclusions LOUDLY. Otherwise tvA silently wins and a run with thetero / tthr
% set would quietly become a plain tvA run: bvar_posterior_tvA does not carry
% the thetero shift, nor the tthr state-dependent-lag / Gamma-prior pass-through
% (those live in bvar_posterior / bvar_posterior_thetero -> rfvar3). build_tthr_irf
% is likewise not defined for a regime-specific A0.
if tvA && thetero
    error('eval_posterior: thetero and tvA are mutually exclusive (set one to 0).');
end
if tvA && tthr
    error(['eval_posterior: tthr (state-dependent lags + Gamma prior) is not ' ...
           'supported with tvA. Set options.tthr = 0 or options.tvA = 0.']);
end

if tvA
    if any(lcLmd(:))
        [lh, likelihood, A0_prior, lambda_prior, var] = ...
            bvar_posterior_tvA(x, y, lags, lcA0, lcLmd, Tsigbrk, delta0, options);
    else
        [lh, likelihood, A0_prior, lambda_prior, var] = ...
            bvar_posterior_tvA(x, y, lags, lcA0, lcLmd, Tsigbrk, [], options);
    end
elseif thetero
    if any(lcLmd(:))
        [lh, likelihood, A0_prior, lambda_prior, var] = ...
            bvar_posterior_thetero(x, y, lags, lcA0, lcLmd, Tsigbrk, delta0, options);
    else
        [lh, likelihood, A0_prior, lambda_prior, var] = ...
            bvar_posterior_thetero(x, y, lags, lcA0, lcLmd, Tsigbrk, [], options);
    end
elseif any(lcLmd(:))
    [lh, likelihood, A0_prior, lambda_prior, var] = ...
        bvar_posterior(x, y, lags, lcA0, lcLmd, Tsigbrk, delta0, options);
else
    [lh, likelihood, A0_prior, lambda_prior, var] = ...
        bvar_posterior(x, y, lags, lcA0, lcLmd, Tsigbrk, [], options);
end
end