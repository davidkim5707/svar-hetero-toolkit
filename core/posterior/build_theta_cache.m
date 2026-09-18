function [cache, stable] = build_theta_cache(x, var, lcA0, lcLmd, ny, lags, ...
    constant, fix_first_regime, accepted_shock, impact_signs, SR_parsed, ...
    upstream_stability_check, horizon, NarrativeRegimeRestrictions, stability_tol, ...
    lcA0_tv, SR_parsed_byregime, impact_signs_byregime)
%BUILD_THETA_CACHE  Precompute all theta-dependent quantities in ONE place.
%
%   Partial-tvA aware:
%   New trailing arg lcA0_tv (subset of lcA0 that is regime-specific):
%     []          -> legacy: full tvA / non-tvA inferred from numel(x)
%     given       -> partial tvA: A0 = [common once] + [tv per regime]
%
%   A0 is unpacked through the shared helper unpackA0_tvA so every consumer
%   reconstructs the regime-stacked A0 identically. nReg is derived from
%   lcA0_tv together with numel(x); under tvA A0, A0inv, By and Psi all become
%   regime-specific (extra trailing dimension). nReg == 1 reproduces the
%   legacy single-A0 behaviour byte-for-byte.
%
%   Under tvA the regime count comes from the A0 stacking (cache.nRegimes),
%   NOT from lambda. cache.lambda is expanded to ny x nRegimes and
%   cache.nlambda set to nRegimes so existing 1:nlambda regime loops iterate
%   over the A0 regimes (a no-heteroskedasticity tvA run gets lambda = ones).
%
%   Homoskedastic-shock handling (local_extractLambda):
%   In the free-regimes (non fix_first) branch the lambda matrix used to be
%   zeroed wholesale before filling the free cells. That left any shock whose
%   lcLmd row is all-false (a homoskedastic shock) at 0 across regimes, which
%   (i) is degenerate for the cached IRF scaling and (ii) is inconsistent with
%   the posterior, which inits lambda = ones. We now zero ONLY the rows that
%   carry free lambda; homoskedastic rows keep the ones-init and reconstruct
%   to a flat row of 1s (constant relative variance), matching the posterior.
%
%   Regime-specific cache fields under tvA:
%     cache.A0     : ny x ny x nRegimes
%     cache.A0inv  : ny x ny x nRegimes
%     cache.By     : ny x ny x lags x nRegimes
%     cache.Psi    : ny x ny x horizon x nRegimes   (or [] if horizon == 0)
%   Plus cache.tvA (logical) and cache.nRegimes.
%
%   Upstream stability (per regime under tvA): rejects the draw if ANY
%   regime companion has a root of modulus >= 1 - stability_tol.
%
%   Requires (standalone): unpackA0_tvA.m

    if nargin < 11, SR_parsed = [];                end
    if nargin < 12, upstream_stability_check = false; end
    if nargin < 13, horizon = 0;                   end
    if nargin < 14, NarrativeRegimeRestrictions = []; end
    if nargin < 15, stability_tol = 1e-6;          end
    if nargin < 16, lcA0_tv = [];                  end   % [] -> legacy detection
    if nargin < 17, SR_parsed_byregime = [];       end   % [] -> no regime gate
    if nargin < 18, impact_signs_byregime = [];    end   % [] -> no per-regime impact gate

    % --- Guard: eval_posterior returned an invalid draw (var = []) ---
    if isempty(var) || ~isstruct(var) || ~isfield(var, 'Bdraw') || isempty(var.Bdraw)
        cache  = struct();
        stable = false;
        return
    end

    % --- Detect regime count + partition (length + lcA0_tv) ---
    nA        = sum(lcA0(:));
    nLmd_free = sum(lcLmd(:));

    if isempty(lcA0_tv)
        % legacy: full tvA or non-tvA inferred from length
        nA_blocks = round((numel(x) - nLmd_free) / nA);
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
            nReg = round((numel(x) - nLmd_free - nA_cm) / nA_tv);
        else
            nReg = 1;             % common-only -> single A0
        end
        if nReg < 1, nReg = 1; end
    end
    tvA = nReg > 1;

    cache.tvA      = tvA;
    cache.nRegimes = nReg;

    % --- A0 and its inverse (shared unpacker; per regime under tvA) ---
    A0stack = unpackA0_tvA(x, ny, lcA0, lcA0_tv, nReg);   % ny x ny x nReg
    if tvA
        cache.A0    = A0stack;
        cache.A0inv = zeros(ny, ny, nReg);
        for r = 1:nReg
            cache.A0inv(:, :, r) = A0stack(:, :, r) \ eye(ny);
        end
    else
        cache.A0    = A0stack(:, :, 1);
        cache.A0inv = cache.A0 \ eye(ny);
    end
    A0offset = nA_cm + nA_tv * nReg;     % lambda block starts after all A0 params

    % --- By (reduced-form lag matrices): By_r = A0_r^{-1} * Aplus ---
    %     Aplus = var.Bdraw(1:lags*ny, :) is the COMMON structural lag block
    %     (rfvar3 estimates one B across regimes; A0 carries the regime).
    Btrim = var.Bdraw(1:lags*ny, :);
    if tvA
        cache.By = zeros(ny, ny, lags, nReg);
        for r = 1:nReg
            for iLag = 1:lags
                rows = (iLag - 1) * ny + (1:ny);
                cache.By(:, :, iLag, r) = cache.A0(:, :, r) \ Btrim(rows, :)';
            end
        end
    else
        cache.By = zeros(ny, ny, lags);
        for iLag = 1:lags
            rows = (iLag - 1) * ny + (1:ny);
            cache.By(:, :, iLag) = cache.A0 \ Btrim(rows, :)';
        end
    end

    % --- Optional upstream stability check (per regime under tvA) ---
    if upstream_stability_check
        if tvA
            for r = 1:nReg
                if local_max_root(cache.By(:, :, :, r), ny, lags) >= 1 - stability_tol
                    cache = struct(); stable = false; return
                end
            end
        else
            if local_max_root(cache.By, ny, lags) >= 1 - stability_tol
                cache = struct(); stable = false; return
            end
        end
    end

    % --- IRF propagator basis Psi^h (per regime under tvA) ---
    %     ir(:, h, :) = Psi(:, :, h) * smat  for h = 1, ..., horizon.
    if horizon > 0
        if tvA
            cache.Psi = zeros(ny, ny, horizon, nReg);
            for r = 1:nReg
                cache.Psi(:, :, :, r) = local_psi(cache.By(:, :, :, r), ny, lags, horizon);
            end
        else
            cache.Psi = local_psi(cache.By, ny, lags, horizon);
        end
    else
        cache.Psi = [];                              % fallback: use impulsdtrf
    end

    % --- Lambda (expanded to nReg columns under tvA) ---
    [lam, nlam] = local_extractLambda(x, A0offset, lcLmd, ny, fix_first_regime);
    if tvA && nlam < nReg
        % no-heteroskedasticity tvA run: regimes defined by A0 only
        lam  = ones(ny, nReg);
        nlam = nReg;
    end
    cache.lambda  = lam;
    cache.nlambda = nlam;

    % --- Constants ---
    cache.accepted_shock              = accepted_shock;
    cache.impact_signs                = impact_signs;
    cache.fix_first                   = fix_first_regime;
    cache.SR_parsed                   = SR_parsed;
    cache.SR_parsed_byregime          = SR_parsed_byregime;      % [] -> no gate
    cache.impact_signs_byregime       = impact_signs_byregime;   % [] -> no gate
    cache.NarrativeRegimeRestrictions = NarrativeRegimeRestrictions;

    stable = true;
end


% =====================================================================
% LOCAL HELPERS  (unchanged)
% =====================================================================

function Psi = local_psi(By, ny, lags, horizon)
    Psi = zeros(ny, ny, horizon);
    Psi(:, :, 1) = eye(ny);                          % tau = 0 (impact)
    for tau = 1:horizon - 1
        Psi_tau = zeros(ny);
        for k = 1:min(tau, lags)
            Psi_tau = Psi_tau + By(:, :, k) * Psi(:, :, tau - k + 1);
        end
        Psi(:, :, tau + 1) = Psi_tau;
    end
end


function mr = local_max_root(By, ny, lags)
    companion = zeros(ny*lags, ny*lags);
    for iLag = 1:lags
        companion(1:ny, (iLag-1)*ny+1 : iLag*ny) = By(:, :, iLag);
    end
    if lags > 1
        companion(ny+1:ny*lags, 1:ny*(lags-1)) = eye(ny*(lags-1));
    end
    mr = max(abs(eig(companion)));
end


function [lambda, nlambda] = local_extractLambda(x1, offA0, lcLmd, ny, fix_first)
% offA0 is passed explicitly: nA_cm + nA_tv*nReg (partial/full tvA) or nA
% (non-tvA), so the lambda block is read from the correct position in x.
    if all(lcLmd(:) == 0)
        lambda  = ones(ny, 1);
        nlambda = 1;
        return
    end

    nlambda = size(lcLmd, 2);
    lambda  = ones(ny, nlambda);

    if fix_first
        lcLmd(:, 1)  = false;
        lambda(:, 1) = 1;
    else
        % HOMOSKEDASTIC-SHOCK FIX: zero ONLY rows that carry free lambda.
        % A homoskedastic shock has an all-false lcLmd row; it must keep the
        % ones-init so it reconstructs to a flat row of 1s below (matching
        % constructLambda in the posterior). The old "lambda(:,:) = 0" left
        % such rows at 0 (degenerate + inconsistent with the posterior).
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