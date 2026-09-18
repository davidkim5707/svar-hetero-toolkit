function out = check_bf_identification(x, lcA0, lcA0_tv, lcLmd, ny, nRegimes, varargin)
% CHECK_BF_IDENTIFICATION
%   Bacchiocchi & Fanelli (2015, OBES, doi:10.1111/obes.12092) local-
%   identification rank check for the heteroskedastic + (partial-)tvA
%   structural block of the sign+hetero SVAR.
%
%   Tests whether the regime covariances {Sigma_s}, s = 1..nRegimes, locally
%   point-identify the free structural parameters
%       theta = x = [ a_cm (nA_cm) ; a_tv (nA_tv x nRegimes, regime-major) ;
%                     lmd_free (ny*(nRegimes-1)) ]
%   i.e. exactly the seedx layout produced in SECTION 4 of
%   svar_run_sign_elliptical. This is B&F Proposition 1 / Rothenberg (1971):
%   full column rank of the n(n+1)/2 * nRegimes  x  a  Jacobian
%       J = d vech(Sigma_s) / d theta'   stacked over regimes,
%   evaluated at a regular point (the posterior mode, or random draws).
%
%   B&F order condition (necessary): a <= nRegimes * ny(ny+1)/2.
%   B&F rank  condition (necessary & sufficient, LOCAL): rank(J) = a.
%
%   SCOPE. This is the POINT-ID condition for {A0_s, lambda_s} from second
%   moments only. It does NOT cover the sign/narrative + Q-rotation layer
%   (that is SET identification: RWZ 2010 / Giacomini-Kitagawa-Read). Read it
%   as: if rank(J)=a, the volatility regimes point-identify A0_s & lambda_s and
%   the sign restrictions merely label/orient the shocks; if rank(J)<a, the
%   deficient directions are pinned only by sign/narrative + prior (+shrinkage).
%
%   LAMBDA CONVENTION. Sigma_s is built as A0_s^{-1} * diag(lmd(:,s)) * A0_s^{-T}
%   with the average-1 row-sum constraint from constructLambda's free-regimes
%   branch (lmd(:,end) = nSig - sum(lmd(:,1:end-1),2)). If your run uses
%   fix_first_regime, or if lmd is a PRECISION rather than a VARIANCE in your
%   rfvar3 path, edit build_lmd / build_g accordingly. VALIDATE the convention
%   once (see the README block printed at the end): Sigma_s(theta_hat) should
%   match the sample covariance of the reduced-form residuals inside regime s.
%
%   DEPENDENCY. unpackA0_tvA.m (your standalone partial-tvA unpacker), used
%   exactly as in SECTION 12/14 so A0_s is reconstructed identically.
%
%   USAGE
%       % at the posterior mode (last stored draw, or your xh):
%       out = check_bf_identification(output.draw_x0(:,end), ...
%                 output.lcA0, output.lcA0_tv, output.lcLmd, ...
%                 options.ny, numel(output.lrange));
%
%       % robustness: loop over several draws / random regular points
%       for j = 1:size(output.draw_x0,2)
%           r(j) = check_bf_identification(output.draw_x0(:,j), ...).rank;
%       end
%
%   OUTPUT (struct)
%       .dim_theta .n_moments .order_ok .order_slack
%       .rank .full_rank .sing_values .smallest_sv .sv_gap_ratio
%       .null_dirs   (only if rank-deficient; columns are unidentified
%                     parameter combinations in seedx coordinates)

    p = inputParser;
    addParameter(p, 'tol', 1e-7);              % relative svd cutoff for rank
    addParameter(p, 'h',   1e-6);              % central finite-diff step
    addParameter(p, 'nSig', nRegimes);         % row-sum target in constructLambda
    addParameter(p, 'verbose', true);
    parse(p, varargin{:});
    tol = p.Results.tol; h = p.Results.h; nSig = p.Results.nSig; vb = p.Results.verbose;

    x = x(:);
    a = numel(x);                              % # free structural params = dim(theta)

    lcA0    = logical(lcA0);
    lcA0_tv = logical(lcA0_tv);
    lcLmd   = logical(lcLmd);
    nA_cm   = sum(lcA0(:) & ~lcA0_tv(:));
    nA_tv   = sum(lcA0_tv(:));
    nA0     = nA_cm + nA_tv * nRegimes;        % length of the A0 block in seedx
    nLmd    = sum(lcLmd(:));

    if a ~= nA0 + nLmd
        warning(['length(x)=%d but nA0+nLmd=%d. Check that x is a tvA seedx ' ...
                 'and that lcA0_tv/lcLmd match the run.'], a, nA0 + nLmd);
    end

    g0 = build_g(x);
    m  = numel(g0);                            % # moments = nRegimes * ny(ny+1)/2

    % ---- central finite-difference Jacobian (m x a) ----
    J = zeros(m, a);
    for k = 1:a
        xp = x; xp(k) = xp(k) + h;
        xm = x; xm(k) = xm(k) - h;
        J(:, k) = (build_g(xp) - build_g(xm)) / (2 * h);
    end

    s   = svd(J);
    rnk = sum(s > tol * s(1));

    out.dim_theta    = a;
    out.n_moments    = m;
    out.order_ok     = (a <= m);
    out.order_slack  = m - a;
    out.rank         = rnk;
    out.full_rank    = (rnk == a);
    out.sing_values  = s;
    out.smallest_sv  = s(end);
    out.sv_gap_ratio = s(min(a, numel(s))) / s(1);

    if ~out.full_rank
        [~, ~, V] = svd(J, 'econ');
        out.null_dirs = V(:, rnk+1:end);       % unidentified directions (seedx coords)
    else
        out.null_dirs = [];
    end

    if vb
        fprintf('\n--- B&F (2015) local-ID rank check ---\n');
        fprintf('  free params  a            : %d\n', a);
        fprintf('  moments  S*ny(ny+1)/2     : %d\n', m);
        fprintf('  order condition a<=m      : %s (slack %d)\n', tern(out.order_ok), out.order_slack);
        fprintf('  rank(J)                   : %d / %d  ->  %s\n', rnk, a, ...
            ternstr(out.full_rank, 'LOCALLY IDENTIFIED', 'RANK-DEFICIENT'));
        fprintf('  a-th singular value       : %.3e\n', s(min(a, numel(s))));
        fprintf('  ratio (a-th / largest)    : %.3e   (tiny => weak/under-identified)\n', out.sv_gap_ratio);
        if ~out.full_rank
            fprintf('  -> %d unidentified direction(s) in out.null_dirs (seedx coords:\n', a - rnk);
            fprintf('     [1:%d]=common A0, [%d:%d]=tv A0 (regime-major), [%d:%d]=free lmd)\n', ...
                nA_cm, nA_cm+1, nA0, nA0+1, a);
        end
        fprintf('  VALIDATE convention once: compare Sigma_s(theta_hat) below to\n');
        fprintf('  cov() of reduced-form residuals within each regime window.\n');
        fprintf('--------------------------------------\n');
    end

    % ============================ nested builders ============================
    function g = build_g(xx)
        A0s = unpackA0_tvA(xx, ny, lcA0, lcA0_tv, nRegimes);   % ny x ny x nRegimes
        lmd = build_lmd(xx);                                   % ny x nRegimes (variance)
        nv  = ny * (ny + 1) / 2;
        g   = zeros(nRegimes * nv, 1);
        for ss = 1:nRegimes
            A0inv = A0s(:, :, ss) \ eye(ny);
            Sg    = A0inv * diag(lmd(:, ss)) * A0inv';
            g((ss-1)*nv + (1:nv)) = Sg(tril(true(ny)));        % vech
        end
    end

    function lmd = build_lmd(xx)
        lmd = ones(ny, nRegimes);
        lmd(lcLmd) = xx(nA0 + (1:nLmd));
        % free-regimes constraint (matches applyFreeRegimes):
        rowSums = sum(lmd(:, 1:(nRegimes-1)), 2);
        lmd(:, nRegimes) = nSig - rowSums;
    end
end

function o = tern(c),            if c, o='OK'; else, o='FAIL'; end, end
function o = ternstr(c, a, b),   if c, o=a;    else, o=b;     end, end