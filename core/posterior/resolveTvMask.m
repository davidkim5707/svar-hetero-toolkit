function lcA0_tv = resolveTvMask(lcA0, ny, options)
%RESOLVETVMASK  Which FREE A0 elements are regime-specific (time-varying).
%   The complement within lcA0 is COMMON (shared across regimes). This is the
%   single source of truth for the partial-tvA partition and is called by
%   bvar_posterior_tvA, the sampler (SECTION 4/12/14), build_theta_cache and
%   SignRestrictionCheck so every site agrees on the split.
%
%   Priority:
%     options.lcA0_tv   : explicit ny x ny logical mask
%     options.tvA_rows  : vector of row indices -> all free elements in them
%     (neither)         : full tvA (every free element regime-specific) = legacy
%
%   Returns lcA0_tv as a subset of lcA0 (intersection guarantees only free
%   elements can be time-varying).
%
%   Diagonal normalization:
%   After resolving, the A0 DIAGONAL is forced to be COMMON across regimes
%   (removed from lcA0_tv). A regime-specific A0 diagonal is observationally
%   equivalent to the regime variance lambda_{i,s}: scaling row i of A0_s by c
%   and lambda_{i,s} by c^2 leaves Sigma_s = A0_s^{-1} Lambda_s A0_s^{-T}
%   unchanged, so a fully-tv row leaves (nSig-1) unidentified scale directions.
%   Holding A0_ii common anchors each tv row's scale (the Bacchiocchi-Fanelli
%   2015 Example-2 q33=0 analog) and loses NO identifiable variation, since the
%   regime-specific own-scale is already carried by lambda_{i,s}.
%
%   Controlled by options.tvA_fix_diagonal (default TRUE). Set it to FALSE to
%   reproduce the exact legacy full-tvA mask (lcA0_tv == lcA0).
%--------------------------------------------------------------------------
    if isfield(options, 'lcA0_tv') && ~isempty(options.lcA0_tv)
        m = logical(options.lcA0_tv);
        if ~isequal(size(m), [ny ny])
            error('resolveTvMask: options.lcA0_tv must be %dx%d logical.', ny, ny);
        end
        lcA0_tv = m & lcA0;
    elseif isfield(options, 'tvA_rows') && ~isempty(options.tvA_rows)
        m = false(ny, ny);
        m(options.tvA_rows(:), :) = true;
        lcA0_tv = m & lcA0;
    else
        lcA0_tv = lcA0;                 % legacy: full tvA
    end

    % --- Hold tv-row DIAGONALS common across regimes (scale normalization) ---
    fix_diag = true;
    if isfield(options, 'tvA_fix_diagonal') && ~isempty(options.tvA_fix_diagonal)
        fix_diag = logical(options.tvA_fix_diagonal);
    end
    if fix_diag
        lcA0_tv(1:ny+1:end) = false;   % strip the diagonal from the tv set
    end
end