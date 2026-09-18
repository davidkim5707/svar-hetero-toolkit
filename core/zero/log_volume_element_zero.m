function logw = log_volume_element_zero(Qdraw, zc, ny)
%LOG_VOLUME_ELEMENT_ZERO  ARRW (2018) construction volume element (storevegBS).
%
%   logw = LOG_VOLUME_ELEMENT_ZERO(Qdraw, zc, ny) returns the log importance
%   weight that reweights null-space-construction draws of the rotation to the
%   uniform (conditional-Haar) measure on the zero-restricted manifold.
%
%   Qdraw is the rotation in the sampler's convention Q = P', i.e. the shock
%   columns p_j are the ROWS of Qdraw (p_j = Qdraw(j,:)' = P(:,j)); this matches
%   the output of gaussian_to_Q_zero / draw_Q_zero_columnwise.
%
%   This is ARRW's storevegBS = LogVolumeElement(f, q_tilde, f_restrictions):
%   the volume element of the spheres-to-Q construction map.  For a sampler that
%   draws (A0, Lambda, B) DIRECTLY in structural coordinates (Sstar prior) -- as
%   opposed to ARRW's NIW reduced-form draw + transform -- this is the ONLY
%   weight term.  ARRW's f_h volume elements (storevefh, storevefhZ) are constant
%   per (A0,Lambda,B) draw and cancel, so they do not enter here.
%
%   Apply as: w = exp(logw); reweight the stored draw by w (then normalize the
%   weights across draws).
%
%   Returns 0 when there are no zero restrictions, or when all zeros lie on a
%   single shock/column (zc.single_column).  The construction then coincides
%   with conditional-Haar and the volume element is constant, so no reweighting
%   is needed -- proven, and confirmed numerically (constant to ~1e-10).
%
%   Implementation calls ARRW's own helpers on the path (verified to reproduce
%   ARRW's exact volume element per draw to ~1e-10): SpheresToQ, SpheresRestriction,
%   LogVolumeElement, perp, NumericalDerivative, LogAbsDet.  The col_order is
%   handled by reordering the zero constraints into build order (Zb) and calling
%   ARRW's natural-order SpheresToQ on Zb -- this IS the col_order construction,
%   up to a fixed permutation of output columns that leaves the volume element
%   unchanged.

    if ~zc.has_zero || zc.single_column
        logw = 0;
        return
    end

    co = zc.col_order(:).';

    % zero constraints reordered into BUILD order: Zb{k} = rows for the k-th built shock
    Zb = cell(1, ny);
    for k = 1:ny
        j = co(k);
        if j <= numel(zc.Z) && ~isempty(zc.Z{j})
            Zb{k} = zc.Z{j};
        else
            Zb{k} = zeros(0, ny);
        end
    end

    % generic (full-rank) sphere block sizes: s(k) = ny - (k-1) - z_{co(k)}
    s = zeros(1, ny);
    for k = 1:ny
        s(k) = ny - (k-1) - size(Zb{k}, 1);
    end
    if any(s < 1)
        error('log_volume_element_zero:infeasible', ...
              'Construction step with non-positive dimension; check col_order/feasibility.');
    end

    % recover sphere coordinates q_tilde from the draw (build order; same perp
    % convention as SpheresToQ, so it round-trips Qdraw regardless of how Qdraw
    % was drawn -- the volume element is a function of the rotation alone)
    P = Qdraw.';                              % P(:,j) = p_j
    q_tilde = zeros(sum(s), 1);
    idx = 0;
    Pb = zeros(ny, 0);
    for k = 1:ny
        Np = perp([Pb, Zb{k}.']);
        q_tilde(idx+1:idx+s(k)) = Np(:, 1:s(k)).' * P(:, co(k));
        idx = idx + s(k);
        Pb = [Pb, P(:, co(k))];
    end

    % ARRW construction volume element (exact ARRW helpers)
    logw = LogVolumeElement(@(x) SpheresToQ(x, Zb, ny), q_tilde, ...
                            @(x) SpheresRestriction(x, Zb, ny));
end