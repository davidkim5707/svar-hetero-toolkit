function Q = gaussian_to_Q_zero(G, zc)
%GAUSSIAN_TO_Q_ZERO  Map a Gaussian matrix G to an orthonormal rotation Q (= P')
%   by ZERO-AWARE Gram-Schmidt in zc.col_order, so the linear zero restrictions
%   in zc hold EXACTLY. This is the manifold-confined replacement for
%   gaussian_to_Q used by the elliptical slice sampler in the zero path.
%
%   Q = gaussian_to_Q_zero(G, zc)
%
%   WHY a replacement is needed: gaussian_to_Q does a plain QR in the fixed
%   column order 1..n, which orthonormalizes the COLUMNS of Q. The zero
%   restrictions, however, are linear constraints on the ROWS of Q (the columns
%   p_j of P = Q'; see build_zero_constraint). Perturbing a single column of the
%   QR input therefore does NOT keep a row of Q on its constraint. This map
%   instead builds the columns p_j of P directly, in zc.col_order, projecting
%   each onto null([previously built columns ; Z{j}]) -- exactly the ARRW (2018)
%   Algorithm-2 construction -- so every p_j stays in its admissible subspace.
%
%   BASIS INVARIANCE: each step applies the ORTHOGONAL PROJECTOR onto
%   null(Mj), computed as K*(K'*g). The projector is unique for a given
%   subspace, so the output is IDENTICAL for every orthonormal basis K that
%   null() may return, and the map is well defined and a.e. continuous
%   wherever rank(Mj) is locally constant. No fixed-basis (e.g. fixed-W_j QR)
%   construction of the null space is required.
%
%   COLUMN FEED: column k of G feeds construction step k, i.e. P(:, col_order(k)).
%   Because step k uses only G(:,1:k) (via the previous columns), perturbing
%   G(:,k) changes P(:, col_order(k:end)) and leaves earlier (already more
%   restricted) columns -- hence all their zeros -- untouched. This is what makes
%   the column-wise elliptical slice sampler stay on the zero manifold.
%
%   SEED: gaussian_to_Q_zero(P(:, col_order), zc) == P' = Q, so the
%   columnwise initializer returns X = P(:, col_order) as the ESS latent.
%
%   No zeros (zc empty / zc.has_zero == false): delegates to gaussian_to_Q, so
%   the no-zero path is unchanged.
%
%   Requires: gaussian_to_Q.m (no-zero delegation).

    if nargin < 2 || isempty(zc) || ~isfield(zc, 'has_zero') || ~zc.has_zero
        Q = gaussian_to_Q(G);
        return
    end

    ny        = size(G, 1);
    col_order = zc.col_order;
    P         = zeros(ny, ny);

    for k = 1:ny
        j = col_order(k);
        g = G(:, k);

        if k > 1
            P_built = P(:, col_order(1:k-1));    % ny x (k-1), already placed
        else
            P_built = zeros(ny, 0);
        end

        Zj = zc.Z{j};                            % z_j x ny (rows c'), [] if none
        Mj = [P_built.'; Zj];                    % (k-1 + z_j) x ny

        if size(Mj, 1) == 0
            K = eye(ny);
        else
            K = null(Mj);                        % ny x dim(null); basis choice
        end                                      %   is irrelevant (see header)

        % Rank check. dim N_j must equal ny - (k-1) - z_j, i.e. the
        % stacked M_j must have full row rank (the fixed-theta analogue of
        % ARRW's regularity condition on F). A violation is an exact,
        % measure-zero alignment between the constraint rows and the built
        % columns; stop rather than build a rotation off the manifold.
        % Zero extra cost: null() is already computed.
        if size(K, 2) ~= ny - size(Mj, 1)
            error('gaussian_to_Q_zero:rankDeficientMj', ...
                ['Construction step %d (shock %d): dim null(M_j) = %d but ' ...
                 'full row rank implies %d. M_j is rank deficient ' ...
                 '(measure-zero alignment); reject this theta draw.'], ...
                k, j, size(K, 2), ny - size(Mj, 1));
        end

        proj = K * (K.' * g);                    % orthogonal projector applied
        nrm  = norm(proj);                       %   to g (basis-invariant)

        if nrm < 1e-12
            % Numerical guard for a MEASURE-ZERO event: g is (numerically)
            % orthogonal to the admissible subspace, so the projection cannot
            % be normalized. Substitute a deterministic unit vector inside the
            % subspace so Q stays orthonormal and the zeros stay satisfied.
            % Under Gaussian / elliptical inputs this branch fires with
            % probability zero and leaves the sampler's target untouched.
            proj = K(:, 1);
            nrm  = 1;
        end

        P(:, j) = proj / nrm;
    end

    Q = P.';
end