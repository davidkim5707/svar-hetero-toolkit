function [X, Q, success, n_outer] = draw_Q_zero_columnwise( ...
    ny, check_full_fn, A0inv, Lambda_all, impact_signs, ...
    sign_regime_dependent, zc, max_outer, max_inner)
%DRAW_Q_ZERO_COLUMNWISE  Column-by-column rotation draw with ZERO restrictions
%   (ARRW 2018, Algorithm 2 construction). Generalizes
%   find_admissible_Q_columnwise: builds the columns of P (= Q') in zc.col_order
%   (most-restricted shock first), drawing each column in the null space of
%       M_j = [ columns built before j ; Z{j} ],
%   so the linear zero restrictions Z{j} hold EXACTLY by construction (no
%   rejection on the zeros, only on the impact signs). When zc has no zeros it
%   delegates to find_admissible_Q_columnwise, reproducing the existing
%   behaviour byte-for-byte.
%
%   [X, Q, success, n_outer] = draw_Q_zero_columnwise( ...
%       ny, check_full_fn, A0inv, Lambda_all, impact_signs, ...
%       sign_regime_dependent, zc, max_outer, max_inner)
%
% ----------------------------------------------------------------------------
%   WHY THIS, NOT REJECTION: the set of rotations satisfying the zeros is of
%   measure zero in O(n), so a random/Haar draw hits it with probability 0.
%   Each zero is a linear restriction c' p_{j0} = 0 on a single column of P
%   (see build_zero_constraint); the null-space construction draws p_j directly
%   in {orthogonal to previous columns} AND {orthogonal to the c's}, which is
%   exactly ARRW's Algorithm 2.
%
%   BUILD ORDER & LABELS: columns are *filled* in zc.col_order, but each goes
%   into its LABELED position, P(:, col_order(k)) = (k-th built column), so
%   column j of the returned P/Q is always shock j. Building most-restricted-
%   first (a) keeps every null space non-empty (zc.feasible) and (b) makes the
%   single-restricted-column case reproduce conditional-Haar exactly, so NO
%   volume-element weight is needed there (zc.single_column == true). When zeros
%   hit >= 2 different shocks the construction measure differs from Haar by the
%   volume element; reweight downstream (separate module). This function does
%   NOT compute that weight.
%
%   CONVENTION (matches check_sr_pass_cached / find_admissible_Q_columnwise):
%     sign_regime_dependent == 0:  smat(:,j) = sqrt(lambda_{m,j}) * A0inv_m * pj
%     sign_regime_dependent == 1:  smat(:,j) = A0inv_m * (sqrt_lambda_m .* pj)
%   with pj = P(:,j) = Q(j,:)'. A0inv may be ny x ny or ny x ny x nRegimes (tvA).
%   The impact-sign SCREEN below is identical to find_admissible_Q_columnwise:
%   it requires the impact signs to hold under EVERY regime.
%
%   IMPORTANT: A0inv and Lambda_all MUST be the same matrices as in the cache
%   used to build zc (build_zero_constraint), so that the constraint vectors in
%   zc.Z and the sign screen here refer to the same draw.
%
%   OUTPUT X (the ESS latent): for the ZERO case X = P(:, col_order), the
%   gaussian_to_Q_zero latent, so gaussian_to_Q_zero(X, zc) == Q and the chain
%   continues with ess_draw_Q_zero_columnwise. For the no-zero case X is the
%   usual gaussian_to_Q latent returned by find_admissible_Q_columnwise. The
%   sampler must map the latent with gaussian_to_Q_zero(., zc) (which reduces to
%   gaussian_to_Q when zc has no zeros), NOT plain gaussian_to_Q.
%
%   Requires: find_admissible_Q_columnwise.m (no-zero delegation).
% ----------------------------------------------------------------------------

    if nargin < 8 || isempty(max_outer), max_outer = 200;   end
    if nargin < 9 || isempty(max_inner), max_inner = 10000; end

    % ---- No-zero fast path: delegate to find_admissible_Q_columnwise -----
    if nargin < 7 || isempty(zc) || ~isfield(zc, 'has_zero') || ~zc.has_zero
        [X, Q, success, n_outer] = find_admissible_Q_columnwise( ...
            ny, check_full_fn, A0inv, Lambda_all, impact_signs, ...
            sign_regime_dependent, max_outer, max_inner);
        return
    end

    % ---- Infeasible zero pattern: no rotation can satisfy it --------------
    if isfield(zc, 'feasible') && ~zc.feasible
        error('draw_Q_zero_columnwise:infeasible', ...
              ['Zero pattern over-restricts (zc.feasible == false): no ' ...
               'orthogonal matrix satisfies all zeros. Fix options.ZeroRestrictions.']);
    end

    nRegimes = size(Lambda_all, 2);
    A0inv_tv = (ndims(A0inv) == 3);          % 3-D => regime-specific A0inv (tvA)
    success  = false;
    X = zeros(ny, ny);
    Q = eye(ny);

    % Precompute sqrt(Lambda) per regime
    sqrtLam = zeros(ny, nRegimes);
    for m = 1:nRegimes
        sqrtLam(:, m) = sqrt(max(Lambda_all(:, m), 0));
    end

    has_restriction = any(impact_signs ~= 0, 1);   % 1 x ny
    col_order = zc.col_order;

    for n_outer = 1:max_outer

        P = zeros(ny, ny);          % columns filled in labeled positions
        all_cols_ok = true;

        for k = 1:ny
            j = col_order(k);                       % shock / labeled column

            % --- admissible subspace: null( [built-before-j ; Z{j}] ) ------
            if k > 1
                P_built = P(:, col_order(1:k-1));    % ny x (k-1)
            else
                P_built = zeros(ny, 0);
            end
            Zj = zc.Z{j};                            % z_j x ny (rows c'), [] if none
            Mj = [P_built.'; Zj];                    % (k-1 + z_j) x ny

            if size(Mj, 1) == 0
                Kj = eye(ny);
            else
                Kj = null(Mj);                       % ny x (ny - rank(Mj))
            end
            d = size(Kj, 2);
            if d < 1
                % over-constrained at this draw (rank dropped) -> rebuild
                all_cols_ok = false; break
            end

            % --- inner loop: draw in the subspace, screen impact signs -----
            found_col = false;
            for inner = 1:max_inner
                w  = randn(d, 1);
                pj = Kj * w;
                nrm = norm(pj);
                if nrm < 1e-12, continue; end
                pj = pj / nrm;                       % uniform on the null-space sphere

                if ~has_restriction(j)
                    found_col = true; break
                end

                col_ok = true;
                for m = 1:nRegimes
                    if A0inv_tv, A0inv_m = A0inv(:, :, m); else, A0inv_m = A0inv; end
                    if sign_regime_dependent == 0
                        sj = sqrtLam(j, m) * (A0inv_m * pj);
                    else
                        sj = A0inv_m * (sqrtLam(:, m) .* pj);
                    end
                    for i = 1:ny
                        if impact_signs(i, j) ~= 0 && impact_signs(i, j) * sj(i) <= 0
                            col_ok = false; break
                        end
                    end
                    if ~col_ok, break; end
                end

                if col_ok, found_col = true; break; end
            end  % inner

            if ~found_col
                all_cols_ok = false; break
            end
            P(:, j) = pj;                            % place into labeled position
        end  % column

        if ~all_cols_ok
            continue
        end

        % ---- Full validation (multi-horizon / narrative / stability) ------
        % Zeros already hold EXACTLY by construction; check_full_fn carries the
        % sign + narrative + stability checks. (A zero-residual safety check
        % belongs in the storage-stage SignRestrictionCheck, not here.)
        Q_candidate = P.';
        if check_full_fn(Q_candidate)
            Q = Q_candidate;
            X = P(:, col_order);                     % ESS latent: gaussian_to_Q_zero(X, zc) == Q
            success = true;
            return
        end
    end  % outer
end