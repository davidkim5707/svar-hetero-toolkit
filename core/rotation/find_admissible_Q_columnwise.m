function [X, Q, success, n_outer] = find_admissible_Q_columnwise( ...
    ny, check_full_fn, A0inv, Lambda_all, impact_signs, ...
    sign_regime_dependent, max_outer, max_inner)
%FIND_ADMISSIBLE_Q_COLUMNWISE  Column-by-column Q initialization (Algorithm 4).
%
%   tvA-aware:
%   A0inv may be ny-by-ny (single A0) or ny-by-ny-by-nRegimes (time-varying
%   A0). When 3-D, the per-column impact screen uses the regime-specific
%   A0inv(:,:,m) for each regime m, so the screen requires the impact signs
%   to hold under EVERY regime's contemporaneous structure. Everything else
%   is unchanged; the signature is identical, so the caller can pass
%   cache.A0inv directly whether or not tvA is on.
%
%   IMPORTANT CONVENTION (matches check_sr_pass_cached / SignRestrictionCheck):
%     sign_regime_dependent == 0:  smat(:,j) = sqrt(lambda_{m,j}) * A0inv_m * pj
%     sign_regime_dependent == 1:  smat(:,j) = A0inv_m * (sqrt_lambda_m .* pj)
%   with pj = Q(j,:)' = P(:,j) and A0inv_m the regime-m inverse.
%
%   INPUTS / OUTPUTS as before; A0inv now optionally 3-D.

if nargin < 7 || isempty(max_outer), max_outer = 200;   end
if nargin < 8 || isempty(max_inner), max_inner = 10000; end

nRegimes = size(Lambda_all, 2);
A0inv_tv = (ndims(A0inv) == 3);       % 3-D => regime-specific A0inv (tvA)
success  = false;
X = zeros(ny, ny);
Q = eye(ny);

% Precompute sqrt(Lambda) for each regime
sqrtLam = zeros(ny, nRegimes);
for m = 1:nRegimes
    sqrtLam(:, m) = sqrt(max(Lambda_all(:, m), 0));
end

% Identify which shocks have impact restrictions
has_restriction = any(impact_signs ~= 0, 1);  % 1-by-ny logical

for n_outer = 1:max_outer

    % ================================================================
    % Build P = Q' column by column.
    %   P(:,j) = j-th column of P = j-th row of Q = Q(j,:)'
    % ================================================================
    P_cols = zeros(ny, 0);
    all_cols_ok = true;

    for j = 1:ny

        found_col = false;

        for inner = 1:max_inner
            % ---- Step 1: draw from N(0, I_n) ----
            xj = randn(ny, 1);

            % ---- Step 2: project onto orthogonal complement ----
            if j > 1
                xj = xj - P_cols * (P_cols' * xj);
            end

            % ---- Step 3: normalize ----
            nrm = norm(xj);
            if nrm < 1e-12
                continue
            end
            pj = xj / nrm;

            % ---- Step 4: per-column sign check (impact only) ----
            if ~has_restriction(j)
                found_col = true;
                P_cols = [P_cols, pj]; %#ok<AGROW>
                break
            end

            % Check impact signs across every regime (regime-specific A0inv)
            col_ok = true;
            for m = 1:nRegimes
                if A0inv_tv
                    A0inv_m = A0inv(:, :, m);
                else
                    A0inv_m = A0inv;
                end

                if sign_regime_dependent == 0
                    sj = sqrtLam(j, m) * (A0inv_m * pj);
                else
                    sj = A0inv_m * (sqrtLam(:, m) .* pj);
                end

                for i = 1:ny
                    if impact_signs(i, j) ~= 0
                        if impact_signs(i, j) * sj(i) <= 0
                            col_ok = false;
                            break
                        end
                    end
                end
                if ~col_ok, break; end
            end

            if col_ok
                found_col = true;
                P_cols = [P_cols, pj]; %#ok<AGROW>
                break
            end
        end  % inner loop

        if ~found_col
            all_cols_ok = false;
            break
        end
    end  % column loop

    if ~all_cols_ok
        continue
    end

    % ---- Assemble Q from P  (P = Q', so Q = P') ----
    Q_candidate = P_cols';

    % ---- Full validation (multi-horizon / narrative / elasticity / stability) ----
    if check_full_fn(Q_candidate)
        Q = Q_candidate;

        [~, R_seed] = qr(randn(ny, ny));
        X = Q * diag(abs(diag(R_seed)));

        success = true;
        return
    end
end

end