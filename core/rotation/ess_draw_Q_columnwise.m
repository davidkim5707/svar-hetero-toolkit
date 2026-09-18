function [X_new, Q_new, n_evals_total] = ess_draw_Q_columnwise(X_current, check_fn, max_shrinks)
%ESS_DRAW_Q_COLUMNWISE  Column-by-column elliptical slice sampling for Q.
%
%   [X_new, Q_new, n_evals_total] = ess_draw_Q_columnwise(X_current, check_fn, max_shrinks)
%
%   Instead of updating all n^2 elements of X simultaneously on a single
%   ellipse (as in ess_draw_Q), this function updates each column x_j of X
%   sequentially via ESS in R^n.  Each column update proposes along a 1D
%   ellipse in n-dimensional space, reconstructs Q = gaussian_to_Q(X),
%   and checks admissibility.
%
%   This reduces the per-step ESS search dimension from n^2 to n, making
%   it much easier for the 1D ellipse to intersect the admissible set in
%   higher-dimensional models.
%
%   The sequential column updates form a valid Gibbs sampler on X:
%     p(x_j | x_{-j}, SR > 0) ∝ [SR(gamma(X)) > 0] * N(x_j; 0, I_n)
%   Since each conditional has a Gaussian kernel, ESS applies.
%
%   INPUTS
%     X_current   : n-by-n real matrix (current Gaussian state)
%     check_fn    : function handle  check_fn(Q) -> true/false
%     max_shrinks : max bracket shrinks per column (default 200)
%
%   OUTPUTS
%     X_new        : n-by-n Gaussian matrix (updated state)
%     Q_new        : orthogonal matrix satisfying check_fn
%     n_evals_total: total check_fn evaluations across all columns

    if nargin < 3 || isempty(max_shrinks)
        max_shrinks = 200;
    end

    n = size(X_current, 1);
    X_new = X_current;
    n_evals_total = 0;

    for j = 1:n
        x_j = X_new(:, j);

        % Draw auxiliary noise for this column
        nu_j = randn(n, 1);

        % Initial angle and bracket
        theta     = 2 * pi * rand();
        theta_min = theta - 2 * pi;
        theta_max = theta;

        found = false;

        for iter = 1:max_shrinks
            n_evals_total = n_evals_total + 1;

            % Elliptical proposal for column j only (R^n, not R^{n^2})
            x_star = nu_j * sin(theta) + x_j * cos(theta);

            % Replace column j and reconstruct full Q
            X_candidate = X_new;
            X_candidate(:, j) = x_star;
            Q_candidate = gaussian_to_Q(X_candidate);

            % Check all sign restrictions under the new Q
            if check_fn(Q_candidate)
                X_new(:, j) = x_star;
                found = true;
                break
            end

            % Shrink bracket
            if theta < 0
                theta_min = theta;
            else
                theta_max = theta;
            end
            theta = theta_min + (theta_max - theta_min) * rand();
        end

        % If not found, keep current column unchanged.
        % This is valid: the old X_new still satisfies SR.
        % The chain simply doesn't move in this column's direction.
    end

    Q_new = gaussian_to_Q(X_new);
end