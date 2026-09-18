function [X_new, Q_new, n_evals_total] = ess_draw_Q_zero_columnwise( ...
    X_current, check_fn, zc, max_shrinks)
%ESS_DRAW_Q_ZERO_COLUMNWISE  Column-wise elliptical slice sampler for the
%   rotation, CONFINED to the zero manifold. Identical in structure to
%   ess_draw_Q_columnwise -- a per-column elliptical slice on the Gaussian
%   latent X -- except the latent is mapped to Q through gaussian_to_Q_zero
%   instead of gaussian_to_Q, so the zero restrictions in zc hold EXACTLY at
%   every proposal and in the returned draw. Only the sign / narrative
%   restrictions are screened by check_fn; the zeros are never rejected.
%
%   [X_new, Q_new, n_evals_total] = ess_draw_Q_zero_columnwise( ...
%       X_current, check_fn, zc, max_shrinks)
%
%   LATENT CONVENTION: X is the gaussian_to_Q_zero latent, i.e. column k of X
%   feeds construction step k -> P(:, col_order(k)), with gaussian_to_Q_zero(X,
%   zc) == Q. The columnwise initializer draw_Q_zero_columnwise returns exactly
%   this X for the zero case. (For the no-zero case, X is the usual
%   gaussian_to_Q latent and this function delegates to ess_draw_Q_columnwise.)
%
%   WHY ZEROS SURVIVE THE SWEEP: perturbing column k of X changes only
%   P(:, col_order(k:end)); the earlier (more restricted) columns and all their
%   zeros are untouched, and the perturbed and later columns are re-projected
%   into their own null spaces by gaussian_to_Q_zero. So each accepted column
%   move keeps Q on the manifold. (Verified: orthonormality and zero residuals
%   ~ 1e-16 under arbitrary single-column elliptical perturbation.)
%
%   This sampler leaves invariant [null-space construction measure] x
%   1[sign/narrative], the SAME confined measure draw_Q_zero_columnwise draws
%   from. When zeros span >= 2 shocks that measure differs from conditional-Haar
%   by the volume element; reweight at storage (separate module). Single
%   restricted shock => the measure already IS conditional-Haar, no weight.
%
%   Requires: gaussian_to_Q_zero.m; ess_draw_Q_columnwise.m (no-zero delegation).

    if nargin < 4 || isempty(max_shrinks), max_shrinks = 20; end

    % ---- No-zero fast path: delegate to the plain column-wise ESS ----------
    if nargin < 3 || isempty(zc) || ~isfield(zc, 'has_zero') || ~zc.has_zero
        [X_new, Q_new, n_evals_total] = ess_draw_Q_columnwise( ...
            X_current, check_fn, max_shrinks);
        return
    end

    n             = size(X_current, 1);
    X_new         = X_current;
    n_evals_total = 0;

    % Sweep the latent columns in feed order (column k -> construction step k).
    for k = 1:n
        x_k     = X_new(:, k);
        nu_k    = randn(n, 1);                   % auxiliary Gaussian (full R^n)
        theta     = 2*pi*rand();
        theta_min = theta - 2*pi;
        theta_max = theta;                       % bracket always contains theta = 0

        for s = 1:max_shrinks
            n_evals_total = n_evals_total + 1;

            x_star          = x_k*cos(theta) + nu_k*sin(theta);   % elliptical proposal
            X_cand          = X_new;
            X_cand(:, k)    = x_star;
            Q_cand          = gaussian_to_Q_zero(X_cand, zc);     % zeros hold by construction

            if check_fn(Q_cand)
                X_new(:, k) = x_star;            % accept this column
                break
            end

            % shrink the slice bracket toward theta = 0 (current state, valid)
            if theta < 0, theta_min = theta; else, theta_max = theta; end
            theta = theta_min + (theta_max - theta_min)*rand();
        end
        % if max_shrinks is exhausted, x_k is kept (theta -> 0 limit is valid).
    end

    Q_new = gaussian_to_Q_zero(X_new, zc);
end