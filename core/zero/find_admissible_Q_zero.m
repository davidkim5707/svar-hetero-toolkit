function [X, Q, success, n_tried] = find_admissible_Q_zero( ...
    check_fn, zc, ny, max_attempts, prescreen_fn)
%FIND_ADMISSIBLE_Q_ZERO  Independent construction-measure rotation repair on
%   the zero manifold. Draws fresh ny x ny standard-Gaussian matrices G, maps
%   each through gaussian_to_Q_zero (the null-space construction, so every
%   zero restriction holds exactly), and keeps the FIRST draw whose rotation
%   passes the full admissibility check.
%
%   [X, Q, success, n_tried] = find_admissible_Q_zero(check_fn, zc, ny, ...
%       max_attempts, prescreen_fn)
%
%   WHY THIS EXISTS: the production global rotation repair must inject a draw
%   from the CORRECT conditional distribution, i.e. the construction measure
%   restricted to the admissible set. Because the proposals here are i.i.d.
%   from the construction measure (independent Gaussian inputs through the
%   same map the chain uses), the first admissible draw has exactly that
%   restricted distribution -- the on-manifold analogue of the "first
%   admissible independent Haar draw" argument used on the no-zero path. A
%   column-by-column search with per-column screening does NOT have this
%   property: conditioning each column on passing its own screen strips
%   state-dependent acceptance factors and distorts the injected law. That
%   search remains valid for INITIALIZATION and seed rescue only.
%
%   LATENT: on success, X is the RAW Gaussian input G, so
%   gaussian_to_Q_zero(X, zc) == Q and (X, Q) is a valid state for the
%   column-wise ESS chain. Do NOT substitute the built columns for X here;
%   the raw input is the exact draw of the Gaussian coordinates.
%
%   OPTIONAL PRESCREEN: prescreen_fn(Q) may implement any cheap NECESSARY
%   condition of check_fn (e.g. an impact-sign screen). Rejecting a WHOLE
%   independent draw on a necessary condition only short-circuits the full
%   check and leaves the first-admissible distribution unchanged, because a
%   draw failing the prescreen would also fail check_fn. It must never pass
%   draws that check_fn would reject more often -- it must be implied by
%   check_fn.
%
%   No zeros (zc empty / zc.has_zero == false): gaussian_to_Q_zero delegates
%   to gaussian_to_Q, so this reproduces independent Haar proposals; the
%   dedicated no-zero routine (find_admissible_Q) remains the default there.
%
%   Requires: gaussian_to_Q_zero.m.

    if nargin < 4 || isempty(max_attempts), max_attempts = 2000; end
    if nargin < 5, prescreen_fn = []; end
    use_prescreen = ~isempty(prescreen_fn);

    success = false;
    X = zeros(ny, ny);
    Q = eye(ny);

    for n_tried = 1:max_attempts
        G  = randn(ny, ny);
        Qc = gaussian_to_Q_zero(G, zc);

        if use_prescreen && ~prescreen_fn(Qc)
            continue
        end

        if check_fn(Qc)
            X = G;
            Q = Qc;
            success = true;
            return
        end
    end
end