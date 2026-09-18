function P = wpercentile(X, w, q)
% Weighted percentiles along the draw axis (columns), midpoint convention
% (matches draw_irfs.ipynb wpercentile; NOT the sampler's Type-7 quantile).
% X: (H, N); w: (N,1); q: quantiles in [0,1]. Returns (H, numel(q)).
    [H, ~] = size(X); q = q(:)';
    P = zeros(H, numel(q));
    w = w(:);
    for h = 1:H
        [xs, ord] = sort(X(h, :));
        ws = w(ord);
        cw = cumsum(ws);
        pp = (cw - 0.5 * ws) / cw(end);
        pp = pp(:)'; xs = xs(:)';
        % drop duplicate plotting positions to keep interp1 sample points unique
        [pp_u, iu] = unique(pp, 'stable');
        xs_u = xs(iu);
        qc = min(max(q, pp_u(1)), pp_u(end));   % clamp -> flat extrapolation
        P(h, :) = interp1(pp_u, xs_u, qc, 'linear');
    end
end
