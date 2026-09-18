function fig = plot_irf_bands(irfs, weights, varnames, shocknames, varargin)
%PLOT_IRF_BANDS  Posterior median and credible bands of impulse responses.
%
%   fig = plot_irf_bands(irfs, weights, varnames, shocknames)
%   fig = plot_irf_bands(..., 'Shocks', idx, 'Horizon', hmax, 'XLabel', txt)
%
%   irfs        ny x H x K x N array of impulse-response draws, as stored in
%               output.sign_irfs_temp (variable, horizon, shock, draw).
%   weights     N x 1 draw weights, or [] for equal weights. Use
%               output.draw_weight after hars and output.draw_weight_final
%               after harsz.
%   varnames    ny names, one per row of panels.
%   shocknames  K names. Only the entries in 'Shocks' are used.
%
%   'Shocks'    shocks to plot, one column of panels each (default 1:K).
%   'Horizon'   last horizon plotted, counting impact as 0 (default H-1).
%   'XLabel'    label of the horizontal axis (default 'Horizon').
%
%   Each panel shows the weighted posterior median (solid), the 68 percent
%   band (shaded), and the 90 percent band (dashed).

    [ny, H, K, N] = size(irfs);
    p = inputParser;
    addParameter(p, 'Shocks',  1:K);
    addParameter(p, 'Horizon', H - 1);
    addParameter(p, 'XLabel',  'Horizon');
    parse(p, varargin{:});
    shocks = p.Results.Shocks(:)';
    hmax   = min(p.Results.Horizon, H - 1);

    if isempty(weights), weights = ones(N, 1); end
    weights   = weights(:) / sum(weights);
    varnames  = cellstr(varnames);
    shocknames = cellstr(shocknames);

    QUANTILES  = [0.05 0.16 0.50 0.84 0.95];
    LINE_COLOR = [0.12 0.31 0.47];
    x = 0:hmax;

    fig = figure('Units', 'inches', 'Position', [1 1 4.2 * numel(shocks), 1.9 * ny + 0.4]);
    for iv = 1:ny
        for is = 1:numel(shocks)
            draws = squeeze(irfs(iv, 1:hmax + 1, shocks(is), :));          % (hmax+1) x N
            q     = wpercentile(reshape(draws, hmax + 1, N), weights, QUANTILES);

            ax = subplot(ny, numel(shocks), (iv - 1) * numel(shocks) + is);
            hold(ax, 'on');
            ax.Toolbar = [];
            fill(ax, [x fliplr(x)], [q(:, 2)' fliplr(q(:, 4)')], LINE_COLOR, ...
                 'FaceAlpha', 0.13, 'EdgeColor', 'none');
            plot(ax, x, q(:, 1), '--', 'Color', LINE_COLOR, 'LineWidth', 1.3);
            plot(ax, x, q(:, 5), '--', 'Color', LINE_COLOR, 'LineWidth', 1.3);
            plot(ax, x, q(:, 3), '-',  'Color', LINE_COLOR, 'LineWidth', 2.4);
            yline(ax, 0, 'k--', 'LineWidth', 1);
            hold(ax, 'off');

            xlim(ax, [0 hmax]);
            set(ax, 'FontSize', 10);
            title(ax, {varnames{iv}, ['to ' shocknames{shocks(is)}]}, ...
                  'FontWeight', 'bold', 'FontSize', 11, 'Interpreter', 'none');
            ylabel(ax, 'Response', 'FontSize', 10);
            if iv == ny, xlabel(ax, p.Results.XLabel, 'FontSize', 10); end
        end
    end
end
