function save_slim(fname, output, options, run_total_sec)
% Keep only what the figures/tables need. Default save format (no -v7.3;
% ~30 MB). To add a field later (FEVD, Lambda, ...), add one S.<field> line.
    S = struct();
    S.sign_irfs_temp = output.sign_irfs_temp;          % (ny x H x K x N) IRF draws

    N = size(output.sign_irfs_temp, 4);
    if isfield(output, 'draw_weight') && ~isempty(output.draw_weight)
        S.draw_weight = output.draw_weight(:);
    else
        warning('save_slim:noWeights', ...
            'output.draw_weight missing/empty -> using uniform ones(%d,1).', N);
        S.draw_weight = ones(N, 1);
    end

    if isfield(output, 'prof') && isfield(output.prof, 't_wall')
        S.t_wall = output.prof.t_wall;                 % sampling-loop seconds (burn-in excluded)
    else
        warning('save_slim:noTwall', 'output.prof.t_wall missing -> storing NaN.');
        S.t_wall = NaN;
    end

    S.run_total_sec  = run_total_sec;                  % end-to-end wall time (s)
    S.ess_diagnostics = output.ess_diagnostics;        % stale-Q breakdown, repair counters
    S.options        = options;                        % full config snapshot
    if isfield(output, 'plot_keep_mask')               % fiscal figures filter draws by this
        S.plot_keep_mask = output.plot_keep_mask;
    end

    save(fname, '-struct', 'S');
    fprintf('[save_slim] wrote %s  (N = %d draws)\n', fname, N);
end
