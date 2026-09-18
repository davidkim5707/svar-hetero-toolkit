function pass = check_narrative_regime_avg(Q, cache, var, NRR, fsign_cell, fsign, ir_cell, xi_override)
%CHECK_NARRATIVE_REGIME_AVG  Regime/window-average narrative restrictions.
%
%   [XI PATCH 07/2026] The labeled shocks are eps_tilde = Q * xi with
%   xi = Lambda^{-1/2} * eps (standardized structural shocks). udraw_true
%   stores eps on the ORIGINAL scale (regime variance lambda), so this
%   function standardizes BY REGIME before rotating:
%
%       xi(win,:) = udraw_true(win,:) ./ sqrt(lambda(:, r))'
%       v(win,:)  = xi(win,:) * Q'          (then the fsign flip)
%
%   With the regime IRFs carrying Dhalf = sqrt(lambda) (local_regime_smat),
%   this restores the forecast-error identity sum_k HD_t(k) = u_{v,t} and
%   makes the checked statistic exactly the paper's H_{k->v,t}.
%
%   Optional 8th argument xi_override: a T x ny matrix of STANDARDIZED
%   shocks that replaces the realized xi rows (used by estimate_nrr_omega;
%   rows outside the restricted dates may be NaN). With nargin < 8 the
%   realized shocks are used.
%
%   Checks per NRR entry (each independently toggleable):
%     (1) FOCAL SIGN: mean of the (flipped) focal labeled shock over the
%         window has the specified sign.
%     (2) FOCAL DOMINANCE (max / min / strict) on mean |HD| over the window.
%     (3) ANTI DOMINANCE (min) on mean |HD| over the window.

    if isempty(NRR), pass = true; return; end

    nlambda = cache.nlambda;
    u_true  = var.udraw_true;                    % T x ny, original scale
    [T_u, ny] = size(u_true);
    use_override = (nargin >= 8) && ~isempty(xi_override);

    for k = 1:numel(NRR)
        res = NRR(k);
        win = res.window(:)';
        r   = res.regime_idx;

        % --- Validate window and regime ---
        if isempty(win) || ~isfinite(r) || r < 1 || r > max(nlambda, 1)
            pass = false; return
        end
        if max(win) > T_u || min(win) < 1
            pass = false; return
        end

        % --- Regime-specific sign flip ---
        if nlambda > 1
            if r > numel(fsign_cell) || isempty(fsign_cell{r})
                pass = false; return
            end
            fs = fsign_cell{r};
        else
            fs = fsign;
        end
        if isempty(fs), pass = false; return; end

        % --- Standardized labeled shocks on the window ---
        lam_r = sqrt(max(cache.lambda(:, r), 0))';   % 1 x ny
        if any(lam_r <= 0), pass = false; return; end
        if use_override
            xi_win = xi_override(win, :);
        else
            xi_win = u_true(win, :) ./ lam_r;
        end
        v_win = (xi_win * Q') .* reshape(fs, 1, []);   % window rows only

        % --- Determine which checks are active ---
        has_focal = isfield(res, 'shock_idx')      && ~isempty(res.shock_idx);
        has_anti  = isfield(res, 'anti_shock_idx') && ~isempty(res.anti_shock_idx);

        if has_focal && isfield(res, 'dominance') && ~isempty(res.dominance)
            focal_mode = lower(res.dominance);
        else
            focal_mode = 'off';
        end
        do_focal_dom = has_focal && ~strcmp(focal_mode, 'off');
        do_anti_dom  = has_anti  && isfield(res, 'anti_dominance') && strcmpi(res.anti_dominance, 'on');

        % [OPTION] nrr_dom window aggregation of |HD|: 'mean_abs' (default) | 'cum_abs'
        %   'mean_abs' : per-month |.| then average -> (1/|W|) sum_t |H_{j,t}|
        %   'cum_abs'  : net over W then |.|        -> |sum_t H_{j,t}|  (ADRR interval)
        if isfield(res, 'dom_agg') && ~isempty(res.dom_agg)
            dom_agg = lower(res.dom_agg);
        else
            dom_agg = 'mean_abs';
        end

        % --- (1) Sign check on focal labeled shock ---
        if has_focal
            j   = res.shock_idx;
            sgn = res.sign;
            if (sgn * mean(v_win(:, j))) <= 0
                pass = false; return
            end
        end

        % --- (2) HD-based dominance / anti-dominance ---
        if do_focal_dom || do_anti_dom
            if r > numel(ir_cell) || isempty(ir_cell{r})
                pass = false; return
            end
            if ~isfield(res, 'variable_idx') || isempty(res.variable_idx)
                error('check_narrative_regime_avg: dominance requires variable_idx.');
            end
            v_idx = res.variable_idx;
            ir_r  = ir_cell{r};

            switch dom_agg
                case {'mean_abs', 'cum_abs'}
                    % one-step (impact-date) contribution per month:
                    %   H_t = Theta_0(v,j) * eps_tilde_{j,t}
                    T_win  = size(v_win, 1);
                    HD_win = zeros(T_win, ny);           % SIGNED per-month HD
                    for ii = 1:T_win
                        HD_t = getHDs_fast(ir_r, v_win(ii, :), v_idx);
                        HD_win(ii, :) = HD_t(:)';
                    end
                    if strcmp(dom_agg, 'mean_abs')       % |.| per month, then mean
                        mean_abs_HD = mean(abs(HD_win), 1);      % 1 x ny
                    else                                 % 'cum_abs': net one-step over W, then |.|
                        mean_abs_HD = abs(sum(HD_win, 1));       % 1 x ny
                    end
                case 'interval_hd'
                    % ADRR multi-horizon interval HD: cumulative contribution of
                    % each shock to the window ENDPOINT, accumulating the delayed
                    % dynamics (Theta_0, Theta_1, ...) of every shock in the window.
                    % getHDs_fast on the full block returns Sum_h Theta_h * eps_{end-h}.
                    % Accumulation depth is capped at the horizons stored in ir_r.
                    Hir    = size(ir_r, 2);              % horizons h = 0..Hir-1
                    dep    = min(size(v_win, 1), Hir);   % accumulation depth
                    v_used = v_win(end-dep+1:end, :);    % last `dep` window months
                    HD_int = getHDs_fast(ir_r, v_used, v_idx);   % ny x 1
                    mean_abs_HD = abs(HD_int(:))';               % 1 x ny
                otherwise
                    error('check_narrative_regime_avg: unknown dom_agg "%s".', dom_agg);
            end

            if do_focal_dom
                j = res.shock_idx;
                if isfield(res, 'other_idx') && ~isempty(res.other_idx)
                    others = res.other_idx;
                else
                    others = setdiff(1:ny, j);
                end
                switch focal_mode
                    case {'on', 'max'}
                        if any(mean_abs_HD(j) <= mean_abs_HD(others)), pass = false; return; end
                    case 'min'
                        if any(mean_abs_HD(j) >= mean_abs_HD(others)), pass = false; return; end
                    case 'strict'
                        if mean_abs_HD(j) <= sum(mean_abs_HD(others)), pass = false; return; end
                    otherwise
                        error('check_narrative_regime_avg: unknown dominance mode "%s".', focal_mode);
                end
            end

            if do_anti_dom
                a = res.anti_shock_idx;
                if isfield(res, 'anti_other_idx') && ~isempty(res.anti_other_idx)
                    anti_others = res.anti_other_idx;
                else
                    anti_others = setdiff(1:ny, a);
                end
                if any(mean_abs_HD(a) >= mean_abs_HD(anti_others)), pass = false; return; end
            end
        end
    end

    pass = true;
end