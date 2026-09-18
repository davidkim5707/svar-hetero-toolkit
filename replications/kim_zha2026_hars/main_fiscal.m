%==========================================================================
% main_fiscal.m -- one-click replication, fiscal application (HARS)
%
% Pipeline:
%   (1) HARS homoskedastic run (REGIME_ON=0, N_theta=100) -> results_fiscal_homosk.mat
%   (2) HARS 3-regime baseline (REGIME_ON=1, N_theta=100)  -> results_fiscal.mat
%   (3) HARS 3-regime ablation (REGIME_ON=1, N_theta=0)    -> results_fiscal_nostrucredraw.mat
%   (4) Paper figures (Figure 3 net-tax, Figure 4 spending) -> output/figures/
%   (5) Efficiency/diagnostics tables                       -> output/
%
% Three HARS runs execute every time. This file also runs standalone.
%==========================================================================
clc; clear; close all;
addpath(genpath(fullfile('..','..','core')), genpath(fullfile('..','..','third_party')));

%% ---- Run 1: homoskedastic (REGIME_ON = 0, N_theta = 100) ---------------
REGIME_ON = 0; NTHETA = 100;
run('run_fiscal.m');
save_slim(fullfile('output', 'results_fiscal_homosk.mat'), output, options, run_total_sec);

%% ---- Run 2: 3-regime baseline (REGIME_ON = 1, N_theta = 100) -----------
clear; REGIME_ON = 1; NTHETA = 100;
addpath(genpath(fullfile('..','..','core')), genpath(fullfile('..','..','third_party')));
run('run_fiscal.m');
save_slim(fullfile('output', 'results_fiscal.mat'), output, options, run_total_sec);

%% ---- Run 3: 3-regime ablation (REGIME_ON = 1, N_theta = 0) -------------
clear; REGIME_ON = 1; NTHETA = 0;
addpath(genpath(fullfile('..','..','core')), genpath(fullfile('..','..','third_party')));
run('run_fiscal.m');
save_slim(fullfile('output', 'results_fiscal_nostrucredraw.mat'), output, options, run_total_sec);

%% ---- Figures + tables --------------------------------------------------
clear; close all;
addpath(genpath(fullfile('..','..','core')), genpath(fullfile('..','..','third_party')));
Rh = load(fullfile('output', 'results_fiscal_homosk.mat'));
R  = load(fullfile('output', 'results_fiscal.mat'));
R0 = load(fullfile('output', 'results_fiscal_nostrucredraw.mat'));
make_fiscal_figures(Rh, R, fullfile('output', 'figures'));
make_fiscal_tables(R, R0, 'output');


%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function make_fiscal_figures(Rh, R, fig_dir)
% MATLAB port of fiscal cells 6 (Figure 3, net-tax shock col 1) and 7
% (Figure 4, spending shock col 2). Both series weighted (the homoskedastic
% baseline is also an ESS run), both filtered by plot_keep_mask, x100 scaling,
% bands [5 50 95] (90% CI). Blue = homoskedastic, red = 3-regime.
    if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
    fiscal_shock_figure(Rh, R, 1, fig_dir, 'figure_3.pdf');   % net-tax shock
    fiscal_shock_figure(Rh, R, 2, fig_dir, 'figure_4.pdf');   % spending shock
    fprintf('[make_fiscal_figures] wrote figure_3.pdf, figure_4.pdf\n');
end


function fiscal_shock_figure(Rh, R, shock_col, fig_dir, fname)
    varnames = {'Net taxes','Gov spending','GDP','3m T-bill','GDP deflator'};
    bands = [0.05 0.50 0.95];
    Hn = 20; x = 0:(Hn-1);
    [bi, bw] = fiscal_filter(Rh);     % homoskedastic baseline (blue)
    [ei, ew] = fiscal_filter(R);      % 3-regime (red)

    base_pct = cell(5,1); ess_pct = cell(5,1); ylims = zeros(5,2);
    for i = 1:5
        bp = 100 * wpercentile(squeeze(bi(i, 1:Hn, shock_col, :)), bw, bands);
        ep = 100 * wpercentile(squeeze(ei(i, 1:Hn, shock_col, :)), ew, bands);
        base_pct{i} = bp; ess_pct{i} = ep;
        lo = min(min(bp(:,1)), min(ep(:,1)));
        hi = max(max(bp(:,3)), max(ep(:,3)));
        m = (hi - lo) * 0.1; ylims(i,:) = [lo-m hi+m];
    end

    fig = figure('Name', fname, 'Units','inches','Position',[1 1 13 10]);
    hB = []; hR = [];
    for i = 1:5
        ax = subplot(3, 2, i); hold(ax, 'on'); ax.Toolbar = [];
        bp = base_pct{i}; ep = ess_pct{i};
        plot(ax, x, bp(:,1), '--', 'Color',[0 0 1 0.7], 'LineWidth',2);
        plot(ax, x, bp(:,3), '--', 'Color',[0 0 1 0.7], 'LineWidth',2);
        hB = plot(ax, x, bp(:,2), '-', 'Color','b', 'LineWidth',2);
        fill(ax, [x fliplr(x)], [ep(:,1)' fliplr(ep(:,3)')], 'r', ...
             'FaceAlpha',0.25, 'EdgeColor','none');
        hR = plot(ax, x, ep(:,2), '-', 'Color','r', 'LineWidth',2.5);
        yline(ax, 0, 'k--', 'LineWidth',1);
        xlim(ax, [0 Hn-1]); ylim(ax, ylims(i,:));
        set(ax, 'XTick', 0:4:Hn, 'FontSize', 10);
        title(ax, varnames{i}, 'FontWeight','bold', 'FontSize',13);
        hold(ax, 'off');
    end
    ax6 = subplot(3, 2, 6); set(ax6, 'Visible', 'off');   % hide unused 6th panel
    add_top_legend([hB hR], {'Homoskedastic model','Heteroskedastic model'});
    exportgraphics(fig, fullfile(fig_dir, fname), 'ContentType','vector');
end


function [irfs, w] = fiscal_filter(R)
% Apply plot_keep_mask (per the notebook loader) and normalize weights.
% irfs: (5, H, 2, N_kept); w: (N_kept, 1) summing to 1.
    irfs = R.sign_irfs_temp;
    N = size(irfs, 4);
    if isfield(R, 'plot_keep_mask') && numel(R.plot_keep_mask) == N
        keep = logical(R.plot_keep_mask(:));
    else
        keep = true(N, 1);
    end
    irfs = irfs(:, :, :, keep);
    w = R.draw_weight(:); w = w(keep); w = w / sum(w);
end


function make_fiscal_tables(R, R0, out_dir)
% mESS (VFJ 2019) over vars 1:5 x horizons 1:12 x shock column 2 (p = 60),
% per compute_irf_mess_fiscal.m. Rows: HARS(100), HARS(0). The CMT row FAILS
% (1 admissible draw in 200,000 candidates; acceptance 0.0005%; no timing).
    vars_use = 1:5; horizons_use = 1:12; shock_use = 2;
    rows = struct('name',{},'kind',{},'N',{},'mess',{},'mess_per_N',{}, ...
                  'sec_eff',{},'hr1000',{},'tsrc',{},'stale',{},'p_post',{},'p_rot',{},'disc',{});
    rows(end+1) = fiscal_hars_row('HARS (N_theta = 100)', R,  vars_use, horizons_use, shock_use);
    rows(end+1) = fiscal_hars_row('HARS (N_theta = 0)',   R0, vars_use, horizons_use, shock_use);

    CMT_NOTE = ['CMT (Carriero, 3-reg): FAILS -- 1 admissible draw in 200,000 candidates ', ...
                '(acceptance 0.0005%); no timing reported.'];
    NOTE = ['Note: HARS rows are timed on this machine (sampling-loop wall time). ', CMT_NOTE];

    txt_file = fullfile(out_dir, 'efficiency_fiscal.txt');
    fid = fopen(txt_file, 'w');
    pr = @(varargin) print_both(fid, varargin{:});
    pr('\nEfficiency / diagnostics -- fiscal application (HARS)\n');
    pr('Vats-Flegal-Jones (2019) mESS; subset vars=1..5, horizons=1..12, shock col 2 (p=60).\n');
    pr('\n=== (A) Efficiency: mESS and time per effective draw ===\n');
    pr('%-40s %8s %10s %10s %12s %16s %7s\n','Sampler','N','mESS','mESS/N','sec/eff','hr per 1000 eff','tsrc');
    pr('%s\n', repmat('-',1,106));
    for k = 1:numel(rows)
        r = rows(k);
        pr('%-40s %8d %10.1f %10.4f %12.3f %16.3f %7s\n', ...
           r.name, r.N, r.mess, r.mess_per_N, r.sec_eff, r.hr1000, r.tsrc);
    end
    pr('%-40s %s\n', 'CMT (Carriero, 3-reg)', 'FAILS (1 admissible draw / 200,000; acceptance 0.0005%)');
    pr('\n=== (B) Stale-Q repair diagnostics ((A0,Lambda) channel; HARS only) ===\n');
    pr('%-40s %10s %18s %16s %12s\n','Sampler','stale-Q','by posterior(%)','by rotation(%)','discarded');
    pr('%s\n', repmat('-',1,100));
    for k = 1:numel(rows)
        r = rows(k);
        pr('%-40s %10d %18.1f %16.1f %12d\n', r.name, r.stale, r.p_post, r.p_rot, r.disc);
    end
    pr('\n%s\n', NOTE);
    fclose(fid);
    type_file(txt_file);

    csv_file = fullfile(out_dir, 'efficiency_fiscal.csv');
    fc = fopen(csv_file, 'w');
    fprintf(fc, ['sampler,N,mESS,mESS_per_N,sec_per_eff_draw,hours_per_1000_eff,', ...
                 'stale_q_events,resolved_posterior_pct,resolved_rotation_pct,discarded_transitions,time_source\n']);
    for k = 1:numel(rows)
        r = rows(k);
        fprintf(fc, '%s,%d,%.4f,%.6f,%.4f,%.6f,%d,%.1f,%.1f,%d,%s\n', ...
            r.name, r.N, r.mess, r.mess_per_N, r.sec_eff, r.hr1000, r.stale, r.p_post, r.p_rot, r.disc, r.tsrc);
    end
    fprintf(fc, 'CMT (Carriero 3-reg),,,,,,,,,,FAILS\n');
    fclose(fc);
    fprintf('[tables] wrote %s\n', csv_file);

    write_fiscal_tex(fullfile(out_dir, 'table_S5_S6_fiscal.tex'), rows, CMT_NOTE);
end


function r = fiscal_hars_row(name, S, vars_use, horizons_use, shock_use)
    [m, N] = mess_of_irfs(S.sign_irfs_temp, vars_use, horizons_use, shock_use);
    t = S.t_wall; dg = S.ess_diagnostics;
    stale  = nzf(dg,'stale_q_events');  byPost = nzf(dg,'stale_resolved_structural');
    byRot  = nzf(dg,'stale_resolved_fallback');  disc = nzf(dg,'stale_unresolved');
    r = struct('name',name,'kind','hars','N',N,'mess',m,'mess_per_N',m/N, ...
               'sec_eff',t/m,'hr1000',(t/m)*1000/3600,'tsrc','[prof]', ...
               'stale',stale,'p_post',pctv(byPost,stale),'p_rot',pctv(byRot,stale),'disc',disc);
end


function [m, N] = mess_of_irfs(irfs, vars, horizons, shock)
    [ny, Ht, Kt, N] = size(irfs);
    vars = vars(vars <= ny); horizons = horizons(horizons <= Ht); shock = shock(shock <= Kt);
    Y = reshape(irfs(vars, horizons, shock, :), [], N)';
    m = mess_vfj(Y);
end


function v = nzf(d, fld)
    if isfield(d, fld) && ~isempty(d.(fld)) && isfinite(d.(fld)(1))
        v = double(d.(fld)(1));
    else
        v = 0;
    end
end

function p = pctv(num, den)
    if den > 0, p = 100 * num / den; else, p = NaN; end
end


function add_top_legend(handles, labels)
    lg = legend(handles, labels, 'Orientation','horizontal', 'Box','off', 'FontSize',12);
    lg.Units = 'normalized'; p = lg.Position;
    lg.Position = [0.5 - p(3)/2, 0.965, p(3), p(4)];
end


function print_both(fid, varargin)
    fprintf(varargin{:}); fprintf(fid, varargin{:});
end

function type_file(fname)
    fid = fopen(fname,'r'); if fid<0, return; end
    while true, ln = fgetl(fid); if ~ischar(ln), break; end, fprintf('%s\n', ln); end
    fclose(fid);
end


function write_fiscal_tex(fname, rows, CMT_NOTE)
    getrow = @(pred) find(arrayfun(pred, rows), 1);
    iH100 = getrow(@(r) strcmp(r.name,'HARS (N_theta = 100)'));
    iH0   = getrow(@(r) strcmp(r.name,'HARS (N_theta = 0)'));
    fid = fopen(fname,'w');
    fprintf(fid, '%% Auto-generated by main_fiscal.m. Requires \\usepackage{booktabs}.\n\n');
    fprintf(fid, '\\begin{table}[htbp]\\centering\n');
    fprintf(fid, '\\caption{Sampling efficiency: hours per 1{,}000 effective draws (fiscal).}\n');
    fprintf(fid, '\\label{tab:S5_fiscal}\n');
    fprintf(fid, '\\begin{tabular}{lrr}\n\\toprule\n');
    fprintf(fid, 'Sampler & $N$ & Hours / 1{,}000 eff. \\\\\n\\midrule\n');
    if ~isempty(iH100), r=rows(iH100); fprintf(fid,'HARS ($N_\\theta=100$) & %d & %.3f \\\\\n', r.N, r.hr1000); end
    fprintf(fid, 'CMT (Carriero, 3-reg) & \\multicolumn{2}{c}{fails} \\\\\n');
    fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
    fprintf(fid, '\\par\\footnotesize %s\n', tex_escape(CMT_NOTE));
    fprintf(fid, '\\end{table}\n\n');
    fprintf(fid, '\\begin{table}[htbp]\\centering\n');
    fprintf(fid, '\\caption{Structural-redraw ablation ($N_\\theta$): repair channel (fiscal).}\n');
    fprintf(fid, '\\label{tab:S6_fiscal}\n');
    fprintf(fid, '\\begin{tabular}{lrrrrr}\n\\toprule\n');
    fprintf(fid, '$N_\\theta$ & Hrs / 1{,}000 & Discarded & Stale-Q & by posterior (\\%%) & by rotation (\\%%) \\\\\n\\midrule\n');
    for ii = [iH100 iH0]
        if isempty(ii), continue; end
        r = rows(ii); ntheta = 100; if contains(r.name,'= 0'), ntheta = 0; end
        fprintf(fid, '%d & %.3f & %d & %d & %.1f & %.1f \\\\\n', ntheta, r.hr1000, r.disc, r.stale, r.p_post, r.p_rot);
    end
    fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n');
    fclose(fid);
    fprintf('[tables] wrote %s\n', fname);
end


function s = tex_escape(s)
    s = strrep(s, '%', '\%'); s = strrep(s, '&', '\&'); s = strrep(s, '_', '\_');
end
