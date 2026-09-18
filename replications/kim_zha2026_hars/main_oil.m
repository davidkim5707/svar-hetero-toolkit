%==========================================================================
% main_oil.m -- one-click replication, oil application (HARS)
%
% Pipeline:
%   (1) HARS baseline run (N_theta = 100) -> output/results_oil.mat (slim)
%   (2) HARS ablation run (N_theta = 0)   -> output/results_oil_nostrucredraw.mat (slim)
%   (3) Paper figures  (Figure 2, Figure S2) -> output/figures/
%   (4) Efficiency/diagnostics tables        -> output/
%
% Both runs execute every time. Expect roughly 2x the single-run wall clock.
% This file also runs standalone.
%==========================================================================
clc; clear; close all;
addpath(genpath(fullfile('..','..','code','core')), genpath(fullfile('..','..','code','third_party')));

%% ---- Run 1: HARS baseline (N_theta = 100) ------------------------------
NTHETA = 100;
run('run_oil.m');
save_slim(fullfile('output', 'results_oil.mat'), output, options, run_total_sec);

%% ---- Run 2: HARS ablation (N_theta = 0) --------------------------------
clear; NTHETA = 0;
addpath(genpath(fullfile('..','..','code','core')), genpath(fullfile('..','..','code','third_party')));
run('run_oil.m');
save_slim(fullfile('output', 'results_oil_nostrucredraw.mat'), output, options, run_total_sec);

%% ---- Figures + tables --------------------------------------------------
clear; close all;
addpath(genpath(fullfile('..','..','code','core')), genpath(fullfile('..','..','code','third_party')));
R  = load(fullfile('output', 'results_oil.mat'));
R0 = load(fullfile('output', 'results_oil_nostrucredraw.mat'));
make_oil_figures(R, 'baselines', fullfile('output', 'figures'));
make_oil_tables(R, R0, 'baselines', 'output');


%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function make_oil_figures(R, baselines_dir, fig_dir)
% MATLAB port of oil cells 3 (Figure 2) and 4 (Figure S2) of draw_irfs.ipynb.
% Conventions: bands [5 50 95] (90% CI, per the notebook); NO x100 scaling;
% oil production (var 1) is cumulated for the ESS series (and CMT in Fig S2);
% weighted midpoint wpercentile for HARS; Type-7 unweighted for the baselines.
    if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
    Hn = 18;                                    % horizons plotted (0..17)
    bands = [0.05 0.50 0.95];
    x = 0:(Hn-1);
    varnames   = {'Oil production','Economic activity index','Real oil price'};
    shocknames = {'oil supply shock','aggregate demand shock','oil-specific demand shock'};

    % ---- HARS (red, weighted) : sign_irfs_temp is (3, H, 3, N) --------------
    ess = R.sign_irfs_temp(:, 1:Hn, :, :);      % (3,18,3,N)
    w   = R.draw_weight(:) / sum(R.draw_weight(:));

    % ---- KM baseline (blue) : Draws_IRFs is (vars, shocks, hor, draws) ------
    km_file = fullfile(baselines_dir, 'irfs_kilian_oil.mat');
    have_km = exist(km_file, 'file') == 2;
    if have_km
        K = load(km_file); km = K.Draws_IRFs(:, :, 1:Hn, :);   % (3,3,18,draws)
    else
        fprintf('[SKIP] %s not found -> Figure 2 drawn HARS-only.\n', km_file);
    end

    % ---- CMT baseline (green) for Figure S2 --------------------------------
    cmt_file = fullfile(baselines_dir, 'irfs_carriero_oil_3regimes.mat');
    have_cmt = exist(cmt_file, 'file') == 2;
    if have_cmt
        C = load(cmt_file);
        cmt = cmt_scale_oil(C.RESPONSE, C.inv_A0_all);   % (ndraw, 3, 3, nhor) rescaled
    else
        fprintf('[SKIP] %s not found -> Figure S2 drawn HARS-only.\n', cmt_file);
    end

    % ================= Figure 2 -- KM (blue) vs HARS (red) =================
    fig2 = figure('Name','Figure 2: KM vs HARS (oil, 3-reg)', ...
                  'Units','inches','Position',[1 1 15 7.5]);
    hHomo = []; hHetero = [];
    for si = 1:3
        for vi = 1:3
            ax = subplot(3, 3, (si-1)*3 + vi); hold(ax, 'on'); ax.Toolbar = [];
            eirf = squeeze(ess(vi, :, si, :));          % (18, N)
            if vi == 1, eirf = cumsum(eirf, 1); end      % oil production -> level
            ep = wpercentile(eirf, w, bands);            % (18,3)
            if have_km
                kirf = squeeze(km(vi, si, :, :));        % (18, draws)
                kp = pctile7(kirf, bands);               % (18,3) Type-7
                plot(ax, x, kp(:,1), '--', 'Color',[0 0 1 0.6], 'LineWidth',2);
                plot(ax, x, kp(:,3), '--', 'Color',[0 0 1 0.6], 'LineWidth',2);
                hHomo = plot(ax, x, kp(:,2), '-', 'Color','b', 'LineWidth',2);
            end
            fill(ax, [x fliplr(x)], [ep(:,1)' fliplr(ep(:,3)')], 'r', ...
                 'FaceAlpha',0.25, 'EdgeColor','none');
            hHetero = plot(ax, x, ep(:,2), '-', 'Color','r', 'LineWidth',3);
            yline(ax, 0, 'k--', 'LineWidth',1);
            xlim(ax, [0 Hn-1]); set(ax, 'FontSize', 11);
            title(ax, sprintf('%s to %s', varnames{vi}, shocknames{si}), ...
                  'FontWeight','bold', 'FontSize',13);
            hold(ax, 'off');
        end
    end
    if have_km
        add_top_legend([hHomo hHetero], {'Homoskedastic model','Heteroskedastic model'});
    else
        add_top_legend(hHetero, {'Heteroskedastic model'});
    end
    exportgraphics(fig2, fullfile(fig_dir, 'figure_2.pdf'), 'ContentType','vector');

    % ================= Figure S2 -- HARS (red) vs CMT (green) ==============
    figS2 = figure('Name','Figure S2: HARS vs CMT (oil, 3-reg)', ...
                   'Units','inches','Position',[1 1 15 7.5]);
    hH = []; hC = [];
    for si = 1:3
        for vi = 1:3
            ax = subplot(3, 3, (si-1)*3 + vi); hold(ax, 'on'); ax.Toolbar = [];
            eirf = squeeze(ess(vi, :, si, :));          % (18, N)
            if vi == 1, eirf = cumsum(eirf, 1); end
            ep = wpercentile(eirf, w, bands);
            fill(ax, [x fliplr(x)], [ep(:,1)' fliplr(ep(:,3)')], 'r', ...
                 'FaceAlpha',0.25, 'EdgeColor','none');
            hH = plot(ax, x, ep(:,2), '-', 'Color','r', 'LineWidth',3);
            if have_cmt
                cirf = squeeze(cmt(:, vi, si, :));       % (ndraw, nhor)
                if vi == 1, cirf = cumsum(cirf, 2); end  % cumulate over horizon
                cp = pctile7(cirf(:, 1:Hn)', bands);     % (18,3) Type-7 over draws
                fill(ax, [x fliplr(x)], [cp(:,1)' fliplr(cp(:,3)')], 'g', ...
                     'FaceAlpha',0.25, 'EdgeColor','none');
                hC = plot(ax, x, cp(:,2), '-', 'Color','g', 'LineWidth',3);
            end
            yline(ax, 0, 'k--', 'LineWidth',1);
            xlim(ax, [0 Hn-1]); set(ax, 'FontSize', 11);
            title(ax, sprintf('%s to %s', varnames{vi}, shocknames{si}), ...
                  'FontWeight','bold', 'FontSize',13);
            hold(ax, 'off');
        end
    end
    if have_cmt
        add_top_legend([hH hC], {'HARS','CMT'});
    else
        add_top_legend(hH, {'HARS'});
    end
    exportgraphics(figS2, fullfile(fig_dir, 'figure_S2.pdf'), 'ContentType','vector');
    fprintf('[make_oil_figures] wrote figure_2.pdf, figure_S2.pdf\n');
end


function out = cmt_scale_oil(RESPONSE, inv_A0_all)
% Per-draw unnormalization of the Carriero oil RESPONSE (notebook cell 4):
% scale each shock column by c = [-invA0(:,1,1); invA0(:,2,2); invA0(:,3,3)].
% RESPONSE: (ndraw, nvar, nshock, nhor); inv_A0_all: (ndraw, n, n).
    nd = size(RESPONSE, 1);
    c = zeros(nd, 3);
    c(:,1) = -inv_A0_all(:, 1, 1);
    c(:,2) =  inv_A0_all(:, 2, 2);
    c(:,3) =  inv_A0_all(:, 3, 3);
    out = RESPONSE .* reshape(c, [nd, 1, 3, 1]);
end


function make_oil_tables(R, R0, baselines_dir, out_dir)
% mESS (VFJ 2019) + HARS repair diagnostics. Subset: vars 1:3 x horizons 1:12
% x shocks 1:3 (p = 108), per compute_irf_mess_oil.m. Rows: HARS(100), HARS(0),
% CMT oil (manual 5112.00 s, 1,000-draw full run, measured).
    vars_use = 1:3; horizons_use = 1:12; shocks_use = 1:3;

    rows = struct('name',{},'kind',{},'N',{},'mess',{},'mess_per_N',{}, ...
                  'sec_eff',{},'hr1000',{},'tsrc',{},'stale',{},'p_post',{},'p_rot',{},'disc',{});
    rows(end+1) = oil_hars_row('HARS (N_theta = 100)', R,  vars_use, horizons_use, shocks_use);
    rows(end+1) = oil_hars_row('HARS (N_theta = 0)',   R0, vars_use, horizons_use, shocks_use);

    cmt_file = fullfile(baselines_dir, 'irfs_carriero_oil_3regimes.mat');
    if exist(cmt_file, 'file') == 2
        C = load(cmt_file);
        irfs = permute(C.RESPONSE, [2 4 3 1]);       % (N,ny,K,H) -> (ny,H,K,N)
        rows(end+1) = oil_baseline_row('CMT (Carriero, 3-reg, 1000 draws)', irfs, 5112.00, ...
                                       vars_use, horizons_use, shocks_use);
    else
        fprintf('[SKIP] %s not found -> CMT oil row dropped.\n', cmt_file);
    end

    NOTE = ['Note: the CMT timing is measured on the authors'' hardware (5112.00 s, ', ...
            '1,000-draw full run); HARS rows are timed on this machine.'];
    print_and_write_tables(rows, out_dir, 'oil', NOTE, 'S3', 'S4');
end


function r = oil_hars_row(name, S, vars_use, horizons_use, shocks_use)
    [m, N] = mess_of_irfs(S.sign_irfs_temp, vars_use, horizons_use, shocks_use);
    t = S.t_wall; dg = S.ess_diagnostics;
    stale  = nzf(dg,'stale_q_events');  byPost = nzf(dg,'stale_resolved_structural');
    byRot  = nzf(dg,'stale_resolved_fallback');  disc = nzf(dg,'stale_unresolved');
    r = struct('name',name,'kind','hars','N',N,'mess',m,'mess_per_N',m/N, ...
               'sec_eff',t/m,'hr1000',(t/m)*1000/3600,'tsrc','[prof]', ...
               'stale',stale,'p_post',pctv(byPost,stale),'p_rot',pctv(byRot,stale),'disc',disc);
end


function r = oil_baseline_row(name, irfs, manual_sec, vars_use, horizons_use, shocks_use)
    [m, N] = mess_of_irfs(irfs, vars_use, horizons_use, shocks_use);
    r = struct('name',name,'kind','baseline','N',N,'mess',m,'mess_per_N',m/N, ...
               'sec_eff',manual_sec/m,'hr1000',(manual_sec/m)*1000/3600,'tsrc','[man]', ...
               'stale',NaN,'p_post',NaN,'p_rot',NaN,'disc',NaN);
end


function [m, N] = mess_of_irfs(irfs, vars, horizons, shock)
% irfs: (ny, H, K, N). Reshape the requested subset to (N x p), mESS.
    [ny, Ht, Kt, N] = size(irfs);
    vars = vars(vars <= ny); horizons = horizons(horizons <= Ht); shock = shock(shock <= Kt);
    Y = reshape(irfs(vars, horizons, shock, :), [], N)';
    m = mess_vfj(Y);
end


function P = pctile7(X, q)
% numpy Type-7 unweighted percentiles along the draw axis (columns).
% X: (H, N); q in [0,1]. Returns (H, numel(q)).
    [H, N] = size(X); q = q(:)';
    pp = (0:N-1) / (N - 1);
    P = zeros(H, numel(q));
    for h = 1:H
        xs = sort(X(h, :));
        P(h, :) = interp1(pp, xs, q, 'linear');
    end
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


function print_and_write_tables(rows, out_dir, tag, NOTE, s_hours, s_ablation)
% Console + .txt + .csv + .tex (Table <s_hours> hours/1000, Table <s_ablation> ablation).
    txt_file = fullfile(out_dir, sprintf('efficiency_%s.txt', tag));
    fid = fopen(txt_file, 'w');
    pr = @(varargin) print_both(fid, varargin{:});
    pr('\nEfficiency / diagnostics -- %s application (HARS)\n', tag);
    pr('Vats-Flegal-Jones (2019) mESS.\n');
    pr('Runtime source: [prof]=stored sampling-loop time (burn-in excluded); [man]=authors'' hardware.\n');
    pr('\n=== (A) Efficiency: mESS and time per effective draw ===\n');
    pr('%-40s %8s %10s %10s %12s %16s %7s\n','Sampler','N','mESS','mESS/N','sec/eff','hr per 1000 eff','tsrc');
    pr('%s\n', repmat('-',1,106));
    for k = 1:numel(rows)
        r = rows(k);
        pr('%-40s %8d %10.1f %10.4f %12.3f %16.3f %7s\n', ...
           r.name, r.N, r.mess, r.mess_per_N, r.sec_eff, r.hr1000, r.tsrc);
    end
    pr('\n=== (B) Stale-Q repair diagnostics ((A0,Lambda) channel; HARS only) ===\n');
    pr('%-40s %10s %18s %16s %12s\n','Sampler','stale-Q','by posterior(%)','by rotation(%)','discarded');
    pr('%s\n', repmat('-',1,100));
    for k = 1:numel(rows)
        r = rows(k);
        if strcmp(r.kind,'hars')
            pr('%-40s %10d %18.1f %16.1f %12d\n', r.name, r.stale, r.p_post, r.p_rot, r.disc);
        else
            pr('%-40s %10s %18s %16s %12s\n', r.name, '--','--','--','--');
        end
    end
    pr('\n%s\n', NOTE);
    fclose(fid);
    type_file(txt_file);

    csv_file = fullfile(out_dir, sprintf('efficiency_%s.csv', tag));
    fc = fopen(csv_file, 'w');
    fprintf(fc, ['sampler,N,mESS,mESS_per_N,sec_per_eff_draw,hours_per_1000_eff,', ...
                 'stale_q_events,resolved_posterior_pct,resolved_rotation_pct,discarded_transitions,time_source\n']);
    for k = 1:numel(rows)
        r = rows(k);
        if strcmp(r.kind,'hars')
            fprintf(fc, '%s,%d,%.4f,%.6f,%.4f,%.6f,%d,%.1f,%.1f,%d,%s\n', ...
                cq(r.name), r.N, r.mess, r.mess_per_N, r.sec_eff, r.hr1000, r.stale, r.p_post, r.p_rot, r.disc, r.tsrc);
        else
            fprintf(fc, '%s,%d,%.4f,%.6f,%.4f,%.6f,,,,,%s\n', ...
                cq(r.name), r.N, r.mess, r.mess_per_N, r.sec_eff, r.hr1000, r.tsrc);
        end
    end
    fclose(fc);
    fprintf('[tables] wrote %s\n', csv_file);

    write_tex(fullfile(out_dir, sprintf('table_%s_%s_%s.tex', s_hours, s_ablation, tag)), rows, NOTE, tag, s_hours, s_ablation);
end


function print_both(fid, varargin)
    fprintf(varargin{:}); fprintf(fid, varargin{:});
end

function type_file(fname)
    fid = fopen(fname,'r'); if fid<0, return; end
    while true, ln = fgetl(fid); if ~ischar(ln), break; end, fprintf('%s\n', ln); end
    fclose(fid);
end

function s = cq(s)
    if contains(s, ','), s = ['"' s '"']; end
end


function write_tex(fname, rows, NOTE, tag, s_hours, s_ablation)
    getrow = @(pred) find(arrayfun(pred, rows), 1);
    iH100 = getrow(@(r) strcmp(r.name,'HARS (N_theta = 100)'));
    iH0   = getrow(@(r) strcmp(r.name,'HARS (N_theta = 0)'));
    iCMT  = getrow(@(r) startsWith(r.name,'CMT'));
    fid = fopen(fname,'w');
    fprintf(fid, '%% Auto-generated by main_%s.m. Requires \\usepackage{booktabs}.\n\n', tag);
    fprintf(fid, '\\begin{table}[htbp]\\centering\n');
    fprintf(fid, '\\caption{Sampling efficiency: hours per 1{,}000 effective draws (%s).}\n', tag);
    fprintf(fid, '\\label{tab:%s_%s}\n', s_hours, tag);
    fprintf(fid, '\\begin{tabular}{lrr}\n\\toprule\n');
    fprintf(fid, 'Sampler & $N$ & Hours / 1{,}000 eff. \\\\\n\\midrule\n');
    if ~isempty(iH100), r=rows(iH100); fprintf(fid,'HARS ($N_\\theta=100$) & %d & %.3f \\\\\n', r.N, r.hr1000); end
    if ~isempty(iCMT),  r=rows(iCMT);  fprintf(fid,'CMT & %d & %.3f \\\\\n', r.N, r.hr1000); end
    fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
    fprintf(fid, '\\par\\footnotesize %s\n', tex_escape(NOTE));
    fprintf(fid, '\\end{table}\n\n');
    fprintf(fid, '\\begin{table}[htbp]\\centering\n');
    fprintf(fid, '\\caption{Structural-redraw ablation ($N_\\theta$): repair channel (%s).}\n', tag);
    fprintf(fid, '\\label{tab:%s_%s}\n', s_ablation, tag);
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
