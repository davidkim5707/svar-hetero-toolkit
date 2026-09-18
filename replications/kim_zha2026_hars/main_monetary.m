%==========================================================================
% main_monetary.m -- one-click replication, monetary application (HARS)
%
% Pipeline:
%   (1) HARS baseline run (N_theta = 100) -> output/results_monetary.mat (slim)
%   (2) HARS ablation run (N_theta = 0)   -> output/results_monetary_nostrucredraw.mat (slim)
%   (3) Paper figures  (displayed + saved to output/figures/)
%   (4) Efficiency/diagnostics tables (printed + saved to output/)
%
% Both runs execute every time: the hours-per-1,000-effective-draws
% comparison is only valid when both HARS runs are timed on the same
% machine. Expect roughly 2x the single-run wall clock.
%==========================================================================
clc; clear; close all;
addpath(genpath(fullfile('..','..','core')), genpath(fullfile('..','..','third_party')));

%% ---- Run 1: HARS baseline (N_theta = 100) ------------------------------
NTHETA = 100;
run('run_monetary.m');
save_slim(fullfile('output', 'results_monetary.mat'), output, options, run_total_sec);

%% ---- Run 2: HARS ablation (N_theta = 0) --------------------------------
clear; NTHETA = 0;
addpath(genpath(fullfile('..','..','core')), genpath(fullfile('..','..','third_party')));
run('run_monetary.m');
save_slim(fullfile('output', 'results_monetary_nostrucredraw.mat'), output, options, run_total_sec);

%% ---- Figures + tables --------------------------------------------------
clear; close all;
addpath(genpath(fullfile('..','..','core')), genpath(fullfile('..','..','third_party')));
R  = load(fullfile('output', 'results_monetary.mat'));
R0 = load(fullfile('output', 'results_monetary_nostrucredraw.mat'));
make_figures(R, 'baselines', fullfile('output', 'figures'));
make_efficiency_tables(R, R0, 'baselines', 'output');


%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function make_figures(R, baselines_dir, fig_dir)
% MATLAB port of the monetary cells of irfs/draw_irfs.ipynb.
% Cosmetic approximations (documented): '--' dash approximates the notebook's
% (0,(2.5,3)) pattern; default (tex) title interpreter + bold weight instead of
% usetex+newtx. Data content matches the notebook exactly (Type-7 percentiles
% for the unweighted baselines, midpoint weighted percentiles for HARS, x100).
    if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

    varnames = {'GDP','GDP deflator','Commodity price index', ...
                'Total reserve','Non-borrowed reserves','Fed funds rate'};
    bands = [0.05 0.50 0.95];
    H = 61; x = 0:(H-1);

    % ---- HARS (red, weighted) -------------------------------------------
    ess = squeeze(R.sign_irfs_temp);                 % (6,61,N)
    w   = R.draw_weight(:) / sum(R.draw_weight(:));
    ess_share = 1 / (numel(w) * sum(w.^2));
    fprintf('[make_figures] HARS draws N = %d, ESS share = %.4f\n', numel(w), ess_share);
    ess_pct = cell(6,1);
    for i = 1:6
        ess_pct{i} = 100 * wpercentile(squeeze(ess(i,:,:)), w, bands);   % (61,3)
    end

    % ---- AR 2018 baseline (blue, unweighted) ----------------------------
    ar_file = fullfile(baselines_dir, 'irfs_ramirez_monetary.mat');
    have_ar = exist(ar_file, 'file') == 2;
    ar_pct = cell(6,1);
    if have_ar
        A = load(ar_file);
        ar = squeeze(A.Draws_IRFs_narrative(:, 1, 1:H, :));   % (6,61,N)
        for i = 1:6
            ar_pct{i} = 100 * pctile7(squeeze(ar(i,:,:)), bands);        % (61,3) Type-7
        end
    else
        fprintf('[SKIP] %s not found -> Figure A drawn HARS-only (no blue baseline).\n', ar_file);
    end

    % ---- CMT baseline (green, unweighted, valid-draw filter) ------------
    cmt_file = fullfile(baselines_dir, 'irfs_carriero_monetary_3regimes.mat');
    have_cmt = exist(cmt_file, 'file') == 2;
    cmt_pct = cell(6,1);
    if have_cmt
        C = load(cmt_file);
        RESP = C.RESPONSE;                                    % (ndraw,6,61)
        valid = sum(abs(reshape(RESP, size(RESP,1), [])), 2) > 0;
        RESP = RESP(valid, :, :);
        fprintf('[make_figures] CMT valid draws = %d\n', size(RESP,1));
        for i = 1:6
            Xi = squeeze(RESP(:, i, :))';                     % (61, Nvalid)
            cmt_pct{i} = 100 * pctile7(Xi, bands);            % (61,3) Type-7
        end
    else
        fprintf('[SKIP] %s not found -> Figures B/C drawn HARS-only (no green CMT).\n', cmt_file);
    end

    % ================= Figure A -- paper Figure 1 =========================
    figA = figure('Name','Figure 1: baseline vs HARS (3-reg)', ...
                  'Units','inches','Position',[1 1 13 10]);
    hHomo = []; hHetero = [];
    for i = 1:6
        ax = subplot(3, 2, i); hold(ax, 'on');
        if have_ar
            bp = ar_pct{i};
            plot(ax, x, bp(:,1), '--', 'Color',[0 0 1 0.7], 'LineWidth',2);
            plot(ax, x, bp(:,3), '--', 'Color',[0 0 1 0.7], 'LineWidth',2);
            hHomo = plot(ax, x, bp(:,2), '-', 'Color','b', 'LineWidth',2);
        end
        ep = ess_pct{i};
        fill(ax, [x fliplr(x)], [ep(:,1)' fliplr(ep(:,3)')], 'r', ...
             'FaceAlpha',0.25, 'EdgeColor','none');
        hHetero = plot(ax, x, ep(:,2), '-', 'Color','r', 'LineWidth',2.5);
        yline(ax, 0, 'k--', 'LineWidth',1);
        % joint ylim of both series' bands, +/-10% margin
        lo = min(ep(:,1)); hi = max(ep(:,3));
        if have_ar, lo = min(lo, min(bp(:,1))); hi = max(hi, max(bp(:,3))); end
        m = (hi - lo) * 0.1; ylim(ax, [lo-m hi+m]);
        xlim(ax, [0 59]);
        set(ax, 'FontSize', 10);
        title(ax, varnames{i}, 'FontWeight','bold', 'FontSize',13);
        hold(ax, 'off');
    end
    if have_ar
        add_top_legend([hHomo hHetero], {'Homoskedastic model','Heteroskedastic model'});
    else
        add_top_legend(hHetero, {'Heteroskedastic model'});
    end
    exportgraphics(figA, fullfile(fig_dir, 'figure_1.pdf'), ...
                   'ContentType','vector');

    % ================= Figure B -- paper Figure S1 ========================
    figB = figure('Name','Figure S1: HARS vs CMT (3-reg)', ...
                  'Units','inches','Position',[1 1 13 10]);
    hH = []; hC = [];
    for i = 1:6
        ax = subplot(3, 2, i); hold(ax, 'on');
        ep = ess_pct{i};
        fill(ax, [x fliplr(x)], [ep(:,1)' fliplr(ep(:,3)')], 'r', ...
             'FaceAlpha',0.25, 'EdgeColor','none');
        hH = plot(ax, x, ep(:,2), '-', 'Color','r', 'LineWidth',2.5);
        lo = min(ep(:,1)); hi = max(ep(:,3));
        if have_cmt
            cp = cmt_pct{i};
            fill(ax, [x fliplr(x)], [cp(:,1)' fliplr(cp(:,3)')], 'g', ...
                 'FaceAlpha',0.25, 'EdgeColor','none');
            hC = plot(ax, x, cp(:,2), '-', 'Color','g', 'LineWidth',2.5);
            lo = min(lo, min(cp(:,1))); hi = max(hi, max(cp(:,3)));
        end
        yline(ax, 0, 'k--', 'LineWidth',1);
        m = (hi - lo) * 0.1; ylim(ax, [lo-m hi+m]);
        xlim(ax, [0 59]);
        set(ax, 'FontSize', 10);
        title(ax, varnames{i}, 'FontWeight','bold', 'FontSize',13);
        hold(ax, 'off');
    end
    if have_cmt
        add_top_legend([hH hC], {'HARS','CMT'});
    else
        add_top_legend(hH, {'HARS'});
    end
    exportgraphics(figB, fullfile(fig_dir, 'figure_S1.pdf'), ...
                   'ContentType','vector');

    fprintf('[make_figures] wrote 2 PDFs to %s\n', fig_dir);
end


function add_top_legend(handles, labels)
% One top-center horizontal legend, Box off.
    lg = legend(handles, labels, 'Orientation','horizontal', 'Box','off', ...
                'FontSize',12);
    lg.Units = 'normalized';
    p = lg.Position;
    lg.Position = [0.5 - p(3)/2, 0.965, p(3), p(4)];
end


function make_efficiency_tables(R, R0, baselines_dir, out_dir)
% mESS (Vats-Flegal-Jones 2019) + HARS repair diagnostics, four rows:
%   1. HARS (N_theta = 100)   from R    time = R.t_wall     [prof]
%   2. HARS (N_theta = 0)      from R0   time = R0.t_wall    [prof]
%   3. CMT (Carriero, 3-reg)   baseline  time = 72435.88 s   [man]  (valid-draw filter)
%   4. AR 2018 homoskedastic   baseline  time =  5663.59 s   [man]  (informational)
    vars_use = 1:6; horizons_use = 1:12; shock_use = 1;

    rows = struct('name',{},'kind',{},'found',{},'N',{},'mess',{},'mess_per_N',{}, ...
                  'sec_eff',{},'hr1000',{},'tsrc',{}, ...
                  'stale',{},'p_post',{},'p_rot',{},'disc',{});

    % ---- Row 1/2: HARS ---------------------------------------------------
    rows(end+1) = hars_row('HARS (N_theta = 100)', R,  vars_use, horizons_use, shock_use);
    rows(end+1) = hars_row('HARS (N_theta = 0)',   R0, vars_use, horizons_use, shock_use);

    % ---- Row 3: CMT (Carriero, 3-reg, 500 draws), valid-draw filter ------
    cmt_file = fullfile(baselines_dir, 'irfs_carriero_monetary_3regimes.mat');
    if exist(cmt_file, 'file') == 2
        C = load(cmt_file);
        RESP = C.RESPONSE;                              % (ndraw,6,61)
        valid = sum(abs(reshape(RESP, size(RESP,1), [])), 2) > 0;
        RESP = RESP(valid, :, :);
        irfs = permute(RESP, [2 3 1]);                  % (6,61,Nvalid)
        irfs = reshape(irfs, size(irfs,1), size(irfs,2), 1, size(irfs,3));
        rows(end+1) = baseline_row('CMT (Carriero, 3-reg, 500 draws)', irfs, 72435.88, ...
                                    vars_use, horizons_use, shock_use);
    else
        fprintf('[SKIP] %s not found -> CMT table row dropped.\n', cmt_file);
    end

    % ---- Row 4: AR 2018 homoskedastic (informational) --------------------
    ar_file = fullfile(baselines_dir, 'irfs_ramirez_monetary.mat');
    if exist(ar_file, 'file') == 2
        A = load(ar_file);
        irfs = permute(A.Draws_IRFs_narrative, [1 3 2 4]);   % (ny,K,H,N) -> (ny,H,K,N)
        rows(end+1) = baseline_row('AR 2018 homoskedastic [informational]', irfs, 5663.59, ...
                                   vars_use, horizons_use, shock_use);
    else
        fprintf('[SKIP] %s not found -> AR 2018 table row dropped.\n', ar_file);
    end

    NOTE = ['Note: CMT and AR 2018 timings are measured on the authors'' hardware ', ...
            '(fixed seconds); HARS rows are timed on this machine. The CMT ', ...
            'comparison is order-of-magnitude (paper Table S1 note).'];

    % ---- Console + .txt --------------------------------------------------
    txt_file = fullfile(out_dir, 'efficiency_monetary.txt');
    fid = fopen(txt_file, 'w');
    pr = @(varargin) print_both(fid, varargin{:});

    pr('\nEfficiency / diagnostics -- monetary application (HARS)\n');
    pr('Vats-Flegal-Jones (2019) mESS; subset vars=1..6, horizons=1..12, shock=1 (p=72).\n');
    pr('Runtime source: [prof]=stored sampling-loop time (burn-in excluded); [man]=authors'' hardware.\n');

    pr('\n=== (A) Efficiency: mESS and time per effective draw ===\n');
    pr('%-40s %8s %10s %10s %12s %16s %7s\n', ...
       'Sampler','N','mESS','mESS/N','sec/eff','hr per 1000 eff','tsrc');
    pr('%s\n', repmat('-', 1, 106));
    for k = 1:numel(rows)
        r = rows(k);
        pr('%-40s %8d %10.1f %10.4f %12.3f %16.3f %7s\n', ...
           r.name, r.N, r.mess, r.mess_per_N, r.sec_eff, r.hr1000, r.tsrc);
    end

    pr('\n=== (B) Stale-Q repair diagnostics ((A0,Lambda) channel; HARS only) ===\n');
    pr('%-40s %10s %18s %16s %12s\n', ...
       'Sampler','stale-Q','by posterior(%)','by rotation(%)','discarded');
    pr('%s\n', repmat('-', 1, 100));
    for k = 1:numel(rows)
        r = rows(k);
        if strcmp(r.kind, 'hars')
            pr('%-40s %10d %18.1f %16.1f %12d\n', ...
               r.name, r.stale, r.p_post, r.p_rot, r.disc);
        else
            pr('%-40s %10s %18s %16s %12s\n', r.name, '--','--','--','--');
        end
    end
    pr(['  stale-Q = MCMC-accepted structural moves that stranded Q; ', ...
        'by posterior = fixed by (A0,Lambda) redraw; by rotation = fixed by global Q-repair; ', ...
        'discarded = unresolved stale-Q events.\n']);
    pr('\n%s\n', NOTE);
    fclose(fid);
    % echo the file we just built to the console too
    type_file(txt_file);

    % ---- CSV (flat) ------------------------------------------------------
    csv_file = fullfile(out_dir, 'efficiency_monetary.csv');
    fc = fopen(csv_file, 'w');
    fprintf(fc, ['sampler,N,mESS,mESS_per_N,sec_per_eff_draw,hours_per_1000_eff,', ...
                 'stale_q_events,resolved_posterior_pct,resolved_rotation_pct,discarded_transitions,time_source\n']);
    for k = 1:numel(rows)
        r = rows(k);
        if strcmp(r.kind, 'hars')
            fprintf(fc, '%s,%d,%.4f,%.6f,%.4f,%.6f,%d,%.1f,%.1f,%d,%s\n', ...
                cq(r.name), r.N, r.mess, r.mess_per_N, r.sec_eff, r.hr1000, ...
                r.stale, r.p_post, r.p_rot, r.disc, r.tsrc);
        else
            fprintf(fc, '%s,%d,%.4f,%.6f,%.4f,%.6f,,,,,%s\n', ...
                cq(r.name), r.N, r.mess, r.mess_per_N, r.sec_eff, r.hr1000, r.tsrc);
        end
    end
    fclose(fc);
    fprintf('[make_efficiency_tables] wrote %s\n', csv_file);

    % ---- LaTeX (booktabs): Table S1 + Table S2 ---------------------------
    write_tex(fullfile(out_dir, 'table_S1_S2_monetary.tex'), rows, NOTE);
end


function r = hars_row(name, S, vars_use, horizons_use, shock_use)
    [m, N] = mess_of_irfs(S.sign_irfs_temp, vars_use, horizons_use, shock_use);
    t = S.t_wall;
    dg = S.ess_diagnostics;
    stale  = nzf(dg, 'stale_q_events');
    byPost = nzf(dg, 'stale_resolved_structural');
    byRot  = nzf(dg, 'stale_resolved_fallback');
    disc   = nzf(dg, 'stale_unresolved');
    r = struct('name',name,'kind','hars','found',true,'N',N,'mess',m, ...
               'mess_per_N',m/N,'sec_eff',t/m,'hr1000',(t/m)*1000/3600,'tsrc','[prof]', ...
               'stale',stale, ...
               'p_post', pctv(byPost,stale), 'p_rot', pctv(byRot,stale), 'disc',disc);
end


function r = baseline_row(name, irfs, manual_sec, vars_use, horizons_use, shock_use)
    [m, N] = mess_of_irfs(irfs, vars_use, horizons_use, shock_use);
    r = struct('name',name,'kind','baseline','found',true,'N',N,'mess',m, ...
               'mess_per_N',m/N,'sec_eff',manual_sec/m,'hr1000',(manual_sec/m)*1000/3600, ...
               'tsrc','[man]','stale',NaN,'p_post',NaN,'p_rot',NaN,'disc',NaN);
end


function [m, N] = mess_of_irfs(irfs, vars, horizons, shock)
% irfs: (ny, H, K_shocks, N). Reshape the requested subset to (N x p), mESS.
    [ny, Ht, Kt, N] = size(irfs);
    vars = vars(vars <= ny); horizons = horizons(horizons <= Ht); shock = shock(shock <= Kt);
    Y = reshape(irfs(vars, horizons, shock, :), [], N)';   % N x p
    m = mess_vfj(Y);
end

function P = pctile7(X, q)
% numpy-compatible (Type-7, 'linear') unweighted percentiles along the draw
% axis (columns). X: (H, N); q: quantiles in [0,1]. Returns (H, numel(q)).
    [H, N] = size(X); q = q(:)';
    pp = (0:N-1) / (N - 1);
    P = zeros(H, numel(q));
    for h = 1:H
        xs = sort(X(h, :));
        P(h, :) = interp1(pp, xs, q, 'linear');
    end
end


function v = nzf(d, fld)
% Safe scalar fetch: field value, or 0 if missing/empty/NaN.
    if isfield(d, fld) && ~isempty(d.(fld)) && isfinite(d.(fld)(1))
        v = double(d.(fld)(1));
    else
        v = 0;
    end
end


function p = pctv(num, den)
    if den > 0, p = 100 * num / den; else, p = NaN; end
end


function print_both(fid, varargin)
% fprintf to both the console and the open file handle fid.
    fprintf(varargin{:});
    fprintf(fid, varargin{:});
end


function type_file(fname)
    fid = fopen(fname, 'r');
    if fid < 0, return; end
    while true
        ln = fgetl(fid);
        if ~ischar(ln), break; end
        fprintf('%s\n', ln);
    end
    fclose(fid);
end


function s = cq(s)
% Quote a CSV field if it contains a comma.
    if contains(s, ','), s = ['"' s '"']; end
end


function write_tex(fname, rows, NOTE)
% Two booktabs tables:
%   S1 = CMT vs HARS(100) hours per 1,000 effective draws.
%   S2 = N_theta 100 vs 0: Hrs/1,000, Discarded, Stale-Q, by posterior-draw %, by rotation %.
% 3 decimals for hours, 1 decimal for %.
    getrow = @(pred) find(arrayfun(pred, rows), 1);
    iH100 = getrow(@(r) strcmp(r.name, 'HARS (N_theta = 100)'));
    iH0   = getrow(@(r) strcmp(r.name, 'HARS (N_theta = 0)'));
    iCMT  = getrow(@(r) startsWith(r.name, 'CMT'));

    fid = fopen(fname, 'w');
    fprintf(fid, '%% Auto-generated by main.m (make_efficiency_tables).\n');
    fprintf(fid, '%% Requires \\usepackage{booktabs}.\n\n');

    % ---- Table S1 -------------------------------------------------------
    fprintf(fid, '\\begin{table}[htbp]\\centering\n');
    fprintf(fid, '\\caption{Sampling efficiency: hours per 1{,}000 effective draws (monetary).}\n');
    fprintf(fid, '\\label{tab:S1_monetary}\n');
    fprintf(fid, '\\begin{tabular}{lrr}\n\\toprule\n');
    fprintf(fid, 'Sampler & $N$ & Hours / 1{,}000 eff. \\\\\n\\midrule\n');
    if ~isempty(iH100)
        r = rows(iH100);
        fprintf(fid, 'HARS ($N_\\theta=100$) & %d & %.3f \\\\\n', r.N, r.hr1000);
    end
    if ~isempty(iCMT)
        r = rows(iCMT);
        fprintf(fid, 'CMT (Carriero, 3-reg) & %d & %.3f \\\\\n', r.N, r.hr1000);
    end
    fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
    fprintf(fid, '\\par\\footnotesize %s\n', tex_escape(NOTE));
    fprintf(fid, '\\end{table}\n\n');

    % ---- Table S2 -------------------------------------------------------
    fprintf(fid, '\\begin{table}[htbp]\\centering\n');
    fprintf(fid, '\\caption{Structural-redraw ablation ($N_\\theta$): repair channel (monetary).}\n');
    fprintf(fid, '\\label{tab:S2_monetary}\n');
    fprintf(fid, '\\begin{tabular}{lrrrrr}\n\\toprule\n');
    fprintf(fid, ['$N_\\theta$ & Hrs / 1{,}000 & Discarded & Stale-Q & ', ...
                  'by posterior (\\%%) & by rotation (\\%%) \\\\\n\\midrule\n']);
    for ii = [iH100 iH0]
        if isempty(ii), continue; end
        r = rows(ii);
        ntheta = 100; if contains(r.name, '= 0'), ntheta = 0; end
        fprintf(fid, '%d & %.3f & %d & %d & %.1f & %.1f \\\\\n', ...
            ntheta, r.hr1000, r.disc, r.stale, r.p_post, r.p_rot);
    end
    fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
    fprintf(fid, '\\end{table}\n');
    fclose(fid);
    fprintf('[make_efficiency_tables] wrote %s\n', fname);
end


function s = tex_escape(s)
    s = strrep(s, '%', '\%');
    s = strrep(s, '&', '\&');
    s = strrep(s, '_', '\_');
end
