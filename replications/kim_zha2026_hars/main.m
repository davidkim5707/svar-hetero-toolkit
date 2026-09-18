%==========================================================================
% main.m -- master one-click replication for all three applications (HARS).
%
% Total pipeline = 7 HARS runs:
%   monetary : baseline (N_theta=100) + ablation (N_theta=0)
%   oil      : baseline (N_theta=100) + ablation (N_theta=0)
%   fiscal   : homoskedastic + 3-regime baseline (N_theta=100) + ablation (N_theta=0)
% Each main_<app>.m also runs standalone. All runs execute fresh every time;
% expect roughly 7x a single-run wall clock.
%
% After the three apps, this script builds the cross-application outputs:
%   Figure 5   -> output/figures/figure_5.pdf  (monetary + oil selected panels)
%   Table 7    -> output/table_7.tex           (repair diagnostics, all apps)
%   Table S7   -> output/table_S7.tex          (CMT vs HARS hours/1,000 + speed-ups)
%==========================================================================
run('main_monetary.m');
run('main_oil.m');
run('main_fiscal.m');

clear; close all;
addpath(genpath(fullfile('..','..','core')), genpath(fullfile('..','..','third_party')));

Rm  = load(fullfile('output','results_monetary.mat'));
Rm0 = load(fullfile('output','results_monetary_nostrucredraw.mat'));
Ro  = load(fullfile('output','results_oil.mat'));
Ro0 = load(fullfile('output','results_oil_nostrucredraw.mat'));
Rf  = load(fullfile('output','results_fiscal.mat'));
Rf0 = load(fullfile('output','results_fiscal_nostrucredraw.mat'));

make_figure5(Rm, Ro, 'baselines', fullfile('output','figures'));
make_table7(Rm, Rm0, Ro, Ro0, Rf, Rf0, 'output');
make_tableS7(Rm, Ro, Rf, 'baselines', 'output');


%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function make_figure5(Rm, Ro, baselines_dir, fig_dir)
% Combined paper Figure 5 (notebook cell 5): 1x2.
%   left  = monetary GDP deflator (var 2), HARS vs CMT, x100, H=61
%   right = oil economic activity (var 2) to oil-specific demand (shock 3),
%           HARS vs CMT, no x100, H=18
    if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
    bands = [0.05 0.50 0.95];

    % ---- left: monetary GDP deflator ----
    em = squeeze(Rm.sign_irfs_temp);                 % (6,61,N)
    wm = Rm.draw_weight(:) / sum(Rm.draw_weight(:));
    ess_mon = 100 * wpercentile(squeeze(em(2,:,:)), wm, bands);   % (61,3)
    cmt_mon = [];
    cmf = fullfile(baselines_dir,'irfs_carriero_monetary_3regimes.mat');
    if exist(cmf,'file')==2
        C = load(cmf); RESP = C.RESPONSE;                         % (ndraw,6,61)
        valid = sum(abs(reshape(RESP,size(RESP,1),[])),2) > 0; RESP = RESP(valid,:,:);
        cmt_mon = 100 * pctile7(squeeze(RESP(:,2,:))', bands);    % (61,3)
    end

    % ---- right: oil economic activity to oil-specific demand ----
    Hn = 18;
    eo = Ro.sign_irfs_temp(:,1:Hn,:,:);              % (3,18,3,N)
    wo = Ro.draw_weight(:) / sum(Ro.draw_weight(:));
    ess_oil = wpercentile(squeeze(eo(2,:,3,:)), wo, bands);       % (18,3), var2 shock3
    cmt_oil = [];
    cof = fullfile(baselines_dir,'irfs_carriero_oil_3regimes.mat');
    if exist(cof,'file')==2
        C = load(cof); cmt = cmt_scale_oil(C.RESPONSE, C.inv_A0_all);  % (nd,3,3,nhor)
        co = squeeze(cmt(:,2,3,:));                  % (nd, nhor)
        cmt_oil = pctile7(co(:,1:Hn)', bands);       % (18,3)
    end

    fig = figure('Name','Figure 5: selected panels','Units','inches','Position',[1 1 13 4.5]);
    % left panel
    axL = subplot(1,2,1); hold(axL,'on'); axL.Toolbar = [];
    xL = 0:60;
    fill(axL,[xL fliplr(xL)],[ess_mon(:,1)' fliplr(ess_mon(:,3)')],'r','FaceAlpha',0.25,'EdgeColor','none');
    hHm = plot(axL, xL, ess_mon(:,2), '-', 'Color','r', 'LineWidth',2.5);
    hCm = [];
    if ~isempty(cmt_mon)
        fill(axL,[xL fliplr(xL)],[cmt_mon(:,1)' fliplr(cmt_mon(:,3)')],'g','FaceAlpha',0.25,'EdgeColor','none');
        hCm = plot(axL, xL, cmt_mon(:,2), '-', 'Color','g', 'LineWidth',2.5);
    end
    yline(axL,0,'k--','LineWidth',1); xlim(axL,[0 59]); set(axL,'FontSize',10);
    title(axL,'GDP deflator to a contractionary monetary policy shock','FontWeight','bold','FontSize',13);
    hold(axL,'off');
    % right panel
    axR = subplot(1,2,2); hold(axR,'on'); axR.Toolbar = [];
    xR = 0:(Hn-1);
    fill(axR,[xR fliplr(xR)],[ess_oil(:,1)' fliplr(ess_oil(:,3)')],'r','FaceAlpha',0.25,'EdgeColor','none');
    plot(axR, xR, ess_oil(:,2), '-', 'Color','r', 'LineWidth',2.5);
    if ~isempty(cmt_oil)
        fill(axR,[xR fliplr(xR)],[cmt_oil(:,1)' fliplr(cmt_oil(:,3)')],'g','FaceAlpha',0.25,'EdgeColor','none');
        plot(axR, xR, cmt_oil(:,2), '-', 'Color','g', 'LineWidth',2.5);
    end
    yline(axR,0,'k--','LineWidth',1); xlim(axR,[0 Hn-1]); set(axR,'FontSize',10);
    title(axR,'Economic activity index to an oil-specific demand shock','FontWeight','bold','FontSize',13);
    hold(axR,'off');

    set(axL, 'Position', [0.07 0.14 0.40 0.72]);   % leave headroom for the legend
    set(axR, 'Position', [0.57 0.14 0.40 0.72]);
    if ~isempty(hCm)
        lg = legend([hHm hCm], {'HARS','CMT'}, 'Orientation','horizontal','Box','off','FontSize',12);
    else
        lg = legend(hHm, {'HARS'}, 'Box','off','FontSize',12);
    end
    lg.Units='normalized'; p=lg.Position; lg.Position=[0.5-p(3)/2, 0.94, p(3), p(4)];
    exportgraphics(fig, fullfile(fig_dir,'figure_5.pdf'), 'ContentType','vector');
    fprintf('[make_figure5] wrote figure_5.pdf\n');
end


function out = cmt_scale_oil(RESPONSE, inv_A0_all)
    nd = size(RESPONSE,1); c = zeros(nd,3);
    c(:,1) = -inv_A0_all(:,1,1); c(:,2) = inv_A0_all(:,2,2); c(:,3) = inv_A0_all(:,3,3);
    out = RESPONSE .* reshape(c,[nd,1,3,1]);
end


function make_table7(Rm, Rm0, Ro, Ro0, Rf, Rf0, out_dir)
% Table 7 (main text): structural-redraw repair diagnostics across the three
% applications, rows N_theta=100/0 per app.
    apps = { 'Monetary', Rm, Rm0; 'Oil', Ro, Ro0; 'Fiscal', Rf, Rf0 };
    rows = struct('app',{},'ntheta',{},'hr1000',{},'disc',{},'stale',{},'p_post',{},'p_rot',{});
    subsets = struct('Monetary',{{1:6,1:12,1}}, 'Oil',{{1:3,1:12,1:3}}, 'Fiscal',{{1:5,1:12,2}});
    for a = 1:3
        nm = apps{a,1}; sub = subsets.(nm);
        for col = 2:3
            S = apps{a,col}; nt = 100*(col==2);
            [m,~] = mess_of_irfs(S.sign_irfs_temp, sub{1}, sub{2}, sub{3});
            hr = (S.t_wall/m)*1000/3600; dg = S.ess_diagnostics;
            st = nzf(dg,'stale_q_events');
            rows(end+1) = struct('app',nm,'ntheta',nt,'hr1000',hr, ...
                'disc',nzf(dg,'stale_unresolved'),'stale',st, ...
                'p_post',pctv(nzf(dg,'stale_resolved_structural'),st), ...
                'p_rot', pctv(nzf(dg,'stale_resolved_fallback'),st)); %#ok<AGROW>
        end
    end
    % console
    fprintf('\n=== Table 7: structural-redraw repair diagnostics (all applications) ===\n');
    fprintf('%-10s %8s %12s %12s %10s %16s %14s\n','App','N_theta','Hrs/1000','Discarded','Stale-Q','byPosterior(%)','byRotation(%)');
    fprintf('%s\n', repmat('-',1,86));
    for k=1:numel(rows)
        r=rows(k);
        fprintf('%-10s %8d %12.3f %12d %10d %16.1f %14.1f\n', r.app, r.ntheta, r.hr1000, r.disc, r.stale, r.p_post, r.p_rot);
    end
    % tex
    fid = fopen(fullfile(out_dir,'table_7.tex'),'w');
    fprintf(fid,'%% Auto-generated by main.m. Requires \\usepackage{booktabs}.\n\n');
    fprintf(fid,'\\begin{table}[htbp]\\centering\n');
    fprintf(fid,'\\caption{Structural-redraw repair diagnostics across applications.}\n\\label{tab:7}\n');
    fprintf(fid,'\\begin{tabular}{llrrrrr}\n\\toprule\n');
    fprintf(fid,'Application & $N_\\theta$ & Hrs / 1{,}000 & Discarded & Stale-Q & by posterior (\\%%) & by rotation (\\%%) \\\\\n\\midrule\n');
    for k=1:numel(rows)
        r=rows(k);
        fprintf(fid,'%s & %d & %.3f & %d & %d & %.1f & %.1f \\\\\n', r.app, r.ntheta, r.hr1000, r.disc, r.stale, r.p_post, r.p_rot);
    end
    fprintf(fid,'\\bottomrule\n\\end{tabular}\n\\end{table}\n');
    fclose(fid);
    fprintf('[make_table7] wrote %s\n', fullfile(out_dir,'table_7.tex'));
end


function make_tableS7(Rm, Ro, Rf, baselines_dir, out_dir)
% Table S7: CMT vs HARS hours per 1,000 effective draws, with speed-ups.
% CMT seconds are fixed constants from the compute_irf_mess reference scripts
% (monetary 72435.88 s; oil 5112.00 s; fiscal CMT fails). HARS hours are
% computed from this machine's slims.
    CMT_MON_SEC = 72435.88; CMT_OIL_SEC = 5112.00;

    % HARS hours per app
    hr_mon = hars_hr(Rm, 1:6, 1:12, 1);
    hr_oil = hars_hr(Ro, 1:3, 1:12, 1:3);
    hr_fis = hars_hr(Rf, 1:5, 1:12, 2);

    % CMT hours per app (from baselines)
    cmt_mon = cmt_hr_monetary(fullfile(baselines_dir,'irfs_carriero_monetary_3regimes.mat'), CMT_MON_SEC);
    cmt_oil = cmt_hr_oil(fullfile(baselines_dir,'irfs_carriero_oil_3regimes.mat'), CMT_OIL_SEC);

    fprintf('\n=== Table S7: CMT vs HARS hours per 1,000 effective draws ===\n');
    fprintf('%-10s %14s %12s %10s\n','App','CMT Hrs/1000','HARS Hrs/1000','Speed-up');
    fprintf('%s\n', repmat('-',1,50));
    fprintf('%-10s %14.3f %12.3f %9.0fx   (order-of-magnitude; extrapolation note)\n','Monetary', cmt_mon, hr_mon, cmt_mon/hr_mon);
    fprintf('%-10s %14.3f %12.3f %9.0fx\n','Oil', cmt_oil, hr_oil, cmt_oil/hr_oil);
    fprintf('%-10s %14s %12.3f %10s\n','Fiscal','fails', hr_fis, '--');

    fid = fopen(fullfile(out_dir,'table_S7.tex'),'w');
    fprintf(fid,'%% Auto-generated by main.m. Requires \\usepackage{booktabs}.\n\n');
    fprintf(fid,'\\begin{table}[htbp]\\centering\n');
    fprintf(fid,'\\caption{CMT vs HARS sampling cost: hours per 1{,}000 effective draws.}\n\\label{tab:S7}\n');
    fprintf(fid,'\\begin{tabular}{lrrr}\n\\toprule\n');
    fprintf(fid,'Application & CMT Hrs / 1{,}000 & HARS Hrs / 1{,}000 & Speed-up \\\\\n\\midrule\n');
    fprintf(fid,'Monetary & %.3f & %.3f & $\\approx$%.0f$\\times$ \\\\\n', cmt_mon, hr_mon, cmt_mon/hr_mon);
    fprintf(fid,'Oil & %.3f & %.3f & $\\approx$%.0f$\\times$ \\\\\n', cmt_oil, hr_oil, cmt_oil/hr_oil);
    fprintf(fid,'Fiscal & \\multicolumn{1}{c}{fails} & %.3f & --- \\\\\n', hr_fis);
    fprintf(fid,'\\bottomrule\n\\end{tabular}\n');
    fprintf(fid,'\\par\\footnotesize The monetary CMT comparison is order-of-magnitude (see the Table S1 note); the CMT hours rest on an extrapolation from the reported run. The fiscal CMT sampler fails (1 admissible draw in 200{,}000 candidates).\n');
    fprintf(fid,'\\end{table}\n');
    fclose(fid);
    fprintf('[make_tableS7] wrote %s\n', fullfile(out_dir,'table_S7.tex'));
end


function hr = hars_hr(S, vars, hor, shock)
    [m,~] = mess_of_irfs(S.sign_irfs_temp, vars, hor, shock);
    hr = (S.t_wall/m)*1000/3600;
end

function hr = cmt_hr_monetary(f, sec)
    hr = NaN; if exist(f,'file')~=2, return; end
    C = load(f); RESP = C.RESPONSE; valid = sum(abs(reshape(RESP,size(RESP,1),[])),2)>0; RESP=RESP(valid,:,:);
    irfs = reshape(permute(RESP,[2 3 1]), size(RESP,2), size(RESP,3), 1, sum(valid));
    [m,~] = mess_of_irfs(irfs, 1:6, 1:12, 1); hr = (sec/m)*1000/3600;
end

function hr = cmt_hr_oil(f, sec)
    hr = NaN; if exist(f,'file')~=2, return; end
    C = load(f); irfs = permute(C.RESPONSE, [2 4 3 1]);
    [m,~] = mess_of_irfs(irfs, 1:3, 1:12, 1:3); hr = (sec/m)*1000/3600;
end


function [m, N] = mess_of_irfs(irfs, vars, horizons, shock)
    [ny, Ht, Kt, N] = size(irfs);
    vars = vars(vars <= ny); horizons = horizons(horizons <= Ht); shock = shock(shock <= Kt);
    Y = reshape(irfs(vars, horizons, shock, :), [], N)';
    m = mess_vfj(Y);
end

function P = pctile7(X, q)
    [H, N] = size(X); q = q(:)';
    pp = (0:N-1) / (N - 1); P = zeros(H, numel(q));
    for h = 1:H, xs = sort(X(h,:)); P(h,:) = interp1(pp, xs, q, 'linear'); end
end

function v = nzf(d, fld)
    if isfield(d, fld) && ~isempty(d.(fld)) && isfinite(d.(fld)(1)), v = double(d.(fld)(1)); else, v = 0; end
end

function p = pctv(num, den)
    if den > 0, p = 100 * num / den; else, p = NaN; end
end
