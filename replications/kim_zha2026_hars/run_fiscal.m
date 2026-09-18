%==========================================================================
% run_fiscal.m -- model, options, data, and sampler call for the fiscal
% application (Kim & Zha, HARS). Builds the 5-variable quarterly VAR (Caldara
% & Kamps 1950-2006 sample, 4 lags, regime breaks 1971Q2 and 1997Q1), the
% joint two-shock revenue+spending sign scheme with narrative anchors, and the
% HARS/ESS options, then times and runs the sampler (hars.m). main_fiscal.m
% sets REGIME_ON and NTHETA per run.
%==========================================================================
%********************************************************
%% Suppress numerical warnings during sampling
%********************************************************
warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:rankDeficientMatrix');
warning('off', 'MATLAB:illConditionedMatrix');
warning('off', 'MATLAB:logOfZero');
warning('off', 'MATLAB:divideByZero');

%********************************************************
%% Data loading and setting
%********************************************************

addpath(genpath(fullfile('..','..','core')), genpath(fullfile('..','..','third_party')));

% === DATASET VARIANT ===
% 'mu_reduced_ver2' -> MU-style 5-var: net taxes, G, GDP, 3-month T-bill (level), GDP deflator
data_variant = 'mu_reduced_ver2';

data_file = sprintf('data/Caldara_DATASET_extended_%s.mat', data_variant);
data_tag  = ['_', data_variant];
load(data_file);
fprintf('Dataset: %s\n', data_file);

if ~isfield(DATASET, 'MAP') || isempty(DATASET.MAP)
    names = cellstr(DATASET.NAMES);
    DATASET.MAP = containers.Map(names, num2cell(1:numel(names)));
end

% Variable selection depends on the nominal block of the chosen variant.
% mu_reduced / mu_reduced_ver2 replace the T-bill + CPI block with an interest
% rate + GDP deflator; ver2 swaps FedFunds for the 3-month T-bill (TB3MS).
% ver3 keeps ver2's NAMES but TAX_S/G_S are Arias-style ratios (flow_t / GDP_{t-1}).
if strcmp(data_variant, 'mu_reduced_ver2')
    select_vars = {'TAX_S', 'G_S', 'GDP_S', 'TB3MS', 'GDP_DEFL'};
    varnames    = {'Net_Taxes', 'Gov_Spending', 'GDP', 'Tbill_3m', 'GDP_Deflator'};
else
    error('run_fiscal:config', 'Only the paper configuration is included in this replication package.');
end

col_idx = zeros(1, length(select_vars));
for i = 1:length(select_vars)
    col_idx(i) = DATASET.MAP(select_vars{i});
end
data_raw = DATASET.TSERIES(:, col_idx);

T0 = 0;
y = data_raw(T0+1:end, :);

filtered_data = struct();
startStr   = char(DATASET.START);
start_date = datetime(str2double(startStr(1:4)), 3*(str2double(startStr(end)) - 1) + 1, 1);
filtered_data.dates = start_date + calmonths(3*(0:size(data_raw,1)-1)');
filtered_data.data = data_raw;

%********************************************************
%% Restrict to the original Caldara-Kamps sample (optional)
%********************************************************
use_original_sample = 1;

if use_original_sample
    sample_start = datetime(1947, 1, 1);
    sample_end   = datetime(2007, 12, 1);
    keep = (filtered_data.dates >= sample_start) & (filtered_data.dates <= sample_end);

    data_raw            = data_raw(keep, :);
    filtered_data.dates = filtered_data.dates(keep);
    filtered_data.data  = data_raw;
    y                   = data_raw;

    sample_tag = '_ck1950_2006';
    fprintf('Sample restricted to %s - %s (%d obs; Caldara-Kamps original)\n', ...
        datestr(sample_start, 'yyyyQQ'), datestr(sample_end, 'yyyyQQ'), size(y, 1));
else
    sample_tag = '_full';
end

fprintf('Data loaded successfully!\n');
fprintf('Sample start: %s\n', datestr(filtered_data.dates(1), 'yyyyQQ'));
fprintf('Sample end:   %s\n', datestr(filtered_data.dates(end), 'yyyyQQ'));
fprintf('Number of observations: %d\n', size(y, 1));
fprintf('Number of variables: %d\n', size(y, 2));

fprintf('\nVariable names:\n');
for i = 1:length(varnames)
    fprintf('%d. %s\n', i, varnames{i});
end

%********************************************************
%% USER SETTINGS
%********************************************************
% --------------------------------------------------------------
% MU-style fiscal identification. This package ships the joint mode.
%
% identification_mode = 'joint':
%   single rotation Q. BC [+ monetary] + revenue + spending share that Q, so
%   the two fiscal columns are mutually orthogonal. (Stronger than MU/CK.)
%
% The ASYMMETRIC sign scheme: the revenue column restricts
% Net Taxes (G free), the spending column restricts G (Net Taxes free), and the
% two carry OPPOSITE GDP signs (rev: GDP<0; spend: GDP>0). GDP sign is IMPOSED,
% so the cumulative multiplier below is partially SIGN-FORCED and is an
% IRF/identification cross-check ONLY -- the substantive multiplier comes from
% the GDP-UNRESTRICTED runs (_nogdpsign).
%
% NON-FISCAL BLOCK IS TOGGLEABLE: the business-cycle and monetary columns are
% each switched on/off via include_business_cycle / include_monetary. The two
% fiscal columns (revenue, spending) always sit AFTER the non-fiscal block, so
% rev_col/spend_col shift automatically with the block size n_nonfiscal.
%
% BUSINESS-CYCLE COLUMN: the standard MU/CK "anchor" non-fiscal shock; drop it
% (include_business_cycle = 0) only when you deliberately want a fiscal-only (or
% fiscal + monetary) identification, otherwise keep it on.
%
% MONETARY COLUMN: absorbs post-2022 interest-rate variation in net taxes; drop
% it (include_monetary = 0) only on the CK sample (1950-2006); keep it on the
% extended sample.
%
% References: ADRR (2018); Blanchard-Perotti (2002); BGSS (2021);
%   Caldara-Kamps (2017); Carriero et al. (2024); Jorgensen-Ravn (2022);
%   Mertens-Ravn (2013); Mountford-Uhlig (2009); Ramey (2011);
%   Ramey-Zubairy (2018); Romer-Romer (2010); Sims (2020).
% --------------------------------------------------------------
if ~exist('REGIME_ON','var'), REGIME_ON = 1; end
regime_on     = REGIME_ON;    % main_fiscal sets 0 for the homoskedastic run
prior_on      = 0;    % 1 = Minnesota prior
svar_svsr_on  = 1;    % 1 = KZ absorb-then-rotate
penalty_offdiagonal_on = 0;
narrative_on  = 1;    % 1 = impose narrative anchors; 0 = sign-only

% === IDENTIFICATION MODE & SHOCK STRUCTURE ===
include_business_cycle = 0;          % 1 -> add a business-cycle column to the non-fiscal block
include_monetary       = 0;          % 1 -> add a monetary column to the non-fiscal block
identification_mode    = 'joint';    % 'joint' | 'sequential'
fiscal_target          = 'spending';  % sequential only: 'revenue' | 'spending'

is_sequential = strcmp(identification_mode, 'sequential');

% --- Non-fiscal block layout (BC first, monetary next, then fiscal) ---
% bc_col / mon_col are [] when the corresponding shock is switched off; the
% fiscal columns always follow the active non-fiscal block.
n_nonfiscal = include_business_cycle + include_monetary;
next_col    = 1;
if include_business_cycle
    bc_col  = next_col;   next_col = next_col + 1;
else
    bc_col  = [];
end
if include_monetary
    mon_col = next_col;   next_col = next_col + 1;
else
    mon_col = [];
end

if strcmp(identification_mode, 'joint')
    do_rev   = true;
    do_spend = true;
    n_shocks = n_nonfiscal + 2;          % [BC] [+ MP] + revenue + spending
    rev_col   = next_col;                % first fiscal column
    spend_col = next_col + 1;            % second fiscal column
else
    error('run_fiscal:config', 'Only the paper configuration is included in this replication package.');
end

assert(n_shocks >= 1, 'No shocks selected: enable at least one fiscal target.');

% === MODEL STRUCTURE ===
options.ny            = size(y, 2);
options.lags          = 4;
options.irf_horizon   = 60;
options.constant      = 1;
options.timetrend     = 0;
options.nexogenous    = 0;
options.exogenous     = [];
options.non_explosive_ = 0;

% === SAMPLING SETTINGS ===
options.firstobs      = options.lags + 1;
options.presample     = 0;
options.ndraw         = 10000;
options.nsep          = 1;
options.nburn         = 1000;
options.max_compute   = 1;
options.nit           = 100;
options.hsnscale      = 0.05;

% === DATA DIMENSIONS ===
options.nobs          = size(y, 1) - options.firstobs + 1;
options.nunits        = 1;

% === PRIOR SETTINGS ===
if prior_on == 1
    options.dummy               = 1;
    options.minn_prior_tau      = 0.2;
    options.minn_prior_decay    = 2.5;
    options.minn_prior_lambda   = 0.5;
    options.minn_prior_mu       = 0.5;
    options.minn_prior_omega    = 0;
    options.unitroot            = ones(options.ny, 1);
else
    options.dummy               = [];
    options.minn_prior_tau      = [];
    options.minn_prior_decay    = [];
    options.minn_prior_lambda   = [];
    options.minn_prior_mu       = [];
    options.minn_prior_omega    = [];
    options.unitroot            = [];
    options.vprior_sig          = [];
end

%********************************************************
%% REGIME SETTINGS (Heteroskedasticity)
%********************************************************
% Original (1950Q1-2006Q4): 4-regime fiscal-history grid.
%   1971Q2 -- end of Bretton Woods (Carriero / BGSS)
%   1979Q4 -- fiscal-volatility decline (Perotti 2005)
%   1997Q1 -- deficit cuts begin (Carriero)
% Extended sample: 5-regime grid (1971Q2, 1997Q1, 2008Q1, 2020Q1, 2022Q1).
% --------------------------------------------------------------
if regime_on == 1
    % Convention: a break date is the LAST obs of the PREVIOUS regime (rfvar3
    % assigns obs t to the new regime iff t > breakInd). Set one quarter earlier
    % so the labelled quarter is the FIRST obs of the new regime.
    if use_original_sample
        options.regimes = [
            datetime(1971, 1, 1), ...   % Bretton Woods end 1971Q2 = first quarter of new regime
            datetime(1996, 10, 1)       % 1997Q1 = first quarter of new regime
        ];
    else
        options.regimes = [
            datetime(1971, 1, 1), ...   % 1971Q2 = first quarter of new regime
            datetime(1996, 10, 1), ...  % 1997Q1 = first quarter of new regime
            datetime(2007, 10, 1)       % 2008Q1 = first quarter of new regime
        ];
    end
    assert(all(options.regimes > filtered_data.dates(1)) && ...
           all(options.regimes < filtered_data.dates(end)), ...
           'A regime break date falls outside the estimation sample.');
else
    options.regimes = [];
end
        % options.regimes = [
        %     datetime(1971, 1, 1), ...   % 1971Q2 = first quarter of new regime
        %     datetime(1996, 10, 1)       % 1997Q1 = first quarter of new regime
        % ];
%********************************************************
%% STRUCTURAL IDENTIFICATION SETTINGS
%********************************************************
options.tvA                   = 0;
options.noLmd                 = 0;
options.A0_restriction        = true(options.ny);
options.a0_prior_setting      = 1;   % flat p(A0) ~ |det(A0)|^{-n}
options.tparam                = 5.703233415698 / 2;
options.tscale                = 1;

%********************************************************
%% SIGN RESTRICTIONS (single Q per run)
%********************************************************
% Variable order: 1 Net Taxes, 2 Gov Spending, 3 GDP, 4 rate, 5 deflator/CPI.
% Asymmetric fiscal scheme (MU 2009): revenue col restricts Net Taxes
% (G free); spending col restricts G (Net Taxes free); opposite GDP signs.
% The optional BC/MP block is appended first; in sequential mode only the
% active fiscal block follows, so the run identifies [BC] [+ MP] + a single
% fiscal shock.
% --------------------------------------------------------------
options.sign_horizon = 4;

sr_common = {};

% Revenue shock (col rev_col): Net Taxes>0, GDP<0, deflator<0 (G UNRESTRICTED)
sr_rev = {
    sprintf('y(1,1:%d,%d) > 0', options.sign_horizon, rev_col), ...
    sprintf('y(3,1:%d,%d) < 0', options.sign_horizon, rev_col), ...
    % sprintf('y(5,1:%d,%d) < 0', options.sign_horizon, rev_col), ...
};

% Spending shock (col spend_col): G>0, GDP>0 (Net Taxes UNRESTRICTED)
sr_spend = {
    sprintf('y(2,1:%d,%d) > 0', options.sign_horizon, spend_col), ...
    sprintf('y(3,1:%d,%d) > 0', options.sign_horizon, spend_col), ...
};

options.SignRestrictions = sr_common;
if do_rev
    options.SignRestrictions = [options.SignRestrictions, sr_rev];
end
if do_spend
    options.SignRestrictions = [options.SignRestrictions, sr_spend];
end

% Fast impact-screen (must match impact-period signs EXACTLY)
options.impact_signs       = zeros(options.ny, options.ny);
if do_rev
    options.impact_signs(1, rev_col)   = +1;   % rev: Net Taxes>0
    options.impact_signs(3, rev_col)   = -1;   % rev: GDP<0
end
if do_spend
    options.impact_signs(2, spend_col) = +1;   % spend: G>0
    options.impact_signs(3, spend_col) = +1;   % spend: GDP>0
end
impact_signs_byregime       = [];     % [] -> no per-regime impact gate

options.sign_regime_dependent  = (svar_svsr_on == 1);
options.penalty_offdiagonal_on = (penalty_offdiagonal_on == 1);
options.store_all_valid_Q      = 1;

% record the identification layout for the sampler / post-processing
options.identification_mode    = identification_mode;
options.include_business_cycle = include_business_cycle;
options.include_monetary       = include_monetary;
options.bc_col                 = bc_col;
options.mon_col                = mon_col;
options.fiscal_cols_active     = [rev_col(do_rev), spend_col(do_spend)];

%********************************************************
%% ESS-SPECIFIC SETTINGS
%********************************************************
options.sign_inneriteration     = 2000;
options.outeriteration          = 2;
options.ess_premix_rounds       = 3;
options.ess_repair_budget       = 2000;
options.upstream_stability_check = false;
options.enforce_stability        = false;
options.use_local_pert_fallback  = true;
options.seed_rescue_rounds = 25;

% --- HARS sampler parameters (Algorithm, Sec 3) ---
% N_ess: local ESS rotations RETAINED per structural draw (Sec 3 step 3).
% >1 amortises the structural move over several Q-draws; 1 = one per sweep.
options.ess_store_rounds    = 5;

% N_theta: max JOINT (A0,Lambda,A+) structural redraws to restore the current Q
% after an MCMC-accepted completed draw strands it (Sec 3 step 1 / Prop 1).
% 0 disables -> every stale event falls to the global Q-repair (avoid). Keep > 0.
if ~exist('NTHETA','var'), NTHETA = 100; end
options.theta_redraw_budget = NTHETA;

%********************************************************
%% NARRATIVE RESTRICTIONS (active fiscal column only)
%********************************************************
% Two distinct signs (ADRR): the IRF normalization (own-variable > 0) is fixed
% by the sign block; the realized-shock sign below is the sign of the DRAWN
% shock on the anchor date (tax cut = negative, spending up = positive).
%
% rev_anchor_set   : 'kj_obra1993'  (KJ 1964 cut -1, OBRA 1993 increase +1)
% spend_anchor_set : 'all_four'     (Korea + Vietnam + Carter-Reagan + 9/11)
% --------------------------------------------------------------
rev_anchor_set   = 'kj_obra1993';
spend_anchor_set = 'all_four';

if narrative_on == 1
    target_dates    = [];
    dominance_types = {};
    shock_indices   = [];
    sign_indices    = [];
    variable_idx    = [];

    % --- revenue anchors (col rev_col, var 1 Net Taxes) ---
    % Each anchor carries its OWN realized-shock sign (ADRR Type A):
    %   tax CUT      -> -1 (net taxes fell)   : kj 1964, egtrra 2001
    %   tax INCREASE -> +1 (net taxes rose)   : obra1993
    % All three are R&R EXOGENOUS changes (kj/egtrra = long-run growth cuts;
    % obra1993 = deficit-driven increase). NB: egtrra's 2001Q3 quarter is the
    % rebate, which R&R apportion as partly endogenous -- see notes.
    if do_rev
        switch rev_anchor_set
            case 'kj_obra1993'
                rev_dates = [datetime(1964, 4, 1), datetime(1993, 10, 1)];
                rev_signs = [-1, +1];
            otherwise
                error('run_fiscal:config', 'Only the paper configuration is included in this replication package.');
        end
        n_rev = numel(rev_dates);
        rev_domin = repmat({'dominance'}, 1, n_rev);
        target_dates    = [target_dates, rev_dates];
        dominance_types = [dominance_types, rev_domin];
        shock_indices   = [shock_indices, rev_col * ones(1, n_rev)];
        sign_indices    = [sign_indices, rev_signs];                 % per-anchor: cut -1 / increase +1
        variable_idx    = [variable_idx, 1 * ones(1, n_rev)];        % Net Taxes
    end

    % --- spending anchors (col spend_col, var 2 Gov Spending, positive draws) ---
    % Ramey-Shapiro military news dates (1950Q3, 1965Q1, 1980Q1, 2001Q3),
    % each shifted to the FIRST quarter in which observed outlays begin to
    % respond. Outlay-response anchors:
    %   Korea 1951Q1 | Vietnam 1965Q3 | Carter-Reagan 1980Q1 | 9/11 2001Q4
    if do_spend
        switch spend_anchor_set
            case 'all_four'   % Korea + Vietnam + Carter-Reagan + 9/11 (CK sample only)
                spend_dates = [datetime(1951, 1, 1), datetime(1965, 7, 1), ...
                               datetime(1980, 1, 1), datetime(2001, 10, 1)];
                spend_domin = {'dominance', 'dominance', 'dominance', 'dominance'};
            otherwise
                error('run_fiscal:config', 'Only the paper configuration is included in this replication package.');
        end
        n_spend = numel(spend_dates);
        target_dates    = [target_dates, spend_dates];
        dominance_types = [dominance_types, spend_domin];
        shock_indices   = [shock_indices, spend_col * ones(1, n_spend)];
        sign_indices    = [sign_indices, +1 * ones(1, n_spend)];      % spending up positive
        variable_idx    = [variable_idx, 2 * ones(1, n_spend)];       % Gov Spending
    end

    num_narr = numel(target_dates);
    assert(numel(dominance_types) == num_narr, ...
        'dominance_types must have one entry per narrative anchor.');
    assert(num_narr > 0, 'No narrative anchors selected for the active fiscal target.');

    % Guard: ARP (2021Q1) requires the extended sample.
    if use_original_sample && any(target_dates == datetime(2021, 1, 1))
        error('The ARP anchor (2021Q1) requires the extended sample: set use_original_sample = 0.');
    end

    full_sample_idx = arrayfun(@(d) find(filtered_data.dates == d), target_dates);
    if any(arrayfun(@(c) isempty(c{1}), num2cell(full_sample_idx)))
        error('At least one narrative anchor date is out of sample. Check dataset range.');
    end

    effective_idx = full_sample_idx - options.lags;
    if any(effective_idx < 1)
        error('Narrative anchor falls within the first %d lags (effective_idx = %d).', ...
              options.lags, min(effective_idx));
    end

    options.NarrativeRestrictions.time_index     = effective_idx;
    options.NarrativeRestrictions.shock_index    = shock_indices;
    options.NarrativeRestrictions.sign           = sign_indices;
    options.NarrativeRestrictions.variable_index = variable_idx;
    options.NarrativeRestrictions.type           = dominance_types;
else
    options.NarrativeRestrictions = [];
end

%********************************************************
%% NRR OMEGA IMPORTANCE WEIGHTS (post-hoc, storage stage)
%********************************************************
options.use_nrr_omega = false;    % AD-RR 1/omega_hat importance weighting (off)

% === DATA STRUCT ===
data.y              = y;
data.varnames       = varnames;
data.filtered_data  = filtered_data;

%********************************************************
%% Run banner
%********************************************************
% build shock label list for the banner
shock_lbls = {};
if do_rev,   shock_lbls{end+1} = 'net-tax';  end
if do_spend, shock_lbls{end+1} = 'spending'; end
shock_line = '';
for s = 1:numel(shock_lbls)
    shock_line = [shock_line, sprintf(' %d %s |', s, shock_lbls{s})];
end
shock_line = shock_line(1:end-1);   % drop trailing '|'

fprintf('\n========================================\n');
fprintf('  Caldara-Kamps %s %d-shock + KZ SVAR-H + Narrative\n', upper(identification_mode), n_shocks);
fprintf('  Shocks:%s\n', shock_line);
if is_sequential
    fprintf('  Sequential: identifying the %s shock ONLY (run the other target separately)\n', fiscal_target);
end
fprintf('  Sampler: hybrid Gibbs-MH-ESS\n');
fprintf('========================================\n');
fprintf('Dataset variant: %s\n', ternary_local(isempty(data_variant), 'gross (log units)', data_variant));
fprintf('Variables: %s\n', strjoin(varnames, ', '));
fprintf('Sample: %s - %s (%d obs)%s\n', ...
    datestr(filtered_data.dates(1), 'yyyyQQ'), ...
    datestr(filtered_data.dates(end), 'yyyyQQ'), size(y, 1), ...
    ternary_local(use_original_sample, '  [Caldara-Kamps original]', ''));
fprintf('Non-fiscal block: BC=%d, monetary=%d (n_nonfiscal=%d)\n', ...
    include_business_cycle, include_monetary, n_nonfiscal);
fprintf('Lags: %d | Regimes: %d\n', options.lags, length(options.regimes) + 1);
fprintf('A0 prior: %d (1=flat |det(A0)|^-n)\n', options.a0_prior_setting);
fprintf('Sign restriction horizon: %d quarters (GDP IMPOSED on active fiscal col)\n', options.sign_horizon);
fprintf('Sign restrictions: %d (across %d identified shocks)\n', numel(options.SignRestrictions), n_shocks);
if narrative_on == 1
    for na = 1:num_narr
        fprintf('Narrative anchor %d/%d: date %s | sign %+d | var %d | obs %d (eff %d) | dominance %s | on shock %d\n', ...
            na, num_narr, datestr(target_dates(na), 'yyyyQQ'), ...
            sign_indices(na), variable_idx(na), full_sample_idx(na), ...
            effective_idx(na), dominance_types{na}, shock_indices(na));
    end
else
    fprintf('Narrative: disabled (sign-only)\n');
end
fprintf('Draws: %d (burn-in %d)\n', options.ndraw, options.nburn);
fprintf('========================================\n\n');

%********************************************************
%% Run sampler
%********************************************************
sampler_script = 'hars';
if exist([sampler_script '.m'], 'file') ~= 2
    error('run_fiscal:samplerMissing', ...
        'Sampler script "%s.m" not found on the path.', sampler_script);
end

clear estimate_nrr_omega_fast

run_start_time = tic;
run(sampler_script);                       % runs hars.m
run_total_sec  = toc(run_start_time);

fprintf('End-to-end run time: %.2f s (%.2f min)\n', run_total_sec, run_total_sec/60);
warning('on', 'all');

%% Local helper
function out = ternary_local(cond, a, b)
    if cond, out = a; else, out = b; end
end