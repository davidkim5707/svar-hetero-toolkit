%==========================================================================
% ex01_harsz_balance_sheet.m -- HARS-Z on an eight-variable monthly U.S. SVAR.
%
% Identifies five labeled shocks (demand, supply, monetary policy, and two
% balance-sheet shocks) with sign, zero, and narrative restrictions and four
% volatility regimes. This is the baseline specification of Kim, "How the
% Financing of Balance-Sheet Policy Shapes Its Effects".
%
% Reads   examples/data/balance_sheet_monthly.csv   2008:M11 to 2026:M6
% Writes  examples/output/ex01_harsz_balance_sheet.mat   (output, options, data)
% Seed    none. The draws depend on the state of the generator at the start.
% Runtime about 1 hour 40 minutes for NDRAW = 10000 on the author's desktop
%         with an eight-worker pool. Set NDRAW = 300 for a run of about 35 minutes that
%         exercises every step. Keep NBURN = 1000. After a short burn-in the
%         sampler can fail to find an admissible state to start sampling from.
%
% Variables (VAR order)
%   1 one-year Treasury yield        5 log real GDP
%   2 log SOMA holdings              6 log GDP deflator
%   3 log(reserves + ON RRP)         7 term spread
%   4 ON RRP / (reserves + ON RRP)   8 excess bond premium
%
% Shocks (columns of the impact matrix)
%   1 demand   2 supply   3 monetary policy
%   4 RFBS, a balance-sheet expansion financed by reserves
%   5 OFBS, a balance-sheet expansion financed by ON RRP
%   6 to 8 unlabeled
%==========================================================================
clear; close all;

NDRAW = 10000;   % recorded draws
NBURN = 1000;    % burn-in draws

example_dir = fileparts(mfilename('fullpath'));
run(fullfile(example_dir, '..', 'setup.m'));
data_file = fullfile(example_dir, 'data', 'balance_sheet_monthly.csv');
save_file = fullfile(example_dir, 'output', 'ex01_harsz_balance_sheet.mat');

%% Data
csv_columns = {'1y_yield', 'log_soma_nb_mbs_long', 'log_reserves_plus_rrp', ...
               'share_onrrp', 'log_gdp', 'log_gdpdef', 'TermSpread', 'EBP'};
varnames    = ["RATE1Y", "SOMA", "log_R_RRP", "ONRRP_share", ...
               "GDP", "GDPDEF", "TermSpread", "EBP"];

data_raw      = readtable(data_file, 'VariableNamingRule', 'preserve');
filtered_data = data_raw(:, [{'dates'}, csv_columns]);
if ~isdatetime(filtered_data.dates)
    filtered_data.dates = datetime(filtered_data.dates, 'InputFormat', 'yyyy-MM-dd');
end
y  = table2array(filtered_data(:, 2:end));
nv = size(y, 2);

i_rate = 1;  i_soma = 2;  i_liab = 3;  i_rrp = 4;
i_gdp  = 5;  i_def  = 6;  i_ts   = 7;

%% VAR and sampler size
options.ny          = nv;
options.lags        = 6;
options.irf_horizon = 60;
options.constant    = 1;
options.timetrend   = 0;
options.exogenous   = [];

options.firstobs    = options.lags + 1;
options.presample   = 0;
options.nobs        = size(y, 1) - options.firstobs + 1;
options.ndraw       = NDRAW;
options.nburn       = NBURN;
options.nsep        = 1;
options.max_compute = 1;      % posterior mode by csminwel
options.nit         = 100;    % csminwel iterations
options.hsnscale    = 0.05;   % scale of the MH proposal covariance

%% Structural model
options.tvA                  = 0;               % A0 common across regimes
options.noLmd                = 0;               % regime-specific shock variances
options.tvA_rows             = [];
options.homoskedastic_shocks = [];              % every shock heteroskedastic
options.A0_restriction       = true(options.ny);
options.a0_prior_setting     = 0;
options.fevd_enable          = false;

% Student-t structural errors. The degrees of freedom come from a two-step
% fit on the regime-standardized residuals of the unrestricted model.
options.tparam = 5.840451331958 / 2;
options.tscale = 1;

%% Minnesota prior
options.dummy             = 1;
options.minn_prior_tau    = 3;
options.minn_prior_decay  = 0.5;
options.minn_prior_lambda = 5;
options.minn_prior_mu     = 1;
options.minn_prior_omega  = 0;
%                           1y  SOMA R+RRP share GDP DEF TS  EBP
options.unitroot          = [1;  1;   1;    0;    1;  1;  0;  0];

%% Volatility regimes
% The sampler assigns a break month to the regime that ends there. Each date
% below is the last month of a regime, and the next regime starts one month
% later, in 2014:M12, 2020:M3, and 2022:M6.
options.regimes = [datetime(2014, 11, 1), ...   % last month of the QE3 purchases
                   datetime(2020,  2, 1), ...   % last month before the COVID purchases
                   datetime(2022,  5, 1)];      % last month before the QT2 runoff

%% Sign restrictions, horizons 1 to 6
options.SignRestrictions = {
    sprintf('y(%d,1:6,1) > 0', i_gdp),  ...   % demand: GDP up
    sprintf('y(%d,1:6,1) > 0', i_def),  ...   % demand: deflator up
    sprintf('y(%d,1:6,1) > 0', i_rate), ...   % demand: yield up
    sprintf('y(%d,1:6,2) > 0', i_gdp),  ...   % supply: GDP up
    sprintf('y(%d,1:6,2) < 0', i_def),  ...   % supply: deflator down
    sprintf('y(%d,1:6,3) > 0', i_rate), ...   % MP: yield up
    sprintf('y(%d,1:6,3) < 0', i_gdp),  ...   % MP: GDP down
    sprintf('y(%d,1:6,3) < 0', i_def),  ...   % MP: deflator down
    sprintf('y(%d,1:6,4) > 0', i_soma), ...   % RFBS: SOMA up
    sprintf('y(%d,1:6,4) > 0', i_liab), ...   % RFBS: reserves + ON RRP up
    sprintf('y(%d,1:6,4) < 0', i_rrp),  ...   % RFBS: ON RRP share down
    sprintf('y(%d,1:6,4) < 0', i_ts),   ...   % RFBS: term spread down
    sprintf('y(%d,1:6,5) > 0', i_soma), ...   % OFBS: SOMA up
    sprintf('y(%d,1:6,5) > 0', i_liab), ...   % OFBS: reserves + ON RRP up
    sprintf('y(%d,1:6,5) > 0', i_rrp),  ...   % OFBS: ON RRP share up
    sprintf('y(%d,1:6,5) < 0', i_ts)          % OFBS: term spread down
};

% The same signs at impact, rows = variables and columns = shocks. The sampler
% uses this matrix as a fast screen before the full check.
impact_signs = zeros(nv, nv);
impact_signs(i_gdp,  1:3) = [+1, +1, -1];
impact_signs(i_def,  1:3) = [+1, -1, -1];
impact_signs(i_rate, [1 3]) = [+1, +1];
impact_signs(i_soma, 4:5) = [+1, +1];
impact_signs(i_liab, 4:5) = [+1, +1];
impact_signs(i_rrp,  4:5) = [-1, +1];
impact_signs(i_ts,   4:5) = [-1, -1];
options.impact_signs           = impact_signs;
options.sign_regime_dependent  = true;
options.penalty_offdiagonal_on = false;
impact_signs_byregime          = [];   % workspace variable the sampler reads; [] sets no per-regime gate

%% Zero restrictions, all at impact
% Demand, supply, and MP shocks do not move SOMA within the month.
% The two balance-sheet shocks do not move GDP or the deflator within the month.
% zero_on_average imposes the zeros on the average transmission, the response
% evaluated at the regime average of the relative variances, so regime_idx
% is a placeholder.
options.zero_on_average  = true;
options.ZeroRestrictions = struct( ...
    'variable_idx', {i_soma, i_soma, i_soma, i_gdp, i_def, i_gdp, i_def}, ...
    'shock_idx',    {1,      2,      3,      4,     4,     5,     5}, ...
    'horizon',      {0,      0,      0,      0,     0,     0,     0}, ...
    'regime_idx',   {1,      1,      1,      1,     1,     1,     1});

%% Narrative restrictions on two windows
% Over each window the focal shock has the stated sign, and its contribution to
% log(reserves + ON RRP) exceeds in absolute value the sum of the
% contributions of the other labeled shocks ('strict').
% 'interval_hd' measures the contribution as the historical decomposition
% accumulated over the window.
row_of = @(yr, mo) find(filtered_data.dates == datetime(yr, mo, 1), 1) - options.lags;
window_covid_qe  = row_of(2020, 4) : row_of(2021, 2);   % inside regime 3
window_rrp_drain = row_of(2023, 9) : row_of(2024, 7);   % inside regime 4

SHOCK_RFBS = 4;
SHOCK_OFBS = 5;
LABELED    = 1:5;
options.nrr_dom_aggregation = 'interval_hd';
options.NarrativeRegimeRestrictions = struct( ...
    'window',         {window_covid_qe,              window_rrp_drain}, ...
    'regime_idx',     {3,                            4}, ...
    'shock_idx',      {SHOCK_RFBS,                   SHOCK_OFBS}, ...
    'variable_idx',   {i_liab,                       i_liab}, ...
    'sign',           {+1,                           -1}, ...
    'other_idx',      {setdiff(LABELED, SHOCK_RFBS), setdiff(LABELED, SHOCK_OFBS)}, ...
    'dominance',      {'strict',                     'strict'}, ...
    'dom_agg',        {'interval_hd',                'interval_hd'}, ...
    'anti_shock_idx', {[],                           []}, ...
    'anti_other_idx', {[],                           []}, ...
    'anti_dominance', {[],                           []});
options.use_nrr_omega = false;   % narrative restrictions enter through SR only

%% HARS-Z settings
options.ess_premix_rounds        = 5;      % pre-mix ESS rounds before the structural move
options.ess_store_rounds         = 5;      % N_ess, rotations retained per theta
options.ess_repair_budget        = 200;    % N_Q, global rotation repair attempts
options.theta_redraw_budget      = 300;    % N_theta, joint structural redraws
options.seed_rescue_rounds       = 25;     % initialization only
options.use_local_pert_fallback  = false;
options.upstream_stability_check = false;
options.enforce_stability        = false;

%% Run
data.y             = y;
data.varnames      = varnames;
data.filtered_data = filtered_data;

harsz;

save(save_file, 'output', 'options', 'data', '-v7.3');
fprintf('Results saved to %s\n', save_file);
