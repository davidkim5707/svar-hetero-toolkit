%==========================================================================
% run_monetary.m -- model, options, data, and sampler call for the monetary
% application (Kim & Zha, HARS). Builds the 6-variable monthly VAR (12 lags,
% regime breaks 1979:10 and 1990:01), the sign and narrative restrictions on
% the monetary policy shock, and the HARS/ESS options, then times and runs
% the sampler (hars.m). main.m sets NTHETA per run (100 baseline, 0 ablation).
%==========================================================================
%% Data loading and setting

% Add code/core and code/third_party to the MATLAB path
addpath(genpath(fullfile('..','..','code','core')), genpath(fullfile('..','..','code','third_party')));

% Load data and variable names
data_raw = load('data/Uhlig_Data_Updated.mat');
varnames = data_raw.varNames;

% Convert the date column to datetime
filtered_data = data_raw;
filtered_data.dates = datetime(data_raw.dates, 'ConvertFrom', 'datenum');

% Endogenous variables
y = data_raw.data;

sampler_script = 'hars';   % HARS sampler (renamed from svar_run_sign_elliptical_final.m)
if exist([sampler_script '.m'], 'file') ~= 2
    error('run_monetary:samplerMissing', ...
        'Sampler script "%s.m" not found on the path.', sampler_script);
end

%********************************************************
%% Define Options Struct for SVAR
%********************************************************

% === USER SETTINGS ===
% --------------------------------------------------------------
% Uhlig (2005):             regime_on = 0, prior_on = 0, svar_svsr_on = 0
% Carriero et al. (2024):   regime_on = 1, prior_on = 0, svar_svsr_on = 0
% Dawis' method:            regime_on = 1, prior_on = 0, svar_svsr_on = 1
% --------------------------------------------------------------
regime_on     = 1;    % 1 = enable regime-switching 
prior_on      = 0;    % 1 = use Minnesota prior
svar_svsr_on  = 1;    % 1 = Dawis (rotational invariance)
penalty_offdiagonal_on = 0; % 1 = Use penalty function to check sign-restrictions

% === MODEL STRUCTURE ===
options.ny            = size(y, 2);       % Number of obs ervable variables
options.lags          = 12;               % Number of lags
options.irf_horizon   = 61;               % Horizon for IRFs
options.constant      = 0;                % Include constant in VAR
options.timetrend     = 0;                % Include time trend
options.nexogenous    = 0;                % No exogenous variables
options.exogenous     = [];               % Empty exogenous matrix
options.non_explosive_ = 0;               % Do not enforce non-explosiveness

% === SAMPLING SETTINGS ===
options.firstobs      = options.lags + 1;                     % First observation for estimation
options.presample     = 0;                                    % Number of presample periods
options.ndraw         = 10000;                                % Total recorded draws
options.nsep          = 1;                                    % Draws per stored draw
options.nburn         = 1000;                                 % Burn-in
options.max_compute   = 1;                                    % Optimization algorithm (1 = csminwel)
options.nit           = 100;                                  % Max iterations
options.hsnscale      = 0.05;                                 % Hessian scaling factor

% === DATA DIMENSIONS ===
options.nobs          = size(y, 1) - options.firstobs + 1;    % Usable observations
options.nunits        = size(y, 3);                           % Number of units (1 if time series)

% === PRIOR SETTINGS ===
if prior_on == 1
    options.dummy               = 1;            % Use Minnesota dummy prior
    options.minn_prior_tau      = 3;            % Overall tightness
    options.minn_prior_decay    = 0.5;          % Lag decay
    options.minn_prior_lambda   = 5;            % Sum-of-coefficient prior
    options.minn_prior_mu       = 1;            % Co-persistence
    options.minn_prior_omega    = 0;            % Prior on shock variances
    options.unitroot            = ones(options.ny, 1);  % Treat all as persistent
else
    options.dummy               = [];
    options.minn_prior_tau      = [];
    options.minn_prior_decay    = [];
    options.minn_prior_lambda   = [];
    options.minn_prior_mu       = []; 
    options.minn_prior_omega    = [];
    options.unitroot            = [];
end

% === REGIME SETTINGS ===
if regime_on == 1
    % Convention: a break date is the LAST obs of the PREVIOUS regime (rfvar3
    % assigns obs t to the new regime iff t > breakInd). To make the labelled
    % event month the FIRST obs of the new regime, set the break one month prior.
    options.regimes = [
        datetime(1979, 9, 1), ...   % Volcker 1979:10 = first month of new (monetarist) regime
        datetime(1989, 12, 1)       % 1990:01 = first month of new regime
    ];
else
    options.regimes = [];
end

% === STRUCTURAL IDENTIFICATION SETTINGS ===
options.tvA                   = 0;  % No time-varying A0
options.noLmd                 = 0;  % Estimate shock variances (homoskedastic)
options.A0_restriction = true(options.ny);  % free all A0 elements (no exclusion restrictions)
options.a0_prior_setting = 1; % 0: Brunnermeier et al. (2021); 1: p(A0) proportional to |det A0|^(-n); 2: Carriero et al. (2024)
options.tparam         = 5.703233415698 / 2;  % Scale for A0 prior
options.tscale         = 1;

% === SIGN RESTRICTIONS ===
options.sign_horizon = 6;
options.SignRestrictions = {
    sprintf('y(2,1:%d,1) < 0', options.sign_horizon), ...
    sprintf('y(3,1:%d,1) < 0', options.sign_horizon), ...
    sprintf('y(5,1:%d,1) < 0', options.sign_horizon), ...  % Interest rate ↑
    sprintf('y(6,1:%d,1) > 0', options.sign_horizon)
};

options.sign_regime_dependent = (svar_svsr_on == 1);  % Rotational invariance
options.sign_inneriteration = 2000;
options.outeriteration = 2;
options.penalty_offdiagonal_on = (penalty_offdiagonal_on == 1); 
options.store_all_valid_Q = 1; 

%********************************************************
% ESS-SPECIFIC SETTINGS
%********************************************************
% The upstream stability filter is off in every application
% (options.upstream_stability_check = false), so the posterior is not
% truncated for stationarity.
options.ess_premix_rounds        = 3;        % pre-mixing rounds before the structural move
options.ess_repair_budget        = 2000;     % global Q-repair budget (N_Q) -> repair_k
options.upstream_stability_check = false;
options.enforce_stability        = false;
options.use_local_pert_fallback  = true;

% --- HARS sampler parameters (Algorithm, Sec 3) ---
% N_ess: local ESS rotations RETAINED per structural draw (Sec 3 step 3).
% >1 amortises the structural move over several Q-draws; 1 = one per sweep.
options.ess_store_rounds    = 5;

% N_theta: max JOINT (A0,Lambda,A+) structural redraws to restore the current Q
% after an MCMC-accepted completed draw strands it (Sec 3 step 1 / Prop 1).
% 0 disables -> every stale event falls to the global Q-repair (avoid). Keep > 0.
if ~exist('NTHETA', 'var'), NTHETA = 100; end   % main.m sets NTHETA per run (100 = baseline, 0 = ablation)
options.theta_redraw_budget = NTHETA;

% Impact sign matrix for ESS column-by-column initialization
options.impact_signs = zeros(options.ny, options.ny);
options.impact_signs(2, 1) = -1;   % variable 2 responds negatively to shock 1
options.impact_signs(3, 1) = -1;   % variable 3 responds negatively to shock 1
options.impact_signs(5, 1) = -1;   % variable 5 responds negatively to shock 1
options.impact_signs(6, 1) = +1;   % variable 6 responds positively to shock 1
impact_signs_byregime       = [];     % [] -> no per-regime impact gate

% === NARRATIVE RESTRICTIONS ===
% Uhlig data, monetary shock (col 1), Oct 1979 Volcker anchor on the federal
% funds rate (var 6). Sign: realized shock POSITIVE on the anchor date.
% type = 'dominance' is the AD-RR Type B (overwhelming-contributor) check,
% matching the AD-RR headline specification (NSR-5, fn.16).
options.NarrativeRestrictions.time_index     = 166;     % time (observation index)
options.NarrativeRestrictions.shock_index    = 1;       % shock (column in structural shocks)
options.NarrativeRestrictions.sign           = 1;   
options.NarrativeRestrictions.variable_index = 6;
options.NarrativeRestrictions.type           = {'dominance'};   % Type B (AD-RR "overwhelming")

%********************************************************
%% NRR OMEGA IMPORTANCE WEIGHTS (post-hoc, storage stage)
%********************************************************
options.use_nrr_omega = false;    % AD-RR 1/omega_hat importance weighting (off)

% === DATA STRUCT ===
data.y              = y;
data.varnames       = varnames;
data.filtered_data  = filtered_data;

% === RUN SVAR ===
clear estimate_nrr_omega_fast

% End-to-end wall clock (burn-in + sampling + post-processing) for the paper.
run_start_time = tic;
run(sampler_script);                          % runs hars.m
run_total_sec  = toc(run_start_time);

fprintf('End-to-end run time: %.2f s (%.2f min)\n', run_total_sec, run_total_sec/60);