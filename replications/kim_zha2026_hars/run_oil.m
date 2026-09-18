%==========================================================================
% run_oil.m -- model, options, data, and sampler call for the oil application
% (Kim & Zha, HARS). Builds the 3-variable monthly VAR (24 lags, regime breaks
% 1990:01 and 2008:01), the sign and elasticity restrictions on the three oil
% shocks, and the HARS/ESS options, then times and runs the sampler (hars.m).
% main_oil.m sets NTHETA per run (100 baseline, 0 ablation).
%==========================================================================
%% Data loading and setting

% Add code/core and code/third_party to the MATLAB path
addpath(genpath(fullfile('..','..','code','core')), genpath(fullfile('..','..','code','third_party')));

% Load data from CSV file
data = load('data/Kilian_Data_Updated.mat');

% Extract variable names from column 2 to the last column
varnames = data.varNames;

% Convert date column to datetime format
filtered_data = data;
filtered_data.data = data.data;
filtered_data.dates = datetime(data.dates, 'ConvertFrom', 'datenum');

% Drop the 'date' column and define y
y = data.data;

%********************************************************
%% Define Options Struct for SVAR - Kilian-style Model
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
narrative_on  = 0;    % 1 = apply narrative restrictions (and tag filename), 0 = none
penalty_offdiagonal_on = 0;

% === BASIC MODEL SETTINGS ===
options.ny            = size(y, 2);
options.lags          = 24;
options.irf_horizon   = 60;
options.constant      = 1;
options.timetrend     = 0;
options.firstobs      = options.lags + 1;

% === SAMPLING AND COMPUTATION SETTINGS ===
options.presample     = 0;
options.ndraw         = 10000;
options.nsep          = 1;
options.nburn         = 1000;
options.max_compute   = 1;
options.nit           = 500;
options.hsnscale      = 0.05;

% === DATA DIMENSIONS ===
options.nobs          = size(y, 1) - options.firstobs + 1;
options.nunits        = size(y, 3);

% === PRIOR SETTINGS ===
if prior_on == 1
    options.dummy               = 1;
    options.minn_prior_tau      = 3;
    options.minn_prior_decay    = 0.5;
    options.minn_prior_lambda   = 5;
    options.minn_prior_mu       = 1;
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
end

% === REGIME SETTINGS ===
if regime_on == 1
   % Convention: a break date is the LAST obs of the PREVIOUS regime (rfvar3
   % assigns obs t to the new regime iff t > breakInd). Set one month earlier so
   % the labelled month is the FIRST obs of the new regime.
   options.regimes = [
       datetime(1989, 12, 1), ...  % 1990:01 = first month of new regime
       datetime(2007, 12, 1)       % 2008:01 = first month of new regime
   ];
else
    options.regimes = [];
end

% === STRUCTURAL IDENTIFICATION SETTINGS ===
options.tvA                = 0;
options.noLmd              = 0;
options.A0_restriction     = true(options.ny);
options.a0_prior_setting   = 1;
options.tparam             = 5.703233415698 / 2;
options.tscale             = 1;

% === SIGN RESTRICTIONS (Kilian Style) ===
options.SignRestrictions = {
    'y(1,1:1,1) < 0', ...  % Shock 1: variable 1 negative
    'y(2,1:1,1) < 0', ...  % Shock 1: variable 2 negative
    'y(3,1:1,1) > 0', ...  % Shock 1: variable 3 positive
    'y(1,1:1,2) > 0', ...  % Shock 2: variable 1 positive
    'y(2,1:1,2) > 0', ...  % Shock 2: variable 2 positive
    'y(3,1:1,2) > 0', ...  % Shock 2: variable 3 positive
    'y(1,1:1,2)/y(3,1:1,2) >= 0', ...      % elasticity lower bound
    'y(1,1:1,2)/y(3,1:1,2) <= 0.0258', ... % elasticity upper bound
    'y(1,1:1,3) > 0', ...  % Shock 3: variable 1 positive
    'y(2,1:1,3) < 0', ...  % Shock 3: variable 2 negative
    'y(3,1:1,3) > 0', ...  % Shock 3: variable 3 positive
    'y(1,1:1,3)/y(3,1:1,3) >= 0', ...      % elasticity lower bound
    'y(1,1:1,3)/y(3,1:1,3) <= 0.0258'      % elasticity upper bound
};
options.sign_regime_dependent  = (svar_svsr_on == 1);
options.sign_inneriteration    = 2000;        % not used by ESS, kept for compatibility
options.outeriteration         = 2;           % not used by ESS, kept for compatibility
options.penalty_offdiagonal_on = (penalty_offdiagonal_on == 1);
options.store_all_valid_Q      = 1;

%********************************************************
% ESS-SPECIFIC SETTINGS
%********************************************************
options.ess_premix_rounds        = 3;        % pre-mixing rounds before the structural move
options.ess_repair_budget        = 2000;     % global Q-repair budget (N_Q) -> repair_k
options.upstream_stability_check = false;
options.enforce_stability        = false;
options.use_local_pert_fallback  = false;

% --- HARS sampler parameters (Algorithm, Sec 3) ---
% N_ess: local ESS rotations RETAINED per structural draw (Sec 3 step 3).
% >1 amortises the structural move over several Q-draws; 1 = one per sweep.
options.ess_store_rounds    = 5;

% N_theta: max JOINT (A0,Lambda,A+) structural redraws to restore the current Q
% after an MCMC-accepted completed draw strands it (Sec 3 step 1 / Prop 1).
% 0 disables -> every stale event falls to the global Q-repair (avoid). Keep > 0.
if ~exist('NTHETA','var'), NTHETA = 100; end
options.theta_redraw_budget = NTHETA;

% === IMPACT SIGN MATRIX (for Algorithm 4 initialization) ===
%                         Supply  AggDem  OilDem
options.impact_signs = [   -1      +1      +1;    % oil production
                           -1      +1      -1;    % real activity
                           +1      +1      +1];   % real oil price
impact_signs_byregime       = [];     % [] -> no per-regime impact gate

% === NARRATIVE RESTRICTIONS (Antolin-Diaz & Rubio-Ramirez 2018, oil) ===
if narrative_on == 1
    tidx = @(yy,mm) find(filtered_data.dates == datetime(yy,mm,1)) - options.lags;

    % Alternative NSR3, Aug 1990 only (AD-RR p.2820 + Online Appendix A):
    % at the Gulf War onset, aggregate demand contributes the least to the
    % unexpected movement in the real oil price (Type A dominance_min, no sign restriction).
    options.NarrativeRestrictions.time_index     = tidx(1990,8);
    options.NarrativeRestrictions.shock_index    = 2;   % aggregate demand
    options.NarrativeRestrictions.variable_index = 3;   % real oil price
    options.NarrativeRestrictions.sign           = 0;   % no sign restriction
    options.NarrativeRestrictions.type           = {'dominance_min'};
else
    options.NarrativeRestrictions = [];
end

% === DATA STRUCT ===
data.y              = y;
data.varnames       = varnames;
data.filtered_data  = filtered_data;

%********************************************************
%% Run sampler
%********************************************************
sampler_script = 'hars';
if exist([sampler_script '.m'], 'file') ~= 2
    error('run_oil:samplerMissing', ...
        'Sampler script "%s.m" not found on the path.', sampler_script);
end

% End-to-end wall clock (burn-in + sampling + post-processing) for the paper.
delete(gcp('nocreate'));
run_start_time = tic;
run(sampler_script);                       % runs hars.m
run_total_sec  = toc(run_start_time);

fprintf('End-to-end run time: %.2f s (%.2f min)\n', run_total_sec, run_total_sec/60);