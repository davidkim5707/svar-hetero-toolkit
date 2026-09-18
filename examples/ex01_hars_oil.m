%==========================================================================
% ex01_hars_oil.m -- HARS on a three-variable monthly model of the oil market.
%
% Identifies an oil supply shock, an aggregate demand shock, and an
% oil-specific demand shock with sign restrictions and bounds on the impact
% supply elasticity, under three volatility regimes. This is the oil
% application of Kim and Zha, "Sharpening Economic Interpretation with HARS".
% replications/kim_zha2026_hars/ builds the paper's figures and tables from
% the same model.
%
% Reads   examples/data/oil_monthly.mat   1971:M1 to 2015:M12
% Writes  examples/output/ex01_hars_oil.mat   (output, options, data)
%         examples/output/ex01_hars_oil.pdf   impulse responses
% Seed    none. The draws depend on the state of the generator at the start.
% Runtime about 2 minutes for NDRAW = 1000 on the author's desktop with an
%         eight-worker pool. The paper uses NDRAW = 10000.
%
% Variables (VAR order)
%   1 oil production growth   2 economic activity index   3 real oil price
%
% Shocks (columns of the impact matrix)
%   1 oil supply   2 aggregate demand   3 oil-specific demand
%==========================================================================
clear; close all;

NDRAW = 1000;    % recorded draws
NBURN = 1000;    % burn-in draws

example_dir = fileparts(mfilename('fullpath'));
run(fullfile(example_dir, '..', 'setup.m'));
data_file   = fullfile(example_dir, 'data', 'oil_monthly.mat');
save_file   = fullfile(example_dir, 'output', 'ex01_hars_oil.mat');
figure_file = fullfile(example_dir, 'output', 'ex01_hars_oil.pdf');

%% Data
data_raw            = load(data_file);
y                   = data_raw.data;
varnames            = data_raw.varNames;
filtered_data       = data_raw;
filtered_data.dates = datetime(data_raw.dates, 'ConvertFrom', 'datenum');

%% VAR and sampler size
options.ny          = size(y, 2);
options.lags        = 24;
options.irf_horizon = 60;
options.constant    = 1;
options.timetrend   = 0;

options.nobs        = size(y, 1) - options.lags;   % observations after the initial lags
options.ndraw       = NDRAW;
options.nburn       = NBURN;
options.nsep        = 1;
options.nit         = 500;    % iterations of csminwel in the posterior mode search
options.hsnscale    = 0.05;   % scale of the MH proposal covariance

%% Structural model
options.tvA              = 0;                 % A0 common across regimes
options.noLmd            = 0;                 % the field must exist; its value is not used
options.A0_restriction   = true(options.ny);
options.a0_prior_setting = 1;                 % p(A0) proportional to |det A0|^(-n)

% Student-t structural errors, degrees of freedom fixed.
options.tparam = 5.703233415698 / 2;
options.tscale = 1;

%% Prior on the reduced form
% No Minnesota dummy observations. The fields must exist and may be empty.
options.minn_prior_tau    = [];
options.minn_prior_decay  = [];
options.minn_prior_lambda = [];
options.minn_prior_mu     = [];
options.minn_prior_omega  = [];
options.unitroot          = [];

%% Volatility regimes
% The sampler assigns a break month to the regime that ends there. Each date
% below is the last month of a regime, and the next regime starts one month
% later, in 1990:M1 and 2008:M1.
options.regimes = [datetime(1989, 12, 1), datetime(2007, 12, 1)];

%% Sign restrictions at impact
% A string y(i, h1:h2, j) > 0 restricts the response of variable i to shock j
% over horizons h1 to h2, with impact counted as horizon 1. A ratio of two
% responses bounds an elasticity.
options.SignRestrictions = {
    'y(1,1:1,1) < 0', ...                     % supply: production down
    'y(2,1:1,1) < 0', ...                     % supply: activity down
    'y(3,1:1,1) > 0', ...                     % supply: price up
    'y(1,1:1,2) > 0', ...                     % aggregate demand: production up
    'y(2,1:1,2) > 0', ...                     % aggregate demand: activity up
    'y(3,1:1,2) > 0', ...                     % aggregate demand: price up
    'y(1,1:1,2)/y(3,1:1,2) >= 0', ...         % production over price, lower bound
    'y(1,1:1,2)/y(3,1:1,2) <= 0.0258', ...    % production over price, upper bound
    'y(1,1:1,3) > 0', ...                     % oil-specific demand: production up
    'y(2,1:1,3) < 0', ...                     % oil-specific demand: activity down
    'y(3,1:1,3) > 0', ...                     % oil-specific demand: price up
    'y(1,1:1,3)/y(3,1:1,3) >= 0', ...         % production over price, lower bound
    'y(1,1:1,3)/y(3,1:1,3) <= 0.0258'         % production over price, upper bound
};

% The same signs at impact, rows = variables and columns = shocks. The sampler
% uses this matrix to initialize the rotation and as a fast screen.
%                        supply  agg. demand  oil demand
options.impact_signs = [   -1        +1          +1;      % oil production
                           -1        +1          -1;      % economic activity
                           +1        +1          +1];     % real oil price
options.sign_regime_dependent  = true;
options.penalty_offdiagonal_on = false;
impact_signs_byregime          = [];   % workspace variable the sampler reads; [] sets no per-regime gate

options.NarrativeRestrictions = [];    % no narrative restrictions in this example

%% HARS settings
options.ess_premix_rounds        = 3;      % pre-mix ESS rounds before the structural move
options.ess_store_rounds         = 5;      % N_ess, rotations retained per theta
options.ess_repair_budget        = 2000;   % N_Q, global rotation repair attempts
options.theta_redraw_budget      = 100;    % N_theta, joint structural redraws
options.use_local_pert_fallback  = false;
options.upstream_stability_check = false;
options.enforce_stability        = false;

%% Run
data.y             = y;
data.varnames      = varnames;
data.filtered_data = filtered_data;

hars;

save(save_file, 'output', 'options', 'data');
fprintf('Results saved to %s\n', save_file);

%% Impulse responses
shocknames = {'oil supply shock', 'aggregate demand shock', 'oil-specific demand shock'};
fig = plot_irf_bands(output.sign_irfs_temp, output.draw_weight, varnames, shocknames, ...
                     'Horizon', 17, 'XLabel', 'Horizon (months)');
exportgraphics(fig, figure_file, 'ContentType', 'vector');
fprintf('Figure saved to %s\n', figure_file);
