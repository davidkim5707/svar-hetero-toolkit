# `options` reference

Both samplers run as scripts. The caller puts two structs in the workspace, `data` and `options`, and the sampler leaves a struct `output`. This file lists every field of `options` that the code under `code/core/` and `code/third_party/` reads. Helper files are cited by file name.

Conventions used in the tables.

- **Default** is the value the code uses when the field is absent or empty. "required" means the code indexes the field unconditionally, and a missing field stops the run with an error. "none" means the field is skipped when it is absent or empty.
- `get_opt(options, name, default)` returns the default when the field is absent **or empty**. Setting a field to `[]` is the same as not setting it for every option read through `get_opt`.
- **Used by** says whether `hars.m`, `harsz.m`, or a helper reachable from them reads the field.
- The driver columns are **ex02** for `examples/ex02_harsz_balance_sheet.m` (calls `harsz`), and **mon**, **oil**, **fis** for `run_monetary.m`, `run_oil.m`, `run_fiscal.m` under `replications/kim_zha2026_hars/` (the last three call `hars`). `examples/ex01_hars_oil.m` uses the values of **oil**, with `ndraw` 1000. "not set" means the driver does not assign the field.

The `data` struct has three fields, all required. They are `data.y` (T x ny matrix), `data.varnames` (used for printing and copied to `output.varnames`), and `data.filtered_data` (a table or struct whose field `dates` is a datetime vector with one entry per row of `data.y`). `filtered_data.dates` is compared with `options.regimes` to locate the regime breaks.

The tables below document 87 top-level option fields. 73 are read in a standard run, 9 sit behind the `thetero` and `tthr` switches, and 5 are read only on the `tvA` path, which both samplers refuse. A further 21 sub-fields of the struct-valued options are documented in sections 6 and 7.

---

## 1. Data and VAR

| Name | Default | Type / allowed values | Read in | Used by | Meaning | Drivers |
|---|---|---|---|---|---|---|
| `ny` | required | positive integer scalar | `hars.m`, `harsz.m` | both | Number of endogenous variables, which must equal the number of columns of `data.y`. | ex02 `size(y,2)` (8). mon `size(y,2)` (6). oil `size(y,2)` (3). fis `size(y,2)` (5). |
| `lags` | required | positive integer scalar | `hars.m`, `harsz.m` | both | Number of VAR lags. It is also added to `NarrativeRestrictions.time_index` to locate the regime of each narrative date. | ex02 6. mon 12. oil 24. fis 4. |
| `constant` | required | 0 or 1 | `hars.m`, `harsz.m`, `bvar_posterior.m` | both | When nonzero a column of ones is added to the exogenous regressors. | ex02 1. mon 0. oil 1. fis 1. |
| `nobs` | required | positive integer scalar | `hars.m`, `harsz.m` | both | Number of usable observations. It closes the last regime in `regime_vec = [0, breakInd - lags, nobs]` (`harsz.m`) and sizes the latent log-scale matrix `delta0` (`harsz.m`). | ex02 `size(y,1) - options.lags`. mon, oil, fis `size(y,1) - options.firstobs + 1`, the same number. |
| `irf_horizon` | required | positive integer scalar | `hars.m`, `harsz.m` | both | Number of horizons stored for every impulse response, with index 1 equal to impact. | ex02 60. mon 61. oil 60. fis 60. |
| `timetrend` | `0` | 0, 1, or 2 | `hars.m`, `harsz.m`, `bvar_posterior.m` | both | 1 adds a linear trend `(1:T)'/T` to the exogenous regressors, and 2 adds its square as well. | all four: 0. |
| `exogenous` | none | T x k matrix | `hars.m`, `harsz.m`, `bvar_posterior.m` | both | Extra exogenous regressors appended after the constant and trend columns. | ex02 `[]`. mon `[]`. oil not set. fis `[]`. |

## 2. Prior

All Minnesota fields below must **exist** in `options`, because `bvar_posterior.m` and `bvar_posterior.m` index them unconditionally. They may be empty. When `minn_prior_tau`, `minn_prior_decay`, `minn_prior_lambda`, and `minn_prior_mu` are all empty and `minn_prior_omega` is empty or 0, no dummy observations are built (`bvar_posterior.m`). That is the setting of the three `hars` drivers.

| Name | Default | Type / allowed values | Read in | Used by | Meaning | Drivers |
|---|---|---|---|---|---|---|
| `a0_prior_setting` | required | 0, 1, or any other value | `bvar_posterior.m` | both | Selects the prior on A0. 0 is an independent Normal with mean `100*I` and standard deviation 200 on every element. 1 is the improper prior proportional to `abs(det(A0))^(-ny)`. Any other value is a Normal prior with zero mean and row variance `Sstar(i,i)`, labelled Carriero et al. (2024) in the code. | ex02 0. mon 1. oil 1. fis 1. |
| `minn_prior_tau` | required field, may be `[]` | scalar or `[]` | `bvar_posterior.m` | both | Overall tightness of the Minnesota dummy observations (`mnprior.tight` in `varprior.m`). | ex02 3. mon, oil, fis `[]`. |
| `minn_prior_decay` | required field, may be `[]` | scalar or `[]` | `bvar_posterior.m` | both | Lag decay of the Minnesota prior. The dummy for lag l is scaled by `l^decay` (`varprior.m`). | ex02 0.5. mon, oil, fis `[]`. |
| `minn_prior_lambda` | required field, may be `[]` | scalar or `[]` | `hars.m`, `harsz.m`, `bvar_posterior.m` | both | Weight on the co-persistence dummy observation. It is passed to `rfvar3` for the seed reduced-form model and to `varprior` as `urprior.lambda`. | ex02 5. mon, oil, fis `[]`. |
| `minn_prior_mu` | required field, may be `[]` | scalar or `[]` | `hars.m`, `harsz.m`, `bvar_posterior.m` | both | Weight on the own-persistence dummy observations. It is passed to `rfvar3` for the seed model and to `varprior` as `urprior.mu`. | ex02 1. mon, oil, fis `[]`. |
| `minn_prior_omega` | required field, may be `[]` | nonnegative scalar or `[]` | `bvar_posterior.m` | both | Weight `vprior.w` on the dummy observations for the residual covariance. 0 or empty adds none (`varprior.m`). | ex02 0. mon, oil, fis `[]`. |
| `unitroot` | field required. `ones(ny,1)` when empty | ny x 1 vector of 0 and 1 | `bvar_posterior.m` | both | Entry i equal to 1 treats variable i as persistent. It scales the first-lag Minnesota dummy (`varprior.m`) and switches the own-persistence dummy on for that variable (`varprior.m`). | ex02 `[1;1;1;0;1;1;0;0]`. mon, oil, fis `[]`. |
| `vprior_sig` | none | ny x 1 vector | `bvar_posterior.m` | both | Prior residual scale per variable for the Minnesota dummies. It overrides the `Sstar` scale when set, and the code errors if it does not have ny entries. | ex02, mon, oil not set. fis `[]`. |
| `Sstar` | overwritten by the sampler | ny x ny matrix | written at `hars.m`, `harsz.m`. Read at `bvar_posterior.m` | both | The sampler sets this field to `compute_Sstar(y, lags)`, the variances of univariate AR(lags) residuals, so a value supplied by the caller is replaced. | not set by any driver. |
| `use_Sstar_for_vprior` | overwritten by the sampler with `true` | logical | written at `hars.m`, `harsz.m`. Read at `bvar_posterior.m` | both | When true and `vprior_sig` is empty, the dummy scale is `sqrt(diag(Sstar))`. Because the sampler always sets it to true, the AR-residual fallback at `bvar_posterior.m` is not reached from `hars` or `harsz`. | not set by any driver. |

## 3. Posterior mode and MH

| Name | Default | Type / allowed values | Read in | Used by | Meaning | Drivers |
|---|---|---|---|---|---|---|
| `ndraw` | required | positive integer | `hars.m`, `harsz.m` | both | Number of stored draws. The main loop runs until this many draws are collected (`harsz.m`). | ex02 `NDRAW` = 10000. mon, oil, fis 10000. |
| `nburn` | required | nonnegative integer | `hars.m`, `harsz.m` | both | Number of burn-in sweeps. Burn-in is skipped when it is 0 (`harsz.m`). | all four: 1000. |
| `nsep` | required | positive integer | `hars.m`, `harsz.m` | both | Number of structural MH sweeps between two storage stages (`harsz.m`). | all four: 1. |
| `nit` | required | positive integer | `hars.m`, `harsz.m` | both | Maximum number of `csminwel` iterations in the posterior mode search. | ex02 100. mon 100. oil 500. fis 100. |
| `hsnscale` | required | positive scalar | `hars.m`, `harsz.m` | both | Scale applied to the inverse Hessian from `csminwel` to form the random-walk MH proposal covariance `Sigma = hsnscale * H`. | all four: 0.05. |

## 4. Structural model (A0, Lambda, regimes, error distribution)

| Name | Default | Type / allowed values | Read in | Used by | Meaning | Drivers |
|---|---|---|---|---|---|---|
| `A0_restriction` | required field. `[]` gives `tril(true(ny))` | ny x ny logical or numeric matrix, or `[]` | `hars.m`, `harsz.m` | both | Nonzero entries mark the free elements of A0 (`lcA0 = A0_restriction ~= 0`). An empty value makes A0 lower triangular. | all four: `true(options.ny)`. |
| `regimes` | required field, may be `[]` | datetime row vector, or `[]` | `hars.m`, `harsz.m` | both | Volatility regime break dates. Each date is matched against `data.filtered_data.dates`, and the matched observation is the last one of the regime that ends there. An empty value gives a single regime. | ex02 2014-11-01, 2020-02-01, 2022-05-01. mon 1979-09-01, 1989-12-01. oil 1989-12-01, 2007-12-01. fis 1971-01-01, 1996-10-01 (`[]` when `REGIME_ON = 0`). |
| `tparam` | field required by `hars.m` and `harsz.m`, may be `[]` | positive scalar or `[]` | `hars.m`, `harsz.m` | both | When non-empty and `regimes` is non-empty, the latent log scales are drawn with `drawt(eout, tparam, tscale^2*tparam)`, an inverse-gamma update with shape `tparam`. The drivers set it to half the Student-t degrees of freedom. `2*tparam` is also the default of `nrr_omega_nu`. | ex02 `5.840451331958/2`. mon, oil, fis `5.703233415698/2`. |
| `tscale` | required when `tparam` and `regimes` are non-empty | positive scalar | `hars.m`, `harsz.m` | both | Enters the rate of the `drawt` update as `tscale^2 * tparam`. | all four: 1. |
| `alpha` | `[]` | scalar or vector | `hars.m`, `harsz.m` | both | Read only when `tparam` is empty. It is passed to `drawdelta` as the outlier probability of a normal-mixture error law. | not set by any driver. |
| `K` | `[]` | scalar or vector | `hars.m`, `harsz.m` | both | Read only when `tparam` is empty. It is passed to `drawdelta` as the outlier scale factor. An empty `K` is the Gaussian case and keeps the latent log scales at zero (`drawdelta.m`). | not set by any driver. |
| `homoskedastic_shocks` | none | vector of row indices of Lambda in 1..ny | `hars.m`, `harsz.m` | both | The relative variances of the listed rows are fixed at 1 in every regime. The rows are removed from the free Lambda mask `lcLmd` and their seed value is set to 1. | ex02 `[]`. others not set. |
| `fix_first_regime` | none (treated as off) | 1 switches it on | `hars.m`, `harsz.m`, `bvar_posterior.m` | both | When equal to 1 the relative variances of regime 1 are fixed at 1 and removed from the free parameters, and the Dirichlet prior applies to the remaining regimes (`bvar_posterior.m`). | not set by any driver. |
| `tvA` | treated as false | logical or 0/1 | `hars.m`, `harsz.m`, `eval_posterior.m` (error message only) | both | A true value stops the run with an error. Both samplers then set the internal `tvA = 0` and keep A0 common across regimes. | all four: 0. |
| `noLmd` | required | any | `hars.m`, `harsz.m` | both | The field is copied to a local variable that no later line reads, so its value has no effect, but a missing field stops the run. | all four: 0. |

## 5. Sign restrictions

| Name | Default | Type / allowed values | Read in | Used by | Meaning | Drivers |
|---|---|---|---|---|---|---|
| `SignRestrictions` | none | cell array of strings, grammar below | `hars.m`, `harsz.m` | both | Sign and elasticity restrictions on impulse responses. They are parsed once by `parse_sign_restrictions.m` and checked under the average transmission and in every regime (`check_sr_pass_cached.m`). | ex02 16 strings on shocks 1 to 5, horizons 1:6. mon 4 strings on shock 1, horizons 1:6. oil 13 strings on shocks 1 to 3, impact only, 4 of them elasticity bounds. fis 4 strings on 2 shocks, horizons 1:4. |
| `sign_regime_dependent` | required when any of `SignRestrictions`, `NarrativeRestrictions`, `NarrativeRegimeRestrictions`, or (in `harsz`) `ZeroRestrictions` is non-empty | logical or 0/1 | `hars.m`, `harsz.m` | both | Order of scaling and rotation in the regime impact matrix. 0 gives `(A0inv*Q')*sqrt(Lambda_r)`. Any other value gives `(A0inv*sqrt(Lambda_r))*Q'` (`check_sr_pass_cached.m`, `harsz.m`). | ex02 `true`. mon, oil, fis `(svar_svsr_on == 1)`, which is true. |
| `penalty_offdiagonal_on` | required under the same condition as `sign_regime_dependent` | logical | `hars.m`, `harsz.m` | both | When true a draw is rejected in a regime if, for some restricted shock j, `abs(QLQ(j,j)) < 3*max(abs(QLQ(j,k)))` over k different from j, where `QLQ = Q*diag(lambda_r)*Q'` (`check_sr_pass_cached.m`). | ex02 `false`. mon, oil, fis `(penalty_offdiagonal_on == 1)`, which is false. |
| `impact_signs` | `zeros(ny, ny)`. Required when `sign_active_regimes` is set (`hars.m`, `harsz.m`) | ny x ny matrix with entries +1, -1, 0. Rows are variables, columns are shocks | `hars.m`, `harsz.m` | both | Impact signs used as a fast screen before the full check (`check_sr_pass_cached.m`) and to guide the column-by-column search for the first rotation. An all-zero matrix switches the screen off. | ex02 8 x 8, signs of the 16 restrictions at impact. mon 6 x 6 with four entries in column 1. oil full 3 x 3. fis 5 x 5 with two entries in each fiscal column. |
| `sign_active_regimes` | none | `containers.Map` with numeric shock ids as keys and vectors of regime indices as values | `hars.m`, `harsz.m` | both | For every shock that is a key of the map, its sign restrictions and its column of `impact_signs` are enforced only in the listed regimes. The sampler builds `SR_parsed_byregime` and `impact_signs_byregime` from it. | not set by any driver. |

### Grammar of `SignRestrictions`

`y(i, h, s)` is the response of variable `i` at horizon index `h` to shock `s`. Horizon index 1 is impact, because the check indexes the second dimension of the IRF array directly (`checkrestrictions_per_shock_fast.m`) and that array starts at impact (`check_sr_pass_cached.m`). Each string is trimmed and matched against four patterns in this order (`parse_sign_restrictions.m`). White space is allowed around every token.

| Form | Pattern | Operators | Right-hand side |
|---|---|---|---|
| sign, one horizon | `y(i,h,s) OP 0` | `>`, `>=`, `<`, `<=` | the literal `0` only |
| sign, horizon range | `y(i,h1:h2,s) OP 0` | `>`, `>=`, `<`, `<=` | the literal `0` only |
| elasticity, one horizon | `y(i,h,s)/y(k,h,s) OP c` | `>=`, `<=` | unsigned number `c`, decimals and exponents allowed |
| elasticity, horizon range | `y(i,h1:h2,s)/y(k,h1:h2,s) OP c` | `>=`, `<=` | unsigned number `c` |

Notes taken from the code.

- In the elasticity forms the horizon and the shock of the denominator must repeat those of the numerator (back-references `\2`, `\3`, `\4` in the patterns).
- The constant `c` cannot carry a sign, because the pattern is `[0-9]*\.?[0-9]+([eE][+-]?\d+)?`.
- A string that matches no pattern produces the warning `parse_sign_restrictions:unknown` and is dropped (`parse_sign_restrictions.m`).
- A range restriction must hold at every horizon in `h1:h2`. Strict operators use a tolerance of `1e-12` (`checkrestrictions_per_shock_fast.m`).
- For each shock the checker tries the column as drawn and then the column multiplied by -1, and records the orientation in `fsign` (`checkrestrictions_per_shock_fast.m`).
- The sampler reads the shock id of each string with the regular expression `,\s*(\d+)\)` (`hars.m`, `harsz.m`), and it raises the cached horizon only from range strings, with the expression `y\(\d+,\s*\d+:(\d+)` (`hars.m`, `harsz.m`). A single-horizon string with `h > 1` does not raise the cached horizon.

Examples accepted by the parser, taken from the drivers.

```matlab
'y(5,1:6,1) > 0'                 % ex02, sprintf('y(%d,1:6,1) > 0', i_gdp)
'y(2,1:6,1) < 0'                 % run_monetary.m
'y(1,1:1,2)/y(3,1:1,2) <= 0.0258' % run_oil.m
```

## 6. Zero restrictions (`harsz.m` only)

`hars.m` contains no read of either field. Its header states that `options.ZeroRestrictions` is ignored (`hars.m`).

| Name | Default | Type / allowed values | Read in | Used by | Meaning | Drivers |
|---|---|---|---|---|---|---|
| `ZeroRestrictions` | none | struct array, sub-fields below | `harsz.m`, passed to `build_zero_constraint.m` | harsz | Each element sets one impulse response to zero exactly. Rotations are built column by column in the null space of the constraints, and when zeros are present each stored draw carries the volume element in `output.draw_weight`. | ex02 7 elements, all at horizon 0: shocks 1, 2, 3 on SOMA and shocks 4, 5 on GDP and the deflator. |
| `zero_on_average` | `true` | logical | `harsz.m` | harsz | True imposes the zeros on the stored average IRF, which uses `A0inv` with Lambda set to 1, so one constraint covers all regimes. False imposes each zero on the transmission of the regime given by `regime_idx` (`build_zero_constraint.m`). | ex02 `true`. |

Sub-fields of `ZeroRestrictions(k)`.

| Sub-field | Default | Type | Read in | Meaning |
|---|---|---|---|---|
| `variable_idx` | required | integer in 1..ny | `build_zero_constraint.m` | Variable whose response is set to zero. |
| `shock_idx` | required | integer in 1..ny | `harsz.m`, `build_zero_constraint.m` | Shock (column of the impact matrix) that carries the zero. The sampler also adds these shocks to the list of stored shocks (`harsz.m`). |
| `horizon` | field required by `harsz.m`. `0` inside `build_zero_constraint.m` | integer >= 0 | `harsz.m`, `build_zero_constraint.m` | Economic horizon of the zero, with 0 equal to impact. The sampler raises the cached horizon to `max(horizon) + 1`. |
| `regime_idx` | required when there is more than one regime, otherwise 1 | integer in 1..number of Lambda columns | `build_zero_constraint.m` | Regime whose transmission carries the zero. With `zero_on_average = true` the value is validated but does not change the constraint. |

## 7. Narrative restrictions

| Name | Default | Type / allowed values | Read in | Used by | Meaning | Drivers |
|---|---|---|---|---|---|---|
| `NarrativeRestrictions` | none | scalar struct with vector sub-fields, one entry per anchor date | `hars.m`, `harsz.m` | both | Point-date narrative restrictions: a sign on the realized labelled shock at each date, and optionally a condition on its historical-decomposition contribution at that date. | ex02 not set. mon one anchor. oil `[]` (the block is behind `narrative_on = 0`). fis 6 anchors, 2 revenue and 4 spending. |
| `NarrativeRegimeRestrictions` | none | struct array, one element per window | `hars.m`, `harsz.m` | both | Window narrative restrictions: a sign on the window mean of the labelled shock, and optionally dominance conditions on its contribution over the window (`check_narrative_regime_avg.m`). | ex02 2 windows. others not set. |
| `nrr_dom_aggregation` | none | string. Only `'interval_hd'` has an effect | `hars.m`, `harsz.m` | both | When equal to `'interval_hd'` and window restrictions exist, the cached IRF horizon is raised to the longest window length. The aggregation rule itself is read from the sub-field `dom_agg` of each window, not from this option. | ex02 `'interval_hd'`. others not set. |
| `use_nrr_omega` | `true` | logical | `hars.m`, `harsz.m` | both | True estimates `omega_hat` for every stored draw and applies `1/omega_hat` importance weights to the summaries when narrative restrictions exist. False leaves the narrative restrictions as an accept or reject check only. | ex02 `false`. mon `false`. oil not set (no narrative restriction is active, so the default has no effect). fis `false`. |
| `nrr_omega_M` | `2000` | positive integer | `hars.m`, `harsz.m` | both | Number of simulated shock draws used to estimate `omega_hat` for each stored draw. | not set by any driver. |
| `nrr_omega_law` | `'student'` | `'student'` or `'gaussian'` (`estimate_nrr_omega.m`). Any string other than `'student'` is treated as Gaussian (`estimate_nrr_omega.m`) | `hars.m`, `harsz.m` | both | Distribution of the simulated standardized shocks in the `omega_hat` estimate. | not set by any driver. |
| `nrr_omega_nu` | `2 * options.tparam` | positive scalar | `hars.m`, `harsz.m` | both | Degrees of freedom of the Student-t law used in the `omega_hat` simulation. | not set by any driver. |

The three `nrr_omega_*` options are assembled into an internal struct `nrr_omega_opts` (`hars.m`, `harsz.m`) that is passed to `check_sr_pass_cached.m` and on to `estimate_nrr_omega_fast.m`. The fields `fast`, `validate`, `validate_K`, `product_form`, `hoist_xi`, and `moment_check` of that internal struct are not fields of `options`.

Sub-fields of `NarrativeRestrictions` (entry i of each vector describes anchor i).

| Sub-field | Default | Type | Read in | Meaning |
|---|---|---|---|---|
| `time_index` | required. When the field is absent or empty the regime lookup at `harsz.m` is skipped | vector of row indices into the residual matrix, that is the row of `data.y` minus `lags` | `hars.m`, `harsz.m`, `check_sr_pass_cached.m` | Dates of the anchors. The sampler adds `options.lags` back to find the regime of each date. |
| `shock_index` | required | vector of shock ids | `hars.m`, `harsz.m`, `check_sr_pass_cached.m`, `SignRestrictionCheck.m` | Labelled shock restricted at each date. |
| `sign` | required | vector with entries +1, -1, 0 | `check_sr_pass_cached.m` | Required sign of the realized labelled shock at the date. 0 imposes no sign. |
| `variable_index` | `1` in `check_sr_pass_cached.m` | vector of variable ids | `check_sr_pass_cached.m`, `SignRestrictionCheck.m` | Variable whose historical decomposition is used in the dominance check. |
| `type` | `'dominance'` | cell array of strings: `'sign'`, `'dominance_max'`, `'dominance_min'`, or anything else (treated as `'dominance'`) | `check_sr_pass_cached.m` | `'sign'` skips the dominance check. `'dominance_max'` requires the contribution of the focal shock to exceed the largest other contribution in absolute value. `'dominance_min'` requires it to be smaller than the smallest other contribution. The default case requires it to exceed the sum of all other absolute contributions. |
| `regime_indices`, `regime_index` | written by the sampler | vector, scalar | written at `hars.m`, `harsz.m`. Read at `check_sr_pass_cached.m` | Regime of each anchor date. The sampler fills these fields from `time_index`, so the caller does not set them. |

Driver values of `NarrativeRestrictions`. mon `time_index = 166`, `shock_index = 1`, `sign = 1`, `variable_index = 6`, `type = {'dominance'}` (`run_monetary.m`). fis `time_index` = six dates minus `lags`, `shock_index` = revenue or spending column, `sign = [-1 +1 +1 +1 +1 +1]`, `variable_index = [1 1 2 2 2 2]`, `type` all `'dominance'` (`run_fiscal.m`). The inactive oil block would set `sign = 0` and `type = {'dominance_min'}` (`run_oil.m`).

Sub-fields of `NarrativeRegimeRestrictions(k)`.

| Sub-field | Default | Type | Read in | Meaning |
|---|---|---|---|---|
| `window` | required | vector of row indices into the residual matrix | `hars.m`, `harsz.m`, `check_narrative_regime_avg.m` | Dates of the window. An empty window, or one outside the sample, makes the check fail. |
| `regime_idx` | required | integer in 1..number of regimes | `harsz.m`, `check_narrative_regime_avg.m` | Regime whose Lambda and IRFs are used for the window. |
| `shock_idx` | field required by the start-up printout (`hars.m`, `harsz.m`). Empty switches the focal checks off | shock id | `hars.m`, `harsz.m`, `check_narrative_regime_avg.m` | Focal labelled shock. |
| `sign` | required when `shock_idx` is set | +1 or -1 | `check_narrative_regime_avg.m` | Required sign of the window mean of the focal shock. The check is `sign * mean > 0`, so 0 always fails. |
| `variable_idx` | required for any dominance check (error at `check_narrative_regime_avg.m`). The field must exist for the printout | variable id | `check_narrative_regime_avg.m` | Variable whose historical decomposition is compared. |
| `dominance` | field required by the printout. Empty or `'off'` switches focal dominance off | `'on'`, `'max'`, `'min'`, `'strict'`, `'off'` | `check_narrative_regime_avg.m` | `'on'` and `'max'` require the focal contribution to exceed every other one. `'min'` requires it to be below every other one. `'strict'` requires it to exceed the sum of the others. |
| `dom_agg` | `'mean_abs'` | `'mean_abs'`, `'cum_abs'`, `'interval_hd'` | `check_narrative_regime_avg.m` | How contributions are aggregated over the window. `'mean_abs'` is the mean of monthly absolute impact contributions. `'cum_abs'` is the absolute value of their sum. `'interval_hd'` is the absolute cumulative contribution to the window endpoint, including delayed effects up to the cached horizon. |
| `other_idx` | `setdiff(1:ny, shock_idx)` | vector of shock ids | `check_narrative_regime_avg.m` | Shocks the focal shock is compared with. |
| `anti_shock_idx` | none | shock id | `check_narrative_regime_avg.m` | Shock that must contribute the least. It is active only with `anti_dominance = 'on'`. |
| `anti_dominance` | off | `'on'` | `check_narrative_regime_avg.m` | Switches the anti-dominance check on. |
| `anti_other_idx` | `setdiff(1:ny, anti_shock_idx)` | vector of shock ids | `check_narrative_regime_avg.m` | Comparison set for the anti-dominance check. |

Driver values, ex02 only. Windows 2020:M4 to 2021:M2 and 2023:M9 to 2024:M7, `regime_idx` 3 and 4, `shock_idx` 4 and 5, `variable_idx` 3 for both, `sign` +1 and -1, `other_idx` the other labelled shocks among 1:5, `dominance = 'strict'`, `dom_agg = 'interval_hd'`, and the three `anti_*` fields empty.

## 8. HARS / ESS sampler settings

| Name | Default | Type / allowed values | Read in | Used by | Meaning | Drivers |
|---|---|---|---|---|---|---|
| `ess_store_rounds` | value of `ess_premix_rounds` | positive integer | `hars.m`, `harsz.m` | both | N_ess. Number of elliptical slice updates of the rotation run after each structural move, all of which are stored as draws under the same structural parameters. | all four: 5. |
| `ess_premix_rounds` | `3` | nonnegative integer | `hars.m`, `harsz.m` | both | Number of elliptical slice updates of the rotation run before each structural move, in burn-in and in sampling. | ex02 5. mon, oil, fis 3. |
| `theta_redraw_budget` | `100` | nonnegative integer | `hars.m`, `harsz.m` | both | N_theta. Maximum number of joint redraws of (A0, Lambda, A+) tried when a completed structural draw leaves the current rotation inadmissible. 0 sends every such event to the global rotation repair. | ex02 300. mon, oil, fis `NTHETA`, which the driver sets to 100 when the caller has not defined it. |
| `ess_repair_budget` | `2000` | positive integer | `hars.m`, `harsz.m` | both | N_Q. Number of independent rotation draws tried in the global repair. In `hars` and in `harsz` without zeros the draws are Haar (`find_admissible_Q`, with the batch size equal to the same number). In `harsz` with zeros they come from `find_admissible_Q_zero`. It is also the inner budget of the columnwise searches in the seed-rescue block (`harsz.m`). | ex02 200. mon, oil, fis 2000. |
| `enforce_stability` | `true` | logical | `hars.m`, `harsz.m` | both | True makes `SignRestrictionCheck` reject a draw at the storage stage when any companion eigenvalue has modulus 1 or more. False keeps such draws. | all four: `false`. |
| `upstream_stability_check` | `false` | logical | `hars.m`, `harsz.m` | both | True makes `build_theta_cache` mark a structural draw unstable when any companion eigenvalue has modulus at or above `1 - upstream_stability_tol`, which the sampler counts as a sampler failure. | all four: `false`. |
| `upstream_stability_tol` | `1e-6` | nonnegative scalar | `hars.m`, `harsz.m` | both | Buffer of the upstream stability check. 0 reproduces the rule of `SignRestrictionCheck`. | not set by any driver. |
| `seed_rescue_rounds` | `25` | positive integer | `hars.m`, `harsz.m` | both | Number of rescue rounds tried before the main loop when the rotation left by burn-in is not admissible at the first sampling state. Initialization only. | ex02 25. fis 25. mon, oil not set. |
| `B_redraw_budget` | `100` | positive integer | `hars.m`, `harsz.m` | both | Maximum number of A+ redraws in each seed-rescue round. Initialization only. | not set by any driver. |
| `use_local_pert_fallback` | `false` | logical | `hars.m`, `harsz.m` | both | True enables a fallback in the search for the first rotation and in the seed rescue: the structural parameters are perturbed locally and the columnwise search is repeated. | ex02 `false`. mon `true`. oil `false`. fis `true`. |
| `cw_outer` | `2000` | positive integer | `hars.m`, `harsz.m` | both | Outer budget of the primary column-by-column search for the first rotation. | not set by any driver. |
| `cw_inner` | `1000` | positive integer | `hars.m`, `harsz.m` | both | Inner budget of the same search. | not set by any driver. |
| `haar_attempts` | `20000` | positive integer | `hars.m`, `harsz.m` | both | Total number of Haar candidates in the brute-force fallback for the first rotation (`find_admissible_Q` argument `max_tries`). In `harsz` this fallback is skipped when zeros are present (`harsz.m`). | not set by any driver. |
| `haar_budget` | `2000` | positive integer | `hars.m`, `harsz.m` | both | Number of candidates per `parfor` batch in the same fallback (`find_admissible_Q` argument `batch_size`). | not set by any driver. |
| `local_pert_scale` | `0.1` | positive scalar | `hars.m`, `harsz.m` | both | The local perturbation is drawn from `N(0, local_pert_scale * Sigma)`. Used only with `use_local_pert_fallback = true`. | not set by any driver. |
| `local_pert_tries` | `500` | positive integer | `hars.m`, `harsz.m` | both | Number of local perturbations tried. | not set by any driver. |
| `cw_outer_local` | `5000` | positive integer | `hars.m`, `harsz.m` | both | Outer budget of the columnwise search at each perturbed point. | not set by any driver. |
| `cw_inner_local` | `10000` | positive integer | `hars.m`, `harsz.m` | both | Inner budget of the columnwise search at each perturbed point. | not set by any driver. |
| `max_theta_retries` | `2000` | positive integer | `hars.m`, `harsz.m` | both | Number of alternative starting values, drawn from `N(xh, H)` around the posterior mode, tried before the run stops with "Could not find an admissible Q". | not set by any driver. |
| `report_every` | `10` | positive integer | `hars.m`, `harsz.m` | both | Progress lines in the initialization fallbacks are printed every this many iterations. | not set by any driver. |

## 9. Output and diagnostics

| Name | Default | Type / allowed values | Read in | Used by | Meaning | Drivers |
|---|---|---|---|---|---|---|
| `report_stability_diag` | `~enforce_stability` | logical | `hars.m`, `harsz.m` | both | True computes the largest companion eigenvalue modulus of every stored draw and writes `output.draw_max_root`, `output.plot_keep_mask`, and `output.trim_threshold`. It does not change the draws. | not set by any driver. |
| `trim_threshold` | `1.02` | positive scalar | `hars.m`, `harsz.m` | both | `output.plot_keep_mask` is true for draws whose largest root is at or below this value. | not set by any driver. |
| `bf_robust_stride` | `1` | positive integer | `hars.m`, `harsz.m` | both | The rank check of `check_bf_identification` is repeated on every `bf_robust_stride`-th stored draw. | not set by any driver. |
| `fevd_enable` | `false` | logical | `hars.m`, `harsz.m` | both | True computes the forecast error variance decomposition by regime and writes the `output.fevd_*` fields. It needs restrictions and more than one regime. | ex02 `false`. others not set. |
| `fevd_QE_regimes` | `[1 3]` | vector of regime indices | `hars.m`, `harsz.m` | both | Regimes averaged into `output.fevd_QE`, after intersection with the regimes that exist. | not set by any driver. |
| `fevd_QT_regimes` | `[2 4]` | vector of regime indices | `hars.m`, `harsz.m` | both | Regimes averaged into `output.fevd_QT`. | not set by any driver. |
| `fevd_bs_shocks` | `[4 5]` | vector of shock ids | `hars.m`, `harsz.m` | both | Shocks for which the console table of FEVD shares is printed. It does not change any stored field. | not set by any driver. |

## 10. Other

These fields are read by helpers that receive `options` from `eval_posterior`, which both samplers call with the full struct (for example `hars.m`, `harsz.m`). No driver sets any of them.

### 10a. Threshold switches reachable through `eval_posterior`

`hars.m` and `harsz.m` contain no code specific to these switches. In particular they build the parameter vector as A0 elements followed by Lambda elements (`harsz.m`), without the extra `delta` block that `bvar_posterior_thetero.m` expects. Whether a run with these switches on works is unclear from code.

| Name | Default | Type / allowed values | Read in | Used by | Meaning |
|---|---|---|---|---|---|
| `thetero` | false | logical | `eval_posterior.m` | both | True routes the posterior evaluation to `bvar_posterior_thetero.m`, which adds a shift `delta_i * d_t` to the structural log variances. |
| `d_t` | required when `thetero` is true | nobs x 1 vector of 0 and 1 | `bvar_posterior_thetero.m` | both | Pre-determined stress indicator that multiplies the variance shift. |
| `delta_shocks` | required when `thetero` is true | vector of shock ids | `bvar_posterior_thetero.m` | both | Shocks that carry a free variance shift. |
| `delta_prior_var` | required when `thetero` is true | positive scalar | `bvar_posterior_thetero.m` | both | Variance of the Normal(0, c) prior on each shift. |
| `tthr` | false | logical | `eval_posterior.m`, `bvar_posterior.m`, `bvar_posterior_thetero.m` | both | True passes `d_inf` and `sd_lags` to `rfvar3`, which interacts the listed lags with the indicator. |
| `d_inf` | required when `tthr` is true | vector of 0 and 1 | `bvar_posterior.m`, `bvar_posterior_thetero.m` | both | State indicator that interacts with the state-dependent lags. |
| `sd_lags` | required when `tthr` is true | vector of lag numbers | `bvar_posterior.m`, `bvar_posterior_thetero.m` | both | Lags whose coefficients differ across the two states. |
| `gamma_prior_sd` | none | positive scalar | `bvar_posterior.m`, `bvar_posterior_thetero.m` | both | Standard deviation of a Normal shrinkage prior on the state-dependent coefficient increments. Empty or infinite means no shrinkage. |
| `sd_eqs` | none (all equations) | vector of equation indices | `bvar_posterior.m` | both | Equations allowed to have state-dependent coefficients. It is read only in `bvar_posterior.m`, not in `bvar_posterior_thetero.m`. |

### 10b. Fields read only on the `tvA` path, which `hars.m` and `harsz.m` refuse

`eval_posterior` calls `bvar_posterior_tvA.m` only when its `tvA` argument is true. Both samplers stop with an error when `options.tvA` is true (`hars.m`, `harsz.m`) and pass `tvA = 0`. No run of `hars` or `harsz` reads these fields.

| Name | Default | Type / allowed values | Read in | Used by | Meaning |
|---|---|---|---|---|---|
| `lcA0_tv` | none | ny x ny logical | `resolveTvMask.m` | neither | Explicit mask of the free A0 elements that differ across regimes. |
| `tvA_rows` | none | vector of row indices | `resolveTvMask.m` | neither | Rows of A0 whose free elements differ across regimes. |
| `tvA_fix_diagonal` | `true` | logical | `resolveTvMask.m` | neither | True keeps the diagonal of A0 common across regimes. |
| `tvA_shrink_tau` | `Inf` | positive scalar | `bvar_posterior_tvA.m` | neither | Strength of the shrinkage of regime-specific A0 elements toward their cross-regime mean. `Inf` switches it off. |
| `tvA_shrink_scale` | `sqrt(diag(Sstar))`, else ones | scalar or ny x 1 vector | `bvar_posterior_tvA.m` | neither | Row scale that standardizes the deviations in the shrinkage term. |

---

## A. Workspace variables other than `options` and `data`

| Variable | Read at | What the code does |
|---|---|---|
| `impact_signs_byregime` | `hars.m`, `harsz.m` (captured by the `btc` wrapper around `build_theta_cache`) | The sampler defines this variable only when `options.sign_active_regimes` is set (`hars.m`, `harsz.m`). In every other case the name must already exist in the caller's workspace. All four drivers set `impact_signs_byregime = []` before the call (`ex02_harsz_balance_sheet.m`, `run_monetary.m`, `run_oil.m`, `run_fiscal.m`). `[]` means no per-regime impact gate (`build_theta_cache.m`, `check_sr_pass_cached.m`). |
| `output` | `hars.m`, `harsz.m` | The samplers never clear `output`. They assign fields into whatever `output` exists, and they test `isfield(output,'bf_check')` and `isfield(output,'sign_irf_regime')`. `output.sign_irf_regime(:,:,:,:,iReg)` is an indexed assignment (`hars.m`, `harsz.m`). A stale `output` from an earlier run in the same workspace therefore leaks into the new one. Both examples start with `clear`. The three `run_*.m` drivers do not, and `main_*.m` clears the workspace between runs. |
| `init_ntried`, `use_local_pert_fallback`, `xA_all`, `draw_Phi_mat` | `exist(...)` guards at `hars.m`, `harsz.m` | Each guard sits in a branch in which the sampler has already assigned the variable, so a stale caller value is overwritten before it is read. |
| `NTHETA` | `run_monetary.m`, `run_oil.m`, `run_fiscal.m` | Read by the three `run_*.m` drivers, not by the samplers. The drivers set it to 100 when it does not exist and copy it to `options.theta_redraw_budget`. |
| `REGIME_ON` | `run_fiscal.m` | Read by `run_fiscal.m` only. It defaults to 1, and 0 sets `options.regimes = []`. |

Because the samplers are scripts, every internal variable they create (`y`, `lags`, `x0`, `Q_current`, `draw_*_all`, and so on) stays in the caller's workspace after the run.

## B. Fields the replication drivers set that no code reads

Checked with a search of `code/core/` and `code/third_party/` for each name. The two examples set none of them, except `noLmd`, which must exist.

| Field | Set in | Note |
|---|---|---|
| `firstobs` | mon, oil, fis | Used only inside the drivers to compute `options.nobs`. |
| `presample` | mon, oil, fis | `bvar_posterior.m` hard-codes `presample = 0`. |
| `max_compute` | mon, oil, fis | The samplers hard-code `max_compute = 1` (`hars.m`, `harsz.m`). |
| `dummy` | mon, oil, fis | No read. |
| `nexogenous` | mon, fis | No read. |
| `non_explosive_` | mon, fis | No read. |
| `nunits` | mon, oil, fis | No read. |
| `sign_horizon` | mon, fis | Used only inside the drivers to build the restriction strings. |
| `sign_inneriteration` | mon, oil, fis | No read. `run_oil.m` marks it "not used by ESS". |
| `outeriteration` | mon, oil, fis | No read. |
| `store_all_valid_Q` | mon, oil, fis | No read. |
| `identification_mode`, `include_business_cycle`, `include_monetary`, `bc_col`, `mon_col`, `fiscal_cols_active` | fis | No read. Stored for bookkeeping (`run_fiscal.m`). |
| `noLmd` | all four | Read into a local variable that is never used (`hars.m`, `harsz.m`). The field must exist, but its value has no effect. |
