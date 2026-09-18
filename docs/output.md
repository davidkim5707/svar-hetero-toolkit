# `output` reference

Both samplers leave their results in the workspace struct `output`. This file lists every field that `code/core/hars.m` and `code/core/harsz.m` write.

Symbols used for dimensions.

| Symbol | Meaning | Source |
|---|---|---|
| `ny` | number of variables | `options.ny` |
| `H` | `options.irf_horizon`. Index 1 is impact | `harsz.m` |
| `nD` | number of stored draws, equal to `output.collected_draws` | `harsz.m` |
| `nReg` | number of volatility regimes, `numel(options.regimes) + 1` | `harsz.m` |
| `K` | number of stored shocks, `numel(accepted_shock)` | `harsz.m` |
| `np` | number of free structural parameters, `sum(lcA0(:)) + sum(lcLmd(:))` | `harsz.m` |
| `nx` | number of exogenous columns (constant, trend terms, `options.exogenous`) | `bvar_posterior.m` |
| `T` | number of rows of `data.y` | `harsz.m` |

Conditions used in the tables.

- **restricted run** means `use_signrestrictions` is true, which holds when any of `SignRestrictions`, `NarrativeRestrictions`, `NarrativeRegimeRestrictions` is non-empty (`hars.m`), or in `harsz.m` also `ZeroRestrictions` (`harsz.m`).
- **hetero** means `any(lcLmd(:))`, that is at least one free relative variance (`harsz.m`). With `options.regimes` empty there is one regime and no free variance.
- **`use_draw_weights`** differs by sampler. In `hars.m` it is `use_nrr_omega && (has_narr_regime || has_narr)` (`hars.m`). In `harsz.m` it is `has_zero_restr || (use_nrr_omega && (has_narr_regime || has_narr))` (`harsz.m`).

## Start here

| What you need | Field | Notes |
|---|---|---|
| IRF draws under the average transmission | `output.sign_irfs_temp` | `ny x H x K x nD`. Shock k of the third dimension is `output.sign_accepted_shock(k)` in the four drivers (see the caveat in section 1). |
| IRF draws by regime | `output.sign_irf_regime` | `ny x H x K x nD x nReg`. Written only when there is more than one regime. |
| Which shock is in which slice | `output.sign_accepted_shock` | Vector of shock ids, sorted. |
| Draw counts | `output.collected_draws`, `output.total_iterations` | Stored draws and outer loop iterations. |
| Weights | `output.use_draw_weights`, then `output.draw_weight_final` (`harsz.m`) or `output.draw_weight` (`hars.m`) | When `use_draw_weights` is false all summaries are unweighted. When it is true, use `output.wmedian(X, 4)` and `output.wquantile(X, p, 4)` on the IRF arrays. |
| Timing | `output.prof.t_wall` | Wall time of the sampling loop in seconds. Burn-in time is printed to the console and is not stored. The end-to-end time `run_total_sec` of the `run_*.m` drivers is a workspace variable, not a field of `output`. |
| Sampler health | `output.ess_diagnostics` | Counts of MH proposals, acceptances, repairs, and failures. |
| Variable names and regime list | `output.varnames`, `output.lrange` | Copied from `data.varnames`. `lrange = 1:nReg`. |

A plain median band for an unweighted run is `median(output.sign_irfs_temp, 4)` and `quantile(output.sign_irfs_temp, [0.16 0.84], 4)`.

The samplers do not clear `output` before writing. If a struct named `output` from an earlier run is still in the workspace, its fields remain unless the new run overwrites them, and `output.sign_irf_regime(:,:,:,:,iReg) = ...` (`hars.m`, `harsz.m`) assigns into the existing array.

---

## 1. Impulse responses

| Field | Dimensions | Written by | Condition | Meaning |
|---|---|---|---|--- |
| `sign_irfs_temp` | `ny x H x K x nD` (variable, horizon, stored shock, draw) | both | restricted run | IRF of each stored shock under the average transmission: reduced-form dynamics `By`, impact matrix `A0 \ Q'` with Lambda set to 1, multiplied by the orientation `fsign` of the shock (`SignRestrictionCheck.m`). In `harsz.m` with `zero_on_average = true` the zero restrictions hold exactly in this array. |
| `sign_irf_regime` | `ny x H x K x nD x nReg` | both | restricted run and `nReg > 1` | IRF of each stored shock in each regime, rebuilt after sampling with impact matrix `(A0 \ Q') * sqrt(Lambda_r)` when `sign_regime_dependent == 0` and `(A0 \ sqrt(Lambda_r)) * Q'` otherwise, multiplied by the regime orientation in `draw_fsign_cell`. |
| `sign_accepted_shock` | vector of length `K` | both | restricted run | Sorted ids of the stored shocks: the shocks named in `SignRestrictions` (or, without sign strings, the narrative shocks), plus in `harsz.m` the shocks that carry a zero restriction (`harsz.m`). |
| `hetero_irf_temp` | `ny x H x ny x nD` | both | no restrictions and `nReg > 1` | IRF with impact matrix `inv(A0)` and no rotation, for every shock. |

Caveat on the third dimension of `sign_irfs_temp`. Its shock list is computed inside `SignRestrictionCheck.m` from the sign strings and `NarrativeRestrictions.shock_index` only. The list in `sign_accepted_shock`, which indexes `sign_irf_regime`, is computed in the sampler and in `harsz.m` also contains the zero-restricted shocks. The two lists coincide in all four drivers. They can differ when a zero-restricted shock has no sign string, or when sign strings and point narrative restrictions name different shocks.

## 2. Draw counts

| Field | Dimensions | Written by | Condition | Meaning |
|---|---|---|---|--- |
| `collected_draws` | scalar | both | always | Number of stored draws, equal to `options.ndraw` at the end of a complete run. |
| `total_iterations` | scalar | both | always | Number of outer loop iterations. Each iteration runs `nsep` structural sweeps and stores up to `ess_store_rounds` draws. |
| `total_parameter_draws_discarded` | scalar | both | restricted run | Number of storage attempts skipped because `SignRestrictionCheck` rejected the state or returned no rotated parameter vector. |
| `total_redraws` | scalar | both | restricted run | Number of stale-rotation events, equal to `ess_diagnostics.total_repair_attempts`. |
| `avg_redraws_per_iteration` | scalar | both | restricted run | `total_redraws / max(total_iterations, 1)`. |

## 3. Weights

| Field | Dimensions | Written by | Condition | Meaning |
|---|---|---|---|--- |
| `use_draw_weights` | logical scalar | both | restricted run | True when the summaries inside the sampler (variance ratios, FEVD) are weighted. The definition differs between the two samplers, see the top of this file. |
| `draw_weight` | `nD x 1` | both | restricted run | **Differs by sampler.** In `hars.m` it is `1/omega_hat` for each draw, and a vector of ones when the omega weighting is off (`hars.m`). In `harsz.m` it is the raw volume element `v` of the zero construction, equal to 1 when there are no zeros, and it is not divided by `omega_hat` (`harsz.m`). |
| `draw_block` | `nD x 1` | harsz | restricted run | Index of the structural draw (theta block) that each stored draw belongs to. Draws stored in the same storage stage share one value. |
| `draw_weightA` | `nD x 1` | harsz | restricted run | `draw_weight` normalized to sum to 1 within each theta block (`harsz.m`). A block whose weights do not sum to a positive finite number gets uniform weights and a warning. |
| `draw_weight_final` | `nD x 1` | harsz | restricted run | `draw_weightA ./ draw_omega_nrr`. These are the weights used by every weighted summary in `harsz.m`. With `use_nrr_omega = false` they equal `draw_weightA`. |
| `draw_omega_nrr` | `nD x 1` | both | restricted run | Estimate `omega_hat` for each draw. It is a vector of ones when `use_nrr_omega` is false or there are no narrative restrictions. |
| `om_diag` | `nD x 1` cell | both | restricted run | One diagnostic struct per draw from `estimate_nrr_omega_fast.m` (fields such as `group_omegas`, `n_groups`, `floored`, `n_pass`, and for the first draw `xi_std`, `xi_kurtosis`, `n_rows`). The cells are empty when the omega weighting is off. |
| `use_nrr_omega` | logical scalar | both | restricted run | Copy of the resolved `options.use_nrr_omega` (default true). |
| `nrr_omega_opts` | struct | both | restricted run | Settings of the omega simulation: `M`, `shock_law`, `nu`, `moment_check`, and `validate` when the simulation ran. |
| `weight_ess_share` | scalar | both | `use_draw_weights` | `(sum w)^2 / (nD * sum w^2)` for the weights the sampler uses (`draw_weight` in `hars.m`, `draw_weight_final` in `harsz.m`). 1 means uniform weights. |
| `weight_ess_share_global` | scalar | harsz | `use_draw_weights` | The same share for `draw_weight ./ draw_omega_nrr` without the within-block normalization. The code marks it as a reference diagnostic that no estimate uses. |
| `wquantile` | function handle `@(X, p, dim)` | both | `use_draw_weights` | Weighted p-quantile of `X` along dimension `dim`, with the sampler's weights captured in the handle. Use `dim = 4` for both IRF arrays. |
| `wmedian` | function handle `@(X, dim)` | both | `use_draw_weights` | Weighted median, equal to `wquantile` with `p = 0.5`. |
| `has_zero_restr` | logical scalar | harsz | restricted run | True when `options.ZeroRestrictions` is non-empty. |
| `ZeroRestrictions` | struct array | harsz | restricted run with zeros | Copy of `options.ZeroRestrictions`. |

## 4. Parameter draws

Every array has the draw index in its last dimension.

| Field | Dimensions | Written by | Condition | Meaning |
|---|---|---|---|--- |
| `draw_x0` | `(np * nD) x 1` | both | always | Free structural parameters of each draw, A0 elements in the order of `lcA0` followed by Lambda elements in the order of `lcLmd`. The code applies `cell2mat` to an `nD x 1` cell of `np x 1` vectors (`hars.m`, `harsz.m`), so the draws are stacked in one column. `reshape(output.draw_x0, [], output.collected_draws)` gives `np x nD`. The usage note in `check_bf_identification.m` indexes it as `draw_x0(:,end)`, which assumes the reshaped form. |
| `draw_Qdraw` | `ny x ny x nD` | both | restricted run | Rotation `Q` of each draw. The impact matrix uses its transpose, `A0 \ Q'`. |
| `draw_Lambda` | `ny x nReg x nD` (`ny x 1 x nD` when no variance is free) | both | restricted run | Relative shock variances by regime, with the last regime filled from the normalization that each row sums to `nReg` (`SignRestrictionCheck.m`). |
| `draw_By` | `ny x ny x lags x nD` | both | restricted run | Reduced-form lag matrices `By(:,:,l) = A0 \ Aplus(:,:,l)`. |
| `draw_Phi` | `(ny*lags + nx) x ny x nD` | both | always | Coefficient draw `Bdraw` from `rfvar3`. Rows `(l-1)*ny + (1:ny)` hold lag l, and their transpose is the structural lag matrix `Aplus(:,:,l)` (`harsz.m`). The remaining `nx` rows are the exogenous terms. |
| `draw_fsign` | `1 x ny x nD` | both | restricted run | Orientation (+1 or -1) of each shock under the average transmission. |
| `draw_fsign_cell` | `ny x nReg x nD` | both | restricted run | Orientation of each shock in each regime. Entries are NaN when `SignRestrictionCheck` returns no per-regime orientation (`harsz.m`). |
| `draw_udraw` | `(T - lags) x ny x nD` | both | always | Structural residuals `var.udraw` at the coefficient draw, on the weighted scale of the regression, where the weight `wt = exp(0.5*lmdseries)` carries the regime variances and latent scales (`linreg.m`). |
| `draw_udraw_true` | `(T - lags) x ny x nD` | both | always | The same residuals divided by `wt`, which returns them to the unweighted scale (`linreg.m`). The narrative checks standardize these by `sqrt(lambda_r)` and rotate them by `Q` to obtain the labelled shocks (`check_narrative_regime_avg.m`). |
| `draw_dout` | `nobs x ny x nD` | both | hetero | Latent log scales `delta` of the error distribution, drawn by `drawt` when `options.tparam` is set and by `drawdelta` otherwise. |
| `draw_var` | `nD x 1` cell | both | always | In a run without restrictions each cell is a struct with `Bdraw` and `udraw_true` (`harsz.m`). In a restricted run the cells are never filled and stay empty. |

## 5. Posterior values per draw

`eval_posterior` returns the **negative** log posterior as its first output (`bvar_posterior.m`). The other three outputs are log values.

| Field | Dimensions | Written by | Condition | Meaning |
|---|---|---|---|--- |
| `draw_lh` | `nD x 1` | both | always | Negative log posterior at the stored structural parameters, before rotation. |
| `draw_likelihood` | `nD x 1` | both | always | Log likelihood term (`log_dnsty - log_prior` in `bvar_posterior.m`). |
| `draw_A0prior` | `nD x 1` | both | always | Log prior of A0. |
| `draw_lambdaprior` | `nD x 1` | both | always | Log prior of Lambda. |
| `draw_lh_wQ` | `nD x 1` | both | restricted run | Negative log posterior evaluated at the rotated parameter vector, in which `Q*A0` replaces A0 (`SignRestrictionCheck.m`, `harsz.m`). |
| `draw_likelihood_wQ` | `nD x 1` | both | restricted run | Log likelihood term at the rotated parameter vector. |
| `draw_A0prior_wQ` | `nD x 1` | both | restricted run | Log prior of the rotated A0. |
| `draw_lambdaprior_wQ` | `nD x 1` | both | restricted run | Log prior of Lambda at the rotated parameter vector. |

## 6. Model description

| Field | Dimensions | Written by | Condition | Meaning |
|---|---|---|---|--- |
| `lcA0` | `ny x ny` logical | both | always | Mask of the free elements of A0, from `options.A0_restriction`. |
| `lcLmd` | `ny x nReg` logical | both | always | Mask of the free relative variances. The last regime is never free, and rows listed in `options.homoskedastic_shocks` are not free. |
| `lcA0_tv` | `ny x ny` logical, all false | both | always | Mask of regime-specific A0 elements. It is all false because both samplers keep A0 common across regimes. |
| `lrange` | `1 x nReg` | both | always | `1:nReg`. |
| `varnames` | as supplied | both | always | Copy of `data.varnames`. |
| `tvA` | scalar, 0 | both | always | Always 0 in these samplers. |

## 7. Sampler diagnostics (`output.ess_diagnostics`)

All fields are scalars unless noted, and all are written in every run. The counters cover the sampling loop only, not burn-in (`total_B_redraws` covers the seed-rescue block before the loop). In a run without restrictions every counter stays at 0 and the budget fields hold the resolved option values.

| Field | Written by | Meaning |
|---|---|--- |
| `total_mh_proposals` | both | Number of MH proposals of (A0, Lambda). |
| `total_mh_accepted` | both | Number of accepted MH proposals. |
| `total_sr_direct` | both | Number of completed structural draws under which the current rotation stayed admissible, so no repair was needed. |
| `total_repair_attempts` | both | Number of stale-rotation events: completed draws, accepted or rejected by MH, that left the current rotation inadmissible. |
| `total_repair_success` | both | Stale events resolved by the global rotation repair. |
| `total_theta_redraws` | both | Number of joint (A0, Lambda, A+) redraw proposals. |
| `total_theta_repair_success` | both | Stale events resolved by a joint redraw, with the rotation kept. |
| `theta_redraw_budget` | both | Resolved value of `options.theta_redraw_budget`. |
| `stale_q_events` | both | Same value as `total_repair_attempts`. |
| `stale_resolved_structural` | both | Same value as `total_theta_repair_success`. |
| `stale_resolved_fallback` | both | Same value as `total_repair_success`. |
| `stale_unresolved` | both | Stale events resolved by neither route. The chain stays at its current state in these cases. |
| `fallback_share_of_theta` | both | `total_repair_success / max(total_sr_direct + total_repair_success + total_theta_repair_success, 1)`. |
| `total_sampler_failures` | both | Completed draws flagged unstable by `build_theta_cache`, plus stale events for which the global repair also failed. |
| `total_ess_evals` | both | Total number of restriction checks made inside the elliptical slice updates. |
| `total_B_redraws` | both | Number of A+ redraws in the seed-rescue block before the main loop. |
| `B_redraw_budget` | both | Resolved value of `options.B_redraw_budget`. |
| `premix_rounds` | both | Resolved value of `options.ess_premix_rounds`. |
| `ness_store` | both | Resolved value of `options.ess_store_rounds`. |
| `repair_budget` | both | Resolved value of `options.ess_repair_budget`. |
| `enforce_stability` | both | Resolved value of `options.enforce_stability`. |
| `repair_order` | both | String. `'Haar_only'` in `hars.m` and `'independent_only'` in `harsz.m`. |
| `use_nrr_omega` | both | Resolved value of `options.use_nrr_omega`. |
| `has_zero_restr` | harsz | True when zero restrictions are present. |

## 8. Timing profile

| Field | Dimensions | Written by | Condition | Meaning |
|---|---|---|---|--- |
| `prof` | struct | both | always | Timers for the sampling loop. `t_eval`, `t_cache`, `t_ess`, `t_repair`, `t_src`, `t_delta` are seconds spent in `eval_posterior`, `build_theta_cache`, the elliptical slice updates, the restriction check and repair, `SignRestrictionCheck`, and `drawt` or `drawdelta`. The matching `n_*` fields are call counts. `t_wall` is the wall time of the sampling loop and `t_sum` is the sum of the six timers. The timers are filled only in a restricted run. |

## 9. Identification check (Bacchiocchi and Fanelli rank condition)

| Field | Dimensions | Written by | Condition | Meaning |
|---|---|---|---|--- |
| `bf_check` | struct | both | `nReg > 1` and hetero | Rank check of `check_bf_identification.m` at the posterior mode. Fields: `dim_theta`, `n_moments`, `order_ok`, `order_slack`, `rank`, `full_rank`, `sing_values`, `smallest_sv`, `sv_gap_ratio`, `null_dirs` (`check_bf_identification.m`). A warning is printed when `full_rank` is false. |
| `bf_rank_draws` | `nbf x 1` | both | restricted run, `nReg > 1`, hetero, and `bf_check` exists | Rank of the same Jacobian at each checked draw. `nbf` is the number of draws kept by `plot_keep_mask` (all draws when that mask does not exist), thinned by `options.bf_robust_stride`. |
| `bf_svgap_draws` | `nbf x 1` | both | same | `sv_gap_ratio` at each checked draw. |
| `bf_full_rank_share` | scalar | both | same | Share of checked draws whose rank equals `bf_check.dim_theta`. |
| `bf_robust_idx` | `1 x nbf` or `nbf x 1` | both | same | Indices of the checked draws. |

## 10. Stability diagnostic

| Field | Dimensions | Written by | Condition | Meaning |
|---|---|---|---|--- |
| `draw_max_root` | `nD x 1` | both | restricted run and `report_stability_diag` (default `~enforce_stability`) | Largest modulus of the companion eigenvalues of `draw_By` for each draw. |
| `plot_keep_mask` | `nD x 1` logical | both | same | True when `draw_max_root <= trim_threshold`. The stored IRFs are not trimmed. |
| `trim_threshold` | scalar | both | same | Resolved value of `options.trim_threshold` (default 1.02). |

## 11. Variance ratios

Medians and quantiles are taken over draws. They are weighted when `use_draw_weights` is true.

| Field | Dimensions | Written by | Condition | Meaning |
|---|---|---|---|--- |
| `Lambda_median` | `ny x nReg` | both | restricted run and `nReg > 1` | Median of `draw_Lambda`. |
| `omega_median` | `ny x nReg` | both | same | Median of the variance ratio `lambda(i,s) / lambda(i,1)`. |
| `omega_16` | `ny x nReg` | both | same | 16th percentile of the same ratio. |
| `omega_84` | `ny x nReg` | both | same | 84th percentile of the same ratio. |

The 16th and 84th percentiles of Lambda are computed (`harsz.m`) but not stored.

## 12. Forecast error variance decomposition

All fields in this section require `options.fevd_enable = true`, a restricted run, `nReg > 1`, and an existing `output.sign_irf_regime` (`hars.m`, `harsz.m`). `nH` is the number of reported horizons. The last dimension of size 3 holds the median, the 16th percentile, and the 84th percentile over draws, weighted when `use_draw_weights` is true.

| Field | Dimensions | Written by | Condition | Meaning |
|---|---|---|---|--- |
| `fevd_horizons` | `1 x nH` | both | see above | Horizons reported, `[0 1 2 5 10 15 20 25 30 35 40 45 50 55 59]` restricted to values below `irf_horizon`. 0 is impact. |
| `fevd_share_regime` | `ny x nH x K x nReg x 3` | both | see above | Share of each stored shock in the forecast error variance of each variable, by regime. The denominator is the total forecast MSE from all `ny` shocks. |
| `fevd_labeled_total_regime` | `ny x nH x nReg x 3` | both | see above | Sum of the shares over the `K` stored shocks. |
| `fevd_accepted_shock` | vector of length `K` | both | see above | Shock ids of the third dimension of `fevd_share_regime`. |
| `fevd_QE_regimes` | vector | both | see above | `options.fevd_QE_regimes` (default `[1 3]`) intersected with `1:nReg`. |
| `fevd_QT_regimes` | vector | both | see above | `options.fevd_QT_regimes` (default `[2 4]`) intersected with `1:nReg`. |
| `fevd_QE` | `ny x nH x K x 3` | both | see above, and `fevd_QE_regimes` non-empty | Shares averaged over the regimes in `fevd_QE_regimes` for each draw, then summarized over draws. |
| `fevd_QT` | `ny x nH x K x 3` | both | see above, and `fevd_QT_regimes` non-empty | The same for `fevd_QT_regimes`. |

---

## Differences between `hars.m` and `harsz.m`

| Item | `hars.m` | `harsz.m` |
|---|---|---|
| `draw_weight` | `1/omega_hat` (`hars.m`) | raw volume element `v` (`harsz.m`) |
| Weights used in weighted summaries and in `wquantile`, `wmedian` | `draw_weight` (`hars.m`) | `draw_weight_final` (`harsz.m`) |
| `use_draw_weights` | narrative omega weighting only (`hars.m`) | also true whenever zeros are present (`harsz.m`) |
| Fields written only by `harsz.m` | | `draw_block`, `draw_weightA`, `draw_weight_final`, `has_zero_restr`, `ZeroRestrictions`, `weight_ess_share_global`, `ess_diagnostics.has_zero_restr` |
| `ess_diagnostics.repair_order` | `'Haar_only'` | `'independent_only'` |
| `sign_accepted_shock` | shocks from sign strings, else narrative shocks | the same plus zero-restricted shocks (`harsz.m`) |
