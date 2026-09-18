# Code

The library. Nothing in this folder loads data or draws a figure for a paper. The examples and the replications call it.

| Folder | Contents |
|---|---|
| `core/` | The samplers `hars.m` and `harsz.m` and the routines they call, written for this project. |
| `third_party/` | Code from external sources. `NOTICE.md` at the repository root lists each file and its origin. |
| `python/` | Post-processing in Python. |

`setup.m` at the repository root adds `core/` and `third_party/` with their subfolders to the MATLAB path.

## How a sampler is called

`hars.m` and `harsz.m` run as scripts, not as functions. The caller defines three things in the workspace and then writes the name of the sampler on a line of its own.

| Workspace variable | Contents |
|---|---|
| `data` | A struct with `y` (T by n observations), `varnames`, and `filtered_data` (a table or struct with a `dates` field in `datetime` format). |
| `options` | A struct with the VAR, prior, restriction, and sampler settings. `docs/options.md` documents every field. |
| `impact_signs_byregime` | Set it to `[]`. The samplers read this variable, and an undefined value stops the run. |

The sampler leaves its results in the struct `output`. `docs/output.md` documents every field. The impulse-response draws are in `output.sign_irfs_temp`, an array of size variables by horizons by shocks by draws, with impact in the first horizon. The draw weights are in `output.draw_weight` after `hars` and in `output.draw_weight_final` after `harsz`.

A sampler reads and writes many workspace variables while it runs. Start a driver with `clear`. Without it, values from an earlier run carry over, and `docs/options.md` describes how a stale `output` enters the new results.

## `core/`

| Subfolder | What it is for |
|---|---|
| `core/rotation/` | Drawing the rotation `Q`. `gaussian_to_Q` maps a Gaussian matrix to a rotation by QR. `ess_draw_Q_columnwise` updates `Q` column by column by elliptical slice sampling. `find_admissible_Q` and `find_admissible_Q_columnwise` search for a rotation that satisfies the restrictions. |
| `core/zero/` | The same steps under zero restrictions. `build_zero_constraint` turns `options.ZeroRestrictions` into linear constraints on the columns of `Q`. `gaussian_to_Q_zero`, `draw_Q_zero_columnwise`, `ess_draw_Q_zero_columnwise`, and `find_admissible_Q_zero` build and update rotations that satisfy the zeros exactly. `log_volume_element_zero` computes the volume element of the construction. Each routine delegates to its counterpart in `core/rotation/` when no zero is set. |
| `core/restrictions/` | Checking restrictions. `parse_sign_restrictions`, `SignRestrictionCheck`, `check_sr_pass_cached`, and the two `checkrestrictions_per_shock` routines handle sign restrictions. `check_narrative_regime_avg`, `estimate_nrr_omega`, and `estimate_nrr_omega_fast` handle narrative restrictions. `check_bf_identification` checks the local identification condition of Bacchiocchi and Fanelli (2015). |
| `core/posterior/` | Posterior evaluation and reduced-form tools. `eval_posterior` evaluates the log posterior through `bvar_posterior`. It routes to `bvar_posterior_thetero` when `options.thetero` is set and to `bvar_posterior_tvA` when `A0` varies across regimes. `hars` and `harsz` keep `A0` common across regimes and use the first route. `build_theta_cache`, `compute_Sstar`, `unpackA0_tvA`, and `resolveTvMask` support them. |
| `core/util/` | Shared helpers. `get_opt` reads an option with a default. `getHDs_fast` computes historical decompositions. `wpercentile` computes weighted quantiles, `mess_vfj` the multivariate effective sample size of Vats, Flegal, and Jones (2019), and `plot_irf_bands` draws posterior medians with 68 and 90 percent bands. `save_slim` saves the subset of `output` that the replication figures read. |

## `third_party/`

| Subfolder | Source |
|---|---|
| `third_party/sims_bpss/` | The optimizer of Christopher Sims (`csminwel`, `csminit`, `bfgsi`, `numgrad`) and VAR tools converted to MATLAB from the R code of the replication package of Brunnermeier, Palia, Sastry, and Sims (2021). |
| `third_party/arrw2018/` | The volume-element routines of the replication package of Arias, Rubio-Ramirez, and Waggoner (2018). |

## `python/`

`hars_z_weights.py` rebuilds the HARS-Z estimator weights from a saved `output` struct and checks them against the weights the sampler stored. It needs `numpy`, and `scipy` or `h5py` to load the file. `ex01` saves in the default format, which `scipy.io.loadmat` reads. `ex02` saves with `-v7.3`, which `h5py` reads.

```python
import h5py
from hars_z_weights import estimator_weights

with h5py.File("examples/output/ex02_harsz_balance_sheet.mat", "r") as f:
    weights = estimator_weights(f["output"])
```
