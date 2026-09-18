# Notice

`LICENSE` holds the MIT text and nothing else, which lets GitHub and other tooling identify it as MIT. This file records what the license covers.

## Scope

The license covers the code in this repository written by the copyright holders. That is `setup.m`, everything under `core/`, `examples/`, `python/`, and `tests/`, and the drivers under `replications/`. Citing the paper is expected scholarly practice and not a licensing condition.

## Not covered

Third-party code keeps its original license and authorship. Every file under `third_party/` comes from an external source. The files under `third_party/sims_bpss/` keep the original attribution in their headers. BPSS refers to the replication package of Brunnermeier, Palia, Sastry, and Sims (2021).

| File | Source | Use |
|---|---|---|
| `csminwel.m` | Christopher Sims optimizer, with Dynare Team additions | verbatim |
| `csminit.m` | Christopher Sims optimizer, with Dynare Team additions | verbatim |
| `bfgsi.m` | Christopher Sims optimizer | verbatim |
| `numgrad.m` | Christopher Sims optimizer | lightly adapted, gradient step changed |
| `rfvar3.m` | Christopher Sims VAR tools, via the BPSS replication R code | converted to MATLAB, lightly modified |
| `varprior.m` | Christopher Sims VAR tools, via the BPSS replication R code | converted to MATLAB, lightly modified |
| `impulsdtrf.m` | Christopher Sims, via the BPSS replication R code | converted to MATLAB, lightly modified |
| `linreg.m` | BPSS replication R code | converted to MATLAB, lightly modified |
| `drawt.m` | BPSS replication R code | converted to MATLAB, lightly modified |
| `drawdelta.m` | BPSS replication R code | converted to MATLAB, lightly modified, with an explicit Gaussian case |

The six files under `third_party/arrw2018/` come from the replication package of Arias, Rubio-Ramirez, and Waggoner (2018), "Inference Based on Structural Vector Autoregressions Identified with Sign and Zero Restrictions: Theory and Applications", *Econometrica*. They are `LogVolumeElement.m`, `SpheresToQ.m`, `SpheresRestriction.m`, `perp.m`, `NumericalDerivative.m`, and `LogAbsDet.m`. `core/zero/log_volume_element_zero.m` calls them for the volume element.

The file `core/posterior/bvar_posterior.m` is a MATLAB reimplementation of R code from the same BPSS replication package. It keeps that attribution in its header.

The license does not cover the datasets under `replications/*/data/` and `examples/data/` or the comparison outputs under `replications/*/baselines/`. The papers give their sources.
