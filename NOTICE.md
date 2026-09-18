# Notice

`LICENSE` holds the MIT text and nothing else, which lets GitHub and other tooling identify it as MIT. This file records what the license covers.

## Scope

The license covers the code in this repository written by the copyright holders. That is `setup.m`, everything under `core/`, and the drivers under `replications/`. Citing the paper is expected scholarly practice and not a licensing condition.

## Not covered

Third-party code keeps its original license and authorship. Every file under `third_party/` comes from an external source, and its header keeps the original attribution. BPSS refers to the replication package of Brunnermeier, Palia, Sastry, and Sims (2021).

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
| `drawdelta.m` | BPSS replication R code | converted to MATLAB, lightly modified |

The file `core/posterior/bvar_posterior.m` is a MATLAB reimplementation of R code from the same BPSS replication package. It keeps that attribution in its header.

The license does not cover the datasets under `replications/*/data/` or the comparison outputs under `replications/*/baselines/`. The paper gives their sources.
