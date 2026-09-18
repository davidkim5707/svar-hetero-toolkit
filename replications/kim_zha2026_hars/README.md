# Replication package for Kim and Zha (2026)

Dawis Kim and Tao Zha, "Sharpening Economic Interpretation with HARS" ([NBER Working Paper w35483](https://www.nber.org/papers/w35483)).

This package replicates the three empirical applications (monetary, oil, and fiscal) with the HARS sampler. One command produces every figure and table the package covers.

## How to run

Open MATLAB, set the working directory to this folder, and run one command.

```matlab
main                        % all three applications end to end
```

To run a single application, use `main_monetary`, `main_oil`, or `main_fiscal`. Each runs standalone.

`main` executes seven HARS runs in total. These are two monetary runs (baseline and ablation), two oil runs (baseline and ablation), and three fiscal runs (homoskedastic, 3-regime baseline, and 3-regime ablation). Every run executes fresh on each call, with no caching. The comparison of hours per 1,000 effective draws holds only when all HARS timings come from the same machine. Expect roughly seven times a single run wall clock.

## Code map

`main.m` is the master driver. It runs the three applications. After they finish, it builds the cross-application outputs Figure 5, Table 7, and Table S7.

`main_monetary.m`, `main_oil.m`, and `main_fiscal.m` are the per-application pipelines. Each runs the HARS baseline and the ablation, and the fiscal pipeline adds the homoskedastic run. Each saves slim results, draws that application's figures, and writes its tables.

`run_monetary.m`, `run_oil.m`, and `run_fiscal.m` hold the model definition for one application. Each loads the data, sets the VAR options, the sign and narrative restrictions, and the sampler settings, and calls `hars`. The caller can set `NTHETA` in every script and `REGIME_ON` in the fiscal script.

The sampler `hars.m` and the routines it calls live in `code/` at the repository root. Every driver adds `code/core` and `code/third_party` to the path with a path relative to this folder, which is why the working directory matters.

`data/` holds the three datasets, one per application.

`baselines/` holds comparison outputs from other samplers.

`output/` holds everything `main` writes. These are the slim results, the figures named by paper number, and the efficiency tables. The repository ships it empty.

## Outputs

`main` produces Figures 1 through 5 together with the appendix figures S1 and S2. It writes the appendix Tables S1 through S7 and the monetary, oil, and fiscal rows of Table 7.

| Driver | Figures | Tables |
|---|---|---|
| `main_monetary` | 1, S1 | S1, S2 |
| `main_oil` | 2, S2 | S3, S4 |
| `main_fiscal` | 3, 4 | S5, S6 |
| `main` | 5 | 7, S7 |

## Baselines

The figures draw their comparison bands and timing rows from four files in `baselines/`.

| File | Role |
|---|---|
| `irfs_ramirez_monetary.mat` | AR (2018) homoskedastic IRFs, the monetary Figure 1 baseline. |
| `irfs_carriero_monetary_3regimes.mat` | CMT monetary IRFs for Figures S1 and 5 and the timing row. |
| `irfs_kilian_oil.mat` | Kilian and Murphy (2012) homoskedastic IRFs, the oil Figure 2 baseline. |
| `irfs_carriero_oil_3regimes.mat` | CMT oil IRFs for Figures S2 and 5 and the timing row. |

`irfs_ramirez_monetary.mat` is 348 MB, over GitHub's 100 MB file limit, and does not ship with the repository. Download it from the [replication folder on Dropbox](https://www.dropbox.com/scl/fo/i2fxistcz7r6slpc58j54/ALYBgJbM-d0D_x5k8_MiIHc?rlkey=uzhibvhrvz41keckbjjmx45bg&st=yfmw6879&dl=0) and place it in `baselines/`. Without it, `main_monetary` still runs and leaves the AR (2018) comparison out of Figure 1 and out of the timing table.

The CMT and AR (2018) timing rows use seconds measured on the authors' hardware rather than this machine's clock. The fiscal application ships no baselines. The CMT fiscal sampler fails on this specification, with one admissible draw in 200,000 candidates. The package rebuilds the fiscal homoskedastic baseline internally.

## Settings

Each run uses 10,000 recorded draws and 1,000 burn-in draws. The ESS settings retain `N_ess` equal to 5 rotations, use a global repair budget `N_Q` of 2,000, and apply three pre-mix rounds. `N_theta` is 100 for each baseline and 0 for each ablation.

The package sets no random seed. The draws depend on the state of MATLAB's random number generator at the start of the run. With the Parallel Computing Toolbox, `find_admissible_Q` also draws its candidate rotations on the pool workers, and the assignment of candidates to workers changes from run to run. Replicator results match the paper's posterior medians and credible bands up to Monte Carlo error rather than draw for draw.

The optimizer `bfgsi` writes a temporary `H.dat` in the working directory, and `csminwel` writes `g1.mat`, `g2.mat`, and `g3.mat`. `.gitignore` excludes all four.
