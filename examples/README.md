# Examples

One script per model. Each script states the whole model in one place, calls a sampler, saves the draws, and draws one figure of impulse responses. Read an example from top to bottom to see what a sampler needs, and copy the nearest one to start a model of your own.

| Script | Sampler | Model | Restrictions | Runtime |
|---|---|---|---|---|
| `ex01_hars_oil.m` | `hars` | Three-variable monthly model of the oil market, 1971:M1 to 2015:M12, 24 lags, three volatility regimes | 9 sign restrictions and 4 elasticity bounds at impact | about 2 minutes for 1,000 recorded draws |
| `ex02_harsz_balance_sheet.m` | `harsz` | Eight-variable monthly U.S. model, 2008:M11 to 2026:M6, 6 lags, four volatility regimes, Minnesota prior | 16 sign restrictions on the first six months, 7 zero restrictions at impact, 2 narrative restrictions on windows | about 35 minutes for 300 recorded draws, about 1 hour 40 minutes for 10,000 |

Runtimes are from the author's desktop with an eight-worker pool. Most of the time in `ex02` goes to the 1,000 burn-in sweeps, which do not shrink with the number of recorded draws.

```matlab
run examples/ex01_hars_oil.m
```

Each script calls `setup.m` itself and locates its data from its own path. The working directory does not matter.

## What an example writes

Everything goes to `examples/output/`, which the repository ships empty.

| File | Contents |
|---|---|
| `<script>.mat` | `output`, `options`, and `data`. `docs/output.md` documents `output`. |
| `<script>.pdf` | Posterior median, 68 percent band (shaded), and 90 percent band (dashed) of the impulse responses, drawn by `plot_irf_bands`. |

`ex01` plots every variable against every shock. `ex02` plots every variable against the two balance-sheet shocks. The figure plots `output.sign_irfs_temp` as stored, the responses under the average transmission in the units of the data. The papers rescale and reorient some responses, and the examples do not.

## The layout of a script

Both scripts follow the same order, and each block sets fields of `options`. `docs/options.md` gives the default, the type, and the meaning of every field.

1. Constants at the top, `NDRAW` and `NBURN`, then the paths.
2. Data. The sampler needs `y`, the variable names, and a `dates` field in `datetime` format.
3. VAR and sampler size.
4. Structural model, which covers the prior on `A0`, the volatility regimes, and the error distribution.
5. Prior on the reduced form.
6. Sign restrictions, as strings, and the matrix of impact signs.
7. Zero restrictions, in `ex02` only.
8. Narrative restrictions, in `ex02` only.
9. Sampler settings.
10. The call, the save, and the figure.

## Adapting an example

- A sign restriction is a string such as `'y(2,1:6,1) < 0'`, the response of variable 2 to shock 1 over horizon indices 1 to 6. Horizon index 1 is impact.
- A volatility break is the last month of the regime that ends, not the first month of the regime that starts.
- The three structs `SignRestrictions`, `ZeroRestrictions`, and `NarrativeRegimeRestrictions` and the matrix `impact_signs` must agree with one another.
- Keep `impact_signs_byregime = [];` and the opening `clear`. `code/README.md` explains both.
- In `ex02`, keep `NBURN = 1000`. After a short burn-in the sampler can fail to find an admissible state to start sampling from.
- A model with tight restrictions can stop at initialization with the message that no admissible rotation was found. Section 8 of `docs/options.md` lists the budgets of that search.

## Data

| File | Contents | Used by |
|---|---|---|
| `data/oil_monthly.mat` | `data` (540 by 3), `dates`, `varNames`. Oil production growth, economic activity index, real oil price. The same file as `replications/kim_zha2026_hars/data/Kilian_Data_Updated.mat`. | `ex01` |
| `data/balance_sheet_monthly.csv` | 212 monthly rows and 8 series. One-year Treasury yield, log SOMA holdings, log of reserves plus ON RRP, ON RRP share, log real GDP, log GDP deflator, term spread, excess bond premium. | `ex02` |

The papers give the sources of the series.
