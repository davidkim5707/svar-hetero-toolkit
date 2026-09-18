# svar-hetero-toolkit

[![tests](https://github.com/davidkim5707/svar-hetero-toolkit/actions/workflows/tests.yml/badge.svg)](https://github.com/davidkim5707/svar-hetero-toolkit/actions/workflows/tests.yml)

MATLAB code for Bayesian structural VARs with heteroskedastic shocks, identified with sign, zero, and narrative restrictions.

| Sampler | Restrictions | Paper |
|---|---|---|
| **HARS**, `code/core/hars.m` | sign and narrative | Dawis Kim and [Tao Zha](http://www.tzha.net/), "Sharpening Economic Interpretation with HARS", [NBER Working Paper w35483](https://www.nber.org/papers/w35483) |
| **HARS-Z**, `code/core/harsz.m` | sign, narrative, and zero | Dawis Kim, "How the Financing of Balance-Sheet Policy Shapes Its Effects", [SSRN](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7214338) |

Both samplers keep the heteroskedastic likelihood, draw the structural parameters by Metropolis-Hastings, and update the rotation by elliptical slice sampling. HARS-Z builds every rotation on the set that satisfies the zero restrictions, and the zeros hold exactly along the chain. With no zero restriction set, `harsz.m` runs the HARS path. The papers describe the models and the theory. This repository documents the code.

## Quick start

```matlab
run setup.m                          % adds code/ to the path
run examples/ex01_hars_oil.m         % HARS on a three-variable oil model, about 2 minutes
run tests/run_tests.m                % self-tests of the zero-restriction routines
```

Each example saves its draws and a figure of impulse responses under `examples/output/`.

## Repository layout

```
code/                 the library
  core/               samplers and their routines
  third_party/        code from external sources
  python/             post-processing in Python
examples/             one script per model, written to be read and adapted
replications/         one folder per paper, written to reproduce its figures and tables
tests/                self-tests
docs/                 reference for the options and output structs
setup.m               adds code/ to the MATLAB path
CHANGELOG.md          what each release contains
```

Every folder has a `README.md` of its own. On every push, GitHub Actions runs the self-tests and `examples/ex01_hars_oil.m` on a Linux runner.

| To do this | Go to |
|---|---|
| Learn how a sampler is called, or adapt one to another model | [`examples/`](examples/README.md) |
| Reproduce the figures and tables of a paper | [`replications/`](replications/README.md) |
| Find a routine, or see how the library is organized | [`code/`](code/README.md) |
| Look up a field of `options` or `output` | [`docs/options.md`](docs/options.md), [`docs/output.md`](docs/output.md) |

### Examples and replications

An example and a replication can estimate the same model. They differ in purpose.

| | Example | Replication |
|---|---|---|
| Purpose | show how to call a sampler and how to state restrictions | reproduce a paper |
| Form | one script with every setting in view | a driver per application, a master driver, and the builders of figures and tables |
| Settings | the paper's, with the draw count at the top of the file | the paper's, fixed |
| Output | the draws and one figure of impulse responses | every figure and table by its number in the paper |
| Comparison with other samplers | none | ships the comparison outputs |

## Requirements

MATLAB R2020a or newer, for `exportgraphics`, with the Statistics and Machine Learning Toolbox. The Parallel Computing Toolbox is optional. The `parfor` loops use a pool when one is available. The posterior mode search uses `csminwel`, which ships under `code/third_party/`, and needs no Optimization Toolbox.

`code/python/hars_z_weights.py` needs Python 3 with `numpy`, and `scipy` or `h5py` to load a saved `.mat`.

## Random numbers

No script sets a random seed. The draws depend on the state of MATLAB's random number generator at the start of the run. With the Parallel Computing Toolbox, `find_admissible_Q` also draws its candidate rotations on the pool workers, and the assignment of candidates to workers changes from run to run. Two runs agree on posterior medians and credible bands up to Monte Carlo error and not draw for draw.

## Citation

Cite the paper whose sampler you use. `CITATION.cff` holds the machine-readable record, which GitHub's "Cite this repository" button reads.

> Kim, D. and T. Zha (2026). *Sharpening Economic Interpretation with HARS*. NBER Working Paper w35483.

> Kim, D. (2026). *How the Financing of Balance-Sheet Policy Shapes Its Effects*. SSRN Working Paper 7214338. https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7214338

## License

MIT. See `LICENSE`, with the scope recorded separately in `NOTICE.md`. Third-party code under `code/third_party/` keeps its own license.
