# Changelog

## 1.0.0

First public release.

- `code/core/hars.m`, the HARS sampler of Kim and Zha (2026), for sign and narrative restrictions under heteroskedastic shocks.
- `code/core/harsz.m`, the HARS-Z sampler of Kim (2026), which adds exact zero restrictions.
- `examples/`, one script per sampler, each with its data and a figure of impulse responses.
- `replications/kim_zha2026_hars/`, the replication package of Kim and Zha (2026).
- `docs/options.md` and `docs/output.md`, a reference for every field of the `options` and `output` structs.
- `tests/`, self-tests of the zero-restriction routines.
