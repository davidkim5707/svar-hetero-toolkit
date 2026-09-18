# Replications

Each folder here reproduces the figures and tables of one paper with one command. A replication fixes every setting at the value the paper used and ships the data and the comparison outputs it needs. To learn how to call a sampler or to adapt one to another model, start from `examples/` instead.

| Folder | Paper | Sampler | Command |
|---|---|---|---|
| `kim_zha2026_hars/` | Kim and Zha (2026), "Sharpening Economic Interpretation with HARS" | `hars` | `main` |

Each folder has its own `README.md` with the run order, the list of outputs by figure and table number, the settings, and the runtime.

## Conventions

A replication folder is named `<authors><year>_<sampler>`. It holds the drivers at its top level and three subfolders.

| Subfolder | Contents |
|---|---|
| `data/` | The datasets the drivers read. |
| `baselines/` | Outputs of other samplers that the figures and tables compare against. A file over GitHub's 100 MB limit is attached to a release, and `download_baselines.m` fetches it. |
| `output/` | Everything the drivers write. The repository ships it empty. |

The drivers run from inside their folder and add `code/core` and `code/third_party` to the path themselves.
