# Validation scripts

Run from a scratch directory; each writes beneath the paths it is given.

- `regression_run.py <dir>` and `compare_runs.py <ref_dir> <new_dir>`:
  seeded, shortened runs of the mass- and SMA-function pipelines (using the
  mass repo's defaults). Run once on `occurrence` master and once on
  `ecc_added`; every `.npy`, `.npz`, and `.csv` output must match.
- `mock_validation.py <ecc_run_dir> <mock_dir> <sigma_e> <seed>`: draws a
  mock population (e from Beta(1.5, 4), 0.30 companions per star in the
  region), detects it with each host's eccentricity-sliced completeness,
  gives detected companions Gaussian posteriors of width `sigma_e` in e,
  and fits the eccentricity distribution. Reuses the Tier 1 and Tier 2
  maps of an existing eccentricity run (`<ecc_run_dir>/mtrue/...`).
- `mock_figure.py`: recovered Beta parameters of the mocks vs. the truth.
- `diagnostics.py <ecc_run_dir> <figure_dir> [<reference_run_dir>]`: the
  sliced average completeness maps, completeness vs. e, and injection
  counts per slice.
- `companion_completeness.py <tier2> <companion> [<results_dir>]`: for one
  companion of an eccentricity run (e.g. `highMstar 120066_0`), the
  completeness its draws get in each e slice from the host map, the
  subsample average, and the all-star average; how many stars and host
  injections support the average there; and a figure of those maps around
  the draws. Use it when a fit reports a companion with zero
  completeness/prior.
