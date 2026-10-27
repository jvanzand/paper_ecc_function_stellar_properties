# Stellar-properties eccentricity workflow

Preliminary work toward a possible paper on whether the eccentricity
distribution of CLS giant companions varies with stellar mass, metallicity,
and age. It mirrors `paper_mass_function_stellar_properties` and
`paper_sma_function_stellar_properties` (occurrence runs, post-fit products,
and paper-item collection).

That code work happens on the `ecc_added` branch of `occurrence` and follows
[ECCENTRICITY_PLAN.md](ECCENTRICITY_PLAN.md), which lists each step, the
default marginal distributions, and the facts about the injections and
posteriors it relies on.

## Workflow

1. `python run_occurrence.py`: runs every entry of `RUN_CONFIGURATIONS`.
   Each run sets `dimensions`: eccentricity is fit in five bins
   (`ECC_EDGES`), and semi-major axis (0.1-30 AU) and companion mass
   (0.4-50 MJ) are marginalized with the flat (uniform per dex) defaults
   of `occurrence.marginal_distributions`. That assumes a and m are
   independent of e. Models: histogram (`piecewise`), Beta, Rayleigh, and
   a two-Rayleigh mixture (`ecc_mixture`).
2. `python run_post_fit_analysis.py`: LaTeX variables, tables, companion
   plots, stellar histograms, and model CDF comparisons.
3. `python collect_paper_items.py`: gathers the paper's figures and tables
   into `paper_items/`.

Runs need `occurrence` checked out on the `ecc_added` branch beside this
repo. Tests: `python -m pytest tests`.

## Contents

- `ECCENTRICITY_PLAN.md`: the implementation plan.
- `run_occurrence.py`, `run_post_fit_analysis.py`, `collect_paper_items.py`:
  the three workflow steps above; `tests/` covers their configuration
  checks.
- `config_dict.py`: the named stellar subsamples (shared with the other
  paper repos).
- `data/derived/star_catalog.json`: the CLS stellar catalog cache.
- `local_paths.example.json`: template for this machine's untracked
  `local_paths.json` (injection-recovery files and companion posteriors).
