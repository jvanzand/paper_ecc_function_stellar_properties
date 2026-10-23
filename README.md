# Stellar-properties eccentricity workflow

Preliminary work toward a possible paper on whether the eccentricity
distribution of CLS giant companions varies with stellar mass, metallicity,
and age. It will mirror `paper_mass_function_stellar_properties` and
`paper_sma_function_stellar_properties` (occurrence runs, post-fit products,
and paper-item collection) once `occurrence` can fit eccentricity.

That code work happens on the `ecc_added` branch of `occurrence` and follows
[ECCENTRICITY_PLAN.md](ECCENTRICITY_PLAN.md), which lists each step, the
default marginal distributions, and the facts about the injections and
posteriors it relies on.

## Contents so far

- `ECCENTRICITY_PLAN.md`: the implementation plan.
- `config_dict.py`: the named stellar subsamples (shared with the other
  paper repos).
- `data/derived/star_catalog.json`: the CLS stellar catalog cache.
- `local_paths.example.json`: template for this machine's untracked
  `local_paths.json` (injection-recovery files and companion posteriors).
