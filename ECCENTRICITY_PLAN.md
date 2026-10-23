# Plan: eccentricity in the occurrence code

Working plan for adding eccentricity to `occurrence` so this repo can measure
the eccentricity distribution of CLS giant companions and test whether it
varies with stellar mass, metallicity, and age. It is the template for the
code work: each step lists the goal, the current state, the design, the files
to change, tests, and what "done" means. Update it as decisions change.

- **Code branch:** every `occurrence` change goes on `ecc_added`, branched
  from `master` (after `parametric_refactor` was merged into it). Merge back
  into `master` only after step 9 passes.
- **This repo** gets its run, post-fit, and collection scripts (mirroring
  `paper_mass_function_stellar_properties` and
  `paper_sma_function_stellar_properties`) once steps 1-5 make an
  eccentricity run possible (step 10).

## Guiding decisions

1. **No higher-dimensional fits.** A run fits one or two dimensions and
   marginalizes over the rest. The fitted and marginalized dimensions are
   chosen per run from semi-major axis (`a`), companion mass (`m`, true or
   M sin i) or mass ratio (`q`), and eccentricity (`e`).
2. **Most general completeness.** Each companion's term uses its own host's
   completeness, built from that host's injections over the fitted
   dimensions (as the mass paper does for `a` and `m`). No pooled
   eccentricity shape is borrowed across hosts. The expected-detection term
   uses the sample-averaged completeness, as now.
3. **No separable joint model.** The occurrence rate is not factored as
   lambda(a, M) f(e). See the caveat below for the assumption marginalizing
   still makes.
4. **Injection sampling.** A marginalized dimension's completeness is averaged
   over the *population* distribution in that dimension; a fitted dimension
   wants even coverage. Injections are reweighted to whichever the run needs
   (step 3), so one set of injections serves every role.

### Caveat: what marginalizing assumes

Marginalizing over a dimension with a fixed distribution assumes the
companions follow that distribution there, independently of the fitted
dimensions, inside the region of interest. For example, an eccentricity fit
that marginalizes over mass assumes the mass distribution does not depend on
eccentricity. Today's code already assumes this for eccentricity: its
completeness maps average over Beta-distributed injections. The assumption is
weakest where known correlations exist (tidal circularization at
a < ~0.1 AU; possibly higher eccentricity for more massive giants), so step 9
checks it by splitting on a marginalized dimension.

**Keep this visible in the code.** Every function that marginalizes (the
injection weights, the marginal completeness maps, and the likelihood terms
for marginalized dimensions) states the independence assumption in its
docstring, and the run option that sets `marginal_distributions` repeats it
in a comment, so a reader cannot miss it.

## Default marginal distributions

Used whenever a dimension is marginalized, unless a run overrides them. Each
is normalized over its range; only its shape matters.

| Dimension | Distribution | Range | Code form |
|---|---|---|---|
| `e` | Beta(0.867, 3.03) (Kipping 2013) | 0 <= e < 1 | `scipy.stats.beta(0.867, 3.03).pdf(e)` |
| `a` | Broken power law, C a^beta [1 - exp(-(a/a0)^gamma)], with C = 350, beta = -0.86, a0 = 3.6 AU, gamma = 1.59 (Fulton et al. 2021), per unit log a | 0.1-30 AU | `likelihood.broken_powerlaw_density([350, log10(3.6), -0.86, 1.59], a)` |
| `m` | Log-normal in log10 M (M_Jup): A exp[-(log10 M - mu)^2 / (2 sigma^2)], with A = 0.09, mu = 0.23, sigma = 0.48 dex (the mass paper's full-sample logG fit), per unit log M | 0.4-50 M_Jup | `likelihood.log_gaussian_density([0.09, 0.23, 0.48], M)` |

Densities for `a` and `m` are per unit log10, matching the code's
occurrence-rate densities; convert to per unit linear coordinate only when a
dimension is fitted on a linear axis. These defaults live in each run's
configuration (`marginal_distributions`, step 1) so a run records what it
assumed.

## Facts about the inputs

- **Injections** (`<star>_recoveries.csv`, about 1000 per star): period and
  semi-amplitude K are drawn log-uniformly (P = 2-10^6 d, K = 0.1-1000 m/s);
  eccentricity from Beta(0.867, 3.03); `inj_au` and `inj_msini` are derived.
  Because log K is linear in log M sin i and log P, the injections are
  uniform in (log a, log M sin i) inside a K window that shifts with period,
  stellar mass, and e. The injection density needed for reweighting (step 3)
  must account for that window's edges.
- **True-mass maps** are made by drawing 10 isotropic inclinations per
  injection (`completeness_utils.recs_msini_converter`); each realization
  keeps the injection's recovery outcome.
- **Companion posteriors** (orvara, `resampled_posteriors_1ksamples/*.h5`,
  1000 draws per system) hold `msec`, `sau`, `esino`, `ecoso` per companion,
  except GL 876, HD 219134, and HD 28185 (10 companions), whose files have
  only `sau` and `msec`. Those need eccentricity posteriors or must be
  excluded from eccentricity fits. orvara samples
  sqrt(e) sin(omega) and sqrt(e) cos(omega) uniformly over the unit disk, so
  its prior is uniform in e, as well as log-uniform in a and M.

## Steps

### 1. Choose fitted and marginalized dimensions per run

- **Goal:** a run states which dimensions it fits (1 or 2), which it
  marginalizes, the bounds of each, and the marginal distributions.
- **Now:** Tier 1 labels (`mtrue`, `msini`, `qtrue`, `qsini`) fix the mass
  coordinate, `a_edges`/`m_edges` fix a 2D grid, and `stack_dim` picks the
  fitted axis; eccentricity is implicitly marginalized over the injection
  Beta distribution.
- **Design:** new run keys, e.g.
  `"dimensions": {"a": {"role": "fit", "edges": [...], "scale": "log"},
  "e": {"role": "fit", "edges": [...], "scale": "linear"},
  "m": {"role": "marginalize", "bounds": (0.4, 50), "kind": "true"}}` and
  `"marginal_distributions": {...}` (defaults above). The existing keys
  become a shorthand that `run.py` translates into this form, so the mass and
  SMA repos keep working unchanged. Validation rejects a dimension listed
  twice, more than two fitted dimensions, or a missing bound.
- **Files:** `run.py` (configuration parsing, Tier labels), the paper repos'
  `run_occurrence.py` (`RUN_DEFAULTS`, validation).
- **Tests:** shorthand and explicit forms resolve to the same configuration;
  invalid role sets are rejected.
- **Done when:** the mass and SMA runs resolve identically to today, and an
  (e fitted; a, m marginalized) run validates.

### 2. Record dimension roles and refuse mismatched reuse

- **Goal:** results from one choice of roles can never be read as another
  (e.g. an (a, e) map feeding a mass function).
- **Design:** write a `dimensions` record (roles, bounds, scales, marginal
  distributions, injection weighting) into Tier 1 maps, Tier 2 sample files,
  `fit_data.npz`, and every chain file. Readers compare it with the run and
  stop with a message naming the difference. Tier 1/2/3 folder names gain a
  short token for the roles (e.g. `fit-a-e_marg-m`) so different choices do
  not share folders.
- **Files:** `main.py`, `completeness_utils.py`, `fit_utils.py`, `mcmc.py`,
  `run.py`, `post_fit_analysis.py` (readers), `convergence.py` (labels).
- **Tests:** a mismatched reuse raises; matching reuse works; old files
  without the record are treated as today's (a, m fitted; e marginalized
  over the injection Beta distribution) and still load.
- **Done when:** every saved product carries its roles and every reader
  checks them.

### 3. Weight injections for each role

- **Goal:** compute completeness under any role assignment from the same
  injections.
- **Design:** each injection (and each inclination realization) gets weight
  w = prod over dimensions of target(x) / injection_density(x), where
  target is the marginal distribution for a marginalized dimension and
  uniform in the chosen scale for a fitted one. Injection densities: uniform
  in log P and log K (with the K-window edges), Beta(0.867, 3.03) in e,
  isotropic in inclination. Completeness in a cell is the weighted fraction
  recovered. Report the effective number of injections per cell,
  (sum w)^2 / sum w^2, and flag cells below a threshold.
- **Files:** new `injection_weights.py`; `completeness_utils.py` (binning
  with weights).
- **Edges of the injection window.** The K window (0.1-1000 m/s) leaves no
  injections at high mass and short period (e.g. 50 M_Jup inside ~0.1 AU)
  and shifts with stellar mass and e, through K proportional to
  M sin i P^(-1/3) (1 - e^2)^(-1/2). Three ways this can go wrong:
  - *Averaging over missing coverage.* If a marginalized dimension's
    completeness is the weighted fraction of the injections that exist, the
    unsampled region simply drops out of the average. That region is mostly
    high-K, where completeness is near 1, so the marginal completeness would
    be biased low. Fill unsampled regions first, with the same physical rule
    as `fill_completeness_nans` (1 above, 0 below the sampled range), then
    average with the target weights over the full range.
  - *Weights that blow up.* The injection density falls toward the window
    edges, so target/density becomes large there and a few injections
    dominate a cell. Cap the weights, and use the effective count to flag
    cells.
  - *True-mass edges.* Inclination draws spread each injection to higher
    true mass, so the injection density in log M_true is the M sin i density
    convolved with the sin i distribution and is not flat near the edges.
    Compute weights with that convolved density, not a constant.
- **Tests:** with uniform weights the maps equal today's; reweighting a
  synthetic set drawn from one distribution reproduces maps from another;
  effective counts match analytic cases; marginalizing over a range that
  extends past the injection window matches the filled analytic result.
- **Done when:** maps can be built for any role assignment, with
  effective-count diagnostics saved alongside them.

### 4. Build completeness over the fitted dimensions

- **Goal:** per-host and sample-averaged completeness maps over the fitted
  dimensions only, marginalized over the rest with step-3 weights.
- **Design:** generalize `single_map_maker`, `average_map`, and
  `build_interpolators` from fixed (a, M sin i or M) axes to the run's
  fitted dimensions, with linear or log axes. Per-host maps keep each host's
  own injections (guiding decision 2); sparse cells are filled the way
  `fill_completeness_nans` does now, or masked if the effective count is too
  low. Extend the inclination conversion to carry e and weights.
- **Files:** `completeness_utils.py`, `rvsearch_borrowed.py`
  (`Completeness` grids and interpolation), `main.py` (`prep_maps`,
  `make_average_map`).
- **Tests:** the (a, m) case reproduces today's maps; an (a, e) map from a
  synthetic set matches its analytic completeness.
- **Done when:** Tier 1 and Tier 2 produce maps for any role assignment.

### 5. Carry eccentricity in the companion samples

- **Goal:** each posterior draw carries (a, M, e), jointly, plus its interim
  prior in every dimension.
- **Design:** `post_sampler2` also reads `esino{k}` and `ecoso{k}` and adds
  e = esino^2 + ecoso^2 from the same draws. Replace the 2xN arrays and the
  `x_samples`/`y_samples` pair in `fit_utils.CompanionSamples` and
  `fit_data.npz` with samples keyed by dimension. The interim prior is
  log-uniform in a and M and uniform in e; express it in each dimension's
  chosen scale (uniform e is p(log e) proportional to e). Marginalized
  dimensions enter each companion's likelihood term through
  target(x)/interim(x) and the ROI bounds, so a companion's weight reflects
  the marginal distribution it is averaged over.
- **Files:** `sampling_utils.py`, `fit_utils.py`, `main.py`
  (`prep_post_draws`, `prep_fit_materials`), `likelihood.py` (caches).
- **Tests:** e computed from esino/ecoso matches orvara's; old `fit_data.npz`
  loads as (a, m); the interim prior integrates correctly in each scale.
- **Done when:** fit inputs exist for any role assignment, with e available.

### 6. Support linear axes

- **Goal:** fit e (bounded, naturally linear) alongside the log axes.
- **Design:** a per-dimension `scale` ("log" or "linear") used by the grids,
  bin widths (dex or linear), occurrence-rate density units, the GP kernel's
  distance, tick formatters, and axis labels.
- **Files:** `likelihood.py`, `mcmc.py`, `fit_utils.py`, `plotting_utils.py`,
  `post_fit_analysis.py`.
- **Tests:** log axes give today's results; a linear-axis histogram fit of a
  synthetic uniform e distribution recovers it.
- **Done when:** a 1D fit in e runs end to end on synthetic data.

### 7. Eccentricity models

- **Goal:** parametric models for f(e).
- **Design:** register Beta(alpha, beta) (Kipping 2013), Rayleigh(sigma),
  and a two-component mixture (e.g. a near-circular Rayleigh plus a
  dynamically hot Rayleigh or Beta), each with density, priors,
  initialization, parameter labels, colors, table labels, and macro names in
  the existing registries (`mcmc_powerlaw.MODEL_REGISTRY`,
  `post_fit_analysis._MODEL_PARAMETER_MACROS` and
  `_MODEL_PARAMETER_LABELS`, `convergence.SMOOTH_PARAMETER_NAMES`).
  Densities are normalized over the fitted e range.
- **Files:** `likelihood.py`, `mcmc_powerlaw.py`, `mcmc.py`,
  `post_fit_analysis.py`, `convergence.py`.
- **Tests:** densities integrate to their amplitude; fits to synthetic draws
  recover their parameters.
- **Done when:** the three models fit synthetic catalogs.

### 8. Post-fit products for any fitted coordinate

- **Goal:** ORD plots, CDF comparisons, parameter tables, one- and
  two-parameter tables, and `variables.tex` work when the fitted coordinate
  is e (or any other).
- **Design:** read the fitted coordinate from the step-2 record instead of
  assuming mass or SMA (the CDF x-label already does this); add e tick
  formatting and labels; name macros by coordinate.
- **Files:** `post_fit_analysis.py`, `plotting_utils.py`, `mcmc.py`
  (`save_model_figures`).
- **Tests:** each product runs on a synthetic e fit; the mass and SMA
  products are unchanged.
- **Done when:** this repo's post-fit products can be defined (step 10).

### 9. Validate

- **Goal:** show the method recovers known eccentricity distributions with
  CLS-like completeness and priors, and check the marginalization caveat.
- **Design:** mock catalogs drawn from a known f(e) (and from a, M
  distributions), "observed" through CLS per-host completeness, with
  orvara-like posteriors (widths and the uniform-e prior). Check recovery of
  f(e), the size of the bias toward high e at low signal-to-noise
  (Lucy-Sweeney), and sensitivity to the marginal distributions (vary them
  within their uncertainties). On the CLS itself, fit f(e) in two
  companion-mass bins and two separation bins to test whether e depends on a
  marginalized dimension.
- **Files:** new `validation/` scripts in this repo; small helpers in
  `occurrence` if needed.
- **Done when:** recovery is unbiased within its uncertainty for the planned
  sample sizes, and the split test is reported.

### 10. Set up this repo's workflow

- **Goal:** `run_occurrence.py`, `run_post_fit_analysis.py`,
  `collect_paper_items.py`, and tests, mirroring the SMA repo.
- **Design:** first runs: e fitted (0-1, linear) with a and m marginalized
  (0.1-30 AU, 0.4-50 M_Jup) for all stars and the high/low mass,
  metallicity, and age samples; models: histogram, Beta, Rayleigh, mixture.
  Exclude a < 0.1 AU (tides).
- **Done when:** the runs validate and the post-fit products resolve.

## Open questions

- Fit e in its own scale, or in sqrt(e) or log e for better behaved
  posteriors near e = 0?
- Minimum effective injection count per cell before a per-host cell is
  masked, and how masked cells enter the likelihood.
- Whether the expected-detection term should also use per-host completeness
  (summed over hosts) rather than the sample average once maps carry
  weights.

## Implementation notes (ecc_added, first pass)

Decisions made while implementing, and how they differ from the steps
above. Each is documented in the code where it applies.

- **Completeness in e: slices, not per-injection weights (steps 3-4).**
  Each star's injections are split into eccentricity slices (edges 0, 0.1,
  0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1; `completeness_utils.ECC_SLICE_EDGES`), and
  each slice gets an ordinary (a, m) map from the original moving-average
  estimator and NaN-fill rules on a grid shared by all slices
  (`sliced_map_maker`). Completeness is piecewise constant in e within a
  slice. Because injections are log-uniform in (a, M sin i) inside the K
  window, binning in (a, m) needs no reweighting there, and the high-mass
  fill rule handles the window's empty corner, so the edge bias described
  in step 3 does not arise. A single slice spanning [0, 1] reproduces the
  original map exactly (tested). Marginalizing integrates the sliced maps
  over the marginalized dimensions with their distributions' exact cell
  probabilities (`engine_completeness`), and stops with an error if a
  marginalized range contains a point with no estimate.
- **Linear axes without changing the engine (step 6).** Eccentricity
  enters the engine as u = 10^e, so the engine's log10 axes are linear in e
  (a bin uniform in e is uniform in log10 u, and orvara's uniform-e prior is
  uniform in log10 u). `dimensions.DimensionSpec.to_engine`/`from_engine`
  convert; plots label ticks with e (`plotting_utils.axis_style`).
- **Engine slots.** A run fits one or two of a, m, e. The engine's x slot
  holds the first fitted dimension (order a, m, e) and its y slot the
  second; for a 1D fit, y holds a marginalized dimension (preference m, a,
  e) as one bin, and the engine stacks over y.
- **Companion weights.** Each draw's completeness-over-prior is multiplied
  by the product of the marginalized dimensions' population densities (zero
  outside their bounds); orvara's interim priors are flat in each
  dimension's coordinate, so no further division is needed
  (`main.weigh_dimension_samples`).
- **Host completeness fallback.** Host values that are NaN *or exactly
  zero* fall back to the sample average at the same point. Zeros were
  added because slicing leaves few injections per box, so "0 of N
  recovered" occurs where true completeness is small but positive; a
  detected companion cannot have zero completeness. In the first all-star
  run about 21% of draws fell back.
- **Sliced average.** Points where no star has an estimate in a slice take
  the value of the next lower-e slice (1,114 points over all slices for all
  stars, mostly in the top slice).
- **Expected detections** use the sample-average completeness only
  (`use_average_completeness=True`); per-star exposure is not implemented
  for explicit roles.
- **Mass-ratio Tier 1 folders** (qtrue, qsini) are not supported with
  explicit roles yet: the default mass distribution is in M_Jup.
- **Missing posterior eccentricities.** GL 876 (b, c, d), HD 219134
  (b-f), and HD 28185 (b, c) use literature values, drawn from split
  normals truncated to [0, 1) independently of a and m
  (`sampling_utils.LITERATURE_ECCENTRICITIES`, matched to the posterior
  files' companion indices by semi-major axis).
- **Eccentricity models (step 7):** `beta` (A, alpha, beta), `rayleigh`
  (A, sigma_e), and `ecc_mixture` (A, f1, sigma1, sigma2; sigma1 < sigma2),
  each A f(e) with f normalized over the fitted e range. Priors are uniform
  in the physical parameters: alpha, beta in [0.1, 20]; sigmas in
  [0.005, 1]; f1 in [0, 1]. Like every model, A is per unit of the stacked
  coordinate (per dex of mass when mass is the one-bin y slot); the
  integrated occurrence the fits report multiplies by the stack width.
- **Speed.** `completeness_utils.moving_average_grid` reproduces the
  original moving-average estimator exactly and is about 500 times faster;
  the sliced maps use it.
- **Workflow (step 10).** Run configurations set `dimensions` in place of
  `a_edges`/`m_edges` (`run_occurrence.ECC_DIMENSIONS`); `run_multiple`
  now accepts either. Runs: `ecc_bounds` (all stars and the high/low mass,
  metallicity, and activity samples; histogram, Beta, Rayleigh, and
  two-Rayleigh models) and `stellar2params` (histogram only). Post-fit
  products stack over m and use the eccentricity models.

## First results (2026-10-22, `ecc_added` 4c9c3ab, not production runs)

- **Regression.** Seeded, shortened mass- and SMA-function runs on master
  and on `ecc_added` match in all 6,515 output arrays and tables
  (`validation/regression_run.py`, `compare_runs.py`).
- **Mocks** (`validation/mock_validation.py`; truth Beta(1.5, 4), 0.30
  companions per star; three seeds each). With e uncertainties of 0.05,
  the rate, alpha, and beta are recovered within about 1 sigma, scattered
  on both sides. With 0.15, rates are recovered but alpha and beta tend
  low (means 1.33 and 3.42), a slightly broader f(e); each seed is within
  about 1.3 sigma, so more seeds are needed before calling it a bias.
- **All-star CLS fit** (0.1-30 AU, 0.4-50 MJ; 104 companions mostly in
  the region): 0.27-0.29 companions per star. Beta alpha = 0.87 +/- 0.08,
  beta = 1.9 +/- 0.3; Rayleigh sigma = 0.26 +/- 0.02. Delta BIC against a
  flat f(e): Beta +26, two-Rayleigh +26, Rayleigh -47 (positive favors the
  model). The two-Rayleigh narrow component sits at its prior floor
  (sigma1 = 0.005, f1 = 0.27): about ten near-circular companions, mostly
  warm Jupiters at 0.1-0.2 AU plus GL 876 b and c.
- **Cost.** The all-star run (four models, 3000-step chains) took 1.4 h and
  0.36 GB; the two-Rayleigh and histogram chains need roughly 4 and 1.5
  times more steps to reach 50 autocorrelation times.

## Open issues from the first results

- The e = 0.8-1 slice has about 60 injections per star and ~70% of each
  star's cells have no estimate, so companions there use the average map.
  Wider top slices, or more injections, would help.
- Tidal circularization: should the inner bound move from 0.1 AU (the
  near-circular warm Jupiters drive the narrow mixture component)?
- The two-Rayleigh sigma1 prior floor (0.005) sets that component's width.
