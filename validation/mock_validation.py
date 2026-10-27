"""Mock-catalog check that eccentricity fits recover a known distribution.

Usage: python mock_validation.py <ecc_run_dir> <mock_dir> <sigma_e> <seed>

Draws companions around the CLS stars from a known population (a and m from
the run's marginal distributions, flat by default, e from Beta(1.5, 4)), detects each with
its host's eccentricity-sliced completeness, gives detected companions
Gaussian "posteriors" (flat priors: log-uniform in a and m, uniform in e,
like orvara), and weights them with the same code as the real catalog. Then
fits the eccentricity distribution and reports the recovered Beta
parameters. Reuses the Tier 1 and Tier 2 maps of <ecc_run_dir>.
"""
import json
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

GITHUB = Path("/Users/judahvz/research/code/GitHub")
ECC_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(GITHUB))
sys.path.insert(0, str(ECC_REPO))

import run_occurrence as ro  # noqa: E402
from occurrence import completeness_utils as cu  # noqa: E402
from occurrence import dimensions as dm  # noqa: E402
from occurrence import main, run  # noqa: E402
from occurrence import marginal_distributions as md  # noqa: E402

source = Path(sys.argv[1]).resolve()
mock = Path(sys.argv[2]).resolve()
sigma_e = float(sys.argv[3])
seed = int(sys.argv[4])
rng = np.random.RandomState(seed)

TRUE_ALPHA, TRUE_BETA, TRUE_RATE = 1.5, 4.0, 0.30   # rate per star in the ROI
SIGMA_LOG_A, SIGMA_LOG_M, N_SAMPLES = 0.02, 0.08, 500
DIMENSIONS = ro.ECC_DIMENSIONS
spec = dm.parse_dimensions(DIMENSIONS)

# Reuse the real run's Tier 1 (all maps) and Tier 2 sliced average
tier1 = mock / "mtrue"
tier2 = tier1 / "allstars"
tier2.mkdir(parents=True, exist_ok=True)
for name in ("mtrue_recoveries", "saved_maps_mtrue", "saved_maps_mtrue_eslices"):
    if not (tier1 / name).exists():
        os.symlink(source / "mtrue" / name, tier1 / name)
for name in ("avg_map_eslices", f"avg_map_{spec.token}"):
    if not (tier2 / name).exists():
        shutil.copytree(source / "mtrue" / "allstars" / name, tier2 / name)

star_df = ro.load_or_make_star_df()
average = cu.load_sliced(tier2 / "avg_map_eslices")
maps_dir = Path(main.sliced_maps_dir(tier1, "mtrue"))


def draw_from(distribution, size):
    """Inverse-CDF draws from a marginal distribution (fine grid)."""
    low, high = distribution.coordinate(distribution.bounds)
    grid = np.linspace(low, high, 20001)
    values = 10**grid if distribution.scale == "log" else grid
    cdf = np.concatenate([[0], np.cumsum(distribution.cell_probabilities(values))])
    coordinate = np.interp(rng.uniform(size=size), cdf, grid)
    return 10**coordinate if distribution.scale == "log" else coordinate


n_true = rng.poisson(TRUE_RATE*len(star_df))
hosts = rng.randint(0, len(star_df), n_true)
a_true = draw_from(spec.dimensions["a"].distribution, n_true)
m_true = draw_from(spec.dimensions["m"].distribution, n_true)
e_true = rng.beta(TRUE_ALPHA, TRUE_BETA, n_true)

samples, comp_lists = {}, [[] for _ in range(len(star_df))]
detected_e = []
for index in range(n_true):
    star = star_df.star_name.iloc[hosts[index]]
    host = cu.load_sliced(maps_dir / star)
    point = (np.array([a_true[index]]), np.array([m_true[index]]),
             np.array([e_true[index]]))
    q = host(point)[0]
    if not np.isfinite(q):
        q = average(point)[0]
    if rng.uniform() >= q:
        continue
    # Observed values scatter around the truth; with flat priors the
    # posterior is a Gaussian around the observed value (truncated for e)
    log_a_obs = np.log10(a_true[index]) + rng.normal(0, SIGMA_LOG_A)
    log_m_obs = np.log10(m_true[index]) + rng.normal(0, SIGMA_LOG_M)
    e_obs = e_true[index] + rng.normal(0, sigma_e)
    a_post = 10**(log_a_obs + rng.normal(0, SIGMA_LOG_A, N_SAMPLES))
    m_post = 10**(log_m_obs + rng.normal(0, SIGMA_LOG_M, N_SAMPLES))
    e_post = np.empty(0)
    while e_post.size < N_SAMPLES:
        draws = e_obs + rng.normal(0, sigma_e, 4*N_SAMPLES)
        e_post = np.concatenate([e_post, draws[(draws >= 0) & (draws < 1)]])
    name = f"mock_{index}"
    samples[name] = [a_post, m_post, e_post[:N_SAMPLES]]
    comp_lists[hosts[index]].append(name)
    detected_e.append(e_true[index])

mock_stars = star_df.copy()
mock_stars["comp_list"] = comp_lists
catalog, _ = main.weigh_dimension_samples(spec, samples, mock_stars, maps_dir,
                                          average)
np.savez(main.dimension_tier2_paths(tier1, "allstars", spec)["catalog"],
         **catalog)
print(f"{n_true} true companions, {len(samples)} detected; "
      f"median detected e = {np.median(detected_e):.3f}, true median "
      f"{np.median(e_true):.3f}")

arguments = dict(ro.RUN_DEFAULTS)
arguments.update({
    "star_df": mock_stars, "tier2_df_cuts_dict": ro.tier2_df_cuts_dict,
    "output_dir": mock, "tier1_list": ["mtrue"], "tier2_list": ["allstars"],
    "tier3_list": ["mock_fit"],
    "dimensions": DIMENSIONS, "run_models_list": ["piecewise", "beta"],
    "plot_models_list": ["piecewise", "beta"], "parallel_fits": False,
    "nsteps": 2000, "burnin": 1000, "piecewise_nsteps": 4000,
    "piecewise_burnin": 1500, "random_seed": seed,
    "prepare_missing": False,
})
run.run_multiple(**arguments)

chains = np.load(mock / "mtrue/allstars/mock_fit/saved_chains/chains_beta_bin0.npz")
flat = chains["flat_chains"]
alpha, beta = flat[:, 1], flat[:, 2]
summary = {
    "sigma_e": sigma_e, "seed": seed, "n_true": int(n_true),
    "n_detected": len(samples),
    "alpha": np.percentile(alpha, [16, 50, 84]).tolist(),
    "beta": np.percentile(beta, [16, 50, 84]).tolist(),
    "true": [TRUE_ALPHA, TRUE_BETA],
}
(mock / "mock_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary))
