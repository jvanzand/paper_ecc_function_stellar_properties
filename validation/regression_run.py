"""Seeded end-to-end runs used as a regression reference for ecc_added.

Usage: python regression_run.py <output_dir>

Runs Tier 1-3 serially with fixed seeds (np.random for the inclination draws
and posterior sampling, random_seed for emcee), so two runs of the same code
give identical arrays. Shortened chains keep it to a few minutes.
"""
import sys
import time
from pathlib import Path

import numpy as np

GITHUB = Path("/Users/judahvz/research/code/GitHub")
MASS_REPO = GITHUB / "paper_mass_function_stellar_properties"
sys.path.insert(0, str(GITHUB))
sys.path.insert(0, str(MASS_REPO))

import run_occurrence as ro  # noqa: E402  (mass repo: catalog, paths, defaults)
from occurrence import run, sampling_utils as su  # noqa: E402

OUTPUT = Path(sys.argv[1]).resolve()
OUTPUT.mkdir(parents=True, exist_ok=True)
star_df = ro.load_or_make_star_df()
paths = ro.load_local_paths()

COMMON = dict(ro.RUN_DEFAULTS)
COMMON.update({
    "star_df": star_df, "tier2_df_cuts_dict": ro.tier2_df_cuts_dict,
    "output_dir": OUTPUT, "recoveries_dir": paths["recoveries_dir"],
    "comp_post_dir": paths["posteriors_dir"], "sampling_func": su.post_sampler2,
    "parallel_fits": False, "make_plots": False,
    "nsteps": 600, "burnin": 300, "piecewise_nsteps": 1200,
    "piecewise_burnin": 400, "random_seed": 1234,
})
RUNS = {
    # Mass function, as in the mass paper
    "reg_mass": dict(
        tier1_list=["mtrue", "qtrue"], tier2_list=["allstars", "highMstar"],
        tier3_list=["reg_mass"], a_edges=[0.1, 10.0],
        m_edges=[0.4, 0.8, 1.6, 3.2, 6.4, 13.0, 26.0, 50.0], stack_dim="a",
        run_models_list=["piecewise", "sigmoid", "logG"],
        plot_models_list=["piecewise"],
    ),
    # SMA function, as in the SMA repo
    "reg_sma": dict(
        tier1_list=["mtrue"], tier2_list=["allstars", "lowMstar"],
        tier3_list=["reg_sma"],
        a_edges=[float(edge) for edge in np.geomspace(0.1, 30.0, 10)],
        m_edges=[0.4, 50.0], stack_dim="m",
        run_models_list=["piecewise", "bpl", "logG"],
        plot_models_list=["piecewise"],
    ),
}

start = time.perf_counter()
np.random.seed(20261022)
for name, settings in RUNS.items():
    arguments = dict(COMMON)
    arguments.update(settings)
    run.run_multiple(**arguments)
print(f"Regression run finished in {time.perf_counter() - start:.0f} s")
