"""Run all occurrence calculations for the eccentricity-function paper.

Post-fit tables, LaTeX variables, and summary plots built from these results
live in ``run_post_fit_analysis.py``.

Edit ``RUN_CONFIGURATIONS`` for scientific choices and ``RUNS_TO_DO`` to
choose which of them run, or call :func:`main` from Python to select runs
and execution behavior. Paths are resolved from
this file, so launching it from another working directory is safe.
"""

import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parent
OCCURRENCE_DIR = PROJECT_DIR.parent / "occurrence"
RESULTS_DIR = PROJECT_DIR / "results"
STAR_CATALOG = PROJECT_DIR / "data" / "derived" / "star_catalog.json"

# Machine-specific input locations live in an untracked JSON file, because the
# recoveries and posteriors are too large to keep in the repository.
LOCAL_PATHS_FILE = PROJECT_DIR / "local_paths.json"
LOCAL_PATHS_EXAMPLE = PROJECT_DIR / "local_paths.example.json"
LOCAL_PATH_KEYS = ("recoveries_dir", "posteriors_dir")
# Full results of the latest main() call, written beneath the output directory
RUN_SUMMARY_NAME = "run_summary.json"

# occurrence is currently a sibling checkout rather than an installed package.
if str(OCCURRENCE_DIR.parent) not in sys.path:
    sys.path.insert(0, str(OCCURRENCE_DIR.parent))

from occurrence import dimensions as occurrence_dimensions  # noqa: E402
from occurrence import run as occurrence_run  # noqa: E402
from occurrence import sampling_utils as su  # noqa: E402

from config_dict import tier2_df_cuts_dict  # noqa: E402


# Every run fits eccentricity in five bins of width 0.2 and marginalizes
# semi-major axis (0.1-30 AU) and companion mass (0.4-50 MJ) with the flat
# defaults of occurrence.marginal_distributions (uniform per dex, the same
# weighting as one stacked bin). The literature distributions are opt-in
# through occurrence.marginal_distributions.population_marginal. Marginalizing
# assumes a and m are independent of e (see that module's docstring).
ECC_EDGES = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
ECC_DIMENSIONS = {
    "e": {"role": "fit", "edges": ECC_EDGES},
    "a": {"role": "marginalize", "bounds": (0.1, 30.0)},
    "m": {"role": "marginalize", "bounds": (0.4, 50.0)},
}
ECC_MODELS = ["piecewise", "beta", "rayleigh", "ecc_mixture"]

# Each entry becomes one call to occurrence.run.run_multiple(). Keep separate
# experiments here instead of copying or commenting out call blocks.
RUN_CONFIGURATIONS = {

    # Main run: eccentricity functions of the full sample and each stellar
    # subsample
    "ecc_bounds": {
        "tier1_list": ["mtrue"],
        "tier2_list": [
            "allstars", "highMstar", "lowMstar", "highFeH", "lowFeH",
            "highAct", "lowAct",
        ],
        "tier3_list": ["ecc_bounds"],
        "dimensions": ECC_DIMENSIONS,
        "run_models_list": ECC_MODELS,
        "plot_models_list": ECC_MODELS,

        "run_fits": True,
        "make_plots": True,
    },

    # Mass-metallicity subsets for the two-parameter tables
    "stellar2params": {
        "tier1_list": ["mtrue"],
        "tier2_list": ["highMstarhighFeH", "highMstarlowFeH",
                       "lowMstarhighFeH", "lowMstarlowFeH"],
        "tier3_list": ["stellar2params"],
        "dimensions": ECC_DIMENSIONS,
        "run_models_list": ["piecewise"],
        "plot_models_list": ["piecewise"],

        "run_fits": True,
        "make_plots": True,
    },

}


# Names from RUN_CONFIGURATIONS to execute when main() is called without
# run_names (including ``python run_occurrence.py``). Set to None to run every entry.
RUNS_TO_DO = ['ecc_bounds']


# Keys a run with ``reuse_fits_from`` may set. Everything else, including bin
# edges and fit settings, is inherited from the source run so the plots always
# describe the fits they reuse.
PLOT_OPTION_KEYS = {
    "make_plots", "plot_models_list", "model_plot_style",
    "n_posterior_draws", "plot_random_seed", "plot_occurrence",
    "plot_density", "plot_cumulative", "plot_corner", "plot_catalog_roi",
    "plot_roi_occurrence", "plot_uncorrected_occurrence_mle",
    "occurrence_legend_loc",
}
DERIVED_RUN_KEYS = PLOT_OPTION_KEYS | {
    "reuse_fits_from", "tier1_list", "tier2_list",
}


RUN_DEFAULTS = {
    "recoveries_mtype": "msini",
    "run_fits": True,
    "make_plots": True,
    "prepare_missing": True,
    "avg_map_only": False,
    "fill_single_nan_with_average": True,
    "run_models_list": ["piecewise"],
    "plot_models_list": ["piecewise"],
    "model_plot_style": "credible",
    "n_posterior_draws": 100,
    "plot_random_seed": 2222,
    "stack_dim": "a",
    "m_unit": "jupiter",
    "plot_occurrence": True,
    "plot_density": True,
    "plot_cumulative": True,
    "plot_corner": True,
    "plot_catalog_roi": True,
    "plot_roi_occurrence": True,
    "plot_uncorrected_occurrence_mle": True,
    "occurrence_legend_loc": "upper right",
    "completeness_type": "single",
    "plot_tier1_maps": False,
    "integration_resolution": (100, 100),
    "use_average_completeness": True,
    "piecewise_parameterization": "independent",
    "piecewise_gp_amplitude": 1.0,
    "piecewise_gp_length_scale_a": None,
    "piecewise_gp_length_scale_m": None,
    "piecewise_gp_jitter": 1e-8,
    "piecewise_gp_infer_hyperparameters": True,
    "piecewise_gp_amplitude_bounds": (0.05, 5.0),
    "piecewise_gp_length_scale_a_bounds": None,
    "piecewise_gp_length_scale_m_bounds": None,
    "max_integrated_occurrence": 1.0,
    "logg_amplitude_bounds": (1e-6, 10.0),
    "logg_sigma_bounds": None,
    "escarpment_amplitude_bounds": (1e-6, 10.0),
    "sigmoid_amplitude_bounds": (1e-6, 10.0),
    "sigmoid_width_bounds": None,
    "bpl_amplitude_bounds": (1e-6, 10.0),
    "bpl_slope_bounds": (-4.0, 4.0),
    "ecc_amplitude_bounds": (1e-6, 10.0),
    "loglinear_amplitude_bounds": (1e-6, 10.0),
    "nwalkers": 50,
    "nsteps": 500,
    "burnin": 200,
    # Production and burn-in steps for the piecewise model alone; None uses
    # nsteps and burnin.
    "piecewise_nsteps": 3000,
    "piecewise_burnin": 500,
    # Fits whose chains span fewer autocorrelation times than this are named
    # in a warning after each run; every experiment also gets
    # saved_chains/convergence.txt with per-parameter details.
    "convergence_target": 50,
    "random_seed": 1234,
    "parallel_fits": True,
    "parallel_mcmc": False,
}


def _catalog_sources(occurrence_dir=OCCURRENCE_DIR):
    cls_dir = Path(occurrence_dir) / "cls_files"
    return {
        "stars": cls_dir / "clsI_tab2.csv",
        "activity": cls_dir / "clsV_tab2.csv",
        "companions": cls_dir / "cls_substellar_comps.csv",
    }


def _file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _catalog_metadata(occurrence_dir=OCCURRENCE_DIR):
    return {
        "schema_version": 1,
        "source_sha256": {
            name: _file_sha256(path)
            for name, path in _catalog_sources(occurrence_dir).items()
        },
    }


def make_star_df(occurrence_dir=OCCURRENCE_DIR):
    """Build the complete CLS stellar catalog used by the paper."""
    sources = _catalog_sources(occurrence_dir)
    cls_i = pd.read_csv(sources["stars"])
    cls_v = pd.read_csv(sources["activity"])
    combined = cls_i.merge(cls_v, on="cps_identifier", how="left")

    columns = [
        "cps_identifier", "Mstar", "Mstar_err", "feh", "feh_err",
        "logrhk", "age",
    ]
    star_df = combined[columns].copy().rename(
        columns={"cps_identifier": "star_name"}
    )
    companions = pd.read_csv(sources["companions"])

    companion_lists = []
    for star_name in star_df["star_name"]:
        rows = companions.loc[companions["cps_identifier"] == star_name]
        companion_names = (
            rows["cps_identifier"] + "_" +
            rows.groupby("cps_identifier").cumcount().astype(str)
        ).tolist()
        # The published sample intentionally excludes the spurious outer
        # companion assigned to HD 167215.
        if "167215" in star_name and companion_names:
            companion_names = companion_names[:1]
        companion_lists.append(companion_names)

    star_df["comp_list"] = companion_lists
    star_df["comp_num"] = star_df["comp_list"].map(len)
    validate_star_df(star_df)
    return star_df


def validate_star_df(star_df, sample_names=None):
    """Validate the catalog and the sample queries selected for a run."""
    required = {"star_name", "Mstar", "feh", "logrhk", "comp_list"}
    missing = required - set(star_df.columns)
    if missing:
        raise ValueError("star catalog is missing columns: {}".format(
            sorted(missing)
        ))
    if star_df.empty:
        raise ValueError("star catalog is empty")
    if star_df["star_name"].isna().any() or star_df["star_name"].duplicated().any():
        raise ValueError("star_name values must be non-null and unique")
    if star_df["Mstar"].isna().any() or (star_df["Mstar"] <= 0).any():
        raise ValueError("Mstar values must be non-null and positive")
    if not star_df["comp_list"].map(lambda value: isinstance(value, list)).all():
        raise ValueError("every comp_list value must be a list")

    sample_names = list(sample_names or [])
    unknown = sorted(set(sample_names) - set(tier2_df_cuts_dict))
    if unknown:
        raise ValueError("unknown stellar samples: {}".format(unknown))
    for name in sample_names:
        query = tier2_df_cuts_dict[name][0]["star_df_query"]
        try:
            selected = star_df.query(query) if query else star_df
        except Exception as error:
            raise ValueError(
                "invalid query for sample {!r}: {!r}".format(name, query)
            ) from error
        if selected.empty:
            raise ValueError("sample {!r} selects no stars".format(name))


def load_or_make_star_df(cache_path=STAR_CATALOG, rebuild=False):
    """Load the frozen catalog, or rebuild and cache it from the CLS tables."""
    cache_path = Path(cache_path)
    if cache_path.is_file() and not rebuild:
        with cache_path.open() as stream:
            cached = json.load(stream)
        if not isinstance(cached, dict) or "records" not in cached:
            raise ValueError(
                "star catalog cache uses an obsolete format; rebuild it with "
                "rebuild_star_catalog=True"
            )
        star_df = pd.DataFrame(
            cached["records"], columns=cached.get("columns")
        )
        validate_star_df(star_df)
        return star_df

    star_df = make_star_df()
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    records = json.loads(
        star_df.to_json(orient="records", double_precision=15)
    )
    cached = _catalog_metadata()
    cached["columns"] = list(star_df.columns)
    cached["records"] = records
    with cache_path.open("w") as stream:
        json.dump(cached, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return star_df


def load_local_paths(path=LOCAL_PATHS_FILE):
    """Read this machine's input directories from ``local_paths.json``."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(
            "{} not found; copy {} to {} and set this machine's "
            "paths".format(path, LOCAL_PATHS_EXAMPLE.name, path.name)
        )
    with path.open() as stream:
        configured = json.load(stream)
    missing = [key for key in LOCAL_PATH_KEYS if not configured.get(key)]
    if missing:
        raise ValueError("{} is missing entries: {}".format(path, missing))
    return {
        key: (PROJECT_DIR / Path(configured[key]).expanduser()).resolve()
        for key in LOCAL_PATH_KEYS
    }


def _selected_run_names(run_names):
    if run_names is None:
        run_names = RUNS_TO_DO
    names = list(RUN_CONFIGURATIONS) if run_names is None else list(run_names)
    unknown = sorted(set(names) - set(RUN_CONFIGURATIONS))
    if unknown:
        raise ValueError("unknown runs: {}".format(unknown))
    if not names:
        raise ValueError("run_names must select at least one run")
    return names


def resolve_run_configuration(name):
    """Return the full configuration for a run, expanding ``reuse_fits_from``.

    A derived run starts from its source run's configuration, overrides only
    plot options and tier subsets, writes to a Tier 3 folder named after
    itself, and never fits.
    """
    entry = RUN_CONFIGURATIONS[name]
    source_name = entry.get("reuse_fits_from")
    if source_name is None:
        if list(entry["tier3_list"]) != [name]:
            raise ValueError(
                "{} must write to a Tier 3 folder named after itself; set "
                "\"tier3_list\": [{!r}]".format(name, name)
            )
        return dict(entry)

    if source_name not in RUN_CONFIGURATIONS:
        raise ValueError("{} reuses fits from unknown run {!r}".format(
            name, source_name
        ))
    source = RUN_CONFIGURATIONS[source_name]
    if "reuse_fits_from" in source:
        raise ValueError("{} cannot reuse fits from derived run {!r}".format(
            name, source_name
        ))
    if len(source["tier3_list"]) != 1:
        raise ValueError(
            "{} must reuse a run with exactly one Tier 3 folder".format(name)
        )
    disallowed = sorted(set(entry) - DERIVED_RUN_KEYS)
    if disallowed:
        raise ValueError(
            "{} reuses fits from {!r} and may only set plot options and tier "
            "subsets; remove {}".format(name, source_name, disallowed)
        )
    for tier_key in ("tier1_list", "tier2_list"):
        extra = sorted(set(entry.get(tier_key, [])) - set(source[tier_key]))
        if extra:
            raise ValueError("{} {} entries are not in {!r}: {}".format(
                name, tier_key, source_name, extra
            ))

    configuration = dict(source)
    configuration.update(entry)
    configuration["tier3_list"] = [name]
    configuration["source_tier3"] = source["tier3_list"][0]
    configuration["run_fits"] = False
    return configuration


def _reused_fit_paths(name, configuration, output_dir):
    """Map each derived Tier 3 folder to the source folder it reuses."""
    source_tier3 = configuration["source_tier3"]
    paths = []
    for tier1 in configuration["tier1_list"]:
        for tier2 in configuration["tier2_list"]:
            tier2_dir = Path(output_dir) / tier1 / tier2
            paths.append((tier2_dir / name, tier2_dir / source_tier3))
    return paths


def _check_reused_fits(name, configuration, output_dir):
    """Require saved fits for every plotted model before anything runs."""
    missing = []
    for _, source_dir in _reused_fit_paths(name, configuration, output_dir):
        chain_dir = source_dir / "saved_chains"
        required = [source_dir / "saved_dicts" / "fit_data.npz"]
        for model_name in configuration["plot_models_list"]:
            if model_name == "piecewise":
                required.append(chain_dir / "chains_piecewise.npz")
            elif not list(chain_dir.glob("chains_{}_bin*.npz".format(
                    model_name))):
                required.append(chain_dir / "chains_{}_bin*.npz".format(
                    model_name))
        missing.extend(str(path) for path in required if not path.exists())
    if missing:
        raise FileNotFoundError(
            "{} reuses fits that do not exist yet; run {!r} first. "
            "Missing: {}".format(
                name, configuration["reuse_fits_from"], missing
            )
        )


def _relative_symlink(link, target):
    """Point ``link`` at ``target`` with a relative path, reusing a match."""
    relative_target = os.path.relpath(target, link.parent)
    if link.is_symlink():
        if os.readlink(link) == relative_target:
            return
        raise FileExistsError(
            "{} already links to {}, not {}".format(
                link, os.readlink(link), relative_target
            )
        )
    if link.exists():
        raise FileExistsError(
            "{} exists and is not a link to {}".format(link, target)
        )
    link.symlink_to(relative_target)


def link_reused_fits(name, configuration, output_dir):
    """Link a derived run's folders to the saved fits of its source run.

    ``saved_chains`` is linked whole. ``saved_dicts`` is a real folder holding
    a link to ``fit_data.npz``, because plotting rewrites the piecewise summary
    there and must not overwrite the source's copy.
    """
    paths = _reused_fit_paths(name, configuration, output_dir)
    for derived_dir, source_dir in paths:
        (derived_dir / "saved_dicts").mkdir(parents=True, exist_ok=True)
        _relative_symlink(
            derived_dir / "saved_chains", source_dir / "saved_chains"
        )
        _relative_symlink(
            derived_dir / "saved_dicts" / "fit_data.npz",
            source_dir / "saved_dicts" / "fit_data.npz",
        )
    return [str(derived_dir) for derived_dir, _ in paths]


def _validate_run_configuration(name, configuration):
    """Check a run before anything is computed.

    A run sets either ``dimensions`` (which dimensions are fit or
    marginalized; occurrence.dimensions.parse_dimensions checks it) or the
    original ``a_edges`` and ``m_edges``, not both.
    """
    if "dimensions" in configuration:
        if "a_edges" in configuration or "m_edges" in configuration:
            raise ValueError(
                "{} sets dimensions, which replace a_edges and m_edges; "
                "remove the edges".format(name)
            )
        try:
            occurrence_dimensions.parse_dimensions(configuration["dimensions"])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError("{} dimensions: {}".format(name, error)) from error
        # The default mass distribution is in Jupiter masses, so mass-ratio
        # Tier 1 spaces cannot be marginalized or fit with it yet
        ratio_tiers = sorted(set(configuration.get("tier1_list", []))
                             & {"qtrue", "qsini"})
        if ratio_tiers:
            raise ValueError(
                "{} sets dimensions, which do not support mass-ratio Tier 1 "
                "spaces yet: {}".format(name, ratio_tiers)
            )
        edge_names = ()
    else:
        edge_names = ("a_edges", "m_edges")
    for edge_name in edge_names:
        edges = configuration[edge_name]
        if len(edges) < 2 or any(
                right <= left for left, right in zip(edges, edges[1:])):
            raise ValueError(
                "{} {} must be strictly increasing".format(name, edge_name)
            )
    settings = dict(RUN_DEFAULTS)
    settings.update(configuration)
    if settings["burnin"] >= settings["nsteps"]:
        raise ValueError("{} burnin must be smaller than nsteps".format(name))
    piecewise_nsteps = settings["piecewise_nsteps"] or settings["nsteps"]
    piecewise_burnin = (settings["burnin"] if settings["piecewise_burnin"] is None
                        else settings["piecewise_burnin"])
    if piecewise_burnin >= piecewise_nsteps:
        raise ValueError(
            "{} piecewise_burnin must be smaller than piecewise_nsteps".format(name)
        )


def _print_plan(configurations, star_df, output_dir):
    print("Catalog: {} stars".format(len(star_df)))
    print("Results: {}".format(output_dir))
    for name, config in configurations.items():
        if "reuse_fits_from" in config:
            print("{} (plots only, reusing {} fits):".format(
                name, config["reuse_fits_from"]
            ))
        else:
            print("{}:".format(name))
        for sample_name in config["tier2_list"]:
            query = tier2_df_cuts_dict[sample_name][0]["star_df_query"]
            count = len(star_df.query(query)) if query else len(star_df)
            print("  {}: {} stars".format(sample_name, count))


def main(
        run_names=None,
        validate_only=False,
        dry_run=False,
        plots_only=False,
        run_fits=None,
        make_plots=None,
        rebuild_star_catalog=False,
        output_dir=RESULTS_DIR,
        recoveries_dir=None,
        comp_post_dir=None):
    """Validate and execute the selected paper calculations.

    ``run_names`` selects keys from ``RUN_CONFIGURATIONS``; when omitted,
    ``RUNS_TO_DO`` is used, and ``None`` there means every run. ``plots_only``
    reuses saved chains instead of fitting. ``run_fits=True`` always reruns
    the Tier 3 fits and overwrites their saved products. Set
    ``make_plots=False`` to skip plots.
    ``dry_run`` reports expanded sample sizes and paths without creating
    results. ``validate_only`` also verifies required input directories, then
    stops. ``recoveries_dir`` and ``comp_post_dir`` default to the
    ``recoveries_dir`` and ``posteriors_dir`` entries in ``local_paths.json``.
    Runs with ``reuse_fits_from`` always skip fitting, whatever ``run_fits``
    says, and require the source run's saved fits to exist.
    """
    names = _selected_run_names(run_names)
    configurations = {
        name: resolve_run_configuration(name) for name in names
    }
    for name, configuration in configurations.items():
        _validate_run_configuration(name, configuration)
    selected_samples = [
        sample
        for configuration in configurations.values()
        for sample in configuration["tier2_list"]
    ]
    star_df = load_or_make_star_df(rebuild=rebuild_star_catalog)
    validate_star_df(star_df, selected_samples)

    if recoveries_dir is None or comp_post_dir is None:
        local_paths = load_local_paths()
        if recoveries_dir is None:
            recoveries_dir = local_paths["recoveries_dir"]
        if comp_post_dir is None:
            comp_post_dir = local_paths["posteriors_dir"]
    recoveries_dir = Path(recoveries_dir)
    comp_post_dir = Path(comp_post_dir)
    output_dir = Path(output_dir)
    missing_inputs = [
        path for path in (recoveries_dir, comp_post_dir) if not path.is_dir()
    ]
    if missing_inputs:
        raise FileNotFoundError(
            "required input directories do not exist: {}".format(
                [str(path) for path in missing_inputs]
            )
        )

    _print_plan(configurations, star_df, output_dir)
    if validate_only:
        for name, configuration in configurations.items():
            if "reuse_fits_from" in configuration:
                _check_reused_fits(name, configuration, output_dir)
    if validate_only or dry_run:
        return []

    if plots_only and run_fits is True:
        raise ValueError("plots_only=True conflicts with run_fits=True")

    results = []
    for name, configuration in configurations.items():
        derived = "reuse_fits_from" in configuration
        if derived:
            _check_reused_fits(name, configuration, output_dir)
            link_reused_fits(name, configuration, output_dir)
        arguments = dict(RUN_DEFAULTS)
        arguments.update(configuration)
        for key in ("reuse_fits_from", "source_tier3"):
            arguments.pop(key, None)
        arguments.update({
            "star_df": star_df,
            "tier2_df_cuts_dict": tier2_df_cuts_dict,
            "output_dir": output_dir,
            "recoveries_dir": recoveries_dir,
            "comp_post_dir": comp_post_dir,
            "sampling_func": su.post_sampler2,
        })
        if plots_only:
            arguments["run_fits"] = False
            arguments["make_plots"] = True
        if run_fits is not None:
            arguments["run_fits"] = run_fits
        if make_plots is not None:
            arguments["make_plots"] = make_plots
        if derived:
            arguments["run_fits"] = False
        results.extend(occurrence_run.run_multiple(**arguments))

    write_run_summary(results, output_dir / RUN_SUMMARY_NAME)
    return results


def write_run_summary(results, path):
    """Save every experiment's full results (paths, chains, convergence).

    The console shows only progress lines and a convergence table; this file
    keeps the details. It is rewritten by each call of :func:`main`.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as stream:
        json.dump(results, stream, indent=2, default=str)
        stream.write("\n")
    print("Full results: {}".format(path))
    return path


if __name__ == "__main__":
    main(run_names=RUNS_TO_DO)
    
    
    
    
    
    
    
    
