"""Gather every figure, table, and variables file the paper uses.

Running this module rebuilds ``paper_items/`` with the same layout as the
LaTeX project: tables and ``variables.tex`` at the top level and figures in
``Figures/``. Copy the folder's contents over the LaTeX project to update the
paper, for example::

    rsync -a paper_items/ /path/to/latex_project/

Figures come from experiment ``plots`` folders written by
``run_occurrence.py``, and whole-sample figures (the stellar histograms) from
``results/paper_figures/``; tables and variables come from
``results/paper_tables/``. Both results folders are written by
``run_post_fit_analysis.py``. Every source is checked before
anything is copied, so a collection is never left half-updated, and every
results macro a collected table uses must be defined in ``variables.tex``.
"""

import re
import shutil
from pathlib import Path

from run_occurrence import PROJECT_DIR, RESULTS_DIR, resolve_run_configuration

from occurrence.post_fit_analysis import (
    PAPER_FIGURES_DIRNAME, PAPER_TABLES_DIRNAME,
)


PAPER_ITEMS_DIR = PROJECT_DIR / "paper_items"
FIGURES_DIRNAME = "Figures"

# Short labels used in figure names, keyed by the file each experiment's
# ``plots`` folder holds. A figure is named "<label>_<tier1>_<tier2>_<run>.png".
PLOT_LABELS = {
    "occurrence_ORD.png": "ORD",
    "companions_by_mstar.png": "companions_Mstar",
    "companions_by_feh.png": "companions_FeH",
    "companions_by_age.png": "companions_Age",
}


def _figures(plot, run, tier1, tier2_list):
    return [
        {"plot": plot, "run": run, "tier1": tier1, "tier2": tier2}
        for tier2 in tier2_list
    ]


# Figures the paper includes, each from one experiment's plots folder. An
# entry may set "name" to override the generated "<label>_<tier1>_<tier2>_<run>"
# file name, e.g. for figures that combine several experiments.
PAPER_FIGURES = [
    # Results: Eccentricity functions for the full sample and each stellar subset
    *_figures("occurrence_ORD.png", "ecc_bounds", "mtrue", [
        "allstars", "highMstar", "lowMstar", "highFeH", "lowFeH",
        "highAct", "lowAct",
    ]),
    # Sample: companions colored by host properties
    *_figures("companions_by_mstar.png", "ecc_bounds", "mtrue", ["allstars"]),
    *_figures("companions_by_feh.png", "ecc_bounds", "mtrue", ["allstars"]),
    *_figures("companions_by_age.png", "ecc_bounds", "mtrue", ["allstars"]),
    # Discussion: normalized model CDFs for every stellar sample pair, made by
    # run_post_fit_analysis.py's cdf_comparison product and saved with the
    # catalog plots
    {"plot": "cdf_comparison.png", "run": "ecc_bounds", "tier1": "mtrue",
     "tier2": "allstars", "name": "CDF_comparison_mtrue.png"},
]

# Figures of the whole sample, not tied to any run: file in
# results/paper_figures/ (written by run_post_fit_analysis.py's
# stellar_histograms product) mapped to its name in Figures/.
PAPER_SAMPLE_FIGURES = {
    "stellar_hist_mstar.png": "hist_Mstar.png",
    "stellar_hist_feh.png": "hist_FeH.png",
    "stellar_hist_age.png": "hist_Age.png",
}

# Tables and variables the paper inputs: file in results/paper_tables/ mapped
# to the name the LaTeX project uses.
PAPER_TABLES = {
    "variables.tex": "variables.tex",
    "model_params_mtrue_allstars_ecc_bounds.tex": "model_params_table.tex",
    "model_params_appendix_ecc_bounds.tex": "model_params_table_appendix.tex",
    "one_parameter_OR_mtrue_ecc_bounds.tex": "one_param_OR_table.tex",
    "two_parameter_OR_mtrue_stellar2params.tex": "two_param_OR_table.tex",
    "two_parameter_OR_reordered_mtrue_stellar2params.tex":
        "two_param_OR_table_reordered.tex",
}


VARIABLES_NAME = "variables.tex"
# Results macros start with a Tier 1 prefix (Mc, Msini, Q, Qsini), e.g.
# \McAllstarsPaperBoundsNstars. Lowercase continuations are matched too, so
# stale names such as \McallstarsNstars are caught.
RESULTS_MACRO = re.compile(r"\\((?:Mc|Msini|Qsini|Q)[A-Za-z]+)")
DEFINED_MACRO = re.compile(r"\\newcommand\{\\([A-Za-z]+)\}")


def undefined_table_macros(plan):
    """Return {table: [macros]} for results macros variables.tex lacks."""
    defined = set(DEFINED_MACRO.findall(
        plan[Path(VARIABLES_NAME)].read_text()
    ))
    undefined = {}
    for destination, source in sorted(plan.items()):
        if destination.suffix != ".tex" or destination.name == VARIABLES_NAME:
            continue
        used = set(RESULTS_MACRO.findall(source.read_text()))
        if used - defined:
            undefined[str(destination)] = sorted(used - defined)
    return undefined


def figure_name(figure):
    """Return the paper filename for a figure entry.

    This is the entry's "name" when given, and otherwise the uniform
    "<label>_<tier1>_<tier2>_<run>.png".
    """
    if "name" in figure:
        return figure["name"]
    if figure["plot"] not in PLOT_LABELS:
        raise ValueError("no label for plot {!r}; add it to PLOT_LABELS or "
                         "give the figure a name".format(figure["plot"]))
    return "{}_{}_{}_{}.png".format(
        PLOT_LABELS[figure["plot"]], figure["tier1"], figure["tier2"],
        figure["run"],
    )


def _figure_source(figure, results_dir):
    """Return the plot file for a figure, checking it against its run."""
    try:
        configuration = resolve_run_configuration(figure["run"])
    except KeyError:
        raise ValueError("figure reads unknown run {!r}".format(
            figure["run"]
        )) from None
    for tier_key, tier in (("tier1_list", figure["tier1"]),
                           ("tier2_list", figure["tier2"])):
        if tier not in configuration[tier_key]:
            raise ValueError("run {!r} has no {} {!r}".format(
                figure["run"], tier_key[:-5], tier
            ))
    tier3 = configuration["tier3_list"][0]
    return (Path(results_dir) / figure["tier1"] / figure["tier2"] / tier3 /
            "plots" / figure["plot"])


def collection_plan(results_dir=RESULTS_DIR):
    """Map each destination path (relative to paper_items/) to its source."""
    plan = {}
    for figure in PAPER_FIGURES:
        destination = Path(FIGURES_DIRNAME) / figure_name(figure)
        if destination in plan:
            raise ValueError("two figures are named {}".format(destination))
        plan[destination] = _figure_source(figure, results_dir)
    sample_figures_dir = Path(results_dir) / PAPER_FIGURES_DIRNAME
    for source_name, destination_name in PAPER_SAMPLE_FIGURES.items():
        destination = Path(FIGURES_DIRNAME) / destination_name
        if destination in plan:
            raise ValueError("two figures are named {}".format(destination))
        plan[destination] = sample_figures_dir / source_name
    tables_dir = Path(results_dir) / PAPER_TABLES_DIRNAME
    for source_name, destination_name in PAPER_TABLES.items():
        destination = Path(destination_name)
        if destination in plan:
            raise ValueError("two items are named {}".format(destination))
        plan[destination] = tables_dir / source_name
    return plan


def main(results_dir=RESULTS_DIR, paper_items_dir=PAPER_ITEMS_DIR,
         dry_run=False):
    """Rebuild ``paper_items/`` from the current results.

    The folder is emptied first so it holds exactly the current paper items.
    Nothing is removed or copied unless every source exists. ``dry_run``
    stops after that check.
    """
    plan = collection_plan(results_dir)
    missing = sorted(str(source) for source in plan.values()
                     if not source.is_file())
    if missing:
        raise FileNotFoundError(
            "{} paper item sources are missing; run run_occurrence.py and "
            "run_post_fit_analysis.py for them first:\n  {}".format(
                len(missing), "\n  ".join(missing)
            )
        )
    undefined = undefined_table_macros(plan)
    if undefined:
        raise ValueError(
            "tables use macros that variables.tex does not define; rerun "
            "run_post_fit_analysis.py for all of them together:\n  {}".format(
                "\n  ".join("{}: {}".format(table, ", ".join(macros))
                            for table, macros in undefined.items())
            )
        )
    print("Collecting {} items into {}".format(len(plan), paper_items_dir))
    if dry_run:
        return plan

    paper_items_dir = Path(paper_items_dir)
    if paper_items_dir.name != PAPER_ITEMS_DIR.name:
        raise ValueError("refusing to empty {}; the collection folder must be "
                         "named {}".format(paper_items_dir, PAPER_ITEMS_DIR.name))
    if paper_items_dir.exists():
        shutil.rmtree(paper_items_dir)
    for destination, source in sorted(plan.items()):
        target = paper_items_dir / destination
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        print("  {} <- {}".format(destination, source.relative_to(
            Path(results_dir)
        )))
    return {paper_items_dir / destination: source
            for destination, source in plan.items()}


if __name__ == "__main__":
    main()
