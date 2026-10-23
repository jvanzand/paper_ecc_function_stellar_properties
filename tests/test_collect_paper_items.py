"""Checks that paper items are named uniformly and collected safely."""

import pytest

import collect_paper_items as collect


def _write_sources(results_dir, plan):
    for source in plan.values():
        source.parent.mkdir(parents=True, exist_ok=True)
        source.write_text(source.name)


def test_figures_follow_the_naming_scheme():
    figure = {"plot": "occurrence_ORD.png", "run": "ecc_bounds",
              "tier1": "mtrue", "tier2": "highFeH"}
    assert collect.figure_name(figure) == "ORD_mtrue_highFeH_ecc_bounds.png"


def test_plan_reads_figures_from_each_run_folder(tmp_path):
    plan = collect.collection_plan(tmp_path)

    old_hosts = plan[collect.Path(
        "Figures/ORD_mtrue_lowAct_ecc_bounds.png"
    )]
    assert old_hosts == (tmp_path / "mtrue" / "lowAct" / "ecc_bounds" /
                         "plots" / "occurrence_ORD.png")
    assert plan[collect.Path("two_param_OR_table.tex")] == (
        tmp_path / "paper_tables" /
        "two_parameter_OR_mtrue_stellar2params.tex"
    )
    assert plan[collect.Path("one_param_OR_table.tex")] == (
        tmp_path / "paper_tables" / "one_parameter_OR_mtrue_ecc_bounds.tex"
    )


def test_figures_must_belong_to_their_run(tmp_path, monkeypatch):
    monkeypatch.setattr(collect, "PAPER_FIGURES", [{
        "plot": "occurrence_ORD.png", "run": "stellar2params",
        "tier1": "mtrue", "tier2": "allstars",
    }])
    with pytest.raises(ValueError, match="allstars"):
        collect.collection_plan(tmp_path)


def test_main_copies_nothing_when_a_source_is_missing(tmp_path):
    results_dir = tmp_path / "results"
    paper_items = tmp_path / "paper_items"
    paper_items.mkdir()
    (paper_items / "keep.txt").write_text("previous collection")
    plan = collect.collection_plan(results_dir)
    _write_sources(results_dir, plan)
    next(iter(plan.values())).unlink()

    with pytest.raises(FileNotFoundError, match="1 paper item sources"):
        collect.main(results_dir, paper_items)
    assert (paper_items / "keep.txt").exists()


def test_main_rebuilds_the_folder_in_the_latex_layout(tmp_path):
    results_dir = tmp_path / "results"
    paper_items = tmp_path / "paper_items"
    paper_items.mkdir()
    (paper_items / "stale.png").write_text("old")
    plan = collect.collection_plan(results_dir)
    _write_sources(results_dir, plan)

    collect.main(results_dir, paper_items)

    assert not (paper_items / "stale.png").exists()
    assert (paper_items / "variables.tex").read_text() == "variables.tex"
    assert (paper_items / "Figures" /
            "ORD_mtrue_lowMstar_ecc_bounds.png").is_file()
    collected = {path.relative_to(paper_items)
                 for path in paper_items.rglob("*") if path.is_file()}
    assert collected == set(plan)


def test_main_only_empties_a_paper_items_folder(tmp_path):
    results_dir = tmp_path / "results"
    _write_sources(results_dir, collect.collection_plan(results_dir))
    with pytest.raises(ValueError, match="refusing"):
        collect.main(results_dir, tmp_path / "Figures")


def test_main_rejects_tables_using_undefined_macros(tmp_path):
    results_dir = tmp_path / "results"
    paper_items = tmp_path / "paper_items"
    plan = collect.collection_plan(results_dir)
    _write_sources(results_dir, plan)
    plan[collect.Path("variables.tex")].write_text(
        r"\newcommand{\McAllstarsSmaBoundsNstars}{\ensuremath{719}}"
    )
    plan[collect.Path("model_params_table.tex")].write_text(
        r"\tablecaption{x} \McAllstarsSmaBoundsNstars & \McallstarsNstars"
    )

    with pytest.raises(ValueError, match="McallstarsNstars"):
        collect.main(results_dir, paper_items)
    assert not paper_items.exists()


def test_tables_may_use_any_defined_macro(tmp_path):
    plan = collect.collection_plan(tmp_path)
    _write_sources(tmp_path, plan)
    plan[collect.Path("variables.tex")].write_text(
        r"\newcommand{\McLowMstarSmaBoundsNeff}{\ensuremath{3.1}}"
    )
    plan[collect.Path("two_param_OR_table.tex")].write_text(
        r"\McLowMstarSmaBoundsNeff"
    )
    assert collect.undefined_table_macros(plan) == {}


def test_figures_may_set_their_own_name(tmp_path):
    plan = collect.collection_plan(tmp_path)
    assert plan[collect.Path("Figures/CDF_comparison_mtrue.png")] == (
        tmp_path / "mtrue" / "allstars" / "ecc_bounds" / "plots" /
        "cdf_comparison.png"
    )


def test_unnamed_figures_need_a_plot_label(tmp_path, monkeypatch):
    monkeypatch.setattr(collect, "PAPER_FIGURES", [{
        "plot": "unlabeled.png", "run": "ecc_bounds",
        "tier1": "mtrue", "tier2": "allstars",
    }])
    with pytest.raises(ValueError, match="PLOT_LABELS"):
        collect.collection_plan(tmp_path)
