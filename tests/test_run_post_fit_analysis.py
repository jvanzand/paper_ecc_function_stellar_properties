"""Checks that post-fit products read the runs they name."""

import pytest

import run_post_fit_analysis as post_fit


def _use_products(monkeypatch, **products):
    monkeypatch.setattr(post_fit, "POST_FIT_PRODUCTS", products)


def test_variables_take_tier3_folders_from_runs(tmp_path, monkeypatch,
                                               extra_runs):
    _use_products(monkeypatch, variables={
        "function": "make_variables",
        "runs": ["ecc_bounds"],
        "three_parameter_runs": ["three_params", "three_params_alt"],
        "tier1_dirs": ["mtrue"],
        "tier2_types": ["allstars", "Mstar"],
    })

    function_name, arguments, folders = post_fit.resolve_product(
        "variables", tmp_path
    )

    assert function_name == "make_variables"
    assert arguments["tier3_dirs"] == ["ecc_bounds"]
    assert arguments["three_parameter_t3"] == [
        "three_params", "three_params_alt",
    ]
    assert arguments["results_dir"] == tmp_path
    assert tmp_path / "mtrue" / "highMstar" / "ecc_bounds" in folders
    assert (tmp_path / "mtrue" / "lowMstarlowFeHhighAct" /
            "three_params") in folders
    assert (tmp_path / "mtrue" / "lowMstarlowFeHhighAct" /
            "three_params_alt") in folders


def test_variables_without_three_parameter_run_skip_it(tmp_path, monkeypatch):
    _use_products(monkeypatch, variables={
        "function": "make_variables", "runs": ["ecc_bounds"],
        "tier1_dirs": ["mtrue"], "tier2_types": ["allstars"],
    })
    _, arguments, _ = post_fit.resolve_product("variables", tmp_path)
    assert arguments["three_parameter_t3"] == []


def test_subset_tables_read_the_run_tier2_folders(tmp_path, monkeypatch):
    _use_products(monkeypatch, tables={
        "function": "make_two_parameter_tables",
        "run": "stellar2params",
        "t1": "mtrue",
    })

    _, arguments, folders = post_fit.resolve_product("tables", tmp_path)

    assert arguments["t3"] == "stellar2params"
    assert arguments["tier2_dirs"] == [
        "highMstarhighFeH", "highMstarlowFeH",
        "lowMstarhighFeH", "lowMstarlowFeH",
    ]
    assert len(folders) == 4


def test_products_reject_tiers_their_run_lacks(tmp_path, monkeypatch):
    _use_products(monkeypatch, table={
        "function": "make_parameter_table", "run": "stellar2params",
        "t1": "mtrue", "t2": "allstars", "models": ["rayleigh"],
    })
    with pytest.raises(ValueError, match="allstars"):
        post_fit.resolve_product("table", tmp_path)


def test_products_cannot_set_tier3_directly(tmp_path, monkeypatch):
    _use_products(monkeypatch, table={
        "function": "make_two_parameter_tables", "run": "stellar2params",
        "t3": "stellar2params",
    })
    with pytest.raises(ValueError, match="run configuration"):
        post_fit.resolve_product("table", tmp_path)


def test_products_reject_unknown_runs(tmp_path, monkeypatch):
    _use_products(monkeypatch, table={
        "function": "make_two_parameter_tables", "run": "stellar_2params",
    })
    with pytest.raises(ValueError, match="unknown run"):
        post_fit.resolve_product("table", tmp_path)


def test_main_requires_results_before_making_anything(tmp_path, monkeypatch):
    _use_products(monkeypatch, tables={
        "function": "make_two_parameter_tables",
        "run": "stellar2params", "t1": "mtrue",
    })
    monkeypatch.setattr(
        post_fit.pfa, "make_two_parameter_tables",
        lambda **kwargs: pytest.fail("made a product without its results"),
    )
    with pytest.raises(FileNotFoundError, match="run_occurrence.py first"):
        post_fit.main(["tables"], results_dir=tmp_path)


def test_main_passes_resolved_arguments(tmp_path, monkeypatch):
    _use_products(monkeypatch, tables={
        "function": "make_two_parameter_tables",
        "run": "stellar2params", "t1": "mtrue",
        "use_latex_variables": False,
    })
    for tier2 in ("highMstarhighFeH", "highMstarlowFeH",
                  "lowMstarhighFeH", "lowMstarlowFeH"):
        (tmp_path / "mtrue" / tier2 / "stellar2params").mkdir(parents=True)
    calls = []
    monkeypatch.setattr(
        post_fit.pfa, "make_two_parameter_tables",
        lambda **kwargs: calls.append(kwargs) or "made",
    )

    outputs = post_fit.main(["tables"], results_dir=tmp_path)

    assert outputs == {"tables": "made"}
    assert calls[0]["t3"] == "stellar2params"
    assert calls[0]["use_latex_variables"] is False


def test_configured_products_resolve(tmp_path):
    for name in post_fit.POST_FIT_PRODUCTS:
        post_fit.resolve_product(name, tmp_path)
    post_fit._selected_product_names(None)


def test_variables_include_each_runs_own_samples(tmp_path, monkeypatch,
                                                extra_runs):
    _use_products(monkeypatch, variables={
        "function": "make_variables",
        "runs": ["ecc_bounds", "full_only"],
        "tier1_dirs": ["mtrue", "msini"],
        "tier2_types": ["allstars", "Mstar"],
    })

    _, arguments, folders = post_fit.resolve_product("variables", tmp_path)

    assert arguments["tier3_dirs"] == ["ecc_bounds", "full_only"]
    full_only = sorted(folder.relative_to(tmp_path) for folder in folders
                       if folder.name == "full_only")
    assert full_only == [
        post_fit.Path("msini/allstars/full_only"),
        post_fit.Path("mtrue/allstars/full_only"),
    ]
    # ecc_bounds has only mtrue
    assert len(folders) == 3 + 2


def test_variables_need_every_run_to_contribute(tmp_path, monkeypatch,
                                                extra_runs):
    _use_products(monkeypatch, variables={
        "function": "make_variables",
        "runs": ["ecc_bounds", "full_only"],
        "tier1_dirs": ["mtrue"], "tier2_types": ["FeH"],
    })
    with pytest.raises(ValueError, match="none of run 'full_only'"):
        post_fit.resolve_product("variables", tmp_path)


def test_variables_reject_samples_no_run_has(tmp_path, monkeypatch,
                                            extra_runs):
    _use_products(monkeypatch, variables={
        "function": "make_variables",
        "runs": ["ecc_bounds", "full_only"],
        "tier1_dirs": ["qsini"], "tier2_types": ["allstars"],
    })
    with pytest.raises(ValueError, match="qsini"):
        post_fit.resolve_product("variables", tmp_path)


def test_variables_use_runs_and_tables_use_run(tmp_path, monkeypatch):
    _use_products(monkeypatch,
                  variables={"function": "make_variables",
                             "run": "ecc_bounds"},
                  table={"function": "make_two_parameter_tables",
                         "runs": ["stellar2params"]})
    for name in ("variables", "table"):
        with pytest.raises(ValueError, match="must name what it reads"):
            post_fit.resolve_product(name, tmp_path)


def _cdf_product(**overrides):
    product = {
        "function": "plot_model_cdf_comparison",
        "models": ["beta", "rayleigh"],
        "sample_pairs": ["Mstar", "Act"],
        "run": "ecc_bounds",
        "tier1": "mtrue",
        "credible": 0.95,
    }
    product.update(overrides)
    return product


def test_cdf_rows_are_high_and_low_sample_pairs(tmp_path, monkeypatch):
    _use_products(monkeypatch, cdf=_cdf_product(labels={"highAct": "Young"}))

    function_name, arguments, paths = post_fit.resolve_product("cdf",
                                                               tmp_path)

    assert function_name == "plot_model_cdf_comparison"
    # The file name does not follow the product name.
    assert arguments["name"] == "cdf_comparison"
    assert "filename" not in arguments
    assert arguments["models"] == ["beta", "rayleigh"]
    assert arguments["credible"] == 0.95
    assert [[curve["t2"] for curve in row] for row in arguments["rows"]] == [
        ["highMstar", "lowMstar"], ["highAct", "lowAct"],
    ]
    assert arguments["rows"][0][0] == {
        "label": post_fit.tier2_df_cuts_dict["highMstar"][1],
        "t1": "mtrue", "t2": "highMstar", "t3": "ecc_bounds",
    }
    assert [curve["label"] for curve in arguments["rows"][1]] == [
        "Young", post_fit.tier2_df_cuts_dict["lowAct"][1],
    ]
    for key in ("sample_pairs", "run", "tier1", "labels"):
        assert key not in arguments


def test_cdf_products_check_every_chain_before_plotting(tmp_path,
                                                        monkeypatch):
    _use_products(monkeypatch, cdf=_cdf_product())
    _, _, paths = post_fit.resolve_product("cdf", tmp_path)
    chains = {path.relative_to(tmp_path) for path in paths
              if path.suffix == ".npz"}
    assert len(chains) == 4*2
    assert tmp_path / "mtrue" / "allstars" / "ecc_bounds" / "plots" in paths
    assert post_fit.Path(
        "mtrue/lowAct/ecc_bounds/saved_chains/chains_rayleigh_bin0.npz"
    ) in chains

    for path in paths:
        if path.suffix != ".npz":
            path.mkdir(parents=True, exist_ok=True)
        elif "rayleigh" not in path.name:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.touch()
    monkeypatch.setattr(
        post_fit.pfa, "plot_model_cdf_comparison",
        lambda **kwargs: pytest.fail("plotted without every chain"),
    )
    with pytest.raises(FileNotFoundError, match="chains_rayleigh_bin0"):
        post_fit.main(["cdf"], results_dir=tmp_path)


def test_cdf_sample_pairs_must_belong_to_the_run(tmp_path, monkeypatch,
                                                 extra_runs):
    _use_products(monkeypatch, cdf=_cdf_product(run="full_only"))
    with pytest.raises(ValueError, match="highMstar"):
        post_fit.resolve_product("cdf", tmp_path)


def test_cdf_labels_must_name_plotted_samples(tmp_path, monkeypatch):
    _use_products(monkeypatch, cdf=_cdf_product(labels={"highFeH": "x"}))
    with pytest.raises(ValueError, match="highFeH"):
        post_fit.resolve_product("cdf", tmp_path)


@pytest.mark.parametrize("change", [
    {"name": "other"}, {"rows": []},
])
def test_cdf_products_cannot_set_derived_arguments(tmp_path, monkeypatch,
                                                   change):
    _use_products(monkeypatch, cdf=_cdf_product(**change))
    with pytest.raises(ValueError, match="sample_pairs"):
        post_fit.resolve_product("cdf", tmp_path)


def test_cdf_products_need_models_and_pairs(tmp_path, monkeypatch):
    product = _cdf_product()
    del product["sample_pairs"]
    _use_products(monkeypatch, cdf=product)
    with pytest.raises(ValueError, match="sample_pairs"):
        post_fit.resolve_product("cdf", tmp_path)


def test_cdf_products_need_a_full_sample_folder_to_save_into(tmp_path,
                                                             monkeypatch):
    import run_occurrence
    configurations = dict(run_occurrence.RUN_CONFIGURATIONS)
    configurations["mass_only"] = dict(
        configurations["ecc_bounds"], tier2_list=["highMstar", "lowMstar"],
        tier3_list=["mass_only"],
    )
    monkeypatch.setattr(run_occurrence, "RUN_CONFIGURATIONS", configurations)
    _use_products(monkeypatch, cdf=_cdf_product(run="mass_only",
                                                sample_pairs=["Mstar"]))
    with pytest.raises(ValueError, match="allstars"):
        post_fit.resolve_product("cdf", tmp_path)


def test_variables_include_standalone_samples(tmp_path, monkeypatch,
                                             extra_runs):
    _use_products(monkeypatch, variables={
        "function": "make_variables",
        "runs": ["ecc_bounds", "standalone"],
        "tier1_dirs": ["mtrue"], "tier2_types": ["allstars"],
        "standalone_tier2_dirs": ["narrow_cuts"],
        "tier3_stack_dims": {"standalone": "a"},
    })
    _, arguments, folders = post_fit.resolve_product("variables", tmp_path)
    assert tmp_path / "mtrue" / "narrow_cuts" / "standalone" in folders
    assert tmp_path / "mtrue" / "allstars" / "ecc_bounds" in folders
    assert arguments["standalone_tier2_dirs"] == ["narrow_cuts"]


def test_cdf_rows_are_labeled_by_stellar_parameter(tmp_path, monkeypatch):
    _use_products(monkeypatch, cdf=_cdf_product(row_labels={"Act": "Stellar age"}))
    _, arguments, _ = post_fit.resolve_product("cdf", tmp_path)
    assert arguments["row_labels"] == ["Stellar Mass", "Stellar age"]
    _use_products(monkeypatch, cdf=_cdf_product(row_labels={"FeH": "x"}))
    with pytest.raises(ValueError, match="FeH"):
        post_fit.resolve_product("cdf", tmp_path)


def test_variables_pass_two_parameter_runs_to_occurrence():
    _, arguments, folders = post_fit.resolve_product("variables")
    assert arguments["two_parameter_t3"] == ["stellar2params"]
    assert "two_parameter_runs" not in arguments
    assert any(folder.name == "stellar2params" for folder in folders)


def test_one_parameter_table_reads_the_one_dimensional_fits():
    function, arguments, folders = post_fit.resolve_product("one_param_table")
    assert function == "make_one_parameter_table"
    assert arguments["t3"] == "ecc_bounds"
    assert arguments["tier2_types"] == ["Mstar", "FeH"]
    assert sorted(folder.parent.name for folder in folders) == [
        "highFeH", "highMstar", "lowFeH", "lowMstar",
    ]


def test_cdf_filename_sets_the_figure_name(tmp_path, monkeypatch):
    _use_products(monkeypatch, cdf=_cdf_product(filename="cdf_q.png"),
                  bare=_cdf_product(filename="cdf_bare"))
    assert post_fit.resolve_product("cdf", tmp_path)[1]["name"] == "cdf_q"
    assert post_fit.resolve_product("bare", tmp_path)[1]["name"] == "cdf_bare"
    for bad in ("plots/cdf.png", "cdf.pdf", ""):
        _use_products(monkeypatch, cdf=_cdf_product(filename=bad))
        with pytest.raises(ValueError, match="filename"):
            post_fit.resolve_product("cdf", tmp_path)


def test_appendix_table_takes_models_from_other_runs(tmp_path, monkeypatch,
                                                    extra_runs):
    _use_products(monkeypatch, table={
        "function": "make_appendix_parameter_table", "run": "ecc_bounds",
        "tier1_dirs": ["mtrue"], "tier2_types": ["allstars", "Mstar"],
        "model_runs": {"escarpment": "full_only"},
    })
    _, arguments, folders = post_fit.resolve_product("table", tmp_path)
    assert arguments["model_tier3"] == {"escarpment": "full_only"}
    assert "model_runs" not in arguments
    borrowed = sorted((folder.parent.parent.name, folder.parent.name)
                      for folder in folders if folder.name == "full_only")
    assert borrowed == [("mtrue", "allstars")]


def test_main_prints_one_line_per_product_and_collects_warnings(
        tmp_path, monkeypatch, capsys):
    _use_products(monkeypatch, variables={
        "function": "make_variables", "runs": ["stellar2params"],
        "tier1_dirs": ["mtrue"], "tier2_types": [],
        "standalone_tier2_dirs": ["highMstarhighFeH"],
    })
    (tmp_path / "mtrue" / "highMstarhighFeH" / "stellar2params").mkdir(
        parents=True)
    tables = tmp_path / "paper_tables"
    tables.mkdir()
    (tables / "variables.tex").write_text(
        "\\newcommand{\\McX}{\\ensuremath{0.30}}\n"
        "\\newcommand{\\McGone}{\\ensuremath{1}}\n"
    )

    def fake_make_variables(**kwargs):
        print("post_fit_analysis.make_variables: skipped experiments: x")
        print("Corner plot saved to: somewhere.png")
        path = kwargs["results_dir"] / "paper_tables" / "variables.tex"
        path.write_text("\\newcommand{\\McX}{\\ensuremath{0.45}}\n"
                        "\\newcommand{\\McNew}{\\ensuremath{2}}\n")
        return path

    monkeypatch.setattr(post_fit.pfa, "make_variables", fake_make_variables)
    post_fit.main(["variables"], results_dir=tmp_path)
    output = capsys.readouterr().out

    product_line = next(line for line in output.splitlines()
                        if line.startswith("[1/1] variables"))
    assert product_line.endswith("  paper_tables/variables.tex (2 commands)")
    assert "saved to" not in output
    assert ("variables.tex: 2 commands; 1 new, 1 removed, 1 values changed"
            in output)
    assert "\\McX: 0.30 -> 0.45" in output
    assert "Warnings:\n  variables: skipped experiments: x" in output
    assert "Made 1 products in" in output and "(1 warning)" in output
    summary = post_fit.json.loads(
        (tables / post_fit.POST_FIT_SUMMARY_NAME).read_text())
    assert summary["variables_changes"]["added"] == ["McNew"]
    assert summary["outputs"]["variables"].endswith("variables.tex")


def test_outputs_are_described_relative_to_the_results(tmp_path):
    plots = tmp_path / "mtrue" / "allstars" / "run" / "plots"
    described = post_fit._describe_output(
        {"a": plots / "one.png", "b": {"c": plots / "two.png"}}, tmp_path)
    assert described == "2 files in mtrue/allstars/run/plots/"
    assert post_fit._describe_output(None, tmp_path) == "done"
