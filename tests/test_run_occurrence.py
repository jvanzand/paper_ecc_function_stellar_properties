"""Small regression suite for the paper orchestration layer."""

import json
import os

import pytest

import run_occurrence as run


def test_catalog_and_active_samples_are_valid():
    star_df = run.make_star_df()
    active_run = next(iter(run.RUN_CONFIGURATIONS.values()))

    run.validate_star_df(
        star_df,
        active_run["tier2_list"],
    )
    assert len(star_df) == 719


def test_catalog_cache_preserves_lists_and_source_checksums(tmp_path):
    cache_path = tmp_path / "catalog.json"
    built = run.load_or_make_star_df(cache_path=cache_path, rebuild=True)
    loaded = run.load_or_make_star_df(cache_path=cache_path)

    assert built.equals(loaded)
    assert isinstance(loaded.iloc[0]["comp_list"], list)
    with cache_path.open() as stream:
        cached = json.load(stream)
    assert set(cached["source_sha256"]) == {"stars", "activity", "companions"}


def test_validate_only_does_not_call_occurrence_runner(monkeypatch):
    monkeypatch.setattr(
        run.occurrence_run,
        "run_multiple",
        lambda **kwargs: (_ for _ in ()).throw(AssertionError("runner called")),
    )
    # Runs that reuse fits need their source run's results, which this
    # check does not depend on
    monkeypatch.setattr(run, "_check_reused_fits", lambda *args: None)

    assert run.main(validate_only=True) == []


def test_mcmc_configuration_keeps_samples_after_burnin():
    assert run.RUN_DEFAULTS["nsteps"] > run.RUN_DEFAULTS["burnin"]


def test_runs_to_do_selects_default_runs(monkeypatch):
    monkeypatch.setattr(run, "RUNS_TO_DO", ["stellar2params"])
    assert run._selected_run_names(None) == ["stellar2params"]


def test_runs_to_do_none_selects_every_run(monkeypatch):
    monkeypatch.setattr(run, "RUNS_TO_DO", None)
    assert run._selected_run_names(None) == list(run.RUN_CONFIGURATIONS)


def test_explicit_run_names_override_runs_to_do(monkeypatch):
    monkeypatch.setattr(run, "RUNS_TO_DO", ["stellar2params"])
    assert run._selected_run_names(["ecc_bounds"]) == ["ecc_bounds"]


def test_runs_to_do_rejects_unknown_names(monkeypatch):
    monkeypatch.setattr(run, "RUNS_TO_DO", ["not_a_run"])
    with pytest.raises(ValueError, match="unknown runs"):
        run._selected_run_names(None)


def test_load_local_paths_reads_both_directories(tmp_path):
    path = tmp_path / "local_paths.json"
    path.write_text(json.dumps({
        "recoveries_dir": str(tmp_path / "recoveries"),
        "posteriors_dir": str(tmp_path / "posteriors"),
    }))

    paths = run.load_local_paths(path)

    assert paths == {
        "recoveries_dir": (tmp_path / "recoveries").resolve(),
        "posteriors_dir": (tmp_path / "posteriors").resolve(),
    }


def test_load_local_paths_explains_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError, match="local_paths.example.json"):
        run.load_local_paths(tmp_path / "local_paths.json")


def test_load_local_paths_rejects_missing_entries(tmp_path):
    path = tmp_path / "local_paths.json"
    path.write_text(json.dumps({"recoveries_dir": "/data/recoveries"}))

    with pytest.raises(ValueError, match="posteriors_dir"):
        run.load_local_paths(path)


def test_example_local_paths_lists_every_key():
    with run.LOCAL_PATHS_EXAMPLE.open() as stream:
        example = json.load(stream)
    assert set(example) == set(run.LOCAL_PATH_KEYS)


def _source_and_derived_runs(**derived_overrides):
    source = {
        "tier1_list": ["mtrue", "qtrue"],
        "tier2_list": ["allstars", "highMstar"],
        "tier3_list": ["source_run"],
        "a_edges": [0.1, 10.0],
        "m_edges": [0.4, 13.0, 50.0],
        "run_models_list": ["piecewise", "loglinear"],
        "plot_models_list": ["piecewise"],
        "run_fits": True,
        "plot_corner": True,
    }
    derived = {
        "reuse_fits_from": "source_run",
        "tier1_list": ["mtrue"],
        "tier2_list": ["allstars"],
        "plot_models_list": ["piecewise", "loglinear"],
        "plot_corner": False,
        **derived_overrides,
    }
    return {"source_run": source, "derived_run": derived}


def test_derived_run_inherits_source_and_never_fits(monkeypatch):
    monkeypatch.setattr(
        run, "RUN_CONFIGURATIONS", _source_and_derived_runs()
    )

    configuration = run.resolve_run_configuration("derived_run")

    assert configuration["m_edges"] == [0.4, 13.0, 50.0]
    assert configuration["tier1_list"] == ["mtrue"]
    assert configuration["tier3_list"] == ["derived_run"]
    assert configuration["source_tier3"] == "source_run"
    assert configuration["plot_models_list"] == ["piecewise", "loglinear"]
    assert configuration["plot_corner"] is False
    assert configuration["run_fits"] is False


def test_derived_run_rejects_fit_settings(monkeypatch):
    monkeypatch.setattr(
        run, "RUN_CONFIGURATIONS",
        _source_and_derived_runs(m_edges=[0.4, 50.0]),
    )
    with pytest.raises(ValueError, match="m_edges"):
        run.resolve_run_configuration("derived_run")


def test_derived_run_rejects_tiers_missing_from_source(monkeypatch):
    monkeypatch.setattr(
        run, "RUN_CONFIGURATIONS",
        _source_and_derived_runs(tier2_list=["lowFeH"]),
    )
    with pytest.raises(ValueError, match="lowFeH"):
        run.resolve_run_configuration("derived_run")


def test_derived_run_cannot_reuse_another_derived_run(monkeypatch):
    configurations = _source_and_derived_runs()
    configurations["second_derived"] = {"reuse_fits_from": "derived_run"}
    monkeypatch.setattr(run, "RUN_CONFIGURATIONS", configurations)
    with pytest.raises(ValueError, match="derived run"):
        run.resolve_run_configuration("second_derived")


def _write_source_fits(output_dir, models=("piecewise", "loglinear")):
    source_dir = output_dir / "mtrue" / "allstars" / "source_run"
    (source_dir / "saved_chains").mkdir(parents=True)
    (source_dir / "saved_dicts").mkdir()
    (source_dir / "saved_dicts" / "fit_data.npz").write_bytes(b"fit")
    for model_name in models:
        name = (
            "chains_piecewise.npz" if model_name == "piecewise"
            else "chains_{}_bin0.npz".format(model_name)
        )
        (source_dir / "saved_chains" / name).write_bytes(b"chain")
    return source_dir


def test_link_reused_fits_uses_relative_links(tmp_path, monkeypatch):
    monkeypatch.setattr(
        run, "RUN_CONFIGURATIONS", _source_and_derived_runs()
    )
    configuration = run.resolve_run_configuration("derived_run")
    source_dir = _write_source_fits(tmp_path)

    run._check_reused_fits("derived_run", configuration, tmp_path)
    run.link_reused_fits("derived_run", configuration, tmp_path)
    run.link_reused_fits("derived_run", configuration, tmp_path)

    derived_dir = tmp_path / "mtrue" / "allstars" / "derived_run"
    assert os.readlink(derived_dir / "saved_chains") == (
        "../source_run/saved_chains"
    )
    fit_link = derived_dir / "saved_dicts" / "fit_data.npz"
    assert os.readlink(fit_link) == "../../source_run/saved_dicts/fit_data.npz"
    assert fit_link.read_bytes() == b"fit"
    assert not (derived_dir / "saved_dicts").is_symlink()
    assert (derived_dir / "saved_chains" / "chains_loglinear_bin0.npz").samefile(
        source_dir / "saved_chains" / "chains_loglinear_bin0.npz"
    )


def test_link_reused_fits_refuses_existing_folders(tmp_path, monkeypatch):
    monkeypatch.setattr(
        run, "RUN_CONFIGURATIONS", _source_and_derived_runs()
    )
    configuration = run.resolve_run_configuration("derived_run")
    _write_source_fits(tmp_path)
    (tmp_path / "mtrue" / "allstars" / "derived_run" / "saved_chains").mkdir(
        parents=True
    )
    with pytest.raises(FileExistsError, match="not a link"):
        run.link_reused_fits("derived_run", configuration, tmp_path)


def test_check_reused_fits_reports_missing_model_chains(tmp_path, monkeypatch):
    monkeypatch.setattr(
        run, "RUN_CONFIGURATIONS", _source_and_derived_runs()
    )
    configuration = run.resolve_run_configuration("derived_run")
    _write_source_fits(tmp_path, models=("piecewise",))
    with pytest.raises(FileNotFoundError, match="chains_loglinear_bin"):
        run._check_reused_fits("derived_run", configuration, tmp_path)


def test_configured_derived_runs_resolve():
    for name, entry in run.RUN_CONFIGURATIONS.items():
        if "reuse_fits_from" in entry:
            run.resolve_run_configuration(name)


def test_runs_write_to_a_tier3_folder_named_after_themselves(monkeypatch):
    configurations = _source_and_derived_runs()
    configurations["source_run"]["tier3_list"] = ["other_folder"]
    monkeypatch.setattr(run, "RUN_CONFIGURATIONS", configurations)
    with pytest.raises(ValueError, match="named after itself"):
        run.resolve_run_configuration("source_run")


def test_configured_runs_are_named_after_their_folders():
    for name in run.RUN_CONFIGURATIONS:
        assert run.resolve_run_configuration(name)["tier3_list"] == [name]


def test_piecewise_steps_are_checked_like_the_shared_steps(monkeypatch):
    configuration = dict(run.resolve_run_configuration("stellar2params"))
    run._validate_run_configuration("ok", dict(
        configuration, piecewise_nsteps=15000, piecewise_burnin=5000))
    with pytest.raises(ValueError, match="piecewise_burnin"):
        run._validate_run_configuration("bad", dict(
            configuration, piecewise_nsteps=1000, piecewise_burnin=2000))


def test_run_summary_keeps_full_results(tmp_path, capsys):
    results = [{"tier1": "mtrue", "chains": {"piecewise": tmp_path / "c.npz"},
                "convergence": [{"parameters": [("bin 0", 10.0, 60.0)]}]}]
    path = run.write_run_summary(results, tmp_path / "out" / "summary.json")
    saved = json.loads(path.read_text())
    assert saved[0]["chains"]["piecewise"] == str(tmp_path / "c.npz")
    assert saved[0]["convergence"][0]["parameters"] == [["bin 0", 10.0, 60.0]]
    assert "Full results:" in capsys.readouterr().out


def test_configured_runs_fit_eccentricity_with_a_and_m_marginalized():
    for name, configuration in run.RUN_CONFIGURATIONS.items():
        run._validate_run_configuration(name, configuration)
        spec = run.occurrence_dimensions.parse_dimensions(
            configuration["dimensions"])
        assert spec.fitted == ("e",)
        assert sorted(spec.marginalized) == ["a", "m"]


def test_dimensions_replace_the_original_edges():
    configuration = dict(run.RUN_CONFIGURATIONS["ecc_bounds"],
                         a_edges=[0.1, 30.0])
    with pytest.raises(ValueError, match="replace a_edges"):
        run._validate_run_configuration("bad", configuration)


def test_invalid_dimensions_are_rejected_before_running():
    dimensions = dict(run.ECC_DIMENSIONS)
    dimensions.pop("m")
    configuration = dict(run.RUN_CONFIGURATIONS["ecc_bounds"],
                         dimensions=dimensions)
    with pytest.raises(ValueError, match="bad dimensions"):
        run._validate_run_configuration("bad", configuration)
