"""Shared test fixtures."""

import pytest

import run_occurrence


THREE_PARAMETER_SAMPLES = [
    "highMstarhighFeHhighAct", "highMstarhighFeHlowAct",
    "highMstarlowFeHhighAct", "highMstarlowFeHlowAct",
    "lowMstarhighFeHhighAct", "lowMstarhighFeHlowAct",
    "lowMstarlowFeHhighAct", "lowMstarlowFeHlowAct",
]


@pytest.fixture
def extra_runs(monkeypatch):
    """Add stub runs that exercise machinery the paper's runs do not use.

    ``three_params`` and ``three_params_alt`` fit the three-parameter
    subsets, ``full_only`` fits only the full sample in two Tier 1 spaces,
    and ``standalone`` fits a single sample outside the high/low pairs.
    """
    base = run_occurrence.RUN_CONFIGURATIONS["ecc_bounds"]
    configurations = dict(run_occurrence.RUN_CONFIGURATIONS)
    configurations.update({
        "three_params": dict(base, tier2_list=THREE_PARAMETER_SAMPLES,
                             tier3_list=["three_params"]),
        "three_params_alt": dict(base, tier2_list=THREE_PARAMETER_SAMPLES,
                                 tier3_list=["three_params_alt"]),
        "full_only": dict(base, tier1_list=["mtrue", "qtrue"],
                          tier2_list=["allstars"], tier3_list=["full_only"]),
        "standalone": dict(base, tier2_list=["narrow_cuts"],
                           tier3_list=["standalone"]),
    })
    monkeypatch.setattr(run_occurrence, "RUN_CONFIGURATIONS", configurations)
    return configurations
