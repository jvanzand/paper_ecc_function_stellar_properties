"""Mock recovery summary: recovered Beta parameters vs. the truth.

Usage: python mock_figure.py <work_dir> <output_png>

Reads the mocks written by mock_validation.py in <work_dir> (the folder
names below) and plots each mock's Beta(alpha, beta) posterior median and
68% interval (left) and the median f(e) (right) against the truth.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy import stats

matplotlib.rc_file("/Users/judahvz/research/code/GitHub/occurrence/matplotlibrc")
# The rc file leaves legend text at matplotlib's 10 pt default and doubles
# legend markers; set both here so the legends match the 20-24 pt axes
LEGEND = {"fontsize": 20, "markerscale": 1.0}

TRUE_ALPHA, TRUE_BETA = 1.5, 4.0
# Mock folder -> e uncertainty of its companions
MOCKS = {
    "mock_s005": 0.05, "mock_s005_seed2": 0.05, "mock_s005_seed3": 0.05,
    "mock_s015": 0.15, "mock_s015_seed2": 0.15, "mock_s015_seed3": 0.15,
}
COLORS = {0.05: "#0072B2", 0.15: "#D55E00"}

work_dir, output = Path(sys.argv[1]), Path(sys.argv[2])
figure, (left, right) = plt.subplots(1, 2, figsize=(17, 7))
e = np.linspace(0.001, 0.999, 400)
for name, sigma_e in MOCKS.items():
    chain_path = (work_dir / name / "mtrue/allstars/mock_fit/saved_chains" /
                  "chains_beta_bin0.npz")
    flat = np.load(chain_path)["flat_chains"]
    alpha, beta = flat[:, 1], flat[:, 2]
    a_low, a_mid, a_high = np.percentile(alpha, [16, 50, 84])
    b_low, b_mid, b_high = np.percentile(beta, [16, 50, 84])
    left.errorbar(a_mid, b_mid, xerr=[[a_mid - a_low], [a_high - a_mid]],
                  yerr=[[b_mid - b_low], [b_high - b_mid]], fmt="o",
                  color=COLORS[sigma_e], ms=11, capsize=4, lw=1.5)
    right.plot(e, stats.beta.pdf(e, a_mid, b_mid), color=COLORS[sigma_e],
               lw=2, alpha=0.8)
left.plot(TRUE_ALPHA, TRUE_BETA, "k*", ms=24, zorder=5)
right.plot(e, stats.beta.pdf(e, TRUE_ALPHA, TRUE_BETA), "k--", lw=3)

# Proxy handles: one entry per e uncertainty, without error bars, so the
# legend symbols do not overlap
left.legend(handles=[
    Line2D([], [], color=COLORS[sigma_e], marker="o", ls="none", ms=11,
           label=rf"Mocks, $\sigma_e = {sigma_e}$")
    for sigma_e in sorted(COLORS)
] + [Line2D([], [], color="k", marker="*", ls="none", ms=20, label="Truth")],
    loc="upper left", **LEGEND)
right.legend(handles=[
    Line2D([], [], color=COLORS[sigma_e], lw=2,
           label=rf"Mocks, $\sigma_e = {sigma_e}$")
    for sigma_e in sorted(COLORS)
] + [Line2D([], [], color="k", ls="--", lw=3,
            label=rf"Truth: Beta({TRUE_ALPHA:g}, {TRUE_BETA:g})")],
    loc="upper right", **LEGEND)

left.set_xlabel(r"$\alpha$")
left.set_ylabel(r"$\beta$")
right.set_xlabel("Eccentricity")
right.set_ylabel("Normalized f(e)")
right.set_xlim(0, 1)
right.set_ylim(0, None)
figure.tight_layout()
figure.savefig(output, dpi=150)
print("Saved", output)
