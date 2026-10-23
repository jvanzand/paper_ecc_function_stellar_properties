"""Mock recovery summary: Beta posteriors per mock vs. the truth."""
import json
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats
matplotlib.rc_file("/Users/judahvz/research/code/GitHub/occurrence/matplotlibrc")
runs = {"mock_s005": ("0.05", 1), "mock_s005_seed2": ("0.05", 2),
        "mock_s005_seed3": ("0.05", 3), "mock_s015": ("0.15", 1),
        "mock_s015_seed2": ("0.15", 2), "mock_s015_seed3": ("0.15", 3)}
colors = {"0.05": "#0072B2", "0.15": "#D55E00"}
figure, (left, right) = plt.subplots(1, 2, figsize=(16, 6.5))
e = np.linspace(0.001, 0.999, 400)
for name, (sigma, seed) in runs.items():
    flat = np.load(f"{name}/mtrue/allstars/mock_fit/saved_chains/chains_beta_bin0.npz")["flat_chains"]
    alpha, beta = flat[:, 1], flat[:, 2]
    a, b = np.percentile(alpha, [16, 50, 84]), np.percentile(beta, [16, 50, 84])
    label = rf"$\sigma_e$ = {sigma}" if seed == 1 else None
    left.errorbar(a[1], b[1], xerr=[[a[1]-a[0]], [a[2]-a[1]]], yerr=[[b[1]-b[0]], [b[2]-b[1]]],
                  fmt="o", color=colors[sigma], ms=10, capsize=4, label=label)
    right.plot(e, stats.beta.pdf(e, np.median(alpha), np.median(beta)), color=colors[sigma], lw=2, alpha=0.8, label=label)
left.plot(1.5, 4.0, "k*", ms=22, label="Truth")
left.set_xlabel(r"$\alpha$"); left.set_ylabel(r"$\beta$"); left.legend(fontsize=18, markerscale=0.6, loc="lower right")
right.plot(e, stats.beta.pdf(e, 1.5, 4.0), "k--", lw=3, label="Truth: Beta(1.5, 4)")
right.set_xlabel("Eccentricity"); right.set_ylabel("Normalized f(e)"); right.set_xlim(0, 1); right.legend(fontsize=18)
figure.tight_layout(); figure.savefig("figures/mock_recovery.png", dpi=150)
