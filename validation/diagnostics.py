"""Diagnostic figures for the eccentricity work.

Usage: python diagnostics.py <ecc_run_dir> <figure_dir> [<reference_run_dir>]

1. sliced_average_maps.png: the all-star average completeness in (a, m) for
   each eccentricity slice.
2. completeness_vs_e.png: completeness marginalized over a and m (default
   population distributions) as a function of e, which is what the
   eccentricity fit's expected-detection term uses, for all stars.
3. slice_counts.png: injections per star in each e slice (after the 10
   inclination draws), with the share of NaN cells in each star's slice map.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import numpy as np

GITHUB = Path("/Users/judahvz/research/code/GitHub")
sys.path.insert(0, str(GITHUB))
from occurrence import completeness_utils as cu  # noqa: E402
from occurrence import dimensions as dm  # noqa: E402

run_dir, figure_dir = Path(sys.argv[1]), Path(sys.argv[2])
figure_dir.mkdir(parents=True, exist_ok=True)
matplotlib.rc_file(GITHUB / "occurrence" / "matplotlibrc")

average = cu.load_sliced(run_dir / "mtrue/allstars/avg_map_eslices")
edges = average.e_edges

# 1. Sliced average maps
figure, axes = plt.subplots(2, 4, figsize=(20, 9), sharex=True, sharey=True)
for index, axis in enumerate(axes.flat):
    mesh = axis.pcolormesh(average.xgrid, average.ygrid, average.zgrids[index],
                           vmin=0, vmax=1, cmap="viridis", shading="auto")
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlim(0.1, 30)
    axis.set_ylim(0.1, 100)
    axis.set_title(f"e = {edges[index]:g}-{edges[index + 1]:g}")
    axis.contour(average.xgrid, average.ygrid, average.zgrids[index], [0.5],
                 colors="w", linewidths=1.5)
for axis in axes[-1]:
    axis.set_xlabel("a [AU]")
for axis in axes[:, 0]:
    axis.set_ylabel(r"$M_c$ [$M_{Jup}$]")
figure.colorbar(mesh, ax=axes, label="Average completeness (white: 50%)",
                shrink=0.8)
figure.savefig(figure_dir / "sliced_average_maps.png", dpi=150,
               bbox_inches="tight")
plt.close(figure)

# 2. Completeness vs e, marginalized over a and m
spec = dm.parse_dimensions({
    "e": {"role": "fit", "edges": [0.0, 1.0]},
    "a": {"role": "marginalize", "bounds": (0.1, 30.0)},
    "m": {"role": "marginalize", "bounds": (0.4, 50.0)},
})
engine = cu.engine_completeness(spec, average)
e = np.linspace(0, 0.999, 400)
figure, axis = plt.subplots(figsize=(8, 5))
axis.plot(e, engine((10**e, np.full_like(e, 2.0))), color="k", lw=2.5,
          label="All stars, e slices")
# The original all-e average map gives one value, the same weighting of a
# and m with e averaged over the injections' Beta distribution
legacy_dir = run_dir / "mtrue/allstars"
original = Path(sys.argv[3]) / "mtrue/allstars/avg_map" if len(sys.argv) > 3 else None
if original is not None and (original / "parent_zgrid.npy").exists():
    one_slice = cu.SlicedCompleteness(
        np.load(original / "parent_xgrid.npy"),
        np.load(original / "parent_ygrid.npy"),
        np.load(original / "parent_zgrid.npy")[None], (0.0, 1.0),
    )
    value = float(cu.engine_completeness(spec, one_slice)((10**0.3, 2.0)))
    axis.axhline(value, color="tab:blue", ls="--", lw=2,
                 label="Original maps (e averaged over injections)")
for boundary in edges[1:-1]:
    axis.axvline(boundary, color="0.8", lw=0.8, zorder=0)
axis.set_xlabel("Eccentricity")
axis.set_ylabel("Mean completeness\nover a and M")
axis.set_xlim(0, 1)
axis.set_ylim(0, None)
axis.legend()
figure.tight_layout()
figure.savefig(figure_dir / "completeness_vs_e.png", dpi=150)
plt.close(figure)

# 3. Injections per slice and NaN share of the stars' slice maps
counts_file = np.load(run_dir / "mtrue/saved_maps_mtrue_eslices/slice_counts.npz")
counts = counts_file["counts"]
nan_share = []
for star in counts_file["stars"]:
    sliced = cu.load_sliced(run_dir / "mtrue/saved_maps_mtrue_eslices" / str(star))
    nan_share.append(np.isnan(sliced.zgrids).mean(axis=(1, 2)))
nan_share = np.array(nan_share)
labels = [f"{edges[i]:g}-{edges[i + 1]:g}" for i in range(len(edges) - 1)]
figure, (left, right) = plt.subplots(1, 2, figsize=(16, 6))
left.boxplot(counts, labels=labels)
left.set_yscale("log")
left.set_xlabel("e slice")
left.set_ylabel("Injections per star")
right.boxplot(nan_share, labels=labels)
right.set_xlabel("e slice")
right.set_ylabel("Share of cells\nwith no estimate")
for axis in (left, right):
    plt.setp(axis.get_xticklabels(), rotation=45, ha="right")
figure.tight_layout()
figure.savefig(figure_dir / "slice_counts.png", dpi=150)
plt.close(figure)
print("Saved diagnostics to", figure_dir)
