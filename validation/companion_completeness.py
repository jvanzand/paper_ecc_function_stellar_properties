"""Check the completeness one companion's posterior draws get, slice by slice.

Usage: python companion_completeness.py <tier2> <companion> [<results_dir>]

For a companion of an eccentricity run (e.g. ``highMstar 120066_0``), and
for each completeness-map e slice its draws fall in (inside the marginalized
a and m ranges), prints:

* how many draws have zero, missing, or positive completeness in the host's
  own sliced map, the subsample's average map, and the all-star average map,
  and the value the catalog finally stores (host, falling back to the
  subsample average where the host is missing or zero);
* at the draws' median (a, M): the subsample average, how many of the
  subsample's stars have an estimate there and how many of those are
  positive, and the host's injections in the moving-average box there.

It also saves ``<companion>_<tier2>_completeness.png``: for each slice, the
host map, the subsample average, and the all-star average around the draws.
<results_dir> defaults to the run script's results/ folder.
"""
import copy
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

ECC_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ECC_REPO))

import run_occurrence as ro  # noqa: E402  (also puts occurrence on the path)
from config_dict import tier2_df_cuts_dict  # noqa: E402
from occurrence import completeness_utils as cu  # noqa: E402
from occurrence import dimensions as dm  # noqa: E402
from occurrence import main  # noqa: E402

tier2, companion = sys.argv[1], sys.argv[2]
results = Path(sys.argv[3]) if len(sys.argv) > 3 else ro.RESULTS_DIR
tier1 = results / "mtrue"
spec = dm.parse_dimensions(ro.ECC_DIMENSIONS)
paths = main.dimension_tier2_paths(str(tier1), tier2, spec)
# The moving-average box of the maps, in dex (moving_average_grid)
BOX_DEX = 0.5

catalog = np.load(paths["catalog"])
if companion not in catalog.files:
    sys.exit(f"{companion} is not in {paths['catalog']}")
a, m, e = (catalog[f"{companion}::{key}"] for key in "ame")
stored = catalog[companion][3]
inside = catalog[f"{companion}::marginal_factor"] > 0

star_df = ro.load_or_make_star_df()
host = star_df[star_df.comp_list.apply(lambda names: companion in names)]
host = host.star_name.iloc[0]
query = tier2_df_cuts_dict[tier2][0].get("star_df_query")
subsample = (star_df.query(query) if query else star_df).star_name.tolist()

maps_dir = Path(main.sliced_maps_dir(str(tier1), "mtrue"))
host_map = cu.load_sliced(maps_dir / host)
average = cu.load_sliced(paths["sliced_average"])
allstars_dir = tier1 / "allstars" / Path(paths["sliced_average"]).name
allstars = cu.load_sliced(allstars_dir) if allstars_dir.exists() else None
stars_per_point = np.load(Path(paths["sliced_average"]) /
                          cu.SLICED_MAP_FILE)["stars_per_point"]
star_maps = [cu.load_sliced(maps_dir / star) for star in subsample
             if (maps_dir / star / cu.SLICED_INTERP_FILE).exists()]
recoveries = pd.read_csv(tier1 / "mtrue_recoveries" /
                         f"{host}_recoveries.csv")
recoveries["slice"] = cu.eccentricity_slice_index(recoveries.inj_e,
                                                  average.e_edges)


def tally(values):
    return (f"{np.sum(values == 0):4d} zero, "
            f"{np.sum(~np.isfinite(values)):4d} missing, "
            f"{np.sum(values > 0):4d} positive")


points = (a, m, e)
sources = [("host map", host_map(points)), ("subsample avg", average(points))]
if allstars is not None:
    sources.append(("all-star avg", allstars(points)))
sources.append(("stored (used)", stored))

print(f"{companion} (host {host}) in {tier2}: {inside.sum()} of {a.size} "
      f"draws inside the marginalized a and m ranges")
if not inside.any():
    sys.exit("None of its draws are inside the marginalized ranges, so the "
             "fit does not use its completeness")
draw_slices = cu.eccentricity_slice_index(e, average.e_edges)
used_slices = sorted(set(draw_slices[inside]))
for index in used_slices:
    selected = inside & (draw_slices == index)
    low, high = average.e_edges[index], average.e_edges[index + 1]
    a_mid, m_mid = np.median(a[selected]), np.median(m[selected])
    print(f"\nSlice e = {low:g}-{high:g}: {selected.sum()} draws, median "
          f"a = {a_mid:.3g} AU, M = {m_mid:.3g} MJ")
    for label, values in sources:
        print(f"  {label:14s} {tally(values[selected])}")
    x_index = np.argmin(np.abs(np.log10(average.xgrid / a_mid)))
    y_index = np.argmin(np.abs(np.log10(average.ygrid / m_mid)))
    point = (np.array([a_mid]), np.array([m_mid]))
    star_values = np.array([star.slice_interpolator(index)(point)[0]
                            for star in star_maps])
    print(f"  At the median: subsample avg = "
          f"{average.slice_interpolator(index)(point)[0]:.3g}; "
          f"{int(stars_per_point[index, y_index, x_index])} of "
          f"{len(subsample)} stars have an estimate there, "
          f"{np.sum(star_values > 0)} of them positive")
    box = recoveries[
        (recoveries.slice == index) &
        (np.abs(np.log10(recoveries.inj_au / a_mid)) <= BOX_DEX/2) &
        (np.abs(np.log10(recoveries.inj_mtrue / m_mid)) <= BOX_DEX/2)
    ]
    orbits = box.drop_duplicates(["inj_au", "inj_e"])
    print(f"  Host injections in the {BOX_DEX}-dex box there: {len(box)} rows "
          f"from {len(orbits)} injected orbits, {int(orbits.recovered.sum())} "
          "recovered")

# Maps around the draws, one column per slice
maps = [("Host " + host, host_map), (f"{tier2} average", average)]
if allstars is not None:
    maps.append(("allstars average", allstars))
a_range = (np.min(a[inside]) / 10, np.max(a[inside]) * 10)
m_range = (np.min(m[inside]) / 10, np.max(m[inside]) * 10)
cmap = copy.copy(plt.get_cmap("viridis"))
cmap.set_bad("lightgray")
fig, axes = plt.subplots(len(maps), len(used_slices),
                         figsize=(6*len(used_slices), 5*len(maps)),
                         squeeze=False, sharex=True, sharey=True,
                         constrained_layout=True)
for column, index in enumerate(used_slices):
    selected = inside & (draw_slices == index)
    for row, (label, sliced) in enumerate(maps):
        ax = axes[row, column]
        mesh = ax.pcolormesh(sliced.xgrid, sliced.ygrid,
                             np.ma.masked_invalid(sliced.zgrids[index]),
                             vmin=0, vmax=1, cmap=cmap, shading="auto")
        ax.scatter(a[selected], m[selected], s=4, color="red",
                   label=f"{companion} draws in this slice")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(*a_range)
        ax.set_ylim(*m_range)
        low, high = average.e_edges[index], average.e_edges[index + 1]
        ax.set_title(f"{label}, e = {low:g}-{high:g}")
        if row == len(maps) - 1:
            ax.set_xlabel("Semi-major axis [AU]")
        if column == 0:
            ax.set_ylabel(r"$M_c$ [$M_{\rm Jup}$]")
axes[0, 0].legend(loc="lower right", fontsize=12)
fig.colorbar(mesh, ax=axes, shrink=0.8).set_label(
    "Completeness (gray: no estimate)"
)
fig.suptitle(f"Completeness at {companion}'s posterior draws ({tier2})",
             fontsize=20)
figure_path = Path(f"{companion}_{tier2}_completeness.png")
fig.savefig(figure_path, dpi=150)
print(f"\nFigure: {figure_path.resolve()}")
