# Copyright (c) 2025-2026, The Regents of the University of Michigan
# This file is from the parsnip project, released under the BSD 3-Clause License.

"""Generate a benchmark plot from cProfile .prof files."""

import argparse
import pstats
from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

FONT_PATH = Path.home() / "Library/Fonts/Afacad[wght].ttf"
if FONT_PATH.exists():
    fm.fontManager.addfont(str(FONT_PATH))
    FONT = fm.FontProperties(fname=str(FONT_PATH))
else:
    FONT = fm.FontProperties()

plt.rcParams["svg.fonttype"] = "path"

STRUCTURE_SET = 105
PACKAGES = ["ase", "pymatgen", "pycifrw", "pycifrw-fast", "parsnip", "gemmi"]
PROF_FILES = {pkg: f"{pkg}_{STRUCTURE_SET}.prof" for pkg in PACKAGES}

LABELS = {
    "ase": "ASE",
    "pymatgen": "pymatgen",
    "pycifrw": "PyCifRW",
    "pycifrw-fast": "PyCifRW (flex)",
    "parsnip": "parsnip",
    "gemmi": "gemmi",
}

CUTOFF = 60  # Axis limit, in seconds. Without this we need a log scale.

parser = argparse.ArgumentParser()
parser.add_argument("--out", default=f"benchmark_{STRUCTURE_SET}.svg")
args = parser.parse_args()

# Collect data
data = {}
for pkg in PACKAGES:
    path = Path(PROF_FILES[pkg])
    if not path.exists():
        continue
    stats = pstats.Stats(str(path))
    data[LABELS[pkg]] = stats.total_tt

# Parsnip first, then fastest to slowest
parsnip_time = data.pop("parsnip")
others = sorted(data.items(), key=lambda kv: kv[1])
names = ["parsnip"] + [n for n, _ in others]
data["parsnip"] = parsnip_time
times = [data[n] for n in names]
y = list(range(len(names)))


def _color(name, t):
    """Get the color for each bar.

    Slow bars get light grey. `parsnip` gets our logo color, PyCIFRW gets the purple
    color from the Australian Nuclear Science and Technology Organisation logo, and
    gemmi gets an orange that matches the other two.
    """
    if t > CUTOFF:
        return "#B0B0B0"
    if "parsnip" in name:
        return "#006C60"
    if "gemmi" in name:
        return "#FF9B7D"
    return "#203D79"


# Set up the main plot
sns.set_style("ticks")
fig, ax = plt.subplots(figsize=(7, 3.2))

for i, (n, t) in enumerate(zip(names, times, strict=False)):
    c = _color(n, t)
    bar_end = min(t, CUTOFF)
    if t > CUTOFF:
        # Hatched bar (based on https://nanobind.readthedocs.io/en/latest/benchmark.html)
        ax.barh(
            i,
            bar_end,
            height=0.55,
            color=c,
            edgecolor="white",
            linewidth=0.5,
            hatch="///",
        )
        ax.annotate(
            f"  {t:.0f} s",
            xy=(bar_end, i),
            va="center",
            ha="left",
            fontsize=9,
            fontproperties=FONT,
            color="#555555",
        )
    else:
        ax.barh(i, bar_end, height=0.55, color=c, edgecolor="white", linewidth=0.5)
        txt = f"{t:.1f} s" if t >= 1 else f"{t * 1000:.0f} ms"
        x = t + CUTOFF * 0.01
        ax.text(
            x, i, f"  {txt}", va="center", ha="left", fontsize=9, fontproperties=FONT
        )

# Clean up the figure aesthetics

ax.set_yticks(y)
ax.set_yticklabels(names)
ax.invert_yaxis()
ax.set_xlim(0, CUTOFF)
ax.set_xlabel("Wall-clock time (s)", fontproperties=FONT)
ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
ax.set_title("ltalirz/cif-parsing-benchmark", fontproperties=FONT, fontsize=12)
for label in ax.get_yticklabels() + ax.get_xticklabels():
    label.set_fontproperties(FONT)

sns.despine(left=True, bottom=False)
ax.tick_params(left=False)

fig.tight_layout()
fig.savefig(args.out)
print(f"Saved to {args.out}")
