#!/usr/bin/env python
# Copyright (c) 2025-2026, The Regents of the University of Michigan
# This file is from the parsnip project, released under the BSD 3-Clause License.

"""Benchmark parsnip parse modes for build_unit_cell."""

import time
from contextlib import nullcontext
from importlib.util import find_spec
from pathlib import Path
from unittest.mock import patch

import numpy as np

from parsnip import CifFile

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_CIF = REPO_ROOT / "tests/sample_data/AFLOW_mC24.cif"

MODES = [
    # (label, parse_mode, use_cfractions)
    ("python_float", "python_float", True),
    ("cfractions", "rational", True),
    ("default", "rational", False),
]


def _hide_cfractions(name, package=None):
    if name == "cfractions":
        return None
    return find_spec(name, package)


def bench(cif_path: Path, n_runs: int = 200, n_warmup: int = 10):
    """Run a benchmark recording the per-atom performance of unit cell construction."""
    cif = CifFile(cif_path)
    n_atoms = cif.build_unit_cell().shape[0]
    print(f"  {n_atoms} atoms in unit cell")

    results = {}
    for label, parse_mode, use_cfractions in MODES:
        if not use_cfractions:
            ctx = patch("importlib.util.find_spec", side_effect=_hide_cfractions)
        else:
            ctx = nullcontext()

        with ctx:
            for _ in range(n_warmup):
                cif.build_unit_cell(parse_mode=parse_mode)

            times = np.empty(n_runs)
            for i in range(n_runs):
                t0 = time.perf_counter()
                cif.build_unit_cell(parse_mode=parse_mode)
                times[i] = time.perf_counter() - t0

        results[label] = times / n_atoms
        mean_us = times.mean() / n_atoms * 1e6
        std_us = times.std() / n_atoms * 1e6
        print(f"  {label:>15s}: {mean_us:.2f} +/- {std_us:.2f} μs/atom")

    return results


def plot(results, out_path: str):
    """Generate a plot of the parse mode benchmark results."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.font_manager as fm
    import matplotlib.pyplot as plt
    import seaborn as sns

    font_path = Path.home() / "Library/Fonts/Afacad[wght].ttf"
    if font_path.exists():
        fm.fontManager.addfont(str(font_path))
        font = fm.FontProperties(fname=str(font_path))
    else:
        font = fm.FontProperties()

    plt.rcParams["svg.fonttype"] = "path"

    colors = {
        "python_float": "#53FBDD",
        "cfractions": "#00C9A4",
        "default": "#006C60",
    }

    labels = list(results.keys())
    means = [results[label].mean() * 1e6 for label in labels]
    stds = [results[label].std() * 1e6 for label in labels]
    display_names = {
        "python_float": "python_float",
        "cfractions": "rational (cfractions)",
        "default": "rational (stdlib)",
    }

    sns.set_style("ticks")
    fig, ax = plt.subplots(figsize=(7, 2.2))

    y = list(range(len(labels)))
    for i, (lab, m, s) in enumerate(zip(labels, means, stds, strict=False)):
        ax.barh(i, m, height=0.55, color=colors[lab], edgecolor="white", linewidth=0.5)
        ax.text(
            max(means) * 0.02,
            i,
            f"{m:.2f} ± {s:.2f} μs/atom  ",
            va="center",
            ha="left",
            fontsize=9,
            fontproperties=font,
        )

    ax.set_yticks(y)
    ax.set_yticklabels([display_names[lab] for lab in labels])
    ax.invert_yaxis()
    ax.set_xlim(0, max(means) * 1.35)
    ax.set_xlabel("Wall-clock time per atom (μs)", fontproperties=font)
    ax.set_title(
        "build_unit_cell parse mode performance", fontproperties=font, fontsize=12
    )
    for label in ax.get_yticklabels() + ax.get_xticklabels():
        label.set_fontproperties(font)

    sns.despine(left=True, bottom=False)
    ax.tick_params(left=False)

    fig.tight_layout()
    fig.savefig(out_path)
    print(f"Saved to {out_path}")


if __name__ == "__main__":
    out = Path(__file__).parent / "benchmark_parse_modes.svg"
    runs = 5_000

    print(f"Benchmarking {DEFAULT_CIF} ({runs} runs per mode)")
    results = bench(DEFAULT_CIF, n_runs=runs)
    plot(results, str(out))
