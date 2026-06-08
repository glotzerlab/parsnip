#!/bin/bash
# Copyright (c) 2025-2026, The Regents of the University of Michigan
# This file is from the parsnip project, released under the BSD 3-Clause License.

set -euo pipefail

REPO_ROOT="$(pwd)"
PLOT_DIR="$REPO_ROOT/doc/generate_benchmark_plots"
BENCH_DIR="$(mktemp -d)"
trap 'rm -rf "$BENCH_DIR"' EXIT

echo "Cloning cif-parsing-benchmark into $BENCH_DIR"
git clone git@github.com:janbridley/cif-parsing-benchmark.git "$BENCH_DIR"
cd "$BENCH_DIR"

uv venv .venv
source .venv/bin/activate

uv sync

tar -xf structures.tar.xz
echo "Running benchmarks..." && bash benchmark.sh

python "$PLOT_DIR/cif_parsing_benchmark_plot.py" --out "$PLOT_DIR/benchmark_105.svg"

echo "... saved to $PLOT_DIR/benchmark_105.svg"
