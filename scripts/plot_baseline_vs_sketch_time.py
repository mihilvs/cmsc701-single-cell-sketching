#!/usr/bin/env python3
"""
Compare baseline runtime vs sketch-path runtimes (stacked sketch + evaluate).

Example (Allen 100k, SRP):
  python scripts/plot_baseline_vs_sketch_time.py \\
    --baseline_csv runs/allen_100k/results/baseline_times.csv \\
    --runs_csv runs/allen_100k/results/runs.csv \\
    --out runs/allen_100k/results/baseline_vs_sketch_time.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline_csv", type=Path, required=True)
    parser.add_argument("--runs_csv", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--method", default="srp", choices=["srp", "dense_rp"])
    args = parser.parse_args()

    baseline_s = float(pd.read_csv(args.baseline_csv)["baseline_seconds"].iloc[0])

    runs = (
        pd.read_csv(args.runs_csv)
        .query("method == @args.method")
        .drop_duplicates(subset=["k"], keep="first")
        .sort_values("k")
    )
    ks = runs["k"].astype(int).tolist()
    sketch_s = runs["sketch_seconds"].astype(float).to_numpy()
    eval_s = runs["eval_seconds"].astype(float).to_numpy()

    fig, ax = plt.subplots(figsize=(8, 5))

    # Baseline: one bar at x=0
    ax.bar(0, baseline_s, width=0.6, color="tab:blue", label="Baseline (PCA + nn + Leiden)")

    # Per k: stacked sketch + evaluate
    x = np.arange(1, len(ks) + 1)
    ax.bar(x, sketch_s, width=0.6, color="tab:orange", label="Sketching")
    ax.bar(
        x,
        eval_s,
        width=0.6,
        bottom=sketch_s,
        color="tab:gray",
        alpha=0.85,
        label="evaluate.py (sketch cluster + kNN overlap)",
    )

    ax.axhline(baseline_s, color="tab:blue", linestyle="--", linewidth=1, alpha=0.6)

    labels = ["Baseline"] + [str(k) for k in ks]
    ax.set_xticks([0, *x])
    ax.set_xticklabels(labels)
    ax.set_xlabel("Sketch dimension k (evaluate columns)")
    ax.set_ylabel("Wall time (seconds)")
    ax.set_yscale("log")
    ax.set_title(f"Baseline vs sketch-path runtime ({args.method}, log scale)")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(True, axis="y", alpha=0.35, which="both")
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=150)
    plt.close(fig)
    print("Wrote", args.out)


if __name__ == "__main__":
    main()
