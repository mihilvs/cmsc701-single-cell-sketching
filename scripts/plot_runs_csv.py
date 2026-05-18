#!/usr/bin/env python3
"""Plot ARI / kNN overlap / ARI-vs-time from a runs.csv (e.g. allen_100k/results/runs.csv)."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def _k_axis_ticks(
    k_values: list[int] | list[float],
    step: int = 100,
    xtick_start: int | None = None,
) -> list[int]:
    """Regular x-axis ticks; optional xtick_start (e.g. 100 for 100k plots, not 50)."""
    kmin, kmax = int(min(k_values)), int(max(k_values))
    if xtick_start is not None:
        start = xtick_start
    elif kmin <= 50:
        start = 50
    else:
        start = (kmin // step) * step
    ticks = list(range(start, kmax + 1, step))
    if xtick_start is None and kmin not in ticks:
        ticks = sorted({kmin, *ticks})
    if kmax not in ticks:
        ticks.append(kmax)
    return sorted(set(ticks))


def _style_k_axis(
    ax,
    k_values: list[int] | list[float],
    xtick_start: int | None = None,
    xtick_step: int = 100,
) -> None:
    ticks = _k_axis_ticks(k_values, step=xtick_step, xtick_start=xtick_start)
    ax.set_xticks(ticks)
    kmin, kmax = int(min(k_values)), int(max(k_values))
    left = min(kmin - 25, ticks[0] - 25)
    ax.set_xlim(left, ticks[-1] + 25)


def plot_from_csv(
    csv_path: Path,
    out_dir: Path,
    xtick_start: int | None = None,
    xtick_step: int = 100,
) -> None:
    df = pd.read_csv(csv_path)
    df = df.drop_duplicates(subset=["method", "k"], keep="first")
    out_dir.mkdir(parents=True, exist_ok=True)
    k_vals = sorted(df["k"].unique())

    fig, ax = plt.subplots(figsize=(7, 5))
    for method, sty in [("srp", "o-"), ("dense_rp", "s-")]:
        sub = df[df["method"] == method].sort_values("k")
        lab = "SRP" if method == "srp" else "Dense RP"
        ax.plot(sub["k"], sub["ari"], sty, label=lab, linewidth=2, markersize=8)
    ax.set_xlabel("Sketch dimension k")
    ax.set_ylabel("ARI vs. baseline Leiden")
    _style_k_axis(ax, k_vals, xtick_start=xtick_start, xtick_step=xtick_step)
    ax.legend()
    ax.grid(True, alpha=0.35)
    fig.tight_layout()
    fig.savefig(out_dir / "ari_vs_k.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5))
    for method, sty in [("srp", "o-"), ("dense_rp", "s-")]:
        sub = df[df["method"] == method].sort_values("k")
        lab = "SRP" if method == "srp" else "Dense RP"
        ax.plot(sub["k"], sub["knn_overlap_pct"], sty, label=lab, linewidth=2, markersize=8)
    ax.set_xlabel("Sketch dimension k")
    ax.set_ylabel("kNN overlap (%)  (K=15)")
    _style_k_axis(ax, k_vals, xtick_start=xtick_start, xtick_step=xtick_step)
    ax.legend()
    ax.grid(True, alpha=0.35)
    fig.tight_layout()
    fig.savefig(out_dir / "knn_overlap_vs_k.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5))
    for method, c, m in [("srp", "tab:orange", "o"), ("dense_rp", "tab:blue", "s")]:
        sub = df[df["method"] == method].sort_values("k")
        lab = "SRP" if method == "srp" else "Dense RP"
        ax.scatter(
            sub["sketch_seconds"],
            sub["ari"],
            s=100,
            c=c,
            marker=m,
            label=lab,
            zorder=3,
        )
        for _, r in sub.iterrows():
            ax.annotate(
                str(int(r["k"])),
                (float(r["sketch_seconds"]), float(r["ari"])),
                textcoords="offset points",
                xytext=(5, 4),
                fontsize=9,
            )
    ax.set_xlabel("Logged sketch wall time (s)")
    ax.set_ylabel("ARI")
    ax.legend()
    ax.grid(True, alpha=0.35)
    fig.tight_layout()
    fig.savefig(out_dir / "ari_vs_sketch_time.png", dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs_csv", type=Path, required=True)
    parser.add_argument("--out_dir", type=Path, required=True)
    parser.add_argument(
        "--xtick_start",
        type=int,
        default=None,
        help="First x-axis tick for sketch dimension k (e.g. 100 for Allen 100k). Default: 50 if data includes k=50.",
    )
    parser.add_argument("--xtick_step", type=int, default=100)
    args = parser.parse_args()
    plot_from_csv(
        args.runs_csv,
        args.out_dir,
        xtick_start=args.xtick_start,
        xtick_step=args.xtick_step,
    )


if __name__ == "__main__":
    main()
