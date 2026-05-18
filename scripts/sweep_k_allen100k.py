#!/usr/bin/env python3
"""
Sweep sketch dimensions for Allen 100k: sketch (SRP and/or dense) + evaluate, append runs.csv.

Example:
  python scripts/sweep_k_allen100k.py --k_list 100,200,300,400,500,600,700,800
  python scripts/sweep_k_allen100k.py --k_list 100,200,300 --methods srp --skip_existing
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from sketch_engine import sketch_h5_filename  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--preprocessed",
        default="runs/allen_100k/data/allen_brain_100k_preprocessed.h5ad",
    )
    parser.add_argument("--output_dir", default="runs/allen_100k/data")
    parser.add_argument(
        "--baseline",
        default="runs/allen_100k/data/allen_brain_100k_baseline_clusters.csv",
    )
    parser.add_argument("--runs_csv", default="runs/allen_100k/results/runs.csv")
    parser.add_argument("--dataset", default="allen_brain_100k")
    parser.add_argument("--k_list", default="100,200,300,400,500,600,700,800")
    parser.add_argument("--methods", default="srp,dense_rp", help="srp, dense_rp, or both")
    parser.add_argument("--chunk_size", type=int, default=5000)
    parser.add_argument("--random_state", type=int, default=42)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--skip_existing", action="store_true")
    parser.add_argument(
        "--eval_only",
        action="store_true",
        help="Skip sketching if .h5 exists; only run evaluate and update CSV.",
    )
    args = parser.parse_args()

    repo = REPO
    py = args.python
    pre = Path(args.preprocessed)
    out_dir = Path(args.output_dir)
    base_csv = Path(args.baseline)
    runs_csv = Path(args.runs_csv)
    ks = [int(x.strip()) for x in args.k_list.split(",") if x.strip()]
    methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    mode_map = {"srp": "sparse", "dense_rp": "dense"}

    if not pre.is_file():
        sys.exit(f"Missing preprocessed file: {pre}")
    if not base_csv.is_file():
        sys.exit(f"Missing baseline CSV: {base_csv}")

    runs_csv.parent.mkdir(parents=True, exist_ok=True)
    existing: set[tuple[str, int]] = set()
    if runs_csv.is_file():
        with open(runs_csv) as fh:
            for row in csv.DictReader(fh):
                existing.add((row["method"], int(row["k"])))

    new_rows: list[dict] = []
    timing_dir = runs_csv.parent / "timing_json"
    timing_dir.mkdir(parents=True, exist_ok=True)

    for method in methods:
        if method not in mode_map:
            sys.exit(f"Unknown method {method!r}")
        mode = mode_map[method]
        for k in ks:
            if args.skip_existing and (method, k) in existing:
                print(f"SKIP (in CSV): {method} k={k}")
                continue

            sketch_name = sketch_h5_filename(str(pre), k, args.random_state, mode=mode)
            sketch_path = out_dir / sketch_name
            metrics_json = timing_dir / f"{method}_k{k}_metrics.json"

            sketch_seconds = 0.0
            if not args.eval_only:
                if args.skip_existing and sketch_path.is_file():
                    print(f"SKIP sketch (file exists): {sketch_path}")
                else:
                    t0 = time.perf_counter()
                    subprocess.run(
                        [
                            py,
                            "sketch_engine.py",
                            "--input",
                            str(pre),
                            "--output_dir",
                            str(out_dir),
                            "--sketch_dim",
                            str(k),
                            "--chunk_size",
                            str(args.chunk_size),
                            "--random_state",
                            str(args.random_state),
                            "--mode",
                            mode,
                        ],
                        cwd=repo,
                        check=True,
                    )
                    sketch_seconds = time.perf_counter() - t0
                    print(f"Sketch {method} k={k}: {sketch_seconds:.1f}s")

            if not sketch_path.is_file():
                sys.exit(f"Missing sketch file: {sketch_path}")

            t0 = time.perf_counter()
            subprocess.run(
                [
                    py,
                    "evaluate.py",
                    "--exact",
                    str(pre),
                    "--sketch",
                    str(sketch_path),
                    "--baseline",
                    str(base_csv),
                    "--timing_json",
                    str(metrics_json),
                ],
                cwd=repo,
                check=True,
            )
            eval_seconds = time.perf_counter() - t0

            with open(metrics_json) as fh:
                metrics = json.load(fh)

            row = {
                "dataset": args.dataset,
                "method": method,
                "k": k,
                "ari": f"{metrics['ari']:.6f}",
                "knn_overlap_pct": f"{metrics['knn_overlap_pct']:.4f}",
                "sketch_seconds": f"{sketch_seconds:.4f}",
                "eval_seconds": f"{eval_seconds:.4f}",
            }
            new_rows.append(row)
            print(f"Done {method} k={k}: ARI={row['ari']} overlap={row['knn_overlap_pct']}%")

    if not new_rows:
        print("No new rows to write.")
        return

    write_header = not runs_csv.is_file()
    with open(runs_csv, "a", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=[
                "dataset",
                "method",
                "k",
                "ari",
                "knn_overlap_pct",
                "sketch_seconds",
                "eval_seconds",
            ],
        )
        if write_header:
            w.writeheader()
        w.writerows(new_rows)
    print(f"Appended {len(new_rows)} rows to {runs_csv}")

    subprocess.run(
        [
            py,
            "scripts/plot_runs_csv.py",
            "--runs_csv",
            str(runs_csv),
            "--out_dir",
            str(runs_csv.parent),
        ],
        cwd=repo,
        check=True,
    )


if __name__ == "__main__":
    main()
