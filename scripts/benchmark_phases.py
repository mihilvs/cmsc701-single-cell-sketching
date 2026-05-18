#!/usr/bin/env python3
"""
Measure wall time and peak RSS for each pipeline phase (preprocess, baseline, sketch, evaluate).

Runs each step as a subprocess and polls `ps -o rss= -p <pid>` (macOS/Linux) for peak resident
set size. RSS units from ps are typically KB on macOS.

Usage (from repo root):
  python scripts/benchmark_phases.py --dataset pbmc3k --output_dir data
  python scripts/benchmark_phases.py --dataset allen_brain --output_dir runs/allen_large/data
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import threading
import time
from dataclasses import dataclass, asdict
from pathlib import Path


@dataclass
class PhaseResult:
    phase: str
    wall_seconds: float
    peak_rss_kb: int
    command: str


def _read_rss_kb(pid: int) -> int:
    try:
        out = subprocess.check_output(
            ["ps", "-o", "rss=", "-p", str(pid)],
            stderr=subprocess.DEVNULL,
            text=True,
        )
        return int(out.strip() or 0)
    except (subprocess.CalledProcessError, ValueError, FileNotFoundError):
        return 0


def run_monitored(cmd: list[str], cwd: Path, poll_interval: float = 0.2) -> tuple[float, int]:
    """Run cmd; return (wall_seconds, peak_rss_kb for main child pid)."""
    t0 = time.perf_counter()
    proc = subprocess.Popen(
        cmd,
        cwd=str(cwd),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    peak = 0

    def poller() -> None:
        nonlocal peak
        while proc.poll() is None:
            rss = _read_rss_kb(proc.pid)
            if rss > peak:
                peak = rss
            time.sleep(poll_interval)
        # one last sample after exit
        rss = _read_rss_kb(proc.pid)
        if rss > peak:
            peak = rss

    th = threading.Thread(target=poller, daemon=True)
    th.start()
    stdout, _ = proc.communicate()
    th.join(timeout=2.0)
    wall = time.perf_counter() - t0
    if proc.returncode != 0:
        sys.stdout.buffer.write(stdout or b"")
        raise subprocess.CalledProcessError(proc.returncode, cmd)
    if stdout:
        sys.stdout.buffer.write(stdout)
    return wall, peak


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", default="pbmc3k", help="e.g. pbmc3k, allen_brain")
    parser.add_argument(
        "--output_dir",
        default="data",
        help="Where preprocessed h5ad and outputs are written (matches data_prep --output_dir)",
    )
    parser.add_argument(
        "--k_list",
        default="50,200,800",
        help="Comma-separated sketch dimensions to benchmark (after preprocess+baseline)",
    )
    parser.add_argument(
        "--seeds",
        default="42",
        help="Comma-separated SparseRandomProjection seeds (non-42 writes *_rs<seed>.h5)",
    )
    parser.add_argument("--python", default=sys.executable, help="Python interpreter to use")
    parser.add_argument(
        "--skip_preprocess",
        action="store_true",
        help="Skip data_prep if preprocessed file already exists",
    )
    parser.add_argument(
        "--skip_baseline",
        action="store_true",
        help="Skip baselines if baseline csv already exists",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=None,
        help="Append one row per phase to this CSV (creates with header if missing)",
    )
    parser.add_argument(
        "--json_out",
        type=Path,
        default=None,
        help="Write full results list as JSON",
    )
    args = parser.parse_args()

    repo = Path(__file__).resolve().parent.parent
    os.chdir(repo)
    cwd = repo

    py = args.python
    dataset = args.dataset
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pre = out_dir / f"{dataset}_preprocessed.h5ad"
    base_csv = out_dir / f"{dataset}_baseline_clusters.csv"

    results: list[PhaseResult] = []

    def run_phase(name: str, cmd: list[str]) -> None:
        print(f"\n=== {name} ===\n{' '.join(cmd)}", flush=True)
        wall, peak_kb = run_monitored(cmd, cwd)
        print(f"  wall: {wall:.3f}s  peak RSS: {peak_kb / 1024:.1f} MB (ps rss field, KB basis)", flush=True)
        results.append(PhaseResult(phase=name, wall_seconds=wall, peak_rss_kb=peak_kb, command=" ".join(cmd)))

    if not args.skip_preprocess:
        run_phase(
            "preprocess",
            [py, "data_prep.py", "--dataset", dataset, "--output_dir", str(out_dir)],
        )
    else:
        if not pre.is_file():
            print("ERROR: --skip_preprocess but missing", pre, file=sys.stderr)
            sys.exit(1)

    if not args.skip_baseline:
        run_phase(
            "baseline",
            [py, "baselines.py", "--input", str(pre), "--output_dir", str(out_dir)],
        )
    else:
        if not base_csv.is_file():
            print("ERROR: --skip_baseline but missing", base_csv, file=sys.stderr)
            sys.exit(1)

    def sketch_h5_basename(k: int, random_state: int) -> str:
        rs = "" if random_state == 42 else f"_rs{random_state}"
        return f"{dataset}_sketched_k{k}{rs}.h5"

    ks = [int(x.strip()) for x in args.k_list.split(",") if x.strip()]
    seeds = [int(x.strip()) for x in args.seeds.split(",") if x.strip()]
    if not seeds:
        print("ERROR: --seeds must list at least one integer.", file=sys.stderr)
        sys.exit(1)

    for seed in seeds:
        for k in ks:
            sketch_h5 = out_dir / sketch_h5_basename(k, seed)
            run_phase(
                f"sketch_k{k}_seed{seed}",
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
                    "3000",
                    "--random_state",
                    str(seed),
                ],
            )
            run_phase(
                f"evaluate_k{k}_seed{seed}",
                [
                    py,
                    "evaluate.py",
                    "--exact",
                    str(pre),
                    "--sketch",
                    str(sketch_h5),
                    "--baseline",
                    str(base_csv),
                ],
            )

    core = [r for r in results if not r.phase.startswith("TOTAL")]
    sketch_eval = [r for r in core if r.phase.startswith(("sketch_", "evaluate_"))]
    sketch_eval_wall = sum(r.wall_seconds for r in sketch_eval)
    sketch_eval_peak = max((r.peak_rss_kb for r in sketch_eval), default=0)
    results.append(
        PhaseResult(
            phase="TOTAL_sketch_and_eval_only",
            wall_seconds=sketch_eval_wall,
            peak_rss_kb=sketch_eval_peak,
            command="(sum wall of sketch_* + evaluate_*; max peak RSS among those phases)",
        )
    )

    total_wall = sum(r.wall_seconds for r in core)
    total_peak = max((r.peak_rss_kb for r in core), default=0)
    results.append(
        PhaseResult(
            phase="TOTAL_all_phases_above",
            wall_seconds=total_wall,
            peak_rss_kb=total_peak,
            command="(sum wall of all non-TOTAL phases; max peak RSS among them)",
        )
    )

    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        with open(args.json_out, "w") as fh:
            json.dump([asdict(r) for r in results], fh, indent=2)
        print("\nWrote", args.json_out)

    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        new_file = not args.csv.is_file()
        with open(args.csv, "a", newline="") as fh:
            w = csv.writer(fh)
            if new_file:
                w.writerow(
                    [
                        "timestamp",
                        "dataset",
                        "output_dir",
                        "phase",
                        "wall_seconds",
                        "peak_rss_mb",
                        "command",
                    ]
                )
            ts = time.strftime("%Y-%m-%dT%H:%M:%S")
            for r in results:
                w.writerow(
                    [
                        ts,
                        dataset,
                        str(out_dir),
                        r.phase,
                        f"{r.wall_seconds:.4f}",
                        f"{r.peak_rss_kb / 1024:.2f}",
                        r.command[:500],
                    ]
                )
        print("Appended", len(results), "rows to", args.csv)


if __name__ == "__main__":
    main()
