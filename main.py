import argparse
import os
import subprocess
import sys
from pathlib import Path

from sketch_engine import sketch_h5_filename

# Run from repo root so module paths resolve.
os.chdir(Path(__file__).resolve().parent)


def run_cmd(cmd: list[str]) -> None:
    print(f"\n{'>' * 50}\nEXECUTING: {' '.join(cmd)}\n{'>' * 50}\n")
    subprocess.run(cmd, check=True)


def run_pipeline(
    dataset: str,
    output_dir: str,
    seeds: list[int],
    methods: list[str],
    k_list: list[int],
    chunk_size: int,
) -> None:
    print(f"\n{'=' * 50}\nSTARTING PIPELINE: {dataset} -> {output_dir}\n{'=' * 50}")

    os.makedirs(output_dir, exist_ok=True)
    pre = os.path.join(output_dir, f"{dataset}_preprocessed.h5ad")
    baseline_csv = os.path.join(output_dir, f"{dataset}_baseline_clusters.csv")
    py = sys.executable

    run_cmd([py, "data_prep.py", "--dataset", dataset, "--output_dir", output_dir])
    run_cmd([py, "baselines.py", "--input", pre, "--output_dir", output_dir])

    mode_map = {"srp": "sparse", "dense_rp": "dense"}
    for seed in seeds:
        print(f"\n--- Projection random_state={seed} ---\n")
        for method in methods:
            mode = mode_map[method]
            for k in k_list:
                sketch_name = sketch_h5_filename(pre, k, seed, mode=mode)
                sketch_path = os.path.join(output_dir, sketch_name)
                run_cmd(
                    [
                        py,
                        "sketch_engine.py",
                        "--input",
                        pre,
                        "--output_dir",
                        output_dir,
                        "--sketch_dim",
                        str(k),
                        "--chunk_size",
                        str(chunk_size),
                        "--random_state",
                        str(seed),
                        "--mode",
                        mode,
                    ]
                )
                run_cmd(
                    [
                        py,
                        "evaluate.py",
                        "--exact",
                        pre,
                        "--sketch",
                        sketch_path,
                        "--baseline",
                        baseline_csv,
                    ]
                )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the full sketching pipeline.")
    parser.add_argument(
        "--scale",
        choices=["small", "large", "all"],
        default="small",
        help="small=PBMC3k in data/; large=Allen in data/ (requires raw download).",
    )
    parser.add_argument(
        "--output_dir",
        default=None,
        help="Override output directory (default: data/ or runs/<name>/data).",
    )
    parser.add_argument("--seeds", default="42", help="Comma-separated projection seeds.")
    parser.add_argument(
        "--methods",
        default="srp",
        help="Comma-separated: srp, dense_rp.",
    )
    parser.add_argument("--k_list", default="50,200,800", help="Comma-separated sketch dims.")
    parser.add_argument("--chunk_size", type=int, default=5000)
    args = parser.parse_args()

    seeds = [int(s.strip()) for s in args.seeds.split(",") if s.strip()]
    methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    k_list = [int(k.strip()) for k in args.k_list.split(",") if k.strip()]

    if args.scale in ("small", "all"):
        out = args.output_dir or "data"
        run_pipeline("pbmc3k", out, seeds, methods, k_list, args.chunk_size)

    if args.scale in ("large", "all"):
        out = args.output_dir or "data"
        run_pipeline("allen_brain", out, seeds, methods, k_list, args.chunk_size)
