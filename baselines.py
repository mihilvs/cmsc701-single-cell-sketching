import json
import os
import argparse
import time

import pandas as pd
import scanpy as sc


def generate_baselines(input_path, output_dir, timing_json: str | None = None):
    """
    Computes exact PCA, kNN, and Leiden clustering on preprocessed data to serve as a baseline.
    """
    sc.settings.verbosity = 3
    
    # 1. Load the Data
    print(f"\n--- Loading Preprocessed Data ---")
    print(f"Reading from: {input_path}")
    try:
        adata = sc.read_h5ad(input_path)
    except FileNotFoundError:
        print(f"Error: Could not find {input_path}. Did you run data_prep.py first?")
        return

    timings: dict[str, float] = {}

    # 2. Exact PCA (The Dense Matrix Math)
    print("\n--- Step 1: Computing Exact PCA ---")
    t0 = time.perf_counter()
    sc.tl.pca(adata, svd_solver='arpack', n_comps=50)
    timings["pca_seconds"] = time.perf_counter() - t0

    # 3. Exact Nearest Neighbors
    print("\n--- Step 2: Computing Exact k-Nearest Neighbors ---")
    t0 = time.perf_counter()
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
    timings["neighbors_seconds"] = time.perf_counter() - t0

    # 4. Baseline Clustering
    print("\n--- Step 3: Running Leiden Clustering ---")
    t0 = time.perf_counter()
    sc.tl.leiden(
        adata,
        resolution=1.0,
        key_added='baseline_leiden',
        flavor='igraph',
        n_iterations=2,
        directed=False,
    )
    timings["leiden_seconds"] = time.perf_counter() - t0
    timings["baseline_total_seconds"] = (
        timings["pca_seconds"] + timings["neighbors_seconds"] + timings["leiden_seconds"]
    )

    # 5. Checkpointing the Results
    print("\n--- Step 4: Saving baseline cluster labels ---")
    os.makedirs(output_dir, exist_ok=True)
    
    # We don't need to save the whole giant matrix again. 
    # We just need a simple table linking the cell's barcode (ID) to its exact cluster.
    results_df = pd.DataFrame(
        {'baseline_cluster': adata.obs['baseline_leiden']}, 
        index=adata.obs_names
    )
    
    output_filename = os.path.basename(input_path).replace("_preprocessed.h5ad", "_baseline_clusters.csv")
    out_file = os.path.join(output_dir, output_filename)
    
    results_df.to_csv(out_file)
    print(f"Success! Baseline clusters saved to: {out_file}")

    if timing_json:
        os.makedirs(os.path.dirname(timing_json) or ".", exist_ok=True)
        with open(timing_json, "w") as fh:
            json.dump(timings, fh, indent=2)
        print(
            f"Baseline timings (s): PCA {timings['pca_seconds']:.3f}, "
            f"neighbors {timings['neighbors_seconds']:.3f}, "
            f"Leiden {timings['leiden_seconds']:.3f}, "
            f"total {timings['baseline_total_seconds']:.3f}"
        )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate exact baselines for scRNA-seq data.")
    parser.add_argument(
        "--input", 
        type=str, 
        default="data/pbmc3k_preprocessed.h5ad", 
        help="Path to the preprocessed .h5ad file."
    )
    parser.add_argument(
        "--output_dir", 
        type=str, 
        default="data", 
        help="Directory to save the resulting baseline CSV."
    )
    parser.add_argument(
        "--timing_json",
        type=str,
        default=None,
        help="Optional path to write per-step wall times (pca, neighbors, leiden).",
    )

    args = parser.parse_args()
    generate_baselines(args.input, args.output_dir, timing_json=args.timing_json)
