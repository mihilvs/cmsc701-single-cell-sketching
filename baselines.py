import os
import argparse
import pandas as pd
import scanpy as sc

def generate_baselines(input_path, output_dir):
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

    # 2. Exact PCA (The Dense Matrix Math)
    print("\n--- Step 1: Computing Exact PCA ---")
    # We calculate the first 50 Principal Components. This is the standard 
    # dense calculation that we want to avoid with our sketching engine later.
    sc.tl.pca(adata, svd_solver='arpack', n_comps=50)

    # 3. Exact Nearest Neighbors
    print("\n--- Step 2: Computing Exact k-Nearest Neighbors ---")
    # This builds the exact graph of which cells are most similar to each other.
    # It runs on the 50 dimensions from the PCA we just calculated.
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)

    # 4. Baseline Clustering
    print("\n--- Step 3: Running Leiden Clustering ---")
    # Leiden is the industry standard community detection algorithm for scRNA-seq.
    # We save these cluster labels into a specific column named 'baseline_leiden'.
    sc.tl.leiden(adata, resolution=1.0, key_added='baseline_leiden')

    # 5. Checkpointing the Results
    print("\n--- Step 4: Saving Ground Truth ---")
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
    print(f"Success! Ground truth clusters saved to: {out_file}")

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
    
    args = parser.parse_args()
    generate_baselines(args.input, args.output_dir)
