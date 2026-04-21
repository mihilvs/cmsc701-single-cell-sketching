import argparse
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from sklearn.metrics import adjusted_rand_score
from sklearn.neighbors import NearestNeighbors

def evaluate_pipeline(exact_h5ad, sketch_h5, baseline_csv):
    """
    Grades the biological fidelity of a sketched scRNA-seq matrix against exact baselines.
    """
    print("\n--- Initializing Evaluation Harness ---")
    sc.settings.verbosity = 0  # Turn off scanpy logs for a clean output report
    
    # 1. Load the Baseline Ground Truth
    print("Loading exact baseline clusters...")
    baseline_df = pd.read_csv(baseline_csv, index_col=0)
    baseline_clusters = baseline_df['baseline_cluster'].values

    # 2. Load the Sketched Data
    print("Loading sketched out-of-core matrix...")
    with h5py.File(sketch_h5, 'r') as f:
        sketched_matrix = f['sketched_matrix'][:]
        
    n_cells, sketch_dim = sketched_matrix.shape
    
    # 3. Compute Sketched Clustering (for ARI)
    print("\n--- Computing Sketched Clusters ---")
    # Wrap the dense sketched array back into an AnnData object
    adata_sketch = ad.AnnData(X=sketched_matrix)
    
    # Compute neighbors and cluster directly on the sketch (no PCA needed)
    sc.pp.neighbors(adata_sketch, n_neighbors=15, use_rep='X')
    sc.tl.leiden(adata_sketch, resolution=1.0, key_added='sketch_leiden')
    
    sketch_clusters = adata_sketch.obs['sketch_leiden'].values
    
    # Grade the clusters
    ari_score = adjusted_rand_score(baseline_clusters, sketch_clusters)
    
    # 4. Compute Geometric kNN Overlap
    print("\n--- Computing kNN Geometric Overlap ---")
    # We need the exact data to find the true high-dimensional neighbors
    adata_exact = sc.read_h5ad(exact_h5ad)
    
    K = 15
    print(f"Extracting exact top-{K} neighbors (this may take a moment)...")
    # sklearn handles the sparse CSR matrix natively here
    nn_exact = NearestNeighbors(n_neighbors=K, algorithm='auto').fit(adata_exact.X)
    exact_indices = nn_exact.kneighbors(return_distance=False)
    
    print(f"Extracting sketched top-{K} neighbors...")
    nn_sketch = NearestNeighbors(n_neighbors=K, algorithm='auto').fit(sketched_matrix)
    sketch_indices = nn_sketch.kneighbors(return_distance=False)
    
    # Calculate the intersection of neighbors cell by cell
    total_overlap = 0
    for i in range(n_cells):
        intersection = np.intersect1d(exact_indices[i], sketch_indices[i])
        total_overlap += len(intersection)
        
    max_possible_overlap = n_cells * K
    overlap_percentage = (total_overlap / max_possible_overlap) * 100
    
    # 5. The Final Report
    print("\n" + "="*40)
    print("          FINAL GRADE REPORT")
    print("="*40)
    print(f"Cells Processed:   {n_cells}")
    print(f"Sketch Dimension:  {sketch_dim} (Target $k$)")
    print(f"Clustering ARI:    {ari_score:.4f}  (1.0 is a perfect match)")
    print(f"kNN Overlap:       {overlap_percentage:.2f}% (100% is perfect geometry)")
    print("="*40 + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate sketch accuracy against exact baselines.")
    parser.add_argument(
        "--exact", 
        type=str, 
        default="data/pbmc3k_preprocessed.h5ad", 
        help="Path to exact preprocessed data."
    )
    parser.add_argument(
        "--sketch", 
        type=str, 
        default="data/pbmc3k_sketched_k500.h5", 
        help="Path to the sketched matrix."
    )
    parser.add_argument(
        "--baseline", 
        type=str, 
        default="data/pbmc3k_baseline_clusters.csv", 
        help="Path to the exact baseline clusters."
    )
    
    args = parser.parse_args()
    evaluate_pipeline(args.exact, args.sketch, args.baseline)
