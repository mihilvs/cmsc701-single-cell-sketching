import os
import argparse
import scanpy as sc

def preprocess_data(dataset_name, output_dir):
    """
    Downloads (or loads) and preprocesses single-cell RNA-seq data.
    """
    sc.settings.verbosity = 3  
    
    # 1. Ingestion
    print(f"\n--- Ingesting {dataset_name} ---")
    if dataset_name == "pbmc3k":
        adata = sc.datasets.pbmc3k()
    elif dataset_name == "pbmc68k_reduced":
        adata = sc.datasets.pbmc68k_reduced()
    elif dataset_name == "allen_brain":
        # Pointing to the raw file you just downloaded via the Census API
        raw_path = "raw_data/allen_mouse_brain.h5ad"
        if not os.path.exists(raw_path):
            raise FileNotFoundError(
                f"\n[ERROR] Missing {raw_path}.\n"
                "Please run fetch_allen_data.py to download the raw data first."
            )
        print(f"Loading large dataset from {raw_path}...")
        adata = sc.read_h5ad(raw_path)
    else:
        raise ValueError(f"Dataset {dataset_name} not currently supported.")
        
    print(f"Original shape: {adata.n_obs} cells x {adata.n_vars} genes")

    # 2. Quality Control & Filtering
    print("\n--- Filtering Dead Cells & Empty Genes ---")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    # 3. Normalization
    print("\n--- Normalizing Counts (Target Sum 10k) & Log1p ---")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # 4. Feature Selection
    print("\n--- Identifying Highly Variable Genes (HVGs) ---")
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    
    adata = adata[:, adata.var.highly_variable]
    print(f"Final Preprocessed shape: {adata.n_obs} cells x {adata.n_vars} genes")

    # 5. Checkpointing
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f"{dataset_name}_preprocessed.h5ad")
    
    print(f"\n--- Saving to Disk ---")
    adata.write(output_path)
    print(f"Success! Data saved to: {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Preprocess scRNA-seq datasets.")
    parser.add_argument(
        "--dataset", 
        type=str, 
        default="pbmc3k", 
        # THE FIX: Added 'allen_brain' to the allowed choices here
        choices=["pbmc3k", "pbmc68k_reduced", "allen_brain"],
        help="The name of the dataset to load and preprocess."
    )
    parser.add_argument(
        "--output_dir", 
        type=str, 
        default="data", 
        help="Directory to save the preprocessed .h5ad file."
    )
    
    args = parser.parse_args()
    preprocess_data(args.dataset, args.output_dir)
