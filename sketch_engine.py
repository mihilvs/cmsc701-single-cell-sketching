import os
import argparse
import h5py
import numpy as np
import anndata as ad
import scipy.sparse as sp
from sklearn.random_projection import SparseRandomProjection

def run_sketching(input_path, output_dir, sketch_dim, chunk_size):
    """
    Executes an out-of-core Sparse Random Projection on a high-dimensional scRNA-seq matrix.
    """
    print(f"\n--- Initializing Out-of-Core Sketch Engine ---")
    
    # 1. Backed Data Loading
    # 'backed=r' ensures the massive matrix remains on disk. 
    # We only load metadata into RAM.
    try:
        adata = ad.read_h5ad(input_path, backed='r')
    except FileNotFoundError:
        print(f"Error: Could not find {input_path}.")
        return

    n_cells, n_genes = adata.shape
    print(f"Dataset detected: {n_cells} cells x {n_genes} genes.")
    print(f"Target sketch dimension (k): {sketch_dim}")
    print(f"Processing in blocks of: {chunk_size} cells")

    # 2. Generating the Achlioptas Matrix
    print("\n--- Constructing Sparse Projection Matrix ---")
    # density='auto' defaults to the standard Achlioptas optimal density (1/sqrt(n_features) or 1/3)
    # We fix the random_state so your experiments are reproducible.
    projector = SparseRandomProjection(n_components=sketch_dim, density='auto', random_state=42)
    
    # The projector needs to know the feature dimension (n_genes) to build the matrix.
    # We fit it on a dummy sparse row of zeros to instantiate the matrix R without loading real data.
    dummy_row = sp.csr_matrix((1, n_genes))
    projector.fit(dummy_row)
    print(f"Projection matrix initialized successfully.")

    # 3. Preparing the Output Sink
    os.makedirs(output_dir, exist_ok=True)
    output_filename = os.path.basename(input_path).replace("_preprocessed.h5ad", f"_sketched_k{sketch_dim}.h5")
    output_path = os.path.join(output_dir, output_filename)
    
    print(f"\n--- Starting Block-wise Projection ---")
    # We open an HDF5 file to write the sketched data block-by-block.
    with h5py.File(output_path, 'w') as f:
        # Pre-allocate a dataset on disk of size (N, k)
        out_dataset = f.create_dataset("sketched_matrix", shape=(n_cells, sketch_dim), dtype=np.float32)
        
        # 4. The Out-of-Core Loop
        for i in range(0, n_cells, chunk_size):
            end_idx = min(i + chunk_size, n_cells)
            
            # I/O: Pull a specific block of rows from the hard drive into RAM
            block = adata.X[i:end_idx]
            
            # Math: Multiply the sparse block by the Achlioptas matrix
            sketched_block = projector.transform(block)
            
            # THE FIX: Cast the sparse matrix to a dense NumPy array
            if sp.issparse(sketched_block):
                sketched_block = sketched_block.toarray()
            
            # I/O: Write the dense sketched block directly back to the hard drive
            out_dataset[i:end_idx, :] = sketched_block
            
            print(f"Processed cells {i} to {end_idx} ({(end_idx/n_cells)*100:.1f}%)")

    print(f"\nSuccess! Fully compressed matrix saved to: {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run out-of-core sparse random projection.")
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
        help="Directory to save the sketched .h5 file."
    )
    parser.add_argument(
        "--sketch_dim", 
        type=int, 
        default=500, 
        help="The target number of dimensions (k) for the sketch."
    )
    parser.add_argument(
        "--chunk_size", 
        type=int, 
        default=1000, 
        help="Number of cells to load into memory per block."
    )
    
    args = parser.parse_args()
    run_sketching(args.input, args.output_dir, args.sketch_dim, args.chunk_size)
