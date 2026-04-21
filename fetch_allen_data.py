import os
import cellxgene_census

def fetch_allen_brain_slice(output_path, max_cells=250000):
    """
    Downloads a large slice of the Allen Mouse Brain dataset using the CELLxGENE Census.
    """
    print("--- Connecting to CELLxGENE Census ---")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    with cellxgene_census.open_soma() as census:
        print("Querying for Mus musculus (Mouse) brain tissue...")
        
        # We query the census for mouse data, specifically looking for brain tissue.
        # We limit the download to 'max_cells' to keep the file size manageable for a local test.
        adata = cellxgene_census.get_anndata(
            census = census,
            organism = "Mus musculus",
            obs_value_filter = "tissue == 'brain'",
            column_names = {"obs": ["tissue", "cell_type"]},
        )
        
        print(f"Dataset retrieved: {adata.n_obs} cells x {adata.n_vars} genes.")
        
        # Slice it down if it's too large for a local test
        if adata.n_obs > max_cells:
            print(f"Subsampling to {max_cells} cells for local stress test...")
            sc.pp.subsample(adata, n_obs=max_cells, random_state=42)
            
        print(f"Saving raw data to {output_path} (This may take a few minutes)...")
        adata.write(output_path)
        print("Download complete!")

if __name__ == "__main__":
    fetch_allen_brain_slice("raw_data/allen_mouse_brain.h5ad")
