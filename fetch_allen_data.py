import os

import cellxgene_census
import scanpy as sc


def fetch_allen_brain_slice(output_path, max_cells=250000):
    """
    Download a subsampled Allen mouse brain slice via the CELLxGENE Census.
    """
    print("--- Connecting to CELLxGENE Census ---")
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

    with cellxgene_census.open_soma() as census:
        print("Querying for Mus musculus (Mouse) brain tissue...")
        adata = cellxgene_census.get_anndata(
            census=census,
            organism="Mus musculus",
            obs_value_filter="tissue == 'brain'",
            column_names={"obs": ["tissue", "cell_type"]},
        )

        print(f"Dataset retrieved: {adata.n_obs} cells x {adata.n_vars} genes.")

        if adata.n_obs > max_cells:
            print(f"Subsampling to {max_cells} cells for local stress test...")
            sc.pp.subsample(adata, n_obs=max_cells, random_state=42)

        print(f"Saving raw data to {output_path} (this may take several minutes)...")
        adata.write(output_path)
        print("Download complete!")


if __name__ == "__main__":
    fetch_allen_brain_slice("raw_data/allen_mouse_brain.h5ad")
