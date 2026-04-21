# Sketching for Large-Scale Single-Cell RNA-seq

This repository contains an out-of-core, sublinear sketching pipeline designed to compress high-dimensional single-cell RNA sequencing (scRNA-seq) datasets. It implements Sparse Random Projections using an Achlioptas-style random matrix to bypass the $\mathcal{O}(N^2)$ memory and time constraints of standard dense distance calculations (e.g., exact PCA and exact kNN).

By utilizing a block-wise processing architecture, this system is capable of scaling to massive datasets (e.g., the $10^6+$ cell Allen Mouse Brain Atlas) on standard hardware without exceeding available RAM.

## Project Architecture

- **`data_prep.py`**: Downloads and preprocesses raw biological data (low-count filtering, log1p normalization, Highly Variable Gene selection).
- **`fetch_allen_data.py`**: Queries the CELLxGENE Census API to download memory-safe, subsampled slices of the massive Allen Mouse Brain dataset.
- **`baselines.py`**: Computes Exact PCA (50 components), exact k-Nearest Neighbors (k=15), and Leiden clustering to establish a mathematically exact ground-truth baseline.
- **`sketch_engine.py`**: The core out-of-core compression loop. Reads sparse data from disk in blocks, multiplies the chunks by the sparse Achlioptas projection matrix, and streams the dense sketches directly back to disk.
- **`evaluate.py`**: Grades the sketched matrix against the exact baseline. Outputs the Adjusted Rand Index (ARI) for clustering fidelity and geometric kNN overlap for structural preservation.
- **`main.py`**: The master orchestrator that executes the full ETL and compression pipeline end-to-end.

## Installation

This pipeline requires Python 3.9+ and relies heavily on `scanpy` and `scikit-learn`.

1. Clone the repository:
   ```bash
   git clone [https://github.com/YOUR_USERNAME/cmsc701-single-cell-sketching.git](https://github.com/YOUR_USERNAME/cmsc701-single-cell-sketching.git)
   cd cmsc701-single-cell-sketching
   ```

2. Create an isolated virtual environment and activate it:
   ```bash
   python3 -m venv sc_env
   source sc_env/bin/activate  # On Windows use: sc_env\Scripts\activate
   ```

3. Install the required dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Pipeline Execution

The pipeline is designed to be executed via the `main.py` orchestrator, which sweeps through multiple target sketch dimensions (k = 50, 200, 800) to map the Pareto frontier of computational efficiency versus biological accuracy.

### 1. Small-Scale Validation (10x PBMC 3k)
For rapid prototyping and mathematical validation, run the pipeline on the built-in 10x Genomics PBMC 3k dataset. This dataset is small enough to process in seconds.

```bash
python main.py --scale small
```

**What this does under the hood:**
1. Downloads the raw PBMC matrix via `scanpy`.
2. Preprocesses the data and saves a checkpoint to `data/pbmc3k_preprocessed.h5ad`.
3. Runs the exact O(N^2) baseline calculations and saves the ground-truth clusters to `data/pbmc3k_baseline_clusters.csv`.
4. Executes the out-of-core sketch engine iteratively for k=50, 200, and 800, saving the compressed tensors to the `data/` directory.
5. Evaluates each sketch and prints the final ARI and kNN Overlap metrics to the terminal.

### 2. Large-Scale Out-of-Core Stress Test (Allen Mouse Brain)
To validate the memory limits of the block-wise architecture, run the pipeline on a massive dataset. 

First, fetch a raw, memory-safe slice (250,000 cells) of the Allen Mouse Brain Atlas using the CELLxGENE API:

```bash
python fetch_allen_data.py
```
*(Note: This downloads several gigabytes of raw data into the `raw_data/` directory and may take 5-10 minutes depending on network bandwidth).*

Once the raw data is downloaded, execute the large-scale pipeline:

```bash
python main.py --scale large
```

During this run, `sketch_engine.py` will actively constrain memory by loading the sparse matrix from the hard drive in isolated blocks (default: 5,000 cells per block), compressing them, and writing them to the output `.h5` file sequentially.


### 3. Running Modules Individually
If you need to tune specific hyperparameters (e.g., changing the chunk size or the target variance), you can bypass `main.py` and invoke the modules individually:

**Preprocess custom data:**
```bash
python data_prep.py --dataset pbmc3k --output_dir custom_data/
```

**Generate a specific sketch size with a custom block size:**
```bash
python sketch_engine.py --input data/pbmc3k_preprocessed.h5ad --sketch_dim 500 --chunk_size 2000
```

**Evaluate a specific output:**
```bash
python evaluate.py --exact data/pbmc3k_preprocessed.h5ad --sketch data/pbmc3k_sketched_k500.h5 --baseline data/pbmc3k_baseline_clusters.csv
```

## Output & Data Structures

To ensure Git operations remain fast, the `.gitignore` explicitly ignores the `/data` and `/raw_data` directories. Upon a successful end-to-end run, your local directory structure will look like this:

```text
├── raw_data/
│   └── allen_mouse_brain.h5ad         # Raw API download
├── data/
│   ├── pbmc3k_preprocessed.h5ad       # Cleaned CSR sparse matrix
│   ├── pbmc3k_baseline_clusters.csv   # Ground truth Leiden labels
│   ├── pbmc3k_sketched_k50.h5         # Dense compressed sketch (k=50)
│   ├── pbmc3k_sketched_k200.h5        # Dense compressed sketch (k=200)
│   └── pbmc3k_sketched_k800.h5        # Dense compressed sketch (k=800)
```

## Evaluation Metrics

The pipeline relies on two primary metrics to judge the efficacy of the Achlioptas matrix projection:

* **Adjusted Rand Index (ARI):** Evaluates how closely the Leiden community detection on the compressed sketch matches the exact high-dimensional graph. Scales from -1.0 (independent) to 1.0 (identical).
* **Geometric kNN Overlap:** A strict structural metric measuring the percentage of true k-nearest neighbors that remain neighbors in the heavily compressed latent space, verifying the Johnson-Lindenstrauss distance preservation bounds.

---
**Authors:** Mihil Sreenilayam & Shreeya Venkatesh
