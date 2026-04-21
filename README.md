# Sketching for Large-Scale Single-Cell RNA-seq

This repository contains an out-of-core, sublinear sketching pipeline designed to compress high-dimensional single-cell RNA sequencing (scRNA-seq) datasets. It implements Sparse Random Projections using an Achlioptas matrix to bypass the memory constraints of standard dense distance calculations.

## Project Architecture
* **`data_prep.py`**: Downloads and preprocesses raw biological data (filtering, log1p normalization, HVG selection).
* **`fetch_allen_data.py`**: Queries the CELLxGENE Census API to download memory-safe slices of the massive Allen Mouse Brain dataset.
* **`baselines.py`**: Computes Exact PCA, exact k-Nearest Neighbors, and Leiden clustering to establish a ground-truth baseline.
* **`sketch_engine.py`**: The core out-of-core compression loop. Reads sparse data from disk in blocks, multiplies by the projection matrix, and writes the dense sketch back to disk.
* **`evaluate.py`**: Grades the sketched matrix against the baseline using Adjusted Rand Index (ARI) and geometric kNN overlap.
* **`main.py`**: The orchestrator that executes the full pipeline end-to-end.

## Installation
Create a virtual environment and install the required dependencies:
```bash
python3 -m venv sc_env
source sc_env/bin/activate
pip install -r requirements.txt
