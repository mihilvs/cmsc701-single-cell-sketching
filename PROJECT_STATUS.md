# Project Status Report: Single-Cell Sketching Pipeline

## Current Development Status (Phase 1 Complete)
The core infrastructure and small-scale prototype (PBMC 3k) are fully implemented, tested, and passing all baseline comparisons. The out-of-core architecture is stable, and memory constraints have been successfully managed for local execution.

### Completed Milestones:
* **Data Pipeline (`data_prep.py`):** Automated Scanpy preprocessing (filtering, log1p normalization, HVG selection) is fully operational.
* **Out-of-Core Engine (`sketch_engine.py`):** The block-wise chunking logic and Achlioptas Sparse Random Projection are functioning correctly and writing directly to disk without memory spikes.
* **Evaluation Harness (`evaluate.py`):** Automated calculation of kNN geometric overlap and clustering Adjusted Rand Index (ARI) against exact PCA baselines is complete.
* **Automation (`main.py`):** The orchestrator is successfully sweeping through multiple target dimensions (k=50, 200, 800) to generate the empirical data needed for the Pareto frontier analysis.

## Next Steps & Remaining Action Items (Phase 2)
With the systems engineering foundation locked in, the focus now shifts to large-scale empirical validation and the theoretical write-up for the final paper.

1. **Large-Scale Data Validation (Allen Mouse Brain):**
   * Execute `python fetch_allen_data.py` to pull a 250k-cell slice via the CELLxGENE Census API. 
   * Run `python main.py --scale large` to stress-test the chunking engine on the massive dataset. *(Note: The block size in `main.py` has been pre-configured to 5,000 to optimize for this specific run).*
2. **Theoretical Proofs:**
   * Draft the formal mathematical analysis for the paper, specifically proving that the Achlioptas matrix construction satisfies the Johnson-Lindenstrauss distortion bounds for our selected target dimensions.
3. **Final Paper Compilation:**
   * Synthesize the empirical results (Pareto frontier of ARI/kNN overlap vs. computational efficiency) and the theoretical proofs into the final 7-9 page conference-style paper.
