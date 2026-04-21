import subprocess
import sys
import argparse

def run_cmd(cmd):
    """Utility to run a shell command and handle errors."""
    print(f"\n{'>'*50}")
    print(f"EXECUTING: {cmd}")
    print(f"{'>'*50}\n")
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"\n[ERROR] Command failed: {cmd}")
        sys.exit(1)

def run_pipeline(dataset):
    print(f"\n{'='*50}")
    print(f"STARTING PIPELINE FOR: {dataset.upper()}")
    print(f"{'='*50}")
    
    data_file = f"data/{dataset}_preprocessed.h5ad"
    baseline_file = f"data/{dataset}_baseline_clusters.csv"
    
    # 1. Data Prep & Baselines
    run_cmd(f"python data_prep.py --dataset {dataset}")
    run_cmd(f"python baselines.py --input {data_file}")
    
    # 2. The Parameter Sweep
    sketch_dims = [50, 200, 800]
    for k in sketch_dims:
        sketch_file = f"data/{dataset}_sketched_k{k}.h5"
        
        # Run the sketching engine (Block size 5,000 handles large data better)
        run_cmd(f"python sketch_engine.py --input {data_file} --sketch_dim {k} --chunk_size 5000")
        
        # Run the grader
        run_cmd(f"python evaluate.py --exact {data_file} --sketch {sketch_file} --baseline {baseline_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the full sketching pipeline.")
    parser.add_argument(
        "--scale", 
        type=str, 
        default="small", 
        choices=["small", "large", "all"],
        help="Which dataset to run: 'small' (PBMC 3k), 'large' (Allen Brain), or 'all'."
    )
    args = parser.parse_args()

    if args.scale in ["small", "all"]:
        run_pipeline("pbmc3k")
        
    if args.scale in ["large", "all"]:
        # Make sure you replaced data_prep.py with the updated version from the previous step!
        run_pipeline("allen_brain")
