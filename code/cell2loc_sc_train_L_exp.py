import time
import os
import logging
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import torch
import cell2location

# Set up logging
log_file = 'procedure_log.txt'
logging.basicConfig(
    filename=log_file,
    filemode='w',  # Overwrite log file each run
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
console = logging.StreamHandler()
console.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
console.setFormatter(formatter)
logging.getLogger('').addHandler(console)

# Start overall timer
start_time = time.time()
logging.info("Script started.")

logging.info("Set environment.")
results_folder = '.'
ref_run_name = f'{results_folder}/reference_signatures_v3_cutoff2_0.03_L'

if not os.path.exists(results_folder):
    os.mkdir(results_folder)
    logging.info(f"Created directory: {results_folder}")

# Load data
logging.info("Loading data.")
t0 = time.time()
adata_ref = sc.read_h5ad('/fs/cbsuvlaminck2/workdir/jp2626/chicken/chicken_qc_processed_rohitadd.h5ad')
adata_ref = adata_ref[adata_ref.obs['site'] == 'LAL'].copy()
logging.info(f"Data loaded in {time.time() - t0:.2f} seconds.")

# sc Reference Filtering
logging.info("Filtering sc reference.")
adata_ref.X = adata_ref.layers['count']

from cell2location.utils.filtering import filter_genes
t0 = time.time()
selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
adata_ref = adata_ref[:, selected].copy()
logging.info(f"Filtering completed in {time.time() - t0:.2f} seconds.")

# sc Reference Modeling
logging.info("Setting up AnnData for regression model.")
t0 = time.time()
cell2location.models.RegressionModel.setup_anndata(
    adata=adata_ref, batch_key='sample', labels_key='sub_celltype_v3_fin'
)
logging.info(f"AnnData setup completed in {time.time() - t0:.2f} seconds.")

# Create the regression model
logging.info("Creating regression model.")
t0 = time.time()
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)
logging.info(f"Model creation completed in {time.time() - t0:.2f} seconds.")

# Train the model
logging.info("Training model.")
t0 = time.time()
mod.train(max_epochs=500, accelerator="gpu")
logging.info(f"Model training completed in {time.time() - t0:.2f} seconds.")

# Save model
logging.info("Saving model.")
t0 = time.time()
mod.save(f"{ref_run_name}", overwrite=True)
logging.info(f"Model saved in {time.time() - t0:.2f} seconds.")

# Export posterior estimates
logging.info("Exporting posterior estimates.")
t0 = time.time()
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500}
)
logging.info(f"Posterior export completed in {time.time() - t0:.2f} seconds.")

# Save anndata object with results
logging.info("Saving AnnData object.")
t0 = time.time()
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref.write(adata_file)
logging.info(f"AnnData object saved in {time.time() - t0:.2f} seconds.")

# Total script runtime
total_time = time.time() - start_time
logging.info(f"Total script runtime: {total_time:.2f} seconds.")
