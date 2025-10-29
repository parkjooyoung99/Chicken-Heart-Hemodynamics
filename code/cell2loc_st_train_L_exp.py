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
import re

# Setup logging
logging.basicConfig(filename='cell2location_log_L.txt', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


start_time = time.time()
logging.info("Script started.")

import glob
adata_files = glob.glob("/fs/cbsuvlaminck5/workdir/jp2626/chickenheart/jy/data/raw_ST_*.h5ad", recursive=True)
adata_files = [   f for f in adata_files if '_L.h5ad' in f and 'L_2' not in f ]

for f in adata_files:
    samplename = os.path.basename(f).replace('raw_ST_', '').replace('.h5ad', '')
    print(samplename)
    print('Set dir')
    logging.info(samplename)
    logging.info("Loading data...")
    results_folder = './' + samplename
    os.makedirs(results_folder, exist_ok=True)
    
    ref_run_name = './reference_signatures_v3_cutoff2_0.03_L'
    run_name = f'{results_folder}/cell2location_map_v3_cutoff2_0.03_L'

    adata_file = f"{ref_run_name}/sc.h5ad"
    adata_ref = sc.read_h5ad(adata_file)
    mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)

    print('Load data')
    adata_st = sc.read_h5ad(f)

    print('Export means_per_cl')
    logging.info("Exporting means per cluster...")
    if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
        inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                        for i in adata_ref.uns['mod']['factor_names']]].copy()
    else:
        inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                        for i in adata_ref.uns['mod']['factor_names']]].copy()
    inf_aver.columns = adata_ref.uns['mod']['factor_names']

    print('Cell2loc')
    logging.info("Preparing Cell2Location model...")
    intersect = np.intersect1d(adata_st.var_names, inf_aver.index)
    adata_st = adata_st[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()

    cell2location.models.Cell2location.setup_anndata(adata=adata_st)

    mod = cell2location.models.Cell2location(
        adata_st, cell_state_df=inf_aver,
        N_cells_per_location=1,
        detection_alpha=20
    )

    logging.info("Training Cell2Location model...")
    train_start_time = time.time()
    mod.train(max_epochs=3000, batch_size=5000, train_size=1, accelerator="gpu")
    logging.info(f"Training completed in {time.time() - train_start_time:.2f} seconds.")

    mod.plot_history(500)
    plt.legend(labels=['full data training'])
    plt.savefig(results_folder + "training_history.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Save data
    logging.info("Exporting posterior estimates...")
    adata_st = mod.export_posterior(
        adata_st, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs}
    )

    print('Save data')
    logging.info("Saving trained model and results...")
    mod.save(f"{run_name}", overwrite=True)
    adata_file = f"{run_name}/sp.h5ad"
    adata_st.write(adata_file)

    end_time = time.time()
    logging.info(f"Total script runtime: {end_time - start_time:.2f} seconds.")
    print(f"Total script runtime: {end_time - start_time:.2f} seconds")

    mod.plot_QC()
    plt.legend(labels=['full data training'])
    plt.savefig(results_folder + "training_QC.png", dpi=300, bbox_inches='tight')
    plt.close()