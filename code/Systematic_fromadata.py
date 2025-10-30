import os
import torch
import pandas as pd
import scanpy as sc
import numpy as np
from sklearn import metrics
import multiprocessing as mp
from GraphST import GraphST
from GraphST.utils import clustering
from sklearn.metrics import adjusted_rand_score
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import random
import math


# Define the function
def detect_domain(adata, file_path, cutoff_point, R_HOME, tools, threshold=0.65, box_size=4,  datatype='Slide', device=None, number_cluster=None, seed = None):
    """
    Function to train the GraphST model on the input data, perform clustering, and save the results.

    Parameters:
    - adata: AnnData object
    - file_path: Directory to save the files (CSV and h5ad)
    - threshold: Threshold value for the clustering (default 0.65, 0.65 is the minimum of moderate recovery)
    - box_size: Size of the box for clustering (default 4)
    - tools: Clustering method (default mclust.  'mclust', 'leiden', 'louvain')
    - datatype: Type of data used ('Slide', default)
    - device: Device to use for computation ('cuda:3' or 'cpu')
    - R_HOME: Path to R environment (default will be used if None)
    - number_cluster: Number of clusters to generate. If None, it will be set based on the overall ARI score systematically (default will be used if None)
    - seed: seed (int)
    """

    if seed is None:
        random.seed(1999) 
    else:
        random.seed(seed)

        
    if device is None:
        device = torch.device('cuda:3' if torch.cuda.is_available() else 'cpu')
    
    # Set R_HOME if provided or use the default - Not working in andreas's
    if R_HOME is not None:
        os.environ['R_HOME'] = R_HOME
    else:
        os.environ['R_HOME'] = os.environ['R_HOME']

    print(f"Using device: {device}")

    if os.path.exists(f'{file_path}/GraphST.h5ad'):
        print('Reading the GraphST result')
        adata = sc.read_h5ad(f'{file_path}/GraphST.h5ad')
    else:
    # Train the model on the provided data
        print('Start modeling')
        model = GraphST.GraphST(adata, datatype=datatype, device=device)
        adata = model.train()

        print('Saving the GraphST result')
        # Save the results to a file
        adata.write_h5ad(f'{file_path}/GraphST.h5ad')

        # Read the result file
        print('Reading the GraphST result')
        adata = sc.read_h5ad(f'{file_path}/GraphST.h5ad')

    # Detect clusters using ARI cluster detection
    if number_cluster is None:
        print('Finding Optimal number of clusters with {tools}' )  # Only printed when number_cluster is None
        if tools == 'mclust':
            for i in range(1, 21):  # Looping from 1 to 20 inclusive
                print(i)
                clustering(adata, n_clusters=i, method=tools)
                adata.obs[f'optimalnum_{i}'] = adata.obs['domain']
        elif tools in ['leiden', 'louvain']:
            for i in range(1, 21):
                print(i)
                clustering(adata, n_clusters=i, method=tools, start=0.1, end=2.0, increment=0.01)
                adata.obs[f'optimalnum_{i}'] = adata.obs['domain']

        # Save clustering results to CSV
        adata.obs.filter(like="optimalnum_").to_csv(f'{file_path}/ARI_{tools}_clust.csv')

        number_cluster = ari_cluster_detection(file_path, cutoff_point,tools)
        print(f"Detected cluster number: {number_cluster}")
    else:
        print(f"Use cluster number: {number_cluster}")

    print(f'Get Domain of {number_cluster} clusters with {tools}')
    adata.obs['domain'] = adata.obs[f'optimalnum_{number_cluster}']
    
    # Remove previous test result
    adata.obs = adata.obs.drop(columns=adata.obs.filter(like="optimalnum_").columns, errors='ignore')    

    # Save final clustered results
    adata.write_h5ad(f'{file_path}/GraphST_ARI_{tools}_cluster{number_cluster}.h5ad')

    return adata, number_cluster


def ari_cluster_detection(file_path, cutoff_point, tools, threshold=0.65, box_size=4):
    import numpy as np
    
    # Read data
    ari = pd.read_csv(f'{file_path}/ARI_{tools}_clust.csv', index_col=0)  # Ensure index column is correctly read
    
    # Initialize combination matrix
    combinations = [(col1, col2) for idx, col1 in enumerate(ari.columns) for col2 in ari.columns[idx+1:]]
    
    rand_score = {}
    
    # Calculate adjusted rand score for all combinations
    for comp1, comp2 in combinations:
        score = adjusted_rand_score(ari[comp1], ari[comp2])
        rand_score[f"{comp1}-{comp2}"] = score
    
    # Create an empty matrix for the rand scores
    rand_matrix = pd.DataFrame(np.zeros((len(ari.columns), len(ari.columns))), index=ari.columns, columns=ari.columns)
    
    # Fill the rand score matrix
    for comp1, comp2 in combinations:
        score = rand_score[f"{comp1}-{comp2}"]
        rand_matrix.loc[comp1, comp2] = rand_matrix.loc[comp2, comp1] = score
    
    # Set diagonal to 1
    np.fill_diagonal(rand_matrix.values, 1)
    
    # Create heatmap
    col_fun = LinearSegmentedColormap.from_list("blue_white_red", [(0.0, "blue"), (0.65, "white"), (1.0, "red")])
    plt.figure(figsize=(10, 8))
    sns.heatmap(rand_matrix, annot=False, cmap=col_fun, cbar_kws={'label': 'Adjusted Rand Index'}, linewidths=0)
    plt.title('ARI Heatmap')
    plt.savefig(f'{file_path}/ARI_{tools}_heatmap.png' )
    plt.show()
    plt.close()

    # Start cluster detection
    boxes = []
    number_cluster = 0
    
    for start_index in rand_matrix.index:
        values = rand_matrix.loc[start_index, :]
        one_index = np.where(values == 1)[0][0]  # Get the index of the first 1
        values = values.iloc[one_index:]   # Get after values
        
        valid_names = values.index
        
        if len(valid_names) == 0:
            continue  # Skip if no valid block
        
        end_index = start_index
        for name in valid_names:
            if values[name] < threshold:
                break
            end_index = name
        
        box_members = valid_names[valid_names.get_loc(start_index):valid_names.get_loc(end_index) + 1]
       
        if len(box_members) >= box_size and number_cluster == 0:
            cluster_name = start_index if cutoff_point == "Start" else end_index
            number_cluster = int(cluster_name.split('_')[1]) # Extract number from domain
        
        #if len(box_members) >= box_size and number_cluster == 0:
         #   number_cluster = int(end_index.split('_')[1])  # Extract number from domain
        
        boxes.append({
            "Box_Name": f"box_{len(boxes) + 1}st",
            "Start": start_index,
            "End": end_index,
            "Size": len(box_members)
        })
    
    # Save the boxes data to CSV
    boxes_df = pd.DataFrame(boxes)
    boxes_df.to_csv(f'{file_path}/ARI_{tools}_cutoff_df.csv', index=False)
    
    return number_cluster


def subsetdomain_umifilt(adata, domain_value, GENE_COUNT_THRESHOLD = 50):
    """
    Subsets the AnnData object for a specific domain and filters genes by raw UMI counts.

    Parameters:
    adata (AnnData): The input AnnData object.
    domain_value (int): The domain value to filter the data.
    GENE_COUNT_THRESHOLD (int, default 50): The minimum total raw UMI count required for a gene to be retained.

    Returns:
    AnnData: A filtered AnnData object.
    """
    # Subset based on the specified domain
    adataa = adata[adata.obs['domain'] == domain_value].copy()
    
    # Replace X with raw counts
    adataa.X = adataa.layers['raw_X']
    
    # Compute total raw counts per gene
    adataa.var['total_raw_counts'] = np.array(np.sum(adataa.X, axis=0)).flatten()
    
    # Filter genes based on UMI count threshold
    adataa = adataa[:, adataa.var['total_raw_counts'] >= GENE_COUNT_THRESHOLD]
    
    # Replace X with unsmoothed counts
    adataa.X = adataa.layers['unsmoothed_X']
    
    return adataa

def filter_and_compute_expression(node_label_df, sm_adata, output, smdataoutput, min_genes=3):
    """
    Filters clusters with more than `min_genes` genes and computes the mean expression 
    across selected genes for each cluster in the AnnData object.
    
    Parameters:
    - node_label_df (pd.DataFrame): DataFrame with 'community_label' column and gene annotations.
    - sm_adata (AnnData): AnnData object containing gene expression data.
    - min_genes (int): Minimum number of genes for a cluster to be kept.
    - output_csv (str): Path to save the filtered DataFrame.

    Returns:
    - sm_adata (AnnData): Updated AnnData object with mean expression per cluster.
    - node_label_df_filt (pd.DataFrame): Filtered DataFrame with clusters having > min_genes genes.
    """
    # Step 1: Filter clusters with more than `min_genes` genes
    cluster_counts = node_label_df['community_label'].value_counts()
    clusters_to_keep = cluster_counts[cluster_counts > min_genes].index
    node_label_df_filt = node_label_df[node_label_df['community_label'].isin(clusters_to_keep)]

    # Save filtered DataFrame to a CSV file
    node_label_df_filt.to_csv(output, index=False)

    # Step 2: Compute mean expression for each cluster
    for i in clusters_to_keep:
        # Filter genes based on the community label
        genes = node_label_df_filt.loc[node_label_df_filt['community_label'] == i, 'name']

        # Subset the AnnData object to include only the selected genes
        tmp = sm_adata[:, sm_adata.var.index.isin(genes)]

        # Compute the mean expression across the selected genes for each observation
        sm_adata.obs[f'module_{i}'] = np.mean(tmp.X, axis=1)
        sm_adata.write_h5ad(smdataoutput)

    return sm_adata, node_label_df_filt


def plot_modules(adata, file_path, ncol=4):
    """
    Function to plot modules on spatial plots.

    Parameters:
    - adata: AnnData object containing the spatial data
    - ncol: Number of columns for the subplot grid (default is 4)
    - file_path: Path to save the resulting plot as a PDF (default is 'figures/module_plot.pdf')
    """
    # Sort module names
    num_modules = sorted(adata.obs.filter(like="module_").columns, key=lambda x: int(x.split('_')[1]))

    # Define the color map
    reds = "Reds"

    # Set up number of rows and columns for subplots
    nrows = math.ceil(len(num_modules) / ncol)  # Adjusted to use ncol

    figsize = 5
    wspace = 0

    # Set up the figure and axes with custom figsize and spacing
    fig, axs = plt.subplots(
        nrows=nrows,
        ncols=ncol,
        figsize=(ncol * figsize + figsize * wspace * (ncol - 1), nrows * figsize),
    )

    # Adjust the spacing between subplots
    plt.subplots_adjust(wspace=wspace)

    # List of modules to plot
    modules = num_modules

    # Plot each module on its own subplot
    for i, module in enumerate(modules):
        sc.pl.spatial(adata, color=module, spot_size=30, cmap=reds, ax=axs[i // ncol, i % ncol], show=False)

    # Hide any unused subplots
    for j in range(i + 1, len(axs.flat)):  # Corrected to use axs.flat (flattened axis)
        fig.delaxes(axs.flat[j])

    plt.tight_layout()

    # Save the plot to a PDF file with the desired width and height
    plt.savefig(f'{file_path}/module_plot.png' , format='pdf', dpi=300)

    # If you want to show the plot too:
    # plt.show()

    print(f"Modules plot saved to {file_path}/module_plot.png")
