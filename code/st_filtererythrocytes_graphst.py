import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
import random
import os
import sys
import time
import csv
import re
import pandas as pd
import scanpy as sc
import numpy as np
import scipy.spatial as scisp
from scipy.sparse import coo_matrix, csr_matrix, csc_matrix
import math
import anndata as ad
import igraph as ig
import plotly.graph_objects as go
import scanpy.external as sce
import scipy.sparse as sp
from statsmodels.nonparametric.smoothers_lowess import lowess
from sklearn.metrics import r2_score
from scipy.interpolate import interp1d
import seaborn as sns
import os
from copy import copy
import matplotlib as mpl
import torch
from sklearn import metrics
import multiprocessing as mp
from GraphST import GraphST


sys.path.append('/workdir/jp2626/graphst/Systematic/ARI/')
from Systematic_fromadata import *

reds = copy(mpl.cm.Reds)
reds.set_under("lightgray")
sc.settings.verbosity = 2 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=100, facecolor="white", frameon=True, figsize=(5, 9))
sc.settings.n_jobs=20
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 500)

import glob

adata_files = glob.glob("/workdir/jp2626/chickenheart/jy/data/smoothie*", recursive=True)
print(adata_files)


for i in adata_files:
    print(f"Processing file: {i}") 
    sixth_element = i.split('/')[6]
    sample = '_'.join(sixth_element.split('_')[1:4]).replace('.h5ad', '')
    start_time = time.time()
    print('Get data')
    meta = pd.read_csv('/fs/cbsuvlaminck3/workdir/jp2626/chickenheart_cell2loc/cell2loc_all_meta.csv', index_col='Unnamed: 0')
    adata = sc.read_h5ad(i)
    meta_filt = meta[meta['sample'].isin(adata.obs['sample'].unique())]
    meta_filt = meta_filt.reindex(adata.obs.index)
    adata.obs['max_pred_celltype'] = meta_filt['max_pred_celltype']
    adata = adata[adata.obs['max_pred_celltype'] != 'Erythrocytes'].copy()
    
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=3000)
    print('Start GraphST ARI')
    file_path = '/workdir/jp2626/chickenheart/jy/graphst_res/noery/' + sample 
    isExist = os.path.exists(file_path)
    print(isExist)
    if not isExist:
        os.mkdir(file_path)
        print('made directory')
        
    print(file_path)
    threshold = 0.65
    device = torch.device('cuda:3' if torch.cuda.is_available() else 'cpu')
    R_HOME = "/home/jp2626/miniconda3/envs/GraphST/lib/R"
    cutoff_point = 'Start'
    tools = 'mclust'
    adata, number_cluster = detect_domain(adata, file_path, cutoff_point ,R_HOME, tools = 'mclust' )