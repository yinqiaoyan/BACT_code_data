##############################################################################
##                            DLPFC 151507                                  ##
##############################################################################

##############################################################################
## This code file is developed by
##
## Kangning Dong and Shihua Zhang.
## Deciphering spatial domains from spatially resolved transcriptomics with 
##    an adaptive graph attention auto-encoder.
##############################################################################

########################## Note ##############################
## Please set the working directory to the source file 
## location.
##############################################################

########################## Note ##############################
## Please first install required R and Python packages.
## R: System_preparation.R
## Python: pip install -r requirements.txt
##         pip install -r requirements_spagcn.txt  (PyTorch)
##         pip install -r requirements_stagate.txt (TensorFlow)
##############################################################

import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import os
import sys
import igraph
import scipy
import anndata

import STAGATE

os.environ['R_HOME'] = '/Library/Frameworks/R.framework/Resources'
os.environ['R_USER'] = '/Users/user_name/opt/anaconda3/envs/envir_name/lib/python3.6/site-packages/rpy2'  # Path of the python package "rpy2"

coor_df = pd.read_csv("../input_data/coordinates_151507.csv")
coor_df.columns = ['coord_x','coord_y']

# generate adata
gene_counts = pd.read_csv("../input_data/gene_matrix_raw_151507.csv")
gene_counts = gene_counts.iloc[:, 1:]
gene_counts = gene_counts.astype(int)

adata = sc.AnnData(gene_counts)
adata.var_names_make_unique()
adata.obsm["spatial"] = coor_df.to_numpy()
# sc.pp.calculate_qc_metrics(adata, inplace=True)

# Normalization
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=5000)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.var['highly_variable'] = True

# Constructing the spatial network
STAGATE.Cal_Spatial_Net(adata, k_cutoff=6, model='KNN')
STAGATE.Stats_Spatial_Net(adata)

# Running STAGATE
adata = STAGATE.train_STAGATE(adata, alpha=0, random_seed=204) 

## mclust
adata = STAGATE.mclust_R(adata, used_obsm='STAGATE', num_cluster=7)
obs_df = adata.obs.dropna()

adata.obs['x'] = adata.obsm['spatial'][:,0]
adata.obs['y'] = adata.obsm['spatial'][:,1]
adata.obs["mclust"].index = coor_df["coord_x"].index

SaveFile = pd.concat([coor_df["coord_x"], coor_df["coord_y"], adata.obs["mclust"]], axis=1)
SaveFile.to_csv("../result_data/dlpfc151507_result_STAGATE.csv", index=False)
