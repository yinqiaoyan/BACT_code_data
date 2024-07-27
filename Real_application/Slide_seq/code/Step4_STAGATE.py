##############################################################################
##                              Slide_seq                                   ##
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

coor_df = pd.read_csv("../input_data/coordinates.csv")
coor_df.columns = ['x', 'y']

geneExpr = pd.read_csv("../input_data/processed_gene_data_5000HVGs.csv")

adata = sc.AnnData(geneExpr.transpose())
adata.var_names_make_unique()
adata.obsm["spatial"] = coor_df.to_numpy()
# sc.pp.calculate_qc_metrics(adata, inplace=True)

# Normalization
# sc.pp.log1p(adata)
adata.var['highly_variable'] = True

# Constructing the spatial network
STAGATE.Cal_Spatial_Net(adata, k_cutoff=6, model='KNN')
STAGATE.Stats_Spatial_Net(adata)

# Running STAGATE
adata = STAGATE.train_STAGATE(adata, alpha=0, n_epochs=2000, random_seed=202)

sc.pp.neighbors(adata, use_rep='STAGATE')
sc.tl.umap(adata)

## mclust
adata = STAGATE.mclust_R(adata, used_obsm='STAGATE', num_cluster=4)
obs_df = adata.obs.dropna()

adata.obs['x'] = adata.obsm['spatial'][:,0]
adata.obs['y'] = adata.obsm['spatial'][:,1]

# Save results
coord_x = pd.Series(data=np.array(coor_df["x"]), index=adata.obs["mclust"].index, name="coord_x", dtype='category')
coord_y = pd.Series(data=np.array(coor_df["y"]), index=adata.obs["mclust"].index, name="coord_y", dtype='category')
SaveFile = pd.concat([coord_x, coord_y, adata.obs["mclust"]], axis=1)
SaveFile.to_csv("../result_data/slideseq_result_STAGATE.csv", index=False)
