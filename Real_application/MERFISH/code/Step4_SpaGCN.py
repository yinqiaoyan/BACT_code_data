##############################################################################
##                            MERFISH_0.19                                  ##
##############################################################################

##############################################################################
## This code file is developed by
##
## Jian Hu, Xiangjie Li, Kyle Coleman, Amelia Schroeder, Nan Ma, David J. Irwin, 
##    Edward B. Lee, Russell T. Shinohara and Mingyao Li.
## SpaGCN: Integrating gene expression, spatial location and histology to
##    identify spatial domains and spatially variable genes by graph 
##    convolutional network.
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

import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import math
from scipy.sparse import issparse
import random, torch
import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import SpaGCN as spg

import scipy
import anndata
#In order to read in image data, we need to install some package. Here we recommend package "opencv"
#inatll opencv in python
#!pip3 install opencv-python
import cv2
import time



coord = pd.read_csv("../input_data/coordinates_0.19.csv")
coord.columns = ['coord_x','coord_y']

# generate adata
geneExpr = pd.read_csv("../input_data/gene_matrix_raw_19.csv")
geneExpr = geneExpr.iloc[:, 1:]

# Set coordinates
x_array=coord["coord_x"]
y_array=coord["coord_y"]

geneExpr_sp = scipy.sparse.csr_matrix(geneExpr.values)
adata = anndata.AnnData(X = geneExpr_sp, obs = coord)

# Normalization
sc.pp.log1p(adata)

# Run SpaGCN
n_clusters = 15
r_seed=t_seed=n_seed=100

adj=spg.calculate_adj_matrix(x=x_array,y=y_array, x_pixel=x_array, y_pixel=y_array, image=0, beta=49, alpha=1, histology=False)
spg.prefilter_genes(adata,min_cells=3) # avoiding all genes are zeros
spg.prefilter_specialgenes(adata)

l=spg.search_l(0.5, adj, start=0.01, end=1000, tol=0.01, max_run=100)
res=spg.search_res(adata, adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)
clf=spg.SpaGCN()
clf.set_l(l)
random.seed(r_seed)
torch.manual_seed(t_seed)
np.random.seed(n_seed)
clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
y_pred, prob=clf.predict()

adata.obs["pred"]= y_pred
adata.obs["pred"]= adata.obs["pred"].astype('category')

# Save results
SaveFile = pd.concat([adata.obs["coord_x"], adata.obs["coord_y"], adata.obs["pred"]], axis=1)
SaveFile.to_csv("../result_data/merfish_0.19_result_SpaGCN.csv", index=False)
