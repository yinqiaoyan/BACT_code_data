##############################################################################
##                              STARmap*                                    ##
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



coord = pd.read_csv("../input_data/coordinates.csv")
coord.columns = ['coord_x','coord_y']

# generate adata
gene_counts = pd.read_csv("../input_data/gene_matrix_raw.csv")
gene_counts = gene_counts.iloc[:, 1:]
gene_counts = gene_counts.astype(int)
gene_counts_sp = scipy.sparse.csr_matrix(gene_counts)
adata = anndata.AnnData(X = gene_counts_sp, obs = coord)

# Set coordinates
x_array=adata.obs["coord_x"].tolist()
y_array=adata.obs["coord_y"].tolist()

# Run SpaGCN
n_clusters = 16  # cell class number
r_seed=t_seed=n_seed=10

adata.obs["pred"]= spg.detect_spatial_domains_ez_mode(adata, 0, x_array, y_array, x_array, y_array, 
                                                      n_clusters=n_clusters, histology=False, s=1, b=49, p=0.5, 
                                                      r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)
adata.obs["pred"]=adata.obs["pred"].astype('category')

# Save results
SaveFile = pd.concat([adata.obs["coord_x"], adata.obs["coord_y"], adata.obs["pred"]], axis=1)
SaveFile.to_csv("../result_data/starmap_result_SpaGCN.csv", index=False)
