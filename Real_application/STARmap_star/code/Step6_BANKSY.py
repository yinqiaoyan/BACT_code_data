##############################################################################
##                              STARmap*                                    ##
##############################################################################

##############################################################################
## This code file is developed by
##
## Vipul Singhal, Nigel Chou, Joseph Lee, Yifei Yue, Jinyue Liu, Wan Kee Chock, 
##    Li Lin, Yun-Ching Chang, Erica Mei Ling Teo, Jonathan Aow, Hwee Kuan Lee, 
##    Kok Hao Chen and Shyam Prabhakar.
## BANKSY unifies cell typing and tissue domain segmentation for scalable 
##    spatial omics data analysis
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

import warnings
warnings.filterwarnings("ignore") 
import os, time
import pandas as pd
import numpy as np
import random
import scipy.sparse as sparse
from scipy.sparse import csr_matrix, issparse
import scanpy as sc
import sklearn

from banksy.initialize_banksy import initialize_banksy
from banksy.run_banksy import run_banksy_multiparam
from banksy_utils.color_lists import spagcn_color

start = time.perf_counter_ns()
random_seed = 234
cluster_algorithm = 'mclust' # 'mclust' 'leiden'
np.random.seed(random_seed)
random.seed(random_seed)

from banksy_utils.load_data import load_adata

# Main tunable variation in running the BANKSY algorithm

coord_keys = ('X', 'Y', 'spatial')
num_clusters = 16  # cell type number

# Read
adata = sc.read_h5ad("../input_data/STARmap_star.h5ad")
adata.var_names_make_unique()

truth_labels = pd.read_csv("../input_data/cell_annotation_file.txt", sep='\s+')
adata.obs['Annotation'] = truth_labels['Annotation'].astype("category")
adata.obs['Annotation']

annotation_key = 'manual_annotations'
adata.obs[annotation_key] = adata.obs["Annotation"]
adata.obs

# Specifying parameters for BANKSY
resolutions = [.9] # clustering resolution for Leiden clustering
pca_dims = [50] # number of dimensions to keep after PCA
lambda_list = [.8] # lambda
k_geom = 6 # spatial neighbours
max_m = 1 # use AGF
nbr_weight_decay = "scaled_gaussian" # can also be "reciprocal", "uniform" or "ranked"

# Initialize Banksy Object
banksy_dict = initialize_banksy(
    adata,
    coord_keys,
    k_geom,
    nbr_weight_decay=nbr_weight_decay,
    max_m=max_m,
    plt_edge_hist=True,
    plt_nbr_weights=True,
    plt_agf_angles=False,
    plt_theta=False,
)

# Run BANKSY using defined parameters
results_df = run_banksy_multiparam(
    adata,
    banksy_dict,
    lambda_list,
    resolutions,
    color_list = spagcn_color,
    max_m = max_m,
    filepath = output_folder,
    key = coord_keys,
    pca_dims = pca_dims,
    annotation_key = annotation_key,
    max_labels = num_clusters,
    cluster_algorithm = cluster_algorithm,
    match_labels = False,
    savefig = False,
    add_nonspatial = False,
    variance_balance = False,
)

SaveFile = pd.concat([adata.obs["X"], adata.obs["Y"], results_df.loc[results_df.index[0], 'labels']], axis=1)
SaveFile.columns = ['coord_x', 'coord_y', 'mclust']
SaveFile.to_csv("../result_data/starmap_result_BANKSY.csv")

