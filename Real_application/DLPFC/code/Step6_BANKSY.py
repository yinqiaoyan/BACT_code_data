##############################################################################
##                            DLPFC 151507                                  ##
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
random_seed = 1231
cluster_algorithm = 'mclust' # 'mclust' 'leiden'
np.random.seed(random_seed)
random.seed(random_seed)

from banksy_utils.load_data import load_adata

# Main tunable variation in running the BANKSY algorithm

save_fig = False
coord_keys = ('array_row', 'array_col', 'spatial')
num_clusters = 7

# Read
adata = sc.read_h5ad("../input_data/151507.h5ad")
adata.var_names_make_unique()

annotation_key = 'manual_annotations'
adata.obs[annotation_key] = adata.obs["ground_truth"]
adata.obs["array_row"] = adata.obsm["spatial"][:,0]
adata.obs["array_col"] = adata.obsm["spatial"][:,1]

num_clusters = 7

resolutions = [.5] # clustering resolution for Leiden clustering
pca_dims = [50] # number of dimensions to keep after PCA
lambda_list = [.8] # lambda
k_geom = 6 # 15 spatial neighbours
max_m = 1 # use AGF
nbr_weight_decay = "scaled_gaussian" # can also be "reciprocal", "uniform" or "ranked"

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

results_df = run_banksy_multiparam(
    adata,
    banksy_dict,
    lambda_list,
    resolutions,  # if cluster_algorithm == 'mclust', then resolutions is NOT used !!!
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


SaveFile = pd.concat([adata.obs["array_row"], adata.obs["array_col"], results_df.loc[results_df.index[0], 'labels'], truth_labels], axis=1)
SaveFile.to_csv("../result_data/dlpfc151507_result_BANKSY.csv")

