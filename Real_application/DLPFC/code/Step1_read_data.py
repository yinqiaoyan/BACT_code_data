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
import scanpy as sc

# MERFISH_0.19 is used for visualizing the cell typing performances of all methods.
adata = sc.read_h5ad("../input_data/151507.h5ad")

df = adata.to_df()
df.to_csv("../input_data/gene_matrix_raw_151507.csv", sep=",")
adata.obs[['array_row', 'array_col']].to_csv("../input_data/coordinates_151507.csv", index=False)
adata.obs.to_csv("../input_data/ground_truth_151507.csv", index=False)

