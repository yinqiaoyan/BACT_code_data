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
adata = sc.read_h5ad("../input_data/MERFISH_0.19.h5ad")
adata.obs['x'] = adata.obsm['spatial'][:,0]
adata.obs['y'] = adata.obsm['spatial'][:,1]

df = adata.to_df()
df.to_csv("../input_data/gene_matrix_raw_19.csv", sep=",")
pd.core.frame.DataFrame(adata.obsm['spatial']).to_csv("../input_data/coordinates_0.19.csv", index=False)
adata.obs.to_csv("../input_data/adata_obs.csv", index=True)

