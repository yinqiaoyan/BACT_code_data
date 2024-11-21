##############################################################
#################       DLPFC 151507       ###################
##############################################################

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

library(SingleCellExperiment)
library(BACT)

## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())


## Read raw data
coord = read.csv("../input_data/coordinates_151507.csv")
geneData_raw = read.csv("../input_data/gene_matrix_raw_151507.csv")

spotLoc_char = geneData_raw[, 1]
geneData_raw_noName = as.matrix(geneData_raw[, -1])
rownames(geneData_raw_noName) = spotLoc_char
# Dim: cells * genes =  4226 * 33538
colnames(coord) = c("coord_x", "coord_y")


## Build SingleCellExperiment object
tmpx = t(geneData_raw_noName)
sce <- SingleCellExperiment(assays=list(counts=tmpx),
                            colData = coord)


## Preprocess the raw data
# 1. For BACT
# Use R function DataPreprocess in the package BACT
# norm.type="logNorm"
sce = DataPreprocess(sce, n.PCs=50, norm.type="logNorm", 
                     select.hvg=TRUE, n.HVGs = 5000)
## Get the processed data after conducting PCA
gene_data_pc = t(reducedDim(sce, "PCA"))


## Save preprocessed data
save(coord, gene_data_pc, file = "../input_data/coord_and_pc_151507.RData")

