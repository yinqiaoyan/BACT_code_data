##############################################################
#################        Slide-seq         ###################
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
library(scran)  # R==4.1.3
# Bioconductor version 3.14 (BiocManager 1.30.16), R 4.1.3
# remotes::install_version("matrixStats", version="1.1.0")
library(scater)
library(BiocSingular)
library(BACT)

## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())


## Read raw data
coord = read.csv("../input_data/BeadLocationsForR.csv")
geneData_raw = read.csv("../input_data/MappedDGEForR.csv")

gene_names = geneData_raw[, 1]
geneData_raw_noName = as.matrix(geneData_raw[, -1])
rownames(geneData_raw_noName) = gene_names
# Dim: genes * cells =  18906 * 24847
cell_name = coord[, 1]
coord = coord[, -1]
colnames(coord) = c("coord_x", "coord_y")


### Randomly select 8000 cells for downstream analysis
set.seed(15)
ids_subset = sample(1:nrow(coord), 8000, replace = F)

coord_subset = coord[ids_subset, ]
data_subset = geneData_raw_noName[, ids_subset]

## Save subset of raw data
write.csv(coord_subset, file = "../input_data/coordinates.csv", row.names = F)
write.csv(data_subset, file = "../input_data/gene_raw_counts.csv")



## Preprocess the raw data
# 1. For BACT
# Use R function DataPreprocess in the package BACT
# norm.type="logNorm"
sce <- SingleCellExperiment(assays=list(counts=data_subset),
                            colData = coord_subset)
sce = DataPreprocess(sce, n.PCs=50, norm.type="logNorm", 
                     select.hvg=TRUE, n.HVGs = 5000)
## Get the processed data after conducting PCA
gene_data_pc = t(reducedDim(sce, "PCA"))


## Save preprocessed data
coord = coord_subset
save(coord, gene_data_pc, file = "../input_data/coord_and_pc.RData")


# 2. For competing methods
# Take log-normalization, and select 5000 HVGs
sce <- SingleCellExperiment(assays=list(counts=data_subset),
                            colData = coord_subset)
sce <- logNormCounts(sce)
dec <- modelGeneVar(sce, assay.type="logcounts")
top <- getTopHVGs(dec, n=5000)
rowData(sce)[["is.HVG"]] <- (rownames(sce) %in% top)
lognorm_data_hvg = sce@assays@data@listData[["logcounts"]]
lognorm_data_hvg = lognorm_data_hvg[rowData(sce)[["is.HVG"]], ]

write.csv(lognorm_data_hvg, file = "../input_data/processed_gene_data_5000HVGs.csv", row.names = F)







