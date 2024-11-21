##############################################################
################## System preparation ########################
##############################################################
rm(list=ls())

### Install basic packages, if necessary
if (!require("aricode", quietly = TRUE)) 
  install.packages("aricode", version="1.0.0")
if (!require("dplyr", quietly = TRUE)) 
  install.packages("dplyr", version="1.0.8")
if (!require("RColorBrewer", quietly = TRUE)) 
  install.packages("RColorBrewer", version="1.1-3")
if (!require("tibble", quietly = TRUE)) 
  install.packages("tibble", version="3.1.6")
if (!require("tidyr", quietly = TRUE)) 
  install.packages("tidyr", version="1.2.0")
if (!require("ggplot2", quietly = TRUE)) 
  install.packages("ggplot2", version="3.3.5")
if (!require("mclust", quietly = TRUE)) 
  install.packages("mclust", version="5.4.9")
if (!require("scales", quietly = TRUE)) 
  install.packages("scales", version="1.1.1")
if (!require("pheatmap", quietly = TRUE)) 
  install.packages("pheatmap", version="1.0.12")
if (!require("latex2exp", quietly = TRUE)) 
  install.packages("latex2exp", version="0.9.4")
if (!require("ggbreak", quietly = TRUE)) 
  install.packages("ggbreak", version="0.1.0")
if (!require("patchwork", quietly = TRUE)) 
  install.packages("patchwork", version="1.1.2")
if (!require("mvtnorm", quietly = TRUE)) 
  install.packages("mvtnorm", version="1.1.3")
if (!require("MCMCpack", quietly = TRUE)) 
  install.packages("MCMCpack", version="1.6.3")
if (!require("coda", quietly = TRUE)) 
  install.packages("coda", version="0.19.4")
if (!require("MASS", quietly = TRUE)) 
  install.packages("MASS", version="7.3.55")
if (!require("Rcpp", quietly = TRUE)) 
  install.packages("Rcpp", version="1.0.10")
if (!require("RcppArmadillo", quietly = TRUE)) 
  install.packages("RcppArmadillo", version="0.12")
if (!require("rBeta2009", quietly = TRUE)) 
  install.packages("rBeta2009", version="1.0")
if (!require("truncnorm", quietly = TRUE)) 
  install.packages("truncnorm", version="1.0.8")
if (!require("edgeR", quietly = TRUE)) 
  install.packages("edgeR", version="3.36.0")

# install the Bioconductor package manager, if necessary
if (!require("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

if (!require("SingleCellExperiment", quietly = TRUE)) 
  BiocManager::install("SingleCellExperiment")
if (!require("scran", quietly = TRUE)) 
  BiocManager::install("scran")
if (!require("scater", quietly = TRUE)) 
  BiocManager::install("scater")
if (!require("BiocSingular", quietly = TRUE)) 
  BiocManager::install("BiocSingular")


### Attention!
# For Windows users, please note that the version of Rtools 
# needs to be compatible with the version of R!


### Install R package "BASS"
devtools::install_github("zhengli09/BASS")


### Install R package "BACT" (our proposed model)
install.packages("PKG_PATH/BACT_1.0.tar.gz", repos = NULL, type="source")
# PKG_PATH: the file path where "BACT_1.0.tar.gz" locates.


