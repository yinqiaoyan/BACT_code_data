##############################################################
###################       STARmap*       #####################
##############################################################

########################## Note ##############################
## Please set the working directory to the source file 
## location.
##############################################################

library(BASS)

## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())


## Read raw data
coord = read.csv("../input_data/coordinates.csv")
geneData_raw = read.csv("../input_data/gene_matrix_raw.csv")

spotLoc_char = geneData_raw[, 1]
geneData_raw_noName = as.matrix(geneData_raw[, -1])
rownames(geneData_raw_noName) = spotLoc_char
dim(geneData_raw_noName)
# Dim: cells * genes =  1207 * 1020
colnames(coord) = c("coord_x", "coord_y")


data_use = t(geneData_raw_noName)
dim(data_use)
rownames(coord) = colnames(data_use)

## Auxiliary function for computing posterior mode
getmode <- function(v) {
  uniqv <- unique(v)
  res <- uniqv[which.max(tabulate(match(v, uniqv)))]
  return(res)
}


## Conduct BASS
set.seed(120)
BASSObject = createBASSObject(X = list(data_use), xy = list(coord), k = 6,
                              psi0=0.5, n0=0.5, alpha0=0.5, step_size=0.15,
                              C = 16, R = 7, init_method = "mclust")

BASSObject_proc = BASS.preprocess(BASSObject, nHVG = 5000, nPC = 50,
                                  geneSelect = "hvgs")
# Remark:
# 1. The number of top PCs is also set to 50 in BACT.
# 2. The step of selecting HVGs is only executed when BASS@P > nHVG. 
# Here, BASS@P < nHVG=5000, so this step was not performed in practice.

BASSObject_res = BASS.run(BASSObject_proc)

BASSObject_res_adj = BASS.postprocess(BASSObject_res, adjustLS = TRUE)
post_samples = BASSObject_res_adj@samples

cIds_mode = apply(post_samples$c, 1, getmode)  # estimated cell clusters
zIds_mode = apply(post_samples$z, 1, getmode)  # estimated region clusters

## Save results
res2 = data.frame(x=coord[,1], y=coord[,2], c=cIds_mode, z=zIds_mode)
write.csv(x = res2, file = "../result_data/starmap_result_BASS.csv")





