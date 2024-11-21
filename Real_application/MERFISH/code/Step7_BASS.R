##############################################################
#################       MERFISH_0.19       ###################
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
section_id = "19"

## Read raw data
coord = read.csv(paste0("../input_data/coordinates_0.", section_id, ".csv"))
geneData_raw = read.csv(paste0("../input_data/gene_matrix_raw_", section_id, ".csv"))

spotLoc_char = geneData_raw[, 1]
geneData_raw_noName = as.matrix(geneData_raw[, -1])
rownames(geneData_raw_noName) = spotLoc_char
dim(geneData_raw_noName)

colnames(coord) = c("coord_x", "coord_y")


data_use = t(geneData_raw_noName)
# The merfish data has been normalized, so only "log1p" is conducted.
data_use = log1p(data_use)
dim(data_use)
rownames(coord) = colnames(data_use)
dim(coord)

## Auxiliary function for computing posterior mode
getmode <- function(v) {
  uniqv <- unique(v)
  res <- uniqv[which.max(tabulate(match(v, uniqv)))]
  return(res)
}

set.seed(120)
BASSObject = createBASSObject(X = list(data_use), xy = list(coord), k = 6,
                              C = 15, R = 8, init_method = "mclust")
# doLogNormalize = F
BASSObject_proc = BASS.preprocess(BASSObject, doLogNormalize = F,
                                  nHVG = 5000, nPC = 50,
                                  geneSelect = "hvgs")
# Remark:
# 1. The number of top PCs is also set to 50 in BACT.
# 2. The step of selecting HVGs is only executed when BASS@P > nHVG. 
# Here, BASS@P < nHVG=5000, so this step was not performed in practice.

BASSObject_res = BASS.run(BASSObject_proc)

BASSObject_res_adj = BASS.postprocess(BASSObject_res, adjustLS = TRUE)
post_samples = BASSObject_res_adj@samples

cIds_mode = apply(post_samples$c, 1, getmode)
zIds_mode = apply(post_samples$z, 1, getmode)


## Save results
res2 = data.frame(x=coord[,1], y=coord[,2], c=cIds_mode, z=zIds_mode)
write.csv(x = res2, file = paste0("../result_data/merfish_0.", section_id, "_result_BASS.csv"))





