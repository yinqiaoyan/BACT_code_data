##############################################################
#################        Slide-seq         ###################
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
geneData_raw = read.csv("../input_data/gene_raw_counts.csv")
dim(geneData_raw)

gene_name = geneData_raw[, 1]
geneData_raw_noName = as.matrix(geneData_raw[, -1])
rownames(geneData_raw_noName) = gene_name
dim(geneData_raw_noName)

colnames(coord) = c("coord_x", "coord_y")


# data_use = t(geneData_raw_noName)
data_use = geneData_raw_noName
dim(data_use)
rownames(coord) = colnames(data_use)
dim(coord)

## Auxiliary function for computing posterior mode
getmode <- function(v) {
  uniqv <- unique(v)
  res <- uniqv[which.max(tabulate(match(v, uniqv)))]
  return(res)
}

set.seed(123)

BASSObject = createBASSObject(X = list(data_use), xy = list(coord), k = 6,
                              C = 4, R = 4, init_method = "mclust")
BASSObject_proc = BASS.preprocess(BASSObject, nHVG = 5000, nPC = 50,
                                  geneSelect = "hvgs")
BASSObject_res = BASS.run(BASSObject_proc)
BASSObject_res_adj = BASS.postprocess(BASSObject_res, adjustLS = TRUE)
post_samples = BASSObject_res_adj@samples

cIds_mode = apply(post_samples$c, 1, getmode)
zIds_mode = apply(post_samples$z, 1, getmode)


## Save results
res2 = data.frame(x=coord[,1], y=coord[,2], c=cIds_mode, z=zIds_mode)
write.csv(x = res2, file = "../result_data/slideseq_result_BASS.csv")



