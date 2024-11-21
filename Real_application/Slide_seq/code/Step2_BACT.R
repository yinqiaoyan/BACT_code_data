##############################################################
#################        Slide-seq         ###################
##############################################################

########################## Note ##############################
## Please set the working directory to the source file 
## location.
##############################################################

library(BACT)

## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())


## Read data
load("../input_data/coord_and_pc.RData")


## Run BACT
# Total execution time is about 1.5 hours
# on a MacBook Pro with Intel Core i5 CPU at 2GHz and 16GB of RAM.
res_list = BACT(gene_data_pc = gene_data_pc, coord = coord, platform = "sc",
                num_init = 4, num_nei = 6,
                a_eta=0, b_eta=1.5, IGkappa=2, IGtau=10, dpAlpha=1,
                a_beta=1, tau_beta=1, tau0=0.01, tau1=0.05, M0=50,
                numOfMCMC=8000, burnIn=4000,
                Is_beta_zero=FALSE, Is_warm_start=TRUE,
                Is_kmeans_use_mean_sd=TRUE,
                Is_print=TRUE, print_gap=500,
                Is_random_seed=TRUE, random_seed=9)


## Auxiliary function for computing posterior mode
getmode <- function(v) {
  uniqv <- unique(v)
  res <- uniqv[which.max(tabulate(match(v, uniqv)))]
  return(res)
}

clIds_mode = apply(res_list$clIds_mcmc, 2, getmode)
save_res = data.frame(coord_x = coord[, 1], coord_y = coord[, 2], labels = clIds_mode)
write.csv(x = save_res, file = "../result_data/slideseq_result_BACT.csv")



