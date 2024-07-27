##############################################################
#################        Slide-seq         ###################
##############################################################

########################## Note ##############################
## Please set the working directory to the source file 
## location.
##############################################################

### edgeR DE analysis ###
library(edgeR)

## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())

## Read data
data_subset = read.csv("../input_data/gene_raw_counts_2.csv")
cell_names = data_subset[, 1]
data_subset = data_subset[, -1]
rownames(data_subset) = cell_names
data_subset = as.matrix(data_subset)  # convert data.frame to matrix

tmpres = read.csv("../result_data/slideseq_result_BACT.csv")
clIds_res = tmpres$labels
group <- rep("control", 8000)
group[clIds_res == 4] = "treat"


### --- data preprocessing --- ###
# Build DGEList object
dgelist <- DGEList(counts = data_subset, group = group)
# dim: 18906 * 8000

# Approach to filter low count data, such as CPM normalization (recommended)
keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
# dim: 13960 * 8000

# Standardization, for example, TMM standardization
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
# dim: 13960 * 8000

### --- DE gene analysis --- ###
# Construct the experimental design matrix based on the grouping information. 
# Notice that the control group comes before the treatment group in the grouping information.
design <- model.matrix(~group)

# Estimate the dispersion of gene expression values
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)

# Model fitting: 
# negative binomial generalized log-linear model
fit <- glmFit(dge, design, robust = TRUE)
lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))


### --- select DE genes --- ###
#读取上述输出的差异倍数计算结果
# gene_diff <- read.delim('/Users/yyq/Documents/MyDocument/project/image (Bioinformatics)/code/校服务器共享集群/sc_resol_st/Puck_180430_1/run-code/control_treat.glmLRT.txt', row.names = 1, sep = '\t', check.names = FALSE)
gene_diff <- lrt$table

# Sort the table by FDR value in ascending order. 
# If the FDR value is the same, continue to sort by log2FC in descending order.
gene_diff <- gene_diff[order(gene_diff$FDR, gene_diff$logFC, decreasing = c(FALSE, TRUE)), ]

# log2FC≥1 & FDR<0.01 mark "up", denote the up-regulated genes
# log2FC≤-1 & FDR<0.01 mark "down", denote the down-regulated genes
# else mark "none", denote the non-DE genes
gene_diff[which(gene_diff$logFC >= 1 & gene_diff$FDR < 0.01),'sig'] <- 'up'
gene_diff[which(gene_diff$logFC <= -1 & gene_diff$FDR < 0.01),'sig'] <- 'down'
gene_diff[which(abs(gene_diff$logFC) <= 1 | gene_diff$FDR >= 0.01),'sig'] <- 'none'

# Output the table of DE genes
gene_diff[which(gene_diff$FDR < 0.01),'onlyFDR'] <- 'Sig'
gene_diff[which(gene_diff$FDR >= 0.01),'onlyFDR'] <- 'none'
gene_diff_select_FDR <- subset(gene_diff, onlyFDR %in% 'Sig')

write.table(gene_diff_select_FDR, file = '../result_data/control_treat.glmLRT.select_FDR.txt', sep = '\t', col.names = NA, quote = FALSE)









