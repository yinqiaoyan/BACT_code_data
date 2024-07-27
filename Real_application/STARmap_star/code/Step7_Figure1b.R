##############################################################
###################       STARmap*       #####################
##############################################################

########################## Note ##############################
## Please set the working directory to the source file 
## location.
##############################################################

library(aricode)
library(ggplot2)

## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())

## True labels
truth_labels <- read.table("../input_data/cell_annotation_file.txt", header = TRUE, sep = "")
truth_labels = truth_labels$Annotation

coord = read.csv("../input_data/coordinates.csv")
plot_color=c("#ff6466", "#ffb610", "#c599f3", "#52c084", "#7b92ce", "#d2d1d0", 
             "#6b1499", "#138320", "#3185eb", "#9d766e", "#b2c7e5", "#a8dc93", 
             "#f29d99", "#FFE0C1", "#EFCDF7", "#d57fbe", "#b6bc6d", "#be9e96")



##############################################################
#####              (1) Cell type annotation              #####
##############################################################
ppdata = data.frame(x = coord[, 1], y = coord[, 2], c = truth_labels)
pp = ggplot(data = ppdata, aes(x=x, y=y)) +
  geom_point(aes(color=factor(c)), size = 6) +
  theme(panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.text = element_text(size = 35),
        legend.title = element_blank(),
        legend.key.size = unit(2, 'cm')) +
  scale_color_manual(values=plot_color) +
  guides(color = guide_legend(override.aes = list(size = 12), ncol = 2))
# pp

ggsave("../figures/starmap_truth_cell.png", pp, width = 20, height = 10, dpi = 100)



##############################################################
#####                     (2) BACT                       #####
##############################################################
## Result of BACT
tmpres = read.csv("../result_data/starmap_result_BACT.csv")
tmpc = tmpres$labels

## Compute ARI
cat("ARI value:", ARI(tmpc, truth_labels))

## Swap label numbers only for better visualization
tmpc2 = tmpc
tmpc_vec = sort(unique(tmpc))
tmpc2_vec = c(13,4,7,2,5,3,11,8,10,18,1,6,9,12,14,15,16,17)  ## BACT
for (ii in 1:length(unique(tmpc))) {
  tmpc2[tmpc == tmpc_vec[ii]] = tmpc2_vec[[ii]]
}
tmpc = tmpc2

ppdata = data.frame(x = coord[, 1], y = coord[, 2], c = tmpc)
pp = ggplot(data = ppdata, aes(x=x, y=y)) +
  geom_point(aes(color=factor(c)), size = 6) +
  theme(panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.text = element_text(size = 35),
        legend.title = element_blank(),
        legend.key.size = unit(2, 'cm')) +
  scale_color_manual(values=plot_color,
                     labels = paste0("C", 1:18)) +
  guides(color = guide_legend(override.aes = list(size = 12), ncol = 2))
# pp

ggsave("../figures/starmap_bact.png", pp, width = 18, height = 10, dpi = 100)



##############################################################
#####                    (3) SpaGCN                      #####
##############################################################
# Result of SpaGCN
tmpres = read.csv("../result_data/starmap_result_SpaGCN.csv")
tmpc = tmpres$pred

## Compute ARI
cat("ARI value:", ARI(tmpc, truth_labels))

## Swap label numbers only for better visualization
tmpc2 = tmpc
tmpc_vec = sort(unique(tmpc))
tmpc2_vec = c(8,7,15,5,6,9,16,13,10,2,4,11,14,1,12,3)  ## SpaGCN
for (ii in 1:length(unique(tmpc))) {
  tmpc2[tmpc == tmpc_vec[ii]] = tmpc2_vec[[ii]]
}
tmpc = tmpc2

ppdata = data.frame(x = coord[, 1], y = coord[, 2], c = tmpc)
pp = ggplot(data = ppdata, aes(x=x, y=y)) +
  geom_point(aes(color=factor(c)), size = 6) +
  theme(panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.text = element_text(size = 35),
        legend.title = element_blank(),
        legend.key.size = unit(2, 'cm')) +
  scale_color_manual(values=plot_color,
                     labels = paste0("C", 1:18)) +
  guides(color = guide_legend(override.aes = list(size = 12), ncol = 2))
# pp

ggsave("../figures/starmap_spagcn.png", pp, width = 18, height = 10, dpi = 100)



##############################################################
#####                    (4) STAGATE                     #####
##############################################################
# Result of STAGATE
tmpres = read.csv("../result_data/starmap_result_STAGATE.csv")
tmpc = tmpres$mclust

## Compute ARI
cat("ARI value:", ARI(tmpc, truth_labels))

## Swap label numbers only for better visualization
tmpc2 = tmpc
tmpc_vec = sort(unique(tmpc))
tmpc2_vec = c(7,4,3,1,9,6,2,8,11,10,5,seq(12,16))  ## STAGATE
for (ii in 1:length(unique(tmpc))) {
  tmpc2[tmpc == tmpc_vec[ii]] = tmpc2_vec[[ii]]
}
tmpc = tmpc2

ppdata = data.frame(x = coord[, 1], y = coord[, 2], c = tmpc)
pp = ggplot(data = ppdata, aes(x=x, y=y)) +
  geom_point(aes(color=factor(c)), size = 6) +
  theme(panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.text = element_text(size = 35),
        legend.title = element_blank(),
        legend.key.size = unit(2, 'cm')) +
  scale_color_manual(values=plot_color,
                     labels = paste0("C", 1:18)) +
  guides(color = guide_legend(override.aes = list(size = 12), ncol = 2))
# pp

ggsave("../figures/starmap_stagate.png", pp, width = 18, height = 10, dpi = 100)



##############################################################
#####                    (5) BANKSY                      #####
##############################################################
# Result of STAGATE
tmpres = read.csv("../result_data/starmap_result_BANKSY.csv")
tmpc = tmpres$mclust

## Compute ARI
cat("ARI value:", ARI(tmpc, truth_labels))

## Swap label numbers only for better visualization
tmpc2 = tmpc
tmpc_vec = sort(unique(tmpc))
tmpc2_vec = c(8,3,2,11,7,4,1,6,10,9,5,seq(12,16))  ## BANKSY
for (ii in 1:length(unique(tmpc))) {
  tmpc2[tmpc == tmpc_vec[ii]] = tmpc2_vec[[ii]]
}
tmpc = tmpc2

ppdata = data.frame(x = coord[, 1], y = coord[, 2], c = tmpc)
pp = ggplot(data = ppdata, aes(x=x, y=y)) +
  geom_point(aes(color=factor(c)), size = 6) +
  theme(panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.text = element_text(size = 35),
        legend.title = element_blank(),
        legend.key.size = unit(2, 'cm')) +
  scale_color_manual(values=plot_color,
                     labels = paste0("C", 1:18)) +
  guides(color = guide_legend(override.aes = list(size = 12), ncol = 2))
# pp

ggsave("../figures/starmap_banksy.png", pp, width = 18, height = 10, dpi = 100)









