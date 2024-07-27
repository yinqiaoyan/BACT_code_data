##############################################################
#################       MERFISH_0.19       ###################
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
adata_obs <- read.csv("../input_data/adata_obs.csv")
truth_labels = adata_obs$cell_class

coord = read.csv("../input_data/coordinates.csv")
plot_color=c("#ff6466", "#f29d99", "#d2d1d0", "#a8dc93", "#6b1499",
             "#c599f3", "#138320", "#3185eb", "#9d766e", "#b2c7e5", 
             "#52c084", "#ffb610", "#FFE0C1", "#EFCDF7", "#d57fbe") 



##############################################################
#####              (1) Cell type annotation              #####
##############################################################
ppdata = data.frame(x = coord[, 1], y = coord[, 2], c = truth_labels)
pp = ggplot(data = ppdata, aes(x=x, y=y)) +
  geom_point(aes(color=factor(c)), size = 3) +
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

ggsave("../figures/merfish_0.19_truth_cell.png", pp, width = 22, height = 11, dpi = 100)



##############################################################
#####                     (2) BACT                       #####
##############################################################
## Result of BACT
tmpres = read.csv("../result_data/merfish_0.19_result_BACT.csv")
tmpc = tmpres$labels

## Compute ARI
cat("ARI value:", ARI(tmpc, truth_labels))

## Swap label numbers only for better visualization
tmpc2 = tmpc
tmpc_vec = sort(unique(tmpc))
tmpc2_vec = c(6,12,3,1,7,4,5,2,8,9,10,11,seq(13,15))  ## BACT
for (ii in 1:length(unique(tmpc))) {
  tmpc2[tmpc == tmpc_vec[ii]] = tmpc2_vec[[ii]]
}
tmpc = tmpc2

ppdata = data.frame(x = coord[, 1], y = coord[, 2], c = tmpc)
pp = ggplot(data = ppdata, aes(x=x, y=y)) +
  geom_point(aes(color=factor(c)), size = 3) +
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
                     labels = paste0("C", 1:15)) +
  guides(color = guide_legend(override.aes = list(size = 12), ncol = 2))
# pp

ggsave("../figures/merfish_0.19_bact.png", pp, width = 18, height = 10, dpi = 100)



##############################################################
#####                    (3) SpaGCN                      #####
##############################################################
# Result of SpaGCN
tmpres = read.csv("../result_data/merfish_0.19_result_SpaGCN.csv")
tmpc = tmpres$pred

## Compute ARI
cat("ARI value:", ARI(tmpc, truth_labels))

## Swap label numbers only for better visualization
tmpc2 = tmpc
tmpc_vec = sort(unique(tmpc))
tmpc2_vec = c(1,7,3,4,8,2,6,9,10,11,12,5,seq(13,15))  ## SpaGCN
for (ii in 1:length(unique(tmpc))) {
  tmpc2[tmpc == tmpc_vec[ii]] = tmpc2_vec[[ii]]
}
tmpc = tmpc2

ppdata = data.frame(x = coord[, 1], y = coord[, 2], c = tmpc)
pp = ggplot(data = ppdata, aes(x=x, y=y)) +
  geom_point(aes(color=factor(c)), size = 3) +
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
                     labels = paste0("C", 1:15)) +
  guides(color = guide_legend(override.aes = list(size = 12), ncol = 2))
# pp

ggsave("../figures/merfish_0.19_spagcn.png", pp, width = 18, height = 10, dpi = 100)



##############################################################
#####                    (4) STAGATE                     #####
##############################################################
# Result of STAGATE
tmpres = read.csv("../result_data/merfish_0.19_result_STAGATE.csv")
tmpc = tmpres$mclust

## Compute ARI
cat("ARI value:", ARI(tmpc, truth_labels))

## Swap label numbers only for better visualization
tmpc2 = tmpc
tmpc_vec = sort(unique(tmpc))
tmpc2_vec = c(1,3,13,8,9,2,6,4,10,5,7,11,12,14,15)  ## STAGATE
for (ii in 1:length(unique(tmpc))) {
  tmpc2[tmpc == tmpc_vec[ii]] = tmpc2_vec[[ii]]
}
tmpc = tmpc2

ppdata = data.frame(x = coord[, 1], y = coord[, 2], c = tmpc)
pp = ggplot(data = ppdata, aes(x=x, y=y)) +
  geom_point(aes(color=factor(c)), size = 3) +
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
                     labels = paste0("C", 1:15)) +
  guides(color = guide_legend(override.aes = list(size = 12), ncol = 2))
# pp

ggsave("../figures/merfish_0.19_stagate.png", pp, width = 18, height = 10, dpi = 100)



##############################################################
#####                    (5) BANKSY                      #####
##############################################################
# Result of STAGATE
tmpres = read.csv("../result_data/merfish_0.19_result_BANKSY.csv")
tmpc = tmpres$mclust

## Compute ARI
cat("ARI value:", ARI(tmpc, truth_labels))

## Swap label numbers only for better visualization
tmpc2 = tmpc
tmpc_vec = sort(unique(tmpc))
tmpc2_vec = c(10,15,3,2,8,5,13,9,7,12,6,1,11,14,4)  ## BANKSY
for (ii in 1:length(unique(tmpc))) {
  tmpc2[tmpc == tmpc_vec[ii]] = tmpc2_vec[[ii]]
}
tmpc = tmpc2

ppdata = data.frame(x = coord[, 1], y = coord[, 2], c = tmpc)
pp = ggplot(data = ppdata, aes(x=x, y=y)) +
  geom_point(aes(color=factor(c)), size = 3) +
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
                     labels = paste0("C", 1:15)) +
  guides(color = guide_legend(override.aes = list(size = 12), ncol = 2))
# pp

ggsave("../figures/merfish_0.19_banksy.png", pp, width = 18, height = 10, dpi = 100)



##############################################################
#####                     Figure1d                       #####
##############################################################
## Results of all methods for the five MERFISH datasets
ari_bact = c(0.382, 0.272, 0.366, 0.448, 0.496) 
ari_spagcn = c(0.274, 0.270, 0.242, 0.264, 0.298)
ari_stagate = c(0.054, 0.076, 0.077, 0.080, 0.073)
ari_banksy = c(0.119, 0.078, 0.094, 0.092, 0.069)


df = data.frame(ARI = c(ari_bact, ari_spagcn, ari_stagate, ari_banksy), 
                Method = rep(c("BACT", "SpaGCN", "STAGATE", "BANKSY"), each = 5))

df$Method <- factor(df$Method, levels = c("BACT", "SpaGCN", "STAGATE", "BANKSY"))
p <- ggplot(df, aes(x = Method, y = ARI)) + 
  stat_boxplot(geom = "errorbar", width=0.2, cex=0.5, coef = 500) +
  geom_boxplot(width=0.3, coef = 500) +
  ylim(c(0,0.6)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 37),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 33))

ggsave("../figures/merfish_boxplot.png", p, width = 12, height = 8, dpi = 100)






