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
#####              (b) Cell type annotation              #####
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
#####                     (c) BACT                       #####
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
#####                    (d) SpaGCN                      #####
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
                     labels = paste0("C", 1:16)) +
  guides(color = guide_legend(override.aes = list(size = 12), ncol = 2))
# pp

ggsave("../figures/starmap_spagcn.png", pp, width = 18, height = 10, dpi = 100)



##############################################################
#####                    (e) STAGATE                     #####
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
                     labels = paste0("C", 1:16)) +
  guides(color = guide_legend(override.aes = list(size = 12), ncol = 2))
# pp

ggsave("../figures/starmap_stagate.png", pp, width = 18, height = 10, dpi = 100)



##############################################################
#####                    (f) BANKSY                      #####
##############################################################
# Result of BANKSY
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
                     labels = paste0("C", 1:16)) +
  guides(color = guide_legend(override.aes = list(size = 12), ncol = 2))
# pp

ggsave("../figures/starmap_banksy.png", pp, width = 18, height = 10, dpi = 100)



##############################################################
#####                     (g) BASS                       #####
##############################################################
# Result of BASS
tmpres = read.csv("../result_data/starmap_result_BASS.csv")
tmpc = tmpres$c

## Compute ARI
cat("ARI value:", ARI(tmpc, truth_labels))

## Swap label numbers only for better visualization
tmpc2 = tmpc
tmpc_vec = sort(unique(tmpc))
tmpc2_vec = c(9,11,13,8,5,4,6,7,3,1,12,15,10,2,14,16)  ## BASS
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
                     labels = paste0("C", 1:16)) +
  guides(color = guide_legend(override.aes = list(size = 12), ncol = 2))
# pp

ggsave("../figures/starmap_bass.png", pp, width = 18, height = 10, dpi = 100)



##############################################################
#####                 (h) ARI boxplots                   #####
##############################################################
## Results of all methods based on the five random repetitions
ari_bact = c(0.561, 0.482, 0.553, 0.629, 0.539)  
ari_spagcn = c(0.432, 0.425, 0.470, 0.418, 0.435) 
ari_stagate = c(0.307, 0.312, 0.287, 0.302, 0.296)  
ari_banksy = c(0.128, 0.100, 0.111, 0.110, 0.117) 
ari_bass = c(0.295, 0.277, 0.220, 0.222, 0.278) 

df = data.frame(ARI = c(ari_bact, ari_spagcn, ari_stagate, ari_banksy, ari_bass), 
                Method = rep(c("BACT", "SpaGCN", "STAGATE", "BANKSY", "BASS"), each = 5))

df$Method <- factor(df$Method, levels = c("BACT", "SpaGCN", "STAGATE", "BANKSY", "BASS"))
p <- ggplot(df, aes(x = Method, y = ARI)) + 
  stat_boxplot(geom = "errorbar", width=0.2, cex=0.5, coef = 500) +
  geom_boxplot(width=0.3, coef = 500) +
  ylim(c(0,0.65)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 37),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 33))

ggsave("../figures/starmap_boxplots.png", p, width = 14, height = 8, dpi = 100)









