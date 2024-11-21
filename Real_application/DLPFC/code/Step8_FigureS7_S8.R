##############################################################
#################       DLPFC 151507       ###################
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
truth_labels <- read.csv("../input_data/ground_truth_151507.csv")
truth_labels = truth_labels$ground_truth

coord = read.csv("../input_data/coordinates_151507.csv")
plot_color=c("#3d7bb0", "#ef8333", "#519f3e", "#e43b35", "#896bb0", "#945a52", 
             "#d27cbb", "#b2c7e5", "#a8dc93", "#f29d99", "#FFE0C1", "#EFCDF7", 
             "#b6bc6d", "#be9e96", "#24DC2E", "#006633", "#FF9933", "#3366CC", 
             "#211DC5", "#A47D83")



##############################################################
#####                     (b) BACT                       #####
##############################################################
## Result of BACT
tmpres = read.csv("../result_data/dlpfc151507_result_BACT.csv")
tmpc = tmpres$labels

## Compute ARI
cat("ARI value:", ARI(tmpc, truth_labels))

## Swap label numbers only for better visualization
tmpc2 = tmpc
tmpc_vec = sort(unique(tmpc))
tmpc2_vec = c(5,3,1,4,6,2,7,seq(8,13))  
for (ii in 1:length(unique(tmpc))) {
  tmpc2[tmpc == tmpc_vec[ii]] = tmpc2_vec[[ii]]
}
tmpc = tmpc2

ppdata = data.frame(x = coord[,2], y = -coord[,1], c = tmpc)
pp = ggplot(data = ppdata, aes(x=x, y=y)) +
  geom_point(aes(color=factor(c)), size = 3.2) +
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
                     labels = paste0("C", 1:length(unique(tmpc)))) +
  guides(color = guide_legend(override.aes = list(size = 12), ncol = 2))
# pp

ggsave("../figures/dlpfc151507_bact.png", pp, width = 11, height = 8, dpi = 100)



##############################################################
#####                    (c) SpaGCN                      #####
##############################################################
# Result of SpaGCN
tmpres = read.csv("../result_data/dlpfc151507_result_SpaGCN.csv")
tmpc = tmpres$pred

## Compute ARI
cat("ARI value:", ARI(tmpc, truth_labels))

## Swap label numbers only for better visualization
tmpc2 = tmpc
tmpc_vec = sort(unique(tmpc))
tmpc2_vec = c(5,2,4,7,1,6,3)
for (ii in 1:length(unique(tmpc))) {
  tmpc2[tmpc == tmpc_vec[ii]] = tmpc2_vec[[ii]]
}
tmpc = tmpc2

ppdata = data.frame(x = coord[, 2], y = -coord[, 1], c = tmpc)
pp = ggplot(data = ppdata, aes(x=x, y=y)) +
  geom_point(aes(color=factor(c)), size = 3.3) +
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
                     labels = paste0("C", 1:7)) +
  guides(color = guide_legend(override.aes = list(size = 12), ncol = 1))
# pp

ggsave("../figures/dlpfc151507_spagcn.png", pp, width = 9, height = 8, dpi = 100)



##############################################################
#####                    (d) STAGATE                     #####
##############################################################
# Result of STAGATE
tmpres = read.csv("../result_data/dlpfc151507_result_STAGATE.csv")
tmpc = tmpres$mclust

## Compute ARI
cat("ARI value:", ARI(tmpc, truth_labels))

## Swap label numbers only for better visualization
tmpc2 = tmpc
tmpc_vec = sort(unique(tmpc))
tmpc2_vec = c(4,1,2,3,7,6,5)
for (ii in 1:length(unique(tmpc))) {
  tmpc2[tmpc == tmpc_vec[ii]] = tmpc2_vec[[ii]]
}
tmpc = tmpc2

ppdata = data.frame(x = coord[, 2], y = -coord[, 1], c = tmpc)
pp = ggplot(data = ppdata, aes(x=x, y=y)) +
  geom_point(aes(color=factor(c)), size = 3.3) +
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
                     labels = paste0("C", 1:7)) +
  guides(color = guide_legend(override.aes = list(size = 12), ncol = 1))
# pp

ggsave("../figures/dlpfc151507_stagate.png", pp, width = 9, height = 8, dpi = 100)



##############################################################
#####                    (e) BANKSY                      #####
##############################################################
# Result of STAGATE
tmpres = read.csv("../result_data/dlpfc151507_result_BANKSY.csv")
tmpc = tmpres$mclust

## Compute ARI
cat("ARI value:", ARI(tmpc, truth_labels))

## Swap label numbers only for better visualization
tmpc2 = tmpc
tmpc_vec = sort(unique(tmpc))
tmpc2_vec = c(6,2,1,3,4,5,7)
for (ii in 1:length(unique(tmpc))) {
  tmpc2[tmpc == tmpc_vec[ii]] = tmpc2_vec[[ii]]
}
tmpc = tmpc2

ppdata = data.frame(x = coord[, 2], y = -coord[, 1], c = tmpc)
pp = ggplot(data = ppdata, aes(x=x, y=y)) +
  geom_point(aes(color=factor(c)), size = 3.3) +
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
                     labels = paste0("C", 1:7)) +
  guides(color = guide_legend(override.aes = list(size = 12), ncol = 1))
# pp

ggsave("../figures/dlpfc151507_banksy.png", pp, width = 9, height = 8, dpi = 100)



##############################################################
#####                     (f) BASS                       #####
##############################################################
# Result of STAGATE
tmpres = read.csv("../result_data/dlpfc151507_result_BASS.csv")
tmpc = tmpres$c

## Compute ARI
cat("ARI value:", ARI(tmpc, truth_labels))

## Swap label numbers only for better visualization
tmpc2 = tmpc
tmpc_vec = sort(unique(tmpc))
tmpc2_vec = c(5,1,3,4,7,6,2)
for (ii in 1:length(unique(tmpc))) {
  tmpc2[tmpc == tmpc_vec[ii]] = tmpc2_vec[[ii]]
}
tmpc = tmpc2

ppdata = data.frame(x = coord[, 2], y = -coord[, 1], c = tmpc)
pp = ggplot(data = ppdata, aes(x=x, y=y)) +
  geom_point(aes(color=factor(c)), size = 3.3) +
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
                     labels = paste0("C", 1:7)) +
  guides(color = guide_legend(override.aes = list(size = 12), ncol = 1))
# pp

ggsave("../figures/dlpfc151507_bass.png", pp, width = 9, height = 8, dpi = 100)



##############################################################
#####                     FigureS8                       #####
##############################################################
## Results of all methods based on five random repetitions
ari_bact = c(0.385, 0.429, 0.418, 0.402, 0.403)  
ari_spagcn = c(0.312, 0.427, 0.410, 0.413, 0.401)  
ari_stagate = c(0.501, 0.498, 0.520, 0.480, 0.505) 
ari_banksy = c(0.399, 0.303, 0.294, 0.396, 0.252) 
ari_bass = c(0.409, 0.356, 0.360, 0.385, 0.324) 


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

ggsave("../figures/FigureS8.png", p, width = 14, height = 8, dpi = 100)








