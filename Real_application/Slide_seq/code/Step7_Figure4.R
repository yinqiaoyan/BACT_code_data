##############################################################
#################        Slide-seq         ###################
##############################################################

########################## Note ##############################
## Please set the working directory to the source file 
## location.
##############################################################

library(ggplot2)

## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())

coord = read.csv("../input_data/coordinates.csv")
plot_color=c("#78A6DE", "#ffb610", "#ff6466", "#6b1499", "#7b92ce",
             "#d2d1d0", "#138320", "#3185eb", "#9d766e", "#b2c7e5", 
             "#a8dc93", "#f29d99", "#FFE0C1", "#EFCDF7")



##############################################################
#####                     (b) BACT                       #####
##############################################################
## Result of BACT
tmpres = read.csv("../result_data/slideseq_result_BACT.csv")
tmpc = tmpres$labels

ppdata = data.frame(x = coord[, 1], y = coord[, 2], c = tmpc)
pp = ggplot(data = ppdata, aes(x=x, y=y)) +
  geom_point(aes(color=factor(c)), size = 2.2) +
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
                     labels = paste0("C", 1:14)) +
  guides(color = guide_legend(override.aes = list(size = 12), ncol = 2))
# pp

ggsave("../figures/slideseq_bact.png", pp, width = 18, height = 10, dpi = 100)



##############################################################
#####                    (c) SpaGCN                      #####
##############################################################
# Result of SpaGCN
tmpres = read.csv("../result_data/slideseq_result_SpaGCN.csv")
tmpc = tmpres$pred

## Swap label numbers only for better visualization
tmpc2 = tmpc
tmpc_vec = sort(unique(tmpc))
tmpc2_vec = c(1,3,2,4)
for (ii in 1:length(unique(tmpc))) {
  tmpc2[tmpc == tmpc_vec[ii]] = tmpc2_vec[[ii]]
}
tmpc = tmpc2

ppdata = data.frame(x = coord[, 1], y = coord[, 2], c = tmpc)
pp = ggplot(data = ppdata, aes(x=x, y=y)) +
  geom_point(aes(color=factor(c)), size = 2.2) +
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
                     labels = paste0("C", 1:4)) +
  guides(color = guide_legend(override.aes = list(size = 12)))
# pp

ggsave("../figures/slideseq_spagcn.png", pp, width = 16.5, height = 10, dpi = 100)



##############################################################
#####                    (d) STAGATE                     #####
##############################################################
# Result of STAGATE
tmpres = read.csv("../result_data/slideseq_result_STAGATE.csv")
tmpc = tmpres$mclust

## Swap label numbers only for better visualization
tmpc2 = tmpc
tmpc_vec = sort(unique(tmpc))
tmpc2_vec = c(1,2,4,3)
for (ii in 1:length(unique(tmpc))) {
  tmpc2[tmpc == tmpc_vec[ii]] = tmpc2_vec[[ii]]
}
tmpc = tmpc2

ppdata = data.frame(x = coord[, 1], y = coord[, 2], c = tmpc)
pp = ggplot(data = ppdata, aes(x=x, y=y)) +
  geom_point(aes(color=factor(c)), size = 2.2) +
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
                     labels = paste0("C", 1:4)) +
  guides(color = guide_legend(override.aes = list(size = 12)))
# pp

ggsave("../figures/slideseq_stagate.png", pp, width = 16.5, height = 10, dpi = 100)



##############################################################
#####                    (e) BANKSY                      #####
##############################################################
# Result of STAGATE
tmpres = read.csv("../result_data/slideseq_result_BANKSY.csv")
tmpc = tmpres$mclust

## Swap label numbers only for better visualization
tmpc2 = tmpc
tmpc_vec = sort(unique(tmpc))
tmpc2_vec = c(3,1,2,4)
for (ii in 1:length(unique(tmpc))) {
  tmpc2[tmpc == tmpc_vec[ii]] = tmpc2_vec[[ii]]
}
tmpc = tmpc2

ppdata = data.frame(x = coord[, 1], y = coord[, 2], c = tmpc)
pp = ggplot(data = ppdata, aes(x=x, y=y)) +
  geom_point(aes(color=factor(c)), size = 2.2) +
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
                     labels = paste0("C", 1:4)) +
  guides(color = guide_legend(override.aes = list(size = 12)))
# pp

ggsave("../figures/slideseq_banksy.png", pp, width = 16.5, height = 10, dpi = 100)



##############################################################
#####                     (f) BASS                       #####
##############################################################
# Result of STAGATE
tmpres = read.csv("../result_data/slideseq_result_BASS.csv")
tmpc = tmpres$c

## Swap label numbers only for better visualization
tmpc2 = tmpc
tmpc_vec = sort(unique(tmpc))
tmpc2_vec = c(1,3,2,4)
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
                     labels = paste0("C", 1:4)) +
  guides(color = guide_legend(override.aes = list(size = 12)))
# pp

ggsave("../figures/slideseq_bass.png", pp, width = 18, height = 10, dpi = 100)












