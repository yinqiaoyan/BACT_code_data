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
ari_spagcn = c(0.274, 0.326, 0.242, 0.264, 0.298)
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



##############################################################
#####                     FigureS1                       #####
##############################################################
data_id = c("04", "09", "14", "24")
plot_color=c("#ff6466", "#f29d99", "#d2d1d0", "#a8dc93", "#6b1499",
             "#c599f3", "#138320", "#3185eb", "#9d766e", "#b2c7e5", 
             "#52c084", "#ffb610", "#FFE0C1", "#EFCDF7", "#d57fbe", "#b6bc6d") 

tmpc2_list = vector("list", 4)
for (iii in 1:4) {
  tmpc2_list[[iii]] = vector("list", 4)
}

## MERFISH_0.04
tmpc2_list[[1]][[1]] = c(7,6,9,5,12,14,11,1,4,10,13,2,3,8)
tmpc2_list[[1]][[2]] = c(1,7,6,11,14,2,12,4,5,3,8,10,9,13)
tmpc2_list[[1]][[3]] = c(7,11,2,6,12,8,3,10,1,13,4,5,9,14)
tmpc2_list[[1]][[4]] = c(2,10,11,7,4,3,5,8,13,12,14,6,1,9)

## MERFISH_0.09
tmpc2_list[[2]][[1]] = c(11,7,1,5,12,2,6,4,3,8,9,10,13)
tmpc2_list[[2]][[2]] = c(7,1,14,6,4,2,5,3,8,12,10,15,9,11,13)
tmpc2_list[[2]][[3]] = c(2,3,7,14,8,12,6,5,1,10,11,9,4,13,15)
tmpc2_list[[2]][[4]] = c(6,5,7,2,3,9,8,4,15,12,10,1,11,13,14)

## MERFISH_0.14
tmpc2_list[[3]][[1]] = c(5,3,7,4,14,6,8,12,2,9,15,13,1,10,11,16)
tmpc2_list[[3]][[2]] = c(7,8,1,6,4,2,14,12,15,5,9,3,10,13,11)
tmpc2_list[[3]][[3]] = c(6,1,15,8,14,4,5,7,2,12,11,3,9,10,13)
tmpc2_list[[3]][[4]] = c(2,12,5,13,9,15,6,10,4,7,11,8,1,3,14)

## MERFISH_0.24
tmpc2_list[[4]][[1]] = c(8,7,1,12,5,6,14,4,2,11,3,9,10,13)
tmpc2_list[[4]][[2]] = c(1,7,2,11,6,3,14,12,9,4,8,5,13,15,10)
tmpc2_list[[4]][[3]] = c(7,2,3,14,9,10,5,11,12,8,13,1,15,6,4)
tmpc2_list[[4]][[4]] = c(8,2,3,11,6,9,5,4,10,13,7,14,12,1,15)



for (id in 1:4) {
  ### --- True labels --- ###
  coord = read.csv(paste0("../input_data/coordinates_0.", data_id[id], ".csv"))
  truth_labels = read.csv(paste0("../input_data/cell_type_annotation_0.", data_id[id], ".csv"))
  truth_labels = truth_labels$cell_class
  
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
  
  ggsave(paste0("../figures/merfish_0.", data_id[id], "_truth_cell.png"), pp, width = 22, height = 11, dpi = 100)
  
  
  
  ### --- BACT --- ###
  tmpres = read.csv(paste0("../result_data/merfish_0.", data_id[id], "_result_BACT.csv"))
  tmpc = tmpres$labels
  
  ## Compute ARI
  cat("ARI value:", ARI(tmpc, truth_labels))
  
  ## Swap label numbers only for better visualization
  tmpc2 = tmpc
  tmpc_vec = sort(unique(tmpc))
  tmpc2_vec = tmpc2_list[[id]][[1]]
  for (ii in 1:length(unique(tmpc))) {
    tmpc2[tmpc == tmpc_vec[ii]] = tmpc2_vec[ii]
  }
  tmpc = tmpc2
  
  ppdata = data.frame(x = coord[, 1], y = coord[, 2], c = tmpc)
  if (id == 4) {
    plot_color_tmp=c("#ff6466", "#f29d99", "#d2d1d0", "#a8dc93", "#6b1499",
                     "#c599f3", "#138320", "#3185eb", "#9d766e", "#b2c7e5", 
                     "#52c084", "#ffb610", "#FFE0C1", "#d57fbe") 
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
      scale_color_manual(values=plot_color_tmp,
                         labels = paste0("C", 1:16)) +
      guides(color = guide_legend(override.aes = list(size = 12), ncol = 2))
  } else {
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
                         labels = paste0("C", 1:16)) +
      guides(color = guide_legend(override.aes = list(size = 12), ncol = 2))
  }
  
  
  ggsave(paste0("../figures/merfish_0.", data_id[id], "_bact.png"), pp, width = 18, height = 10, dpi = 100)
  
  
  
  ### --- SpaGCN --- ###
  tmpres = read.csv(paste0("../result_data/merfish_0.", data_id[id], "_result_SpaGCN.csv"))
  tmpc = tmpres$pred
  
  ## Compute ARI
  cat("ARI value:", ARI(tmpc, truth_labels))
  
  ## Swap label numbers only for better visualization
  tmpc2 = tmpc
  tmpc_vec = sort(unique(tmpc))
  tmpc2_vec = tmpc2_list[[id]][[2]]
  for (ii in 1:length(unique(tmpc))) {
    tmpc2[tmpc == tmpc_vec[ii]] = tmpc2_vec[ii]
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
  
  ggsave(paste0("../figures/merfish_0.", data_id[id], "_spagcn.png"), pp, width = 18, height = 10, dpi = 100)
  
  
  
  ### --- STAGATE --- ###
  tmpres = read.csv(paste0("../result_data/merfish_0.", data_id[id], "_result_STAGATE.csv"))
  tmpc = tmpres$mclust
  
  ## Compute ARI
  cat("ARI value:", ARI(tmpc, truth_labels))
  
  ## Swap label numbers only for better visualization
  tmpc2 = tmpc
  tmpc_vec = sort(unique(tmpc))
  tmpc2_vec = tmpc2_list[[id]][[3]]
  for (ii in 1:length(unique(tmpc))) {
    tmpc2[tmpc == tmpc_vec[ii]] = tmpc2_vec[ii]
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
  
  ggsave(paste0("../figures/merfish_0.", data_id[id], "_stagate.png"), pp, width = 18, height = 10, dpi = 100)
  
  
  
  ### --- BANKSY --- ###
  tmpres = read.csv(paste0("../result_data/merfish_0.", data_id[id], "_result_BANKSY.csv"))
  tmpc = tmpres$mclust
  
  ## Compute ARI
  cat("ARI value:", ARI(tmpc, truth_labels))
  
  ## Swap label numbers only for better visualization
  tmpc2 = tmpc
  tmpc_vec = sort(unique(tmpc))
  tmpc2_vec = tmpc2_list[[id]][[4]]
  for (ii in 1:length(unique(tmpc))) {
    tmpc2[tmpc == tmpc_vec[ii]] = tmpc2_vec[ii]
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
  
  ggsave(paste0("../figures/merfish_0.", data_id[id], "_banksy.png"), pp, width = 18, height = 10, dpi = 100)
  
}





