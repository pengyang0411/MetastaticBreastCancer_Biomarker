##-------------------------------
## Figure Panel S3
## Author: Peng Yang
##-------------------------------
library(ggplot2)
library(dplyr)
library(Seurat)
library(ggpubr)
library(tibble)
library(fgsea)
library(viridis)
library(RColorBrewer)
library(hrbrthemes)
library(DoMultiBarHeatmap)
library("ComplexHeatmap")
library(rstatix)
library(dittoSeq)
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
library(EnhancedVolcano)

load('/Data/combined_nonTumor_final.RData')

##-------------------------------
## Figure S3A-C cellchat
##-------------------------------

for(s in 1:3){
  
  load('CellChat/cell_chat_matched_', c('Baseline', 'Early Progressor', 'Late Progressor')[s], '.RData')
  
  pdf(paste0('FigureS3/FigureS3', s, '.pdf'), width = 8, height = 8)
  print(netVisual_chord_gene(cellchat, sources.use = NULL, targets.use = NULL, signaling = pathways.show.all, legend.pos.x = 8,
                             scale = T, thresh = 0.05,
                             title.name = s))
  dev.off()
  
}

##-------------------------------
## Figure S3D: Umap plot across 
## different sample status
##-------------------------------

Idents(combined_nontumor) <- combined_nontumor$CT_L2
combined_nontumor$CT_L2.num <- as.factor(paste(as.numeric(as.factor(combined_nontumor$CT_L2))))

pt3d <- DimPlot(combined_nontumor, reduction = 'umap.harmony', group.by = 'CT_L2', 
                split.by = 'Status_new', 
                label = F, repel = T, label.size = 5) + 
  NoAxes() + ggtitle(' ') + 
  scale_x_reverse() +
  theme_void() + 
  # scale_color_manual(# breaks = unique(merged_combined$orig.ident.num),
  #   labels = paste(1:length(unique(combined_nontumor$CT_L1)), as.character(combined_nontumor$CT_L2)), 
  #   values = colorRampPalette(33)) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 6),
                              label.hjust = 0,
                              label.theme = element_text())) + #NoLegend() +
  theme(strip.text.x = element_blank())

pt3d

pdf('FigureS3/FigureS3d.pdf', width = 16, height = 5)
print(pt3d)
dev.off()



##-------------------------------
## Figure 3E: Marker gene plot
##-------------------------------

Idents(combined_nontumor) <- combined_nontumor$CT_L2

if(file.exists('FigureS3/Data/cell_markers.RData')){
  load('FigureS3/Data/cell_markers.RData')
}else{
  merged_combined.maker <- FindAllMarkers(combined_nontumor, only.pos = TRUE, 
                                          min.pct = 0.10, logfc.threshold = 0.25)
  save(merged_combined.maker, file = 'FigureS3/Data/cell_markers.RData')
  
  write.csv(merged_combined.maker, file = 'FigureS3/Data/cell_markers.csv')
  
}

merged_combined.maker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

## Make the heatmap for marker genes
merged_combined.maker %>%
  group_by(cluster) %>%
  top_n(n = 1, wt = avg_log2FC) -> top1



FigureS3e <- DotPlot(combined_nontumor, features = unique(top1$gene)) + 
  theme(axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1, size = 10)) + 
  ylab(' ') + xlab(' ') + 
  # scale_y_discrete(labels = c('LP', 'EP', 'BL'),
  #                  breaks = c('Late Progressor', 'Early Progressor', 'Baseline')) +
  scale_color_gradientn(
    colors = c("#5DBCFF", "#6DCCFF", "white", "#F494BE", "#F484AE")) +
  theme(
    text = element_text(size = 10),
    panel.grid.major = element_line(colour = "grey90", size=0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
    # axis.text.y=element_text(angle = -30, vjust = 1, hjust = 1),
    legend.position="top",
    legend.title=element_text(size=10)) +
  guides(size = guide_legend(title.position="top", 
                             title = 'Percentage of Expression',
                             title.hjust = 0.3,
                             ncol = 6,
                             byrow = T,
                             override.aes = list(stroke = 0.4, ncol = 3)),
         color = guide_colourbar(title.position = "left", title = 'Average Expression',
                                 title.vjust = 0.8,
                                 barwidth = 10))


pdf(paste0('FigureS3/Figure3e.pdf'), width = 8, height = 8)
print(FigureS3e)
dev.off()


##-------------------------------
## Figure 3g: Relative other T
## cell distribution
##-------------------------------

## Create the dataframe
df_nontumor_L2 <- c()
for(sample in unique(combined_nontumor$orig.ident)){
  df_nontumor_L2 <- rbind(df_nontumor_L2, 
                          cbind(cbind(table(combined_nontumor$CT_L2[combined_nontumor$orig.ident == sample]), sample),
                                as.character(combined_nontumor$Status_new[combined_nontumor$orig.ident == sample][1])))
}

df_nontumor_L2 <- cbind(df_nontumor_L2, rownames(df_nontumor_L2))

colnames(df_nontumor_L2) <- c('Count', 'Sample', 'Status', 'CellTypes')
df_nontumor_L2 <- as.data.frame(df_nontumor_L2)

## Make factors
df_nontumor_L2$Count <- as.numeric(df_nontumor_L2$Count)
df_nontumor_L2$Sample <- as.factor(df_nontumor_L2$Sample)
df_nontumor_L2$Status <- as.factor(df_nontumor_L2$Status)
df_nontumor_L2$CellTypes <- as.factor(df_nontumor_L2$CellTypes)

cells_2 <- c('gdT', 'dnT', 'MAIT')

df_nontumor_L2_frac <- df_nontumor_L2
for(i in unique(df_nontumor_L2$Sample)){
  
  Indx <- which(df_nontumor_L2$Sample == i)
  
  df_nontumor_L2_frac[Indx, 'Count'] = df_nontumor_L2[Indx, 'Count'] / sum(df_nontumor_L2[Indx, 'Count'])
  
}

df_nontumor_L2_frac_sub <- df_nontumor_L2_frac[df_nontumor_L2_frac$CellTypes %in% cells_2, ]

category <- rep(NA, nrow(df_nontumor_L2_frac_sub))
category[df_nontumor_L2_frac_sub$CellTypes %in% c('gdT', 'dnT', 'MAIT')] <- 'other T'

df_nontumor_L2_frac_sub$category <- as.factor(category)

str(df_nontumor_L2_frac_sub)
## Adding stat test
stat.test <- df_nontumor_L2_frac_sub %>%
  group_by(category, CellTypes) %>%
  t_test(Count ~ Status) %>%
  # adjust_pvalue(method = "BH") %>%
  add_significance("p")
stat.test$p.adj.signif <- stat.test$p.signif
stat.test

# Add p-values onto the box plots
stat.test <- stat.test %>%
  add_xy_position(x = "CellTypes", dodge = 0.5)

dodge <- position_dodge(width = 0.5)
pt_nontumor_L2_frac_sub <- ggplot(df_nontumor_L2_frac_sub) + 
  geom_boxplot(size = 0.3, outlier.shape = NA, width = 0.618, lwd = 0.5,
               aes(y = Count, x = CellTypes, color = Status)) +
  facet_wrap(~category, scales = 'free') + 
  geom_jitter(shape = 16, position = position_jitterdodge(), 
              size = 2,  alpha = 0.9, aes(fill = Status, y = Count, x = CellTypes, color = Status)) +
  theme_classic() + 
  ggtitle(" ") +
  xlab(" ") + ylab("Cell Fraction")  +
  # geom_pwc(
  #   aes(group = Status), tip.length = 0,
  #   method = "wilcox_test", label = "p.adj.format",
  #   bracket.nudge.y = -0.08,
  # )  +
  # stat_compare_means(aes(group = Status), comparisons = my_comparisons) +
  # stat_compare_means(aes(group = Status), label = "p.signif", method = "t.test") +
  theme(legend.position="bottom",
        # legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        panel.border = element_blank(),
        strip.background = element_blank(),
        # strip.text.x = element_blank(),
        # legend.title = element_text(color = "transparent"),
        # legend.text = element_text(color = "transparent"),
  ) + guides(color=guide_legend(title = 'Sample Status'),
             fill = FALSE)   +
  scale_fill_manual(values=c('#cccc66', '#66cccc', '#ff9966'))  +
  scale_color_manual(values=c('#cccc66', '#66cccc', '#ff9966')) #+
  # stat_pvalue_manual(
  #   stat.test,  label = "p.signif", tip.length = 0.01, hide.ns = T,
  # )

pdf('FigureS3/Figure3g.pdf', width = 10, height = 8)
print(pt_nontumor_L2_frac_sub)
dev.off()




