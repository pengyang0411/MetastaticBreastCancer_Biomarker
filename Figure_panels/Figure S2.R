##-------------------------------
## Figure Panel S2
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

library(Signac)
library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(EnhancedVolcano)
load('/Data/merged_data_Harmony_new.RData')
load('/Users/py11/Box Sync/CDK4/Data/merged_data_Harmony_new.RData')

## Subset of tumor cells
Idents(merged_combined) <- merged_combined$Tumor_cells
merged_combined <- subset(merged_combined, idents = 'tumor cells')


##-------------------------------
## Figure S2A: Heatmap compares
## sample status
##-------------------------------
Idents(merged_combined) <- merged_combined$metastatic
# merged_combined <- ScaleData(object = merged_combined, features = rownames(merged_combined))

## For heatmap
merged_combined.maker <- FindAllMarkers(merged_combined, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
merged_combined.maker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
head(merged_combined.maker)

## Make the heatmap for marker genes
merged_combined.maker %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

pdf(paste0('FigureS2/FigureS2a.pdf'), width = 8, height = 6)
DoMultiBarHeatmap(merged_combined, features = top10$gene, group.by = 'Status_new',
                  additional.group.by = c('orig.ident', 'metastatic'), label = F)  + scale_fill_viridis()
dev.off()



##-------------------------------
## Figure S2B: Volcano plot
## between pleural and liver
## for baseline samples
##-------------------------------
Idents(merged_combined) <- merged_combined$Status_new
combined_BL <-  subset(merged_combined, idents = 'Baseline')

Idents(combined_BL) <- combined_BL$metastatic

if(file.exists(paste0('FigureS2/Data/DE_tumor_BL_Pleural_Liver.csv'))){
  thisDE <- read.csv(paste0('FigureS2/Data/DE_tumor_BL_Pleural_Liver.csv'))
}else{
  thisDE <- FindMarkers(object = combined_BL, 
                        ident.1 = 'Liver', ident.2 = 'Pleural effusion', 
                        min.pct = 0, test.use = "MAST")
  thisDE$gene <- rownames(thisDE)
  write.csv(thisDE, file = paste0('FigureS2/Data/DE_tumor_BL_Pleural_Liver.csv'))
}

# thisDE$avg_log2FC <- thisDE$avg_log2FC + rnorm(length(thisDE$avg_log2FC ), 0, 0.1)

FigureS2b <- EnhancedVolcano(thisDE,
                             lab = thisDE$gene,
                             # selectLab = selected_genes,
                             x = 'avg_log2FC',
                             y = 'p_val_adj',
                             title = paste0('Pleural effusion (Left) vs Liver (Right)'),
                             pCutoff = max(thisDE$p_val_adj),
                             FCcutoff = 1,
                             pointSize = 2.0,
                             labSize = 4.0,
                             subtitle = NULL,
                             # colors = c('#FD9AA0', '#6DCCFD'),
                             drawConnectors = TRUE,
                             legendPosition = 'right',
                             xlim = c(-3, 3)) + theme_bw() + 
  theme(
    text = element_text(size = 12),
    panel.grid.major = element_line(colour = "grey90", size=0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position="none",
    legend.title=element_text(size=10)) + 
  guides(size = guide_legend(title.position="bottom", #title = 'P values',
                             title.hjust = 0.5,
                             ncol = 4,
                             byrow = T,
                             override.aes = list(stroke = 0.4, ncol = 3)),
         fill = guide_colourbar(title.position = "bottom", title.hjust = 0.5)) + 
  scale_color_manual(values = c('white', '#6DCCFD', '#FD9AA0'))


pdf(paste0('FigureS2/FigureS2b.pdf'), width = 5, height = 4)
print(FigureS2b)
dev.off()



##-------------------------------
## Figure S2C: Volcano plot
## between pleural and liver
## for late progressor samples
##-------------------------------
Idents(merged_combined) <- merged_combined$Status_new
combined_LP <-  subset(merged_combined, idents = 'Late Progressor')

Idents(combined_LP) <- combined_LP$metastatic

if(file.exists(paste0('FigureS2/Data/DE_tumor_LP_Pleural_Liver.csv'))){
  thisDE <- read.csv(paste0('FigureS2/Data/DE_tumor_LP_Pleural_Liver.csv'))
}else{
  thisDE <- FindMarkers(object = combined_LP, 
                        ident.1 = 'Liver', ident.2 = 'Pleural effusion', 
                        min.pct = 0, test.use = "MAST")
  thisDE$gene <- rownames(thisDE)
  write.csv(thisDE, file = paste0('FigureS2/Data/DE_tumor_LP_Pleural_Liver.csv'))
}

# thisDE$avg_log2FC <- thisDE$avg_log2FC + rnorm(length(thisDE$avg_log2FC ), 0, 0.1)

FigureS2c <- EnhancedVolcano(thisDE,
                             lab = thisDE$gene,
                             # selectLab = selected_genes,
                             x = 'avg_log2FC',
                             y = 'p_val_adj',
                             title = paste0('Pleural effusion (Left) vs Liver (Right)'),
                             pCutoff = max(thisDE$p_val_adj),
                             FCcutoff = 1,
                             pointSize = 2.0,
                             labSize = 4.0,
                             subtitle = NULL,
                             # colors = c('#FD9AA0', '#6DCCFD'),
                             drawConnectors = TRUE,
                             legendPosition = 'right',
                             xlim = c(-3, 3)) + theme_bw() + 
  theme(
    text = element_text(size = 12),
    panel.grid.major = element_line(colour = "grey90", size=0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position="none",
    legend.title=element_text(size=10)) + 
  guides(size = guide_legend(title.position="bottom", #title = 'P values',
                             title.hjust = 0.5,
                             ncol = 4,
                             byrow = T,
                             override.aes = list(stroke = 0.4, ncol = 3)),
         fill = guide_colourbar(title.position = "bottom", title.hjust = 0.5)) + 
  scale_color_manual(values = c('white', 'white', '#6DCCFD', '#FD9AA0'))


pdf(paste0('FigureS2/FigureS2c.pdf'), width = 5, height = 4)
print(FigureS2c)
dev.off()




##-------------------------------
## Figure S2C: Heatmap compares
## all
##-------------------------------
Idents(merged_combined) <- merged_combined$metastatic
combined_Pleural <-  subset(merged_combined, idents = 'Pleural effusion')


Idents(combined_Pleural) <- combined_Pleural$Status_new
## For heatmap
if(file.exists(paste0('FigureS2/Data/DE_tumor_Pleural_all.RData'))){
  load('FigureS2/Data/DE_tumor_Pleural_all.RData')
}else{
  combined_Pleural.maker <- FindAllMarkers(combined_Pleural, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
  save(combined_Pleural.maker, file = 'FigureS2/Data/DE_tumor_Pleural_all.RData')
}
combined_Pleural.maker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
head(combined_Pleural.maker)

## Make the heatmap for marker genes
combined_Pleural.maker %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

# pdf(paste0('Figure2/Figure2b.pdf'), width = 8, height = 6)
# DoMultiBarHeatmap(merged_combined, features = top10$gene, group.by = 'Status_new', 
#                   additional.group.by = c('orig.ident', 'metastatic'), label = F)  + scale_fill_viridis() 
# dev.off()

Idents(combined_Pleural) <- combined_Pleural$Status_new
FigureS2d <- DotPlot(combined_Pleural, features = top10$gene) + 
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
                             # title.hjust = 0.3, 
                             ncol = 6,
                             byrow = T,
                             override.aes = list(stroke = 0.4, ncol = 3)),
         color = guide_colourbar(title.position = "left", title = 'Average Expression',
                                 title.vjust = 0.8,
                                 barwidth = 10))

pdf(paste0('FigureS2/FigureS2d.pdf'), width = 8, height = 3.5)
print(FigureS2d)
dev.off()

