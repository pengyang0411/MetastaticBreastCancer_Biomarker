##-------------------------------
## Figure Panel 5
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
library(scales)
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
library(EnhancedVolcano)

load('/Data/combined_nonTumor_final.RData')

##-------------------------------------
## Figure S5A: Tsne plot
##-------------------------------------


Idents(combined_nontumor) <- combined_nontumor$CT_L1
combined_Myeloid <- subset(combined_nontumor, idents = c('Mono', 'DC'))


Idents(combined_Myeloid) <- combined_Myeloid$CT_L2
combined_Myeloid <- subset(combined_Myeloid, 
                           idents = c('CD16 Mono', 'cDC1', 'cDC2', 'pDC',
                                      'Macrophage'))

combined_Myeloid <- RunTSNE(combined_Myeloid, reduction = 'umap.pca', dims.use = 1:10, perplexity = 200)


ptS5A <- DimPlot(combined_Myeloid, reduction = 'tsne', group.by = 'RNA_snn_res.0.05', 
                 # split.by = 'Status_new', 
                 label = T, repel = T, label.size = 5) + 
  NoAxes() + ggtitle(' ') + 
  scale_x_reverse() +
  theme_void() + 
  # scale_color_manual(# breaks = unique(merged_combined$orig.ident.num),
  #   labels = paste(1:length(unique(combined_Tcells$CT_L2)), T_cells),
  #   values = umapColor[1:13]) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 6),
                              label.hjust = 0,
                              label.theme = element_text())) + #NoLegend() +
  theme(strip.text.x = element_blank())

pdf('FigureS5/FigureS5A.pdf', width = 7, height = 5)
print(ptS5A)
dev.off()

##-------------------------------------
## Figure S5A: Heatmap
##-------------------------------------

Myeloid_markders <- c(readxl::read_excel('FigureS5/Data/Myeloid_markders.xlsx', sheet = 1)$Gene,
                      readxl::read_excel('FigureS5/Data/Myeloid_markders.xlsx', sheet = 2)$Gene)

Idents(combined_Myeloid) <- combined_Myeloid$RNA_snn_res.0.05


FigureS5B <-  DotPlot(combined_Myeloid, features = unique(Myeloid_markders)) + 
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


pdf(paste0('FigureS5B/FigureS5B.pdf'), width = 10, height = 5)
print(FigureS5B)
dev.off()


##-------------------------------------
## Figure S5C: Heatmap of DE genes 
## Across different metastatic site
## for baseline samples
##-------------------------------------

Idents(combined_Myeloid) <- combined_Myeloid$RNA_snn_res.0.05

if(file.exists('FigureS5/Data/cell_markers_tsne.RData')){
  load('FigureS5/Data/cell_markers_tsne.RData')
}else{
  merged_combined.maker <- FindAllMarkers(combined_Myeloid, only.pos = TRUE, 
                                          min.pct = 0.10, logfc.threshold = 0.25)
  save(merged_combined.maker, file = 'FigureS5/Data/cell_markers_tsne.RData')
  
}

merged_combined.maker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

## Make the heatmap for marker genes
merged_combined.maker %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top1

identical(as.character(top1$cluster), cell_types)


FigureS5c <- DotPlot(combined_Myeloid, features = unique(top1$gene)) + 
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


pdf(paste0('FigureS5/FigureS5c.pdf'), width = 15, height = 5)
print(FigureS5c)
dev.off()


##------------------------------------------
## Figure S5D: Tsne of B cells
##------------------------------------------
Idents(combined_nontumor) <- combined_nontumor$CT_L1
combined_B <- subset(combined_nontumor, idents = 'B')

Idents(combined_B) <- combined_B$CT_L2
combined_B <- subset(combined_B, idents = c('B intermediate', 'B memory', 'B naive',
                                            'Plasmablast'))

combined_B <- RunTSNE(combined_B, reduction = 'umap.pca', dims.use = 1:10, perplexity = 300)

ptS5d <- DimPlot(combined_B, reduction = 'tsne', group.by = 'CT_L2', 
                 # split.by = 'Status_new', 
                 label = T, repel = T, label.size = 5) + 
  NoAxes() + ggtitle(' ') + 
  scale_x_reverse() +
  theme_void() + 
  guides(color = guide_legend(override.aes = list(shape = 15, size = 6),
                              label.hjust = 0,
                              label.theme = element_text())) + #NoLegend() +
  theme(strip.text.x = element_blank()) + 
  theme(legend.position = 'bottom')

pdf('FigureS5/FigureS5d.pdf', height = 5, width = 5)
print(ptS5d)
dev.off()

##------------------------------------------
## Figure S5E: B cell fractions
##------------------------------------------

# Create the dataframe
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



cells_2 <- c('B intermediate', 'B memory', 'B naive',
             'CD16 Mono', 'Macrophage',
             'CD4 CTL', 'CD4 Naive', 'CD4 Proliferating', 'CD4 TCM', 'CD4 TEM',
             'CD8 Naive', 'CD8 Proliferating', 'CD8 TCM', 'CD8 TEM',
             'Treg', 'gdT', 'dnT', 'MAIT', 
             'NK', 'NK Proliferating', 'NK_CD56bright', 'ILC',
             'CAFs', 'PVL', 'cDC1', 'cDC2', 'pDC', 
             'Endothelial', 'Eryth', 'HSPC', 'Plasmablast', 'Platelet')

df_nontumor_L2_frac <- df_nontumor_L2
for(i in unique(df_nontumor_L2$Sample)){
  
  Indx <- which(df_nontumor_L2$Sample == i)
  
  df_nontumor_L2_frac[Indx, 'Count'] = df_nontumor_L2[Indx, 'Count'] / sum(df_nontumor_L2[Indx, 'Count'])
  
}

df_nontumor_L2_frac_sub <- df_nontumor_L2_frac[df_nontumor_L2_frac$CellTypes %in% cells_2, ]

category <- rep(NA, nrow(df_nontumor_L2_frac_sub))
category[df_nontumor_L2_frac_sub$CellTypes %in% c('B naive', 'CD4 Naive', 'CD8 Naive')] <- 'Naive cells'
category[df_nontumor_L2_frac_sub$CellTypes %in% c('CD4 Proliferating', 'CD8 Proliferating')] <- 'Proliferating cells'
category[df_nontumor_L2_frac_sub$CellTypes %in% c('B intermediate', 'B memory')] <- 'B cells'
category[df_nontumor_L2_frac_sub$CellTypes %in% c('CD16 Mono', 'Macrophage')] <- 'Macrophage'
category[df_nontumor_L2_frac_sub$CellTypes %in% c('CD4 CTL', 'CD4 TCM', 'CD4 TEM')] <- 'CD4'
category[df_nontumor_L2_frac_sub$CellTypes %in% c('CD8 TCM', 'CD8 TEM')] <- 'CD8'
category[df_nontumor_L2_frac_sub$CellTypes %in% c('Treg')] <- 'Treg'
category[df_nontumor_L2_frac_sub$CellTypes %in% c('gdT', 'dnT', 'MAIT')] <- 'other T'
category[df_nontumor_L2_frac_sub$CellTypes %in% c('NK', 'NK Proliferating', 'NK_CD56bright', 'ILC')] <- 'NK'
category[df_nontumor_L2_frac_sub$CellTypes %in% c('CAFs', 'PVL')] <- 'CAFs'
category[df_nontumor_L2_frac_sub$CellTypes %in% c('cDC1', 'cDC2', 'pDC')] <- 'DC'
category[df_nontumor_L2_frac_sub$CellTypes %in% c('Endothelial', 'Eryth', 'HSPC', 'Plasmablast', 'Platelet')] <- 'Endothelial'

df_nontumor_L2_frac_sub$category <- as.factor(category)



df_B_cells_L2_frac_sub <- df_nontumor_L2_frac_sub[df_nontumor_L2_frac_sub$CellTypes %in% c('B intermediate', 'B memory', 'B naive', 'Plasmablast'), ]
df_B_cells_L2_frac_sub$CellTypes <- factor(as.character(df_B_cells_L2_frac_sub$CellTypes))

ptS5e <- ggplot(df_B_cells_L2_frac_sub) + 
  geom_boxplot(size = 0.3, outlier.shape = NA, width = 0.618, lwd = 0.5,
               aes(y = Count, x = CellTypes, color = Status)) +
  geom_jitter(shape = 16, position = position_jitterdodge(), 
              size = 2,  alpha = 0.9, aes(fill = Status, y = Count, x = CellTypes, color = Status)) +
  theme_classic() + 
  ggtitle(" ") +
  xlab(" ") + ylab("Cell Fraction")  +
  theme(legend.position="top",
        # legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        panel.border = element_blank(),
        strip.background = element_blank(),
  ) + guides(color=guide_legend(title = ' '),
             fill = FALSE)   +
  scale_fill_manual(values=c('#cccc66', '#66cccc', '#ff9966'))  +
  scale_color_manual(values=c('#cccc66', '#66cccc', '#ff9966')) #+
# stat_pvalue_manual(
#   stat.test,  label = "p.signif", tip.length = 0.01, hide.ns = T,
# )

pdf('FigureS5/FigureS5e.pdf', width = 5, height = 4)
print(ptS5e)
dev.off()


##------------------------------------------
## Figure S5F: violin plot of CSF1R and CCL2
##------------------------------------------

my_comparisons<- list( c("Baseline", "Early Progressor"), 
                       c("Baseline", "Late Progressor"),
                       c("Early Progressor", "Late Progressor"))
Idents(combined_Myeloid) <- combined_Myeloid$Status_new

Figure_CSF1R <- VlnPlot(combined_Myeloid, features = 'CSF1R', pt.size = 0) +
  geom_boxplot(color = "black", outlier.shape = NA, width = 0.2) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             label="p.signif", method="wilcox.test") + ylim(c(0, 6.5)) + NoLegend() +
  xlab(' ') + ggtitle('Expression level of CSF1R') + ylab(' ') +
  scale_fill_manual(values = c('#cccc66', '#66cccc', '#ff9966'))  +  
  scale_x_discrete(labels = c('LP', 'EP', 'BL'),
                   breaks = c('Late Progressor', 'Early Progressor', 'Baseline')) 
pdf(paste0('FigureS5/FigureS5f1.pdf'), width = 3.5, height = 4)
print(Figure_CSF1R)
dev.off()

Figure_CCL2 <- VlnPlot(combined_Myeloid, features = 'CCL2', pt.size = 0) +
  geom_boxplot(color = "black", outlier.shape = NA, width = 0.2) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             label="p.signif", method="wilcox.test") + ylim(c(0, 6.5)) + NoLegend() +
  xlab(' ') + ggtitle('Expression level of CCL2') + ylab(' ') +
  scale_fill_manual(values = c('#cccc66', '#66cccc', '#ff9966'))  +  
  scale_x_discrete(labels = c('LP', 'EP', 'BL'),
                   breaks = c('Late Progressor', 'Early Progressor', 'Baseline')) 
pdf(paste0('FigureS5/FigureS5f2.pdf'), width = 3.5, height = 4)
print(Figure_CCL2)
dev.off()

