##-------------------------------
## Figure Panel S4
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

library(pheatmap)
library(RColorBrewer)
library(wesanderson)
library(viridis)
library(ggsci)
library(tidytext)
library(ggpubr)
library(cowplot)
library(facetscales)
library(latex2exp)
library(ggstatsplot)
library(scales)

load('/Data/combined_Tcells.RData')


##-------------------------------
## Figure S4A: Heatmap of metastatistic
## for CD4T cells
##-------------------------------

CD4T_cells <- c('CD4 CTL', 'CD4 Naive', 'CD4 Proliferating', 'CD4 TCM', 'CD4 TEM',
                'Treg', 'gdT')
Idents(combined_Tcells) <- combined_Tcells$CT_L2
combined_CD4Tcells <- subset(combined_Tcells, idents = CD4T_cells)

Idents(combined_CD4Tcells) <- combined_CD4Tcells$Status_new
CD4T <- readxl::read_excel('/Users/yangpeng/Box Sync/CDK4/Tcell_function/41591_2023_2371_MOESM3_ESM.xlsx', sheet = 6, col_names = T)
df_CD4T <- data.frame(cluster = CD4T$`Supplementary Table 5. Top 50 DEGs for 12 CD4 T cell clusters`[-1],
                      gene    = CD4T$...2[-1])



FigureS4b <- DotPlot(combined_CD4Tcells, features =  df_CD4T$gene[df_CD4T$cluster == 'CD4_c4_Tstr']) + 
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


pdf(file = 'FigureS4/FigureS4a.pdf', width = 12, height = 4)
print(FigureS4a)
dev.off()


##-------------------------------
## Figure S4b: Adding module score
## for CD4T cells
##-------------------------------
##--------------------------
## Curated gene list
##--------------------------
MP_meta <- readxl::read_excel('/Tcell_function/41591_2023_2371_MOESM3_ESM.xlsx', sheet = 7, col_names = T, skip = 1)
MP_meta <- as.data.frame(MP_meta)
org_MP_Names <- names(MP_meta)
# MP_meta <- list(stress = MP_meta$`Stress response`)
names(MP_meta) = paste0('Curated', 1:ncol(MP_meta))

combined_CD4Tcells <- AddModuleScore(
  object = combined_CD4Tcells,
  features = MP_meta,
  ctrl = 5,
  name = 'Curated'
)

colorForStatus <- c('#cccc66', '#66cccc', '#ff9966')
my_comparisons <- list( c("Baseline", "Early Progressor"), 
                        c("Baseline", "Late Progressor"), 
                        c("Early Progressor", "Late Progressor") )
pdf('FigureS4/FigureS4b.pdf', width = 5, height = 5)
VlnPlot(combined_CD4Tcells, features = paste0('Curated', which(org_MP_Names == 'Stress response')), pt.size = 0.001) + 
  stat_compare_means(comparisons = my_comparisons) + ylim(c(-0.1, 0.35)) +
  scale_fill_manual(
    values = colorForStatus) + NoLegend() + ggtitle(' ') + xlab(' ') + 
  stat_summary(fun.y = median, geom='point', size = 20, colour = "white", shape = 95)
dev.off()



##-------------------------------
## Figure S4C: Heatmap of metastatistic
## for CD8T cells
##-------------------------------

CD8T_cells <- c('CD8 Naive', 'CD8 Proliferating', 'CD8 TCM', 'CD8 TEM', 
                'Treg', 'gdT')
Idents(combined_Tcells) <- combined_Tcells$CT_L2
combined_CD8Tcells <- subset(combined_Tcells, idents = CD8T_cells)

Idents(combined_CD8Tcells) <- combined_CD8Tcells$Status_new
# if(file.exists('FigureS4/Data/CD8T_cell_markers.RData')){
#   load('FigureS4/Data/CD8T_cell_markers.RData')
# }else{
#   combined_CD8Tcells.maker <- FindAllMarkers(combined_CD8Tcells, only.pos = TRUE, 
#                                              min.pct = 0.10, logfc.threshold = 0.25)
#   save(combined_CD8Tcells.maker, file = 'FigureS4/Data/CD8T_cell_markers.RData')
#   
# }

CD8T <- readxl::read_excel('/Users/yangpeng/Box Sync/CDK4/Tcell_function/41591_2023_2371_MOESM3_ESM.xlsx', sheet = 4, col_names = T)
df_CD8T <- data.frame(cluster = CD8T$`Supplementary Table 3. Top 50 differentially expressed genes (DEGs) for 14 CD8 T cell clusters`[-1],
                      gene    = CD8T$...2[-1])
# df_CD4T$gene[df_CD4T$cluster == 'CD4_c4_Tstr']
Idents(combined_CD8Tcells) <- combined_CD8Tcells$Status_new


FigureS4c <- DotPlot(combined_CD8Tcells, features = df_CD8T$gene[df_CD8T$cluster == 'CD8_c4_Tstr']) + 
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


pdf(file = 'FigureS4/FigureS4c.pdf', width = 12, height = 4)
print(FigureS4c)
dev.off()






##-------------------------------
## Figure S4D: Adding module score
## for CD8T cells
##-------------------------------
##--------------------------
## Curated gene list
##--------------------------
MP_meta <- readxl::read_excel('/Tcell_function/41591_2023_2371_MOESM3_ESM.xlsx', sheet = 5, col_names = T, skip = 1)
# MP_meta <- readxl::read_excel('/Users/py11/Box Sync/CDK4/Tcell_function/41591_2023_2371_MOESM3_ESM.xlsx', sheet = 5, col_names = T, skip = 1)
MP_meta <- as.data.frame(MP_meta)
org_MP_Names <- names(MP_meta)
# MP_meta <- list(stress = MP_meta$`Stress response`)
# MP_meta <- list(MP_meta)
names(MP_meta) = paste0('Curated', 1:ncol(MP_meta))


combined_CD8Tcells <- AddModuleScore(
  object = combined_CD8Tcells,
  features = MP_meta,
  ctrl = 5,
  name = 'Curated',
)

my_comparisons <- list( c("Baseline", "Early Progressor"), 
                        c("Baseline", "Late Progressor"), 
                        c("Early Progressor", "Late Progressor") )

Idents(combined_CD8Tcells) <- combined_CD8Tcells$Status_new
pdf('FigureS4/FigureS4d.pdf', width = 5, height = 5)
VlnPlot(combined_CD8Tcells, features = paste0('Curated', which(org_MP_Names == 'Stress response')), pt.size = 0.001) + 
  stat_compare_means(comparisons = my_comparisons) + ylim(c(-0.2, 0.5)) +
  scale_fill_manual(
    values = colorForStatus) + NoLegend() + ggtitle(' ') + xlab(' ') + 
  stat_summary(fun.y = median, geom='point', size = 20, colour = "white", shape = 95)
dev.off()


##-------------------------------------
## Figure S4F: Cell fractions stratified
## by metastatistic sites
##-------------------------------------

umapColor <- c("#b3cde0", "#C8CADF", "#9ACBDE", "#6DCCDD","#6497b1", "#EDC6DD","#F092B1","#F27FA5",
               "#F47892", "#F6A395","#F8AD77","#E7B080","#C9AFA2","#ABADC5", "#C4DA5D")

df_metastatic_L2 <- c()

for(sample in unique(combined_Tcells$orig.ident)){
  df_metastatic_L2 <- rbind(df_metastatic_L2, 
                            cbind(cbind(cbind(table(combined_Tcells$CT_L2[combined_Tcells$orig.ident == sample]), sample),
                                        as.character(combined_Tcells$Status_new[combined_Tcells$orig.ident == sample][1])),
                                  as.character(combined_Tcells$metastatic[combined_Tcells$orig.ident == sample][1])) )
}

df_metastatic_L2 <- cbind(df_metastatic_L2, rownames(df_metastatic_L2))

colnames(df_metastatic_L2) <- c('Count', 'Sample', 'Status', 'Metastatic', 'CellTypes')
df_metastatic_L2 <- as.data.frame(df_metastatic_L2)

df_metastatic_L2 <- df_metastatic_L2[df_metastatic_L2$Count != 0, ]
## Make factors

df_metastatic_L2$Count <- as.numeric(df_metastatic_L2$Count)
df_metastatic_L2$Sample <- as.factor(df_metastatic_L2$Sample)
df_metastatic_L2$Status <- as.factor(df_metastatic_L2$Status)
df_metastatic_L2$Metastatic <- as.factor(df_metastatic_L2$Metastatic)
df_metastatic_L2$CellTypes <- factor(df_metastatic_L2$CellTypes,
                                     levels = T_cells)


ptS4f <- ggplot(df_metastatic_L2, aes(fill = CellTypes, y = Count, x = Sample)) +
  geom_bar(position="fill", stat="identity") +
  facet_grid(~ Metastatic, drop = T, scales = 'free', space = 'free') +
  theme_classic() +
  scale_fill_manual(values = umapColor[1:13]) +
  xlab(" ") + ylab(" ")  +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
    strip.background = element_blank(),
    strip.text.x = element_blank()
    # strip.text.x = element_text(size = 10)
  ) + NoLegend()

pdf(paste0('FigureS4/FigureS4f.pdf'), width = 7, height = 3)
print(ptS4f)
dev.off()



##-------------------------------------
## Figure S4g: Heatmap of DE genes 
## Across different metastatic site
## for baseline samples
##-------------------------------------


Idents(combined_Tcells) <- combined_Tcells$Status_new
combined_BL <- subset(combined_Tcells, idents = 'Baseline')
Idents(combined_BL) <- combined_BL$metastatic


if(file.exists('FigureS4/Data/cell_markers_BL.RData')){
  load('FigureS4/Data/cell_markers_BL.RData')
}else{
  merged_combined.maker <- FindAllMarkers(combined_BL, only.pos = TRUE, 
                                          min.pct = 0.10, logfc.threshold = 0.25)
  save(merged_combined.maker, file = 'FigureS4/Data/cell_markers_BL.RData')
  
}

merged_combined.maker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

## Make the heatmap for marker genes
merged_combined.maker %>%
  group_by(cluster) %>%
  top_n(n = 15, wt = avg_log2FC) -> top1


FigureS4g <- DotPlot(combined_BL, features = unique(top1$gene)) + 
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


pdf(paste0('FigureS4/FigureS4g.pdf'), width = 8, height = 3)
print(FigureS4g)
dev.off()


##-------------------------------------
## Figure S4H: Heatmap of DE genes 
## Across different metastatic site
## for Late Progressor samples
##-------------------------------------


Idents(combined_Tcells) <- combined_Tcells$Status_new
combined_LP <- subset(combined_Tcells, idents = 'Late Progressor')
Idents(combined_LP) <- combined_LP$metastatic
combined_LP <- subset(combined_LP, idents = c('Liver', 'Pleural effusion'))

if(file.exists('FigureS4/Data/cell_markers_LP.RData')){
  load('FigureS4/Data/cell_markers_LP.RData')
}else{
  merged_combined.maker <- FindAllMarkers(combined_LP, only.pos = TRUE, 
                                          min.pct = 0.10, logfc.threshold = 0.25)
  save(merged_combined.maker, file = 'FigureS4/Data/cell_markers_LP.RData')
  
}

merged_combined.maker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

## Make the heatmap for marker genes
merged_combined.maker %>%
  group_by(cluster) %>%
  top_n(n = 15, wt = avg_log2FC) -> top1


FigureS5b <- DotPlot(combined_LP, features = unique(top1$gene)) + 
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


pdf(paste0('FigureS5/FigureS5b.pdf'), width = 8, height = 3)
print(FigureS4H)
dev.off()
