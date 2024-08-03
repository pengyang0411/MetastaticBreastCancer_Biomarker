##-------------------------------
## Figure Panel 4
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

## cell type colors

umapColor <- c("#b3cde0", "#C8CADF", "#9ACBDE", "#6DCCDD","#6497b1", "#EDC6DD","#F092B1","#F27FA5",
               "#F47892", "#F6A395","#F8AD77","#E7B080","#C9AFA2","#ABADC5", "#C4DA5D")


##-------------------------------
## Figure 4A: Umap plot across 
## different sample status
##-------------------------------
combined_Tcells$CT_L2 <- factor(as.character(combined_Tcells$CT_L2),
                                levels = T_cells)

Idents(combined_Tcells) <- combined_Tcells$CT_L2
combined_Tcells$CT_L2.num <- as.factor(as.numeric(as.factor(combined_Tcells$CT_L2)))


pt4a <- DimPlot(combined_Tcells, reduction = 'tsne', group.by = 'CT_L2.num', 
                # split.by = 'metastatic',
                label = F, repel = T, label.size = 5) + 
  NoAxes() + ggtitle(' ') + 
  scale_x_reverse() +
  theme_void() + 
  scale_color_manual(# breaks = unique(merged_combined$orig.ident.num),
    labels = paste(1:length(unique(combined_Tcells$CT_L2)), T_cells),
    values = umapColor[1:13]) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 6),
                              label.hjust = 0,
                              label.theme = element_text())) + NoLegend() #+
  # theme(strip.text.x = element_blank())

pt4a

pdf('Figure4/Figure4a.pdf', width = 5, height = 5)
print(pt4a)
dev.off()



##-------------------------------
## Figure 4B: Marker gene plot
##-------------------------------

Idents(combined_Tcells) <- combined_Tcells$CT_L2

if(file.exists('Figure4/Data/cell_markers.RData')){
  load('Figure4/Data/cell_markers.RData')
}else{
  merged_combined.maker <- FindAllMarkers(combined_Tcells, only.pos = TRUE, 
                                          min.pct = 0.10, logfc.threshold = 0.25)
  save(merged_combined.maker, file = 'Figure4/Data/cell_markers.RData')
  
}

merged_combined.maker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

## Make the heatmap for marker genes
merged_combined.maker %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC) -> top1

identical(as.character(top1$cluster), cell_types)


Figure4b <- DotPlot(combined_Tcells, features = unique(top1$gene)) + 
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


pdf(paste0('Figure4/Figure4b.pdf'), width = 8, height = 5)
print(Figure4b)
dev.off()

##-------------------------------
## Figure 4C: T cell distribution
##-------------------------------


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
df_metastatic_L2$CellTypes <- factor(df_metastatic_L2$CellTypes,
                                     levels = T_cells)


pt4c <- ggplot(df_metastatic_L2, aes(fill = CellTypes, y = Count, x = Sample)) +
  geom_bar(position="fill", stat="identity") +
  facet_grid(~ Status, drop = T, scales = 'free', space = 'free') +
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

pdf(paste0('Figure4/Figure4c.pdf'), width = 7, height = 3)
print(pt4c)
dev.off()


##-------------------------------
## Figure 4D: Heatmap plot of 
## DE genes for metatastics from
## CD8 T cells
##-------------------------------

## Create the dataframe
df_nontumor_L2 <- c()
for(sample in unique(combined_Tcells$orig.ident)){
  df_nontumor_L2 <- rbind(df_nontumor_L2, 
                          cbind(cbind(table(combined_Tcells$CT_L2[combined_Tcells$orig.ident == sample]), sample),
                                as.character(combined_Tcells$Status_new[combined_Tcells$orig.ident == sample][1])))
}

df_nontumor_L2 <- cbind(df_nontumor_L2, rownames(df_nontumor_L2))

colnames(df_nontumor_L2) <- c('Count', 'Sample', 'Status', 'CellTypes')
df_nontumor_L2 <- as.data.frame(df_nontumor_L2)

## Make factors
df_nontumor_L2$Count <- as.numeric(df_nontumor_L2$Count)
df_nontumor_L2$Sample <- as.factor(df_nontumor_L2$Sample)
df_nontumor_L2$Status <- as.factor(df_nontumor_L2$Status)
df_nontumor_L2$CellTypes <- as.factor(df_nontumor_L2$CellTypes)

cells_2 <- c('CD4 CTL', 'CD4 TCM', 'CD4 TEM',
             'CD8 TCM', 'CD8 TEM')

df_nontumor_L2_frac <- df_nontumor_L2[df_nontumor_L2$CellTypes %in% cells_2, ]

for(i in unique(df_nontumor_L2$Sample)){
  
  Indx <- which(df_nontumor_L2$Sample == i)
  
  df_nontumor_L2_frac[Indx, 'Count'] = df_nontumor_L2[Indx, 'Count'] / sum(df_nontumor_L2[Indx, 'Count'])
  
}

df_nontumor_L2_frac_sub <- df_nontumor_L2_frac[df_nontumor_L2_frac$CellTypes %in% cells_2, ]

category <- rep(NA, nrow(df_nontumor_L2_frac_sub))
category[df_nontumor_L2_frac_sub$CellTypes %in% c('CD4 CTL', 'CD4 TCM', 'CD4 TEM')] <- 'CD4+T'
category[df_nontumor_L2_frac_sub$CellTypes %in% c('CD8 TCM', 'CD8 TEM')] <- 'CD8+T'

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
Figure4d <- ggplot(df_nontumor_L2_frac_sub) + 
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

pdf(paste0('Figure4/Figure4d.pdf'), width = 8, height = 3)
print(Figure4d)
dev.off()


##-------------------------------
## Figure 4E: Functional plot
## from CD8 T cells
##-------------------------------
CD8T <- readxl::read_excel('/Tcell_function/41591_2023_2371_MOESM3_ESM.xlsx', sheet = 5, col_names = T)
marker.list <- list()

for(i in 1:length(CD8T)){
  
  FunctionName <- CD8T[[i]][1]
  
  marker.list[[FunctionName]] <- CD8T[[i]][2:c(which.max(is.na(CD8T[[i]])) - 1)]
  
}

combined_CD8Tcells <-  AddModuleScore(combined_CD8Tcells,
                                      features = marker.list,
                                      ctrl = 5,
                                      name = "FunctionScore")
## Adding to meta data
for(i in 1:length(marker.list)){
  colnames(combined_CD8Tcells@meta.data)[colnames(combined_CD8Tcells@meta.data) == paste0("FunctionScore", i)] <- names(marker.list)[i]
}
Idents(combined_CD8Tcells) <- combined_CD8Tcells$CT_L2
Differentiation <- c("Naïve", "Activation:Effector function", "Exhaustion")
Function <- c("TCR Signaling", "Cytotoxicity", "Cytokine/Cytokine receptor",
              "Chemokine/Chemokine receptor", "Senescence", "Anergy",
              "NFKB Signaling", "Stress response", "MAPK Signaling", "Adhesion",
              "IFN Response")
Metabolism <- c("Oxidative phosphorylation", "Glycolysis", "Fatty acid metabolism")
Apoptosis <- c("Pro-apoptosis", "Anti-apoptosis")
MarkerNameVector <- c(Differentiation, Function, Metabolism, Apoptosis)
FunctionScoreMatrix <- matrix(0,
                              ncol = length(unique(combined_CD8Tcells$CT_L2)),
                              nrow = length(marker.list))
colnames(FunctionScoreMatrix) <- unique(combined_CD8Tcells$CT_L2)
rownames(FunctionScoreMatrix) <- MarkerNameVector
for(ci in 1:ncol(FunctionScoreMatrix)){
  for(ri in 1:nrow(FunctionScoreMatrix)){
    FunctionVec <- as_tibble(combined_CD8Tcells@meta.data) %>% dplyr::pull(MarkerNameVector[ri])
    fv <- mean(FunctionVec[combined_CD8Tcells$CT_L2 == colnames(FunctionScoreMatrix)[ci]])
    FunctionScoreMatrix[ri, ci] <- fv
  }
}
FunctionScoreMatrix <- t(apply(FunctionScoreMatrix, 1, rescale, to=c(-1, 1)))
# orderC = c("CD8_c13", "CD8_c3", "CD8_c6", "CD8_c0", "CD8_c11", "CD8_c9", "CD8_c10", "CD8_c12", "CD8_c8", "CD8_c2", "CD8_c7", "CD8_c4", "CD8_c5", "CD8_c1")
FunctionScoreMatrix <- FunctionScoreMatrix[,CD8T_cells]
my.breaks <- c(seq(-1, 0, by=0.1), seq(0.1, 1, by=0.1))
my.colors <- c(
  colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks)/2),
  colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks)/2))
## cellType_col <- data.frame(cell.type = CD8_Obj_CellType)
## rownames(cellType_col) <- colnames(FunctionScoreMatrix)
signatureType_row <- data.frame(Signature.type = c(
  rep("Differentiation", length(Differentiation)),
  rep("Function", length(Function)),
  rep("Metabolism", length(Metabolism)),
  rep("Apoptosis", length(Apoptosis))))
rownames(signatureType_row) <- MarkerNameVector

figurePath <- 'Figure4/'
pheatmap(FunctionScoreMatrix,
         show_colnames = T,
         show_rownames = T,
         ## annotation_col = cellType_col,
         annotation_row = signatureType_row,
         gaps_row = c(3, 14, 17),
         cluster_rows = F,
         cluster_cols = F,
         breaks = my.breaks,
         color = my.colors,
         border_color = "NA",
         fontsize = 8,
         width = 5,
         height = 3.8,
         filename = file.path(figurePath, paste0("Figure4e.pdf")))


##-------------------------------
## Figure 4F Boxplot for HSP90AA
## of CD4+T cells
##-------------------------------

my_comparisons<- list( c("Baseline", "Early Progressor"), 
                       c("Baseline", "Late Progressor"),
                       c("Early Progressor", "Late Progressor"))
Idents(combined_CD4Tcells) <- combined_CD4Tcells$Status_new
## Making the violin plot for CD8T cells with UBC and UBB
Figure4f <- VlnPlot(combined_CD4Tcells, features = 'HSP90AA1', pt.size = 0) +
  geom_boxplot(color = "black", outlier.shape = NA, width = 0.2) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             label="p.signif", method="wilcox.test") + ylim(c(0, 6.5)) + NoLegend() +
  xlab(' ') + ggtitle('Expression level of HSP90AA1') + ylab(' ') +
  scale_fill_manual(values = c('#cccc66', '#66cccc', '#ff9966'))  +  
  scale_x_discrete(labels = c('LP', 'EP', 'BL'),
                 breaks = c('Late Progressor', 'Early Progressor', 'Baseline')) 
pdf(paste0('Figure4/Figure4f.pdf'), width = 3.5, height = 4)
print(Figure4f)
dev.off()



##-------------------------------
## Figure 4G Boxplot for HSP90AB1
## of CD4+T cells
##-------------------------------

my_comparisons<- list( c("Baseline", "Early Progressor"), 
                       c("Baseline", "Late Progressor"),
                       c("Early Progressor", "Late Progressor"))
Idents(combined_CD4Tcells) <- combined_CD4Tcells$Status_new
## Making the violin plot for CD8T cells with UBC and UBB
Figure4g <- VlnPlot(combined_CD4Tcells, features = 'HSP90AB1', pt.size = 0) +
  geom_boxplot(color = "black", outlier.shape = NA, width = 0.2) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             label="p.signif", method="wilcox.test") + ylim(c(0, 6.5)) + NoLegend() +
  xlab(' ') + ggtitle('Expression level of HSP90AB1') + ylab(' ') +
  scale_fill_manual(values = c('#cccc66', '#66cccc', '#ff9966'))  +  
  scale_x_discrete(labels = c('LP', 'EP', 'BL'),
                   breaks = c('Late Progressor', 'Early Progressor', 'Baseline')) 
pdf(paste0('Figure4/Figure4g.pdf'), width = 3.5, height = 4)
print(Figure4g)
dev.off()

##-------------------------------
## Figure 4H Boxplot for HSPA8
## of CD4+T cells
##-------------------------------

my_comparisons<- list( c("Baseline", "Early Progressor"), 
                       c("Baseline", "Late Progressor"),
                       c("Early Progressor", "Late Progressor"))
Idents(combined_CD4Tcells) <- combined_CD4Tcells$Status_new
## Making the violin plot for CD8T cells with UBC and UBB
Figure4h <- VlnPlot(combined_CD4Tcells, features = 'HSPA8', pt.size = 0) +
  geom_boxplot(color = "black", outlier.shape = NA, width = 0.2) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             label="p.signif", method="wilcox.test") + ylim(c(0, 6.5)) + NoLegend() +
  xlab(' ') + ggtitle('Expression level of HSPA8') + ylab(' ') +
  scale_fill_manual(values = c('#cccc66', '#66cccc', '#ff9966'))  +  
  scale_x_discrete(labels = c('LP', 'EP', 'BL'),
                   breaks = c('Late Progressor', 'Early Progressor', 'Baseline')) 
pdf(paste0('Figure4/Figure4h.pdf'), width = 3.5, height = 4)
print(Figure4h)
dev.off()

##-------------------------------
## Figure 4I: Functional plot
## from CD4 T cells
##-------------------------------


CD4T <- readxl::read_excel('/Tcell_function/41591_2023_2371_MOESM3_ESM.xlsx', sheet = 7, col_names = T)
marker.list <- list()

for(i in 1:length(CD4T)){
  
  FunctionName <- CD4T[[i]][1]
  
  marker.list[[FunctionName]] <- CD4T[[i]][2:c(which.max(is.na(CD4T[[i]])) - 1)]
  
}

combined_CD4Tcells <-  AddModuleScore(combined_CD4Tcells,
                                      features = marker.list,
                                      ctrl = 5,
                                      name = "FunctionScore")
## Adding to meta data
for(i in 1:length(marker.list)){
  colnames(combined_CD4Tcells@meta.data)[colnames(combined_CD4Tcells@meta.data) == paste0("FunctionScore", i)] <- names(marker.list)[i]
}
Idents(combined_CD4Tcells) <- combined_CD4Tcells$CT_L2
Differentiation <- c("Naïve", "Activation/Effector function", "Exhaustion")
Function <- c("TCR signaling", "Cytotoxicity", "Cytokine/Cytokine receptor",
              "Chemokine/Chemokine receptor", "Stress response", "Adhesion",
              "IFN response", "Treg signature", "Costimulatory molecules")
Metabolism <- c("OXPHOS", "Glycolysis", "Lipid metabolism")
Apoptosis <- c("Pro-apoptosis", "Anti-apoptosis")
MarkerNameVector <- c(Differentiation, Function, Metabolism, Apoptosis)
FunctionScoreMatrix <- matrix(0,
                              ncol = length(unique(combined_CD4Tcells$CT_L2)),
                              nrow = length(marker.list))
colnames(FunctionScoreMatrix) <- unique(combined_CD4Tcells$CT_L2)
rownames(FunctionScoreMatrix) <- MarkerNameVector
for(ci in 1:ncol(FunctionScoreMatrix)){
  for(ri in 1:nrow(FunctionScoreMatrix)){
    FunctionVec <- as_tibble(combined_CD4Tcells@meta.data) %>% dplyr::pull(MarkerNameVector[ri])
    fv <- mean(FunctionVec[combined_CD4Tcells$CT_L2 == colnames(FunctionScoreMatrix)[ci]])
    FunctionScoreMatrix[ri, ci] <- fv
  }
}
FunctionScoreMatrix <- t(apply(FunctionScoreMatrix, 1, rescale, to=c(-1, 1)))
# orderC = c("CD8_c13", "CD8_c3", "CD8_c6", "CD8_c0", "CD8_c11", "CD8_c9", "CD8_c10", "CD8_c12", "CD8_c8", "CD8_c2", "CD8_c7", "CD8_c4", "CD8_c5", "CD8_c1")
FunctionScoreMatrix <- FunctionScoreMatrix[,CD4T_cells]
my.breaks <- c(seq(-1, 0, by=0.1), seq(0.1, 1, by=0.1))
my.colors <- c(
  colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks)/2),
  colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks)/2))
## cellType_col <- data.frame(cell.type = CD8_Obj_CellType)
## rownames(cellType_col) <- colnames(FunctionScoreMatrix)
signatureType_row <- data.frame(Signature.type = c(
  rep("Differentiation", length(Differentiation)),
  rep("Function", length(Function)),
  rep("Metabolism", length(Metabolism)),
  rep("Apoptosis", length(Apoptosis))))
rownames(signatureType_row) <- MarkerNameVector

figurePath <- 'Figure4/'
pheatmap(FunctionScoreMatrix,
         show_colnames = T,
         show_rownames = T,
         ## annotation_col = cellType_col,
         annotation_row = signatureType_row,
         gaps_row = c(3, 12, 15),
         cluster_rows = F,
         cluster_cols = F,
         breaks = my.breaks,
         color = my.colors,
         border_color = "NA",
         fontsize = 8,
         width = 5,
         height = 3.8,
         filename = file.path(figurePath, paste0("Figure4i.pdf")))

##-------------------------------
## Figure 4J Boxplot for HSP90AA
## of CD8+T cells
##-------------------------------

my_comparisons<- list( c("Baseline", "Early Progressor"), 
                       c("Baseline", "Late Progressor"),
                       c("Early Progressor", "Late Progressor"))
Idents(combined_CD8Tcells) <- combined_CD8Tcells$Status_new
## Making the violin plot for CD8T cells with UBC and UBB
Figure4j <- VlnPlot(combined_CD8Tcells, features = 'HSP90AA1', pt.size = 0) +
  geom_boxplot(color = "black", outlier.shape = NA, width = 0.2) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             label="p.signif", method="wilcox.test") + ylim(c(0, 6.5)) + NoLegend() +
  xlab(' ') + ggtitle('Expression level of HSP90AA1') + ylab(' ') +
  scale_fill_manual(values = c('#cccc66', '#66cccc', '#ff9966'))  +  
  scale_x_discrete(labels = c('LP', 'EP', 'BL'),
                   breaks = c('Late Progressor', 'Early Progressor', 'Baseline')) 
pdf(paste0('Figure4/Figure4j.pdf'), width = 3.5, height = 4)
print(Figure4j)
dev.off()



##-------------------------------
## Figure 4L Boxplot for HSP90AB1
## of CD8+T cells
##-------------------------------

my_comparisons<- list( c("Baseline", "Early Progressor"), 
                       c("Baseline", "Late Progressor"),
                       c("Early Progressor", "Late Progressor"))
Idents(combined_CD8Tcells) <- combined_CD8Tcells$Status_new
## Making the violin plot for CD8T cells with UBC and UBB
Figure4k <- VlnPlot(combined_CD8Tcells, features = 'HSP90AB1', pt.size = 0) +
  geom_boxplot(color = "black", outlier.shape = NA, width = 0.2) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             label="p.signif", method="wilcox.test") + ylim(c(0, 6.5)) + NoLegend() +
  xlab(' ') + ggtitle('Expression level of HSP90AB1') + ylab(' ') +
  scale_fill_manual(values = c('#cccc66', '#66cccc', '#ff9966'))  +  
  scale_x_discrete(labels = c('LP', 'EP', 'BL'),
                   breaks = c('Late Progressor', 'Early Progressor', 'Baseline')) 
pdf(paste0('Figure4/Figure4k.pdf'), width = 3.5, height = 4)
print(Figure4k)
dev.off()

##-------------------------------
## Figure 4L Boxplot for HSPA8
## of CD8+T cells
##-------------------------------

my_comparisons<- list( c("Baseline", "Early Progressor"), 
                       c("Baseline", "Late Progressor"),
                       c("Early Progressor", "Late Progressor"))
Idents(combined_CD8Tcells) <- combined_CD8Tcells$Status_new
## Making the violin plot for CD8T cells with UBC and UBB
Figure4l <- VlnPlot(combined_CD8Tcells, features = 'HSPA8', pt.size = 0) +
  geom_boxplot(color = "black", outlier.shape = NA, width = 0.2) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             label="p.signif", method="wilcox.test") + ylim(c(0, 6.5)) + NoLegend() +
  xlab(' ') + ggtitle('Expression level of HSPA8') + ylab(' ') +
  scale_fill_manual(values = c('#cccc66', '#66cccc', '#ff9966'))  +  
  scale_x_discrete(labels = c('LP', 'EP', 'BL'),
                   breaks = c('Late Progressor', 'Early Progressor', 'Baseline')) 
pdf(paste0('Figure4/Figure4l.pdf'), width = 3.5, height = 4)
print(Figure4l)
dev.off()


##----------------------------------
## Figure 4M Heatmap for exhaustive
## genest
##----------------------------------
Exhausted_T <- c('PD1',
                 'CTLA4',
                 'CD38',
                 'CXCL13',
                 'MKI67',
                 'GZMB',
                 'BATF4',
                 'TOX',
                 'IRF4',
                 'CCL3',
                 'CXCR6',
                 'CSF1',
                 'IL13',
                 'GNLY',
                 'FASLG',
                 'TNF',
                 'LAG3',
                 'TIGIT')

Exhausted_T_new <- c('CD38', 'MKI67', 'GZMB', 'TOX', 'CCL3', 'CXCR6', 'CSF1',
                     'GNLY', 'FASLG',  'TNF', 'LAG3', 'TIGIT')

## For all the T cells
Idents(combined_Tcells) <- combined_Tcells$Status_new
Figure4m <- DotPlot(combined_Tcells, features = Exhausted_T) + 
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


pdf(paste0('Figure4/Figure4m.pdf'), width = 8, height = 3)
print(Figure4m)
dev.off()

