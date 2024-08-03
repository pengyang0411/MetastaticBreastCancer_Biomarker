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

##------------------------------------------
## Figure 5A: Umap of NK cells
##------------------------------------------

Idents(combined_nontumor) <- combined_nontumor$CT_L1
combined_NK <- subset(combined_nontumor, idents = 'NK')

Idents(combined_NK) <- combined_NK$CT_L2
combined_NK <- subset(combined_NK, 
                      idents = c('NK', 'NK Proliferating', 'NK_CD56bright'))

combined_NK <- RunTSNE(combined_NK, reduction = 'umap.pca', dims.use = 1:10, perplexity = 300)

pt5a <- DimPlot(combined_NK, reduction = 'tsne', group.by = 'CT_L2', 
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
  theme(strip.text.x = element_blank()) + theme(legend.position = 'bottom')



pdf('Figure5/Figure5a.pdf', height = 5, width = 5)
print(pt5a)
dev.off()


##-------------------------------------
## Figure 5B: Cell fraction across
## NK cells
##-------------------------------------


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

cells_NK <- c('ILC', 'NK', 'NK Proliferating', 'NK_CD56bright')

df_nontumor_L2_frac <- df_nontumor_L2

for(i in unique(df_nontumor_L2$Sample)){
  
  Indx <- which(df_nontumor_L2$Sample == i)
  
  df_nontumor_L2_frac[Indx, 'Count'] = df_nontumor_L2[Indx, 'Count'] / sum(df_nontumor_L2[Indx, 'Count'])
  
}

df_nontumor_NK_sub <- df_nontumor_L2_frac[df_nontumor_L2_frac$CellTypes %in% cells_NK, ]

## Adding stat test
stat.test <- df_nontumor_NK_sub %>%
  group_by(CellTypes) %>%
  t_test(Count ~ Status) %>%
  # adjust_pvalue(method = "BH") %>%
  add_significance("p")
stat.test$p.adj.signif <- stat.test$p.signif
stat.test

# Add p-values onto the box plots
stat.test <- stat.test %>%
  add_xy_position(x = "CellTypes", dodge = 0.5)

dodge <- position_dodge(width = 0.5)
Figure5b <- ggplot(df_nontumor_NK_sub) + 
  geom_boxplot(size = 0.3, outlier.shape = NA, width = 0.618, lwd = 0.5,
               aes(y = Count, x = CellTypes, color = Status)) +
  facet_wrap(~category, scales = 'free') + 
  geom_jitter(shape = 16, position = position_jitterdodge(), 
              size = 2,  alpha = 0.9, aes(fill = Status, y = Count, x = CellTypes, color = Status)) +
  theme_classic() + 
  ggtitle(" ") +
  xlab(" ") + ylab("Cell Fraction")  +
  theme(legend.position="bottom",
        # legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        panel.border = element_blank(),
        strip.background = element_blank(),
  ) + guides(color=guide_legend(title = 'Sample Status'),
             fill = FALSE)   +
  scale_fill_manual(values=c('#cccc66', '#66cccc', '#ff9966'))  +
  scale_color_manual(values=c('#cccc66', '#66cccc', '#ff9966')) #+
# stat_pvalue_manual(
#   stat.test,  label = "p.signif", tip.length = 0.01, hide.ns = T,
# )

pdf(paste0('Figure5/Figure5b.pdf'), width = 8, height = 3)
print(Figure5b)
dev.off()


##-------------------------------------
## Figure 5C: Heatmap of NK cells
##-------------------------------------


Idents(combined_NK) <- combined_NK$CT_L2
if(file.exists('Figure5/Data/cell_markers_NK_L2.RData')){
  load('Figure5/Data/cell_markers_NK_L2.RData')
}else{
  merged_combined.maker <- FindAllMarkers(combined_NK, only.pos = TRUE, 
                                          min.pct = 0.10, logfc.threshold = 0.25)
  save(merged_combined.maker, file = 'Figure5/Data/cell_markers_NK_L2.RData')
  
}

merged_combined.maker %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

## Make the heatmap for marker genes
merged_combined.maker %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top1

Figure5c <- DotPlot(combined_NK, features = unique(top1$gene)) + 
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
    legend.position="none",
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


pdf(paste0('Figure5/Figure5c.pdf'), width = 8, height = 2.2)
print(Figure5c)
dev.off()


##-------------------------------------
## Figure 5D: MP anlaysis of 
## Bcells, NK cells, M1 and M2 cells
## from plerual effusion
##-------------------------------------


##-----------B Cells-----------##
Idents(combined_nontumor) <- combined_nontumor$CT_L1
combined_B <- subset(combined_nontumor, idents = 'B')

Idents(combined_B) <- combined_B$CT_L2
combined_B <- subset(combined_B, idents = c('B intermediate', 'B memory', 'B naive',
                                            'Plasmablast'))

##-----------M1 Cells-----------##
Idents(combined_nontumor) <- combined_nontumor$CT_L1
combined_Myeloid <- subset(combined_nontumor, idents = c('Mono', 'DC'))

Idents(combined_Myeloid) <- combined_Myeloid$CT_L2
combined_Myeloid <- subset(combined_Myeloid, 
                           idents = c('CD16 Mono', 'cDC1', 'cDC2', 'pDC',
                                      'Macrophage'))

Idents(combined_Myeloid) <- combined_Myeloid$CT_L2
combined_M1 <- subset(combined_Myeloid, 
                      idents = c('CD16 Mono', 'Macrophage'))

##-----------M1 Cells-----------##
Idents(combined_Myeloid) <- combined_Myeloid$CT_L2
combined_M2 <- subset(combined_Myeloid, 
                      idents = c('cDC1', 'cDC2', 'pDC'))



MP_meta <- readxl::read_excel('/Users/yangpeng/Box Sync/CDK4/Data_nonTumor/T_cell MPs.xlsx')
MP_meta <- as.data.frame(MP_meta)
org_MP_Names <- names(MP_meta)
names(MP_meta) = paste0('MP', 1:ncol(MP_meta))

##----------- B cells -----------##
Idents(combined_B) <- combined_B$metastatic
pleural_combined <- subset(combined_B, idents = 'Pleural effusion')

Idents(pleural_combined) <- pleural_combined$Status_new


pleural_combined <- AddModuleScore(
  object = pleural_combined,
  features = as.list(MP_meta),
  ctrl = 5,
  name = 'MP'
)

res_MP <- array(0, c(ncol(MP_meta) * 2, 4))

for(j in 1:ncol(MP_meta)){
  
  ## meta data
  res_MP[j,1] <- 'LP vs BL' #paste(pairs[i,], collapse = '_vs_')
  res_MP[j,2] <- org_MP_Names[j]
  
  ## Calculate the median difference
  res_MP[j,3] <- median(pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Late Progressor',1]) - 
    median(pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Baseline',1])
  
  res_MP[j,4] <- t.test(pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Late Progressor',1],
                        pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Baseline',1])$p.value
  
  i = j + ncol(MP_meta)
  
  ## meta data
  res_MP[i,1] <- 'LP vs EP' #paste(pairs[i,], collapse = '_vs_')
  res_MP[i,2] <- org_MP_Names[j]
  
  ## Calculate the median difference
  res_MP[i,3] <- median(pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Late Progressor',1]) - 
    median(pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Early Progressor',1])
  
  res_MP[i,4] <- t.test(pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Late Progressor',1],
                        pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Early Progressor',1])$p.value
}

res_MP_df <- data.frame(comparison = as.factor(res_MP[,1]),
                        MP = as.factor(res_MP[,2]),
                        changes = as.numeric(res_MP[,3]),
                        p.value = as.numeric(res_MP[j,4]))
library(dplyr)
res_MP_df <- res_MP_df %>%
  as_tibble() %>%
  dplyr::arrange(desc(changes))

df_MP_B <- as.data.frame(res_MP_df)
df_MP_B$celltype = 'B cells'

##----------- NK cells -----------##
Idents(combined_NK) <- combined_NK$metastatic
pleural_combined <- subset(combined_NK, idents = 'Pleural effusion')

Idents(pleural_combined) <- pleural_combined$Status_new


pleural_combined <- AddModuleScore(
  object = pleural_combined,
  features = as.list(MP_meta),
  ctrl = 5,
  name = 'MP'
)

res_MP <- array(0, c(ncol(MP_meta) * 2, 4))

for(j in 1:ncol(MP_meta)){
  
  ## meta data
  res_MP[j,1] <- 'LP vs BL' #paste(pairs[i,], collapse = '_vs_')
  res_MP[j,2] <- org_MP_Names[j]
  
  ## Calculate the median difference
  res_MP[j,3] <- median(pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Late Progressor',1]) - 
    median(pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Baseline',1])
  
  res_MP[j,4] <- t.test(pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Late Progressor',1],
                        pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Baseline',1])$p.value
  
  i = j + ncol(MP_meta)
  
  ## meta data
  res_MP[i,1] <- 'LP vs EP' #paste(pairs[i,], collapse = '_vs_')
  res_MP[i,2] <- org_MP_Names[j]
  
  ## Calculate the median difference
  res_MP[i,3] <- median(pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Late Progressor',1]) - 
    median(pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Early Progressor',1])
  
  res_MP[i,4] <- t.test(pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Late Progressor',1],
                        pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Early Progressor',1])$p.value
}

res_MP_df <- data.frame(comparison = as.factor(res_MP[,1]),
                        MP = as.factor(res_MP[,2]),
                        changes = as.numeric(res_MP[,3]),
                        p.value = as.numeric(res_MP[j,4]))
library(dplyr)
res_MP_df <- res_MP_df %>%
  as_tibble() %>%
  dplyr::arrange(desc(changes))

df_MP_NK <- as.data.frame(res_MP_df)
df_MP_NK$celltype = 'NK cells'

##----------- M1 cells -----------##
Idents(combined_M1) <- combined_M1$metastatic
pleural_combined <- subset(combined_M1, idents = 'Pleural effusion')

Idents(pleural_combined) <- pleural_combined$Status_new


pleural_combined <- AddModuleScore(
  object = pleural_combined,
  features = as.list(MP_meta),
  ctrl = 5,
  name = 'MP'
)

res_MP <- array(0, c(ncol(MP_meta) * 2, 4))

for(j in 1:ncol(MP_meta)){
  
  ## meta data
  res_MP[j,1] <- 'LP vs BL' #paste(pairs[i,], collapse = '_vs_')
  res_MP[j,2] <- org_MP_Names[j]
  
  ## Calculate the median difference
  res_MP[j,3] <- median(pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Late Progressor',1]) - 
    median(pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Baseline',1])
  
  res_MP[j,4] <- t.test(pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Late Progressor',1],
                        pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Baseline',1])$p.value
  
  i = j + ncol(MP_meta)
  
  ## meta data
  res_MP[i,1] <- 'LP vs EP' #paste(pairs[i,], collapse = '_vs_')
  res_MP[i,2] <- org_MP_Names[j]
  
  ## Calculate the median difference
  res_MP[i,3] <- median(pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Late Progressor',1]) - 
    median(pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Early Progressor',1])
  
  res_MP[i,4] <- t.test(pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Late Progressor',1],
                        pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Early Progressor',1])$p.value
}

res_MP_df <- data.frame(comparison = as.factor(res_MP[,1]),
                        MP = as.factor(res_MP[,2]),
                        changes = as.numeric(res_MP[,3]),
                        p.value = as.numeric(res_MP[j,4]))
library(dplyr)
res_MP_df <- res_MP_df %>%
  as_tibble() %>%
  dplyr::arrange(desc(changes))

df_MP_M1 <- as.data.frame(res_MP_df)
df_MP_M1$celltype = 'M1 cells'

##----------- M2 cells -----------##
Idents(combined_M2) <- combined_M2$metastatic
pleural_combined <- subset(combined_M2, idents = 'Pleural effusion')

Idents(pleural_combined) <- pleural_combined$Status_new


pleural_combined <- AddModuleScore(
  object = pleural_combined,
  features = as.list(MP_meta),
  ctrl = 5,
  name = 'MP'
)

res_MP <- array(0, c(ncol(MP_meta) * 2, 4))

for(j in 1:ncol(MP_meta)){
  
  ## meta data
  res_MP[j,1] <- 'LP vs BL' #paste(pairs[i,], collapse = '_vs_')
  res_MP[j,2] <- org_MP_Names[j]
  
  ## Calculate the median difference
  res_MP[j,3] <- median(pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Late Progressor',1]) - 
    median(pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Baseline',1])
  
  res_MP[j,4] <- t.test(pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Late Progressor',1],
                        pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Baseline',1])$p.value
  
  i = j + ncol(MP_meta)
  
  ## meta data
  res_MP[i,1] <- 'LP vs EP' #paste(pairs[i,], collapse = '_vs_')
  res_MP[i,2] <- org_MP_Names[j]
  
  ## Calculate the median difference
  res_MP[i,3] <- median(pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Late Progressor',1]) - 
    median(pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Early Progressor',1])
  
  res_MP[i,4] <- t.test(pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Late Progressor',1],
                        pleural_combined[[paste0('MP', j)]][pleural_combined$Status_new == 'Early Progressor',1])$p.value
}

res_MP_df <- data.frame(comparison = as.factor(res_MP[,1]),
                        MP = as.factor(res_MP[,2]),
                        changes = as.numeric(res_MP[,3]),
                        p.value = as.numeric(res_MP[j,4]))
library(dplyr)
res_MP_df <- res_MP_df %>%
  as_tibble() %>%
  dplyr::arrange(desc(changes))

df_MP_M2 <- as.data.frame(res_MP_df)
df_MP_M2$celltype = 'M2 cells'


## Making it to dataframe
df_MP <- rbind(df_MP_B, df_MP_NK, df_MP_M1, df_MP_M2)
df_MP <- as.data.frame(df_MP)
df_MP$celltype <- factor(df_MP$celltype, 
                         levels = c('B cells', 'NK cells', 'M1 cells', 'M2 cells'))

## generate the plot
pt5d <- ggplot(df_MP, aes(comparison, reorder(MP, changes), fill = changes)) +
  geom_tile() +
  ggtitle('MP Analysis') +
  facet_wrap(~celltype, nrow = 1) + 
  # scale_fill_viridis(discrete = FALSE) +
  scale_fill_gradientn(
    ## limits = c(-1.5, 1.5),
    colors = c("#5DBCFF", "#6DCCFF", "white", "#F494BE", "#F484AE")) +
  xlab(' ') + ylab(' ') +
  guides(fill=guide_legend(title="Median difference in module score")) +
  theme(legend.position="right",
        axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1, size = 10),
        # panel.border = element_rect(fill = 'transparent'),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        # strip.text.x = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10)
  ) +  scale_colour_gradient(high = 'dark orange', low = 'dark blue')# +
# scale_x_discrete(labels=c('LP vs BL', 'LP vs EP')) #+
# scale_y_discrete(labels=rev(orders))
pdf('Figure5/Figure5d.pdf', height = 6, width = 6.5)
print(pt5d)
dev.off()

##------------------------------------------
## Figure 5e: Umap of Myeloid cells
##------------------------------------------

combined_Myeloid <- RunTSNE(combined_Myeloid, reduction = 'umap.pca', dims.use = 1:10, perplexity = 200)

pt5e <- DimPlot(combined_Myeloid, reduction = 'tsne', group.by = 'CT_L2', 
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
  theme(strip.text.x = element_blank()) + theme(legend.position = 'bottom')

pdf('Figure5/Figure5e.pdf', height = 5, width = 5)
print(pt5e)
dev.off()


##------------------------------------------
## Figure 5f: cell fraction of M1 cells
##------------------------------------------

cells_M1 <- c('CD16 Mono', 'Marcophage')

df_nontumor_M1_sub <- df_nontumor_L2_frac[df_nontumor_L2_frac$CellTypes %in% cells_M1, ]

## Adding stat test
stat.test <- df_nontumor_M1_sub %>%
  group_by(CellTypes) %>%
  t_test(Count ~ Status) %>%
  # adjust_pvalue(method = "BH") %>%
  add_significance("p")
stat.test$p.adj.signif <- stat.test$p.signif
stat.test

# Add p-values onto the box plots
stat.test <- stat.test %>%
  add_xy_position(x = "CellTypes", dodge = 0.5)

dodge <- position_dodge(width = 0.5)
Figure5f <- ggplot(df_nontumor_M1_sub) + 
  geom_boxplot(size = 0.3, outlier.shape = NA, width = 0.618, lwd = 0.5,
               aes(y = Count, x = CellTypes, color = Status)) +
  facet_wrap(~category, scales = 'free') + 
  geom_jitter(shape = 16, position = position_jitterdodge(), 
              size = 2,  alpha = 0.9, aes(fill = Status, y = Count, x = CellTypes, color = Status)) +
  theme_classic() + 
  ggtitle(" ") +
  xlab(" ") + ylab("Cell Fraction")  +
  theme(legend.position="bottom",
        # legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        panel.border = element_blank(),
        strip.background = element_blank(),
  ) + guides(color=guide_legend(title = 'Sample Status'),
             fill = FALSE)   +
  scale_fill_manual(values=c('#cccc66', '#66cccc', '#ff9966'))  +
  scale_color_manual(values=c('#cccc66', '#66cccc', '#ff9966')) #+
# stat_pvalue_manual(
#   stat.test,  label = "p.signif", tip.length = 0.01, hide.ns = T,
# )

pdf(paste0('Figure5/Figure5f.pdf'), width = 8, height = 3)
print(Figure5f)
dev.off()



##------------------------------------------
## Figure 5G: cell fraction of M2 cells
##------------------------------------------

cells_M2 <- c('cDC1', 'cDC2', 'pDC')

df_nontumor_M2_sub <- df_nontumor_L2_frac[df_nontumor_L2_frac$CellTypes %in% cells_M2, ]

## Adding stat test
stat.test <- df_nontumor_M2_sub %>%
  group_by(CellTypes) %>%
  t_test(Count ~ Status) %>%
  # adjust_pvalue(method = "BH") %>%
  add_significance("p")
stat.test$p.adj.signif <- stat.test$p.signif
stat.test

# Add p-values onto the box plots
stat.test <- stat.test %>%
  add_xy_position(x = "CellTypes", dodge = 0.5)

dodge <- position_dodge(width = 0.5)
Figure5g <- ggplot(df_nontumor_M2_sub) + 
  geom_boxplot(size = 0.3, outlier.shape = NA, width = 0.618, lwd = 0.5,
               aes(y = Count, x = CellTypes, color = Status)) +
  facet_wrap(~category, scales = 'free') + 
  geom_jitter(shape = 16, position = position_jitterdodge(), 
              size = 2,  alpha = 0.9, aes(fill = Status, y = Count, x = CellTypes, color = Status)) +
  theme_classic() + 
  ggtitle(" ") +
  xlab(" ") + ylab("Cell Fraction")  +
  theme(legend.position="bottom",
        # legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        panel.border = element_blank(),
        strip.background = element_blank(),
  ) + guides(color=guide_legend(title = 'Sample Status'),
             fill = FALSE)   +
  scale_fill_manual(values=c('#cccc66', '#66cccc', '#ff9966'))  +
  scale_color_manual(values=c('#cccc66', '#66cccc', '#ff9966')) #+
# stat_pvalue_manual(
#   stat.test,  label = "p.signif", tip.length = 0.01, hide.ns = T,
# )

pdf(paste0('Figure5/Figure5g.pdf'), width = 8, height = 3)
print(Figure5g)
dev.off()


##------------------------------------------
## Figure 5H: Tsen plot of Myeloid cells
##------------------------------------------

combined_Myeloid <- RunTSNE(combined_Myeloid, reduction = 'umap.pca', dims.use = 1:10, perplexity = 200)

pt5h <- DimPlot(combined_Myeloid, reduction = 'tsne', group.by = 'CT_L2', 
                 split.by = 'Status_new',
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
  theme(strip.text.x = element_blank()) + theme(legend.position = 'bottom')

pdf('Figure5/Figureh.pdf', height = 5, width = 15)
print(pt5h)
dev.off()


##------------------------------------------
## Figure 5I: Heatmap of Myeloid cells
##------------------------------------------
Idents(combined_M1) <- combined_M1$CT_L2
if(file.exists('Figure5/Data/cell_markers_M1_L2.RData')){
  load('Figure5/Data/cell_markers_M1_L2.RData')
}else{
  merged_combined.maker <- FindAllMarkers(combined_M1, only.pos = TRUE, 
                                          min.pct = 0.10, logfc.threshold = 0.25)
  save(merged_combined.maker, file = 'Figure5/Data/cell_markers_M1_L2.RData')
  
}

merged_combined.maker %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

## Make the heatmap for marker genes
merged_combined.maker %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> topM1

Idents(combined_M2) <- combined_M2$CT_L2
if(file.exists('Figure5/Data/cell_markers_M2_L2.RData')){
  load('Figure5/Data/cell_markers_M2_L2.RData')
}else{
  merged_combined.maker <- FindAllMarkers(combined_M2, only.pos = TRUE, 
                                          min.pct = 0.10, logfc.threshold = 0.25)
  save(merged_combined.maker, file = 'Figure5/Data/cell_markers_M2_L2.RData')
  
}

merged_combined.maker %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

## Make the heatmap for marker genes
merged_combined.maker %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> topM2

Figure5i <- DotPlot(combined_Myeloid, features = unique(topM1$gene, topM2$gene)) + 
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
    legend.position="none",
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


pdf(paste0('Figure5/Figure5i.pdf'), width = 8, height = 3)
print(Figure5i)
dev.off()


##------------------------------------------
## Figure 5J: Violin plot of M1 cells
##------------------------------------------

## Run MP analysis
Myeloid_markders <- c(readxl::read_excel('Figure5/Data/Myeloid_markders.xlsx', sheet = 1)$Gene,
                      readxl::read_excel('Figure5/Data/Myeloid_markders.xlsx', sheet = 2)$Gene)


MP_meta_list <- list()
for(i in 1:2){
  
  MP_meta_list[[i]] <- readxl::read_excel('Figure6/Data/Myeloid_markders.xlsx', sheet = 1)$Gene
  
}
org_MP_Names <- c('M1', 'M2')


combined_M1 <- AddModuleScore(
  object = combined_M1,
  features = MP_meta_list,
  ctrl = 5,
  name = 'M'
)

colorForStatus <- c('#cccc66', '#66cccc', '#ff9966')
my_comparisons <- list( c("Baseline", "Early Progressor"), 
                        c("Baseline", "Late Progressor"), 
                        c("Early Progressor", "Late Progressor") )
Idents(combined_M1) <- combined_M1$Status_new

pdf('Figure5/Figure5j.pdf', width = 5, height = 5)
VlnPlot(combined_M1, features = 'M1', pt.size = 0.001) + 
  stat_compare_means(comparisons = my_comparisons) + ylim(c(-0.6, 1.5)) +
  scale_fill_manual(
    values = colorForStatus) + NoLegend() + ggtitle('M1') + xlab(' ') + 
  stat_summary(fun.y = median, geom='point', size = 20, colour = "white", shape = 95)
dev.off()

##------------------------------------------
## Figure 5K: Violin plot of M2 cells
##------------------------------------------

pdf('Figure6/Figure6_M2.pdf', width = 5, height = 5)
VlnPlot(combined_M1, features = 'M2', pt.size = 0.001) + 
  stat_compare_means(comparisons = my_comparisons) + ylim(c(-0.6, 1.5)) +
  scale_fill_manual(
    values = colorForStatus) + NoLegend() + ggtitle('M2') + xlab(' ') + 
  stat_summary(fun.y = median, geom='point', size = 20, colour = "white", shape = 95)
dev.off()


