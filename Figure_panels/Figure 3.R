##-------------------------------
## Figure Panel 3
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

library(Signac)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
library(EnhancedVolcano)
load('Data/combined_nonTumor_final.RData')


cell_colors <- c(
                 "#66cccc", # CD4 T
                 "#F9B26C", # CD8 T
                 "#cccc66", # other T
                 "#ffcccc", # B
                 "#F494BE", # Mono
                 "#ff9966", # DC
                 "#ff99cc", # NK
                 "#6DCCDD", # CAFs
                 "#9999cc", # PVL
                 "#EDCAE0", # Endothelial
                 "grey") # other

cell_types <- c( "CD4 T", "CD8 T", "other T", "B", "Mono", "DC", "NK","CAFs", 
                 "PVL", "Endothelial","other")

combined_nontumor$CT_L1 <- factor(as.character(combined_nontumor$CT_L1), 
                                  levels = cell_types)

##-------------------------------
## Figure 3A: Umap plot across 
## different sample status
##-------------------------------

Idents(combined_nontumor) <- combined_nontumor$CT_L1
combined_nontumor$CT_L1.num <- as.factor(as.numeric(as.factor(combined_nontumor$CT_L1)))

pt3a <- DimPlot(combined_nontumor, reduction = 'umap.harmony', group.by = 'CT_L1.num', 
                split.by = 'Status_new', 
                label = T, repel = T, label.size = 5) + 
  NoAxes() + ggtitle(' ') + 
  scale_x_reverse() +
  theme_void() + 
  scale_color_manual(# breaks = unique(merged_combined$orig.ident.num),
    labels = paste(1:length(unique(combined_nontumor$CT_L1)), cell_types),
    values = cell_colors) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 6),
                              label.hjust = 0,
                              label.theme = element_text())) + NoLegend() +
  theme(strip.text.x = element_blank())

pt3a

pdf('Figure3/Figure3a.pdf', width = 12, height = 4)
print(pt3a)
dev.off()

##-------------------------------
## Figure 3B: Umap plot of genes
## across different cell types
##-------------------------------

Idents(combined_nontumor) <- combined_nontumor$CT_L1

if(file.exists('Figure3/Data/cell_markers.RData')){
  load('Figure3/Data/cell_markers.RData')
}else{
  merged_combined.maker <- FindAllMarkers(combined_nontumor, only.pos = TRUE, 
                                          min.pct = 0.10, logfc.threshold = 0.25)
  save(merged_combined.maker, file = 'Figure3/Data/cell_markers.RData')
  
}

merged_combined.maker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

## Make the heatmap for marker genes
merged_combined.maker %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC) -> top1

identical(as.character(top1$cluster), cell_types)


Figure3b <- DotPlot(combined_nontumor, features = unique(top1$gene)) + 
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


pdf(paste0('Figure3/Figure3b.pdf'), width = 8, height = 5)
print(Figure3b)
dev.off()

##-------------------------------
## Figure 3C: Relative non-tumor
## cell abundance across samples
##-------------------------------


## Create the dataframe
df_nontumor_L1 <- c()
for(sample in unique(combined_nontumor$orig.ident)){
  df_nontumor_L1 <- rbind(df_nontumor_L1, 
                          cbind(cbind(table(combined_nontumor$CT_L1[combined_nontumor$orig.ident == sample]), sample),
                                as.character(combined_nontumor$Status_new[combined_nontumor$orig.ident == sample][1])))
}

df_nontumor_L1 <- cbind(df_nontumor_L1, rownames(df_nontumor_L1))

colnames(df_nontumor_L1) <- c('Count', 'Sample', 'Status', 'CellTypes')
df_nontumor_L1 <- as.data.frame(df_nontumor_L1)

## Make factors
df_nontumor_L1$Count <- as.numeric(df_nontumor_L1$Count)
df_nontumor_L1$Sample <- as.factor(df_nontumor_L1$Sample)
df_nontumor_L1$Status <- as.factor(df_nontumor_L1$Status)
df_nontumor_L1$CellTypes <- factor(df_nontumor_L1$CellTypes, 
                                   cell_types)



## make the plot
pt3c <- ggplot(df_nontumor_L1, aes(fill = CellTypes, y = Count, x = Sample)) +
  geom_bar(position="fill", stat="identity") +
  facet_grid(~ Status, drop = T, scales = 'free', space = 'free') +
  theme_classic() +
  scale_fill_manual(values = cell_colors) +
  xlab(" ") + ylab(" ")  +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
    strip.background = element_blank(),
    strip.text.x = element_blank()
    # strip.text.x = element_text(size = 10)
  ) + NoLegend()

pdf('Figure3/Figure3c.pdf', width = 7, height = 3)
print(pt3c)
dev.off()


##-------------------------------
## Figure 3D: Relative non-tumor
## cell distribution
##-------------------------------

cells_1 <- c('CD4 T', 'CD8 T', 'B',  'Mono', 'DC', 'NK', 'CAFs')
df_nontumor_L1_frac <- df_nontumor_L1
for(i in unique(df_nontumor_L1$Sample)){
  
  Indx <- which(df_nontumor_L1$Sample == i)
  
  df_nontumor_L1_frac[Indx, 'Count'] = df_nontumor_L1[Indx, 'Count'] / sum(df_nontumor_L1[Indx, 'Count'])
  
}

df_nontumor_L1_frac_sub <- df_nontumor_L1_frac[df_nontumor_L1_frac$CellTypes %in% cells_1, ]
df_nontumor_L1_frac_sub$CellTypes <- factor(as.character(df_nontumor_L1_frac_sub$CellTypes), levels = cells_1)

## Adding stat test
stat.test <- df_nontumor_L1_frac_sub %>%
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
Figure3d <- ggplot(df_nontumor_L1_frac_sub) + 
  geom_boxplot(size = 0.3, outlier.shape = NA, width = 0.618, lwd = 0.5,
               aes(y = Count, x = CellTypes, color = Status)) +
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
        # axis.text.x = element_text(angle = 45, hjust=1),
        panel.border = element_blank(),
        # strip.background = element_blank(),
        # strip.text.x = element_blank(),
        # legend.title = element_text(color = "transparent"),
        # legend.text = element_text(color = "transparent"),
  ) + guides(color=guide_legend(title = 'Sample Status'),
             fill = FALSE)  + 
  scale_fill_manual(values=c('#cccc66', '#66cccc', '#ff9966'))  +  
  scale_color_manual(values=c('#cccc66', '#66cccc', '#ff9966')) + 
  stat_pvalue_manual(
    stat.test,  label = "p.signif", tip.length = 0.01, hide.ns = T,
  )

pdf('Figure3/Figure3d.pdf', width = 7, height = 3.5)
print(Figure3d)
dev.off()

##-------------------------------
## Figure 3E: MP Analysis for
## T cells
##-------------------------------


df_MP <- c()


pairs <- array(0, c(3,2))
pairs[1, ] <- c('Late Progressor', 'Early Progressor')
pairs[2, ] <- c('Late Progressor',  'Baseline')
pairs[3, ] <- c('Early Progressor', 'Baseline')


for(j in 1:3){
  
  tmp      <- read.csv(paste0('Data/T_cell_MP_', pairs[j,1], '_', pairs[j,2],'.csv'))
  tmp$pair <- paste0(pairs[j,1], ' vs ', pairs[j,2])
  
  df_MP <- rbind(df_MP, tmp)
  
}


df_MP$comparison <- as.factor(df_MP$comparison)
df_MP$MP         <- as.factor(df_MP$MP)
df_MP$changes    <- as.numeric(df_MP$changes)
df_MP$pair       <- as.factor(df_MP$pair)


Figure3e <- ggplot(df_MP, aes(pair, reorder(MP, changes), fill = changes)) +
  geom_tile() +
  ggtitle('MP Analysis') +
  # scale_fill_viridis(discrete = FALSE) +
  scale_fill_gradientn(
    limits = c(-0.3, 0.3),
    colors = c("#5DBCFF", "#6DCCFF", "white", "#F494BE", "#F484AE")) +
  xlab(' ') + ylab(' ') +
  guides(fill=guide_legend(title="Median difference in module score")) +
  theme(legend.position="right",
        axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1, size = 10),
        # panel.border = element_rect(fill = 'transparent'),
        panel.border = element_blank(),
        panel.background = element_blank(),
        # strip.background = element_blank(),
        # strip.text.x = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10)
  ) +  scale_colour_gradient(high = 'dark orange', low = 'dark blue') +
  scale_x_discrete(labels=c('EP vs BL', 'LP vs BL', 'LP vs EP')) #+
# scale_y_discrete(labels=rev(orders))



pdf(paste0('Figure3/Figure3e.pdf'), width = 5.5, height = 5)
print(Figure3e)
dev.off()


##-------------------------------
## Figure 3F: Ligand-receptor 
## interaction analysis 
##-------------------------------
for(s in c('Baseline', 'Early Progressor', 'Late Progressor')){
  
  pdf(paste0('Figure3/Figure3f_cellchat_chord_', s, '.pdf'), width = 8, height = 8)
  print(netVisual_chord_gene(cellchat, sources.use = NULL, targets.use = NULL, signaling = pathways.show.all, legend.pos.x = 8,
                             scale = T, thresh = 0.05,
                             title.name = s))
  dev.off()
  
}


