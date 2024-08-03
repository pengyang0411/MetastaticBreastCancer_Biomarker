##-------------------------------
## Figure Panel S1
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
library(Matrix)
library(ggplot2)
library(patchwork)
library(EnhancedVolcano)

load('/Data/merged_data_Harmony_new.RData')

##-------------------------------
## Figure S1B: UMAP PCA across 
## samples
##-------------------------------

umapColor <- c("#6DCCDD","#9ACBDE","#C8CADF","#EDC6DD","#F092B1","#F27FA5",
               "#F47892", "#F6A395","#F8AD77","#E7B080","#C9AFA2","#ABADC5",
               "#AEB9AC","#B9C984", "#C4DA5D")
SampleNames <- c("PA3",   "PA5",   "PA45",  "PA95",  "PA120", "PA110", "PA131", "PA153",
                 'PA11', 'PA46', 'PA144',
                 'PA1', 'PA3 #1', 'PA3 #2',
                 'PA126', 'PA130', 'PA139', 'PA99', 'PA165')

merged_combined$orig.ident <- factor(merged_combined$orig.ident, 
                                     levels = SampleNames)


Idents(merged_combined) <- merged_combined$orig.ident
merged_combined$orig.ident.num <- as.factor(as.numeric(as.factor(merged_combined$orig.ident)))
# merged_combined <- subset(merged_combined, idents = 'PA130', invert = T)
ptS1b <- DimPlot(merged_combined, reduction = 'umap.pca', split.by = 'orig.ident', 
                  ncol = 6,
                label = F, repel = F, label.size = 5) + 
  NoAxes() + ggtitle(' ') + 
  scale_x_reverse() +
  theme_void() + 
  scale_color_manual(# breaks = unique(merged_combined$orig.ident.num),
    labels = paste(1:length(unique(merged_combined$orig.ident)), SampleNames), 
    values = colorRampPalette(umapColor)(length(unique(merged_combined$orig.ident)))) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 6),
                              label.hjust = 0,
                              label.theme = element_text())) + NoLegend()

# pt1b

pdf('FigureS1/FigureS1b.pdf', width = 8, height = 4)
print(ptS1b)
dev.off()


## Subset of tumor cells
Idents(merged_combined) <- merged_combined$Tumor_cells
merged_combined <- subset(merged_combined, idents = 'tumor cells')


##-------------------------------
## Figure S1C: Volcano plot
## between Baseline vs Early Progressor
##-------------------------------


if(file.exists(paste0('FigureS1/Data/DE_tumor_EP_BL.csv'))){
  thisDE <- read.csv(paste0('FigureS1/Data/DE_tumor_EP_BL.csv'))
}else{
  thisDE <- FindMarkers(object = merged_combined, 
                        ident.1 = 'Early Progressor', ident.2 = 'Baseline', 
                        min.pct = 0, test.use = "MAST")
  thisDE$gene <- rownames(thisDE)
  write.csv(thisDE, file = paste0('FigureS1/Data/DE_tumor_EP_BL.csv'))
}

# thisDE$avg_log2FC <- thisDE$avg_log2FC + rnorm(length(thisDE$avg_log2FC ), 0, 0.1)

FigureS1c <- EnhancedVolcano(thisDE,
                             lab = thisDE$gene,
                             # selectLab = selected_genes,
                             x = 'avg_log2FC',
                             y = 'p_val_adj',
                             title = paste0('Baseline (Left) vs Early Progressor (Right)'),
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


pdf(paste0('FigureS1/FigureS1c.pdf'), width = 5, height = 4)
print(FigureS1c)
dev.off()



##-------------------------------
## Figure S1D: Volcano plot
## between Baseline vs LP
##-------------------------------


Idents(merged_combined) <- merged_combined$Status_new

if(file.exists(paste0('FigureS1/Data/DE_tumor_LP_BL.csv'))){
  thisDE <- read.csv(paste0('FigureS1/Data/DE_tumor_LP_BL.csv'))
}else{
  thisDE <- FindMarkers(object = merged_combined, 
                        ident.1 = 'Late Progressor', ident.2 = 'Baseline', 
                        min.pct = 0, test.use = "MAST")
  thisDE$gene <- rownames(thisDE)
  write.csv(thisDE, file = paste0('FigureS1/Data/DE_tumor_LP_BL.csv'))
}

# thisDE$avg_log2FC <- thisDE$avg_log2FC + rnorm(length(thisDE$avg_log2FC ), 0, 0.1)

FigureS1d <- EnhancedVolcano(thisDE,
                            lab = thisDE$gene,
                            # selectLab = selected_genes,
                            x = 'avg_log2FC',
                            y = 'p_val_adj',
                            title = paste0('Baseline (Left) vs Late Progressor (Right)'),
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


pdf(paste0('FigureS1/FigureS1d.pdf'), width = 5, height = 4)
print(FigureS1d)
dev.off()



##-------------------------------
## Figure S1E: Volcano plot
## between Early vs Late Progressor
##-------------------------------


if(file.exists(paste0('FigureS1/Data/DE_tumor_LP_EP.csv'))){
  thisDE <- read.csv(paste0('FigureS1/Data/DE_tumor_LP_EP.csv'))
}else{
  thisDE <- FindMarkers(object = merged_combined, 
                        ident.1 = 'Late Progressor', ident.2 = 'Early Progressor', 
                        min.pct = 0, test.use = "MAST")
  thisDE$gene <- rownames(thisDE)
  write.csv(thisDE, file = paste0('FigureS1/Data/DE_tumor_LP_EP.csv'))
}

# thisDE$avg_log2FC <- thisDE$avg_log2FC + rnorm(length(thisDE$avg_log2FC ), 0, 0.1)

FigureS1e <- EnhancedVolcano(thisDE,
                            lab = thisDE$gene,
                            # selectLab = selected_genes,
                            x = 'avg_log2FC',
                            y = 'p_val_adj',
                            title = paste0('Early (Left) vs Late Progressor (Right)'),
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


pdf(paste0('FigureS1/FigureS1e.pdf'), width = 5, height = 4)
print(FigureS1e)
dev.off()



##-------------------------------
## Figure S1F: MP Analysis
##-------------------------------


df_MP <- c()

pairs <- array(0, c(3,2))
pairs[1, ] <- c('Late Progression', 'Early Progression')
pairs[2, ] <- c('Baseline', 'Early Progression')
pairs[3, ] <- c('Baseline', 'Late Progression')

Labels <- array(0, c(3,2))
Labels[1, ] <- c('Late Progressor', 'Early Progressor')
Labels[2, ] <- c('Early Progressor', 'Baseline')
Labels[3, ] <- c('Late Progressor', 'Baseline')

for(j in 1:3){
  
  tmp      <- read.csv(paste0('/Data_Tumor/Tumor_MP_', pairs[j,1], '_', pairs[j,2],'.csv'))#[,c('pathway', 'NES')]
  tmp$pair <- paste0(Labels[j,1], ' vs ', Labels[j,2])
  
  df_MP <- rbind(df_MP, tmp)
  
}


df_MP$comparison <- as.factor(df_MP$pair)
df_MP$MP         <- as.factor(df_MP$MP)
df_MP$changes    <- as.numeric(df_MP$changes)
df_MP$padj       <- p.adjust(df_MP$p.value, method = 'BH')



df_MP$pair[df_MP$pair == 'Early Progressor vs Baseline'] <- 'EP vs BL'
df_MP$pair[df_MP$pair == 'Late Progressor vs Baseline'] <- 'LP vs BL'
df_MP$pair[df_MP$pair == 'Late Progressor vs Early Progressor'] <- 'LP vs EP'
df_MP$pair <- factor(df_MP$pair, levels = c('EP vs BL', 'LP vs BL', 'LP vs EP'))
MPmat <- xtabs(changes~MP+pair,df_MP[,c('MP','changes','pair')])
# orders <- c('MP14 EMT-III', 'MP8 Proteasomal degradation', 'MP6 Hypoxia',
#             'MP17 Interferon/MHC-II (I)', 'MP20 MYC', 'MP11 Translation initiation',
#             'MP22 Secreted I', 'MP40 Unassigned-I', 'MP34 Platelet-activation',
#             'MP2  Cell Cycle - G1/S', 'MP1  Cell Cycle - G2/M', 
#             'MP18 Interferon/MHC-II (II)', 'MP10 Protein maturation',
#             'MP38 Glutathione', 'MP19 Epithelial Senescence', 
#             'MP39 Metal-response', 'MP5 Stress', 'MP30 PDAC-classical')

ptS1f <- ggplot(df_MP, aes(pair, reorder(MP, changes), fill = changes)) +
  geom_tile() +
  ggtitle('MP Analysis') +
  # scale_fill_viridis(discrete = FALSE) +
  scale_fill_gradientn(
    limits = c(-0.31, 0.31),
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



pdf(paste0('FigureS1/FigureS1f.pdf'), width = 6, height = 5)
print(ptS1f)
dev.off()


