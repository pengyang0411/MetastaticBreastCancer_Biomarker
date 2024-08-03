##-------------------------------
## Figure Panel 2
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
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
library(EnhancedVolcano)
load('/merged_data_Harmony_new.RData')

##---------------------------------
## Figure 2A: pie chart
##---------------------------------

df_BL <- data.frame(n = c(0, 6, 0, 2),
                    perc = c(0, 6, 0, 2) / 8,
                    label = factor(c('Ascities', 'Liver', 'Pelvic bone', 'Pleural effusion'),
                                   levels = c('Ascities', 'Liver', 'Pelvic bone', 'Pleural effusion')))


pt2chart_a <- ggplot(df_BL, aes(x="", y = perc, fill = label)) +
  # facet_wrap(~Status) + 
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) + 
  xlab(" ") + ylab(" ")  +
  theme_void() + NoLegend()


df_EP <- data.frame(n = c(0, 1, 0, 2),
                    perc = c(0, 1, 0, 2) / 3,
                    label = factor(c('Ascities', 'Liver', 'Pelvic bone', 'Pleural effusion'),
                                   levels = c('Ascities', 'Liver', 'Pelvic bone', 'Pleural effusion')))

pt2chart_b <- ggplot(df_EP, aes(x="", y = perc, fill = label)) +
  # facet_wrap(~Status) + 
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) + 
  xlab(" ") + ylab(" ")  +
  theme_void() + NoLegend()

df_LP <- data.frame(n = c(2, 2, 1, 2),
                    perc = c(2, 2, 1, 2) / 7,
                    label = factor(c('Ascities', 'Liver', 'Pelvic bone', 'Pleural effusion'),
                                   levels = c('Ascities', 'Liver', 'Pelvic bone', 'Pleural effusion')))


pt2chart_c <- ggplot(df_LP, aes(x="", y = perc, fill = label)) +
  # facet_wrap(~Status) + 
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) + 
  xlab(" ") + ylab(" ")  +
  theme_void() + NoLegend()

pdf('Figure2/Figure2chat_a.pdf', width = 4, height = 4)
print(pt2chart_a)
dev.off()

pdf('Figure2/Figure2chat_b.pdf', width = 4, height = 4)
print(pt2chart_b)
dev.off()

pdf('Figure2/Figure2chat_c.pdf', width = 4, height = 4)
print(pt2chart_c)
dev.off()

pdf('Figure2/Figure2chat_legend.pdf', width = 5, height = 4)
ggplot(df_LP, aes(x="", y = perc, fill = label)) +
  # facet_wrap(~Status) + 
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) + 
  xlab(" ") + ylab(" ")  + theme(legend.position = 'right')
dev.off()



##---------------------------------
## Figure 2B
##---------------------------------
## Create the dataframe
df_tumor <- c()
for(sample in unique(merged_combined$orig.ident)){
  df_tumor <- rbind(df_tumor, c(sample, as.character(merged_combined$metastatic[merged_combined$orig.ident == sample][1]),
                                sum(merged_combined$Tumor_cells[merged_combined$orig.ident == sample] == 'tumor cells'), 'Tumor Cell'))
  df_tumor <- rbind(df_tumor, c(sample, as.character(merged_combined$metastatic[merged_combined$orig.ident == sample][1]),
                                sum(merged_combined$Tumor_cells[merged_combined$orig.ident == sample] != 'tumor cells'), 'Non-Tumor Cell'))
}
colnames(df_tumor) <- c('Sample', 'Status', 'Count', 'Condition')
df_tumor <- as.data.frame(df_tumor)
df_tumor$Status <- as.factor(df_tumor$Status)
df_tumor$Sample <- as.factor(df_tumor$Sample)
df_tumor$Count <- as.numeric(df_tumor$Count)
df_tumor$Condition <- factor(df_tumor$Condition, levels = c('Tumor Cell', 'Non-Tumor Cell'))


## make the plot
pt2b <- ggplot(df_tumor, aes(fill = Condition, y = Count, x = Sample)) +
  geom_bar(position="fill", stat="identity") +
  facet_grid(~ Status, drop = T, scales = 'free', space = 'free') +
  theme_classic() +
  scale_fill_manual(values = c('#FD9AA0', '#6DCCFD')) +
  xlab(" ") + ylab(" ")  +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
    strip.background = element_blank(),
    strip.text.x = element_blank()
    # strip.text.x = element_text(size = 10)
  ) + NoLegend()

pdf('Figure2/Figure2b.pdf', width = 6, height = 2.5)
print(pt2b)
dev.off()


## Subset of tumor cells
Idents(merged_combined) <- merged_combined$Tumor_cells
merged_combined <- subset(merged_combined, idents = 'tumor cells')

##-------------------------------
## Figure 2C: Heatmap across
## sample status, metastatic site, etc
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
  top_n(n = 7, wt = avg_log2FC) -> top10


Idents(merged_combined) <- merged_combined$metastatic
Figure2c <- DotPlot(merged_combined, features = top10$gene) + 
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

pdf(paste0('Figure2/Figure2c.pdf'), width = 8, height = 3.5)
print(Figure2c)
dev.off()

##-------------------------------
## Figure 2D: Tumor cells for
## liver and pleural effusion
##-------------------------------
Idents(merged_combined) <- merged_combined$metastatic
liver_combined <- subset(merged_combined, idents = 'Liver')
colorForStatus <- c('#cccc66', '#66cccc', '#ff9966')
## Making umap plot to seperate the status
pt2d1 <- DimPlot(liver_combined, reduction = 'umap.harmony', group.by = 'Status_new', 
        split.by = 'Status_new', ncol = 2,
        label = F, repel = T, label.size = 5) + 
  NoAxes() + ggtitle(' ') + 
  scale_x_reverse() +
  theme_void() + 
  scale_color_manual(
    values = colorForStatus) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 6),
                              # label.hjust = 10, label.vjust = 2,  title.vjust = 2,
                              label.theme = element_text())) + #NoLegend() +
  theme(strip.text.x = element_blank(),
        legend.position = 'none',
        legend.spacing.y = unit(2.0, 'cm'))

pdf('Figure2/Figure2d1.pdf', width = 10, height = 8)
print(pt2d1)
dev.off()


Idents(merged_combined) <- merged_combined$metastatic
pleural_combined <- subset(merged_combined, idents = 'Pleural effusion')

## Making umap plot to seperate the status
pt2d2 <- DimPlot(pleural_combined, reduction = 'umap.harmony', group.by = 'Status_new', 
                 split.by = 'Status_new', ncol = 2,
                 label = F, repel = T, label.size = 5) + 
  NoAxes() + ggtitle(' ') + 
  scale_x_reverse() +
  theme_void() + 
  scale_color_manual(
    values = colorForStatus) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 6),
                              label.hjust = 0, 
                              label.theme = element_text())) + #NoLegend() +
  theme(strip.text.x = element_blank(),
        legend.position = 'none')

pdf('Figure2/Figure2c2.pdf', width = 10, height = 8)
print(pt2d2)
dev.off()

##-------------------------------
## Figure 2D: Trajectory analysis 
## for liver and pleural effusion
##-------------------------------


## Trajectory analysis for liver
names(liver_combined@reductions)[4] = 'UMAP'

liver_combined <- as.cell_data_set(liver_combined)
liver_combined <- cluster_cells(cds = liver_combined, reduction_method = "UMAP")
liver_combined <- learn_graph(liver_combined, use_partition = TRUE)


set.seed(55)
RootCells <- rownames(merged_combined@meta.data)[merged_combined@meta.data$Status_new == 'Baseline' & merged_combined@meta.data$metastatic == 'Liver']
liver_combined <- order_cells(liver_combined, root_cells = sample(RootCells, replace = F, 1))


figure2e1 <- 
plot_cells(cds = liver_combined,
           color_cells_by = "pseudotime",
           show_trajectory_graph = T,
           trajectory_graph_color = "white",
           trajectory_graph_segment_size = 1,
           graph_label_size = 2,
           cell_size = 2,
           label_cell_groups = F,
           label_groups_by_cluster = F,
           label_branch_points = F,
           label_roots = T,
           label_leaves = T) +
  scale_color_gradientn(colours =
                          c("#0000DD", "#AAAAFF", "#EEEEEE", "#FFAAAA", "#AA0000"))+
  scale_x_reverse()+
  theme_void() + theme(text = element_text(size = 6),
                       legend.title = element_text(size = 10),
                       # legend.key.size = 10,
                       legend.position = "top") + 
  guides(color=guide_legend(title="Pseudotime")) + 
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10))


pdf(paste0('Figure2/Figure2e1.pdf'), width = 6, height = 6)
print(figure2e1)
dev.off()

## Trajectory analysis for pleural effusion
names(pleural_combined@reductions)[4] = 'UMAP'

pleural_combined <- as.cell_data_set(pleural_combined)
pleural_combined <- cluster_cells(cds = pleural_combined, reduction_method = "UMAP")
pleural_combined <- learn_graph(pleural_combined, use_partition = TRUE)

set.seed(55)
RootCells <- rownames(merged_combined@meta.data)[merged_combined@meta.data$Status_new == 'Baseline' & merged_combined@meta.data$metastatic == 'Pleural effusion']
pleural_combined <- order_cells(pleural_combined, root_cells = sample(RootCells, 1))


figure2e2 <- 
  plot_cells(cds = pleural_combined,
             color_cells_by = "pseudotime",
             show_trajectory_graph = T,
             trajectory_graph_color = "white",
             trajectory_graph_segment_size = 1,
             graph_label_size = 2,
             cell_size = 2,
             label_cell_groups = F,
             label_groups_by_cluster = F,
             label_branch_points = F,
             label_roots = T,
             label_leaves = T) +
  scale_color_gradientn(colours =
                          c("#0000DD", "#AAAAFF", "#EEEEEE", "#FFAAAA", "#AA0000"))+
  scale_x_reverse()+
  theme_void() + theme(text = element_text(size = 6),
                       legend.title = element_text(size = 10),
                       # legend.key.size = 10,
                       legend.position = "none") + 
  guides(color=guide_legend(title="Pseudotime")) + 
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10))



pdf(paste0('Figure2/Figure2e2.pdf'), width = 6, height = 6)
print(figure2e2)
dev.off()




##-------------------------------
## Figure 2F: Volcano plot for 
## Pleural effusion BL vs LP
##-------------------------------

Idents(merged_combined) <- merged_combined$metastatic
pleural_combined <- subset(merged_combined, idents = 'Pleural effusion')

Idents(pleural_combined) <- pleural_combined$Status_new


## Late Progressor vs Baseline
if(file.exists(paste0('Figure2/Data/DE_Pleural_LP_BL.csv'))){
  thisDE <- read.csv(paste0('Figure2/Data/DE_Pleural_LP_BL.csv'))
}else{
  thisDE <- FindMarkers(object = pleural_combined, 
                        ident.1 = 'Late Progressor', ident.2 = 'Baseline', 
                        min.pct = 0, test.use = "MAST")
  thisDE$gene <- rownames(thisDE)
  write.csv(thisDE, file = paste0('Figure2/Data/DE_Pleural_LP_BL.csv'))
}

# thisDE$avg_log2FC <- thisDE$avg_log2FC + rnorm(length(thisDE$avg_log2FC ), 0, 0.1)

Figure2f <- EnhancedVolcano(thisDE,
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
                            ylim = c(0, 80),
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


pdf(paste0('Figure2/Figure2f.pdf'), width = 5, height = 4)
print(Figure2f)
dev.off()

##-------------------------------
## Figure 2G: Volcano plot for 
## Pleural effusion BL vs LP
##-------------------------------

## Early Progressor versus Baseline
if(file.exists(paste0('Figure2/Data/DE_Pleural_EP_LP.csv'))){
  thisDE <- read.csv(paste0('Figure2/Data/DE_Pleural_EP_LP.csv'))
}else{
  thisDE <- FindMarkers(object = pleural_combined, 
                        ident.1 = 'Late Progressor', ident.2 = 'Early Progressor', 
                        min.pct = 0, test.use = "MAST")
  thisDE$gene <- rownames(thisDE)
  write.csv(thisDE, file = paste0('Figure2/Data/DE_Pleural_EP_LP.csv'))
}

# thisDE$avg_log2FC <- thisDE$avg_log2FC + rnorm(length(thisDE$avg_log2FC ), 0, 0.1)

Figure2g <- EnhancedVolcano(thisDE,
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
  scale_color_manual(values = c('white', 'white', '#6DCCFD', '#FD9AA0'))


pdf(paste0('Figure2/Figure2g.pdf'), width = 5, height = 4)
print(Figure2g)
dev.off()


##----------------------------
## Figure 2H: Hallmark pathway 
##----------------------------
## Hallmark pathway
pathways_list <- c('hallmark')
pathways_path <- c('Data/msigdb_v2022.1.Hs_GMTs/h.all.v2022.1.Hs.symbols.gmt')
## Load pathways gmt file
pathways <- gmtPathways(gmt.file = pathways_path[1])
DE_files <- c('Figure2/Data/DE_Pleural_EP_LP.csv',
              'Figure2/Data/DE_Pleural_LP_BL.csv')
df_hallmark <- c()
for(j in 1:2){
  
  thisDE <- read.csv(paste0(DE_files[j]))
  
  ## Run fGSEA
  TESTX<-thisDE
  TESTX$gene<-TESTX$gene
  T1.mkers<-TESTX
  res.t1 <- T1.mkers %>%
    dplyr::filter(p_val_adj < 0.01) %>%
    dplyr::select(gene, avg_log2FC) %>%
    tibble::deframe()
  fgseaRes.t1 <- fgsea(pathways=pathways, stats=res.t1, nperm=1000)
  fgseaResTidy.t1 <- fgseaRes.t1 %>%
    as_tibble() %>%
    dplyr::arrange(desc(NES))
  Indx <- c(order(fgseaResTidy.t1$NES)[1:5],
            order(fgseaResTidy.t1$NES, decreasing = T)[1:5])
  
  if(j == 1) fgseaResTidy.t1$pair <- 'LP vs EP'
  if(j == 2) fgseaResTidy.t1$pair <- 'LP vs BL'
  
  df_hallmark <- rbind(df_hallmark, fgseaResTidy.t1[Indx, ])
  
}

df_hallmark <- as.data.frame(df_hallmark)
df_hallmark <- df_hallmark[df_hallmark$pathway != 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', ]
df_hallmark$pathway <- as.factor(unlist(lapply(as.character(df_hallmark$pathway), function(x) strsplit(x, split = 'HALLMARK_')[[1]][2])))



pt2h <- ggplot(df_hallmark, aes(x = pair, y = reorder(pathway, NES))) +        ## global aes
  geom_point(aes(fill = NES, size =  -log(pval)),
             color = 'black',
             shape = 21,
             stroke = 0.01)  +  
  ggtitle('Hallmark Pathway') + 
  # scale_x_discrete(labels=c('EP vs BL', 'LP vs BL', 'LP vs EP')) + 
  xlab("") + ylab("") +
  labs(size='-log(P-values)') +
  # scale_fill_gradient(high = 'dark orange', low = 'dark blue') + 
  scale_fill_gradientn(
    ## limits = c(-1.5, 1.5),
    colors = c("#5DBCFF", "#6DCCFF", "white", "#F494BE", "#F484AE"))+
  # scale_size(range = c(0, 3.2), limits = c(0, 100), breaks = c(0,20,40,60,80,100))+             ## to tune the size of circles
  theme(
    text = element_text(size = 10), 
    plot.title = element_text(hjust = 1),
    panel.grid.major = element_line(colour = "grey90", size=0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1, size = 10),
    # axis.text.y=element_text(angle = -30, vjust = 1, hjust = 1),
    legend.position="none",
    legend.title=element_text(size=10)) +
  guides(size = guide_legend(title.position="top",
                             title.hjust = 0.5,
                             ncol = 1,
                             byrow = T,
                             override.aes = list(stroke = 0.4)),
         fill = guide_colourbar(title.position = "top", title.hjust = 0.5))


pdf(paste0('Figure2/Figure2h.pdf'), width = 3.8, height = 5)
print(pt2h)
dev.off()




##-------------------------------
## Figure 2I: Volcano plot for 
## liver BL vs LP
##-------------------------------

Idents(merged_combined) <- merged_combined$metastatic
liver_combined <- subset(merged_combined, idents = 'Liver')

Idents(liver_combined) <- liver_combined$Status_new


if(file.exists(paste0('Figure2/Data/DE_Liver_LP_BL.csv'))){
  thisDE <- read.csv(paste0('Figure2/Data/DE_Liver_LP_BL.csv'))
}else{
  thisDE <- FindMarkers(object = liver_combined, 
                        ident.1 = 'Late Progressor', ident.2 = 'Baseline', 
                        min.pct = 0, test.use = "MAST")
  thisDE$gene <- rownames(thisDE)
  write.csv(thisDE, file = paste0('Figure2/Data/DE_Liver_LP_BL.csv'))
}

# thisDE$avg_log2FC <- thisDE$avg_log2FC + rnorm(length(thisDE$avg_log2FC ), 0, 0.1)

Figure2i <- EnhancedVolcano(thisDE,
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
                       selectLab = c('S100A6',
                                     'EIF3H',
                                     'LDHA',
                                     'PPA1',
                                     'ENO1',
                                     'RPL4',
                                     'RPL5',
                                     'RPS2',
                                     'RPS13',
                                     'TUBA1A',
                                     'TUBB',
                                     'MGP', 'MDK', 'MTRNR2LB', 'PAK1', 
                                     'MALAT1', 'TPT1', 'APOE', 'APOC1', 
                                     'APOA1', 'IGHA1', 'IGHG1', 'FTL',
                                     'IGHM'),
                       # colors = c('#FD9AA0', '#6DCCFD'),
                       drawConnectors = TRUE,
                       legendPosition = 'right',
                       # ylim = c(0, 330),
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


pdf(paste0('Figure2/Figure2i.pdf'), width = 5, height = 4)
print(Figure2i)
dev.off()




##----------------------------
## Hallmark pathway and GOBP 
## pathway analysis
##----------------------------
## Hallmark pathway
pathways_list <- c('hallmark', 'gobp')
pathways_path <- c('Data/msigdb_v2022.1.Hs_GMTs/h.all.v2022.1.Hs.symbols.gmt')

## Load pathways gmt file
pathways <- gmtPathways(gmt.file = pathways_path[1])
## Run fGSEA
TESTX<-thisDE
TESTX$gene<-TESTX$gene
T1.mkers<-TESTX
res.t1 <- T1.mkers %>%
  dplyr::filter(p_val_adj < 0.01) %>%
  dplyr::select(gene, avg_log2FC) %>%
  tibble::deframe()
fgseaRes.t1 <- fgsea(pathways=pathways, stats=res.t1, nperm=1000)
fgseaResTidy.t1 <- fgseaRes.t1 %>%
  as_tibble() %>%
  dplyr::arrange(desc(NES))
Indx <- c(order(fgseaResTidy.t1$NES)[1:5],
          order(fgseaResTidy.t1$NES, decreasing = T)[1:5])

df_hallmark <- as.data.frame(fgseaResTidy.t1[Indx, ])
df_hallmark$pathway <- as.factor(unlist(lapply(as.character(df_hallmark$pathway), function(x) strsplit(x, split = 'HALLMARK_')[[1]][2])))
df_hallmark$pair <- factor('LP vs BL')

pt2j <- ggplot(df_hallmark, aes(x = pair, y = reorder(pathway, NES))) +        ## global aes
  geom_point(aes(fill = NES, size =  -log(pval)),
             color = 'black',
             shape = 21,
             stroke = 0.01)  +  
  ggtitle('Hallmark Pathway') + 
  # scale_x_discrete(labels=c('EP vs BL', 'LP vs BL', 'LP vs EP')) + 
  xlab("") + ylab("") +
  labs(size='-log(P-values)') +
  # scale_fill_gradient(high = 'dark orange', low = 'dark blue') + 
  scale_fill_gradientn(
    ## limits = c(-1.5, 1.5),
    colors = c("#5DBCFF", "#6DCCFF", "white", "#F494BE", "#F484AE"))+
  # scale_size(range = c(0, 3.2), limits = c(0, 100), breaks = c(0,20,40,60,80,100))+             ## to tune the size of circles
  theme(
    text = element_text(size = 10), 
    plot.title = element_text(hjust = 1),
    panel.grid.major = element_line(colour = "grey90", size=0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1, size = 10),
    # axis.text.y=element_text(angle = -30, vjust = 1, hjust = 1),
    legend.position="none",
    legend.title=element_text(size=10)) +
  guides(size = guide_legend(title.position="top",
                             title.hjust = 0.5,
                             ncol = 6,
                             byrow = T,
                             override.aes = list(stroke = 0.4)),
         fill = guide_colourbar(title.position = "top", title.hjust = 0.5))


pdf(paste0('Figure2/Figure2j.pdf'), width = 3, height = 5)
print(pt2j)
dev.off()




