##-------------------------------
## Figure Panel 6
## Author: Peng Yang
##-------------------------------

library(ggpubr)
library(Seurat)
library(SeuratData)
library(ggalluvial)
## load harmony integration
load('/Data/merged_data_Harmony_new.RData')
## Subset the samples
Idents(merged_combined) <- merged_combined$orig.ident
merged_combined <- subset(merged_combined, 
                          idents = c('PA3', 'PA3 #1', 'PA3 #2'))


##-------------------------------
## Figure 7A: Umap plot across 
## different samples and cell
## types
##-------------------------------

Idents(merged_combined) <- merged_combined$Tumor_cells

pt7b1 <- DimPlot(merged_combined, reduction = 'umap.pca', #group.by = 'CT_L2.num', 
                 # split.by = 'Status_new', 
                 label = F, repel = T, label.size = 5) + 
  NoAxes() + ggtitle(' ') + 
  scale_x_reverse() +
  theme_void() +
  scale_color_manual(labels = c('Tumor Cell', 'Non-tumor Cell'),
                     # values = c('#FD9AA0', '#6DCCFD')) +
                     values = c('#FD9AA0', '#6DCCFD')) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 13), 
                              title.position = 'bottom', 
                              nrow = 1,
                              label.hjust = 0,
                              label.theme = element_text())) + #NoLegend() +
  theme(strip.text.x = element_blank(), legend.position = 'bottom',
        legend.text = element_text(size = 13))

pt7b1

pdf('Figure6/Figure6b1.pdf', width = 5, height = 6)
print(pt7b1)
dev.off()


pt7b2 <- DimPlot(merged_combined, reduction = 'umap.harmony', #group.by = 'CT_L2.num', 
                 # split.by = 'Status_new', 
                 label = F, repel = T, label.size = 5) + 
  NoAxes() + ggtitle(' ') + 
  scale_x_reverse() +
  theme_void() + 
  scale_color_manual(labels = c('Tumor Cell', 'Non-tumor Cell'),
                     # values = c('#FD9AA0', '#6DCCFD')) +
                     values = c('#FD9AA0', '#6DCCFD')) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 6), 
                              title.position = 'bottom', 
                              nrow = 1,
                              label.hjust = 0,
                              label.theme = element_text())) + NoLegend() +
  theme(strip.text.x = element_blank())

pt7b2

pdf('Figure6/Figure6b2.pdf', width = 5, height = 5)
print(pt7b2)
dev.off()


Idents(merged_combined) <- merged_combined$orig.ident

pt7b3 <- DimPlot(merged_combined, reduction = 'umap.pca', #group.by = 'CT_L2.num', 
                 # split.by = 'Status_new', 
                 label = F, repel = T, label.size = 5) + 
  NoAxes() + ggtitle(' ') + 
  scale_x_reverse() +
  theme_void() + 
  scale_color_manual(# breaks = unique(merged_combined$orig.ident.num),
    # labels = paste(1:length(unique(combined_Tcells$CT_L2)), T_cells),
    values = c("#6DCCDD", "#EDCAE0", "#F9B26C")) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 13), 
                              title.position = 'bottom', 
                              nrow = 1,
                              label.hjust = 0,
                              label.theme = element_text())) + #NoLegend() +
  theme(strip.text.x = element_blank(), legend.position = 'bottom',
        legend.text = element_text(size = 13))


pdf('Figure6/Figure7b3.pdf', width = 5, height = 6)
print(pt7b3)
dev.off()


pt7b4 <- DimPlot(merged_combined, reduction = 'umap.harmony', #group.by = 'CT_L2.num', 
                 # split.by = 'Status_new', 
                 label = F, repel = T, label.size = 5) + 
  NoAxes() + ggtitle(' ') + 
  scale_x_reverse() +
  theme_void() + 
  scale_color_manual(# breaks = unique(merged_combined$orig.ident.num),
    # labels = paste(1:length(unique(combined_Tcells$CT_L2)), T_cells),
    values = c("#6DCCDD", "#EDCAE0", "#F9B26C")) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 6), 
                              title.position = 'bottom', 
                              nrow = 1,
                              label.hjust = 0,
                              label.theme = element_text())) + NoLegend() +
  theme(strip.text.x = element_blank())

pdf('Figure6/Figure6b4.pdf', width = 5, height = 5)
print(pt7b4)
dev.off()


##-------------------------------
## Figure 6c: Dotplot for DE
## Analysis
##-------------------------------
Idents(merged_combined) <- merged_combined$Tumor_cells
merged_combined_tumor <- subset(merged_combined, idents = 'tumor cells')

Idents(merged_combined) <- merged_combined$orig.ident

if(file.exists('Figure6/Data/DE_integrated_sample_marker.csv')){
  load('Figure6/Data/DE_integrated_sample_marker.RData')
}else{
  merged_combined.maker <- FindAllMarkers(merged_combined, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
  merged_combined.maker %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
  head(merged_combined.maker)
  save(merged_combined.maker, file = "Figure7/Data/DE_integrated_sample_marker.RData")
}



## Make the heatmap for marker genes
merged_combined.maker %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10


Figure6c <- DotPlot(merged_combined, features = top10$gene) + 
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
                             title.hjust = 0.5,
                             ncol = 4,
                             byrow = T,
                             override.aes = list(stroke = 0.4, ncol = 3)),
         fill = guide_colourbar(title.position = "top", title.hjust = 0.5))

pdf(paste0('Figure6/Figure7c.pdf'), width = 8, height = 3)
print(Figure6c)
dev.off()


##-------------------------------
## Figure 7D: Tumor cell dist
##-------------------------------
df_tumor <- c()
for(sample in unique(merged_combined$orig.ident)){
  df_tumor <- rbind(df_tumor, c(sample, 
                                sum(merged_combined$Tumor_cells[merged_combined$orig.ident == sample] == 'tumor cells'), 'Tumor cells'))
  df_tumor <- rbind(df_tumor, c(sample, 
                                sum(merged_combined$Tumor_cells[merged_combined$orig.ident == sample] != 'tumor cells'), 'Non-Tumor cells'))
}
colnames(df_tumor) <- c('Sample', 'Count', 'Condition')
df_tumor <- as.data.frame(df_tumor)


## Adding status informations
Status <- rep('Baseline', nrow(df_tumor))
Status[which(df_tumor$Sample %in% c('PA3 #1', 'PA3 #2')) ] = 'Late Progressor'
df_tumor$Status <- as.factor(Status)
df_tumor$Sample <- as.factor(df_tumor$Sample)
df_tumor$Count <- as.numeric(df_tumor$Count)
df_tumor$Condition <- factor(df_tumor$Condition, levels = c('Tumor cells', 'Non-Tumor cells'))


pt7d <- ggplot(df_tumor, aes(fill = Condition, y = Count, x = Sample)) +
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

pdf('Figure6/Figure6d.pdf', width = 4, height = 3)
print(pt7d)
dev.off()


##-------------------------------
## Figure 7E: Tumor hallmark
##-------------------------------
# df_hallmark <- read.csv(paste0('Figure7/Data/Tumor_pathway_hallmark_Baseline_Progression.csv'))
comparision <- c('PA3 #1', 'PA3 #2')
pathways_path <- c('/Data/msigdb_v2022.1.Hs_GMTs/h.all.v2022.1.Hs.symbols.gmt')
df_hallmark <- c()
for(i in 1:2){
  thisDE <- FindMarkers(object = merged_combined_tumor, 
                        ident.1 = comparision[i], ident.2 = "PA3", 
                        min.pct = 0, test.use = "MAST")
  thisDE$gene <- rownames(thisDE)
  
  pathways <- gmtPathways(gmt.file = pathways_path)
  
  TESTX<-thisDE
  # TESTX$gene<-TESTX$X
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
  
  df_hallmark <- rbind(df_hallmark, data.frame(fgseaResTidy.t1[Indx, 1:7]))
  
}

df_hallmark$pathway <- as.factor(unlist(lapply(as.character(df_hallmark$pathway), function(x) strsplit(x, split = 'HALLMARK_')[[1]][2])))
df_hallmark$pair <- factor(rep(c('PA3 #1 vs PA3', 'PA3 #2 vs PA3'), each = 10))

pt6e <- ggplot(df_hallmark, aes(x = pair, y = reorder(pathway, NES))) +        ## global aes
  geom_point(aes(fill = NES, size =  -log(pval)),
             color = 'black',
             shape = 21,
             stroke = 0.01)  +  
  ggtitle('Hallmark Pathway') + 
  # scale_x_discrete(labels=c('LP vs BL')) + 
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
    legend.position="right",
    legend.title=element_text(size=10)) +
  guides(size = guide_legend(title.position="top",
                             title.hjust = 0.5,
                             ncol = 1,
                             byrow = T,
                             override.aes = list(stroke = 0.4)),
         fill = guide_colourbar(title.position = "top", title.hjust = 0.5))


pdf(paste0('Figure6/Figure6e.pdf'), width = 4.6, height = 5)
print(pt6e)
dev.off()

##-------------------------------
## Figure 6f: Ligand and receptor
## interaction analysis
##-------------------------------

load('CellChat/cell_chat_matched_LP.RData')
pdf(paste0('Figures/Figure6f.pdf'), width = 7, height = 6.5)
print(netVisual_bubble(cellchat, remove.isolate = FALSE, title.name = s))
dev.off()


##-------------------------------
## Figure 6g: Umap plot of non
## tumor cells
##-------------------------------
Idents(combined_nontumor) <- combined_nontumor$CT_L1
combined_nontumor$CT_L1.num <- as.factor(as.numeric(as.factor(combined_nontumor$CT_L1)))

pt7g <- DimPlot(combined_nontumor, reduction = 'umap.harmony', group.by = 'CT_L1.num', 
                 # split.by = 'Status_new', 
                 label = T, repel = T, label.size = 5) + 
  NoAxes() + ggtitle(' ') + 
  scale_x_reverse() +
  theme_void() + 
  scale_color_manual(# breaks = unique(merged_combined$orig.ident.num),
    labels = paste(1:length(unique(combined_nontumor$CT_L1)), cell_types),
    values = cell_colors) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 6), 
                              title.position = 'bottom', 
                              nrow = 1,
                              label.hjust = 0,
                              label.theme = element_text())) + NoLegend() +
  theme(strip.text.x = element_blank())


pdf('Figure6/Figure6g.pdf', width = 5, height = 5)
print(pt7g)
dev.off()


##-------------------------------
## Figure 6h: Non-tumor cell
## distribution
##-------------------------------

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
df_nontumor_L1       <- df_nontumor_L1[df_nontumor_L1$Count != 0, ]
df_nontumor_L1$Sample <- as.factor(df_nontumor_L1$Sample)
df_nontumor_L1$Status <- as.factor(df_nontumor_L1$Status)
df_nontumor_L1$CellTypes <- factor(df_nontumor_L1$CellTypes, levels = cell_types)


df_nontumor_L1_frac <- df_nontumor_L1[df_nontumor_L1$CellTypes %in% cell_types, ]
for(i in unique(df_nontumor_L1_frac$Sample)){
  
  Indx <- which(df_nontumor_L1_frac$Sample == i)
  
  df_nontumor_L1_frac[Indx, 'Count'] = df_nontumor_L1_frac[Indx, 'Count'] / sum(df_nontumor_L1_frac[Indx, 'Count'])
  
}


## Pluvial plot
df_nontumor_L1_frac$subject <- as.numeric(df_nontumor_L1_frac$CellTypes)

pt6h <- ggplot(df_nontumor_L1_frac,
               aes(x = Sample, stratum = CellTypes, alluvium = subject,
                   y = Count,
                   fill = CellTypes, label = CellTypes)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow() + ylab("Fractions") + xlab(" ") + 
  geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 3, min.y = 0.02) + theme_void() +
  theme(legend.position = "none",
        axis.title.x = element_text(),
        axis.text.x = element_text())  +
  ggtitle(" ") + #scale_fill_brewer(palette="PiYG") + 
  scale_fill_manual(values = cell_colors)

pdf('Figure6/Figure6h.pdf', width = 3.5, height = 4.5)
print(pt6h)
dev.off()


##-------------------------------
## Figure 7I: MP Analysis
##-------------------------------

## Subsetting CD8T cells
CD8T_cells <- c('CD8 Naive', 'CD8 Proliferating', 'CD8 TCM', 'CD8 TEM', 
                'Treg', 'gdT')
Idents(merged_combined) <- merged_combined$CT_L2
combined_CD8T_cells <- subset(merged_combined, idents = CD8T_cells)

Idents(combined_CD8T_cells) <- combined_CD8T_cells$Status_new

## Subsetting CD4T cells
CD4T_cells <- c('CD4 CTL', 'CD4 Naive', 'CD4 Proliferating', 'CD4 TCM', 'CD4 TEM',
                'Treg', 'gdT')
Idents(merged_combined) <- merged_combined$CT_L2
combined_CD4T_cells <- subset(merged_combined, idents = CD4T_cells)

Idents(combined_CD4T_cells) <- combined_CD4T_cells$Status_new

## For CD8T cells
MP_meta <- readxl::read_excel('/Tcell_function/41591_2023_2371_MOESM3_ESM.xlsx', sheet = 5, col_names = T, skip = 1)
MP_meta <- as.data.frame(MP_meta)
org_MP_Names <- names(MP_meta)
names(MP_meta) = paste0('MP', 1:ncol(MP_meta))

combined_CD8T_cells <- AddModuleScore(
  object = combined_CD8T_cells,
  features = as.list(MP_meta),
  ctrl = 5,
  name = 'CD8T'
)

res_MP_CD8T <- array(0, c(ncol(MP_meta)*2, 4))

for(j in 1:ncol(MP_meta)){
  
  ## meta data
  res_MP_CD8T[j,1] <- 'PA3 #1 vs PA3' #paste(pairs[i,], collapse = '_vs_')
  res_MP_CD8T[j,2] <- org_MP_Names[j]
  
  ## Calculate the median difference
  res_MP_CD8T[j,3] <- median(combined_CD8T_cells[[paste0('CD8T', j)]][combined_CD8T_cells$orig.ident == 'PA3 #1',1]) - 
    median(combined_CD8T_cells[[paste0('CD8T', j)]][combined_CD8T_cells$orig.ident == 'PA3',1])
  
  res_MP_CD8T[j,4] <- t.test(combined_CD8T_cells[[paste0('CD8T', j)]][combined_CD8T_cells$orig.ident == 'PA3 #1',1],
                             combined_CD8T_cells[[paste0('CD8T', j)]][combined_CD8T_cells$orig.ident == 'PA3',1])$p.value
  
  i = j + ncol(MP_meta)
  
  ## meta data
  res_MP_CD8T[i,1] <- 'PA3 #2 vs PA3' #paste(pairs[i,], collapse = '_vs_')
  res_MP_CD8T[i,2] <- org_MP_Names[j]
  
  ## Calculate the median difference
  res_MP_CD8T[i,3] <- median(combined_CD8T_cells[[paste0('CD8T', j)]][combined_CD8T_cells$orig.ident == 'PA3 #2',1]) - 
    median(combined_CD8T_cells[[paste0('CD8T', j)]][combined_CD8T_cells$orig.ident == 'PA3',1])
  
  res_MP_CD8T[i,4] <- t.test(combined_CD8T_cells[[paste0('CD8T', j)]][combined_CD8T_cells$orig.ident == 'PA3 #2',1],
                             combined_CD8T_cells[[paste0('CD8T', j)]][combined_CD8T_cells$orig.ident == 'PA3',1])$p.value
}

res_MP_CD8T_df <- data.frame(comparison = as.factor(res_MP_CD8T[,1]),
                             MP = as.factor(res_MP_CD8T[,2]),
                             changes = as.numeric(res_MP_CD8T[,3]),
                             p.value = as.numeric(res_MP_CD8T[,4]),
                             celltypes = 'CD8T')

## For CD4T cells
MP_meta <- readxl::read_excel('/CDK4/Tcell_function/41591_2023_2371_MOESM3_ESM.xlsx', sheet = 7, col_names = T, skip = 1)
MP_meta <- as.data.frame(MP_meta)
org_MP_Names <- names(MP_meta)
names(MP_meta) = paste0('MP', 1:ncol(MP_meta))

combined_CD4T_cells <- AddModuleScore(
  object = combined_CD4T_cells,
  features = as.list(MP_meta),
  ctrl = 5,
  name = 'CD4T'
)

res_MP_CD4T <- array(0, c(ncol(MP_meta)*2, 4))

for(j in 1:ncol(MP_meta)){
  
  ## meta data
  res_MP_CD4T[j,1] <- 'PA3 #1 vs PA3' #paste(pairs[i,], collapse = '_vs_')
  res_MP_CD4T[j,2] <- org_MP_Names[j]
  
  ## Calculate the median difference
  res_MP_CD4T[j,3] <- median(combined_CD4T_cells[[paste0('CD4T', j)]][combined_CD4T_cells$orig.ident == 'PA3 #1',1]) - 
    median(combined_CD4T_cells[[paste0('CD4T', j)]][combined_CD4T_cells$orig.ident == 'PA3',1])
  
  res_MP_CD4T[j,4] <- t.test(combined_CD4T_cells[[paste0('CD4T', j)]][combined_CD4T_cells$orig.ident == 'PA3 #1',1],
                             combined_CD4T_cells[[paste0('CD4T', j)]][combined_CD4T_cells$orig.ident == 'PA3',1])$p.value
  
  i = j + ncol(MP_meta)
  
  ## meta data
  res_MP_CD4T[i,1] <- 'PA3 #2 vs PA3' #paste(pairs[i,], collapse = '_vs_')
  res_MP_CD4T[i,2] <- org_MP_Names[j]
  
  ## Calculate the median difference
  res_MP_CD4T[i,3] <- median(combined_CD4T_cells[[paste0('CD4T', j)]][combined_CD4T_cells$orig.ident == 'PA3 #2',1]) - 
    median(combined_CD4T_cells[[paste0('CD4T', j)]][combined_CD4T_cells$orig.ident == 'PA3',1])
  
  res_MP_CD4T[i,4] <- t.test(combined_CD4T_cells[[paste0('CD4T', j)]][combined_CD4T_cells$orig.ident == 'PA3 #2',1],
                             combined_CD4T_cells[[paste0('CD4T', j)]][combined_CD4T_cells$orig.ident == 'PA3',1])$p.value
}

res_MP_CD4T_df <- data.frame(comparison = as.factor(res_MP_CD4T[,1]),
                             MP = as.factor(res_MP_CD4T[,2]),
                             changes = as.numeric(res_MP_CD4T[,3]),
                             p.value = as.numeric(res_MP_CD4T[,4]),
                             celltypes = 'CD4T')

## Merge the CD4T and CD8T cells
res_MP_df <- rbind(res_MP_CD8T_df, res_MP_CD4T_df)
res_MP_df$celltypes <- as.factor(res_MP_df$celltypes)


pt6i <- ggplot(res_MP_df, aes(comparison, reorder(MP, changes), fill = changes)) +
  geom_tile() + facet_wrap(~celltypes) + 
  ggtitle('MP Analysis') +
  # scale_fill_viridis(discrete = FALSE) +
  scale_fill_gradientn(
    limits = c(-.25, .25),
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
  )

pdf(paste0('Figure6/Figure6i.pdf'), width = 6, height = 6)
print(pt6i)
dev.off()


