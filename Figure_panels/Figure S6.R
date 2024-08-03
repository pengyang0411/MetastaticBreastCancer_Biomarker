##----------------------------------
## Figure S6
## Author: Peng Yang
##----------------------------------

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
load('/Data/combined_nonTumor_final.RData')

## Subset the samples
Idents(combined_nontumor) <- combined_nontumor$orig.ident
combined_nontumor <- subset(combined_nontumor, 
                            idents = c('PA3', 'PA3 #1', 'PA3 #2'))

load('/Data/merged_data_Harmony_new.RData')
## Subset the samples
Idents(merged_combined) <- merged_combined$orig.ident
merged_combined <- subset(merged_combined, 
                          idents = c('PA3', 'PA3 #1', 'PA3 #2'))

##-------------------------------
## Figure S6A: MP Analysis
##-------------------------------


Idents(merged_combined) <- merged_combined$Tumor_cells
merged_combined_tumor <- subset(merged_combined, idents = 'tumor cells')

MP_meta <- readxl::read_excel('/Data/inferCNV/MP signatures.xlsx')
MP_meta <- as.data.frame(MP_meta)
org_MP_Names <- names(MP_meta)
names(MP_meta) = paste0('MP', 1:ncol(MP_meta))


merged_combined_tumor <- AddModuleScore(
  object = merged_combined,
  features = as.list(MP_meta),
  ctrl = 5,
  name = 'MP'
)

res_MP <- array(0, c(ncol(MP_meta)*2, 4))

for(j in 1:ncol(MP_meta)){
  
  ## meta data
  res_MP[j,1] <- 'PA3 #1 vs PA3' #paste(pairs[i,], collapse = '_vs_')
  res_MP[j,2] <- org_MP_Names[j]
  
  ## Calculate the median difference
  res_MP[j,3] <- median(merged_combined_tumor[[paste0('MP', j)]][merged_combined_tumor$orig.ident == 'PA3 #1',1]) - 
    median(merged_combined_tumor[[paste0('MP', j)]][merged_combined_tumor$orig.ident == 'PA3',1])
  
  res_MP[j,4] <- t.test(merged_combined_tumor[[paste0('MP', j)]][merged_combined_tumor$orig.ident == 'PA3 #1',1],
                        merged_combined_tumor[[paste0('MP', j)]][merged_combined_tumor$orig.ident == 'PA3',1])$p.value
  
  i = j + ncol(MP_meta)
  
  ## meta data
  res_MP[i,1] <- 'PA3 #2 vs PA3' #paste(pairs[i,], collapse = '_vs_')
  res_MP[i,2] <- org_MP_Names[j]
  
  ## Calculate the median difference
  res_MP[i,3] <- median(merged_combined_tumor[[paste0('MP', j)]][merged_combined_tumor$orig.ident == 'PA3 #2',1]) - 
    median(merged_combined_tumor[[paste0('MP', j)]][merged_combined_tumor$orig.ident == 'PA3',1])
  
  res_MP[i,4] <- t.test(merged_combined_tumor[[paste0('MP', j)]][merged_combined_tumor$orig.ident == 'PA3 #2',1],
                        merged_combined_tumor[[paste0('MP', j)]][merged_combined_tumor$orig.ident == 'PA3',1])$p.value
}

res_MP_df <- data.frame(comparison = as.factor(res_MP[,1]),
                        MP = as.factor(res_MP[,2]),
                        changes = as.numeric(res_MP[,3]),
                        p.value = as.numeric(res_MP[j,4]))
library(dplyr)
res_MP_df <- res_MP_df %>%
  as_tibble() %>%
  dplyr::arrange(desc(changes))

df_MP <- as.data.frame(res_MP_df)

# df_MP <- read.csv('Figure7/Data/Tumor_MP_Baseline_Progression.csv')
# df_MP$pair <- 'LP vs BL'

ptS6a <- ggplot(df_MP, aes(comparison, reorder(MP, changes), fill = changes)) +
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
  ) 
# scale_y_discrete(labels=rev(orders))



pdf(paste0('FigureS6/FigureS6a.pdf'), width = 5.2, height = 5)
print(ptS6a)
dev.off()



##--------------------------
## Figure S6B-D Ligand-receptor
## analysis
##--------------------------

for(s in 1:2){
  
  load(file = paste0('CellChat/cell_chat_matched_', c('Baseline', 'Late Progressor')[s], '.RData'))
  
  if(s == 1){
    pdf(paste0('FigureS6/FigureS6B.pdf'), width = 8, height = 8)
    print(netVisual_bubble(cellchat, remove.isolate = FALSE, title.name = s))
    dev.off()
  }
  
  pdf(paste0('FigureS6/FigureS6C-D', s, '.pdf'), width = 8, height = 8)
  print(netVisual_chord_gene(cellchat, sources.use = NULL, targets.use = NULL, signaling = pathways.show.all, legend.pos.x = 8,
                             scale = T, thresh = 0.05,
                             title.name = s))
  dev.off()
  
}


##--------------------------
## Figure S6E: Umap
##--------------------------

Idents(combined_nontumor) <- combined_nontumor$CT_L2
combined_nontumor$CT_L2.num <- as.factor(as.numeric(as.factor(combined_nontumor$CT_L2)))

ptS6e1 <- DimPlot(combined_nontumor, reduction = 'umap.pca', group.by = 'CT_L2.num', 
                  # split.by = 'Status_new', 
                  label = F, repel = T, label.size = 5) + 
  NoAxes() + ggtitle(' ') + 
  scale_x_reverse() +
  theme_void() + 
  # scale_color_manual(# breaks = unique(merged_combined$orig.ident.num),
  #   labels = paste(1:length(unique(combined_nontumor$CT_L1)), cell_types),
  #   values = cell_colors) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 6), 
                              title.position = 'bottom', 
                              nrow = 1,
                              label.hjust = 0,
                              label.theme = element_text())) + NoLegend() +
  theme(strip.text.x = element_blank())



pdf('FigureS6/FigureS6e1.pdf', width = 5, height = 5)
print(ptS6e1)
dev.off()


ptS6e2 <- DimPlot(combined_nontumor, reduction = 'umap.harmony', group.by = 'CT_L2', 
                  # split.by = 'Status_new', 
                  label = F, repel = T, label.size = 5) + 
  NoAxes() + ggtitle(' ') + 
  scale_x_reverse() +
  theme_void() + 
  # scale_color_manual(# breaks = unique(merged_combined$orig.ident.num),
  #   labels = paste(1:length(unique(combined_nontumor$CT_L1)), cell_types),
  #   values = cell_colors) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 6), 
                              title.position = 'bottom', 
                              # ncol = 1,
                              label.hjust = 0,
                              label.theme = element_text())) + #NoLegend() +
  theme(strip.text.x = element_blank(), legend.position = 'right')



pdf('FigureS6/FigureS6e2.pdf', width = 8, height = 5)
print(ptS6e1)
dev.off()

##--------------------------
## Figure S6F cell fraction
##--------------------------

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
df_nontumor_L2       <- df_nontumor_L2[df_nontumor_L2$Count != 0, ]
df_nontumor_L2$Sample <- as.factor(df_nontumor_L2$Sample)
df_nontumor_L2$Status <- as.factor(df_nontumor_L2$Status)
df_nontumor_L2$CellTypes <- factor(df_nontumor_L2$CellTypes, levels = cell_types)

df_nontumor_L2_frac <- df_nontumor_L2#[df_nontumor_L2$CellTypes %in% cells_2, ]
for(i in unique(df_nontumor_L2_frac$Sample)){
  
  Indx <- which(df_nontumor_L2_frac$Sample == i)
  
  df_nontumor_L2_frac[Indx, 'Count'] = df_nontumor_L2_frac[Indx, 'Count'] / sum(df_nontumor_L2_frac[Indx, 'Count'])
  
}
df_nontumor_L2_frac$subject <- as.numeric(df_nontumor_L2_frac$CellTypes)
ptS6f <- ggplot(df_nontumor_L2_frac,
                aes(x = Sample, stratum = CellTypes, alluvium = subject,
                    y = Count,
                    fill = CellTypes, label = CellTypes)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow() + ylab("Fractions") + xlab(" ") + 
  geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 3, min.y = 0.1) + theme_void() +
  theme(legend.position = "none",
        axis.title.x = element_text(),
        axis.text.x = element_text())  +
  ggtitle(" ") #+ scale_fill_brewer(palette="Set2")

pdf('FigureS6/FigureS6f.pdf', width = 4, height = 6)
print(ptS6f)
dev.off()


