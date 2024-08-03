##------------------------------------
## Cellchat analysis 
## Author: Peng Yang
##-----------------------------------
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
library(CellChat)
library(patchwork)
##------------------------------------------
## Load the data
##------------------------------------------


Status <- c('Baseline', 'Early Progressor', 'Late Progressor')
df_all <- c()
for(s in Status){
  
  if(file.exists(paste0('CellChat/cell_chat_',s, '.RData'))){
    load(paste0('CellChat/cell_chat_',s, '.RData'))
  }else{
    
    load('/Data/combined_cellchat.RData')
    
    Idents(merged_combined) <- merged_combined$Status_new
    merged_combined <- subset(merged_combined, ident = s)
    
    # This is a combined data from two biological conditions: normal and diseases
    data.input = normalizeData(merged_combined@assays$RNA@counts) # normalized data matrix
    meta = data.frame(labels = as.character(merged_combined$labels),
                      status = merged_combined$Status_new,
                      id = merged_combined$orig.ident) # a dataframe with rownames containing cell mata data
    # rownames(meta) <- colnames(data.input)
    # 
    # cell.use = rownames(meta)[meta$status == s] # extract the cell names from disease data
    # 
    # # Prepare input data for CelChat analysis
    # data.input = data.input[, cell.use]
    # meta = meta[cell.use, ]
    # # meta = data.frame(labels = meta$labels[cell.use], row.names = colnames(data.input)) # manually create a dataframe consisting of the cell labels
    unique(meta$labels) # check the cell labels
    levels(meta$labels) <-   unique(meta$labels)
    
    rm(merged_combined)
    
    ##------------------------------------------
    ## Create cellchat object
    ##------------------------------------------
    cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
    
    ##------------------------------------------
    ## Set the ligand-receptor interaction database
    ##------------------------------------------
    CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
    # showDatabaseCategory(CellChatDB)
    
    
    # use a subset of CellChatDB for cell-cell communication analysis
    CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
    # use all CellChatDB for cell-cell communication analysis
    # CellChatDB.use <- CellChatDB # simply use the default CellChatDB
    
    # set the used database in the object
    cellchat@DB <- CellChatDB.use
    ##------------------------------------------
    ## Preprocessing the expression data for 
    ## cell-cell communication analysis
    ##------------------------------------------
    
    # subset the expression data of signaling genes for saving computation cost
    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
    future::plan("multisession", workers = 4) # do parallel
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    
    
    ##------------------------------------------
    ## Compute the communication probability and
    ## infer cellular communication network
    ##------------------------------------------
    cellchat <- computeCommunProb(cellchat, type = "triMean")
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    
    ##------------------------------------------
    ## Infer the cell-cell communication at a 
    ## signaling pathway level
    ##------------------------------------------
    
    cellchat <- computeCommunProbPathway(cellchat)
    
    ##------------------------------------------
    ## Calculate the aggregated cell-cell 
    ## communication network
    ##------------------------------------------
    cellchat <- aggregateNet(cellchat)
    
    save(cellchat, file = paste0('CellChat/cell_chat_', s, '.RData'))
  }
  ##------------------------------------------
  ## Figures
  ##------------------------------------------
  pathways.show.all <- cellchat@netP$pathways
  
  pt_bubble <- netVisual_bubble(cellchat, remove.isolate = FALSE, title.name = s,
                                return.data = T)
  
  df_all <- rbind(df_all, cbind(pt_bubble$communication, s))

  # pdf(paste0('CellChat/cellchat_bubble_', s, '.pdf'), width = 7, height = 6.5)
  # print(pt_bubble$gg.obj)
  # dev.off()
  # 
  # pdf(paste0('CellChat/cellchat_chord_', s, '.pdf'), width = 8, height = 8)
  # print(netVisual_chord_gene(cellchat, sources.use = NULL, targets.use = NULL, signaling = pathways.show.all, legend.pos.x = 8,
  #                            scale = T, thresh = 0.05,
  #                            title.name = s))
  # dev.off()
  
}

names(df_all)[14] = 'Status'
df_all$Status <- factor(df_all$Status, levels = c('Baseline', 'Early Progressor', 'Late Progressor'))

color.heatmap = c("Spectral", 
                  "viridis"); 
n.colors = 10; direction = -1; thresh = 0.05; 

comparison = NULL; group = NULL; remove.isolate = FALSE; 

max.dataset = NULL; min.dataset = NULL; min.quantile = 0; 

max.quantile = 1; line.on = TRUE; line.size = 0.2; color.text.use = TRUE; 

color.text = NULL; dot.size.min = NULL; dot.size.max = NULL; 

title.name = NULL; font.size = 10; font.size.title = 10; 

show.legend = TRUE; grid.on = TRUE; color.grid = "grey90"; 

angle.x = 90; vjust.x = NULL; hjust.x = NULL; 
angle.x = 90; vjust.x = NULL; hjust.x = NULL; 
values <- c(1, 2, 3)
names(values) <- c("p > 0.05", "0.01 < p < 0.05", "p < 0.01")
dot.size.max = max(df_all$pval)
dot.size.min = min(df_all$pval)
g <- ggplot(df_all, aes(x = source.target, y = interaction_name_2, 
               color = prob, size = pval)) + geom_point(pch = 16) + 
  theme_linedraw() + theme(panel.grid.major = element_blank()) +
  facet_grid(~Status) + 
  scale_radius(range = c(dot.size.min, dot.size.max), 
                 breaks = sort(unique(df_all$pval)), labels = names(values)[values %in% 
                                                                          sort(unique(df_all$pval))], name = "p-value") + 
  theme(axis.text.x = element_text(angle = angle.x, hjust = hjust.x,
                                   vjust = vjust.x), axis.title.x = element_blank(),
        axis.title.y = element_blank()) + scale_x_discrete(position = "bottom")

color.heatmap = c("Spectral", 
                  "viridis")
color.use <- RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap[1])
color.use <- rev(color.use)
g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), 
                                na.value = "white", limits = c(quantile(df_all$prob, 
                                                                        0, na.rm = T), quantile(df_all$prob, 1, na.rm = T)), 
                                breaks = c(quantile(df_all$prob, 0, na.rm = T), quantile(df_all$prob, 
                                                                                         1, na.rm = T)), labels = c("min", "max")) + guides(color = guide_colourbar(barwidth = 0.5, 
                                                                                                                                                                    title = "Commun. Prob."))
g <- g + theme(text = element_text(size = font.size), plot.title = element_text(size = font.size.title)) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 6))
if (grid.on) {
  if (length(unique(df_all$source.target)) > 1) {
    g <- g + geom_vline(xintercept = seq(1.5, length(unique(df_all$source.target)) - 
                                           0.5, 1), lwd = 0.1, colour = color.grid)
  }
  if (length(unique(df_all$interaction_name_2)) > 1) {
    g <- g + geom_hline(yintercept = seq(1.5, length(unique(df_all$interaction_name_2)) - 
                                           0.5, 1), lwd = 0.1, colour = color.grid)
  }
}

g <- g + theme(text = element_text(size = font.size), 
               plot.title = element_text(size = font.size.title)) + 
  theme(legend.title = element_text(size = 18), legend.text = element_text(size = 18)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA),
        strip.text = element_text(colour = 'black', size = 18))

pdf(paste0('CellChat/cellchat_bubble_all.pdf'), width = 20, height = 7)
print(g)
dev.off()
# dim(cellchat@)



## For the matched sample

Status <- c('Baseline', 'Late Progressor')

for(s in Status){
  
  if(file.exists(paste0('CellChat/cell_chat_matched_',s, '.RData'))){
    load(paste0('CellChat/cell_chat_matched_',s, '.RData'))
  }else{
    
    load('/Users/yangpeng/Box Sync/CDK4/Data/combined_cellchat.RData')
    
    Idents(merged_combined) <- merged_combined$orig.ident
    if(s == 'Baseline'){
      merged_combined <- subset(merged_combined, ident = 'PA3')
    }else{
      merged_combined <- subset(merged_combined, ident = c('PA3 #1', 'PA3 #2'))
    }
    
    
    # This is a combined data from two biological conditions: normal and diseases
    data.input = normalizeData(merged_combined@assays$RNA@counts) # normalized data matrix
    meta = data.frame(labels = as.character(merged_combined$labels),
                      status = merged_combined$Status_new,
                      id = merged_combined$orig.ident) # a dataframe with rownames containing cell mata data
    # rownames(meta) <- colnames(data.input)
    # 
    # cell.use = rownames(meta)[meta$status == s] # extract the cell names from disease data
    # 
    # # Prepare input data for CelChat analysis
    # data.input = data.input[, cell.use]
    # meta = meta[cell.use, ]
    # # meta = data.frame(labels = meta$labels[cell.use], row.names = colnames(data.input)) # manually create a dataframe consisting of the cell labels
    unique(meta$labels) # check the cell labels
    levels(meta$labels) <-   unique(meta$labels)
    
    rm(merged_combined)
    
    ##------------------------------------------
    ## Create cellchat object
    ##------------------------------------------
    cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
    
    ##------------------------------------------
    ## Set the ligand-receptor interaction database
    ##------------------------------------------
    CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
    # showDatabaseCategory(CellChatDB)
    
    
    # use a subset of CellChatDB for cell-cell communication analysis
    CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
    # use all CellChatDB for cell-cell communication analysis
    # CellChatDB.use <- CellChatDB # simply use the default CellChatDB
    
    # set the used database in the object
    cellchat@DB <- CellChatDB.use
    ##------------------------------------------
    ## Preprocessing the expression data for 
    ## cell-cell communication analysis
    ##------------------------------------------
    
    # subset the expression data of signaling genes for saving computation cost
    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
    future::plan("multisession", workers = 4) # do parallel
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    
    
    ##------------------------------------------
    ## Compute the communication probability and
    ## infer cellular communication network
    ##------------------------------------------
    cellchat <- computeCommunProb(cellchat, type = "triMean")
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    
    ##------------------------------------------
    ## Infer the cell-cell communication at a 
    ## signaling pathway level
    ##------------------------------------------
    
    cellchat <- computeCommunProbPathway(cellchat)
    
    ##------------------------------------------
    ## Calculate the aggregated cell-cell 
    ## communication network
    ##------------------------------------------
    cellchat <- aggregateNet(cellchat)
    
    save(cellchat, file = paste0('CellChat/cell_chat_matched_',s, '.RData'))
  }
  ##------------------------------------------
  ## Figures
  ##------------------------------------------
  pathways.show.all <- cellchat@netP$pathways
  
  
  pdf(paste0('CellChat/cellchat_matched_bubble_', s, '.pdf'), width = 7, height = 6.5)
  print(netVisual_bubble(cellchat, remove.isolate = FALSE, title.name = s))
  dev.off()
  
  pdf(paste0('CellChat/cellchat_matched_chord_', s, '.pdf'), width = 8, height = 8)
  print(netVisual_chord_gene(cellchat, sources.use = NULL, targets.use = NULL, signaling = pathways.show.all, legend.pos.x = 8,
                             scale = T, thresh = 0.05,
                             title.name = s))
  dev.off()
  
}

