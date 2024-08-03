##-------------------------------
## Figure Panel 1
## Author: Peng Yang
##-------------------------------

## Load the merged data with Harmony integration
load('Data/merged_data_Harmony_new.RData')

##-------------------------------
## Figure 1C: UMAP PCA across 
## samples
##-------------------------------

umapColor <- c("#6DCCDD","#9ACBDE","#C8CADF","#EDC6DD","#F092B1","#F27FA5",
               "#F47892", "#F6A395","#F8AD77","#E7B080","#C9AFA2","#ABADC5",
               "#AEB9AC","#B9C984", "#C4DA5D")
SampleNames <- c("PA3",   "PA5",   "PA45",  "PA95",  "PA120", "PA110", "PA131", "PA153",
                 'PA11', 'PA46', 'PA144',
                 'PA1', 'PA3 #1', 'PA3 #2',
                 'PA126', 'PA139', 'PA99', 'PA165')

merged_combined$orig.ident <- factor(merged_combined$orig.ident, 
                                     levels = SampleNames)


pt1c <- DimPlot(merged_combined, reduction = 'umap.pca', group.by = 'orig.ident.num', 
                label = T, repel = F, label.size = 5) + 
  NoAxes() + ggtitle(' ') + 
  scale_x_reverse() +
  theme_void() + 
  scale_color_manual(# breaks = unique(merged_combined$orig.ident.num),
    labels = paste(1:length(unique(merged_combined$orig.ident)), SampleNames), 
    values = colorRampPalette(umapColor)(length(unique(merged_combined$orig.ident))))  + NoLegend() 

pdf('Figure1/Figure1c.pdf', width = 8, height = 8)
print(pt1c)
dev.off()


##-------------------------------
## Figure 1D: Tumor cell fraction
##-------------------------------
## Create the dataframe
df_tumor <- c()
for(sample in unique(merged_combined$orig.ident)){
  df_tumor <- rbind(df_tumor, c(sample, as.character(merged_combined$Status_new[merged_combined$orig.ident == sample][1]),
                                sum(merged_combined$Tumor_cells[merged_combined$orig.ident == sample] == 'tumor cells'), 'Tumor Cell'))
  df_tumor <- rbind(df_tumor, c(sample, as.character(merged_combined$Status_new[merged_combined$orig.ident == sample][1]),
                                sum(merged_combined$Tumor_cells[merged_combined$orig.ident == sample] != 'tumor cells'), 'Non-Tumor Cell'))
}
colnames(df_tumor) <- c('Sample', 'Status', 'Count', 'Condition')
df_tumor <- as.data.frame(df_tumor)
df_tumor$Status <- as.factor(df_tumor$Status)
df_tumor$Sample <- as.factor(df_tumor$Sample)
df_tumor$Count <- as.numeric(df_tumor$Count)
df_tumor$Condition <- factor(df_tumor$Condition, levels = c('Tumor Cell', 'Non-Tumor Cell'))
## make the plot
pt1d <- ggplot(df_tumor, aes(fill = Condition, y = Count, x = Sample)) +
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

pdf('Figure1/Figure1d.pdf', width = 6, height = 2.5)
print(pt1d)
dev.off()

##-------------------------------
## Figure 1E: UMAP of tumor cells
## versus non-tumor cells
##-------------------------------

pt1e <- DimPlot(merged_combined, reduction = 'umap.pca', group.by = 'Tumor_cells', 
                # split.by = 'Tumor_cells', cols = 1, 
                label = F, repel = F) + NoAxes() + ggtitle(' ') +  #+ NoLegend()
  scale_x_reverse() + 
  guides(color = guide_legend(override.aes = list(size=8), ncol=2) ) +
  theme(legend.position = 'bottom') + 
  scale_color_manual(labels = c('Tumor Cell', 'Non-tumor Cell'),
                     # values = c('#FD9AA0', '#6DCCFD')) +
                     values = c('#FD9AA0', 'light grey')) +
  scale_alpha_continuous(range = c(1,0.2)) + 
  guides(color = guide_legend(override.aes = list(shape = 15, size = 6),
                              label.hjust = 0,
                              label.theme = element_text())) + NoLegend()

pdf('Figure1/Figure1e_tumor.pdf', width = 6, height = 6)
print(pt1e)
dev.off()

pt1e <- DimPlot(merged_combined, reduction = 'umap.pca', group.by = 'Tumor_cells', 
                # split.by = 'Tumor_cells', cols = 1, 
                label = F, repel = F) + NoAxes() + ggtitle(' ') +  #+ NoLegend()
  scale_x_reverse() + 
  guides(color = guide_legend(override.aes = list(size=8), ncol=2) ) +
  theme(legend.position = 'bottom') + 
  scale_color_manual(labels = c('Tumor Cell', 'Non-tumor Cell'),
                     # values = c('#FD9AA0', '#6DCCFD')) +
                     values = c('light grey', '#6DCCFD')) +
  scale_alpha_continuous(range = c(1,0.2)) + 
  guides(color = guide_legend(override.aes = list(shape = 15, size = 6),
                              label.hjust = 0,
                              label.theme = element_text())) + NoLegend()

pdf('Figure1/Figure1e_nontumor.pdf', width = 6, height = 6)
print(pt1e)
dev.off()


##-------------------------------
## Figure 1F: UMAP of samples
## in different status
##-------------------------------
colorForStatus <- c('#cccc66', '#66cccc', '#ff9966')

pt1f <- DimPlot(merged_combined, reduction = 'umap.pca', group.by = 'Status_new', 
                label = F, repel = TRUE) + NoAxes() + ggtitle(' ') + # + scale_color_brewer(palette="Set2") + 
  scale_color_manual(values = c(colorForStatus[1], 'light grey', 'light grey')) + 
  scale_x_reverse() + 
  guides(color = guide_legend(override.aes = list(size=8), ncol=3) )+
  theme(legend.position = 'bottom') + NoLegend()

pdf('Figure1/Figure1f_BL.pdf', width = 6, height = 6)
print(pt1f)
dev.off()

pt1f <- DimPlot(merged_combined, reduction = 'umap.pca', group.by = 'Status_new', 
                label = F, repel = TRUE) + NoAxes() + ggtitle(' ') + # + scale_color_brewer(palette="Set2") + 
  scale_color_manual(values = c('light grey', colorForStatus[2], 'light grey')) + 
  scale_x_reverse() + 
  guides(color = guide_legend(override.aes = list(size=8), ncol=3) )+
  theme(legend.position = 'bottom') + NoLegend()

pdf('Figure1/Figure1f_EP.pdf', width = 6, height = 6)
print(pt1f)
dev.off()

pt1e <- DimPlot(merged_combined, reduction = 'umap.pca', group.by = 'Status_new', 
                label = F, repel = TRUE) + NoAxes() + ggtitle(' ') + # + scale_color_brewer(palette="Set2") + 
  scale_color_manual(values = c('light grey', 'light grey', colorForStatus[3])) + 
  scale_x_reverse() + 
  guides(color = guide_legend(override.aes = list(size=8), ncol=3) )+
  theme(legend.position = 'bottom') + NoLegend()

pdf('Figure1/Figure1f_LP.pdf', width = 6, height = 6)
print(pt1f)
dev.off()

##-------------------------------
## Figure 1G: Heatmap of samples
## in different status
##-------------------------------
## Subset of tumor cells
Idents(merged_combined) <- merged_combined$Tumor_cells
merged_combined <- subset(merged_combined, idents = 'tumor cells')
merged_combined <- ScaleData(object = merged_combined, features = rownames(merged_combined))
Idents(merged_combined) <- merged_combined$Status_new
## For heatmap
merged_combined.maker <- FindAllMarkers(merged_combined, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
merged_combined.maker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
head(merged_combined.maker)

## Make the heatmap for marker genes
merged_combined.maker %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
# library(DoMultiBarHeatmap)
pdf(paste0('Figure1/Figure1g.pdf'), width = 8, height = 7)
print(DoHeatmap(merged_combined, features = top10$gene,
                group.colors = brewer.pal(n = 3, name = "Set2")) +
        scale_fill_viridis() +
        theme(text = element_text(size = 10),
              legend.text = element_text(size=10)) +
        guides(fill=guide_legend(title=" "))
      )
dev.off()



##-------------------------------
## Figure 1H: Hallmark pathway
##-------------------------------
source('/CDK4/source.R')
pathways_list <- c('hallmark', 'gobp')
df_hallmark <- read_pathway_tumor(i = 1)
df_hallmark <- df_hallmark[df_hallmark$pathway != 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', ]
df_hallmark$pathway <- as.factor(unlist(lapply(as.character(df_hallmark$pathway), function(x) strsplit(x, split = 'HALLMARK_')[[1]][2])))

pt1j <- ggplot(df_hallmark, aes(x = pair, y = reorder(pathway, NES))) +        ## global aes
  geom_point(aes(fill = NES, size =  -log(pval)),
             color = 'black',
             shape = 21,
             stroke = 0.01)  +  
  ggtitle('Hallmark Pathway') + 
  scale_x_discrete(labels=c('EP vs BL', 'LP vs BL', 'LP vs EP')) + 
  xlab("") + ylab("") +
  labs(size='-log(P-values)') +
  # scale_fill_gradient(high = 'dark orange', low = 'dark blue') + 
  scale_fill_gradientn(
    ## limits = c(-1.5, 1.5),
    colors = c("#5DBCFF", "#6DCCFF", "white", "#F494BE", "#F484AE"))+
  # scale_size(range = c(0, 3.2), limits = c(0, 100), breaks = c(0,20,40,60,80,100))+             ## to tune the size of circles
  theme(
    text = element_text(size = 10),
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


pdf(paste0('Figure1/Figure1h.pdf'), width = 5, height = 5)
print(pt1h)
dev.off()


